#!/usr/bin/env python3
from pathlib import Path
import math
import urllib.request
import json
import numpy as np
import pandas as pd

MATRIX = Path("data/public/GSE156336/GSE156336_Feature_Counts_RNA_seq_CPM.txt.gz")
MANIFEST = Path("data/public/GSE156336/manifest.tsv")
OUT = Path("results/biological_interpretation")
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
TOP_N = 100

def cohens_d(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    va = np.var(a, ddof=1)
    vb = np.var(b, ddof=1)
    df = len(a) + len(b) - 2
    if df <= 0:
        return np.nan
    pooled = math.sqrt(((len(a)-1)*va + (len(b)-1)*vb) / df)
    if pooled == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / pooled

def main():
    print("=" * 70)
    print("GSE156336 BIOLOGICAL INTERPRETATION")
    print("=" * 70)

    if not MATRIX.exists():
        raise SystemExit(f"MATRIX_NOT_FOUND: {MATRIX}")
    if not MANIFEST.exists():
        raise SystemExit(f"MANIFEST_NOT_FOUND: {MANIFEST}")

    manifest = pd.read_csv(MANIFEST, sep="\t")

    # GEO matrix: first column is gene identifier; remaining columns are samples.
    # GEO matrix has 22 sample names in the header, while each data row
    # contains one gene identifier followed by the 22 sample values.
    header_only = pd.read_csv(
        MATRIX,
        sep=r"\s+",
        engine="python",
        compression="gzip",
        nrows=0
    )
    sample_cols = [
        str(x).strip().strip('"')
        for x in header_only.columns
    ]

    df = pd.read_csv(
        MATRIX,
        sep=r"\s+",
        engine="python",
        compression="gzip",
        header=None,
        skiprows=1,
        names=["gene_id"] + sample_cols
    )

    print(f"MATRIX_ROWS: {len(df)}")
    print(f"MATRIX_COLUMNS: {len(df.columns)}")
    print(f"MANIFEST_SAMPLES: {len(manifest)}")

    gene_col = df.columns[0]

    required = manifest["matrix_column"].tolist()
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise SystemExit(
            "MATRIX_SAMPLE_MISMATCH: " + ", ".join(missing)
        )

    print("MATRIX_SAMPLE_RECONCILIATION: PASS")

    expected = {
        "CALR": 3,
        "JAK2": 3,
        "Triple_negative": 12,
        "Healthy": 4,
    }

    observed = manifest["group"].value_counts().to_dict()
    print("GROUP_COUNTS:", observed)

    for group, n in expected.items():
        if observed.get(group, 0) != n:
            raise SystemExit(
                f"GROUP_COUNT_ERROR: {group}: expected {n}, found {observed.get(group,0)}"
            )

    # Numeric expression matrix.
    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Remove duplicate/empty identifiers.
    df[gene_col] = df[gene_col].astype(str).str.strip().str.strip('"')
    df = df[df[gene_col] != ""].copy()

    # CALR-vs-non-CALR ET effect landscape.
    calr_cols = manifest.loc[
        manifest["group"] == "CALR", "matrix_column"
    ].tolist()

    noncalr_cols = manifest.loc[
        manifest["group"].isin(["JAK2", "Triple_negative"]),
        "matrix_column"
    ].tolist()

    rows = []

    for _, row in df.iterrows():
        a = row[calr_cols].dropna().to_numpy(dtype=float)
        b = row[noncalr_cols].dropna().to_numpy(dtype=float)

        if len(a) == 0 or len(b) == 0:
            continue

        calr_mean = float(np.mean(a))
        noncalr_mean = float(np.mean(b))
        log2diff = float(
            np.log2(calr_mean + 1) - np.log2(noncalr_mean + 1)
        )

        d = cohens_d(a, b)

        rows.append({
            "gene_id": row[gene_col],
            "CALR_n": len(a),
            "nonCALR_ET_n": len(b),
            "CALR_mean_CPM": calr_mean,
            "nonCALR_ET_mean_CPM": noncalr_mean,
            "mean_difference_CPM": calr_mean - noncalr_mean,
            "log2_pseudocount_ratio": log2diff,
            "Cohens_d": d,
            "abs_Cohens_d": abs(d) if not np.isnan(d) else np.nan,
            "direction": (
                "Higher_in_CALR" if log2diff > 0
                else "Lower_in_CALR"
            ),
        })

    effects = pd.DataFrame(rows)
    effects = effects.sort_values(
        "abs_Cohens_d",
        ascending=False
    )

    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    effects.to_csv(
        TABLES / "GSE156336_CALR_vs_nonCALR_ET_effects.tsv",
        sep="\t",
        index=False
    )

    higher = effects[
        effects["direction"] == "Higher_in_CALR"
    ].head(TOP_N)

    lower = effects[
        effects["direction"] == "Lower_in_CALR"
    ].sort_values(
        "log2_pseudocount_ratio"
    ).head(TOP_N)

    higher.to_csv(
        TABLES / "GSE156336_top_CALR_higher_genes.tsv",
        sep="\t",
        index=False
    )

    lower.to_csv(
        TABLES / "GSE156336_top_CALR_lower_genes.tsv",
        sep="\t",
        index=False
    )

    # PCA.
    samples = manifest["matrix_column"].tolist()
    X = np.log2(
        df[samples].to_numpy(dtype=float) + 1
    )

    keep = np.var(X, axis=1) > 0
    X = X[keep]

    X = X - X.mean(axis=1, keepdims=True)

    U, S, Vt = np.linalg.svd(
        X,
        full_matrices=False
    )

    scores = Vt[:2].T * S[:2]

    pca = manifest[
        ["sample_accession", "matrix_column", "group"]
    ].copy()

    pca["PC1"] = scores[:, 0]
    pca["PC2"] = scores[:, 1]

    explained = (S ** 2) / np.sum(S ** 2)

    pca.to_csv(
        TABLES / "GSE156336_PCA_coordinates.tsv",
        sep="\t",
        index=False
    )

    # CALR expression.
    calr_rows = df[
        df[gene_col].astype(str).isin(["811", "811.0"])
    ]

    if len(calr_rows) == 1:
        calr_row = calr_rows.iloc[0]

        calr_expression = manifest[
            ["sample_accession", "matrix_column", "group", "description"]
        ].copy()

        calr_expression["CALR_expression"] = [
            float(calr_row[x])
            for x in calr_expression["matrix_column"]
        ]

        calr_expression.to_csv(
            TABLES / "GSE156336_CALR_expression.tsv",
            sep="\t",
            index=False
        )

    # Figures.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8,6))

        for group in ["CALR", "JAK2", "Triple_negative", "Healthy"]:
            s = pca[pca["group"] == group]
            plt.scatter(
                s["PC1"],
                s["PC2"],
                s=70,
                label=group
            )

            for _, r in s.iterrows():
                plt.annotate(
                    r["matrix_column"],
                    (r["PC1"], r["PC2"]),
                    xytext=(4,4),
                    textcoords="offset points",
                    fontsize=8
                )

        plt.xlabel(f"PC1 ({explained[0]*100:.1f}% variance)")
        plt.ylabel(f"PC2 ({explained[1]*100:.1f}% variance)")
        plt.title("GSE156336 exploratory PCA")
        plt.legend()
        plt.tight_layout()
        plt.savefig(
            FIGURES / "GSE156336_PCA.png",
            dpi=300
        )
        plt.close()

        top = effects.head(20).sort_values("Cohens_d")

        plt.figure(figsize=(10,7))
        plt.barh(
            range(len(top)),
            top["Cohens_d"]
        )
        plt.yticks(
            range(len(top)),
            top["gene_id"].astype(str),
            fontsize=8
        )
        plt.xlabel("Cohen's d")
        plt.title("Top exploratory CALR-associated gene effects")
        plt.tight_layout()
        plt.savefig(
            FIGURES / "GSE156336_CALR_top_effects.png",
            dpi=300
        )
        plt.close()

        print("FIGURES: CREATED")

    except Exception as exc:
        print("FIGURES: SKIPPED:", exc)

    # Report.
    report = []
    report.append("# GSE156336 Biological Interpretation Report")
    report.append("")
    report.append("## Scope")
    report.append("")
    report.append(
        "Exploratory analysis of the GEO-provided CPM matrix from "
        "GSE156336."
    )
    report.append("")
    report.append(
        "**Important limitation:** the CALR group contains only "
        "3 samples. Gene-level effect sizes and PCA are therefore "
        "hypothesis-generating and require independent validation."
    )
    report.append("")
    report.append("## Cohort")
    report.append("")
    for group in ["CALR", "JAK2", "Triple_negative", "Healthy"]:
        report.append(
            f"- {group}: {observed.get(group, 0)} samples"
        )

    report.append("")
    report.append("## PCA")
    report.append("")
    report.append(
        f"- PC1 explained variance: {explained[0]*100:.2f}%"
    )
    report.append(
        f"- PC2 explained variance: {explained[1]*100:.2f}%"
    )

    report.append("")
    report.append("## Top exploratory CALR-vs-non-CALR effects")
    report.append("")
    report.append(
        "| Gene ID | Mean difference | Cohen's d | Direction |"
    )
    report.append(
        "|---|---:|---:|---|"
    )

    for _, r in effects.head(20).iterrows():
        report.append(
            f"| {r['gene_id']} | "
            f"{r['mean_difference_CPM']:.4f} | "
            f"{r['Cohens_d']:.4f} | "
            f"{r['direction']} |"
        )

    report.append("")
    report.append("## Scientific interpretation")
    report.append("")
    report.append(
        "This analysis identifies expression features that differ "
        "between the three CALR-mutated samples and the combined "
        "JAK2/triple-negative ET samples."
    )
    report.append("")
    report.append(
        "Because CALR n=3, these results should not be presented as "
        "definitive differential-expression findings. They are useful "
        "for generating biological hypotheses and selecting candidates "
        "for independent validation."
    )
    report.append("")
    report.append("## Reproducibility")
    report.append("")
    report.append(
        f"- Input: `{MATRIX}`"
    )
    report.append(
        "- Matrix type: GEO-provided CPM"
    )
    report.append(
        "- CALR matrix identifier: 811"
    )
    report.append(
        f"- Matrix rows analysed: {len(df)}"
    )
    report.append(
        f"- Samples analysed: {len(manifest)}"
    )

    (OUT / "BIOLOGICAL_INTERPRETATION_REPORT.md").write_text(
        "\n".join(report) + "\n"
    )

    print("")
    print("=" * 70)
    print("BIOLOGICAL INTERPRETATION COMPLETE")
    print("=" * 70)
    print("REPORT:", OUT / "BIOLOGICAL_INTERPRETATION_REPORT.md")
    print("EFFECTS:", TABLES / "GSE156336_CALR_vs_nonCALR_ET_effects.tsv")
    print("PCA:", TABLES / "GSE156336_PCA_coordinates.tsv")
    print("")
    print("TOP 15 CALR-vs-non-CALR ET EFFECTS")
    print(
        effects[
            [
                "gene_id",
                "mean_difference_CPM",
                "Cohens_d",
                "direction"
            ]
        ].head(15).to_string(index=False)
    )

if __name__ == "__main__":
    main()
