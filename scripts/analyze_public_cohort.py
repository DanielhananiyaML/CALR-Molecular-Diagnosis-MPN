#!/usr/bin/env python3
"""Reproducible CALR expression analysis for GSE156336.

Accepts the NCBI GEO processed CPM matrix or raw feature-count matrix.
The script resolves CALR using its Entrez Gene ID (811), reconciles
GEO matrix columns against the repository manifest, checks for duplicate
samples and missing CALR values, and produces reproducible group summaries.
"""
from __future__ import annotations

import argparse
import gzip
import pathlib

import pandas as pd


CALR_IDS = {"CALR", "ENSG00000179218", "811"}


def resolve_calr(index):
    """Return matrix index entries corresponding to CALR."""
    hits = []

    for value in map(str, index):
        normalized = value.strip().strip('"').strip("'")
        base = normalized.split(".")[0].upper()

        if (
            base in CALR_IDS
            or base.startswith("CALR|")
            or base.endswith("|CALR")
        ):
            hits.append(value)

    return hits


def read_matrix(path):
    """Read a GEO whitespace-delimited expression matrix."""
    opener = gzip.open if str(path).endswith(".gz") else open

    with opener(path, "rt") as fh:
        matrix = pd.read_csv(
            fh,
            sep=r"\s+",
            engine="python",
            index_col=0,
        )

    matrix.index = (
        matrix.index.astype(str)
        .str.strip()
        .str.strip('"')
        .str.strip("'")
    )

    matrix.columns = (
        matrix.columns.astype(str)
        .str.strip()
        .str.strip('"')
        .str.strip("'")
    )

    return matrix


def mean_difference_ci(x, y, confidence=0.95):
    """Mean difference (x - y) with a Welch 95% CI.

    Uses the Welch-Satterthwaite degrees of freedom and Student's
    t critical value. With very small groups this remains explicitly
    exploratory and should not be interpreted as definitive.
    """
    import math

    x = pd.Series(x, dtype="float64").dropna()
    y = pd.Series(y, dtype="float64").dropna()

    difference = float(x.mean() - y.mean())

    nx = len(x)
    ny = len(y)

    if nx < 2 or ny < 2:
        return difference, float("nan"), float("nan")

    vx = float(x.var(ddof=1))
    vy = float(y.var(ddof=1))

    standard_error = math.sqrt((vx / nx) + (vy / ny))

    if standard_error == 0:
        return difference, difference, difference

    # Welch-Satterthwaite degrees of freedom.
    variance_x = vx / nx
    variance_y = vy / ny

    numerator = (variance_x + variance_y) ** 2
    denominator = (
        (variance_x ** 2) / (nx - 1)
        + (variance_y ** 2) / (ny - 1)
    )

    degrees_of_freedom = numerator / denominator

    # Student's t critical value for a two-sided confidence interval.
    from scipy.stats import t

    alpha = 1.0 - confidence
    critical = float(
        t.ppf(1.0 - alpha / 2.0, degrees_of_freedom)
    )

    lower = difference - critical * standard_error
    upper = difference + critical * standard_error

    return difference, lower, upper


def cohens_d(x, y):
    """Cohen's d for two independent groups."""
    import math

    x = pd.Series(x, dtype="float64").dropna()
    y = pd.Series(y, dtype="float64").dropna()

    nx = len(x)
    ny = len(y)

    if nx < 2 or ny < 2:
        return float("nan")

    vx = float(x.var(ddof=1))
    vy = float(y.var(ddof=1))

    denominator = nx + ny - 2

    if denominator <= 0:
        return float("nan")

    pooled_variance = (
        ((nx - 1) * vx) +
        ((ny - 1) * vy)
    ) / denominator

    pooled_sd = math.sqrt(pooled_variance)

    if pooled_sd == 0:
        return float("nan")

    return float((x.mean() - y.mean()) / pooled_sd)


def exploratory_comparisons(result, reference_group="CALR"):
    """Calculate exploratory CALR-group effect sizes and uncertainty."""
    comparisons = []

    x = result.loc[
        result["group"] == reference_group,
        "CALR_expression",
    ]

    for group in sorted(result["group"].unique()):
        if group == reference_group:
            continue

        y = result.loc[
            result["group"] == group,
            "CALR_expression",
        ]

        difference, ci_low, ci_high = mean_difference_ci(x, y)

        comparisons.append(
            {
                "reference_group": reference_group,
                "comparison_group": group,
                "reference_n": len(x),
                "comparison_n": len(y),
                "reference_mean": float(x.mean()),
                "comparison_mean": float(y.mean()),
                "mean_difference": difference,
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "cohens_d": cohens_d(x, y),
            }
        )

    return pd.DataFrame(comparisons)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument(
        "--manifest",
        default="data/public/GSE156336/manifest.tsv",
    )
    ap.add_argument(
        "--out",
        default="results/tables/GSE156336_CALR_expression.tsv",
    )
    args = ap.parse_args()

    # ------------------------------------------------------------
    # 1. Read manifest
    # ------------------------------------------------------------
    manifest = pd.read_csv(args.manifest, sep="\t")

    required_columns = {
        "sample_accession",
        "matrix_column",
        "group",
        "description",
    }

    missing_columns = required_columns - set(manifest.columns)

    if missing_columns:
        raise SystemExit(
            "MANIFEST_SCHEMA_FAILED: missing columns "
            f"{sorted(missing_columns)}"
        )

    # ------------------------------------------------------------
    # 2. Check duplicate sample accessions / matrix columns
    # ------------------------------------------------------------
    if manifest["sample_accession"].duplicated().any():
        duplicates = manifest.loc[
            manifest["sample_accession"].duplicated(),
            "sample_accession",
        ].tolist()

        raise SystemExit(
            "DUPLICATE_SAMPLE_ACCESSIONS: "
            f"{duplicates}"
        )

    if manifest["matrix_column"].duplicated().any():
        duplicates = manifest.loc[
            manifest["matrix_column"].duplicated(),
            "matrix_column",
        ].tolist()

        raise SystemExit(
            "DUPLICATE_MATRIX_COLUMNS: "
            f"{duplicates}"
        )

    # ------------------------------------------------------------
    # 3. Read and normalize GEO matrix
    # ------------------------------------------------------------
    matrix = read_matrix(args.matrix)

    print(f"MATRIX_ROWS: {matrix.shape[0]}")
    print(f"MATRIX_SAMPLES: {matrix.shape[1]}")

    # ------------------------------------------------------------
    # 4. Verify every manifest sample exists in matrix
    # ------------------------------------------------------------
    missing = [
        column
        for column in manifest["matrix_column"]
        if column not in matrix.columns
    ]

    if missing:
        raise SystemExit(
            "SAMPLE_RECONCILIATION_FAILED: "
            f"missing {len(missing)} matrix columns: {missing}"
        )

    # ------------------------------------------------------------
    # 5. Verify matrix sample count
    # ------------------------------------------------------------
    if len(manifest) != len(matrix.columns):
        raise SystemExit(
            "SAMPLE_COUNT_MISMATCH: "
            f"manifest={len(manifest)}, "
            f"matrix={len(matrix.columns)}"
        )

    # ------------------------------------------------------------
    # 6. Resolve CALR
    # ------------------------------------------------------------
    hits = resolve_calr(matrix.index)

    if len(hits) != 1:
        raise SystemExit(
            "CALR_RESOLUTION_FAILED: "
            f"expected exactly one CALR row, found {hits}"
        )

    calr_index = hits[0]

    print(f"CALR_MATRIX_ID: {calr_index}")

    # ------------------------------------------------------------
    # 7. Extract CALR expression using explicit matrix mapping
    # ------------------------------------------------------------
    expression = pd.to_numeric(
        matrix.loc[
            calr_index,
            manifest["matrix_column"],
        ],
        errors="coerce",
    )

    if expression.isna().any():
        missing_samples = manifest.loc[
            expression.isna(),
            "sample_accession",
        ].tolist()

        raise SystemExit(
            "MISSING_EXPRESSION_VALUES: "
            f"CALR contains missing/non-numeric values for "
            f"{missing_samples}"
        )

    # ------------------------------------------------------------
    # 8. Build reproducible expression table
    # ------------------------------------------------------------
    result = manifest.copy()
    result["CALR_expression"] = expression.to_numpy()

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    result.to_csv(
        out_path,
        sep="\t",
        index=False,
    )

    # ------------------------------------------------------------
    # 9. Group-level exploratory summary
    # ------------------------------------------------------------
    summary = (
        result
        .groupby("group")["CALR_expression"]
        .agg(
            count="count",
            mean="mean",
            median="median",
            std="std",
        )
    )

    summary_path = out_path.with_name(
        "GSE156336_CALR_group_summary.tsv"
    )

    summary.to_csv(
        summary_path,
        sep="\t",
    )

    # ------------------------------------------------------------
    # 10. Exploratory effect sizes and uncertainty
    # ------------------------------------------------------------
    comparisons = exploratory_comparisons(result)

    comparison_path = out_path.with_name(
        "GSE156336_CALR_exploratory_comparisons.tsv"
    )

    comparisons.to_csv(
        comparison_path,
        sep="\t",
        index=False,
    )

    print("\nCALR GROUP SUMMARY")
    print("==================")
    print(summary.to_string())

    print("\nOUTPUTS")
    print("=======")
    print(f"Expression table: {out_path}")
    print(f"Group summary:    {summary_path}")
    print(f"Comparisons:      {comparison_path}")

    print("\nEXPLORATORY EFFECT SIZES")
    print("========================")
    print(comparisons.to_string(index=False))


if __name__ == "__main__":
    main()
