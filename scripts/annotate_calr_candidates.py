from pathlib import Path
import time

import pandas as pd
import requests

BASE = Path("results/biological_interpretation/tables")

INPUT = BASE / "GSE156336_CALR_candidate_genes_filtered.tsv"
OUT_ALL = BASE / "GSE156336_CALR_candidate_genes_annotated.tsv"
OUT_HIGH = BASE / "GSE156336_CALR_candidate_genes_annotated_higher.tsv"
OUT_LOW = BASE / "GSE156336_CALR_candidate_genes_annotated_lower.tsv"
OUT_SUMMARY = BASE / "GSE156336_gene_annotation_summary.txt"

# Load candidate table
df = pd.read_csv(
    INPUT,
    sep="\t",
    dtype={"gene_id": str}
)

df["gene_id"] = (
    df["gene_id"]
    .astype(str)
    .str.strip()
    .str.strip('"')
)

gene_ids = (
    df["gene_id"]
    .dropna()
    .drop_duplicates()
    .tolist()
)

print(f"CANDIDATE_GENES: {len(gene_ids)}")

# NCBI Gene E-utilities
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

annotations = []
batch_size = 100

for start in range(0, len(gene_ids), batch_size):

    batch = gene_ids[start:start + batch_size]

    print(
        f"ANNOTATING: {start + 1}-"
        f"{min(start + batch_size, len(gene_ids))}"
    )

    params = {
        "db": "gene",
        "id": ",".join(batch),
        "retmode": "json",
        "tool": "CALR_Molecular_Diagnosis_MPN",
        "email": "daniellhananiya@gmail.com",
    }

    response = requests.get(
        URL,
        params=params,
        timeout=60,
    )

    response.raise_for_status()

    result = response.json().get("result", {})

    for gene_id in batch:

        record = result.get(str(gene_id), {})

        organism = record.get("organism", {})

        annotations.append(
            {
                "gene_id": str(gene_id),
                "gene_symbol": record.get("name", ""),
                "gene_description": record.get(
                    "description", ""
                ),
                "organism": organism.get(
                    "scientificname", ""
                ),
                "chromosome": record.get(
                    "chromosome", ""
                ),
                "map_location": record.get(
                    "maplocation", ""
                ),
                "gene_type": record.get(
                    "geneticsource", ""
                ),
                "annotation_status": (
                    "annotated"
                    if record
                    else "not_found"
                ),
            }
        )

    time.sleep(0.34)

# Annotation table
annotation_df = pd.DataFrame(annotations)

# Merge with biological results
df = df.merge(
    annotation_df,
    on="gene_id",
    how="left",
    validate="one_to_one",
)

# Save complete annotated table
df.to_csv(
    OUT_ALL,
    sep="\t",
    index=False,
)

# Split by direction
higher = df[
    df["direction"] == "Higher_in_CALR"
].copy()

lower = df[
    df["direction"] == "Lower_in_CALR"
].copy()

higher.to_csv(
    OUT_HIGH,
    sep="\t",
    index=False,
)

lower.to_csv(
    OUT_LOW,
    sep="\t",
    index=False,
)

# Summary
annotated = (
    df["annotation_status"] == "annotated"
).sum()

not_found = (
    df["annotation_status"] == "not_found"
).sum()

summary = f"""
GSE156336 CALR CANDIDATE GENE ANNOTATION SUMMARY
================================================

Total candidate genes: {len(df)}

Higher in CALR: {len(higher)}

Lower in CALR: {len(lower)}

Successfully annotated: {annotated}

Not found: {not_found}

Annotation source:
NCBI Gene / Entrez E-utilities

Scientific limitation:
These candidates originate from an exploratory comparison
containing only 3 CALR samples. They are not confirmed
biomarkers or definitive differentially expressed genes
without appropriate statistical and independent validation.
"""

OUT_SUMMARY.write_text(
    summary.strip() + "\n",
    encoding="utf-8",
)

# Display results
print()
print("ANNOTATION COMPLETE")
print("===================")
print(f"TOTAL CANDIDATES: {len(df)}")
print(f"ANNOTATED: {annotated}")
print(f"NOT_FOUND: {not_found}")
print(f"HIGHER_IN_CALR: {len(higher)}")
print(f"LOWER_IN_CALR: {len(lower)}")

print()
print("TOP 20 CALR-HIGHER CANDIDATES")
print("==============================")

columns = [
    "gene_id",
    "gene_symbol",
    "gene_description",
    "CALR_mean_CPM",
    "nonCALR_ET_mean_CPM",
    "mean_difference_CPM",
    "Cohens_d",
    "CALR_detected_n",
]

print(
    higher[columns]
    .sort_values(
        "Cohens_d",
        ascending=False,
    )
    .head(20)
    .to_string(index=False)
)

print()
print("OUTPUT FILES")
print("============")
print(OUT_ALL)
print(OUT_HIGH)
print(OUT_LOW)
print(OUT_SUMMARY)
