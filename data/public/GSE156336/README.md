# Public cohort: GSE156336

This project uses **GSE156336** as the first public transcriptomic cohort for CALR-focused analysis.

NCBI GEO describes GSE156336 as RNA-seq of 18 essential-thrombocythemia patient samples and 4 healthy controls, including CALR-mutated, JAK2V617F, triple-negative and healthy groups. Processed raw-count and CPM files are provided by GEO.

Source: NCBI GEO accession GSE156336.

## Download

Run:

```bash
python scripts/download_public_cohort.py --accession GSE156336
```

The script downloads the GEO processed raw-count file and records the retrieval metadata. The downloaded data are intentionally **not committed** to this repository by default.

## Reproducibility rule

Do not manually copy or alter the downloaded matrix. The analysis should begin from the GEO-provided file and preserve the accession, URL and SHA-256 checksum in the generated provenance file.
