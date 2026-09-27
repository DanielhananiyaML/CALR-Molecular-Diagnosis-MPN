# Public cohort: GSE156336

This project uses **GSE156336** as the first public transcriptomic cohort for CALR-focused analysis.

NCBI GEO describes GSE156336 as an RNA-seq cohort containing 18 essential-thrombocythemia patient samples and 4 healthy controls, including CALR-mutated, JAK2V617F, triple-negative and healthy groups.

For the current analysis, the **GEO-provided CPM matrix** is used:

`GSE156336_Feature_Counts_RNA_seq_CPM.txt.gz`

Source: NCBI GEO accession GSE156336.

## Download

The repository download script supports both the GEO-provided CPM and raw feature-count matrices.

To download the CPM matrix used by the current analysis:

```bash
python scripts/download_public_cohort.py \
  --dataset GSE156336_CPM \
  --output data/public/GSE156336/GSE156336_Feature_Counts_RNA_seq_CPM.txt.gz
```

To download the GEO-provided raw feature-count matrix instead:

```bash
python scripts/download_public_cohort.py \
  --dataset GSE156336_RAW \
  --output data/public/GSE156336/GSE156336_Raw_Feature_Counts_RNA_seq.txt.gz
```

Downloaded data are intentionally **not committed** to this repository by default.

## Reproducibility rule

Do not manually copy, edit, or transform the GEO-provided matrix before analysis.

The current CALR expression analysis begins from the GEO-provided CPM matrix and preserves the accession, source URL, filename, file size and SHA-256 checksum in the associated provenance record.

The current CPM matrix checksum is:

`725c334c6d7aa04ea6c27a39a430d4e42169806d83fbd18eb78f68b0a75ddb8a`

The analysis uses the repository manifest to reconcile GEO sample identifiers with the matrix columns and resolves CALR using its matrix identifier.
