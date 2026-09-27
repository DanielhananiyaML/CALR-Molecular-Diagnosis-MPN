# Public Cohort Results — Analysis-Ready

## Dataset selection

GSE156336 was selected as the first public transcriptomic cohort because NCBI GEO identifies it as an RNA-seq study of 18 essential-thrombocythemia patients and 4 healthy controls, including CALR-mutated, JAK2V617F and triple-negative groups.

## Current status

The repository contains the cohort manifest, provenance metadata and reproducible download/analysis scripts. The processed count matrix is **not bundled** because the current execution environment cannot resolve the NCBI FTP host. Therefore, no expression values are reported here and no biological inference is made from an unverified local copy.

## Reproducibility gate

Before quantitative interpretation, the downloaded matrix must pass:

- accession verification;
- gzip integrity check;
- SHA-256 capture;
- sample-column/manifest reconciliation;
- CALR identifier discovery;
- missing-value and duplicate-sample checks.

Only after these checks pass should the CALR expression table and downstream statistics be generated.
