# Public-data milestone execution — v0.6.1

## Objective
Execute the quantitative CALR expression milestone using the NCBI GEO GSE156336 processed matrix.

## Verified dataset
NCBI GEO identifies GSE156336 as an RNA-seq study with 22 samples: 18 essential-thrombocythemia patients and 4 healthy controls. The sample list contains 3 CALR, 3 JAK2, 12 triple-negative and 4 healthy-control samples.

## Execution
The NCBI record and supplementary matrix URL were verified. Runtime acquisition was attempted from the execution environment but failed because external NCBI hostname resolution was unavailable. Therefore the quantitative matrix could not be loaded locally.

## Scientific decision
No CALR expression values, fold changes, p-values, or biological conclusions are reported. The analysis remains gated until the exact NCBI matrix is downloaded and passes integrity, sample-reconciliation and CALR-identifier checks.

## Reproducible completion command
```bash
python scripts/download_public_cohort.py --dataset GSE156336_CPM --output data/public/GSE156336/GSE156336_Feature_Counts_RNA_seq_CPM.txt.gz
python scripts/analyze_public_cohort.py --matrix data/public/GSE156336/GSE156336_Feature_Counts_RNA_seq_CPM.txt.gz
```
