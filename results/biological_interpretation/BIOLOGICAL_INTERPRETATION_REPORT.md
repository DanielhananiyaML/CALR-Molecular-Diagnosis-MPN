# GSE156336 Biological Interpretation Report

## Scope

Exploratory analysis of the GEO-provided CPM matrix from GSE156336.

**Important limitation:** the CALR group contains only 3 samples. Gene-level effect sizes and PCA are therefore hypothesis-generating and require independent validation.

## Cohort

- CALR: 3 samples
- JAK2: 3 samples
- Triple_negative: 12 samples
- Healthy: 4 samples

## PCA

- PC1 explained variance: 31.51%
- PC2 explained variance: 14.94%

## Top exploratory CALR-vs-non-CALR effects

| Gene ID | Mean difference | Cohen's d | Direction |
|---|---:|---:|---|
| 730005 | 0.1051 | 7.8708 | Higher_in_CALR |
| 7365 | 0.7090 | 7.6402 | Higher_in_CALR |
| 7364 | 0.7907 | 6.3848 | Higher_in_CALR |
| 83715 | 0.8798 | 4.2348 | Higher_in_CALR |
| 100423020 | 0.1446 | 3.3645 | Higher_in_CALR |
| 648691 | 0.0853 | 3.2662 | Higher_in_CALR |
| 3050 | 0.2417 | 3.2432 | Higher_in_CALR |
| 391766 | 0.0743 | 3.1994 | Higher_in_CALR |
| 100506405 | 0.0629 | 3.0356 | Higher_in_CALR |
| 728461 | 0.0288 | 3.0094 | Higher_in_CALR |
| 100422882 | 0.0288 | 3.0094 | Higher_in_CALR |
| 28510 | 0.0288 | 3.0094 | Higher_in_CALR |
| 406996 | 0.0288 | 3.0094 | Higher_in_CALR |
| 441009 | 0.2100 | 3.0088 | Higher_in_CALR |
| 606551 | 0.5365 | 3.0039 | Higher_in_CALR |
| 26189 | 0.0317 | 2.8752 | Higher_in_CALR |
| 767565 | 0.0317 | 2.8752 | Higher_in_CALR |
| 794 | 0.0971 | 2.8729 | Higher_in_CALR |
| 4082 | 281.1568 | 2.8615 | Higher_in_CALR |
| 654780 | 0.2164 | 2.8327 | Higher_in_CALR |

## Scientific interpretation

This analysis identifies expression features that differ between the three CALR-mutated samples and the combined JAK2/triple-negative ET samples.

Because CALR n=3, these results should not be presented as definitive differential-expression findings. They are useful for generating biological hypotheses and selecting candidates for independent validation.

## Reproducibility

- Input: `data/public/GSE156336/GSE156336_Feature_Counts_RNA_seq_CPM.txt.gz`
- Matrix type: GEO-provided CPM
- CALR matrix identifier: 811
- Matrix rows analysed: 25702
- Samples analysed: 22
