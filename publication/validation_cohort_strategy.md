# Public-cohort validation strategy

## Rationale

The primary RNA-seq cohort, GSE156336, is small and is therefore best treated as an exploratory dataset. Independent public cohorts provide a way to test whether CALR-associated transcriptional patterns are directionally consistent across sample types and technologies.

## Cohorts

1. **GSE156336** — PBMC RNA-seq; 18 ET patients and 4 healthy controls. The GEO record identifies three CALR, three JAK2 and four healthy samples among the explicitly highlighted subset used in later validation work.
2. **GSE103237** — expression profiling in PV and ET; 7 CALR-mutated ET and 17 JAK2V617F ET samples, with normal-donor expression data represented through the related normal-donor series.
3. **GSE54644** — broad MPN expression cohort used only for pathway/context comparison.

## Planned analysis

- harmonize gene identifiers to HGNC symbols where possible;
- extract CALR expression and a predefined JAK-STAT/MAPK/TNF-NFκB pathway panel;
- compare CALR-mutated ET with healthy controls where the design permits;
- compare CALR-mutated ET with JAK2-mutated ET in cohorts where both groups are present;
- report effect sizes and confidence intervals alongside p-values;
- avoid pooling incompatible platforms without explicit batch/platform modelling;
- treat small groups as exploratory and avoid clinical prediction claims.

## Current status

The registry and analysis plan are complete. Quantitative validation is **pending acquisition of the original GEO matrices in an environment with external-data access**. No numerical validation result is asserted by this repository version.
