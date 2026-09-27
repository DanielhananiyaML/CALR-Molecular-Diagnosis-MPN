# Public Cohort Analysis Plan

## Cohort

The first public-data layer uses NCBI GEO **GSE156336**, an RNA-seq cohort containing 18 essential-thrombocythemia patient samples and 4 healthy controls. The series includes CALR-mutated, JAK2V617F and triple-negative groups.

## Planned analyses

1. Verify matrix dimensions and sample identifiers.
2. Verify sample-group mapping against the GEO manifest.
3. Extract CALR expression using the matrix's actual identifier convention.
4. Perform exploratory CALR expression summaries by mutation group.
5. Where sample size permits, compare CALR-mutated ET with JAK2-mutated, triple-negative and healthy groups.
6. Report effect sizes and uncertainty rather than treating a small cohort as definitive.
7. Keep mutation classification and expression association separate from clinical diagnostic claims.

## Important limitation

GSE156336 is an ET transcriptomic cohort. It is **not** a CALR mutation prevalence cohort and cannot by itself establish population-level CALR variant frequencies, prognosis or diagnostic sensitivity.

## Extension

A second cohort may be added later for independent validation, preferably one with explicit CALR mutation subtype information or a larger MPN cohort.
