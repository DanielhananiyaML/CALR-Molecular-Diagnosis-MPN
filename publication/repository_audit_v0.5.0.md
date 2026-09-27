# Repository Audit — v0.5.0 candidate

## Scientific scope

The repository now contains a public-data analysis layer based on NCBI GEO GSE156336. The scope remains exploratory computational research and does not constitute a clinical diagnostic validation study.

## Verified

- Canonical CALR Type 1 and Type 2 protein reconstruction remains intact.
- Protein-level reconstruction is explicitly distinguished from nucleotide-level live reconstruction.
- 20 existing scientific tests pass before the public-cohort test layer.
- GSE156336 is documented as 22 samples: 18 ET patients and 4 healthy controls.
- Sample manifest contains 3 CALR, 3 JAK2, 12 triple-negative and 4 healthy samples.
- Download script records source URL and SHA-256 checksum.
- Quantitative results are withheld until the actual GEO matrix is retrieved and validated.

## Gate before release

The project should not claim cohort-derived expression findings until the matrix is downloaded, checksum-recorded, parsed and analyzed successfully.

## Recommended release

After successful cohort execution, package the public-data layer as **v0.5.0**. Do not label the full project v1.0.0 until the manuscript, provenance, independent validation and final audit are complete.
