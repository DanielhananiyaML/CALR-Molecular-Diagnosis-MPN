# CALR Molecular Diagnosis in BCR::ABL1-Negative Myeloproliferative Neoplasms

**Daniel Hananiya — Medical Laboratory Scientist | Molecular Biology Researcher | Bioinformatics Enthusiast**

A reproducible computational and literature-grounded project for characterising **CALR exon 9 insertions/deletions (indels)** in BCR::ABL1-negative myeloproliferative neoplasms (MPNs), with emphasis on Type 1, Type 2 and non-classical/type-like variants.

## Why this project matters

CALR exon 9 indels are an important molecular class in essential thrombocythemia (ET) and primary myelofibrosis (PMF), particularly among cases lacking JAK2 V617F. The two canonical variants are the 52-bp deletion (Type 1) and 5-bp insertion (Type 2); together they account for most CALR-mutated MPNs. CALR indels are heterogeneous, so fragment size alone can be insufficient for complete molecular interpretation.

This repository deliberately separates **reference facts**, **computational classification**, **public evidence**, and **clinical interpretation**. It does not present simulated data as patient evidence and is not a clinically validated diagnostic assay.

## Molecular reference standard

The project uses **CALR MANE Select transcript NM_004343.4 / ENST00000316448.10** as the primary reference. NCBI and ClinGen identify NM_004343.4 as the MANE Select transcript. A key versioning issue is documented for the historical Type 1 notation: the commonly cited legacy `NM_004343.3:c.1092_1143del52` maps to the current MANE Select `NM_004343.4:c.1099_1150del` representation. The biological event is the same 52-bp exon 9 deletion.

## Core canonical variants

| Class | Historical/common HGVS | Current MANE Select HGVS | Protein consequence | Length change |
|---|---|---|---|---:|
| Type 1 | `NM_004343.3:c.1092_1143del52` | `NM_004343.4:c.1099_1150del` | `p.Leu367fs` / commonly `p.L367fs*46` | -52 bp |
| Type 2 | `NM_004343.4:c.1154_1155insTTGTC` | same | `p.Lys385fs` / commonly `p.K385fs*47` | +5 bp |

Both are exon 9 frameshift events that alter the CALR C-terminal sequence. Published work and current diagnostic literature emphasise that non-classical indels require careful characterisation rather than classification by fragment length alone.

## Computational workflow

`Reference definition → HGVS parsing → indel length/frame calculation → canonical Type 1/Type 2 matching → type-like classification flag → protein consequence check → evidence table → diagnostic-method interpretation → report`

## Repository structure

```text
CALR-Molecular-Diagnosis-MPN/
├── README.md
├── CITATION.cff
├── LICENSE
├── analysis/
│   └── CALR_exon9_analysis.ipynb
├── scripts/
│   ├── calr_indel_analysis.py
│   ├── calr_variant_annotation.py
│   ├── calr_sequence_analysis.py
│   ├── calr_reference.py
│   ├── calr_consequence.py
│   ├── calr_genbank_reference.py
│   └── calr_reconstruct.py
├── data/
│   ├── raw/
│   ├── processed/
│   └── metadata/project_scope.yaml
├── results/
│   ├── tables/
│   └── figures/
├── figures/
├── references/README.md
├── publication/README.md
├── tests/test_calr_engine.py
└── .github/workflows/ci.yml
```

## Current computational engine

The first engine release is intentionally transparent and testable. It can:

- parse common CALR exon 9 HGVS insertion/deletion strings;
- calculate net nucleotide change and reading-frame effect;
- identify the canonical Type 1 and Type 2 patterns;
- flag non-canonical frameshifts for deeper review;
- distinguish legacy Type 1 notation from the current MANE Select representation;
- generate a structured TSV report suitable for downstream analysis.

The engine **does not infer pathogenicity from frame shift alone** and does not substitute for clinical review.

## Sequence-aware reference layer

The repository now includes a runtime reference retriever for **NCBI RefSeq NM_004343.4** and a sequence-aware consequence engine. The reference is fetched at runtime rather than silently vendored, and the retrieval script records the accession, URL and SHA-256 checksum. NCBI currently identifies NM_004343.4 as the reviewed CALR mRNA paired with NP_004334.1. citeturn0search6turn0search10

The consequence engine applies supported HGVS indels to a supplied coding sequence, verifies the net frame change, and translates the altered coding sequence. It is deliberately a **sequence-analysis engine, not a clinical interpretation engine**.

## Diagnostic interpretation layer

The project will compare:

- PCR fragment analysis / capillary electrophoresis;
- high-resolution melting (HRM);
- Sanger sequencing;
- targeted NGS.

Evidence indicates that CALR exon 9 heterogeneity creates important limitations for assays that report only fragment length. Non-standard profiles should be reflexed to sequence-level characterisation where appropriate.

## Scientific guardrails

1. Use current transcript/version identifiers wherever possible.
2. Preserve historical HGVS aliases when they explain older literature or database records.
3. Treat CALR exon 9 as an indel-characterisation problem, not merely a 52-bp-vs-5-bp screening problem.
4. Separate somatic MPN evidence from germline ClinVar classifications.
5. Do not label a variant clinically actionable solely because it is a frameshift.
6. Do not create patient-level results from simulated or illustrative data.
7. Release only after runtime tests, provenance checks and manuscript/repository audit pass.

## Key evidence

- Klampfl et al. and Nangalia et al. established somatic CALR exon 9 frameshift mutations as major drivers in JAK2-nonmutated MPNs.
- CALR mutation characterisation literature describes Type 1, Type 2 and type-like variants and stresses precise mutation characterisation.
- Diagnostic-method comparisons show that fragment analysis, HRM, Sanger and NGS have different detection characteristics.
- Recent diagnostic literature further highlights pitfalls with non-standard CALR indels and the need for sequence-level confirmation.

## Status

**Phase 1 — Molecular foundation and computational engine: completed.**

**Phase 2 — Sequence-aware analysis, evidence integration and diagnostic workflow: active.**

The current build includes a GenBank-aware runtime reference parser and a reconstruction utility for the canonical Type 1 and Type 2 indels. The reconstruction step is deliberately kept separate from CI because it depends on live NCBI retrieval; it records source and CDS checksums when run.

Current automated test status: **23 tests passed locally**.

The CI workflow uses `pytest` so both the original engine tests and the sequence-aware parser tests are executed.

The v0.3.3 build adds a public-evidence register and an educational molecular workflow figure. The project is not a clinical diagnostic pipeline and should not be released as such.

## v0.3.5 — protein-level reconstruction milestone

A new protein-validation layer has been added using the authoritative human CALR reference protein **NP_004334.1**, cross-referenced to the reviewed UniProt CALR entry **P27797**. The reference is 417 amino acids long. The local FASTA checksum is recorded in `data/metadata/protein_reference_provenance.yaml`.

The new reconstruction produces reproducible protein-level representations of the canonical exon 9 consequences:

- **Type 1:** `NM_004343.4:c.1099_1150del` → `p.Leu367ThrfsTer46 / p.L367fs*46`, reconstructed length 411 aa.
- **Type 2:** `NM_004343.4:c.1154_1155insTTGTC` → `p.Lys385AsnfsTer47 / p.K385fs*47`, reconstructed length 430 aa.

The novel C-terminal sequences reproduce the published Type 1 and Type 2 benchmark tails. This is a **protein-level validation bridge**: it validates the established protein consequences against the reviewed reference protein, but it does not replace nucleotide-level reconstruction directly from `NM_004343.4`.

Current automated test status: **23 tests passed locally**.


## v0.4.0 scientific interpretation layer

The current development branch adds:
- sequence-derived Results and Discussion documents;
- quantitative CALR protein composition table;
- an evidence claim register;
- a research-only Type 1-like / Type 2-like annotation framework;
- a repository audit documenting scientific and clinical guardrails.

The current release recommendation is **not yet v1.0.0**. The next planned extension is a public ET/PMF dataset analysis with transparent provenance.

## Public-data phase

The project now includes a reproducible public-data layer using **NCBI GEO GSE156336**, an RNA-seq ET cohort with 18 patient samples and 4 healthy controls, including CALR-mutated, JAK2V617F and triple-negative groups. The cohort is used for exploratory transcriptomic analysis, not for clinical validation or population prevalence estimation.

The dataset is not bundled. Use `scripts/download_public_cohort.py` to retrieve the GEO-provided processed count matrix and generate a SHA-256 provenance record. Quantitative analysis should only be interpreted after the matrix passes the repository's reproducibility checks.

## v0.6.1 — independent public-cohort validation framework

An independent public-cohort registry has been added using GEO **GSE103237** and **GSE54644** to complement the primary GSE156336 RNA-seq cohort. GSE103237 contains 7 CALR-mutated and 17 JAK2V617F ET samples within a broader PV/ET expression study; GSE54644 provides broader MPN/JAK-STAT pathway context. These cohorts are explicitly separated by analytical role and are not pooled automatically.

The repository currently makes **no numerical claim from these validation cohorts** because the original matrices have not yet been acquired in the current execution environment. This is intentional: accession metadata are evidence; quantitative results require the underlying files and checksum provenance.


## Public-data milestone

GSE156336 is verified against the NCBI GEO record. Quantitative CALR expression analysis is implemented and integrity-gated; biological inference is withheld until the NCBI processed matrix is acquired locally. See `publication/milestone_v0.6.1.md`.
