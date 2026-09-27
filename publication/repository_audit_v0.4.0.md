# v0.4.0 Repository Audit

## Scope

Audit covers scientific claims, reference-version consistency, computational outputs, test status, provenance, manuscript alignment and clinical-interpretation safeguards.

## Findings

- Reference transcript: `NM_004343.4` — consistently represented.
- Reference protein: `NP_004334.1` — consistently represented in the protein layer.
- Canonical Type 1 and Type 2 HGVS representations — documented with historical/current distinction.
- Protein lengths — represented in results tables.
- Charge analysis — explicitly labeled as a residue-count proxy.
- Clinical claims — restricted; no patient-level diagnostic claims.
- Non-standard CALR indels — explicitly addressed.
- Germline/somatic distinction — explicitly documented.
- Diagnostic method limitations — incorporated.
- Automated test suite — retained as the release gate.
- Figures — documented as research/portfolio outputs rather than clinical validation figures.

## Release recommendation

**Do not release v1.0.0 yet.**

Recommended next addition: a public ET/PMF dataset with transparent provenance and a small, reproducible CALR-focused cohort analysis. After that analysis, perform a second audit and decide whether the repository is ready for v1.0.0.
