# Raw reference data policy

The CALR reference is retrieved from NCBI RefSeq at runtime rather than silently redistributed in the repository.

Primary reference:
- `NM_004343.4` — reviewed CALR mRNA
- paired protein: `NP_004334.1`
- source: NCBI RefSeq

Runtime retrieval scripts:
- `scripts/calr_reference.py` — FASTA retrieval
- `scripts/calr_genbank_reference.py` — GenBank retrieval with feature/CDS parsing

The GenBank workflow records a SHA-256 checksum for both the raw record and the extracted CDS. This makes future reruns auditable if the external record changes. NCBI currently lists NM_004343.4 as the reviewed CALR RefSeq mRNA and NP_004334.1 as its paired protein. citeturn0search0turn0search1
