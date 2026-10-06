# Embedded bcftools consensus

bcftools revision: `92c31be0319d296b21fec9ddf70686c7ef06d6f3`.
The files in `upstream/` are unchanged imports with their original licenses.
`ConsensusEmbedded.c` adapts allele application and forward chains for resident
sequence construction. STAR selects donor alleles explicitly before applying them.
`ConsensusSymbols.h` makes the bcftools region-index symbols private to this adapter.
All consumers link STAR's one bundled HTSlib; no bcftools executable is invoked.
Unused standalone CLI sections are removed by the linker. Linux was validated;
other platforms and linkers need their corresponding build and output comparison.
