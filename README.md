# alphafold
sandbox for secondary structure modeling
- [x] `Z_final` (circle) computed N different ways by going around circle should agree
- [x] `bpp` computation based on circle,  compared to analytical for small sequences
- [x] `Z_linear` compared to analytical for small sequences
- [x] `Z_linear` gives same answer computed N different ways
- [x] Check strands that can form >= 1 base pair (`Z_final` computed N ways agree)
- [ ] Are there conventions that simplify `l`, `l_BP`, `C_init`, `C_std`?
- [ ] Check `bpp` computation against Vienna/Mccaskill style
- [ ] Work out `Z` derivatives, compare analytical vs. numerical
- [ ] Work out `bpp` derivatives, compare analytical vs. numerical
- [ ] Add AU pairs & base pair stack motifs
- [ ] Compute MFE
- [ ] Coaxial stacks
- [ ] run on tRNAs --> calculate `F` metric
- [ ] 5' and 3' dangles
- [ ] "log-linear" loop model
- [ ] noncanonical BPs

