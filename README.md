# RvED
Incorporating Effect Directions Boosts Statistical Power  of Region-Based Association Tests

A preliminary version of package can be found in RvED_0.1.tar.gz. More details including data simulation, model formulation and effect direction datasets would be updated soon. 

A very useful tool (ImpG-Summary) to impute summary statistics can be found in https://github.com/huwenboshi/ImpG

Suggested:
- Step 1: Download trait-specific summary statistics from public dataset such as GWAS-catalog, GWASdb and et al.
- Step 2: Use ImpG-Summary and the reference panel (1000 Genome or genotype data of researchers) to impute summary statistics of all variants in reference panel.
- Step 3: Conduct RvED method for any genome region by incorporating variants effect directions obtained in Step 2.
