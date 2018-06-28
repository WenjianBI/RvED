# RvED
Incorporating Effect Directions Boosts Statistical Power  of Region-Based Association Tests

A preliminary version of package can be found in RvED_0.1.tar.gz. More details including data simulation, model formulation and effect direction datasets would be updated soon. 

A very useful tool (ImpG-Summary) to impute summary statistics can be found in https://github.com/huwenboshi/ImpG

## Suggested:
### Step 1: 
- Download trait-specific summary statistics from public dataset such as GWAS-catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics).
- Download VCF files of 1000 Genome in http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/
- Download sample panels of 1000 Genome in http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/?C=D;O=A

### Step 2: 
Use ImpG-Summary to impute summary statistics.
- Generate *Haplotype file* and *SNP mapping file* from downloaded *1000 Genome files*.

### Step 3: 
Conduct RvED method for any genome region by incorporating variants effect directions obtained in Step 2.


