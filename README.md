# RvED
Incorporating Effect Directions Boosts Statistical Power  of Region-Based Association Tests

## R Package RvED_0.1.tar.gz
A preliminary version of package can be found in RvED_0.1.tar.gz. More details including data simulation, model formulation and effect direction datasets would be updated soon. 

NOTE: please install R package of kinship before installing R package of RvED.

## Some bioinformatics tools 
- A tool of ImpG-Summary to impute summary statistics can be found in https://github.com/huwenboshi/ImpG
- A tool of vcf2beagle.jar to transform VCF file to Beagle file can be found in  https://faculty.washington.edu/browning/beagle_utilities/utilities.html#vcf2beagle

## Suggestion about variants effect directions assignment:
Apply ImpG-Summary to impute summary statistics
### Pre-imputation 
- Download trait-specific summary statistics from public dataset such as GWAS-catalog (https://www.ebi.ac.uk/gwas/downloads/summary-statistics).
- Download VCF files of 1000 Genome in http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/
- Use vcf2beagle.jar to generate *bgl.gz* file and *markers* file for 1000 Genome data.
- Download sample panels of 1000 Genome in http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/?C=D;O=A
### NOTE about imputation
- All files of .gz should be extracted first.

