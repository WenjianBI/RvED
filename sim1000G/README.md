## This folder is for simulations 

- Step1: download 1000 Genome reference panels and only retain EUR population.
- Step2: extract alleles between chr4:77356278 to chr4:77703432 (as the example in 'sim1000G' R package).

### The bash code
```
cd /home/wbi1/one_sided_SKAT/sim1000G
vcf_file=/home/wbi1/one_sided_SKAT/summary-statistics-ImpG/chr4.1kg.phase3.v5a.vcf.gz
zcat $vcf_file | awk '{if(NR<=5 || ($2>=77356278 && $2<=77703432)) print $0}' > region.bwj.vcf

cat integrated_call_samples_v3.20130502.ALL.panel | awk '$3 ~ /EUR/ {print $1}' > EUR.panel
module load vcftools
vcftools --vcf region.bwj.vcf --out region.EUR --recode --keep EUR.panel --maf 0.0001
cat region.EUR.recode.vcf | java -jar vcf2beagle.jar . region.EUR
```
