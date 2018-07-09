## This folder is for simulations 

- Step 1: download 1000 Genome reference panels and only retain EUR population.
- Step 2: extract SNP alleles between chr4:77356278 to chr4:77703432 (as the example in 'sim1000G' R package).
- Step 3: remove SNP alleles without polymorphism in EUR population. After filtering, there remains 671 rare variants (RVs, MAF < 0.05) and 826 common variants (CVs, MAF > 0.05)
- Step 4: simulate genotype based on 'sim1000G' R package, select 10 RVs and randomly assign directions of 1 or -1 to generate continuous phenotype. Effect size is given based on minor allele frequency in 1000G dataset. 
- Step 5: use linear regression directly to calculate Z scores for all CVs in this region and use ImpG to impute the Z scores for all RVs in this region.
- Step 6: evaluate the direction accuracy of the 10 RVs
- Step 7: repeat the processes of 'Step 5' and 'Step 6' for 1000 times to get an average accuracy rate.

### The bash code
```
cd /home/wbi1/one_sided_SKAT/sim1000G
vcf_file=/home/wbi1/one_sided_SKAT/summary-statistics-ImpG/chr4.1kg.phase3.v5a.vcf.gz
zcat $vcf_file | awk '{if(NR<=5 || ($2>=77356278 && $2<=77703432)) print $0}' > region.bwj.vcf

cat integrated_call_samples_v3.20130502.ALL.panel | awk '$3 ~ /EUR/ {print $1}' > EUR.panel
module load vcftools
vcftools --vcf region.bwj.vcf --out region.EUR --recode --keep EUR.panel --maf 0.0001 --remove-indels
cat region.EUR.recode.vcf | java -jar vcf2beagle.jar . region.EUR
```

### The R code
```
library(sim1000G)

vcf_file="Y:/one_sided_SKAT/sim1000G/region.EUR.recode.vcf"
vcf = readVCF(vcf_file, min_maf = NA, max_maf = NA,maxNumberOfVariants = 3500)
af = apply(vcf$gt1+vcf$gt2,1,mean)/2
maf = ifelse(af<0.5,af,1-af)
region.EUR <- read.delim("Y:/one_sided_SKAT/sim1000G/region.EUR.markers", header=FALSE)
combine=cbind(region.EUR,af,maf)
colnames(combine)=c("SNP","pos","REF","ALT","ALT.AF","MAF")

### Simulate genotype data
startSimulation(vcf, totalNumberOfIndividuals = 2000)
ids = generateUnrelatedIndividuals(1000)

genotype = retrieveGenotypes(ids)

### Simulate phenotype data
pos.CVs=which(combine$MAF>0.05)
pos.RVs=which(combine$MAF<0.05)

n.causalRVs=10
c=0.2   # an efficient for effect size assignment
pos.cRVs=sample(pos.CVs,n.causalRVs)
dir.cRVs=sample(c(1,-1),n.causalRVs,replace = T)
maf.cRVs=combine$MAF[pos.cRVs]
eff.cRVs=c*abs(log10(maf.cRVs))

phenotype = colSums(t(genotype[,pos.cRVs])*dir.cRVs*eff.cRVs)

### Calculate Z-score for common variants
for(i in pos.CVs){
  basic=combine[i,]
  g=genotype[,i]
  lm.res=glm(phenotype~g)
  lm.coef=summary(lm.res)$coefficients
  z=lm.coef["g","Estimate"]/lm.coef["g","Std. Error"]
  if(i==pos.CVs[1]) Z.CVs=data.frame(basic,z)
  else Z.CVs=rbind(Z.CVs,data.frame(basic,z))
}
```

### The bash code
```
cp -R /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/backup/ImpG-master/ImpG-Bins /home/wbi1/one_sided_SKAT/sim1000G;
cd  /home/wbi1/one_sided_SKAT/sim1000G/ImpG-Bins;
module load gcc;
make;
python /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/bwj.py \
  -p /home/wbi1/one_sided_SKAT/summary-statistics-ImpG/integrated_call_samples_v3.20130502.ALL.panel \
  -o /home/wbi1/one_sided_SKAT/sim1000G/output \
  -m /home/wbi1/one_sided_SKAT/sim1000G/region.EUR.markers \
  -b /home/wbi1/one_sided_SKAT/sim1000G/region.EUR.bgl \
  -t /home/wbi1/one_sided_SKAT/sim1000G/Zscore.txt \
  --pop EUR --bin /home/wbi1/one_sided_SKAT/sim1000G/ImpG-Bins --maf 0.0001 --lambd 0.1
```
