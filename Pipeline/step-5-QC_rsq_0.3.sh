#! /bin/bash

# how to submit it to cluster
# bsub  -R 'rusage[mem=4000]' -J wbi1[1-22] -q biostat -P 1234 '/home/wbi1/database/Gkang-dbGaP/phs000276/RootStudyConsent_phs000276.NHLBI_NFBC66.v2.p1.c1.GRU/Data_Process/codes/step-5-QC_rsq_0.3.sh $LSB_JOBINDEX'

# This file is to remove imputated variants whose R-sq less than 0.3 and maf greater than 0.05
chr=$1

# load module
module load vcflib/051916
# set working directory
cd /home/wbi1/database/Gkang-dbGaP/phs000276/RootStudyConsent_phs000276.NHLBI_NFBC66.v2.p1.c1.GRU/Data_Process

vcffilter -f "R2 > 0.3" -f "MAF < 0.05" step_4_after_imputation/chr$chr.dose.vcf.gz > step_5_QC/chr$chr.maf.0.05.R2.0.3.dose.vcf  
/home/wbi1/Pipeline_of_BWJ/ANNOTATION_hg19.sh step_6_ANNO /home/wbi1/database/Gkang-dbGaP/phs000276/RootStudyConsent_phs000276.NHLBI_NFBC66.v2.p1.c1.GRU/Data_Process/step_5_QC/chr$chr.maf.0.05.R2.0.3.dose.vcf chr$chr.maf.0.05.R2.0.3.dose 



