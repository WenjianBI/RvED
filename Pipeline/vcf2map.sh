#! /bin/bash

# This file is to generate a mapping file with columns of SNP_id, SNP_pos and region_id based on VCF file and annotation generated by ANNOVAR software. Here is an example to extract SNP (without indel) in exon region.

if (($#!=4))
then
echo "vcf2R.sh vcf_file anno_file output_dir output_prefix"
exit
fi

VCF=$1      # VCF file. NOTE that VCF file and annotation file should follow the same reference panel.
ANNO=$2     # annotation file (by ANNOVAR) with columns of 'Func.refGene', 'Gene.refGene', 'avsnp147'
DIR=$3      # Output directory
PREFIX=$4   # Prefix for output files

module load vcftools/0.1.15

mkdir ${DIR} -p
cd ${DIR}

# Use ANNO file to select variants and use PED file to select subjects 
cat ${ANNO} | awk -F "\t" '{if($6~/\yexonic/)print $0}' > ${PREFIX}.anno
cat ${PREFIX}.anno | awk -F "\t" '{print $1"\t"$14}' > ${PREFIX}.anno.pos 

# Use vcftools to generate numeric matrix for genotype
vcftools --vcf ${VCF} --keep ${SUBJ} --out ${PREFIX} --012 --positions ${PREFIX}.anno.pos
vcftools --vcf ${VCF} --keep ${SUBJ} --out ${PREFIX} --freq --positions ${PREFIX}.anno.pos

## Example
vcf_file=/home/wbi1/one_sided_SKAT/Update_NFBC66_2018_07_10/chr21.maf.0.05.R2.0.3.dose.vcf
anno_file=/home/wbi1/one_sided_SKAT/Update_NFBC66_2018_07_10/chr21.maf.0.05.R2.0.3.dose.anno.hg19_multianno.txt
output_dir=/home/wbi1/one_sided_SKAT/Update_NFBC66_2018_07_10
output_prefix=chr21.maf.0.05.R2.0.3.dose
# VCF=$vcf_file;ANNO=$anno_file;DIR=$output_dir;PREFIX=$output_prefix
vcf2R.sh vcf_file anno_file output_dir output_prefix
