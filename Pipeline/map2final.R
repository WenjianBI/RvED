
rm(list=ls())
library(RvED)
options(stringsAsFactors = F)
anno <- read.delim("Y:/one_sided_SKAT/Update_NFBC66_2018_07_10/Results/chr21.maf.0.05.R2.0.3.dose.reformat.anno", header=FALSE)
colnames(anno)=c("chr","pos","REF","ALT","rsID","Gene","Func")
raw <- read.csv("Y:/one_sided_SKAT/Update_NFBC66_2018_07_10/Results/chr21.maf.0.05.R2.0.3.dose.raw", sep="")
raw = raw[,c(-1,-3,-4,-5,-6)]  # only retain IID and genotype
frq <- read.csv("Y:/one_sided_SKAT/Update_NFBC66_2018_07_10/Results/chr21.maf.0.05.R2.0.3.dose.frq", sep="")
z <- read.csv("Y:/one_sided_SKAT/summary-statistics-ImpG/output/chr21/all.impz", sep="")
pheno1 <- read.csv("Y:/one_sided_SKAT/Update_NFBC66_2018_07_10/pheno1.csv")

colnames(pheno1)[2]="IID"

## NOTE: SNPs in raw and frq files MUST exist in anno file

# some SNPs may are assigned to multiple regions (e.g. genes)
post.anno = function(anno){
  update.anno=c()
  for(i in 1:nrow(anno)){
    genes=unlist(strsplit(anno$Gene[i],split=",|;"))
    tmp=cbind(CHR=anno$chr[i],
              pos=anno$pos[i],
              SNP_id=anno$rsID[i],
              Gene=genes)
    update.anno=rbind(update.anno,tmp)
  }
  update.anno=as.data.frame(update.anno)
  update.anno$SNP=paste0(update.anno$CHR,":",update.anno$pos)
  return(update.anno)
}

update.anno=post.anno(anno)
merge.anno=merge(update.anno,frq,all = T)
merge.anno=merge(merge.anno,z,all.x = T)
# case 1: same ref and alt alleles: keep the orignal z score
# case 2: reverse ref and alt alleles: opposite of z score
# case 3: NA for other cases
merge.anno$update.Z.score=with(merge.anno,ifelse(Ref_allele==A2&Alt_Allele==A1,
                                                 Z.score,
                                                 ifelse(Ref_allele==A1&Alt_Allele==A2,
                                                        -1*Z.score,
                                                        NA)))
merge.anno$rawID=with(merge.anno,paste0("X",CHR,".",pos,"_",A1))
merge.gap = merge(pheno1,raw)


Genes=unique(update.anno$Gene)
g=Genes[1]
g.anno=subset(merge.anno,Gene==g)

pheno.name="FS_KOL_H"
Phen=merge.gap[,pheno.name]
Geno=as.matrix(merge.gap[,g.anno$rawID])
Cova=as.matrix(merge.gap$SEX,nrow(merge.gap),1)
Wt=dbeta(g.anno$MAF,1,25)^2
Nu=sign(g.anno$update.Z.score)
Nu=ifelse(is.na(Nu),0,Nu)
Nu.SKAT=Nu*Wt

names(Phen)=rownames(Geno)=rownames(Cova)=merge.gap$IID

SKATs(Phen,Geno,Cova,Wt)
Burden.d(Phen,Geno,Cova,Wt,Nu)
SKAT.d(Phen,Geno,Cova,Wt,Nu.SKAT)

subset(merge.gap[,c("FS_KOL_H","BMI","SBP","DBP")],is.na(BMI))
