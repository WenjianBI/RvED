
# The following up function 

# Input:
# output_dir/output_prefix: should be the same as in vcf2map.sh. In output_dir, there exists files of 1) .reformat.anno, 2) .raw 3) .frq  
# z: results from software of ImpG
# pheno: R matrix of n*(m+1) where n is sample size and m is number of phenotypes. Column 1 is "IID" (the same as .raw data) and other columns are for different phenotypes. Column names are required.
# covar: R matrix of n*(k+1) where n is sample size and m is number of covariates. Column 1 is "IID" (the same as .raw data) and other columns are for different covariates. Column names are required.

map2final = function(output_dir,
                     output_prefix,
                     z,
                     pheno,
                     covar,
                     add.suffix="",
                     min.nSNPs=3,
                     min.nMACs=15)
{
  ## Specify files in output_dir
  anno.file=paste0(output_dir,"/",output_prefix,".reformat.anno")
  raw.file=paste0(output_dir,"/",output_prefix,".raw")
  frq.file=paste0(output_dir,"/",output_prefix,".frq")
  ## Read in files
  anno <- read.delim(anno.file, header=FALSE)
  colnames(anno)=c("chr","pos","REF","ALT","rsID","Gene","Func")
  raw <- read.csv(raw.file, sep="")
  frq <- read.csv(frq.file, sep="")
  ###
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
  if(any(intersect(colnames(covar),colnames(raw))!="IID")) stop("Except 'IID', no other column names can be shared by covar and raw.")
  if(any(intersect(colnames(pheno),colnames(raw))!="IID")) stop("Except 'IID', no other column names can be shared by pheno and raw.")
  merge.gap = merge(covar,raw)
  merge.gap = merge(pheno,merge.gap)
  
  Genes=unique(update.anno$Gene)
  
  outidx=F  # indicate if output exists
  for(p in colnames(pheno)[-1]){
    print(paste("Analyzing phenotype of",p,"......."))
    
    Phen=merge.gap[,p]
    Cova=as.matrix(merge.gap[,colnames(covar)[-1],drop=F])
    
    for(g in Genes){
      g.anno=subset(merge.anno,Gene==g)
      Geno=as.matrix(merge.gap[,g.anno$rawID,drop=F])
      if(fail.check(Geno,min.nSNPs,min.nMACs)) next;
      Wt=dbeta(g.anno$MAF,1,25)^2
      Nu=sign(g.anno$update.Z.score)
      Nu=ifelse(is.na(Nu),0,Nu)
      Nu.SKAT=Nu*Wt
      names(Phen)=rownames(Geno)=rownames(Cova)=merge.gap$IID
      p.SKATs=SKATs(Phen,Geno,Cova,Wt)
      p.BT.d=Burden.d(Phen,Geno,Cova,Wt,Nu)
      p.SKAT.d=SKAT.d(Phen,Geno,Cova,Wt,Nu.SKAT)$pvalue
      
      tmp=data.frame(phenotype=p,
                     gene=g,
                     SKAT=p.SKATs$p.SKAT,BT=p.SKATs$p.Burden,SKAT.O=p.SKATs$p.SKAT.O,
                     SKAT.d=p.SKAT.d,BT.d=p.BT.d)
      if(outidx==F){
        output=tmp
        outidx=T
      }else{
        output=rbind.data.frame(output,tmp)
      }
    }
  }
  output.file=paste0(output_dir,"/",output_prefix,add.suffix,".final.csv")
  write.csv(output,file=output.file,row.names = F)
}

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

# check if total number of SNPs exceed min.nSNPs and total number of minor alleles exceed min.nMACs
fail.check=function(Geno,min.nSNPs,min.nMACs){
  if(ncol(Geno)<min.nSNPs) return(TRUE)
  AC1=apply(Geno, 2, sum, na.rm=T)
  AC2=apply(2-Geno, 2, sum, na.rm=T)
  MAC=pmin(AC1,AC2)
  if(sum(MAC)<min.nMACs) return(TRUE)
  return(FALSE)
}

###### Examples:
# options(stringsAsFactors = F)
# library(RvED)
# source("map2final.R")
# output_dir="Y:/one_sided_SKAT/Update_NFBC66_2018_07_10/Results"
# output_prefix="chr21.maf.0.05.R2.0.3.dose"
# z <- read.csv("Y:/one_sided_SKAT/summary-statistics-ImpG/output/chr21/all.impz", sep="")
# Phen.all <- read.csv("Y:/one_sided_SKAT/Update_NFBC66_2018_07_10/pheno1.csv")
# pheno=Phen.all[,c("SUBJID","FS_KOL_H","FS_KOL_L","FS_TRIGL")]
# covar=Phen.all[,c("SUBJID","SEX")]
# colnames(pheno)[1]=colnames(covar)[1]="IID"
# colnames(covar)[2]="GENDER"
# min.nSNPs=3
# min.nMACs=15
# 
# map2final(output_dir,output_prefix,z,pheno,covar)
