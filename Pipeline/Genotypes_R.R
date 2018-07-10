
rm(list=ls())

# Further transformation to faciliate the usage of R codes
# First, run the function of Geno.read
# Then, run the following codes

Dir.G="Y:/one_sided_SKAT/NFBC66_2017_09_18/exon/Genotypes"
Dir.R="Y:/one_sided_SKAT/NFBC66_2017_09_18/exon/Genotypes_R"
dir.create(Dir.R)

for(chr in 1:22)
{
  print(chr)
  prefix=paste0("chr",chr)
  Geno.read(Dir.G,Dir.R,prefix)
}


Geno.read=function(Dir.G,
                   Dir.R,
                   prefix)
{
  setwd(Dir.G)
  # Read in genotype matrix and assign rownames/colnames based on SNPs and Subjects
  mat=read.table(paste0(prefix,".012"),row.names = 1)
  ind=read.table(paste0(prefix,".012.indv"))
  pos=read.table(paste0(prefix,".012.pos"))
  rownames(mat)=ind$V1
  colnames(mat)=paste0("chr",pos$V1,"_",pos$V2)
  # Read in annotation file and transform data
  anno.tmp=read.delim(paste0(prefix,".anno"), header=FALSE)
  anno.tmp$pos=paste0("chr",anno.tmp$V1,"_",anno.tmp$V14)
  which.dup=grep(",|;",anno.tmp$V7)
  which.uniq=setdiff(1:nrow(anno.tmp),which.dup)
  
  anno=anno.tmp[which.uniq,c(7,17)]
  for(i in which.dup)
  {
    genes=unlist(strsplit(anno.tmp[i,7],split=",|;"))
    pos=anno.tmp[i,17]
    for(j in genes)
    {
      anno=rbind(anno,c(j,pos))
    }
  }
  # Read in maf data
  freq=read.table(paste0(prefix,".frq"),skip = 1)
  freq$pos=paste0("chr",freq$V1,"_",freq$V2)
  freq$maf=as.numeric(gsub(".*:","",freq$V6))
  freq$maf=ifelse(freq$maf<0.05,freq$maf,1-freq$maf)
  maf=freq[,c(7,8)]
  # 
  output=paste0(Dir.R,"/",prefix,".RData")
  save(mat,anno,maf,file=output)
}

