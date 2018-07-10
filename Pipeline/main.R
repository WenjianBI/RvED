
# bsub -R 'rusage[mem=500]' -J wbi1[1-66] -q short -P 1234 -app R-3.3.1 'R --no-save --arg %I --$LSB_JOBINDEX </home/wbi1/one_sided_SKAT/NFBC66_2017_09_18/exon/main.R'

args.orig <- commandArgs()
print(args.orig)
n.cpu = as.numeric(substr(args.orig[5], 3, 10)) #filenum is the # of CPU;

source('/home/wbi1/one_sided_SKAT/subfunc_v1_2017_06_26/subfunc-v0.8.R')
dir="/home/wbi1/one_sided_SKAT/NFBC66_2017_09_18"                                          # overall directory
DIR.d="/home/wbi1/one_sided_SKAT/BEN_2017_09_21/Directions-predication/combined-results"   # directory of directions

setwd(dir)
pheno1=read.csv("pheno1.csv")
Covars=as.matrix(pheno1[,"SEX",drop=F])
rownames(Covars)=pheno1$SUBJID

n.cpu.table=data.frame(n.cpu=1:66,
                       chr=rep(1:22,3),
                       phe=rep(8:10,each=22))

chr=n.cpu.table$chr[n.cpu.table$n.cpu==n.cpu]
phe=n.cpu.table$phe[n.cpu.table$n.cpu==n.cpu]
# The below 3 rows are for preparation of variants directions
if(phe==8) pheno="hdl"
if(phe==9) pheno="ldl"
if(phe==10) pheno="trig"

setwd(DIR.d)
direction=read.csv(paste0(pheno,"_chr",chr,".csv"))
colnames(direction)[2]="pos1"

setwd(dir)
geno.file=paste0("exon/Genotypes_R/chr",chr,".RData")
load(geno.file)
Phen=pheno1[,phe]
Phen=log(Phen)  # make the phenotype look like normal distribution
Phen.name=colnames(pheno1)[phe]
genes=unique(anno$V7)
i=genes[1]
maf$pos1=as.numeric(gsub("chr.*_","",maf$pos))
maf=merge(maf,direction,all.x=T,by="pos1")
# direction predication
maf$d=ifelse(maf$Estimate>0,1,-1)
maf$d[is.na(maf$d)]=0


DIR.out=paste0("/home/wbi1/one_sided_SKAT/NFBC66_2017_09_18/exon/Results/",Phen.name)
dir.create(DIR.out,recursive = T)
setwd(DIR.out)

pvals.output=c()
genes.output=c()
for(i in genes)
{
  print(i)
  pos=anno$pos[anno$V7==i]
  if(length(pos)>=3)
  {
    SNPs.pos=which(is.element(maf$pos,pos))  # positions in "mat" and "maf" are the same
    Geno=as.matrix(mat[,SNPs.pos])
    if(sum(c(Geno))<15) next
    d=maf$d[SNPs.pos]
    # summary of directions
    n0=sum(d==0);
    np=sum(d==1)   # number of positive
    nn=sum(d==-1)  # number of negative
    #
    pvals=SKATs(Phen,Geno,Covars)  # original methods including SKAT, Burden and SKAT-O methods
    # The below are to output 7 p-values
    Wt=Wt.func(Geno)
    Nu.SKAT.d=d*Wt
    Nu.BT.d=d
    pval.SKAT.d=SKAT.d(Phen,Geno,Covars,Wt,Nu.SKAT.d)$pvalue
    pval.BT.d=Burden.d(Phen,Geno,Covars,Wt,Nu.BT.d)
    # The below is to combine 4 p-values together as output
    pvals.tot=c(pvals$p.SKAT,
                pvals$p.Burden,
                pvals$p.SKAT.O,
                pval.SKAT.d,
                pval.BT.d,
                n0,np,nn)
    genes.output=c(genes.output,i)
    pvals.output=rbind(pvals.output,
                       pvals.tot)
  }
}
rownames(pvals.output)=genes.output
colnames(pvals.output)=c("p.SKAT","p.Burden","p.SKAT.O",
                         "pval.SKAT.d","pval.BT.d",
                         "n.0","n.pos","n.neg")
write.csv(pvals.output,file = paste0("chr",chr,".csv"))
