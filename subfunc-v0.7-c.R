lib.loc="/home/wbi1/R/x86_64-pc-linux-gnu-library/3.3" 
library(CompQuadForm,lib.loc = lib.loc)
library(kinship,lib.loc = lib.loc)
library(SKAT,lib.loc = lib.loc)
# module load gcc/4.8.1

# original SKAT package function, could be used to compute burden-w, SKAT and SKAT-O results
SKATs=function(Phen,       # vector of n samples with name of subjects.
               Geno,       # matrix of n*m, where n is sample size and m is SNPs number. rownames(Geno) should be subjects. Missing values should be assigned as NA.
               Cova=NULL,  # matrix of n*k, where n is sample size and k is covariates number. rownames(Cova) should be subjects.
               Wt=NULL,         # Weight vector of m for m SNPs. Should be the same order as nrow(Geno)
               trait="c")  # "c" or "b" for "continuous" or "binary".
{
  if(is.null(Wt)) Wt=Wt.func(Geno)
  if(trait=="c") obj=SKAT_Null_Model(Phen~Cova,out_type = "C",Adjustment = FALSE)
  if(trait=="b") obj=SKAT_Null_Model(Phen~Cova,out_type = "D",Adjustment = FALSE)
  p.SKAT=SKAT(Geno,obj,r.corr=0,weights=sqrt(Wt))$p.value   # the weigting is actually the square root of weights in (Wu, et al). 
  p.Burden=SKAT(Geno,obj,r.corr=1,weights=sqrt(Wt))$p.value
  p.SKAT.O=SKAT(Geno,obj,method="SKATO",weights=sqrt(Wt))$p.value
  return(list(p.SKAT=p.SKAT,
              p.Burden=p.Burden,
              p.SKAT.O=p.SKAT.O))
}


Burden.d=function(Phen,       # vector of n samples with name of subjects.
                  Geno,       # matrix of n*m, where n is sample size and m is SNPs number. rownames(Geno) should be subjects. Missing values should be assigned as NA.
                  Cova=NULL,  # matrix of n*k, where n is sample size and k is covariates number. rownames(Cova) should be subjects.
                  Wt=NULL,    # Weight vector of m for m SNPs. Should be the same order as nrow(Geno)
                  Nu=NULL,    # Direction vector of m for m SNPs. Only 1, -1 or 0(equivalent to the removal of the corresponding SNP) is accepted.
                  trait="c",  # "c" or "b" for "continuous" or "binary".
                  one.tail=T) # if one tailed or not
{
  if(is.null(Wt)) Wt=Wt.func(Geno)
  if(is.null(Nu))
  {
    warning("There is no assignment of Nu, we still use tradtional burden test")
    Nu=rep(1,length(Wt));
    one.tail=F
  }
  Geno.b=apply(t(Geno)*Wt*Nu,2,sum)    # collapse multiple rare variants into one single variable
  if(trait=="c") res=glm(Phen~Geno.b+Cova)
  if(trait=="b") res=glm(Phen~Geno.b+Cova,family = binomial)
  pval=summary(res)$coefficients[2,4]
  if(one.tail)
  {
    coef=summary(res)$coefficients[2,1]
    pval=ifelse(coef>0,pval/2,1-pval/2)
  }
  return(pval)
}



# compute p-value based on SKAT-d method. 
SKAT.d=function(Phen,       # vector of n samples with name of subjects.
                Geno,       # matrix of n*m, where n is sample size and m is SNPs number. rownames(Geno) should be subjects. Missing values should be assigned as NA.
                Cova=NULL,  # matrix of n*k, where n is sample size and k is covariates number. rownames(Cova) should be subjects.
                Wt=NULL,         # vector of m for m SNPs. Should be the same order as nrow(Geno)
                Nu=NULL,         # vector of m for m SNPs. Should be the same order as Geno
                trait="c",  # "c" or "b" for "continuous" or "binary". "fam.c" for family-based "continuous" trait) 
                kin=NULL,   # if trait = "fam.c", one of kin, ped and ibd should be specified. kin is output of makekinship or bdsmatrix.ibd.
                ped=NULL,   # a dataframe with four columns named with c("famid","id","father.id","mother.id"). Used as input of makekinship.
                ibd=NULL,   # a dataframe with three columns named with c("id1","id2","x"). Used as input of bdsmatrix.ibd. Should contain 
                acc=1e-7,
                lim=1e5)
{
  if(is.null(Wt)) Wt=Wt.func(Geno)
  if(is.null(Nu))
  {
    warning("There is no assignment of Nu, we use default setting of Wt.")
    Nu=Wt
  }
  check.res=Input.check(Phen,Geno,Cova,Wt,Nu,trait,kin,ped,ibd)
  Phen=check.res$Phen
  Geno=check.res$Geno
  X=check.res$X   # updated covariate matriax
  Wt=check.res$Wt
  Nu=check.res$Nu
  kin=check.res$kin
  
  mu=Nu*1/Wt  # W=diag(Wt); mu=solve(W)%*%Nu
  res=obs.func(Phen,X,trait,Geno,Wt,mu,kin)
  Stats=res$Stats     # observable statistics
  Stats.SKAT=res$Stats.SKAT
  Sigma=res$Sigma     # covariance matrix of t(G)
  mu1=mu*sqrt(Wt)     # sqrt(W)%*%mu
  
  null.pars=null.para(Sigma,mu1)
  lambda=null.pars$lambda          # parameter "lambda" for weighted sum of chi-square variables (null hypothesis)
  delta=null.pars$delta            # parameter "delta" for weighted sum of chi-square variables (null hypothesis)
  
  # Statistics and p-values
  pvalue=davies(Stats,lambda,delta=delta,acc=acc,lim=lim)$Qq # /davies(sum(Wt),eig,sigma=sigma)$Qq
  pvalue.SKAT=davies(Stats.SKAT,lambda,acc=acc,lim=lim)$Qq # /davies(sum(Wt),eig,sigma=sigma)$Qq
  return(list(pvalue=pvalue,
              pvalue.SKAT=pvalue.SKAT,
              Stats=Stats,
              Stats.SKAT=Stats.SKAT))
}

# compute covariate matrix (under null hypothesis) and statistics based on observation
obs.func=function(Phen,X,trait,Geno,Wt,mu,kin)
{
  if(trait=="c")
  {
    res.glm=glm(Phen~X)
    est.sd=sigma(res.glm)
    sqr.W.t.G=t(Geno)*sqrt(Wt)   # sqrt(W)%*%t(G)
    sqr.W.t.G.X=sqr.W.t.G%*%X    # sqrt(W)%*%t(G)%*%X
    Sigma=(sqr.W.t.G%*%t(sqr.W.t.G)-sqr.W.t.G.X%*%solve(t(X)%*%X)%*%t(sqr.W.t.G.X))/est.sd^2  # covariance matrix
    sqrt.Stats=(t(Geno)%*%res.glm$residuals/est.sd^2+mu)*sqrt(Wt)
    sqrt.Stats.SKAT=(t(Geno)%*%res.glm$residuals/est.sd^2)*sqrt(Wt)
    Stats=t(sqrt.Stats)%*%sqrt.Stats
    Stats.SKAT=t(sqrt.Stats.SKAT)%*%sqrt.Stats.SKAT
  }
  if(trait=="b")
  { 
    res.glm=glm(Phen~X,family = binomial)
    est.mu = res.glm$fitted.values
    v=est.mu*(1-est.mu)
    
    VG=Geno*sqrt(v) # sqrt(diag(v))%*%Geno 
    VX=X*sqrt(v)    # sqrt(diag(v))%*%X: X is n*(1+k) matrix where k is number of covariates
    temp=t(VG)%*%VX
    GPG=t(VG)%*%VG-temp%*%solve(t(VX)%*%VX)%*%t(temp)
    Sigma=t(t(GPG*sqrt(Wt))*sqrt(Wt))   # diag(sqrt(W))%*%GPG%*%diag(sqrt(W))
    sqrt.Stats=(t(Geno)%*%(Phen-est.mu)+mu)*sqrt(Wt)
    sqrt.Stats.SKAT=(t(Geno)%*%(Phen-est.mu))*sqrt(Wt)
    Stats=t(sqrt.Stats)%*%sqrt.Stats
    Stats.SKAT=t(sqrt.Stats.SKAT)%*%sqrt.Stats.SKAT
  }
  if(trait=="fam.c")
  {
    ID=rownames(kin)
    Data=data.frame(Phen,ID,X)
    if(ncol(X)==1) exprs<-"Phen ~ 1"
    else exprs<-paste("Phen ~", paste(names(Data)[-c(1,2,3)],collapse=" + "))
    # res.lme=lmekin(as.formula(exprs),Data,random=~1|ID,varlist=list(as.matrix(kin)))
    # the above line also output the same result but the process is very slow
    res.lme=lmekin(as.formula(exprs),Data,random=~1|ID,varlist=kin)
    Sigma.lme=res.lme$theta[1]*2*as.matrix(kin)+res.lme$theta[2]*diag(nrow(kin))    # Simga matrix from linear mixed effect model
    S=solve(Sigma.lme)                                      # inverse matrix
    sqr.W.t.G=t(Geno)*sqrt(Wt)   # sqrt(W)%*%t(G) where S is solve(Sigma.lme)
    sqr.W.t.G.X=sqr.W.t.G%*%S%*%X    # sqrt(W)%*%t(G)%*%S%*%X
    # covariance matrix: sqrt(W)*t(G)*S*t(S)*G*sqrt(W)
    Sigma=(sqr.W.t.G%*%S%*%t(sqr.W.t.G)-sqr.W.t.G.X%*%solve(t(X)%*%S%*%X)%*%t(sqr.W.t.G.X))  
  
    resi=res.lme$residuals
    sqrt.Stats=(t(Geno)%*%(S%*%resi)+mu)*sqrt(Wt)
    sqrt.Stats.SKAT=(t(Geno)%*%(S%*%resi))*sqrt(Wt)
    Stats=t(sqrt.Stats)%*%sqrt.Stats
    Stats.SKAT=t(sqrt.Stats.SKAT)%*%sqrt.Stats.SKAT
  }
  return(list(Sigma=Sigma,        # t(Q)%*%P%*%Q where P is covariance matrix of "first derivation". Used to compute null distribution
              Stats=Stats,        # estimated "first derivation" based on model. Used to compute statistics
              Stats.SKAT=Stats.SKAT))  
}
      
null.para=function(Sigma,mu1)  
{
  res=eigen(Sigma,symmetric = T)
  lambda=res$values
  
  lamb.pos=which(lambda>0)
  ok=which(lambda>mean(lambda[lamb.pos])/1e5)

  lambda=lambda[ok]
  V=res$vector[,ok]
  
  delta=(t(V)%*%mu1)^2/lambda
  return(list(lambda=lambda,delta=delta))
}

Input.check=function(Phen,Geno,Cova,Wt,Nu,trait,kin,ped,ibd)
{
  if(!is.matrix(Geno)) stop('Input of "Geno" should be a matrix.')
  if(length(Phen)!=nrow(Geno)) stop('length(Phen) should be the same with nrow(Geno).')
  if(ncol(Geno)!=length(Wt)) stop('ncol(Geno) should be the same with length(Wt).')
  if(ncol(Geno)!=length(Nu)) stop('ncol(Geno) should be the same with length(Nu).')
  if(any(names(Phen)!=rownames(Geno))) stop('names(Phen) should be the same with rownames(Geno).')
  if(!is.null(Cova)){
    if(!is.matrix(Cova)) stop('Input of "Cova" should be a matrix.')
    if(length(Phen)!=nrow(Cova)) stop('length(Phen) should be the same with nrow(Cova).')
  }
  
  Geno=Geno.impute(Geno)
  # remove SNPs with no mutation and check if covariates have linear dependence
  null.SNPs=which(apply(Geno,2,sum)==0)
  if(length(null.SNPs)>0){
    if(length(null.SNPs)==ncol(Geno)) stop("All SNPs within the set have no mutation.")
    Geno=Geno[,-1*null.SNPs,drop=F];
    Wt=Wt[-1*null.SNPs];
    Nu=Nu[-1*null.SNPs]
  }
  # check the linear independence of Covariate matrix
  X=cbind(rep(1,nrow(Cova)),Cova)
  if(qr(X)$rank!=ncol(X)) stop("Please double check the linear dependence of all covariates.")
  # compute kinship matrix and reorder Phenotype, Genotype and/or Covaraite matrix to be consistent with kinship matrix
  if(trait=="fam.c")
  {
    kin=kin.check(kin,ped,ibd)    # here we defaultly treat subjects without kinship matrix infomation as unrelated subjects
    IDs.p=names(Phen)
    IDs.k=dimnames(kin)[[1]]
    pos.k=match(IDs.k, IDs.p)     # position of each subject of kinship in phenotype matrix
    tmpidx=!is.na(pos.k)  
    tmpkin=kin[tmpidx, tmpidx]    # remove those subjects without phenotype   
    pos.unr=setdiff(1:length(Phen),pos.k)
    n.unr=length(pos.unr)
    odr=c(pos.k[tmpidx],pos.unr)  # updated order of subjects
    kin=bdsmatrix(c(tmpkin@blocksize,rep(1,n.unr)),
                  c(tmpkin@blocks,rep(0.5,n.unr)),
                  tmpkin@rmat,
                  list(names(Phen)[odr],names(Phen)[odr]))
    Geno=Geno[odr,];
    Phen=Phen[odr];
    X=X[odr,]
  }
  return(list(Geno=Geno,
              Phen=Phen,
              X=X,
              Wt=Wt,
              Nu=Nu,
              kin=kin))
}


kin.check=function(kin,ped,ibd)
{
  # The priority is "kin" followed by "ped" followed by "ibd"
  if(!is.null(kin))
  {
    if(!is.bdsmatrix(kin)) stop('Input of "kin" should be of class "bdsmatrix"')
    else return(kin)
  } 
  else if(!is.null(ped))
  {
    if(!is.data.frame(ped)|colnames(ped)!=c("famid","id","father.id","mother.id")) 
      stop('Input of "ped" should be a dataframe with 4 columns of c("famid","id","father.id","mother.id")')
    else {
      kin=makekinship(ped$famid,ped$id,ped$father.id,ped$mother.id)
      return(kin) 
    }
  } 
  else if(!is.null(ibd)) 
  {
    if(!is.data.frame(ibd)|any(colnames(ibd)!=c("id1","id2","x"))) 
      stop('Input of "ibd" should be a dataframe with 3 columns of c("id1","id2","x")')
    else {
      kin=bdsmatrix.ibd(ibd$id1,ibd$id2,ibd$x,diagonal = 0.5)
      return(kin)
    }
  }else stop("One of kin, ped and ibd should be specified for family-based analysis.")
}


Geno.impute=function(Geno)
{
  MAF=apply(Geno, 2, mean,na.rm=T)/2
  for(pp in 1:ncol(Geno)) {
    IDX<-which(is.na(Geno[,pp]))
    if(length(IDX)>0) {
      Geno[IDX,pp]=rbinom(length(IDX), 2, MAF[pp])
    }
  }
  return(Geno)
}

Wt.func=function(Geno)
{
  MAF=apply(Geno, 2, mean,na.rm=T)/2
  Wt=dbeta(MAF,1,25)^2
  return(Wt)
}
