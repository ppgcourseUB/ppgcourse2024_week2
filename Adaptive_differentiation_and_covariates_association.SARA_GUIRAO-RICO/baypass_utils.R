###############################################
###function: Simulate under the inference model
##############################################
simulate.baypass <- function(omega.mat,
                             nsnp=1000,
                             beta.coef=NA,
                             beta.pi=c(1,1),
                             pop.trait=0,
                             sample.size=100,
                             pi.maf=0.05,
                             suffix="sim",
                             print.sim.params.values=FALSE,
                             remove.fixed.loci=FALSE,
                             output.bayenv.format=FALSE,
                             coverage=NA){
#sample.size=matrix with npop colums or vector of length npop with count. If matrix and poolseq data, the population sample size (i.e., haploid size of the pools for Pool-Seq data) is set to the colMax of sample.size
#coverage = matrix with npop colums or vector of length npop with coverage => activate poolseq data

 require(mvtnorm)
 if(!is.matrix(omega.mat)){stop("omega.mat must be a matrix (e.g., as obtained from *mat_omega.out BayPass output)")}
 npop=nrow(omega.mat)
 
 if(sum(is.na(beta.coef))==0){
  if(length(pop.trait)!=npop){stop("Trait dimension must have the same size as the rank of omega.mat")}
  simu.cov=TRUE
  nsnp.asso=length(beta.coef)
  if(nsnp>0){beta.coef=c(rep(0,nsnp),beta.coef)}
  nsnp=nsnp + nsnp.asso
  alpha.cov = beta.coef %*% t(pop.trait)
 }else{simu.cov=FALSE}
 
 if(length(sample.size)==1){
  NN=matrix(sample.size,nsnp,npop)
  poolsize=rep(sample.size,npop)
  }else{
  sample.size=as.matrix(sample.size)
  if(ncol(sample.size)==1){#c'est un vecteur
    if(nrow(sample.size)!=npop){stop("Sample size dimension must be of length 1 or have the same size as the rank of omega.mat or must be a matrix with npop columns")}
    NN=matrix(rep(sample.size,nsnp),nsnp,npop,byrow=TRUE)
    poolsize=as.numeric(sample.size)
  }else{
   tmp.snp=sample(1:nrow(sample.size),nsnp,replace=TRUE)
   NN=sample.size[tmp.snp,]
   poolsize=apply(sample.size,2,max)
  }}

 if(length(coverage)==1){
  if(is.na(coverage)){
  poolseq=FALSE
  }else{poolseq=TRUE ; NN.coverage=matrix(coverage,nsnp,npop)}
  }else{
   poolseq=TRUE
   coverage=as.matrix(coverage)
   if(ncol(coverage)==1){#c'est un vecteur
    if(nrow(coverage)!=npop){
     stop("Coverage dimension must be of length 1, or have the same size as the rank of omega.mat or must be a matrix with npop columns")
    }
   NN.coverage=matrix(rep(coverage,nsnp),nsnp,npop,byrow=TRUE)
  }else{
   tmp.snp=sample(1:nrow(coverage),nsnp,replace=TRUE)
   NN.coverage=coverage[tmp.snp,]
  }
 }
 
 if(poolseq){
  YY.coverage=NN.coverage*0
  NN=matrix(rep(poolsize,nsnp),nsnp,npop,byrow=TRUE)
  }

 if(length(beta.pi)!=2){stop("beta.pi must be of length 2")}else{
  if(sum(beta.pi==1)==2){Pi=runif(nsnp,pi.maf,1-pi.maf)}else{
   Pi=rbeta(nsnp,beta.pi[1],beta.pi[2])
   Pi[Pi<pi.maf]=pi.maf ; Pi[Pi>1-pi.maf]=1-pi.maf
  }
 }
 
 ALPHA=YY=matrix(0,nsnp,npop)
 for(i in 1:nsnp){
   mat=Pi[i]*(1-Pi[i])*omega.mat
   ALPHA[i,]=rmvnorm(1,rep(Pi[i],npop),mat)
 }
 if(simu.cov){ALPHA=ALPHA + alpha.cov}
 ALPHA_tr=ALPHA ; ALPHA_tr[ALPHA>1]=1 ; ALPHA_tr[ALPHA<0]=0 
 for(i in 1:nsnp){
  for(j in 1:npop){
    YY[i,j]=rbinom(1,size=NN[i,j],prob=ALPHA_tr[i,j])
    if(poolseq){
     YY.coverage[i,j]=rbinom(1,size=NN.coverage[i,j],prob=YY[i,j]/NN[i,j])
    }
  }
  if(i%%(nsnp/10)==0){cat(i,"SNP simulated out of",nsnp,"\n")}
  }
  
 if(output.bayenv.format){remove.fixed.loci=TRUE}#Bayenv doesn't accept fixed loci
 
 if(remove.fixed.loci){
 if(poolseq){tmp.freq=rowSums(YY.coverage)/rowSums(NN.coverage)}else{tmp.freq=rowSums(YY)/rowSums(NN)}
 snp.sel=tmp.freq>0 & tmp.freq<1
 nsnp=sum(snp.sel)
 YY=YY[snp.sel,] ; NN=NN[snp.sel,] ; Pi=Pi[snp.sel] ; ALPHA=ALPHA[snp.sel,] 
 if(poolseq){YY.coverage=YY.coverage[snp.sel,];NN.coverage=NN.coverage[snp.sel,]}
 if(simu.cov){beta.coef=beta.coef[snp.sel]}
 cat("Number of SNPs removed: ",sum(!snp.sel),"\n")
 }

 if(print.sim.params.values){
  write.table(file=paste("pi.",suffix,sep=""),Pi,quote=F,col.names=F,row.names=F)
  write.table(file=paste("alpha.",suffix,sep=""),ALPHA,quote=F,col.names=F,row.names=F)

  if(simu.cov){
   write.table(file=paste("betacoef.",suffix,sep=""),beta.coef,quote=F,col.names=F,row.names=F)
   write.table(file=paste("pheno.",suffix,sep=""),t(pop.trait),quote=F,col.names=F,row.names=F)
  }
 }

 mat_baypass=cbind(YY,NN)
 all2.pos=2*(1:npop)
 mat_baypass[,all2.pos-1]=YY ; mat_baypass[,all2.pos]=NN-YY
 
 if(!poolseq | print.sim.params.values){
  write.table(file=paste("G.",suffix,sep="") ,mat_baypass,quote=F,col.names=F,row.names=F)
  if(output.bayenv.format){
   mat_bayenv=rbind(YY,NN)
   all2.pos=2*(1:nsnp)
   mat_bayenv[all2.pos-1,]=YY ; mat_bayenv[all2.pos,]=NN-YY
   mat_bayenv=cbind(mat_bayenv,rep("",nsnp))
   write.table(file=paste("bayenv_freq.",suffix,sep=""),mat_bayenv,sep="\t",quote=F,col.names=F,row.names=F)
  }
 }
  
 if(poolseq){
   mat_baypass=cbind(YY.coverage,NN.coverage)
   all2.pos=2*(1:npop)
   mat_baypass[,all2.pos-1]=YY.coverage ; mat_baypass[,all2.pos]=NN.coverage-YY.coverage
   write.table(file=paste("Gpool.",suffix,sep="") ,mat_baypass,quote=F,col.names=F,row.names=F)

  if(output.bayenv.format){
   mat_bayenv=rbind(YY.coverage,NN.coverage)
   all2.pos=2*(1:nsnp)
   mat_bayenv[all2.pos-1,]=YY.coverage ; mat_bayenv[all2.pos,]=NN.coverage-YY.coverage
   mat_bayenv=cbind(mat_bayenv,rep("",nsnp))
   write.table(file=paste("bayenv_freq_pool.",suffix,sep=""),mat_bayenv,sep="\t",quote=F,col.names=F,row.names=F)
  }
  write.table(file=paste("poolsize.",suffix,sep=""),t(poolsize),quote=F,col.names=F,row.names=F)
 }

 if(simu.cov){
  if(poolseq){
  list(Y.pool=YY.coverage,N.pool=NN.coverage,Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat,betacoef.sim=beta.coef)
  }else{
  list(Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat,betacoef.sim=beta.coef)
  }
 }else{
  if(poolseq){
   list(Y.pool=YY.coverage,N.pool=NN.coverage,Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat)
  }else{
   list(Y.sim=YY,N.sim=NN,Pi.sim=Pi,alpha.sim=ALPHA,omega.sim=omega.mat)
  }
}
}

###############################################
###function: transform geno input files into YY and NN
##############################################

geno2YN<-function(genofile){
 data=read.table(genofile)
 npop=ncol(data)/2 ; all1=seq(1,2*npop,2) ; all2=all1+1
 YY=data[,all1]
 NN=YY+data[,all2]
 list(YY=YY,NN=NN)
}

###############################################
###function: Compute FMD distance (Forstner and Moonen, 2003) between two covariance matrices
##############################################
fmd.dist<-function(mat1,mat2){
 require(geigen)
 return(sqrt(sum(log(geigen(mat1,mat2)$values)**2)))
}

###############################################
###function: Simulate covariate value correlated to a given Omega PC
##############################################
simulate.PCcorrelated.covariate <- function(omega,axis=1,targeted.rho=0.1,tol=0.01){
 npops=nrow(omega)
 om.svd=svd(omega)
 PC=om.svd$u[,axis]
 PC.scaled=scale(PC)
 ###hybrid ALGO between brut force (since PC not necessarily gaussian) and projection
 #See https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable 
 tol=tol*abs(targeted.rho)
 target.min=max(-1,targeted.rho-tol)
 target.max=min(1,targeted.rho+tol)
 C = matrix(targeted.rho,2,2)
 diag(C) = 1
 C=chol(C)
 cc=100 ; cnt=0
 while(cc<target.min | cc>target.max){
   nn=(cbind(PC.scaled,rnorm(npops))%*% C)[,2] 
   cc=cor(PC,nn)
   cnt=cnt+1
 }
 cat(targeted.rho," found in ",cnt," iterations\n")   
 return(nn)
}

###############################################
###function: Plot Omega after spectral decomposition: Omega=UDU' with D=diagonal matrix of singular values
##############################################
plot.omega <- function(omega,PC=c(1,2),pop.names=paste0("Pop",1:nrow(omega)),main=expression("SVD of "*Omega),col=rainbow(nrow(omega)),pch=16,pos=2){
 om.svd=svd(omega)
 eig=om.svd$d
 pcent.var=100*eig/sum(eig)
 plot(om.svd$u[,PC],main=main,pch=pch,col=col,
      xlab=paste0("PC",PC[1]," (",round(pcent.var[PC[1]],2),"%)"),
      ylab=paste0("PC",PC[2]," (",round(pcent.var[PC[2]],2),"%)")
      )
 text(om.svd$u[,PC[1]],om.svd$u[,PC[2]],pop.names,col=col,pos=pos)
 list(PC=om.svd$u,eig=eig,pcent.var=pcent.var)
 }

#################################################
#####fonction: compute genetic offset from estimated regression coefficient
#################################################
compute_genetic_offset<-function(beta.coef=NULL,regfile="summary_betai.out",candidate.snp=NULL,
                                 covfile="cov.baypass",newenv=NULL,refenv=NULL,scalecov=TRUE,compute.rona=FALSE){
  #adapted from genetic.gap function of the LEA package to compute genetic gap (Gain and Francois, 2023)
  #beta.coef=matrix of nsnp X ncovariates coefficient. If NULL, user may provide a BayPass output file
  #regfile = BayPass output file (either [outprefix_]summary_betai_reg.out or [outprefix_]summary_betai.out if IS or MC estimates respectively)
  #covfile = BayPass covariate file (-efile) 
  #newenv = covariate values for the new environment(s): either a vector of ncovariates is only a single new environment or a matrix of with covariate values for several new environments in column (i.e., ncovariates X nenvironments)
  #refenv = covariable values for the ref environment(s) to be used instead of covfile (i.e.: default: ref environment = those in cov files specified by each row of covariable value). Same format as newenv 
  require(data.table)
  ##reading and checking covariate
  cov=as.matrix(read.table(covfile))
  n.cov=nrow(cov) ; n.pop=ncol(cov)
  ##reading and checking newenv
  if(is.null(newenv)){
    stop("Please Provide a vector or a matrix (n_covariates X n_newenv) of covariate values for the target environment(s)\n")
  }else{newenv=as.matrix(newenv)}
  if(nrow(newenv)!=n.cov){
    stop("Check newenv argument: The number of covariates for the new environment is not equal to the number of covariates in the original covariate file\n")
  }
  n.newenv=ncol(newenv)
  ##Scaling pop and new covariates 
  if(scalecov){
    m.cov=rowMeans(cov) ; sd.cov=apply(cov,1,sd)
    cov=(cov-m.cov)/sd.cov
    newenv=(newenv-m.cov)/sd.cov
  }
  ##reading and checking refenv
  if(is.null(refenv)){##refenv is set to cov if not provided (default)
    cat("Genetic Offset statistics will be estimated between each of the",n.pop,"population environment (specified by",n.cov,"covariables)\n provided in the original covariate file and the",n.newenv,"target environments provided with newenv argument.\n
  The resulting GO matrix will contain",n.pop,"reference environment (row) times",n.newenv,"target environments (column) entries\n")
    refenv=cov ; n.refenv=n.pop
  }else{
    refenv=as.matrix(refenv)
    if(nrow(refenv)!=n.cov){
      stop("Check refenv argument: the number of covariates for the reference environment is not equal to the number of covariates in the original covariate file\n")
    }
    n.refenv=ncol(refenv)
    ##Scaling pop and new covariates 
    if(scalecov){refenv=(refenv-m.cov)/sd.cov}
    cat("Genetic Offset statistics will be estimated between each of the",n.refenv,"reference environment (specified by",n.cov,"covariables)\n provided with refenv argument and the",n.newenv,"target environments provided with newenv argument.\n
  The resulting GO matrix will contain",n.refenv,"reference environment (row) times",n.newenv,"target environments (column) entries\n")
  }
  ##reading regression coefficient
  if(is.null(beta.coef)){
    type=NULL
    if(grepl("summary_betai_reg.out",regfile)){
      type="IS"
      beta.coef=matrix(fread(regfile,data.table=FALSE,header=T)$Beta_is,ncol=n.cov) #nsnp,ncov
    }
    if(grepl("summary_betai.out",regfile)){
      type="MC" 
      beta.coef=matrix(fread(regfile,data.table=FALSE,header=T)$M_Beta,ncol=n.cov) #nsnp,ncov
    }
    if(is.null(type)){
      stop("Check regfile argument: No proper BayPass output files with estimates of regression coefficient (i.e., either (either [outprefix_]summary_betai_reg.out or [outprefix_]summary_betai.out if IS or MC estimates respectively) could be found.\nPlease provide a valid BayPass output file (regfile argument) OR a matrix of estimated regression coefficient (beta.coef argument)")
    }
  }else{
    beta.coef=as.matrix(beta.coef)
    if(ncol(beta.coef)!=n.cov){stop("ERROR: The number of covariates (columns) of the beta.coef matrix is not equal to the number of covariates in the original covariate file\n")} 
  }
  
  ###################
  ggap=matrix(0,n.refenv,n.newenv)
  if(!is.null(candidate.snp)){beta.coef=beta.coef[candidate.snp,]}
  BtB=t(beta.coef)%*%beta.coef/nrow(beta.coef) # cov(beta.coef)
  eig <- eigen(BtB,symmetric = TRUE) #eigendecompositon of BtB=(beta.coef*t(beta.coef))/nsnp=UDU' car symetrique
  
  for(i in 1:n.refenv){#usually more target than ref
    #    M=t(cov-newenv[,i])%*%t(beta.coef) #sum of reg. coeff for each covariable weighted by the difference between the current and new env. <=> sum. of expected difference in predicted and current allele frequency! 
    # E(f_cur_i - f_new_i) = (Pi_i + S_j(beta_ij*Ycur_j)) - (Pi_i + S_j(beta_ij*Ynew_j)) = S_j(beta_ij*(Ycur_j-Ynew_j))     
    #    ggap[,i] = rowMeans(M**2) #diag(M %*% t(M))/nrow(beta.coef)
    #    mean.diff.freq[,i]=rowMeans(M)
    diff=newenv-refenv[,i]
    ggap[i,] = colSums(eig$values*(t(eig$vectors)%*%diff)**2)  #GG=(x-x*)%*%BtB%*%(x-x*)'=((x-x*)%*%U)%*%D%*%((x-x*)%*%U)'
  }
  #covariable importance
  covimp=(eig$vectors**2)%*%(eig$values)
  out=list(go = ggap,BtB.eigenvalues = eig$values, BtB.eigenvectors = eig$vectors,covimp=as.vector(covimp))
  
  if(compute.rona){
    out$rona=matrix(0,n.refenv,n.newenv)
    for(i in 1:n.refenv){out$rona[i,]=rowMeans(abs(t(newenv-refenv[,i])%*%t(beta.coef)))}
  } 
  
  return(out)
}

###############################################
###function: concatenate sub-datasets results
############################################## 
concatenate_res<-function(dir="./",anaprefix="ana",extension="",nsubsets=2,
                          snpdet_prefix="./detsnp.sub",retrieve_pi_xtx=TRUE,retrieve_bfis=TRUE,retrieve_c2=FALSE){
  #extension should be the same for snpdet files and baypass output files (i.e., all files should be compressed the same way or not)
  #snp_det_prefix should include the path (since it may be different than the directory containing BayPass output files)
  require(data.table)
  corepref=""
  if(nchar(anaprefix)>0){corepref=paste0(anaprefix,"_")}
  if(nsubsets<1){stop("Check nsubsets: at least one subset is needed\n")}
  
  for(i in 1:nsubsets){
    cat("Processing run",i,"out of",nsubsets,"\n")
    tmp.snpdet=fread((paste0(snpdet_prefix,i,extension)),data.table = F)[,1:2]
    colnames(tmp.snpdet)=c("CHR","POS")
    tmp.nsnps=nrow(tmp.snpdet)
    if(retrieve_pi_xtx){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_pi_xtx.out",extension),data.table=F)[,c("M_P","M_XtX","XtXst")]
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(retrieve_bfis){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_betai_reg.out",extension),data.table=F)$"BF(dB)"
      tmp=matrix(as.numeric(tmp),tmp.nsnps)      
      colnames(tmp)=paste0("BFis_cov_",1:ncol(tmp))
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(retrieve_c2){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_contrast.out",extension),data.table=F)$"C2"
      tmp=matrix(tmp,tmp.nsnps)
      colnames(tmp)=paste0("C2_",1:ncol(tmp))
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(i==1){res=tmp.snpdet}else{res=rbind(res,tmp.snpdet)}
  }
  res=res[order(res$CHR,res$POS),]
  return(res)
}
