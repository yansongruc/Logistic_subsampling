################   DKDO   ################
library(mvtnorm)
library(ggplot2)
library(foreach)
library(doParallel)
source("GeneData")

D=GeneData(1e4,c(0.5,0.5),"mzNormal")
X=D$X
Y=D$Y

########   without interception   ########
DKDO_Knn=function(X,Y,K,sub){
  NN=nabor::knn(X,X,K)$nn.idx
  phat=rep(0,length(Y))
  for(i in 1:length(Y)){
    phat[i]=sum(Y[NN[i,]])/K
  }
  what=phat*(1-phat)
  
  Z=X
  for(i in 1:nrow(X)){
    Z[i,]=sqrt(what[i])*X[i,]
  }
  k=round(sub/2/ncol(Z))
  t=rank(Z[,1],ties.method ="first")
  idx=c(which(t<=k),which(t>=length(t)-k+1))
  for(j in 2:ncol(Z)){
    idz=c(1:nrow(Z))[-idx]
    z=Z[-idx,]
    t=rank(z[,j],ties.method ="first")
    idt=idz[c(which(t<=k),which(t>=length(t)-k+1))]
    idx=c(idx,idt)
  }
  
  parhat=getMLE(X[idx,],Y[idx],1)$par
  return(list(ssp.idx=idx,par=parhat))
}


DKDO_pilot=function(X,Y,r0,r1){
  # obtain a pilot estimator
  n0=length(Y)-sum(Y)
  pprop=rep(1/(2*n0),length(Y))
  pprop[Y==1]=1/(2*sum(Y))
  idx0=sample(1:length(Y),r0,prob = pprop)
  beta0=getMLE(X[idx0,],Y[idx0],1/pprop[idx0])$par
  if(is.na(beta0[1])){
    return(list(ssp.idx=NA,par=NA))
  }
  else{
    # 
    phat=1-1/(1+exp(X%*%beta0))
    what=phat*(1-phat)
    # if(intercept=="Yes"){
    #   XX=X[,-1]
    # }
    # else if(intercept=="No"){
    #   XX=X
    # }
    Z=X
    for(i in 1:nrow(X)){
      Z[i,]=sqrt(what[i])*X[i,]
    }
    k=round(r1/2/ncol(Z))
    loc=Z[idx0,]
    idx=idx0
    for(j in 1:ncol(Z)){
      idz=c(1:nrow(Z))[-idx]
      z=Z[-idx,]
      t=rank(z[,j],ties.method ="first")
      idt=idz[c(which(t<=k),which(t>=length(t)-k+1))]
      idx=c(idx,idt)
    }
    parhat=getMLE(X[idx,],Y[idx],1)$par
    return(list(ssp.idx=idx,par=parhat))
  }
  # if(intercept=="Yes"){
  #   parhat=getMLE(XX[idx,],Y[idx],1)$par
  #   return(list(ssp.idx=c(idx,idx0),par=c(beta0[1],parhat)))
  # }
  # else if(intercept=="No"){
  #   parhat=getMLE(XX[c(idx,idx0),],Y[c(idx,idx0)],1)$par
  #   return(list(ssp.idx=c(idx,idx0),par=parhat))
  # }
}

