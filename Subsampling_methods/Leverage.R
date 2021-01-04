################   Leverage   ################
library(nabor)
library(mvtnorm)
library(ggplot2)
library(foreach)
library(doParallel)
source("GeneData")

### without intercept
LevDeter_Knn=function(X,Y,K,sub){
  NN=nabor::knn(X,X,K)$nn.idx
  phat=rep(0,length(Y))
  for(i in 1:length(Y)){
    phat[i]=sum(Y[NN[i,]])/K
  }
  what=phat*(1-phat)
  WX=diag(sqrt(what))%*%X
  H=svd(WX)$u
  h=apply(H^2,1,sum)
  t=rank(-h,ties.method ="first")
  idx=which(t<=sub)
  
  parhat=getMLE(X[idx,],Y[idx],1)$par
  return(list(ssp.idx=idx,par=parhat))
}


LevDeter_Pilot=function(X,Y,r0,r1){
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
    #   XX[,1]=X[,2]+beta0[1]/beta0[2]
    # }
    # else if(intercept=="No"){
    #   XX=X
    # }
    WX=diag(c(sqrt(what)))%*%X
    H=svd(WX)$u
    h=apply(H^2,1,sum)
    t=rank(-h,ties.method ="first")
    idx=which(t<=r1)
    parhat=getMLE(X[c(idx,idx0),],Y[c(idx,idx0)],1)$par
    return(list(ssp.idx=c(idx,idx0),par=parhat))
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

LevDeter_Pilotb=function(X,Y,r0,r1){
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
    #   XX[,1]=X[,2]+beta0[1]/beta0[2]
    # }
    # else if(intercept=="No"){
    #   XX=X
    # }
    WX=X
    for(i in 1:nrow(X)){
      WX[i,]=sqrt(what[i])*X[i,]
    }
    MT=solve(t(WX[idx0,])%*%WX[idx0,])
    h=rep(0,length(Y))
    for(i in 1:length(Y)){
      h[i]=t(WX[i,])%*%MT%*%WX[i,]
    }
    t=rank(-h,ties.method ="first")
    idx=which(t<=r1)
    parhat=getMLE(X[c(idx,idx0),],Y[c(idx,idx0)],1)$par
    return(list(ssp.idx=c(idx,idx0),par=parhat))
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



