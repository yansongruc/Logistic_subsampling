################   MED   ################ 
library(mvtnorm)
library(nabor)
library(mined)
library(ggplot2)
library(foreach)
library(doParallel)
source("GeneData")
source("test7_1")
source("test5_1")

D=GeneData(1e4,rep(0.5,2),"mixNormal")
X=D$X
Y=D$Y

MEDkw=function(X,Y,K,sub,intercept=c("Yes","No")){
  #Step 1: obtain the pilot \hat pi
  NN=nabor::knn(X,X,K)$nn.idx
  phat=rep(0,length(Y))
  for(i in 1:length(Y)){
    phat[i]=sum(Y[NN[i,]])/K
  }
  
  # Step 2: obtain the more efficient estimator \hat\beta_{w}
  what=phat*(1-phat)
  What=solve(t(X)%*%(X*c(what)))
  PI.mMSE=sqrt((2*phat*(1-phat))^2*rowSums((X%*%What)^2))
  PI.mMSE=PI.mMSE/sum(PI.mMSE)
  
  call=match.call()
  intercept=match.arg(intercept)
  if(intercept=="Yes"){
    loc=SelectMinED(X[,-1],log(PI.mMSE),sub,1,2)$points
    idx=nabor::knn(X[,-1],loc,1)$nn.idx
  }
  else if(intercept=="No"){
    loc=SelectMinED(X,log(PI.mMSE),sub,1,2)$points
    idx=nabor::knn(X,loc,1)$nn.idx
  }
  beta=getMLE(X[idx,],Y[idx],1/PI.mMSE[idx])$par
  return(list(ssp.idx=idx,par=beta))
}




MEDpw=function(X,Y,r0,r1){
  #Step 1: obtain the pilot estimator
  n0=length(Y)-sum(Y)
  pprop=rep(1/(2*n0),length(Y))
  pprop[Y==1]=1/(2*sum(Y))
  set.seed(0)
  idx0=sample(1:length(Y),r0,prob = pprop)
  beta0=getMLE(X[idx0,],Y[idx0],1/pprop[idx0])$par
  
  # Step 2: obtain the more efficient estimator \hat\beta_{w}
  phat=1-1/(1+exp(X%*%beta0))
  what=phat*(1-phat)
  What=solve(t(X)%*%(X*c(what)))
  PI.mMSE=sqrt((2*phat*(1-phat))^2*rowSums((X%*%What)^2))
  PI.mMSE=PI.mMSE/sum(PI.mMSE)
  loc=SelectMinED(X,log(PI.mMSE),r1,1,2)$points
  idx=nabor::knn(X,loc,1)$nn.idx
  #beta=getMLE(X[idx,],Y[idx],1/PI.mMSE[idx])$par
  beta=getMLE(X[c(idx,idx0),],Y[c(idx,idx0)],1/c(PI.mMSE[idx],pprop[idx0]))$par
  return(list(ssp.idx=idx,par=beta))
}

## not good enough
MEDpuw=function(X,Y,r0,r1){
  #Step 1: obtain the pilot \hat beta
  n0=length(Y)-sum(Y)
  pprop=rep(1/(2*n0),length(Y))
  pprop[Y==1]=1/(2*sum(Y))
  loc0=SelectMinED(X,log(pprop),r0,1,2)$points
  idx0=nabor::knn(X,loc0,1)$nn.idx
  beta0=beta0t=getMLE(X[idx0,],Y[idx0],1)$par
  beta0[1]=beta0[1]+log(sum(Y)/n0)
  
  # Step 2: obtain the more efficient estimator \hat\beta_{uw}
  phat=1-1/(1+exp(X%*%beta0))
  what=phat*(1-phat)
  What=length(Y)*solve(t(X[idx0,])%*%(X[idx0,]*c(what[idx0,])))
  PI.mMSE=sqrt((2*phat*(1-phat))^2*rowSums((X%*%What)^2))
  PI.mMSE=PI.mMSE/sum(PI.mMSE)
  loc=SelectMinED(X,log(PI.mMSE),r1,1,2)$points
  idx=nabor::knn(X,loc,1)$nn.idx
  betat=getMLE(X[idx,],Y[idx],1)$par
  beta=betat+beta0
  
  # Step 3: combine the two estimators \hat\beta_0 and \hat\beta_p
  what1=(1-1/(1+exp(X%*%beta0t)))*(1/(1+exp(X%*%beta0t)))
  M1=t(X[idx0,])%*%(X[idx0,]*c(what1[idx0,]))
  what2=(1-1/(1+exp(X%*%betat)))*(1/(1+exp(X%*%betat)))
  M2=t(X[idx,])%*%(X[idx,]*c(what2[idx,]))
  betaf=solve(M1+M2)%*%(M1%*%beta0+M2%*%beta)
  return(list(ssp.idx=idx,par=beta))
}
