################   More efficient estimation based on Poisson subsampling   ################ 
library(mvtnorm)
library(ggplot2)
library(foreach)
library(doParallel)


ME=function(X,Y,r0,r1){
  #Step 1: obtain the pilot \hat\beta_0
  n0=length(Y)-sum(Y)
  pprop=rep(1/(2*n0),length(Y))
  pprop[Y==1]=1/(2*sum(Y))
  u0=runif(length(Y))
  idx0=which(u0<r0*pprop)
  wei0=rep(1,length(idx0))
  for(i in 1:length(idx0)){
    if(r1*pprop[idx0[i]]>1)
      wei0[i]=r1*pprop[idx0[i]]
  }
  betap0=beta0=getMLE(X[idx0,],Y[idx0],wei0)$par
  if(is.na(beta0[1])){
    return(list(ssp.idx=NA,par=NA))
  }
  else{
    beta0[1]=beta0[1]+log(sum(Y)/n0)
    
    # Step 2: obtain the more efficient estimator \hat\beta_{uw}
    phat=1-1/(1+exp(X%*%beta0))
    what=phat*(1-phat)
    wei=rep(1,length(idx0))
    for(i in 1:length(idx0)){
      if(r0*pprop[idx0[i]]<1)
        wei[i]=r0*pprop[idx0[i]]
    }
    W0=length(Y)*solve(t(X[idx0,])%*%(X[idx0,]*c(what[idx0,]/wei)))
    PI.mMSE=sqrt((Y-phat)^2*rowSums((X%*%W0)^2))
    phi0=sum(PI.mMSE[idx0]/wei)
    PI.p=PI.mMSE/phi0
    u=runif(length(Y))
    idx=which(u<r1*PI.p)
    wei1=rep(1,length(idx))
    for(i in 1:length(idx)){
      if(r1*PI.p[idx[i]]>1)
        wei1[i]=r1*PI.p[idx[i]]
    }
    betap=getMLE(X[idx,],Y[idx],wei1)$par
    beta=betap+beta0
    
    # Step 3: combine the two estimators \hat\beta_0 and \hat\beta_p
    what1=(1-1/(1+exp(X%*%betap0)))*(1/(1+exp(X%*%betap0)))
    M1=t(X[idx0,])%*%(X[idx0,]*c(what1[idx0,]))
    what2=(1-1/(1+exp(X%*%betap)))*(1/(1+exp(X%*%betap)))
    M2=t(X[idx,])%*%(X[idx,]*c(what2[idx,]))
    betaf=solve(M1+M2)%*%(M1%*%beta0+M2%*%beta)
    
    return(list(ssp.idx=c(idx0,idx),par=betaf))
  }
  # gg=0
  # while (gg<=1) {
  #   gt=getMLE(X[idx0,],Y[idx0],wei0)
  #   if(gt$message=="Successful convergence"){
  #     betap0=beta0=gt$par
  #     gg=gg+1
  #   }
  # }
  
}
