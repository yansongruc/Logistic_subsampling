################   IBOSS   ################ 
library(mvtnorm)
library(nabor)
library(ggplot2)
library(foreach)
library(doParallel)
source("GeneData")

D=GeneData(1e4,rep(0.5,2),"mzNormal")
X=D$X
Y=D$Y


# delta=2.5 when X follows mzNormal, nzNormal, mixNormal;
# delta=0.5 when X follows T3
IBOSSL=function(X,Y,r0,r1,delta,intercept=c("Yes","No")){
  call=match.call()
  intercept=match.arg(intercept)
  # Stage 1
  n0=length(Y)-sum(Y)
  pprop=rep(1/(2*n0),length(Y))
  pprop[Y==1]=1/(2*sum(Y))
  idx0=sample(1:length(Y),r0,prob = pprop)
  if(intercept=="Yes"){
    beta0=getMLE(cbind(1,X)[idx0,],Y[idx0],1/pprop[idx0])$par
    if(is.na(beta0[1])){
      return(list(ssp.idx=NA,par=NA))
    }
    else{
      c=cbind(1,X)%*%beta0
      f=c*c*(1-1/(1+exp(c)))^(ncol(X)+1)*(1/(1+exp(c)))^(ncol(X)+1)
      cstar=c[which.max(f)]
      idxB=rep(0,length(Y))
      for(i in 1:length(Y)){
        if(min(abs(c[i]-cstar),abs(c[i]+cstar))<=delta)
          idxB[i]=1
      }
      B=which(idxB==1)
      if(length(B)<=r1){
        XB=X[B,]
        idx=1:length(B)
      }
      else{
        # stage 2
        XB=X[B,]
        k=round(r1/2/ncol(X))
        t=rank(XB[,1],ties.method ="first")
        idx=c(which(t<=k),which(t>=length(t)-k+1))
        for(j in 2:ncol(X)){
          idxb=c(1:nrow(XB))[-idx]
          xb=XB[-idx,]
          t=rank(xb[,j],ties.method ="first")
          idt=idxb[c(which(t<=k),which(t>=length(t)-k+1))]
          idx=c(idx,idt)
        }
      }
      beta=getMLE(cbind(1,rbind(X[idx0,],XB[idx,])),c(Y[idx0],Y[B][idx]),1)$par
      return(list(ssp.idx=c(idx0,B[idx]),par=beta))
    }
  }
  else if(intercept=="No"){
    beta0=getMLE(X[idx0,],Y[idx0],1/pprop[idx0])$par
    if(is.na(beta0[1])){
      return(list(ssp.idx=NA,par=rep(NA,ncol(X))))
    }
    else{
      c=X%*%beta0
      f=c*c*(1-1/(1+exp(c)))^(ncol(X)+1)*(1/(1+exp(c)))^(ncol(X)+1)
      cstar=c[which.max(f)]
      idxB=rep(0,length(Y))
      for(i in 1:length(Y)){
        if(min(abs(c[i]-cstar),abs(c[i]+cstar))<=delta)
          idxB[i]=1
      }
      B=which(idxB==1)
      if(length(B)<=r1){
        XB=X[B,]
        idx=1:length(B)
      }
      else{
        # stage 2
        XB=X[B,]
        k=round(r1/2/ncol(X))
        t=rank(XB[,1],ties.method ="first")
        idx=c(which(t<=k),which(t>=length(t)-k+1))
        for(j in 2:ncol(X)){
          idxb=c(1:nrow(XB))[-idx]
          xb=XB[-idx,]
          t=rank(xb[,j],ties.method ="first")
          idt=idxb[c(which(t<=k),which(t>=length(t)-k+1))]
          idx=c(idx,idt)
        }
      }
      beta=getMLE(rbind(X[idx0,],XB[idx,]),c(Y[idx0],Y[B][idx]),1)$par
      return(list(ssp.idx=c(idx0,B[idx]),par=beta))
    }
  }
}

