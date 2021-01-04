################   Information versus k   ################
library(nabor)
#library(mined)
library(mvtnorm)
#library(OSMAC)
library(ggplot2)
library(foreach)
library(doParallel)
source("GeneData")
source("more_efficient")
#source("MED")
source("IBOSSL")
source("Leverage")
source("DKDO")
source("getMLE")
source("twostep")

######### Conditional on full data #########
set.seed(666)
n=1e4
beta=rep(0.5,7)
D=GeneData(n,beta,"Exp")
X=D$X
#X=log(X)
Y=D$Y

### uniform
unif_inf=function(X,beta){
  w=(1-1/(1+exp(X%*%beta)))*(1/(1+exp(X%*%beta)))
  DI=matrix(0,500,n)
  for(j in 1:500){
    id=sample(1:n,n,replace = FALSE)
    XX=X[id,]
    ww=w[id]
    I=matrix(0,length(beta),length(beta))
    for(i in 1:n){
      I=I+ww[i]*(XX[i,]%*%t(XX[i,]))
      DI[j,i]=sum(diag(I))
    }
  }
  return(apply(DI,2,mean))
}
Unif_I=unif_inf(X,beta)

### MORE(A-opt)
Aopt_inf=function(X,Y,beta){
  p=1-1/(1+exp(X%*%beta))
  w=p*(1-p)
  W=solve(t(X)%*%(X*c(w)))
  pi=abs(Y-p)*sqrt(apply((X%*%W)^2,1,sum))
  pi=pi/sum(pi)
  DI=matrix(0,500,n)
  for(i in 1:500){
    id=sample(1:n,n,FALSE,pi)
    XX=X[id,]
    ww=w[id]
    I=matrix(0,length(beta),length(beta))
    for(j in 1:n){
      I=I+ww[j]*(XX[j,]%*%t(XX[j,]))
      DI[i,j]=sum(diag(I))
    }
  }
  return(apply(DI,2,mean))
}
Aopt_I=Aopt_inf(X,Y,beta)

### MORE(L-opt)
Lopt_inf=function(X,Y,beta){
  p=1-1/(1+exp(X%*%beta))
  w=p*(1-p)
  pi=abs(Y-p)*sqrt(apply(X^2,1,sum))
  pi=pi/sum(pi)
  DI=matrix(0,500,n)
  for(i in 1:500){
    id=sample(1:n,n,FALSE,pi)
    XX=X[id,]
    ww=w[id]
    I=matrix(0,length(beta),length(beta))
    for(j in 1:n){
      I=I+ww[j]*(XX[j,]%*%t(XX[j,]))
      DI[i,j]=sum(diag(I))
    }
  }
  return(apply(DI,2,mean))
}
Lopt_I=Lopt_inf(X,Y,beta)

### leverage
lev_inf=function(X,beta){
  w=(1-1/(1+exp(X%*%beta)))*(1/(1+exp(X%*%beta)))
  WX=c(sqrt(w))*X
  H=(svd(WX)$u)%*%t(svd(WX)$u)
  pi=diag(H)
  t=rank(-pi,ties.method = "first")
  dI=rep(0,n)
  I=matrix(0,length(beta),length(beta))
  for(i in 1:n){
    I=I+w[t==i]*(X[t==i,]%*%t(X[t==i,]))
    dI[i]=sum(diag(I))
  }
  return(dI)
}
Lev_I=lev_inf(X,beta)


### DKDO
dkdo_inf=function(X,beta){
  w=(1-1/(1+exp(X%*%beta)))*(1/(1+exp(X%*%beta)))
  WX=c(sqrt(w))*X
  # k=round(n/2/ncol(WX))
  # t=rank(WX[,1],ties.method ="first")
  # idx=c(which(t<=k),which(t>=length(t)-k+1))
  # for(j in 2:ncol(WX)){
  #   idz=c(1:nrow(WX))[-idx]
  #   z=WX[-idx,]
  #   t=rank(z[,j],ties.method ="first")
  #   idt=idz[c(which(t<=k),which(t>=length(t)-k+1))]
  #   idx=c(idx,idt)
  # }
  # idr=setdiff(1:n,idx)
  # XX=X[c(idx,idr),]
  # ww=w[c(idx,idr)]
  # I=matrix(0,length(beta),length(beta))
  # for(j in 1:n){
  #   I=I+ww[j]*(XX[j,]%*%t(XX[j,]))
  #   dI[j]=sum(diag(I))
  # }
  WXm=apply(WX,2,mean)
  WXd=apply((WX-rep(1,n)%*%t(WXm))^2,1,sum)
  t=rank(-WXd,ties.method = "first")
  dI=rep(0,n)
  I=matrix(0,length(beta),length(beta))
  for(i in 1:n){
    I=I+w[t==i]*(X[t==i,]%*%t(X[t==i,]))
    dI[i]=sum(diag(I))
  }
  return(dI)
}
Dkdo_I=dkdo_inf(X,beta)

Df_data=data.frame(trI=c(Unif_I,Aopt_I,Lopt_I,Lev_I,Dkdo_I),
                   method=rep(c("UNIF","MORE(A-opt)","MORE(L-opt)","LEV","DKDO"),each=n),
                   k=rep(1:n,times=5))
write.csv(Df_data,"Logistic/Main/Information/lExpInformation.csv")

# P=ggplot(Df_data,aes(x=k,y=trI,group=method,colour=method))+
#   geom_line(aes(linetype=method))+
#   scale_linetype_manual(values=c(1,2,1,2,1))+
#   scale_color_manual(values=c("#4DAF4A","#E41A1C","#984EA3","#984EA3","orange"))+
#   xlab("k")+ylab("Information")+theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y = element_text(angle=90),
#         legend.justification=c(0,1),
#         legend.position = c(0,1),
#         legend.background = element_rect(colour = "black"))
P=ggplot(Df_data,aes(x=k,y=trI,group=method,colour=method))+
  geom_line(aes(linetype=method))+
  scale_linetype_manual(values=c(1,2,1,2,1))+
  scale_color_manual(values=c("#4DAF4A","#E41A1C","#984EA3","#984EA3","orange"))+
  xlab("k")+ylab("Information")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(angle=90),
        legend.position = "none",
        legend.background = element_rect(colour = "black"))
P


######### Unconditional #########
detI=function(n,beta,seed){
  ### generate the data
  set.seed(seed)
  D=GeneData(n,beta,"nzNormal")
  X=D$X
  Y=D$Y
  p=1-1/(1+exp(X%*%beta))
  w=p*(1-p)
  DI=matrix(0,5,n)
  
  ### uniform
  I1=matrix(0,length(beta),length(beta))
  for(i in 1:n){
    I1=I1+w[i]*(X[i,]%*%t(X[i,]))
    DI[1,i]=sum(diag(I1))
  }
  
  ### MORE(A-opt)
  W=solve(t(X)%*%(X*c(w)))
  pi=abs(Y-p)*sqrt(apply((X%*%W)^2,1,sum))
  pi=pi/sum(pi)
  id=sample(1:n,n,FALSE,pi)
  XX=X[id,]
  ww=w[id]
  I2=matrix(0,length(beta),length(beta))
  for(i in 1:n){
    I2=I2+ww[i]*(XX[i,]%*%t(XX[i,]))
    DI[2,i]=sum(diag(I2))
  }
  
  ### MORE(L-opt)
  pi=abs(Y-p)*sqrt(apply(X^2,1,sum))
  pi=pi/sum(pi)
  id=sample(1:n,n,FALSE,pi)
  XX=X[id,]
  ww=w[id]
  I3=matrix(0,length(beta),length(beta))
  for(i in 1:n){
    I3=I3+ww[i]*(XX[i,]%*%t(XX[i,]))
    DI[3,i]=sum(diag(I3))
  }
  
  ### leverage
  WX=c(sqrt(w))*X
  H=(svd(WX)$u)%*%t(svd(WX)$u)
  pil=diag(H)
  t=rank(-pil,ties.method = "first")
  I4=matrix(0,length(beta),length(beta))
  for(i in 1:n){
    I4=I4+w[t==i]*(X[t==i,]%*%t(X[t==i,]))
    DI[4,i]=sum(diag(I4))
  }
  
  ### DKDO
  WXm=apply(WX,2,mean)
  WXd=apply((WX-rep(1,n)%*%t(WXm))^2,1,sum)
  t=rank(-WXd,ties.method = "first")
  I5=matrix(0,length(beta),length(beta))
  for(i in 1:n){
    I5=I5+w[t==i]*(X[t==i,]%*%t(X[t==i,]))
    DI[5,i]=sum(diag(I5))
  }
  
  return(DI)
}

cl<- makeCluster(10) 
registerDoParallel(cl) 
ResultNzUInf= foreach(i=1:100,
                  .combine=rbind,
                  .packages=c("mvtnorm")) %dopar% detI(1e4,rep(0.5,7),666*i+18)
stopCluster(cl)

IU=IA=IL=Il=ID=matrix(0,100,1e4)
for(i in 1:100){
  IU[i,]=ResultNzUInf[5*i-4,]
  IA[i,]=ResultNzUInf[5*i-3,]
  IL[i,]=ResultNzUInf[5*i-2,]
  Il[i,]=ResultNzUInf[5*i-1,]
  ID[i,]=ResultNzUInf[5*i,]
}

Df_Data=data.frame(trI=c(apply(IU,2,mean),apply(IA,2,mean),apply(IL,2,mean),apply(Il,2,mean),apply(ID,2,mean)),
                   method=rep(c("UNIF","MORE(A-opt)","MORE(L-opt)","LEV","DKDO"),each=n),
                   k=rep(1:n,times=5))
write.csv(Df_Data,"Logistic/Main/Information/NzUInf.csv")
# P=ggplot(Df_data,aes(x=k,y=trI,group=method,colour=method))+
#   geom_line(aes(linetype=method))+
#   scale_linetype_manual(values=c(1,2,1,2,1))+
#   scale_color_manual(values=c("#4DAF4A","#E41A1C","#984EA3","#984EA3","orange"))+
#   xlab("k")+ylab("Information")+theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y = element_text(angle=90),
#         legend.justification=c(0,1),
#         legend.position = c(0,1),
#         legend.background = element_rect(colour = "black"))
P=ggplot(Df_Data,aes(x=k,y=trI,group=method,colour=method))+
  geom_line(aes(linetype=method))+
  scale_linetype_manual(values=c(1,2,1,2,1))+
  scale_color_manual(values=c("#4DAF4A","#E41A1C","#984EA3","#984EA3","orange"))+
  xlab("k")+ylab("Information")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(angle=90),
        legend.position = "none",
        legend.background = element_rect(colour = "black"))
P


