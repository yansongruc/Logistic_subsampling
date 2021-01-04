################   Information versus n   ################
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
Inform2=function(n){
  set.seed(100)
  beta=rep(0.5,7)
  D=GeneData(n,beta,"Exp")
  X=D$X
  X=log(X)
  Y=D$Y
  p=1-1/(1+exp(X%*%beta))
  w=p*(1-p)
  res=rep(0,5)
  
  # uniform
  r1=rep(0,500)
  for(i in 1:500){
    id=sample(1:n,500,replace = FALSE)
    I=t(X[id,])%*%(X[id,]*c(w[id]))
    r1[i]=sum(diag(I))
  }
  res[1]=mean(r1)
  
  # Aopt
  W=solve(t(X)%*%(X*c(w)))
  pi=abs(Y-p)*sqrt(apply((X%*%W)^2,1,sum))
  pi=pi/sum(pi)
  r2=rep(0,500)
  for(i in 1:500){
    id=sample(1:n,500,FALSE,pi)
    I=t(X[id,])%*%(X[id,]*c(w[id]))
    r2[i]=sum(diag(I))
  }
  res[2]=mean(r2)
  
  # Lopt
  pi=abs(Y-p)*sqrt(apply(X^2,1,sum))
  pi=pi/sum(pi)
  r3=rep(0,500)
  for(i in 1:500){
    id=sample(1:n,500,FALSE,pi)
    I=t(X[id,])%*%(X[id,]*c(w[id]))
    r3[i]=sum(diag(I))
  }
  res[3]=mean(r3)
 
  # LEV
  WX=c(sqrt(w))*X
  H=(svd(WX)$u)%*%t(svd(WX)$u)
  pi=diag(H)
  t=rank(-pi,ties.method = "first")
  id=which(t<=500)
  I=t(X[id,])%*%(X[id,]*c(w[id]))
  res[4]=sum(diag(I))
  
  # DKDO
  WXm=apply(WX,2,mean)
  WXd=apply((WX-rep(1,n)%*%t(WXm))^2,1,sum)
  t=rank(-WXd,ties.method = "first")
  id=which(t<=500)
  I=t(X[id,])%*%(X[id,]*c(w[id]))
  res[5]=sum(diag(I))
  
  return(res)
}

nseq=c(600,800,c(2:10)*500)
A=sapply(nseq,Inform2)
lExpInfn=data.frame(trI=c(A),
                  method=rep(c("UNIF","MORE(A-opt)","MORE(L-opt)","LEV","DKDO"),times=length(nseq)),
                  n=rep(nseq,each=5))
write.csv(lExpInfn,"Logistic/Main/Information/lExpInfn.csv")

# P=ggplot(MzInfn,aes(x=n,y=trI,group=method,colour=method))+
#   geom_line(aes(linetype=method))+
#   geom_point(aes(shape=method),size=2)+
#   scale_shape_manual(values=c(1,2,3,4,5))+
#   scale_linetype_manual(values=c(1,2,1,2,1))+
#   scale_color_manual(values=c("#4DAF4A","#E41A1C","#984EA3","#984EA3","orange"))+
#   xlab("n")+ylab("Information")+theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y = element_text(angle=90),
#         legend.justification=c(0,1),
#         legend.position = c(0,1),
#         legend.background = element_rect(colour = "black"))
P=ggplot(lExpInfn,aes(x=n,y=trI,group=method,colour=method))+
  geom_line(aes(linetype=method))+
  geom_point(aes(shape=method),size=2)+
  scale_shape_manual(values=c(1,2,3,4,5))+
  scale_linetype_manual(values=c(1,2,1,2,1))+
  scale_color_manual(values=c("#4DAF4A","#E41A1C","#984EA3","#984EA3","orange"))+
  xlab("n")+ylab("Information")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(angle=90),
        legend.position = "none",
        legend.background = element_rect(colour = "black"))
P
