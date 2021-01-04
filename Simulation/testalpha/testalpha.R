################   alpha value   ################
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


testalpha=function(X,Y,beta,p,Z,rr,r){
  n=nrow(X)
  res=rep(0,6)
  
  ### UNIF
  id=sample(1:n,r*n)
  res[1]=(svd(Z[id,])$d[ncol(X)])^2/r
  
  ### MORE(A-opt)
  W=solve(t(Z)%*%Z)
  pi=abs(Y-p)*sqrt(apply((X%*%W)^2,1,sum))
  pi=pi/sum(pi)
  id=sample(1:n,r*n,FALSE,pi)
  res[2]=(svd(Z[id,])$d[ncol(X)])^2/r
  
  ### MORE(L-opt)
  pi=abs(Y-p)*sqrt(apply(X^2,1,sum))
  pi=pi/sum(pi)
  id=sample(1:n,r*n,FALSE,pi)
  res[3]=(svd(Z[id,])$d[ncol(X)])^2/r
  
  ### IBOSS
  c=X%*%beta
  f=c*c*(1-1/(1+exp(c)))^(ncol(X)+1)*(1/(1+exp(c)))^(ncol(X)+1)
  cstar=c[which.max(f)]
  idxB=rep(0,n)
  for(i in 1:n){
    if(min(abs(c[i]-cstar),abs(c[i]+cstar))<=2.5)  # 0.5 for T3
      idxB[i]=1
  }
  B=which(idxB==1)
  if(length(B)<=r*n){
    XB=X[B,]
    id=1:length(B)
  }
  else{
    XB=X[B,]
    k=round(r*n/2/ncol(X))
    t=rank(XB[,1],ties.method ="first")
    id=c(which(t<=k),which(t>=length(t)-k+1))
    for(j in 2:ncol(X)){
      idxb=c(1:nrow(XB))[-id]
      xb=XB[-id,]
      t=rank(xb[,j],ties.method ="first")
      idt=idxb[c(which(t<=k),which(t>=length(t)-k+1))]
      id=c(id,idt)
    }}
  res[4]=(svd(Z[id,])$d[ncol(X)])^2/r
    
  
  ### LEV
  H=(rr$u)%*%t(rr$u)
  pi=diag(H)
  t=rank(-pi,ties.method = "first")
  id=which(t<=r*n)
  res[5]=(svd(Z[id,])$d[ncol(X)])^2/r
  
  ### DKDO
  k=round(r*n/2/ncol(Z))
  t=rank(Z[,1],ties.method ="first")
  id=c(which(t<=k),which(t>=length(t)-k+1))
  for(j in 2:ncol(Z)){
    idz=c(1:n)[-id]
    z=Z[-id,]
    t=rank(z[,j],ties.method ="first")
    idt=idz[c(which(t<=k),which(t>=length(t)-k+1))]
    id=c(id,idt)
  }
  res[6]=(svd(Z[id,])$d[ncol(X)])^2/r
  
  return(res)
}

rseq=c(0.05,c(1:9)/10)
alphamain=function(n,beta,rseq,seed){
  ### generate the data
  set.seed(seed)
  D=GeneData(n,beta,"mzExp")
  X=D$X
  #X=log(X)
  Y=D$Y
  p=1-1/(1+exp(X%*%beta))
  w=p*(1-p)
  Z=c(sqrt(w))*X
  rr=svd(Z)
  sigmaz2=(rr$d[ncol(Z)])^2
  
  testalpha1=function(r){return(testalpha(X,Y,beta,p,Z,rr,r))}
  A=sapply(rseq,testalpha1)
  out=log(0.1*sigmaz2/A)/log(n)+1
  #out=log(0.5*sigmaz2/A)/log(n)+1

  return(out)
}

cl<- makeCluster(15) 
registerDoParallel(cl) 
ResultlExpalpha= foreach(i=1:500,
                      .combine=rbind,
                      .packages=c("mvtnorm")) %dopar% alphamain(1e4,rep(0.5,7),rseq,666*i+18)
stopCluster(cl)

AU=AA=AL=AI=Al=AD=matrix(0,500,10)
for(i in 1:500){
  t=ResultlExpalpha
  AU[i,]=t[6*i-5,]
  AA[i,]=t[6*i-4,]
  AL[i,]=t[6*i-3,]
  AI[i,]=t[6*i-2,]
  Al[i,]=t[6*i-1,]
  AD[i,]=t[6*i,]
}

ResultlExpA=data.frame(alpha=c(apply(AU,2,mean),apply(AA,2,mean),apply(AL,2,mean),
                             apply(AI,2,mean),apply(Al,2,mean),apply(AD,2,mean)),
                     r=rep(rseq,times=6),
                     method=rep(c("UNIF","MORE(A-opt)","MORE(L-opt)","IBOSS","LEV","DKDO"),each=length(rseq)))
write.csv(ResultlExpA,"Logistic/Main/Simulation/testalpha/ResultlExpA.csv")

# P=ggplot(ResultMzA,aes(x=r,y=alpha,group=method,colour=method))+
#   geom_line(aes(linetype=method))+
#   geom_point(aes(shape=method),size=2)+
#   scale_shape_manual(values=c(1,3,4,5,5,2))+
#   scale_linetype_manual(values=c(1,1,1,1,2,1))+
#   scale_color_manual(values=c("#4DAF4A","#A65628","#E41A1C","#984EA3","#984EA3","orange"))+
#   xlab("r")+ylab(expression(alpha))+theme_bw()+
#   scale_x_continuous(breaks = rseq[-1])+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y = element_text(angle=90),
#         legend.justification=c(1,0),
#         legend.position = c(1,0),
#         legend.background = element_rect(colour = "black"),
#         legend.key.size = unit(0.18, "inches"))
P=ggplot(ResultlExpA,aes(x=r,y=alpha,group=method,colour=method))+
  geom_line(aes(linetype=method))+
  geom_point(aes(shape=method),size=2)+
  scale_shape_manual(values=c(1,3,4,5,5,2))+
  scale_linetype_manual(values=c(1,1,1,1,2,1))+
  scale_color_manual(values=c("#4DAF4A","#A65628","#E41A1C","#984EA3","#984EA3","orange"))+
  xlab("r")+ylab(expression(alpha))+theme_bw()+
  scale_x_continuous(breaks = rseq[-1])+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(angle=90),
        legend.position ="none")

P




n=1e5
Rseq=c(c(1:5)*2/1000,c(1:4)*2/100)
Alphamain=function(n,beta,rseq,seed){
  ### generate the data
  set.seed(seed)
  D=GeneData(n,beta,"mzNormal")
  X=D$X
  #X=log(X)
  Y=D$Y
  p=1-1/(1+exp(X%*%beta))
  w=p*(1-p)
  Z=c(sqrt(w))*X
  rr=svd(Z)
  sigmaz2=(rr$d[ncol(Z)])^2
  
  testalpha1=function(r){return(testalpha(X,Y,beta,p,Z,rr,r))}
  A=sapply(rseq,testalpha1)
  out=log(0.5*sigmaz2/A)/log(n)+1
  
  return(out)
}

cl<- makeCluster(15) 
registerDoParallel(cl) 
ResultMzAlpha= foreach(i=1:500,
                         .combine=rbind,
                         .packages=c("mvtnorm")) %dopar% alphamain(1e5,rep(0.5,7),Rseq,666*i+18)
stopCluster(cl)

AU=AA=AL=AI=Al=AD=matrix(0,500,10)
for(i in 1:500){
  t=ResultlExpalpha
  AU[i,]=t[6*i-5,]
  AA[i,]=t[6*i-4,]
  AL[i,]=t[6*i-3,]
  AI[i,]=t[6*i-2,]
  Al[i,]=t[6*i-1,]
  AD[i,]=t[6*i,]
}

ResultlExpA=data.frame(alpha=c(apply(AU,2,mean),apply(AA,2,mean),apply(AL,2,mean),
                               apply(AI,2,mean),apply(Al,2,mean),apply(AD,2,mean)),
                       r=rep(rseq,times=6),
                       method=rep(c("UNIF","MORE(A-opt)","MORE(L-opt)","IBOSS","LEV","DKDO"),each=length(rseq)))
write.csv(ResultlExpA,"Logistic/Main/Simulation/testalpha/ResultlExpA.csv")

# P=ggplot(ResultMzA,aes(x=r,y=alpha,group=method,colour=method))+
#   geom_line(aes(linetype=method))+
#   geom_point(aes(shape=method),size=2)+
#   scale_shape_manual(values=c(1,3,4,5,5,2))+
#   scale_linetype_manual(values=c(1,1,1,1,2,1))+
#   scale_color_manual(values=c("#4DAF4A","#A65628","#E41A1C","#984EA3","#984EA3","orange"))+
#   xlab("r")+ylab(expression(alpha))+theme_bw()+
#   scale_x_continuous(breaks = rseq[-1])+
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.y = element_text(angle=90),
#         legend.justification=c(1,0),
#         legend.position = c(1,0),
#         legend.background = element_rect(colour = "black"),
#         legend.key.size = unit(0.18, "inches"))
P=ggplot(ResultlExpA,aes(x=r,y=alpha,group=method,colour=method))+
  geom_line(aes(linetype=method))+
  geom_point(aes(shape=method),size=2)+
  scale_shape_manual(values=c(1,3,4,5,5,2))+
  scale_linetype_manual(values=c(1,1,1,1,2,1))+
  scale_color_manual(values=c("#4DAF4A","#A65628","#E41A1C","#984EA3","#984EA3","orange"))+
  xlab("r")+ylab(expression(alpha))+theme_bw()+
  scale_x_continuous(breaks = rseq[-1])+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(angle=90),
        legend.position ="none")

P

