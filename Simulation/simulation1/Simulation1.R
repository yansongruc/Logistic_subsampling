################   Simulation 1: MSE and prediction error v.s. k   ################
library(nabor)
#library(mined)
library(mvtnorm)
library(OSMAC)
library(ggplot2)
library(foreach)
library(doParallel)
source("GeneData")
source("more_efficient")
#source("MED")
source("IBOSSL")
source("Leverage")
source("DKDO")

###############################################
# Study change of MSE virus k                 #
# Fix full sample size n=1e4                  #
# k=300,400,600,800,1000,1200                 #
###############################################
Fixn=function(X,Y,Xt,Yt,sub){
  d=ncol(X)
  res=matrix(0,5,d+1)
  aa=twostep(X,Y,200,sub-200,"mmse")
  if(is.na(aa[1])){
    res[1,2:(d+1)]=rep(NA,d)
  }
  else{
    res[1,2:(d+1)]=aa$par
  }
  res[2,]=ME(cbind(1,X),Y,200,sub-200)$par
  res[3,2:(d+1)]=IBOSSL(X,Y,200,sub,2.5,"No")$par   # 0.5 for T3
  res[4,2:(d+1)]=LevDeter_Knn(X,Y,50,sub)$par
  res[5,2:(d+1)]=DKDO_Knn(X,Y,50,sub)$par
  mseres=apply((res-rep(1,5)%*%t(c(0,beta)))^2,1,sum)
  p=1/(1+exp(-cbind(1,Xt)%*%t(res)))
  yhat=round(p)
  err=apply((yhat-Yt%*%t(rep(1,5)))^2,2,mean)
  res=cbind(res,mseres)
  res=cbind(res,err)
  return(res)
}

n=1e4
beta=rep(0.5,7)
subseq=c(3e2,4e2,6e2,8e2,1e3,1.2e3)
Fixn_k=function(seed){
  set.seed(seed)
  D=GeneData(n,beta,"mzExp")
  X=D$X
  #X=log(X)    # only for Exp
  Y=D$Y
  Dt=GeneData(n,beta,"mzExp")
  Xt=Dt$X
  #Xt=log(Xt)    # only for Exp
  Yt=Dt$Y
  
  fullres=getMLE(X,Y,1)$par
  fullmse=sum((fullres-beta)^2)
  fullp=1/(1+exp(-cbind(1,Xt)%*%c(0,fullres)))
  fully=round(fullp)
  fullerr=mean((fully-Yt)^2)
  Fixn1=function(sub){return(Fixn(X,Y,Xt,Yt,sub))}
  out=sapply(subseq,Fixn1,simplify = FALSE)
  return(list(full=c(fullres,fullmse,fullerr),sub=out))
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultMEk=foreach(i=1:1000,
                   .combine=cbind,
                   .packages=c("mvtnorm","nabor","mined","OSMAC")) %dopar% Fixn_k(666*i+18)
stopCluster(cl)

fullmse=fullerr=rep(0,1000)
for(i in 1:1000){
  fullmse[i]=ResultMEk[,i]$full[8]
  fullerr[i]=ResultMEk[,i]$full[9]
}
fullmse=mean(fullmse,na.rm=TRUE)
fullerr=mean(fullerr,na.rm=TRUE)

osmmse=osmerr=mormse=morerr=ibomse=iboerr=levmse=leverr=dkdmse=dkderr=matrix(0,6,1000)
for(j in 1:1000){
  t=ResultMEk[,j]$sub
  for(i in 1:6){
    osmmse[i,j]=t[[i]][1,9]
    osmerr[i,j]=t[[i]][1,10]
    mormse[i,j]=t[[i]][2,9]
    morerr[i,j]=t[[i]][2,10]
    ibomse[i,j]=t[[i]][3,9]
    iboerr[i,j]=t[[i]][3,10]
    levmse[i,j]=t[[i]][4,9]
    leverr[i,j]=t[[i]][4,10]
    dkdmse[i,j]=t[[i]][5,9]
    dkderr[i,j]=t[[i]][5,10]
  }
}
mean1=function(x){
  return(mean(x,na.rm=TRUE))
}
osmmse=apply(osmmse,1,mean1)
osmerr=apply(osmerr,1,mean1)
mormse=apply(mormse,1,mean1)
morerr=apply(morerr,1,mean1)
ibomse=apply(ibomse,1,mean1)
iboerr=apply(iboerr,1,mean1)
levmse=apply(levmse,1,mean1)
leverr=apply(leverr,1,mean1)
dkdmse=apply(dkdmse,1,mean1)
dkderr=apply(dkderr,1,mean1)

MseMEk=data.frame(mse=c(osmmse,mormse,ibomse,levmse,dkdmse,rep(fullmse,6)),
                  err=c(osmerr,morerr,iboerr,leverr,dkderr,rep(fullerr,6)),
                  Methods=as.factor(rep(c("OSMAC","MORE","IBOSS","LEV","DKDO","FULL"),each=6)),
                  k=rep(subseq,times=6))
write.csv(MseMEk,"MseMEk.csv")

PMEkmse=ggplot(MseMEk,aes(x=k,y=mse,group=Methods,colour=Methods))+
  geom_line(linetype=1)+theme_bw()+
  geom_point(aes(shape=Methods),size=2)+
  scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("k")+ylab("MSE")+
  theme(axis.text.y = element_text(angle=90),
        legend.justification=c(1,1),
        legend.position = c(1,1),
        legend.background = element_rect(colour = "black"))
PMEkmse

PMEkerr=ggplot(MseMEk,aes(x=k,y=err,group=Methods,colour=Methods))+
  geom_line(linetype=1)+theme_bw()+
  geom_point(aes(shape=Methods),size=2)+
  scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("k")+ylab("Prediction error")+
  theme(axis.text.y = element_text(angle=90),
        legend.justification=c(1,1),
        legend.position = c(1,1),
        legend.background = element_rect(colour = "black"))
PMEkerr



