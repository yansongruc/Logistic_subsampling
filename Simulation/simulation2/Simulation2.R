################   Simulation 2: MSE and prediction error v.s. n   ################
library(nabor)
library(mined)
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
# Study change of MSE virus n                 #
# Fix subsample size k=500                    #
# n=2.5e3, 5e3, 1e4, 1.5e4, 2e4               #
###############################################
Fixk=function(n,beta,Xt,Yt,seed){
  set.seed(seed)
  D=GeneData(n,beta,"mzExp")
  X=D$X
  #X=log(X)  # only for Exp
  Y=D$Y
  d=length(beta)
  
  res=matrix(0,6,d+1)
  res[1,2:(d+1)]=twostep(X,Y,200,300,"mmse")$par
  res[2,]=ME(cbind(1,X),Y,200,300)$par
  res[3,2:(d+1)]=IBOSSL(X,Y,200,500,2.5,"No")$par   # 0.5 for T3
  res[4,2:(d+1)]=LevDeter_Knn(X,Y,50,500)$par
  res[5,2:(d+1)]=DKDO_Knn(X,Y,50,500)$par
  res[6,2:(d+1)]=getMLE(X,Y,1)$par
  mseres=apply((res-rep(1,6)%*%t(c(0,beta)))^2,1,sum)
  p=1/(1+exp(-cbind(1,Xt)%*%t(res)))
  yhat=round(p)
  err=apply((yhat-Yt%*%t(rep(1,6)))^2,2,mean)
  res=cbind(res,mseres)
  res=cbind(res,err)
  return(res)
}

nseq=c(2.5e3,5e3,1e4,1.5e4,2e4)
beta=rep(0.5,7)
Fixk_n=function(seed){
  set.seed(seed+111)
  Dt=GeneData(1e4,beta,"mzExp")
  Xt=Dt$X
  #Xt=log(Xt)  # only for Exp
  Yt=Dt$Y
  
  Fixk1=function(n){return(Fixk(n,rep(0.5,7),Xt,Yt,seed))}
  out=sapply(nseq,Fixk1,simplify = FALSE)
  return(out)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultMEn=foreach(i=1:1000,
                  .combine=cbind,
                  .packages=c("mvtnorm","nabor","mined","OSMAC")) %dopar% Fixk_n(666*i+18)
stopCluster(cl)

osmmse=osmerr=mormse=morerr=ibomse=iboerr=levmse=leverr=dkdmse=dkderr=fulmse=fulerr=matrix(0,5,1000)
for(j in 1:1000){
  t=ResultMEn[,j]
  for(i in 1:5){
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
    fulmse[i,j]=t[[i]][6,9]
    fulerr[i,j]=t[[i]][6,10]
  }
}
mean1=function(x){
  return(mean(x,na.rm=TRUE))
}
osmmse=apply(osmmse,1,mean)
osmerr=apply(osmerr,1,mean)
mormse=apply(mormse,1,mean)
morerr=apply(morerr,1,mean)
ibomse=apply(ibomse,1,mean)
iboerr=apply(iboerr,1,mean)
levmse=apply(levmse,1,mean)
leverr=apply(leverr,1,mean)
dkdmse=apply(dkdmse,1,mean)
dkderr=apply(dkderr,1,mean)
fulmse=apply(fulmse,1,mean)
fulerr=apply(fulerr,1,mean)

MseMEn=data.frame(mse=c(osmmse,mormse,ibomse,levmse,dkdmse,fulmse),
                    err=c(osmerr,morerr,iboerr,leverr,dkderr,fulerr),
                    Methods=as.factor(rep(c("OSMAC","MORE","IBOSS","LEV","DKDO","FULL"),each=5)),
                    n=rep(nseq,times=6))
write.csv(MseMEn,"MseMEn.csv")

PMEnmse=ggplot(MseMEn,aes(x=log(n),y=log(mse),group=Methods,colour=Methods))+
  geom_line(linetype=1)+theme_bw()+
  geom_point(aes(shape=Methods),size=2)+
  scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("log(n)")+ylab("log(MSE)")+
  theme(axis.text.y = element_text(angle=90),
        legend.justification=c(1,1),
        legend.position = "none",
        legend.background = element_rect(colour = "black"))
PMEnmse

PMEnerr=ggplot(MseMEn,aes(x=log(n),y=err,group=Methods,colour=Methods))+
  geom_line(linetype=1)+theme_bw()+
  geom_point(aes(shape=Methods),size=2)+
  scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("log(n)")+ylab("Prediction error")+
  theme(axis.text.y = element_text(angle=90),
        legend.justification=c(1,1),
        legend.position ="none",
        legend.background = element_rect(colour = "black"))
PMEnerr


