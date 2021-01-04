################   Real data: phoneme   ################
library(nabor)
library(mined)
library(mvtnorm)
library(OSMAC)
library(ggplot2)
library(foreach)
library(doParallel)
library(readr)
library(splines)
source("more_efficient")
#source("MED")
source("IBOSSL")
source("Leverage")
source("DKDO")

### data treatment
phoneme=read_csv("Logistic/Main/phoneme/phoneme.data")
Y=phoneme$speaker
X=phoneme[,-c(1,258)]
idaa=which(Y=="aa")
iddcl=which(Y=="dcl")
idsh=which(Y=="sh")
idiy=which(Y=="iy")
#idao=which(Y=="ao")
DataX=X[c(idaa,iddcl,idsh,idiy),]
DataY=Y[c(idaa,iddcl,idsh,idiy)]

idi=c(1:2,696,697,2325:2327,1453:1455)
Dataidi=data.frame(freq=rep(c(1:256),each=10),
                   logp=c(as.matrix(DataX[idi,])),
                   class=factor(rep(c("aa","aa","dcl","dcl","iy","iy","iy","sh","sh","sh"),times=256)))
Pidi=ggplot(Dataidi,aes(x=freq,y=logp,colour=class))+
  geom_line(aes(linetype=class))+
  scale_linetype_manual(values=c(1,2,1,2))+
  scale_color_manual(values=c("#4DAF4A","#4DAF4A","#E41A1C","#E41A1C"))+
  xlab("Frequency")+ylab("log-periodogram")+theme_bw()+
  theme(axis.text.y = element_text(angle=90),
        legend.justification=c(1,1),
        legend.position = c(1,1),
        legend.background = element_rect(colour = "black"))

Pidi

H=ns(1:256,df=8)
HX=t(t(H)%*%t(DataX))
HY=c(rep(0,length(idaa)+length(iddcl)),rep(1,length(idsh)+length(idiy)))
id=sample(1:nrow(HX),nrow(HX),replace = FALSE)
Xtrain=HX[id,]
Xtrain=t(t(Xtrain) / apply(Xtrain, 2, sd))
Ytrain=HY[id]

Dtrain=as.data.frame(cbind(Xtrain,Ytrain))
a=glm(Ytrain~.,family = binomial,data=Dtrain)
summary(a)
# 6 is not 
Xtrain=Xtrain[,-6]
Dtrain=as.data.frame(cbind(Xtrain,Ytrain))
a=glm(Ytrain~.,family = binomial,data=Dtrain)
summary(a)


###### comparison
subseq=c(300,500,700,1000,1200,1500)
Main=function(Xtrain,Ytrain,subseq,seed){
  set.seed(seed)
  it=sample(1:nrow(Xtrain),2500)
  X=Xtrain[it,]
  Y=Ytrain[it]
  Xt=Xtrain[-it,]
  Yt=Ytrain[-it]
  
  ### full data
  t1=proc.time()
  full=getMLE(cbind(1,X),Y,1)$par
  tfull=(proc.time()-t1)[[3]]
  pfull=1/(1+exp(-cbind(1,Xt)%*%full))
  yfull=round(pfull)
  errfull=mean((Yt-yfull)^2)
  
  ### OSMAC
  Osm=function(sub){
    t1=proc.time()
    r=twostep(cbind(1,X),Y,200,sub-200,"mmse")$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))}
  osm=sapply(subseq,Osm,simplify = FALSE)
  errosm=tosm=rep(NA,6)
  parosm=matrix(NA,6,8)
  for(i in 1:6){
    if(is.null(osm[[i]]$res[1])){
    }
    else{
      parosm[i,]=osm[[i]]$res
      posm=1/(1+exp(-cbind(1,Xt)%*%osm[[i]]$res))
      yosm=round(posm)
      errosm[i]=mean((Yt-yosm)^2)
      tosm[i]=osm[[i]]$time
    }
  }
  
  ### MORE
  Mor=function(sub){
    t1=proc.time()
    r=ME(cbind(1,X),Y,200,sub-200)$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))}
  mor=sapply(subseq,Mor,simplify = FALSE)
  errmor=tmor=rep(NA,6)
  parmor=matrix(NA,6,8)
  for(i in 1:6){
    if(is.na(mor[[i]]$res[1])){
    }
    else{
      parmor[i,]=mor[[i]]$res
      pmor=1/(1+exp(-cbind(1,Xt)%*%mor[[i]]$res))
      ymor=round(pmor)
      errmor[i]=mean((Yt-ymor)^2)
      tmor[i]=mor[[i]]$time
    }
  }
  
  ### IBOSS
  Ibo=function(sub){
    t1=proc.time()
    r=IBOSSL(X,Y,200,sub-200,2.5,"Yes")$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))}
  ibo=sapply(subseq,Ibo,simplify = FALSE)
  erribo=tibo=rep(NA,6)
  paribo=matrix(NA,6,8)
  for(i in 1:6){
    if(is.na(ibo[[i]]$res[1])){
    }
    else{
      paribo[i,]=ibo[[i]]$res
      pibo=1/(1+exp(-cbind(1,Xt)%*%ibo[[i]]$res))
      yibo=round(pibo)
      erribo[i]=mean((Yt-yibo)^2)
      tibo[i]=ibo[[i]]$time
    }
  }
  
  ### Leverage: Pilot
  Lev=function(sub){
    t1=proc.time()
    r=LevDeter_Pilot(cbind(1,X),Y,200,sub-200)$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))}
  lev=sapply(subseq,Lev,simplify = FALSE)
  errlev=tlev=rep(NA,6)
  parlev=matrix(NA,6,8)
  for(i in 1:6){
    if(is.na(lev[[i]]$res[1])){
    }
    else{
      parlev[i,]=lev[[i]]$res
      plev=1/(1+exp(-cbind(1,Xt)%*%lev[[i]]$res))
      ylev=round(plev)
      errlev[i]=mean((Yt-ylev)^2)
      tlev[i]=lev[[i]]$time
    }
  }

  ### DKDO: pilot
  Dkd=function(sub){
    t1=proc.time()
    r=DKDO_pilot(cbind(1,X),Y,200,sub-200)$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))}
  dkd=sapply(subseq,Dkd,simplify = FALSE)
  errdkd=tdkd=rep(NA,6)
  pardkd=matrix(NA,6,8)
  for(i in 1:6){
    if(is.na(dkd[[i]]$res[1])){
    }
    else{
      pardkd[i,]=dkd[[i]]$res
      pdkd=1/(1+exp(-cbind(1,Xt)%*%dkd[[i]]$res))
      ydkd=round(pdkd)
      errdkd[i]=mean((Yt-ydkd)^2)
      tdkd[i]=dkd[[i]]$time
    }
  }

  return(list(Fullpar=full,Fullerr=errfull,Fullt=tfull,Osmpar=parosm,Osmerr=errosm,Osmt=tosm,
              Morpar=parmor,Morerr=errmor,Mort=tmor,Ibopar=paribo,Iboerr=erribo,Ibot=tibo,
              Levpar=parlev,Leverr=errlev,Levt=tlev,Dkdpar=pardkd,Dkderr=errdkd,Dkdt=tdkd))
}

cl<- makeCluster(15) 
registerDoParallel(cl) 
Resultpm=foreach(i=1:1000,
                  .combine=cbind,
                  .packages=c("mvtnorm","nabor","mined")) %dopar% Main(Xtrain,Ytrain,subseq,666*i+18)
stopCluster(cl)

fullpar=matrix(0,1000,8)
fullerr=fullt=rep(NULL,1000)
for(i in 1:1000){
  fullpar[i,]=Resultpm[,i]$Fullpar
  fullerr[i]=Resultpm[,i]$Fullerr
  fullt[i]=Resultpm[,i]$Fullt
}
Fullpar=apply(fullpar,2,mean)
#[1] -12.601923  16.442459  -3.021905   4.867152  -3.910880   2.733905  -9.483273   6.170050
Fullerr=mean(fullerr)
#[1] 0.02580547

osmpar=osmerr=osmt=morpar=morerr=osmt=ibopar=iboerr=ibot=matrix(NA,1000,6)
levpar=leverr=levt=dkdpar=dkderr=dkdpar=matrix(NA,1000,6)
for(i in 1:1000){
  t=Resultpm[,i]
  osmerr[i,]=t$Osmerr
  morerr[i,]=t$Morerr
  iboerr[i,]=t$Iboerr
  leverr[i,]=t$Leverr
  dkderr[i,]=t$Dkderr
  osmt[i,]=t$Osmt
  mort[i,]=t$Mort
  ibot[i,]=t$Ibot
  levt[i,]=t$Levt
  dkdt[i,]=t$Dkdt
  
  for(j in 1:6){
    if(is.na(t$Osmt[j])){
    }
    else{
      osmpar[i,j]=sum((t$Osmpar[j,]-Fullpar)^2)
    }
    
    if(is.na(t$Mort[j])){
    }
    else{
      morpar[i,j]=sum((t$Morpar[j,]-Fullpar)^2)
    }
    
    if(is.na(t$Ibot[j])){
    }
    else{
      ibopar[i,j]=sum((t$Ibopar[j,]-Fullpar)^2)
    }
    
    if(is.na(t$Levt[j])){
    }
    else{
      levpar[i,j]=sum((t$Levpar[j,]-Fullpar)^2)
    }
    
    if(is.na(t$Dkdt[j])){
    }
    else{
      dkdpar[i,j]=sum((t$Dkdpar[j,]-Fullpar)^2)
    }
  }
}

mean1=function(x){return(mean(x,na.rm=TRUE))}
sd1=function(x){return(sd(x,na.rm=TRUE))}

data_pm=data.frame(nNA=c(apply(apply(osmpar,2,is.na),2,mean),apply(apply(morpar,2,is.na),2,mean),apply(apply(ibopar,2,is.na),2,mean),
                         apply(apply(levpar,2,is.na),2,mean),apply(apply(dkdpar,2,is.na),2,mean),rep(0,6)),
                   parmse=c(apply(osmpar,2,mean1),apply(morpar,2,mean1),apply(ibopar,2,mean1),
                            apply(levpar,2,mean1),apply(dkdpar,2,mean1),rep(mean(apply((fullpar-rep(1,1000)%*%t(Fullpar))^2,1,sum)),6)),
                   parsd=c(apply(osmpar,2,sd1),apply(morpar,2,sd1),apply(ibopar,2,sd1),
                           apply(levpar,2,sd1),apply(dkdpar,2,sd1),rep(sd(apply((fullpar-rep(1,1000)%*%t(Fullpar))^2,1,sum)),6)),
                   errmean=c(apply(osmerr,2,mean1),apply(morerr,2,mean1),apply(iboerr,2,mean1),
                             apply(leverr,2,mean1),apply(dkderr,2,mean1),rep(Fullerr,6)),
                   errsd=c(apply(osmerr,2,sd1),apply(morerr,2,sd1),apply(iboerr,2,sd1),
                           apply(leverr,2,sd1),apply(dkderr,2,sd1),rep(sd(fullerr),6)),
                   tmean=c(apply(osmt,2,mean1),apply(mort,2,mean1),apply(ibot,2,mean1),
                           apply(levt,2,mean1),apply(dkdt,2,mean1),rep(Fullt,6)),
                   Method=factor(rep(c("OSMAC","MORE","IBOSS","LEV","DKDO","FULL"),each=6)),
                   k=rep(c(300,500,700,1000,1200,1500),times=6))
write.csv(data_pm,"Logistic/Main/phoneme/data_pm.csv")

idd=c(7:18)
p_pmpar=ggplot(data_pm[-idd,],aes(x=k,y=log(parmse),group=Method,colour=Method))+
  theme_bw()+theme(axis.text.y = element_text(angle=90),
                   legend.justification=c(1,1),
                   legend.position = c(1,1),
                   legend.background = element_rect(colour = "black"))+
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method),size=2)+scale_shape_manual(values=c(1,2,4,6))+
  scale_linetype_manual(values=c(1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","orange","#E41A1C","#377EB8"))+
  xlab("k")+ylab("log(MSE)")
p_pmpar

p_pmerr=ggplot(data_pm,aes(x=k,y=errmean,colour=Method))+
  theme_bw()+theme(axis.text.y = element_text(angle=90),
                   legend.justification=c(1,1),
                   legend.position = c(1,1),
                   legend.background = element_rect(colour = "black"))+
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method),size=2)+scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_linetype_manual(values=c(1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("k")+ylab("PER")
p_pmerr


p_pmpar=ggplot(data_pm[-idd,],aes(x=k,y=parmse,group=Method,colour=Method))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_rect(colour = "black"),
                      axis.text=element_text(size=12),
                      axis.text.y = element_text(angle=90,size=12),
                      axis.title=element_text(size=14),
                      legend.justification=c(1,1),
                      legend.position = c(1,1),
                      legend.title = element_blank(),
                      legend.background = element_rect(colour = "black"),
                      legend.key.width=unit(2,"line"),
                      legend.key.height=unit(1,"line"))+
  geom_line(aes(linetype=Method), position=pd)+
  geom_point(position=pd,size=1)+scale_shape_manual(values=c(1,1,1,1,1,1,1,1,1,1))+
  scale_linetype_manual(values=c(1,2,3,1,1,1,2,3,1,1))+
  scale_size_manual(values=c(1,1,1,1,1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","#4DAF4A","#4DAF4A","black","#A65628","#E41A1C","#E41A1C","#E41A1C","#377EB8","#984EA3"))+
  xlab("k")+ylab("MSE")
p_pmpar

p_pmerr=ggplot(data_pm,aes(x=k,y=errmean,colour=Method))+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text=element_text(size=12),
        axis.text.y = element_text(angle=90,size=12),
        axis.title=element_text(size=14),
        legend.justification=c(1,1),
        legend.position = c(1,1),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.width=unit(2,"line"),
        legend.key.height=unit(1,"line"))+
  geom_line(aes(linetype=Method), position=pd)+
  geom_point(position=pd,size=1)+scale_shape_manual(values=c(1,1,1,1,1,1,1,1,1,1))+
  scale_linetype_manual(values=c(1,2,3,1,1,1,2,3,1,1))+
  scale_size_manual(values=c(1,1,1,1,1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","#4DAF4A","#4DAF4A","black","#A65628","#E41A1C","#E41A1C","#E41A1C","#377EB8","#984EA3"))+
  xlab("k")+ylab("MSE")
p_pmerr

