################   Real data: covertype   ################
library(nabor)
library(mined)
library(mvtnorm)
#library(MaxPro)
#library(OSMAC)
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
source("getMLE")
source("twostep")

### data treatment
covtype=read.csv("~/Logistic/Main/covtype/covtype.data", header=FALSE)
type1idx=which(covtype$V55==1)
type2idx=which(covtype$V55==2)  
type7idx=which(covtype$V55==7)   # 1 and 7 is for imbalance data

###### comparison for covertype 1 & 2
idx=c(type1idx,type2idx)
idx=sample(idx,length(idx),replace = FALSE)
Data=covtype[idx,]
Data$V55=Data$V55-1
X=Data[,1:10]
x1=x2=rep(0,length(idx))
x1[which(Data$V11==1)]=1
x1[which(Data$V12==1)]=2
x1[which(Data$V13==1)]=3
x1[which(Data$V14==1)]=4
X=cbind(X,x1)
for(i in 1:40){
  x2[which(Data[,14+i]==1)]=i
}
X=cbind(X,x2)
X=t(t(X)/apply(X, 2, sd))
Y=Data$V55
D=as.data.frame(cbind(X,Y))
a=glm(Y~.-1,family = binomial,data=D)
summary(a)

Xtrain=X[,-10]
Xtrain[,4]=log(Xtrain[,4]+1)
Xtrain[,5]=log(Xtrain[,5]-min(Xtrain[,5])+1)
Xtrain[,7]=log(Xtrain[,7])
#Xtrain[,8]=log(Xtrain[,8])
ida=which(Xtrain[,7]==-Inf)
Xtrain=Xtrain[-ida,]
Ytrain=Y[-ida]

D=as.data.frame(cbind(Xtrain,Ytrain))
a=glm(Ytrain~.-1,family = binomial,data=D)
summary(a)


###### comparison for covertypr 1 & 2
subseq=c(500,1000,1500,2000,2500,3000,4000)
Main1=function(Xtrain,Ytrain,subseq,seed){
  set.seed(seed)
  it=sample(1:nrow(Xtrain),4e5,replace=FALSE)
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
  parosm=matrix(NA,7,12)
  errosm=tosm=rep(NA,7)
  for(i in 1:7){
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
  parmor=matrix(NA,7,12)
  errmor=tmor=rep(NA,7)
  for(i in 1:7){
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
  erribo=tibo=rep(NA,7)
  paribo=matrix(NA,7,12)
  for(i in 1:7){
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
    idd=sample(1:4e5,sub*2,replace=FALSE)
    Xd=X[idd,]
    Yd=Y[idd]
    r=LevDeter_Pilot(cbind(1,Xd),Yd,200,sub-200)$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))
  }
  lev=sapply(subseq,Lev,simplify = FALSE)
  errlev=tlev=rep(NA,7)
  parlev=matrix(NA,7,12)
  for(i in 1:7){
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
    idd=sample(1:4e5,sub*2,replace=FALSE)
    #idd=sample(1:4e5,sub/5*9,replace=FALSE)
    Xd=X[idd,]
    Yd=Y[idd]
    r=DKDO_pilot(cbind(1,Xd),Yd,200,sub-200)$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))
  }
  dkd=sapply(subseq,Dkd,simplify = FALSE)
  errdkd=tdkd=rep(NA,7)
  pardkd=matrix(NA,7,12)
  for(i in 1:7){
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
  return(list(Fullpar=full,Fullerr=errfull,Fullt=tfull,Osmpar=parosm,Osmerr=errosm,Osmt=tosm,Morpar=parmor,Morerr=errmor,Mort=tmor,
              Ibopar=paribo,Iboerr=erribo,Ibot=tibo,Levpar=parlev,Leverr=errlev,Levt=tlev,Dkdpar=pardkd,Dkderr=errdkd,Dkdt=tdkd))
}

cl<- makeCluster(15) 
registerDoParallel(cl) 
Resultct=foreach(i=1:1000,
                .combine=cbind,
                .packages=c("mvtnorm","nabor","mined")) %dopar% Main1(Xtrain,Ytrain,subseq,666*i+18)
stopCluster(cl)

fullpar=matrix(NA,1000,12)
fullerr=fullt=rep(NA,1000)
for(i in 1:1000){
  fullpar[i,]=Resultct[,i]$Fullpar
  fullerr[i]=Resultct[,i]$Fullerr
  fullt[i]=Resultct[,i]$Fullt
}
Fullpar=apply(fullpar,2,mean)
Fullerr=mean(fullerr)
Fullt=mean(fullt)

osmpar=osmerr=osmt=morpar=morerr=mort=matrix(NA,1000,7)
ibopar=iboerr=ibot=levpar=leverr=levt=dkdpar=dkderr=dkdt=matrix(NA,1000,7)
for(i in 1:1000){
  t=Resultct[,i]
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
  
  for(j in 1:7){
    osmpar[i,j]=sum((t$Osmpar[j,]-Fullpar)^2)
    morpar[i,j]=sum((t$Morpar[j,]-Fullpar)^2)
    ibopar[i,j]=sum((t$Ibopar[j,]-Fullpar)^2)
    levpar[i,j]=sum((t$Levpar[j,]-Fullpar)^2)
    dkdpar[i,j]=sum((t$Dkdpar[j,]-Fullpar)^2)
  }
}

data_ct=data.frame(parmse=c(apply(osmpar,2,mean),apply(morpar,2,mean),
                             apply(ibopar,2,mean),apply(levpar,2,mean),apply(dkdpar,2,mean),
                             rep(mean(apply((fullpar-rep(1,1000)%*%t(Fullpar))^2,1,sum)),7)),
                    parsd=c(apply(osmpar,2,sd),apply(morpar,2,sd),
                            apply(ibopar,2,sd),apply(levpar,2,sd),apply(dkdpar,2,sd),
                            rep(sd(apply((fullpar-rep(1,1000)%*%t(Fullpar))^2,1,sum)),7)),
                    errmean=c(apply(osmerr,2,mean),apply(morerr,2,mean),
                              apply(iboerr,2,mean),apply(leverr,2,mean),apply(dkderr,2,mean),
                              rep(Fullerr,7)),
                    errsd=c(apply(osmerr,2,sd),apply(morerr,2,sd),
                            apply(iboerr,2,sd),apply(leverr,2,sd),apply(dkderr,2,sd),
                            rep(sd(fullerr),7)),
                    tmean=c(apply(osmt,2,mean),apply(mort,2,mean),
                            apply(ibot,2,mean),apply(levt,2,mean),apply(dkdt,2,mean),
                            rep(Fullt,7)),
                    Method=factor(rep(c("OSMAC","MORE","IBOSS","LEV","DKDO","FULL"),each=7)),
                    k=rep(c(500,1000,1500,2000,2500,3000,4000),times=6))
write.csv(data_ct,"data_ct.csv")

pd=position_dodge(10)
#brewer.pal(8,'Set1')
p_ctpar=ggplot(data_ct,aes(x=k,y=parmse,group=Method,colour=Method))+
  theme_bw()+theme(panel.border = element_rect(),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position = c(1,1),
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method),size=2)+scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_linetype_manual(values=c(1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("k")+ylab("MSE")
p_ctpar

p_cterr=ggplot(data_ct,aes(x=k,y=errmean,colour=Method))+theme_bw()+
  theme(panel.border = element_rect(),
        axis.text=element_text(size=12),
        axis.text.y = element_text(angle=90,size=12),
        axis.title=element_text(size=14),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.width=unit(2,"line"),
        legend.key.height=unit(1,"line"))+
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method),size=2)+scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_linetype_manual(values=c(1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("k")+ylab("PER")
p_cterr

p_ctt=ggplot(data_ct,aes(x=k,y=tmean,colour=Method))+theme_bw()+
  theme(panel.border = element_rect(),
        axis.text=element_text(size=12),
        axis.text.y = element_text(angle=90,size=12),
        axis.title=element_text(size=14),
        legend.justification=c(1,1),
        legend.position = c(1,1),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.width=unit(2,"line"),
        legend.key.height=unit(1,"line"))+
  geom_line(aes(linetype=Method))+
  geom_point(size=1)+scale_shape_manual(values=c(1,1,1,1,1,1))+
  scale_linetype_manual(values=c(1,1,1,1,1,1))+
  scale_size_manual(values=c(1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("computational time")+ylab("prediction error")
p_ctt


###### comparison for covertypr 1 & 7, which is imbalance data
idx=c(type1idx,type7idx)
idx=sample(idx,length(idx),replace = FALSE)
Data=covtype[idx,]
id7=which(Data$V55==7)
Data$V55[id7]=Data$V55[id7]-7
X=Data[,1:10]
x1=x2=rep(0,length(idx))
x1[which(Data$V11==1)]=1
x1[which(Data$V12==1)]=2
x1[which(Data$V13==1)]=3
x1[which(Data$V14==1)]=4
X=cbind(X,x1)
for(i in 1:40){
  x2[which(Data[,14+i]==1)]=i
}
X=cbind(X,x2)
X=t(t(X)/apply(X, 2, sd))
Y=Data$V55
D=as.data.frame(cbind(X,Y))
a=glm(Y~.,family = binomial,data=D)
summary(a) # all significant

Xtrain=X 
Xtrain[,4]=log(Xtrain[,4]+1)
Xtrain[,5]=log(Xtrain[,5]-min(Xtrain[,5])+1)
Xtrain[,7]=log(Xtrain[,7])
ida=which(Xtrain[,7]==-Inf)
Xtrain=Xtrain[-ida,]
Ytrain=Y[-ida]

D=as.data.frame(cbind(Xtrain,Ytrain))
a=glm(Ytrain~.,family = binomial,data=D)
summary(a)

###### comparison for covertype 1 & 7
subseq=c(500,1000,1500,2000,2500,3000,4000)
Main2=function(Xtrain,Ytrain,subseq,seed){
  set.seed(seed)
  it=sample(1:nrow(Xtrain),1.9e5,replace=FALSE)
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
  parosm=matrix(NA,7,13)
  errosm=tosm=rep(NA,7)
  for(i in 1:7){
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
  parmor=matrix(NA,7,13)
  errmor=tmor=rep(NA,7)
  for(i in 1:7){
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
  erribo=tibo=rep(NA,7)
  paribo=matrix(NA,7,13)
  for(i in 1:7){
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
    idd=sample(1:1.9e5,sub*3,replace=FALSE)
    Xd=X[idd,]
    Yd=Y[idd]
    r=LevDeter_Pilot(cbind(1,Xd),Yd,200,sub-200)$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))
  }
  lev=sapply(subseq,Lev,simplify = FALSE)
  errlev=tlev=rep(NA,7)
  parlev=matrix(NA,7,13)
  for(i in 1:7){
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
    idd=sample(1:1.9e5,sub*3,replace=FALSE)
    Xd=X[idd,]
    Yd=Y[idd]
    r=DKDO_pilot(cbind(1,Xd),Yd,200,sub-200)$par
    t=(proc.time()-t1)[[3]]
    return(list(res=r,time=t))
  }
  dkd=sapply(subseq,Dkd,simplify = FALSE)
  errdkd=tdkd=rep(NA,7)
  pardkd=matrix(NA,7,13)
  for(i in 1:7){
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
  
  return(list(Fullpar=full,Fullerr=errfull,Fullt=tfull,Osmpar=parosm,Osmerr=errosm,Osmt=tosm,Morpar=parmor,Morerr=errmor,Mort=tmor,
              Ibopar=paribo,Iboerr=erribo,Ibot=tibo,Levpar=parlev,Leverr=errlev,Levt=tlev,Dkdpar=pardkd,Dkderr=errdkd,Dkdt=tdkd))
}

cl<- makeCluster(15) 
registerDoParallel(cl) 
Resultcti=foreach(i=1:1000,
                 .combine=cbind,
                 .packages=c("mvtnorm","nabor","mined")) %dopar% Main2(Xtrain,Ytrain,subseq,666*i+18)
stopCluster(cl)

fullpar=matrix(NA,1000,13)
fullerr=fullt=rep(NA,1000)
for(i in 1:1000){
  fullpar[i,]=Resultcti[,i]$Fullpar
  fullerr[i]=Resultcti[,i]$Fullerr
  fullt[i]=Resultcti[,i]$Fullt
}
Fullpar=apply(fullpar,2,mean)
Fullerr=mean(fullerr)
Fullt=mean(fullt)

osmpar=osmerr=osmt=morpar=morerr=mort=matrix(NA,1000,7)
ibopar=iboerr=ibot=levpar=leverr=levt=dkdpar=dkderr=dkdt=matrix(NA,1000,7)
for(i in 1:1000){
  t=Resultcti[,i]
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
  
  for(j in 1:7){
    osmpar[i,j]=sum((t$Osmpar[j,]-Fullpar)^2)
    morpar[i,j]=sum((t$Morpar[j,]-Fullpar)^2)
    ibopar[i,j]=sum((t$Ibopar[j,]-Fullpar)^2)
    levpar[i,j]=sum((t$Levpar[j,]-Fullpar)^2)
    dkdpar[i,j]=sum((t$Dkdpar[j,]-Fullpar)^2)
  }
}

data_cti=data.frame(parmse=c(apply(osmpar,2,mean),apply(morpar,2,mean),
                             apply(ibopar,2,mean),apply(levpar,2,mean),apply(dkdpar,2,mean),
                            rep(mean(apply((fullpar-rep(1,1000)%*%t(Fullpar))^2,1,sum)),7)),
                   parsd=c(apply(osmpar,2,sd),apply(morpar,2,sd),
                           apply(ibopar,2,sd),apply(levpar,2,sd),apply(dkdpar,2,sd),
                           rep(sd(apply((fullpar-rep(1,1000)%*%t(Fullpar))^2,1,sum)),7)),
                   errmean=c(apply(osmerr,2,mean),apply(morerr,2,mean),
                             apply(iboerr,2,mean),apply(leverr,2,mean),apply(dkderr,2,mean),
                             rep(Fullerr,7)),
                   errsd=c(apply(osmerr,2,sd),apply(morerr,2,sd),
                           apply(iboerr,2,sd),apply(leverr,2,sd),apply(dkderr,2,sd),
                           rep(sd(fullerr),7)),
                   tmean=c(apply(osmt,2,mean),apply(mort,2,mean),
                           apply(ibot,2,mean),apply(levt,2,mean),apply(dkdt,2,mean),
                           rep(Fullt,7)),
                   Method=factor(rep(c("OSMAC","MORE","IBOSS","LEV","DKDO","FULL"),each=7)),
                   k=rep(c(500,1000,1500,2000,2500,3000,4000),times=6))
write.csv(data_cti,"data_cti.csv")

pd=position_dodge(10)
#brewer.pal(8,'Set1')
p_ctipar=ggplot(data_cti,aes(x=k,y=parmse,group=Method,colour=Method))+
  theme_bw()+theme(panel.border = element_rect(),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position = c(1,1),
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method),size=2)+scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_linetype_manual(values=c(1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("k")+ylab("MSE")
p_ctipar

p_ctierr=ggplot(data_cti,aes(x=k,y=errmean,colour=Method))+theme_bw()+
  theme(panel.border = element_rect(),
        axis.text=element_text(size=12),
        axis.text.y = element_text(angle=90,size=12),
        axis.title=element_text(size=14),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.width=unit(2,"line"),
        legend.key.height=unit(1,"line"))+
  geom_line(aes(linetype=Method))+
  geom_point(aes(shape=Method),size=2)+scale_shape_manual(values=c(1,2,3,4,5,6))+
  scale_linetype_manual(values=c(1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","orange","#A65628","#E41A1C","#984EA3","#377EB8"))+
  xlab("k")+ylab("PER")
p_ctierr

# p_ctit=ggplot(data_cti[-c(36:42),],aes(x=errmean,y=tmean,colour=Method))+theme_bw()+
#   theme(panel.border = element_rect(),
#         axis.text=element_text(size=12),
#         axis.text.y = element_text(angle=90,size=12),
#         axis.title=element_text(size=14),
#         legend.justification=c(1,1),
#         legend.position = c(1,1),
#         legend.title = element_blank(),
#         legend.background = element_rect(colour = "black"),
#         legend.key.width=unit(2,"line"),
#         legend.key.height=unit(1,"line"))+
#   geom_line(aes(linetype=Method))+
#   geom_vline(xintercept = data_su$errmean[29],linetype=2,color="orange")+
#   geom_hline(yintercept = data_su$tmean[29],linetype=2,color="orange")+
#   geom_point(size=1)+scale_shape_manual(values=c(1,1,1,1))+
#   scale_linetype_manual(values=c(1,1,1,1))+
#   scale_size_manual(values=c(1,1,1,1))+
#   scale_color_manual(values=c("#4DAF4A","#E41A1C","#984EA3","#377EB8"))+
#   xlab("prediction error")+ylab("computational time")
# p_ctit

