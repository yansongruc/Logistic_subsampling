################   Test the quantile of deterministic methods in nonuniform   ################
library(nabor)
library(mined)
library(mvtnorm)
library(OSMAC)
library(ggplot2)
library(foreach)
library(doParallel)
source("GeneData")
source("MED")
source("IBOSSL")
source("Leverage")
source("more_efficient")
source("test7_1")



####### Given full data ######
QD=function(n,beta,sub){
  set.seed(100)
  D=GeneData(n,beta,"mzNormal")
  X=D$X
  Y=D$Y
  
  # subsampling methods
  Ro=Rm=rep(0,300)
  for(i in 1:300){
    Ro[i]=sum((twostep(X,Y,200,sub-200,"mmse")$par-beta)^2)
    Rm[i]=sum((ME(X,Y,200,sub-200)-beta)^2)
  }
  Dmed=sum((MEDkw(X,Y,50,sub)$par-beta)^2)
  #Di=sum((IBOSSL(X,Y,200,sub,"Yes")$par-beta)^2) #only for T3
  Di=sum((IBOSSL(X,Y,200,sub,"No")$par-beta)^2)
  Dl=sum((LevDeter_Knn(X,Y,50,sub)$par-beta)^2)
  Dd=sum((InfDeter_Knn(X,Y,50,sub,"Dopt")$par-beta)^2)
  
  return(list(Osmac=Ro,More=Rm,Med=Dmed,IBOSS=Di,Lev=Dl,Dkdo=Dd))
}

ResultQD_mz=QD(1e4,rep(0.5,2),400)
DataQD_T3=data.frame(res=c(ResultQD_T3$More,ResultQD_T3$Osmac,
                           ResultQD_T3$Med,ResultQD_T3$IBOSS,ResultQD_T3$Lev,ResultQD_T3$Dkdo),
                     methods=factor(c(rep(c("MORE","OSMAC"),each=300),"MED","IBOSS","LEVERAGE","DKDO")))

#library(RColorBrewer)
#display.brewer.pal(8,'Set1')
#brewer.pal(8,'Set1') "#377EB8" "#984EA3"
PQD_T3=ggplot() +
  geom_density(aes(x=res,fill=methods),data=DataQD_T3[1:600,],alpha=0.3)+
  geom_vline(aes(xintercept=DataQD_T3[601,1]),linetype=1,color="#A65628")+
  geom_vline(aes(xintercept=DataQD_T3[602,1]),linetype=1,color="orange")+
  geom_vline(aes(xintercept=DataQD_T3[603,1]),linetype=1,color="#E41A1C")+
  geom_vline(aes(xintercept=DataQD_T3[604,1]),linetype=1,color="#4DAF4A")+
  theme(axis.text.y = element_text(angle=90),
        legend.justification=c(1,1),
        legend.position = c(1,1),
        legend.background = element_rect(colour = "black"))+
  xlab("MSE")+ylab("Empirical density")
PQD_T3

#mz- seed 100

#nz, mix, T3 -seed 1000







####### Conditional on population ######
QD=function(seed){
  set.seed(seed)
  n=1e4
  beta=c(0.5,0.5)
  sub=400
  #D=GeneData(n,beta,"Exp")
  X1=rnorm(n)
  X2=rexp(n,rate=2)
  X=cbind(X1,X2)
  X=D$X
  Y=D$Y
  
  # subsampling methods
  Ro=Rm=rep(0,500)
  for(i in 1:500){
    Ro[i]=sum((twostep(X,Y,200,sub-200,"mmse")$par-beta)^2)
    Rm[i]=sum((ME(X,Y,200,sub-200)$par-beta)^2)
  }
  Dmed=sum((MEDkw(X,Y,50,sub,"No")$par-beta)^2)
  Di=sum((IBOSSL(X,Y,200,sub-200,2.5,"No")$par-beta)^2) #only for T3
  #Di=sum((IBOSSL(X,Y,200,sub,"No")$par-beta)^2)
  Dl=sum((LevDeter_Knn(X,Y,50,sub)$par-beta)^2)
  Dd=sum((InfDeter_Knn(X,Y,50,sub,"Dopt")$par-beta)^2)
  
  return(list(Osmac=Ro,More=Rm,Med=Dmed,IBOSS=Di,Lev=Dl,Dkdo=Dd))
}

cl<- makeCluster(15) 
registerDoParallel(cl) 
ResultQD_me=foreach(i=1:500,
                   .combine=rbind,
                   .packages=c("mvtnorm","nabor","mined")) %dopar% QD(666*i+818)
stopCluster(cl)

ro=rm=matrix(0,500,500)
dm=di=dl=dd=rep(0,500)
for(i in 1:500){
  t=ResultQD_me[i,]
  ro[i,]=t$Osmac
  rm[i,]=t$More
  dm[i]=t$Med
  di[i]=t$IBOSS
  dl[i]=t$Lev
  dd[i]=t$Dkdo
}

DataQD_ME=data.frame(res=c(c(ro),c(rm),mean(dm),mean(di),mean(dl),mean(dd)),
                  methods=factor(c(rep(c("OSMAC","MORE"),each=25e4),"MED","IBOSS","LEVERAGE","DKDO")))

write.csv(DataQD_ME,"Logistic/Main/Figure1/DataQD_ME.csv")


# PQD_Mz=ggplot() +geom_boxplot(aes(x=methods,y=log(res),fill=methods),data=DataQD_Mz[1:5e5,],alpha=0.5)+
#   scale_fill_manual(values = c("#984EA3","#377EB8"),guide=FALSE)+theme_bw()+
#   geom_hline(aes(yintercept=log(res),colour="MED"),data=DataQD_Mz[5e5+1,],linetype=2)+
#   geom_hline(aes(yintercept=log(res),colour="IBOSS"),data=DataQD_Mz[5e5+2,],linetype=1)+
#   geom_hline(aes(yintercept=log(res),colour="LEV"),data=DataQD_Mz[5e5+3,],linetype=1)+
#   geom_hline(aes(yintercept=log(res),colour="DKDO"),data=DataQD_Mz[5e5+4,],linetype=1)+
#   scale_colour_manual(values=c("#4DAF4A","#A65628","#E41A1C","#FF7F00"))+
#   theme(axis.text.y = element_text(angle=90),
#         legend.title = element_blank(),
#         legend.justification=c(0.5,0),
#         legend.position = c(0.5,0),
#         legend.background = element_rect(colour = "black"))+
#   ylab("log(MSE)")
# PQD_Mz



PQD_ME=ggplot() +geom_boxplot(aes(x=methods,y=log(res),fill=methods),data=DataQD_ME[1:5e5,],alpha=0.5)+
  scale_fill_manual(values = c("#984EA3","#377EB8"),guide=FALSE)+theme_bw()+
  geom_hline(aes(yintercept=log(res),colour="MED"),data=DataQD_ME[5e5+1,],linetype=2)+
  geom_hline(aes(yintercept=log(res),colour="IBOSS"),data=DataQD_ME[5e5+2,],linetype=1)+
  geom_hline(aes(yintercept=log(res),colour="LEV"),data=DataQD_ME[5e5+3,],linetype=1)+
  geom_hline(aes(yintercept=log(res),colour="DKDO"),data=DataQD_ME[5e5+4,],linetype=1)+
  scale_colour_manual(values=c("#4DAF4A","#A65628","#E41A1C","#FF7F00"))+
  theme(axis.text.y = element_text(angle=90),
        legend.title = element_blank(),
        legend.position = "none")+
  ylab("log(MSE)")
PQD_ME







