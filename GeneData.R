##### generate various types of data
GeneData=function(n,beta,Type=c("mzNormal","nzNormal","ueNormal","mixNormal","T3","Exp","mzExp")){
  Sigma=0.5*matrix(1,length(beta),length(beta))+diag(0.5,length(beta))
  call=match.call()
  Type=match.arg(Type)
  if(Type=="mzNormal"){
    X=rmvnorm(n,rep(0,length(beta)),Sigma)
  }
  else if(Type=="nzNormal"){
    X=rmvnorm(n,rep(1.5,length(beta)),Sigma)
  }
  else if(Type=="ueNormal"){
    for(i in 1:length(beta)){
      Sigma[i,i]=i^2
    }
    X=rmvnorm(n,rep(0,length(beta)),Sigma)
  }
  else if(Type=="mixNormal"){
    X=0.5*rmvnorm(n,rep(1,length(beta)),Sigma)+0.5*rmvnorm(n,rep(-1,length(beta)),Sigma)
  }
  else if(Type=="T3"){
    X=0.1*rmvt(n,sigma=Sigma,df=3)
  }
  else if(Type=="Exp"){
    X=matrix(rexp(length(beta)*n,rate=2),nrow=n,ncol=length(beta))
  }
  else if(Type=="mzExp"){
    dd=round(length(beta)/2)
    if(dd==1){
      X1=rnorm(n)
    }
    else{
      X1=rmvnorm(n,rep(0,dd),Sigma[1:dd,1:dd])
    }
    X2=matrix(rexp((length(beta)-dd)*n,rate=2),nrow=n,ncol=length(beta)-dd)
    X=cbind(X1,X2)
  }
  P=1/(1+exp(-X%*%beta))
  Y=rep(NULL,n)
  for (i in 1:n) {
    Y[i]=rbinom(1,1,P[i])
  }
  return(list(X=X,Y=Y))
}


