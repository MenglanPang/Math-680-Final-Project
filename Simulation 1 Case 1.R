###Simulation 1

source("cv tweedie solpath.R")


f1<-function(x)   x
f2<-function(x)  (3*x^2-1)/6
f3<-function(x)  (5*x^3-3*x)/10

library(MASS)
library(tweedie)

p<-8
n<-1000
rho<-1.5
set.seed(680)

Simulation_Case1 <- function(sim) { 
RESULTS<-NULL

for (omega in c(0,0.5)){

  Cov<-matrix(omega,nrow=p,ncol=p)
  diag(Cov)<-1
  
  X<-mvrnorm(n,rep(0,p),Cov) 
  
  B1<-f1(X)[,1:3]
  B2<-f2(X)[,1:3]
  B3<-f3(X)[,1:3]
  
  coef<-matrix(c(rep(0.5,3),rep(0.2,3),rep(0.5,3))*c(1,-1,1),ncol=1)
  
  logu<-0.3+cbind(B1,B2,B3)%*%coef
  
  
  Y<-rep(NA,n)
  for (i in 1:n){
    Y[i]<-rtweedie(1,xi=1.5,mu=exp(logu[i]),phi=1)
  }
  
  Data<-cbind(f1(X),f2(X),f3(X))
  
  Data.trainX<-Data[1:(n/2),as.vector(matrix(1:24,nrow=3,byrow=TRUE))]
  Data.testX<-Data[(n/2+1):n,as.vector(matrix(1:24,nrow=3,byrow=TRUE))]
  
  trainY<-Y[1:(n/2)]
  testY<-Y[(n/2+1):n]
  
  for (meth in 1:2){
    if (meth==1){
      alpha<-1
      g<-rep(1,24)
    }
    if (meth==2){
      alpha<-1
      g <-rep(3,8)
    }      
      
      best.lam<-CV.tweedie(X=Data.trainX,y=trainY,m=10,k=5,rho=rho,alpha=alpha,g=g,b0.int=0.3,b.int=c(0.5,0.2,0.5,-0.5,-0.2,-0.5,0.5,0.2,0.5,rep(0,15)))$best.lam

      test.Fit<-irlsbmd(X=Data.testX,y=testY,v=1/nrow(Data.testX),w=sqrt(g),tau=alpha,rho=rho,g=g,lam=best.lam,b0=0.3,b=c(0.5,0.2,0.5,-0.5,-0.2,-0.5,0.5,0.2,0.5,rep(0,15))) 
      
      Beta<-c(test.Fit$b0,test.Fit$b)
      
      beta.coef<-ifelse(test.Fit$b==0,0,1)
      C.coef<-sum(beta.coef[1:9]==1)
      IC.coef<-sum(beta.coef[10:24]==1)
      
      beta.block<-unlist(lapply(1:8,function(i) sum(beta.coef[(i-1)*3+1:3])))
      
      C.block<-sum(beta.block[1:3]>=1)
      IC.block<-sum(beta.block[4:8]>=1)
      
      v<-rep(1/nrow(Data.testX),nrow(Data.testX))
      Data.testX_1<-cbind(1,Data.testX)
      neg.loglik<-sum(v*(testY*exp(-(rho-1)*Data.testX_1%*%Beta)/(rho-1)+exp((2-rho)*Data.testX_1%*%Beta)/(2-rho)))
      
      Results<-c(C.block,IC.block,C.coef,IC.coef,neg.loglik)

    
    RESULTS<-rbind(RESULTS,Results)
  }
  
  meth<-3
  g <- rep(3,8)
  Fit.meth3.lam<-NULL
  Fit.meth3.dev<-NULL
  alpha_v<-seq(0.1,1,0.2)
  for (alpha in alpha_v){
    Fit.meth3<-CV.tweedie(X=Data.trainX,y=trainY,m=10,k=5,rho=rho,alpha=alpha,g=g,b0.int=0.3,b.int=c(0.5,0.2,0.5,-0.5,-0.2,-0.5,0.5,0.2,0.5,rep(0,15)))
    Fit.meth3.lam<-c(Fit.meth3.lam,Fit.meth3$best.lam)
    Fit.meth3.dev<-c(Fit.meth3.dev,Fit.meth3$min.Dev)
  }
  best.alpha<-alpha_v[which.min(Fit.meth3.dev)]
  best.lam<-Fit.meth3.lam[which.min(Fit.meth3.dev)]  

  test.Fit<-irlsbmd(X=Data.testX,y=testY,v=1/nrow(Data.testX),w=sqrt(g),tau=best.alpha,rho=rho,g=g,lam=best.lam,b0=0.3,b=c(0.5,0.2,0.5,-0.5,-0.2,-0.5,0.5,0.2,0.5,rep(0,15))) 
    
  Beta<-c(test.Fit$b0,test.Fit$b)
  
  beta.coef<-ifelse(test.Fit$b==0,0,1)
  C.coef<-sum(beta.coef[1:9]==1)
  IC.coef<-sum(beta.coef[10:24]==1)
  
  beta.block<-unlist(lapply(1:8,function(i) sum(beta.coef[(i-1)*3+1:3])))
  
  C.block<-sum(beta.block[1:3]>=1)
  IC.block<-sum(beta.block[4:8]>=1)
  
  v<-rep(1/nrow(Data.testX),nrow(Data.testX))
  Data.testX_1<-cbind(1,Data.testX)
  neg.loglik<-sum(v*(testY*exp(-(rho-1)*Data.testX_1%*%Beta)/(rho-1)+exp((2-rho)*Data.testX_1%*%Beta)/(2-rho)))
  
  Results<-c(C.block,IC.block,C.coef,IC.coef,neg.loglik)
    
  RESULTS<-rbind(RESULTS,Results)
   
}
write.csv(RESULTS,paste("Case1_",sim,".csv",sep=""))
return(RESULTS)
}

# must load method library in order to use mc.reset.stream
library(methods)

# parallel library is defautly installed by R, just need to load it
library(parallel)
n.sim<-50
res <- mclapply(seq(n.sim),Simulation_Case1,mc.cores = getOption("mc.cores", 50L))
Output<-Reduce("+", res)/length(res) 
write.csv(Output,"Sim_case1.csv")



