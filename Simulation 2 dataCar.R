library(insuranceData)
library(dummies)
library(tableone)
library(tweedie)
library(HDtweedie)
library(MASS)
source("cv tweedie solpath.R")

data(dataCar)
Data<-dataCar
var<-c("veh_value","veh_body","veh_age","gender","area","agecat")
catvar<-c("veh_body","gender","area","agecat")

CreateTableOne(vars=var,factorVars=catvar,data=Data)

value<-Data$veh_value
body<-dummy(Data$veh_body)[,-1]
veh_age<-Data$veh_age
male<-ifelse(Data$gender=="M",1,0)
area<-dummy(Data$area)[,-1]
age<-dummy(Data$agecat)[,-1]

X<-cbind(value,body,veh_age,male,area,age)

beta.star<-round(glm(Data$clm~X,family=binomial)$coef,2)

size<-nrow(X)
dataX<-X  
set.seed(500)

SimCar <- function(sim) { 
#Index<-sort(sample(1:nrow(X),size=size,prob=rep(1/nrow(X),nrow(X))))

#dataX<-X[Index,]


logu<-cbind(1,dataX)%*%beta.star

Y<-rep(NA,nrow(dataX))
for (i in 1:nrow(dataX)){
  Y[i]<-rtweedie(1,xi=1.5,mu=exp(logu[i]),phi=1)
}

noise<-mvrnorm(nrow(dataX),rep(0,10),diag(1,10))
colnames(noise)<-paste("N",1:10,sep="")

DataX<-cbind(dataX,noise)
DataY<-Y

Data.trainX<-DataX[1:(size/2),]
Data.testX<-DataX[(size/2+1):size,]
  
trainY<-DataY[1:(size/2)]
testY<-DataY[(size/2+1):size]

beta.est<-matrix(NA,nrow=ncol(DataX),ncol=3)

group1<-rep(1,35)

best.lam<-CV.tweedie(X=Data.trainX,y=trainY,m=10,k=5,rho=1.5,alpha=1,g=group1,b0.int=beta.star[1],b.int=c(beta.star[-1],rep(0,10)))$best.lam


beta.est[,1]<-irlsbmd(X=Data.testX,y=testY,v=1/nrow(Data.testX),w=sqrt(group1),tau=1,rho=1.5,g=group1,lam=best.lam,b0=beta.star[1],b=c(beta.star[-1],rep(0,10)))$b
 
      
group2<-c(1,12,1,1,5,5,rep(1,10))

best.lam<-CV.tweedie(X=Data.trainX,y=trainY,m=10,k=5,rho=1.5,alpha=1,g=group2,b0.int=beta.star[1],b.int=c(beta.star[-1],rep(0,10)))$best.lam

      
beta.est[,2]<-irlsbmd(X=Data.testX,y=testY,v=1/nrow(Data.testX),w=sqrt(group2),tau=1,rho=1.5,g=group2,lam=best.lam,b0=beta.star[1],b=c(beta.star[-1],rep(0,10)))$b
 
 
group3 <- group2
Fit.meth3.lam<-NULL
Fit.meth3.dev<-NULL
alpha_v<-seq(0.1,1,0.2)
for (alpha in alpha_v){
  Fit.meth3<-CV.tweedie(X=Data.trainX,y=trainY,m=10,k=5,rho=1.5,alpha=alpha,g=group3,b0.int=beta.star[1],b.int=c(beta.star[-1],rep(0,10)))
  Fit.meth3.lam<-c(Fit.meth3.lam,Fit.meth3$best.lam)
  Fit.meth3.dev<-c(Fit.meth3.dev,Fit.meth3$min.Dev)
}
best.alpha<-alpha_v[which.min(Fit.meth3.dev)]
best.lam<-Fit.meth3.lam[which.min(Fit.meth3.dev)]  

beta.est[,3]<-irlsbmd(X=Data.testX,y=testY,v=1/nrow(Data.testX),w=sqrt(group3),tau=best.alpah,rho=1.5,g=group3,lam=best.lam,b0=beta.star[1],b=c(beta.star[-1],rep(0,10)))$b


beta.coef<-ifelse(beta.est==0,0,1)

C.coef<-apply(beta.coef[1:25,]==1,2,sum)
IC.coef<-apply(beta.coef[26:35,]==1,2,sum)

beta.block<-rbind(beta.coef[1,],apply(beta.coef[2:13,],2,sum),beta.coef[14:15,],apply(beta.coef[16:20,],2,sum),apply(beta.coef[21:25,],2,sum),beta.coef[26:35,])

C.block<-apply(beta.block[1:6,]>=1,2,sum)
IC.block<-apply(beta.block[7:16,]>=1,2,sum)

Results<-rbind(beta.est,C.coef,IC.coef,C.block,IC.block)
write.csv(Results,paste("SimCar_",sim,".csv",sep=""))
return(Results)
}

library(methods)

# parallel library is defautly installed by R, just need to load it
library(parallel)
n.sim<-50
res <- mclapply(seq(n.sim),SimCar,mc.cores = getOption("mc.cores", 50L))
#res <- mclapply(seq(n.sim),SimCar,mc.cores = 1)
Output<-Reduce("+", res)/length(res) 
Bias<-Output[1:25,]-beta.star[-1]

SE<-apply(array(unlist(res), c(39, 3, n.sim)), c(1,2), sd)[1:25,]

Var<-apply(array(unlist(res), c(39, 3, n.sim)), c(1,2), var)[1:25,]

rMSE<-sqrt(Bias^2+Var)

Out<-rbind(Output,Bias,SE,Var,rMSE)
write.csv(Out,"SimCar.csv")

