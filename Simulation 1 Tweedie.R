###Simulation Case 1
library(MASS)
set.seed(680)
w<-0.5
p<-8
n<-500
Cov<-matrix(w,nrow=p,ncol=p)
diag(Cov)<-1
X<-mvrnorm(n,rep(0,p),Cov) 


f1<-function(x)   x
f2<-function(x)  (3*x^2-1)/6
f3<-function(x)  (5*x^3-3*x)/10

B1<-f1(X)[,1:3]
B2<-f2(X)[,1:3]
B3<-f3(X)[,1:3]


coef<-matrix(c(rep(0.5,3),rep(0.2,3),rep(0.5,3))*c(1,-1,1),ncol=1)

logu<-0.3+cbind(B1,B2,B3)%*%coef

ro<-1.5
phi<-1

Y<-rep(NA,n)
for (i in 1:n){
Y[i]<-rtweedie(1,xi=1.5,mu=exp(logu[i]),phi=1)
}

Data_case1<-cbind(f1(X),f2(X),f3(X))

Data_case1.b<-Data_case1[,as.vector(matrix(1:24,nrow=3,byrow=TRUE))]
  
###Simulation Case 2
set.seed(680)
w<-0.5
p<-6
n<-500
Cov<-matrix(w,nrow=p,ncol=p)
diag(Cov)<-1
Z<-mvrnorm(n,rep(0,p),Cov)

X<-cbind(Z[,1]+matrix(rnorm(n*3,0,0.1),nrow=n,ncol=3),Z[,2:6])

B1<-f1(X)[,1:3]
B2<-f2(X)[,1:3]
B3<-f3(X)[,1:3]

coef<-matrix(c(rep(0.5,3),rep(0.2,3),rep(0.5,3))*c(1,-1,1),ncol=1)

logu<-0.3+cbind(B1,B2,B3)%*%coef

ro<-1.5
phi<-1

Y<-rep(NA,n)
for (i in 1:n){
  Y[i]<-rtweedie(1,xi=1.5,mu=exp(logu[i]),phi=1)
}

Data_case2<-cbind(f1(X),f2(X),f3(X))

Data_case2.b<-Data_case2[,as.vector(matrix(1:24,nrow=3,byrow=TRUE))]

###Simulation Case 3
set.seed(680)
w<-0.5
p<-8
n<-500
Cov<-matrix(w,nrow=p,ncol=p)
diag(Cov)<-1
X<-mvrnorm(n,rep(0,p),Cov) 

B1_6<-f1(X)[,1:6]

coef<-matrix(rep(c(1,-1),3),ncol=1)

logu<-0.3+B1_6%*%coef

ro<-1.5
phi<-1

Y<-rep(NA,n)
for (i in 1:n){
  Y[i]<-rtweedie(1,xi=1.5,mu=exp(logu[i]),phi=1)
}

Data_case3.b<-Data_case3[,as.vector(matrix(1:24,nrow=3,byrow=TRUE))]




