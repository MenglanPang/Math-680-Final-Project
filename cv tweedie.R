source("tweedieAlgo2.R")

CV.tweedie<-function(X,y,m,k,rho,alpha,g,b0.int,b.int){ 
  n<-nrow(X)
  p<-ncol(X)  
  n.block<-length(g)
  X.block<-lapply(seq(n.block),
                  function(i) X[,(sum(g[1:i])-g[i]+1):sum(g[1:i])])
  tau<-alpha
  
  initial<-calculatelam1(X,y,v=1/nrow(X),w=sqrt(g),tau,rho,g)
  
  lam1<-initial$lam1
  
  if(n>=p){
    lamm<-0.001*lam1
  }else{
    lamm<-0.05*lam1
  }
  lam.vec<-exp(rev(c(log(lamm),(seq(m-1)/(m-1))*(log(lam1)-log(lamm))+log(lamm))))
  
  #Create k equally size folds 
  folds <- cut(seq(1,nrow(X)),breaks=k,labels=FALSE)
  Dev<-NULL 
  
  #Perform k fold cross validation 
  for(i in 1:k)
  {   
    Indexes <- which(folds==i,arr.ind=TRUE)   
    testX <- X[Indexes, ]   
    trainX <- X[-Indexes, ] 
    testy<- y[Indexes] 
    trainy <- y[-Indexes]
    
    Fit<-solpath(X=trainX,y=trainy,v=1/nrow(trainX),w=sqrt(g),g=g,tau=alpha,rho=rho,lam.vec=lam.vec,b0.int=b0.int,b.int=b.int)
    Beta.hat.M<-rbind(Fit$b0,Fit$beta)
    
    v<-rep(1/nrow(testX),nrow(testX))
    testX_1<-cbind(1,testX)
    dev<-2*v*(testy*exp(-(rho-1)*testX_1%*%Beta.hat.M)/(rho-1)+exp((2-rho)*testX_1%*%Beta.hat.M)/(2-rho))
    
    Dev<-rbind(Dev,dev)   
  } 
  cv.Dev<-apply(Dev,2,sum) 
  best.lam<-lam.vec[which.min(cv.Dev)] 
  min.Dev<-min(cv.Dev)
  return(list(best.lam=best.lam,min.Dev=min.Dev)) 
}

# alpha<-1
# g<-rep(3,8)
# CV.fit<-CV.tweedie(Data.trainX,trainY,m=10,k=5,rho=rho,alpha=alpha,g=g,b0.int=0.3,b.int=c(0.5,0.2,0.5,-0.5,-0.2,-0.5,0.5,0.2,0.5,rep(0,15)))
#   
