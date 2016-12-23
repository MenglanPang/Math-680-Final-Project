source('tweedieAlgo1.R')
# This function calculate lambda1 for a given data set
calculatelam1<-function(X,y,v,w,tau,rho,g,tol=1e-8){
  n<-nrow(X)
  p<-ncol(X)
  n.block<-length(g)
  X.block<-lapply(seq(n.block),function(i) X[,(sum(g[1:i])-g[i]+1):sum(g[1:i])])
  
  b0<-0
  iterating=TRUE
  while(iterating){
    coeff.v<-(rho-1)*y*exp(-(rho-1)*b0)+
      (2-rho)*exp((2-rho)*b0)
    v.update<-v*coeff.v
    y.update<-b0+v/v.update*
      (y*exp(-(rho-1)*b0)-exp((2-rho)*b0))
    U0<-sum(-v.update*(y.update-b0))
    gamma0<-sum(v.update)
    b0.new<-b0-U0/gamma0
    if((abs(b0.new-b0)<tol)){
      iterating=FALSE
    }
    b0<-b0.new
  }
  U<-lapply(X.block,function(x) {
    x<-as.matrix(x)
    apply((-y*exp(-(rho-1)*b0)+exp((2-rho)*b0))*v*x,2,sum) }
            )
  
  U.lam<-sapply(seq(n.block),function(i) sqrt(crossprod(U[[i]]))/(tau*w[i]))
  lam1<-max(U.lam)
  return(list(lam1=lam1,b0=b0,U.lam=U.lam))
}

# Computing the solution path
solpath<-function(X,y,m=20,v,w,tau,rho,g,lam.vec=NULL,tol=1e-8,maxit=1e20,maxit.in=3e8,b0.int,b.int){
  n<-nrow(X)
  p<-ncol(X)
  n.block<-length(g)
  X.block<-lapply(seq(n.block),
                  function(i) X[,(sum(g[1:i])-g[i]+1):sum(g[1:i])])

  if (is.null(lam.vec)==0){
    initial<-irlsbmd(X,y,b0.int,b.int,v,w,tau,rho,g,
                     lam.vec[1],tol,maxit,maxit.in)
    m<-length(lam.vec)
    b0.vec<-rep(NA,m)
    b.vec<-matrix(NA,p,m)
    b0.vec[1]<-initial$b0
    b.vec[,1]<-initial$b
  }else{
    initial<-calculatelam1(X,y,v,w,tau,rho,g)
    lam1<-initial$lam1
    if(n>=p){
      lamm<-0.001*lam1
    }else{
      lamm<-0.05*lam1
    }
    lam.vec<-exp(rev(c(log(lamm),(seq(m-1)/(m-1))*(log(lam1)-log(lamm))+log(lamm))))
    b0.vec<-rep(NA,m)
    b.vec<-matrix(NA,p,m)
    b0.vec[1]<-initial$b0
    b.vec[,1]<-rep(0,p)
  }
  
  U.lam<-initial$U.lam
  Sc<-NULL
  for(i in 2:m){
    b0.initial<-b0.vec[i-1]
    b.initial<-b.vec[,i-1]
    b.block<-lapply(seq(n.block),
                    function(i) b.initial[(sum(g[1:i-1])+1):sum(g[1:i])])
    Sc<-union(which(U.lam<2*lam.vec[i]-lam.vec[i-1]),Sc)
    U.lam[Sc]<-NA
    
    iterating=TRUE
    while(iterating){
      X.block.new<-X.block
      X.block.new[Sc]<-NULL
      X.new<-do.call(cbind,X.block.new)
      if (length(Sc)==0){
        g.new<-g
      }else{
        g.new<-g[-Sc]
      }
      b.block.new<-b.block
      b.block.new[Sc]<-NULL
      b.initial<-unlist(b.block.new)
      res<-irlsbmd(X.new,y,b0.initial,
                   b.initial,v,w,tau,rho,g.new,
                   lam.vec[i],tol,
                   maxit,maxit.in)
      b0.update<-res$b0
      b.update<-res$b
      
      V<-intersect(Sc,which(U.lam>lam.vec[i]))
      Sc<-setdiff(Sc,V)
      
      b.initial<-rep(NA,p)
      for(j in Sc){
        b.initial[(sum(g[1:j-1])+1):sum(g[1:j])]<-0
      }
      b.initial[is.na(b.initial)]<-b.update
      if(length(V)==0){
        iterating=FALSE
        U.lam[!is.na(U.lam)]<-res$U.lam
        b0.vec[i]<-b0.update
        b.vec[,i]<-b.initial
      }else{
        b0.initial<-b0.update
      }
    }
  }
  return(list(b0=b0.vec,beta=b.vec,lam=lam.vec))
}

