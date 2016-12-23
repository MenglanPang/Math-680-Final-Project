
positive.part<-function(x){
  if(x>0){
    return(x)
  }else{
    return(0)
  }
}

# g is a vector indicating length of each block e.g. if g=c(1,3,2), then g represents
# a collection of 3 groups of length 1, 3, 2, respectively
# X is required to be arranged in the same way as g
# This function computes the grouped elastic net penalized solution for tweedie model
irlsbmd<-function(X,y,b0,b,v,w,tau,rho,g,lam,tol=1e-8,maxit=1e20,maxit.in=3e8){
  n.block<-length(g)
  X.block<-lapply(seq(n.block),function(i) X[,(sum(g[1:i-1])+1):sum(g[1:i])])
  Xb<-as.vector(crossprod(t(X),b))
  k<-0
  iterating.outer = TRUE
  #outer layer
  while(iterating.outer){
    k<-k+1
    b.block<-lapply(seq(n.block),function(i) b[(sum(g[1:i-1])+1):sum(g[1:i])])
    coeff.v<-(rho-1)*y*exp(-(rho-1)*(b0+Xb))+
      (2-rho)*exp((2-rho)*(b0+Xb))
    v.update<-as.vector(v*coeff.v)
    y.update<-as.vector(b0+Xb+(v/v.update)*(y*exp(-(rho-1)*(b0+Xb))-exp((2-rho)*(b0+Xb))))
    H<-lapply(X.block, function(x) crossprod(x,v.update*x))
    gamma<-sapply(H, function(x) max(eigen(x)$values))
    H0<-gamma0<-sum(v.update)
    
    iterating.inner = TRUE
    j<-0
    b.update<-b
    b0.update<-b0
    Xb.update<-Xb
    #inner layer
    while(iterating.inner){
      j<-j+1
      b.block.update<-lapply(seq(n.block),function(i) b.update[(sum(g[1:i-1])+1):sum(g[1:i])])
      U<-lapply(X.block,function(x) {x<-as.matrix(x) 
                                     apply(-v.update*(y.update-b0.update-Xb.update)*x,2,sum) })
      b.block.ls<-lapply(seq(n.block),function(i) (gamma[i]*b.block.update[[i]]-U[[i]])*
        positive.part(1-lam*tau*w[i]/sqrt(crossprod(gamma[i]*b.block.update[[i]]-U[[i]])))/
        (gamma[i]+lam*(1-tau)))
      b.ls<-unlist(b.block.ls)
      Xb.ls<-as.vector(crossprod(t(X),b.ls))
      U0<-sum(-v.update*(y.update-b0.update-Xb.ls))
      b0.ls<-b0.update-U0/gamma0

      if((j  >= maxit.in)||(sqrt(crossprod(c(b0.update-b0.ls,b.update-b.ls))) < tol)){
        iterating.inner = FALSE
      }
      b.update<-b.ls
      b0.update<-b0.ls
      Xb.update<-Xb.ls
    }
    if((k >= maxit)||(sqrt(crossprod(c(b0-b0.update,b-b.update))) < tol)){
      iterating.outer = FALSE
    }
    b<-b.ls
    b0<-b0.ls
    Xb<-Xb.update
  }
  U<-lapply(X.block,function(x) {x<-as.matrix(x)
                                 apply((-y*exp(-(rho-1)*(b0+Xb))+
                                         exp((2-rho)*(b0+Xb)))*v*x,2,sum)})
  # U.lam is U/tau*w, it is the value compared to lambda in KKT condition checking
  U.lam<-sapply(seq(n.block),function(i) sqrt(crossprod(U[[i]]))/(tau*w[i]))
  return(list(b0=b0,b=b,U.lam=U.lam))
}
