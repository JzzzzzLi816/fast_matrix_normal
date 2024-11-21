sample_matrix_normal <- function(n = 10,p=5, M=NULL, U=NULL, V=NULL){
 
  if(!is.null(M)){
    n <- nrow(M)
    p <- ncol(M)
  } 
  if(!is.null(U)){
    n <- ncol(U)
  } 
  if(!is.null(V)){
    p <- ncol(V)
  } 
  if(is.null(U)){
    U <- diag(n)
  }
  if(is.null(V)){
    V <- diag(p)
  }
  if(is.null(M)){
    M = matrix(0,nrow=n,ncol=p)
  }
  Ru <- chol(U)
  Rv <- chol(V)
  Z = matrix(rnorm(n*p),nrow=n,ncol=p)
  Y <- M + crossprod(Ru,Z)%*%Rv
  return(Y)
}

n = 100
p = 100
U <- matrix(0.5,nrow=n,ncol=n) + 0.5*diag(n)
V <- matrix(0.8,nrow=p,ncol=p) + 0.2*diag(p)

Y<- sample_matrix_normal(U=U, V=V)



#cor(Y)
#cor(t(Y))
