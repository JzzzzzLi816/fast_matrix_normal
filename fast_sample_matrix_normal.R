sample_matrix_normal <- function(num_samp = 1, n = 10,p=5, M=NULL, U_cov=NULL, V_cov=NULL, U_prec = NULL, V_prec = NULL, useCov = TRUE){
  if(!is.null(M)){
    n <- nrow(M)
    p <- ncol(M)
  }
  if(is.null(M)){
    M = matrix(0,nrow=n,ncol=p)
  }
  result <- array(rnorm(n * p * num_samp), dim = c(n, p, num_samp))
  if (useCov) {
    if(is.null(U_cov)) {
      U_cov <- diag(n)
    }
    if (is.null(V_cov)) {
      V_cov <- diag(p)
    }
    Ru = chol(U_cov)
    Rv = chol(V_cov)
    return(sapply(1:num_samp, function(i) M + crossprod(Ru, result[,,i]) %*% Rv))
  } else {
    if(is.null(U_prec)) {
      U_prec <- diag(n)
    }
    if (is.null(V_prec)) {
      V_prec <- diag(p)
    }
    Ru = chol(U_prec)
    Rv = chol(V_prec)
    Ru_inv = backsolve(Ru, diag(nrow(Ru)))
    Rv_inv = chol(Rv, diag(nrow(Rv)))
    return(sapply(1:num_samp, function(i) M + tcrossprod(Ru_inv %*% result[,,i], Rv_inv)))
 }
}

n = 100
p = 100
U <- matrix(0.5,nrow=n,ncol=n) + 0.5*diag(n)
V <- matrix(0.8,nrow=p,ncol=p) + 0.2*diag(p)

Y<- sample_matrix_normal(U=U, V=V)



#cor(Y)
#cor(t(Y))
