library(Rcpp)
sourceCpp("~/Desktop/Michigan/Biostat/BIOSTAT615/fast_matrix_normal/fast_matrixNormal.cpp")
source("~/Desktop/Michigan/Biostat/BIOSTAT615/fast_matrix_normal/helper.R")
fast_dmatnorm <- function(Z, M, U, V, log=TRUE, Precision=FALSE, tol=1e-8) {
  dc = dim_check(M,U,V)
  M <- dc[1]
  U <- dc[2]
  V <- dc[3]
  n <- dc[4]
  p <- dc[5]
  # U,v square
  #  M nxp
  # cholesky check - pos definite
   check_matnorm(Z,M,U,V,tol)

    Ru <- chol(U)
    Rv <- chol(V)
    
    # solving for A
    D <- Z - M 

    # calculate log density first
    # RCPP for sum(A^2)
  if (!Precision) {
    A2 = forwardsolve(Ru, D, upper.tri=TRUE,transpose=TRUE)
    A = t(forwardsolve(t(Rv), t(A2)))
    # use Rcpp function
    numerator = -0.5 * sum_of_squares(A)#sum(A^2)
    denom = (n*p/2)*log(2*pi) + sum(n*log(diag(Rv))) + sum(p*log(diag(Ru)))
  } else {
    B = Ru %*% D %*% t(Rv)
    # use Rcpp function
    numerator = -0.5 * sum_of_squares(B)#sum(B^2)
    denom = (n*p/2)*log(2*pi) - sum(n*log(diag(Rv))) - sum(p*log(diag(Ru)))
  }
    
    log.dens = numerator - denom
  
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}











