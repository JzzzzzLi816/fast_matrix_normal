fast_dmatnorm <- function(Z, M, U, V, log=TRUE, Precision=FALSE, tol=1e-8) {
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
  # U,v square
  #  M nxp
  # check for nas
  # cholesky check - pos definite
  check_matnorm(M,U,V,tol)

    Ru <- chol(U)
    Rv <- chol(V)
    
    # solving for A
    D <- Z - M 

    # calculate log density first
    # RCPP for sum(A^2)
  if (!Precision) {
    A2 = forwardsolve(Ru, D, upper.tri=TRUE,transpose=TRUE)
    A = t(forwardsolve(t(Rv), t(A2)))
    numerator = -0.5 * sum(A^2)
    denom = (n*p/2)*log(2*pi) + sum(n*log(diag(Rv))) + sum(p*log(diag(Ru)))
  } else {
    B = Ru %*% D %*% t(Rv)
    numerator = -0.5 * sum(B^2)
    denom = (n*p/2)*log(2*pi) - sum(n*log(diag(Rv))) - sum(p*log(diag(Ru)))
  }
    
    log.dens = numerator - denom
  
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}


set.seed(123)
A <- datasets::CO2[1:10, 4:5]
M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11))
V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)
#' V # Right covariance matrix (2 x 2), say the covariance between parameters.
U <- diag(10) # Block of left-covariance matrix ( 84 x 84), say the covariance between subjects.

#system.time(
fast_dmatnorm(A, M, U, V)
#)
fast_dmatnorm(A, M, U, V, log = FALSE)

# original matrixnorm
#system.time(
dmatnorm(A, M, U, V)
#)
dmatnorm(A, M, U, V, log = FALSE)


is.sym <- function(X, tol=1e-8) {
  if (dim(X)[1] != dim(X)[2]) {FALSE}
  
  else if (all(abs(X-t(X)) < tol)) {TRUE}
  else {FALSE}
}



is.pos.def <- function(X, tol=1e-8) {
  eigenvalues <- eigen(X, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) > tol) {TRUE}
  else {FALSE}
}


# MATRIX CHECK
check_matnorm <- function(M, U, V, tol=1e-8) {
  if (anyNA(M)) {
    stop("M contains missing values.", call. = FALSE)
  }
  if (anyNA(U)) {
    stop("U contains missing values.")
  }
  if (anyNA(V)) {
    stop("V contains missing values.")
  }
  if (nrow(M) != nrow(U)) {
    stop("The mean matrix M has different sample size than the scale sample size
         matrix U. M has ", dim(M)[[1]], "rows, and U has ", dim(U)[[1]], ".")
  }
  if (ncol(M) != nrow(V)) {
    stop("The mean matrix M has different number of parameters than scale
         parameter matrix V: M  -- ", dim(M)[2], "; V -- ", dim(V)[1], ".")
  }
  if (!is.sym(U, tol)) {
    stop("U is not symmetric.")
  }
  if (!is.sym(V, tol)) {
    stop("V is not symmetric.")
  }
  
  if (!is.pos.def(U, tol)) {
    stop("U is not positive definite. Calculation may not be accurate.
         Possibly lower tolerance.")
  }
  if (!is.pos.def(V, tol)) {
    stop("V is not positive definite. Calculation may not be accurate.
         Possibly lower tolerance.")
  }
  return(invisible())
}
