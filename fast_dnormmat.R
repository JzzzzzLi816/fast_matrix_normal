#fast_dmatnorm <- function(Z, M, U, V, log=TRUE, Precision=FALSE) {
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
  
  #cholesky check?
  
  Ru <- chol(U)
  Rv <- chol(V)
  
  # solving for A
  D <- Z - M 
  #print(dim(D))
  A2 = forwardsolve(Ru, D, upper.tri=TRUE,transpose=TRUE)
  
  #A = backsolve(Rv, A2)
  A = A2 %*% chol2inv(Rv)
  # calculate log density first
  # rcpp?
  numerator = -0.5 * sum(A^2)
  denom = (n*p/2)*log(2*pi) + sum(n*log(diag(Rv))) + sum(p*log(diag(Ru)))
  
  log.dens = numerator - denom
  
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
#}


set.seed(123)
A <- datasets::CO2[1:10, 4:5]
M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11))
V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)
#' V # Right covariance matrix (2 x 2), say the covariance between parameters.
U <- diag(10) # Block of left-covariance matrix ( 84 x 84), say the covariance between subjects.

fast_dmatnorm(A, M, U, V)
fast_dmatnorm(A, M, U, V, log = FALSE)


dmatnorm(A, M, U, V)
dmatnorm(A, M, U, V, log = FALSE)

# matrix norm code
dmatnorm <- function(
    A,
    M,
    U,
    V,
    tol = .Machine$double.eps^0.5,
    log = TRUE
) {
  n <- nrow(A)
  p <- ncol(A)
  
  # Checks
  if (is.data.frame(A)) A <- as.matrix(A)
  #if (sum(dim(A) == dim(M)) != 2) stop("M must have same dimensions as A.")
  #check_matnorm(s = 1, M, U, V, tol)
  
  # The Log Density
  log.dens <- (-n * p / 2) * log(2 * pi) - p / 2 * log(det(U)) - n / 2 * log(det(V)) +
    -1 / 2 * tr(solve(U) %*% (A - M) %*% solve(V) %*% t(A - M))
  
  # Return
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}
