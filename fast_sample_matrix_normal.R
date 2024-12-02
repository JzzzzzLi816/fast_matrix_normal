fast_rmatnorm <- function(num_samp = 1, n = 10, p = 5, M = NULL, U_cov = NULL, V_cov = NULL, U_prec = NULL, V_prec = NULL, useCov = TRUE) {
  if (!is.null(M)) {
    n <- nrow(M)
    p <- ncol(M)
  }
  if (is.null(M)) {
    M <- matrix(0, nrow = n, ncol = p)
  }
  
  if (useCov) {
    # Default covariance matrices
    if (is.null(U_cov)) {
      U_cov <- diag(n)
    }
    if (is.null(V_cov)) {
      V_cov <- diag(p)
    }
    
    # Cholesky decomposition of covariance matrices
    Ru <- chol(U_cov)
    Rv <- chol(V_cov)
    if (num_samp == 1) {
      Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
      Y <- M + crossprod(Ru, Z) %*% Rv
      return(Y)
    } else {
      return_res <- array(NA, dim = c(n, p, num_samp))
      Z <- array(rnorm(n * p * num_samp), dim = c(n, p, num_samp))
      for (i in 1:num_samp) {
        return_res[,,i] <- M + crossprod(Ru, Z[,,i]) %*% Rv
      }
      return(return_res)
    }
    
  } else {
    # Default precision matrices
    if (is.null(U_prec)) {
      U_prec <- diag(n)
    }
    if (is.null(V_prec)) {
      V_prec <- diag(p)
    }
    
    # Cholesky decomposition and inversion of precision matrices
    Ru <- chol(U_prec)
    Rv <- chol(V_prec)
    Ru_inv <- backsolve(Ru, diag(nrow(Ru)))
    Rv_inv <- backsolve(Rv, diag(nrow(Rv)))
    if (num_samp == 1) {
      Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
      Y <- M + tcrossprod(Ru_inv %*% Z, Rv_inv)
      return(Y)
    } else {
      return_res <- array(NA, dim = c(n, p, num_samp))
      Z <- array(rnorm(n * p * num_samp), dim = c(n, p, num_samp))
      for (i in 1:num_samp) {
        return_res[,,i] <- M + tcrossprod(Ru_inv %*% Z[,,i], Rv_inv)
      }
      return(return_res)
    }
  }
}
