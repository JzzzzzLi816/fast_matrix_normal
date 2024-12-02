dim_check <- function(M, U, V) {
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
  return(M,U,V, n, p)
}


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
check_matnorm <- function(Z = NULL, # For sampling function, there's no Z to check
                          M, 
                          U, 
                          V, 
                          tol=1e-8) {
  if (!is.null(Z) & anyNA(Z)) {
    stop("Z contains missing values.", call. = FALSE)
  }
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