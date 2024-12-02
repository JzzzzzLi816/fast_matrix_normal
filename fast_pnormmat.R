

fast_pnormmat <-function(Z, 
                         M, 
                         U, 
                         V, 
                         log=TRUE, 
                         Precision=FALSE, 
                         tol=1e-8) {
  dc = dim_check(M,U,V)
  M <- dc$M
  U <- dc$U
  V <- dc$V
  n <- dc$n
  p <- dc$p
  
  check_matnorm(Z,M,U,V,tol)
  
  # If the matrix is large, we use monte carlo integration;
  # If the matrix is smaller, we use vectorization and Kronecker method.
  if(n*p > 1000) {
    
  } else {
    
  }
}