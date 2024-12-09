
fast_pnormmat <-function(Lower = -Inf, # Lower bound matrix
                         Upper = Inf, # Upper bound matrix
                         M = NULL, # Mean matrix
                         U_cov = NULL, 
                         V_cov = NULL, 
                         U_prec = NULL, 
                         V_prec = NULL, 
                         useCov = TRUE,
                         log=TRUE, 
                         method = "naive_monte_carlo", # naive_monte_carlo or vectorized
                         N = NULL, # optional number of samples for monte carlo 
                         tol=1e-8,
                         algorithm = mvtnorm::GenzBretz() # default algorithm for mvtnorm::pmvnorm
                         ) {
  
  # Matrix dimension checks
  dc <- if (useCov) {
    dim_check(M = M, U = U_cov, V = V_cov)
  } else {
    dim_check(M = M, U = U_prec, V = V_prec)
  }
  M <- dc[[1]]; n <- dc[[4]]; p <- dc[[5]]
  
  # Check Lower and Upper bounds using vectorized_check
  bounds <- vectorized_check(Lower, Upper, n, p)
  Lower <- bounds$Lower
  Upper <- bounds$Upper
  
  # Validate covariance or precision matrices
  if (useCov) {
    check_matnorm(Z = NULL, M = M, U = U_cov, V = V_cov, tol = tol)
    U <- U_cov; V <- V_cov
  } else {
    check_matnorm(Z = NULL, M = M, U = U_prec, V = V_prec, tol = tol)
    U <- solve(U_prec); V <- solve(V_prec)
  }
  
  if (method == "vectorized") {
    # Vectorized method using mvtnorm
    cdf <- mvtnorm::pmvnorm(
      lower = Lower,
      upper = Upper,
      mean = as.vector(M),
      sigma = kronecker(V, U),
      algorithm = algorithm
    )
    method <- "mvnorm computation"
  } else {
    # Monte Carlo method
    if (is.null(N)) {
      N <- max(n*p*100)  # Adaptive sample size
    }
    # Use our sampler
    samples <- fast_rmatnorm(num_samp = N, 
                             n = n, 
                             p = p, 
                             M = M, 
                             U_cov = U, 
                             V_cov = V, 
                             U_prec = U, 
                             V_prec = V, 
                             useCov = useCov)
    within_bounds <- sweep(samples, c(1, 2), Lower, `>=`) & sweep(samples, c(1, 2), Upper, `<=`)
    valid_samples <- apply(within_bounds, 3, all)
    count <- sum(valid_samples)
    cdf <- count / N
    method <- "monte carlo simulation"
  }
  
  df <- data.frame(method = method, cdf = cdf, log_cdf = log(cdf))
  return(df)
}
