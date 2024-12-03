
source("fast_dnormmat.R")
source("fast_sample_matrix_normal.R")
source("helper.R")
fast_pnormmat <-function(Lower = -Inf, # Lower bound matrix
                         Upper = Inf, # Upper bound matrix
                         M = NULL, # Mean matrix
                         U_cov = NULL, 
                         V_cov = NULL, 
                         U_prec = NULL, 
                         V_prec = NULL, 
                         useCov = TRUE,
                         log=TRUE, 
                         method = NULL, # naive_monte_carlo or vectorized
                         N = NULL, # optional number of samples for monte carlo 
                         tol=1e-8,
                         algorithm = mvtnorm::GenzBretz() # default algorithm for mvtnorm::pmvnorm
                         ) {
  
  
  # If the user specified to vectorized to calculate precise result, otherwise we use monte carlo method.
  if(useCov) {
    dc <- dim_check(M = M, U = U_cov, V = V_cov)
    M <- dc$M
    U <- dc$U
    V <- dc$V
    n <- dc$n
    p <- dc$p
    check_matnorm(Z = NULL, M = M,U = U_cov,V = V_cov, tol = tol)
    if(method == "vectorized") {
      vectorized_check(Lower, Upper, n, p)
      
      # Compute the probability using mvtnorm::pmvnorm
      cdf <- mvtnorm::pmvnorm(
        lower = Lower,
        upper = Upper,
        mean = as.vector(M),
        sigma = kronecker(V, U),
        algorithm = algorithm,
        ...
      )
      method <- "mvnorm computation"
      return(cdf)
    } else { # naive monte carlo method
      if(is.null(N)) {
        N = n*p*10
      }
      samples <- fast_rmatnorm(num_samp = N, n = n, p = p, M = M,U = U_cov,V = V_cov, useCov = TRUE)
      # compare across num_samp samples
      within_bounds <- (samples >= Lower & samples <= Upper)
      
      # Check if all elements of each sample are within bounds
      valid_samples <- apply(within_bounds, 3, all)  # Applies across the 3rd dimension
      
      # Estimate the CDF
      count <- sum(valid_samples)
      cdf <- count / N
      method <- "monte carlor simulation"
      return(cdf)
    }
  } 
  # ignore precision matrix usage for now
  # else { # Use precision matrix
  #   check_matnorm(Z = NULL, M = M,U = U_prec,V = V_prec, tol = tol)
  #   if(method == "vectorized") {
  #     vectorized_check(Lower, Upper, n, p)
  #     
  #   } else { # naive monte carlo method
  #     if(is.null(N)) {
  #       N = n*p*10
  #     }
  #     samples <- fast_rmatnorm(num_samp = N, n = n, p = p, M = M,U = U_prec,V = V_prec, useCov = FALSE)
  #     
  #     within_bounds <- (samples >= Lower & samples <= Upper)
  #     valid_samples <- apply(within_bounds, 3, all)  
  #     
  #     count <- sum(valid_samples)
  #     cdf <- count / N
  #   }
  # }
  
  df <- data.frame(method = method, cdf = cdf, log_cdf = log(cdf))
  return(df)
}