#' Approximate the CDF of a Matrix Normal Distribution
#'
#' This function calculates the cumulative distribution function (CDF) of a matrix normal distribution
#' It supports three methods: naive Monte Carlo approximation, Sobol sequence-based quasi-Monte Carlo, and direct computation using mvtnorm::pmvnorm.
#'
#' @param Lower A matrix of lower bounds for the integration region. Defaults to \code{-Inf}.
#' @param Upper A matrix of upper bounds for the integration region. Defaults to \code{Inf}.
#' @param M The mean of the matrix normal distribution.
#' @param U_cov The row covariance matrix. Used when \code{useCov = TRUE}.
#' @param V_cov The column covariance matrix. Used when \code{useCov = TRUE}.
#' @param U_prec A matrix specifying the row precision matrix (inverse of covariance). Used when \code{useCov = FALSE}.
#' @param V_prec A matrix specifying the column precision matrix (inverse of covariance). Used when \code{useCov = FALSE}.
#' @param useCov A logic indicating whether to use covariance matrices (\code{TRUE}) or precision matrices (\code{FALSE}). Defaults to \code{TRUE}.
#' @param method A character string specifying the method to use for the CDF approximation. Defaults to \code{"naive_monte_carlo"}. Options are: 
#' \itemize{
#'   \item \code{"naive_monte_carlo"}: Uses naive Monte Carlo integration.
#'   \item \code{"sobol"}: Uses Sobol sequence-based quasi-Monte Carlo integration (requires the `randtoolbox` package).
#'   \item \code{"pmvnorm"}: Uses the `mvtnorm::pmvnorm` function for direct computation.
#' }
#' @param N An integer specifying the number of samples to draw for Monte Carlo or Sobol methods. Defaults to \code{NULL}, in which case an adaptive sample size is used.
#' @param tol A numeric value specifying the tolerance for precision checks in covariance/precision matrices. Defaults to \code{1e-8}.
#' @param max_iter An integer specifying the maximum number of iterations for iterative methods. Defaults to \code{1000}.
#' @param algorithm A function specifying the integration algorithm to use with mvtnorm::pmvnorm. Defaults to \code{mvtnorm::GenzBretz()}.
#' 
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{method}: The method used for the computation (e.g., \code{"naive monte carlo"}, \code{"sobol monte carlo"}, \code{"mvnorm computation"}).
#'   \item \code{cdf}: The estimated CDF value.
#'   \item \code{log_cdf}: The log of the CDF value.
#' }
#' 
#' @examples
#' # Define the dimensions and M, U, V
#' n <- 3
#' p <- 2
#' # Create the mean matrix
#' M <- matrix(0, nrow = n, ncol = p)
#' matrixU <- matrix(c(1, 2, 0, 2, -1, -2, 1, 3, 0), nrow = n, ncol = n)
#' matrixV <- matrix(c(2, 0, 1, 4), nrow = p, ncol = p)
#' U_cov <- crossprod(matrixU)  # Row covariance matrix
#' V_cov <- crossprod(matrixV)  # Column covariance matrix
#'
#' # Set the lower and upper bounds for the integration
#' Lower <- matrix(-10, nrow = n, ncol = p)
#' Upper <- matrix(10, nrow = n, ncol = p)
#'
#' # Compute the CDF using the naive Monte Carlo method
#' fast_pnormmat(Lower = Lower, Upper = Upper, M = M, 
#'               U_cov = U_cov, V_cov = V_cov, method = "naive_monte_carlo", N = 1000)
#'
#' # Compute the CDF using the Sobol method
#' fast_pnormmat(Lower = Lower, Upper = Upper, M = M, 
#'               U_cov = U_cov, V_cov = V_cov, method = "sobol", N = 1000)
#'
#' # Compute the CDF using the mvtnorm package
#' fast_pnormmat(Lower = Lower, Upper = Upper, M = M, 
#'               U_cov = U_cov, V_cov = V_cov, method = "pmvnorm")
#' 
#' @export
fast_pnormmat <-function(Lower = -Inf, # Lower bound matrix
                         Upper = Inf, # Upper bound matrix
                         M = NULL, # Mean matrix
                         U_cov = NULL, 
                         V_cov = NULL, 
                         U_prec = NULL, 
                         V_prec = NULL, 
                         useCov = TRUE,
                         method = "naive_monte_carlo", # naive_monte_carlo, sobol, pmvnorm
                         N = NULL, # optional number of samples for monte carlo 
                         tol=1e-8,
                         max_iter = 1000,
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
  
  if (method == "pmvnorm") {
    # Vectorized method using mvtnorm
    cdf <- mvtnorm::pmvnorm(
      lower = Lower,
      upper = Upper,
      mean = as.vector(M),
      sigma = kronecker(V, U),
      algorithm = algorithm
    )
    method <- "mvnorm computation"
  } else if (method == "sobol") {
    
    if (!requireNamespace("randtoolbox", quietly = TRUE)) {
      stop("The 'randtoolbox' package is required for the Sobol method. Please install it using install.packages('randtoolbox').")
    }
    if (is.null(N)) {
      N <- max(2000, 10*n*p)  # Adaptive sample size
    }
    sobol_points <- randtoolbox::sobol(n = N, dim = n * p)
    
    # Transform Sobol to normal distribution (vectorized) because that's how sobol works
    Z <- qnorm(sobol_points)  
    
    # Reshape Z to n x p x N
    Z <- array(Z, dim = c(n, p, N))
    
    Ru <- chol(U)
    Rv <- chol(V)
    samples <- array(0, dim = c(n, p, N))
    for (i in 1:N) { #back in matrix form
      samples[,,i] <- M + crossprod(Ru, Z[,,i]) %*% Rv
    }
    
    # Check bounds
    within_bounds <- sweep(samples, c(1, 2), Lower, `>=`) & sweep(samples, c(1, 2), Upper, `<=`)
    valid_samples <- apply(within_bounds, 3, all)
    count <- sum(valid_samples)
    cdf <- count / N
    method <- "sobol monte carlo"
    
  } else if (method == "naive_monte_carlo"){
    # Naive monte Carlo method
    if (is.null(N)) {
      N <- max(2000, 10*n*p)  # Adaptive sample size
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
    method <- "naive monte carlo"
  } else {
    stop("Invalid method. Please choose one of 'naive_monte_carlo', 'sobol', or 'pmvnorm'.")
  }
  
  df <- data.frame(method = method, cdf = cdf, log_cdf = log(cdf))
  return(df)
}

