set.seed(123)
sim_func_par = function(n, p, mean_val = 0, sd_val = 1) {
  M = matrix(rnorm(n * p, mean = mean_val, sd = sd_val), nrow = n, ncol = p)
  random_matrixU <- matrix(rnorm(n^2), nrow = n, ncol = n)
  random_matrixV <- matrix(rnorm(p^2), nrow = p, ncol = p)
  U = crossprod(random_matrixU)
  V = crossprod(random_matrixV)
  return(list(M = M, U = U, V = V))
}
# check output matches below
# have to import bc functions below have checks removed
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
MSE = function(true_val, est_val) {
  mean((true_val - est_val)^2)
}

n_list = 10:200
p = 30
fast_time = c()
matrixNormal_time = c()
for (n in n_list) {
  n_MUV = sim_func_par(n, p)
  M = n_MUV$M
  U = n_MUV$U
  V = n_MUV$V
  fast_time = c(fast_time, system.time(fast_rmatnorm(M = M, U_cov = U, V_cov = V))[["elapsed"]])
  matrixNormal_time = c(matrixNormal_time, system.time(matrixNormal::rmatnorm(M = M, U = U, V = V))[["elapsed"]])
}

plot(x = n_list, y = fast_time, type = "l", col = "red", ylim = c(min(matrixNormal_time), max(matrixNormal_time)), xlab = "n", ylab = "Time consumed in second")
lines(x = n_list, y = matrixNormal_time, type = "l", col = "blue")
legend("topleft", legend = c("fast matrix normal sampling", "original matrix normal sampling"), col = c("red", "blue"), lty = 1)


# Input parameters
M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11))  # Mean matrix
V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)  # Right covariance
U <- diag(10)  # Left covariance
n <- nrow(M)  # Rows of M
p <- ncol(M)  # Columns of M

res_gen = function(num_samp, true_M, true_U, true_V, multi_fast) {
  n = nrow(true_U)
  p = nrow(true_V)
  sampled_mean = apply(multi_fast, c(1, 2), mean)
  mse_mean = MSE(true_M, sampled_mean)
  flattened_samples = matrix(NA, nrow = num_samp, ncol = n * p)
  for (i in 1:num_samp) {
    flattened_samples[i, ] = as.vector(multi_fast[,,i])
  }
  empirical_cov = cov(flattened_samples)
  theoretical_cov = kronecker(V, U)
  mse_cov = MSE(empirical_cov, theoretical_cov)
  return(list(mse_mean = mse_mean, mse_cov = mse_cov))
}


# cov_case 
multi_fast_1000 = fast_rmatnorm(num_samp = 1000, M = M, U_cov = U, V_cov = V)
multi_fast_10000 = fast_rmatnorm(num_samp = 10000, M = M, U_cov = U, V_cov = V)

res_gen(1000, M, U, V, multi_fast_1000)
res_gen(10000, M, U, V, multi_fast_10000)
# precision case
U_prec = solve(U)
V_prec = solve(V)
multi_fast_1000_prec = fast_rmatnorm(num_samp = 1000, M = M, U_prec = U_prec, V_prec = V_prec, useCov = FALSE)
multi_fast_10000_prec = fast_rmatnorm(num_samp = 10000, M = M, U_prec = U_prec, V_prec = V_prec, useCov = FALSE)
res_gen(1000, M, U, V, multi_fast_1000_prec)
res_gen(10000, M, U, V, multi_fast_10000_prec)

