library(Rcpp)
sourceCpp("~/Desktop/Michigan/Biostat/BIOSTAT615/fast_matrix_normal/fast_matrixNormal.cpp")
set.seed(123)
# 2 different test sets
# to make sure the output is the same
A <- datasets::CO2[1:10, 4:5]
#A <- matrix(rnorm(10000*2, 5, 2), nrow = 10000, ncol = 2)
M <- cbind(stats::rnorm(10, 435, 296), stats::rnorm(10, 27, 11))
#V <- matrix(rnorm(1000*1000))
V <- matrix(c(87, 13, 13, 112), nrow = 2, ncol = 2, byrow = TRUE)
#' V # Right covariance matrix (2 x 2), say the covariance between parameters.
U <- diag(10) # Block of left-covariance matrix, say the covariance between subjects.

# check output matches below
# have to import bc functions below have checks removed
fast_dmatnorm(A, M, U, V)
fast_dmatnorm(A, M, U, V, log = FALSE)

# original matrixnorm
dmatnorm(A, M, U, V)
dmatnorm(A, M, U, V, log = FALSE)


# Random matrix
# checking speed of Rcpp
# checked 1000x100 matrix too 
B <- matrix(rnorm(1000*1000, 5, 2), nrow = 1000, ncol = 1000)

library(microbenchmark)
# Benchmark
microbenchmark(
  R = sum(B^2),
  Rcpp = sum_of_squares(B),
  times = 100
)


#library(matrixNormal)

# this is dmatnorm from matrixNormal without the checks
dmatnorm <- function(
    A,
    M,
    U,
    V,
    tol = .Machine$double.eps^0.5,
    log = TRUE) {
  n <- nrow(A)
  p <- ncol(A)
  
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

# fast dmatnorm with no checks
fast_dmatnorm <- function(Z, M, U, V, log=TRUE, Precision=FALSE, tol=1e-8) {
  Ru <- chol(U)
  Rv <- chol(V)
  
  # solving for A
  D <- Z - M 
  
  # calculate log density first
  if (!Precision) {
    A2 = forwardsolve(Ru, D, upper.tri=TRUE,transpose=TRUE)
    A = t(forwardsolve(t(Rv), t(A2)))
    # use Rcpp function
    numerator = -0.5 * sum_of_squares(A)
    denom = (n*p/2)*log(2*pi) + sum(n*log(diag(Rv))) + sum(p*log(diag(Ru)))
  } else {
    B = Ru %*% D %*% t(Rv)
    # use Rcpp function
    numerator = -0.5 * sum_of_squares(B)
    denom = (n*p/2)*log(2*pi) - sum(n*log(diag(Rv))) - sum(p*log(diag(Ru)))
  }
  
  log.dens = numerator - denom
  
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}

# bigger speed check vs matrixNormal
# try n=[10,100,1000]
create_mats <- function(n, p) {
  random_matrixU <- matrix(rnorm(n^2), nrow = n, ncol = n)
  random_matrixV <- matrix(rnorm(p^2), nrow = p, ncol = p)
  
  D <- matrix(rnorm(n * p, mean = 2, sd = 6), nrow = n, ncol = p)
  M <- matrix(rnorm(n * p, mean = 1, sd = 2), nrow = n, ncol = p)
  U <- crossprod(random_matrixU)
  V <- crossprod(random_matrixV)
  
  return(list(D = D, M = M, U = U, V = V))
}
  

res1 = create_mats(10,3)
res_10.3 <- microbenchmark(
  fast_matnorm = fast_dmatnorm(res1$D, res1$M, res1$U, res1$V),
  matnorm = dmatnorm(res1$D, res1$M, res1$U, res1$V),
  times = 100
)

res2 = create_mats(100,30)
res_100.30 <- microbenchmark(
  fast_matnorm = fast_dmatnorm(res2$D, res2$M, res2$U, res2$V),
  matnorm = dmatnorm(res2$D, res2$M, res2$U, res2$V),
  times = 100
)

res3 = create_mats(1000,30)
res_1000.30 <- microbenchmark(
  fast_matnorm = fast_dmatnorm(res3$D, res3$M, res3$U, res3$V),
  matnorm = dmatnorm(res3$D, res3$M, res3$U, res3$V),
  times = 100
)

library(ggplot2)
#bms = c(res_10.3, res_100.30, res_1000.30)
autoplot(res_10.3) + 
  labs(
    title = "Execution time for n=10, p=3"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  geom_violin(fill = "darkblue")

autoplot(res_100.30) + 
  labs(
    title = "Execution time for n=100, p=30"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  geom_violin(fill = "orange")


autoplot(res_1000.30) + 
  labs(
    title = "Execution time for n=1000, p=30"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  geom_violin(fill = "darkgreen")
