# Quantile function of the Marchenko-Pastur distribution with
# ratio parameter gam (0 < gam <= 1) and unit scale.
#
# input:
# q: probability
# gam: ratio parameter

qMarchenkoPastur <- function(q, gam)
{
a <- (1 - sqrt(gam))^2
b <- (1 + sqrt(gam))^2
dens <- function(x) ifelse(x > 0, sqrt(pmax((b - x) * (x - a), 0)) / (2 * pi * gam * x), 0)
cdf <- function(t) if (t <= a) 0 else if (t >= b) 1 else integrate(dens, a, t)$value
uniroot(function(t) cdf(t) - q, c(a, b), tol = 1e-10)$root
}

# Helper function that returns the squared norms of the augmented
# part of the first K+1 eigenvectors of the covariance matrix of the
# augmented data.
#
# input:
# X.C: centered original data
# n: sample size (number of rows of X.C)
# p: data dimension (number of columns of X.C)
# naug: number of components which should be augmented
# sigma: noise standard deviation
# K: largest dimension considered
#
# If p + naug > n, the eigenvectors are obtained via the n x n
# Gram matrix as only the first K+1 eigenvectors are needed.

PCAaugHDNorm <- function(X.C, n, p, naug, sigma, K)
{
Aug <- matrix(rnorm(n * naug, mean = 0, sd = sigma), ncol = naug)
Aug.C <- sweep(Aug, 2, colMeans(Aug), "-")
if (p + naug <= n) {
  Xaug <- cbind(X.C, Aug.C)
  EVaug <- eigen(crossprod(Xaug), symmetric = TRUE)$vectors[, 1:(K + 1), drop = FALSE]
  res <- colSums(EVaug[-(1:p), , drop = FALSE]^2)
} else {
  G <- tcrossprod(X.C) + tcrossprod(Aug.C)
  EVG <- eigen(G, symmetric = TRUE)
  V <- EVG$vectors[, 1:(K + 1), drop = FALSE]
  res <- colSums(crossprod(Aug.C, V)^2) / EVG$values[1:(K + 1)]
}
res
}

# Main function for high-dimensional PCA augmentation
# (Radojicic and Virta, 2025)
# Input:
# X: n x p numeric data matrix
#   - naug: number of augmented components, if NULL
#           max(p, min(ceiling(n/2), 1000))
#   - nrep: number of times the augmentation is done
#   - sigma2: noise variance, if NULL estimated using the
#             Marchenko-Pastur corrected median eigenvalue
#   - K: largest dimension considered, if NULL min(20, min(p, n-1) - 1)
PCAaugHD <- function(X, naug = NULL, nrep = 1, sigma2 = NULL, K = NULL)
{
  data.name <- deparse(substitute(X))
  method <- "HDPA"
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  m <- min(p, n - 1)
  if (m < 3) stop("'X' must have at least 3 columns and 4 rows")
  if (is.null(naug)) naug <- max(p, min(ceiling(n / 2), 1000))
  checkPosInt(naug, "naug")
  checkPosInt(nrep, "nrep")
  naug <- as.integer(naug)
  nrep <- as.integer(nrep)
  if (is.null(K)) K <- min(20, m - 1)
  checkPosInt(K, "K")
  K <- as.integer(K)
  if (K > m - 1) {
    K <- m - 1L
    warning("'K' was reduced to ", K)
  }
  gam <- p / n
  rho <- naug / n

  MEAN <- colMeans(X)
  X.C <- sweep(X, 2, MEAN, "-")
  EV.Mdata <- eigen(crossprod(X.C) / n, symmetric = TRUE)
  EV.vec <- EV.Mdata$vectors
  EV.val <- EV.Mdata$values

  if (is.null(sigma2)) {
    LambdaMed <- EV.val[ceiling(m / 2)]
    sigma2 <- if (gam <= 1) LambdaMed / qMarchenkoPastur(0.5, gam) else
      LambdaMed / (gam * qMarchenkoPastur(0.5, 1 / gam))
  } else {
    if (!is.numeric(sigma2) || length(sigma2) != 1 || sigma2 <= 0) stop("'sigma2' must be a positive number")
  }

  NORMS <- replicate(nrep, PCAaugHDNorm(X.C, n, p, naug, sqrt(sigma2), K))
  NORMS <- matrix(NORMS, ncol = nrep)
  aveNORMS <- rowMeans(NORMS)

  # debiased eigenvalues
  a <- EV.val[1:(K + 1)] - sigma2 * (1 + gam)
  LAMBDA <- a / 2 + sqrt(pmax(a^2 - 4 * sigma2^2 * gam, 0)) / 2
  # rescaled norms of the augmented parts
  hn <- LAMBDA * (LAMBDA + (gam + rho) * sigma2) / (LAMBDA + sigma2) * aveNORMS
  # criterion: differences of consecutive values
  gn <- diff(hn)
  est.k <- which.min(gn)

  S <- X.C %*% EV.vec
  colnames(S) <- paste0("PC.", 1:p)
  RES <- list(method = method, k = est.k, fn = hn, phin = LAMBDA,
              gn = gn, lambda = EV.val, W = EV.vec,
              S = S, MU = MEAN, data.name = data.name, sigma2 = sigma2,
              index.start = 1)
  class(RES) <- "ladle"
  RES
}
