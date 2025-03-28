# Helper function that returns the norm of the augmented 
# noise eigenvector part
#
# input:
# X.C: centered original data
# n: sample size (number of rows of X.C)
# p: data dimension (number of columns of X.C)
# naug: number of components which should be augmented
# sigma: estimate of the noise standard deviation

PCAaugNorm <- function(X.C,n,p,naug,sigma)
{
Aug <- matrix(rnorm(n*naug, mean=0, sd=sigma), ncol=naug)
MEAN.Aug <- colMeans(Aug)
Aug.C <- sweep(Aug, 2, MEAN.Aug, "-")
Xaug <- cbind(X.C, Aug.C)
Maug <- crossprod(Xaug)/(n - 1)
MaugEV <- eigen(Maug, symmetric = TRUE)
MaugEVj <- MaugEV$vectors[-(1:p), , drop = FALSE]
colSums(MaugEVj^2)
}

# Main function for PCA augmentation
# Input:
# X: n x p numeric data matrix
# noise: rule to estimate the noise variance options are:
#   - median: median of eigenvalues
#   - quantile: mean of the  <= than the alpha 
#               quantile of the eigenvalues
#   - last: last eigenvalue
#   - known: if it is known
#   - naug: number of augmented components
#   - nrep: number of times the augmentation is done
#   - sigma2: if noise="known", this gives the noise variance
#   - alpha: the quantile to be used if noise="quantile"
PCAaug <- function (X, noise = "median", naug = 1, nrep = 1, sigma2=NULL, alpha=NULL) 
{
  data.name <- deparse(substitute(X))
  noise <- match.arg(noise, c("median", "last", "quantile", "known"))
  method <- "PCA"
  p <- ncol(X)
  n <- nrow(X)
  MEAN <- colMeans(X)
  X.C <- sweep(X, 2, MEAN, "-")
  Mdata <- crossprod(X.C)/(n - 1)
  EV.Mdata <- eigen(Mdata, symmetric = TRUE)
  EV.vec <- EV.Mdata$vectors
  EV.val <- EV.Mdata$values
  
  SigHat2 <- switch(noise,
                    median = median(EV.val),
                    last = min(EV.val),
                    quantile = { if (is.null(alpha)) stop("'alpha' must be numeric")
                      cut <- quantile(EV.val, probs = alpha)
                      mean(EV.val[EV.val <=  cut]) },
                    known = ifelse(is.null(sigma2), 
                                   stop("'sigma2' must be numeric"), 
                                   sigma2)
  )
  
  sigma <- sqrt(SigHat2)
  
  LAMBDA <- (EV.val-SigHat2)
  LAMBDA[LAMBDA < 0] <- 0
  LAMBDA2 <- c(LAMBDA,0)
  
  NORMS <- replicate(nrep, PCAaugNorm(X.C,n,p,naug,sigma))
  
  aveNORMS <- rowMeans(NORMS)
  
  # fn: measure of variation of augmented eigenvectors
  fn <- (cumsum(c(0, aveNORMS[1:p])))
  # phin: normalized eigenvalues
  phin <- (LAMBDA2/(cumsum(LAMBDA2)+1))
  # combined info
  gn <- fn + phin
  est.k <- which.min(gn) - 1
  W <- EV.vec
  S <- X.C %*% W
  colnames(S) <- paste0("PC.", 1:p)
  RES <- list(method = method, k = est.k, fn = fn, phin = phin, 
              gn = gn, lambda = EV.val, W = W, 
              S = S, MU = MEAN, data.name = data.name, sigma2 = SigHat2)
  class(RES) <- "ladle"
  RES
}
 
