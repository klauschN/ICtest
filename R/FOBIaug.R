# Helper function that returns the norm of the augmented
# eigenvector part for FOBI
#
# input:
# X.C: centered original data
# n: sample size (number of rows of X.C)
# p: data dimension (number of columns of X.C)
# naug: number of components which should be augmented
#
# As FOBI is affine equivariant the variance of the augmented
# Gaussian components does not matter, they are therefore
# simulated as standard normal.

FOBIaugNorm <- function(X.C, n, p, naug)
{
Aug <- matrix(rnorm(n*naug), ncol=naug)
Xaug <- cbind(X.C, Aug)
Maug <- MFobi(Xaug)
MaugEV <- eigen(Maug, symmetric = TRUE)
MaugEVj <- MaugEV$vectors[-(1:p), , drop = FALSE]
colSums(MaugEVj^2)
}

# Main function for FOBI augmentation
# Input:
# X: n x p numeric data matrix
#   - naug: number of augmented components
#   - nrep: number of times the augmentation is done
FOBIaug <- function (X, naug = 1, nrep = 1)
{
  data.name <- deparse(substitute(X))
  method <- "FOBI"
  X <- as.matrix(X)
  checkPosInt(naug, "naug")
  checkPosInt(nrep, "nrep")
  p <- ncol(X)
  n <- nrow(X)
  MEAN <- colMeans(X)
  X.C <- sweep(X, 2, MEAN, "-")
  COV <- crossprod(X.C)/(n - 1)
  EVD.COV <- eigen(COV, symmetric = TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% tcrossprod(diag((1/EVD.COV$values)^0.5),
                                                 EVD.COV$vectors)

  Mdata <- MFobi(X.C)
  EV.Mdata <- eigen(Mdata, symmetric = TRUE)
  EV.vec <- EV.Mdata$vectors
  EV.val <- EV.Mdata$values

  LAMBDA2 <- c(EV.val, 0)

  NORMS <- replicate(nrep, FOBIaugNorm(X.C, n, p, naug))
  NORMS <- matrix(NORMS, ncol = nrep)

  aveNORMS <- rowMeans(NORMS)

  # fn: measure of variation of augmented eigenvectors
  fn <- (cumsum(c(0, aveNORMS[1:p])))
  # phin: normalized eigenvalues
  phin <- (LAMBDA2/(cumsum(LAMBDA2)+1))
  # combined info
  gn <- fn + phin
  est.k <- which.min(gn) - 1
  W <- crossprod(EV.vec, COV.inv.sqrt)
  S <- tcrossprod(X.C, W)
  colnames(S) <- paste0("IC.", 1:p)
  RES <- list(method = method, k = est.k, fn = fn, phin = phin,
              gn = gn, lambda = EV.val, W = W,
              S = S, MU = MEAN, data.name = data.name)
  class(RES) <- "ladle"
  RES
}
