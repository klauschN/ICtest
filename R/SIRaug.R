# Helper function that returns the norm of the augmented
# eigenvector part for SIR
#
# input:
# X.C: centered original data
# y: response vector
# h: number of slices
# n: sample size (number of rows of X.C)
# p: data dimension (number of columns of X.C)
# naug: number of components which should be augmented
# ...: passed on to covSIR
#
# As SIR is affine equivariant the variance of the augmented
# components does not matter, they are therefore simulated
# as standard normal independent of y.

SIRaugNorm <- function(X.C, y, h, n, p, naug, ...)
{
Aug <- matrix(rnorm(n*naug), ncol=naug)
Xaug <- cbind(X.C, Aug)
Maug <- MSIR(Xaug, y = y, h = h, ...)
MaugEV <- eigen(Maug, symmetric = TRUE)
MaugEVj <- MaugEV$vectors[-(1:p), , drop = FALSE]
colSums(MaugEVj^2)
}

# Main function for SIR augmentation
# Input:
# X: n x p numeric data matrix
#   - y: response vector
#   - h: number of slices
#   - naug: number of augmented components
#   - nrep: number of times the augmentation is done
#   - ...: passed on to covSIR
SIRaug <- function (X, y, h = 10, naug = 1, nrep = 1, ...)
{
  data.name <- deparse(substitute(X))
  method <- "SIR"
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

  Mdata <- MSIR(X.C, y = y, h = h, ...)
  EV.Mdata <- eigen(Mdata, symmetric = TRUE)
  EV.vec <- EV.Mdata$vectors
  EV.val <- EV.Mdata$values

  LAMBDA2 <- c(EV.val, 0)

  NORMS <- replicate(nrep, SIRaugNorm(X.C, y, h, n, p, naug, ...))
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
  colnames(S) <- paste0("SIC.", 1:p)
  RES <- list(method = method, k = est.k, fn = fn, phin = phin,
              gn = gn, lambda = EV.val, W = W,
              S = S, MU = MEAN, data.name = data.name)
  class(RES) <- "ladle"
  RES
}
