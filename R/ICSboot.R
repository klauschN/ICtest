# input of the following function is 
# x: vector of eigenvalues (ordered)
# k: the number of successive elements for which the
#    varinace is computed 1 < k <= length(x)

# Output
# - rolling variances of order k for x
# - order vector which permutes x so that the elements with the smallest
#   variance are at the end
# - index of the k elements which gave smallest variance
# - orderd x
# - ordered variances 

fixOrder <- function(x, k)
{
  P <- length(x)
  Index <- 1:P
  RollVarX <-  roll_var(x,k)
  Start <- which.min(RollVarX)
  End <- Start + k - 1
  Smallest <- Start:End
  Order <- c(Index[-Smallest], Index[Smallest])
  Ox <- x[Order]
  RollVarOX <- roll_var(Ox,k)
  
  RES <- list(RollVarX = RollVarX, Selected=Smallest, 
              Order=Order, xOrdered=Ox, RollVarOX=RollVarOX)
  RES 
}



#### Teststatistic:, VARIANCE of last p-k eigenvalues!

ICS_boot_teststatistic <- function (X, k, S1=cov, S2=cov4, S1args=NULL, S2args=NULL) 
{
  n <- nrow(X)
  p <- ncol(X)
  MEAN <- colMeans(X)
  Xc <- sweep(X, 2, MEAN, "-")
  COV <- do.call("S1", c(list(Xc), S1args))
  EVD.COV <- eigen(COV, symmetric = TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*% 
    t(EVD.COV$vectors)
  Y <- tcrossprod(Xc, COV.inv.sqrt)
  COV4 <- do.call("S2", c(list(Y), S2args))
  EVD.COV4 <- eigen(COV4, symmetric = TRUE, only.values=TRUE)
  D <- EVD.COV4$values
  orderD <- fixOrder(D, p-k)
  Dord <- orderD$xOrdered
  TEST.STATISTIC.X <- sum((Dord[(k + 1):p] - mean(Dord[(k + 1):p]))^2)
  return(TEST.STATISTIC.X)
}

# Subfunction 1, NGCA model

ICS_boot_normal_S1 <- function (Z1, Z2, Winv, MEAN, k, S1=cov, S2=cov4, S1args=NULL, S2args=NULL) 
{
  n <- nrow(Z1)
  p <- ncol(Winv)
  covZ2 <- cov(Z2)
  Z2tilde <- rmvnorm(n, sigma=covZ2)
  Z1tilde <- Z1[sample(1:n, n, replace = TRUE), , drop = FALSE]
  Zstar <- cbind(Z1tilde, Z2tilde)
  Xstar <- tcrossprod(Zstar, Winv)
  Xstar <- sweep(Xstar, 2, MEAN, "+")
  TEST.STATISTIC.Xstar <- ICS_boot_teststatistic(Xstar, k, S1=S1, S2=S2, S1args=S1args, S2args=S2args)
  TEST.STATISTIC.Xstar
}

# Subfunction 1, NGICA model

ICS_boot_normal_S2 <- function (Z1, Z2, Winv, MEAN, k, S1=cov, S2=cov4, S1args=NULL, S2args=NULL) 
{
  n <- nrow(Z1)
  p <- ncol(Winv)
  covZ2 <- cov(Z2)
  Z2tilde <- rmvnorm(n, sigma=covZ2)
  Z1tilde <- apply(Z1, 2, SampFunction, n = n)
  Zstar <- cbind(Z1tilde, Z2tilde)
  Xstar <- Xstar2 <- tcrossprod(Zstar, Winv)
  Xstar <- sweep(Xstar, 2, MEAN, "+")
  TEST.STATISTIC.Xstar <- ICS_boot_teststatistic(Xstar, k, S1=S1, S2=S2, S1args=S1args, S2args=S2args)
  TEST.STATISTIC.Xstar
}

# Internal function NGCA model

ICS_subnormal_boot_S1 <- function (X, k, S1=cov, S2=cov4, S1args=NULL, S2args=NULL, n.boot = 200) 
{
  n <- nrow(X)
  p <- ncol(X)
  MEAN <- colMeans(X)
  Xc <- sweep(X, 2, MEAN, "-")
  COV <- do.call("S1", c(list(Xc), S1args))
  EVD.COV <- eigen(COV, symmetric = TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*% 
    t(EVD.COV$vectors)
  Y <- tcrossprod(Xc, COV.inv.sqrt)
  COV4 <- do.call("S2", c(list(Y), S2args))
  EVD.COV4 <- eigen(COV4, symmetric = TRUE)
  D <- EVD.COV4$values
  orderD <- fixOrder(D, p-k)
  Dord <- orderD$xOrdered
  W <- crossprod(EVD.COV4$vectors, COV.inv.sqrt)[orderD$Order, ]
  TEST.STATISTIC.X <- sum((Dord[(k + 1):p] - mean(Dord[(k + 1):p]))^2)
  names(TEST.STATISTIC.X) = "T"
  PARAMETER <- n.boot
  names(PARAMETER) <- c("replications")
  Z <- tcrossprod(Xc, W)
  Z1 <- Z[, seq_len(k), drop = FALSE]
  Z2 <- Z[, -(0:k), drop = FALSE]
  Winv <- solve(W)
  TEST.STATISTICS.Xstar <- t(replicate(n.boot, ICS_boot_normal_S1(Z1, Z2,
                                                                  Winv, MEAN, k, S1=S1, S2=S2, S1args=S1args, S2args=S2args)))
  names(Z) <- paste0("IC.", 1:p)
  PVAL <- (sum(TEST.STATISTIC.X < TEST.STATISTICS.Xstar) + 
             1)/(n.boot + 1)
  ALTERNATIVE <- paste0("the last ", p - k, " components are not gaussian")
  RES <- list(statistic = n * TEST.STATISTIC.X, p.value = PVAL, 
              parameter = PARAMETER, alternative = ALTERNATIVE, k = k, 
              W = W, S = Z, D = D, MU = MEAN)
  RES
}

# Internal function NGICA model

ICS_subnormal_boot_S2 <- function (X, k, S1=cov, S2=cov4, S1args=NULL, S2args=NULL, n.boot = 200) 
{
  n <- nrow(X)
  p <- ncol(X)
  MEAN <- colMeans(X)
  Xc <- sweep(X, 2, MEAN, "-")
  COV <- do.call("S1", c(list(Xc), S1args))
  EVD.COV <- eigen(COV, symmetric = TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*% 
    t(EVD.COV$vectors)
  Y <- tcrossprod(Xc, COV.inv.sqrt)
  COV4 <- do.call("S2", c(list(Y), S2args))
  EVD.COV4 <- eigen(COV4, symmetric = TRUE)
  D <- EVD.COV4$values
  orderD <- fixOrder(D, p-k)
  Dord <- orderD$xOrdered
  W <- crossprod(EVD.COV4$vectors, COV.inv.sqrt)[orderD$Order, ]
  TEST.STATISTIC.X <- sum((Dord[(k + 1):p] - mean(Dord[(k + 1):p]))^2)
  names(TEST.STATISTIC.X) = "T"
  PARAMETER <- n.boot
  names(PARAMETER) <- c("replications")
  Z <- tcrossprod(Xc, W)
  colnames(Z) <- paste0("IC.", 1:p)
  Z1 <- Z[, seq_len(k), drop = FALSE]
  Z2 <- Z[, -(0:k), drop = FALSE]
  Winv <- solve(W)
  TEST.STATISTICS.Xstar <- t(replicate(n.boot, ICS_boot_normal_S2(Z1,Z2, 
                                                                  Winv, MEAN, k, S1=S1, S2=S2, S1args=S1args, S2args=S2args)))
  PVAL <- (sum(TEST.STATISTIC.X < TEST.STATISTICS.Xstar) + 
             1)/(n.boot + 1)
  ALTERNATIVE <- paste0("the last ", p - k, " components are not gaussian")
  RES <- list(statistic = n * TEST.STATISTIC.X, p.value = PVAL, 
              parameter = PARAMETER, alternative = ALTERNATIVE, k = k, 
              W = W, S = Z, D = D, MU = MEAN)
  RES
}

# Main toplayer function

ICSboot <- function (X, k, S1=cov, S2=cov4, S1args=NULL, S2args=NULL, n.boot = 200, s.boot = "B1")
  {
    DNAME <- deparse(substitute(X))
    S1NAME <- deparse(substitute(S1))
    S2NAME <- deparse(substitute(S2))
    X <- as.matrix(X)
    s.boot <- match.arg(s.boot, c("B1", "B2"))
    if (s.boot == "B1") {
      RES <- ICS_subnormal_boot_S1(X = X, k = k, S1=S1, S2=S2, S1args=S1args, S2args=S2args, n.boot=n.boot)
      METHOD <- paste("ICA subgaussianity bootstrapping test using", S1NAME, "-", S2NAME, "and strategy B1")
    }
    else {
      RES <- ICS_subnormal_boot_S2(X = X, k = k, S1=S1, S2=S2, S1args=S1args, S2args=S2args, n.boot=n.boot)
      METHOD <- paste("ICA subgaussianity bootstrapping test using", S1NAME, "-", S2NAME, "and strategy B2")
    }
    RES <- c(RES, method = METHOD, data.name = DNAME, s.boot = s.boot)
    RES <- RES[c(1:3, 10:11, 4:9, 12)]
    class(RES) <- c("ictest", "htest")
    RES
  }




