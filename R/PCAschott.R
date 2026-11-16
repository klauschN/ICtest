# Schott's (2006) high-dimensional test for subsphericity
# Input:
# X: n x p numeric data matrix, p can be larger than n
# k: number of eigenvalues larger than the equal ones, 0 <= k <= p-2
PCAschott <- function(X, k)
    {
    DNAME <- deparse(substitute(X))
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    checkK(k, p, p - 2)

    MEAN <- colMeans(X)
    Xc <- sweep(X,2,MEAN,"-")
    COV <- crossprod(Xc)/(n-1)

    METHOD <- "PCA subsphericity test using Schott's test"

    EV <- eigen(COV, symmetric = TRUE)
    EV.VALUES <- EV$values
    ind <- (k+1):p
    q <- p-k

    # generalization of the sphericity statistic of Ledoit and Wolf (2002)
    U <- mean(EV.VALUES[ind]^2) / mean(EV.VALUES[ind])^2 - 1
    # Schott (2006), Section 5: n - 1 are the degrees of freedom of COV
    TESTSTATISTIC <- (n-1-k) * q * U / 2
    names(TESTSTATISTIC) = "T"
    PARAMETER <- (q*(q+1)/2)-1
    names(PARAMETER) <- c("df")
    PVAL <- pchisq(TESTSTATISTIC, df=PARAMETER, lower.tail = FALSE)
    ALTERNATIVE <- paste0("the last ",q, " eigenvalues are not equal")

    W <- EV$vectors
    S <- Xc %*% W
    colnames(S) <- paste0("PC.",1:p)

    RES <- list(statistic = TESTSTATISTIC, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=S, D=EV.VALUES,
    MU=MEAN, SCATTER=COV)
    class(RES) <- c("ictest", "htest")
    RES
    }
