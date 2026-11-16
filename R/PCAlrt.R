# Likelihood ratio test for subsphericity with Bartlett-Lawley correction
# Input:
# X: n x p numeric data matrix, n > p
# k: number of eigenvalues larger than the equal ones, 0 <= k <= p-2
PCAlrt <- function(X, k)
    {
    DNAME <- deparse(substitute(X))
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    checkK(k, p, p - 2)
    if (n <= p) stop("the likelihood ratio test requires more observations than variables, use PCAschott instead")

    MEAN <- colMeans(X)
    Xc <- sweep(X,2,MEAN,"-")
    COV <- crossprod(Xc)/(n-1)

    METHOD <- "PCA subsphericity test using the likelihood ratio test"

    EV <- eigen(COV, symmetric = TRUE)
    EV.VALUES <- EV$values
    ind <- (k+1):p
    q <- p-k
    MEANev <- mean(EV.VALUES[ind])

    # Bartlett-Lawley correction factor, see Lawley (1956) and Schott (2006), eq. (2)
    cn <- (n-1-k) - (2*q + 1 + 2/q)/6 + MEANev^2 * sum(1/(EV.VALUES[-ind]-MEANev)^2)

    TESTSTATISTIC <- cn * (q*log(MEANev)- sum(log(EV.VALUES[ind])))
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
