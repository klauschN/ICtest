SIRasymp <- function(X, y, k, h=10, ...)
    {
    DNAME <- deparse(substitute(X))
    X <- as.matrix(X)
    checkK(k, ncol(X))
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n

    EVD.COV <- eigen(COV, symmetric=TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)

    Y <- tcrossprod(Xc, COV.inv.sqrt)
    S2 <- covSIR(Y, y, h=h, ...)
    EVD.S2 <- eigen(S2, symmetric = TRUE)

    D <- EVD.S2$values
    W <- crossprod(EVD.S2$vectors, COV.inv.sqrt)
    S <- tcrossprod(Xc, W)

    TESTSTATISTIC <-  n*sum(D[(k+1):p])
    names(TESTSTATISTIC) = "T"

    # number of slices actually used (smaller than h for a discrete y)
    H <- attr(S2, "nslices")
    if (k >= H - 1) stop("'k' must be smaller than the number of slices minus 1 (", H - 1, ")")
    PARAMETER <- (p-k)*(H-k-1)
    names(PARAMETER) <- "df"
    PVAL <- pchisq(TESTSTATISTIC, df=PARAMETER, lower.tail = FALSE)
    METHOD <- "SIR test for subspace dimension"
    ALTERNATIVE <- paste0("the last ",p-k, " eigenvalues are not zero")
      
    colnames(S) <- paste0("SIC.",1:p)
    
    RES <- list(statistic = TESTSTATISTIC, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=S, D=D,
    MU=MEAN)
    class(RES) <- c("ictest", "htest")
    RES
    }
