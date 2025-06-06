FOBIasymp <- function(X, k, type="S3", model="NGCA" , method="satterthwaite")
    {
    DNAME <- deparse(substitute(X))
    type <-  match.arg(type, c("S1", "S2", "S3"))
    model <-  match.arg(model, c("NGCA", "ICA"))
    method <- match.arg(method, c("satterthwaite","integration","saddlepoint"))
    p <- ncol(X)
    n <- nrow(X)
    X <- as.matrix(X)

    MU <- colMeans(X)
    X.C <- sweep(X, 2, MU, "-")

    COV <- crossprod(X.C) / n
    EVD.COV <- eigen(COV, symmetric = TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% tcrossprod(diag((1/EVD.COV$values)^0.5),
                                                   EVD.COV$vectors)
    X.C.W <- tcrossprod(X.C, COV.inv.sqrt)
    r <- sqrt(rowSums(X.C.W^2))
    X.C.W.r <- r * X.C.W
    COV4 <- crossprod(X.C.W.r)/n
    EVD.COV4 <- eigen(COV4, symmetric = TRUE)
    W <- crossprod(EVD.COV4$vectors, COV.inv.sqrt)

    ORDER <- order((EVD.COV4$values-p-2)^2, decreasing = TRUE)

    D <- EVD.COV4$values[ORDER]
    W <- W[ORDER, ]
    Z <- tcrossprod(X.C, W)
    
    if (model=="ICA"){
        Z4.k <- colMeans(Z[,0:k, drop=FALSE]^4)
        Z4.pk <- rep(3, p-k)
        SumZ4 <- sum(c(Z4.k, Z4.pk))

        Sigma1 <- SumZ4 - p + 8
        Sigma2 <- 4
        MTEXT <- "in an ICA model"
        }
    
    if (model=="NGCA"){
        Z2 <- Z^2
        rsZ2 <-rowSums(Z2)
        MeanZ2 <- mean(rsZ2)
        VarZ2 <- mean((rsZ2-MeanZ2)^2)

        Sigma1 <- VarZ2 + 8
        Sigma2 <- 4
        MTEXT <- "in an NGCA model"
        }
    
    colnames(Z) <- paste0("IC.",1:p)
    
    RES <- switch(type, S1 = {
        Tk <- n*sum((D[(k+1):p]-(p+2))^2) / (p-k)
        names(Tk) <- "T"
        WEIGHTS <- c(2*Sigma1 / (p-k), 2*Sigma1/(p-k) + Sigma2)
        DFs <- c((p-k-1)*(p-k+2)/2,1)
        PARAMETER <- c(WEIGHTS[1], DFs[1], WEIGHTS[2], DFs[2])
        names(PARAMETER) <- c("w1", "df1", "w2", "df2")
        PVAL <- 1 - pchisqsum(Tk, DFs, WEIGHTS, method=method)
        METHOD <- c("FOBI subgaussianty test using a weighted sum of chisquare test", MTEXT)
        ALTERNATIVE <- paste0("there are fewer than ",p-k, " gaussian components")
        list(statistic = Tk, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=Z, D=D, MU=MU, sigma1 = 
        Sigma1, sigma2 = Sigma2, type=type)
        }, S2 = {
        Dk <- D[(k+1):p]
        nk <- length(Dk)
        m1R22 <- sum(Dk)/nk
        s2R22 <- sum((Dk-m1R22)^2)/nk
        PART1 <- n*(p-k)/(2*Sigma1) *  s2R22
        PART2 <- n*(m1R22-(p+2))^2 / (2*Sigma1 / (p-k) + Sigma2)
        DF1 <- (p-k-1)*(p-k+2)/2
        DF2 <- DF1 + 1
        Tk <- PART1 + PART2
        names(Tk) <- "T"
        PVAL <- 1- pchisq(Tk, df=DF2)
        PARAMETER <-  DF2
        names(PARAMETER) <- "df"
        METHOD <- c("FOBI subgaussianty test using a chisquare test (Variant I)", MTEXT)
        ALTERNATIVE <- paste0("there are fewer than ",p-k, " gaussian components")
        list(statistic = Tk, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=Z, D=D, MU=MU, sigma1 = 
        Sigma1, sigma2 = Sigma2, type=type)
        }, S3 = {
        Dk <- D[(k+1):p]
        nk <- length(Dk)
        m1R22 <- sum(Dk)/nk
        s2R22 <- sum((Dk-m1R22)^2)/nk
        Tk <- n*(p-k)/(2*Sigma1) *  s2R22
        names(Tk) <- "T"
        DF1 <- (p-k-1)*(p-k+2)/2
        PVAL <- 1- pchisq(Tk, df=DF1)
        PARAMETER <-  DF1
        names(PARAMETER) <- "df"
        METHOD <- c("FOBI subgaussianty test using a chisquare test (Variant II)",MTEXT)
        ALTERNATIVE <- paste0("there are fewer than ",p-k, " gaussian components")
        list(statistic = Tk, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=Z, D=D, MU=MU, sigma1 = 
        Sigma1, sigma2 = Sigma2, type=type, model=model)
        }
        )

    class(RES) <- c("ictest", "htest")
    RES
    }
