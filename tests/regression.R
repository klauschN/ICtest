# Regression tests for bugs fixed in version 0.3-8
library(ICtest)
set.seed(1)

n <- 200
p <- 5
X <- cbind(rexp(n), runif(n), matrix(rnorm(n * (p - 2)), n))
Xdf <- as.data.frame(X)
expectError <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))

# PCAlrt: k = 0 gives a finite statistic, S equals Xc %*% W
res <- PCAlrt(X, 0)
stopifnot(is.finite(res$statistic))
Xc <- sweep(X, 2, colMeans(X))
stopifnot(isTRUE(all.equal(unname(res$S), unname(Xc %*% res$W))))
# the components agree with those of PCAasymp up to the signs of the columns,
# which depend on the LAPACK implementation
S1 <- unname(PCAlrt(X, 2)$S)
S2 <- unname(PCAasymp(X, 2)$S)
stopifnot(isTRUE(all.equal(S1, sweep(S2, 2, sign(colSums(S1 * S2)), "*"))))
expectError(PCAlrt(X, p - 1))
expectError(PCAlrt(matrix(rnorm(40 * 50), 40), 2))

# PCAschott: statistic (n-1-k) q U / 2, scale invariant, works for p > n
res <- PCAschott(X, 2)
l <- eigen(cov(X), only.values = TRUE)$values[3:p]
stopifnot(isTRUE(all.equal(unname(res$statistic), (n - 3) * (p - 2) * (mean(l^2) / mean(l)^2 - 1) / 2)))
stopifnot(isTRUE(all.equal(res$statistic, PCAschott(10 * X, 2)$statistic)))
stopifnot(is.finite(PCAschott(X, 0)$statistic))
expectError(PCAschott(X, p - 1))
res <- PCAschott(matrix(rnorm(40 * 50), 40), 2)
stopifnot(is.finite(res$statistic), res$parameter == 48 * 49 / 2 - 1)

# NGPP: tanh and gauss objectives are about zero for Gaussian data
Z <- matrix(rnorm(20000), ncol = 2)
stopifnot(all(abs(ICtest:::computeObj_C(Z, 3, 1)) < 1e-3))
stopifnot(all(abs(ICtest:::computeObj_C(Z, 4, 1)) < 1e-3))
# invalid arguments give errors, three non-linearities with two weights work
expectError(NGPP(X, 1, nl = "foo"))
expectError(NGPP(X, 1, method = "deflation"))
res <- NGPP(X, 2, nl = c("skew", "pow3", "tanh"), alpha = c(0.4, 0.3))
stopifnot(ncol(res$S) == 2)
stopifnot(NGPPsim(X, 1, N = 20)$data.name == "X")

# kSearch: backward search without rejection returns the smallest k
Y <- matrix(rnorm(n * 4), n)
res <- kSearch(Y, PCAasymp, search = "backward", early_stop = FALSE, alpha = 1e-8)
stopifnot(res$k == 0)
res2 <- kSearch(Y, ICtest::PCAasymp, search = "backward", early_stop = FALSE, alpha = 1e-8)
stopifnot(res2$k == 0)

# bootstrap tests: column names of S
stopifnot(identical(colnames(FOBIboot(X, 2, n.boot = 5)$S), paste0("IC.", 1:p)))
stopifnot(identical(colnames(ICSboot(X, 2, n.boot = 5)$S), paste0("IC.", 1:p)))
stopifnot(is.null(names(FOBIboot(X, 2, n.boot = 5)$S)))

# k checks
expectError(FOBIasymp(X, k = p))
expectError(ICSboot(X, k = p - 1, n.boot = 5))
expectError(PCAboot(X, k = p - 1, n.boot = 5))
expectError(PCAasymp(X, k = -1))
expectError(PCAasymp(X, k = 1.5))
stopifnot(is.finite(FOBIasymp(X, k = p - 1)$statistic))

# data frames are accepted
yc <- X[, 1] + rnorm(n, sd = 0.5)
invisible(PCAladle(Xdf, n.boot = 5))
invisible(FOBIladle(Xdf, n.boot = 5))
invisible(SIRladle(Xdf, yc, n.boot = 5))
invisible(PCAaug(Xdf))
invisible(SIRboot(Xdf, yc, k = 1, n.boot = 5))
expectError(PCAladle(X, n.boot = 5, ncomp = p))

# SIR with a discrete response
yd <- rbinom(n, 1, plogis(X[, 1]))
M <- covSIR(X, yd)
stopifnot(attr(M, "nslices") == 2)
res <- SIRasymp(X, yd, k = 0)
stopifnot(res$parameter == p)
expectError(SIRasymp(X, yd, k = 1))
invisible(SIRladle(X, factor(yd), n.boot = 5, ncomp = 2))
# ties in a continuous response
yt <- round(yc)
invisible(SIRasymp(X, yt, k = 1))

# PCAaugHD
res <- PCAaugHD(X)
stopifnot(res$k >= 1, res$k <= length(res$gn), res$index.start == 1)

# checks added after the 0.3-8 code review
# positive integer arguments
expectError(PCAaug(X, naug = 0))
expectError(FOBIaug(X, nrep = -1))
expectError(SIRaug(X, yc, naug = 1.5))
expectError(PCAaugHD(X, nrep = 1.5))
expectError(PCAboot(X, 1, n.boot = 0))
expectError(FOBIboot(X, 1, n.boot = 2.5))
expectError(ICSboot(X, 1, n.boot = NA))
expectError(SIRboot(X, yc, 1, n.boot = -1))
# NGPP: k and the weights alpha
expectError(NGPP(X, 0))
expectError(NGPP(X, p + 1))
expectError(NGPPsim(X, p, N = 5))
expectError(NGPP(X, 1, alpha = 1.5))
expectError(NGPP(X, 1, alpha = c(0.5, -0.1)))
invisible(NGPP(X, 1, nl = "pow3"))
# simulated p-values are never zero
stopifnot(NGPPsim(X, 0, N = 5)$p.value >= 1/6)
# log cosh does not overflow
stopifnot(is.finite(ICtest:::computeObj_C(matrix(c(800, rnorm(99))), 3, 1)))
# general ladle objects have no components
lad <- ladle(X, function(x) cov(x), n.boots = 5)
expectError(plot(lad))
expectError(components(lad))
# kSearch: empty range and unsuccessful search
expectError(kSearch(X, PCAasymp, min_dim = 3, max_dim = 2))
ks <- suppressWarnings(kSearch(X, PCAasymp, alpha = 1))
stopifnot(is.na(ks$k), inherits(ks, "kSearch"), !inherits(ks, "htest"))
invisible(capture.output(print(ks), summary(ks)))
