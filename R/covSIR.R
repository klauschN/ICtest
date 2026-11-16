covSIR <- function(X, y, h = 10, ...) {
  n <- nrow(X)
  if (is.factor(y) || length(unique(y)) <= h) {
    # discrete response: one slice for each observed value
    slices <- as.integer(factor(y))
  } else {
    # repeated quantiles (ties in y) would give non-unique breaks
    breaks <- unique(quantile(y, probs = seq(0, 1, by = 1/h), ...))
    slices <- cut(y, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  }
  slicemean <- do.call(rbind, by(X, slices, colMeans))
  nh <- as.numeric(table(slices))
  COV <- cov.wt(slicemean, wt=nh/n, method="ML")$cov
  attr(COV, "nslices") <- length(nh)
  COV
}
