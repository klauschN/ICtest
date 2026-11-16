#####
#
# general internal functions for ladle
#
########

# Function to compute for a specific iteration fn
# input
# EVboot - eigenvectors for a bootstrapping M matrix
# EVdata - eigenvectors for the M matrix for the observed data
# r - the value for the number of components to consider

fi <- function(EVboot, EVdata, r)
  {
  fni <- numeric(r)
  for (ii in 1:r) {
    fni[ii] <- det(crossprod(EVdata[,1:ii],EVboot[,1:ii]))
  }
  1-abs(fni)
  }

# print method for ladle object.
# prints everything but "S"
print.ladle <- function(x, ...)
  {
  print.listof(x[!(names(x) %in% c("S", "index.start"))], ...)
  }

# plot method for ladle object.
# ploting function is the same as plot.ictest
plot.ladle <- function (x, which = "all", ...)
  {
    if (is.null(x$S)) stop("the ladle object contains no components, plot the criterion with ladleplot() instead")
    which <- match.arg(which, c("all", "k"))
    if (which == "all") {
      S <- x$S
    }
    else {
      S <- x$S[, 0:x$k, drop = FALSE]
    }
    if (ncol(S) <= 2) {
      plot(S, ...)
    }
    else {
      if (any(class(S) %in% c("mts", "xts", "zoo"))) {plot(S, ...)} else {pairs(S, ...)}
    }
  }

summary.ladle <- function(object, ...)
{
  cat("\n")
  cat(object$method, "Ladle estimator for", object$data.name,"\n")
  cat("k: ", object$k)
  cat("\n")
}


# extracting the components is the same as for ictest
components.ladle <- function (x, which = "all", ...)
{
    if (is.null(x$S)) stop("the ladle object contains no components")
    which <- match.arg(which, c("all", "k"))
    if (which == "all") {
        S <- x$S
    }
    else {
        S <- x$S[, 0:x$k, drop = FALSE]
    }
    S
}

# check that k is a single integer with 0 <= k <= kmax
# (kmax defaults to p - 1)
checkK <- function(k, p, kmax = p - 1)
  {
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != round(k) || k < 0 || k > kmax)
    stop("'k' must be an integer between 0 and ", kmax)
  invisible(TRUE)
  }

# check that x is a single positive integer (used for naug, nrep, n.boot, N)
checkPosInt <- function(x, name)
  {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != round(x) || x < 1)
    stop("'", name, "' must be a positive integer")
  invisible(TRUE)
  }
