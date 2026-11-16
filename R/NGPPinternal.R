# Internal helper for NGPP, NGPPsim and NGPPest:
# checks the non-linearities and the weights and returns
# the non-linearities as integers and the full weight vector.
# If alpha has one element less than nl, the last weight is 1 - sum(alpha).

parseNL <- function(nl, alpha)
{
  nl <- match.arg(nl, c("skew", "pow3", "tanh", "gauss"), several.ok = TRUE)
  if (!is.numeric(alpha)) stop("'alpha' must be numeric")
  if (length(alpha) == length(nl) - 1) {
    alpha <- c(alpha, 1 - sum(alpha))
  } else if (length(alpha) != length(nl)) {
    stop("'alpha' must have length length(nl) or length(nl) - 1")
  }
  # weights given in full are used as such (as documented), but all weights,
  # including a completed last one, must lie in [0, 1]
  if (any(!is.finite(alpha)) || any(alpha < 0) || any(alpha > 1))
    stop("the weights 'alpha' must lie between 0 and 1")
  list(nl_int = match(nl, c("skew", "pow3", "tanh", "gauss")), alpha = alpha)
}
