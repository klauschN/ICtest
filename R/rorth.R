rorth <- function(k) { 
  temp <- qr(matrix(rnorm(k * k), k))
  sign_vec <- sign(diag(qr.R(temp)))
  sweep(qr.Q(temp), 2, sign_vec, "*")
  }
