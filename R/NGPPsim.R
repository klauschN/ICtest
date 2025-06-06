NGPPsim <-
function(X, k, nl = c("skew", "pow3"), alpha = 0.8, N = 1000, eps = 1e-6, verbose = FALSE, maxiter = 100){
  
  X_name <- names(X)
  
  n <- nrow(X)
  
  # Parse non-linearities to integers
  nl_int <- match(nl, c("skew", "pow3", "tanh", "gauss"))

  if(length(alpha) != length(nl_int)){
    alpha <- c(alpha, 1 - sum(alpha))
  }
  
  # nl_int <- c(1, 2)
  # alpha <- 0.8
  
  # Whiten data
  MU <- colMeans(X)
  p <- ncol(X)
  cov_x <- cov(X)
  X <- whiten(X)
  
  # Estimate the initial k + 1 rows of U with FOBI
  FOBI_U <- FOBI(X)$W
  FOBI_S <- X%*%t(FOBI_U)
  FOBI_obj <- computeObj_C(FOBI_S, nl_int, alpha)
  U <- FOBI_U[(order(FOBI_obj, decreasing=TRUE)[1:(k + 1)]), , drop=FALSE]
  
  # Deflation-based estimation
    
  orth <- diag(p)
  
  for(i in 1:(k+1)){
    if(i > 1){
      orth <- orth - U[(i-1), ]%*%t(U[(i-1), ])
      
    }
    u <- U[i, ]
    
    crit <- 1
    iter <- 0
    while(crit > eps){
      
      # Orthogonalization
      Tvec <- computeTVec_C(u, X, nl_int, alpha)
      u2 <- orth%*%Tvec
      u2 <- u2/sqrt(sum(u2^2))
      
      # Criterion evaluation
      crit <- sqrt(min(sum((u2 - u)^2), sum((u2 + u)^2)))
      u <- u2
      
      iter <- iter + 1
      if(iter > maxiter){
        # stop("too many iterations")
        break
      } 
    }
    
    if(verbose == TRUE) print(iter)
    U[i, ] <- t(u)
  }
  
  # Compute the final objective function values
  S <- X%*%t(U)
  colnames(S) <- sapply(1:(k+1), function(i) paste("IC.", i, sep=""))
  obj <- computeObj_C(S, nl_int, alpha)
  
  
  # Simulated test of normality
  obj_gauss <- rep(0, N)
  
  for(i in 1:N){
    X_gauss <- matrix(rnorm((p - k)*n), nrow = n)
    
    # Whiten data
    X_gauss <- whiten(X_gauss)
    
    # Estimate the initial row u with FOBI
    FOBI_U <- FOBI(X_gauss)$W
    FOBI_S <- X_gauss%*%t(FOBI_U)
    FOBI_obj <- computeObj_C(FOBI_S, nl_int, alpha)
    u <- FOBI_U[(order(FOBI_obj, decreasing=TRUE)[1]), ]
    
    crit <- 1
    iter <- 0
    while(crit > eps){
      
      # Orthogonalization
      u2 <- computeTVec_C(u, X_gauss, nl_int, alpha)
      u2 <- u2/sqrt(sum(u2^2))
      
      # Criterion evaluation
      crit <- sqrt(min(sum((u2 - u)^2), sum((u2 + u)^2)))
      u <- u2
      
      iter <- iter + 1
      if(iter > maxiter){
       # stop("too many iterations")
        break
      } 
    }
    
    S_gauss <- X_gauss%*%u
    obj_gauss[i] <- computeObj_C(S_gauss, nl_int, alpha)
  }
  
  
  res_statistic <- obj[k + 1]
  names(res_statistic) <- "Objective function value"
  
  res_parameter <- N
  names(res_parameter) <- "Repetitions"
  
  res <- list(statistic = res_statistic,
              p.value = mean(res_statistic < obj_gauss),
              parameter = res_parameter,
              method = "Testing for sub-Gaussianity using NGPP",
              data.name = X_name,
              alternative = paste("There are less than", p - k, "Gaussian components"),
              k = k,
              W = U%*%symmetricPower_C(cov_x, -0.5),
              S = S,
              D = c(obj),
              MU = MU)
  
  class(res) <- c("ictest","htest")
  
  return(res)
}
