#' @title 
#' Relevant Component Estimation via Iterative Search
#' 
#' @author Katariina Perkonoja
#'
#' @description 
#' Performs dimension estimation using various search strategies (forward, backward, binary)
#' and hypothesis testing methods that evaluate whether the true dimension \eqn{d \leq k}.
#' The search continues until a p-value exceeds the specified alpha threshold, which indicates
#' that the null hypothesis is not rejected. The search can optionally continue to test all values, as there may not be a global optimum 
#' in the tested range.
#'
#' @param X A data matrix with \eqn{p>1} columns.
#' @param method A function that performs a hypothesis test given \code{X} and dimension \code{k}. 
#' For more details on supported tests, see the
#' "Details" section below.
#' 
#' @details
#' This function is designed to work with the following tests:
#' \itemize{
#'   \item \code{\link[ICtest]{FOBIasymp}}, \code{\link[ICtest]{FOBIboot}}
#'   \item \code{\link[ICtest]{ICSboot}}
#'   \item \code{\link[ICtest]{NGPPsim}}
#'   \item \code{\link[ICtest]{PCAasymp}}, \code{\link[ICtest]{PCAboot}}, \code{\link[ICtest]{PCAschott}}
#'   \item \code{\link[ICtest]{SIRasymp}}, \code{\link[ICtest]{SIRboot}}
#' }
#' These tests evaluate the null hypothesis that the true dimension \eqn{d \leq k}, and return
#' a p-value indicating whether the null is rejected or not.
#'
#' @param search A character string specifying the search strategy. 
#'    Options are \code{forward} (default), \code{backward} and \code{binary}. For more details on search techniques, see the
#' "Details" section below.
#' 
#' @details
#' The search work as follows:
#' \itemize{
#'   \item \code{"forward"} (default) performs a linear search starting from the smallest possible value of \code{k} and incrementing upward. Returns the smallest \code{k} where the null is not rejected.
#'   \item \code{"backward"} performs a linear search starting from the largest possible value of \code{k} and decrementing downward. Returns \code{k}+1, where \code{k} is the largest dimension where the null is rejected.
#'   \item \code{"binary"} splits the possible range of \code{k}s in half and performs a binary search. This search may skip testing some possible \code{k} values, but it is typically faster than linear search.
#'   }
#' 
#' @param alpha A numeric significance level (default = 0.05) used as the threshold for rejecting the null hypothesis.
#' @param early_stop Logical or \code{NULL}. Controls whether the search stops early once the null hypothesis is not rejected. 
#' If \code{NULL} (default), the behavior is determined by the selected search strategy: \code{TRUE} for \code{"forward"} and \code{"backward"}, 
#' and \code{FALSE} for \code{"binary"}. If explicitly set to \code{TRUE} or \code{FALSE}, this overrides the strategy's default behavior.
#' @param min_dim Optional integer to limit the minimum dimension to search. By default, this is set according to the \code{method} being used.
#' @param max_dim Optional integer to limit the maximum dimension to search. By default, this is set according to the \code{method} being used.
#' @param ... Additional arguments passed to the testing function \code{method}.
#' 
#' @return An object of class \code{ictest} inheriting from \code{htest}, with additional elements:
#' \describe{
#'   \item{\code{tested.ks}}{An integer vector of all tested values of \code{k} during the search.}
#'   \item{\code{tested.ks.pvals}}{A numeric vector of p-values corresponding to each tested \code{k}.}
#' }
#'
#' @examples
#' # Applying forward search with PCAasymp while evaluating 
#' # all valid values of k (early_stop = FALSE)
#' n <- 200
#' S <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5),  
#'            rnorm(n), rnorm(n), rnorm(n))
#' A <- rorth(5)
#' X <- S %*% t(A)
#' result <- kSearch(
#'   X = X, 
#'   method = PCAasymp, 
#'   alpha = 0.05, 
#'   search = "forward", 
#'   early_stop = FALSE 
#' )
#' 
#' result
#'
#' @export

kSearch <- function(X, 
                    method, 
                    search = c("forward", "backward", "binary"), 
                    alpha = 0.05, 
                    early_stop = NULL,
                    min_dim = NULL,
                    max_dim = NULL,
                    ...) {
  
  # Match the search argument
  search <- match.arg(search)
  
  # Number of columns (dimensions)
  p <- ncol(X)
  
  # Define the ranges for k based on the method
  k_range <- switch(deparse(substitute(method)),
                    "FOBIasymp" = 0:(p - 1),
                    "FOBIboot" = 0:(p - 1),
                    "ICSboot" = 1:(p - 2),
                    "NGPPsim" = 1:(p-1),
                    "PCAasymp" = 0:(p - 2),
                    "PCAboot" = 0:(p - 2),
                    "PCAschott" = 0:(p - 2),
                    "SIRasymp" = 0:(p - 1),
                    "SIRboot" = 0:(p - 1),
                    stop("Unknown test function")
  )

  # Limit k_range to min_dim and max_dim if necessary
  min_dim <- if (is.null(min_dim)) min(k_range) else max(min(k_range), min_dim)
  max_dim <- if (is.null(max_dim)) p else min(p, max_dim)
  k_range <- k_range[k_range >= min_dim & k_range <= max_dim]
  
  # Default early_stop logic based on search method if not provided
  if (is.null(early_stop)) {
    early_stop <- switch(search,
                         forward = TRUE,
                         backward = TRUE,
                         binary = FALSE)
  }
  
  # Initialize empty progress bar
  #library(progress)
  pb <- NULL
  
  # Function to execute the test for a given k
  test_k <- function(k, tested_ks, pvals, test_results) {
    test_result <- tryCatch({
      method(X, k, ...)
    }, error = function(e) {
      # Print the error message and prompt the user
      message("\nERROR: The test function encountered an issue.")
      message("Please check if the test function is used correctly\nand that all necessary parameters are passed.\n")
      stop(e)  # Stop the function with the error
    })
  
    tested_ks <- c(tested_ks, k)
    pvals <- c(pvals, test_result$p.value)
    test_results[[as.character(k)]] <- test_result
    
    list(
      pval = test_result$p.value,
      tested_ks = tested_ks,
      pvals = pvals,
      test_results = test_results
    )
  }
  
  # Forward search with progress bar
  forward_search <- function() {
    
    tested_ks <- c()
    pvals <- c()
    test_results <- list()
    
    if (early_stop == FALSE) {
      
      pb <- progress_bar$new(
        format = "  Searching [:bar] :percent :elapsed",
        total = length(k_range), clear = FALSE, width = 60,
        show_after = 0
      )
    } 
    
    for (k in k_range) {
      
      res <- test_k(k, tested_ks, pvals, test_results)
      tested_ks <- res$tested_ks
      pvals <- res$pvals
      test_results <- res$test_results
      
      if (!is.null(pb)) pb$tick()
      
    
      if (res$pval > alpha && early_stop){
        
        return(list(est_dim = k, 
                    tested_ks = tested_ks, 
                    pvals = pvals, 
                    test_results = test_results))
        
      }
    }
    
    return(list(
      est_dim = head(tested_ks[pvals > alpha], 1),
      tested_ks = tested_ks,
      pvals = pvals,
      test_results = test_results
    ))
  }
  
  # Backward search with progress bar
  backward_search <- function() {
    
    tested_ks <- pvals <- c()
    test_results <- list()
    
    if (early_stop == FALSE) {
      
      pb <- progress_bar$new(
        format = "  Searching [:bar] :percent :elapsed",
        total = length(k_range), clear = FALSE, width = 60,
        show_after = 0
      )
    } 
    
    for (k in rev(k_range)) {
      
      res <- test_k(k, tested_ks, pvals, test_results)
      tested_ks <- res$tested_ks
      pvals <- res$pvals
      test_results <- res$test_results
      
      if (!is.null(pb)) pb$tick()
      
      if (res$pval <= alpha && early_stop) {
        
        return(list(est_dim = min(k + 1, max(k_range)), 
                    tested_ks = tested_ks,
                    pvals = pvals,
                    test_results = test_results
                    ))
        
      }
    }
    
    return(list(
      est_dim = min(head(tested_ks[pvals <= alpha], 1) + 1, max(k_range)),
      tested_ks = tested_ks,
      pvals = pvals,
      test_results = test_results
    ))
  }
  
  # Binary search
  binary_search <- function() {
    
    tested_ks <- pvals <- c()
    test_results <- list()
    
    if (early_stop == FALSE) {
      
      pb <- progress_bar$new(
        format = "  Searching [:bar] :percent :elapsed",
        total = length(log2(max(k_range))), clear = FALSE, width = 60,
        show_after = 0
        )
      }
    
    low <- min(k_range)
    high <- max(k_range)
    best_k <- NULL
    
    while (low <= high) {
      
      mid <- floor((low + high) / 2)
      res <- test_k(mid, tested_ks, pvals, test_results)
      tested_ks <- res$tested_ks
      pvals <- res$pvals
      test_results <- res$test_results
      
      if (!is.null(pb) && !pb$finished) {
        
        pb$tick()
      
        }
      
      if (res$pval > alpha) {
        
        best_k <- mid
        
        if (early_stop) {
          
          return(list(est_dim = best_k, 
                      tested_ks = tested_ks,
                      pvals = pvals,
                      test_results = test_results))
          
        }
        
        high <- mid - 1 # keep searching left for possibly smaller valid k
        
      } else {
        
        low <- mid + 1
        
      }
    }
    
    return(list(est_dim = best_k,  
                tested_ks = tested_ks,
                pvals = pvals,
                test_results = test_results))
  }
  
  search_result <- switch(search,
                          forward = forward_search(),
                          backward = backward_search(),
                          binary = binary_search())
  
  est_dim <- search_result$est_dim
  tested_ks <- search_result$tested_ks
  pvals <- search_result$pvals
  test_results <- search_result$test_results
  
  
  # Check if no valid k was found
  if (length(est_dim) == 0 || is.null(est_dim)) {
    warning("No value of k found where the null hypothesis is not rejected.")
    
    # Return result with only tested_ks and p_values, but still assign class
    result <- list(
      tested.ks = tested_ks,
      tested.ks.pvals = pvals
    )
  } else {
    # If a valid k was found, include the test results
    result <- c(
      test_results[[as.character(est_dim)]],
      list(
        tested.ks = tested_ks,
        tested.ks.pvals = pvals
      )
    )
  }
  
  # Assign the 'kSearch' class to the result, regardless of the outcome
  class(result) <- c("kSearch", "ictest", "htest")
  
  # Return the result
  cat("\n")
  return(result)
}

#' @title 
#' Summarizing an Object of Class kSearch
#'
#' @description 
#' Summarizes a \code{kSearch} object.
#'
#' @usage 
#' \method{summary}{kSearch}(object, ...)
#'
#' @param object An object of class \code{kSearch}.
#' @param ... Further arguments.
#'
#' @author 
#' Katariina Perkonoja
#'
#' @seealso 
#' \code{\link{print.kSearch}}, \code{\link{kSearch}}
#'
#' @examples 
#' n <- 500
#' S <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n))
#' A <- rorth(4)
#' X <- S %*% t(A)
#' result <- kSearch(X, method = PCAasymp)
#' summary(result)
#'
#' @keywords methods
#' @keywords print
#'
#' @method summary kSearch
#' @exportS3Method summary kSearch
summary.kSearch <- function(object, ...) {
  cat(paste0("kSearch estimator using \"", object$method, "\" for ", object$data.name, "\n"))
  cat("k:", object$k, "\n")
}

#' @title 
#' Printing an Object of Class kSearch
#'
#' @description 
#' Basic printing of an object of class \code{kSearch}. Prints the estimated dimension \code{k} along with
#' all tested values and corresponding p-values.
#'
#' @usage 
#' \method{print}{kSearch}(x, ...)
#'
#' @param x An object of class \code{kSearch}.
#' @param digits Number of significant digits to use when printing p-values. Defaults to \code{getOption("digits")}.
#' @param ... Further arguments.
#'
#' @author 
#' Katariina Perkonoja
#'
#' @seealso 
#' \code{\link{summary.kSearch}}, \code{\link{kSearch}}
#'
#' @examples 
#' n <- 500
#' S <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n))
#' A <- rorth(4)
#' X <- S %*% t(A)
#' result <- kSearch(X, method = PCAasymp)
#' result
#'
#' @keywords methods
#' @keywords print
#'
#' @method print kSearch
#' @exportS3Method print kSearch
print.kSearch <- function(x, digits = getOption("digits"), ...) {
  # Inform about the estimated k being printed
  if (!is.null(x$k)) {
    cat("Results for estimated k = ", x$k, ":\n", sep = "")
    # Call the inherited print method from htest
    NextMethod()
  }
  
  #  Add the additional information for kSearch
  if (!is.null(x$tested.ks) && !is.null(x$tested.ks.pvals)) {
    cat("Tested k values and respective p-values:\n")
    tested_ks <- x$tested.ks
    p_values <- x$tested.ks.pvals
    for (i in seq_along(tested_ks)) {
      cat("k = ", tested_ks[i], ", p-value = ", format.pval(p_values[i], digits = max(1L, digits - 3L)), "\n", sep ="")
    }
  }
  
  cat("\n")
  invisible(x)
}

