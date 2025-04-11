#' @title 
#' Relevant Component Estimation via Iterative Search
#'
#' @description 
#' Performs dimension estimation using various search strategies (forward, backward, binary)
#' and hypothesis testing methods that evaluate whether the true dimension \eqn{d \leq k}.
#' The search continues until a p-value exceeds the specified alpha threshold, which indicates
#' that the null hypothesis is not rejected. The search can optionally continue to test all values, as there may not be a global optimum 
#' in the tested range.
#'
#' @param X A data matrix with \eqn{p>1} columns.
#' @param test_func A function that performs a hypothesis test given \code{X} and dimension \code{k}. 
#' For more details on supported tests, see the
#' "Details" section below.
#' 
#' @details
#' This function is designed to work with the following tests:
#' \itemize{
#'   \item \code{\link[ICtest]{FOBIasymp}}, \code{\link[ICtest]{FOBIboot}}
#'   \item \code{\link[ICtest]{ICSboot}}, 
#'   \item \code{\link[ICtest]{NGPPsim}}
#'   \item \code{\link[ICtest]{PCAasymp}}, \code{\link[ICtest]{PCAboot}}, 
#'   \item \code{\link[ICtest]{PCAschott}}, \item \code{\link[ICtest]{SIRasymp}}, 
#'   \item \code{\link[ICtest]{SIRboot}}
#' }
#' These tests evaluate the null hypothesis that the true dimension \eqn{d \leq k}, and return
#' a p-value indicating whether the null is rejected or not.
#'
#' @param method A character string specifying the search strategy. 
#'    Options are \code{forward} (default), \code{backward} and \code{binary}. For more details on search techniques, see the
#' "Details" section below.
#' 
#' @details
#' The search methods work as follows:
#' \itemize{
#'   \item \code{"forward"} (default) performs a linear search starting from the smallest possible value of \code{k} and incrementing upward.
#'   \item \code{"backward"} performs a linear search starting from the largest possible value of \code{k} and decrementing downward.
#'   \item \code{"binary"} first splits the possible range of \code{k}s in half and performs a binary search in the lower half. If no suitable \code{k} is found in the lower half, it continues the search in the upper half. This method may skip testing some possible \code{k} values, but it is typically faster than exhaustive linear search. The result from the binary search can be further refined using \code{"forward"} or \code{"backward"} search, if desired. Note that with this method the progress bar does not finish at 100 \%.
#'   }
#' 
#' @param alpha A numeric significance level (default = 0.05) used as the threshold for accepting the null hypothesis.
#' @param early_stop Logical. If \code{TRUE}, the search stops at the first \code{k} where the null hypothesis cannot be rejected
#'   (p-value > \code{alpha}). Note that in this case, the progress bar does not necessarily finish at 100\%. If \code{FALSE}, it continues through all plausible values of \code{k}.
#' @param max_dim Optional integer to limit the maximum dimension to search (defaults to \code{ncol(X)}).
#' @param ... Additional arguments passed to the testing function \code{test_func}.
#' 
#' @return A list with the following components:
#' \describe{
#'   \item{\code{test_name}}{A character string giving the name of the test function used.}
#'   \item{\code{estimated_k}}{An integer indicating the estimated number of relevant components \code{k} where the null hypothesis \eqn{d \leq k} is first accepted (or the smallest if \code{early_stop = FALSE}).}
#'   \item{\code{test_result_est_k}}{The test result for the estimated number of relevant components.}
#'   \item{\code{tested_ks}}{An integer vector of all tested values of \code{k} during the search.}
#'   \item{\code{p_values}}{A numeric vector of p-values corresponding to each tested \code{k}.}
#' }
#'
#' @examples
#' # Example 1: Using FOBIboot for forward search
#' set.seed(123)
#' X <- matrix(rnorm(1000), ncol = 10)  # Create a random data matrix
#' result <- kSearch(
#'   X = X, 
#'   test_func = FOBIboot, 
#'   alpha = 0.05, 
#'   method = "forward", 
#'   early_stop = TRUE, 
#'   n.boot = 1000
#' )
#' print(result)
#' 
#' # Example 2: Using PCAasymp for forward search without early stop
#' set.seed(1)
#' n <- 200
#' S <- cbind(rnorm(n, sd = 2), rnorm(n, sd = 1.5), rnorm(n), rnorm(n), rnorm(n))
#' A <- rorth(5)
#' X <- S %*% t(A)
#' result <- kSearch(
#'   X = X, 
#'   test_func = PCAasymp, 
#'   alpha = 0.05, 
#'   method = "forward", 
#'   early_stop = FALSE, 
#'   n.boot = 1000
#' )
#' print(result)
#'
#' @export

kSearch <- function(X, 
                    test_func, 
                    method = c("forward", "backward", "binary"), 
                    alpha = 0.05, 
                    early_stop = TRUE, 
                    max_dim = NULL,
                    ...) {
  
  # Match the method argument
  method <- match.arg(method)
  
  # Number of columns (dimensions)
  p <- ncol(X)
  max_dim <- if (is.null(max_dim)) p else min(p, max_dim)
  
  # Define the ranges for k based on the test_func
  k_range <- switch(deparse(substitute(test_func)),
                    "FOBIasymp" = 0:(p - 1),
                    "FOBIboot" = 0:(p - 1),
                    "ICSboot" = 1:(p - 2),
                    # NGPPsim should be checked, help says k <= p, but
                    # it throws an error so using p-1 for now:
                    # while (crit > eps) { : missing value where TRUE/FALSE needed
                    "NGPPsim" = 1:(p-1),
                    "PCAasymp" = 0:(p - 2),
                    "PCAboot" = 0:(p - 2),
                    "PCAschott" = 0:(p - 2),
                    "SIRasymp" = 0:(p - 1),
                    "SIRboot" = 0:(p - 1),
                    stop("Unknown test function")
  )

  # Limit k_range to max_dim if necessary
  k_range <- k_range[k_range <= max_dim]  
  
  # Load the necessary library for progress bar
  library(progress)
  
  pb <- progress_bar$new(
    format = "  Searching [:bar] :percent :elapsed",
    total = length(k_range), clear = FALSE, width = 60,
    show_after = 0
  )
  
  pb$tick(0)
  flush.console()
  
  
  # Initialize variables to store the results
  tested_ks <- c()
  pvals <- c()
  found_k <- NULL
  test_results <- list()
  
  # Function to execute the test for a given k
  test_k <- function(k) {
    test_result <- tryCatch({
      test_func(X, k, ...)
    }, error = function(e) {
      # Print the error message and prompt the user
      message("\nERROR: The test function encountered an issue.")
      message("Please check if the test function is used correctly\nand that all necessary parameters are passed.\n")
      stop(e)  # Stop the function with the error
    })
    
    tested_ks <<- c(tested_ks, k)
    pvals <<- c(pvals, test_result$p.value)
    
    # Store the test result for the last valid dimension
    if (test_result$p.value > alpha) {
      test_results[[as.character(k)]] <<- test_result
    }
    
    return(test_result$p.value)
  }
  
  # Forward search with progress bar
  forward_search <- function() {

    for (k in k_range) {
      pval <- test_k(k)
      pb$tick()  # Update progress bar
      if (pval > alpha && early_stop) {
        return(k)
      }
    }
    return(head(tested_ks[pvals > alpha], 1))  # First accepted if no early stop
  }
  
  # Backward search with progress bar
  backward_search <- function() {

    for (k in rev(k_range)) {
      pval <- test_k(k)
      pb$tick()  # Update progress bar
      if (pval > alpha && early_stop) {
        return(k)
      }
    }
    return(tail(tested_ks[pvals > alpha], 1))  # Last accepted if no early stop
  }
  
  # Binary search with progress bar
  binary_search <- function() {
    low <- min(k_range)
    high <- max(k_range)
    best_k <- NULL

    # Helper function to search a range
    search_half <- function(low, high) {
      local_best <- NULL
      while (low <= high) {
        mid <- floor((low + high) / 2)
        pval <- test_k(mid)
        pb$tick()
        if (pval > alpha) {
          local_best <- min(mid, local_best)
          if (early_stop){
            return(local_best)
          } 
          high <- mid - 1
        } else {
          low <- mid + 1
        }
      }
      return(local_best)
    }
    
    # Split k_range into halves
    mid_point <- floor((min(k_range) + max(k_range)) / 2)
    
    # First prefer searching lower half
    lower_k <- search_half(min(k_range), mid_point)
    
    if (!is.null(lower_k)) {
      best_k <- lower_k
    } else {
      # Only search upper half if lower half had no valid k
      upper_k <- search_half(mid_point + 1, max(k_range))
      if (!is.null(upper_k)) {
        best_k <- upper_k
        }
    }
    
    return(best_k)
  }
  
  # Execute the chosen search method
  est_dim <- switch(method,
                    forward = forward_search(),
                    backward = backward_search(),
                    binary = binary_search()
  )
  
  # Line break before returning result
  cat("\n")  
  return(list(
    test_name = deparse(substitute(test_func)),
    estimated_k = est_dim,
    test_result_est_k = test_results[[as.character(est_dim)]],
    tested_ks = tested_ks,
    p_values = pvals
  ))
  
}
