#' Bin environmental data
#'
#' Bin each environment variable and count the number of observation of all combinations of bins. Is then used to filter out rare combinations of variables.
#'
#' @param X matrix (or data.frame) of environmental data.
#' @param n number of bins for each variable in the data (just one number).
#' @param precisions the precision used to cut each column of the data into bins (a vector, as long as the number of columns of X); overrides \code{n}.
#'
#' @return A data.frame with all observed combinations of binned variables and a column \code{n} containing the number of observations of this combination. Each bin is identified by the average value of each variable in the bin. Therefore, it changes from bin to bin, i.e. even with a precision of 0.1, binned values are not all 0.05, 0.15, etc.
#'
#' @export
#'
#' @examples
#' X <- data.frame(a=runif(100), b=runif(100)*10)
#' rasterize(X, n=2)
#' rasterize(X, precisions=c(0.5, 5))
rasterize <- function(X, n=10, precisions=NULL) {
  X_binned <- X
  
  if (is.null(precisions)) {
    # safety check
    if (n < 2) {
      stop("n should be >= 2.")
    }
    # cut based on n
    for (j in 1:ncol(X)) {
      X_binned[,j] <- cut(X[,j], breaks=n)
    }
    
  } else {
    # safety check
    if (length(precisions) != ncol(X)) {
      stop("The vector of precisions has ", length(precisions), " elements when if should have as many elements as there are columns in X (", ncol(X), ").")
    }
    # round columns to the given precision
    for (j in 1:ncol(X)) {
      X_binned[,j] <- round_any(X[,j], accuracy=precisions[j])
    }
  }
    
  # compute average value inside each bin
  names(X_binned) <- paste0("b", names(X_binned))
  XX <- cbind(X_binned, X)
  bins <- dplyr::summarise_each_(dplyr::group_by_(XX, .dots=names(X_binned)), funs="mean", vars=names(X)) 

  # compute number of occurrences in each bin
  counts <- dplyr::count_(X_binned, vars=names(X_binned))

  # combine bin values and counts
  C <- data.frame(bins[,names(X)], n=counts$n)
  
  return(C)
}

# allow British spelling
rasterise <- rasterize

# duplicate round_any function from plyr to avoid depending on it
round_any <- function(x, accuracy) {
  round(x/accuracy) * accuracy
}
