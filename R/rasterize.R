#' Bin environmental data
#'
#' Bin each environment variable and count the number of observation of all combinations of bins. Is then used to filter out rare combinations of variables.
#'
#' @param X matrix (or data.frame) of environmental data.
#' @param n number of bins for each variable in the data (either one number or a vector with one number per column of \code{X}).
#' @param precision the precision used to cut each column of the data into bins (either one number or a vector with one number per column of \code{X}); overrides \code{n}.
#' @param breaks list of vectors of cut points for each variable (should have as many elements as columns in \code{X}); overrides \code{n} and \code{precision}.
#'
#' @return A data.frame with all observed combinations of binned variables and a column \code{n} containing the number of observations of this combination. Each bin is identified by the average value of each variable in the bin. Therefore, it changes from bin to bin, i.e. even with a precision of 0.1, binned values are not all 0.05, 0.15, etc.
#'
#' @export
#'
#' @examples
#' X <- data.frame(a=runif(100), b=runif(100)*10)
#' rasterize(X, n=2)
#' rasterize(X, n=c(2,4))
#' rasterize(X, precision=1)
#' rasterize(X, precision=c(0.5, 5))
#' rasterize(X, breaks=list(c(0, 0.5, 1), c(0, 2, 4, 10)))
rasterize <- function(X, n=10, precision=NULL, breaks=NULL) {
  X_binned <- X
  nX <- ncol(X)

  if (!is.null(breaks)) {
    # safety check
    if (length(breaks) != nX) {
      stop("`breaks` has ", length(breaks), " elements. It should have as many as there are columns in X (", nX, ").")
    }
    # cut X based on breaks
    for (j in 1:nX) {
      X_binned[,j] <- cut(X[,j], breaks=breaks[[j]])
    }

  } else if (!is.null(precision)) {
    # safety check
    if (length(precision) == 1) {
      precision <- rep(precision, times=nX)
    } else if (length(precision) != nX) {
      stop("`precision` has ", length(precision), " elements. It should have as many as there are columns in X (", nX, ").")
    }
    # round X to the given precision
    for (j in 1:nX) {
      X_binned[,j] <- round_any(X[,j], accuracy=precision[j])
    }

  } else {
    # safety check
    if (length(n) == 1) {
      n <- rep(n, times=nX)
    } else if (length(n) != nX) {
      stop("`n` has ", length(n), " elements. It should have as many as there are columns in X (", nX, ").")
    }
    if (any(n < 2)) {
      stop("`n` should always be >= 2.")
    }
    # cut X based on n
    for (j in 1:nX) {
      X_binned[,j] <- cut(X[,j], breaks=n[j])
    }
  }

  # compute average value inside each bin
  Xb <- stats::aggregate(X, by=X_binned, FUN="mean", na.rm=TRUE)

  # compute number of occurrences in each bin
  Xn <- as.data.frame(table(X_binned))
  # remove bins with 0 occurences which are not counted in the aggregate call above
  Xn <- Xn[which(Xn$Freq > 0),]

  # combine bin values and counts
  C <- data.frame(Xb[,-(1:ncol(X))], n=Xn$Freq)

  return(C)
}

# allow British spelling
rasterise <- rasterize

# duplicate round_any function from plyr to avoid depending on it
round_any <- function(x, accuracy) {
  round(x/accuracy) * accuracy
}
