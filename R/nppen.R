#' Non-Parametric Probabilistic Ecological Niche
#'
#' Compute the probability of presence of a taxon based on the environment.
#
#' @param X matrix or data.frame containing the environmental data at locations of presence (possibly rasterized using \code{\link{rasterize}}).
#' @param Y matrix or data.frame containing the environmental data at the points where the probability of presence needs to be predicted (i.e., usually on a grid).
#' @param fast when TRUE, a single total covariance matrix is used for each element of Y, instead of one matrix per permutation (i.e. per line of X). This is immensely faster and leads to only marginal (<2\%) changes in the predicted probabilities as long as X is large and points of Y are distributed within the range covered by X (which they should be for the model to make sense anyway).
#'
#' @return A vector of probabilities of length \code{nrow(Y)}.
#'
#' @references
#' Beaugrand, G., Lenoir, S., Ibañez, F., and Manté, C. (2011) _A new model to assess the probability of occurrence of a species, based on presence-only data_. Marine Ecology Progress Series, *424*, 175-190. \url{http://www.int-res.com/abstracts/meps/v424/p175-190/}
#'
#' @encoding UTF-8
#' @export
#'
#' @examples
#' # define environmental data
#' set.seed(1)
#' X <- data.frame(temp=rnorm(200, 15, 3), sal=rnorm(200, 37, 1))
#' # define target points: one well inside the niche, one on the border
#' Y <- data.frame(temp=c(15, 18.7), sal=c(37, 38.7))
#' # represent both in environmental space
#' plot(X$temp, X$sal)
#' points(Y$temp, Y$sal, col="red")
#' # compute probability of presence
#' nppen(X, Y, fast=FALSE)
#' nppen(X, Y, fast=TRUE)
#' # higher probability for the first point, as expected
#'
#' # rasterize X to make is smaller (hence faster to compute) and avoid spatial
#' # bias in the underlying data (e.g. preferential sampling in a give area)
#' X_binned <- rasterize(X)
#' nrow(X)
#' nrow(X_binned)
#' plot(X_binned$temp, X_binned$sal)
#' points(Y$temp, Y$sal, col="red")
#' # compute probability of presence
#' nppen(X_binned[,-ncol(X_binned)], Y, fast=FALSE) # NB: remove column `n`
#' # simply binning loose some information regarding where the centre of the
#' # niche is and the niche appear more spread (hence higher probability for
#' # the second point here).
#'
#' # only keep common observations (= observed more than once)
#' X_binned_common <- X_binned[X_binned$n > 1,]
#' points(X_binned_common$temp, X_binned_common$sal, pch=16)
#' nppen(X_binned_common[,-ncol(X_binned)], Y, fast=FALSE)
#' # reducing to common observations helps better represent the initial niche
nppen <- function(X, Y, fast=TRUE) {
  # convert to matrices, for speed
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  # sizes
  n <- nrow(X)
  m <- nrow(Y)

  # compute reference distance for each point of Y:
  #   the actual Mahalanobis distance to X
  # NB: mahalanobis() has
  #   sweep(Y, 2, colMeans(X))
  # but this
  Yc <- t(t(Y)-colMeans(X))
  # is faster
  invCov <- solve(stats::cov(X))
  d2 <- rowSums((Yc %*% invCov) * Yc)

  # compute the probabilities of occurrence of those distances:
  #   for a given element of Y (y_j), swap it with an element of X (x_i)
  #   and compute the distance from x_i to X - x_i + y_j
  #   then count how many of these permutations lead to a distance
  #   larger than distance from y_j to X
  p <- rep(-1, m)  # NB: pre-allocating the vector is a bit faster
  if (fast) {

    for (iy in 1:m) {
      y <- Y[iy,]

      # use a global matrix which is y + X instead of swapping y for
      # each element of X
      Z <- rbind(y, X)
      inCov <- solve(stats::cov(Z))
      Xc <- t(t(X) - colMeans(Z))
      d2m <- rowSums((Xc %*% invCov) * Xc)

      # compute the probability
      p[iy] <- sum(d2m > d2[iy]) / n
      # NB: it is actually faster to do this in the loop than to store
      #     all values in a matrix and compute the probabilities on the
      #     matrix. I don't understand why, but...
    }

  } else {

    for (iy in 1:m) {
      y <- Y[iy,]

      d2m <- rep(-1, n)
      for (ix in 1:n) {
        # create the matrix X - x_i + y_j
        Z <- rbind(y, X[-ix,])

        # compute mahalanobis' distance
        Xc <- X[ix,] - colMeans(Z)
        invCov <- solve(stats::cov(Z))
        d2m[ix] <- c(Xc %*% invCov %*% Xc)
        # NB: since X[ix,] is a vector, we don't have to worry about the transpositions and the rowSums
      }

      # compute the probability
      p[iy] <- sum(d2m > d2[iy]) / n
    }

  }

  return(p)
}

# # Reference implementation (slow but accurate)
# nppen_ref <- function(X, Y) {
#   X <- as.matrix(X)
#   Y <- as.matrix(Y)
#
#   n <- nrow(X)
#   m <- nrow(Y)
#
#   # compute reference distance for each point of Y
#   d2 <- stats::mahalanobis(Y, colMeans(X), stats::cov(X))
#
#   # compute the probabilities of occurrence of those distances
#   p <- c()
#   for (iy in 1:m) {
#     y <- Y[iy,]
#     d2m <- c()
#     for (ix in 1:n) {
#       Z <- rbind(y, X[-ix,])
#       d2m[ix] <- stats::mahalanobis(X[ix,], colMeans(Z), stats::cov(Z))
#     }
#     p[iy] <- sum(d2m > d2[iy]) / n
#   }
#
#   return(p)
# }
