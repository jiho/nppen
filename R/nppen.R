#
#     Non-Parametric Probabilistic Ecological Niche
#
#     Beaugrand, G, Lenoir, S, Ibañez, F, and Manté, C
#     A new model to assess the probability of occurrence of a species,
#     based on presence-only data
#     Marine Ecology Progress Series, 424:175--190, 2011
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

nppen <- function(X, Y, fast=TRUE)
#
#   Compute the probability of presence of a given taxon based on the
#   environmental descriptors
#
#   X   matrix or data.frame containing the environmental data at locations of
#       presence (possibly rasterized)
#
#   Y   matrix or data.frame containing the environmental data at the points
#       where the probability of presence needs to be predicted (i.e., usually
#       on a grid)
#
#   fast    when TRUE, a single total covariance matrix is used for each
#           element of Y, instead of one matrix per permutation (i.e. per line
#           of X). This is immensely faster and leads to only marginal (<2%)
#           changes in the predicted probabilities as long as X is large and
#           points of Y as distributed within the range covered by X
#
{
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
    invCov <- solve(cov(X))
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
            inCov <- solve(cov(Z))
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
                invCov <- solve(cov(Z))
                d2m[ix] <- c(Xc %*% invCov %*% Xc)
                # NB: since X[ix,] is a vector, we don't have to worry about the transpositions and the rowSums
            }

            # compute the probability
            p[iy] <- sum(d2m > d2[iy]) / n
        }

    }

    return(p)
}

# Reference implementation (slow but accurate)
nppen.ref <- function(X, Y) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    n <- nrow(X)
    m <- nrow(Y)

    # compute reference distance for each point of Y
    d2 <- mahalanobis(Y, colMeans(X), cov(X))

    # compute the probabilities of occurrence of those distances
    p <- c()
    for (iy in 1:m) {
        y <- Y[iy,]
        d2m <- c()
        for (ix in 1:n) {
            Z <- rbind(y, X[-ix,])
            d2m[ix] <- mahalanobis(X[ix,], colMeans(Z), cov(Z))
        }
        p[iy] <- sum(d2m > d2[iy]) / n
    }

    return(p)
}
