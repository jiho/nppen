#
#     Non-Parametric Probabilitic Ecological Niche
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

rasterize <- function(X, n=10, precisions=sapply(lapply(X,range),diff)/n)
#
#   Reduce the precision of all columns of a data.frame to bins
#   and count the number of occurence of each combination of bin values
#   This is a bit like reducing the precise information on locations
#   in a 2D plane to pixels of a given grey level, hence the name
#
#   X           matrix of data
#
#   n           scalar, number of bins for all columns of the data
#
#   precisions  vector, overrides n, the precision which is used to
#               cut each column of the data into bins
#
{
    # Checks
    if (length(precisions) != ncol(X)) {
        stop("The vector of precisions does not have as many elements as there are columns in the data X")
    }

    # Round columns to the given precision
    require("plyr")
    for (j in 1:ncol(X)) {
        X[,j] = round_any(X[,j], accuracy=precisions[j])
    }

   # Compute how many times each combination of values is present
   X = count(X, vars=names(X))

   return(X)
}

nppen <- function(X, Y)
#
#   Compute the probability of presence of a given taxon based on
#   the environmental descriptors
#
#   X	matrix or data.frame containing the environmental data at
#   	locations of presence (possibly rasterized)
#
#   Y	matrix or data.frame containing the environmental data at
#   	the points where the probability of presence needs to be
#   	predicted (i.e., usually on a grid)
#
{

	# vector sizes
	n = nrow(X)
	m = nrow(Y)

	# compute reference distance for each point of Y:
	# the actual Mahalanobis distance to X
	d2 = mahalanobis(Y, colMeans(X), cov(X))

	# compute the probabilities of occurrence of those distances
	# for each element of Y (y_j), swap it with an element of X (x_i)
	# and compute the distance from x_i to X - x_i + y_j
	# then count how many of these permutations lead to a distance
	# larger than the reference one
	p = c()
	for (iy in 1:m) {
		y = Y[iy,]
		d2m = c()
		for (ix in 1:n) {
			Z = rbind(y, X[-ix,])
			d2m[ix] = mahalanobis(X[ix,], colMeans(Z), cov(Z))
		}
		p[iy] = sum(d2m > d2[iy]) / n
	}

	return(p)
}
