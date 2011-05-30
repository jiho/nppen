#
#     Non-Parametric Probabilitic Ecological Niche
#
#   Compute the probability of presence of a given taxon based on 
#   the environmental parameters
#
#   X	matrix or data.frame containing the environmental data at
#   	locations of presence
#   	
#   Y	matrix or data.frame containing the environmental data at
#   	the points where the probability of presence needs to be
#   	predicted (i.e., usually on a grid)
#
#   Beaugrand, G, Lenoir, S, Ibañez, F, and Manté, C
#   A new model to assess the probability of occurrence of a species,
#   based on presence-only data
#   Marine Ecology Progress Series, 424:175--190, 2011
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

nppen <- function(X, Y) {

	# vector sizes
	n = nrow(X)
	m = nrow(Y)

	# compute distance between all points of Y and the centroid of X
	d2 = mahalanobis(Y, mean(X), cov(X))

	# compute probabilities of occurrence of those distances
	# to do this compute the distance between each element of X and the rest of X+each element of Y. For each element of Y count how many of these permutations lead to a distance larger than the one computed.
	p = c()
	for (iy in 1:m) {
		y = Y[iy,]
		d2m = c()
		for (ix in 1:n) {
			Z = rbind(y, X[-ix,])
			d2m[ix] = mahalanobis(X[ix,], mean(Z), cov(Z))
		}
		p[iy] = sum(d2[iy] > d2m)/n
	}
	
	return(p)
}

nppen.ibanez <- function(X, Y) {
	# note about mahalanobis computation
	#   mahalanobis(Z[1,], mean(Z), cov(Z))
	# is equivalent to
	#   Zs = as.data.frame(scale(Z))
	#   k = as.numeric(Zs[1,])
	#   k %*% solve(cov(Zs)) %*% k
	# which is a form close to what Ibanez uses

	# vector sizes
	n = nrow(X)
	m = nrow(Y)

	p = c()
	for (iy in 1:m) {
		# distance from the current element of Y to centroid of X
		Z = rbind(Y[iy,], X)
		Zs = scale(Z)
		k = Zs[1,] - apply(Zs[-1,], 2, mean)
		invCov = solve(cov(Zs))
		e0 = as.numeric(k %*% invCov %*% k)

		d2m = c()
		for (ix in 1:n) {
			# distance from the current element of X and the centoid of the rest of X + the current y
			k = Zs[1+ix,] - apply(Zs[-(1+ix),], 2, mean)
			d2m[ix] = k %*% invCov %*% k	
		}
		p[iy] = sum(e0 > d2m)/n
	}
	return(p)
}
