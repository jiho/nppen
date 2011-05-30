#
#     Test the NPPEN method from Beaugrand et al 2011
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

## Simulated data                                           {

# reference matrix
n = 100
X = data.frame(
	T = runif(n)*5,
	S = 36+runif(n),
	bathy = 10+runif(n)*300
)

# prediction matrix
m = 50
Y = data.frame(
	T = 1+runif(m)*3,
	S = 36+runif(m)/2,
	bathy = 20+runif(m)*200
)

source("nppen.R")

system.time(nppen(X, Y))
system.time(nppen.ibanez(X, Y))
system.time(nppen.vector(X, Y))

nppen(X, Y)
nppen.ibanez(X, Y)
nppen.vector(X, Y)


# }
