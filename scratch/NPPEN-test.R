#
#     Test the NPPEN method from Beaugrand et al 2011
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

pdf("nppen_plots.pdf")

set.seed(1234)

## Simulate data                                            {

# reference matrix
n = 1000
X = data.frame(
	T = rnorm(n, mean=5),
	S = 36+rnorm(n)#,
    # bathy = rnorm(n, mean=300, sd=30)
)
# add some data in a given location of the state space
X2 = data.frame(
	T = rnorm(n*3, mean=5),
	S = 36+rnorm(n)#,
    # bathy = rnorm(n, mean=300, sd=30)
)
X2 <- X2[X2$T>5.5 & X2$T<6.5 & X2$S<37 & X2$S>34.5,]
X <- rbind(X, X2)

# prediction matrix
m = 200
Y = data.frame(
	T = rnorm(m, mean=5, sd=0.4),
	S = 36+rnorm(m, sd=0.5)#,
    # bathy = rnorm(m, mean=280, sd=10)
)

# Visualize the relative range of X and Y
library("ggplot2")
library("reshape")
# put all data in the same data.frame
Z <- rbind(data.frame(X, source="X"), data.frame(Y, source="Y"))
ggplot(Z) + geom_point(aes(x=T, y=S, colour=source), alpha=0.5)
# -> OK, Y is distributed within the distribution of X

# }


source("nppen.R")

## Test rasterization                                       {

Xr <- rasterize(X, precision=c(0.2, 0.3))
(nrow(Xr))

qplot(T, S, data=X)
last_plot() %+% Xr
qplot(T, S, size=freq, data=Xr)


# }

## Compute probabilities of occurence with nppen            {

# reference implementation
system.time(pref <- nppen.ref(Xr[,-ncol(Xr)], Y))

# my own matrix computation
system.time(pm <- nppen(Xr[,-ncol(Xr)], Y, fast=FALSE))
# ~ 30% less time

# the fast, possibly inacurate, way
r <- 10
system.time(for(i in 1:r){ p <- nppen(Xr[,-ncol(Xr)], Y) })
# wow, faaaast

sum(pref != pm)
# -> OK, my implementation is equivalent to the reference one
sum(pref != p)
# - > slight differences

# difference between the reference implementation and the fast one
error <- (abs(pref - p) / pref) * 100
error[is.nan(error)] <- 0   # in case of division by 0
summary(error)


# plot the results
Yp <- data.frame(Y, p=p, pref=pref, error)

# plot probas
pX <- ggplot(mapping=aes(x=T, y=S)) + geom_point(data=Xr, size=1)
pX + geom_point(aes(colour=pref), data=Yp, size=2, alpha=0.8)
pX + geom_point(aes(colour=p), data=Yp, size=2, alpha=0.8)
# -> OK that seems to be logical

# plot the error
pX + geom_point(aes(colour=abs(pref-p), alpha=abs(pref-p)), data=Yp, size=2)
pX + geom_point(aes(colour=error, alpha=error), data=Yp, size=2) + scale_colour_continuous("% error") + scale_alpha_continuous("% error")
# -> maximal error for the outer most points


# estimate the confidence in the interpolation
#   use the mahalanobis distance from the full (non rasterized) could
#   indeed, if some values are observed frequently, they should count more in
#   the proba estimation but they should boost our confidence in the result
confid <- sqrt(mahalanobis(Y, colMeans(X), cov(X)))
# rescale in the correct direction, between 0 and 1
confid <- (max(confid) - confid) / max(confid)

ggplot(mapping=aes(x=T, y=S)) + geom_point(data=X, alpha=0.2)  + geom_density2d(data=X, alpha=0.7, colour="black") + geom_point(aes(colour=confid), data=Yp, size=2, alpha=0.9) + scale_colour_continuous("confid")
# replot p on the same graph to compare
ggplot(mapping=aes(x=T, y=S)) + geom_point(data=X, alpha=0.2)  + geom_density2d(data=X, alpha=0.7, colour="black") + geom_point(aes(colour=p), data=Yp, size=2, alpha=0.9) + scale_colour_continuous("p")
# -> p and confid are a bit different, this is what's expected


# estimate the confidence in the interpolation
#   another way is to use the number of observation per bin in X as a measure
#   of confidence
ggplot(mapping=aes(x=T, y=S)) +
    geom_tile(aes(fill=freq), data=Xr) +
    scale_fill_continuous(high="#B2FF7E", low="#FB7B80") +
    geom_point(aes(colour=p), data=Yp, size=2, alpha=0.9) +
    scale_colour_continuous(high="black", low="white") +
    opts(panel.background=theme_blank(), panel.grid.minor=theme_blank(), panel.grid.major=theme_blank(), panel.border=theme_rect(colour="grey50"))

# }

dev.off()

# # Use parallel computation to speed things up
# # split X in pieces
# # nPieces <- 4
# # Y$piece <- as.numeric(cut(1:nrow(Y), breaks=quantile(1:nrow(Y), seq(0,1,length=nPieces+1)), include.lowest=T))
# # unique(Y$piece)
# # compute
# library("multicore")
# library("doMC")
# library("foreach")
# registerDoMC(cores=4)
# system.time(prp <- foreach(iter(Y, by="row", chunksize=1000), .combine="c") %dopar% nppen(Xr, Y))


# write the results for comparison with MATLAB
write.table(X, file="X.txt", sep="\t", row.names=F)
write.table(Xr, file="X_raster.txt", sep="\t", row.names=F)
write.table(Y, file="Y.txt", sep="\t", row.names=F)
write.table(Yp, file="Y_pred.txt", sep="\t", row.names=F)
