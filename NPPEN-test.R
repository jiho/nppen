#
#     Test the NPPEN method from Beaugrand et al 2011
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

pdf("nppen_plots.pdf")

set.seed(123)

# Simulate data
# reference matrix
n = 100
X = data.frame(
	T = rnorm(n, mean=5),
	S = 36+rnorm(n)#,
    # bathy = rnorm(n, mean=300, sd=30)
)

# prediction matrix
m = 50
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
qplot(T, S, colour=source, data=Z)
# -> OK, Y is distributed within the distribution of X

# Test rasterization
source("nppen.R")
qplot(T, S, data=X)
last_plot() %+% rasterize(X, n=10)

# Compute probabilities of occurence with nppen (with rasterization or not)
Xr <- rasterize(X, n=20)
Xr <- Xr[,-ncol(Xr)]
system.time(p <- nppen(X, Y))
system.time(pr <- nppen(Xr, Y))



Yp <- data.frame(Y, p=p, pr=pr)
# verify that those match the distribution in X
ggplot(mapping=aes(x=T, y=S)) + geom_point(data=Xr) + geom_point(aes(colour=p), data=Yp)
# -> OK that works

dev.off()

# write the results for comparison with MATLAB
write.table(X, file="X.txt", sep="\t", row.names=F)
write.table(Xr, file="X_raster.txt", sep="\t", row.names=F)
write.table(Y, file="Y.txt", sep="\t", row.names=F)
write.table(Yp, file="Y_pred.txt", sep="\t", row.names=F)
