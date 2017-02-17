#
# Test the various R implementation and the MATLAB implementation
#
# (c) 2017 Jean-Olivier Irisson, GNU General Public License v3
#

library("tidyverse")

set.seed(1234)

## Simulate data ----

n <- 1000

# reference matrix
X <- data.frame(V1=rnorm(n), V2=rnorm(n))

# test matrix
Y <- data.frame(V1=rnorm(n), V2=rnorm(n))
Y <- Y[1:300,]

# put all data in the same data.frame
p <- ggplot() + geom_point(aes(V1, V2), alpha=0.5) + coord_fixed()
p %+% X
p %+% Y

write_delim(X, path="X.txt")
write_delim(Y, path="Y.txt")

## Compare of implementations ----

devtools::load_all()

# reference implementation vs my matrix computation
system.time(p_ref <- nppen::nppen_ref(X, Y[1:20,]))
#  user  system elapsed
# 3.678   0.138   3.825
system.time(p_mine <- nppen(X, Y[1:20,], fast=FALSE))
#  user  system elapsed
# 2.372   0.167   2.550
all(p_ref == p_mine)
# [1] TRUE
# -> faster and exactly equal

# R implementation vs MATLAB
system.time(Y$p_slow <- nppen(X, Y[,1:2], fast=FALSE))
# read result from MATLAB
Y$p_matlab <- read.table("Yprobs_NPPEN_Matlab.csv")[,1]
all(Y$p_slow == Y$p_matlab)
# [1] TRUE
# -> exactly equal

# exact computation vs fast approximation
# code = cf above
#   user  system elapsed
# 34.285   2.776  37.154
system.time(Y$p_fast <- nppen(X, Y[,1:2], fast=TRUE))
#  user  system elapsed
# 0.052   0.005   0.057

Y$error <- abs(Y$p_slow - Y$p_fast)
Y$percent_error <- Y$error / Y$p_slow * 100
summary(Y$error)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.00000 0.00100 0.00076 0.00100 0.00400
qplot(Y$error, binwidth=0.001)
# -> orders of magnitude faster and relatively accurate when close to the centre of the niche

# what is the source of error and where are the errors
ggplot() +
  geom_point(aes(V1, V2), data=X, alpha=0.1) +
  geom_point(aes(V1, V2, size=error), data=Y, alpha=0.5) +
  scale_size(range=c(0.1,10)) +
  coord_equal()
ggplot() +
  geom_point(aes(V1, V2), data=X, alpha=0.3) +
  geom_point(aes(V1, V2, colour=p_slow, size=percent_error), data=Y, alpha=0.8) +
  scale_colour_distiller(palette="Spectral") +
  coord_equal()
qplot(p_slow, error, data=Y, geom=c("point", "smooth"))
qplot(p_slow, percent_error, data=Y, geom=c("point", "smooth"))
# -> errors are maximal in the middle of the probability range (=when prediction is uncertain) and relative error seems to be higher at the limit of the niche (in environmental space) (i.e. where proba is low => percent error gets higher, mathematically)

# evaluate how error relates to distance from the core of the reference data or from the reference data points
Y$dist_center <- sqrt((Y$V1 - mean(X$V1))^2 + (Y$V2 - mean(X$V2))^2)
Y$dist_closest <- apply(Y[,1:2], 1, function(x) {
  min(sqrt((X$V1 - x[1])^2 + (X$V2 - x[2])^2))
})
Y$dist_close <- apply(Y[,1:2], 1, function(x) {
  mean(sort(sqrt((X$V1 - x[1])^2 + (X$V2 - x[2])^2))[1:10])
})

qplot(dist_center, error, data=Y)
qplot(dist_center, percent_error, data=Y)
qplot(dist_closest, error, data=Y)
qplot(dist_closest, percent_error, data=Y)
qplot(dist_close, error, data=Y)
qplot(dist_close, percent_error, data=Y)
qplot(dist_center+dist_close*5, error, data=Y)
# -> nothing obvious, a bit difficult to predict (but overall, the errors are small)
