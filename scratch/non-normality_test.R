#
# Test how sensitive NPPEN is to non multinormal data (in particular to multimodal data)
#
# (c) 2017 Jean-Olivier Irisson, GNU General Public License v3
#

library("tidyverse")

## Simulate data
set.seed(123)

n <- 2000

# reference, multinormal data
Xn <- data.frame(V1=rnorm(n), V2=rnorm(n))
Xn$source <- "n"

# skewed data
p <- (0:100)/100
qplot(p, dbeta(p, 8, 4))
Xs <- data.frame(V1=(rbeta(n, 8, 4)-0.5)*6, V2=(rbeta(n, 9, 2)-0.5)*6)
Xs$source <- "s"

# multimodal data
sd <- 0.7
m <- 1
Xm <- rbind(
  data.frame(V1=rnorm(n, -m, sd), V2=rnorm(n, -m, sd)),
  data.frame(V1=rnorm(n,  m, sd), V2=rnorm(n,  m, sd))
)
Xm$source <- "m"

# test matrix
grid <- seq(-4, 4, length.out=60)
Y <- expand.grid(V1=grid, V2=grid)

# put all data in the same data.frame
p <- ggplot() + geom_point(aes(V1, V2), alpha=0.5) + coord_fixed()
p %+% Xn
p %+% Xs
p %+% Xm
p %+% Y

devtools::load_all()

## Map niches ----

Y$n <- nppen(Xn[,1:2], Y[,1:2])
Y$s <- nppen(Xs[,1:2], Y[,1:2])
Y$m <- nppen(Xm[,1:2], Y[,1:2])

Yt <- gather(Y, "source", "proba", -starts_with("V"))

Xt <- rbind(Xn, Xs, Xm)

ggplot(Yt) +
  geom_raster(aes(V1, V2, fill=proba)) +
  geom_contour(aes(V1, V2, z=proba), colour="black") +
  geom_point(aes(V1, V2), data=Xt, size=0.2) +
  facet_wrap(~source) + scale_fill_distiller(palette="Blues") +
  coord_fixed() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
ggplot() +
  geom_density2d(aes(V1, V2), data=Xt, colour="red", n=200, h=c(1.5, 1.5)) +
  geom_contour(aes(V1, V2, z=proba), data=Yt, colour="black") +
  facet_wrap(~source)
  coord_fixed() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
# -> Multimodal is horrid
#    Skewed is made symmetrical, not ideal
#    Normal is re-centered which is actually better than the estimation from the data


