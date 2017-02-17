#
# Test how much rasterisation helps remove the effect of sampling bias
#
# (c) 2017 Jean-Olivier Irisson, GNU General Public License v3
#

library("tidyverse")

## Simulate data
set.seed(123)

n <- 3000

# simple reference matrix
# = subset of the possibles
Xr <- data.frame(V1=rnorm(n), V2=rnorm(n))

# biased reference matrix
# = subset of the possibles + oversampling in some regions of the state space
Xb <- rbind(Xr, filter(data.frame(V1=rnorm(n), V2=rnorm(n)), V1>0, V2<2.5, V2>-2.5))
n <- 1000

# test matrix
grid <- seq(-2.5, 2.5, length.out=50)
Y <- expand.grid(V1=grid, V2=grid)

# put all data in the same data.frame
p <- ggplot() + geom_point(aes(V1, V2), alpha=0.5) + coord_fixed()
p %+% Xr
p %+% Xb
p %+% Y

devtools::load_all()

## Rasterise ----

XrR <- rasterise(Xr, n=30)
p %+% Xr
p %+% XrR

XbR <- rasterise(Xb, n=30)
p %+% Xb
p %+% XbR


## Test effect of rasterisation on proba ----

Y$r  <- nppen(Xr, Y[,1:2])
Y$rR <- nppen(XrR[,1:2], Y[,1:2])
Y$b  <- nppen(Xb, Y[,1:2])
Y$bR <- nppen(XbR[,1:2], Y[,1:2])

Yt <- gather(Y, "source", "proba", -starts_with("V"))

Xr$source  <- "r"
XrR$source <- "rR"
Xb$source  <- "b"
XbR$source <- "bR"
Xt <- bind_rows(Xr, XrR, Xb, XbR)

ggplot(Xt) + geom_point(aes(V1, V2), alpha=0.1) + geom_vline(aes(xintercept=0), colour="black") + facet_wrap(~source) + coord_fixed()
ggsave("raster_points.png")

ggplot(Yt) + geom_raster(aes(V1, V2, fill=proba)) + geom_contour(aes(V1, V2, z=proba), colour="black") + geom_vline(aes(xintercept=0), colour="black") + facet_wrap(~source) + scale_fill_distiller(palette="Blues") + coord_fixed() + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
ggsave("raster_proba.png")
# -> not perfect but helps a lot
