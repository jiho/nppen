# Non-Parametric Probabilistic Ecological Niche

Basic implementation of:

Beaugrand, G., Lenoir, S., Ibañez, F., and Manté, C. (2011) [_A new model to assess the probability of occurrence of a species, based on presence-only data_](http://www.int-res.com/abstracts/meps/v424/p175-190/). Marine Ecology Progress Series, *424*, 175-190.

in R.

This implementation contains an approximation that makes it very fast in the common conditions of application and is turned on by default.

The code is not on CRAN (and probably will never be, by itself, since it is just two functions). To install:

```r
devtools::install_github("jiho/nppen")
# NB: install package `devtools` if you do not have it already.
```


