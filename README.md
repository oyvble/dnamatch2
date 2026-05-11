
# dnamatch2

A fast large scale DNA database search algorithm to search between trace samples (also mixture profiles) and between trace samples and reference samples. 

## Requiring EuroForMix

It is recommended to have euroformix R-package installed (at least v4.0)
before installing dnamatch2. 

## Installation of other required packages

``` r
install.packages(c('gWidgets2tcltk','igraph')) #required
```

## Installation of dnamatch2

Installation directly from source through GitHub (does not require
compilation):

``` r
install.packages("remotes")
remotes::install_github("oyvble/dnamatch2")
```

## Run dnamatch2

``` r
library(dnamatch2);gui() 
```
