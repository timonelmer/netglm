
<!-- README.md is generated from README.Rmd. Please edit that file -->

# netglm

<!-- badges: start -->
<!--![GitHub release (latest by date)](https://img.shields.io/github/v/release/timonelmer/netglm)
![GitHub Release Date](https://img.shields.io/github/release-date/timonelmer/netglm) -->

![GitHub
issues](https://img.shields.io/github/issues-raw/timonelmer/netglm)
![GitHub All
Releases](https://img.shields.io/github/downloads/timonelmer/netglm/total)
<!-- [![Codecov test coverage](https://codecov.io/gh/timonelmer/dnetglm/branch/master/graph/badge.svg)](https://codecov.io/gh/timonelmer/netglm?branch=master) -->

<!-- badges: end -->

The goal of `netglm` is to provide a software which allows to model
(multi-group) network data in a general linear regression framework.

## Installation

You can install the development version of netglm from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("timonelmer/netglm")
```

## Example

This is a basic example which shows you how to estimate a multi-group
QAP model. But first, let us create some test data (inspired by the
example in function `sna::netlm`) we will imagine that we have two
independently measured networks (dependent variables) with four
predictor variables each (i.e. four networks). For example, friendship
networks in two different school classes as dependent variables.

``` r

# create two sets of random networks as independent variables
ivnet1<-sna::rgraph(20,4) 
ivnet2<-sna::rgraph(20,4)

# create dependent variables as functions of the independent variables
dv1<-ivnet1[1,,]+4*ivnet1[2,,]+2*ivnet1[3,,]   # Note that the fourth graph  is unrelated
dv1 <- dv1 + rnorm(400,mean = 1, sd = 1)

dv2 <- 2*ivnet2[1,,]+3*ivnet2[2,,]+3*ivnet2[3,,]
dv2 <- dv2 + rnorm(400,mean = 1, sd = 1)
dvs <- list(dv1, dv2) # put the two dependent variables in a list

iv1 <- list(ivnet1[1,,],ivnet1[2,,],ivnet1[3,,], ivnet1[4,,]) #put all independent variables in a list
iv2 <- list(ivnet2[1,,],ivnet2[2,,],ivnet2[3,,], ivnet2[4,,])
ivs <- list(iv1, iv2) # put the two lists of independent variables in a list
```

Now that we have created some test data, we can estimate a multi-group
QAP model with the `QAP.MG` function:

``` r
library(netglm)
fit <- QAP.MG(dvs, ivs, iv.names = c("intercept",paste0("IV",1:4)), samples = 1000)
#> 
#>  replace diagonal values with NAs
#>  get the observed estimates
#> Loading required package: foreach
#> Loading required package: doParallel
#> Loading required package: iterators
#> Loading required package: parallel
#> 
#>  estimating permuted networks with mode yQAP
```

We can inspect the output as follows:

``` r
fit$output # for fixed effects estimates and their distribution
#>           Estimates p(1sided) abs(p)    adj.d    Exp.V Exp.V.sd  2.5th P
#> intercept   0.97029     0.000  0.000 21.20541  4.72336  0.17699  4.35717
#> IV1         1.39495     1.000  0.000  8.00505 -0.00250  0.17457 -0.35509
#> IV2         3.55058     1.000  0.000 20.15557 -0.00841  0.17658 -0.35319
#> IV3         2.55761     1.000  0.000 13.94993  0.01096  0.18256 -0.34899
#> IV4         0.10572     0.731  0.269  0.60746 -0.00525  0.18267 -0.37624
#>           97.5th P significance
#> intercept  5.05518          ***
#> IV1        0.33972          ***
#> IV2        0.33642          ***
#> IV3        0.34768          ***
#> IV4        0.36077

fit$r.squared # for r.squared measures
#>     r.squared adj.r.squared 
#>     0.8110922     0.8100914
```
