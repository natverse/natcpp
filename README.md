
# natcpp

<!-- badges: start -->
[![natverse](https://img.shields.io/badge/natverse-Part%20of%20the%20natverse-a241b6)](https://natverse.org)
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](https://natverse.org/natcpp/reference/)
[![R-CMD-check](https://github.com/natverse/natcpp/workflows/R-CMD-check/badge.svg)](https://github.com/natverse/natcpp/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of *natcpp* is to provide accelerated routines through compiled C++ 
code for basic functions to support the [natverse](https://natverse.org) computational neuroanatomy ecosystem.
*natcpp* is not designed to be used by end users but instead the 
[nat](https://natverse.org/nat/) package will use it when available to enable
large speed-ups for some primitive operations. End users should therefore not count on the  *natcpp* functions/interface being stable, especially during its
initial development.

Our plan is to get this onto
CRAN quite quickly so that binary packages will be available to Mac/Windows users who do not have a full C compiler toolchain installed. Once on CRAN, we will try to keep the package simple and stable, with infrequent
updates. A development version will be available via the relevant github
repository (see below). 
We do not intend to specify a GitHub remote so these updates will not be 
installed by default i.e. we will prefer CRAN binaries unless the user
intervenes.

## Installation

You will be able to install the released version of natcpp from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("natcpp")
```

For now please do:

``` r
install.packages("natmanager")
natmanager::install(pkgs = 'natverse/natcpp')
```


## Example

As already noted, the intention is that end users will only use the nat package
(and that will automagically use natcpp when available). However to prove that
everything is set up properly, you could do:

``` r
library(nat)
library(natcpp)
topntail(as.seglist(Cell07PNs[[1]]))

```

.
