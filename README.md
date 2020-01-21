
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dsb

<!-- badges: start -->

<!-- badges: end -->

The goal of dsb is: see biorxiv preprint <DOI:Here> we developed a
method specifically for normalizing CITEseq data,accounting for
protein-specific background noise caused by unbound antibody captured
and sequenced in droplets as well as correcting for the technical
component of variation in protein library size. We account for these
sources of experimental noisewith the empty droplets containing
background antibody reads and by defining a covariate of non-staining
proteins and isotype control antibodies respectively.

## Installation

You can install the released version of dsb from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dsb")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(dsb)
## basic example code
```
