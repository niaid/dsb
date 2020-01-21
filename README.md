
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dsb

<!-- badges: start -->

<!-- badges: end -->

This package was developed at [John Tsang’s
Lab](https://www.niaid.nih.gov/research/john-tsang-phd) by Matt Mulè and
Andrew Martins. The package implements our normalization and denoising
method for CITEseq data. Technical discussion of how the method works
can be found in [the biorxiv preprint](https://biorxiv.org) We utilized
the dsb package to normalize CITEseq data reported in this paper
[](https://)

In [the biorxiv preprint](https://biorxiv.org), comparing unstained
control cells and empty droplets we found the major contribotor to
background noise in CITEseq data is unbound antibody captured and
sequenced in droplets. DSB corrects for this background by leveraging
empty droplets which serve as a “built in” noise measurement in any
droplet capture single cell platform (e.g. 10X, dropseq, indrop).

In addition we define a per-cell denoising covariate to account for the
technical component of library size differences between cells which
removes spurious cluster formation derived from globally dimcells
clustering together.

## installation

You can install the released version of dsb in your R session with the
command
below

``` r
# this is analagous to install.packages("package), you need the package devtools to install a package from a github repository like this one. 
require(devtools)
#> Loading required package: devtools
#devtools::install_github(repo = 'MattPM/dsb')
```

## quickstart

``` r
# load the example data
library(dsb)
#load(cells_citeseq_mtx)
#load(empty_drop_citeseq_mtx)

# normalize
normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx,
                                        empty_drop_matrix = empty_drop_citeseq_mtx)
```

## The full version (recommended)

By default dsb defines the per cell denoising covariate by fitting a
gaussian mixture model to the log + 10 counts of each cell and defining
the noise ocvariates the mean. We reccomend including the counts from
isotype controls in each cell in the denoising covariates.

``` r

# define a vector of the isotype controls in the data 
isotypes = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT","MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT")

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx,
                                        empty_drop_matrix = empty_drop_citeseq_mtx,
                                        use.isotype.control = TRUE,
                                        isotype.control.name.vec = isotypes)
```

# visualize distributions of CD4 and CD8

plot raw points (overplotted) and points with labeled density
distributions (similar to
flow)

``` r
# plot this and avoid plotting by adding a density gradient like a flowjo plot
# this nice density function is from here: https://slowkow.com/notes/ggplot2-color-by-density/
get_density = function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

library(ggplot2)
data.plot = normalized_matrix %>% t %>%
  as.data.frame() %>% 
  dplyr::select(CD4_PROT, CD8_PROT) 
data.plot = data.plot %>%   dplyr::mutate(density = get_density(data.plot$CD4_PROT, data.plot$CD8_PROT, n = 100)) 

# plot with and without density gradient
p1 = ggplot(data.plot, aes(x = CD8_PROT, y = CD4_PROT, color = density)) +
  geom_point(size = 0.4) +
  geom_vline(xintercept = 0, color = "red", linetype = 2) + 
  geom_hline(yintercept = 0, color = "red", linetype = 2) + 
  viridis::scale_color_viridis(option = "B") +  
  scale_shape_identity() 
p2 = ggplot(data.plot, aes(x = CD8_PROT, y = CD4_PROT)) +
  geom_point(size = 0.4) +
  geom_vline(xintercept = 0, color = "red", linetype = 2) + 
  geom_hline(yintercept = 0, color = "red", linetype = 2) 

cowplot::plot_grid(p1,p2)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

The plots above show the actual protein distributions. There is no
artificial jitter added to points.

# How do I get the empty droplets?

See Vignettes for more info.
