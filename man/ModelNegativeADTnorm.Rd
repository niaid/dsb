% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelNegativeADTnorm.r
\name{ModelNegativeADTnorm}
\alias{ModelNegativeADTnorm}
\title{ModelNegativeADTnorm R function: Normalize single cell antibody derived tag (ADT) protein data.
This function defines the background level for each protein by fitting a 2 component Gaussian
mixture after log transformation. Empty Droplet ADT counts are not supplied. The fitted background
mean of each protein across all cells is subtracted from the log transformed counts. Note this is
distinct from and unrelated to the 2 component mixture used in the second step of
`DSBNormalizeProtein` which is fitted to all proteins of each cell. After this background correction
step, `ModelNegativeADTnorm` then models and removes technical cell to cell variations using the
same step II procedure as in the DSBNormalizeProtein function using identical function arguments.
This is a experimental function that performs well in testing and is motivated by our observation
in Supplementary Fig 1 in the dsb paper showing that the fitted background mean was concordant with
the mean of ambient ADTs in both empty droplets and unstained control cells. We recommend using
`ModelNegativeADTnorm` if empty droplets are not available.
See <https://www.nature.com/articles/s41467-022-29356-8> for details of the algorithm.}
\usage{
ModelNegativeADTnorm(
  cell_protein_matrix,
  denoise.counts = TRUE,
  use.isotype.control = TRUE,
  isotype.control.name.vec = NULL,
  define.pseudocount = FALSE,
  pseudocount.use,
  fast.km = FALSE,
  quantile.clipping = FALSE,
  quantile.clip = c(0.001, 0.9995),
  return.stats = FALSE
)
}
\arguments{
\item{cell_protein_matrix}{Raw protein ADT UMI count data to be normalized. Cells - columns
Proteins (ADTs) - rows.}

\item{denoise.counts}{Recommended function default `denoise.counts = TRUE` and `use.isotype.control = TRUE`.
This runs step II of the dsb algorithm to define and remove cell to cell technical noise.}

\item{use.isotype.control}{Recommended function default `denoise.counts = TRUE` and `use.isotype.control = TRUE`.
This includes isotype controls in defining the dsb technical component.}

\item{isotype.control.name.vec}{A vector of the names of the isotype control proteins in the rows of the cells
and background matrix e.g. isotype.control.name.vec = c('isotype1', 'isotype2').}

\item{define.pseudocount}{FALSE (default) uses the value 10 optimized for protein ADT data.}

\item{pseudocount.use}{Must be defined if `define.pseudocount = TRUE`. This is the pseudocount to be added to
raw ADT UMI counts. Otherwise the default pseudocount used.}

\item{fast.km}{Recommended to set this parameter to `TRUE` for large datasets. If `fast.km = TRUE`, the function defines
cell level background for step II with a a k=2 k-means cluster instead of a 2 component gaussian mixture. Increases speed
~10x on 1 million cell benchmark with minimal impact on results.}

\item{quantile.clipping}{FALSE (default), if outliers or a large range of values for some proteins are observed
(e.g. -50 to 50) these are often from rare outlier cells. re-running the function with `quantile.clipping = TRUE`
will adjust by applying 0.001 and 0.998th quantile value clipping to trim values to those max and min values.}

\item{quantile.clip}{if `quantile.clipping = TRUE`, a vector of the lowest and highest quantiles to clip. These can
be tuned to the dataset size. The default c(0.001, 0.9995) optimized to clip only a few of the most extreme outliers.}

\item{return.stats}{if TRUE, returns a list, element 1 $dsb_normalized_matrix is the normalized adt matrix element 2
$dsb_stats is the internal stats used by dsb during denoising (the background mean, isotype control values, and the
final dsb technical component that is regressed out of the counts)}
}
\value{
Normalized ADT data are returned as a standard R "matrix" of cells (columns), proteins (rows) that can be
added to Seurat, SingleCellExperiment or python anndata object - see vignette. If return.stats = TRUE, function
returns a list: x$dsb_normalized_matrix normalized matrix, x$protein_stats are mean and sd of log transformed cell,
background and the dsb normalized values (as list). x$technical_stats includes the dsb technical component value for
each cell and each variable used to calculate the technical component.
}
\description{
ModelNegativeADTnorm R function: Normalize single cell antibody derived tag (ADT) protein data.
This function defines the background level for each protein by fitting a 2 component Gaussian
mixture after log transformation. Empty Droplet ADT counts are not supplied. The fitted background
mean of each protein across all cells is subtracted from the log transformed counts. Note this is
distinct from and unrelated to the 2 component mixture used in the second step of
`DSBNormalizeProtein` which is fitted to all proteins of each cell. After this background correction
step, `ModelNegativeADTnorm` then models and removes technical cell to cell variations using the
same step II procedure as in the DSBNormalizeProtein function using identical function arguments.
This is a experimental function that performs well in testing and is motivated by our observation
in Supplementary Fig 1 in the dsb paper showing that the fitted background mean was concordant with
the mean of ambient ADTs in both empty droplets and unstained control cells. We recommend using
`ModelNegativeADTnorm` if empty droplets are not available.
See <https://www.nature.com/articles/s41467-022-29356-8> for details of the algorithm.
}
\examples{
library(dsb) # load example data cells_citeseq_mtx and empty_drop_matrix included in package

# use a subset of cells and background droplets from example data
cells_citeseq_mtx = cells_citeseq_mtx[ ,1:400]
empty_drop_matrix = empty_drop_citeseq_mtx[ ,1:400]

# example I
adt_norm = dsb::ModelNegativeADTnorm(
  # step I: remove ambient protein noise modeled by a gaussian mixture
  cell_protein_matrix = cells_citeseq_mtx,

  # recommended step II: model and remove the technical component of each cell's protein data
  denoise.counts = TRUE,
  use.isotype.control = TRUE,
  isotype.control.name.vec = rownames(cells_citeseq_mtx)[67:70]
)

}
\references{
https://www.nature.com/articles/s41467-022-29356-8
}
\author{
Matthew P. Mulè, \email{mattmule@gmail.com}
}
