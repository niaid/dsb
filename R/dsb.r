#' the DSB normalization function
#'
#' @param cell_protein_matrix the raw protein count data to be normalized cells = columns, proteins = rows
#' @param empty_drop_matrix the raw empty droplet protein count data used for background correction
#' @param denoise.counts set to TRUE by default (recommended) - defines and regresses each cell's technical component using background proteins in each cell and itotype controls when
#' @param use.isotype.control set to TRUE by default (recommended) with denoise.counts = TRUE include isotype controls in defining the ?
#' @param isotype.control.name.vec a vector of the names of the isotype control proteins exactly as in the data to be normalized
#' @param define.pseudocount TRUE FALSE : false by default users can supply a pseudocount besides the default 10 which is optimal for CITEseq data.
#' @param pseudocount.use the pseudocount to use if overriding the default pseudocount by setting efine.pseudocount to TRUE
#'
#' @return a normalized R "matrix" of cells by proteins that can be added to any Seurat, SingleCellExperiment or python anndata object - see vignette
#' @export
#'
#' @importFrom limma removeBatchEffect
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom stats prcomp
#' @importFrom stats sd
#' @examples
#' library(dsb) # lazy load example data cells_citeseq_mtx and empty_drop_matrix included in package
#'
#' # use a subset of 40 cells and bacground droplets from example data
#' cells_citeseq_mtx = cells_citeseq_mtx[ ,1:400]
#' empty_drop_matrix = empty_drop_citeseq_mtx[ ,1:400]
#'
#' adt_norm = dsb::DSBNormalizeProtein(
#'
#'   # step I: remove ambient protien noise reflected in counts from empty droplets
#'   cell_protein_matrix = cells_citeseq_mtx,
#'   empty_drop_matrix = empty_drop_matrix,
#'
#'   # recommended step II: model and remove the technical component of each cell's protein library
#'   denoise.counts = TRUE,
#'   use.isotype.control = TRUE,
#'   isotype.control.name.vec = rownames(cells_citeseq_mtx)[67:70]
#' )
#'
#' # example II - experiments without isotype controls
#'
#' adt_norm = dsb::DSBNormalizeProtein(
#'   cell_protein_matrix = cells_citeseq_mtx,
#'   empty_drop_matrix = empty_drop_matrix,
#'   denoise.counts = FALSE
#' )
#'
DSBNormalizeProtein = function(cell_protein_matrix, empty_drop_matrix,
                               denoise.counts = TRUE,  use.isotype.control = TRUE, isotype.control.name.vec = NULL,
                               define.pseudocount = FALSE,pseudocount.use){

  adt = cell_protein_matrix %>% as.matrix()
  adtu = empty_drop_matrix %>% as.matrix()

  if(define.pseudocount == TRUE) {
    adtu_log = log(adtu + pseudocount.use)
    adt_log = log(adt + pseudocount.use)
  } else {
    adtu_log = log(adtu + 10)
    adt_log = log(adt + 10)
  }
  # ambient correction background rescaling
  mu_u = apply(adtu_log, 1 , mean)
  sd_u = apply(adtu_log, 1 , sd)
  norm_adt = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u)
  if(denoise.counts == FALSE){
    return(norm_adt)
  } else {
    # 2 component mixture
    cellwise_background_mean = apply(norm_adt, 2, function(x) {
      g = mclust::Mclust(x, G=2, warn = F , verbose = F)
      return(g$parameters$mean[1])
    })
    gc()
    # using isotypes
    if (use.isotype.control == TRUE) {
      # define technical component for each cell as primary latent component
      noise_matrix = rbind(norm_adt[isotype.control.name.vec, ], cellwise_background_mean)
      get_noise_vector = function(noise_matrix) {
        g = stats::prcomp(t(noise_matrix), scale = TRUE)
        return(g$x[ ,1])
      }
      noise_vector = get_noise_vector(noise_matrix)
      denoised_adt = limma::removeBatchEffect(norm_adt, covariates = noise_vector)
    } else {
      # without isotype controls
      noise_vector = cellwise_background_mean
      denoised_adt = limma::removeBatchEffect(norm_adt, covariates = noise_vector)
    }
    return(denoised_adt)
  }
}




