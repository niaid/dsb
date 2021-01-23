#' the DSB normalization function
#'
#' @param cell_protein_matrix the data to be normalized cells = columns, proteins = rows
#' @param empty_drop_matrix the empty droplets used for background correction
#' @param define.pseudocount TRUE FALSE : false by default users can supply a pseudocount besides the default 10 which is optimal for CITEseq data.
#' @param pseudocount.use the pseudocount to use if overriding the default pseudocount by setting efine.pseudocount to TRUE
#' @param denoise.counts account for a per cell background correction factor by getting the mean of negative proteins in that cell + or -  isotype controls
#' @param use.isotype.control should the denoising covariate include isotype controls - recommended but false by default
#' @param isotype.control.name.vec a vector of the names of the isotype control proteins exactly as in the data to be normalized
#'
#' @return a normalized matrix of cells by proteins that can be added to any Seurat,  SingleCellExperiment or python anndata object - see vignette
#' @export
#'
#' @importFrom limma removeBatchEffect
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom stats prcomp
#' @importFrom stats sd
#' @examples
#' library(dsb) # lazy loads example data cells_citeseq_mtx and empty_drop_matrix included in package
#'
#' adt_norm = DSBNormalizeProtein(
#'   # remove ambient protien noise reflected in counts from empty droplets
#'   cell_protein_matrix = cells_citeseq_mtx, # cell-containing droplet raw protein count matrix
#'   empty_drop_matrix = empty_drop_citeseq_mtx, # empty/background droplet raw protein counts
#'
#'   # recommended step II: model and remove the technical component of each cell's protein library
#'   denoise.counts = TRUE, # model and remove each cell's technical component
#'   use.isotype.control = TRUE, # use isotype controls to define the technical component
#'   isotypecontrol.name.vec = rownames(cells_citeseq_mtx)[67:70] # vector of isotype control names
#' )
DSBNormalizeProtein = function(cell_protein_matrix, empty_drop_matrix, define.pseudocount = FALSE,
                               pseudocount.use, denoise.counts = TRUE, use.isotype.control = TRUE,
                               isotype.control.name.vec = NULL){

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
        g = prcomp(t(noise_matrix), scale = TRUE)
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




