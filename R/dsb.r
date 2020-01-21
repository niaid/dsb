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
#' @return a normalized matrix that can be added to any S4 object e.g. Seurat or SingleCellExperiment Object
#' @export
#'
#' @importFrom limma removeBatchEffect
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom stats prcomp
#' @importFrom stats sd
DSBNormalizeProtein = function(cell_protein_matrix, empty_drop_matrix, define.pseudocount = FALSE,
                               pseudocount.use, denoise.counts = TRUE, use.isotype.control = FALSE,
                               isotype.control.name.vec = NULL){
  # suppressMessages(require(mclust))
  # suppressMessages(require(limma))
  adt = cell_protein_matrix %>% as.matrix()
  adtu = empty_drop_matrix %>% as.matrix()

  if(define.pseudocount == TRUE) {
    adtu_log = log(adtu + pseudocount.use)
    adt_log = log(adt + pseudocount.use)
  } else {
    adtu_log = log(adtu + 10)
    adt_log = log(adt + 10)
  }
  # apply scaling by background (subtract control means from counts and divide by control sd)
  mu_u = apply(adtu_log, 1 , mean)
  sd_u = apply(adtu_log, 1 , sd)
  norm_adt = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u)
  # transpose
  if(denoise.counts == FALSE){
    return(norm_adt)
  } else {
    # per cell denoising: fit mixture of 2 gaussians to each cell's protein data, get background mean (mu 1)
    cellwise_background_mean = apply(norm_adt, 2, function(x) {
      g = mclust::Mclust(x, G=2, warn = F , verbose = F)
      return(g$parameters$mean[1])
    })
    gc()
    # using isotypes: (recommended)
    if (use.isotype.control == TRUE) {
      # define pc1 through isotypes and background protein as a cellwise correction factor
      noise_matrix = rbind(norm_adt[isotype.control.name.vec, ], cellwise_background_mean)
      get_noise_vector = function(noise_matrix) {
        g = prcomp(t(noise_matrix), scale = TRUE)
        return(g$x[ ,1])
      }
      noise_vector = get_noise_vector(noise_matrix)
      denoised_adt = limma::removeBatchEffect(norm_adt, covariates = noise_vector)
    } else {
      noise_vector = cellwise_background_mean
      denoised_adt = limma::removeBatchEffect(norm_adt, covariates = noise_vector)
    }
    return(denoised_adt)
  }
}
