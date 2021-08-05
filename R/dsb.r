#' the DSB normalization function
#'
#' @param cell_protein_matrix the raw protein count data to be normalized cells = columns, proteins = rows
#' @param empty_drop_matrix the raw empty droplet protein count data used for background correction
#' @param denoise.counts set to TRUE by default (recommended) - defines and regresses each cell's technical component using background proteins in each cell and isotype controls when
#' @param use.isotype.control set to TRUE by default (recommended) with denoise.counts = TRUE include isotype controls in defining the ?
#' @param isotype.control.name.vec to be used if denoise.counts = TRUE. a vector of the names of the isotype control proteins in the rows of the cells and background matrix e.g. c('isotype1', 'isotype2') or rownames(cells_citeseq_mtx)[grepl('sotype', rownames(cells_citeseq_mtx))]
#' @param define.pseudocount TRUE FALSE : false by default users can supply a pseudocount besides the default 10 which is optimal for CITEseq data.
#' @param pseudocount.use the pseudocount to use if overriding the default pseudocount by setting define.pseudocount to TRUE
#' @param quantile.clipping default to true, apply 0.001 and 0.998th quantile value clipping to handle low and high magnitude outliers - useful for downstream modeling and visualization
#' @param quantile.clip if using quantile clipping a vector of the lowest and highest quantile to clip, the default c(0.001, 0.9995) optimized to clip only a few of the most extreme outliers.
#' @param return.stats if TRUE, returns a list, element 1 is the normalized adt matrix element 2 is the internal stats used by dsb during denoising (the background mean, isotype control values, and the final dsb technical component that is regressed out of the counts)
#'
#' @return a normalized R "matrix" of cells by proteins that can be added to any Seurat, SingleCellExperiment or python anndata object - see vignette
#' @export
#'
#' @importFrom limma removeBatchEffect
#' @importFrom mclust Mclust mclustBIC
#' @importFrom stats prcomp sd quantile
#' @examples
#' library(dsb) # lazy load example data cells_citeseq_mtx and empty_drop_matrix included in package
#'
#' # use a subset of cells and background droplets from example data
#' cells_citeseq_mtx = cells_citeseq_mtx[ ,1:400]
#' empty_drop_matrix = empty_drop_citeseq_mtx[ ,1:400]
#'
#' # example I
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
#' adt_norm = dsb::DSBNormalizeProtein(
#'   cell_protein_matrix = cells_citeseq_mtx,
#'   empty_drop_matrix = empty_drop_matrix,
#'   denoise.counts = FALSE
#' )
#'
#'# example III - return dsb internal stats used during denoising for each cell
#'# returns a 2 element list - the normalized matrix and the internal stats
#' dsb_object = dsb::DSBNormalizeProtein(
#'    cell_protein_matrix = cells_citeseq_mtx,
#'    empty_drop_matrix = empty_drop_matrix,
#'    isotype.control.name.vec = rownames(cells_citeseq_mtx)[67:70],
#'    return.stats = TRUE
#' )
#'
#' # the dsb normalized matrix to be used in downstream analysis
#' dsb_object$dsb_normalized_matrix
#'
#' # the internal dsb stats; can be examined e.g. for outliers
#' dsb_object$dsb_stats
#'
#'
DSBNormalizeProtein = function(cell_protein_matrix, empty_drop_matrix, denoise.counts = TRUE,
                               use.isotype.control = TRUE, isotype.control.name.vec = NULL,
                               define.pseudocount = FALSE, pseudocount.use,quantile.clipping = FALSE,
                               quantile.clip = c(0.001, 0.9995), return.stats = FALSE){
  # error handling
  if(!is.null(isotype.control.name.vec)) {
    if(!all(isotype.control.name.vec %in% rownames(cell_protein_matrix))) {
      stop("the following isotype control cannot be found in the rownames of `cell_protein_matrix`: ",
           setdiff(isotype.control.name.vec, rownames(cell_protein_matrix)))
    }
  }
  if (isTRUE(return.stats) & isFALSE(denoise.counts)) {
    stop("set return.stats = FALSE to run dsb if denoise.counts = FALSE")
  }
  if (isTRUE(use.isotype.control) & is.null(isotype.control.name.vec)) {
    stop('must specify a vector of isotype control names if use.isotype.control is set to TRUE (recommended)')
  }
  # convert input data to matrix (e.g. for sparse matrix conversion) and log transform
  adt = cell_protein_matrix %>% as.matrix()
  adtu = empty_drop_matrix %>% as.matrix()
  # log transform
  if(isTRUE(define.pseudocount)) {
    adtu_log = log(adtu + pseudocount.use)
    adt_log = log(adt + pseudocount.use)
  } else {
    adtu_log = log(adtu + 10)
    adt_log = log(adt + 10)
  }
  # rescale cells based on expected noise in background drops
  print("correcting ambient protein background noise")
  mu_u = apply(adtu_log, 1 , mean)
  sd_u = apply(adtu_log, 1 , sd)
  norm_adt = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u)
  # step II calculate dsb technical component and regress out of ambient corrected values
  if(isTRUE(denoise.counts)){
    print(paste0('calculating dsb technical component for each cell to remove cell to cell techncial noise'))
    cellwise_background_mean = apply(norm_adt, 2, function(x) {
      g = mclust::Mclust(x, G=2, warn = FALSE, verbose = FALSE)
      return(g$parameters$mean[1])
      })
    gc()
    if (isTRUE(use.isotype.control)) {
      noise_matrix = rbind(norm_adt[isotype.control.name.vec, ], cellwise_background_mean)
      get_noise_vector = function(noise_matrix) {
        g = stats::prcomp(t(noise_matrix), scale = TRUE)
        return(g$x[ ,1])
        }
      noise_vector = get_noise_vector(noise_matrix)
      norm_adt = limma::removeBatchEffect(norm_adt, covariates = noise_vector)
      } else {
        print(paste0('denoising without using isotype controls is not recommended. ',
                     'If the raw ADT data include isotype controls, then recommended usage of dsb is to ',
                     "set `denoise.counts` = TRUE, `use.isotype.control` = TRUE, and define `isotype.control.name.vec`"))
        noise_vector = cellwise_background_mean
        norm_adt = limma::removeBatchEffect(norm_adt, covariates = noise_vector)
      }
    }
  # apply quantile clipping of outliers
  if (isTRUE(quantile.clipping)) {
    ql = apply(norm_adt, 1, FUN = stats::quantile, quantile.clip[1])
    qh = apply(norm_adt, 1, FUN = stats::quantile, quantile.clip[2])
    for (i in 1:nrow(norm_adt)) {
      norm_adt[i, ] = ifelse(norm_adt[i, ] < ql[i], ql[i], norm_adt[i, ])
      norm_adt[i, ] = ifelse(norm_adt[i, ] > qh[i], qh[i], norm_adt[i, ])
      }
    }
  # return dsb internal stats
  if (isTRUE(return.stats)) {
    print('returning a list object; access normalized matrix with x$dsb_normalized_matrix & stats with x$dsb_stats')
    dsb_stats = cbind(t(noise_matrix), dsb_technical_component = noise_vector)
    ret_obj = list('dsb_normalized_matrix' = norm_adt, 'dsb_stats' = dsb_stats)
    return(ret_obj)
    } else {
      return(norm_adt)
      }
}
