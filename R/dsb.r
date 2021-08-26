#' Normalize single cell antibody derived tag (ADT) protein data with the DSBNormalizeProtein function. This single function runs step I (ambient protein background correction) and step II (defining and removing cell to cell technical variation) of the dsb normalization method. See <https://www.biorxiv.org/content/10.1101/2020.02.24.963603v3> for details of the algorithm.
#' @author Matthew P. Mulè, \email{mattmule@@gmail.com}
#' @param cell_protein_matrix Raw protein ADT count data to be normalized with cells as columns and proteins as rows. See vignette, this is defined after quality control outlier cell removal based on the filtered output from Cell Ranger. Any CITE-seq count alignment tool can be used to define this as well.
#' @param empty_drop_matrix Raw empty droplet protein count data used for background correction with cells as columns and proteins as rows. This can easily be defined from the raw output from Cell Ranger (see vignette). Any count alignment tool for CITE-seq can be used to align and define these background drops.
#' @param denoise.counts TRUE (default) recommended to keep this TRUE and use with use.isotype.control = TRUE. This runs step II of the dsb algorithm to define and remove cell to cell technical noise.
#' @param use.isotype.control TRUE (default) recommended to use this with denoise.counts = TRUE. This includes isotype controls in defining the dsb technical component.
#' @param isotype.control.name.vec A vector of the names of the isotype control proteins in the rows of the cells and background matrix e.g. isotype.control.name.vec = c('isotype1', 'isotype2') or rownames(cells_citeseq_mtx)[grepl('sotype', rownames(cells_citeseq_mtx))]
#' @param define.pseudocount FALSE (default) uses the value 10 optimized for protein ADT data.
#' @param pseudocount.use the pseudocount to use if overriding the default pseudocount by setting define.pseudocount = TRUE
#' @param quantile.clipping FALSE (default), if outliers or a large range of values for some proteins is seen (e.g. -50 to 50) re-run with quantile.clipping = TRUE. This applies 0.001 and 0.998th quantile value clipping to handle low and high magnitude outliers.
#' @param quantile.clip if quantile.clipping = TRUE, one can provide a vector of the lowest and highest quantile to clip, these can be tuned to the dataset size. The default c(0.001, 0.9995) optimized to clip only a few of the most extreme outliers.
#' @param return.stats if TRUE, returns a list, element 1 $dsb_normalized_matrix is the normalized adt matrix element 2 $dsb_stats is the internal stats used by dsb during denoising (the background mean, isotype control values, and the final dsb technical component that is regressed out of the counts)
#'
#' @return The normalized ADT data are returned as a standard R "matrix" of cells (columns) by proteins (rows) that can be added to any Seurat, SingleCellExperiment or python anndata object - see vignette.
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
#'   # step I: remove ambient protein noise reflected in counts from empty droplets
#'   cell_protein_matrix = cells_citeseq_mtx,
#'   empty_drop_matrix = empty_drop_matrix,
#'
#'   # recommended step II: model and remove the technical component of each cell's protein data
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
#' # the internal dsb stats; can be examined for outliers see vignette FAQ
#' dsb_object$dsb_stats
#'
#'
DSBNormalizeProtein = function(cell_protein_matrix, empty_drop_matrix, denoise.counts = TRUE,
                               use.isotype.control = TRUE, isotype.control.name.vec = NULL,
                               define.pseudocount = FALSE, pseudocount.use, quantile.clipping = FALSE,
                               quantile.clip = c(0.001, 0.9995), return.stats = FALSE){
  # error handling
  if (!is.null(isotype.control.name.vec)) {
    if(!all(isotype.control.name.vec %in% rownames(cell_protein_matrix))) {
      stop("the following isotype.control.name.vec are not in the rownames of `cell_protein_matrix`: ",
           setdiff(isotype.control.name.vec, rownames(cell_protein_matrix)))
      }
  }
  if (isFALSE(denoise.counts)) {
    print(paste0("running step I ambient correction only, not removing cell to cell technical noise.",
                 " Setting use.isotype.control to FALSE"))
    use.isotype.control = FALSE
    isotype.control.name.vec = NULL
  }
  if (isTRUE(denoise.counts) & isFALSE(use.isotype.control)) {
    warning(
      'denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available.',
      ' If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE',
      ' and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix'
      )
    iso_detect = rownames(cell_protein_matrix)[grepl('sotype|Iso|iso', rownames(cell_protein_matrix))]
    if (length(iso_detect) > 0) {
      print('potential isotype controls detected: ')
      print(iso_detect)
    }
  }
  if(isTRUE(return.stats) & isFALSE(denoise.counts)) {
    stop("set return.stats = FALSE if denoise.counts = FALSE")
    }
  if (isTRUE(use.isotype.control) & is.null(isotype.control.name.vec)) {
    stop('must specify a vector of isotype control names if use.isotype.control is set to TRUE')
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
    dsb_stats = cbind(t(noise_matrix), dsb_technical_component = abs(noise_vector))
    ret_obj = list('dsb_normalized_matrix' = norm_adt, 'dsb_stats' = dsb_stats)
    return(ret_obj)
    } else {
      return(norm_adt)
    }
}
