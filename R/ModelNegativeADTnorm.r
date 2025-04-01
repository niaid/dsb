### Model negative adt norm
#' ModelNegativeADTnorm R function: Normalize single cell antibody derived tag (ADT) protein data.
#' This function defines the background level for each protein by fitting a 2 component Gaussian
#' mixture after log transformation. Empty Droplet ADT counts are not supplied. The fitted background
#' mean of each protein across all cells is subtracted from the log transformed counts. Note this is
#' distinct from and unrelated to the 2 component mixture used in the second step of
#' `DSBNormalizeProtein` which is fitted to all proteins of each cell. After this background correction
#' step, `ModelNegativeADTnorm` then models and removes technical cell to cell variations using the
#' same step II procedure as in the DSBNormalizeProtein function using identical function arguments.
#' This is a experimental function that performs well in testing and is motivated by our observation
#' in Supplementary Fig 1 in the dsb paper showing that the fitted background mean was concordant with
#' the mean of ambient ADTs in both empty droplets and unstained control cells. We recommend using
#' `ModelNegativeADTnorm` if empty droplets are not available.
#' See <https://www.nature.com/articles/s41467-022-29356-8> for details of the algorithm.
#' @author Matthew P. Mulè, \email{mattmule@@gmail.com}
#' @references https://www.nature.com/articles/s41467-022-29356-8
#' @param cell_protein_matrix Raw protein ADT UMI count data to be normalized. Cells - columns
#' Proteins (ADTs) - rows.
#' @param denoise.counts Recommended function default `denoise.counts = TRUE` and `use.isotype.control = TRUE`.
#' This runs step II of the dsb algorithm to define and remove cell to cell technical noise.
#' @param use.isotype.control Recommended function default `denoise.counts = TRUE` and `use.isotype.control = TRUE`.
#' This includes isotype controls in defining the dsb technical component.
#' @param isotype.control.name.vec A vector of the names of the isotype control proteins in the rows of the cells
#' and background matrix e.g. isotype.control.name.vec = c('isotype1', 'isotype2').
#' @param define.pseudocount FALSE (default) uses the value 10 optimized for protein ADT data.
#' @param pseudocount.use Must be defined if `define.pseudocount = TRUE`. This is the pseudocount to be added to
#' raw ADT UMI counts. Otherwise the default pseudocount used.
#' @param quantile.clipping FALSE (default), if outliers or a large range of values for some proteins are observed
#'  (e.g. -50 to 50) these are often from rare outlier cells. re-running the function with `quantile.clipping = TRUE`
#'  will adjust by applying 0.001 and 0.998th quantile value clipping to trim values to those max and min values.
#' @param quantile.clip if `quantile.clipping = TRUE`, a vector of the lowest and highest quantiles to clip. These can
#'  be tuned to the dataset size. The default c(0.001, 0.9995) optimized to clip only a few of the most extreme outliers.
#' @param fast.km Recommended to set this parameter to `TRUE` for large datasets. If `fast.km = TRUE`, the function defines
#' cell level background for step II with a a k=2 k-means cluster instead of a 2 component gaussian mixture. Increases speed
#' ~10x on 1 million cell benchmark with minimal impact on results.
#' @param return.stats if TRUE, returns a list, element 1 $dsb_normalized_matrix is the normalized adt matrix element 2
#'  $dsb_stats is the internal stats used by dsb during denoising (the background mean, isotype control values, and the
#'  final dsb technical component that is regressed out of the counts)
#'
#' @return Normalized ADT data are returned as a standard R "matrix" of cells (columns), proteins (rows) that can be
#' added to Seurat, SingleCellExperiment or python anndata object - see vignette. If return.stats = TRUE, function
#' returns a list: x$dsb_normalized_matrix normalized matrix, x$protein_stats are mean and sd of log transformed cell,
#' background and the dsb normalized values (as list). x$technical_stats includes the dsb technical component value for
#' each cell and each variable used to calculate the technical component.
#' @export
#'
#' @importFrom limma removeBatchEffect
#' @importFrom mclust Mclust mclustBIC
#' @importFrom stats prcomp sd quantile kmeans
#' @examples
#' library(dsb) # load example data cells_citeseq_mtx and empty_drop_matrix included in package
#'
#' # use a subset of cells and background droplets from example data
#' cells_citeseq_mtx = cells_citeseq_mtx[ ,1:400]
#' empty_drop_matrix = empty_drop_citeseq_mtx[ ,1:400]
#'
#' # example I
#' adt_norm = dsb::ModelNegativeADTnorm(
#'   # step I: remove ambient protein noise modeled by a gaussian mixture
#'   cell_protein_matrix = cells_citeseq_mtx,
#'
#'   # recommended step II: model and remove the technical component of each cell's protein data
#'   denoise.counts = TRUE,
#'   use.isotype.control = TRUE,
#'   isotype.control.name.vec = rownames(cells_citeseq_mtx)[67:70]
#' )
#'
ModelNegativeADTnorm = function(cell_protein_matrix,
                                denoise.counts = TRUE,
                                use.isotype.control = TRUE,
                                isotype.control.name.vec = NULL,
                                define.pseudocount = FALSE,
                                pseudocount.use,
                                fast.km = FALSE,
                                quantile.clipping = FALSE,
                                quantile.clip = c(0.001, 0.9995),
                                return.stats = FALSE){

  paste0('using empirical background estimate, ubiquitously expressed proteins',
         ' may be 0-centered;\n recommended usage of dsb: provide background ADT',
         ' matrix from\n  empty droplets as argument to empty_drop_matrix ')

  # formatting checks on step 2 functions
  a = isotype.control.name.vec
  cm = rownames(cell_protein_matrix)
  iso_detect = cm[grepl(pattern = 'sotype|Iso|iso|control|CTRL|ctrl|Ctrl|ontrol', x = cm)]

  # isotype.control.name.vec specified but some isotypes are not in input matrices
  if (!is.null(a) & !isTRUE(all(a %in% cm))){
    stop(paste0("some elements of isotype.control.name.vec are not in input data rownames: \n",
                setdiff(a,cm))
    )
  }
  # step II = FALSE - remind user isotypes unused. set isotype-related args to FALSE, NULL
  if (isFALSE(denoise.counts)) {
    print(paste0("Not running dsb step II (removal of cell to cell technical noise)",
                 " Setting use.isotype.control and isotype.control.name.vec to FALSE and NULL"))
    use.isotype.control = FALSE
    isotype.control.name.vec = NULL
    if (length(iso_detect) > 0) {
      print('potential isotype controls detected: ')
      print(iso_detect)
    }
  }
  # use.isotype.control = TRUE but isotype.control.name.vec = NULL:
  if (isTRUE(use.isotype.control) & is.null(isotype.control.name.vec)) {
    if (length(iso_detect) > 0) {
      print('potential isotype controls detected: ')
      print(iso_detect)
    }
    stop('if use.isotype.control = TRUE, set isotype.control.name.vec to rownames of isotype controls')
  }
  # denoise.counts = TRUE with use.isotype.control = FALSE:
  if (isTRUE(denoise.counts) & isFALSE(use.isotype.control)) {
    warning(paste0(
      '`use.isotype.control` = FALSE is not recommended if setting `denoise.counts` = TRUE',
      ' \nwhen isotype controls are available.\n If data include isotype controls',
      ' set `denoise.counts` = TRUE `use.isotype.control` = TRUE',
      ' \and set `isotype.control.name.vec` to a vector of isotype control rownames'
    ))
    if (length(iso_detect) > 0) {
      print('potential isotype controls detected: ')
      print(iso_detect)
    }
  }

  # STEP I - protein level background correction - create norm_adt
  # if matrices are dgTMatrix coerce into regular matrix
  adt = as.matrix(cell_protein_matrix)
  if(isTRUE(define.pseudocount)) {
    adt_log = log(adt + pseudocount.use)
  } else {
    adt_log = log(adt + 1)
  }
  # empirical estimate of protein background as fitted background mean across cells
  p.model = apply(adt_log, 1, function(x){
    mclust::Mclust(data = x, G = 2, verbose = FALSE, warn = FALSE)
  })
  mu1 = unlist(lapply(p.model, function(x) x$parameters$mean[[1]]))
  if (isFALSE(length(mu1) == length(rownames(adt_log)))) {
    ad_name = setdiff(rownames(adt_log), names(mu1))
    warning(
      paste0('empirical background cound not be fit for ',
             length(ad_name), ' proteins: ',
             paste(ad_name, collapse = ', '),
             '. Value returned will be log transformed without background correction.')
    )
    ad = as.numeric(rep(x = 0, length(ad_name)))
    names(ad) = ad_name
    mu1 = c(mu1 , ad)
    mu1 = mu1[match(rownames(adt_log) , names(mu1) )]
  }
  norm_adt = apply( adt_log, 2, function(x) (x - mu1) )


  # STEP II
  # dsb step II calculate dsb technical component and regress out of ambient corrected values
  if(isTRUE(denoise.counts)){
    print(paste0('fitting models to each cell for dsb technical component and',
                 ' removing cell to cell technical noise'))
    if(isTRUE(fast.km)){
      error_cols <- character(0)
      cellwise_background_mean <- sapply(seq_len(ncol(norm_adt)), function(i) {
        tryCatch({
          km <- stats::kmeans(norm_adt[, i], centers = 2, iter.max = 100, nstart = 5)
          min(km$centers)
        }, error = function(e) {
          error_cols <<- c(error_cols, colnames(norm_adt)[i])
          cat(sprintf("km.fast issue detecting background for cell %s: %s\n",
                      colnames(norm_adt)[i], e$message))
          NULL
        })
      }) # end of apply statement to define cell_background_mean
      # Error handling for km.fast
      if (length(error_cols) > 0) {
        cat("issues detecting background mean occurred in cells:",
            paste(error_cols, collapse = ", "), "\n")
      }
    } else {
      # original per cell background mean definition
      cellwise_background_mean = apply(norm_adt, 2, function(x) {
        g = mclust::Mclust(x, G=2, warn = FALSE, verbose = FALSE)
        return(g$parameters$mean[1])
      })
    }
    gc()
    if (isTRUE(use.isotype.control)) {
      noise_matrix = rbind(norm_adt[isotype.control.name.vec, ], cellwise_background_mean)
      get_noise_vector = function(noise_matrix) {
        g = stats::prcomp(t(noise_matrix), scale = TRUE)
        return(g$x[ ,1])
      }
      # regress out eigenvector of noise matrix
      noise_vector = get_noise_vector(noise_matrix)
      norm_adt = limma::removeBatchEffect(norm_adt, covariates = noise_vector)
    } else {
      # regress out mu1
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
  # write internal stats to R list to return if return.stats = TRUE
  if(isTRUE(return.stats)) {
    print('returning stats and results in a list; access normalized matrix with x$dsb_normalized_matrix')
    protein_stats = list(

      'raw cell matrix stats' = data.frame(
        cell_mean = apply(adt_log, 1 , mean),
        cell_sd = apply(adt_log, 1 , sd)),

      'dsb normalized matrix stats' =
        data.frame(
          dsb_mean = apply(norm_adt, 1 , mean),
          dsb_sd = apply(adt_log, 1 , sd))
    )
    if(isTRUE(denoise.counts)) {
      technical_stats = cbind(t(noise_matrix), dsb_technical_component = noise_vector)
    } else{
      technical_stats = NULL
    }
    ret_obj = list(
      'dsb_normalized_matrix' = norm_adt,
      'technical_stats' = technical_stats,
      'protein_stats' = protein_stats
    )
    return(ret_obj)
  } else {
    return(norm_adt)
  }
}
