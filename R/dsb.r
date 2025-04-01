#' DSBNormalizeProtein R function: Normalize single cell antibody derived tag (ADT) protein data.
#' This function corrects for both protein specific and cell to cell technical noise in antibody derived tag (ADT) data.
#' For datasets without access to empty drops use dsb::ModelNegativeADTnorm.
#' See <https://www.nature.com/articles/s41467-022-29356-8> for details of the algorithm.
#' @author Matthew P. Mulè, \email{mattmule@@gmail.com}
#' @references https://www.nature.com/articles/s41467-022-29356-8
#' @param cell_protein_matrix Raw protein ADT UMI count data to be normalized. Cells = columns, protein antibody = rows.
#' @param empty_drop_matrix Raw empty droplet / background ADT UMI count data used for background correction.
#' Cells - columns and Proteins (ADTs) - rows. See vignettes for how to define background matrix from the
#' raw_feature_bc_matrix output from Cell Ranger or from other alignment tools such as kallisto and Cite-Seq-Count.
#' For datasets without access to empty drops use dsb::ModelNegativeADTnorm.
#' @param scale.factor one of `standardize` or `mean.subtract`.
#' Scale factor specifies how to implement protein level denoising. The recommended default is `standardize` which is the method
#' described in Mulè et al 2022 as "Step I". For each protein, this subtracts the mean and divides by the standard deviation
#' of that protein observed in the empty droplets, making the resulting value interpretable as the number of standard deviations
#' above the average of the background for that protein observed in empty droplets.
#' `mean.subtract`, subtracts the mean without dividing by the standard deviation; can be used in scenarios where low
#' background levels are detected systematically for most proteins in the dataset background standard deviation may be unstable.
#' @param denoise.counts Recommended function default `denoise.counts = TRUE` and `use.isotype.control = TRUE`.
#' This removes remove cell to cell technical noise as described as "step II" in Mulè et al 2022.
#' @param use.isotype.control Recommended function default `denoise.counts = TRUE` and `use.isotype.control = TRUE`.
#' This includes isotype controls in defining the dsb technical component.
#' @param isotype.control.name.vec A vector of the names of the isotype control proteins in the rows of the cells
#' and background matrix e.g. `isotype.control.name.vec = c('isotype1', 'isotype2')`.
#' @param define.pseudocount `FALSE` (default) uses the value 10 optimized for protein ADT data.
#' Any pseudocount can be used by setting this argument to `FALSE` and specifying `pseudocount.use`.
#' @param pseudocount.use Must be defined if `define.pseudocount = TRUE`. This is the pseudocount to be added to
#' raw ADT UMI counts. Otherwise the default pseudocount used.
#' @param quantile.clipping FALSE (default), if outliers or a large range of values for some proteins are observed
#'  (e.g. -50 to 50) these are often from rare outlier cells. re-running the function with `quantile.clipping = TRUE`
#'  will adjust by applying 0.001 and 0.998th quantile value clipping to trim values to those max and min values. If
#'  range of normalized values are still very broad and high (e.g. above 40) try setting `scale.factor = mean.subtract`.
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
#' dsb_output = dsb::DSBNormalizeProtein(
#'    cell_protein_matrix = cells_citeseq_mtx,
#'    empty_drop_matrix = empty_drop_matrix,
#'    isotype.control.name.vec = rownames(cells_citeseq_mtx)[67:70],
#'    return.stats = TRUE
#' )
#'
#' # the dsb normalized matrix to be used in downstream analysis is dsb_output$dsb_normalized_matrix
#' # protein level stats are in dsb_output$protein_stats
#' # cell-level stats are in dsb_output$technical_stats
#'
DSBNormalizeProtein = function(cell_protein_matrix,
                               empty_drop_matrix,
                               denoise.counts = TRUE,
                               use.isotype.control = TRUE,
                               isotype.control.name.vec = NULL,
                               define.pseudocount = FALSE,
                               pseudocount.use,
                               quantile.clipping = FALSE,
                               quantile.clip = c(0.001, 0.9995),
                               fast.km = FALSE,
                               scale.factor = c('standardize', 'mean.subtract')[1],
                               return.stats = FALSE){
  # check background matrix provided
  if (is.null(empty_drop_matrix)) {
    stop("'empty_drop_matrix' was not specified. To normalize without empty drops, use ModelNegativeADTnorm function.")
  }

  # formatting checks on input matrices
  a = isotype.control.name.vec
  b = rownames(empty_drop_matrix)
  cm = rownames(cell_protein_matrix)

  if (!isTRUE(all.equal(cm, b))){
    diff = c(setdiff(cm, b), setdiff(b, cm))

    # stop if nrow of cells and background is not equal
    if(!isTRUE(all.equal( nrow(cell_protein_matrix), nrow(empty_drop_matrix)))){
      stop(paste0(
        'the number of rows (features) in `cell_protein_matrix` and `empty_drop_matrix` do not match check: ',
        diff))
    }

    # stop if the rownames are not the same
    if (length(diff) > 0) {
      stop(paste0('rows of cell and background matrices have mis-matching names: ', diff))
    }

    # reorder rows of empty drop matrix to match the cells if in a different order
    if (!identical(rownames(cell_protein_matrix), rownames(empty_drop_matrix)) &&
        setequal(rownames(cell_protein_matrix), rownames(empty_drop_matrix))) {
      rmatch = match(rownames(cell_protein_matrix), rownames(empty_drop_matrix))
      empty_drop_matrix = empty_drop_matrix[rmatch, , drop = FALSE]
      warning("Reordered 'empty_drop_matrix' rows to match 'cell_protein_matrix'")
    }
  }

  # formatting checks on step 2 functions
  iso_detect = cm[grepl(pattern = 'sotype|Iso|iso|control|CTRL|ctrl|Ctrl|ontrol', x = cm)]

  # isotype.control.name.vec specified but some isotypes are not in input matrices
  if (!is.null(a) & !isTRUE(all(a %in% b)) & !isTRUE(all(a %in% cm))){
    stop(paste0("some elements of isotype.control.name.vec are not in input data rownames: \n",
                'cell_protein_matrix - ', setdiff(a,cm),
                ' \nempty_drop_matrix - ', setdiff(a,b))
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
  } # end step 2 argument checks

  # STEP I - protein level background correction
  # coerce to regular matrix if stored as a sparse matrix
  adt = as.matrix(cell_protein_matrix)
  adtu = as.matrix(empty_drop_matrix)
  if(isTRUE(define.pseudocount)) {
    adtu_log = log(adtu + pseudocount.use)
    adt_log = log(adt + pseudocount.use)
  } else {
    adtu_log = log(adtu + 10)
    adt_log = log(adt + 10)
  }
  print("correcting ambient protein background noise")
  mu_u = rowMeans(adtu_log)
  sd_u = apply(adtu_log, 1, stats::sd)
  if (scale.factor == 'standardize') {
    # print low sd proteins
    if (any(sd_u < 0.05)) {
      print(paste0('proteins below have low background log norm sd < 0.05',
                   ' check raw and normalized distributions. ',
                   ' protein stats can be returned with return.stats = TRUE'
      ))
      print(names(which(sd_u<0.05)))
    }
    # standardize
    norm_adt = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u)
  }
  if(scale.factor == 'mean.subtract'){
    norm_adt = apply(adt_log, 2, function(x) (x  - mu_u))
  }
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
  #### End of step 2

  # apply quantile clipping of outliers
  if (isTRUE(quantile.clipping)) {
    ql = apply(norm_adt, 1, FUN = stats::quantile, quantile.clip[1])
    qh = apply(norm_adt, 1, FUN = stats::quantile, quantile.clip[2])
    for (i in 1:nrow(norm_adt)) {
      norm_adt[i, ] = ifelse(norm_adt[i, ] < ql[i], ql[i], norm_adt[i, ])
      norm_adt[i, ] = ifelse(norm_adt[i, ] > qh[i], qh[i], norm_adt[i, ])
    }
  }
  if(isTRUE(return.stats)) {
    print('returning results in a list; normalized matrix accessed x$dsb_normalized_matrix')
    protein_stats = list(

      'raw cell matrix stats' = data.frame(
        cell_mean = apply(adt_log, 1 , mean),
        cell_sd = apply(adt_log, 1 , sd)
      ),

      'dsb normalized matrix stats' =
        data.frame(
          dsb_mean = apply(norm_adt, 1 , mean),
          dsb_sd = apply(adt_log, 1 , sd)
        )
    )
    if(isTRUE(denoise.counts)) {
      if(isTRUE(use.isotype.control)){
        technical_stats = cbind(t(noise_matrix), dsb_technical_component = noise_vector)
      } else { #when use.isotype.control is FALSE, the noise_matrix does not exist
        technical_stats = data.frame(dsb_technical_component = noise_vector)
      }
    } else{
      technical_stats = NULL
    }
    if(is.null(empty_drop_matrix)){
      background.stats = NULL
    }else{
      background.stats = data.frame(
        background_mean = mu_u,
        background_sd = sd_u
      )
      protein_stats = c(protein_stats, background.stats)
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
