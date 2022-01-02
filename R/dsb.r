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
#' @return Normalized ADT data are returned as a standard R "matrix" of cells (columns), proteins (rows) that can be added to Seurat, SingleCellExperiment or python anndata object - see vignette. If return.stats = TRUE, function returns a list: x$dsb_normalized_matrix normalized matrix, x$protein_stats are mean and sd of log transformed cell, background and the dsb normalized values (as list). x$technical_stats includes the dsb technical component value for each cell and each variable used to calculate the technical component.
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
DSBNormalizeProtein = function(cell_protein_matrix, empty_drop_matrix, denoise.counts = TRUE,
                               use.isotype.control = TRUE, isotype.control.name.vec = NULL,
                               define.pseudocount = FALSE, pseudocount.use, quantile.clipping = FALSE,
                               quantile.clip = c(0.001, 0.9995), return.stats = FALSE){
  # input formatting checks:

  # check format of input cell and background matrices
  # rows (features / proteins) must be of the same length and in the same order
  # detect unequal rows and parse through possible reasons with specific error messages:
  a = isotype.control.name.vec
  b = rownames(empty_drop_matrix)
  c = rownames(cell_protein_matrix)
  if (!isTRUE(all.equal(rownames(cell_protein_matrix), rownames(empty_drop_matrix)))){

    # the number of rows in cell and background matrix are not the same:
    # stop and print base error message
    stopifnot(isTRUE(all.equal(nrow(cell_protein_matrix), nrow(empty_drop_matrix))))

    # some rows have different names in cell and/or background matrix:
    # stop and print difference in both cell and background
    diff = c(setdiff(c,b), setdiff(b,c))
    if (length(diff) > 0) {
      stop(paste0('rows of cell and background matrices have mis-matching names: \n', diff))
    }

    # no difference in elements c,b but rownames are not equal-the rows are the same but in different order:
    # warn and match background matrix rows to cell matrix rows.
    if (length(diff < 0)) {
      warning('rows (proteins) of cell_protein_matrix and empty_drop_matrix are not in the same order')
      rmatch = match(x = rownames(cell_protein_matrix), table = rownames(empty_drop_matrix) )
      empty_drop_matrix = empty_drop_matrix[rmatch, ]
      print('reordered empty_drop_matrix rows to match cell_protein_matrix rows')
    }
  }

  # step II denoising argument checks

  # isotype.control.name.vec is specified but some isotypes are not in input matrices
  # there is a name formatting error: print mismatched elements of isotype.control.name.vec in matrix rows
  if (!is.null(a) & !isTRUE(all(a %in% b)) & !isTRUE(all(a %in% c))){
    stop(paste0("some elements of isotype.control.name.vec are not in input data rownames: \n",
                'cell_protein_matrix - ', setdiff(a,b), ' \nempty_drop_matrix - ', setdiff(a,c))
    )
  }

  # denoise.counts = FALSE:
  # remind that isotypes are not being used, automatically set isotype-related args to FALSE, NULL
  if (isFALSE(denoise.counts)) {
    print(paste0("Running step I ambient correction and log transformation, not running step II removal of cell to cell technical noise.",
                 " Setting use.isotype.control and isotype.control.name.vec to FALSE and NULL"))
    use.isotype.control = FALSE
    isotype.control.name.vec = NULL
  }

  # isotype control related errors
  # Try to detect isotypes in matrix to recommend in error messages.
  iso_detect = rownames(cell_protein_matrix)[grepl('sotype|Iso|iso|control|CTRL|ctrl|Ctrl', rownames(cell_protein_matrix))]

  # use.isotype.control set to TRUE but isotype.control.name.vec was not specified:
  # stop and set isotype.control.name.vec to names of isotype control rows
  if (isTRUE(use.isotype.control) & is.null(isotype.control.name.vec)) {
    stop('if use.isotype.control = TRUE, set isotype.control.name.vec to names of isotype control rows')
    if (length(iso_detect) > 0) {
      print('potential isotype controls detected: ')
      print(iso_detect)
    }
  }

  # denoise.counts = TRUE but use.isotype.control = FALSE:
  # warn. Try to detect isotypes in matrix; recommend isotype usage during step II.
  if (isTRUE(denoise.counts) & isFALSE(use.isotype.control)) {
    warning('denoise.counts = TRUE with use.isotype.control = FALSE not recommended if isotype controls are available.\n',
            ' If data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE\n',
            ' and set `isotype.control.name.vec` to a vector of isotype control rownames from cell_protein_matrix'
    )
    if (length(iso_detect) > 0) {
      print('potential isotype controls detected: ')
      print(iso_detect)
    }
  }

  # Normalization function
  # log transform raw cell and background matrices with default or user defined pseudocount
  adt = cell_protein_matrix %>% as.matrix()
  adtu = empty_drop_matrix %>% as.matrix()
  if(isTRUE(define.pseudocount)) {
    adtu_log = log(adtu + pseudocount.use)
    adt_log = log(adt + pseudocount.use)
  } else {
    adtu_log = log(adtu + 10)
    adt_log = log(adt + 10)
  }
  # dsb step I rescale cells based on expected noise in background drops
  print("correcting ambient protein background noise")
  mu_u = apply(adtu_log, 1 , mean)
  sd_u = apply(adtu_log, 1 , sd)
  norm_adt = apply(adt_log, 2, function(x) (x  - mu_u) / sd_u)

  # dsb step II calculate dsb technical component and regress out of ambient corrected values
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

  # return object formatting if return.stats = TRUE
  # return cell level dsb internal stats and protein level stats with the normalized values
  if(isTRUE(return.stats) & isTRUE(denoise.counts)) {
    print('returning list; access normalized matrix with x$dsb_normalized_matrix, protein stats list with x$protein_stats')

    technical_stats = cbind(t(noise_matrix), dsb_technical_component = noise_vector)
    protein_stats = list('background matrix stats' = data.frame(background_mean = mu_u, background_sd = sd_u),
                         'cell matrix stats' = data.frame(cell_mean = apply(adt_log, 1 , mean), cell_sd = apply(adt_log, 1 , sd)),
                         'dsb normalized stats' = data.frame(dsb_mean = apply(norm_adt, 1 , mean), dsb_sd = apply(adt_log, 1 , sd))
    )
    ret_obj = list(
      'dsb_normalized_matrix' = norm_adt,
      'technical_stats' = technical_stats,
      'protein_stats' = protein_stats
    )
    return(ret_obj)
  }
  if(isTRUE(return.stats) & isFALSE(denoise.counts)) {
    print('returning list; access normalized matrix with x$dsb_normalized_matrix, dsb and protein stats with x$protein_stats')

    protein_stats = list('background matrix stats' = data.frame(background_mean = mu_u, background_sd = sd_u),
                         'cell matrix stats' = data.frame(cell_mean = apply(adt_log, 1 , mean), cell_sd = apply(adt_log, 1 , sd)),
                         'dsb normalized stats' = data.frame(dsb_mean = apply(norm_adt, 1 , mean), dsb_sd = apply(adt_log, 1 , sd))
    )
    ret_obj = list(
      'dsb_normalized_matrix' = norm_adt,
      'protein_stats' = protein_stats
    )
    return(ret_obj)
  }
  # default: return a matrix of the normalized values
  if(isFALSE(return.stats)) {
    return(norm_adt)
  }
}


