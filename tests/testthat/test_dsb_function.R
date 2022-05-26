# unit testing for dsb function
library(testthat)
context("testing dsb norm")

testthat::test_that(desc = "run dsb on example data", code = {
  set.seed(1)

  # test error handling
  cell = dsb::cells_citeseq_mtx
  empty = dsb::empty_drop_citeseq_mtx

  # test mis-specified isotypes
  isotypes = rownames(cell)[grepl(x = rownames(cell), pattern = 'otyp')]
  isotypes[2] = 'mis-specified_isotype'
  testthat::expect_error(
    norm = DSBNormalizeProtein(cell_protein_matrix = cell[ ,1:50],
                               empty_drop_matrix = empty,
                               # input mis specified isotype in function
                               use.isotype.control = TRUE,
                               isotype.control.name.vec = isotypes)
  )

  # test scrambled rownames
  scramble_row = sample(rownames(empty),size = 87, replace = FALSE)
  cell = cell[scramble_row, ]
  testthat::expect_warning(
    DSBNormalizeProtein(cell_protein_matrix = cell[ ,1:50],
                        empty_drop_matrix = empty,
                        use.isotype.control = FALSE)
  )

  #test mismatched rows
  rownames(cell)[5] = 'wrong name'
  testthat::expect_error(
    norm = DSBNormalizeProtein(cell_protein_matrix = cell,
                               empty_drop_matrix = empty,
                               use.isotype.control = FALSE)
  )

  # test a missing row
  cell = cell[-5, ]
  testthat::expect_error(
    norm = DSBNormalizeProtein(cell_protein_matrix = cell,
                               empty_drop_matrix = empty,
                               use.isotype.control = FALSE)
  )

  # test DSBNormalizeProtein without step II
  result = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:100],
                               empty_drop_matrix = empty_drop_citeseq_mtx,
                               denoise.counts = FALSE, use.isotype.control = FALSE)

  # test non default scaling
  r2 =
    DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:100],
                         empty_drop_matrix = empty_drop_citeseq_mtx,
                         define.pseudocount = TRUE,
                         pseudocount.use = 5,
                         use.isotype.control = TRUE,
                         isotype.control.name.vec = rownames(cells_citeseq_mtx)[grepl(
                           rownames(cells_citeseq_mtx), pattern = 'otyp')],
                         quantile.clipping = TRUE, scale.factor = 'mean.subtract',
                         return.stats = FALSE)
  r2mean = mean(rowMeans(r2))
  testthat::expect_lt(object = r2mean, 1)

  # test mis-specified non-default scaling
  testthat::expect_error(
      DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:100],
                          empty_drop_matrix = empty_drop_citeseq_mtx,
                          scale.factor = 'not.an.option')
  )

  # test output return value of function run on example data.
  testthat::expect_gt(object = mean(rowMeans(result)), 1)

  # test output dimensions of returned matrix
  testthat::expect_equal(ncol(result), expected = 100)
  testthat::expect_equal(nrow(result), expected = 87)

  # test full options: defined pseudocount, isotype controls, outlier clip, stats
  result = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:100],
                               empty_drop_matrix = empty_drop_citeseq_mtx,
                               define.pseudocount = TRUE,
                               pseudocount.use = 5,
                               use.isotype.control = TRUE,
                               isotype.control.name.vec = rownames(cells_citeseq_mtx)[grepl(rownames(cells_citeseq_mtx), pattern = 'otyp')],
                               quantile.clipping = TRUE,
                               return.stats = TRUE
  )
  # test returned value is a list with return.stats = TRUE
  testthat::expect_type(object = result, type = 'list')


  # test return of dsb tech componant with denoise counts TRUE but use isotype FALSE and return stats TRUE
  result = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:100],
                               empty_drop_matrix = empty_drop_citeseq_mtx,
                               define.pseudocount = TRUE,
                               pseudocount.use = 5,
                               use.isotype.control = FALSE,
                               quantile.clipping = TRUE,
                               return.stats = TRUE
  )

  # test that a list is returned
  testthat::expect_type(object = result, type = 'list')


  # test returned matrix value and dimensions
  normprot = result$dsb_normalized_matrix
  testthat::expect_equal(ncol(normprot), expected = 100)
  testthat::expect_equal(nrow(normprot), expected = 87)
  testthat::expect_gt(object = mean(rowMeans(normprot)), 1)

  #######################################################
  # test ModelNegativeADTnorm
  result =
    ModelNegativeADTnorm(cell_protein_matrix = cells_citeseq_mtx[ ,1:100],
                       define.pseudocount = TRUE,
                       pseudocount.use = 5,
                       use.isotype.control = TRUE,
                       isotype.control.name.vec = rownames(cells_citeseq_mtx)[grepl(
                         rownames(cells_citeseq_mtx), pattern = 'otyp')],
                       quantile.clipping = TRUE,
                       quantile.clip = c(0.0001, 0.9999),
                       return.stats = TRUE
                       )
  testthat::expect_type(object = result, type = 'list')
  normprot = result$dsb_normalized_matrix
  testthat::expect_equal(ncol(normprot), expected = 100)
  testthat::expect_equal(nrow(normprot), expected = 87)
  testthat::expect_gt(object = mean(rowMeans(normprot)), 0.1)

  # simple implementation of ModelNegativeADTnorm
  result =ModelNegativeADTnorm(cell_protein_matrix = cells_citeseq_mtx[ ,1:100], denoise.counts = FALSE)
  normprot = result
  testthat::expect_equal(ncol(normprot), expected = 100)
  testthat::expect_equal(nrow(normprot), expected = 87)
  testthat::expect_gt(object = mean(rowMeans(normprot)), 0.1)

  # mu1 denoising implementation of ModelNegativeADTnorm
  testthat::expect_warning(
    ModelNegativeADTnorm(cell_protein_matrix = cells_citeseq_mtx[ ,1:50],
                         denoise.counts = TRUE,
                         use.isotype.control = FALSE)
  )



})

