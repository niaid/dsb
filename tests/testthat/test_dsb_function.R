# unit testing for dsb function
context("testing dsb norm")
library(testthat)

testthat::test_that(desc = "run dsb on example data", code = {
  # initialize
  set.seed(1)
  isotype_names = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT",
                    "MouseIgG2akappaisotype_PROT","RatIgG2bkIsotype_PROT")

  # test appropriate errors and warnings
  # no isotype control name vector specified with use.isotype.control = TRUE
  expect_error(
    object = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx,
                                 empty_drop_matrix = empty_drop_citeseq_mtx,
                                 denoise.counts = TRUE,
                                 use.isotype.control = TRUE)

  )
  # incorrectly specified isotype control names (not matching rownames of raw input matrix)
  expect_error(
    object = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx,
                                 empty_drop_matrix = empty_drop_citeseq_mtx,
                                 denoise.counts = TRUE,
                                 use.isotype.control = TRUE,
                                 isotype.control.name.vec = c('incorrect_isotype1','RatIgG2bkIsotype_PROT' )
                                 )
  )
  # error if return stats TRUE but no denoising step run
  expect_error(
    object = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx,
                                 empty_drop_matrix = empty_drop_citeseq_mtx,
                                 denoise.counts = FALSE, return.stats = TRUE)

  )

  # warn if isotypes not specified and denoising
  expect_warning(
    object = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:100],
                                 empty_drop_matrix = empty_drop_citeseq_mtx[ ,1:100],
                                 denoise.counts = TRUE, use.isotype.control = FALSE)
      )

  # test function runs with options specified for 100% code coverage
  result = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx,
                               empty_drop_matrix = empty_drop_citeseq_mtx,
                               denoise.counts = FALSE,
                               define.pseudocount = FALSE,
                               use.isotype.control = FALSE)
  # test result value
  expect_gt(object = mean(rowMeans(result)), 1)

  # test mu1 only denoising
  tester = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:50],
                                 empty_drop_matrix = empty_drop_citeseq_mtx[ ,1:50],
                                 denoise.counts = TRUE,
                                 use.isotype.control = FALSE)

  # test full options on data subset
  result = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:400],
                               empty_drop_matrix = empty_drop_citeseq_mtx[ ,1:400],
                               define.pseudocount = TRUE,
                               pseudocount.use = 0.5,
                               use.isotype.control = TRUE,
                               isotype.control.name.vec = isotype_names,
                               return.stats = TRUE,
                               quantile.clipping = TRUE)

  # check stats return andcorrelaiton of technical component
  dsb_stats = result$dsb_stats
  result = result$dsb_normalized_matrix
  # test vals
  expect_gt(object = mean(rowMeans(result)), 1)
})
