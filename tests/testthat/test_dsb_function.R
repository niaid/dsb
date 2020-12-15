context("testing dsb norm")
library(testthat)

testthat::test_that(desc = "run dsb on example data", code = {
  set.seed(1)

  # test the dsb normalization function on the example data with defaults
  result = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx,
                               empty_drop_matrix = empty_drop_citeseq_mtx,
                               denoise.counts = FALSE,
                               use.isotype.control = FALSE)

  # the mean of the rowmeans is 1.61, make sure returned value is > 1
  expect_gt(object = mean(rowMeans(result)), 1)

  # make sure returned dimensions are correct
  expect_equal(ncol(result), expected = 2872)
  expect_equal(nrow(result), expected = 87)

  # test full options: defined pseudocount and isotype controls.
  result = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx,
                               empty_drop_matrix = empty_drop_citeseq_mtx,
                               define.pseudocount = TRUE,
                               pseudocount.use = 10,
                               use.isotype.control = TRUE,
                               isotype.control.name.vec = c("Mouse IgG2bkIsotype_PROT", "MouseIgG1kappaisotype_PROT",
                                                            "MouseIgG2akappaisotype_PROT", "RatIgG2bkIsotype_PROT"))

  # make sure returned dimensions are correct and mean of rowmeans is greater than 1
  expect_equal(ncol(result), expected = 2872)
  expect_equal(nrow(result), expected = 87)
  expect_gt(object = mean(rowMeans(result)), 1)


})
