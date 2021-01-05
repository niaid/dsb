#' @title small example CITE-seq protein dataset for 87 surface protein in 2872 cells
#'
#' @references Kotliarov et. al. 2020 Nat. Medicine
#'
#' @description A matrix of cells by antibodies Raw CITEseq data used for example scripts of the dsb package. raw data processed with CITE-seq-count https://hoohm.github.io/CITE-seq-Count/
#'
#' @format An R matrix, rows: 87 proteins, columns: 2872 cells
#' \describe{
#' \item{cells_citeseq_mtx}{R matrix of cells by proteins; a random distribution of a maximum of 100 cells per cluster from the 30 clusters reported in Kotliarov et. al. 2020}
#' }
"cells_citeseq_mtx"
