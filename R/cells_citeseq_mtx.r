#' @title small example CITE-seq protein dataset for 87 surface protein in 2872 cells
#'
#' @references Kotliarov et. al. 2020 Nat. Medicine
#'
#' @description A matrix of raw UMI counts for surface proteins for surface proteins measured with CITE-seq antibodies. This data is used for example scripts in the dsb package. Raw data was processed with CITE-seq-Count.
#'
#' @format An R matrix, rows: 87 proteins, columns: 2872 cells
#' \describe{
#' \item{cells_citeseq_mtx}{R matrix of cells by proteins; a random distribution of a maximum of 100 cells per cluster from the 30 clusters reported in Kotliarov et. al. 2020}
#' }
"cells_citeseq_mtx"
