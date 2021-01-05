#' @title small example CITE-seq protein dataset for 87 surface protein in 8005 empty droplets
#'
#' @references Kotliarov et. al. 2020 Nat. Medicine
#'
#' @description A matrix of cells by antibodies Raw CITEseq data used for example scripts of the dsb package. raw data processed with CITE-seq-count https://hoohm.github.io/CITE-seq-Count/
#'
#' @format An R matrix, rows: 87 proteins, columns: 8005 empty droplets.
#' \describe{
#' \item{empty_drop_citeseq_mtx}{R matrix of empty / background droplets from a CITE-seq experiment. Negative drops were called on cell hashing data with Seurat's HTODemux function and cross referencing mRNA in driplets against patient genotypes with Demuxlet. Ambiguous drops, and with less than 80 unique mRNA were removed. This is used for robust estimation of the background distribution of each protein}
#' }
"empty_drop_citeseq_mtx"
