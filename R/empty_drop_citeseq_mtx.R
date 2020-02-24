#' @title sample CITEseq protein data 87 protein by 8005 empty droplets
#'
#' @description A matrix of cells by antibodies Raw CITEseq data used for example scripts of the dsb package. raw data processed with CITE-seq-count https://hoohm.github.io/CITE-seq-Count/
#'
#' @format A matrix 87 protein by 8005 empty droplets confirmed by cell hashing and cross referencing sample genotypes.
#' \describe {
#' \item{empty_drop_citeseq_mtx}{is a matrix of empty droplets from a CITEseq experiment using cell hashing antibodies called by Seurat's HTODemux function as negative cells, called by Demuxlet via referencing patient genotype data as ambiguous drops, and with less than 80 unique genes detected. This is used for robust estimation of the background distribution of each protein.}
#' }
"empty_drop_citeseq_mtx"
