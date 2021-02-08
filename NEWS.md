# dsb 0.1.0

• dsb is now hosted on CRAN: [go to dsb on CRAN](https://cran.r-project.org/web/packages/dsb/index.html)

• Documentation is improved from the beta release based on user feedback with a clearer workflow for defining and performing QC on background droplets and cells defined from raw protein UMI data. Additional code is added for integrating dsb with Seurat, Scanpy in Python, and  Bioconductor's SingleCellExperiment class. Added an FAQ section based on user questions and a workflow for both multiplexing experiments or for non-multiplexing / single lane experiments. See [the updated documentation on CRAN](https://cran.r-project.org/web/packages/dsb/vignettes/dsb_normalizing_CITEseq_data.html)

• The current dsb package defaults are: `denoise.counts = TRUE` and `use.isotype.control = TRUE`. Isotype contorls are not required for normalization and protein panels of any size can be normalized. See package documentation.  



