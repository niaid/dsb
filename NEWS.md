# dsb 0.2.0

### Enhancements  

#### Return dsb internal stats  
Some advanced users may desire to look into the latent noise variables estimated by dsb. Setting `return.stats = TRUE` now returns a list, the first element is the dsb normalized ADT values, the second element is a dataframe of noise variables for each cell, including the derived "dsb technical component" along with its constituent variables, the background mean and the isotype control protein values for each cell.  

#### quantile clipping  
Rarely, a small minority of cells can have outlier values for a single protein, for example expression at level -25 after normalization. These values can now be clipped to be the lowest / highest quantile by default c(0.001, 0.9995). To use this feature set `quantile.clipping = TRUE` and quantile.clip to be a vector of lowest and highest quantile values to use; defaults to c(0.001, 0.9995).  

#### Documentation  
Documentation has been overhauled, the main vignette now includes simpler estimation of cells vs empty drops by using the EmptyDrops method that is now implemented by Cell Ranger by default. A workflow for multimodal (CITE-seq) single cell analysis using the dsb method for normalization is provided including loading raw data, quality control, clustering and multiple versions of the Seurat Weighted Nearest Neighbors joint mRNA and protein clustering method. Added to FAQ section in vignette. 

# dsb 0.1.0

• dsb is now hosted on CRAN: [go to dsb on CRAN](https://CRAN.R-project.org/package=dsb)

• Documentation is improved from the beta release based on user feedback with a clearer workflow for defining and performing QC on background droplets and cells defined from raw protein UMI data. Additional code is added for integrating dsb with Seurat, Scanpy in Python, and  Bioconductor's SingleCellExperiment class. Added an FAQ section based on user questions and a workflow for both multiplexing experiments or for non-multiplexing / single lane experiments. See the updated documentation on CRAN.

• The current dsb package defaults are: `denoise.counts = TRUE` and `use.isotype.control = TRUE`. Isotype controls are not required for normalization. See package documentation.  

