# dsb 1.0.3

This is a patch with no function updates. Startup message changed for linux build "Note". 

# dsb 1.0.2

This is a patch with no function updates. Readme is slimmed down with links to the full vignettes on CRAN. Updated citation information. Updated license to CC0 per CRAN request. 

# dsb 1.0.1 

*this is a patch to fix a vignette rendering issue in the 1.0.0 release*  

**Release notes from dsb 1.0.0:**

This is the feature complete version of dsb being released with the publication of our preprint in Nature Communications.

### Enhancements  

- A new function added to implement a method to normalize ADTs for datasets where empty droplets are not available: `ModelNegativeADTnorm`.  
- updated internal code for `DSBNormalizeProtein` to implement additional error checking and messages / warnings during function run.  

### New vignettes
"Additional Topics - quantile.clipping - scale.factor - Python and Bioc - multiplexing - multi batch - FAQ"  
"Normalizing ADTs for datasets without empty droplets"  
"Understanding how the dsb method works"  

# dsb 1.0.0 

This is the feature complete version of dsb being released with the publication of our preprint in Nature Communications. 

### Enhancements  

- A new function added to implement a method to normalize ADTs for datasets where empty droplets are not available: `ModelNegativeADTnorm`.  
- updated internal code for `DSBNormalizeProtein` to implement additional error checking and messages / warnings during function run.  

### New vignettes
"Additional Topics - quantile.clipping - scale.factor - Python and Bioc - multiplexing - multi batch - FAQ"  
"Normalizing ADTs for datasets without empty droplets"  
"Understanding how the dsb method works"  



# dsb 0.3.0 

### Enhancements  

Additional error checking on input cell and background matrices:  
 - stop if input matrix rows are not equivalent length https://github.com/niaid/dsb/issues/29
 - stop if any names in input matrices are not equivalent  
 - warn if the rows are not in the same order and reorder to match.  

Improve warning and error messages for isotype control name matching issues.

Advanced users can now examine protein (mean sd) and cell level stats to output if `return.stats = TRUE`. If `denoise.counts = FALSE`. The output includes protein level stats. If `denoise.counts = TRUE`, the output contains cell level stats including the dsb technical component and derivative variables used in step II.  

dsb with full options specified: 

```
# full options: defined pseudocount, isotype controls, outlier clip, stats
  result = DSBNormalizeProtein(cell_protein_matrix = cells_citeseq_mtx[ ,1:100],
                               empty_drop_matrix = empty_drop_citeseq_mtx,
                               define.pseudocount = TRUE,
                               pseudocount.use = 5,
                               use.isotype.control = TRUE,
                               isotype.control.name.vec = rownames(cells_citeseq_mtx)[grepl(
                                                          rownames(cells_citeseq_mtx), pattern = 'otyp')],
                               quantile.clipping = TRUE,
                               return.stats = TRUE
  )
  
# normalized data 
result$dsb_normalized_matrix

# protein-level statistics
result$protein_stats

# cell-level statistics
result$technical_stats


```
If `return.stats = FALSE`, the output is a R matrix equivalent to `result$dsb_normalized_matrix` above. 


Add unit tests for changes.

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

