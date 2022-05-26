
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/dsb)](https://CRAN.R-project.org/package=dsb)
<!-- badges: end -->

# <a href='https://CRAN.R-project.org/package=dsb/'><img src='man/figures/sticker2.png' align="right" width="150" /></a> dsb: a method for normalizing and denoising antibody derived tag data from CITE-seq, ASAP-seq, TEA-seq and related assays.

The dsb R package is available on [**CRAN: latest dsb
release**](https://CRAN.R-project.org/package=dsb)  
to install in R: `install.packages('dsb')`

#### Vignettes:

1.  [**Using dsb in an end to end CITE-seq workflow including WNN
    clustering**](https://CRAN.R-project.org/package=dsb/vignettes/end_to_end_workflow.html)  
2.  [**How the dsb method
    works**](https://CRAN.R-project.org/package=dsb/vignettes/understanding_dsb.html)  
3.  [**Normalizing ADTs if empty drops are not
    available**](https://CRAN.R-project.org/package=dsb/vignettes/no_empty_drops.html)  
4.  [**FAQ
    etc.**](https://CRAN.R-project.org/package=dsb/vignettes/additional_topics.html)

[**Mulè, Martins, and Tsang, Nature Communications
(2022)**](https://www.nature.com/articles/s41467-022-29356-8) describes
this method and deconvolution of ADT noise sources.

dsb is also available in [**Python** through
*muon*](https://muon.readthedocs.io/en/latest/omics/citeseq.html)

Check out [**Recent Publications**](#pubications) that used this method
for ADT normalization.

In the first Vignette, we demonstrate an end-to-end basic CITE-seq
analysis starting from UMI count alignment output files from Cell
Ranger. Standard output files from Cell Ranger are perfectly set up to
use dsb. Our method is also compatible with any alignment tool; see:
[**using other alignment tools**](#otheraligners). We load unfiltered
UMI data containing cells and empty droplets, perform QC on cells and
background droplets, normalize with dsb, and demonstrate protein-based
clustering and multimodal RNA+Protein joint clustering using dsb
normalized values with Seurat’s Weighted Nearest Neighbor method.

## Background and motivation <a name="background_motivation"></a>

Protein data derived from sequencing antibody derived tags (ADTs) in
CITE-seq and other related assays has substantial background noise.
[**Our paper**](https://www.nature.com/articles/s41467-022-29356-8)
outlines experiments and analysis designed to dissect sources of noise
in ADT data we used to developed our method. We found all experiments
measuring ADTs capture protein-specific background noise because ADT
reads in empty / background drops (outnumbering cell-containing droplets
\> 10-fold in all experiments) were highly concordant with ADT levels in
unstained spike-in cells. We therefore utilize background droplets which
capture the *ambient component* of protein background noise to correct
values in cells. We also remove technical cell-to-cell variations by
defining each cell’s dsb “technical component”, a conservative
adjustment factor derived by combining isotype control levels with each
cell’s specific background level fitted with a single cell model.

## Installation and quick overview <a name="installation"></a>

The method is carried out in a single step with a call to the
`DSBNormalizeProtein()` function.  
`cells_citeseq_mtx` - a raw ADT count matrix `empty_drop_citeseq_mtx` -
a raw ADT count matrix from non-cell containing empty / background
droplets.  
`denoise.counts = TRUE` - implement step II to define and remove the
‘technical component’ of each cell’s protein library.  
`use.isotype.control = TRUE` - include isotype controls in the modeled
dsb technical component.

``` r
# install.packages('dsb')
library(dsb)

adt_norm = DSBNormalizeProtein(
  cell_protein_matrix = cells_citeseq_mtx, 
  empty_drop_matrix = empty_drop_citeseq_mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = rownames(cells_citeseq_mtx)[67:70]
  )
```

For a full tutorial vignette using 10X genomics data demonstrating how
to quality control cells and empty droplets, normalize ADTs with dsb and
ue these values with joint mRNA and protein clustering to generate the
joint map below, please see [**the main vignette on
CRAN**](https://CRAN.R-project.org/package=dsb/vignettes/end_to_end_workflow.html)

<img src="man/figures/multimodal_heatmap.png" />

### Some recent publications using dsb <a name="pubications"></a>

Publications from outside institutes  
[Jardine et al. *Nature*
(2021)](https://www.nature.com/articles/s41586-021-03929-x)  
[Mimitou et. al *Nature Biotechnology*
(2021)](https://www.nature.com/articles/s41587-021-00927-2)  
[COMBAT consortium et al. *Cell*
(2022)](https://www.cell.com/cell/pdf/S0092-8674(22)00070-8.pdf)

Publications from the Tsang lab  
[Kotliarov et al. *Nature Medicine*
(2020)](https://www.nature.com/articles/s41591-020-0769-8)  
[Liu et al. *Cell*
(2021)](https://www.cell.com/cell/pdf/S0092-8674(21)00168-9.pdf)

### using other alignment algorithms <a name="otheraligners"></a>

dsb was developed prior to 10X Genomics supporting CITE-seq or hashing
data and we routinely use other alignment pipelines.

A note on alignment and how to use dsb with Cell Ranger is detailed in
the main vignette. Cells and empty droplets are used by default by dsb.

<img src="man/figures/readme_cheatsheet.png" />

To use dsb properly with CITE-seq-Count you need to align background.
One way to do this is to set the `-cells` argument to \~ 200000. That
will align the top 200000 barcodes in terms of ADT library size, making
sure you capture the background. Please refer to [**CITE-seq-count
documentation**](https://hoohm.github.io/CITE-seq-Count/Running-the-script/)

``` bash
CITE-seq-Count -R1 TAGS_R1.fastq.gz  -R2 TAGS_R2.fastq.gz \
 -t TAG_LIST.csv -cbf X1 -cbl X2 -umif Y1 -umil Y2 \
  -cells 200000 -o OUTFOLDER
```

If you already aligned your mRNA with Cell Ranger or something else but
wish to use a different tool like kallisto or Cite-seq-count for ADT
alignment, you can provide the latter with whitelist of cell barcodes to
align. A simple way to do this is to extract all barcodes with at least
k mRNA where we set k to a tiny number to retain cells *and* cells
capturing ambient ADT reads:

``` r
library(Seurat)
umi = Read10X(data.dir = 'data/raw_feature_bc_matrix/')
k = 3 
barcode.whitelist = 
  rownames(
    CreateSeuratObject(counts = umi,
                       min.features = k,  # retain all barcodes with at least k raw mRNA
                       min.cells = 800, # this just speeds up the function by removing genes. 
                       )@meta.data 
    )

write.table(barcode.whitelist,
file =paste0(your_save_path,"barcode.whitelist.tsv"), 
sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
```

With the example dataset in the vignette this retains about 150,000
barcodes.

Now you can provide that as an argument to `-wl` in CITE-seq-count to
align the ADTs and then proceed with the dsb analysis example.

``` bash
CITE-seq-Count -R1 TAGS_R1.fastq.gz  -R2 TAGS_R2.fastq.gz \
 -t TAG_LIST.csv -cbf X1 -cbl X2 -umif Y1 -umil Y2 \
  -wl path_to_barcode.whitelist.tsv -o OUTFOLDER
```

This whitelist can also be provided to Kallisto.  
[kallisto bustools
documentation](https://www.kallistobus.tools/tutorials/kb_kite/python/kb_kite/)

``` bash
kb count -i index_file -g gtf_file.t2g -x 10xv3 \
-t n_cores -w path_to_barcode.whitelist.tsv -o output_dir \
input.R1.fastq.gz input.R2.fastq.gz
```

Next one can similarly define cells and background droplets empirically
with protein and mRNA based thresholding as outlined in the main
tutorial.

### A note on Cell Ranger –expect-cells <a name="cellranger"></a>

Note *whether or not you use dsb*, if you want to define cells using the
`filtered_feature_bc_matrix` file, you should make sure to properly set
the Cell Ranger `--expect-cells` argument roughly equal to the estimated
cell recovery per lane based on number of cells you loaded in the
experiment. see [the note from 10X about
this](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview#cell_calling).
The default value of 3000 is relatively low for modern experiments. Note
cells and empty droplets can also be defined directly from the
`raw_feature_bc_matrix` using any method, including simple protein and
mRNA library size based thresholding because this contains all droplets.

Topics covered in other vignettes on CRAN: **Integrating dsb with
Bioconductor, integrating dsb with python/Scanpy, Using dsb with data
lacking isotype controls, integrating dsb with sample multiplexing
experiments, using dsb on data with multiple batches, advanced usage -
using a different scale / standardization based on empty droplet levels,
returning internal stats used by dsb, outlier clipping with the
quantile.clipping argument, other FAQ.**
