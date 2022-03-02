
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/dsb)](https://CRAN.R-project.org/package=dsb)
<!-- badges: end -->

# <a href='https://mattpm.github.io/dsb/'><img src='man/figures/sticker2.png' align="right" width="175" /></a> dsb: a method for normalizing and denoising antibody derived tag data from CITE-seq, ASAP-seq, TEA-seq and related assays.

The dsb R package is available on CRAN: to install:
`install.packages('dsb')`  
[**CRAN: latest dsb release, documentation and
vignettes**](https://CRAN.R-project.org/package=dsb)

dsb is also available for *Python* with muon:  
[**muon
documentation**](https://muon.readthedocs.io/en/latest/omics/citeseq.html)

Please cite the dsb manuscript if you used our software or found the
experimental and modeling derivation of ADT noise sources in our paper
helpful:  
[**dsb
manuscript**](https://www.biorxiv.org/content/10.1101/2020.02.24.963603v3)

Recent publications using the dsb method  
[**Recent Publications**](#pubications)

See news.md for updates.  
[news.md](https://github.com/niaid/dsb/blob/master/NEWS.md)

Below we demonstrate an end-to-end basic CITE-seq analysis starting from
UMI count alignment output files from Cell Ranger. Standard output files
from Cell Ranger are perfectly set up to use dsb. Our method is also
compatible with any alignment tool; see: [using other alignment
tools](#otheraligners). We load unfiltered UMI data containing cells and
empty droplets, perform QC on cells and background droplets, normalize
with dsb, and demonstrate protein-based clustering and multimodal
RNA+Protein joint clustering using dsb normalized values with Seurat’s
Weighted Nearest Neighbor method.

## Table of Contents

1.  [Background and motivation](#background_motivation)
2.  [Installation and quick overview](#installation)
3.  [Download public 10X Genomics data](#tutorial)
4.  [Step 1 A note on alignment of ADTs](#step1)
5.  [Step 2 Load RNA and ADT data and define droplet quality control
    metadata](#step2)
6.  [Step 3 Quality control on cell-containing and background
    droplets](#step3)
7.  [Step 4 Normalize protein data with the DSBNormalizeProtein
    Function](#step4)
8.  [Integrating dsb with Seurat](#seurat)
9.  [Clustering cells based on dsb normalized protein using
    Seurat](#seurat_clustering)
10. [dsb derived cluster interpretation](#interpretation)
11. [Weighted Nearest Neighbor multimodal clustering using dsb
    normalized values with Seurat](#wnn)

Topics covered in other vignettes on CRAN: **Integrating dsb with
Bioconductor, integrating dsb with python/Scanpy, Using dsb with data
lacking isotype controls, integrating dsb with sample multiplexing
experiments, using dsb on data with multiple batches, advanced usage -
using a different scale / standardization based on empty droplet levels,
returning internal stats used by dsb, outlier clipping with the
quantile.clipping argument, other FAQ.**

## Background and motivation <a name="background_motivation"></a>

Protein data derived from sequencing antibody derived tags (ADTs) in
CITE-seq and other related assays suffers from substantial background
noise. We performed experiments designed to dissect protein noise in
CITE-seq data; our method dsb is based on 3 key findings detailed in our
[**paper**](https://www.biorxiv.org/content/10.1101/2020.02.24.963603v3).
Briefly:

1)  Based on unstained cell spike-in experiments, we found a major
    source of protein-specific background noise comes from ambient,
    unbound antibody encapsulated in droplets.

2)  All experiments capture ambient protein specific background noise
    with a “built-in” control–ADT reads in empty / background drops
    (which outnumber cell-containing droplets \> 10-fold in all
    experiments) were highly concordant with ADT levels in unstained
    spike-in cells. These background droplets thus capture the *ambient
    component* of protein background noise. Our method (step I) uses
    background droplets to correct for this protein-specific noise.

3)  Technical cell-to-cell variations (i.e. capture and RT efficiency,
    non-specific binding, sequencing depth) should be normalized across
    cells. Since library size correction is not appropriate for ADTs we
    derived a more conservative technical factor. In step II, we define
    and remove variation associated with each cell’s dsb “technical
    component”. The technical component is defined after ambient
    correction by combining (by default) isotype control levels with the
    cell’s specific background level fitted with a single cell model.
    Combining these variables captured important statistical properties
    not exhibited by using individual variables such as isotype controls
    alone.

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
install.packages('dsb')
library(dsb)

adt_norm = DSBNormalizeProtein(
  cell_protein_matrix = cells_citeseq_mtx, 
  empty_drop_matrix = empty_drop_citeseq_mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = rownames(cells_citeseq_mtx)[67:70]
  )
```

## Download public 10X Genomics data <a name="tutorial"></a>

Download BOTH the Feature / cell matrix (filtered) and Feature / cell
matrix (raw) count matrices here: [public 10X Genomics CITE-seq
data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3_nextgem).

## Step 1 A note on alignment of ADTs <a name="step1"></a>

Here is a visual showing for the workflow we will follow starting from
alignment and proceeding through QC and normalization  
<img src="man/figures/readme_cheatsheet.png" />

The Cell Ranger `raw_feature_bc_matrix` includes every possible cell
barcode (columns) x genes / ADT (rows); about 7 Million barcodes for the
V3 assay. A tiny subset of the columns in the `raw_feature_bc_matrix`
contain your cells–those are also in the `filtered_feature_bc_matrix`
file. Another subset (often ranging 50,000-200,000) of the
`raw_feature_bc_matrix` contain empty droplets capturing ambient, free
floating antibodies, this can be a majority of your ADT sequencing
counts\! Instead of throwing that sequencing data away, we extract empty
droplets capturing ADT using the code below and use them with dsb. This
same workflow also applies to Kallisto and CITE-Seq-Count, see [using
other alignment tools](#otheraligners).

Note *whether or not you use dsb*, if you want to define cells using the
`filtered_feature_bc_matrix` file, you should make sure to properly set
the Cell Ranger `--expect-cells` argument roughly equal to the estimated
cell recovery per lane based on number of cells you loaded in the
experiment. Cell Ranger uses this in part to define cells by using the
[empty drops
method](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y).
The default value is set to **3000-that’s low for most modern
experiments** – see [the note from 10X about
this](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview#cell_calling).
Note cells and empty droplets can also be defined directly from the
`raw_feature_bc_matrix` using any method, including simple protein and
mRNA library size based thresholding because this contains all
droplets.

## Step 2 Load RNA and ADT data and define droplet quality control metadata <a name="step2"></a>

Download and un-compress ‘Feature / cell matrix (filtered)’ and ‘Feature
/ cell matrix (raw)’ which will be automatically renamed
`filtered_feature_bc_matrix` and `raw_feature_bc_matrix`. To follow
along without changing any file paths you can add them to a directory
`data/` in your working directory.

Your working directory should now be structured:  
data/  
  |\_filtered\_feature\_bc\_matrix  
  |\_raw\_feature\_bc\_matrix

Here we use the convenience function from Seurat `Read10X` which will
automatically detect multiple assays and create two element list `Gene
Expression` and `Antibody Capture`.

``` r
library(dsb)

# read raw data using the Seurat function "Read10X" 
raw = Seurat::Read10X("data/raw_feature_bc_matrix/")
cells = Seurat::Read10X("data/filtered_feature_bc_matrix/")

# define cell-containing barcodes and separate cells and empty drops
stained_cells = colnames(cells$`Gene Expression`)
background = setdiff(colnames(raw$`Gene Expression`), stained_cells)

# split the data into separate matrices for RNA and ADT
prot = raw$`Antibody Capture`
rna = raw$`Gene Expression`
```

Now calculate some standard meta data for cells that we will use for
quality control using standard approaches used in scRNAseq analysis
pipelines.

``` r
# create metadata of droplet QC stats used in standard scRNAseq processing
mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below

md = data.frame(
  rna.size = log10(Matrix::colSums(rna)), 
  prot.size = log10(Matrix::colSums(prot)), 
  n.gene = Matrix::colSums(rna > 0), 
  mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
)
# add indicator for barcodes Cell Ranger called as cells
md$drop.class = ifelse(rownames(md) %in% stained_cells, 'cell', 'background')

# remove barcodes with no evidence of capture in the experiment
md = md[md$rna.size > 0 & md$prot.size > 0, ]
```

**we now have 103,075 barcodes in `md` with RNA and protein data** only
about 7,000 contain cells. With dsb we use a majority of this data
instead of discarding everything except the
cells.

## Step 3 Quality control cells and background droplets <a name="step3"></a>

We next visualize the log library size (total reads) of mRNA and protein
coloring values with MT gene proportion to see where in the distribution
potentially low-quality, dying cells are, and where the highest density
(number) of droplets are. We use this to subset the empty drops to
narrow in on the major background distribution – the area of high
density in the leftmost plot. We can see in the 3rd plot from the left
with cells colored by MT gene proportion that some empty drops as
defined by Cell Ranger may be low-quality apoptotic cells (which have
high MT gene content). It is likely inevitable to have some low quality
cells in the ambient matrix but we can filter these out as well by
setting upper and lower cutoffs on the background matrix library size.
The Normalized values are robust to different background thresholds
used, so long as one does not omit this major background population (See
Supplementary Fig. 7 in the paper for more information).
<img src="man/figures/drop_thresholds2.png" />

``` r
background_drops = rownames(
  md[ md$prot.size > 1.5 & 
      md$prot.size < 3 & 
      md$rna.size < 2.5, ]
  ) 
background.adt.mtx = as.matrix(prot[ , background_drops])
```

We next quality control cells in a similar manner with additional
quality control metrics calculated as in any standard scRNAseq analysis,
e.g. see [Luecken *et. al.* 2019 *Mol Syst
Biol*](https://www.embopress.org/doi/full/10.15252/msb.20188746).

``` r
# calculate statistical thresholds for droplet filtering.
cellmd = md[md$drop.class == 'cell', ]

# filter drops with + / - 3 median absolute deviations from the median library size
rna.mult = (3*mad(cellmd$rna.size))
prot.mult = (3*mad(cellmd$prot.size))
rna.lower = median(cellmd$rna.size) - rna.mult
rna.upper = median(cellmd$rna.size) + rna.mult
prot.lower = median(cellmd$prot.size) - prot.mult
prot.upper = median(cellmd$prot.size) + prot.mult

# filter rows based on droplet qualty control metrics
qc_cells = rownames(
  cellmd[cellmd$prot.size > prot.lower & 
         cellmd$prot.size < prot.upper & 
         cellmd$rna.size > rna.lower & 
         cellmd$rna.size < rna.upper & 
         cellmd$mt.prop < 0.14, ]
  )
```

Sanity check: are the number of cells passing QC in line with the
expected recovery from the experiment?

``` r
length(qc_cells)
```

**\[1\] 4096**

Yes. After quality control above we have 4096 cells which is in line
with the ~5000 cells loaded in this experiment.

Now subset the metadata ADT and RNA matrices.

``` r
cell.adt.raw = as.matrix(prot[ , qc_cells])
cell.rna.raw = rna[ ,qc_cells]
cellmd = cellmd[qc_cells, ]
```

## Optional step; remove proteins without staining

Some proteins in an experiment may not work for bioinformatic reasons or
may target a very rare cell population that was absent in the
experiment. Proteins without counts in stained cells can be removed from
both matrices prior to normalization. In this experiment, CD34 has
essentially no data (a maximum raw UMI value of 4 across all cells in
the experiment). We therefore remove it. In many cases, removing
proteins is not necessary, but we recommend checking out your raw data.

``` r
# flter 
pm = sort(apply(cell.adt.raw, 1, max))
pm2 = apply(background.adt.mtx, 1, max)
head(pm)
```

prot pmax  
CD34\_TotalSeqB 4  
CD80\_TotalSeqB 60  
CD274\_TotalSeqB 75  
IgG2b\_control\_TotalSeqB 90

``` r
# remove the non staining protein 
cell.adt.raw  = cell.adt.raw[!rownames(cell.adt.raw) == 'CD34_TotalSeqB', ]
background.adt.mtx = background.adt.mtx[!rownames(background.adt.mtx) == 'CD34_TotalSeqB', ]
```

## Step 4 Normalize protein data with the DSBNormalizeProtein Function <a name="step4"></a>

We are now ready to use dsb to normalize and denoise the ADTs. For data
with isotype control proteins, set `denoise.counts = TRUE` and
`use.isotype.control = TRUE` and provide a vector containing names of
isotype control proteins (the rownames of the protein matrix that are
isotype controls). For data without isotype controls, see the vignette
section *Using dsb with data lacking isotype
controls*.

``` r
#normalize protein data for the cell containing droplets with the dsb method. 
dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot, 
  empty_drop_matrix = negative_mtx_rawprot, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = rownames(cells_mtx_rawprot)[29:31] 
  )
# note: normalization takes ~ 20 seconds
# system.time()
# user  system elapsed 
#  20.799   0.209  21.783 
```

The function returns a matrix of normalized protein values which can be
integrated with any single cell analysis software. We provide an example
with Seurat, Bioconductor and Scanpy below.

Advanced users may want to examine internal stats used in dsb, in that
case use `return.stats = TRUE`. If the range of values after
normalization is large, this is often due to a few low or high magnitude
outliers. The simplest way to address these outliers is to clip the
maximum values by setting `quantile.clipping = TRUE`. Finally a
different pseudocount can be used with `define.pseudocount = TRUE` and
`pseudocount.use`. Please see the vignettes on CRAN to learn more about
different parameters that can be used in dsb

``` r
# dsb with non-standard options 
cell.adt.dsb.2 = DSBNormalizeProtein(
  cell_protein_matrix = cell.adt.raw,
  empty_drop_matrix = background.adt.mtx,
  denoise.counts = TRUE,
  use.isotype.control = TRUE,
  isotype.control.name.vec = rownames(cell.adt.raw)[29:31],
  define.pseudocount = TRUE,
  pseudocount.use = 1
  quantile.clipping = TRUE,
  return.stats = TRUE, 
)
```

## Integrating dsb with Seurat <a name="seurat"></a>

Create a Seurat object. Make sure to add the dsb normalized matrix
`cell.adt.dsb` to the `data` slot, not the `counts` slot.

``` r
# Seurat workflow 
library(Seurat)

# integrating with Seurat
stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))

# create Seurat object note: min.cells is a gene filter, not a cell filter
s = Seurat::CreateSeuratObject(counts = cell.rna.raw, 
                               meta.data = cellmd,
                               assay = "RNA", 
                               min.cells = 20)

# add dsb normalized matrix "dsb_norm_prot" to the "CITE" data (not counts!) slot
s[["CITE"]] = Seurat::CreateAssayObject(data = cell.adt.dsb)
```

## Clustering cells based on dsb normalized protein using Seurat <a name="seurat_clustering"></a>

Here we will cluster the cells and annotate them based on dsb normalized
protein levels. This is similar to the workflows used in our paper
[Kotliarov *et al.* 2020](https://doi.org/10.1038/s41591-020-0769-8). We
first run spectral clustering using Seurat directly on the dsb
normalized protein values **without** reducing dimensionality of the
cells x protein matrix with PCA.

``` r
# define proteins to use in clustering (non-isptype controls)
prots = rownames(s@assays$CITE@data)[1:28]

# cluster and run umap 
s = Seurat::FindNeighbors(object = s, dims = NULL,assay = 'CITE', 
                          features = prots, k.param = 30, 
                          verbose = FALSE)

# direct graph clustering 
s = Seurat::FindClusters(object = s, resolution = 1, 
                         algorithm = 3, 
                         graph.name = 'CITE_snn', 
                         verbose = FALSE)
# umap (optional)
# s = Seurat::RunUMAP(object = s, assay = "CITE", features = prots,
#                     seed.use = 1990, min.dist = 0.2, n.neighbors = 30,
#                     verbose = FALSE)

# make results dataframe 
d = cbind(s@meta.data, 
          as.data.frame(t(s@assays$CITE@data))
          # s@reductions$umap@cell.embeddings)
          )
```

To see if we recovered the expected cell populations, it is often more
informative and interpretable to first look at the summarized (median or
mean) protein expression in each cluster–we recommend doing this before
trying to look at cells in 2-d visualization plots like umap.

## dsb derived cluster interpretation <a name="interpretation"></a>

dsb values are interpretable as the number of standard deviations of
each protein from the expected noise with additional correction for cell
to cell technical variations.

``` r
# calculate the median protein expression separately for each cluster 
adt_plot = d %>% 
  dplyr::group_by(CITE_snn_res.1) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("CITE_snn_res.1") 
# plot a heatmap of the average dsb normalized values for each cluster
pheatmap::pheatmap(t(adt_plot), 
                   color = viridis::viridis(25, option = "B"), 
                   fontsize_row = 8, border_color = NA)
```

<img src="man/figures/prot_heatmap.png" width="400" height="400" />

Now we can cell type based on median dsb normalized protein level per
cluster.

``` r
clusters = c(0:13)
celltype = c("CD4_Tcell_Memory", # 0 
             "CD14_Monocytes", #1
             "CD14_Monocytes_Activated", #2
             "CD4_Naive_Tcell", #3
             "B_Cells", #4
             "NK_Cells", #5
             "CD4_Naive_Tcell_CD62L+", #6
             "CD8_Memory_Tcell", #7
             "DC", #8
             "CD8_Naive_Tcell", #9
             "CD4_Effector", #10
             "CD16_Monocyte", #11
             "DOUBLETS", #12
             "DoubleNegative_Tcell" #13
)

s@meta.data$celltype = plyr::mapvalues(
  x = s@meta.data$CITE_snn_res.1, 
  from = clusters,  to = celltype
  )

# # optional -- dimensionality reduction plot 
# Seurat::DimPlot(s, reduction = 'umap', group.by = 'celltype',
#               label = TRUE, repel = TRUE, label.size = 2.5, pt.size = 0.1) + 
#   theme_bw() + NoLegend() + ggtitle('dsb normalized protein')
```

<img src="man/figures/umap.png" width="400" height="400" />

## Weighted Nearest Neighbor multimodal clustering using dsb normalized values with Seurat <a name="wnn"></a>

Below we demonstrate using Seurat’s weighted nearest neighbors
multimodal clustering method with dsb normalized values as input for the
protein assay. The performance of this algorithm is better on larger
datasets but we demonstrate here on this small dataset as an example.

Below we show a modified version of WNN directly using normalized
protein values as well as the default original implementation of WNN
which uses PCA on protein data.

For a dataset with a smaller number of proteins, we have found that just
using the dsb normalized cells x protein directly rather than
compressing the ADT data into principal components can improve the
resulting clusters and interpretation of the joint embedding. Datasets
generated with recently available pre-titrated panels consisting of more
than 100 or 200 proteins may benefit more from dimensionality reduction
with
PCA.

### Method 1 – use protein + RNA directly without compressing protein into principle components

``` r
## WNN with dsb values 
# use RNA pearson residuals as normalized values for RNA pca 
DefaultAssay(s) = "RNA"
s = NormalizeData(s, verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = 'vst', verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# set up dsb values to use in WNN analysis 
DefaultAssay(s) = "CITE"
#  use normalized protein values as a dimensionality reduction object.
VariableFeatures(s) = prots

# run true pca to initialize dr pca slot for WNN 
s = ScaleData(s, assay = 'CITE', verbose = FALSE)
s = RunPCA(s, reduction.name = 'pdsb', features = VariableFeatures(s), verbose = FALSE)

# make matrix of norm values to add as dr embeddings
pseudo = t(s@assays$CITE@data)[,1:29]
pseudo_colnames = paste('pseudo', 1:29, sep = "_")
colnames(pseudo) = pseudo_colnames
s@reductions$pdsb@cell.embeddings = pseudo

# run WNN directly using dsb normalized values. 
s = FindMultiModalNeighbors(
  object = s,
  reduction.list = list("pca", "pdsb"),
  weighted.nn.name = "dsb_wnn", 
  knn.graph.name = "dsb_knn",
  modality.weight.name = "dsb_weight",
  snn.graph.name = "dsb_snn",
  dims.list = list(1:30, 1:29), 
  verbose = FALSE
)

s = FindClusters(s, graph.name = "dsb_knn", 
                 algorithm = 3, resolution = 1.5,
                 random.seed = 1990,  verbose = FALSE)
```

Now we can visualize the multimodal results with a multimodal heatmap of
mRNA and protein within each cluster.

``` r
# create multimodal heatmap 
vf = VariableFeatures(s,assay = "RNA")

Idents(s) = "dsb_knn_res.1.5"
DefaultAssay(s)  = "RNA"
rnade = FindAllMarkers(s, features = vf, only.pos = TRUE)
gene_plot = rnade %>% filter(avg_log2FC > 1 ) %>%  group_by(cluster) %>% top_n(3) %$% gene %>% unique 

s@meta.data$celltype_subcluster = paste(s@meta.data$celltype, s@meta.data$dsb_knn_res.1.5)

d = cbind(s@meta.data, 
          # protein 
          as.data.frame(t(s@assays$CITE@data)), 
          # mRNA
          as.data.frame(t(as.matrix(s@assays$RNA@data[gene_plot, ]))),
          s@reductions$umap@cell.embeddings)

# combined data 
adt_plot = d %>% 
  dplyr::group_by(dsb_knn_res.1.5) %>% 
  dplyr::summarize_at(.vars = c(prots, gene_plot), .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("dsb_knn_res.1.5") 


# make a combined plot 
suppressMessages(library(ComplexHeatmap))
# protein heatmap 
prot_col = circlize::colorRamp2(breaks = seq(-10,30, by = 2), colors = viridis::viridis(n = 21, option = "B", end = 0.95))
p1 = Heatmap(t(adt_plot)[prots, ], name = "protein",col = prot_col, use_raster = T,
             row_names_gp = gpar(color = "black", fontsize = 5))

# mRNA heatmap 
mrna = t(adt_plot)[gene_plot, ]
rna_col = circlize::colorRamp2(breaks = c(-2,-1,0,1,2), colors = colorspace::diverge_hsv(n = 5))
p2 = Heatmap(t(scale(t(mrna))), name = "mRNA", col = rna_col, use_raster = T, 
             clustering_method_columns = 'average',
             column_names_gp = gpar(color = "black", fontsize = 7), 
             row_names_gp = gpar(color = "black", fontsize = 5))

ht_list = p1 %v% p2
draw(ht_list)
```

<img src="man/figures/multimodal_heatmap.png" />

### Method 2 – the default Seurat WNN using PCA on protein normalized with dsb.

``` r
# use pearson residuals as normalized values for pca 
DefaultAssay(s) = "RNA"
s = NormalizeData(s, verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = 'vst', verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# set up dsb values to use in WNN analysis (do not normalize with CLR, use dsb normalized values)
DefaultAssay(s) = "CITE"
VariableFeatures(s) = prots
s = s %>% ScaleData() %>% RunPCA(reduction.name = 'apca')

# run WNN 
s = FindMultiModalNeighbors(
  s, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), 
  modality.weight.name = "RNA.weight"
)

# cluster 
s <- RunUMAP(s, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
s <- FindClusters(s, graph.name = "wsnn", algorithm = 3, resolution = 1.5, verbose = FALSE, random.seed = 1990)

# we can next proceed as above.
```

This tutorial is a guide. It cam be modified to fit your needs. For
topics see these vignettes on CRAN:

### understanding dsb

a “under the hood” guide with every assumption and modeling step to get
more intuition for how the method works.

### Other topics and FAQ

Integrating dsb with Bioconductor  
Integrating dsb with python/Scanpy  
Using dsb with data lacking isotype controls  
Integrating dsb with sample multiplexing experiments  
outlier handling with quantile clipping  
returning internal stats used by dsb  
Frequently Asked Questions

### Some recent publications using dsb <a name="pubications"></a>

Publications from outside institutes  
[Jardine et al. *Nature*
(2021)](https://www.nature.com/articles/s41586-021-03929-x)  
[Mimitou et. al *Nature Biotechnology*
(2021)](https://www.nature.com/articles/s41587-021-00927-2)  
[COMBAT consortium et al. *Cell*
(2022)](https://www.cell.com/cell/pdf/S0092-8674\(22\)00070-8.pdf)

Publications from the Tsang lab  
[Kotliarov et al. *Nature Medicine*
(2020)](https://www.nature.com/articles/s41591-020-0769-8)  
[Liu et al. *Cell*
(2021)](https://www.cell.com/cell/pdf/S0092-8674\(21\)00168-9.pdf)

### using other alignment algorithms <a name="otheraligners"></a>

dsb was developed prior to 10X Genomics supporting CITE-seq or hashing
data and we routinely use other alignment pipelines.

[**CITE-seq-count
documentation**](https://hoohm.github.io/CITE-seq-Count/Running-the-script/)

To use dsb properly with CITE-seq-Count you need to align background.
One way to do this is to set the `-cells` argument to ~ 200000. That
will align the top 200000 barcodes in terms of ADT library size, making
sure you capture the background.

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
