tidySingleCellExperiment - part of tidytranscriptomics
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidySingleCellExperiment/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/stemangiola/tidySingleCellExperiment/actions)
<!-- badges: end -->

**Brings SingleCellExperiment to the tidyverse!**

Website:
[tidySingleCellExperiment](https://stemangiola.github.io/tidySingleCellExperiment/articles/introduction.html)

Please also have a look at

- [tidySummarizedExperiment](https://stemangiola.github.io/tidySummarizedExperiment/)
  for tidy manipulation of SummarizedExperiment objects)
- [tidyseurat](https://stemangiola.github.io/tidyseurat/) for tidy
  manipulation of Seurat objects
- [tidybulk](https://stemangiola.github.io/tidybulk/) for tidy bulk
  RNA-seq data analysis
- [tidygate](https://github.com/stemangiola/tidygate) for adding custom
  gate information to your tibble
- [tidyHeatmap](https://stemangiola.github.io/tidyHeatmap/) for heatmaps
  produced with tidy principles

# Introduction

tidySingleCellExperiment provides a bridge between Bioconductor
single-cell packages \[@amezquita2019orchestrating\] and the tidyverse
\[@wickham2019welcome\]. It enables viewing the Bioconductor
*SingleCellExperiment* object as a tidyverse tibble, and provides
SingleCellExperiment-compatible *dplyr*, *tidyr*, *ggplot* and *plotly*
functions. This allows users to get the best of both Bioconductor and
tidyverse worlds.

## Functions/utilities available

| SingleCellExperiment-compatible Functions | Description                                                                        |
|-------------------------------------------|------------------------------------------------------------------------------------|
| `all`                                     | After all `tidySingleCellExperiment` is a SingleCellExperiment object, just better |

| tidyverse Packages | Description                                        |
|--------------------|----------------------------------------------------|
| `dplyr`            | All `dplyr` tibble functions (e.g. `select`)       |
| `tidyr`            | All `tidyr` tibble functions (e.g. `pivot_longer`) |
| `ggplot2`          | `ggplot` (`ggplot`)                                |
| `plotly`           | `plot_ly` (`plot_ly`)                              |

| Utilities         | Description                                                      |
|-------------------|------------------------------------------------------------------|
| `as_tibble`       | Convert cell-wise information to a `tbl_df`                      |
| `join_features`   | Add feature-wise information, returns a `tbl_df`                 |
| `aggregate_cells` | Aggregate cell gene-transcription abundance as pseudobulk tissue |

## Installation

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("tidySingleCellExperiment")
```

Load libraries used in this vignette.

``` r
# Bioconductor single-cell packages
library(scater)
library(scran)
library(SingleR)
library(SingleCellSignalR)

# Tidyverse-compatible packages
library(purrr)
library(magrittr)
library(tidyHeatmap)

# Both
library(tidySingleCellExperiment)
```

# Data representation of `tidySingleCellExperiment`

This is a *SingleCellExperiment* object but it is evaluated as a tibble.
So it is compatible both with SingleCellExperiment and tidyverse.

``` r
pbmc_small_tidy <- tidySingleCellExperiment::pbmc_small 
```

**It looks like a tibble**

``` r
pbmc_small_tidy
```

    ## # A SingleCellExperiment-tibble abstraction: 80 × 17
    ## # [90mFeatures=230 | Cells=80 | Assays=counts, logcounts[0m
    ##    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
    ##    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
    ##  1 ATGC… SeuratPro…         70           47 0               A             g2    
    ##  2 CATG… SeuratPro…         85           52 0               A             g1    
    ##  3 GAAC… SeuratPro…         87           50 1               B             g2    
    ##  4 TGAC… SeuratPro…        127           56 0               A             g2    
    ##  5 AGTC… SeuratPro…        173           53 0               A             g2    
    ##  6 TCTG… SeuratPro…         70           48 0               A             g1    
    ##  7 TGGT… SeuratPro…         64           36 0               A             g1    
    ##  8 GCAG… SeuratPro…         72           45 0               A             g1    
    ##  9 GATA… SeuratPro…         52           36 0               A             g1    
    ## 10 AATG… SeuratPro…        100           41 0               A             g1    
    ## # ℹ 70 more rows
    ## # ℹ 10 more variables: RNA_snn_res.1 <fct>, file <chr>, ident <fct>,
    ## #   PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>, PC_4 <dbl>, PC_5 <dbl>, tSNE_1 <dbl>,
    ## #   tSNE_2 <dbl>

**But it is a SingleCellExperiment object after all**

``` r
assay(pbmc_small_tidy, "counts")[1:5, 1:5]
```

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##         ATGCCAGAACGACT CATGGCCTGTGCAT GAACCTGATGAACC TGACTGGATTCTCA
    ## MS4A1                .              .              .              .
    ## CD79B                1              .              .              .
    ## CD79A                .              .              .              .
    ## HLA-DRA              .              1              .              .
    ## TCL1A                .              .              .              .
    ##         AGTCAGACTGCACA
    ## MS4A1                .
    ## CD79B                .
    ## CD79A                .
    ## HLA-DRA              1
    ## TCL1A                .

# Annotation polishing

We may have a column that contains the directory each run was taken
from, such as the “file” column in `pbmc_small_tidy`.

``` r
pbmc_small_tidy$file[1:5]
```

    ## [1] "../data/sample2/outs/filtered_feature_bc_matrix/"
    ## [2] "../data/sample1/outs/filtered_feature_bc_matrix/"
    ## [3] "../data/sample2/outs/filtered_feature_bc_matrix/"
    ## [4] "../data/sample2/outs/filtered_feature_bc_matrix/"
    ## [5] "../data/sample2/outs/filtered_feature_bc_matrix/"

We may want to extract the run/sample name out of it into a separate
column. Tidyverse `extract` can be used to convert a character column
into multiple columns using regular expression groups.

``` r
# Create sample column
pbmc_small_polished <-
    pbmc_small_tidy |>
    extract(file, "sample", "../data/([a-z0-9]+)/outs.+", remove=FALSE)

# Reorder to have sample column up front
pbmc_small_polished |>
    select(sample, everything())
```

    ## # A SingleCellExperiment-tibble abstraction: 80 × 18
    ## # [90mFeatures=230 | Cells=80 | Assays=counts, logcounts[0m
    ##    .cell sample orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents
    ##    <chr> <chr>  <fct>           <dbl>        <int> <fct>           <fct>        
    ##  1 ATGC… sampl… SeuratPro…         70           47 0               A            
    ##  2 CATG… sampl… SeuratPro…         85           52 0               A            
    ##  3 GAAC… sampl… SeuratPro…         87           50 1               B            
    ##  4 TGAC… sampl… SeuratPro…        127           56 0               A            
    ##  5 AGTC… sampl… SeuratPro…        173           53 0               A            
    ##  6 TCTG… sampl… SeuratPro…         70           48 0               A            
    ##  7 TGGT… sampl… SeuratPro…         64           36 0               A            
    ##  8 GCAG… sampl… SeuratPro…         72           45 0               A            
    ##  9 GATA… sampl… SeuratPro…         52           36 0               A            
    ## 10 AATG… sampl… SeuratPro…        100           41 0               A            
    ## # ℹ 70 more rows
    ## # ℹ 11 more variables: groups <chr>, RNA_snn_res.1 <fct>, file <chr>,
    ## #   ident <fct>, PC_1 <dbl>, PC_2 <dbl>, PC_3 <dbl>, PC_4 <dbl>, PC_5 <dbl>,
    ## #   tSNE_1 <dbl>, tSNE_2 <dbl>

# Preliminary plots

Set colours and theme for plots.

``` r
# Use colourblind-friendly colours
friendly_cols <- dittoSeq::dittoColors()

# Set theme
custom_theme <-
    list(
        scale_fill_manual(values=friendly_cols),
        scale_color_manual(values=friendly_cols),
        theme_bw() +
            theme(
                panel.border=element_blank(),
                axis.line=element_line(),
                panel.grid.major=element_line(size=0.2),
                panel.grid.minor=element_line(size=0.1),
                text=element_text(size=12),
                legend.position="bottom",
                aspect.ratio=1,
                strip.background=element_blank(),
                axis.title.x=element_text(margin=margin(t=10, r=10, b=10, l=10)),
                axis.title.y=element_text(margin=margin(t=10, r=10, b=10, l=10))
            )
    )
```

    ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

We can treat `pbmc_small_polished` as a tibble for plotting.

Here we plot number of features per cell.

``` r
pbmc_small_polished |>
    ggplot(aes(nFeature_RNA, fill=groups)) +
    geom_histogram() +
    custom_theme
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](man/figures/plot1-1.png)<!-- -->

Here we plot total features per cell.

``` r
pbmc_small_polished |>
    ggplot(aes(groups, nCount_RNA, fill=groups)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.1) +
    custom_theme
```

![](man/figures/plot2-1.png)<!-- -->

Here we plot abundance of two features for each group.

``` r
pbmc_small_polished |>
    join_features(features=c("HLA-DRA", "LYZ")) |>
    ggplot(aes(groups, .abundance_counts + 1, fill=groups)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(aes(size=nCount_RNA), alpha=0.5, width=0.2) +
    scale_y_log10() +
    custom_theme
```

    ## tidySingleCellExperiment says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.

![](man/figures/unnamed-chunk-10-1.png)<!-- -->

# Preprocess the dataset

We can also treat `pbmc_small_polished` as a *SingleCellExperiment*
object and proceed with data processing with Bioconductor packages, such
as *scran* \[@lun2016pooling\] and *scater* \[@mccarthy2017scater\].

``` r
# Identify variable genes with scran
variable_genes <-
    pbmc_small_polished |>
    modelGeneVar() |>
    getTopHVGs(prop=0.1)

# Perform PCA with scater
pbmc_small_pca <-
    pbmc_small_polished |>
    runPCA(subset_row=variable_genes)
```

    ## Warning in check_numbers(k = k, nu = nu, nv = nv, limit = min(dim(x)) - : more
    ## singular values/vectors requested than available

    ## Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth =
    ## TRUE, : You're computing too large a percentage of total singular values, use a
    ## standard svd instead.

``` r
pbmc_small_pca
```

    ## # A SingleCellExperiment-tibble abstraction: 80 × 18
    ## # [90mFeatures=230 | Cells=80 | Assays=counts, logcounts[0m
    ##    .cell orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents groups
    ##    <chr> <fct>           <dbl>        <int> <fct>           <fct>         <chr> 
    ##  1 ATGC… SeuratPro…         70           47 0               A             g2    
    ##  2 CATG… SeuratPro…         85           52 0               A             g1    
    ##  3 GAAC… SeuratPro…         87           50 1               B             g2    
    ##  4 TGAC… SeuratPro…        127           56 0               A             g2    
    ##  5 AGTC… SeuratPro…        173           53 0               A             g2    
    ##  6 TCTG… SeuratPro…         70           48 0               A             g1    
    ##  7 TGGT… SeuratPro…         64           36 0               A             g1    
    ##  8 GCAG… SeuratPro…         72           45 0               A             g1    
    ##  9 GATA… SeuratPro…         52           36 0               A             g1    
    ## 10 AATG… SeuratPro…        100           41 0               A             g1    
    ## # ℹ 70 more rows
    ## # ℹ 11 more variables: RNA_snn_res.1 <fct>, file <chr>, sample <chr>,
    ## #   ident <fct>, PC1 <dbl>, PC2 <dbl>, PC3 <dbl>, PC4 <dbl>, PC5 <dbl>,
    ## #   tSNE_1 <dbl>, tSNE_2 <dbl>

If a tidyverse-compatible package is not included in the
tidySingleCellExperiment collection, we can use `as_tibble` to
permanently convert `tidySingleCellExperiment` into a tibble.

``` r
# Create pairs plot with GGally
pbmc_small_pca |>
    as_tibble() |>
    select(contains("PC"), everything()) |>
    GGally::ggpairs(columns=1:5, ggplot2::aes(colour=groups)) +
    custom_theme
```

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

![](man/figures/pc_plot-1.png)<!-- -->

# Identify clusters

We can proceed with cluster identification with *scran*.

``` r
pbmc_small_cluster <- pbmc_small_pca

# Assign clusters to the 'colLabels' of the SingleCellExperiment object
colLabels(pbmc_small_cluster) <-
    pbmc_small_pca |>
    buildSNNGraph(use.dimred="PCA") |>
    igraph::cluster_walktrap() %$%
    membership |>
    as.factor()
```

    ## Warning in (function (to_check, X, clust_centers, clust_info, dtype, nn, :
    ## detected tied distances to neighbors, see ?'BiocNeighbors-ties'

``` r
# Reorder columns
pbmc_small_cluster |> select(label, everything())
```

    ## # A SingleCellExperiment-tibble abstraction: 80 × 19
    ## # [90mFeatures=230 | Cells=80 | Assays=counts, logcounts[0m
    ##    .cell  label orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8 letter.idents
    ##    <chr>  <fct> <fct>           <dbl>        <int> <fct>           <fct>        
    ##  1 ATGCC… 2     SeuratPro…         70           47 0               A            
    ##  2 CATGG… 2     SeuratPro…         85           52 0               A            
    ##  3 GAACC… 2     SeuratPro…         87           50 1               B            
    ##  4 TGACT… 1     SeuratPro…        127           56 0               A            
    ##  5 AGTCA… 2     SeuratPro…        173           53 0               A            
    ##  6 TCTGA… 2     SeuratPro…         70           48 0               A            
    ##  7 TGGTA… 1     SeuratPro…         64           36 0               A            
    ##  8 GCAGC… 2     SeuratPro…         72           45 0               A            
    ##  9 GATAT… 2     SeuratPro…         52           36 0               A            
    ## 10 AATGT… 2     SeuratPro…        100           41 0               A            
    ## # ℹ 70 more rows
    ## # ℹ 12 more variables: groups <chr>, RNA_snn_res.1 <fct>, file <chr>,
    ## #   sample <chr>, ident <fct>, PC1 <dbl>, PC2 <dbl>, PC3 <dbl>, PC4 <dbl>,
    ## #   PC5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>

And interrogate the output as if it was a regular tibble.

``` r
# Count number of cells for each cluster per group
pbmc_small_cluster |>
    count(groups, label)
```

    ## tidySingleCellExperiment says: A data frame is returned for independent data analysis.

    ## # A tibble: 8 × 3
    ##   groups label     n
    ##   <chr>  <fct> <int>
    ## 1 g1     1        12
    ## 2 g1     2        14
    ## 3 g1     3        14
    ## 4 g1     4         4
    ## 5 g2     1        10
    ## 6 g2     2        11
    ## 7 g2     3        10
    ## 8 g2     4         5

We can identify and visualise cluster markers combining
SingleCellExperiment, tidyverse functions and tidyHeatmap
\[@mangiola2020tidyheatmap\]

``` r
# Identify top 10 markers per cluster
marker_genes <-
    pbmc_small_cluster |>
    findMarkers(groups=pbmc_small_cluster$label) |>
    as.list() |>
    map(~ .x |>
        head(10) |>
        rownames()) |>
    unlist()

# Plot heatmap
pbmc_small_cluster |>
    join_features(features=marker_genes) |>
    group_by(label) |>
    heatmap(.feature, .cell, .abundance_counts, .scale="column")
```

    ## tidySingleCellExperiment says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis.

    ## tidyHeatmap says: (once per session) from release 1.7.0 the scaling is set to "none" by default. Please use scale = "row", "column" or "both" to apply scaling

    ## Warning: The `.scale` argument of `heatmap()` is deprecated as of tidyHeatmap 1.7.0.
    ## ℹ Please use scale (without dot prefix) instead: heatmap(scale = ...)
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](man/figures/unnamed-chunk-11-1.png)<!-- -->

# Reduce dimensions

We can calculate the first 3 UMAP dimensions using the
SingleCellExperiment framework and *scater*.

``` r
pbmc_small_UMAP <-
    pbmc_small_cluster |>
    runUMAP(ncomponents=3)
```

And we can plot the result in 3D using plotly.

``` r
pbmc_small_UMAP |>
    plot_ly(
        x=~`UMAP1`,
        y=~`UMAP2`,
        z=~`UMAP3`,
        color=~label,
        colors=friendly_cols[1:4]
    )
```

<figure>
<img src="man/figures/plotly.png" alt="plotly screenshot" />
<figcaption aria-hidden="true">plotly screenshot</figcaption>
</figure>

# Cell type prediction

We can infer cell type identities using *SingleR* \[@aran2019reference\]
and manipulate the output using tidyverse.

``` r
# Get cell type reference data
blueprint <- celldex::BlueprintEncodeData()

# Infer cell identities
cell_type_df <-

    assays(pbmc_small_UMAP)$logcounts |>
    Matrix::Matrix(sparse = TRUE) |>
    SingleR::SingleR(
        ref = blueprint,
        labels = blueprint$label.main,
        method = "single"
    ) |>
    as.data.frame() |>
    as_tibble(rownames="cell") |>
    select(cell, first.labels)
```

``` r
# Join UMAP and cell type info
pbmc_small_cell_type <-
    pbmc_small_UMAP |>
    left_join(cell_type_df, by="cell")
```

    ## Warning in is_sample_feature_deprecated_used(x, when(by, !is.null(.) ~ by, :
    ## tidySingleCellExperiment says: from version 1.3.1, the special columns
    ## including cell id (colnames(se)) has changed to ".cell". This dataset is
    ## returned with the old-style vocabulary (cell), however we suggest to update
    ## your workflow to reflect the new vocabulary (.cell)

``` r
# Reorder columns
pbmc_small_cell_type |>
    select(cell, first.labels, everything())
```

    ## Warning in is_sample_feature_deprecated_used(.data, (enquos(..., .ignore_empty
    ## = "all") %>% : tidySingleCellExperiment says: from version 1.3.1, the special
    ## columns including cell id (colnames(se)) has changed to ".cell". This dataset
    ## is returned with the old-style vocabulary (cell), however we suggest to update
    ## your workflow to reflect the new vocabulary (.cell)

    ## # A SingleCellExperiment-tibble abstraction: 80 × 23
    ## # [90mFeatures=230 | Cells=80 | Assays=counts, logcounts[0m
    ##    cell          first.labels orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
    ##    <chr>         <chr>        <fct>           <dbl>        <int> <fct>          
    ##  1 ATGCCAGAACGA… CD4+ T-cells SeuratPro…         70           47 0              
    ##  2 CATGGCCTGTGC… CD8+ T-cells SeuratPro…         85           52 0              
    ##  3 GAACCTGATGAA… CD8+ T-cells SeuratPro…         87           50 1              
    ##  4 TGACTGGATTCT… CD4+ T-cells SeuratPro…        127           56 0              
    ##  5 AGTCAGACTGCA… CD4+ T-cells SeuratPro…        173           53 0              
    ##  6 TCTGATACACGT… CD4+ T-cells SeuratPro…         70           48 0              
    ##  7 TGGTATCTAAAC… CD4+ T-cells SeuratPro…         64           36 0              
    ##  8 GCAGCTCTGTTT… CD4+ T-cells SeuratPro…         72           45 0              
    ##  9 GATATAACACGC… CD4+ T-cells SeuratPro…         52           36 0              
    ## 10 AATGTTGACAGT… CD4+ T-cells SeuratPro…        100           41 0              
    ## # ℹ 70 more rows
    ## # ℹ 17 more variables: letter.idents <fct>, groups <chr>, RNA_snn_res.1 <fct>,
    ## #   file <chr>, sample <chr>, ident <fct>, label <fct>, PC1 <dbl>, PC2 <dbl>,
    ## #   PC3 <dbl>, PC4 <dbl>, PC5 <dbl>, tSNE_1 <dbl>, tSNE_2 <dbl>, UMAP1 <dbl>,
    ## #   UMAP2 <dbl>, UMAP3 <dbl>

We can easily summarise the results. For example, we can see how cell
type classification overlaps with cluster classification.

``` r
# Count number of cells for each cell type per cluster
pbmc_small_cell_type |>
    count(label, first.labels)
```

    ## tidySingleCellExperiment says: A data frame is returned for independent data analysis.

    ## # A tibble: 11 × 3
    ##    label first.labels     n
    ##    <fct> <chr>        <int>
    ##  1 1     CD4+ T-cells     2
    ##  2 1     CD8+ T-cells     8
    ##  3 1     NK cells        12
    ##  4 2     B-cells         10
    ##  5 2     CD4+ T-cells     6
    ##  6 2     CD8+ T-cells     2
    ##  7 2     Macrophages      1
    ##  8 2     Monocytes        6
    ##  9 3     Macrophages      1
    ## 10 3     Monocytes       23
    ## 11 4     Erythrocytes     9

We can easily reshape the data for building information-rich faceted
plots.

``` r
pbmc_small_cell_type |>

    # Reshape and add classifier column
    pivot_longer(
        cols=c(label, first.labels),
        names_to="classifier", values_to="label"
    ) |>

    # UMAP plots for cell type and cluster
    ggplot(aes(UMAP1, UMAP2, color=label)) +
    geom_point() +
    facet_wrap(~classifier) +
    custom_theme
```

    ## tidySingleCellExperiment says: A data frame is returned for independent data analysis.

![](man/figures/unnamed-chunk-15-1.png)<!-- -->

We can easily plot gene correlation per cell category, adding
multi-layer annotations.

``` r
pbmc_small_cell_type |>

    # Add some mitochondrial abundance values
    mutate(mitochondrial=rnorm(dplyr::n())) |>

    # Plot correlation
    join_features(features=c("CST3", "LYZ"), shape="wide") |>
    ggplot(aes(CST3 + 1, LYZ + 1, color=groups, size=mitochondrial)) +
    geom_point() +
    facet_wrap(~first.labels, scales="free") +
    scale_x_log10() +
    scale_y_log10() +
    custom_theme
```

    ## Warning in is_sample_feature_deprecated_used(x, when(by, !is.null(.) ~ by, :
    ## tidySingleCellExperiment says: from version 1.3.1, the special columns
    ## including cell id (colnames(se)) has changed to ".cell". This dataset is
    ## returned with the old-style vocabulary (cell), however we suggest to update
    ## your workflow to reflect the new vocabulary (.cell)

![](man/figures/unnamed-chunk-16-1.png)<!-- -->

# Nested analyses

A powerful tool we can use with tidySingleCellExperiment is tidyverse
`nest`. We can easily perform independent analyses on subsets of the
dataset. First we classify cell types into lymphoid and myeloid, and
then nest based on the new classification.

``` r
pbmc_small_nested <-
    pbmc_small_cell_type |>
    filter(first.labels != "Erythrocytes") |>
    mutate(cell_class=dplyr::if_else(`first.labels` %in% c("Macrophages", "Monocytes"), "myeloid", "lymphoid")) |>
    nest(data=-cell_class)
```

    ## Warning: There were 2 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `data = map(...)`.
    ## Caused by warning in `is_sample_feature_deprecated_used()`:
    ## ! tidySingleCellExperiment says: from version 1.3.1, the special columns including cell id (colnames(se)) has changed to ".cell". This dataset is returned with the old-style vocabulary (cell), however we suggest to update your workflow to reflect the new vocabulary (.cell)
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

``` r
pbmc_small_nested
```

    ## # A tibble: 2 × 2
    ##   cell_class data           
    ##   <chr>      <list>         
    ## 1 lymphoid   <SnglCllE[,40]>
    ## 2 myeloid    <SnglCllE[,31]>

Now we can independently for the lymphoid and myeloid subsets (i) find
variable features, (ii) reduce dimensions, and (iii) cluster using both
tidyverse and SingleCellExperiment seamlessly.

``` r
pbmc_small_nested_reanalysed <-
    pbmc_small_nested |>
    mutate(data=map(
        data, ~ {
            .x <- runPCA(.x, subset_row=variable_genes)

            variable_genes <-
                .x |>
                modelGeneVar() |>
                getTopHVGs(prop=0.3)

            colLabels(.x) <-
                .x |>
                buildSNNGraph(use.dimred="PCA") |>
                igraph::cluster_walktrap() %$%
                membership |>
                as.factor()

            .x |> runUMAP(ncomponents=3)
        }
    ))

pbmc_small_nested_reanalysed
```

    ## # A tibble: 2 × 2
    ##   cell_class data           
    ##   <chr>      <list>         
    ## 1 lymphoid   <SnglCllE[,40]>
    ## 2 myeloid    <SnglCllE[,31]>

We can then unnest and plot the new classification.

``` r
pbmc_small_nested_reanalysed |>

    # Convert to tibble otherwise SingleCellExperiment drops reduced dimensions when unifying data sets.
    mutate(data=map(data, ~ .x |> as_tibble())) |>
    unnest(data) |>

    # Define unique clusters
    unite("cluster", c(cell_class, label), remove=FALSE) |>

    # Plotting
    ggplot(aes(UMAP1, UMAP2, color=cluster)) +
    geom_point() +
    facet_wrap(~cell_class) +
    custom_theme
```

![](man/figures/unnamed-chunk-19-1.png)<!-- -->

We can perform a large number of functional analyses on data subsets.
For example, we can identify intra-sample cell-cell interactions using
*SingleCellSignalR* \[@cabello2020singlecellsignalr\], and then compare
whether interactions are stronger or weaker across conditions. The code
below demonstrates how this analysis could be performed. It won’t work
with this small example dataset as we have just two samples (one for
each condition). But some example output is shown below and you can
imagine how you can use tidyverse on the output to perform t-tests and
visualisation.

``` r
pbmc_small_nested_interactions <-
    pbmc_small_nested_reanalysed |>

    # Unnest based on cell category
    unnest(data) |>

    # Create unambiguous clusters
    mutate(integrated_clusters=first.labels |> as.factor() |> as.integer()) |>

    # Nest based on sample
    nest(data=-sample) |>
    mutate(interactions=map(data, ~ {

        # Produce variables. Yuck!
        cluster <- colData(.x)$integrated_clusters
        data <- data.frame(assays(.x) |> as.list() |> extract2(1) |> as.matrix())

        # Ligand/Receptor analysis using SingleCellSignalR
        data |>
            cell_signaling(genes=rownames(data), cluster=cluster) |>
            inter_network(data=data, signal=_, genes=rownames(data), cluster=cluster) %$%
            `individual-networks` |>
            map_dfr(~ bind_rows(as_tibble(.x)))
    }))

pbmc_small_nested_interactions |>
    select(-data) |>
    unnest(interactions)
```

If the dataset was not so small, and interactions could be identified,
you would see something like below.

``` r
tidySingleCellExperiment::pbmc_small_nested_interactions
```

    ## # A tibble: 100 × 9
    ##    sample  ligand          receptor ligand.name receptor.name origin destination
    ##    <chr>   <chr>           <chr>    <chr>       <chr>         <chr>  <chr>      
    ##  1 sample1 cluster 1.PTMA  cluster… PTMA        VIPR1         clust… cluster 2  
    ##  2 sample1 cluster 1.B2M   cluster… B2M         KLRD1         clust… cluster 2  
    ##  3 sample1 cluster 1.IL16  cluster… IL16        CD4           clust… cluster 2  
    ##  4 sample1 cluster 1.HLA-B cluster… HLA-B       KLRD1         clust… cluster 2  
    ##  5 sample1 cluster 1.CALM1 cluster… CALM1       VIPR1         clust… cluster 2  
    ##  6 sample1 cluster 1.HLA-E cluster… HLA-E       KLRD1         clust… cluster 2  
    ##  7 sample1 cluster 1.GNAS  cluster… GNAS        VIPR1         clust… cluster 2  
    ##  8 sample1 cluster 1.B2M   cluster… B2M         HFE           clust… cluster 2  
    ##  9 sample1 cluster 1.PTMA  cluster… PTMA        VIPR1         clust… cluster 3  
    ## 10 sample1 cluster 1.CALM1 cluster… CALM1       VIPR1         clust… cluster 3  
    ## # ℹ 90 more rows
    ## # ℹ 2 more variables: interaction.type <chr>, LRscore <dbl>

# Aggregating cells

Sometimes, it is necessary to aggregate the gene-transcript abundance
from a group of cells into a single value. For example, when comparing
groups of cells across different samples with fixed-effect models.

In tidySingleCellExperiment, cell aggregation can be achieved using the
`aggregate_cells` function.

``` r
pbmc_small_tidy |>
  aggregate_cells(groups, assays = "counts")
```

    ## class: SummarizedExperiment 
    ## dim: 230 2 
    ## metadata(0):
    ## assays(1): counts
    ## rownames(230): ACAP1 ACRBP ... ZNF330 ZNF76
    ## rowData names(1): feature
    ## colnames(2): g1 g2
    ## colData names(4): groups .aggregated_cells orig.ident file
