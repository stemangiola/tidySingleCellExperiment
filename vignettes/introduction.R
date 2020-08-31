## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set( fig.path = "man/figures/")

## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(knitr)
knitr::opts_chunk$set(cache = TRUE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)

library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(SingleCellExperiment)
library(SingleCellSignalR)
library(scater)
library(scran)
library(igraph)
library(tidyHeatmap)
library(tidySCE)


my_theme = 	
  list(
    scale_fill_brewer(palette="Set1"),
    scale_color_brewer(palette="Set1") ,
    theme_bw() +
  	theme(
  		panel.border = element_blank(),
  		axis.line = element_line(),
  		panel.grid.major = element_line(size = 0.2),
  		panel.grid.minor = element_line(size = 0.1),
  		text = element_text(size=12),
  		legend.position="bottom",
  		aspect.ratio=1,
  		strip.background = element_blank(),
  		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
  		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
  	)
  )


## ----eval=FALSE---------------------------------------------------------------
#  install.packages("tidySCE")

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("stemangiola/tidySCE")

## -----------------------------------------------------------------------------
pbmc_small_tidy = tidySCE::pbmc_small %>% tidy()

## -----------------------------------------------------------------------------
pbmc_small_tidy

## -----------------------------------------------------------------------------
pbmc_small_tidy@assays

## -----------------------------------------------------------------------------
pbmc_small_polished =
  pbmc_small_tidy %>%
  extract(file, "sample", "../data/([a-z0-9]+)/outs.+", remove = F) 

pbmc_small_polished %>%
  select(sample, everything())

## ----plot1, cache=TRUE--------------------------------------------------------

pbmc_small_polished %>%
  tidySCE::ggplot(aes(nFeature_RNA, fill=groups)) + 
  geom_histogram() +
  my_theme


## ----plot2, cache=TRUE--------------------------------------------------------

pbmc_small_polished %>%
  tidySCE::ggplot(aes(groups, nCount_RNA, fill=groups)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  my_theme



## -----------------------------------------------------------------------------
pbmc_small_polished %>% 
  join_transcripts(transcripts = c("HLA-DRA" ,     "LYZ" )) %>%
  ggplot(aes(groups, abundance_counts + 1, fill=groups)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(size=nCount_RNA), alpha=0.5, width = 0.2) + 
  scale_y_log10() + 
  my_theme


## ----preprocess, cache=TRUE---------------------------------------------------

variable_gens = 
  pbmc_small_polished %>%
  modelGeneVar() %>%
  getTopHVGs(prop=0.1)

pbmc_small_pca = 
  pbmc_small_polished %>% 
  runPCA(   subset_row = variable_gens  )

pbmc_small_pca

## ----pc_plot, cache=TRUE------------------------------------------------------

pbmc_small_pca %>%
  as_tibble %>%
  select(contains("PC"), everything()) %>%
	  GGally::ggpairs(columns = 1:5, ggplot2::aes(colour=groups))


## ----cluster, cache=TRUE------------------------------------------------------

pbmc_small_cluster = pbmc_small_pca


# Assigning to the 'colLabels' of the 'sce'.
colLabels(pbmc_small_cluster) =
  pbmc_small_pca %>%
  buildSNNGraph(use.dimred="PCA") %>%
  igraph::cluster_walktrap() %$%
  membership %>%
  as.factor

pbmc_small_cluster %>% select(label, everything())

## ----cluster count, cache=TRUE------------------------------------------------

pbmc_small_cluster %>%
  tidySCE::count(groups, label)


## -----------------------------------------------------------------------------
# Identify markers
marker_genes = 
  pbmc_small_cluster %>%
  findMarkers(groups=pbmc_small_cluster$label) %>% 
  as.list %>% 
  map(~ .x %>% head(10) %>% rownames) %>% 
  unlist

# Plot heatmap
pbmc_small_cluster %>%
  join_transcripts(transcripts = marker_genes) %>% 
  group_by(label) %>% 
  heatmap(transcript, cell, abundance_counts, .scale = "column")
  

## ----umap, cache=TRUE---------------------------------------------------------

pbmc_small_UMAP = 
  pbmc_small_cluster %>%
  runUMAP(ncomponents = 3)


## ----umap plot, eval=FALSE----------------------------------------------------
#  
#  pbmc_small_UMAP %>%
#  	plot_ly(
#  		x = ~`UMAP_1`,
#  		y = ~`UMAP_2`,
#  		z = ~`UMAP_3`,
#  		color = ~ label
#  	)
#  

## ----eval=FALSE---------------------------------------------------------------
#  blueprint = SingleR::BlueprintEncodeData()
#  
#  cell_type_df =
#    pbmc_small_UMAP@assays@data$logcounts %>%
#    Matrix::Matrix(sparse = TRUE) %>%
#     SingleR::SingleR(
#         ref = blueprint ,
#         labels = blueprint$label.main,
#         method = "single"
#     ) %>%
#    as.data.frame() %>%
#    as_tibble(rownames="cell") %>%
#    select(cell, first.labels)

## -----------------------------------------------------------------------------
pbmc_small_cell_type =
  pbmc_small_UMAP %>%
  left_join(cell_type_df, by="cell")

pbmc_small_cell_type %>%
  tidySCE::select(cell, first.labels, everything())

## -----------------------------------------------------------------------------
pbmc_small_cell_type %>%
  count(label, first.labels)

## -----------------------------------------------------------------------------
pbmc_small_cell_type %>%
  
  # Reshaping
  pivot_longer(
    cols=c(label, first.labels), 
    names_to = "classifier", values_to = "label"
  ) %>%
  
  # Plotting
  ggplot(aes(UMAP1, UMAP2, color=label)) +
  geom_point() +
  facet_wrap(~classifier) +
  my_theme


## -----------------------------------------------------------------------------
pbmc_small_cell_type %>% 
  
  # Add mitochondrial abundance
  mutate(mitochondrial = rnorm(n())) %>%
  
  # Plot correlation
  join_transcripts(transcripts = c("CST3" ,     "LYZ" ), shape = "wide") %>%
  ggplot(aes(CST3 +1, LYZ + 1, color=groups, size=mitochondrial)) +
  geom_point() + 
  facet_wrap(~first.labels, scales = "free") +
  scale_x_log10() +
  scale_y_log10() +
  my_theme

## -----------------------------------------------------------------------------
pbmc_small_nested = 
  pbmc_small_cell_type %>%
  filter(first.labels != "Erythrocytes") %>%
  mutate(cell_class = if_else(`first.labels` %in% c("Macrophages", "Monocytes"), "myeloid", "lmphoid")) %>%
  nest(data = -cell_class)

pbmc_small_nested

## -----------------------------------------------------------------------------
pbmc_small_nested_reanalysed = 
  pbmc_small_nested %>%
  mutate(data = map(
    data, ~ {
      
      .x = runPCA(.x, subset_row = variable_gens )
      
      variable_gens = 
        .x %>%
        modelGeneVar() %>%
        getTopHVGs(prop=0.3)
      
      colLabels(.x) =
        .x %>%
        buildSNNGraph(use.dimred="PCA") %>%
        igraph::cluster_walktrap() %$%
        membership %>%
        as.factor
      
      .x %>% runUMAP(ncomponents = 3)
      
    }
  )) 

pbmc_small_nested_reanalysed

## -----------------------------------------------------------------------------
pbmc_small_nested_reanalysed %>%
  
  # Convert to tibble otherwise SingleCellExperiment drops reduced dimensions when unifying data sets.
  mutate(data = map(data, ~ .x %>% as_tibble)) %>%
  unnest(data) %>%

  # Define unique clusters
  unite("cluster", c(cell_class, label), remove=FALSE) %>%
  
  # Plotting
  ggplot(aes(UMAP1, UMAP2, color=cluster)) +
  geom_point() +
  facet_wrap(~cell_class) +
  my_theme



## ---- eval = FALSE------------------------------------------------------------
#  pbmc_small_nested_interactions =
#    pbmc_small_nested_reanalysed %>%
#  
#    # Unnest based on cell category
#    unnest(data) %>%
#  
#    # Create unambiguous clusters
#    mutate(integrated_clusters = first.labels %>% as.factor %>% as.integer) %>%
#  
#    # Nest based on sample
#  	tidySCE::nest(data = -sample) %>%
#  	tidySCE::mutate(interactions = map(data, ~ {
#  		
#  		# Produce variables. Yuck!
#  		cluster = .x@colData$integrated_clusters
#  		data = data.frame(.x@assays@data %>% as.list %>% .[[1]] %>% as.matrix)
#  
#  		# Ligand/Receptor analysis using SingleCellSignalR
#  		data %>%
#  			cell_signaling(genes=rownames(data),cluster=cluster) %>%
#  			inter_network(data = data, signal = ., genes = rownames(data), cluster = cluster) %$%
#  			`individual-networks` %>%
#  			map_dfr(~bind_rows(as_tibble(.x)))
#  	}))
#  
#  pbmc_small_nested_interactions %>%
#    select(-data) %>%
#    unnest(interactions)

## ---- echo = FALSE------------------------------------------------------------
tidySCE::pbmc_small_nested_interactions

