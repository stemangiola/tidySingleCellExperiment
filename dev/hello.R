library(scRNAseq)
sce <- GrunPancreasData()
sce


library(scran)
library(scater)
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats, percent_subsets="altexps_ERCC_percent")
sce <- sce[,!qcfilter$discard]
summary(qcfilter$discard)


clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
summary(sizeFactors(sce))

sce2 <- computeSpikeFactors(sce, "ERCC")
summary(sizeFactors(sce2))

sce <- logNormCounts(sce)

dec <- modelGeneVar(sce)
top.hvgs <- getTopHVGs(dec, prop=0.1)

sce <- runPCA(sce, subset_row=top.hvgs)
reducedDimNames(sce)

dec2 <- modelGeneVarWithSpikes(sce, 'ERCC')
sced <- denoisePCA(sce, dec2, subset.row=getTopHVGs(dec2, prop=0.1))
ncol(reducedDim(sced, "PCA"))
output <- getClusteredPCs(reducedDim(sce))
npcs <- metadata(output)$chosen
reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[,1:npcs,drop=FALSE]
npcs

# In this case, using the PCs that we chose from getClusteredPCs().
g <- buildSNNGraph(sce, use.dimred="PCAsub")
cluster <- igraph::cluster_walktrap(g)$membership

# Assigning to the 'colLabels' of the 'sce'.
colLabels(sce) <- factor(cluster)
table(colLabels(sce))

sce <- runTSNE(sce, dimred="PCAsub")

# PCA
sce@int_colData@listData$reducedDims

# Gene location
sce@rowRanges

# Cell annotation
sce@colData

# normalisations
sce@assays




# Combining:
rbind(sce, sce)
cbind(sce, sce)

# Subsetting:
sce[1:10,]
sce[,1:5]

sce2 <- sce
sce2[1:10,] <- sce[11:20,]

# Can also use subset()
sce$WHEE <- sample(LETTERS, ncol(sce), replace=TRUE)
subset(sce, , WHEE=="A")

# Can also use split()
split(sce, sample(LETTERS, nrow(sce), replace=TRUE))
