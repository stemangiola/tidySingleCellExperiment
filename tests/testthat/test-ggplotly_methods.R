library(MASS, include.only = "rnegbin")
data("pbmc_small")
# Mock up ADT and cell hashing experiments
set.seed(2023-08-29)
# Antibody tags
pos_myus <- sample(x = 7500:20000, size = 5)
int_myus <- sample(x = 100:1500, size = 5)
neg_myus <- sample(x = 1:30, size = 15)
all_myus <- c(pos_myus, int_myus, neg_myus)
all_myus <- sample(x = all_myus, size = length(all_myus))

mat <- list()
for(i in seq_along(all_myus)) {
  mat[[i]] <- rnegbin(n = dim(pbmc_small)[[2]], mu = all_myus[[i]], theta = all_myus[[i]]/500)
}
mat <- Reduce(f = cbind, x = mat)
colnames(mat) <- paste("Ab", seq_along(mat[1,]), sep = "-")
rownames(mat) <- colnames(pbmc_small)

altExps(pbmc_small)[["ADT"]] <- SingleCellExperiment(assays = list(counts = t(mat), logcounts = log10(t(mat) + 1)))

# Cell hashing
HTO_myus <- sample(x = c(100, 100000), size = 6, replace = TRUE)
mat <- list()
for(i in seq_along(HTO_myus)) {
  mat[[i]] <- rnegbin(n = dim(pbmc_small)[[2]], mu = HTO_myus[[i]], theta = HTO_myus[[i]]/500)
}

mat <- Reduce(f = cbind, x = mat)
colnames(mat) <- paste("HTO", seq_along(mat[1,]), sep = "-")
rownames(mat) <- colnames(pbmc_small)

altExps(pbmc_small)[["Hashtag demultiplex"]] <- SingleCellExperiment(assays = list(counts = t(mat), logcounts = log10(t(mat) + 1)))

df <- pbmc_small
df$number <- rnorm(ncol(df))
df$factor <- sample(gl(3, 1, ncol(df)))

test_that("ggplot()", {
    # cell metadata
    p <- ggplot(df, aes(factor, number))
    expect_silent(show(p))
    expect_s3_class(p, "ggplot")
    # assay data
    g <- sample(rownames(df), 1)
    fd <- join_features(df, g, shape="wide", assays = "counts")
    p <- ggplot(fd, aes(factor, .data[[make.names(g)]]))
    expect_silent(show(p))
    expect_s3_class(p, "ggplot")
    # reduced dimensions
    p <- ggplot(df, aes(PC_1, PC_2, col=factor))
    expect_silent(show(p))
    expect_s3_class(p, "ggplot")
})

test_that("plotly()", {
    # cell metadata
    p <- plot_ly(df, x=~factor, y=~number, type="violin")
    expect_silent(show(p))
    expect_s3_class(p, "plotly")
    # assay data
    g <- sample(rownames(df), 1)
    fd <- join_features(df, g, shape="wide", assays = "counts")
    p <- plot_ly(fd, x=~factor, y=g, type="violin")
    expect_silent(show(p))
    expect_s3_class(p, "plotly")
    # reduced dimensions
    p <- plot_ly(fd, x=~PC_1, y=~PC_2, type="scatter", mode="markers")
    expect_silent(show(p))
    expect_s3_class(p, "plotly")
})
