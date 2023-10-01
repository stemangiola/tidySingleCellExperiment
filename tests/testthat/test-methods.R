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
  mat[[i]] <-   MASS::rnegbin(n = dim(pbmc_small)[[2]], mu = all_myus[[i]], theta = all_myus[[i]]/500)
}
mat <- Reduce(f = cbind, x = mat)
colnames(mat) <- paste("Ab", seq_along(mat[1,]), sep = "-")
rownames(mat) <- colnames(pbmc_small)

altExps(pbmc_small)[["ADT"]] <- SingleCellExperiment(assays = list(counts = t(mat), logcounts = log10(t(mat) + 1)))

# Cell hashing
HTO_myus <- sample(x = c(100, 100000), size = 6, replace = TRUE)
mat <- list()
for(i in seq_along(HTO_myus)) {
  mat[[i]] <-   MASS::rnegbin(n = dim(pbmc_small)[[2]], mu = HTO_myus[[i]], theta = HTO_myus[[i]]/500)
}

mat <- Reduce(f = cbind, x = mat)
colnames(mat) <- paste("HTO", seq_along(mat[1,]), sep = "-")
rownames(mat) <- colnames(pbmc_small)

altExps(pbmc_small)[["Hashtag demultiplex"]] <- SingleCellExperiment(assays = list(counts = t(mat), logcounts = log10(t(mat) + 1)))

df <- pbmc_small

test_that("show()", {
  txt <- capture.output(show(df))
  expect_lt(length(txt), 20)
  expect_equal(grep("SingleCellExperiment", txt), 1)
  i <- grep(str <- ".*Features=([0-9]+).*", txt)
  expect_equal(gsub(str, "\\1", txt[i]), paste(nrow(df)))
  i <- grep(str <- ".*Cells=([0-9]+).*", txt)
  expect_equal(gsub(str, "\\1", txt[i]), paste(ncol(df)))
  i <- grep(".*s=.*", txt)
  j <- grep(".cell*", txt) -1
  header_text <- paste(txt[i:j], collapse = "") |> 
    stringr::str_remove_all(pattern = "# ") |> 
    stringr::str_remove_all(pattern = "\033") |> 
    stringr::str_remove_all(pattern = "\\[90m ") |> 
    stringr::str_remove_all(pattern = "\\[90m") |> 
    stringr::str_remove_all(pattern = "\\[0m")
  x <- header_text |> 
    stringr::str_remove(pattern = ".+s=") |> 
    strsplit(split = ", ") |> 
    unlist()
  y <- assayNames(df)
  for (k in seq_along(altExps(pbmc_small))) {
    y <- append(x = y, paste(altExpNames(pbmc_small)[[k]], assayNames(altExps(pbmc_small)[[k]]), sep = "-"))
  }
  for(j in seq_along(y)) {
    expect_contains(x, y[j])
  }
}
)

test_that("join_features()", {
  gs <- sample(rownames(df), 3)
  # long
  fd <- join_features(df, gs, shape="long")
  expect_s3_class(fd, "tbl_df")
  expect_setequal(unique(fd$.feature), gs)
  expect_true(all(table(fd$.feature) == ncol(df)))
  expect_identical(
    matrix(fd$.abundance_counts, nrow=length(gs)),
    as.matrix(unname(counts(df)[fd$.feature[seq_along(gs)], ])))
  # wide
  fd <- join_features(df, gs, shape="wide", assays="counts")
  expect_s4_class(fd, "SingleCellExperiment")
  expect_null(fd$.feature)
  expect_identical(
    unname(t(as.matrix(as_tibble(fd)[, make.names(gs)]))),
    as.matrix(unname(counts(df)[gs, ])))
  
  # Add features from altExp if they exist
  if(!is.null(altExp(df))) {
    gs <- sample(rownames(altExp(df)), 3)
    # long
    fd <- join_features(df, gs, shape="long")
    expect_s3_class(fd, "tbl_df")
    expect_setequal(unique(fd$.feature), gs)
    expect_true(all(table(fd$.feature) == ncol(df)))
    expect_identical(
      matrix(fd |> select(starts_with(".abundance")) |> pull(1), nrow=length(gs)),
      as.matrix(unname(assays(altExp(df))[[1]][fd$.feature[seq_along(gs)], ])))
    # wide
    fd <- join_features(df, gs, shape="wide", assays="ADT-counts")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_null(fd$.feature)
    expect_identical(
      unname(t(as.matrix(as_tibble(fd)[, make.names(gs)]))),
      as.matrix(unname(assays(altExp(df))[[1]][gs, ])))
  }
})

test_that("as_tibble()", {
  fd <- as_tibble(df)
  expect_s3_class(fd, "tbl_df")
  expect_equal(nrow(fd), ncol(df))
  ncd <- ncol(colData(df))
  ndr <- vapply(reducedDims(df), ncol, integer(1))
  expect_equal(ncol(fd), sum(1, ncd, ndr))
  # duplicated PCs
  reducedDim(df, "PCB") <- reducedDim(df, "PCA")
  fd <- as_tibble(mutate(df, abc=1))
  expect_equal(ncol(fd), ncol(as_tibble(df))+1)
})

test_that("aggregate_cells()", {
  df$factor <- sample(gl(3, 1, ncol(df)))
  df$string <- sample(c("a", "b"), ncol(df), TRUE)
  tbl <- distinct(select(df, factor, string))
  fd <- aggregate_cells(df, c(factor, string))
  expect_identical(assayNames(fd), assayNames(df))
  # [HLC: aggregate_cells() currently
  # reorders features alphabetically]
  fd <- fd[rownames(df), ]
  expect_s4_class(fd, "SummarizedExperiment")
  expect_equal(dim(fd), c(nrow(df), nrow(tbl)))
  foo <- mapply(
    f=tbl$factor,
    s=tbl$string,
    \(f, s) {
      expect_identical(
        df |> 
          filter(factor == f, string == s) |>
          assay() |> rowSums() |> as.vector(),
        fd[, fd$factor == f & fd$string == s] |>
          assay() |> as.vector())
    })
  # specified 'assays' are subsetted
  expect_error(aggregate_cells(df, c(factor, string), assays="x"))
  fd <- aggregate_cells(df, c(factor, string), assays="counts")
  expect_identical(assayNames(fd), "counts")
})
