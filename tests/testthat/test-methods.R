
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
  gs <- c(sample(rownames(df), 3), sample(rownames(altExps(df)[[1]]), 3))
  sce_counts_combined <- do.call("rbind", append(lapply(altExps(df), counts), values = list(as.matrix(counts(df))), after = 0))
  # long
  fd <- join_features(df, gs, shape="long")
  expect_s3_class(fd, "tbl_df")
  expect_setequal(unique(fd$.feature), gs)
  expect_true(all(table(fd$.feature) == ncol(df)))
  expect_identical(
    {
      fd <- select(fd,.cell, .feature, starts_with(".abundance"))
      fd <- select(fd,.cell, .feature, ends_with("-counts")|ends_with("_counts"))
      fd <- mutate(fd, counts = case_when(is.na(unlist(pick(3))) ~ unlist(pick(4)), .default = unlist(pick(3))))
      fd <- pivot_wider(
        data=fd,
        id_cols = ".feature",
        names_from = ".cell",
        values_from = "counts")
      fd <- as.data.frame(fd[order(fd$.feature),])
      rownames(fd) <- fd$.feature
      fd <- fd[, -1]
      fd <- select(fd, sort(colnames(df)))
      fd <- as.matrix(fd)
      fd
    },
    as.matrix(sce_counts_combined[sort(gs), sort(colnames(df))]))

  # wide
  fd <- join_features(df, gs, shape="wide", assays="counts")
  expect_s4_class(fd, "SingleCellExperiment")
  expect_null(fd$.feature)
  expect_identical(
    {
      fd <- as.data.frame(select(fd, cell = .cell, make.names(sort(gs))))
      fd <- fd[order(fd$cell), ]
      rownames(fd) <- fd$cell
      fd <- fd[, -1]
      t(as.matrix(fd))
    }, {
      sce_wide_mat <- as.matrix(sce_counts_combined[sort(gs), sort(colnames(df))])
      rownames(sce_wide_mat) <- make.names(rownames(sce_wide_mat))
      sce_wide_mat
    }
  )
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
  fd <- aggregate_cells(df, .sample = c(factor, string), assays = assayNames(df))
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
  # Aggregate when using multiple assays
  assays_to_use <- c("logcounts", "ADT-logcounts")
  fd <- aggregate_cells(df, .sample = c(factor, string), assays = assays_to_use)
  expect_identical(assayNames(fd), assays_to_use)
  fd_all_features <- tidySingleCellExperiment:::get_all_features(df) |>
    filter(assay_id %in% assays_to_use) |> pull(feature) |> sort()
  expect_identical(fd_all_features, sort(rownames(fd)))
})
