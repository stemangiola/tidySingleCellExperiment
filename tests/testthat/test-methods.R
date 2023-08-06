df <- pbmc_small

test_that("show()", {
    txt <- capture.output(show(df))
    expect_lt(length(txt), 20)
    expect_equal(grep("SingleCellExperiment", txt), 1)
    i <- grep(str <- ".*Features=([0-9]+).*", txt)
    expect_equal(gsub(str, "\\1", txt[i]), paste(nrow(df)))
    i <- grep(str <- ".*Cells=([0-9]+).*", txt)
    expect_equal(gsub(str, "\\1", txt[i]), paste(ncol(df)))
})

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
    fd <- join_features(df, gs, shape="wide", assay="counts")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_null(fd$.feature)
    expect_identical(
        unname(t(as.matrix(as_tibble(fd)[, make.names(gs)]))),
        as.matrix(unname(counts(df)[gs, ])))
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
