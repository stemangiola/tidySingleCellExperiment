library(S4Vectors)
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
df$number <- sample(seq(ncol(df)))
df$factor <- sample(
    factor(1:3, labels=paste0("g", 1:3)),
    ncol(df), TRUE, c(0.1, 0.3, 0.6))

test_that("arrange()", {
    expect_identical(
        arrange(df, number), 
        df[, order(df$number)])
    suppressWarnings({
        fd <- df %>%
            scater::logNormCounts() %>% 
            scater::runPCA()
    })
    expect_identical(
        arrange(fd, PC1), 
        fd[, order(reducedDim(fd)[, 1])])
    fd <- df %>%
        mutate(foo=seq(ncol(df))) %>%
        arrange(foo) %>% select(-foo)
    expect_identical(fd, df)
})

test_that("bind_rows()", {
    # warn about duplicated cells names
    expect_warning(fd <- bind_rows(df, df))
    # cell names should be unique after binding
    expect_true(!any(duplicated(pull(fd, .cell))))
})

test_that("bind_cols()", {
    fd <- bind_cols(df, select(df, factor))
    i <- grep("^factor", names(colData(fd)))
    expect_length(i, 2)
    expect_identical(fd[[i[1]]], df$factor)
    expect_identical(fd[[i[2]]], df$factor)
    expect_identical(
        select(fd, -starts_with("factor")), 
        select(df, -factor))
})

test_that("distinct()", {
    fd <- distinct(df, factor)
    expect_equal(nrow(fd), nlevels(df$factor))
    expect_identical(fd[[1]], unique(df$factor))
})

test_that("filter()", {
    fd <- filter(df, factor %in% levels(df$factor))
    expect_identical(df, fd)
    fd <- filter(df, factor == "g1")
    expect_equal(ncol(fd), sum(df$factor == "g1"))
    # missing cell names
    fd <- df; colnames(fd) <- NULL
    expect_silent(filter(df, number == 1))
    expect_message(fd <- filter(fd, number < 10))
    expect_type(pull(fd, .cell), "character")
    expect_null(colnames(fd))
})

test_that("group_by()", {
    fd <- group_by(df, factor)
    expect_equal(n_groups(fd), nlevels(df$factor))
    expect_equal(group_size(fd), tabulate(df$factor))
})

test_that("summaris/ze()", {
    fd <- mutate(df, n=runif(ncol(df)))
    ne <- summarise(fd, a=mean(n))
    mo <- summarize(fd, b=mean(n))
    expect_identical(ne$a, mean(fd$n))
    expect_identical(ne$a, mo$b)
})

test_that("mutate()", {
    fd <- mutate(df, peter="pan")
    expect_true(all(fd$peter == "pan"))
    fd <- mutate(df, number=paste(number))
    expect_identical(fd$number, paste(df$number))
    
    # special columns are blocked
    df |>
      mutate(.cell=1) |>
      expect_error("you are trying to mutate a column that is view only")
    
    df |>
      mutate(PC_10=1) |>
      expect_error("you are trying to mutate a column that is view only")
})

test_that("rename()", {
    fd <- rename(df, num=number, fac=factor)
    expect_identical(fd$num, df$number)
    expect_identical(fd$fac, df$factor)
    
    df |> 
      rename(ne=mo) |> 
      expect_error("Column `mo` doesn't exist")
    
    # special columns are blocked
    # ...'to' cannot be special
    
    df |>
      rename(a=PC_1) |>
      expect_error("you are trying to rename a column that is view only")  
    
    df |> 
      rename(a=.cell) |> 
      expect_error("you are trying to rename a column that is view only")
    # ...'from' cannot be special
    
    df |> 
      rename(PC_1=number) |> 
      expect_error("These names are duplicated")
    
    df |> 
      rename(.cell=number) |> 
      expect_error("These names are duplicated")
})

test_that("left_join()", {
    y <- df |> 
        distinct(factor) |> 
        mutate(string=letters[seq(nlevels(df$factor))])
    fd <- left_join(df, y, by="factor")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(n <- ncol(colData(fd)), ncol(colData(df))+1)
    expect_identical(colData(fd)[-n], colData(df))
})

test_that("left_join(), with DataFrame y", {
    y <- df |> 
        distinct(factor) |> 
        mutate(string=letters[seq(nlevels(df$factor))]) |> 
        DataFrame()
    fd <- left_join(df, y, by="factor")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(n <- ncol(colData(fd)), ncol(colData(df))+1)
    expect_identical(colData(fd)[-n], colData(df))
})

test_that("inner_join()", {
    y <- df |> 
        distinct(factor) |> 
        mutate(string=letters[seq(nlevels(df$factor))]) |> 
        slice(1)
    fd <- inner_join(df, y, by="factor")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(n <- ncol(colData(fd)), ncol(colData(df))+1)
    expect_equal(ncol(fd), sum(df$factor == fd$factor[1]))
})

test_that("inner_join(), with DataFrame y", {
    y <- df |> 
        distinct(factor) |> 
        mutate(string=letters[seq(nlevels(df$factor))]) |> 
        slice(1) |> DataFrame()
    fd <- inner_join(df, y, by="factor")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(n <- ncol(colData(fd)), ncol(colData(df))+1)
    expect_equal(ncol(fd), sum(df$factor == fd$factor[1]))
})

test_that("right_join()", {
    y <- df |>
        distinct(factor) |>
        mutate(string=letters[seq(nlevels(df$factor))]) |>
        slice(1)
    fd <- right_join(df, y, by="factor")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(n <- ncol(colData(fd)), ncol(colData(df))+1)
    expect_equal(ncol(fd), sum(df$factor == fd$factor[1]))
})

test_that("right_join(), with DataFrame y", {
    y <- df |>
        distinct(factor) |>
        mutate(string=letters[seq(nlevels(df$factor))]) |>
        slice(1) |> DataFrame()
    fd <- right_join(df, y, by="factor")
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(n <- ncol(colData(fd)), ncol(colData(df))+1)
    expect_equal(ncol(fd), sum(df$factor == fd$factor[1]))
})

test_that("full_join()", {
    # w/ duplicated cell names
    y <- tibble(factor="g2", other=1:3)
    fd <- expect_message(full_join(df, y, by="factor", relationship="many-to-many"))
    expect_s3_class(fd, "tbl_df")
    expect_true(all(is.na(fd$other[fd$factor != "g2"])))
    expect_true(all(!is.na(fd$other[fd$factor == "g2"])))
    expect_equal(nrow(fd), ncol(df)+2*sum(df$factor == "g2"))
    # w/o duplicates
    y <- tibble(factor="g2", other=1)
    fd <- expect_silent(full_join(df, y, by="factor"))
    expect_s4_class(fd, "SingleCellExperiment")
    expect_identical(
        select(fd, -other), 
        mutate(df, factor=paste(factor)))
})

test_that("full_join(), with DataFrame y", {
    # w/ duplicated cell names
    y <- tibble(factor="g2", other=1:3) |> DataFrame()
    fd <- expect_message(full_join(df, y, by="factor", relationship="many-to-many"))
    expect_s3_class(fd, "tbl_df")
    expect_true(all(is.na(fd$other[fd$factor != "g2"])))
    expect_true(all(!is.na(fd$other[fd$factor == "g2"])))
    expect_equal(nrow(fd), ncol(df)+2*sum(df$factor == "g2"))
    # w/o duplicates
    y <- tibble(factor="g2", other=1) |> DataFrame()
    fd <- expect_silent(full_join(df, y, by="factor"))
    expect_s4_class(fd, "SingleCellExperiment")
    expect_identical(
        select(fd, -other), 
        mutate(df, factor=paste(factor)))
})

test_that("slice()", {
    expect_identical(slice(df), df[, 0])
    expect_identical(slice(df, ncol(df)+1), df[, 0])
    expect_identical(slice(df, 1), df[, 1])
    expect_identical(slice(df, -1), df[, -1])
    i <- sample(ncol(df), 5)
    expect_identical(slice(df, i), df[, i])
    expect_identical(slice(df, -i), df[, -i])
})

test_that("slice_sample()", {
    pbmc_small |>
        slice_sample(n=0) |>
        ncol() |>
        expect_equal(0)
    pbmc_small |>
        slice_sample(n=50) |>
        ncol() |>
        expect_equal(50)
})

test_that("slice_head()", {
    pbmc_small |>
        slice_head(n=0) |>
        ncol() |>
        expect_equal(0)
    pbmc_small |>
        slice_head(n=50) |>
        ncol() |>
        expect_equal(50)
    expect_equal(
        colnames(pbmc_small) |> head(n=50),
        pbmc_small |> slice_head(n=50) |> colnames()
    )
})

test_that("slice_tail()", {
    pbmc_small |>
        slice_tail(n=0) |>
        ncol() |>
        expect_equal(0)
    pbmc_small |>
        slice_tail(n=50) |>
        ncol() |>
        expect_equal(50)
    expect_equal(
        colnames(pbmc_small) |> tail(n=50),
        pbmc_small |> slice_tail(n=50) |> colnames()
    )
})

test_that("slice_min()", {
    pbmc_small |>
        slice_min(nFeature_RNA, n=0) |>
        ncol() |>
        expect_equal(0)
    pbmc_small |>
        slice_min(nFeature_RNA, n=5) |>
        ncol() |>
        expect_equal(5)
    expect_equal(
        pbmc_small |> as_tibble() |>
            arrange(nFeature_RNA) |>
            head(n=5) %>% pull(.cell),
        pbmc_small |> slice_min(nFeature_RNA, n=5) |> colnames()
  )
})

test_that("slice_max()", {
    pbmc_small |>
        slice_max(nFeature_RNA, n=0) |>
        ncol() |>
        expect_equal(0)
    pbmc_small |>
        slice_max(nFeature_RNA, n = 5) |>
        ncol() |>
        expect_equal(5)
    expect_equal(
        pbmc_small |> as_tibble() |>
            arrange(desc(nFeature_RNA)) |>
            head(n=5) %>% pull(.cell),
        pbmc_small |> slice_max(nFeature_RNA, n=5) |> colnames()
  )
})

test_that("slice_min() slice_max() tibble input for order_by", {
  pbmc_small |>
    slice_min(tibble::tibble(nFeature_RNA, nCount_RNA), n=5) |>
    ncol() |>
    expect_equal(5)
  pbmc_small |>
    slice_max(tibble::tibble(nFeature_RNA, nCount_RNA), n=5) |>
    ncol() |>
    expect_equal(5)
})

test_that("select()", {
    fd <- select(df, .cell, number)
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(dim(fd), dim(df))
    fd <- select(df, number)
    expect_s3_class(fd, "tbl_df")
    expect_equal(nrow(fd), ncol(df))
})

test_that("sample_n()", {
    fd <- sample_n(df, n <- 50)
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(nrow(fd), nrow(df))
    expect_equal(ncol(fd), n)
    fd <- sample_n(df, 1e3, TRUE)
    expect_s3_class(fd, "tbl_df")
    expect_equal(nrow(fd), 1e3)
})

test_that("sample_frac()", {
    fd <- sample_frac(df, 0.1)
    expect_s4_class(fd, "SingleCellExperiment")
    expect_equal(nrow(fd), nrow(df))
    expect_equal(ncol(fd), ncol(df)/10)
    fd <- sample_frac(df, 10, TRUE)
    expect_s3_class(fd, "tbl_df")
    expect_equal(nrow(fd), ncol(df)*10)
})

test_that("count()", {
    fd <- count(df, factor)
    expect_s3_class(fd, "tbl_df")
    expect_equal(nrow(fd), nlevels(df$factor))
    expect_identical(fd$n, tabulate(df$factor))
})

test_that("add_count()", {
    fd <- add_count(df, factor)
    expect_identical(select(fd, -n), df)
    expect_identical(fd$n, unname(c(table(df$factor)[df$factor])))
})

test_that("rowwise()", {
    df |> 
    summarise(sum(lys)) |>
    expect_error("object 'lys' not found")
  
    df$lys <- replicate(ncol(df), sample(10, 3), FALSE)
    fd <- df |> rowwise() |> summarise(sum(lys))
    expect_s3_class(fd, "tbl_df")
    expect_equal(dim(fd), c(ncol(df), 1))
    expect_identical(fd[[1]], sapply(df$lys, sum))
})
