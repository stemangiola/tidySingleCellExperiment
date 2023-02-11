context("dplyr test")

library(magrittr)

test_that("arrange", {
    tt_pca_aranged <-
        pbmc_small %>%
        arrange(groups) %>%
        scater::logNormCounts() %>%
        scater::runPCA()

    tt_pca <-
        pbmc_small %>%
        scater::logNormCounts() %>%
        scater::runPCA()



    expect_equal(
        reducedDims(tt_pca_aranged)$PCA[sort(colnames(tt_pca_aranged)), 1:3] %>% abs() %>% head(),
        reducedDims(tt_pca_aranged)$PCA[sort(colnames(tt_pca_aranged)), 1:3] %>% abs() %>% head(),
        tollerance = 1e-3
    )
})

test_that("bind_rows", {
    tt_bind <- bind_rows(pbmc_small, pbmc_small)

    tt_bind %>%
        select(cell) %>%
        tidySingleCellExperiment:::to_tib() %>%
        dplyr::count(cell) %>%
        dplyr::count(n) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("bind_cols", {
    tt_bind <- pbmc_small %>% select(groups)

    pbmc_small %>%
        bind_cols(tt_bind) %>%
        select(groups...7) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("distinct", {
    pbmc_small %>%
        distinct(groups) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("filter", {
    pbmc_small %>%
        filter(groups == "g1") %>%
        ncol() %>%
        expect_equal(44)
})

test_that("group_by", {
    pbmc_small %>%
        group_by(groups) %>%
        nrow() %>%
        expect_equal(80)
})

test_that("summarise", {
    pbmc_small %>%
        summarise(mean(nCount_RNA)) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("mutate", {
    pbmc_small %>%
        mutate(groups = 1) %>%
        distinct(groups) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("rename", {
    pbmc_small %>%
        rename(s_score = groups) %>%
        select(s_score) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("left_join", {
    pbmc_small %>%
        left_join(pbmc_small %>%
                      distinct(groups) %>%
                      mutate(new_column = 1:2)) %>%
        colData() %>%
        ncol() %>%
        expect_equal(10)
})

test_that("inner_join", {
    pbmc_small %>%
        inner_join(pbmc_small %>%
                       distinct(groups) %>%
                       mutate(new_column = 1:2) %>%
                       slice(1)) %>%
        ncol() %>%
        expect_equal(36)
})

test_that("right_join", {
    pbmc_small %>%
        right_join(pbmc_small %>%
                       distinct(groups) %>%
                       mutate(new_column = 1:2) %>%
                       slice(1)) %>%
        ncol() %>%
        expect_equal(36)
})

test_that("full_join", {
    pbmc_small %>%
        full_join(tibble::tibble(groups = "g1", other = 1:4)) %>%
        nrow() %>%
        expect_equal(212)
})

test_that("slice", {
    pbmc_small %>%
        slice(1) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("select", {
    pbmc_small %>%
        select(cell, orig.ident) %>%
        class() %>%
        as.character() %>%
        expect_equal("SingleCellExperiment")

    pbmc_small %>%
        select(orig.ident) %>%
        class() %>%
        as.character() %>%
        .[1] %>%
        expect_equal("tbl_df")
})

test_that("sample_n", {
    pbmc_small %>%
        sample_n(50) %>%
        ncol() %>%
        expect_equal(50)

    expect_equal(   pbmc_small %>% sample_n(500, replace = TRUE) %>% ncol,   31  )
})

test_that("sample_frac", {
    pbmc_small %>%
        sample_frac(0.1) %>%
        ncol() %>%
        expect_equal(8)

    expect_equal(   pbmc_small %>% sample_frac(10, replace = TRUE) %>% ncol,   31  )
})

test_that("count", {
    pbmc_small %>%
        count(groups) %>%
        nrow() %>%
        expect_equal(2)
})

test_that("add count", {
  pbmc_small %>%
    add_count(groups) %>%
    nrow() %>%
    expect_equal(2)
})
