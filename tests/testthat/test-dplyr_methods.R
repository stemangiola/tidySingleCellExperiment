context("dplyr test")

library(magrittr)

tt <- pbmc_small %>% tidy()

test_that("arrange", {
    tt_pca_aranged <-
        tt %>%
        arrange(groups) %>%
        scater::logNormCounts() %>%
        scater::runPCA()

    tt_pca <-
        tt %>%
        scater::logNormCounts() %>%
        scater::runPCA()


    expect_equal(
        tt_pca_aranged@int_colData@listData$reducedDims$PCA[sort(colnames(tt_pca_aranged)), 1:3] %>% abs() %>% head(),
        tt_pca@int_colData@listData$reducedDims$PCA[sort(colnames(tt_pca_aranged)), 1:3] %>% abs() %>% head(),
        tollerance = 1e-3
    )
})

test_that("bind_rows", {
    tt_bind <- bind_rows(tt, tt)

    tt_bind %>%
        select(cell) %>%
        tidySCE:::to_tib() %>%
        dplyr::count(cell) %>%
        dplyr::count(n) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("bind_cols", {
    tt_bind <- tt %>% select(groups)

    tt %>%
        bind_cols(tt_bind) %>%
        select(groups...7) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("distinct", {
    tt %>%
        distinct(groups) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("filter", {
    tt %>%
        filter(groups == "g1") %>%
        ncol() %>%
        expect_equal(44)
})

test_that("group_by", {
    tt %>%
        group_by(groups) %>%
        nrow() %>%
        expect_equal(80)
})

test_that("summarise", {
    tt %>%
        summarise(mean(nCount_RNA)) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("mutate", {
    tt %>%
        mutate(groups = 1) %>%
        distinct(groups) %>%
        nrow() %>%
        expect_equal(1)
})

test_that("rename", {
    tt %>%
        rename(s_score = groups) %>%
        select(s_score) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("left_join", {
    tt %>%
        left_join(tt %>%
                      distinct(groups) %>%
                      mutate(new_column = 1:2)) %>%
        `@`(colData) %>%
        ncol() %>%
        expect_equal(10)
})

test_that("inner_join", {
    tt %>%
        inner_join(tt %>%
                       distinct(groups) %>%
                       mutate(new_column = 1:2) %>%
                       slice(1)) %>%
        ncol() %>%
        expect_equal(36)
})

test_that("right_join", {
    tt %>%
        right_join(tt %>%
                       distinct(groups) %>%
                       mutate(new_column = 1:2) %>%
                       slice(1)) %>%
        ncol() %>%
        expect_equal(36)
})

test_that("full_join", {
    tt %>%
        full_join(tibble::tibble(groups = "g1", other = 1:4)) %>%
        nrow() %>%
        expect_equal(212)
})

test_that("slice", {
    tt %>%
        slice(1) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("select", {
    tt %>%
        select(cell, orig.ident) %>%
        class() %>%
        as.character() %>%
        expect_equal("tidySCE")

    tt %>%
        select(orig.ident) %>%
        class() %>%
        as.character() %>%
        .[1] %>%
        expect_equal("tbl_df")
})

test_that("sample_n", {
    tt %>%
        sample_n(50) %>%
        ncol() %>%
        expect_equal(50)
})

test_that("sample_frac", {
    tt %>%
        sample_frac(0.1) %>%
        ncol() %>%
        expect_equal(8)
})

test_that("count", {
    tt %>%
        count(groups) %>%
        nrow() %>%
        expect_equal(2)
})
