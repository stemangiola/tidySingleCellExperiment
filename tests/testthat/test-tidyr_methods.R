context("tidyr test")

tt <-   pbmc_small %>%   mutate(col2 = "other_col")

test_that("nest_unnest", {
    col_names <- tt %>% colData %>% colnames() %>% c("cell")

    x <- tt %>%
        nest(data = -groups) %>%
        unnest(data) %>%
        scater::logNormCounts() %>%
        scater::runPCA()
    y <- tt %>%
        scater::logNormCounts() %>%
        scater::runPCA()


    expect_equal(
        reducedDims(x)$PCA %>%
            as.data.frame() %>%
            as_tibble(rownames = "cell") %>%
            arrange(cell) %>%
            pull(PC1) %>%
            abs(),
        reducedDims(x)$PCA %>%
            as.data.frame() %>%
            as_tibble(rownames = "cell") %>%
            arrange(cell) %>%
            pull(PC1) %>%
            abs()
    )
})

test_that("unite separate", {
    un <- tt %>% unite("new_col", c(groups, col2), sep = ":")

    un %>%
        select(new_col) %>%
        slice(1) %>%
        pull(new_col) %>%
        expect_equal("g2:other_col")

    se <-
        un %>%
        separate(
            col = new_col,
            into = c("orig.ident", "groups"),
            sep = ":"
        )

    se %>%
        select(orig.ident) %>%
        ncol() %>%
        expect_equal(1)
})

test_that("extract", {
    tt %>%
        extract(groups,
                into = "g",
                regex = "g([0-9])",
                convert = TRUE) %>%
        pull(g) %>%
        class() %>%
        expect_equal("integer")
})

test_that("pivot_longer", {
    tt %>%
        pivot_longer(c(orig.ident, groups),
                     names_to = "name",
                     values_to = "value") %>%
        class() %>%
        .[1] %>%
        expect_equal("tbl_df")
})
