context("tidyr test")

tt <-
  pbmc_small %>%
  tidy() %>%
  mutate(col2 = "other_col")

test_that("nest_unnest", {
  col_names <- colnames(tt@colData) %>% c("cell")

  x <- tt %>%
    nest(data = -groups) %>%
    unnest(data) %>%
    scater::logNormCounts() %>%
    scater::runPCA()
  y <- tt %>%
    scater::logNormCounts() %>%
    scater::runPCA()


  expect_equal(
    x@int_colData@listData$reducedDims$PCA %>% as.data.frame() %>% as_tibble(rownames = "cell") %>% arrange(cell) %>% pull(PC1) %>% abs(),
    y@int_colData@listData$reducedDims$PCA %>% as.data.frame() %>% as_tibble(rownames = "cell") %>% arrange(cell) %>% pull(PC1) %>% abs()
  )
})

test_that("unite separate", {
  un <- tt %>% unite("new_col", c(groups, col2), sep = ":")

  expect_equal(un %>% select(new_col) %>% slice(1) %>% pull(new_col), "g2:other_col")

  se <- un %>% separate(col = new_col, into = c("orig.ident", "groups"), sep = ":")

  expect_equal(se %>% select(orig.ident) %>% ncol(), 1)
})

test_that("extract", {
  expect_equal(
    tt %>% extract(groups, into = "g", regex = "g([0-9])", convert = TRUE) %>% pull(g) %>% class(),
    "integer"
  )
})

test_that("pivot_longer", {
  expect_equal(
    tt %>% pivot_longer(c(orig.ident, groups), names_to = "name", values_to = "value") %>% class() %>% .[1],
    "tbl_df"
  )
})
