context('methods test')

data("pbmc_small")
library(dplyr)
test_that("join_features",{


  pbmc_small %>%
    join_features("CD3D") %>%
    slice(1) %>%
    tidySingleCellExperiment::pull(.abundance_counts) %>%
    expect_equal(4, tolerance=0.1)


})


test_that("duplicated PCA matrices",{

  pbmc_small@int_colData@listData$reducedDims$PCA2 = pbmc_small@int_colData@listData$reducedDims$PCA

  pbmc_small %>%
    mutate(aa = 1) |>
    as_tibble() |>
    ncol() |>
    expect_equal(
      (pbmc_small |> as_tibble() |> ncol()) + 1
    )


})
