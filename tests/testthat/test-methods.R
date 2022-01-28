context('methods test')

data("pbmc_small")
library(dplyr)
test_that("join_features",{


  pbmc_small %>%
    join_features("CD3D") %>%
    slice(1) %>%
    tidySingleCellExperiment::pull(abundance_counts) %>%
    expect_equal(4, tolerance=0.1)


})
