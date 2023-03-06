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

test_that("aggregate_cells() returns expected values", {
  
  # Create pseudo-bulk object for testing
  pbmc_pseudo_bulk <- pbmc_small |>
    tidySingleCellExperiment::aggregate_cells(c(groups, ident), assays = "counts")
  
  # Check row length is unchanged
  pbmc_pseudo_bulk |>
    nrow() |>
    expect_equal(pbmc_small |> nrow())
  
  # Check column length is correctly modified
  pbmc_pseudo_bulk |> 
    ncol() |>
    expect_equal(pbmc_small |>
                   as_tibble() |>
                   select(groups, ident) |>
                   unique() |>
                   nrow()
    )
  
  # Spot check for correctly aggregated count value of ACAP1 gene
  assay(pbmc_pseudo_bulk, "counts")["ACAP1", "g1___0"] |>
    expect_equal(assay(pbmc_small, "counts")["ACAP1", pbmc_small |>
                                               as_tibble() |>
                                               filter(groups == "g1", ident == 0) |>
                                               pull(.cell)] |>
                   sum())
})
