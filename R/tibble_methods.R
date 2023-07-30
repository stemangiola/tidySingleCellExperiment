#' @rdname as_tibble
#' @inherit tibble::as_tibble
#' 
#' @examples
#' pbmc_small |> as_tibble()
#' 
#' @importFrom tibble as_tibble
#' @importFrom pkgconfig get_config
#' @importFrom SummarizedExperiment colData
#' @export
as_tibble.SingleCellExperiment <- function(x, ...,
    .name_repair=c("check_unique", "unique", "universal", "minimal"),
    rownames=pkgconfig::get_config("tibble::rownames", NULL)) {
    colData(x) %>%
        as.data.frame() %>%
        tibble::as_tibble(rownames=c_(x)$name) %>%


        # Attach reduced dimensions
        when(

            # Only if I have reduced dimensions and special datasets
            ncol(x@int_colData@listData$reducedDims) > 0 ~ (.) %>% bind_cols(
              special_datasets_to_tibble(x, ...)
            ),

            # Otherwise skip
            ~ (.)
        )
}

#' @rdname glimpse
#' @inherit pillar::glimpse
#'
#' @examples
#' pbmc_small |> glimpse()
#' 
#' @importFrom tibble glimpse
#' @export
glimpse.tidySingleCellExperiment = function(x, width = NULL, ...){
    x %>%
        as_tibble() %>%
        tibble::glimpse(width = width, ...)
}
