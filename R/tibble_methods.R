#' @name as_tibble
#' @rdname as_tibble
#' @inherit tibble::as_tibble
#' @return `tibble`
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> as_tibble()
#' 
#' @importFrom tibble as_tibble
#' @importFrom pkgconfig get_config
#' @importFrom SummarizedExperiment colData
#' @export
as_tibble.SingleCellExperiment <- function(x, ...,
    .name_repair=c("check_unique", "unique", "universal", "minimal"),
    rownames=pkgconfig::get_config("tibble::rownames", NULL)) {
    df <- colData(x) %>%
        as(Class = "data.frame", strict = FALSE) %>%
        tibble::as_tibble(rownames=c_(x)$name)
    # Attach reduced dimensions only if 
    # there are any and for special datasets
    if (length(reducedDims(x))) {
        fd <- special_datasets_to_tibble(x, ...)
        df <- bind_cols(df, fd)
    }
    return(df)
}

#' @name glimpse
#' @rdname glimpse
#' @inherit pillar::glimpse
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> glimpse()
#' 
#' @importFrom tibble glimpse
#' @export
glimpse.tidySingleCellExperiment <- function(x, width=NULL, ...){
    x %>%
        as_tibble() %>%
        tibble::glimpse(width=width, ...)
}
