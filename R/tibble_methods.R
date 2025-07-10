#' @name as_tibble
#' @rdname as_tibble
#' @inherit tibble::as_tibble
#' @return `tibble`
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> as_tibble()
#' 
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#' 
#' @importFrom tibble as_tibble
#' @importFrom pkgconfig get_config
#' @importFrom SummarizedExperiment colData
#' @export
as_tibble.SingleCellExperiment <- function(x, ...,
    .name_repair=c("check_unique", "unique", "universal", "minimal"),
    rownames=pkgconfig::get_config("tibble::rownames", NULL)) {
    df <- colData(x) %>%
        as.data.frame() %>%
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
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#' 
#' @importFrom tibble glimpse
#' @export
glimpse.tidySingleCellExperiment <- function(x, width=NULL, ...){
    x %>%
        as_tibble() %>%
        tibble::glimpse(width=width, ...)
}
