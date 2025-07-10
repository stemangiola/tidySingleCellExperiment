#' @name ggplot
#' @rdname ggplot
#' @inherit ggplot2::ggplot
#' @title Create a new \code{ggplot} from a \code{tidySingleCellExperiment}
#' @return `ggplot`
#'
#' @examples
#' library(ggplot2)
#' data(pbmc_small)
#' pbmc_small |> 
#'   ggplot(aes(groups, nCount_RNA)) +
#'   geom_boxplot()
#'     
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#'     
#' @importFrom purrr map
#' @importFrom rlang quo_name
#' @importFrom ggplot2 aes ggplot
#' @export
ggplot.SingleCellExperiment <- function(data=NULL, 
    mapping=aes(), ..., environment=parent.frame()) {
    
    # Deprecation of special column names
    .cols <- mapping %>% 
        unlist() %>% map(~ quo_name(.x)) %>% 
        unlist() %>% as.character()
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }
    
    data %>%
        as_tibble() %>%
        ggplot2::ggplot(mapping=mapping)
}

# addressing R CMD CHECK NOTE "no visible global function definition for 'aes'"
globalVariables("aes", "tidySingleCellExperiment")