#' @name %>%
#' @rdname pipe
#' @keywords internal
#' 
#' @title Pipe operator
#' @description See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @usage lhs \%>\% rhs
#' @param lhs A value or the `magrittr` placeholder.
#' @param rhs A function call using the `magrittr` semantics.
#' @return The result of calling `rhs(lhs)`.
#' 
#' @examples
#' `%>%` <- magrittr::`%>%`
#' letters %>% head(n=3)
#' 
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, et al. Welcome to the tidyverse. Journal of Open Source Software. 2019;4(43):1686. https://doi.org/10.21105/joss.01686
#' 
#' @importFrom magrittr %>%
#' @export
NULL
