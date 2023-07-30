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
#' @importFrom magrittr %>%
#' @export
NULL
