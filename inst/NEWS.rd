\name{NEWS}
\title{News for Package \pkg{tidySingleCellExperiment}}

\section{Changes in version 1.4.0, Bioconductor 3.14 Release}{
\itemize{
    \item Improved sample_n, and sample_frac functions.
    \item Add join_features prefix.
    \item Dropped tidy method as never needed.
    \item Add unnest_tidySingleCellExperiment for nested data that was not produce with tidySingleCellExperiment::nest() but rather with tidyr::nest().
}}

\section{Changes in version 1.5.1, Bioconductor 3.15 Release}{
\itemize{
    \item Rely of ttservice package for shared function with tidySingleCellExperiment to avoid clash
    \item Use .cell for cell column name to avoid errors when cell column is defined by the user
}}

\section{Changes in version 1.19.2, Bioconductor 3.22 Release}{
\itemize{
    \item \strong{BREAKING CHANGE}: Changed default shape parameter in \code{join_features()} from "long" to "wide". 
    This means that \code{join_features()} now returns a SingleCellExperiment object by default instead of a tibble. 
    To get the old behavior, explicitly specify \code{shape="long"}.
}}

