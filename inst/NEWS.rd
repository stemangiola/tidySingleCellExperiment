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

