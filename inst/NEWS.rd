\name{NEWS}
\title{News for Package \pkg{tidySingleCellExperiment}}

\section{Changes in version 1.19.2, Bioconductor 3.22 Release}{
\itemize{
    \item Soft deprecated \code{bind_rows()} in favor of \code{append_samples()} from ttservice.
    \item Added \code{append_samples()} method for SingleCellExperiment objects.
    \item \code{bind_rows()} is not a generic method in dplyr and may cause conflicts.
    \item Users are encouraged to use \code{append_samples()} instead.
}}

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

