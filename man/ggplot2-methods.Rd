% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot2_methods.R
\name{ggplot}
\alias{ggplot}
\title{Create a new ggplot from a tidySingleCellExperiment object}
\arguments{
\item{.data}{Default dataset to use for plot. If not already a data.frame,
will be converted to one by \code{\link[=fortify]{fortify()}}. If not specified,
must be supplied in each layer added to the plot.}

\item{mapping}{Default list of aesthetic mappings to use for plot.
If not specified, must be supplied in each layer added to the plot.}

\item{...}{Other arguments passed on to methods. Not currently used.}

\item{environment}{DEPRECATED. Used prior to tidy evaluation.}
}
\value{
A ggplot
}
\description{
\code{ggplot()} initializes a ggplot object. It can be used to
declare the input data frame for a graphic and to specify the
set of plot aesthetics intended to be common throughout all
subsequent layers unless specifically overridden.
}
\details{
\code{ggplot()} is used to construct the initial plot object,
and is almost always followed by \code{+} to add component to the
plot. There are three common ways to invoke \code{ggplot()}:

The first method is recommended if all layers use the same
data and the same set of aesthetics, although this method
can also be used to add a layer using data from another
data frame. See the first example below. The second
method specifies the default data frame to use for the plot,
but no aesthetics are defined up front. This is useful when
one data frame is used predominantly as layers are added,
but the aesthetics may vary from one layer to another. The
third method initializes a skeleton \code{ggplot} object which
is fleshed out as layers are added. This method is useful when
multiple data frames are used to produce different layers, as
is often the case in complex graphics.
}
\examples{

library(ggplot2)

tidySingleCellExperiment::pbmc_small \%>\%
    
    tidySingleCellExperiment::ggplot(aes(groups, nCount_RNA)) +
    geom_boxplot()
}
