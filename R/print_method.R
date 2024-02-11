# This file is a replacement of the unexported functions in the tibble
# package, in order to specify "tibble abstraction in the header"

#' @name tbl_format_header
#' @rdname tbl_format_header
#' @inherit pillar::tbl_format_header
#'
#' @examples
#' # TODO
#'
#' @importFrom rlang names2
#' @importFrom pillar align
#' @importFrom pillar get_extent
#' @importFrom pillar style_subtle
#' @importFrom pillar tbl_format_header
#' @export

tbl_format_header.tidySingleCellExperiment <- function(x, setup, ...) {

    number_of_features <- x |> attr("number_of_features")
    assay_names <- x |> attr("assay_names")


    # Change name
    named_header <- setup$tbl_sum
    names(named_header) <- "A SingleCellExperiment-tibble abstraction"

    if (all(names2(named_header) == "")) {
        header <- named_header
    } else {
        header <- paste0(
            align(paste0(names2(named_header), ":"), space=NBSP),
            " ", named_header) %>%
            # Add further info single-cell

          append(sprintf(
              "\033[90m Features=%s | Cells=%s | Assays=%s\033[39m",
              number_of_features, nrow(x),
              paste(assay_names, collapse=", ")), after=1)

    }
    style_subtle(pillar___format_comment(header, width=setup$width))
}

#' @name formatting
#' @rdname formatting
#' @aliases print
#' @inherit tibble::formatting
#' @param n_extra number of extra lines
#' @return Prints a message to the console describing
#'   the contents of the `tidySingleCellExperiment`.
#'
#' @examples
#' data(pbmc_small)
#' print(pbmc_small)
#'
#' @importFrom vctrs new_data_frame
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleCellExperiment altExpNames
#' @export

print.SingleCellExperiment <- function(x, ..., n = NULL, width = NULL, n_extra = NULL) {
  if (length(names(altExps(x))) > 0) {
    alt_exp_assays <- list()
    assay_names_list <- lapply(altExps(x), assayNames)
    assay_names_df <- stack(assay_names_list)
    assay_names_string <- c(assayNames(x), paste(assay_names_df$ind, assay_names_df$values, sep = "-")) |>
      paste(collapse = ", ")
  } else {
    assay_names_string <- paste(assayNames(x), collapse = ", ")
  }
  x |>
    as_tibble(n_dimensions_to_return = 5) |>
    new_data_frame(class = c("tidySingleCellExperiment", "tbl")) %>%
    add_attr(nrow(x), "number_of_features") %>%
    add_attr(assay_names_string, "assay_names") %>%
    print()

  invisible(x)
}
