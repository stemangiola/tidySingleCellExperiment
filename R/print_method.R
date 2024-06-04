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
                paste(assay_names, collapse=", ")
            ), after=1)
    }
    style_subtle(pillar___format_comment(header, width=setup$width))
}

#' @name formatting
#' @rdname formatting
#' @aliases print
#' @inherit tibble::formatting
#' @return Prints a message to the console describing
#'   the contents of the `tidySingleCellExperiment`.
#'
#' @param x A SingleCellExperiment object.
#' @param n (optional) The maximum number of rows to show.
#' @param width (optional) The maximum width of the output.
#' @examples
#' data(pbmc_small)
#' print(pbmc_small)
#'
#' @importFrom vctrs new_data_frame
#' @importFrom SummarizedExperiment assayNames
#' @importFrom crayon red
#' @importFrom dplyr mutate across all_of
#' @importFrom crayon red
#' @export
print.SingleCellExperiment <- function(x, ..., n=NULL, width=NULL) {
    # Retrieve special view-only columns
    special_cols <- get_special_columns(x)

    # Convert to tibble
    tibble_x <- x |> as_tibble(n_dimensions_to_return=5)

    # Modify column titles for special columns
    col_titles <- names(tibble_x)
    modified_col_titles <- ifelse(col_titles %in% special_cols, 
                                    crayon::red(as.character(col_titles)), 
                                    col_titles)

    names(tibble_x) <- modified_col_titles

    # Print the tibble with modified column titles
    new_data_frame(tibble_x, class=c("tidySingleCellExperiment", "tbl")) %>%
        add_attr(nrow(x), "number_of_features") %>%
        add_attr(assayNames(x), "assay_names") %>%
        print(n = n, width = width)

    invisible(x)
}