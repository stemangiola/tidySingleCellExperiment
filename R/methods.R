


setClass("tidySingleCellExperiment", contains = "SingleCellExperiment")

#' @importFrom methods show
#' @import SingleCellExperiment
#' @importFrom magrittr %>%
setMethod(
    f = "show",
    signature = "tidySingleCellExperiment",
    definition = function(object) {
        object %>%
            print()
    }
)

#' tidy for SingleCellExperiment
#'
#' @param object A SingleCellExperiment object
#'
#' @return A tidySingleCellExperiment object
#'
#' @name tidy
#'
#' @examples
#'
#' tidySingleCellExperiment::pbmc_small %>% tidy()
#' @export
tidy <- function(object) {
    UseMethod("tidy", object)
}

#' @importFrom methods as
#'
#' @param object A SingleCellExperiment object
#'
#' @export
tidy.SingleCellExperiment <- function(object) {
    as(object, "tidySingleCellExperiment")
}

#' Add differential transcription information to a tbl using edgeR.
#'
#' \lifecycle{experimental}
#'
#' @description join_transcripts() extracts and joins information for specific
#'   transcripts
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name join_transcripts
#' @rdname join_transcripts
#'
#' @param .data A tidy SingleCellExperiment object
#' @param transcripts A vector of transcript identifiers to join
#' @param all If TRUE return all
#' @param exclude_zeros If TRUE exclude zero values
#' @param shape Format of the returned table "long" or "wide"
#'
#' @details This function extracts information for specified transcripts and
#'   returns the information in either long or wide format.
#'
#' @return A `tbl` containing the information.for the specified transcripts
#'
#' @examples
#'
#' tidySingleCellExperiment::pbmc_small %>%
#'     tidy() %>%
#'     join_transcripts(transcripts=c("HLA-DRA", "LYZ"))
#' @export
#'
join_transcripts <- function(.data,
                             transcripts = NULL,
                             all = FALSE,
                             exclude_zeros = FALSE,
                             shape = "long") {
    UseMethod("join_transcripts", .data)
}
#' @export
join_transcripts.default <-
    function(.data,
             transcripts = NULL,
             all = FALSE,
             exclude_zeros = FALSE,
             shape = "long") {
        print("This function cannot be applied to this object")
    }
#' @importFrom tidyselect contains
#' @importFrom tidyselect everything
#' @export
join_transcripts.tidySingleCellExperiment <-
    function(.data,
             transcripts = NULL,
             all = FALSE,
             exclude_zeros = FALSE,
             shape = "long") {
        message(data_frame_returned_message)

        my_tibble =
            .data %>%
            as_tibble()

        # Shape is long
        if (shape == "long")
            my_tibble %>%
            left_join(
                get_abundance_sc_long(
                    .data = .data,
                    transcripts = transcripts,
                    all = all,
                    exclude_zeros = exclude_zeros
                ),
                by = "cell"
            ) %>%
            select(cell, transcript, contains("abundance"), everything())

        # Shape if wide
        else
            my_tibble  %>% left_join(get_abundance_sc_wide(
                .data = .data,
                transcripts = transcripts,
                all = all
            ),
            by = "cell")

    }
