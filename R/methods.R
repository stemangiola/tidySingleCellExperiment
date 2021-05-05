setClass("tidySingleCellExperiment", contains = "SingleCellExperiment")

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
#' tidySingleCellExperiment::pbmc_small
#' @export
tidy <- function(object) {
    UseMethod("tidy", object)
}

#' @importFrom methods as
#' @importFrom lifecycle deprecate_warn
#'
#' @param object A SingleCellExperiment object
#'
#' @export
tidy.SingleCellExperiment <- function(object) {

    # DEPRECATE
    deprecate_warn(
        when = "1.1.1",
        what = "tidy()",
        details = "tidySingleCellExperiment says: tidy() is not needed anymore."
    )

    object
}

setMethod(
    f = "show",
    signature = "SingleCellExperiment",
    definition = function(object) {
        if (isTRUE(x = getOption(x = "restore_SingleCellExperiment_show", default = FALSE))) {
            f <- getMethod(
                f = "show",
                signature = "SingleCellExperiment",
                where = asNamespace(ns = "SingleCellExperiment")
            )
            f(object = object)
        } else {
            object %>%
                print()
        }
    }
)

#' Add differential featureion information to a tbl using edgeR.
#'
#' \lifecycle{experimental}
#'
#' @description join_features() extracts and joins information for specific
#'   features
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name join_features
#' @rdname join_features
#'
#' @param .data A tidy SingleCellExperiment object
#' @param features A vector of feature identifiers to join
#' @param all If TRUE return all
#' @param exclude_zeros If TRUE exclude zero values
#' @param shape Format of the returned table "long" or "wide"
#' @param ... Parameters to pass to join wide, i.e. assay name to extract feature abundance from
#'
#' @details This function extracts information for specified features and
#'   returns the information in either long or wide format.
#'
#' @return A `tbl` containing the information.for the specified features
#'
#' @examples
#'
#' tidySingleCellExperiment::pbmc_small %>%
#'
#'     join_features(features=c("HLA-DRA", "LYZ"))
#' @export
#'
join_features <- function(.data,
                             features = NULL,
                             all = FALSE,
                             exclude_zeros = FALSE,
                             shape = "long", ...) {
    UseMethod("join_features", .data)
}
#' @export
join_features.default <-
    function(.data,
             features = NULL,
             all = FALSE,
             exclude_zeros = FALSE,
             shape = "long", ...) {
        print("This function cannot be applied to this object")
    }
#' @importFrom tidyselect contains
#' @importFrom tidyselect everything
#' @export
join_features.SingleCellExperiment <-
    function(.data,
             features = NULL,
             all = FALSE,
             exclude_zeros = FALSE,
             shape = "long", ...) {

        # CRAN Note
        cell = NULL
        feature= NULL

        # Shape is long
        if (shape == "long")
          .data %>%
            left_join(
                get_abundance_sc_long(
                    .data = .data,
                    features = features,
                    all = all,
                    exclude_zeros = exclude_zeros
                ),
                by = "cell"
            ) %>%
            select(cell, feature, contains("abundance"), everything())

        # Shape if wide
        else
          .data  %>% left_join(get_abundance_sc_wide(
                .data = .data,
                features = features,
                all = all, ...
            ),
            by = "cell")

    }
