#' @importFrom tibble as_tibble
#' @importFrom SingleCellExperiment colData
#'
#' @keywords internal
#'
#' @param .data A tidySCE
#'
#' @noRd
to_tib <- function(.data) {
    colData(.data) %>%
        as.data.frame() %>%
        as_tibble(rownames="cell")
}

# Greater than
gt <- function(a, b) {
    a > b
}

# Smaller than
st <- function(a, b) {
    a < b
}

# Negation
not <- function(is) {
    !is
}

# Raise to the power
pow <- function(a, b) {
    a^b
}

# Equals
eq <- function(a, b) {
    a == b
}

prepend <- function(x, values, before=1) {
    n <- length(x)
    stopifnot(before > 0 && before <= n)
    if (before == 1) {
        c(values, x)
    }
    else {
        c(x[seq_len(before - 1)], values, x[before:n])
    }
}
#' Add class to abject
#'
#'
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class <- function(var, name) {
    if (!name %in% class(var)) class(var) <- prepend(class(var), name)

    var
}

#' Remove class to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
#' @keywords internal
drop_class <- function(var, name) {
    class(var) <- class(var)[!class(var) %in% name]
    var
}

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom utils tail
#'
#' @param .data A tidySCE
#' @param transcripts A character
#' @param all A boolean
#'
#'
#' @return A tidySCE object
#'
#'
#' @noRd
get_abundance_sc_wide <- function(.data, transcripts=NULL, all=FALSE) {

    # Solve CRAN warnings
    . <- NULL

    # For SCE there is not filed for variable features
    variable_feature <- c()

    # Check if output would be too big without forcing
    if (
        length(variable_feature) == 0 &
            is.null(transcripts) &
            all == FALSE
    ) {
        stop("
                Your object does not contain variable transcript labels,
                transcript argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of transcript names
                3. set all=TRUE (this will output a very large object, does your computer have enough RAM?)
                ")
    }

    # Get variable features if existing
    if (
        length(variable_feature) > 0 &
            is.null(transcripts) &
            all == FALSE
    ) {
        variable_genes <- variable_feature
    } # Else
    else {
        variable_genes <- NULL
    }

    # Just grub last assay
    assays(.data) %>%
        as.list() %>%
        tail(1) %>%
        .[[1]] %>%
        when(
            variable_genes %>% is.null() %>% `!`() ~ (.)[variable_genes, , drop=FALSE],
            transcripts %>% is.null() %>% `!`() ~ (.)[transcripts, , drop=FALSE],
            ~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
        ) %>%
        as.matrix() %>%
        t() %>%
        as_tibble(rownames="cell")
}

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom purrr when
#' @importFrom purrr map2
#'
#' @param .data A tidySCE
#' @param transcripts A character
#' @param all A boolean
#' @param exclude_zeros A boolean
#'
#' @return A tidySCE object
#'
#'
#' @noRd
get_abundance_sc_long <- function(.data, transcripts=NULL, all=FALSE, exclude_zeros=FALSE) {

    # Solve CRAN warnings
    . <- NULL

    # For SCE there is not filed for variable features
    variable_feature <- c()

    # Check if output would be too big without forcing
    if (
        length(variable_feature) == 0 &
            is.null(transcripts) &
            all == FALSE
    ) {
        stop("
                Your object does not contain variable transcript labels,
                transcript argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of transcript names
                3. set all=TRUE (this will output a very large object, does your computer have enough RAM?)
                ")
    }


    # Get variable features if existing
    if (
        length(variable_feature) > 0 &
            is.null(transcripts) &
            all == FALSE
    ) {
        variable_genes <- variable_feature
    } # Else
    else {
        variable_genes <- NULL
    }

    assay_names <- assays(.data) %>% names()


    assays(.data) %>%
        as.list() %>%

        # Take active assay
        map2(
            assay_names,

            ~ .x %>%
                when(
                    variable_genes %>% is.null() %>% `!`() ~ .x[variable_genes, , drop=FALSE],
                    transcripts %>% is.null() %>% `!`() ~ .x[toupper(rownames(.x)) %in% toupper(transcripts), , drop=FALSE],
                    all ~ .x,
                    ~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
                ) %>%

                # Replace 0 with NA
                when(exclude_zeros ~ (.) %>% {
                    x <- (.)
                    x[x == 0] <- NA
                    x
                }, ~ (.)) %>%
                as.matrix() %>%
                data.frame() %>%
                as_tibble(rownames="transcript") %>%
                tidyr::pivot_longer(
                    cols=-transcript,
                    names_to="cell",
                    values_to="abundance" %>% paste(.y, sep="_"),
                    values_drop_na=TRUE
                )
            # %>%
            # mutate_if(is.character, as.factor) %>%
        ) %>%
        Reduce(function(...) left_join(..., by=c("transcript", "cell")), .)
}

#' @importFrom dplyr select_if
#' @importFrom S4Vectors DataFrame
#'
#' @keywords internal
#'
#' @param .data A tibble
#' @param SingleCellExperiment_object A tidySCE
#'
#' @noRd
as_meta_data <- function(.data, SingleCellExperiment_object) {

    # Solve CRAN warnings
    . <- NULL

    col_to_exclude <- get_special_columns(SingleCellExperiment_object)

    .data %>%
        select_if(!colnames(.) %in% col_to_exclude) %>%
        # select(-one_of(col_to_exclude)) %>%
        data.frame(row.names="cell") %>%
        DataFrame()
}

#' @importFrom purrr map_chr
#'
#' @keywords internal
#'
#' @param SingleCellExperiment_object A tidySCE
#'
#' @noRd
#'
get_special_columns <- function(SingleCellExperiment_object) {
    get_special_datasets(SingleCellExperiment_object) %>%
        map(~ .x %>% colnames()) %>%
        unlist() %>%
        as.character()
}

get_special_datasets <- function(SingleCellExperiment_object) {
    rd <- SingleCellExperiment_object@int_colData@listData$reducedDims

    map2(rd %>% as.list(), names(rd), ~ {
        mat <- .x[, seq_len(min(5, ncol(.x))), drop=FALSE]

        # Set names as SCE is much less constrained and there could be missing names
        if (length(colnames(mat)) == 0) colnames(mat) <- sprintf("%s%s", .y, seq_len(ncol(mat)))

        mat
    })
}

get_needed_columns <- function() {
    # c("cell",  "orig.ident", "nCount_RNA", "nFeature_RNA")
    c("cell")
}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {
    v <- quo_name(quo_squash(v))
    gsub("^c\\(|`|\\)$", "", v) %>%
        strsplit(", ") %>%
        unlist()
}

#' @importFrom purrr when
#' @importFrom dplyr select
#' @importFrom rlang expr
#' @importFrom tidyselect eval_select
select_helper <- function(.data, ...) {
    loc <- tidyselect::eval_select(expr(c(...)), .data)

    dplyr::select(.data, loc)
}

data_frame_returned_message = "tidySCE says: A data frame is returned for independent data analysis."
duplicated_cell_names = "tidySCE says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis."
