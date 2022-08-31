#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment colData
#'
#' @keywords internal
#'
#' @param .data A tidySingleCellExperiment
#'
#' @noRd
to_tib <- function(.data) {
    colData(.data) %>%
        as.data.frame() %>%
        as_tibble(rownames=c_(.data)$name)
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
#' @importFrom SummarizedExperiment assays
#'
#' @param .data A tidySingleCellExperiment
#' @param features A character
#' @param all A boolean
#' @param ... Parameters to pass to join wide, i.e. assay name to extract feature abundance from
#'
#'
#' @return A tidySingleCellExperiment object
#'
#'
#' @noRd
get_abundance_sc_wide <- function(.data, features=NULL, all=FALSE, assay = assays(.data) %>% as.list() %>% tail(1) %>% names,  prefix = "" ) {

    # Solve CRAN warnings
    . <- NULL

    # For SCE there is not filed for variable features
    variable_feature <- c()

    # Check if output would be too big without forcing
    if (
        length(variable_feature) == 0 &
            is.null(features) &
            all == FALSE
    ) {
        stop("
                Your object does not contain variable feature labels,
                feature argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of feature names
                3. set all=TRUE (this will output a very large object, does your computer have enough RAM?)
                ")
    }

    # Get variable features if existing
    if (
        length(variable_feature) > 0 &
            is.null(features) &
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
      .[[assay]] %>%
        when(
            variable_genes %>% is.null() %>% `!`() ~ (.)[variable_genes, , drop=FALSE],
            features %>% is.null() %>% `!`() ~ (.)[features, , drop=FALSE],
            ~ stop("It is not convenient to extract all genes, you should have either variable features or feature list to extract")
        ) %>%
        as.matrix() %>%
        t() %>%
        as_tibble(rownames=c_(.data)$name) %>%

        # Add prefix
        setNames(c(c_(.data)$name, sprintf("%s%s", prefix, colnames(.)[-1])))
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
#' @importFrom SummarizedExperiment assays
#'
#' @param .data A tidySingleCellExperiment
#' @param features A character
#' @param all A boolean
#' @param exclude_zeros A boolean
#'
#' @return A tidySingleCellExperiment object
#'
#'
#' @noRd
get_abundance_sc_long <- function(.data, features=NULL, all=FALSE, exclude_zeros=FALSE) {

    # Solve CRAN warnings
    . <- NULL

    # For SCE there is not filed for variable features
    variable_feature <- c()

    # Check if output would be too big without forcing
    if (
        length(variable_feature) == 0 &
            is.null(features) &
            all == FALSE
    ) {
        stop("
                Your object does not contain variable feature labels,
                feature argument is empty and all arguments are set to FALSE.
                Either:
                1. use detect_variable_features() to select variable feature
                2. pass an array of feature names
                3. set all=TRUE (this will output a very large object, does your computer have enough RAM?)
                ")
    }


    # Get variable features if existing
    if (
        length(variable_feature) > 0 &
            is.null(features) &
            all == FALSE
    ) {
        variable_genes <- variable_feature
    } # Else
    else {
        variable_genes <- NULL
    }

    assay_names <- assays(.data) %>% names()

    # Check that I have assay manes
    if(length(assay_names) == 0)
      stop("tidySingleCellExperiment says: there are no assays names in the source SingleCellExperiment.")

    assays(.data) %>%
        as.list() %>%

        # Take active assay
        map2(
            assay_names,

            ~ .x %>%
                when(
                    variable_genes %>% is.null() %>% `!`() ~ .x[variable_genes, , drop=FALSE],
                    features %>% is.null() %>% `!`() ~ .x[toupper(rownames(.x)) %in% toupper(features), , drop=FALSE],
                    all ~ .x,
                    ~ stop("It is not convenient to extract all genes, you should have either variable features or feature list to extract")
                ) %>%

                # Replace 0 with NA
                when(exclude_zeros ~ (.) %>% {
                    x <- (.)
                    x[x == 0] <- NA
                    x
                }, ~ (.)) %>%
                as.matrix() %>%
                data.frame(check.names = FALSE) %>%
                as_tibble(rownames=".feature") %>%
                tidyr::pivot_longer(
                    cols=- .feature,
                    names_to=c_(.data)$name,
                    values_to=".abundance" %>% paste(.y, sep="_"),
                    values_drop_na=TRUE
                )
            # %>%
            # mutate_if(is.character, as.factor) %>%
        ) %>%
        Reduce(function(...) full_join(..., by=c(".feature", c_(.data)$name)), .)
}

#' @importFrom dplyr select_if
#' @importFrom S4Vectors DataFrame
#'
#' @keywords internal
#'
#' @param .data A tibble
#' @param SingleCellExperiment_object A tidySingleCellExperiment
#'
#' @noRd
as_meta_data <- function(.data, SingleCellExperiment_object) {

    # Solve CRAN warnings
    . <- NULL

    col_to_exclude <- get_special_columns(SingleCellExperiment_object)

    .data_df =
      .data %>%
      select_if(!colnames(.) %in% col_to_exclude) %>%
      data.frame()

    # Set row names in a robust way. the argument row.names of the data.frame function does not work for 1-row data frames
    rownames(.data_df) = .data_df |> pull(!!c_(SingleCellExperiment_object)$symbol)
    .data_df = .data_df |> select(-!!c_(SingleCellExperiment_object)$symbol)

    .data_df %>% DataFrame()

}

#' @importFrom purrr map_chr
#'
#' @keywords internal
#'
#' @param SingleCellExperiment_object A tidySingleCellExperiment
#'
#' @noRd
#'
get_special_columns <- function(SingleCellExperiment_object) {
    get_special_datasets(SingleCellExperiment_object) %>%
        map(~ .x %>% colnames()) %>%
        unlist() %>%
        as.character()
}

get_special_datasets <- function(SingleCellExperiment_object, n_dimensions_to_return = Inf) {
    rd <- SingleCellExperiment_object@int_colData@listData$reducedDims

    map2(rd %>% as.list(), names(rd), ~ {
        mat <- .x[, seq_len(min(n_dimensions_to_return, ncol(.x))), drop=FALSE]

        # Set names as SCE is much less constrained and there could be missing names
        if (length(colnames(mat)) == 0) colnames(mat) <- sprintf("%s%s", .y, seq_len(ncol(mat)))

        mat
    })
}

get_needed_columns <- function(.data) {

  c(c_(.data)$name)
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

data_frame_returned_message = "tidySingleCellExperiment says: A data frame is returned for independent data analysis."
duplicated_cell_names = "tidySingleCellExperiment says: This operation lead to duplicated cell names. A data frame is returned for independent data analysis."

# This function is used for the change of special sample column to .sample
# Check if "sample" is included in the query and is not part of any other existing annotation
#' @importFrom stringr str_detect
#' @importFrom stringr regex
is_sample_feature_deprecated_used = function(.data, user_columns, use_old_special_names = FALSE){

  old_standard_is_used_for_cell =
    (
      ( any(str_detect(user_columns  , regex("\\bcell\\b"))) & !any(str_detect(user_columns  , regex("\\W*(\\.cell)\\W*")))  ) |
        "cell" %in% user_columns
    ) &
    !"cell" %in% colnames(colData(.data))

  old_standard_is_used = old_standard_is_used_for_cell

  if(old_standard_is_used){
    warning("tidySingleCellExperiment says: from version 1.3.1, the special columns including cell id (colnames(se)) has changed to \".cell\". This dataset is returned with the old-style vocabulary (cell), however we suggest to update your workflow to reflect the new vocabulary (.cell)")

    use_old_special_names = TRUE
  }

  use_old_special_names
}

get_special_column_name_symbol = function(name){
  list(name = name, symbol = as.symbol(name))
}

# Key column names
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
ping_old_special_column_into_metadata = function(.data){

  metadata(.data)$cell__ = get_special_column_name_symbol("cell")

  .data
}

get_special_column_name_cell = function(name){
  list(name = name, symbol = as.symbol(name))
}

cell__ = get_special_column_name_symbol(".cell")

#' @importFrom S4Vectors metadata
c_ =  function(x){
  # Check if old deprecated columns are used
  if("cell__" %in% names(metadata(x))) cell__ = metadata(x)$cell__
  return(cell__)
}

#' Add attribute to abject
#'
#' @keywords internal
#' @noRd
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}
