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

# Get assays
get_all_assays <- function(x) {
  assay_names <- names(assays(x))
  alt_exp_assays <- list()
  alt_exp_assay_names_list <- lapply(altExps(x), assayNames)
  names(assay_names) <- rep("Main", length(assay_names))
  alt_exp_assay_names_df <- stack(alt_exp_assay_names_list)
  alt_exp_assay_names <- paste(alt_exp_assay_names_df$ind, alt_exp_assay_names_df$values, sep = "-")
  names(alt_exp_assay_names) <- alt_exp_assay_names_df$ind
  all_assay_names_df <- rbind(stack(assay_names), alt_exp_assay_names_df)
  all_assay_names <- c(assay_names, alt_exp_assay_names)
  all_assay_names_ext_df <- stack(all_assay_names)
  all_assay_names_ext_df <- cbind(all_assay_names_ext_df, all_assay_names_df$values)
  colnames(all_assay_names_ext_df) <- c("assay_id", "exp_id", "assay_name")
  return(all_assay_names_ext_df)
}

# Get list of features
get_all_features <- function(x) {
  all_assay_names_ext_df <- get_all_assays(x)
  features_lookup <- vector("list", length = length(all_assay_names_ext_df$assay_id))
  main_features <- vector("list", length = 1)
  names(main_features) <- "Main"
  main_features[["Main"]] <- rownames(rowData(x))
  temp_funct <- function(x) rownames(rowData(x))
  alt_exp_features <- lapply(altExps(x), temp_funct)
  feature_df <- stack(c(main_features, alt_exp_features))
  colnames(feature_df) <- c("feature", "exp_id")
  feature_df <- merge(feature_df, all_assay_names_ext_df, by = "exp_id")
  return(feature_df)
}

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom utils tail
#' @importFrom SummarizedExperiment assays
#' @importFrom stats setNames
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
get_abundance_sc_wide <- function(.data, features=NULL, all=FALSE, assay = assays(.data) |> as.list() |> tail(1) |> names(),  prefix = "" ) {
  
  # Solve CRAN warnings
  . <- NULL
  
  # For SCE there is no field for variable features
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
  
  # Get selected features
  feature_df <- get_all_features(.data)
  selected_features <- feature_df[(feature_df$feature %in% features), ]
  selected_features_df <- selected_features[(selected_features$assay_id %in% assay),]
  if(!(nrow(selected_features_df) > 0 && all(selected_features_df$assay_id %in% assay))) stop("tidySingleCellExperiment says: Please specify correct assay.")
  selected_features_exp <- unique(selected_features_df$exp_id)
  if(length(selected_features_exp) > 1) stop("Please avoid mixing features from different experiments.")
  selected_features_assay <- unique(selected_features_df$assay_name)
  
  if(selected_features_exp == "Main") {
    assays(.data)[[assay]] %>%
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
  } else {
    assays(altExps(.data)[[selected_features_exp]])[[selected_features_assay]] %>%
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
get_abundance_sc_long <- function(.data, features = NULL, all = FALSE, exclude_zeros = FALSE) {
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
  
  assay_names <- names(assays(.data))
  
  # Check that I have assay names - can you even have an sce object with no assays?
  if (length(assay_names) == 0) {
    stop("tidySingleCellExperiment says: there are no assays names in the source SingleCellExperiment.")
  }
  
    # Get selected features
  feature_df <- get_all_features(.data)
  selected_features <- feature_df[(feature_df$feature %in% features), ]
  selected_features_exp <- unique(selected_features$exp_id)
  if(length(selected_features_exp) > 1) stop("Please avoid mixing features from different experiments.")
  selected_features_assay_names <- unique(selected_features$assay_id)
  
  if (selected_features_exp == "Main") {
    assays(.data) %>%
      as.list() %>%
      # Take active assay
      map2(
        assay_names,
        ~ .x %>%
          when(
            variable_genes %>% is.null() %>% `!`() ~ .x[variable_genes, , drop = FALSE],
            features %>% is.null() %>% `!`() ~ .x[toupper(rownames(.x)) %in% toupper(features), , drop = FALSE],
            all ~ .x,
            ~ stop("It is not convenient to extract all genes, you should have either variable features or feature list to extract")
          ) %>%
          # Replace 0 with NA
          when(exclude_zeros ~ (.) %>%
                 {
                   x <- (.)
                   x[x == 0] <- NA
                   x
                 }, ~ (.)) %>%
          as.matrix() %>%
          data.frame(check.names = FALSE) %>%
          as_tibble(rownames = ".feature") %>%
          tidyr::pivot_longer(
            cols = -.feature,
            names_to = c_(.data)$name,
            values_to = ".abundance" %>% paste(.y, sep = "_"),
            values_drop_na = TRUE
          )
        # %>%
        # mutate_if(is.character, as.factor) %>%
      ) %>%
      Reduce(function(...) full_join(..., by = c(".feature", c_(.data)$name)), .)
  } else {
    assays(altExps(.data)[[selected_features_exp]]) %>%
      as.list() %>%
      # Take active assay
      map2(
        selected_features_assay_names,
        ~ .x %>%
          when(
            variable_genes %>% is.null() %>% `!`() ~ .x[variable_genes, , drop = FALSE],
            features %>% is.null() %>% `!`() ~ .x[toupper(rownames(.x)) %in% toupper(features), , drop = FALSE],
            all ~ .x,
            ~ stop("It is not convenient to extract all genes, you should have either variable features or feature list to extract")
          ) %>%
          # Replace 0 with NA
          when(exclude_zeros ~ (.) %>%
                 {
                   x <- (.)
                   x[x == 0] <- NA
                   x
                 }, ~ (.)) %>%
          as.matrix() %>%
          data.frame(check.names = FALSE) %>%
          as_tibble(rownames = ".feature") %>%
          tidyr::pivot_longer(
            cols = -.feature,
            names_to = c_(.data)$name,
            values_to = ".abundance" %>% paste(.y, sep = "_"),
            values_drop_na = TRUE
          )
        # %>%
        # mutate_if(is.character, as.factor) %>%
      ) %>%
      Reduce(function(...) full_join(..., by = c(".feature", c_(.data)$name)), .)
  }
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
  
  col_to_exclude <-
    
    # special_datasets_to_tibble(SingleCellExperiment_object) |>
    # colnames()
    get_special_columns(SingleCellExperiment_object) |>
    
    
    # I need this in case we have multiple reduced dimension data frames with overlapping names of the columns.
    # For example multiple PCA versions
    vctrs::vec_as_names(repair = "unique") |>
    
    # To avoid name change by the bind_cols of as_tibble
    trick_to_avoid_renaming_of_already_unique_columns_by_dplyr()
  
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

#' @importFrom purrr reduce
#' @importFrom tibble enframe
special_datasets_to_tibble = function(.singleCellExperiment, ...){
  x =
    .singleCellExperiment |>
    get_special_datasets(...) %>%
    map(~ .x %>% when(
      
      # If row == 1 do a trick
      dim(.) %>% is.null() ~ {
        (.) %>%
          tibble::enframe() %>%
          spread(name, value)
      },
      
      # Otherwise continue normally
      ~ as_tibble(.)
    )) %>%
    reduce(bind_cols)
  
  # To avoid name change by the bind_cols of as_tibble
  colnames(x) = colnames(x) |> trick_to_avoid_renaming_of_already_unique_columns_by_dplyr()
  
  x
}

#' @importFrom stringr str_replace_all
trick_to_avoid_renaming_of_already_unique_columns_by_dplyr = function(x){
  x |> str_replace_all("\\.\\.\\.", "___")
}

#' Get specific annotation columns
#'
#' @keywords internal
#' @noRd
#' 
#' @importFrom rlang enquo
#' @importFrom purrr map
#' @importFrom dplyr distinct_at
#' @importFrom magrittr equals
#' @importFrom dplyr vars
#' 
#' @param .data A tibble
#' @param .col A vector of column names
#' 
#' @return A character
get_specific_annotation_columns = function(.data, .col){
  
  # Comply with CRAN NOTES
  . = NULL
  
  # Make col names
  .col = enquo(.col)
  
  # x-annotation df
  n_x = .data %>% distinct_at(vars(!!.col)) %>% nrow
  
  # element wise columns
  .data %>%
    select(-!!.col) %>%
    colnames %>%
    map(
      ~
        .x %>%
        when(
          .data %>%
            distinct_at(vars(!!.col, .x)) %>%
            nrow %>%
            equals(n_x) ~ (.),
          ~ NULL
        )
    ) %>%
    
    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist
}

#' Subset columns
#'
#' @keywords internal
#' @noRd
#' 
#' @importFrom rlang enquo
#' 
#' @param .data A tibble
#' @param .column A vector of column names
#'
#' @return A tibble
subset = function(.data, .column)	{
  
  # Make col names
  .column = enquo(.column)
  
  # Check if column present
  if(quo_names(.column) %in% colnames(.data) %>% all %>% `!`)
    stop("nanny says: some of the .column specified do not exist in the input data frame.")
  
  .data %>%
    
    # Selecting the right columns
    select(	!!.column,	get_specific_annotation_columns(.data, !!.column)	) %>%
    distinct()
}

feature__ = get_special_column_name_symbol(".feature")
sample__ = get_special_column_name_symbol(".sample")
