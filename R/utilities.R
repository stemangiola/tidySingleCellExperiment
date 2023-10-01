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
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class <- function(var, name) {
    if (!name %in% class(var)) 
        class(var) <- prepend(class(var), name)
    return(var)
}

#' Remove class to abject
#'
#' @keywords internal
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

#' get abundance wide
#'
#' @keywords internal
#'
#' @importFrom magrittr "%$%"
#' @importFrom utils tail
#' @importFrom stats setNames
#' @importFrom purrr reduce
#' @importFrom dplyr full_join
#' @importFrom SummarizedExperiment assay assayNames
#'
#' @param .data A `tidySingleCellExperiment`
#' @param features A character
#' @param all A boolean
#' @param ... Parameters to pass to join wide, i.e., 
#'   `assay` to extract feature abundances from
#'
#' @return A tidySingleCellExperiment object
#'
#' @noRd
get_abundance_sc_wide <- function(.data, features=NULL, all=FALSE, prefix="", ...) {
    
  arg_list <- c(mget(ls(environment(), sorted=F)), match.call(expand.dots=F)$...)
  assays_to_use <- eval(arg_list$assays)
  if(is.null(assays_to_use)) stop("Please provide assay names")
  
  
  # Solve CRAN warnings
  . <- NULL
  
  # For SCE there is not filed for variable features
  variable_feature <- c()
  
  # Check if output would be too big without forcing
  if (isFALSE(all) && is.null(features)) {
    if (!length(variable_feature)) {
      stop("Your object does not contain variable feature labels,\n",
           " feature argument is empty and all arguments are set to FALSE.\n",
           " Either:\n",
           " 1. use detect_variable_features() to select variable feature\n",
           " 2. pass an array of feature names\n",
           " 3. set all=TRUE (this will output a very large object;",
           " does your computer have enough RAM?)")
    } else {
      # Get variable features if existing
      variable_genes <- variable_feature
    }
  } else {
    variable_genes <- NULL
  }
  
  if (!is.null(variable_genes)) {
    gs <- variable_genes
  } else if (!is.null(features)) {
    gs <- features
  } else {
    stop("It is not convenient to extract all genes.",
         " You should have either variable features,",
         " or a feature list to extract.")
  }
  # Get selected features and assays
  feature_df <- get_all_features(.data)
  selected_features <- feature_df[(feature_df$feature %in% gs), ]
  selected_features <- selected_features[selected_features$assay_id %in% assays_to_use,]
  selected_experiments_list <- split(x = selected_features, f = as.character(selected_features$exp_id))
  extract_feature_values <- function(exp) {
    selected_features_exp <- as.character(unique(exp$exp_id))
    selected_features_assay <- as.character(unique(exp$assay_name))
    selected_features_assay_names <- as.character(unique(exp$assay_id))
    if(selected_features_exp == "Main") {
      selected_features_from_exp <- rownames(assay(.data, selected_features_assay_names))[(rownames(assay(.data, selected_features_assay_names)) %in% gs)]
      mtx <- assay(.data, selected_features_assay_names)[selected_features_from_exp,]
      if(is.null(dim(mtx))) mtx <- matrix(mtx, byrow = TRUE, nrow = 1, ncol = length(mtx))
      mtx %>%
        `colnames<-`(colnames(.data)) %>%
        as.matrix() %>% t() %>%
        as_tibble(rownames=c_(.data)$name, .name_repair = "minimal") %>%
        setNames(c(c_(.data)$name, sprintf("%s%s", prefix, selected_features_from_exp)))
    } else {
      selected_features_from_exp <- rownames(altExps(.data)[[selected_features_exp]])[(rownames(altExps(.data)[[selected_features_exp]]) %in% gs)]
      mtx <- assay(altExps(.data)[[selected_features_exp]], selected_features_assay)[selected_features_from_exp,]
      if(is.null(dim(mtx))) mtx <- matrix(mtx, byrow = TRUE, nrow = 1, ncol = length(mtx))
      mtx %>%
        `colnames<-`(colnames(.data)) %>%
        as.matrix() %>% t() %>%
        as_tibble(rownames=c_(.data)$name, .name_repair = "minimal") %>%
        setNames(c(c_(.data)$name, sprintf("%s%s", prefix, selected_features_from_exp)))
    }
  }
  suppressMessages({
    lapply(selected_experiments_list, extract_feature_values) |> 
        purrr::reduce(full_join)
    })
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
#' @importFrom purrr reduce
#' @importFrom dplyr full_join
#' @importFrom SummarizedExperiment assays assayNames
#'
#' @param .data A tidySingleCellExperiment
#' @param features A character
#' @param all A boolean
#' @param exclude_zeros A boolean
#'
#' @return A tidySingleCellExperiment object
#'
#' @noRd
get_abundance_sc_long <- function(.data, features = NULL, all = FALSE, exclude_zeros = FALSE, ...) {
  
  arg_list <- c(mget(ls(environment(), sorted=F)), match.call(expand.dots=F)$...)
  assays_to_use <- eval(arg_list$assays)
  if(is.null(assays_to_use)) assays_to_use <- get_all_assays(.data)$assay_id
  
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
  
  # Get selected features and assays
  feature_df <- get_all_features(.data)
  selected_features <- feature_df[(feature_df$feature %in% features), ]
  selected_features <- selected_features[selected_features$assay_id %in% assays_to_use,]
  selected_experiments_list <- split(x = selected_features, f = as.character(selected_features$exp_id))
  
  extract_feature_values <- function(exp) {
    selected_features_exp <- as.character(unique(exp$exp_id))
    selected_features_assay <- as.character(unique(exp$assay_name))
    selected_features_assay_names <- as.character(unique(exp$assay_id))
    if (selected_features_exp == "Main") {
      assays(.data)[selected_features_assay] %>%
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
    } else {
      assays(altExps(.data)[[selected_features_exp]])[selected_features_assay] %>%
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
  suppressMessages({
    lapply(selected_experiments_list, extract_feature_values) |> 
        purrr::reduce(full_join)
    })
}

#' @importFrom dplyr select any_of
#' @importFrom S4Vectors DataFrame
#'
#' @keywords internal
#'
#' @param .data A tibble
#' @param SingleCellExperiment_object A `tidySingleCellExperiment`
#'
#' @noRd
as_meta_data <- function(.data, SingleCellExperiment_object) {

    col_to_exclude <-
        get_special_columns(SingleCellExperiment_object) |>
        # Need this in case we have multiple reduced dimensions 
        # with overlapping column names, e.g., multiple PCAs
        vctrs::vec_as_names(repair="unique") |>
        # To avoid name change by the 'bind_cols()' of 'as_tibble()'
        trick_to_avoid_renaming_of_already_unique_columns_by_dplyr()
    
    .data_df <- .data %>%
        select(-any_of(col_to_exclude)) %>%
        data.frame()
    
    # Set row names in a robust way; the 'row.names' argument 
    # of 'data.frame()' does not work for 1-row 'data.frame's
    sym <- c_(SingleCellExperiment_object)$symbol
    rownames(.data_df) <- pull(.data_df, !!sym)
    
    .data_df <- select(.data_df, -!!sym)
    return(DataFrame(.data_df))
}

#' @importFrom purrr map_chr
#'
#' @keywords internal
#'
#' @param SingleCellExperiment_object A `tidySingleCellExperiment`
#'
#' @noRd
get_special_columns <- function(SingleCellExperiment_object) {
  get_special_datasets(SingleCellExperiment_object) %>%
    map(~ .x %>% colnames()) %>%
    unlist() %>%
    as.character()
}

#' @importFrom SingleCellExperiment reducedDims
get_special_datasets <- function(SingleCellExperiment_object, n_dimensions_to_return=Inf) {
    
    rd <- reducedDims(SingleCellExperiment_object)

    map2(as.list(rd), names(rd), ~ {
        n_dims <- min(n_dimensions_to_return, ncol(.x))
        mat <- .x[, seq_len(n_dims), drop=FALSE]
        # Set names as SCE is much less constrained 
        # and there could be missing names
        if (is.null(colnames(mat))) colnames(mat) <- 
            sprintf("%s%s", .y, seq_len(ncol(mat)))
        return(mat)
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

#' returns variables from an expression
#' @param expression an expression
#' @importFrom rlang enexpr
#' @return list of symbols
return_arguments_of <- function(expression){
    variables <- enexpr(expression) |> as.list()
    if(length(variables) > 1) {
        variables <- variables[-1] # removes first element which is function
    }
    variables
}

#' @importFrom purrr when
#' @importFrom dplyr select
#' @importFrom rlang expr
#' @importFrom tidyselect eval_select
select_helper <- function(.data, ...) {
    loc <- tidyselect::eval_select(expr(c(...)), .data)
    dplyr::select(.data, loc)
}

data_frame_returned_message <- paste(
    "tidySingleCellExperiment says:",
    "A data frame is returned for independent data analysis.")
duplicated_cell_names <- paste(
    "tidySingleCellExperiment says:",
    "This operation lead to duplicated cell names.",
    "A data frame is returned for independent data analysis.")

# This function is used for the change of special sample column to .sample
# Check if "sample" is included in the query and is not part of any other existing annotation
#' @importFrom stringr str_detect
#' @importFrom stringr regex
  
is_sample_feature_deprecated_used <- function(.data, 
    user_columns, use_old_special_names=FALSE) {
    
    cell <- any(str_detect(user_columns, regex("\\bcell\\b")))
    .cell <- any(str_detect(user_columns, regex("\\W*(\\.cell)\\W*")))
    
    old_standard_is_used <- 
        !"cell" %in% colnames(colData(.data)) &&
        ("cell" %in% user_columns || (cell && !.cell))
    
    if (old_standard_is_used) {
        warning("tidySingleCellExperiment says:",
            " from version 1.3.1, the special columns including",
            " cell id (colnames(se)) has changed to \".cell\".",
            " This dataset is returned with the old-style vocabulary (cell),",
            " however, we suggest to update your workflow",
            " to reflect the new vocabulary (.cell).")
        use_old_special_names <- TRUE
    }
    use_old_special_names
}

get_special_column_name_symbol <- function(name) {
    list(name=name, symbol=as.symbol(name))
}

# Key column names
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-

ping_old_special_column_into_metadata <- function(.data) {
    metadata(.data)$cell__ <- get_special_column_name_symbol("cell")
    return(.data)
}

get_special_column_name_cell <- function(name) {
    list(name=name, symbol=as.symbol(name))
}


#' @importFrom S4Vectors metadata
c_ <- function(x) {
    # Check if old deprecated columns are used
    if ("cell__" %in% names(metadata(x)))
        cell__ <- metadata(x)$cell__
    return(cell__)
}

#' Add attribute to abject
#'
#' @keywords internal
#' @noRd
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr <- function(var, attribute, name) {
    attr(var, name) <- attribute
    return(var)
}

#' @importFrom tidyr spread
#' @importFrom tibble enframe
#' @importFrom purrr map reduce
special_datasets_to_tibble <- function(.singleCellExperiment, ...) {
    x <- .singleCellExperiment %>%
        get_special_datasets(...) %>%
        map(~ {
            if (!is.null(dim(.x)))
                return(as_tibble(.x))
            # If row == 1 do a trick
            .x %>%
                tibble::enframe() %>%
                tidyr::spread(name, value)
        }) %>% purrr::reduce(bind_cols)
    
    # To avoid name change by the 'bind_cols()' of 'as_tibble()'
    colnames(x) <- colnames(x) |> 
        trick_to_avoid_renaming_of_already_unique_columns_by_dplyr()
    return(x)
}

#' @importFrom stringr str_replace_all
trick_to_avoid_renaming_of_already_unique_columns_by_dplyr <- function(x) {
    str_replace_all(x, "\\.\\.\\.", "___")
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
get_specific_annotation_columns <- function(.data, .col) {
    
    # Comply with CRAN NOTES
    . <- NULL
    
    # Make col names
    .col <- enquo(.col)
    
    # x-annotation df
    n_x <- .data |> distinct_at(vars(!!.col)) |> nrow()
    
    # Exclude columns that have more values than my .col
    columns_unique_length = .data |> select(-!!.col) |> lapply(function(x) unique(x) |> length())
    columns_unique_length = columns_unique_length[columns_unique_length<=n_x]
    
    .sample = .data |> select(!!.col) |> unite(".sample", !!.col) |> pull(.sample)
    
    # element wise columns
    columns_unique_length |>
      names() |> 
        map(~ {
            n_.x <- .data |> pull(all_of(.x)) |> paste(.sample)  |> unique() |> length()
            if (n_.x == n_x) .x else NULL
        }) %>%
        # Drop NULL
        { (.)[lengths((.)) != 0] } |>
        unlist()
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

subset <- function(.data, .column)	{
    
    # Make col names
    .column <- enquo(.column)
    
    # Check if column present
    if (!all(quo_names(.column) %in% colnames(.data)))
        stop("nanny says: some of the .column specified",
            " do not exist in the input data frame.")
    
    .data |>
        # Selecting the right columns
        select(!!.column, get_specific_annotation_columns(.data, !!.column)) %>%
        distinct()
}


splitColData <- function(x, f) {
  # This is by @jma1991 
  # at https://github.com/drisso/SingleCellExperiment/issues/55
  
  i <- split(seq_along(f), f)
  
  v <- vector(mode = "list", length = length(i))
  
  names(v) <- names(i)
  
  for (n in names(i)) { v[[n]] <- x[, i[[n]]] }
  
  return(v)
  
}

cell__ <- get_special_column_name_symbol(".cell")
feature__ <- get_special_column_name_symbol(".feature")
sample__ <- get_special_column_name_symbol(".sample")
