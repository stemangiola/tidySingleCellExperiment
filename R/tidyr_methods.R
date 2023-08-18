#' @name unnest
#' @rdname unnest
#' @inherit tidyr::unnest
#' @aliases unnest_single_cell_experiment
#' @return `tidySingleCellExperiment`
#' 
#' @examples
#' pbmc_small |> 
#'     nest(data=-groups) |> 
#'     unnest(data)
#'
#' @importFrom tidyr unnest
#' @importFrom purrr when
#' @export
unnest.tidySingleCellExperiment_nested <- function(data, cols, ..., 
    keep_empty=FALSE, ptype=NULL, names_sep=NULL, names_repair="check_unique", 
    .drop, .id, .sep, .preserve) {
    
    cols <- enquo(cols)
    
    unnest_single_cell_experiment(data, !!cols, ..., 
        keep_empty=keep_empty, ptype=ptype,
        names_sep=names_sep, names_repair=names_repair)
}

#' @rdname unnest
#' @importFrom methods is
#' @importFrom tidyr unnest
#' @importFrom rlang quo_name 
#' @importFrom rlang enquo 
#' @importFrom purrr reduce
#' @importFrom purrr when
#' @importFrom purrr imap
#' @export
unnest_single_cell_experiment  <-  function(data, cols, ..., 
    keep_empty=FALSE, ptype=NULL, names_sep=NULL, names_repair="check_unique", 
    .drop, .id, .sep, .preserve) {
    
    # Need this otherwise crashes map
    .data_ <- data
    cols <- enquo(cols)
    
    # If my only column to unnest() is a 'tidySingleCellExperiment'
    # [HLC: comment says 'only', but only the first entry is being checked.
    # is this intentional? or, what happens if, e.g., the 2nd is a tidySCE?]
    .test <- .data_ |> pull(!!cols) |> _[[1]] |> is("SingleCellExperiment")
    if (.test) {
        # Do my trick to unnest()
        .data_ |>
            mutate(!!cols := imap(
                !!cols, ~ .x |>
                    bind_cols_(
                        # Attach back the columns used for nesting
                        .data_ |> 
                            select(-!!cols) |>
                            slice(rep(.y, nrow(as_tibble(.x))))
                    )
            )) |>
            pull(!!cols) |>
            reduce(bind_rows)
    } else {
        # Else do normal stuff
        .data_ |>
            drop_class("tidySingleCellExperiment_nested") |>
            tidyr::unnest(!!cols, ..., keep_empty=keep_empty, 
                ptype=ptype, names_sep=names_sep, names_repair=names_repair) |>
            add_class("tidySingleCellExperiment_nested")
    }
}

#' @name nest
#' @rdname nest
#' @inherit tidyr::nest
#' @return `tidySingleCellExperiment_nested`
#'
#' @examples
#' pbmc_small |> 
#'     nest(data=-groups) |> 
#'     unnest(data)
#'     
#' @importFrom tidyr nest
#' @importFrom rlang enquos
#' @importFrom rlang :=
#' @export
nest.SingleCellExperiment <- function(.data, ..., .names_sep = NULL) {
    cols <- enquos(...)
    col_name_data <- names(cols)
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>%
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    my_data__ <- .data
    cols_sym <- as.symbol(col_name_data)
    cell_sym <- c_(my_data__)$symbol
    
    my_data__ %>%
        # This is needed otherwise nest goes into loop and fails
        to_tib() %>%
        tidyr::nest(...) %>%
        mutate(
            !!cols_sym := map(
                !!cols_sym, ~ {
                    my_data__ %>%
                        # Subset cells
                        filter(!!cell_sym %in% pull(.x, !!cell_sym)) %>%
                        # Subset columns
                        select(colnames(.x))
                }
            )
        ) %>%
        # Coerce to tidySingleCellExperiment_nested for unnesting
        add_class("tidySingleCellExperiment_nested")
}

#' @name extract
#' @rdname extract
#' @inherit tidyr::extract
#' @return `tidySingleCellExperiment`
#' 
#' @examples
#' pbmc_small|>
#'   extract(groups, 
#'     into="g", 
#'     regex="g([0-9])", 
#'     convert=TRUE)
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom tidyr extract
#' @export
extract.SingleCellExperiment <- function(data, col, into, 
    regex="([[:alnum:]]+)", remove=TRUE, convert=FALSE, ...) {
    col <- enquo(col)
    
    # Deprecation of special column names
    .cols <- c(quo_name(col), into)
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }
    
    colData(data) <-
        data %>%
        as_tibble() %>%
        tidyr::extract(col=!!col, into=into, 
            regex=regex, remove=remove, convert=convert, ...) %>%
        as_meta_data(data)
    
    data
}

#' @name pivot_longer
#' @rdname pivot_longer
#' @inherit tidyr::pivot_longer
#' @return `tidySingleCellExperiment`
#' 
#' @export
#' @examples
#' # See vignette("pivot") for examples and explanation
#' pbmc_small |> pivot_longer(
#'   cols=c(orig.ident, groups),
#'   names_to="name", values_to="value")
#' 
#' @importFrom ellipsis check_dots_used
#' @importFrom tidyr pivot_longer
#' @export
pivot_longer.SingleCellExperiment <- function(data,
    cols, ..., cols_vary = "fastest", names_to = "name",
    names_prefix = NULL, names_sep = NULL, names_pattern = NULL,
    names_ptypes = NULL, names_transform = NULL, names_repair = "check_unique",
    values_to = "value", values_drop_na = FALSE, values_ptypes = NULL,
    values_transform = NULL) {
    cols <- enquo(cols)
    
    message(data_frame_returned_message)
    
    # Deprecation of special column names
    .cols <- c(quo_names(cols))
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }
    
    data %>%
        as_tibble() %>%
        tidyr::pivot_longer(!!cols,
            ...,
            cols_vary = cols_vary,
            names_to = names_to,
            names_prefix = names_prefix,
            names_sep = names_sep,
            names_pattern = names_pattern,
            names_ptypes = names_ptypes,
            names_transform = names_transform,
            names_repair = names_repair,
            values_to = values_to,
            values_drop_na = values_drop_na,
            values_ptypes = values_ptypes,
            values_transform = values_transform)
}

#' @name unite
#' @rdname unite
#' @inherit tidyr::unite
#' @return `tidySingleCellExperiment`
#' 
#' @examples
#' pbmc_small |> unite(
#'   col="new_col", 
#'   c(orig.ident, groups))
#'     
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom ellipsis check_dots_unnamed
#' @importFrom tidyr unite
#' @export
unite.SingleCellExperiment <- function(data, col, 
    ..., sep="_", remove=TRUE, na.rm=FALSE) {
    
    # Check that we are not modifying a key column
    cols <- enquo(col)
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }
    
    .view_only_cols <- c(
        get_special_columns(data),
        get_needed_columns(data))
    
    .test <- intersect(
        quo_names(cols), 
        .view_only_cols)
    
    if (remove && length(.test)) {
        stop("tidySingleCellExperiment says:",
            " you are trying to rename a column",
            " that is view only ", 
            paste(.view_only_cols, collapse=", "),
            " (it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one.")
    }
    
    colData(data) <- data %>%
        as_tibble() %>%
        tidyr::unite(!!cols, ..., sep=sep, remove=remove, na.rm=na.rm) %>%
        as_meta_data(data)
    
    data
}

#' @name separate
#' @rdname separate
#' @inherit tidyr::separate
#' @return `tidySingleCellExperiment`
#' 
#' @examples
#' un <- pbmc_small |> unite("new_col", c(orig.ident, groups))
#' un |> separate(new_col, c("orig.ident", "groups"))
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom ellipsis check_dots_used
#' @importFrom tidyr separate
#' @export
separate.SingleCellExperiment <- function(data, col, into, 
    sep="[^[:alnum:]]+", remove=TRUE, convert=FALSE, 
    extra="warn", fill="warn", ...) {
    
    # Check that we are not modifying a key column
    cols <- enquo(col)
    
    # Deprecation of special column names
    .cols <- c(quo_names(cols))
    if (is_sample_feature_deprecated_used(data, .cols)) {
        data <- ping_old_special_column_into_metadata(data)
    }
    
    .view_only_cols <- c(
        get_special_columns(data),
        get_needed_columns(data))
    
    .test <- intersect(
        quo_names(cols), 
        .view_only_cols)
    
    if (remove && length(.test)) {
        stop("tidySingleCellExperiment says:",
            " you are trying to rename a column",
            " that is view only ",
            paste(.view_only_cols, collapse=", "),
            "(it is not present in the colData).",
            " If you want to mutate a view-only column,",
            " make a copy and mutate that one.")
    }
    
    colData(data) <-
        data %>%
        as_tibble() %>%
        tidyr::separate(
            !!cols, into=into, sep=sep, remove=remove, 
            convert=convert, extra=extra, fill=fill, ...) %>%
        as_meta_data(data)
    
    data
}
