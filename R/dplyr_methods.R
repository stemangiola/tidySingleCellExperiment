#' @name arrange
#' @rdname arrange
#' @inherit dplyr::arrange
#' @family single table verbs
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> 
#'     arrange(nFeature_RNA)
#'     
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @export
arrange.SingleCellExperiment <- function(.data, ..., .by_group=FALSE) {
    new_metadata <- 
        .data |> 
        as_tibble() |> 
        dplyr::arrange(..., .by_group=.by_group)
    
    .data[, pull(new_metadata, !!c_(.data)$symbol)]
}

#' @name bind_rows
#' @rdname bind_rows
#' @inherit ttservice::bind_rows
#' 
#' @examples
#' data(pbmc_small)
#' tt <- pbmc_small
#' bind_rows(tt, tt)
#'
#' tt_bind <- tt |> select(nCount_RNA, nFeature_RNA)
#' tt |> bind_cols(tt_bind)
#' 
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom rlang dots_values
#' @importFrom ttservice bind_rows
#' @importFrom SingleCellExperiment cbind
#' @export
bind_rows.SingleCellExperiment <- function(..., .id=NULL, add.cell.ids=NULL) {
    tts <- flatten_if(dots_values(...), is_spliced)
    
    new_obj <- SingleCellExperiment::cbind(tts[[1]], tts[[2]])
    
    # If duplicated cell names
    if (new_obj %>% colnames %>% duplicated %>% which %>% length %>% gt(0))
        warning("tidySingleCellExperiment says:",
            " you have duplicated cell names;",
            " they will be made unique.")
    colnames(new_obj) <- make.unique(colnames(new_obj), sep="_")
    
    new_obj
}

#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' @importFrom rlang dots_values
#' @importFrom ttservice bind_cols
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
bind_cols_ <- function(..., .id=NULL) {
    tts <- tts <- flatten_if(dots_values(...), is_spliced)
    
    colData(tts[[1]]) <- bind_cols(colData(tts[[1]]) %>% as.data.frame(),
        tts[[2]], .id=.id) %>% DataFrame()
    
    tts[[1]]
}

#' @rdname bind_rows
#' @aliases bind_cols
#' @export
bind_cols.SingleCellExperiment <- bind_cols_

#' @name distinct
#' @rdname distinct
#' @inherit dplyr::distinct
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> distinct(groups)
#'
#' @importFrom dplyr distinct
#' @export
distinct.SingleCellExperiment <- function(.data, ..., .keep_all=FALSE) {
    message(data_frame_returned_message)
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    .data %>%
        as_tibble() %>%
        dplyr::distinct(..., .keep_all=.keep_all)
}

#' @name filter
#' @rdname filter
#' @inherit dplyr::filter
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> filter(groups == "g1")
#'
#' # Learn more in ?dplyr_tidy_eval
#' 
#' @importFrom purrr map
#' @importFrom dplyr filter
#' @export
filter.SingleCellExperiment <- function(.data, ..., .preserve=FALSE) {
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    new_meta <- .data |>
        as_tibble() |>
        dplyr::filter(..., .preserve=.preserve)
    
    # Try to solve missing colnames
    .cell <- c_(.data)$symbol
    if (colnames(.data) |> is.null()) {
        message("tidySingleCellExperiment says: ",
            "the input object does not have cell names (colnames(...)).\n",
            "Therefore, the cell column in the filtered tibble abstraction ",
            "will still include an incremental integer vector.")
        new_meta <- new_meta %>% mutate(!!.cell := as.integer(!!.cell))
    }
    
    .data[, pull(new_meta, !!.cell)]
}

#' @name group_by
#' @rdname group_by
#' @inherit dplyr::group_by
#' @seealso \code{}
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> group_by(groups)
#'     
#' @importFrom dplyr group_by_drop_default
#' @importFrom dplyr group_by
#' @export
group_by.SingleCellExperiment <- function(.data, ..., 
    .add=FALSE, .drop=group_by_drop_default(.data)) {
    
    message(data_frame_returned_message)
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    .data %>%
        as_tibble() %>%
        dplyr::group_by(..., .add=.add, .drop=.drop)
}


#' @name summarise
#' @aliases summarize
#' @inherit dplyr::summarise
#' @family single table verbs
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> summarise(mean(nCount_RNA))
#'
#' @importFrom dplyr summarise
#' @importFrom purrr map
#' @export
summarise.SingleCellExperiment <- function(.data, ...) {
    message(data_frame_returned_message)
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    .data %>%
        as_tibble() %>%
        dplyr::summarise(...)
}

#' @name summarise
#' @rdname summarise
#' @importFrom dplyr summarize
#' @export
summarize.SingleCellExperiment <- summarise.SingleCellExperiment

#' @name mutate
#' @rdname mutate
#' @inherit dplyr::mutate
#' @family single table verbs
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> mutate(nFeature_RNA=1)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom rlang enquos
#' @importFrom dplyr mutate
#' @importFrom purrr map
#' @export
mutate.SingleCellExperiment <- function(.data, ...) {
    
    # Check that we are not modifying a key column
    cols <- enquos(...) %>% names()
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    .view_only_cols <- c(
        get_special_columns(.data),
        get_needed_columns(.data))
    
    .test <- cols |>
        intersect(.view_only_cols) |>
        length()
    
    if (.test) {
        stop("tidySingleCellExperiment says:",
            " you are trying to mutate a column that is view only",
            " ", paste(.view_only_cols, collapse=", "),
            " (it is not present in the colData).",
            " If you want to mutate a view-only column, make a copy",
            " (e.g. mutate(new_column=", cols[1], ")) and mutate that one.")
    }
    
    colData(.data) <-
        .data %>%
        as_tibble() %>%
        dplyr::mutate(...) %>%
        as_meta_data(.data)
    
    .data
}

#' @name rename
#' @rdname rename
#' @inherit dplyr::rename
#' @family single table verbs
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> rename(s_score=nFeature_RNA)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom tidyselect eval_select
#' @importFrom dplyr rename
#' @export
rename.SingleCellExperiment <- function(.data, ...) {
    
    # Check that we are not modifying a key column
    read_only_columns <- c(
        get_needed_columns(.data),
        get_special_columns(.data))

    # Small df to be more efficient
    df <- .data[1, 1] |> as_tibble() 
    
    # What columns we are going to create
    cols_from <- tidyselect::eval_select(expr(c(...)), df) |> names()
    
    # What are the columns before renaming
    original_columns <- df |> colnames()
    
    # What the column after renaming would be
    new_colums <- df |> rename(...) |> colnames()
    
    # What column you are impacting
    changed_columns <- original_columns |> setdiff(new_colums)
    
    # Check that you are not impacting any read-only columns
    if (any(changed_columns %in% read_only_columns)) {
        stop("tidySingleCellExperiment says:",
            " you are trying to rename a column that is view only",
            " ", paste(changed_columns, collapse=", "),
            " (it is not present in the colData).",
            " If you want to rename a view-only column, make a copy",
            " (e.g., mutate(", cols_from[1], "=",  changed_columns[1], ")).")
    }
    
    colData(.data) <- 
        colData(.data) |>
        as.data.frame() |>
        dplyr::rename(...) |>
        DataFrame()

    .data
}

#' @name rowwise
#' @rdname rowwise
#' @inherit dplyr::rowwise
#'
#' @examples
#' # TODO
#'
#' @importFrom dplyr rowwise
#' @export
rowwise.SingleCellExperiment <- function(data, ...) {
    message(data_frame_returned_message)
    
    data %>%
        as_tibble() %>%
        dplyr::rowwise(...)
}

#' @name left_join
#' @rdname left_join
#' @inherit dplyr::left_join
#'
#' @examples
#' data(pbmc_small)
#' tt <- pbmc_small
#' tt |> left_join(tt |>  
#'   distinct(groups) |> 
#'   mutate(new_column=1:2))
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr left_join
#' @importFrom dplyr count
#' @export
left_join.SingleCellExperiment <- function(x, y, 
    by=NULL, copy=FALSE, suffix=c(".x", ".y"), ...) {
    
    # Deprecation of special column names
    .cols <- if (!is.null(by)) by else colnames(y)
    if (is_sample_feature_deprecated_used(x, .cols)) {
        x <- ping_old_special_column_into_metadata(x)
    }
    
    z <- x |>
        as_tibble() |>
        dplyr::left_join(y, by=by, copy=copy, suffix=suffix, ...)
    
    # If duplicated cells returns tibble
    if (any(duplicated(z[[c_(x)$name]]))) {
        message(duplicated_cell_names)
        return(z)
    }
    
    # Otherwise return updated tidySingleCellExperiment
    colData(x) <- z |> as_meta_data(x)
    return(x)
}

#' @name inner_join
#' @rdname inner_join
#' @inherit dplyr::inner_join
#'
#' @examples
#' data(pbmc_small)
#' tt <- pbmc_small
#' tt |> inner_join(tt |> 
#'   distinct(groups) |>  
#'   mutate(new_column=1:2) |> 
#'   slice(1))
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr inner_join
#' @importFrom dplyr pull
#' @export
inner_join.SingleCellExperiment <- function(x, y, 
    by=NULL, copy=FALSE, suffix=c(".x", ".y"), ...) {
    
    # Deprecation of special column names
    .cols <- if (!is.null(by)) by else colnames(y)
    if (is_sample_feature_deprecated_used(x, .cols)) {
        x <- ping_old_special_column_into_metadata(x)
    }
    
    z <- x |>
        as_tibble() |>
        dplyr::inner_join(y, by=by, copy=copy, suffix=suffix, ...)
    
    # If duplicated cells returns tibble
    if (any(duplicated(z[[c_(x)$name]]))) {
        message(duplicated_cell_names)
        return(z)
    }
    
    # Otherwise return updated tidySingleCellExperiment
    new_obj <- x[, pull(z, c_(x)$name)]
    colData(new_obj) <- z |> as_meta_data(new_obj)
    return(new_obj)
}

#' @name right_join
#' @rdname right_join
#' @inherit dplyr::right_join
#'
#' @examples
#' data(pbmc_small)
#' tt <- pbmc_small
#' tt |> right_join(tt |> 
#'   distinct(groups) |> 
#'   mutate(new_column=1:2) |> 
#'   slice(1))
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr right_join
#' @importFrom dplyr pull
#' @export
right_join.SingleCellExperiment <- function(x, y, 
    by=NULL, copy=FALSE, suffix=c(".x", ".y"), ...) {
    
    # Deprecation of special column names
    .cols <- if (!is.null(by)) by else colnames(y)
    if (is_sample_feature_deprecated_used(x, .cols)) {
        x <- ping_old_special_column_into_metadata(x)
    }
    
    z <- x |>
        as_tibble() |>
        dplyr::right_join(y, by=by, copy=copy, suffix=suffix, ...)
    
    # If duplicated cells returns tibble
    if (any(duplicated(z[[c_(x)$name]]))) {
        message(duplicated_cell_names)
        return(z)
    } 
    
    # Otherwise return updated tidySingleCellExperiment
    new_obj <- x[, pull(z, c_(x)$name)]
    colData(new_obj) <- z |> as_meta_data(new_obj)
    return(new_obj)
}

#' @name full_join
#' @rdname full_join
#' @inherit dplyr::full_join
#'
#' @examples
#' data(pbmc_small)
#' tt <- pbmc_small
#' tt |> full_join(tibble::tibble(groups="g1", other=1:4))
#'
#' @importFrom dplyr full_join
#' @importFrom dplyr pull
#' @export
full_join.SingleCellExperiment <- function(x, y, 
    by=NULL, copy=FALSE, suffix=c(".x", ".y"), ...) {
    
    # Deprecation of special column names
    .cols <- if (!is.null(by)) by else colnames(y)
    if (is_sample_feature_deprecated_used(x, .cols)) {
        x <- ping_old_special_column_into_metadata(x)
    }
    
    z <- x |>
        as_tibble() |>
        dplyr::full_join(y, by=by, copy=copy, suffix=suffix, ...) 
    
    # If duplicated cells returns tibble
    if (any(duplicated(z[[c_(x)$name]]))) {
        message(duplicated_cell_names)
        return(z)
    }
       
    # Otherwise return updated tidySingleCellExperiment
    new_obj <- x[, pull(z, c_(x)$name)]
    colData(new_obj) <- z |> as_meta_data(x)
    return(new_obj)
}

#' @name slice
#' @rdname slice
#' @aliases slice_head slice_tail 
#'   slice_sample slice_min slice_max
#' @inherit dplyr::slice
#' @family single table verbs
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> slice(1)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr slice
#' @export
slice.SingleCellExperiment <- function(.data, ..., .by=NULL, .preserve=FALSE) {
    new_meta <- .data |>
        colData() |>
        as.data.frame() |>
        dplyr::slice(..., .by=.by, .preserve=.preserve)
    
    .data[, rownames(new_meta)]
}

#' @name slice_sample
#' @rdname slice
#' @inherit dplyr::slice_sample
#' @examples
#' data(pbmc_small)
#' pbmc_small |> slice_sample(n=1)
#' pbmc_small |> slice_sample(prop=0.1)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr slice_sample
#' @export
slice_sample.SingleCellExperiment <- function(.data, ..., n=NULL,
    prop=NULL, by=NULL, weight_by=NULL, replace=FALSE) {
    lifecycle::signal_superseded("1.0.0", "sample_n()", "slice_sample()")

    if (!is.null(n))
        new_meta <-
            .data |>
            colData() |>
            as_tibble(rownames=c_(.data)$name) |>
            select(-everything(), c_(.data)$name, {{ by }}, {{ weight_by }}) |>
            slice_sample(..., n=n, by={{ by }},
                weight_by={{ weight_by }}, replace=replace)
    else if (!is.null(prop))
        new_meta <-
            .data |>
            colData() |>
            as_tibble(rownames=c_(.data)$name) |>
            select(-everything(), c_(.data)$name, {{ by }}, {{ weight_by }}) |>
            slice_sample(..., prop=prop, by={{ by }},
                weight_by={{ weight_by }}, replace=replace)
    else
        stop("tidySingleCellExperiment says:",
            " you should provide `n` or `prop` arguments")

    count_cells <- new_meta %>%
        select(!!c_(.data)$symbol) %>%
        count(!!c_(.data)$symbol)

    .max_cell_count <- ifelse(nrow(count_cells)==0, 0, max(count_cells$n))

    # If repeated cells due to replacement
    if (.max_cell_count |> gt(1)){
        message("tidySingleCellExperiment says: When sampling with replacement",
            " a data frame is returned for independent data analysis.")
        .data |>
            as_tibble()  |>
            right_join(new_meta %>% 
                select(!!c_(.data)$symbol), by=c_(.data)$name)
    } else {
        .data[, pull(new_meta, !!c_(.data)$symbol)]
    }
}

#' @name slice_head
#' @rdname slice
#' @inherit dplyr::slice_head
#' @examples
#' data(pbmc_small)
#' # First rows based on existing order
#' pbmc_small |> slice_head(n=5)
#' 
#' @importFrom dplyr slice_head
#' @importFrom tibble rowid_to_column
#' @export
slice_head.SingleCellExperiment <- function(.data, ..., n, prop, by=NULL) {
    row_number___ <- NULL
    idx <- .data |>
        colData() |>
        as.data.frame() |>
        select(-everything(), {{ by }}) |>
        rowid_to_column(var='row_number___')  |>
        slice_head(..., n=n, prop=prop, by={{ by }}) |>
        pull(row_number___)

    new_obj <- .data[, idx]
    new_obj
}

#' @name slice_tail
#' @rdname slice
#' @inherit dplyr::slice_tail
#' @examples
#' data(pbmc_small)
#' # First rows based on existing order
#' pbmc_small |> slice_tail(n=5)
#' 
#' @importFrom dplyr slice_tail
#' @importFrom tibble rowid_to_column
#' @export
slice_tail.SingleCellExperiment <- function(.data, ..., n, prop, by=NULL) {
    row_number___ <- NULL
    idx <- .data |>
        colData() |>
        as.data.frame() |>
        select(-everything(), {{ by }}) |>
        rowid_to_column(var='row_number___')  |>
        slice_tail(..., n=n, prop=prop, by={{ by }}) |>
        pull(row_number___)

    new_obj <- .data[, idx]
    new_obj
}

#' @name slice_min
#' @rdname slice
#' @inherit dplyr::slice_min
#' @examples
#' data(pbmc_small)
#' 
#' # Rows with minimum and maximum values of a metadata variable
#' pbmc_small |> slice_min(nFeature_RNA, n=5)
#'
#' # slice_min() and slice_max() may return more rows than requested
#' # in the presence of ties.
#' pbmc_small |>  slice_min(nFeature_RNA, n=2)
#'
#' # Use with_ties=FALSE to return exactly n matches
#' pbmc_small |> slice_min(nFeature_RNA, n=2, with_ties=FALSE)
#'
#' # Or use additional variables to break the tie:
#' pbmc_small |> slice_min(tibble::tibble(nFeature_RNA, nCount_RNA), n=2)
#'
#' # Use by for group-wise operations
#' pbmc_small |> slice_min(nFeature_RNA, n=5, by=groups)
#'
#' @importFrom dplyr slice_min
#' @importFrom tibble rowid_to_column
#' @export
slice_min.SingleCellExperiment <- function(.data, order_by, ..., n, prop,
    by=NULL, with_ties=TRUE, na_rm=FALSE) {
    row_number___ <- NULL
    order_by_variables <- return_arguments_of(!!enexpr(order_by))

    idx <- .data |>
        colData() |>
        as.data.frame() |>
        select(-everything(), !!!order_by_variables, {{ by }}) |>
        rowid_to_column(var ='row_number___')  |>
        slice_min(
            order_by={{ order_by }}, ..., n=n, prop=prop, by={{ by }},
            with_ties=with_ties, na_rm=na_rm
        ) |>
        pull(row_number___)

    new_obj <- .data[, idx]
    new_obj
}

#' @name slice_max
#' @rdname slice
#' @inherit dplyr::slice_max
#' @examples
#' data(pbmc_small)
#' # Rows with minimum and maximum values of a metadata variable
#' pbmc_small |> slice_max(nFeature_RNA, n=5)
#' 
#' @importFrom dplyr slice_max
#' @importFrom tibble rowid_to_column
#' @export
slice_max.SingleCellExperiment <- function(.data, order_by, ..., n, prop,
    by=NULL, with_ties=TRUE, na_rm=FALSE) {
    row_number___ <- NULL

    order_by_variables <- return_arguments_of(!!enexpr(order_by))

    idx <- .data |>
        colData() |>
        as.data.frame() |>
        select(-everything(), !!!order_by_variables, {{ by }}) |>
        rowid_to_column(var ='row_number___')  |>
        slice_max(
            order_by={{ order_by }}, ..., n=n, prop=prop, by={{ by }},
            with_ties=with_ties, na_rm=na_rm
        ) |>
        pull(row_number___)

    new_obj <- .data[, idx]
    new_obj
}


#' @name select
#' @rdname select
#' @inherit dplyr::select
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> select(cell, orig.ident)
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr select
#' @export
select.SingleCellExperiment <- function(.data, ...) {
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    
    new_obj <- .data |>
        as_tibble() |>
        select_helper(...)
    
    # If key columns are missing, return tibble
    if (!all(get_needed_columns(.data) %in% colnames(new_obj))) {
        message(
            "tidySingleCellExperiment says: Key columns are missing.",
            " A data frame is returned for independent data analysis.")
        return(new_obj)
    }
    
    # Otherwise return updated tidySingleCellExperiment
    colData(.data) <- new_obj |> as_meta_data(.data)
    return(.data)
}

#' @name sample_n
#' @rdname sample_n
#' @aliases sample_frac
#' @inherit dplyr::sample_n
#' @return `tidySingleCellExperiment`
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> sample_n(50)
#' pbmc_small |> sample_frac(0.1)
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr sample_n
#' @export
sample_n.SingleCellExperiment <- function(tbl, size, 
    replace=FALSE, weight=NULL, .env=NULL, ...) {
    
    lifecycle::signal_superseded("1.0.0", "sample_n()", "slice_sample()")
    
    new_meta <- colData(tbl) %>%
        as.data.frame() %>%
        as_tibble(rownames=c_(tbl)$name) %>%
        dplyr::sample_n( size, replace=replace, weight=weight, .env=.env, ...)
    
    count_cells <- new_meta %>% 
        select(!!c_(tbl)$symbol) %>% 
        count(!!c_(tbl)$symbol)
    
    # If repeted cells
    if (count_cells$n %>% max() %>% gt(1)) {
        message("tidySingleCellExperiment says:",
            " When sampling with replacement a data frame",
            " is returned for independent data analysis.")
        tbl %>%
            as_tibble() %>%
            right_join(new_meta %>% select(!!c_(tbl)$symbol), by=c_(tbl)$name)
    } else {
        new_obj <- tbl[, new_meta %>% pull(!!c_(tbl)$symbol)]
        new_obj
    }
}

#' @rdname sample_n
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr sample_frac
#' @export
sample_frac.SingleCellExperiment <- function(tbl, size=1, 
    replace=FALSE, weight=NULL, .env=NULL, ...) {
    
    lifecycle::signal_superseded("1.0.0", "sample_frac()", "slice_sample()")
    
    new_meta <- colData(tbl) %>%
        as.data.frame() %>%
        as_tibble(rownames=c_(tbl)$name) %>%
        dplyr::sample_frac(size, replace=replace, weight=weight, .env=.env, ...)
    
    count_cells <- new_meta %>% 
        select(!!c_(tbl)$symbol) %>% 
        count(!!c_(tbl)$symbol)
    
    # If repeated cells
    if (count_cells$n %>% max() %>% gt(1)) {
        message("tidySingleCellExperiment says:",
            " When sampling with replacement a data frame",
            " is returned for independent data analysis.")
        tbl %>%
            as_tibble() %>%
            right_join(new_meta %>% select(!!c_(tbl)$symbol), by=c_(tbl)$name)
    } else {
        new_obj <- tbl[, new_meta %>% pull(!!c_(tbl)$symbol)]
        new_obj
    }
}

#' @name count
#' @rdname count
#' @inherit dplyr::count
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> count(groups)
#'     
#' @importFrom dplyr count
#' @export
count.SingleCellExperiment <- function(x, ..., 
    wt=NULL, sort=FALSE, name=NULL, 
    .drop=group_by_drop_default(x)) {
    
    message(data_frame_returned_message)
    
    # Deprecation of special column names
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(x, .cols)) {
        x <- ping_old_special_column_into_metadata(x)
    }
    
    x |>
        as_tibble() |>
        dplyr::count(..., wt=!!enquo(wt), sort=sort, name=name, .drop=.drop)
}

#' @rdname count
#' @aliases add_count
#' @importFrom dplyr add_count
#' @export
add_count.SingleCellExperiment <- function(x, ..., 
    wt=NULL, sort=FALSE, name=NULL) {
    
    # Deprecation of special column names
    .cols <- enquos(..., .ignore_empty="all") %>% 
        map(~ quo_name(.x)) %>% unlist()
    if (is_sample_feature_deprecated_used(x, .cols)) {
        x <- ping_old_special_column_into_metadata(x)
    }
    
    colData(x) <- x |>
        as_tibble() |>
        dplyr::add_count(..., wt=!!enquo(wt), sort=sort, name=name) |>
        as_meta_data(x)
    
    x
}

#' @name pull
#' @rdname pull
#' @inherit dplyr::pull
#' 
#' @examples
#' data(pbmc_small)
#' pbmc_small |> pull(groups)
#'     
#' @importFrom ellipsis check_dots_used
#' @importFrom dplyr pull
#' @export
pull.SingleCellExperiment <- function(.data, var=-1, name=NULL, ...) {
    var <- enquo(var)
    name <- enquo(name)
    
    # Deprecation of special column names
    .cols <- quo_name(var)
    if (is_sample_feature_deprecated_used(.data, .cols)) {
        .data <- ping_old_special_column_into_metadata(.data)
    }
    .data %>%
        as_tibble() %>%
        dplyr::pull(var=!!var, name=!!name, ...)
}

#' @name group_split
#' @rdname group_split
#' @inherit dplyr::group_split title description details sections source
#'
#' @param .data A `tidySingleCellExperiment` object
#' @param ... The grouping variables by which to split the object
#'
#' @export
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small |> group_split(pbmc_small, groups)
group_split.SingleCellExperiment <- function(.data, ...) {
  
  var_list <- enquos(...)
  
  group_df <- .data |> 
      as_tibble()
  
  group_df <- group_df |> 
      unite("group_column___", !!!var_list, remove = FALSE)
    
  group_list <- group_df$group_column___
  
  groups <- group_list |> 
      unique()
  
  v <- vector(mode = "list", length = length(groups))
  
  for (i in seq_along(v)) {
      v[[i]] <- .data[,group_list == groups[[i]]]
  }
  
  v
  
}

