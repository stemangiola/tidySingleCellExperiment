#' @importFrom methods getMethod
setMethod(
    f="show",
    signature="SingleCellExperiment",
    definition=function(object) {
        opt <- getOption("restore_SingleCellExperiment_show", default=FALSE)
        if (isTRUE(opt)) {
            f <- getMethod(
                f="show",
                signature="SummarizedExperiment",
                where=asNamespace(ns="SummarizedExperiment"))
            f(object=object)
        } else { print(object) }
    }
)

setClass("tidySingleCellExperiment", contains="SingleCellExperiment")

#' @name join_features
#' @rdname join_features
#' @inherit ttservice::join_features
#' @aliases join_features,SingleCellExperiment-method
#'
#' @return A `tidySingleCellExperiment` object
#'   containing information for the specified features.
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small %>% join_features(
#'   features=c("HLA-DRA", "LYZ"))
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr contains
#' @importFrom dplyr everything
#' @importFrom ttservice join_features
#' @importFrom stringr str_c
#' @importFrom stringr str_subset
#' @export
setMethod("join_features", "SingleCellExperiment", function(.data,
    features=NULL, all=FALSE, exclude_zeros=FALSE, shape="long", ...) {
    # CRAN Note
    .cell <- NULL
    .feature <- NULL
    arg_list <- c(mget(ls(environment(), sorted=F)), match.call(expand.dots=F)$...)
    all_assays <- get_all_assays(.data)$assay_id
    if(is.null(arg_list$assays)) assays_from_join_call <- all_assays
    # Shape is long
    if (shape == "long") {
        # Suppress generic data frame creation message produced by left_join
        suppressMessages({
          .data <-
            .data %>%
            left_join(
              by=c_(.data)$name,
              get_abundance_sc_long(
                .data=.data,
                features=features,
                all=all,
                exclude_zeros=exclude_zeros, 
                ...)) %>%
            select(!!c_(.data)$symbol, .feature,
                   contains(".abundance"), everything())
        })
      
        # Provide data frame creation and abundance column message
        if (any(class(.data) == "tbl_df")) {
            
            abundance_columns <-
                .data %>%
                colnames() %>%
                stringr::str_subset('.abundance_')
            
            message(stringr::str_c("tidySingleCellExperiment says: join_features produces",
                " duplicate cell names to accomadate the long data format. For this reason, a data", 
                " frame is returned for independent data analysis. Assay feature abundance is", 
                " appended as ", 
                stringr::str_flatten_comma(abundance_columns, last = " and "), "."
            ))
        }
      
        .data
        
    # Shape if wide
    } else if (shape == "wide"){
        if(is.null(arg_list$assays)) stop("Please provide assays")
        .data  %>%
            left_join(
                by=c_(.data)$name,
                get_abundance_sc_wide(
                    .data=.data,
                    features=features,
                    all=all, 
                    ...))
    }
})

#' @name tidy
#' @rdname tidy
#' @title tidy for `SingleCellExperiment`
#'
#' @param object A `SingleCellExperiment` object.
#' @return A `tidySingleCellExperiment` object.
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small
#'
#' @export
tidy <- function(object) {
    UseMethod("tidy", object)
}

#' @rdname tidy
#' @importFrom lifecycle deprecate_warn
#' @export
tidy.SingleCellExperiment <- function(object) {

    # DEPRECATE
    deprecate_warn(
        when="1.1.1",
        what="tidy()",
        details="tidySingleCellExperiment says: tidy() is not needed anymore.")

    return(object)
}

#' @name aggregate_cells
#' @rdname aggregate_cells
#' @inherit ttservice::aggregate_cells
#' @aliases aggregate_cells,SingleCellExperiment-method
#'
#' @examples
#' data(pbmc_small)
#' pbmc_small_pseudo_bulk <- pbmc_small |>
#'   aggregate_cells(c(groups, ident), assays="counts")
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#' @importFrom tibble enframe
#' @importFrom Matrix rowSums
#' @importFrom ttservice aggregate_cells
#' @importFrom SummarizedExperiment assays assays<- assayNames
#' @importFrom S4Vectors split
#' @importFrom stringr str_remove
#' @importFrom dplyr full_join
#' @importFrom dplyr left_join
#' @importFrom dplyr group_by
#' @importFrom dplyr pick
#' @importFrom dplyr group_rows
#' @importFrom dplyr group_keys
#' @importFrom dplyr bind_rows
#' @importFrom dplyr pull
#' @importFrom tidyr unite
#' @importFrom tidyr separate
#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom purrr set_names
#'
#'
#' @export
setMethod("aggregate_cells", "SingleCellExperiment",  function(.data,
         .sample=NULL, assays=NULL,
         aggregation_function=Matrix::rowSums,
         ...) {
  
  # Fix NOTEs
  feature <- NULL
  .sample <- enquo(.sample)

  arg_list <- c(mget(ls(environment(), sorted = F)), match.call(expand.dots = F)$...)
  assays_to_use <- eval(arg_list$assays)
  if (is.null(assays_to_use)) assays_to_use <- tail(names(assays(.data)), n = 1)

  sample_groups <- .data |>
    as_tibble() |>
    group_by(pick({{ .sample }}))

  sample_group_idx <- sample_groups |>
    group_rows()

  sample_group_keys <- sample_groups |>
    group_keys()

  .sample_names <- colnames(sample_group_keys)

  grouping_factor_names <- sample_group_keys |>
    unite(col = "grouping_factor", !!.sample, sep = "___") |>
    pull(grouping_factor)

  sce_split <- map(.x = seq_along(sample_group_idx), .f = \(.num) .data[, sample_group_idx[[.num]]]) |>
    purrr::set_names(grouping_factor_names)

  grouping_factor <-
    .data |>
    colData() |>
    as_tibble() |>
    select(!!.sample) |>
    suppressMessages() |>
    unite("my_id_to_split_by___", !!.sample, sep = "___") |>
    pull(my_id_to_split_by___) |>
    as.factor()

  list_count_cells <- table(grouping_factor) |>
    enframe(name = "grouping_factor", value = ".aggregated_cells") |>
    mutate(.aggregated_cells = as.integer(.aggregated_cells))

  feature_df <- get_all_features(.data)
  selected_features <- feature_df[feature_df$assay_id %in% assays_to_use, ]
  selected_experiments_list <- split(x = selected_features, f = as.character(selected_features$exp_id))
  if ("Main" %in% names(selected_experiments_list)) selected_experiments_list <- selected_experiments_list[c("Main", setdiff(names(selected_experiments_list), "Main"))]

  aggregate_assays_fun <- function(exp) {
    selected_exp <- unique(exp$exp_id)
    selected_assays <- exp |> distinct(assay_name, .keep_all = TRUE)
    if (selected_exp == "Main") {
      aggregate_sce_fun <- function(sce) {
        aggregated_vals <- assays(sce)[selected_assays$assay_name] |>
          as.list() |>
          map(.f = \(.list) aggregation_function(.list))
        map(.x = seq_along(aggregated_vals), \(.num) enframe(x = aggregated_vals[[.num]], name = ".feature", value = selected_assays$assay_id[[.num]])) |>
          suppressMessages(reduce(full_join))
      }
      aggregated_list <- lapply(sce_split, aggregate_sce_fun) |>
        list_transpose() |>
        map(.f = \(.list) .list |> bind_rows(.id = "grouping_factor"))
      interim_res <- map(.x = seq_along(aggregated_list), .f = \(.num) aggregated_list[[.num]] |> 
            separate(col = grouping_factor, into = .sample_names, sep = "___")) |> 
        purrr::set_names(nm = selected_exp)
      map(.x = seq_along(interim_res), .f = \(.num) interim_res[[.num]] |> mutate(data_source = names(interim_res)[[.num]])) |> 
        purrr::reduce(full_join) |> 
        mutate(data_source = ifelse(data_source == "Main", yes = "RNA", no = data_source)) |> 
        select(data_source, everything())
    } else {
      aggregate_sce_fun <- function(sce) {
        aggregated_vals <- assays(altExps(sce)[[selected_exp]])[selected_assays$assay_name] |>
          as.list() |>
          set_names(selected_assays$assay_id) |>
          map(.f = \(.list) aggregation_function(.list))
        map(.x = seq_along(aggregated_vals), \(.num) enframe(x = aggregated_vals[[.num]], name = ".feature", value = selected_assays$assay_id[[.num]])) |>
          suppressMessages(reduce(full_join))
      }
      aggregated_list <- lapply(sce_split, aggregate_sce_fun) |>
        list_transpose() |>
        map(.f = \(.list) .list |> bind_rows(.id = "grouping_factor"))
      interim_res <- map(.x = seq_along(aggregated_list), .f = \(.num) aggregated_list[[.num]] |> 
                           separate(col = grouping_factor, into = .sample_names, sep = "___")) |> 
        purrr::set_names(nm = selected_exp)
      map(.x = seq_along(interim_res), .f = \(.num) interim_res[[.num]] |> 
            mutate(data_source = names(interim_res)[[.num]])) |>
        purrr::reduce(full_join) |> 
        mutate(data_source = ifelse(data_source == "Main", yes = "RNA", no = data_source)) |> 
        select(data_source, everything())
    }
  }
  se <- lapply(selected_experiments_list, aggregate_assays_fun) |> 
    purrr::reduce(full_join)
    se |> 
    tidybulk::as_SummarizedExperiment(
      .sample = .sample_names,
      .transcript = .feature,
      .abundance = setdiff(colnames(se), c("data_source", .sample_names, ".feature")))
})
