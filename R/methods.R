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
#' @importFrom dplyr group_split
#' @importFrom dplyr full_join
#'
#'
#' @export
setMethod("aggregate_cells", "SingleCellExperiment",  function(.data,
                               .sample = NULL,
                               slot = "data",
                               assays = NULL,
                               aggregation_function = Matrix::rowSums) {
  # Fix NOTEs
  feature <- NULL

  .sample <- enquo(.sample)

  # Subset only wanted assays
  if (!is.null(assays)) {
    assay_info <- get_all_assays(.data)
    if (!any(assay_info$assay_id %in% assays)) stop("Please select an appropriate assay name")
    selected_assays <- assay_info[assay_info$assay_id %in% assays, ]
    selected_experiments_list <- split(x = selected_assays, f = as.character(selected_assays$exp_id))
    if ("Main" %in% names(selected_experiments_list)) selected_experiments_list <- selected_experiments_list[c("Main", setdiff(names(selected_experiments_list), "Main"))]

    aggregate_exp <- function(exp) {
      selected_exp <- unique(exp$exp_id)
      if (selected_exp == "Main") {
        .data@assays@data <- .data@assays@data[exp$assay_name]
      } else {
        col_data <- colData(.data)
        .data <- altExps(.data)[[selected_exp]]
        colData(.data) <- col_data
        .data@assays@data <- .data@assays@data[exp$assay_name]
        names(.data@assays@data) <- exp$assay_id
      }
      .data %>%
        nest(data = -!!.sample) %>%
        mutate(.aggregated_cells = as.integer(map(data, ~ ncol(.x)))) %>%
        mutate(data = map(data, ~

          # loop over assays
          map2(
            as.list(assays(.x)), names(.x@assays),

            # Get counts
            ~ .x %>%
              aggregation_function(na.rm = TRUE) %>%
              enframe(
                name  = "feature",
                value = sprintf("%s", .y)
              ) %>%
              mutate(feature = as.character(feature))
          ) %>%
            Reduce(function(...) full_join(..., by = c("feature")), .))) %>%
        left_join(.data %>% as_tibble() %>% subset(!!.sample), by = quo_names(.sample)) %>%
        unnest(data)
    }
    suppressMessages({
      Reduce(f = full_join, x = lapply(selected_experiments_list, aggregate_exp))
    })
  } %>%
    drop_class("tidySingleCellExperiment_nested") %>%
    as_SummarizedExperiment(.sample = !!.sample, .transcript = feature, .abundance = !!as.symbol(names(.data@assays)))
})
