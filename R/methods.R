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
#' @importFrom rlang ensym
#' @importFrom rlang as_name
#' @importFrom magrittr "%>%"
#' @importFrom tibble enframe
#' @importFrom Matrix rowSums
#' @importFrom ttservice aggregate_cells
#' @importFrom SummarizedExperiment assays assays<- assayNames
#' @importFrom S4Vectors split
#' @importFrom stringr str_remove
#' @importFrom dplyr group_split
#' @importFrom dplyr full_join
#' @importFrom purrr reduce
#'
#'
#' @export
setMethod("aggregate_cells", "SingleCellExperiment",  function(.data,
         .sample=NULL, slot="data", assays=NULL,
         aggregation_function=Matrix::rowSums,
         ...) {
  
  # Fix NOTEs
  feature <- NULL
  .sample <- enquo(.sample)
  
  # Subset only wanted assays
  if (!is.null(assays)) {
    assay_info <- get_all_assays(.data)
    if(!any(assay_info$assay_id %in% assays)) stop("Please select an appropriate assay name")
    selected_assays <- assay_info[assay_info$assay_id %in% assays,]
    selected_exp <- unique(selected_assays$exp_id)
    selected_experiments_list <- split(x = selected_assays, f = as.character(selected_assays$exp_id))
    if("Main" %in% selected_exp) selected_experiments_list <- selected_experiments_list[c("Main", setdiff(names(selected_experiments_list), "Main"))]
  }
  
  aggregate_exps <- function(exp) {
    if(unique(exp$exp_id) == "Main") 
    {
      assays(.data) <- assays(.data)[exp$assay_name]
    } else {
      col_data <- colData(.data)
      .data <- altExps(.data)[[unique(exp$exp_id)]]
      colData(.data) <- col_data
      assays(.data) <- assays(.data)[exp$assay_name]
    }
    grouping_factor =
      .data |>
      colData() |>
      as_tibble() |>
      select(!!.sample) |>
      suppressMessages() |>
      unite("my_id_to_split_by___", !!.sample, sep = "___") |>
      pull(my_id_to_split_by___) |>
      as.factor()
    
    list_count_cells = table(grouping_factor) |> as.list()
    
    # New method
    list_assays =
      .data |>
      assays() |>
      as.list() |>
      map(~ .x |> splitColData(grouping_factor)) |>
      unlist(recursive=FALSE)
    
    list_assays =
      list_assays |>
      map2(names(list_assays), ~ {
        # Get counts
        .x %>%
          aggregation_function(na.rm=TRUE) %>%
          enframe(
            name =".feature",
            value="x") %>% # sprintf("%s", .y)) %>%
          
          # In case we don't have rownames
          mutate(.feature=as.character(.feature))
      }) |>
      enframe(name = ".sample") |>
      
      # Clean groups
      mutate(assay_name = assayNames(!!.data) |> rep(each=length(levels(grouping_factor)))) |>
      mutate(.sample = .sample |> str_remove(assay_name) |> str_remove("\\.")) |>
      group_split(.sample) |>
      map(~ .x |>  unnest(value) |> pivot_wider(names_from = assay_name, values_from = x) ) |>
      
      # Add cell count
      map2(
        list_count_cells,
        ~ .x |> mutate(.aggregated_cells = .y)
      )
    
    
    do.call(rbind, list_assays) |>
      
      left_join(
        .data |>
          colData() |>
          as_tibble() |>
          subset(!!.sample) |>
          unite("my_id_to_split_by___", !!.sample, remove=FALSE, sep = "___"),
        by= join_by(".sample" == "my_id_to_split_by___")
      )
  }
  
  lapply(selected_experiments_list, aggregate_exps) |> 
    bind_rows() |>
    as_SummarizedExperiment(
      .sample=.sample,
      .transcript=.feature,
      .abundance=!!as.symbol(names(.data@assays)))
})
