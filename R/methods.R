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
#' @export
setMethod("join_features", "SingleCellExperiment", function(.data,
    features=NULL, all=FALSE, exclude_zeros=FALSE, shape="long", ...) {
    # CRAN Note
    .cell <- NULL
    .feature <- NULL
    
    # Shape is long
    if (shape == "long") {
        .data %>%
        left_join(
            by=c_(.data)$name,
            get_abundance_sc_long(
                .data=.data,
                features=features,
                all=all,
                exclude_zeros=exclude_zeros)) %>%
        select(!!c_(.data)$symbol, .feature, 
            contains(".abundance"), everything())
    # Shape if wide
    } else {
        .data  %>% 
            left_join(
                by=c_(.data)$name,
                get_abundance_sc_wide(
                    .data=.data,
                    features=features,
                    all=all, ...))
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
#' @export
setMethod("aggregate_cells", "SingleCellExperiment", function(.data,
    .sample=NULL, slot="data", assays=NULL, 
    aggregation_function=Matrix::rowSums,
    .metadata_columns_to_keep = NULL) {
    
    # Fix NOTEs
    feature <- NULL
    .sample <- enquo(.sample)
    
    # Subset only wanted assays
    if (!is.null(assays)) {
        assays(.data) <- assays(.data)[assays]
    }
    
    .data %>%
        nest(data=-!!.sample) %>%
        mutate(.aggregated_cells=as.integer(map(data, ~ ncol(.x)))) %>% 
        mutate(
            data=map(data, ~ {
                # Loop over assays
                map2(as.list(assays(.x)), assayNames(.x), ~ {
                    # Get counts
                    .x %>%
                        aggregation_function(na.rm=TRUE) %>%
                        enframe(
                            name ="feature",
                            value=sprintf("%s", .y)) %>%
                        mutate(feature=as.character(feature))
                }) %>% 
                Reduce(function(...) full_join(..., by="feature"), .)
            })
        ) %>%
        left_join(
            .data %>% 
                as_tibble() %>% 
                subset(!!.sample), 
            by=quo_names(.sample)) %>%
        unnest(data) %>%
        drop_class("tidySingleCellExperiment_nested") %>%
        as_SummarizedExperiment(
            .sample=!!.sample, 
            .transcript=feature, 
            .abundance=!!as.symbol(names(.data@assays)))
})
