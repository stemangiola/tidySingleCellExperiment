#' (DEPRECATED) Extract and join information for transcripts.
#'
#'
#' @description join_transcripts() extracts and joins information for specified transcripts
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name join_transcripts
#' @rdname join_transcripts
#'
#' @param .data A tidySingleCellExperiment object
#' @param transcripts A vector of transcript identifiers to join
#' @param all If TRUE return all
#' @param exclude_zeros If TRUE exclude zero values
#' @param shape Format of the returned table "long" or "wide"
#' @param ... Parameters to pass to join wide, i.e. assay name to extract transcript abundance from
#'
#' @details DEPRECATED, please use join_features()
#'
#' @return A `tbl` containing the information.for the specified transcripts
#'
#' @examples
#'
#' print("DEPRECATED")
#'
#'
#' @export
#'
join_transcripts <- 
    function(.data,
    transcripts=NULL,
    all=FALSE,
    exclude_zeros=FALSE,
    shape="long", ...)
    {
        UseMethod("join_transcripts", .data)
    }
#' @export
join_transcripts.default <-
    function(.data,
        transcripts=NULL,
        all=FALSE,
        exclude_zeros=FALSE,
        shape="long", ...)
    {
        print("tidySingleCellExperiment says:",
            " This function cannot be applied to this object")
    }
#' @export
join_transcripts.Seurat <-
    function(.data,
        transcripts=NULL,
        all=FALSE,
        exclude_zeros=FALSE,
        shape="long", ...)
    {
        deprecate_warn(
            "1.1.2", "join_transcripts()", 
            "tidySingleCellExperiment::join_features()")
        
        .data %>%
            join_features(features=transcripts,
                all=all,
                exclude_zeros=exclude_zeros,
                shape=shape, ...)
    }
