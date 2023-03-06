#' This is a generalisation of ifelse that accepts an object and return an objects
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @importFrom purrr as_mapper
#' @importFrom magrittr equals
#'
#' @param .x A tibble
#' @param .p A boolean
#' @param .f1 A function
#' @param .f2 A function
#'
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
  switch(.p %>% not() %>% sum(1),
         as_mapper(.f1)(.x),
         if (.f2 %>% is.null %>% not())
           as_mapper(.f2)(.x)
         else
           .x)
}

#' as_SummarizedExperiment
#'
#' @keywords internal
#' @noRd
#' 
#' @description as_SummarizedExperiment creates a `SummarizedExperiment` object from a `tbl` or `tidybulk` tbl formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#' @importFrom utils data
#' @importFrom tidyr pivot_longer
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @return A `SummarizedExperiment` object
as_SummarizedExperiment = function(.data,
                                    .sample = NULL,
                                    .transcript = NULL,
                                    .abundance = NULL) {
  
  # Get column names
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)
  col_names = get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
  .sample = col_names$.sample
  .transcript = col_names$.transcript
  .abundance = col_names$.abundance
  
  # If present get the scaled abundance
  .abundance_scaled =
    .data %>%
    ifelse_pipe(
      ".abundance_scaled" %in% ((.) %>% get_tt_columns() %>% names) &&
        # .data %>% get_tt_columns() %$% .abundance_scaled %>% is.null %>% not() &&
        quo_name((.) %>% get_tt_columns() %$% .abundance_scaled) %in% ((.) %>% colnames),
      ~ .x %>% get_tt_columns() %$% .abundance_scaled,
      ~ NULL
    )
  
  # Get which columns are sample wise and which are feature wise
  col_direction = get_x_y_annotation_columns(.data,
                                             !!.sample,
                                             !!.transcript,
                                             !!.abundance,
                                             !!.abundance_scaled)
  sample_cols = col_direction$horizontal_cols
  feature_cols = col_direction$vertical_cols
  counts_cols = col_direction$counts_cols
  
  colData =
    .data %>%
    select(!!.sample, sample_cols) %>%
    distinct() %>%
    
    # Unite if multiple sample columns
    tidyr::unite(!!sample__$name, !!.sample, remove = FALSE, sep = "___") |>
    
    arrange(!!sample__$symbol) %>% {
      S4Vectors::DataFrame(
        (.) %>% select(-!!sample__$symbol),
        row.names = (.) %>% pull(!!sample__$symbol)
      )
    }
  
  rowData =
    .data %>%
    select(!!.transcript, feature_cols) %>%
    distinct() %>%
    
    # Unite if multiple sample columns
    tidyr::unite(!!feature__$name, !!.transcript, remove = FALSE, sep = "___") |>
    
    arrange(!!feature__$symbol) %>% {
      S4Vectors::DataFrame(
        (.) %>% select(-!!feature__$symbol),
        row.names = (.) %>% pull(!!feature__$symbol)
      )
    }
  
  my_assays =
    .data %>%
    
    # Unite if multiple sample columns
    tidyr::unite(!!sample__$name, !!.sample, remove = FALSE, sep = "___") |>
    
    # Unite if multiple sample columns
    tidyr::unite(!!feature__$name, !!.transcript, remove = FALSE, sep = "___") |>
    
    select(!!sample__$symbol,
           !!feature__$symbol,
           !!.abundance,
           !!.abundance_scaled,
           counts_cols) %>%
    distinct() %>%
    
    pivot_longer( cols=-c(!!feature__$symbol,!!sample__$symbol), names_to="assay", values_to= ".a") %>%
    tidyr::nest(`data` = -`assay`) %>%
    mutate(`data` = `data` %>%  map(
      ~ .x %>%
        spread(!!sample__$symbol, .a) %>%
        
        # arrange sample
        select(!!feature__$symbol, rownames(colData)) |>
        
        # Arrange symbol
        arrange(!!feature__$symbol) |>
        
        # Convert
        as_matrix(rownames = feature__$name)
    ))
  
  # Build the object
  SummarizedExperiment::SummarizedExperiment(
    assays = my_assays %>% pull(`data`) %>% stats::setNames(my_assays$assay),
    rowData = rowData,
    colData = colData
  )
  
}

#' Get column names either from user or from attributes
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang quo_is_symbol
#' @importFrom rlang quo_is_symbolic
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#'
#' @return A list of column enquo or error
get_sample_transcript_counts = function(.data, .sample, .transcript, .abundance){
  
  if( quo_is_symbolic(.sample) ) .sample = .sample
  else if(".sample" %in% (.data %>% get_tt_columns() %>% names))
    .sample =  get_tt_columns(.data)$.sample
  else my_stop()
  
  if( quo_is_symbolic(.transcript) ) .transcript = .transcript
  else if(".transcript" %in% (.data %>% get_tt_columns() %>% names))
    .transcript =  get_tt_columns(.data)$.transcript
  else my_stop()
  
  if(  quo_is_symbolic(.abundance) ) .abundance = .abundance
  else if(".abundance" %in% (.data %>% get_tt_columns() %>% names))
    .abundance = get_tt_columns(.data)$.abundance
  else my_stop()
  
  list(.sample = .sample, .transcript = .transcript, .abundance = .abundance)
  
}

get_tt_columns = function(.data){
  if(
    .data %>% attr("internals") %>% is.list() &&
    "tt_columns" %in% names(.data %>% attr("internals"))
  ) #& "internals" %in% (.data %>% attr("internals") %>% names()))
    .data %>% attr("internals") %$% tt_columns
  else NULL
}

#' get_x_y_annotation_columns
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_symbol
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .horizontal The name of the column horizontally presented in the heatmap
#' @param .vertical The name of the column vertically presented in the heatmap
#' @param .abundance The name of the transcript/gene abundance column
#' @param .abundance_scaled The name of the transcript/gene scaled abundance column
#'
#' @description This function recognise what are the sample-wise columns and transcrip-wise columns
#'
#' @return A list
get_x_y_annotation_columns = function(.data, .horizontal, .vertical, .abundance, .abundance_scaled){
  
  # Comply with CRAN NOTES
  . = NULL
  
  # Make col names
  .horizontal = enquo(.horizontal)
  .vertical = enquo(.vertical)
  .abundance = enquo(.abundance)
  .abundance_scaled = enquo(.abundance_scaled)
  
  # x-annotation df
  n_x = .data %>% select(!!.horizontal) |> distinct() |> nrow()
  n_y = .data %>% select(!!.vertical) |> distinct() |> nrow()
  
  # Sample wise columns
  horizontal_cols=
    .data %>%
    select(-!!.horizontal, -!!.vertical, -!!.abundance) %>%
    colnames %>%
    map(
      ~
        .x %>%
        when(
          .data %>%
            select(!!.horizontal, !!as.symbol(.x)) %>%
            distinct() |>
            nrow() %>%
            equals(n_x) ~ .x,
          ~ NULL
        )
    ) %>%
    
    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist
  
  # Transcript wise columns
  vertical_cols=
    .data %>%
    select(-!!.horizontal, -!!.vertical, -!!.abundance, -horizontal_cols) %>%
    colnames %>%
    map(
      ~
        .x %>%
        ifelse_pipe(
          .data %>%
            select(!!.vertical, !!as.symbol(.x)) |>
            distinct() |>
            nrow() %>%
            equals(n_y),
          ~ .x,
          ~ NULL
        )
    ) %>%
    
    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist
  
  # Counts wise columns, at the moment scaled counts is treated as special and not accounted for here
  counts_cols =
    .data %>%
    select(-!!.horizontal, -!!.vertical, -!!.abundance) %>%
    
    # Exclude horizontal
    ifelse_pipe(!is.null(horizontal_cols),  ~ .x %>% select(-horizontal_cols)) %>%
    
    # Exclude vertical
    ifelse_pipe(!is.null(vertical_cols),  ~ .x %>% select(-vertical_cols)) %>%
    
    # Exclude scaled counts if exist
    ifelse_pipe(.abundance_scaled %>% quo_is_symbol,  ~ .x %>% select(-!!.abundance_scaled) ) %>%
    
    # Select colnames
    colnames %>%
    
    # select columns
    map(
      ~
        .x %>%
        ifelse_pipe(
          .data %>%
            select(!!.vertical, !!.horizontal, !!as.symbol(.x)) %>%
            distinct() |>
            nrow() %>%
            equals(n_x * n_y),
          ~ .x,
          ~ NULL
        )
    ) %>%
    
    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist
  
  list(  horizontal_cols = horizontal_cols,  vertical_cols = vertical_cols, counts_cols = counts_cols )
}

#' Get matrix from tibble
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames The column name of the input tibble that will become the rownames of the output matrix
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#' tibble(.feature = "CD3G", count=1) |> as_matrix(rownames=.feature)
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  rownames = enquo(rownames)
  tbl %>%
    
    # Through warning if data frame is not numerical beside the rownames column (if present)
    ifelse_pipe(
      do_check &&
        tbl %>%
        # If rownames defined eliminate it from the data frame
        ifelse_pipe(!quo_is_null(rownames), ~ .x[,-1], ~ .x) %>%
        dplyr::summarise_all(class) %>%
        tidyr::gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% not() %>% any(),
      ~ {
        warning("tidybulk says: there are NON-numerical columns, the matrix will NOT be numerical")
        .x
      }
    ) %>%
    as.data.frame() %>%
    
    # Deal with rownames column if present
    ifelse_pipe(
      !quo_is_null(rownames),
      ~ .x %>%
        magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
        select(-1)
    ) %>%
    
    # Convert to matrix
    as.matrix()
}