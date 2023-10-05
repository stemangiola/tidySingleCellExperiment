#' @name plot_ly
#' @rdname plot_ly
#' @aliases plot_ly,SingleCellExperiment-method
#' @inherit plotly::plot_ly
#' @return `plotly`
#' 
#' @examples 
#' data(pbmc_small)
#' pbmc_small |> 
#'     plot_ly(x = ~ nCount_RNA, y = ~ nFeature_RNA)
#' 
#' @importFrom tidySummarizedExperiment plot_ly
#' @export
setMethod("plot_ly", "SingleCellExperiment",
    function(data=data.frame(), ..., 
             type=NULL, name=NULL, 
             color=NULL, colors=NULL, alpha=NULL,
             stroke=NULL, strokes=NULL, alpha_stroke=1,
             size=NULL, sizes=c(10, 100),
             span=NULL, spans=c(1, 20),
             symbol=NULL, symbols=NULL,
             linetype=NULL, linetypes=NULL,
             split=NULL, frame=NULL,
             width=NULL, height=NULL, source="A") {
      
      data %>%
        
          # This is a trick to not loop the call
          as_tibble() %>%
          plot_ly(..., type=type, name=name,
                  color=color, colors=colors, alpha=alpha,
                  stroke=stroke, strokes=strokes, alpha_stroke=alpha_stroke,
                  size=size, sizes=sizes, 
                  span=span, spans=spans,
                  symbol=symbol, symbols=symbols, 
                  linetype=linetype, linetypes=linetypes,
                  split=split, frame=frame,
                  width=width, height=height, source=source)
})
