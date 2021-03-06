#' Initiate a plotly visualization
#'
#'
#' This function maps R objects to [plotly.js](https://plot.ly/javascript/),
#' an (MIT licensed) web-based interactive charting library. It provides
#' abstractions for doing common things (e.g. mapping data values to
#' fill colors (via `color`) or creating [animation]s (via `frame`)) and sets
#' some different defaults to make the interface feel more 'R-like'
#' (i.e., closer to [plot()] and [ggplot2::qplot()]).
#'
#' @details Unless `type` is specified, this function just initiates a plotly
#' object with 'global' attributes that are passed onto downstream uses of
#' [add_trace()] (or similar). A [formula] must always be used when
#' referencing column name(s) in `data` (e.g. `plot_ly(mtcars, x=~wt)`).
#' Formulas are optional when supplying values directly, but they do
#' help inform default axis/scale titles
#' (e.g., `plot_ly(x=mtcars$wt)` vs `plot_ly(x=~mtcars$wt)`)
#'
#' @param data A data frame (optional) or [crosstalk::SharedData] object.
#' @param ... Arguments (i.e., attributes) passed along to the trace `type`.
#' See [schema()] for a list of acceptable attributes for a given trace `type`
#' (by going to `traces` -> `type` -> `attributes`). Note that attributes
#' provided at this level may override other arguments
#' (e.g. `plot_ly(x=1:10, y=1:10, color=I("red"), marker=list(color="blue"))`).
#' @param type A character string specifying the trace type
#'   (e.g. `"scatter"`, `"bar"`, `"box"`, etc).
#' If specified, it *always* creates a trace, otherwise
#' @param name Values mapped to the trace's name attribute. Since a trace can
#' only have one name, this argument acts very much like `split` in that it
#' creates one trace for every unique value.
#' @param color Values mapped to relevant 'fill-color' attribute(s)
#' (e.g. [fillcolor](https://plot.ly/r/reference#scatter-fillcolor),
#' [marker.color](https://plot.ly/r/reference#scatter-marker-color),
#' [textfont.color](https://plot.ly/r/reference/#scatter-textfont-color), etc.).
#' The mapping from data values to color codes may be controlled using
#' `colors` and `alpha`, or avoided altogether via [I()]
#'   (e.g., `color=I("red")`).
#' Any color understood by [grDevices::col2rgb()] may be used in this way.
#' @param colors Either a colorbrewer2.org palette name
#'   (e.g. "YlOrRd" or "Blues"),
#' or a vector of colors to interpolate in hexadecimal "#RRGGBB" format,
#' or a color interpolation function like `colorRamp()`.
#' @param stroke Similar to `color`, but values are mapped to relevant 'stroke-color' attribute(s)
#' (e.g., [marker.line.color](https://plot.ly/r/reference#scatter-marker-line-color)
#'  and [line.color](https://plot.ly/r/reference#scatter-line-color)
#' for filled polygons). If not specified, `stroke` inherits from `color`.
#' @param strokes Similar to `colors`, but controls the `stroke` mapping.
#' @param alpha A number between 0 and 1 specifying the alpha channel applied to `color`.
#' Defaults to 0.5 when mapping to [fillcolor](https://plot.ly/r/reference#scatter-fillcolor) and 1 otherwise.
#' @param alpha_stroke Similar to `alpha`, but applied to `stroke`.
#' @param symbol (Discrete) values mapped to [marker.symbol](https://plot.ly/r/reference#scatter-marker-symbol).
#' The mapping from data values to symbols may be controlled using
#' `symbols`, or avoided altogether via [I()] (e.g., `symbol=I("pentagon")`).
#' Any [pch] value or [symbol name](https://plot.ly/r/reference#scatter-marker-symbol) may be used in this way.
#' @param symbols A character vector of [pch] values or [symbol names](https://plot.ly/r/reference#scatter-marker-symbol).
#' @param linetype (Discrete) values mapped to [line.dash](https://plot.ly/r/reference#scatter-line-dash).
#' The mapping from data values to symbols may be controlled using
#' `linetypes`, or avoided altogether via [I()] (e.g., `linetype=I("dash")`).
#' Any `lty` (see [par]) value or [dash name](https://plot.ly/r/reference#scatter-line-dash) may be used in this way.
#' @param linetypes A character vector of `lty` values or [dash names](https://plot.ly/r/reference#scatter-line-dash)
#' @param size (Numeric) values mapped to relevant 'fill-size' attribute(s)
#' (e.g., [marker.size](https://plot.ly/r/reference#scatter-marker-size),
#' [textfont.size](https://plot.ly/r/reference#scatter-textfont-size),
#' and [error_x.width](https://plot.ly/r/reference#scatter-error_x-width)).
#' The mapping from data values to symbols may be controlled using
#' `sizes`, or avoided altogether via [I()] (e.g., `size=I(30)`).
#' @param sizes A numeric vector of length 2 used to scale `size` to pixels.
#' @param span (Numeric) values mapped to relevant 'stroke-size' attribute(s)
#' (e.g.,
#' [marker.line.width](https://plot.ly/r/reference#scatter-marker-line-width),
#' [line.width](https://plot.ly/r/reference#scatter-line-width) for filled polygons,
#' and [error_x.thickness](https://plot.ly/r/reference#scatter-error_x-thickness))
#' The mapping from data values to symbols may be controlled using
#' `spans`, or avoided altogether via [I()] (e.g., `span=I(30)`).
#' @param spans A numeric vector of length 2 used to scale `span` to pixels.
#' @param split (Discrete) values used to create multiple traces (one trace per value).
#' @param frame (Discrete) values used to create animation frames.
#' @param width Width in pixels (optional, defaults to automatic sizing).
#' @param height Height in pixels (optional, defaults to automatic sizing).
#' @param source a character string of length 1. Match the value of this string
#' with the source argument in [event_data()] to retrieve the
#' event data corresponding to a specific plot (shiny apps can have multiple plots).
#' @author Carson Sievert
#' @references <https://plotly-r.com/overview.html>
#' @seealso \itemize{
#'  \item For initializing a plotly-geo object: [plot_geo()]
#'  \item For initializing a plotly-mapbox object: [plot_mapbox()]
#'  \item For translating a ggplot2 object to a plotly object: [ggplotly()]
#'  \item For modifying any plotly object: [layout()], [add_trace()], [style()]
#'  \item For linked brushing: [highlight()]
#'  \item For arranging multiple plots: [subplot()], [crosstalk::bscols()]
#'  \item For inspecting plotly objects: [plotly_json()]
#'  \item For quick, accurate, and searchable plotly.js reference: [schema()]
#' }
#'
#' @return A plotly
#' @export
#' @examples
#' \dontrun{
#' # plot_ly() tries to create a sensible plot based on the information you
#' # give it. If you don't provide a trace type, plot_ly() will infer one.
#' plot_ly(economics, x=~pop)
#' plot_ly(economics, x=~date, y=~pop)
#' # plot_ly() doesn't require data frame(s), which allows one to take
#' # advantage of trace type(s) designed specifically for numeric matrices
#' plot_ly(z=~volcano)
#' plot_ly(z=~volcano, type="surface")
#'
#' # plotly has a functional interface: every plotly function takes a plotly
#' # object as it's first input argument and returns a modified plotly object
#' add_lines(plot_ly(economics, x=~date, y=~ unemploy / pop))
#'
#' # To make code more readable, plotly imports the pipe operator from magrittr
#' economics %>%
#'     plot_ly(x=~date, y=~ unemploy / pop) %>%
#'     add_lines()
#'
#' # Attributes defined via plot_ly() set 'global' attributes that
#' # are carried onto subsequent traces, but those may be over-written
#' plot_ly(economics, x=~date, color=I("black")) %>%
#'     add_lines(y=~uempmed) %>%
#'     add_lines(y=~psavert, color=I("red"))
#'
#' # Attributes are documented in the figure reference -> https://plot.ly/r/reference
#' # You might notice plot_ly() has named arguments that aren't in this figure
#' # reference. These arguments make it easier to map abstract data values to
#' # visual attributes.
#' p <- plot_ly(iris, x=~Sepal.Width, y=~Sepal.Length)
#' add_markers(p, color=~Petal.Length, size=~Petal.Length)
#' add_markers(p, color=~Species)
#' add_markers(p, color=~Species, colors="Set1")
#' add_markers(p, symbol=~Species)
#' add_paths(p, linetype=~Species)
#' }
#'
plot_ly <- function(data=data.frame(), ..., type=NULL, name=NULL,
    color=NULL, colors=NULL, alpha=NULL,
    stroke=NULL, strokes=NULL, alpha_stroke=1,
    size=NULL, sizes=c(10, 100),
    span=NULL, spans=c(1, 20),
    symbol=NULL, symbols=NULL,
    linetype=NULL, linetypes=NULL,
    split=NULL, frame=NULL,
    width=NULL, height=NULL, source="A") {
    UseMethod("plot_ly")
}

#' @export
#'
plot_ly.default <- function(data=data.frame(), ..., type=NULL, name=NULL,
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
        drop_class("tbl_df") %>%
        plotly::plot_ly(...,
            type=type, name=name,
            color=color, colors=colors, alpha=alpha,
            stroke=stroke, strokes=strokes, alpha_stroke=alpha_stroke,
            size=size, sizes=sizes,
            span=span, spans=spans,
            symbol=symbol, symbols=symbols,
            linetype=linetype, linetypes=linetypes,
            split=split, frame=frame,
            width=width, height=height, source=source
        )
}

#' @importFrom plotly plot_ly
#' @export
plot_ly.SingleCellExperiment <- function(data=data.frame(), ..., type=NULL, name=NULL,
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
        plot_ly(...,
            type=type, name=name,
            color=color, colors=colors, alpha=alpha,
            stroke=stroke, strokes=strokes, alpha_stroke=alpha_stroke,
            size=size, sizes=sizes,
            span=span, spans=spans,
            symbol=symbol, symbols=symbols,
            linetype=linetype, linetypes=linetypes,
            split=split, frame=frame,
            width=width, height=height, source=source
        )
}
