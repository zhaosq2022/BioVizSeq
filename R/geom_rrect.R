#' Rounded rectangle
#' 
#' draws rounded rectangle by using the locations 
#' of the four corners (xmin, xmax, ymin and ymax) like geom_rect().
#' 
#' @title geom_rrect
#' 
#' @eval ggplot2:::rd_aesthetics("geom", "rrect")
#' 
#' @param mapping Set of aesthetic mappings created by [aes()]. If specified and
#'   `inherit.aes = TRUE` (the default), it is combined with the default mapping
#'   at the top level of the plot. You must supply `mapping` if there is no plot
#'   mapping.
#' @param data A data.frame, or other object, will override the plot data. All objects will be
#'             fortified to produce a data frame.
#' @param stat Name of stat to modify data.
#' @param position The position adjustment to use for overlapping points on this layer.
#' @param r The radius of rounded corners. 
#' @param na.rm If "FALSE" (default), missing values are removed with a warning. If "TRUE", missing values are silently removed, logical.
#' @param show.legend Whether to show legend, logical.
#' @param inherit.aes Whether to inherit aesthetic mappings, logical, defaults to "TRUE". 
#' @param ... additional parameter, e.g. color, linewidth, alpha.
#' @export
#' @author Shiqi Zhao
#' @return ggplot object
#' @examples
#' library(ggplot2)
#' df <- data.frame(
#'  xmin = c(1, 2, 3),
#'  xmax = c(2, 3, 4),
#'  ymin = c(1, 2, 3),
#'  ymax = c(2, 3, 4),
#'  category = c("A", "B", "C")
#'  )
#'  
#'  p <- ggplot(df) +
#'    geom_rrect(aes(xmin = xmin, xmax = xmax, 
#'               ymin = ymin, ymax = ymax, fill = category), 
#'              r = 0.4, linewidth = 1, colour = "black") 
#'  
#'  print(p)
#' 
geom_rrect <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       r = 0.2,
                       ...,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRrect,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      r = r,
      na.rm = na.rm,
      ...
    )
  )
}

#' @importFrom ggplot2 ggproto
#' @importFrom grid roundrectGrob
#' @importFrom ggplot2 draw_key_polygon
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 Geom

GeomRrect <- ggproto(
  "GeomRrect", Geom,
  default_aes = aes(
    colour = NA, fill = "grey35", linewidth = 0.5, 
    linetype = 1, alpha = NA
  ),
  required_aes = c("xmin", "xmax", "ymin", "ymax"),
  
  draw_panel = function(self, data, panel_params, coord, r = 0.2) {
    coords <- coord$transform(data, panel_params)
    
    gl <- lapply(seq_along(coords$xmin), function(i) {
      grid::roundrectGrob(
        x = coords$xmin[i], 
        y = coords$ymin[i],
        width = (coords$xmax[i] - coords$xmin[i]),
        height = (coords$ymax[i] - coords$ymin[i]),
        r = grid::unit(r, "native"),
        default.units = "native",
        just = c("left", "bottom"),
        gp = grid::gpar(
          col = coords$colour[i],
          fill = alpha(coords$fill[i], coords$alpha[i]),
          lwd = coords$linewidth[i] * .pt,
          lty = coords$linetype[i],
          lineend = "butt"
        )
      )
    })
    
    grobs <- do.call(grid::gList, gl)
    grob <- grid::grobTree(children = grobs)
    grob$name <- grid::grobName(grob, "geom_rrect")
    grob
  },
  draw_key = ggplot2::draw_key_polygon,
  rename_size = TRUE
)