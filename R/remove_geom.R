#' Remove geom layer

#' @param ggplot2_object ggplot plot object
#'     and the variables involved (i.e. compound or species)
#' @param geom_type string; layer to remove

#' @return ggplot object
#' @export


remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  return(ggplot2_object)
}
