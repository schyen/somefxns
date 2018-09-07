#' pca_biplot
#'
#' uses ggfortify and ggplot2's autoplot to make biplot.
#' uses prcomp to calculate pca
#'
#' @param data dataframe. Independant variable in rows,
#'    dependant variable in columns. Variable names in row and column names
#' @param metadata dataframe. metadata for the independant variables
#'    rownames must match rownames in data
#' @param colour string. column name in metadata to colour code score
#' @param label logic. default TRUE. label score points
#' @param loadings logic. default TRUE. When set, show loadings on biplot
#' @param loadings.label logic. default TRUE. When set, label loadings
#' @param loadings.label.colour string. default 'darkred'
#'
#' @return
#'    list of:
#'         \code{p_biplot} ggplot biplot
#'         \code{pcx} prcomp object.
#'         \code{pcx_summary} dataframe. summary of pca

pca_biplot <- function(data, metadata, colour = NULL, label = TRUE,
                       loadings = TRUE, loadings.label = TRUE,
                       loadings.label.colour = 'darkred') {

  # check inputs----------------------------------------------------------------
  if(is.null(colour)) colour <- 'black'
  else {
    check <- any(colour %in% colnames(metadata))
    print(check)
    if(!check) {
      stop("colour must be a column name in metadata", call. = FALSE)
    }
  }

  check <- identical(sort(rownames(data)), sort(rownames(metadata)))
  if(!check) {
    stop("rownames of data must be identical to rownames of metadata")
  }
  # perform pca-----------------------------------------------------------------
  pcx <- prcomp(data)

  pcx_summary <- summary(pcx)

  message("\nSummary of PCA:\n")
  print(pcx_summary)

  # plottig biplot--------------------------------------------------------------
  p_biplot <- ggplot2::autoplot(pcx, data = metadata,
                                # specifying score labelling
                                colour = colour, label = label,
                                label.repel = TRUE,
                                # specifying loading labelling
                                loadings = loadings,
                                loadings.label = loadings.label,
                                loadings.label.colour = loadings.label.colour)

  p_biplot <- p_biplot +
    ggplot2::scale_shape_manual(values=21:25) +
    ggplot2::theme_bw(10) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(face='bold'),
                   axis.title.y = ggplot2::element_text(face='bold'),
                   panel.grid = ggplot2::element_blank(),
                   legend.key = ggplot2::element_rect(colour='white'),
                   legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                   panel.border = ggplot2::element_rect(size=2, colour='black'),
                   panel.spacing = ggplot2::unit(.05, 'npc'),
                   strip.text.x = ggplot2::element_text(face='bold'))

  print(p_biplot)

  return(list('pcx' = pcx, 'pcx_summary' = pcx_summary, 'p_biplot' = p_biplot))
}
