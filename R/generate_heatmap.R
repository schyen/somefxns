#' explore_heatmap.R
#'
#' Produce heatmap of count data. uses euclidean distances (recommdended) of clr transformed data; zero replacement completed by Bayesian-multiplicative replacement
#'
#' @param data dataframe. long format, with column of xdata, variable
#'     (such as taxonomy, or compounds), and value (such as abundance or concentration)
#' @param x_var string. variable to be heatmap x-axis. example: sample
#'     should be a column name in data
#' @param y_var string. variable to be heatmap y-axis. example: species or compound
#'     should be column name in data
#' @param value_var string. column containing heatmap values
#' @param x_dist string; dist method for x dendrogram;
#'     method for computing phylogenetic distance;
#'     supplied to the dist(); default method: "euclidean"
#'     dist method options: "euclidean", "maximum", "manhattan",
#'     "canberra", "binary" or "minkowski"
#' @param x_hclust  string; hclust method for x dendrogram;
#'     Linkage criteria to determine distance between sets of
#'     observations; to be used in hclust(); default method: "complete"
#'     hclust method options:"ward.D", "ward.D2", "single", "complete",
#'    "average" (= UPGMA), "mcquitty" (= WPGMA),
#'    "median" (= WPGMC) or "centroid" (= UPGMC)
#' @param y_dist string; dist method for ydata dendrogram;
#'     method for computing phylogenetic distance;
#'     supplied to the dist(); default method: "euclidean"
#'     dist method options: "euclidean", "maximum", "manhattan",
#'     "canberra", "binary" or "minkowski"
#' @param y_hclust string; hclust method for ydata dendrogram;
#'     Linkage criteria to determine distance between sets of
#'     observations; to be used in hclust(); default method: "complete"
#'     hclust method options:"ward.D", "ward.D2", "single", "complete",
#'     "average" (= UPGMA), "mcquitty" (= WPGMA),
#'     "median" (= WPGMC) or "centroid" (= UPGMC)
#'     When set to FALSE, xdata listed along y-axis
#' @param x_meta dataframe. default NULL. Any metadata on heatmap x axis
#'     variable (such as sample) to be included
#'     x_var should be the identifier to link metadata
#' @param y_meta dataframe. default NULL. any metadata  on heatmap y axis
#'     variable (such as species) to be included
#'     y_var should be identifier to link metadata
#' @param x_label string. default NULL. Should be column in x_meta.
#'     How x-axis of heatmap should be labelled. When set to NULL, x_var is used
#' @param y_label string. default NULL. Should be a column in y_meta.
#'     How y-axis of heatmap should be labelled. When set to NULL, y_var is used
#' @param heatmap_w numeric; default 20; width (cm) of heatmap panel
#' @param heatmap_h numeric; default 20; height (cm) of heatmap panel
#' @param dendro_x_h numeric; default 8; height (cm) of x-axis dendrogram;
#'     helps with ensuring x-axis dendrogram labels are not cut off
#' @param dendro_y_w numeric; default 8; width(cm) of y-axis dendrogram;
#'     helps with ensuring y-axis dendrogram labels are not cut off
#' @param output = string; path and file name of resulting heatmap png
#'     Example: "/path/to/my_heatmap.png"
#' @return
#'     \code{pdata} dataframe for heatmap
#'     \code{ddata_x} dataframe for x-axis dendrogram
#'     \code{ddata_x} dataframe for y-axis dendrogram
#'     \code{p_hl} ggplot of heatmap and legend
#'     \code{p_x} ggplot of x-axis dendrogram
#'     \code{p_y} ggplot of y-axis dendrogram
#'     \code{g} gtable grob of entire heatmap and dendrograms

generate_heatmap <- function(data, x_var, y_var, value_var,
                         x_dist = 'euclidean', x_hclust = 'complete',
                         y_dist = 'euclidean', y_hclust = 'complete',
                         x_meta = NULL, y_meta = NULL,
                         x_label = NULL, y_label = NULL,
                         heatmap_w = 20, heatmap_h = 20,
                         dendro_x_h = 8, dendro_y_w = 8,
                         out) {

  # formatting data into wide format so can calculate dendrograms
  ## heatmap x-axis data
  form <- as.formula(paste(x_var, y_var, sep = "~"))
  xdata <- reshape2::dcast(data, formula = form, value.var = value_var)

  # setting rownames
  row.names(xdata) <- xdata[ , x_var]
  xdata <- xdata[ ,colnames(xdata) != x_var]

  ## heatmap y-axis data
  form <- as.formula(paste(y_var, x_var, sep = "~"))
  ydata <- reshape2::dcast(data, formula = form, value.var = value_var)
  row.names(ydata) <- ydata[, y_var]
  ydata <- ydata[, colnames(ydata) != y_var]

  # calculating dendrograms------------------------------------------------------
  dd_xdata <- as.dendrogram(hclust(dist(xdata, method=x_dist), # x as row
                                     method=x_hclust))
  ord_xdata <- order.dendrogram(dd_xdata)

  dd_ydata <- as.dendrogram(hclust(dist(ydata, method=y_dist), # ydata as row
                                 method=y_hclust))
  ord_ydata <- order.dendrogram(dd_ydata)

  # extracting dendrogram data
  ddata_x <- ggdendro::dendro_data(dd_xdata)
  ddata_y <- ggdendro::dendro_data(dd_ydata)

  # reading in metadata if necessary--------------------------------------------
  if(!is.null(x_meta)) {
    if(is.null(x_label)) {
      stop("x_label must not be null when x_meta is specified")
    }

    # copying dendrogram labels
    new_x_lab <- ggdendro::label(ddata_x)
    new_x_lab$label <- as.character(new_x_lab$label)
    colnames(new_x_lab)[3] <- x_var

    # adding x_label
    new_x_lab <- dplyr::inner_join(new_x_lab, x_meta, by=x_var)

    # resetting strain dendrogram labels to taxonomy (or strain number)
    ddata_x$label <- new_x_lab[, x_label]
  }

  if(!is.null(y_meta)) {
    if(is.null(y_label)) {
      stop("y_label must not be null when y_meta is specified")
    }

    # copying dendrogram labels
    new_y_lab <- ggdendro::label(ddata_y)
    new_y_lab$label <- as.character(new_y_lab$label)
    colnames(new_y_lab)[3] <- y_var

    # adding y_label
    new_y_lab <- dplyr::inner_join(new_y_lab, y_meta, by=y_var)

    # resetting strain dendrogram labels to taxonomy (or strain number)
    ddata_y$label <- new_y_lab[, y_label]
  }


  # heatmap data-----------------------------------------------------------------
  ddata <- xdata[ord_xdata, ord_ydata]

  # column for sample names, to set order
  ddata$sampleID <- rownames(ddata)
  ddata$sampleID <- with(ddata, factor(sampleID, levels=sampleID,
                                       ordered=TRUE))
  # melting for ggplot
  pdata <- reshape2::melt(ddata, id.vars = x_var, variable.name = y_var,
                          value.name = value_var)

  ## Set up a blank theme--------------------------------------------------------
  theme_none <- function() {
    ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(colour=NA),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
    )
    }

  # creating heatmap--------------------------------------------------------------
  p_h <- ggplot2::ggplot(pdata, ggplot2::aes_string(x=x_var, y=y_var))

  p_hl <- p_h +
    ggplot2::geom_tile(ggplot2::aes_string(fill=value_var)) +
    ggplot2::scale_fill_gradient2(low = "#9A0038", mid = "#FFFFFF", high = "#003AC3",
                         midpoint=0,
                         breaks=c(ceiling(min(pdata[,value_var])), 0,
                                  floor(max(pdata[, value_var]))),
                         ggplot2::guide_legend(title='Normalized Count')) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(0,0,0,0, 'cm'),
          legend.direction='horizontal')

  # remove legend from p_hl
  p_h <- p_hl + ggplot2::theme(legend.position='none')

  ## X axis dendrogram 1----------------------------------------------------------
  # calculating plot limits
  xlim <- range(c(ggdendro::segment(ddata_x)$x, ggdendro::segment(ddata_x)$xend))
  ylim <- range(c(ggdendro::segment(ddata_x)$y, ggdendro::segment(ddata_x)$yend))

  x_name <- as.character(ggdendro::label(ddata_x)$label)
  name_len <- lapply(x_name, function(x) nchar(x))
  name_len <- unlist(name_len)
  name_len <- max(name_len)
  xadjust <- name_len*0.75

  xmin <- floor(xlim[1]) - 0.5
  xmax <- ceiling(xlim[2])+ 0.5
  xlim <- c(xmin, xmax)

  ymin <- floor(ylim[1]) - xadjust
  ymax <- ceiling(ylim[2]) + 1
  ylim <- c(ymin, ymax)

  # plotting x-axis dendrogram
  p_x <- ggplot2::ggplot() +
    ggplot2::geom_segment(data=ddata_x$segment,
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::geom_text(data=ddata_x$labels, ggplot2::aes(x=x, y=y, label=label),
                       angle=-90, hjust=0, vjust=0.3) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits=xlim) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits=ylim) +
    theme_none() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0,0.2,0,0.5), 'cm')) +#t,r,b,l
    ggplot2::labs(x=NULL, y=NULL)

  ## Y axis dendrogram 2----------------------------------------------------------
  # calculating plot limits
  xlim <- range(c(ggdendro::segment(ddata_y)$x, ggdendro::segment(ddata_y)$xend))
  ylim <- range(c(ggdendro::segment(ddata_y)$y, ggdendro::segment(ddata_y)$yend))
  ydata_name <- as.character(ggdendro::label(ddata_y)$label)
  name_len <- lapply(ydata_name, function(x) nchar(x))
  name_len <- unlist(name_len)
  name_len <- max(name_len)

  yadjust <- name_len*1.6

  xmin <- floor(xlim[1]) - 0.5
  xmax <- ceiling(xlim[2]) + 0.5
  xlim <- c(xmin, xmax)

  ymin <- floor(ylim[1]) - yadjust
  ymax <- ceiling(ylim[2])
  ylim <- c(ymin, ymax)

  # plotting y-axis dendrogram
  p_y <- ggplot2::ggplot() +
    ggplot2::geom_segment(data=ddata_y$segment,
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggplot2::geom_text(data=ddata_y$labels, ggplot2::aes(x=x, y=y, label=label),
              hjust=1, vjust=0.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits=xlim) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits=ylim) +
    theme_none() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(-0.1,0,0.6,0), 'cm')) +
    ggplot2::labs(x=NULL, y=NULL)

  # plot caption
  x_dist <- x_dist
  x_hclust <- x_hclust
  y_dist <- y_dist
  y_hclust <- y_hclust

  caption <- stringr::str_c(sprintf("X-Axis Dendrogram:\n     Distance method: %s\n     Linkage method: %s\n\n",x_dist, x_hclust),
                   sprintf("Y-Axis Dendrogram:\n     Distance method: %s\n     Linkage method: %s", y_dist, y_hclust))
  # putting plots together--------------------------------------------------------
  # extracting legend
  g_hl <- ggplot2::ggplotGrob(p_hl)
  legend <- grep("guide", g_hl$layout$name)
  g_l <- g_hl[["grobs"]][[legend]]

  # putting legend and captions together
  g_c <- grid::textGrob(label=caption, x=0.2, y=0.3, just='left')

  g_cl <- gridExtra::grid.arrange(g_c, g_l, ncol=1)

  # converting plots to grobs
  g_h <- ggplot2::ggplotGrob(p_h)
  g_x <- ggplot2::ggplotGrob(p_x)
  g_y <- ggplot2::ggplotGrob(p_y)

  mat <- matrix(list(g_x, g_h, g_cl, g_y), nrow=2)
  g <- gtable::gtable_matrix(name='heatmap',grobs=mat,
                     widths=ggplot2::unit(c(heatmap_w,dendro_y_w),'cm'),
                     heights=ggplot2::unit(c(dendro_x_h, heatmap_h), 'cm'))

  # saving plot
  ggplot2::ggsave(out, g, width=(heatmap_w+dendro_y_w)*1.01,
         height=(heatmap_h+dendro_x_h)*1.01, unit='cm')

  msg <- sprintf("\nSaving heatmap:\n%s", out)
  message(msg)

  return(list('pdata'=pdata, 'ddata_x'=ddata_x,
              'ddata_x'=ddata_x, 'p_hl'=p_hl,'p_x'=p_x,'p_y'=p_y,
              'g' = g))

}
