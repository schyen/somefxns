#' survivalCurve
#'
#' plotting survival curve of MIC plate assay
#'
#' @param pdata dataframe. MIC data with metadata, adjusted absorbance,
#'     and percent survival
#' @param well_include string or vector of string. Default NULL. Which wells to
#'     include in analysis
#' @param x string. Default \code{'minute'}. x axis of growth curve.
#' @param y string. Default \code{'survival'}. y axis of growth curve.
#' @param timepoint Default \code{99}. Which timepoint in the growth curve
#'     should the survival curve be made from. Should correspond to the
#'     \code{Repeat} value in the reading. If multiple timepoints provided,
#'     produces survival curve for each timepoint.
#' @param preview logic. Defulat \code{TRUE}. Preview plot
#'
#' @import ggplot2
#' @import dplyr
#' @return growth curve. ggplot.
#' @export

survivalCurve <- function(pdata, well_include=NULL, timepoint=max(pdata$Repeat),
                          x='wellconc', y='survival', preview=TRUE) {

  scaleFUN <- function(x) signif(x, 4)

  # remove wells
  if(!is.null(well_include)) {
    pdata <- pdata[pdata$Well %in% well_include,]
  }

  pdata[,x] <- as.numeric(pdata[,x])
  pdata <- pdata[pdata$Repeat %in% timepoint,]
  pdata$time_h <- pdata$minute/60

  pdata$label <- sprintf("Incub. time: %0.2f h",pdata$time_h)
  pdata[,x] <- as.numeric(pdata[,x])
  pdata <- pdata %>%
    arrange(!! sym(x))

  y_upper <- 10*ceiling(max(pdata$survival[!is.na(pdata$survival)])/10)

  p <- ggplot(pdata, aes_string(x=x, y=y)) +
    geom_vline(xintercept=0, colour='grey60') +
    geom_hline(yintercept=0, colour='grey60') +
    # geom_hline(yintercept=10, colour='red', linetype='dashed') +
    # geom_smooth() +
    geom_point() +
    scale_y_continuous(breaks=seq(y_upper,0,-20)) +
    scale_x_continuous(trans='log', breaks=unique(pdata[,x]),
                       labels = scaleFUN) +
    theme_bw(10)

  if(length(timepoint)>1) {
    p <- p + facet_wrap(~label)
  }
  else {
    p <- p + labs(subtitle=unique(pdata$label))
  }

  if(preview==TRUE) {
    print(p)
  }

  return(list('pdata' = pdata, 'p' = p))
}
