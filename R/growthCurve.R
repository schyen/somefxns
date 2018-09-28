#' growthCurve
#'
#' plotting growth curve of MIC plate assay
#'
#' @param pdata dataframe. MIC data with metadata, adjusted absorbance,
#'     and percent survival
#' @param well_include string or vector of string. Default NULL. Which wells to
#'     include in analysis
#' @param x string. Default \code{'minute'}. x axis of growth curve.
#' @param y string. Default \code{'adj'}. y axis of growth curve.
#' @param colour string. Default \code{'abxlevel'}. characteristic to colour code by.
#' @param group string. Default \code{'abxlevel'}. characteristic to group curves by
#' @param preview logic. Default \code{TRUE}. Preview plot
#'
#' @import ggplot2
#' @return growth curve. ggplot.
#' @export

growthCurve <- function(pdata, well_include=NULL, x='minute', y='adj',
                        colour='abxlevel', group='abxlevel', preview=TRUE) {

  # remove wells
  pdata <- pdata[pdata$Well %in% well_include,]

  # setting order of abx treatment
  factor_order <- sort(unique(pdata[ ,group]))
  pdata[ ,group] <- factor(pdata[ ,group], levels=factor_order)
  p <- ggplot(pdata, aes_string(x=x, y=y, colour=colour)) +
    # geom_point() +
    geom_vline(xintercept=0, colour='grey60') +
    geom_hline(yintercept=0, colour='grey60') +
    geom_smooth(aes_string(group=group)) +
    theme_bw(10)

  if(preview==TRUE) {
    print(p)
  }

  return(p)

}
