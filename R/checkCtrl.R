#' checkCtrl
#'
#' check negative and positive control of MIC assay
#'
#' @param plateDF dataframe. first sheet of victor
#' @param metadata dataframe. well metadata must have columns:
#'     Well, platerow, platecol, strain, curveID, media, abx, abxlevel
#' @param well_include string or vector of string. Default NULL. Which wells to
#'     include in analysis
#' @param x string. x-axis of plot. default \code{'minute'}
#' @param y string. y-axis of plot. default \code{'abs'}
#' @param colour string. characteristic to colour by on plot. Default 'abxlevel'
#' @param group string. characteristic to colour by on plot. Default 'abxlevel'
#'
#' @import ggplot2
#' @return ggplot. growth curve of the negative and postive control

checkCtrl <- function(plateDF, metadata, well_include=NULL,
                      x='minute', y='abs', colour='abxlevel',
                      group='abxlevel') {

  newcolname = colnames(plateDF)
  newcolname[length(newcolname)] <- 'abs'
  colnames(plateDF) <- newcolname

  # adding metadata
  full <- merge(plateDF, metadata, 'Well')

  # converting time to minute
  full$minute <- (full$Repeat - 1) * 10
  full$hour <- full$minute/60

  # remove wells
  if(!is.null(well_include)) {
    full <- full[full$Well %in% well_include,]
  }

  # setting order of abx treatment
  factor_order <- rev(sort(unique(full[ ,group])))
  full[ ,group] <- factor(full[ ,group], levels=factor_order)

  p <- ggplot(full, aes_string(x=x, y=y, colour=colour)) +
    geom_point() +
    geom_vline(xintercept=0, colour='grey60') +
    geom_hline(yintercept=0, colour='grey60') +
    geom_smooth(aes_string(group=group)) +
    theme_bw(10)
  print(p)

  return(p)

}
