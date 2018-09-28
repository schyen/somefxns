#' plateAdjust
#'
#' normalized absorbance values by adjusting to blank
#'
#' @param plateDF dataframe. first spreadsheet of victor file
#' @param metadata dataframe. well metadata. must have columns:
#'     Well, platerow, platecol, strain, curveID, media, abx, abxlevel
#' @param well_include string or vector of string. Default NULL. Which wells to
#'     include in analysis
#' @param blank_by \code{c('row','col')}. Default \code{'row'}.
#'     Sets orientation of blanks.
#'     \code{'row'} means blanks are in the same row as corresponding wells
#'     \code{'col'} means blanks are in the same column as corresponding wells
#' @param blank_label string. Default \code{'blank'}.
#'     How are blanks identified in metadata
#'
#' @return
#'     full dataframe of with metadata, absorbance, adjusted absorbance in one dataframe
#'
plateAdjust <- function(plateDF, metadata, well_include=NULL, blank_by = 'row',
                        blank_label = 'blank') {

  if(!blank_by %in% c('row','col')) {
    stop("blank_by must be either 'row' or 'col'")
  }

  if (blank_by == 'row') blank_by <- 'platerow'
  else blank_by <- 'platecol'

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
  full <- rename(full, plateno = Plate)

  absAdjust <- function(d) {
    blankval <- d$abs[d$strain == blank_label]
    d$adj <- d$abs - blankval
    return(d)
  }

  # adjust absorbance using blank
  full <- full %>%
    group_by(!! sym(blank_by), Repeat) %>%
    do(absAdjust(.)) %>%
    as.data.frame()

  return(full)

}
