#' plateAdjust
#'
#' normalized absorbance values by adjusting to blank
#'
#' @param plateDF dataframe. first spreadsheet of victor file
#' @param metadata dataframe. well metadata. must have columns:
#'     Well, platerow, platecol, strain, curveID, welltype, media, abx, wellconc
#' @param well_include string or vector of string. Default NULL. Which wells to
#'     include in analysis
#' @param blank_by \code{c('row','col', 'T0')}. Default \code{'row'}.
#'     Sets orientation of blanks.
#'     \code{'row'} means blanks are in the same row as corresponding wells
#'     \code{'col'} means blanks are in the same column as corresponding wells
#'     \code{'timepoint'} means blank by reading at time 0 (aka repeat 1)
#' @param blank_label string. Default \code{'blank'}.
#'     How are blanks identified in metadata
#'
#' @import dplyr
#' @return
#'     full dataframe of with metadata, absorbance, adjusted absorbance in one dataframe
#' @export

plateAdjust <- function(plateDF, metadata, well_include=NULL, blank_by = 'row',
                        blank_label = 'blank') {

  if(!blank_by %in% c('row','col','T0')) {
    stop("blank_by must be either 'row', 'col' or 'T0'")
  }

  if (blank_by == 'row') blank_by <- 'platerow'
  if (blank_by== 'col') blank_by <- 'platecol'

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

  if(blank_by != 'T0') {
    absAdjust <- function(d) {
      blankval <- d$abs[d$welltype == blank_label]
      d$adj <- d$abs - blankval
      return(d)
    }

    # adjust absorbance using blank
    full <- full %>%
      group_by(!! sym(blank_by), Repeat) %>%
      do(absAdjust(.)) %>%
      as.data.frame()
  }
  else {
    absAdjust <- function(d) {
      blankval <- d$abs[d$Repeat == 1]
      d$adj <- d$abs - blankval
      return(d)
    }

    full <- full %>%
      group_by(Well) %>%
      arrange(Repeat) %>%
      do(absAdjust(.)) %>%
      as.data.frame()

    # remove timepoint 0
    full <- filter(full, Repeat != 1)
  }

  return(full)

}
