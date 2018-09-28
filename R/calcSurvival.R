#' calcSurvival
#'
#' Calculate percent survival using growth control well as 100%
#'
#' @param full dataframe. metadata, absorbance, adjusted absorbance
#' @param ctrl_by \code{c('row', 'col')} Default 'row'
#'     \code{'row'} means growth control wells are in the same row as corresponding wells
#'     \code{'col'} means growth control wells are in the same column as corresponding wells
#'
#' @import dplyr
#' @import rlang

#' @return dataframe. full, with percent survival in survival column.

calcSurvival <- function (full, ctrl_by = 'row') {

  if(!blank_by %in% c('row','col')) {
    stop("blank_by must be either 'row' or 'col'")
  }

  if (ctrl_by == 'row') ctrl_by <- 'platerow'
  else ctrl_by <- 'platecol'

  # calculating percent survival
  calcSurvival <- function(d) {
    growthctrl <- abs(d$adj[d$welltype == 'growthctrl'])

    # only calculate survival if growth control well exists in current df
    if (!all(is.na(growthctrl))) {
      growthctrl <- growthctrl[!is.na(growthctrl)]
      d$survival <- round(abs(d$adj),3) / round(growthctrl,3) * 100
    }
    else {
      d$survival <- NA
    }

    return(d)
  }

  # calculate percent survival
  full <- full %>%
    group_by(!! sym(ctrl_by), Repeat) %>%
    do(calcSurvival(.)) %>%
    as.data.frame()

}
