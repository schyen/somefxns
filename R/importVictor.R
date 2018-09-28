#' importVictor
#'
#' import raw data from victor
#'
#' @param file path and file name
#'
#' \code{df} - dataframe of wells, absorbance. is the first sheet of the victor file;
#' @import readxl
#' @export

importVictor <- function(file){

  sheet <- excel_sheets(file)

  # reading in read table
  raw <- as.data.frame(read_excel(file, sheet=sheet[1]))
  return(raw)

}
