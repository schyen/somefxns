#' importVictor
#'
#' import raw data from victor
#'
#' @param file path and file name
#'
#' @return list of:
#'     df - dataframe of wells, absorbance. is the first sheet of the victor file
#'     plate - 2nd sheet of victor file
#'     protocol - 3rd sheet of victor file
#'     errors - 4th sheet of victor file
#'     notes - 5th sheet of victor file
#'

importVictor <- function(file){

  sheet <- excel_sheets(file)

  # reading in read table
  raw <- list('df' = as.data.frame(read_excel(file, sheet=sheet[1])) ,
              'plate' = read_excel(file, sheet=sheet[2]),
              'protocol' = read_excel(file, sheet=sheet[3]),
              'errors' = read_excel(file, sheet=sheet[4]),
              'notes' = read_excel(file, sheet=sheet[5])
  )
  return(raw)

}
