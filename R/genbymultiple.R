#' genbymultiple
#'
#' generate sequence by multiples
#'
#' @param x numeric; starting value
#' @param n numeric; number of values to be in sequence
#' @param multiple numeric; factor by which sequence is generated
#'
#' @return vector of numbers increasing by a multiple
#' @export

genbymultiple <- function(x, n, multiple) {
  out <- c()
  curr <- x
  counter <- 1
  while(counter <= n) {
    out <- c(out, curr)
    curr <- curr*multiple
    counter <- counter + 1
  }
  return(out)
}
