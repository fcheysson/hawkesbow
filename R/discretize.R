#' discrete
#'
#' @param hawkes 
#' @param length 
#' @param binsize 
#'
#' @return
#' @export
#'
#' @examples
discrete <- function(hawkes, length=NULL, binsize=NULL) {
  if (is.null(length) & is.null(binsize))
    stop("One of length or binsize must be specified.")
  if (!is.null(binsize) & !is.null(length))
    cat("Warning: Both length and binsize specified. binsize will be ignored.\n")
  if (!is.null(length)) {
    if (length <= 0)
      stop("length is less than or equal to zero.")
    ti <- seq(0, hawkes$T, length.out=length+1)
    bin <- sapply(1:length, function(i) {
      sum(hawkes$p > ti[i] & hawkes$p < ti[i+1])
    })
  } else {
    if (hawkes$T %% binsize != 0)
      cat("Warning: hawkes$T is not a multiplier of binsize. Last bin will be discarded.\n")
    ti <- seq(0, hawkes$T, by=binsize)
    bin <- sapply(1:(length(ti)-1), function(i) {
      sum(hawkes$p > ti[i] & hawkes$p < ti[i+1])
    })
  }
  return(bin)
}