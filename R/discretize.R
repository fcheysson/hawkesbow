#' Discretizes a Hawkes simulation
#'
#' @param hawkes An object created by the function `hawkes`
#' @param length (Either) The length for the output vector
#' @param binsize (Either) The binsize for the discretization
#'
#' @return
#' The vector of counts
#' @export
#'
#' @examples
#' x = hawkes(100, fun=1, repr=0.5, family="exp", rate=2)
#' y = discrete(x, length=100)
#' z = discrete(x, binsize=1)
#' all(y == z)
discrete <- function(hawkes, length=NULL, binsize=NULL) {
  if (is.null(length) & is.null(binsize))
    stop("One of length or binsize must be specified.")
  if (!is.null(binsize) & !is.null(length))
    message("Warning: Both length and binsize specified. binsize will be ignored.\n")
  if (!is.null(length)) {
    if (length <= 0)
      stop("length is less than or equal to zero.")
    ti <- seq(0, hawkes$end, length.out=length+1)
    bin <- sapply(1:length, function(i) {
      sum(hawkes$p > ti[i] & hawkes$p < ti[i+1])
    })
  } else {
    if (hawkes$end %% binsize != 0)
      message("Warning: hawkes$end is not a multiplier of binsize. Last bin will be discarded.\n")
    ti <- seq(0, hawkes$end, by=binsize)
    bin <- sapply(1:(length(ti)-1), function(i) {
      sum(hawkes$p > ti[i] & hawkes$p < ti[i+1])
    })
  }
  return(bin)
}
