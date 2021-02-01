#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib hawkesbow
#' @importFrom Rcpp sourceCpp loadModule
#' @importFrom graphics arrows axis legend matplot points segments
#' @importFrom methods is
#' @importFrom stats fft integrate optim rexp rpois runif
#' @importFrom utils modifyList tail
## usethis namespace: end
NULL

loadModule("HawkesModule", TRUE)

.onUnload <- function (libpath) {
  library.dynam.unload("hawkesbow", libpath)
}
