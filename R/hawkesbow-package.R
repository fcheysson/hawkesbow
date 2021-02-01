#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib hawkesbow
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp loadModule
## usethis namespace: end
NULL

loadModule("HawkesModule", TRUE)

.onUnload <- function (libpath) {
  library.dynam.unload("hawkesbow", libpath)
}
