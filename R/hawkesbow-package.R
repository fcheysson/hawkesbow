#' @details
#' To be implemented later:
#' - Variance and confidence interval for the estimated parameters
#' - Spectral density based goodness-of-fit tests
#' - Custom built-kernels
#' @examples
#' # Simulate an exponential Hawkes process with baseline intensity 1,
#' # reproduction mean 0.5 and exponential fertility function with rate 2.
#' x <- hawkes(10, fun=1, repr=0.5, family="exp", rate=2)
#'
#' # Plot its conditional intensity function
#' oldpar = par()
#' par(mfrow = c(2, 1), mar = c(4.1, 4.1, 1.1, 2.1))
#' plot(x, intensity = TRUE)
#' # and its poisson cluster representation
#' plot(x, intensity = FALSE)
#' par(oldpar)
#'
#' # Estimate the parameters from the arrival times of `x`
#' # using maximum likelihood estimation
#' opt = mle(x$p, "Exponential", x$end)
#' opt$par                          # Estimated parameters
#' opt$model$ddloglik(x$p, x$end)     # Hessian matrix of the log-likelihood
#'
#' # Estimate the parameters from count data using Whittle's method
#' y = discrete(x, binsize = 1)
#' opt = whittle(y, "Exponential", binsize = 1)
#' opt$par                          # Estimated parameters
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib hawkesbow
#' @importFrom Rcpp sourceCpp loadModule
#' @importFrom graphics plot arrows axis legend matplot points segments
#' @importFrom methods is new
#' @importFrom stats fft integrate optim rexp rpois runif
#' @importFrom utils modifyList tail
## usethis namespace: end
NULL

loadModule("HawkesModule", TRUE)

.onUnload <- function (libpath) {
  library.dynam.unload("hawkesbow", libpath)
}
