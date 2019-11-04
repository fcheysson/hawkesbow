#' Fitting Hawkes processes from continuous data
#'
#' This function fits a Hawkes process to continuous data by minimizing the likelihood.
#'
#' @param model An object of class Model
#' @param data An object of class ContinuousData
#' @param opts (Optional) To be passed to `nloptr`
#' @param ... Additional arguments passed to `nloptr`
#'
#' @return
#' Returns the `optim` output
#' @export
#'
#' @examples
#' x = hawkes(100,1,1,2)
#' model = new(Exponential)
#' data = new(ContinuousData, x$p, x$T)
#' mle(model, data)
mle <- function(model, data, opts = NULL, ...) {
  # Likelihood function (for nloptr)
  nlopt_fn <- function(param) {
    model$param <- param
    return( lapply(X = model$likngrad(), FUN = "-") )
  }

  if (is.null(opts))
    opts <- list("algorithm" = "NLOPT_LD_LBFGS") else
      if (is.null(opts[["algorithm"]]))
        opts <- c(opts, "algorithm" = "NLOPT_LD_LBFGS")

  opt <- nloptr(x0 = model$param, eval_f = nlopt_fn, lb = rep(.0001, 3), ...)

  model$vcov <- solve(-model$hessian())
  model$opt <- opt
}
