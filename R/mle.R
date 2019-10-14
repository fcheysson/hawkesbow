#' Title
#'
#' @param data 
#' @param init 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
mle <- function(data, T, init=c(1,1,2), opts = NULL, ...) {
  model <- new(ExpHawkes)
  model$data <- data
  model$T <- T
  
  # Likelihood function (for nloptr)
  nlopt_fn <- function(param_) {
    model$param <- param_
    return( lapply(X = model$likngrad(), FUN = "-") )
  }
  
  if (is.null(opts))
    opts <- list("algorithm" = "NLOPT_LD_LBFGS") else
      if (is.null(opts[["algorithm"]]))
        opts <- c(opts, "algorithm" = "NLOPT_LD_LBFGS")
  
  opt <- nloptr(x0 = init, eval_f = nlopt_fn, lb = rep(.0001, 3), opts = opts, ...)
  
  model$vcov <- solve(-model$hessian())
  model$opt <- opt
  
  return( model )
}

#' Title
#'
#' @param model 
#' @param opts 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
mle_ <- function(model, opts = NULL, ...) {
  # Likelihood function (for nloptr)
  nlopt_fn <- function(param_) {
    model$param <- param_
    return( lapply(X = model$likngrad(), FUN = "-") )
  }
  
  if (is.null(opts))
    opts <- list("algorithm" = "NLOPT_LD_LBFGS") else
      if (is.null(opts[["algorithm"]]))
        opts <- c(opts, "algorithm" = "NLOPT_LD_LBFGS")
  
  opt <- nloptr(x0 = model$param, eval_f = nlopt_fn, lb = rep(.0001, 3), opts = opts, ...)
  
  model$vcov <- solve(-model$hessian())
  model$opt <- opt
}