#' Title
#'
#' @param hawkes 
#'
#' @return
#' @export
#'
#' @examples
whittle <- function(data, binsize, init=c(1,.5,2), trunc=5, ...) {
  model <- new(ExpHawkes)
  model$ddata <- data
  model$binsize <- binsize
  
  n <- length(data)
  
  # Periodogram
  dft <- fft(data - mean(data))
  I <- Mod(dft)^2 / n

  # Whittle pseudo likelihood function (for optim)
  wlik <- function(param_) {
    model$param <- param_
    return( model$wLik(I, trunc) )
  }
  
  opt <- optim(par=init, fn = wlik, hessian = TRUE, 
               lower = rep(.0001, 3), upper = c(Inf, .9999, Inf), method = "L-BFGS-B", ...)
  
  model$vcov <- solve(opt$hessian)
  model$opt <- opt

  return( model )
}

#' Title
#'
#' @param trunc 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
whittle_ <- function(model, trunc=5, ...) {
  data <- model$ddata
  n <- length(data)
  
  # Periodogram
  dft <- fft(data - mean(data))
  I <- Mod(dft)^2 / n
  
  # Whittle pseudo likelihood function (for optim)
  wlik <- function(param_) {
    model$param <- param_
    return( model$wLik(I, trunc) )
  }
  
  opt <- optim(par=model$param, fn = wlik, hessian = TRUE, 
               lower = rep(.0001, 3), upper = c(Inf, .9999, Inf), method = "L-BFGS-B", ...)
  
  model$vcov <- solve(opt$hessian)
  model$opt <- opt
}
