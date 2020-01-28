#' Fitting Hawkes processes from discrete data
#'
#' This function fits a Hawkes process to discrete data by minimizing the Whittle contrast.
#'
#' @param model An object of class Model
#' @param data An object of class DiscreteData
#' @param trunc (Optional) The number of foldings taken into account due to aliasing
#' @param ... Additional arguments passed to `optim`
#'
#' @return
#' Returns the `optim` output
#' @export
#'
#' @examples
#' x = hawkes(100, fun = 1, repr = .5, family = "exp", rate = 1)
#' y = discrete(x, length = 100)
#' model = new(Exponential)
#' data = new(DiscreteData, y, 1)
#' whittle(model, data)
whittle <- function(model, data, trunc=5L, ...) {

    model$attach(data)
    counts <- data$counts
    n <- length(counts)

    # Periodogram
    dft <- fft(counts - mean(counts))
    I <- Mod(dft)^2 / n
    I <- I[-1]      # remove value corresponding to density in 0 (not well defined for centered processes)

    # Whittle pseudo likelihood function (for optim)
    wlik <- function(param) {
        model$param <- param
        return( model$whittle(I, trunc) )
    }

    optargs = list(hessian = TRUE,
                   lower = rep(.0001, length(model$param)),
                   upper = c(Inf, .9999, rep(Inf, length(model$param)-2)),
                   method = "L-BFGS-B")
    optargs = modifyList(optargs, list(...))
    opt <- do.call(optim, c(list(par = model$param, fn = wlik), optargs))

    # model$param <- opt$
    return(opt);
}
