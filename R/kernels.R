#' @name Model
#' @title C++ abstract class for Hawkes processes
#'
#' @description This is a C++ abstract class for Hawkes processes, which holds
#' methods for the estimation of its parameters.
#'
#' @field param Vector of parameters of the Hawkes process, of the form \eqn{(\eta, \mu, ...)}.
#' @field binsize Bin size for the count sequences.
#' @field new(DerivedClass,(param),(binsize)) Constructor for derived classes; `param` and/or `binsize`
#' can be safely ignored.
#' @field mean() Returns the expected value on \eqn{[0,\mathrm{end}]}.
#' @field dmean() Returns the Jacobian matrix of the expected value on \eqn{[0,\mathrm{end}]}.
#' @field ddmean() Returns the Hessian matrix of the expected value on \eqn{[0,\mathrm{end}]}.
#' @field f(xi) Returns the spectral density function of the time-continuous count sequence. \itemize{
#' \item `xi` A numeric vector of frequencies.
#' }
#' @field f1(xi,trunc) Returns the spectral density function of the discrete time count sequence. \itemize{
#' \item `xi` A numeric vector of frequencies.
#' \item `trunc` The number of foldings to take into account for the aliasing.
#' }
#' @field whittle(I,trunc) Returns the log-spectral likelihood of a discrete time count sequence. \itemize{
#' \item `I` The periodogram of the count sequence.
#' \item `trunc` The number of foldings to take into account for the aliasing.
#' }
#' @field loglik(events,end) Returns the log-likelihood of a sequence of arrival times. \itemize{
#' \item `events` The sequence of arrival times.
#' \item `end` The endpoint of the observation window \eqn{[0,\mathrm{end}]}.
#' }
#' @field dloglik(events,end) Returns the Jacobian matrix of the log-likelihood of a sequence of arrival times. \itemize{
#' \item `events` The sequence of arrival times.
#' \item `end` The endpoint of the observation window \eqn{[0,\mathrm{end}]}.
#' }
#' @field ddloglik(events,end) Returns the Hessian matrix of the log-likelihood of a sequence of arrival times. \itemize{
#' \item `events` The sequence of arrival times.
#' \item `end` The endpoint of the observation window \eqn{[0,\mathrm{end}]}.
#' }
#'
#' @details This serves as a basis for the Hawkes model and its count sequence,
#' with conditional intensity function
#' \deqn{\lambda(t) = \eta + \mu \sum_{T_i < t} h^\ast(t - T_i).}
#' As an abstract class, an object of class `Model` should never be directly
#' instanciated, but rather one of its derived class.
#' The constructor can take no argument, in which case the vector `param` is
#' initialised to sensible values and `binsize` defaults to 1.
#' Alternatively, `param` and/or `binsize` can be specified.
#'
#' @examples
#' # Simulate 1000 exponential Hawkes processes on \eqn{[0, 100]},
#' # and average the periodogram of the count sequences with bin size 1
#' # at each frequency.
#' I = rep(0, 100)
#' for (k in 1:1e3) {
#'     x = hawkes(100, fun = 1, repr = .5, family = "exp", rate = 2)
#'     y = discrete(x, binsize = 1)
#'     I = I + Mod(fft(y - mean(y)))^2 / length(y)
#' }
#'
#' # Check that the averaged periodogram correctly approximates the spectral
#' # density function of the count sequence
#' model = new(Exponential)
#' model$param = c(1, .5, 2)
#' model$binsize = 1
#'
#' z = 2 * pi * 0:99 / 100     # Frequencies of the periodogram
#' plot(z, I / 1e3, type = "l") # Averaged periodogram
#' lines(z, model$f1(xi = z, trunc = 10L), col = "red")
#'
#' @seealso [Exponential]
NULL

#' @name Exponential
#' @title Reproduction kernels for the Hawkes processes
#'
#' @description These classes are derived from the class `Model`, each implementing
#' a different reproduction kernel for the Hawkes process.
#' They inherit all fields from [Model].
#'
#' @details
#' \itemize{
#' \item The kernel `Exponential` has density function
#' \deqn{h^\ast(t) = \beta \exp(-\beta t) 1_{\{t \ge 0\}}.}
#' Its vector of parameters must be of the form \eqn{(\eta, \mu, \beta)}.
#' Both `loglik`, its derivatives, and `whittle` can be used with this reproduction kernel.
#'
#' \item The kernel `SymmetricExponential` has density function
#' \deqn{h^\ast(t) = 0.5 \beta \exp(-\beta |t|).}
#' Its vector of parameters must be of the form \eqn{(\eta, \mu, \beta)}.
#' Only `whittle` can be used with this reproduction kernel.
#'
#' \item The kernel `Gaussian` has density function
#' \deqn{h^\ast(t) = \frac{1}{\sigma \sqrt{2\pi}}\exp\left(-\frac{(t-\nu)^2}{2\sigma^2}\right).}
#' Its vector of parameters must be of the form \eqn{(\eta, \mu, \nu, \sigma^2)}.
#' Only `whittle` is available with this reproduction kernel.
#'
#' \item The kernel `PowerLaw` has density function
#' \deqn{h^\ast(t) = \theta a^\theta (t+a)^{-\theta-1} 1_{\{\theta > 0 \}}.}
#' Its vector of parameters must be of the form \eqn{(\eta, \mu, \theta, a)}.
#' Both `loglik`, its derivatives, and `whittle` can be used with this reproduction kernel.
#'
#' \item The kernels `Pareto3`, `Pareto2` and `Pareto1` have density function
#' \deqn{h_\theta^\ast(t) = \theta a^\theta t^{-\theta - 1} 1_{\{t > a\}},}
#' with \eqn{\theta} = 3, 2 and 1 respectively.
#' Their vectors of parameters must be of the form \eqn{(\eta, \mu, a)}.
#' Only `whittle` is available with this reproduction kernel.
#' }
#'
#' @seealso [Model]
NULL

#' @name SymmetricExponential
#' @rdname Exponential
NULL

#' @name Gaussian
#' @rdname Exponential
NULL

#' @name PowerLaw
#' @rdname Exponential
NULL

#' @name Pareto3
#' @rdname Exponential
NULL

#' @name Pareto2
#' @rdname Exponential
NULL

#' @name Pareto1
#' @rdname Exponential
NULL
