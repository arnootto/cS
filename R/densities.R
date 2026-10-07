# =======================================================================
# Contaminated survival models
# Parameterization: (1 - alpha) * Reference(hazard h)  +  alpha * Delayed(hazard eta * h)
#   alpha in (0,1) : probability of belonging to the DELAYED component
#   eta   in (0,1) : hazard multiplier of that component
# No label swapping: the second component is always the low-hazard (delayed) one.
# =======================================================================

# ---- densities and survival functions ---------------------------------

#' Contaminated exponential (cE) distribution
#'
#' Density, survival function and hazard function of the contaminated
#' exponential distribution
#' \deqn{f_\text{cE}(t;\lambda,\alpha,\eta) = (1-\alpha)\lambda e^{-\lambda t}
#'   + \alpha\eta\lambda e^{-\eta\lambda t}, \quad t > 0,}
#' a mixture of a reference exponential with rate \eqn{\lambda} and a delayed
#' exponential with rate \eqn{\eta\lambda}.
#'
#' @param t      numeric vector of times \eqn{> 0}.
#' @param lambda rate of the reference component \eqn{> 0}.
#' @param alpha  probability of belonging to the delayed component
#'   \eqn{\in (0,1)}. May be a vector of the same length as \code{t}.
#' @param eta    hazard multiplier of the delayed component \eqn{\in (0,1)}.
#' @param log    logical; if \code{TRUE} the log-density is returned
#'   (\code{dcexp} only).
#'
#' @return Numeric vector of (log-)densities (\code{dcexp}), survival
#'   probabilities (\code{scexp}) or hazards (\code{hcexp}).
#' @seealso \code{\link{rcexp}}, \code{\link{ml.cexp}}
#' @export
#' @examples
#' t <- c(0.5, 1, 5)
#' dcexp(t, lambda = 1, alpha = 0.1, eta = 0.5)
#' scexp(t, lambda = 1, alpha = 0.1, eta = 0.5)
#' hcexp(t, lambda = 1, alpha = 0.1, eta = 0.5)
#'
#' # marginal hazard for different eta (Figure 1a of the paper)
#' t <- seq(0.001, 18, length = 500)
#' plot(t, hcexp(t, lambda = 1, alpha = 0.1, eta = 0.5), type = "l",
#'      ylim = c(0, 1.05), xlab = "time", ylab = "hazard")
#' lines(t, hcexp(t, lambda = 1, alpha = 0.1, eta = 0.2), col = 2)
#' abline(h = 1, lty = 3)
#' @importFrom stats dexp
dcexp <- function(t, lambda, alpha, eta, log = FALSE) {
  d <- (1 - alpha) * dexp(t, rate = lambda) +
    alpha * dexp(t, rate = eta * lambda)
  if (log) log(d) else d
}

#' @rdname dcexp
#' @export
scexp <- function(t, lambda, alpha, eta) {
  (1 - alpha) * exp(-lambda * t) + alpha * exp(-eta * lambda * t)
}

#' @rdname dcexp
#' @export
hcexp <- function(t, lambda, alpha, eta) {
  dcexp(t, lambda, alpha, eta) / scexp(t, lambda, alpha, eta)
}

# delayed component has hazard eta * h0, i.e. Weibull scale sigma * eta^(-1/k)
# with eta < 1 this scale is LARGER than sigma, so the times are longer

#' Contaminated Weibull (cW) distribution
#'
#' Density, survival function and hazard function of the contaminated Weibull
#' distribution
#' \deqn{f_\text{cW}(t;\sigma,k,\alpha,\eta) = (1-\alpha)f_\text{W}(t;\sigma,k)
#'   + \alpha f_\text{W}(t;\sigma\eta^{-1/k},k), \quad t > 0,}
#' where \eqn{f_\text{W}(t;\sigma,k)} is the Weibull density with scale
#' \eqn{\sigma} and shape \eqn{k}. The delayed component has hazard
#' \eqn{\eta h_\text{W}(t;\sigma,k)}, i.e. Weibull scale
#' \eqn{\sigma\eta^{-1/k}}, which is larger than \eqn{\sigma} when
#' \eqn{\eta < 1}. \code{h0weib} is the hazard of the reference Weibull,
#' \eqn{h_\text{W}(t;\sigma,k) = (k/\sigma)(t/\sigma)^{k-1}}.
#'
#' @param t     numeric vector of times \eqn{> 0}.
#' @param sigma scale of the reference component \eqn{> 0}.
#' @param k     shape \eqn{> 0}, common to both components.
#' @param alpha probability of belonging to the delayed component
#'   \eqn{\in (0,1)}. May be a vector of the same length as \code{t}.
#' @param eta   hazard multiplier of the delayed component \eqn{\in (0,1)}.
#' @param log   logical; if \code{TRUE} the log-density is returned
#'   (\code{dcweib} only).
#'
#' @return Numeric vector of (log-)densities (\code{dcweib}), survival
#'   probabilities (\code{scweib}) or hazards (\code{hcweib},
#'   \code{h0weib}).
#' @seealso \code{\link{rcweib}}, \code{\link{ml.cweib}}
#' @export
#' @examples
#' t <- c(1, 5, 20)
#' dcweib(t, sigma = 5, k = 1.5, alpha = 0.1, eta = 0.5)
#' scweib(t, sigma = 5, k = 1.5, alpha = 0.1, eta = 0.5)
#' hcweib(t, sigma = 5, k = 1.5, alpha = 0.1, eta = 0.5)
#'
#' # marginal hazard against the reference hazard (Figure 2 of the paper)
#' t <- seq(0.001, 36, length = 500)
#' plot(t, h0weib(t, sigma = 5, k = 1.5), type = "l", lty = 3,
#'      xlab = "time", ylab = "hazard")
#' lines(t, hcweib(t, sigma = 5, k = 1.5, alpha = 0.1, eta = 0.5), col = 2)
#' @importFrom stats dweibull
dcweib <- function(t, sigma, k, alpha, eta, log = FALSE) {
  d <- (1 - alpha) * dweibull(t, shape = k, scale = sigma) +
    alpha * dweibull(t, shape = k, scale = sigma * eta^(-1 / k))
  if (log) log(d) else d
}

#' @rdname dcweib
#' @export
scweib <- function(t, sigma, k, alpha, eta) {
  H <- (t / sigma)^k
  (1 - alpha) * exp(-H) + alpha * exp(-eta * H)
}

#' @rdname dcweib
#' @export
hcweib <- function(t, sigma, k, alpha, eta) {
  dcweib(t, sigma, k, alpha, eta) / scweib(t, sigma, k, alpha, eta)
}

#' @rdname dcweib
#' @export
h0weib <- function(t, sigma, k) (k / sigma) * (t / sigma)^(k - 1)

# proportion of the delayed component among those still at risk at time t

#' Proportion of delayed observations still at risk (cE)
#'
#' Proportion of the delayed component among the observations that have not
#' yet experienced the event by time \eqn{t}, under the contaminated
#' exponential distribution,
#' \deqn{\pi(t) = \frac{\alpha e^{-\eta\lambda t}}{S_\text{cE}(t;\lambda,\alpha,\eta)}.}
#' It equals \eqn{\alpha} at \eqn{t = 0} and tends to \eqn{1} as
#' \eqn{t \to \infty}.
#'
#' @inheritParams dcexp
#' @return Numeric vector of proportions in \eqn{(0,1)}.
#' @export
#' @examples
#' pi.cexp(c(0, 1, 5, 20), lambda = 1, alpha = 0.1, eta = 0.5)
pi.cexp <- function(t, lambda, alpha, eta) {
  alpha * exp(-eta * lambda * t) / scexp(t, lambda, alpha, eta)
}
