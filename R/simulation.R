# ---- simulation (z = 1 means delayed component) -----------------------

#' Random generation from the cE and cW distributions
#'
#' Generates random times from the contaminated exponential
#' (\code{rcexp}) or contaminated Weibull (\code{rcweib}) distribution,
#' together with the component labels. Each observation is first assigned to
#' the delayed component with probability \code{alpha} (\eqn{z = 1}) and its
#' time is then drawn from the corresponding component.
#'
#' @param n      number of observations.
#' @param lambda rate of the reference exponential component \eqn{> 0}.
#' @param sigma  scale of the reference Weibull component \eqn{> 0}.
#' @param k      Weibull shape \eqn{> 0}.
#' @param alpha  probability of belonging to the delayed component
#'   \eqn{\in (0,1)}.
#' @param eta    hazard multiplier of the delayed component \eqn{\in (0,1)}.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{t}}{numeric vector of simulated times.}
#'     \item{\code{z}}{integer vector of component labels: \code{1} for the
#'       delayed component, \code{0} for the reference component.}
#'   }
#' @seealso \code{\link{dcexp}}, \code{\link{dcweib}},
#'   \code{\link{rcexp.cv}}
#' @export
#' @examples
#' set.seed(1)
#' d1 <- rcexp(1000, lambda = 0.5, alpha = 0.15, eta = 0.3)
#' mean(d1$z)
#' d2 <- rcweib(1000, sigma = 10, k = 1.5, alpha = 0.15, eta = 0.25)
#' tapply(d2$t, d2$z, median)
#' @importFrom stats rbinom rexp
rcexp <- function(n, lambda, alpha, eta) {
  z <- rbinom(n, 1, alpha)
  t <- rexp(n, rate = ifelse(z == 1, eta * lambda, lambda))
  list(t = t, z = z)
}

#' @rdname rcexp
#' @export
#' @importFrom stats rweibull
rcweib <- function(n, sigma, k, alpha, eta) {
  z <- rbinom(n, 1, alpha)
  t <- rweibull(n, shape = k, scale = ifelse(z == 1, sigma * eta^(-1 / k), sigma))
  list(t = t, z = z)
}

# ---- simulation from the concomitant models ---------------------------

#' Random generation from the cE and cW models with concomitant variables
#'
#' Generates random times from the contaminated exponential
#' (\code{rcexp.cv}) or contaminated Weibull (\code{rcweib.cv}) model in which
#' the probability of belonging to the delayed component depends on
#' concomitant variables through a logit link,
#' \deqn{\alpha_i = \frac{\exp(\mathbf{v}_i^\top\boldsymbol\gamma)}
#'   {1 + \exp(\mathbf{v}_i^\top\boldsymbol\gamma)}.}
#'
#' @param n      number of observations.
#' @param lambda rate of the reference exponential component \eqn{> 0}.
#' @param sigma  scale of the reference Weibull component \eqn{> 0}.
#' @param k      Weibull shape \eqn{> 0}.
#' @param eta    hazard multiplier of the delayed component \eqn{\in (0,1)}.
#' @param gamma  numeric vector of regression coefficients for the mixing
#'   weight, in the order of the columns of the model matrix of
#'   \code{alpha.formula} (intercept first).
#' @param covars a data frame with \code{n} rows containing the concomitant
#'   variables. If \code{NULL} (default), a data frame with an \code{id}
#'   column is created and only an intercept can be used.
#' @param alpha.formula a one-sided \code{\link[stats]{formula}} for the
#'   mixing weight, e.g. \code{~ x1 + x2}. Default \code{~1}.
#'
#' @return A data frame with columns \code{t} (simulated times), \code{z}
#'   (\code{1} = delayed component), \code{alpha} (the mixing weight of each
#'   observation) and the columns of \code{covars}.
#' @seealso \code{\link{ml.cexp.cv}}, \code{\link{ml.cweib.cv}}
#' @export
#' @examples
#' set.seed(1)
#' n  <- 1000
#' cv <- data.frame(rural = rbinom(n, 1, 0.4), indigenous = rbinom(n, 1, 0.1))
#' g  <- c(log(0.15/0.85), log(2.5), log(3))
#' d1 <- rcexp.cv(n, lambda = 0.1, eta = 0.3, gamma = g,
#'                covars = cv, alpha.formula = ~ rural + indigenous)
#' head(d1)
#' @importFrom stats model.matrix rbinom rexp
rcexp.cv <- function(n, lambda, eta, gamma, covars = NULL, alpha.formula = ~1) {
  if (is.null(covars)) covars <- data.frame(id = 1:n)
  V <- model.matrix(alpha.formula, data = covars)
  alpha <- as.vector(exp(V %*% gamma)/(1 + exp(V %*% gamma)))
  z <- rbinom(n, 1, alpha)
  t <- rexp(n, rate = ifelse(z == 1, eta*lambda, lambda))
  data.frame(t = t, z = z, alpha = alpha, covars)
}

#' @rdname rcexp.cv
#' @export
#' @importFrom stats rweibull
rcweib.cv <- function(n, sigma, k, eta, gamma, covars = NULL, alpha.formula = ~1) {
  if (is.null(covars)) covars <- data.frame(id = 1:n)
  V <- model.matrix(alpha.formula, data = covars)
  alpha <- as.vector(exp(V %*% gamma)/(1 + exp(V %*% gamma)))
  z <- rbinom(n, 1, alpha)
  t <- rweibull(n, shape = k,
                scale = ifelse(z == 1, sigma*eta^(-1/k), sigma))
  data.frame(t = t, z = z, alpha = alpha, covars)
}
