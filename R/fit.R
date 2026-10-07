# ---- weighted complete-data log-likelihoods (v = prob. delayed) -------
# Internal helpers (not exported), used in the M-step of ml.cexp and ml.cweib.

ll.cexp <- function(t, lambda, eta, v) {
  (1 - v) * dexp(t, rate = lambda, log = TRUE) +
    v * dexp(t, rate = eta * lambda, log = TRUE)
}

ll.cweib <- function(t, sigma, k, eta, v) {
  (1 - v) * dweibull(t, shape = k, scale = sigma, log = TRUE) +
    v * dweibull(t, shape = k, scale = sigma * eta^(-1 / k), log = TRUE)
}

# ---- reference (one-component) survival models ------------------------

#' Fit the exponential reference model
#'
#' Fits the one-component exponential distribution by maximum likelihood. The
#' MLE is \eqn{\hat\lambda = 1/\bar{t}}.
#'
#' @details
#' **Standard errors.** Obtained by inverting the negative Hessian of the
#' log-likelihood at the MLE, computed numerically with
#' \code{\link[stats]{optimHess}}.
#'
#' @param t numeric vector of observed times \eqn{> 0}.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{results}}{data frame with the estimate and standard error
#'       of \eqn{\lambda}.}
#'     \item{\code{lambda}}{numeric; estimated rate \eqn{\hat\lambda}.}
#'     \item{\code{loglike}}{numeric; maximised log-likelihood.}
#'     \item{\code{AIC}}{numeric; Akaike information criterion.}
#'     \item{\code{BIC}}{numeric; Bayesian information criterion.}
#'   }
#' @seealso \code{\link{fit.w}}, \code{\link{ml.cexp}}
#' @export
#' @examples
#' set.seed(1)
#' d <- rcexp(1000, lambda = 0.5, alpha = 0.15, eta = 0.3)
#' fit.e(d$t)$results
#' @importFrom stats dexp optimHess
fit.e <- function(t) {
  n <- length(t)

  lambda <- 1 / mean(t)
  ll <- sum(dexp(t, rate = lambda, log = TRUE))

  obs.ll <- function(par, t) sum(dexp(t, rate = par[1], log = TRUE))
  Hess <- optimHess(par = lambda, fn = obs.ll, t = t)
  se <- sqrt(diag(solve(-Hess)))

  results <- data.frame(estimates = lambda, se = se)
  rownames(results) <- "lambda"

  list(results = results,
       lambda = lambda,
       loglike = ll,
       AIC = -2 * ll + 1 * 2,
       BIC = -2 * ll + 1 * log(n))
}

#' Fit the Weibull reference model
#'
#' Fits the one-component Weibull distribution with scale \eqn{\sigma} and
#' shape \eqn{k} by maximum likelihood.
#'
#' @details
#' **Estimation.** For fixed \eqn{k} the MLE of the scale is
#' \eqn{\hat\sigma(k) = \{\overline{t^k}\}^{1/k}}, so the shape is found by
#' maximising the profile log-likelihood over \eqn{k \in (0.05, 50)} with
#' \code{\link[stats]{optimize}}.
#'
#' **Standard errors.** Obtained by inverting the negative Hessian of the
#' log-likelihood at the MLE, computed numerically with
#' \code{\link[stats]{optimHess}}.
#'
#' @param t numeric vector of observed times \eqn{> 0}.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{results}}{data frame with the estimates and standard
#'       errors of \eqn{\sigma} and \eqn{k}.}
#'     \item{\code{sigma}, \code{k}}{numeric; estimated scale and shape.}
#'     \item{\code{loglike}}{numeric; maximised log-likelihood.}
#'     \item{\code{AIC}}{numeric; Akaike information criterion.}
#'     \item{\code{BIC}}{numeric; Bayesian information criterion.}
#'   }
#' @seealso \code{\link{fit.e}}, \code{\link{ml.cweib}}
#' @export
#' @examples
#' set.seed(1)
#' d <- rcweib(1000, sigma = 10, k = 1.5, alpha = 0.15, eta = 0.25)
#' fit.w(d$t)$results
#' @importFrom stats dweibull optimize optimHess
fit.w <- function(t) {
  n <- length(t)

  pl <- function(k) {                       # profile log-likelihood in k
    sigma <- mean(t^k)^(1 / k)
    sum(dweibull(t, shape = k, scale = sigma, log = TRUE))
  }
  k     <- optimize(pl, c(0.05, 50), maximum = TRUE)$maximum
  sigma <- mean(t^k)^(1 / k)
  ll    <- pl(k)

  obs.ll <- function(par, t) sum(dweibull(t, shape = par[2], scale = par[1], log = TRUE))
  Hess <- optimHess(par = c(sigma, k), fn = obs.ll, t = t)
  se <- sqrt(diag(solve(-Hess)))

  results <- data.frame(estimates = c(sigma, k), se = se)
  rownames(results) <- c("sigma", "k")

  list(results = results,
       sigma = sigma, k = k,
       loglike = ll,
       AIC = -2 * ll + 2 * 2,
       BIC = -2 * ll + 2 * log(n))
}

# ---- contaminated exponential -----------------------------------------

#' Contaminated exponential (cE) model: EM estimation
#'
#' Fits the contaminated exponential distribution
#' \eqn{(1-\alpha)\,\mathcal{E}(\lambda) + \alpha\,\mathcal{E}(\eta\lambda)}
#' by maximum likelihood using the EM algorithm.
#'
#' @details
#' **E-step.** Computes the posterior probability that each observation
#' belongs to the delayed component,
#' \deqn{v_i = \frac{\alpha f_\text{E}(t_i;\eta\lambda)}
#'   {f_\text{cE}(t_i;\lambda,\alpha,\eta)}.}
#'
#' **M-step.** \eqn{\alpha} is updated in closed form as \eqn{\bar{v}};
#' \eqn{\lambda} and \eqn{\eta} are updated by maximising the weighted
#' complete-data log-likelihood with \code{\link[stats]{optim}}, using
#' \eqn{\log\lambda} and \eqn{\mathrm{logit}(\eta)} so that \eqn{\eta}
#' stays in \eqn{(0,1)}.
#'
#' **Convergence.** The algorithm stops when the change in the observed
#' log-likelihood is below \code{em.tol} or after \code{em.maxit}
#' iterations.
#'
#' **Standard errors.** Obtained by inverting the negative Hessian of the
#' observed log-likelihood at the MLE on the natural parameter scale,
#' computed numerically with \code{\link[stats]{optimHess}}.
#'
#' @param t        numeric vector of observed times \eqn{> 0}.
#' @param start    optional numeric vector of starting values
#'   \code{c(lambda, alpha, eta)}. Default \code{c(1/mean(t), 0.1, 0.3)}.
#' @param method   optimisation method passed to \code{\link[stats]{optim}}
#'   in the M-step. Default \code{"BFGS"}.
#' @param reltol   relative tolerance passed to \code{\link[stats]{optim}}.
#'   Default \code{1e-12}.
#' @param maxit    maximum number of \code{\link[stats]{optim}} iterations in
#'   each M-step. Default \code{10000}.
#' @param em.tol   convergence threshold on the change in observed
#'   log-likelihood between successive EM iterations. Default \code{1e-8}.
#' @param em.maxit maximum number of EM iterations. Default \code{5000}.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{results}}{data frame with estimates and standard errors of
#'       \eqn{\lambda}, \eqn{\alpha} and \eqn{\eta}.}
#'     \item{\code{lambda}, \code{alpha}, \code{eta}}{numeric; parameter
#'       estimates.}
#'     \item{\code{time.ratio}}{numeric; \eqn{1/\hat\eta}, the ratio of the
#'       delayed mean to the reference mean (\eqn{> 1}).}
#'     \item{\code{v}}{numeric vector; posterior probability that each
#'       observation belongs to the delayed component.}
#'     \item{\code{loglike}}{numeric; maximised log-likelihood.}
#'     \item{\code{AIC}, \code{BIC}}{numeric; information criteria of the cE
#'       model.}
#'     \item{\code{AIC.ref}, \code{BIC.ref}}{numeric; information criteria of
#'       the exponential reference model.}
#'     \item{\code{iter}}{integer; number of EM iterations.}
#'   }
#' @seealso \code{\link{dcexp}}, \code{\link{fit.e}},
#'   \code{\link{ml.cexp.cv}}
#' @export
#' @examples
#' set.seed(1)
#' d1 <- rcexp(2000, lambda = 0.5, alpha = 0.15, eta = 0.3)
#' f1 <- ml.cexp(d1$t)
#' f1$results
#' f1$time.ratio
#' table(true = d1$z, fitted = as.numeric(f1$v > 0.5))
#' @importFrom stats dexp optim optimHess plogis qlogis
ml.cexp <- function(t, start = NULL, method = "BFGS", reltol = 1e-12,
                    maxit = 10000, em.tol = 1e-8, em.maxit = 5000) {
  n <- length(t)

  lambda0    <- 1 / mean(t)
  ref.loglik <- sum(dexp(t, rate = lambda0, log = TRUE))

  if (is.null(start)) {
    start <- c(lambda0, 0.1, 0.3)           # lambda, alpha, eta
  }
  lambda <- start[1]; alpha <- start[2]; eta <- start[3]

  cexp.ml <- function(b.hat, t, v) {
    lambda.hat <- exp(b.hat[1])
    eta.hat    <- plogis(b.hat[2])          # keeps eta in (0,1)
    sum(ll.cexp(t, lambda = lambda.hat, eta = eta.hat, v = v))
  }

  loglik <- NULL
  loglik[1] <- -Inf
  cond <- FALSE
  i <- 1
  while (cond == FALSE) {
    i <- i + 1
    # E-step
    v <- alpha * dexp(t, rate = eta * lambda) / dcexp(t, lambda, alpha, eta)
    alpha <- mean(v)
    # M-step
    fit <- optim(par = c(log(lambda), qlogis(eta)),
                 fn = cexp.ml, t = t, v = v, method = method,
                 control = list(fnscale = -1, maxit = maxit, reltol = reltol))
    lambda <- exp(fit$par[1])
    eta    <- plogis(fit$par[2])
    loglik[i] <- sum(dcexp(t, lambda, alpha, eta, log = TRUE))
    if (abs(loglik[i] - loglik[i - 1]) < em.tol || i >= em.maxit) cond <- TRUE
  }

  obs.ll <- function(par, t) sum(dcexp(t, par[1], par[2], par[3], log = TRUE))
  par  <- c(lambda, alpha, eta)
  Hess <- optimHess(par = par, fn = obs.ll, t = t)
  se   <- sqrt(diag(solve(-Hess)))

  results <- data.frame(estimates = par, se = se)
  rownames(results) <- c("lambda", "alpha", "eta")

  lc <- loglik[i]
  v  <- alpha * dexp(t, rate = eta * lambda) / dcexp(t, lambda, alpha, eta)

  list(results = results,
       lambda = lambda, alpha = alpha, eta = eta,
       time.ratio = 1 / eta,                # delayed mean / reference mean (> 1)
       v = v,                               # P(delayed | t)
       loglike = lc,
       AIC = -2 * lc + 3 * 2,
       BIC = -2 * lc + 3 * log(n),
       AIC.ref = -2 * ref.loglik + 1 * 2,
       BIC.ref = -2 * ref.loglik + 1 * log(n),
       iter = i)
}

# ---- contaminated Weibull ---------------------------------------------

#' Contaminated Weibull (cW) model: EM estimation
#'
#' Fits the contaminated Weibull distribution
#' \eqn{(1-\alpha)\,\mathcal{W}(\sigma,k) + \alpha\,\mathcal{W}(\sigma\eta^{-1/k},k)}
#' by maximum likelihood using the EM algorithm.
#'
#' @details
#' **Initialisation.** If \code{start} is not supplied, \eqn{\sigma} and
#' \eqn{k} are initialised at the MLEs of the Weibull reference model, with
#' \eqn{\alpha^{(0)} = 0.1} and \eqn{\eta^{(0)} = 0.3}.
#'
#' **E-step.** Computes the posterior probability that each observation
#' belongs to the delayed component,
#' \deqn{v_i = \frac{\alpha f_\text{W}(t_i;\sigma\eta^{-1/k},k)}
#'   {f_\text{cW}(t_i;\sigma,k,\alpha,\eta)}.}
#'
#' **M-step.** \eqn{\alpha} is updated in closed form as \eqn{\bar{v}};
#' \eqn{\sigma}, \eqn{k} and \eqn{\eta} are updated by maximising the
#' weighted complete-data log-likelihood with \code{\link[stats]{optim}},
#' using \eqn{\log\sigma}, \eqn{\log k} and \eqn{\mathrm{logit}(\eta)}.
#'
#' **Convergence.** The algorithm stops when the change in the observed
#' log-likelihood is below \code{em.tol} or after \code{em.maxit}
#' iterations.
#'
#' **Standard errors.** Obtained by inverting the negative Hessian of the
#' observed log-likelihood at the MLE on the natural parameter scale,
#' computed numerically with \code{\link[stats]{optimHess}}.
#'
#' @inheritParams ml.cexp
#' @param start optional numeric vector of starting values
#'   \code{c(sigma, k, alpha, eta)}.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{results}}{data frame with estimates and standard errors of
#'       \eqn{\sigma}, \eqn{k}, \eqn{\alpha} and \eqn{\eta}.}
#'     \item{\code{sigma}, \code{k}, \code{alpha}, \code{eta}}{numeric;
#'       parameter estimates.}
#'     \item{\code{time.ratio}}{numeric; \eqn{\hat\eta^{-1/\hat k}}, the
#'       ratio of any quantile of the delayed component to the corresponding
#'       quantile of the reference component (\eqn{> 1}).}
#'     \item{\code{v}}{numeric vector; posterior probability that each
#'       observation belongs to the delayed component.}
#'     \item{\code{loglike}}{numeric; maximised log-likelihood.}
#'     \item{\code{AIC}, \code{BIC}}{numeric; information criteria of the cW
#'       model.}
#'     \item{\code{AIC.ref}, \code{BIC.ref}}{numeric; information criteria of
#'       the Weibull reference model.}
#'     \item{\code{iter}}{integer; number of EM iterations.}
#'   }
#' @seealso \code{\link{dcweib}}, \code{\link{fit.w}},
#'   \code{\link{ml.cweib.cv}}
#' @export
#' @examples
#' set.seed(1)
#' d3 <- rcweib(2000, sigma = 10, k = 1.5, alpha = 0.15, eta = 0.25)
#' f3 <- ml.cweib(d3$t)
#' f3$results
#' f3$time.ratio
#'
#' \donttest{
#' # Table 1 of the paper: January 2024 Ecuadorian births
#' t <- births$lag + 0.5
#' fe <- fit.e(t); fce <- ml.cexp(t)
#' fw <- fit.w(t); fcw <- ml.cweib(t)
#' round(c(E = fe$BIC, cE = fce$BIC, W = fw$BIC, cW = fcw$BIC))
#' fcw$results
#' }
#' @importFrom stats dweibull optim optimHess plogis qlogis
ml.cweib <- function(t, start = NULL, method = "BFGS", reltol = 1e-12,
                     maxit = 10000, em.tol = 1e-8, em.maxit = 5000) {
  n <- length(t)

  ref.ml <- function(b, t) {
    sum(dweibull(t, shape = exp(b[2]), scale = exp(b[1]), log = TRUE))
  }
  ref <- optim(par = c(log(mean(t)), 0), fn = ref.ml, t = t, method = method,
               control = list(fnscale = -1, maxit = maxit, reltol = reltol))
  ref.loglik <- ref$value

  if (is.null(start)) {
    start <- c(exp(ref$par[1]), exp(ref$par[2]), 0.1, 0.3)  # sigma, k, alpha, eta
  }
  sigma <- start[1]; k <- start[2]; alpha <- start[3]; eta <- start[4]

  cweib.ml <- function(b.hat, t, v) {
    sigma.hat <- exp(b.hat[1])
    k.hat     <- exp(b.hat[2])
    eta.hat   <- plogis(b.hat[3])           # keeps eta in (0,1)
    sum(ll.cweib(t, sigma = sigma.hat, k = k.hat, eta = eta.hat, v = v))
  }

  loglik <- NULL
  loglik[1] <- -Inf
  cond <- FALSE
  i <- 1
  while (cond == FALSE) {
    i <- i + 1
    # E-step
    v <- alpha * dweibull(t, shape = k, scale = sigma * eta^(-1 / k)) /
      dcweib(t, sigma, k, alpha, eta)
    alpha <- mean(v)
    # M-step
    fit <- optim(par = c(log(sigma), log(k), qlogis(eta)),
                 fn = cweib.ml, t = t, v = v, method = method,
                 control = list(fnscale = -1, maxit = maxit, reltol = reltol))
    sigma <- exp(fit$par[1])
    k     <- exp(fit$par[2])
    eta   <- plogis(fit$par[3])
    loglik[i] <- sum(dcweib(t, sigma, k, alpha, eta, log = TRUE))
    if (abs(loglik[i] - loglik[i - 1]) < em.tol || i >= em.maxit) cond <- TRUE
  }

  obs.ll <- function(par, t) sum(dcweib(t, par[1], par[2], par[3], par[4], log = TRUE))
  par  <- c(sigma, k, alpha, eta)
  Hess <- optimHess(par = par, fn = obs.ll, t = t)
  se   <- sqrt(diag(solve(-Hess)))

  results <- data.frame(estimates = par, se = se)
  rownames(results) <- c("sigma", "k", "alpha", "eta")

  lc <- loglik[i]
  v  <- alpha * dweibull(t, shape = k, scale = sigma * eta^(-1 / k)) /
    dcweib(t, sigma, k, alpha, eta)

  list(results = results,
       sigma = sigma, k = k, alpha = alpha, eta = eta,
       time.ratio = eta^(-1 / k),           # delayed quantile / reference quantile (> 1)
       v = v,
       loglike = lc,
       AIC = -2 * lc + 4 * 2,
       BIC = -2 * lc + 4 * log(n),
       AIC.ref = -2 * ref.loglik + 2 * 2,
       BIC.ref = -2 * ref.loglik + 2 * log(n),
       iter = i)
}
