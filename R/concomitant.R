# =======================================================================
# Contaminated survival models with a concomitant variable
# alpha_i = P(i delayed | v_i) = exp(v_i'gamma)/(1 + exp(v_i'gamma))
# Labelling: eta in (0,1), so the second component is always the DELAYED one.
#   No swap is needed: the constraint identifies the labelling for every i,
#   whatever value alpha_i takes.
# =======================================================================

# ---- concomitant contaminated exponential -----------------------------

#' Contaminated exponential (cE) model with concomitant variables
#'
#' Fits the contaminated exponential model in which the probability that
#' observation \eqn{i} belongs to the delayed component depends on
#' concomitant variables \eqn{\mathbf{v}_i} through a logit link,
#' \deqn{\alpha_i = \frac{\exp(\mathbf{v}_i^\top\boldsymbol\gamma)}
#'   {1 + \exp(\mathbf{v}_i^\top\boldsymbol\gamma)},}
#' so that \eqn{t_i} follows \eqn{(1-\alpha_i)\,\mathcal{E}(\lambda) +
#' \alpha_i\,\mathcal{E}(\eta\lambda)}. The parameters \eqn{\lambda} and
#' \eqn{\eta} are common to all observations. Estimation is by maximum
#' likelihood using the EM algorithm.
#'
#' @details
#' **Initialisation.** If \code{init} is not supplied,
#' \eqn{\lambda^{(0)} = 1/\bar{t}}, \eqn{\eta^{(0)} = 0.3}, the intercept of
#' \eqn{\boldsymbol\gamma} is \eqn{\mathrm{logit}(0.1)} and the remaining
#' coefficients are zero.
#'
#' **E-step.** Computes the posterior probability that each observation
#' belongs to the delayed component.
#'
#' **M-step.** The complete-data log-likelihood separates into a gating part
#' for \eqn{\boldsymbol\gamma} (a weighted logistic regression) and a survival
#' part for \eqn{(\lambda, \eta)}; both are maximised with
#' \code{\link[stats]{optim}}.
#'
#' **Standard errors.** Obtained by inverting the negative Hessian of the
#' observed log-likelihood on the unconstrained scale
#' \eqn{(\log\lambda, \mathrm{logit}(\eta), \boldsymbol\gamma)}, computed
#' numerically with \code{\link[stats]{optimHess}}.
#'
#' @param formula a \code{\link[stats]{formula}} whose left-hand side is the
#'   observed time, e.g. \code{t ~ 1}. The right-hand side is not used.
#' @param alpha.formula a one-sided \code{\link[stats]{formula}} for the
#'   mixing weight, e.g. \code{~ x1 + x2}. Default \code{~1} (constant
#'   mixing weight).
#' @param data    a data frame containing the variables in \code{formula} and
#'   \code{alpha.formula}.
#' @param init    optional list of starting values with components
#'   \code{lambda}, \code{eta} and \code{gamma}.
#' @param method  optimisation method passed to \code{\link[stats]{optim}}
#'   in the M-step. Default \code{"BFGS"}.
#' @param reltol  relative tolerance passed to \code{\link[stats]{optim}}.
#'   Default \code{1e-12}.
#' @param maxit   maximum number of \code{\link[stats]{optim}} iterations in
#'   each M-step. Default \code{10000}.
#' @param em.tol  convergence threshold on the change in observed
#'   log-likelihood between successive EM iterations. Default \code{1e-8}.
#' @param em.maxit maximum number of EM iterations. Default \code{1000}.
#' @param hessian logical; if \code{TRUE} (default) standard errors,
#'   \eqn{z}-values and \eqn{p}-values are computed.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{results}}{data frame with estimates (and, if
#'       \code{hessian = TRUE}, standard errors, \eqn{z}-values and
#'       \eqn{p}-values) of \eqn{\log\lambda}, \eqn{\mathrm{logit}(\eta)} and
#'       \eqn{\boldsymbol\gamma}.}
#'     \item{\code{lambda}, \code{eta}}{numeric; parameter estimates.}
#'     \item{\code{gamma}}{numeric vector; estimated coefficients of the
#'       mixing weight.}
#'     \item{\code{alpha}}{numeric vector; fitted mixing weight
#'       \eqn{\hat\alpha_i} of each observation.}
#'     \item{\code{alpha.bar}}{numeric; mean of \code{alpha}.}
#'     \item{\code{time.ratio}}{numeric; \eqn{1/\hat\eta}.}
#'     \item{\code{v}}{numeric vector; posterior probability that each
#'       observation belongs to the delayed component.}
#'     \item{\code{V}}{the model matrix of \code{alpha.formula}.}
#'     \item{\code{t}}{the observed times.}
#'     \item{\code{loglike}}{numeric; maximised log-likelihood.}
#'     \item{\code{npar}}{integer; number of estimated parameters.}
#'     \item{\code{converged}}{logical; \code{TRUE} if the EM algorithm
#'       converged.}
#'     \item{\code{iter}}{integer; number of EM iterations.}
#'     \item{\code{AIC}, \code{BIC}}{numeric; information criteria.}
#'   }
#' @seealso \code{\link{ml.cexp}}, \code{\link{ml.cweib.cv}},
#'   \code{\link{lrt}}, \code{\link{rcexp.cv}}
#' @export
#' @examples
#' set.seed(1)
#' n  <- 1000
#' cv <- data.frame(rural = rbinom(n, 1, 0.4), indigenous = rbinom(n, 1, 0.1))
#' fo <- ~ rural + indigenous
#' g  <- c(log(0.15/0.85), log(2.5), log(3))
#' d1 <- rcexp.cv(n, lambda = 0.1, eta = 0.3, gamma = g,
#'                covars = cv, alpha.formula = fo)
#' f0 <- ml.cexp.cv(t ~ 1, ~1, data = d1)
#' f1 <- ml.cexp.cv(t ~ 1, fo, data = d1)
#' round(f1$results, 3)
#' lrt(f0, f1)
#' @importFrom stats model.frame model.response model.matrix optim optimHess
#'   dexp plogis qlogis pnorm
ml.cexp.cv <- function(formula, alpha.formula = ~1, data, init = NULL,
                       method = "BFGS", reltol = 1e-12, maxit = 10000,
                       em.tol = 1e-8, em.maxit = 1000, hessian = TRUE) {

  mf <- model.frame(formula, data)
  t  <- model.response(mf, "numeric")
  V  <- model.matrix(alpha.formula, data = data)
  n  <- length(t)

  # Q1: gating part
  lc1 <- function(gamma.hat, V, w) {
    vg.hat    <- V %*% gamma.hat
    alpha.hat <- exp(vg.hat)/(1 + exp(vg.hat))
    sum((1 - w)*log(1 - alpha.hat) + w*log(alpha.hat))
  }

  # Q2: survival part, eta in (0,1) during estimation
  lc2 <- function(par.hat, t, w) {
    lambda.hat <- exp(par.hat[1])
    eta.hat    <- plogis(par.hat[2])
    sum((1 - w)*dexp(t, rate = lambda.hat, log = TRUE) +
          w*dexp(t, rate = eta.hat*lambda.hat, log = TRUE))
  }

  # observed-data log-likelihood, unconstrained scale (for the Hessian)
  cexp.cv.ml <- function(b.hat, t, V) {
    vg.hat    <- V %*% b.hat[-(1:2)]
    alpha.hat <- exp(vg.hat)/(1 + exp(vg.hat))
    sum(dcexp(t, exp(b.hat[1]), alpha.hat, plogis(b.hat[2]), log = TRUE))
  }

  if (is.null(init)) {
    lambda <- 1/mean(t)
    eta    <- 0.3
    gamma  <- c(log(0.1/0.9), rep(0, ncol(V) - 1))
  } else {
    lambda <- init$lambda; eta <- init$eta; gamma <- init$gamma
  }

  loglik <- -Inf; converged <- FALSE; iter <- 0
  while (!converged && iter < em.maxit) {
    iter <- iter + 1
    old.loglik <- loglik

    # E-step
    alpha <- as.vector(exp(V %*% gamma)/(1 + exp(V %*% gamma)))
    w <- alpha*dexp(t, rate = eta*lambda)/dcexp(t, lambda, alpha, eta)

    # M-step
    Q1 <- optim(par = gamma, fn = lc1, V = V, w = w, method = method,
                control = list(fnscale = -1, maxit = maxit, reltol = reltol))
    Q2 <- optim(par = c(log(lambda), qlogis(eta)), fn = lc2, t = t, w = w,
                method = method,
                control = list(fnscale = -1, maxit = maxit, reltol = reltol))

    gamma  <- Q1$par
    lambda <- exp(Q2$par[1])
    eta    <- plogis(Q2$par[2])

    alpha  <- as.vector(exp(V %*% gamma)/(1 + exp(V %*% gamma)))
    loglik <- sum(dcexp(t, lambda, alpha, eta, log = TRUE))
    if (abs(loglik - old.loglik) < em.tol) converged <- TRUE
  }

  w <- alpha*dexp(t, rate = eta*lambda)/dcexp(t, lambda, alpha, eta)

  pars <- c(log(lambda), qlogis(eta), gamma)
  if (hessian) {
    hess <- optimHess(pars, fn = cexp.cv.ml, t = t, V = V)
    se   <- sqrt(diag(base::solve(-hess)))
    tval <- pars/se
    results <- data.frame(Estimate = pars, Std.Error = se, t.value = tval,
                          P.t = round(2*pnorm(abs(tval), lower.tail = FALSE), 5))
  } else {
    results <- data.frame(Estimate = pars)
  }
  rownames(results) <- c("log_lambda", "logit_eta", paste0("gamma_", colnames(V)))

  npar <- length(pars)
  list(results = results, lambda = lambda, eta = eta, gamma = gamma,
       alpha = alpha, alpha.bar = mean(alpha), time.ratio = 1/eta,
       v = w, V = V, t = t,
       loglike = loglik, npar = npar,
       converged = converged, iter = iter,
       AIC = -2*loglik + 2*npar, BIC = -2*loglik + npar*log(n))
}

# ---- concomitant contaminated Weibull ---------------------------------

#' Contaminated Weibull (cW) model with concomitant variables
#'
#' Fits the contaminated Weibull model in which the probability that
#' observation \eqn{i} belongs to the delayed component depends on
#' concomitant variables \eqn{\mathbf{v}_i} through a logit link,
#' \deqn{\alpha_i = \frac{\exp(\mathbf{v}_i^\top\boldsymbol\gamma)}
#'   {1 + \exp(\mathbf{v}_i^\top\boldsymbol\gamma)},}
#' so that \eqn{t_i} follows \eqn{(1-\alpha_i)\,\mathcal{W}(\sigma,k) +
#' \alpha_i\,\mathcal{W}(\sigma\eta^{-1/k},k)}. The parameters \eqn{\sigma},
#' \eqn{k} and \eqn{\eta} are common to all observations. Estimation is by
#' maximum likelihood using the EM algorithm.
#'
#' @details
#' **Initialisation.** If \code{init} is not supplied,
#' \eqn{\sigma^{(0)} = \bar{t}}, \eqn{k^{(0)} = 1}, \eqn{\eta^{(0)} = 0.3},
#' the intercept of \eqn{\boldsymbol\gamma} is \eqn{\mathrm{logit}(0.1)} and
#' the remaining coefficients are zero.
#'
#' **E-step.** Computes the posterior probability that each observation
#' belongs to the delayed component.
#'
#' **M-step.** The complete-data log-likelihood separates into a gating part
#' for \eqn{\boldsymbol\gamma} (a weighted logistic regression) and a survival
#' part for \eqn{(\sigma, k, \eta)}; both are maximised with
#' \code{\link[stats]{optim}}.
#'
#' **Standard errors.** Obtained by inverting the negative Hessian of the
#' observed log-likelihood on the unconstrained scale
#' \eqn{(\log\sigma, \log k, \mathrm{logit}(\eta), \boldsymbol\gamma)},
#' computed numerically with \code{\link[stats]{optimHess}}.
#'
#' @inheritParams ml.cexp.cv
#' @param init optional list of starting values with components
#'   \code{sigma}, \code{k}, \code{eta} and \code{gamma}.
#'
#' @return A list with the same components as \code{\link{ml.cexp.cv}},
#'   with \code{sigma} and \code{k} in place of \code{lambda}, and
#'   \code{time.ratio} equal to \eqn{\hat\eta^{-1/\hat k}}, the ratio of any
#'   quantile of the delayed component to the corresponding quantile of the
#'   reference component. The rows of \code{results} are \code{log_sigma},
#'   \code{log_k}, \code{logit_eta} and the \code{gamma_} coefficients.
#' @seealso \code{\link{ml.cweib}}, \code{\link{ml.cexp.cv}},
#'   \code{\link{lrt}}, \code{\link{rcweib.cv}}
#' @export
#' @examples
#' set.seed(1)
#' n  <- 1000
#' cv <- data.frame(rural = rbinom(n, 1, 0.4), indigenous = rbinom(n, 1, 0.1))
#' g  <- c(log(0.15/0.85), log(2.5), log(3))
#' d2 <- rcweib.cv(n, sigma = 20, k = 1.4, eta = 0.3, gamma = g,
#'                 covars = cv, alpha.formula = ~ rural + indigenous)
#' w1 <- ml.cweib.cv(t ~ 1, ~ rural + indigenous, data = d2)
#' round(w1$results, 3)
#'
#' \donttest{
#' # Table 2 of the paper: January 2024 Ecuadorian births
#' births$t <- births$lag + 0.5
#' f0.w <- ml.cweib.cv(t ~ 1, ~1, data = births)
#' f1.w <- ml.cweib.cv(t ~ 1, ~ home + rural + indigenous + out_province,
#'                     data = births)
#' lrt(f0.w, f1.w)
#' round(f1.w$results, 4)
#' round(exp(f1.w$gamma[-1]), 3)   # odds ratios for membership
#' }
#' @importFrom stats model.frame model.response model.matrix optim optimHess
#'   dweibull plogis qlogis pnorm
ml.cweib.cv <- function(formula, alpha.formula = ~1, data, init = NULL,
                        method = "BFGS", reltol = 1e-12, maxit = 10000,
                        em.tol = 1e-8, em.maxit = 1000, hessian = TRUE) {

  mf <- model.frame(formula, data)
  t  <- model.response(mf, "numeric")
  V  <- model.matrix(alpha.formula, data = data)
  n  <- length(t)

  lc1 <- function(gamma.hat, V, w) {
    vg.hat    <- V %*% gamma.hat
    alpha.hat <- exp(vg.hat)/(1 + exp(vg.hat))
    sum((1 - w)*log(1 - alpha.hat) + w*log(alpha.hat))
  }

  lc2 <- function(par.hat, t, w) {
    sigma.hat <- exp(par.hat[1])
    k.hat     <- exp(par.hat[2])
    eta.hat   <- plogis(par.hat[3])
    sum((1 - w)*dweibull(t, shape = k.hat, scale = sigma.hat, log = TRUE) +
          w*dweibull(t, shape = k.hat,
                     scale = sigma.hat*eta.hat^(-1/k.hat), log = TRUE))
  }

  cweib.cv.ml <- function(b.hat, t, V) {
    vg.hat    <- V %*% b.hat[-(1:3)]
    alpha.hat <- exp(vg.hat)/(1 + exp(vg.hat))
    sum(dcweib(t, exp(b.hat[1]), exp(b.hat[2]), alpha.hat,
               plogis(b.hat[3]), log = TRUE))
  }

  if (is.null(init)) {
    sigma <- mean(t); k <- 1; eta <- 0.3
    gamma <- c(log(0.1/0.9), rep(0, ncol(V) - 1))
  } else {
    sigma <- init$sigma; k <- init$k; eta <- init$eta; gamma <- init$gamma
  }

  loglik <- -Inf; converged <- FALSE; iter <- 0
  while (!converged && iter < em.maxit) {
    iter <- iter + 1
    old.loglik <- loglik

    # E-step
    alpha <- as.vector(exp(V %*% gamma)/(1 + exp(V %*% gamma)))
    w <- alpha*dweibull(t, shape = k, scale = sigma*eta^(-1/k))/
      dcweib(t, sigma, k, alpha, eta)

    # M-step
    Q1 <- optim(par = gamma, fn = lc1, V = V, w = w, method = method,
                control = list(fnscale = -1, maxit = maxit, reltol = reltol))
    Q2 <- optim(par = c(log(sigma), log(k), qlogis(eta)), fn = lc2,
                t = t, w = w, method = method,
                control = list(fnscale = -1, maxit = maxit, reltol = reltol))

    gamma <- Q1$par
    sigma <- exp(Q2$par[1])
    k     <- exp(Q2$par[2])
    eta   <- plogis(Q2$par[3])

    alpha  <- as.vector(exp(V %*% gamma)/(1 + exp(V %*% gamma)))
    loglik <- sum(dcweib(t, sigma, k, alpha, eta, log = TRUE))
    if (abs(loglik - old.loglik) < em.tol) converged <- TRUE
  }

  w <- alpha*dweibull(t, shape = k, scale = sigma*eta^(-1/k))/
    dcweib(t, sigma, k, alpha, eta)

  pars <- c(log(sigma), log(k), qlogis(eta), gamma)
  if (hessian) {
    hess <- optimHess(pars, fn = cweib.cv.ml, t = t, V = V)
    se   <- sqrt(diag(base::solve(-hess)))
    tval <- pars/se
    results <- data.frame(Estimate = pars, Std.Error = se, t.value = tval,
                          P.t = round(2*pnorm(abs(tval), lower.tail = FALSE), 5))
  } else {
    results <- data.frame(Estimate = pars)
  }
  rownames(results) <- c("log_sigma", "log_k", "logit_eta",
                         paste0("gamma_", colnames(V)))

  npar <- length(pars)
  list(results = results, sigma = sigma, k = k, eta = eta, gamma = gamma,
       alpha = alpha, alpha.bar = mean(alpha), time.ratio = eta^(-1/k),
       v = w, V = V, t = t, loglike = loglik, npar = npar,
       converged = converged, iter = iter,
       AIC = -2*loglik + 2*npar, BIC = -2*loglik + npar*log(n))
}

# ---- likelihood ratio test --------------------------------------------

#' Likelihood ratio test for nested concomitant models
#'
#' Compares two nested fits from \code{\link{ml.cexp.cv}} or
#' \code{\link{ml.cweib.cv}}, typically a model with a constant mixing weight
#' (\code{alpha.formula = ~1}) against one with concomitant variables. The
#' statistic \eqn{D = 2(\ell_1 - \ell_0)} is compared with a
#' \eqn{\chi^2} distribution on \eqn{p_1 - p_0} degrees of freedom.
#'
#' Both fits must contain the components \code{loglike} and \code{npar}, so
#' the test is not intended for comparing a contaminated model with its
#' reference model (the null value \eqn{\alpha = 0} lies on the boundary of
#' the parameter space; use AIC or BIC instead).
#'
#' @param fit0 fit of the restricted (null) model.
#' @param fit1 fit of the larger model.
#'
#' @return A data frame with the two log-likelihoods (\code{loglike0},
#'   \code{loglike1}), the test statistic (\code{stat}), the degrees of
#'   freedom (\code{df}) and the \eqn{p}-value (\code{p}).
#' @seealso \code{\link{ml.cexp.cv}}, \code{\link{ml.cweib.cv}}
#' @export
#' @examples
#' set.seed(1)
#' n  <- 1000
#' cv <- data.frame(rural = rbinom(n, 1, 0.4))
#' d1 <- rcexp.cv(n, lambda = 0.1, eta = 0.3, gamma = c(-1.5, 1),
#'                covars = cv, alpha.formula = ~ rural)
#' f0 <- ml.cexp.cv(t ~ 1, ~1, data = d1)
#' f1 <- ml.cexp.cv(t ~ 1, ~ rural, data = d1)
#' lrt(f0, f1)
#' @importFrom stats pchisq
lrt <- function(fit0, fit1) {
  D  <- 2*(fit1$loglike - fit0$loglike)
  df <- fit1$npar - fit0$npar
  data.frame(loglike0 = fit0$loglike, loglike1 = fit1$loglike,
             stat = D, df = df, p = pchisq(D, df, lower.tail = FALSE))
}
