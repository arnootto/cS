#' cS: Contaminated Survival Models
#'
#' @description
#' Density, survival, hazard and random generation functions, and maximum
#' likelihood estimation, for the contaminated survival (cS) model.
#'
#' The cS model is a two-component mixture in which both components share a
#' reference survival model \eqn{\mathcal{R}(\boldsymbol\theta)} and differ
#' only in the level of the hazard. Its survival function is
#' \deqn{S_\text{C}(t;\boldsymbol\theta,\alpha,\eta)
#'   = (1-\alpha)\,S_\text{R}(t;\boldsymbol\theta)
#'   + \alpha\,S_\text{R}(t;\boldsymbol\theta)^{\eta},}
#' where \eqn{\alpha \in (0,1)} is the probability of belonging to the
#' \emph{delayed} component and \eqn{\eta \in (0,1)} is the hazard multiplier
#' of that component. Because \eqn{\eta < 1}, the second component is always
#' the low-hazard (delayed) one, so no label swapping is needed.
#' Two members are provided:
#' \describe{
#'   \item{\strong{cE}}{Contaminated exponential: the reference model is
#'     exponential with rate \eqn{\lambda}, and the delayed component is
#'     exponential with rate \eqn{\eta\lambda}.}
#'   \item{\strong{cW}}{Contaminated Weibull: the reference model is Weibull
#'     with scale \eqn{\sigma} and shape \eqn{k}, and the delayed component is
#'     Weibull with scale \eqn{\sigma\eta^{-1/k}} and the same shape.}
#' }
#'
#' The mixing weight may also depend on concomitant variables through a logit
#' link, \eqn{\alpha_i = \exp(\mathbf{v}_i^\top\boldsymbol\gamma) /
#' \{1 + \exp(\mathbf{v}_i^\top\boldsymbol\gamma)\}}.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{dcexp}}, \code{\link{scexp}}, \code{\link{hcexp}},
#'         \code{\link{rcexp}}}{Density, survival, hazard and random
#'         generation for the cE distribution.}
#'   \item{\code{\link{dcweib}}, \code{\link{scweib}}, \code{\link{hcweib}},
#'         \code{\link{rcweib}}}{Density, survival, hazard and random
#'         generation for the cW distribution.}
#'   \item{\code{\link{fit.e}}, \code{\link{fit.w}}}{Maximum likelihood fit
#'         of the exponential and Weibull reference models.}
#'   \item{\code{\link{ml.cexp}}, \code{\link{ml.cweib}}}{EM estimation of
#'         the cE and cW models.}
#'   \item{\code{\link{ml.cexp.cv}}, \code{\link{ml.cweib.cv}}}{EM estimation
#'         of the cE and cW models with concomitant variables.}
#'   \item{\code{\link{lrt}}}{Likelihood ratio test for nested concomitant
#'         models.}
#'   \item{\code{\link{births}}}{January 2024 Ecuadorian birth registration
#'         data.}
#' }
#'
#' @references
#' Otto, A. (2026). Who registers late? Identifying delayed birth
#' registrations from two administrative dates.
#'
#' @docType package
#' @name cS-package
#' @aliases cS
"_PACKAGE"
