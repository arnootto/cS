#' Ecuadorian birth registration data, January 2024
#'
#' Birth and registration dates of the 17,144 live births that occurred in
#' Ecuador in January 2024, with four characteristics recorded for each birth.
#' These are the data analysed in Otto (2026).
#'
#' The data were taken from the 2024 live birth statistics (Estadisticas de
#' Nacidos Vivos, ENV) published by the National Institute of Statistics and
#' Census of Ecuador (INEC). Births with a missing or negative registration
#' lag, or with the place of birth not stated, were removed, as were births
#' with a missing value in any of the four indicators. The script that
#' builds the data from the INEC file is in the \code{data-raw} folder of the
#' package source.
#'
#' The registry records calendar dates, so \code{lag} is a whole number of
#' days and can be zero. In the paper the survival time is taken as
#' \code{lag + 0.5}, the midpoint of the day on which registration took
#' place. Every birth in the file has a registration date, so no observation
#' is censored.
#'
#' @format A data frame with 17,144 rows and 7 variables:
#' \describe{
#'   \item{\code{birth_date}}{date of birth (\code{Date}).}
#'   \item{\code{registration_date}}{date of registration at the civil
#'     registry (\code{Date}).}
#'   \item{\code{lag}}{days from birth to registration.}
#'   \item{\code{home}}{\code{1} if the birth occurred at home, \code{0}
#'     otherwise.}
#'   \item{\code{rural}}{\code{1} if the mother's residence is rural,
#'     \code{0} otherwise.}
#'   \item{\code{indigenous}}{\code{1} if the mother self-identifies as
#'     Indigenous, \code{0} otherwise.}
#'   \item{\code{out_province}}{\code{1} if the birth was registered in a
#'     province other than the province of birth, \code{0} otherwise.}
#' }
#'
#' @source Instituto Nacional de Estadistica y Censos (INEC) (2024).
#'   Estadisticas de Nacidos Vivos 2024 (ENV), open data.
#'   \url{https://www.ecuadorencifras.gob.ec}
#'
#' @references
#' Otto, A. (2026). Who registers late? Identifying delayed birth
#' registrations from two administrative dates.
#'
#' @examples
#' data(births)
#' str(births)
#' births$t <- births$lag + 0.5
#' summary(births$t)
"births"
