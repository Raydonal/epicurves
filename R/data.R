#' COVID-19 surveillance data, Paraíba, Brazil (2020)
#'
#' Daily COVID-19 case and death counts for the Brazilian state of Paraíba,
#' covering the first epidemic wave used throughout both companion
#' manuscripts: the SIRD-GA application (daily cases) and the piecewise
#' ensemble application (cumulative/daily deaths).
#'
#' @format A data frame with 289 rows (18 March-31 December 2020) and 6
#'   columns:
#' \describe{
#'   \item{date}{Calendar date.}
#'   \item{epi_week}{Epidemiological week.}
#'   \item{new_cases}{Newly confirmed cases on this date.}
#'   \item{total_cases}{Cumulative confirmed cases.}
#'   \item{new_deaths}{Newly registered COVID-19 deaths on this date (first
#'     non-zero value on 31 March 2020).}
#'   \item{deaths}{Cumulative COVID-19 deaths.}
#' }
#' @source State-level epidemiological surveillance bulletins, as compiled
#'   for the authors' 2023 ensemble precursor study.
"paraiba_covid19"

#' Population of Paraíba, Brazil
#'
#' Fixed population size `N` used throughout the SIRD-GA manuscript's
#' real-data application.
#' @format A numeric scalar, `4.06e6` (2022 census baseline).
"paraiba_population"
