# Prepares the `paraiba_covid19` package dataset from the raw surveillance
# export used in the 2023 ensemble precursor manuscript and in both
# companion papers' real-data applications.
#
# Source: state-level COVID-19 surveillance bulletin data for Paraíba,
# Brazil, 18 March-31 December 2020 (n = 289 daily observations), as
# compiled in data-raw/paraiba_covid19_raw.xlsx.

library(readxl)

raw <- as.data.frame(read_excel("data-raw/paraiba_covid19_raw.xlsx"))

paraiba_covid19 <- data.frame(
  date = as.Date(raw$date),
  epi_week = raw$epi_week,
  new_cases = raw$newCases,
  total_cases = raw$totalCases,
  new_deaths = raw$newDeaths,
  deaths = raw$deaths
)

# Population of Paraíba used throughout both manuscripts (2022 census
# baseline), exported alongside the data for convenience.
paraiba_population <- 4.06e6

usethis::use_data(paraiba_covid19, overwrite = TRUE)
usethis::use_data(paraiba_population, overwrite = TRUE)
