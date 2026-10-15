# COVID-19 surveillance data, Paraíba, Brazil (2020)

Daily COVID-19 case and death counts for the Brazilian state of Paraíba,
covering the first epidemic wave used throughout both companion
manuscripts: the SIRD-GA application (daily cases) and the piecewise
ensemble application (cumulative/daily deaths).

## Usage

``` r
paraiba_covid19
```

## Format

A data frame with 289 rows (18 March-31 December 2020) and 6 columns:

- date:

  Calendar date.

- epi_week:

  Epidemiological week.

- new_cases:

  Newly confirmed cases on this date.

- total_cases:

  Cumulative confirmed cases.

- new_deaths:

  Newly registered COVID-19 deaths on this date (first non-zero value on
  31 March 2020).

- deaths:

  Cumulative COVID-19 deaths.

## Source

State-level epidemiological surveillance bulletins, as compiled for the
authors' 2023 ensemble precursor study.
