# Define an event of a mass-action compartmental model

An event moves one individual from compartment `from` to compartment
`to` at rate `rate * x[from] * x[by]` per unit of time and per capita of
the population (fractions `x`), or `rate * x[from]` when `by` is `NULL`.
This covers infection (`by` = the infectious compartment), progression,
recovery, death, vaccination and waning of immunity.

## Usage

``` r
ga_event(from, to, rate, by = NULL)
```

## Arguments

- from, to:

  Compartment names. `to` may be `NULL` for removal from the population.

- rate:

  Name of the rate parameter.

- by:

  Optional name of a second compartment multiplying the rate (mass
  action).

## Value

An object of class `ga_event`.

## Examples

``` r
ga_event("S", "I", rate = "lambda", by = "I")
#> $from
#> [1] "S"
#> 
#> $to
#> [1] "I"
#> 
#> $rate
#> [1] "lambda"
#> 
#> $by
#> [1] "I"
#> 
#> attr(,"class")
#> [1] "ga_event"
```
