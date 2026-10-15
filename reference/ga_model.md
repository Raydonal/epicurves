# Define a mass-action compartmental model

Define a mass-action compartmental model

## Usage

``` r
ga_model(compartments, events)
```

## Arguments

- compartments:

  Character vector of compartment names. All compartments are explicit;
  the fractions sum to one at the start and the events conserve that sum
  unless `to` is `NULL`.

- events:

  Named list of
  [`ga_event()`](https://raydonal.github.io/epicurves/reference/ga_event.md)
  objects.

## Value

An object of class `ga_model`, used by
[`ga_loglik()`](https://raydonal.github.io/epicurves/reference/ga_loglik.md),
[`ga_simulate()`](https://raydonal.github.io/epicurves/reference/ga_simulate.md)
and
[`ga_fit()`](https://raydonal.github.io/epicurves/reference/ga_fit.md).

## Examples

``` r
seird <- ga_model(c("S", "E", "I", "R", "D"), list(
  infection   = ga_event("S", "E", "lambda", by = "I"),
  progression = ga_event("E", "I", "sigma"),
  recovery    = ga_event("I", "R", "gamma"),
  death       = ga_event("I", "D", "alpha")))
seird
#> Mass-action compartmental model: S, E, I, R, D 
#>   infection    S -> E  at rate lambda * S * I
#>   progression  E -> I  at rate sigma * E
#>   recovery     I -> R  at rate gamma * I
#>   death        I -> D  at rate alpha * I
#> Rate parameters: lambda, sigma, gamma, alpha 
```
