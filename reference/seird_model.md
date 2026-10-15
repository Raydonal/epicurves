# SEIRD compartmental model (susceptible-exposed-infectious-recovered-dead)

Extends SIRD with an exposed (non-infectious, incubating) compartment,
the natural extension for diseases with an appreciable incubation period
such as COVID-19. State: `(s, e, i, r)`, with the death fraction
implicit as `1 - s - e - i - r`. Events: infection (`S+I -> E+I`, rate
`lambda*s*i`), progression (`E -> I`, rate `sigma*e`), recovery
(`I -> R`, rate `gamma*i`) and disease-induced death (`I -> `(implicit
`D`), rate `alpha*i`).

## Usage

``` r
seird_model()
```

## Value

A `jump_model` with compartments `c("s","e","i","r")` and parameters
`lambda`, `sigma`, `gamma`, `alpha`.
