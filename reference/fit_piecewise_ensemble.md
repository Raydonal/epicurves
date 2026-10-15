# Sequential piecewise ensemble model

Reproduces the sequential partition-point identification algorithm of
the manuscript: starting from an initial data window, candidate
partition points are searched by grid search
([`find_partition_point()`](https://raydonal.github.io/epicurves/reference/find_partition_point.md));
the window is then advanced by `update_every` observations and, at each
update, the most recent (still-unsplit) segment is re-examined for a new
structural change, which is confirmed only if the pooled MSE of the best
split is at least `threshold` lower than the MSE of leaving that segment
unsplit.

## Usage

``` r
fit_piecewise_ensemble(
  y,
  t = NULL,
  window = 90L,
  update_every = 20L,
  min_obs = 12L,
  criterion = "mse",
  beta_fixed = 1.5,
  threshold = 0.5,
  confirm_criterion = c("mse_threshold", "aic", "permutation"),
  alpha = 0.05,
  B = 99L
)
```

## Arguments

- y:

  Numeric vector of cumulative counts (the full observed series).

- t:

  Optional numeric time index (defaults to `seq_along(y)`).

- window:

  Size of the initial data window used to search for the first partition
  point (default 90, as in the manuscript).

- update_every:

  Number of observations by which the window advances at each sequential
  update (default 20).

- min_obs:

  Minimum admissible segment length (default 12; use 15 to match
  Simulation Scenario 2 of the manuscript).

- criterion, beta_fixed:

  Passed to
  [`fit_ensemble()`](https://raydonal.github.io/epicurves/reference/fit_ensemble.md).

- threshold:

  Relative MSE reduction required to confirm a new partition point under
  `confirm_criterion = "mse_threshold"` (default 0.5, i.e. the split
  must at least halve the MSE relative to the unsplit fit — this is the
  rule described in the manuscript).

- confirm_criterion:

  How a candidate split is confirmed. `"mse_threshold"` (default)
  reproduces the manuscript's rule exactly. `"aic"` instead compares a
  changepoint-penalised AIC for the split against the unsplit fit —
  splitting doubles the number of growth-curve parameters (three models
  re-fit on each side) and adds one degree of freedom for the
  partition-point location itself, and `"aic"` is the more conservative
  choice when a flat MSE-reduction threshold over-segments noisy real
  series (see the package vignette for a worked comparison on both
  simulated and real surveillance data). `"permutation"` instead
  performs a formal hypothesis test
  ([`permutation_test_partition()`](https://raydonal.github.io/epicurves/reference/permutation_test_partition.md)):
  the candidate split is confirmed only if its MSE reduction is unlikely
  under a null, unsplit growth process, at level `alpha`, properly
  accounting for the multiple comparisons implicit in scanning the whole
  candidate grid. This is the criterion we recommend when statistical
  significance, not just a descriptive rule, is required (e.g. for
  direct comparison against nonparametric permutation-based detectors
  such as
  [`ecp::e.divisive`](https://rdrr.io/pkg/ecp/man/e.divisive.html)); it
  is markedly more expensive computationally than the other two.

- alpha:

  Significance level used only when `confirm_criterion = "permutation"`.

- B:

  Number of permutation-test null replicates, used only when
  `confirm_criterion = "permutation"`; passed to
  [`permutation_test_partition()`](https://raydonal.github.io/epicurves/reference/permutation_test_partition.md).

## Value

An object of class `piecewise_fit`: a list with `cuts` (sorted confirmed
partition points), `segments` (list of `ensemble_fit` objects, one per
final segment, each carrying the original time range in `$t_range`),
`fitted` (full fitted curve aligned with `y`), and `history` (list of
grid-search diagnostics at every update step, suitable for reproducing
the manuscript's MSE-profile figures).
