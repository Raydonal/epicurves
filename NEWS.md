# epicurves 0.0.0.9000

Initial implementation, built to validate and reproduce the results of two
companion manuscripts derived from the thesis.

## Findings made while validating against the thesis

- **Richards growth-curve sign correction.** The thesis's Richards-curve
  equation (Chapter 3), transcribed literally, uses `exp(+beta*gamma*t)`;
  this drives the curve to 0 as `t -> Inf` instead of converging to the
  carrying capacity `K`, inconsistent with the standard Richards (1959)
  form and with the thesis's own (correct) Gompertz/Logistic equations.
  `richards_curve()` implements the corrected form,
  `exp(-beta*gamma*t)`. This needs to be corrected in the Article 2
  manuscript's methods section.
- **Partition-point over-segmentation.** The flat "50% MSE reduction"
  split-confirmation rule described in the thesis can over-segment long,
  structurally stable phases, both in a re-simulation of the thesis's own
  Scenario 1 and in the real Paraíba death series. `fit_piecewise_ensemble()`
  exposes a `confirm_criterion = "aic"` option (a BIC-style,
  parameter-count-penalised alternative) that is markedly more
  conservative while still recovering the Scenario 1 change points to the
  same accuracy as the thesis's own reported values. Recommended for
  Article 2 as a robustness improvement, and as a direct response to the
  2023 referees' novelty concerns about the (non-piecewise) precursor
  paper.
- **SIRD-GA real-data identifiability: two basins.** Fitting SIRD-GA to raw
  daily case counts is sensitive to day-of-week reporting artefacts; a
  7-day centered rolling average (`smooth_series()`) is recommended before
  estimation. On the smoothed Paraíba series, a 5-start sweep
  (`sird_ga_fit_multistart()`) reveals two distinct optima: a
  "fast-dynamics" basin (4/5 starts, loglik ~2947, R0~1.03, mean
  infectious period ~1.4 days, epidemiologically implausible) and a
  "slow-dynamics" basin (1/5 starts, loglik ~2924, R0~1.03, mean
  infectious period ~4.5 days, matching the published estimates and known
  COVID-19 natural history). Within the fast basin, gamma and alpha trade
  off freely while gamma+alpha stays fixed, mirroring the thesis's own
  appendix finding that gamma/alpha are estimated far less precisely than
  lambda or R0. Recommended for Article 1: report the multi-start summary
  and prefer the slow-dynamics basin on epidemiological grounds, rather
  than reporting an unconstrained MLE point estimate. See the SIRD-GA
  vignette for the full table and discussion.
