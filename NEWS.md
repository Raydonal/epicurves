# epicurves (development)

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
- **SIRD-GA real-data identifiability.** Fitting SIRD-GA to raw daily case
  counts is sensitive to day-of-week reporting artefacts; a 7-day
  centered rolling average is recommended before estimation (see the
  SIRD-GA vignette for a worked comparison against the raw series).
