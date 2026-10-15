# Analysis scripts

Scripts that produce the numbers, tables and figures of the manuscripts that
accompany `epicurves`. They are not part of the installed package (the folder
is listed in `.Rbuildignore`).

Each folder holds the scripts and a `results/` folder with the `.rds` and
`.log` files they wrote. The numbers in the manuscripts were taken from
these files, so they can be checked without rerunning the computations.

| Folder | Content |
|---|---|
| `sird-ga/` | SIRD model under the Gaussian approximation: simulation study, profile likelihood of the death rate, Paraíba application |
| `piecewise-ensemble/` | Growth-curve ensembles and segmentation rules: Monte Carlo studies, Paraíba application |
| `rjournal/` | Software article: validation of the general engine, worked examples, bootstrap, timing |

## Running

Install the package and the optional `edetect` package, then run each script
from its own folder. Every script reads and writes paths relative to that
folder.

```r
remotes::install_github("Raydonal/epicurves")
remotes::install_github("castlaboratory/edetect")
```

```bash
cd analysis/sird-ga
Rscript sim_study.R 1 20000 100 4    # scenario 1, N = 20000, 100 replicates, 4 cores
```

The `.sh` files run whole batches and skip results that already exist.
The simulation studies take hours on four cores; launch them detached, for
example `setsid nohup ./resume_all.sh >> results/run_all.log 2>&1 &`.

## Scripts by result

### `sird-ga/`

| Result | Scripts |
|---|---|
| Simulation study tables (bias, RMSE, coverage) | `sim_study.R` (one scenario and population size per call; `run_all.sh`, `resume_all.sh` run the grid), `summarise_sim.R`, `make_sim_tables.R` |
| Wald standard errors against Monte Carlo spread | `sim_se_diagnostic.R` |
| Profile likelihood of the death rate (figure) | `profile_alpha_sim.R`, `profile_alpha_fig.R` |
| Paraíba, constant transmission | `paraiba_fit.R` |
| Paraíba, transmission changing every four weeks (tables and figure) | `paraiba_tv.R`, `paraiba_report.R` |
| Paraíba, profile in the case reporting probability | `paraiba_tv_profile.R` |
| Paraíba, sensitivity to the death reporting probability | `paraiba_tv_q.R` |
| Paraíba, parametric bootstrap | `paraiba_bootstrap.R`, `paraiba_bootstrap_summary.R` |
| Accuracy of the numerical integration | `check_Qk.R` |

### `piecewise-ensemble/`

| Result | Scripts |
|---|---|
| Monte Carlo comparison of segmentation rules (table) | `scenarios.R`, `mc_study.R` (`run_after_article1.sh` runs all scenarios), `summarise_mc.R`, `make_mc_tables.R` |
| Detection rates expected by chance | `mc_chance_baseline.R` |
| Overdispersed series without change points (table) | `mc_negbin.R` (`run_negbin.sh`), `summarise_negbin.R` |
| E-detector run length to false alarm | `mc_arl.R` (`run_arl.sh`), `summarise_arl.R` |
| Paraíba segmentation (table and figure) | `paraiba_rules.R`, `paraiba_rules_fig.R` |

### `rjournal/`

| Result | Scripts |
|---|---|
| Validation of the general engine on SEIRD and SIRDV models | `validate_generic.R` (`run_validation.sh`) |
| SEIRD worked example | `example_seird.R`, `example_seird_sum.R`, `check_Qk.R` |
| Paraíba growth-curve example | `example_paraiba.R` |
| Parametric bootstrap of the SEIRD example | `example_bootstrap.R` (`run_bootstrap.sh`) |
| Cost of one likelihood evaluation | `timing.R` |

`run_examples.sh` runs the examples and the timing after the validation.

The figure scripts write to `figures/` in their folder.
