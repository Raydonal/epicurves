# epicurves: Likelihood Inference for Epidemic Curves with Stochastic Compartmental Models and Growth-Curve Ensembles

Likelihood-based tools for epidemic curves. Stochastic compartmental
models are fitted to surveillance counts through the Gaussian (linear
noise) approximation of the underlying density-dependent Markov jump
process and a Kalman filter: any mass-action model can be declared by
its events, observed through the incidence of events (new cases, new
deaths) or the prevalence of compartments, with under-reporting and
time-varying rates. Logistic, Gompertz and Richards growth curves are
fitted to daily counts by Poisson or negative binomial likelihood,
combined with Akaike weights, and the series is divided into phases by
penalised likelihood or by anytime-valid sequential e-detectors. The
package accompanies two articles derived from Oliveira's doctoral thesis
and reproduces their analyses.

## See also

Useful links:

- <https://github.com/Raydonal/epicurves>

- Report bugs at <https://github.com/Raydonal/epicurves/issues>

## Author

**Maintainer**: Raydonal Ospina <raydonal@castlab.org>

Authors:

- Raydonal Ospina <raydonal@castlab.org>

- Jaciele Oliveira <jaciele.jesus@ufpe.br>

- Cristiano Ferraz <cferraz@castlab.org>

- André Leite
