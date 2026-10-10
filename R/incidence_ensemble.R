#' Daily incidence implied by a growth curve
#'
#' For a cumulative growth curve \eqn{C(u)} with \eqn{C(0)=y_0}, returns the
#' expected count in each unit interval, \eqn{\mu_u = C(u) - C(u-1)},
#' \eqn{u = 1, 2, \ldots}.
#'
#' @param u Integer local times (1 = first day of the segment).
#' @param model One of `"logistic"`, `"gompertz"`, `"richards"`.
#' @param y0,K,gamma Curve parameters (`0 < y0 < K`, `gamma > 0`).
#' @param beta Richards shape parameter.
#' @export
growth_incidence <- function(u, model, y0, K, gamma, beta = 1.5) {
  cf <- switch(model,
    logistic = function(t) logistic_curve(t, y0, K, gamma),
    gompertz = function(t) gompertz_curve(t, y0, K, gamma),
    richards = function(t) richards_curve(t, y0, K, gamma, beta))
  pmax(cf(u) - cf(u - 1), 1e-8)
}

count_loglik <- function(x, mu, family, size) {
  if (family == "poisson") sum(stats::dpois(x, mu, log = TRUE))
  else sum(stats::dnbinom(x, mu = mu, size = size, log = TRUE))
}

#' Fit one growth curve to daily counts by maximum likelihood
#'
#' Poisson or negative-binomial likelihood for the daily counts `x` of a
#' segment, with mean given by [growth_incidence()]. Parameters are
#' optimised on the scale \eqn{(\log y_0, \log(K - y_0), \log\gamma[, \log\phi])}
#' from a small grid of starting values.
#'
#' @param x Daily counts of the segment (local times `1, ..., length(x)`).
#' @param model Growth model.
#' @param family `"poisson"` or `"negbin"`.
#' @param beta Richards shape (held fixed).
#' @param start Optional warm start on the unconstrained scale (e.g. `par`
#'   of a previous fit), tried in addition to the default grid.
#' @return A list with `coef`, `loglik`, `aic`, `bic`, `k` (number of
#'   estimated parameters), `fitted`, `model`, `family`, `beta`, `n`.
#' @export
fit_growth_incidence <- function(x, model = c("logistic", "gompertz", "richards"),
                                 family = c("poisson", "negbin"), beta = 1.5, start = NULL) {
  model <- match.arg(model); family <- match.arg(family)
  n <- length(x); u <- seq_len(n)
  tot <- max(sum(x), 1)
  unpack <- function(v) {
    y0 <- exp(v[1]); K <- y0 + exp(v[2]); g <- exp(v[3])
    list(y0 = y0, K = K, gamma = g, size = if (family == "negbin") exp(v[4]) else Inf)
  }
  nll <- function(v) {
    p <- unpack(v)
    if (!all(is.finite(unlist(p[1:3]))) || p$gamma > 5) return(1e12)
    mu <- growth_incidence(u, model, p$y0, p$K, p$gamma, beta)
    if (!all(is.finite(mu))) return(1e12)
    ll <- count_loglik(x, mu, family, p$size)
    if (!is.finite(ll)) 1e12 else -ll
  }
  grid <- expand.grid(rem = c(1.5, 5), g = c(0.03, 0.12), y0f = c(0.05, 1))
  starts <- lapply(seq_len(nrow(grid)), function(r) {
    v0 <- c(log(grid$y0f[r] * tot + 1), log(grid$rem[r] * tot), log(grid$g[r]))
    if (family == "negbin") v0 <- c(v0, log(10))
    v0
  })
  if (!is.null(start) && length(start) == 3L + (family == "negbin")) starts <- c(list(start), starts)
  best <- NULL
  for (v0 in starts) {
    f <- tryCatch(stats::optim(v0, nll, method = "Nelder-Mead",
                               control = list(maxit = 2000, reltol = 1e-10)),
                  error = function(e) NULL)
    if (!is.null(f) && (is.null(best) || f$value < best$value)) best <- f
  }
  best <- tryCatch({
    f2 <- stats::optim(best$par, nll, method = "BFGS", control = list(maxit = 300))
    if (f2$value < best$value) f2 else best
  }, error = function(e) best)
  p <- unpack(best$par)
  k <- 3L + (family == "negbin")
  vc <- tryCatch({
    H <- stats::optimHess(best$par, nll)
    V <- solve(H)
    if (all(is.finite(V)) && all(diag(V) > 0)) V else NULL
  }, error = function(e) NULL)
  ll <- -best$value
  structure(list(coef = c(y0 = p$y0, K = p$K, gamma = p$gamma, size = p$size),
                 loglik = ll, aic = -2 * ll + 2 * k, bic = -2 * ll + k * log(n), k = k,
                 fitted = growth_incidence(u, model, p$y0, p$K, p$gamma, beta),
                 model = model, family = family, beta = beta, n = n, par = best$par,
                 vcov = vc),
            class = "growth_incidence_fit")
}

#' @export
predict.growth_incidence_fit <- function(object, newdata, ...) {
  cf <- object$coef
  growth_incidence(newdata, object$model, cf[["y0"]], cf[["K"]], cf[["gamma"]], object$beta)
}

#' Forecast band of a growth-incidence fit
#'
#' Delta-method band \eqn{\hat\mu_u \pm z\,\widehat{\mathrm{se}}(\hat\mu_u)}, with the
#' standard error propagated from the inverse observed information of the
#' unconstrained parameters by a numerical gradient.
#'
#' @param object A `growth_incidence_fit`.
#' @param u Local times.
#' @param z Band half-width multiplier.
#' @return Matrix with columns `mu`, `lower`, `upper`.
#' @export
predict_band <- function(object, u, z = 2) {
  mu_of <- function(v) growth_incidence(u, object$model, exp(v[1]), exp(v[1]) + exp(v[2]),
                                        exp(v[3]), object$beta)
  v <- object$par
  mu <- mu_of(v)
  se <- rep(Inf, length(u))
  if (!is.null(object$vcov)) {
    V <- object$vcov[1:3, 1:3, drop = FALSE]
    G <- vapply(1:3, function(j) {
      h <- 1e-5 * max(1, abs(v[j])); e <- rep(0, length(v)); e[j] <- h
      (mu_of(v + e) - mu_of(v - e)) / (2 * h)
    }, numeric(length(u)))
    G <- matrix(G, nrow = length(u))
    se <- sqrt(pmax(rowSums((G %*% V) * G), 0))
  }
  se[!is.finite(se)] <- Inf
  hw <- if (z == 0) 0 else z * se
  cbind(mu = mu, lower = pmax(mu - hw, 1e-6), upper = mu + hw)
}

#' Ensemble of growth curves fitted to daily counts
#'
#' Fits the Logistic, Gompertz and Richards incidence models by maximum
#' likelihood and combines them with Akaike weights
#' \eqn{w_i \propto \exp(-\Delta\mathrm{AIC}_i/2)} (`weights = "akaike"`) or
#' with inverse-MSE weights on the count scale (`weights = "inverse_mse"`).
#'
#' @param x Daily counts.
#' @param family `"poisson"` or `"negbin"`.
#' @param weights Weighting scheme.
#' @param beta Richards shape (held fixed).
#' @param models Models included.
#' @param warm Optional previous `incidence_ensemble` whose estimates are
#'   used as warm starts.
#' @return An object of class `incidence_ensemble`.
#' @export
fit_ensemble_incidence <- function(x, family = c("poisson", "negbin"),
                                   weights = c("akaike", "inverse_mse"), beta = 1.5,
                                   models = c("gompertz", "richards", "logistic"),
                                   warm = NULL) {
  family <- match.arg(family); weights <- match.arg(weights)
  fits <- lapply(models, function(m) {
    st <- if (!is.null(warm) && !is.null(warm$fits[[m]])) warm$fits[[m]]$par
    fit_growth_incidence(x, m, family, beta, start = st)
  })
  names(fits) <- models
  aic <- vapply(fits, `[[`, numeric(1), "aic")
  mse <- vapply(fits, function(f) mean((x - f$fitted)^2), numeric(1))
  w <- if (weights == "akaike") {
    z <- exp(-(aic - min(aic)) / 2); z / sum(z)
  } else {
    z <- 1 / mse; z / sum(z)
  }
  fitted <- as.vector(vapply(fits, `[[`, numeric(length(x)), "fitted") %*% w)
  size <- if (family == "negbin") sum(w * vapply(fits, function(f) f$coef[["size"]], numeric(1))) else Inf
  structure(list(fits = fits, weights = w, aic = aic, mse = mse, fitted = fitted,
                 family = family, size = size, scheme = weights, n = length(x),
                 loglik = count_loglik(x, fitted, family, size)),
            class = "incidence_ensemble")
}

#' @export
predict.incidence_ensemble <- function(object, newdata, ...) {
  as.vector(vapply(object$fits, function(f) predict(f, newdata),
                   numeric(length(newdata))) %*% object$weights)
}

#' Online segmentation of an epidemic curve with e-detectors
#'
#' Sequentially monitors daily counts against a forecast band of a
#' growth-curve ensemble fitted to the current segment, and opens a new
#' segment when an anytime-valid e-detector (package \pkg{edetect}) raises an
#' alarm. Within a segment the ensemble is refitted every `refit_every` days
#' on data observed before the day being monitored, so each band is
#' predictable with respect to the past.
#'
#' The band \eqn{[\ell_t, u_t]} is the envelope, over the models carrying at
#' least `w_min` ensemble weight, of their delta-method bands
#' ([predict_band()]). Two one-sided detectors run in parallel, each with
#' target average run length `2 * arl`: an upward detector tests the
#' composite null "conditional mean at most \eqn{u_t}" and a downward
#' detector the null "at least \eqn{\ell_t}". If the conditional mean lies in
#' the band, the probability of any false alarm by time \eqn{t} is at most
#' \eqn{t/\mathrm{arl}} (union bound over the two detectors). With
#' `chart = "poisson"` the detectors are Poisson u charts with exposure
#' \eqn{u_t} (resp. \eqn{\ell_t}); with `chart = "subgaussian"` they monitor
#' the negative-binomial Pearson residuals \eqn{(x_t-u_t)/s(u_t)} and
#' \eqn{(x_t-\ell_t)/s(\ell_t)}, \eqn{s(\mu)=\sqrt{\mu(1+\mu/\hat\phi)}},
#' declared sub-Gaussian with unit scale (an approximation).
#'
#' @param x Daily counts.
#' @param d_min Minimum segment length before monitoring starts.
#' @param arl Target average run length to a false alarm.
#' @param refit_every Refit interval (days) within a segment.
#' @param chart `"poisson"` or `"subgaussian"`.
#' @param z Half-width multiplier of each model's band (`0` monitors the
#'   point forecast).
#' @param w_min Minimum ensemble weight for a model to enter the envelope.
#' @param weights Ensemble weighting scheme.
#' @param type e-detector type, `"sr"` or `"cusum"`.
#' @param beta Richards shape.
#' @return A list with `cuts` (last day of each closed segment), `alarms`
#'   (alarm times), `delay`, `band` (matrix of monitored forecast, lower and
#'   upper band, `NA` on unmonitored days), `segments` and the final ensemble
#'   of every segment (`fits`).
#' @export
edetect_segment <- function(x, d_min = 20L, arl = 1000, refit_every = 7L,
                            chart = c("poisson", "subgaussian"), z = 2, w_min = 0.05,
                            weights = c("akaike", "inverse_mse"),
                            type = c("sr", "cusum"), beta = 1.5) {
  if (!requireNamespace("edetect", quietly = TRUE)) stop("package 'edetect' is required")
  chart <- match.arg(chart); weights <- match.arg(weights); type <- match.arg(type)
  family <- if (chart == "poisson") "poisson" else "negbin"
  mk <- function(dir) {
    if (chart == "poisson") {
      edetect::edetect_design(arl = 2 * arl, class = "poisson", center = 1, direction = dir, type = type)
    } else {
      edetect::edetect_design(arl = 2 * arl, class = "subgaussian", center = 0, sigma = 1,
                              direction = dir, type = type)
    }
  }
  d_up <- mk("up"); d_dn <- mk("down")
  Tn <- length(x)
  band <- matrix(NA_real_, Tn, 3, dimnames = list(NULL, c("mu", "lower", "upper")))
  cuts <- integer(); alarms <- integer(); fits <- list()
  start <- 1L
  ens_band <- function(ens, u) {
    keep <- names(ens$weights)[ens$weights >= w_min]
    b <- lapply(ens$fits[keep], predict_band, u = u, z = z)
    c(mu = predict(ens, u), lower = min(vapply(b, function(m) m[1, "lower"], 1)),
      upper = max(vapply(b, function(m) m[1, "upper"], 1)))
  }
  while (start + d_min <= Tn) {
    up <- edetect::edetect_init(d_up); dn <- edetect::edetect_init(d_dn)
    t <- start + d_min
    ens <- NULL; last_fit <- -Inf; alarm <- FALSE
    while (t <= Tn) {
      if (is.null(ens) || t - last_fit >= refit_every) {
        ens <- tryCatch(fit_ensemble_incidence(x[start:(t - 1L)], family, weights, beta,
                                               warm = ens),
                        error = function(e) ens)
        last_fit <- t
      }
      bt <- ens_band(ens, t - start + 1L)
      if (!is.finite(bt[["mu"]])) { t <- t + 1L; next }
      if (!is.finite(bt[["upper"]])) bt[["upper"]] <- 100 * max(bt[["mu"]], 1)
      if (!is.finite(bt[["lower"]])) bt[["lower"]] <- 1e-6
      bt[["lower"]] <- max(bt[["lower"]], 1e-6)
      band[t, ] <- bt
      if (chart == "poisson") {
        up <- suppressWarnings(edetect::edetect_update(up, x[t], n = bt[["upper"]]))
        dn <- suppressWarnings(edetect::edetect_update(dn, x[t], n = bt[["lower"]]))
      } else {
        s <- function(m) sqrt(m * (1 + m / ens$size))
        up <- suppressWarnings(edetect::edetect_update(up, (x[t] - bt[["upper"]]) / s(bt[["upper"]])))
        dn <- suppressWarnings(edetect::edetect_update(dn, (x[t] - bt[["lower"]]) / s(bt[["lower"]])))
      }
      hit <- c(up = !is.na(up$alarm_time), dn = !is.na(dn$alarm_time))
      if (any(hit)) {
        alarm <- TRUE
        ce <- if (hit[["up"]]) up$change_estimate else dn$change_estimate
        new_start <- max((start + d_min - 1L) + ce, start + d_min)
        fits[[length(fits) + 1L]] <- ens
        cuts <- c(cuts, new_start - 1L); alarms <- c(alarms, t)
        start <- new_start
        break
      }
      t <- t + 1L
    }
    if (!alarm) { fits[[length(fits) + 1L]] <- ens; break }
  }
  list(cuts = cuts, alarms = alarms, delay = alarms - cuts, band = band,
       segments = data.frame(start = c(1L, cuts + 1L), end = c(cuts, Tn)), fits = fits,
       settings = list(d_min = d_min, arl = arl, refit_every = refit_every, chart = chart,
                       z = z, w_min = w_min, weights = weights, type = type))
}

#' Offline segmentation of daily counts by penalised likelihood
#'
#' Binary segmentation of the daily-count series: a segment is split at the
#' candidate day maximising the summed ensemble log-likelihood of the two
#' halves, and the split is kept only if it lowers the Bayesian information
#' criterion, with the split charged for the parameters of one additional
#' curve and the change location. The parameter count of an ensemble is the
#' weight-averaged count of its members (three for every curve used here):
#' the members are alternative descriptions of a single mean, not additive
#' components. Each half is fitted by [fit_ensemble_incidence()].
#'
#' @param x Daily counts.
#' @param d_min Minimum segment length.
#' @param family,weights,beta Passed to [fit_ensemble_incidence()].
#' @param step Spacing of the coarse candidate grid (days); the optimum is
#'   then refined at unit spacing within `step` days.
#' @param max_cuts Maximum number of change points.
#' @return A list with `cuts` and `segments`.
#' @export
bic_segment <- function(x, d_min = 20L, family = c("poisson", "negbin"),
                        weights = c("akaike", "inverse_mse"), beta = 1.5,
                        step = 2L, max_cuts = 10L) {
  family <- match.arg(family); weights <- match.arg(weights)
  ens_k <- 3L + (family == "negbin")   # Akaike-weighted parameter count of one curve
  ll_of <- function(a, b) {
    e <- tryCatch(fit_ensemble_incidence(x[a:b], family, weights, beta), error = function(e) NULL)
    if (is.null(e)) -Inf else e$loglik
  }
  queue <- list(c(1L, length(x))); cuts <- integer()
  while (length(queue) && length(cuts) < max_cuts) {
    seg <- queue[[1]]; queue <- queue[-1]
    a <- seg[1]; b <- seg[2]; n <- b - a + 1L
    if (n < 2L * d_min) next
    cand <- seq(a + d_min - 1L, b - d_min, by = step)
    ll0 <- ll_of(a, b)
    llc <- vapply(cand, function(d) ll_of(a, d) + ll_of(d + 1L, b), numeric(1))
    if (step > 1L) {                 # refine around the coarse optimum
      d0 <- cand[which.max(llc)]
      fine <- setdiff(seq(max(a + d_min - 1L, d0 - step + 1L), min(b - d_min, d0 + step - 1L)), cand)
      if (length(fine)) {
        cand <- c(cand, fine)
        llc <- c(llc, vapply(fine, function(d) ll_of(a, d) + ll_of(d + 1L, b), numeric(1)))
      }
    }
    j <- which.max(llc)
    bic0 <- -2 * ll0 + ens_k * log(n)
    bic1 <- -2 * llc[j] + (2L * ens_k + 1L) * log(n)
    if (is.finite(bic1) && bic1 < bic0) {
      cuts <- sort(c(cuts, cand[j]))
      queue <- c(queue, list(c(a, cand[j]), c(cand[j] + 1L, b)))
    }
  }
  list(cuts = cuts, segments = data.frame(start = c(1L, cuts + 1L), end = c(cuts, length(x))))
}
