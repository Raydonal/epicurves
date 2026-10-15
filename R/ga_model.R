#' Define an event of a mass-action compartmental model
#'
#' An event moves one individual from compartment `from` to compartment `to`
#' at rate `rate * x[from] * x[by]` per unit of time and per capita of the
#' population (fractions `x`), or `rate * x[from]` when `by` is `NULL`. This
#' covers infection (`by` = the infectious compartment), progression,
#' recovery, death, vaccination and waning of immunity.
#'
#' @param from,to Compartment names. `to` may be `NULL` for removal from the
#'   population.
#' @param rate Name of the rate parameter.
#' @param by Optional name of a second compartment multiplying the rate
#'   (mass action).
#' @return An object of class `ga_event`.
#' @examples
#' ga_event("S", "I", rate = "lambda", by = "I")
#' @export
ga_event <- function(from, to, rate, by = NULL) {
  structure(list(from = from, to = to, rate = rate, by = by), class = "ga_event")
}

#' Define a mass-action compartmental model
#'
#' @param compartments Character vector of compartment names. All
#'   compartments are explicit; the fractions sum to one at the start and
#'   the events conserve that sum unless `to` is `NULL`.
#' @param events Named list of [ga_event()] objects.
#' @return An object of class `ga_model`, used by [ga_loglik()],
#'   [ga_simulate()] and [ga_fit()].
#' @examples
#' seird <- ga_model(c("S", "E", "I", "R", "D"), list(
#'   infection   = ga_event("S", "E", "lambda", by = "I"),
#'   progression = ga_event("E", "I", "sigma"),
#'   recovery    = ga_event("I", "R", "gamma"),
#'   death       = ga_event("I", "D", "alpha")))
#' seird
#' @export
ga_model <- function(compartments, events) {
  stopifnot(is.character(compartments), length(compartments) >= 2L, !anyDuplicated(compartments),
            is.list(events), length(events) >= 1L, !is.null(names(events)), all(nzchar(names(events))),
            all(vapply(events, inherits, logical(1), "ga_event")))
  for (ev in events) {
    nm <- c(ev$from, ev$to, ev$by)
    if (!all(nm %in% compartments)) stop("unknown compartment in event: ", paste(setdiff(nm, compartments), collapse = ", "))
  }
  d <- length(compartments)
  V <- matrix(0, d, length(events), dimnames = list(compartments, names(events)))
  for (e in seq_along(events)) {
    V[events[[e]]$from, e] <- V[events[[e]]$from, e] - 1
    if (!is.null(events[[e]]$to)) V[events[[e]]$to, e] <- V[events[[e]]$to, e] + 1
  }
  params <- unique(vapply(events, `[[`, character(1), "rate"))
  structure(list(compartments = compartments, events = events, V = V, params = params),
            class = "ga_model")
}

#' @export
print.ga_model <- function(x, ...) {
  cat("Mass-action compartmental model:", paste(x$compartments, collapse = ", "), "\n")
  for (nm in names(x$events)) {
    ev <- x$events[[nm]]
    cat(sprintf("  %-12s %s -> %s  at rate %s * %s%s\n", nm, ev$from, if (is.null(ev$to)) "(out)" else ev$to,
                ev$rate, ev$from, if (is.null(ev$by)) "" else paste0(" * ", ev$by)))
  }
  cat("Rate parameters:", paste(x$params, collapse = ", "), "\n")
  invisible(x)
}

#' Observation streams for the Gaussian-approximation likelihood
#'
#' `ga_incidence()` describes counts of an event over each observation
#' interval (new infections, new deaths, new hospital admissions);
#' `ga_prevalence()` describes counts of a compartment at the end of each
#' interval. Reported counts are modelled as a binomial thinning of the
#' latent count with probability `prob` plus Gaussian measurement error with
#' variance `noise` times the latent count, approximated by a Gaussian with
#' the same mean and variance.
#'
#' @param event Name of an event of the model.
#' @param compartment Name of a compartment of the model.
#' @param prob Name of the reporting-probability parameter, or a number in
#'   `(0, 1]`.
#' @param noise Name of the measurement-error parameter, or a non-negative
#'   number.
#' @return An object of class `ga_obs`.
#' @export
ga_incidence <- function(event, prob = 1, noise = 0) {
  structure(list(type = "incidence", target = event, prob = prob, noise = noise), class = "ga_obs")
}

#' @rdname ga_incidence
#' @export
ga_prevalence <- function(compartment, prob = 1, noise = 0) {
  structure(list(type = "prevalence", target = compartment, prob = prob, noise = noise), class = "ga_obs")
}

# Internal: augmented stoichiometry (compartments + one counter per observed event),
# reactant indices and observation descriptors, all 0-based for the C++ engine.
ga_spec <- function(model, obs, rate_names) {
  ev_obs <- unique(vapply(Filter(function(o) o$type == "incidence", obs), `[[`, character(1), "target"))
  if (!all(ev_obs %in% names(model$events))) stop("unknown event in observations: ", paste(setdiff(ev_obs, names(model$events)), collapse = ", "))
  comp_obs <- vapply(Filter(function(o) o$type == "prevalence", obs), `[[`, character(1), "target")
  if (!all(comp_obs %in% model$compartments)) stop("unknown compartment in observations")
  d <- length(model$compartments)
  Vc <- matrix(0, length(ev_obs), ncol(model$V))
  for (j in seq_along(ev_obs)) Vc[j, match(ev_obs[j], names(model$events))] <- 1
  V <- rbind(model$V, Vc)
  idx <- function(nm) if (is.null(nm)) -1L else match(nm, model$compartments) - 1L
  list(V = V, d = d, counters = ev_obs,
       pidx = match(vapply(model$events, `[[`, character(1), "rate"), rate_names) - 1L,
       r1 = vapply(model$events, function(e) idx(e$from), integer(1)),
       r2 = vapply(model$events, function(e) idx(e$by), integer(1)),
       obs_type = vapply(obs, function(o) if (o$type == "incidence") 1L else 2L, integer(1)),
       obs_index = vapply(obs, function(o) if (o$type == "incidence") d + match(o$target, ev_obs) - 1L
                          else match(o$target, model$compartments) - 1L, integer(1)))
}

ga_obs_value <- function(v, theta) if (is.character(v)) theta[[v]] else v

#' Gaussian-approximation log-likelihood of a compartmental model
#'
#' Kalman-filter log-likelihood of counts observed at times
#' `Delta, 2 Delta, ..., n Delta` under the Gaussian (linear noise)
#' approximation of the density-dependent Markov jump process defined by
#' `model`. Over each interval the ordinary differential equation of the
#' mean, its resolvent and the Lyapunov equation of the state-noise
#' covariance are integrated jointly by the classical Runge-Kutta method in
#' compiled code. Incidence observations are handled with event counters
#' and an augmented state that contains the state at the start and the end
#' of each interval.
#'
#' @param model A [ga_model()].
#' @param theta Named list or vector with the rate parameters of `model` and
#'   any parameters named in `obs`.
#' @param init Named vector of initial fractions, one per compartment
#'   (normalised to sum to one).
#' @param data Matrix or data frame of counts, one column per element of
#'   `obs` (in the same order), one row per interval. `NA` marks a missing
#'   count.
#' @param obs List of [ga_incidence()] / [ga_prevalence()] streams.
#' @param N Population size.
#' @param Delta Length of the observation interval.
#' @param nsub Runge-Kutta steps per interval.
#' @param varying Optional named list of time-varying rate parameters: each
#'   element is a vector with one value per interval, overriding the
#'   constant value in `theta`.
#' @param full If `TRUE`, return the list with one-step-ahead predicted
#'   means (`fitted`) and standard deviations (`sd`) of the counts and the
#'   path of the mean (`x`, fractions) instead of the scalar.
#' @return The log-likelihood, or a list if `full = TRUE`.
#' @export
ga_loglik <- function(model, theta, init, data, obs, N, Delta = 1, nsub = 10L,
                      varying = NULL, full = FALSE) {
  theta <- as.list(theta)
  Y <- as.matrix(data); storage.mode(Y) <- "double"
  if (ncol(Y) != length(obs)) stop("data must have one column per observation stream")
  n <- nrow(Y)
  rn <- model$params
  if (!all(rn %in% c(names(theta), names(varying)))) stop("missing rate parameters: ", paste(setdiff(rn, c(names(theta), names(varying))), collapse = ", "))
  th <- matrix(vapply(rn, function(p) if (p %in% names(theta)) as.numeric(theta[[p]]) else NA_real_, 1),
               nrow = 1, dimnames = list(NULL, rn))
  if (length(varying)) {
    th <- th[rep(1, n), , drop = FALSE]
    for (p in names(varying)) { stopifnot(length(varying[[p]]) == n); th[, p] <- varying[[p]] }
  }
  p_obs <- vapply(obs, function(o) ga_obs_value(o$prob, theta), 1)
  t_obs <- vapply(obs, function(o) ga_obs_value(o$noise, theta), 1)
  x0 <- init[model$compartments]
  if (anyNA(x0) || any(x0 < 0) || anyNA(th) || any(th < 0) || any(p_obs <= 0 | p_obs > 1) || any(t_obs < 0))
    return(if (full) list(loglik = -Inf) else -Inf)
  x0 <- x0 / sum(x0)
  sp <- ga_spec(model, obs, rn)
  x0a <- c(x0, numeric(length(sp$counters)))
  out <- ga_loglik_cpp(Y, N, Delta, x0a, th, sp$V, sp$pidx, sp$r1, sp$r2, sp$obs_type, sp$obs_index,
                       p_obs, t_obs, nsub = as.integer(nsub))
  if (!full) return(out$loglik)
  if (!is.null(out$x)) colnames(out$x) <- c(model$compartments, if (length(sp$counters)) paste0("n_", sp$counters))
  out
}

#' Exact simulation of a compartmental model and its surveillance counts
#'
#' Gillespie simulation of the Markov jump process defined by `model`, with
#' counts aggregated over intervals of length `Delta` and reported counts
#' drawn as binomial thinning plus Gaussian measurement error (rounded and
#' truncated at zero), the data-generating process assumed by [ga_loglik()].
#'
#' @inheritParams ga_loglik
#' @param init_counts Named integer vector of initial counts per compartment.
#' @param n Number of intervals.
#' @param varying Optional named list of rate parameters that change between
#'   intervals, each a vector of length `n` (as in [ga_loglik()]); a rate
#'   listed here may be absent from `theta`.
#' @return A list with `data` (matrix of reported counts, one column per
#'   stream), `latent` (matrix of latent counts per stream) and `states`
#'   (compartment counts at the end of each interval, first row the initial
#'   state).
#' @export
ga_simulate <- function(model, theta, init_counts, obs, n, Delta = 1, varying = NULL) {
  theta <- as.list(theta)
  x <- as.numeric(init_counts[model$compartments]); names(x) <- model$compartments
  N <- sum(x)
  E <- length(model$events)
  from <- match(vapply(model$events, `[[`, "", "from"), model$compartments)
  by <- vapply(model$events, function(e) if (is.null(e$by)) NA_integer_ else match(e$by, model$compartments), 1L)
  rate <- vapply(model$events, `[[`, "", "rate")
  for (p in names(varying)) stopifnot(length(varying[[p]]) == n)
  k <- vapply(rate, function(p) as.numeric(if (p %in% names(varying)) varying[[p]][1] else theta[[p]]), 1)
  V <- model$V
  states <- matrix(NA_real_, n + 1, length(x), dimnames = list(NULL, model$compartments))
  states[1, ] <- x
  counts <- matrix(0, n, E, dimnames = list(NULL, names(model$events)))
  t <- 0
  for (j in seq_len(n)) {
    tend <- j * Delta
    for (p in names(varying)) k[rate == p] <- varying[[p]][j]
    repeat {
      r <- k * x[from] * ifelse(is.na(by), 1, x[by] / N)
      tot <- sum(r)
      if (tot <= 0) { t <- tend; break }
      t <- t + stats::rexp(1, tot)
      if (t > tend) { t <- tend; break }
      e <- sample.int(E, 1L, prob = r)
      x <- x + V[, e]
      counts[j, e] <- counts[j, e] + 1
    }
    states[j + 1, ] <- x
  }
  latent <- sapply(obs, function(o) if (o$type == "incidence") counts[, o$target] else states[-1, o$target])
  latent <- matrix(latent, nrow = n)
  rep_counts <- latent
  for (s in seq_along(obs)) {
    p <- ga_obs_value(obs[[s]]$prob, theta); t2 <- ga_obs_value(obs[[s]]$noise, theta)
    o <- stats::rbinom(n, size = latent[, s], prob = p) + stats::rnorm(n, 0, sqrt(t2 * latent[, s]))
    rep_counts[, s] <- pmax(round(o), 0)
  }
  nm <- vapply(obs, function(o) paste0(o$type, "_", o$target), "")
  colnames(latent) <- colnames(rep_counts) <- nm
  list(data = rep_counts, latent = latent, states = states)
}
