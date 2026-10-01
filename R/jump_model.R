#' Define a density-dependent Markov jump compartmental model
#'
#' A general representation of a compartmental epidemic model as a set of
#' elementary events, each with an integer (here, real-valued after
#' normalisation) jump vector and a rate function of the current state and
#' parameters. [sird_drift()], [sird_jacobian()] and [sird_diffusion()]
#' are the SIRD-specific, hand-derived analogues of [jm_drift()],
#' [jm_jacobian()] and [jm_diffusion()] below; a `jump_model` makes the
#' same Gaussian-diffusion-approximation machinery available for
#' **any** compartmental structure (SEIR, SEIRD, SIRDV, age- or
#' risk-stratified models, ...), at the cost of a numerical rather than
#' analytic Jacobian. [general_ga_fit()] and friends use a `jump_model`
#' exactly as [sird_ga_fit()] uses the hand-derived SIRD functions.
#'
#' @param compartments Character vector naming the explicit state
#'   compartments (the model's state space has one *implicit* compartment,
#'   recovered as `1 - sum(explicit compartments)`, exactly as `d` is
#'   implicit in [sird_drift()]).
#' @param events A list of elementary events, each a list with `jump` (a
#'   numeric vector the same length as `compartments`, giving the change
#'   to each explicit compartment when the event fires) and `rate` (a
#'   function of `(state, parms)` returning the event's instantaneous
#'   rate).
#' @return An object of class `jump_model`.
#' @export
new_jump_model <- function(compartments, events) {
  stopifnot(all(vapply(events, function(e) length(e$jump) == length(compartments), logical(1))))
  structure(list(compartments = compartments, events = events), class = "jump_model")
}

#' Drift vector of a jump model
#' @param model A `jump_model`.
#' @param state,parms As in [sird_drift()].
#' @export
jm_drift <- function(model, state, parms) {
  b <- numeric(length(model$compartments))
  for (ev in model$events) b <- b + ev$rate(state, parms) * ev$jump
  b
}

#' Diffusion covariance matrix of a jump model
#' @inheritParams jm_drift
#' @export
jm_diffusion <- function(model, state, parms) {
  d <- length(model$compartments)
  Sigma <- matrix(0, d, d)
  for (ev in model$events) Sigma <- Sigma + ev$rate(state, parms) * outer(ev$jump, ev$jump)
  Sigma
}

#' Numerical Jacobian of a jump model's drift
#'
#' Central-difference Jacobian of [jm_drift()]. Used in place of an
#' analytic Jacobian (as hand-derived for SIRD in [sird_jacobian()]) so
#' that [new_jump_model()] generalises to arbitrary compartmental
#' structures without per-model algebra. The cost is a small, well-
#' controlled numerical error (step `h`), validated in the package tests
#' against the analytic SIRD Jacobian.
#' @inheritParams jm_drift
#' @param h Finite-difference step.
#' @export
jm_jacobian <- function(model, state, parms, h = 1e-6) {
  d <- length(state)
  J <- matrix(0, d, d)
  for (j in seq_len(d)) {
    up <- state; up[j] <- up[j] + h
    dn <- state; dn[j] <- dn[j] - h
    J[, j] <- (jm_drift(model, up, parms) - jm_drift(model, dn, parms)) / (2 * h)
  }
  J
}

#' SEIRD compartmental model (susceptible-exposed-infectious-recovered-dead)
#'
#' Extends SIRD with an exposed (non-infectious, incubating) compartment,
#' the natural extension for diseases with an appreciable incubation
#' period such as COVID-19. State: `(s, e, i, r)`, with the death fraction
#' implicit as `1 - s - e - i - r`. Events: infection (`S+I -> E+I`, rate
#' `lambda*s*i`), progression (`E -> I`, rate `sigma*e`), recovery
#' (`I -> R`, rate `gamma*i`) and disease-induced death (`I -> `(implicit
#' `D`), rate `alpha*i`).
#'
#' @return A `jump_model` with compartments `c("s","e","i","r")` and
#'   parameters `lambda`, `sigma`, `gamma`, `alpha`.
#' @export
seird_model <- function() {
  new_jump_model(
    compartments = c("s", "e", "i", "r"),
    events = list(
      infection    = list(jump = c(-1,  1, 0, 0), rate = function(x, p) p[["lambda"]] * x[1] * x[3]),
      progression  = list(jump = c( 0, -1, 1, 0), rate = function(x, p) p[["sigma"]]  * x[2]),
      recovery     = list(jump = c( 0,  0,-1, 1), rate = function(x, p) p[["gamma"]]  * x[3]),
      death        = list(jump = c( 0,  0,-1, 0), rate = function(x, p) p[["alpha"]]  * x[3])
    )
  )
}

#' SIRD compartmental model, as a `jump_model`
#'
#' The same SIRD model as [sird_drift()]/[sird_jacobian()]/
#' [sird_diffusion()], expressed as a [new_jump_model()], purely so that
#' [general_ga_transition()]/[general_ga_loglik()] can be validated
#' against the hand-derived, analytic-Jacobian SIRD-GA engine (which
#' remains the implementation used for the main SIRD-GA results). See
#' `test-general-engine.R` for the numerical cross-check.
#'
#' @return A `jump_model` with compartments `c("s","i","r")` and
#'   parameters `lambda`, `gamma`, `alpha`.
#' @export
sird_model <- function() {
  new_jump_model(
    compartments = c("s", "i", "r"),
    events = list(
      infection = list(jump = c(-1,  1, 0), rate = function(x, p) p[["lambda"]] * x[1] * x[2]),
      recovery  = list(jump = c( 0, -1, 1), rate = function(x, p) p[["gamma"]]  * x[2]),
      death     = list(jump = c( 0, -1, 0), rate = function(x, p) p[["alpha"]]  * x[2])
    )
  )
}
