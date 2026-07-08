## R/cumulant-c2222.R
## --------------------------------------------------------------------------
## Closed-form (window-truncation-free) piece of the pruned-SS marginal 4th
## cumulant, following the Andreasen-Fernandez-Villaverde-Rubio-Ramirez
## (2018) pruned state-space hierarchy [AFVRR] and Mutschler (2015) [M2015].
##
## BACKGROUND. compute_fourth_cumulant() (R/cumulant-likelihood.R) writes the
## marginal excess kurtosis of an observable as kappa4 = CHAIN + TRACE, where
## for the centered observable y_i - E[y_i] = a'w + w'Mw (w = innovation
## history, M symmetric):
##   kappa4 = 48 * a'S M S M S a      (CHAIN)
##          + 48 * tr((MS)^4)         (TRACE)
## The CHAIN term is ALREADY exact closed form (.fourth_cumulant_chain_closed,
## no truncation) because the "M x2-accumulation" part of M contracts against
## the linear part `a` through the stationary tensor C2211 = cum(x2,x2,x1,x1)
## (solved via a tensor-Lyapunov equation, .fourth_cross_cumulant).
##
## The TRACE term (48*tr((MS)^4)) is still computed by summing a truncated
## innovation-history window (.fourth_cumulant_qform_marginal), because it
## involves M FOUR times, and M = Mc (contemporaneous ghxx/ghxu/ghuu quadratic
## in (x1_{t-1}, e_t)) + Mx (the ghx-weighted x2-accumulation across all past
## lags). Profiling shows this window sum is ~99.8% of compute_fourth_cumulant's
## cost (see .claude memory: wave-2026-07-08-hvp-solution-t2.md).
##
## SCOPE OF THIS FILE (READ BEFORE USE): expanding tr((MS)^4) with M = Mc + Mx
## by multilinearity of the trace gives 6 distinct necklace terms (grouping
## the C(4,k) assignments of Mc/Mx to the 4 trace slots by cyclic symmetry):
##   tr((McS)^4)                                         [0 Mx]  -- FINITE, closed form
##   4*tr(Mx S Mx S Mx S Mc S)                           [1 Mc]  -- needs cum(x2,x2,x2,x1) (C2221)
##   4*tr(Mx S Mx S Mc S Mc S) + 2*tr(Mx S Mc S Mx S Mc S) [2 Mc] -- needs C2211 (adjacent) + cum4(x2) (alternating)
##   4*tr(Mx S Mc S Mc S Mc S)                           [3 Mc]  -- needs cum(x2,x1,x1,x1) (=0, x1 Gaussian, vanishes)
##   tr((MxS)^4)                                         [4 Mc=0] -- needs cum4(x2) = C2222
##
## Deriving C2222 = cum(x2,x2,x2,x2) and C2221 = cum(x2,x2,x2,x1) in closed
## form requires a further tensor-Lyapunov closure whose RHS recurses into a
## GENUINELY NONZERO 5th-order object cum(x2,x2,x2,x1,x1) (three quadratic-in-
## Gaussian x2 "necklace" factors do not Wick-cancel the way the 2-quadratic
## C2211 RHS does). This was investigated previously (2026-07-01, orchestrator
## + a delegated agent) and assessed as a genuinely large, error-prone,
## multi-tensor research build -- NOT attempted here (see .claude memory:
## cumulant-third-cumulant-skewness-bug.md, "the FULL Andreasen-FV-RR 4th-
## moment hierarchy... genuinely a major standalone build, NOT one more
## tensor").
##
## WHAT THIS FILE SHIPS: the CHAIN term (delegated, exact, unchanged) plus a
## CONTEMPORANEOUS-ONLY closed form for TRACE, i.e. just the "[0 Mx]" term
## above, tr((McS)^4). This OMITS the x2-self and all x2/Mc mixed trace terms.
## Historical MC evidence (cumulant-third-cumulant-skewness-bug.md) shows the
## omitted terms can be the DOMINANT part of the trace itself (~93% in one
## persistent-state case) -- but the trace as a WHOLE is only ~0.1%-1.3% of
## total kappa4 (the chain dominates ~99.5%+), so the omission's impact on
## the FULL marginal kappa4 is small and, per the oracle test in
## tests/testthat/test-cumulant-c2222.R, may fall inside Monte-Carlo noise for
## moderately-persistent models. It should NOT be assumed to hold near the
## unit root (large rho(hx)) where the trace's share of kappa4 grows.
##
## Given this, .fourth_cumulant_closed_form() is a PARTIAL / SCOPED closed
## form: it removes the O(n_lag) dense-matrix trace loop entirely (genuine
## O(1) cost after the existing C211/C2211 tensor solves), at the cost of
## dropping the x2-driven trace corrections. compute_fourth_cumulant()'s
## `method = "closed_form"` option routes here; default stays `"window"`
## (see R/cumulant-likelihood.R) because the window-limit oracle (closed
## form vs .fourth_cumulant_qform_marginal at large n_lag) does NOT match --
## by construction, since the omitted terms are real and nonzero.
## --------------------------------------------------------------------------


#' Contemporaneous-only closed form of the 4th-cumulant TRACE term
#'
#' Computes \eqn{48 \cdot \mathrm{tr}((M_{c,i} S_g)^4)} for each observable
#' \eqn{i}, where \eqn{M_{c,i}} is the observable's OWN contemporaneous
#' quadratic form in \eqn{(x^{(1)}_{t-1}, e_t)} (from \code{ghxx}/\code{ghxu}/
#' \code{ghuu}) and \eqn{S_g = \mathrm{blockdiag}(\Sigma_x, \Sigma_e)}.
#'
#' This is the exact closed-form contribution of the "[0 Mx]" necklace term
#' in the tr((Mc+Mx)S)^4 expansion (see file header); it OMITS all terms
#' involving the x2-accumulation operator Mx (i.e. cum4(x2) and its cross
#' terms with Mc), which require the (not implemented) full AFVRR 4th-
#' cumulant hierarchy.
#'
#' @param ghx,ghu,ghxx,ghxu,ghuu Endo-row decision-rule blocks (as in
#'   \code{compute_fourth_cumulant}).
#' @param hx n_s x n_s state transition (only used for n_s/dimension).
#' @param Sigma_x,Sigma_e Stationary state / shock covariance.
#' @return n_endo numeric vector: the contemporaneous-only trace contribution
#'   to kappa4, per observable.
#' @noRd
.fourth_cumulant_trace_contemp_closed <- function(ghx, ghu, ghxx, ghxu, ghuu,
                                                  hx, Sigma_x, Sigma_e) {
  n_endo <- nrow(ghx); n_s <- nrow(hx); n_exo <- nrow(Sigma_e)
  d  <- n_s + n_exo
  Sg <- matrix(0, d, d)
  Sg[seq_len(n_s), seq_len(n_s)] <- Sigma_x
  Sg[(n_s + 1L):d, (n_s + 1L):d] <- Sigma_e
  bi <- seq_len(n_s); ei <- (n_s + 1L):d

  out <- numeric(n_endo)
  for (i in seq_len(n_endo)) {
    Mc <- matrix(0, d, d)
    Mc[bi, bi] <- 0.5 * matrix(ghxx[i, ], n_s, n_s)
    Xui <- matrix(ghxu[i, ], n_exo, n_s, byrow = TRUE)
    Mc[bi, ei] <- 0.5 * t(Xui); Mc[ei, bi] <- 0.5 * Xui
    Mc[ei, ei] <- 0.5 * matrix(ghuu[i, ], n_exo, n_exo, byrow = TRUE)
    Mc <- (Mc + t(Mc)) / 2

    MS  <- Mc %*% Sg
    MS2 <- MS %*% MS
    out[i] <- 48 * sum(MS2 * t(MS2))    # 48 * tr((Mc Sg)^4)
  }
  out
}


#' Partial closed-form marginal 4th cumulant of pruned observables
#'
#' Same overall shape/quantity as \code{compute_fourth_cumulant()} (marginal
#' excess-kurtosis \code{kurtosis_obs} and the diagonal-filled \code{c4_obs}
#' tensor), but computed WITHOUT the innovation-history truncation window
#' for the CHAIN term (exact, delegated to \code{.fourth_cumulant_chain_closed})
#' and with a CONTEMPORANEOUS-ONLY closed form for the TRACE term (see file
#' header for exactly which cross-terms are omitted and why).
#'
#' @param dr2 A \code{DecisionRules2} (order-2) object.
#' @param model Parsed model.
#' @param params Named numeric parameter vector (default \code{model$param_values}).
#' @param obs_vars Unused (kept for interface symmetry with the window path;
#'   the underlying computation is always over all endogenous variables).
#' @return List with \code{kurtosis_obs}, \code{c4_obs} (same shape as
#'   \code{compute_fourth_cumulant}), plus \code{kappa4_chain} (exact),
#'   \code{kappa4_trace_contemp} (partial trace), and a \code{scope_note}
#'   string documenting the omission.
#' @noRd
.fourth_cumulant_closed_form <- function(dr2, model, params = NULL,
                                         obs_vars = NULL) {
  endo      <- dr2$endo_names
  exo       <- dr2$exo_names
  state_idx <- dr2$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)

  if (is.null(params)) params <- model$param_values

  ghx  <- dr2$ghx;  ghu  <- dr2$ghu
  ghxx <- dr2$ghxx %||% NULL
  ghxu <- dr2$ghxu %||% NULL
  ghuu <- dr2$ghuu %||% NULL

  shock_stderr <- .get_shock_stderr(model, exo, params)
  Sigma_e <- diag(shock_stderr^2, n_exo)

  hx  <- ghx [state_idx, , drop = FALSE]
  hu  <- ghu [state_idx, , drop = FALSE]
  hxx <- if (!is.null(ghxx)) ghxx[state_idx, , drop = FALSE] else NULL
  hxu <- if (!is.null(ghxu)) ghxu[state_idx, , drop = FALSE] else NULL
  huu <- if (!is.null(ghuu)) ghuu[state_idx, , drop = FALSE] else NULL

  Sigma_x <- .state_covariance(hx, hu, Sigma_e)

  kurtosis <- rep(0, n_endo); names(kurtosis) <- endo
  c4_obs   <- matrix(0, n_endo, n_endo * n_endo * n_endo)
  kappa4_chain         <- rep(0, n_endo); names(kappa4_chain) <- endo
  kappa4_trace_contemp <- rep(0, n_endo); names(kappa4_trace_contemp) <- endo

  if (n_s > 0L && !is.null(hxx) && !is.null(hxu) && !is.null(huu)) {
    C211_4 <- .solve_third_cross_cumulant(
      hx, .third_cumulant_rhs(hxx, Sigma_x, hxu = hxu, huu = huu,
                              hu = hu, hx = hx, Sigma_e = Sigma_e))
    C2211 <- .fourth_cross_cumulant(hx, hu, hxx, hxu, huu,
                                    Sigma_x, Sigma_e, C211_4)
    chain <- .fourth_cumulant_chain_closed(
      ghx, ghu, ghxx, ghxu, ghuu, hx, Sigma_x, Sigma_e, C211_4, C2211)
    trace_contemp <- .fourth_cumulant_trace_contemp_closed(
      ghx, ghu, ghxx, ghxu, ghuu, hx, Sigma_x, Sigma_e)

    kappa4_vec <- chain + trace_contemp
    kappa4_chain <- chain
    kappa4_trace_contemp <- trace_contemp

    var_total <- diag(compute_moments(dr2, model, params = params)$var_cov)

    for (i in seq_len(n_endo)) {
      k2 <- var_total[i]
      if (k2 > 1e-30) {
        kurtosis[i] <- kappa4_vec[i] / k2^2
        col_idx <- 1 + (i - 1L) + n_endo * (i - 1L) + n_endo^2 * (i - 1L)
        if (col_idx <= ncol(c4_obs)) c4_obs[i, col_idx] <- kappa4_vec[i]
      }
    }
  }

  list(
    kurtosis_obs = kurtosis,
    c4_obs       = c4_obs,
    kappa4_chain = kappa4_chain,
    kappa4_trace_contemp = kappa4_trace_contemp,
    scope_note = paste0(
      "PARTIAL closed form (R/cumulant-c2222.R): chain term is EXACT ",
      "(no truncation, delegated to .fourth_cumulant_chain_closed); trace ",
      "term is CONTEMPORANEOUS-ONLY (48*tr((Mc*Sigma_g)^4)) and OMITS the ",
      "x2-self (cum4(x2), i.e. the full C2222/AFVRR hierarchy) and x2-Mc ",
      "cross trace terms. Not exact at any rho(hx); do not use as a drop-in ",
      "replacement near the unit root. See the file header for the full ",
      "necklace-term accounting.")
  )
}
