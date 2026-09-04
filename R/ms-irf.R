## R/ms-irf.R
## --------------------------------------------------------------------------
## ms_irf() -- impulse responses for Markov-switching DSGE models.
##
## Three objects are produced for each shock:
##
##   1. CONDITIONAL ("staying in regime r"): the regime is held at r for the
##      whole horizon.  This is the textbook regime-r IRF and is computed by
##      the ordinary fixed-regime recursion
##          y_1 = ghu_r eps_r ,  y_h = ghx_r y_{h-1}[state]
##      i.e. exactly what compute_irfs() does on the regime-r model.
##
##   2. EXPECTED ("starting in regime r, regimes then evolve under P"): the
##      Markov-averaged response
##          IRF_h = E[ y_h | s_1 = r ]
##      computed by carrying the JOINT (response, regime) accumulator
##          m_h^{(j)} = E[ y_h 1{s_h = j} | s_1 = r ]
##          m_1^{(j)} = 1{j = r} ghu_r eps_r
##          m_h^{(j)} = ghx_j ( sum_i P[i,j] m_{h-1}^{(i)}[state] ) ,  h >= 2
##          IRF_h     = sum_j m_h^{(j)}
##      This is the linear (order-1) case of the Koop-Pesaran-Potter
##      generalized IRF: because the state-space is linear given the regime
##      path, the expectation over future regime paths is exact -- no
##      simulation is needed, and the recursion above IS that expectation.
##
##   3. ERGODIC: sum_r pi0_r * EXPECTED_r, i.e. the response to a shock hitting
##      at a random date with the regime drawn from pi0.  Note the impact
##      vector itself is regime-dependent (a one-standard-deviation shock is
##      LARGER in a high-volatility regime), which is why the ergodic IRF is
##      not simply the ergodic average of the conditional IRFs.
##
## ORACLES (test-ms-irf.R):
##   * CONDITIONAL[[r]] == compute_irfs() on the regime-r decision rules.
##     Different function, different code path -- a genuine cross-check.
##   * With P = I, EXPECTED[[r]] == CONDITIONAL[[r]].  The expected path runs
##     the full Markov accumulator, so this is NOT construction-circular: it
##     pins that the P[i,j] bookkeeping degenerates correctly.
##   * With P = I, ERGODIC == sum_r pi0_r * CONDITIONAL[[r]].
## --------------------------------------------------------------------------


## Internal: lower-Cholesky impact factor of Sigma_e, matching compute_irfs().
## @noRd
.ms_irf_impact_factor <- function(Sigma_e, n_exo) {
  tryCatch(
    t(chol(Sigma_e)),
    error = function(e)
      diag(sqrt(pmax(diag(Sigma_e), 0)), nrow = n_exo)
  )
}


## Internal: build the per-regime IRF ingredients (ghx, ghu, Sigma_e) for
## either an ms_dsge_spec (shock-variance switching, common decision rules)
## or an MsDecisionRules object (structural switching).
## @noRd
.ms_irf_regimes <- function(dr, model, ms_spec, params) {
  if (inherits(dr, "MsDecisionRules")) {
    h    <- length(dr$dr)
    dr1  <- dr$dr[[1L]]
    exo  <- dr1$exo_names
    Se   <- .get_shock_cov(model, exo, params)
    list(
      h         = h,
      P         = dr$P,
      pi0       = dr$pi0,
      names     = if (!is.null(names(dr$dr))) names(dr$dr)
                  else paste0("regime", seq_len(h)),
      endo      = dr1$endo_names,
      exo       = exo,
      state_idx = dr1$state_idx,
      ghx       = lapply(dr$dr, `[[`, "ghx"),
      ghu       = lapply(dr$dr, `[[`, "ghu"),
      ## Per-regime Sigma_e: the structural solver switches the decision rules,
      ## not the shock block, so the common Sigma_e is used for every regime.
      Sigma_e   = replicate(h, Se, simplify = FALSE),
      sub_dr    = dr$dr
    )
  } else {
    if (!inherits(ms_spec, "ms_dsge_spec"))
      stop("ms_irf: ms_spec must be an ms_dsge_spec object when `dr` is a ",
           "single DecisionRules object. Pass an MsDecisionRules object ",
           "instead for structural (decision-rule) switching.", call. = FALSE)
    exo   <- dr$exo_names
    n_exo <- length(exo)
    Se    <- .get_shock_cov(model, exo, params)
    h     <- ms_spec$n_regimes

    scales <- ms_spec$shock_scales
    sc1    <- scales[[1L]]
    if (!is.null(names(sc1))) {
      if (length(sc1) != n_exo || !identical(sort(names(sc1)), sort(exo)))
        stop(sprintf(
          "ms_irf: shock_scales names (%s) do not match model exo names (%s).",
          paste(names(sc1), collapse = ","), paste(exo, collapse = ",")),
          call. = FALSE)
      scales <- lapply(scales, function(v) v[exo])
    } else if (length(sc1) != n_exo) {
      stop(sprintf("ms_irf: shock_scales length (%d) != n_exo (%d).",
                   length(sc1), n_exo), call. = FALSE)
    }

    Se_list <- lapply(scales, function(sc) Se * outer(sc, sc))
    sub_dr  <- lapply(Se_list, function(Se_r) { d <- dr; d$Sigma_e <- Se_r; d })

    list(
      h         = h,
      P         = ms_spec$transition,
      pi0       = ms_spec$pi0,
      names     = ms_spec$regime_names,
      endo      = dr$endo_names,
      exo       = exo,
      state_idx = dr$state_idx,
      ghx       = replicate(h, dr$ghx, simplify = FALSE),
      ghu       = replicate(h, dr$ghu, simplify = FALSE),
      Sigma_e   = Se_list,
      sub_dr    = sub_dr
    )
  }
}


## Internal: wrap a list of n_periods x n_endo matrices as an IRFCollection.
## @noRd
.ms_irf_collection <- function(mats, endo, exo, n_periods) {
  names(mats) <- exo
  class(mats) <- "IRFCollection"
  attr(mats, "n_periods")  <- n_periods
  attr(mats, "endo_names") <- endo
  attr(mats, "exo_names")  <- exo
  mats
}


#' Regime-conditional and Markov-averaged impulse responses
#'
#' Computes impulse response functions for a Markov-switching DSGE model:
#' the response conditional on \emph{staying} in a given regime, the response
#' \emph{starting} in a given regime with the regime then evolving under the
#' transition matrix, and the ergodic-weighted response.
#'
#' Works for both switching styles:
#' \itemize{
#'   \item shock-variance switching -- pass the single
#'     \code{\link{solve_perturbation}} decision rules as \code{dr} together
#'     with an \code{\link{ms_dsge_spec}};
#'   \item structural switching -- pass an \code{MsDecisionRules} object
#'     (from \code{\link{solve_ms_perturbation}}) as \code{dr}; \code{ms_spec}
#'     is then unused.
#' }
#'
#' The impact vector for shock \eqn{k} in regime \eqn{r} is the \eqn{k}-th
#' column of the lower Cholesky factor of \eqn{\Sigma_e^{(r)}}, scaled by
#' \code{shock_size} -- the same Dynare convention \code{\link{compute_irfs}}
#' uses.  A one-standard-deviation shock is therefore genuinely larger in a
#' high-volatility regime, which is why \code{ergodic} is a weighted average
#' of \code{expected} and not of \code{conditional}.
#'
#' @param dr  A \code{DecisionRules} object (shock-variance switching) or an
#'   \code{MsDecisionRules} object (structural switching).
#' @param model  \code{dynhr_mod} object (for the shock covariance).
#' @param ms_spec  An \code{\link{ms_dsge_spec}}; required when \code{dr} is a
#'   single \code{DecisionRules} object, ignored for \code{MsDecisionRules}.
#' @param n_periods  IRF horizon (default 40).
#' @param shock_size  Shock size in standard-deviation units (default 1).
#' @param params  Named numeric parameter vector (default
#'   \code{model$param_values}).
#'
#' @return A list of class \code{"ms_irf"}:
#'   \describe{
#'     \item{\code{conditional}}{Length-\code{h} named list of
#'       \code{IRFCollection}s: regime \eqn{r} held fixed over the horizon.}
#'     \item{\code{expected}}{Length-\code{h} named list of
#'       \code{IRFCollection}s: impact in regime \eqn{r}, regimes then
#'       evolving under \code{P}.}
#'     \item{\code{ergodic}}{A single \code{IRFCollection}: \code{pi0}-weighted
#'       average of \code{expected}.}
#'     \item{\code{transition}, \code{pi0}, \code{regime_names}}{The MS setup
#'       actually used.}
#'   }
#' @seealso \code{\link{compute_irfs}}, \code{\link{ms_dsge_spec}},
#'   \code{\link{solve_ms_perturbation}}
#' @export
ms_irf <- function(dr, model, ms_spec = NULL, n_periods = 40L,
                    shock_size = 1, params = NULL) {

  n_periods <- as.integer(n_periods)
  if (length(n_periods) != 1L || is.na(n_periods) || n_periods < 1L)
    stop("ms_irf: n_periods must be a positive integer.", call. = FALSE)
  if (is.null(params)) params <- model$param_values

  reg <- .ms_irf_regimes(dr, model, ms_spec, params)

  h         <- reg$h
  P         <- reg$P
  endo      <- reg$endo
  exo       <- reg$exo
  state_idx <- reg$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)

  ## Per-regime impact factors (Dynare convention, matching compute_irfs()).
  L_list <- lapply(reg$Sigma_e, .ms_irf_impact_factor, n_exo = n_exo)

  ## Impact response y_1 for (regime r, shock k).
  impact <- array(0, c(n_endo, n_exo, h))
  for (r in seq_len(h)) {
    for (k in seq_len(n_exo)) {
      eps <- L_list[[r]][, k] * shock_size
      ## Same fallback as compute_irfs(): a shock declared with zero variance
      ## has an all-zero Cholesky column, which would give an identically zero
      ## IRF even though ghu[, k] is perfectly well defined.
      if (all(eps == 0) && shock_size != 0) {
        eps <- numeric(n_exo)
        eps[k] <- shock_size
      }
      impact[, k, r] <- as.numeric(reg$ghu[[r]] %*% eps)
    }
  }

  blank <- function() {
    m <- matrix(0, n_periods, n_endo)
    colnames(m) <- endo
    rownames(m) <- paste0("t", seq_len(n_periods))
    m
  }

  ## ---- 1. CONDITIONAL: regime r held fixed --------------------------------
  conditional <- vector("list", h)
  for (r in seq_len(h)) {
    ghx_r <- reg$ghx[[r]]
    mats  <- vector("list", n_exo)
    for (k in seq_len(n_exo)) {
      m <- blank()
      y <- impact[, k, r]
      m[1L, ] <- y
      if (n_periods > 1L) for (t in seq.int(2L, n_periods)) {
        y <- as.numeric(ghx_r %*% y[state_idx])
        m[t, ] <- y
      }
      mats[[k]] <- m
    }
    conditional[[r]] <- .ms_irf_collection(mats, endo, exo, n_periods)
  }

  ## ---- 2. EXPECTED: Markov accumulator over regime paths ------------------
  expected <- vector("list", h)
  for (r in seq_len(h)) {
    mats <- vector("list", n_exo)
    for (k in seq_len(n_exo)) {
      m <- blank()
      ## mcur[, j] = E[ y_h 1{s_h = j} | s_1 = r ]
      mcur <- matrix(0, n_endo, h)
      mcur[, r] <- impact[, k, r]
      m[1L, ] <- rowSums(mcur)
      if (n_periods > 1L) for (t in seq.int(2L, n_periods)) {
        ## state block of the previous accumulator, pushed through P then ghx_j
        prev_state <- mcur[state_idx, , drop = FALSE]     # n_state x h
        mixed      <- prev_state %*% P                    # n_state x h (col j)
        mnew <- matrix(0, n_endo, h)
        for (j in seq_len(h))
          mnew[, j] <- as.numeric(reg$ghx[[j]] %*% mixed[, j])
        mcur <- mnew
        m[t, ] <- rowSums(mcur)
      }
      mats[[k]] <- m
    }
    expected[[r]] <- .ms_irf_collection(mats, endo, exo, n_periods)
  }

  ## ---- 3. ERGODIC: pi0-weighted average of EXPECTED ------------------------
  pi0 <- reg$pi0
  erg_mats <- vector("list", n_exo)
  for (k in seq_len(n_exo)) {
    acc <- blank()
    for (r in seq_len(h)) acc <- acc + pi0[r] * expected[[r]][[k]]
    erg_mats[[k]] <- acc
  }
  ergodic <- .ms_irf_collection(erg_mats, endo, exo, n_periods)

  names(conditional) <- reg$names
  names(expected)    <- reg$names

  structure(
    list(
      conditional  = conditional,
      expected     = expected,
      ergodic      = ergodic,
      transition   = P,
      pi0          = pi0,
      regime_names = reg$names,
      n_regimes    = h,
      n_periods    = n_periods,
      endo_names   = endo,
      exo_names    = exo
    ),
    class = c("ms_irf", "list")
  )
}


#' @export
#' @noRd
print.ms_irf <- function(x, ...) {
  cat(sprintf("<ms_irf>  %d regimes, %d shocks, %d periods\n",
              x$n_regimes, length(x$exo_names), x$n_periods))
  cat("  Regimes:", paste(x$regime_names, collapse = ", "), "\n")
  cat("  Ergodic weights:", paste(round(x$pi0, 4), collapse = ", "), "\n")
  cat("  Components: $conditional[[r]], $expected[[r]], $ergodic\n")
  invisible(x)
}
