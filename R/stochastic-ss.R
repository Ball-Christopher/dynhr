## R/stochastic-ss.R
## --------------------------------------------------------------------------
## Stochastic steady state (SSS) and Generalized Impulse Response Functions
## (GIRFs) for pruned second-order DSGE models.
##
## Dependencies: .order2_aug_system, .order2_stationary_moments,
##               .order2_conditional_moments, .get_shock_cov, .get_shock_stderr
##               (all in R/stochsimul-monolith.R — no cross-file edits needed)
##
## References:
##   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F. (2018).
##     The pruned state-space system for non-linear DSGE models.
##     Review of Economic Studies, 85(1), 1-49.
##   Koop, G., Pesaran, M. H., & Potter, S. M. (1996). Impulse response
##     analysis in nonlinear multivariate models. Journal of Econometrics,
##     74(1), 119-147.
## --------------------------------------------------------------------------


## ============================================================================
## 1. STOCHASTIC STEADY STATE
## ============================================================================

#' Stochastic (ergodic / risky) steady state for a pruned order-2 DSGE model
#'
#' Returns the fixed point of \eqn{E[y_t]} under the ergodic distribution of
#' the pruned second-order state-space system with Gaussian shocks (Andreasen,
#' Fernandez-Villaverde & Rubio-Ramirez 2018).  This is a closed-form result
#' via the Lyapunov fixed point; no simulation or fixed-point iteration is
#' needed.  For a linear (order-1) model, the stochastic SS equals the
#' deterministic SS (zero risk correction).
#'
#' The SSS solves the fixed point of the augmented pruned-state mean:
#' \deqn{\mu_\xi = (I - T)^{-1}(c + c_u)}
#' where \eqn{T}, \eqn{c}, \eqn{c_u} come from the augmented system
#' (see \code{.order2_aug_system}).  The level SSS is then
#' \deqn{\text{sss} = y^* + D_\xi \mu_\xi + \tfrac{1}{2} g_{ss} + c_v}
#' This equals \code{compute_moments_order2(dr, model)$mean} to machine
#' precision (Gate 3 in the test suite verifies this).
#'
#' Non-Gaussian shocks (PSKF skewness): the Gaussian ergodic mean formula is
#' approximate under skewness; a warning is emitted via
#' \code{.order2_stationary_moments} in that case (inherited verbatim).
#'
#' @param dr DecisionRules2 object returned by \code{solve_perturbation(order=2)}
#' @param model dynhr_mod object
#' @param params Optional named numeric parameter vector (default:
#'   \code{model$param_values})
#' @return Named numeric vector of length \code{n_endo} giving the stochastic
#'   steady state in level units (same variables and ordering as
#'   \code{dr$endo_names}).  The attribute \code{"risk_correction"} holds
#'   \eqn{\text{sss} - y^*}: the deviation from the deterministic steady state,
#'   also known as the precautionary risk premium.
#' @seealso \code{\link{compute_moments_order2}} (returns the same mean),
#'   \code{\link{compute_girf}} (GIRF from the SSS)
#' @export
stochastic_steady_state <- function(dr, model, params = NULL) {
  if (!inherits(dr, "DecisionRules2")) {
    stop(
      "stochastic_steady_state: requires a DecisionRules2 object ",
      "(solve_perturbation(order=2)). Got class: ",
      paste(class(dr), collapse = ", ")
    )
  }
  if (is.null(params)) params <- model$param_values

  Sigma_e <- .get_shock_cov(model, dr$exo_names, params)
  sys     <- .order2_aug_system(dr, Sigma_e)
  st      <- .order2_stationary_moments(sys)

  sss <- st$mean
  names(sss) <- dr$endo_names
  attr(sss, "risk_correction") <- sss - dr$ys
  sss
}


## ============================================================================
## 2. GENERALIZED IMPULSE RESPONSE FUNCTIONS FROM THE SSS
## ============================================================================

#' Generalized Impulse Response Functions (GIRFs) from the stochastic SS
#'
#' Computes the Koop-Pesaran-Potter (1996) GIRF from the stochastic steady
#' state (SSS) of a pruned second-order DSGE model, or from a caller-supplied
#' initial state for state-dependent GIRF exercises:
#' \deqn{\text{GIRF}_h(k) = E[y_{t+h} \mid \varepsilon_t = e_k,\;
#'   x_t = x_0]
#'   - E[y_{t+h} \mid x_t = x_0]}
#' for each shock \eqn{k}, where \eqn{e_k} is the \eqn{k}-th column of the
#' lower Cholesky factor of \eqn{\Sigma_\varepsilon} (scaled by
#' \code{shock_size}) and \eqn{x_0} is the initial state (defaults to the
#' stochastic steady state).
#'
#' Future shocks are integrated out analytically using
#' \code{.order2_conditional_moments} (augmented Lyapunov propagation); no
#' Monte Carlo is needed.  At order 1 (all quadratic terms zero), GIRFs equal
#' standard IRFs from the deterministic SS.
#'
#' Period-indexing: row h of the returned matrix corresponds to h periods
#' after the shock (h = 1 is the impact period, the same convention as
#' \code{compute_irfs()}).  The h = 1 impact is computed explicitly via the
#' output policy functions:
#' \deqn{\text{GIRF}_1(k) = g_{hu} e_k + g_{hxu}(e_k \otimes x^*_1)
#'   + \tfrac{1}{2} g_{huu}(e_k \otimes e_k)}
#' where \eqn{x^*_1} is the first-order state deviation from the deterministic
#' SS (from \code{initial_state}, or the ergodic mean when
#' \code{initial_state = NULL}).
#' Periods h = 2, \ldots, n use analytic Lyapunov propagation from the
#' post-shock state \eqn{x^*_1 + h_u e_k}; the second-order correction
#' \eqn{E[x^{(2)}]} is initialized to 0 in the propagator (see
#' \code{.order2_conditional_moments}), introducing an \eqn{O(\sigma^2)}
#' approximation that is absorbed into the \eqn{c_u / c_v} constants after
#' a few steps.  This error is negligible for typical DSGE calibrations
#' (shock stderr 1--3 \%); full augmented-state initialization is deferred.
#'
#' @param dr DecisionRules2 object returned by \code{solve_perturbation(order=2)}
#' @param model dynhr_mod object
#' @param n_periods Horizon length (default 40)
#' @param shock_size Shock magnitude in units of the shock's standard
#'   deviation (default 1 = a one-standard-deviation impulse). The impact is
#'   \code{shock_size} times the relevant column of \code{chol(Sigma_e)}, so
#'   \code{shock_size = 2} is a 2-sigma impulse, not 2 in level units.
#' @param params Optional named numeric parameter vector (default:
#'   \code{model$param_values})
#' @param initial_state Optional named numeric vector over the state variables
#'   (names matching \code{dr$endo_names[dr$state_idx]}, in level units) from
#'   which the GIRF is computed.  When \code{NULL} (default), the function
#'   starts from the stochastic steady state, reproducing the standard ergodic
#'   GIRF.  Supply a different state to obtain a state-dependent GIRF (the
#'   RBC_state_dependent_GIRF use case): the GIRF then measures the response
#'   relative to the baseline path that starts from the same \code{initial_state}
#'   without a shock.
#' @param shock_name Optional character string naming a single shock
#'   (must appear in \code{dr$exo_names}).  When non-\code{NULL}, only that
#'   shock's GIRF is computed and the returned \code{IRFCollection} contains
#'   a single entry.  \code{NULL} (default) computes all shocks.
#' @param n_replications Optional integer.  \code{NULL} (default) uses the
#'   analytic pruned-state-space GIRF, which integrates future shocks out
#'   exactly (exact under Gaussianity at order 2).  An integer \code{>= 1}
#'   instead estimates the GIRF by Monte Carlo (Koop-Pesaran-Potter): the
#'   pruned simulation is averaged over \code{n_replications} future shock
#'   paths for both the shocked and the (shock-free-impact) baseline
#'   expectation.  The MC path is slower but assumption-light -- use it when
#'   the analytic integration does not apply (order 3+, non-Gaussian shocks,
#'   occasionally-binding-constraint / regime models).
#' @param seed Optional integer seed for the Monte-Carlo path
#'   (\code{n_replications} non-\code{NULL}); makes the MC GIRF reproducible.
#'   Ignored on the analytic path.
#' @return An \code{IRFCollection} object (same structure as
#'   \code{compute_irfs()}): a named list with one entry per (selected) shock;
#'   each entry is a \code{n_periods x n_endo} numeric matrix of GIRF paths
#'   (deviation from the baseline path starting at \code{initial_state}).
#'   Extra attributes: \code{type = "GIRF"},
#'   \code{sss} = the stochastic steady state vector,
#'   \code{initial_state} = the state used as the GIRF starting point,
#'   \code{order = 2L}, \code{n_periods}, \code{endo_names}, \code{exo_names}.
#'   The object is compatible with \code{plot_irfs()} and any downstream code
#'   that accepts \code{IRFCollection}.
#' @seealso \code{\link{stochastic_steady_state}},
#'   \code{\link{compute_irfs}} (standard linear IRF from det-SS),
#'   \code{\link{compute_irfs_order2}} (non-linear IRF from det-SS)
#' @export
compute_girf <- function(dr, model, n_periods = 40L, shock_size = 1,
                         params = NULL, initial_state = NULL,
                         shock_name = NULL, n_replications = NULL,
                         seed = NULL) {
  if (!inherits(dr, "DecisionRules2")) {
    stop(
      "compute_girf: requires a DecisionRules2 object ",
      "(solve_perturbation(order=2)). Got class: ",
      paste(class(dr), collapse = ", ")
    )
  }
  if (is.null(params)) params <- model$param_values
  n_periods <- as.integer(n_periods)

  ## --- n_replications: switch between the analytic and Monte-Carlo paths ---
  ## NULL (default) -> analytic pruned-state-space GIRF (future shocks
  ##   integrated out exactly; exact under Gaussianity at order 2).
  ## integer >= 1   -> Koop-Pesaran-Potter Monte-Carlo GIRF: average the
  ##   pruned simulation over n_replications future shock paths. Needed when
  ##   the analytic integration does not apply (order 3+, non-Gaussian shocks,
  ##   OBC/regime models). Slower but assumption-light.
  use_mc <- !is.null(n_replications)
  if (use_mc) {
    n_replications <- as.integer(n_replications)
    if (length(n_replications) != 1L || is.na(n_replications) ||
        n_replications < 1L) {
      stop("compute_girf: n_replications must be a single integer >= 1 (or NULL).")
    }
  }

  Sigma_e <- .get_shock_cov(model, dr$exo_names, params)
  sys     <- .order2_aug_system(dr, Sigma_e)
  st      <- .order2_stationary_moments(sys)
  sss     <- st$mean
  names(sss) <- dr$endo_names

  state_idx  <- dr$state_idx
  state_names <- dr$endo_names[state_idx]

  ## --- initial_state: resolve the first-order state deviation from det-SS ---
  ## x1_start is a length-n_s vector (deviation from det-SS ys).
  if (is.null(initial_state)) {
    ## Default: stochastic steady state
    x1_start <- (sss - dr$ys)[state_idx]
  } else {
    ## Caller supplied a state in level units.  Accept either a full n_endo
    ## vector (named over endo_names) or a length-n_s vector over state vars.
    if (!is.numeric(initial_state)) {
      stop("compute_girf: initial_state must be a numeric vector.")
    }
    if (length(initial_state) == length(state_idx)) {
      ## Compact form: treat as the state-variable subvector directly.
      if (!is.null(names(initial_state))) {
        if (!all(names(initial_state) %in% state_names)) {
          stop(
            "compute_girf: initial_state names must be state variable names. ",
            "Expected subset of: ", paste(state_names, collapse = ", ")
          )
        }
        ## Reorder to match state_idx ordering.
        x1_start <- initial_state[state_names] - dr$ys[state_idx]
      } else {
        x1_start <- initial_state - dr$ys[state_idx]
      }
    } else if (length(initial_state) == length(dr$endo_names)) {
      ## Full n_endo vector: extract state subvector.
      if (!is.null(names(initial_state))) {
        x1_start <- (initial_state - dr$ys)[state_idx]
      } else {
        x1_start <- (initial_state - dr$ys)[state_idx]
      }
    } else {
      stop(
        "compute_girf: initial_state has length ", length(initial_state),
        " but must have length n_state (", length(state_idx),
        ") or n_endo (", length(dr$endo_names), ")."
      )
    }
  }

  endo   <- dr$endo_names
  exo    <- dr$exo_names
  n_endo <- length(endo)
  n_exo  <- length(exo)
  n_s    <- length(state_idx)

  ## --- shock_name: select shocks to compute ---
  if (!is.null(shock_name)) {
    if (!is.character(shock_name) || length(shock_name) != 1L) {
      stop("compute_girf: shock_name must be a single character string or NULL.")
    }
    if (!(shock_name %in% exo)) {
      stop(
        "compute_girf: shock_name '", shock_name,
        "' not found in dr$exo_names. Available: ",
        paste(exo, collapse = ", ")
      )
    }
    shock_idx <- which(exo == shock_name)
  } else {
    shock_idx <- seq_len(n_exo)
  }

  ## Sub-matrices for state variables only
  hx  <- dr$ghx[state_idx, , drop = FALSE]   # n_s x n_s
  hu  <- dr$ghu[state_idx, , drop = FALSE]   # n_s x n_exo
  hxu <- dr$ghxu[state_idx, , drop = FALSE]  # n_s x n_s*n_exo

  ## Cholesky factor for shock scaling (same convention as compute_irfs):
  ## IRF for shock k uses the k-th column of chol(Sigma_e) (lower triangular).
  L_chol <- tryCatch(
    t(chol(Sigma_e)),
    error = function(e) {
      diag(sqrt(pmax(diag(Sigma_e), 0)), nrow = n_exo)
    }
  )

  irfs <- vector("list", length(shock_idx))
  names(irfs) <- exo[shock_idx]

  if (use_mc) {
    ## ======================================================================
    ## Monte-Carlo (Koop-Pesaran-Potter) GIRF path.
    ## GIRF(h) = E[y_{t+h} | x_t = x_start, eps_t = e_k]
    ##         - E[y_{t+h} | x_t = x_start]            (eps_t random)
    ## with future shocks (h>=2) drawn from N(0, Sigma_e) in both expectations.
    ## Matches the validated Gate-7 oracle convention exactly.
    ## ======================================================================
    if (!is.null(seed)) set.seed(as.integer(seed))

    ## Pruned-recursion state matrices (state rows only).
    hxx <- dr$ghxx[state_idx, , drop = FALSE]
    huu <- dr$ghuu[state_idx, , drop = FALSE]
    hss <- dr$ghss[state_idx]
    ghx <- dr$ghx; ghu <- dr$ghu; ghxx <- dr$ghxx
    ghuu <- dr$ghuu; ghss <- dr$ghss; ghxu <- dr$ghxu

    pruned_output <- function(x1p, x2p, e) {
      as.numeric(ghx %*% x1p + ghu %*% e) +
        as.numeric(ghx %*% x2p) +
        0.5 * as.numeric(ghxx %*% (x1p %x% x1p)) +
        as.numeric(ghxu %*% (e %x% x1p)) +
        0.5 * as.numeric(ghuu %*% (e %x% e)) +
        0.5 * ghss
    }
    pruned_step <- function(x1p, x2p, e) {
      list(
        x1 = as.numeric(hx %*% x1p + hu %*% e),
        x2 = as.numeric(hx %*% x2p +
                          0.5 * hxx %*% (x1p %x% x1p) +
                          hxu %*% (e %x% x1p) +
                          0.5 * huu %*% (e %x% e) +
                          0.5 * hss)
      )
    }
    ## Pre-draw the future shock sequences (h >= 2), shared across the baseline
    ## and impulse paths.  Each replication uses IDENTICAL draws for both paths
    ## (common random numbers) so the common stochastic variation cancels in the
    ## difference, leaving only the low-variance impulse signal.
    ##
    ## IMPACT (h=1) convention — must match the analytic path (below):
    ##   GIRF(h) = E[y_{t+h} | x_t, eps_t = e_k] - E[y_{t+h} | x_t, eps_t = 0].
    ## The baseline impact is the DETERMINISTIC zero shock, NOT a random draw.
    ## Using a random baseline impact (eps_t ~ N) makes the estimator converge to
    ##   analytic - 0.5*ghuu*vec(Sigma_e)
    ## i.e. it spuriously subtracts the unconditional sigma^2 correction (the
    ## ghuu*Sigma term), biasing every GIRF.  With baseline impact = 0 the h=1
    ## difference is exact (deterministic) and h>=2 share the future draws, so the
    ## Monte-Carlo GIRF is an unbiased, low-variance estimator of the analytic one.
    ##
    ## Layout: eps_draws[rep, h, ] is the n_exo-vector for replication rep,
    ## period h.  Only h >= 2 rows are used (shared future shocks).

    eps_draws <- array(
      stats::rnorm(n_replications * n_periods * n_exo),
      dim = c(n_replications, n_periods, n_exo)
    )
    e_zero <- numeric(n_exo)
    ## Scale by Cholesky once: eps_draws[rep, h, ] <- L_chol %*% z
    for (rep in seq_len(n_replications)) {
      for (h in seq_len(n_periods)) {
        eps_draws[rep, h, ] <- as.numeric(L_chol %*% eps_draws[rep, h, ])
      }
    }

    for (ki in seq_along(shock_idx)) {
      k   <- shock_idx[ki]
      e_k <- as.numeric(L_chol[, k]) * shock_size

      ## Accumulate the PAIRED difference impulse - baseline per replication.
      acc_diff <- matrix(0, n_periods, n_endo)
      x2_0 <- numeric(n_s)

      for (rep in seq_len(n_replications)) {
        ## --- Baseline path: zero impact shock (h=1), shared future (h>=2) ---
        x1b <- x1_start; x2b <- x2_0
        ## --- Impulse path: structured impact e_k (h=1), shared future (h>=2) ---
        x1i <- x1_start; x2i <- x2_0

        for (h in seq_len(n_periods)) {
          e_shared <- eps_draws[rep, h, ]              # shared draw (future only)
          e_base   <- if (h == 1L) e_zero else e_shared  # baseline: 0 at impact
          e_imp    <- if (h == 1L) e_k    else e_shared  # impulse: e_k at impact

          acc_diff[h, ] <- acc_diff[h, ] +
            pruned_output(x1i, x2i, e_imp) - pruned_output(x1b, x2b, e_base)

          stb <- pruned_step(x1b, x2b, e_base)
          sti <- pruned_step(x1i, x2i, e_imp)
          x1b <- stb$x1; x2b <- stb$x2
          x1i <- sti$x1; x2i <- sti$x2
        }
      }

      girf_mat <- acc_diff / n_replications
      dimnames(girf_mat) <- list(paste0("t", seq_len(n_periods)), endo)
      irfs[[ki]] <- girf_mat
    }
  } else {
    ## ======================================================================
    ## Analytic pruned-state-space GIRF path (default).
    ## ======================================================================
    n_tail <- n_periods - 1L
    if (n_tail > 0L) {
      cm_base_tail <- .order2_conditional_moments(sys, x1_start, n_tail)
    }

    for (ki in seq_along(shock_idx)) {
      k <- shock_idx[ki]
      ## Shock vector: k-th Cholesky column, scaled
      e_k <- as.numeric(L_chol[, k]) * shock_size

      ## --- h = 1: explicit impact ---
      ## GIRF_1 = output(x1_start + shock) - output(x1_start):
      ##   y1_diff = ghu * e_k               (first-order)
      ##   y2_diff = ghxu*(e_k⊗x1_start) + 0.5*ghuu*(e_k⊗e_k)  (quadratic cross)
      ## ghss and ghx*x2 cancel between shocked and baseline.
      girf_1 <- as.numeric(dr$ghu %*% e_k) +
                as.numeric(dr$ghxu %*% (e_k %x% x1_start)) +
                0.5 * as.numeric(dr$ghuu %*% (e_k %x% e_k))

      girf_mat <- matrix(0, n_periods, n_endo,
                         dimnames = list(paste0("t", seq_len(n_periods)), endo))
      girf_mat[1L, ] <- girf_1

      if (n_tail > 0L) {
        ## --- h = 2..n_periods: analytic propagation from post-shock state ---
        ## x1 after one pruned step: hx*x1_start + hu*e_k
        x1_after <- as.numeric(hx %*% x1_start + hu %*% e_k)

        cm_shock_tail <- .order2_conditional_moments(sys, x1_after, n_tail)

        girf_mat[2L:n_periods, ] <- cm_shock_tail$mean - cm_base_tail$mean
      }

      irfs[[ki]] <- girf_mat
    }
  }

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")     <- n_periods
  attr(irfs, "endo_names")    <- endo
  attr(irfs, "exo_names")     <- exo[shock_idx]
  attr(irfs, "order")         <- 2L
  attr(irfs, "type")          <- "GIRF"
  attr(irfs, "method")        <- if (use_mc) "monte-carlo" else "analytic"
  if (use_mc) attr(irfs, "n_replications") <- n_replications
  attr(irfs, "sss")           <- sss
  attr(irfs, "initial_state") <- x1_start + dr$ys[state_idx]
  irfs
}
