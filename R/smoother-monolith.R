### ===========================================================================
### dynhr_smoother.R -- Kalman smoother + historical decomposition for dynhr
### ===========================================================================
### Provides:
###   build_dsge_state_space()  -- extract compact state-space from m + dr
###   kalman_smoother()         -- RTS backward smoother for DSGE form
###   historical_decomposition() -- per-shock contribution to all endo vars
###
### DSGE state-space form (first-order perturbation):
###   s_t = T s_{t-1} + R eps_t          (n_state x 1)
###   y_t = Z s_{t-1} + D eps_t          (n_obs x 1)
###
### where s_t = [backward; mixed] variables, y_t = observables.
### Note: y_t depends on s_{t-1} and eps_t shares between both equations
### (correlated state/observation noise).
###
### For historical decomposition, all endo vars are recovered via:
###   yall_t = ghx * s_{t-1} + ghu * eps_t
### ===========================================================================

# ---------------------------------------------------------------------------
#' Build compact DSGE state-space from model and decision rules
#'
#' Extracts T, R, Z, D matrices from ghx/ghu using lead_lag_incidence
#' to determine the correct ghx column ordering.
#'
#' State-space form:
#'   \eqn{s_t   = T s_{t-1} + R \varepsilon_t}  (n_state x 1)
#'   \eqn{yall_t = ghx \cdot s_{t-1} + ghu \cdot \varepsilon_t}  (n_endo x 1, all vars)
#'   \eqn{y_t   = Z s_{t-1} + D \varepsilon_t}  (n_obs x 1, observables only)
#'
#' @param m         Parsed model object (from parse_mod)
#' @param dr        Decision rules (from stoch_simul()$dr)
#' @param obs_vars  Character vector of observable names
#' @param verbose   Print diagnostic info (default TRUE)
#' @param params    Named parameter vector used to evaluate the shock
#'   covariance \code{Sigma_e} from the \code{shocks;} block (default:
#'   \code{m$param_values}). Pass the draw-specific vector when \code{dr}
#'   was solved at non-default parameters.
#' @return List with T_mat, R_mat, Z_mat, D_mat, Sigma_e, ghx, ghu,
#'   indices, names, and the observation intercept \code{d =
#'   dr$ys[obs_vars]} (with the full steady state in \code{ys}). The
#'   intercept is what makes a state space self-contained: everything that
#'   consumes one takes observables in LEVELS and subtracts \code{d}.
#' @export
# ---------------------------------------------------------------------------
build_dsge_state_space <- function(m, dr, obs_vars, verbose = TRUE,
                                   params = m$param_values) {
  
  endo_names <- m$var_names
  n_endo     <- length(endo_names)
  n_state    <- ncol(dr$ghx)
  n_shock    <- ncol(dr$ghu)
  
  stopifnot(nrow(dr$ghx) == n_endo,
            nrow(dr$ghu) == n_endo)
  
  ## ---- Step 1: Determine ghx column ordering from lead_lag_incidence ----
  ## lli[1, j] > 0 means variable j appears at t-1 (is a state variable)
  ## The VALUE of lli[1, j] gives the Jacobian column index, which
  ## determines the ordering of ghx columns.
  
  lli <- m$lead_lag_incidence
  ghx_col_to_endo <- NULL
  
  if (!is.null(lli) && is.matrix(lli) && nrow(lli) >= 1) {
    lag_row <- lli[1, ]
    has_lag <- which(lag_row > 0)
    
    if (length(has_lag) == n_state) {
      ## Sort by Jacobian column index to get ghx column ordering
      ghx_col_to_endo <- has_lag[order(lag_row[has_lag])]
      
      if (verbose) {
        state_names_ordered <- endo_names[ghx_col_to_endo]
        cat(sprintf("  ghx column order (from lli): %s\n",
                    paste(state_names_ordered, collapse = ", ")))
      }
    } else {
      warning(sprintf(
        "lead_lag_incidence lag count (%d) != ghx columns (%d). Falling back.",
        length(has_lag), n_state
      ))
    }
  }
  
  ## ---- Fallback: try declaration order of state variables ----
  if (is.null(ghx_col_to_endo)) {
    vc <- m$variable_classification
    back_vars  <- character(0)
    mixed_vars <- character(0)
    
    if (is.list(vc) && !is.data.frame(vc)) {
      if (!is.null(vc$predetermined)) back_vars  <- vc$predetermined
      if (!is.null(vc$backward))      back_vars  <- c(back_vars, vc$backward)
      if (!is.null(vc$mixed))         mixed_vars <- vc$mixed
    }
    
    state_vars <- c(back_vars, mixed_vars)
    state_idx_unsorted <- match(state_vars, endo_names)
    
    ## Try declaration order first (most common for dynhr)
    decl_order <- sort(state_idx_unsorted)
    pred_mixed_order <- state_idx_unsorted   # [pred, mixed] Dynare convention
    
    ## Test both: whichever gives stable T eigenvalues is correct
    T_decl <- dr$ghx[decl_order, , drop = FALSE]
    T_pm   <- dr$ghx[pred_mixed_order, , drop = FALSE]
    
    eig_decl <- max(abs(eigen(T_decl, only.values = TRUE)$values))
    eig_pm   <- max(abs(eigen(T_pm,   only.values = TRUE)$values))
    
    if (eig_decl < 1.0) {
      ghx_col_to_endo <- decl_order
      if (verbose) cat("  ghx column order: declaration order (verified by eigenvalues)\n")
    } else if (eig_pm < 1.0) {
      ghx_col_to_endo <- pred_mixed_order
      if (verbose) cat("  ghx column order: [predetermined, mixed] (verified by eigenvalues)\n")
    } else {
      warning(sprintf(
        "Neither ordering gives stable T. max|eig|: decl=%.4f, pred_mixed=%.4f",
        eig_decl, eig_pm
      ))
      ghx_col_to_endo <- decl_order  # best guess
    }
  }
  
  ## ---- Step 2: Build compact state-space matrices ----
  ## T_mat: state transition (rows = state vars in ghx column order)
  ## R_mat: shock impact on states
  T_mat <- dr$ghx[ghx_col_to_endo, , drop = FALSE]   # n_state x n_state
  R_mat <- dr$ghu[ghx_col_to_endo, , drop = FALSE]    # n_state x n_shock
  
  ## Z_mat: observation equation (obs rows of ghx)
  ## D_mat: direct shock -> obs
  obs_idx <- match(obs_vars, endo_names)
  if (any(is.na(obs_idx))) {
    stop("Observables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))
  }
  Z_mat <- dr$ghx[obs_idx, , drop = FALSE]             # n_obs x n_state
  D_mat <- dr$ghu[obs_idx, , drop = FALSE]              # n_obs x n_shock
  
  ## ---- Step 3: Verify eigenvalue stability ----
  eig_mod <- abs(eigen(T_mat, only.values = TRUE)$values)
  state_names_ordered <- endo_names[ghx_col_to_endo]
  
  if (verbose) {
    cat(sprintf("  State-space: %d states, %d obs, %d shocks\n",
                n_state, length(obs_vars), n_shock))
    cat(sprintf("  State vars: %s\n", paste(state_names_ordered, collapse = ", ")))
    cat(sprintf("  T eigenvalues: [%.4f, %.4f]",
                min(eig_mod), max(eig_mod)))
    if (all(eig_mod < 1.0)) cat(" -- all stable [OK]\n")
    else cat(sprintf(" -- UNSTABLE (max=%.4f) [X]\n", max(eig_mod)))
  }
  
  if (any(eig_mod >= 1.0)) {
    warning(sprintf(
      "T matrix has unstable eigenvalues (max |lambda| = %.4f).",
      max(eig_mod)
    ))
  }
  
  ## Shock covariance from the shocks; block. ghx/ghu are unit-shock
  ## responses (Sigma_e is NOT baked into them), so downstream filters
  ## must use this as Q -- kalman_smoother() defaults to it.
  Sigma_e <- .get_shock_cov(m, m$varexo_names, params)

  ## Observation intercept. Every filtering and smoothing entry point takes
  ## observables in LEVELS and subtracts the model's own steady state; the
  ## state space therefore carries that steady state with it. It used to carry
  ## no `ys` at all, which is why kalman_smoother() -- reachable only through a
  ## pre-built state space -- silently REQUIRED data already in deviations
  ## while kalman_filter() took levels. On the bundled nk_demo (observable
  ## steady states 0.5, 2, 4) that cost 33233 log points with no warning.
  ys_obs <- NULL
  if (!is.null(dr$ys) && !is.null(names(dr$ys))) {
    v <- dr$ys[obs_vars]
    if (length(v) == length(obs_vars) && !anyNA(v)) ys_obs <- as.numeric(v)
  }

  structure(
    list(
      T_mat             = T_mat,
      R_mat             = R_mat,
      Z_mat             = Z_mat,
      D_mat             = D_mat,
      Sigma_e           = Sigma_e,
      d                 = ys_obs,               # obs intercept: ys[obs_vars]
      ys                = dr$ys,                # full steady state (all endo)
      ghx               = dr$ghx,
      ghu               = dr$ghu,
      ghx_col_to_endo   = ghx_col_to_endo,     # maps ghx column i -> endo index
      obs_idx           = obs_idx,
      state_names       = state_names_ordered,  # in ghx column order
      obs_names         = obs_vars,            # observable names (character)
      endo_names        = endo_names,
      shock_names       = m$varexo_names,
      n_state           = n_state,
      n_obs             = length(obs_vars),
      n_shock           = n_shock,
      n_endo            = n_endo,
      timing            = "lagged"              # Convention A: y_t = Z s_{t-1} + D eps_t
    ),
    class = "dsge_ss"
  )
}


#' Which observation components carry information at this period?
#'
#' Skip-aware Cholesky of the innovation covariance. The i-th Cholesky pivot
#' IS the conditional variance of observable i given components 1..i-1, which
#' is exactly the quantity the univariate (Koopman-Durbin) filter tests one
#' observable at a time in \code{.kf_univariate_loop_R()}
#' (R/kalman-filter.R): \code{if (F_star > kalman_tol) \{ update \} else
#' \{ skip \}}. Dropping a component here therefore reproduces the filter's
#' decision, and because a skipped component contributes no row to the factor
#' the later conditional variances are computed WITHOUT conditioning on it --
#' again matching the sequential filter.
#'
#' The kept subset's F is positive definite by construction: its Cholesky
#' pivots are precisely the retained conditional variances, all above `tol`.
#'
#' @param F_t Innovation covariance (symmetric), n x n.
#' @param tol Pivot floor. Absolute, matching the filter's `kalman_tol`
#'   (Dynare's convention), but raised to a relative floor when `F_t` is
#'   badly scaled -- a unit-root smoother initialises P at 1e6 * I, so an
#'   absolute-only threshold is meaningless there.
#' @return Logical vector: TRUE for components that carry information.
#' @noRd
.smoother_informative_obs <- function(F_t, tol = 1e-10) {
  n <- nrow(F_t)
  if (n == 0L) return(logical(0))
  dmax <- max(abs(diag(F_t)), 0, na.rm = TRUE)
  tol  <- max(tol, n * dmax * .Machine$double.eps)
  keep <- logical(n)
  L    <- matrix(0, n, n)
  for (i in seq_len(n)) {
    prev <- seq_len(i - 1L)
    piv  <- F_t[i, i] - if (i > 1L) sum(L[i, prev]^2) else 0
    if (!is.finite(piv) || piv <= tol) next     # zero/negative variance: skip
    keep[i] <- TRUE
    L[i, i] <- sqrt(piv)
    if (i < n) {
      for (j in (i + 1L):n) {
        cr <- F_t[j, i] - if (i > 1L) sum(L[j, prev] * L[i, prev]) else 0
        L[j, i] <- cr / L[i, i]
      }
    }
  }
  keep
}


# ---------------------------------------------------------------------------
#' Kalman smoother for DSGE state-space
#'
#' DSGE form (note timing: \eqn{y_t} depends on \eqn{s_{t-1}}):
#' \preformatted{
#'   s_t = T s_{t-1} + R eps_t
#'   y_t = Z s_{t-1} + D eps_t
#' }
#'
#' State and observation noise are correlated (shared eps_t):
#'   Var(R eps) = R Q R'
#'   Var(D eps) = D Q D'
#'   Cov(R eps, D eps) = R Q D'
#'
#' Forward pass: standard KF with correlated noise.
#' Backward pass: Durbin-Koopman (2012) fixed-interval smoother -- ONE
#' adjoint recursion \eqn{(r_t, N_t)} delivers the smoothed shocks, the
#' smoothed states and the smoothed state covariances.
#'
#' \strong{Why not Rauch-Tung-Striebel.}  RTS assumes \eqn{y_{t+1}} is
#' conditionally independent of \eqn{s_t} given \eqn{s_{t+1}}.  Under this
#' package's lag-1 observation timing that is FALSE: \eqn{y_{t+1} = Z s_t +
#' D \varepsilon_{t+1}} loads \eqn{s_t} directly, so the RTS state pass is
#' not exact here (measured up to 17\% of the state scale on a 2-shock
#' fixture, ~7e-5 relative on the nk_2obc never-binding fixture).  Writing
#' the observation on the LAGGED state instead -- \eqn{\alpha_t := s_{t-1}},
#' \eqn{\alpha_{t+1} = T \alpha_t + R \varepsilon_t}, \eqn{y_t = Z \alpha_t +
#' D \varepsilon_t} -- puts the model in exactly Durbin & Koopman's
#' correlated-noise form (DK 2012 \eqn{\S}4.5, \eqn{\S}6.2), whose backward
#' recursion IS exact:
#' \preformatted{
#'   r_{t-1} = Z_t' F_t^{-1} v_t + L_t' r_t,   L_t = T - K_t Z_t,  r_T = 0
#'   N_{t-1} = Z_t' F_t^{-1} Z_t + L_t' N_t L_t,                   N_T = 0
#'   s_{t-1|T} = s_{t-1|t-1} + P_{t-1|t-1} r_{t-1}
#'   V_{t-1|T} = P_{t-1|t-1} - P_{t-1|t-1} N_{t-1} P_{t-1|t-1}
#' }
#' (the DK prediction pair \eqn{(a_t, P_t)} for \eqn{\alpha_t} is this
#' filter's UPDATED pair \eqn{(s_{t-1|t-1}, P_{t-1|t-1})}; correlated noise
#' changes \eqn{K_t}, and hence \eqn{L_t}, but not the form of the two
#' recursions).  Because no observation loads \eqn{s_T}, the last smoothed
#' state and covariance are the last FILTERED ones,
#' \eqn{s_{T|T}} and \eqn{P_{T|T}}, and the recursion also yields the
#' pre-sample \eqn{s_{0|T}} (returned as \code{smoothed_initial}) for free.
#' The smoothed states produced this way satisfy the transition
#' \eqn{s_{t|T} = T s_{t-1|T} + R \varepsilon_{t|T}} identically (~1e-16),
#' and the smoothed shocks are unchanged by this route.
#'
#' \strong{Exact diffuse initialisation is NOT supported.}  There is no
#' \code{lik_init} argument: \eqn{P_{0|0}} is the unconditional (Lyapunov)
#' covariance, with a large-\eqn{\kappa} diagonal fallback when the Lyapunov
#' solve fails on a unit root.  The smoothed moments are then exact for that
#' proper large-variance prior, not for the exact-diffuse limit; use
#' \code{kalman_filter(lik_init = "diffuse")} when the diffuse likelihood is
#' what is needed.
#'
#' @param data   \code{T x n_obs} matrix of observables in \strong{levels},
#'   in the column order of \code{obs_vars}. The model's own steady state is
#'   subtracted, exactly as \code{\link{kalman_filter}()} does; pass
#'   \code{d = 0} for series already in deviations.
#' @param dr     Decision rules from \code{\link{solve_perturbation}()} --
#'   the same second argument \code{\link{kalman_filter}()} takes.
#' @param model  Parsed model from \code{\link{parse_mod}()}.
#' @param params Named parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param obs_vars Character vector naming the observables, in the column
#'   order of \code{data}.
#' @param me_variance Base measurement-error variance, added to the diagonal
#'   of the innovation covariance for every observable and period -- the same
#'   argument, with the same meaning, as in \code{\link{kalman_filter}()}.
#'   Default \code{0}. Set it to the value you filtered with, or the two will
#'   describe different noise models.
#' @param d      Observation-intercept override, length \code{n_obs}.
#'   \code{NULL} (default) means the model's steady state
#'   \code{dr$ys[obs_vars]}, as in \code{\link{kalman_filter}()}. Pass
#'   \code{d = 0} to smooth data that is already in deviations.
#' @param lik_init Initialisation of \eqn{P_{0|0}}: \code{"auto"} (default)
#'   uses the unconditional Lyapunov covariance and falls back to a
#'   large-diagonal prior on a unit root; \code{"stationary"} demands the
#'   Lyapunov solution and errors if it does not exist; \code{"kappa"} forces
#'   the large-diagonal prior. There is no exact-diffuse option -- use
#'   \code{\link{kalman_filter}(lik_init = "diffuse")} for that.
#' @param kalman_tol Conditional-variance floor below which an observation
#'   component is treated as carrying no information and dropped for that
#'   period, matching \code{\link{kalman_filter}}'s univariate fallback.
#' @param Q      n_shock x n_shock shock covariance. Default \code{NULL}:
#'   use \code{ss$Sigma_e} (the covariance from the \code{shocks;} block),
#'   which makes the forward-pass \code{loglik} identical to
#'   \code{kalman_filter()} on the same data. \code{ghu} holds unit-shock
#'   responses, so the identity is only correct when every shock has
#'   \code{stderr 1}; falls back to the identity (with a warning) only for
#'   hand-built state spaces lacking \code{Sigma_e}.
#' @param me_extra  \code{n_obs x T} matrix of per-period additive
#'   measurement-error variances (from a \code{filter_tunes} block), or
#'   \code{NULL} (no filter tunes).
#' @param shock_scale  \code{n_exo x T} matrix of per-period shock
#'   standard-deviation scale factors (from a \code{heteroskedastic_shocks}
#'   block), or \code{NULL} (constant shock variances).
#' @return List with \code{smoothed_states} (\eqn{T \times n_{state}}, row
#'   \eqn{t} is \eqn{s_{t|T}}), \code{smoothed_shocks}, \code{filtered_states},
#'   \code{filtered_cov} / \code{predicted_cov} / \code{smoothed_cov},
#'   \code{smoothed_initial} (\eqn{s_{0|T}}),
#'   \code{smoothed_initial_cov} (\eqn{V_{0|T}}) and \code{loglik}.
#' @seealso \code{\link{kalman_filter}} for the filtered pass and the exact
#'   diffuse likelihood, \code{\link{build_dsge_state_space}} for the state
#'   space this builds internally, \code{\link{smoother2histval}} to turn a
#'   smoother run into an initial history, and
#'   \code{\link{kf_innovation_diagnostics}} to check the innovations are
#'   white.
#' @examples
#' mod <- system.file("extdata/models/nk_demo.mod", package = "dynhr")
#' m   <- parse_mod(mod, verbose = FALSE)
#' cp  <- compile_model(m, verbose = FALSE)
#' ss  <- solve_steady(cp, m$param_values)
#' dr  <- solve_perturbation(m, cp, ss$values, m$param_values)
#' Y   <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                       package = "dynhr")))
#' obs <- c("ygr", "infl", "intr")
#'
#' ## Filter and smooth take the same arguments.
#' kf <- kalman_filter(Y, dr, m, m$param_values, obs_vars = obs)
#' sm <- kalman_smoother(Y, dr, m, m$param_values, obs_vars = obs)
#'
#' dim(sm$smoothed_states)      # T x n_state, row t is s_{t|T}
#' dim(sm$smoothed_shocks)      # T x n_shock
#'
#' ## With measurement error, pass the SAME me_variance to both.
#' kf2 <- kalman_filter(Y, dr, m, m$param_values, obs_vars = obs,
#'                      me_variance = 0.01)
#' sm2 <- kalman_smoother(Y, dr, m, m$param_values, obs_vars = obs,
#'                        me_variance = 0.01)
#' @export
# ---------------------------------------------------------------------------
kalman_smoother <- function(data, dr, model, params = NULL,
                            obs_vars, me_variance = 0,
                            d = NULL, Q = NULL, me_extra = NULL,
                            shock_scale = NULL,
                            lik_init = c("auto", "stationary", "kappa"),
                            kalman_tol = 1e-10) {
  lik_init <- match.arg(lik_init)

  ## ---- One shape, one data convention -----------------------------------
  ## Up to 0.9.3 this function took a pre-built `dsge_ss` as its second
  ## argument and, because such an object carried no steady state, observables
  ## in DEVIATIONS -- while every sibling entry point took (data, dr, model,
  ## params, obs_vars, ...) and LEVELS. Two entry points to the same recursion
  ## silently requiring different data is a defect, not an interface, so the
  ## state-space convention is now internal (.kalman_smoother_ss()) and the
  ## state space carries its own intercept.
  if (inherits(dr, "dsge_ss") || (is.list(dr) && !is.null(dr$T_mat)))
    stop("kalman_smoother: the second argument is a decision-rule object, ",
         "not a state space. Replace\n",
         "    ss <- build_dsge_state_space(model, dr, obs_vars)\n",
         "    kalman_smoother(data, ss, ...)\n",
         "with\n",
         "    kalman_smoother(data, dr, model, params, obs_vars = obs_vars, ...)\n",
         "NOTE the data convention also changed: this call takes observables ",
         "in LEVELS and subtracts the model's steady state itself, exactly as ",
         "kalman_filter() does. If your series are already in deviations, pass ",
         "`d = 0`.", call. = FALSE)
  if (missing(model) || is.null(model))
    stop("kalman_smoother: `model` is required -- call it the same way you ",
         "call kalman_filter(): kalman_smoother(data, dr, model, params, ",
         "obs_vars = ...).", call. = FALSE)
  if (missing(obs_vars) || is.null(obs_vars))
    stop("kalman_smoother: `obs_vars` is required (which observables the ",
         "columns of `data` correspond to).", call. = FALSE)
  if (is.null(params)) params <- model$param_values

  ss <- build_dsge_state_space(model, dr, obs_vars, verbose = FALSE,
                               params = params)
  .kalman_smoother_ss(data, ss, d = d, me_variance = me_variance, Q = Q,
                      me_extra = me_extra, shock_scale = shock_scale,
                      lik_init = lik_init, kalman_tol = kalman_tol)
}


#' Kalman smoother on a pre-built state space
#'
#' The recursion itself. \code{kalman_smoother()} is the public face of this:
#' it builds the state space from \code{(dr, model, params, obs_vars)} and
#' calls this. Kept as an internal entry point for the in-package consumers
#' that already hold a \code{dsge_ss} -- conditional-forecast.R,
#' diag-orchestrate.R, forecast-eval.R, shock-decomposition.R -- and for
#' hand-built state spaces, which have no model behind them.
#'
#' @param ss A \code{dsge_ss}. Its \code{d} field is the observation
#'   intercept and is subtracted from \code{data}, so \code{data} is in
#'   LEVELS. A hand-built list without a \code{d} field carries no steady
#'   state, so nothing is subtracted and its data is in deviations by
#'   construction -- there is no other reading available to it.
#' @param d Intercept override, length \code{n_obs}; \code{NULL} (default)
#'   means \code{ss$d}. Pass \code{d = 0} for data already in deviations.
#' @inheritParams kalman_smoother
#' @noRd
.kalman_smoother_ss <- function(data, ss, d = NULL, me_variance = 0,
                                Q = NULL, me_extra = NULL, shock_scale = NULL,
                                lik_init = c("auto", "stationary", "kappa"),
                                kalman_tol = 1e-10) {
  lik_init <- match.arg(lik_init)

  ## Observation intercept. kalman_filter() subtracts d = dr$ys[obs_vars] from
  ## the data; so does this, via the state space's own `d` field.
  d_obs <- if (!is.null(d)) rep_len(as.numeric(d), ncol(data))
           else if (!is.null(ss$d) && !anyNA(ss$d)) as.numeric(ss$d)
           else NULL
  if (!is.null(d_obs) && length(d_obs) != ncol(data))
    stop(sprintf(paste0("kalman_smoother: the observation intercept has %d ",
                        "entr%s but the data has %d column(s). A recycled ",
                        "intercept would demean the wrong observables."),
                 length(d_obs), if (length(d_obs) == 1L) "y" else "ies",
                 ncol(data)), call. = FALSE)
  if (!is.null(d_obs) && all(d_obs == 0)) d_obs <- NULL

  ## Counters for the singular-F diagnostic raised after the forward pass.
  n_sing_periods <- 0L
  n_sing_dropped <- 0L

  ## Convert current-state dsge_ss to lagged-state before extracting matrices.
  ## ss_convert_timing() is a no-op when ss$timing == "lagged".
  if (inherits(ss, "dsge_ss") && !is.null(ss$timing) && ss$timing != "lagged")
    ss <- ss_convert_timing(ss)

  TT      <- nrow(data)
  n_s     <- ss$n_state
  n_obs   <- ss$n_obs
  n_shk   <- ss$n_shock

  TT_mat  <- ss$T_mat
  R_mat   <- ss$R_mat
  Z_mat   <- ss$Z_mat
  D_mat   <- ss$D_mat

  if (is.null(Q)) {
    if (!is.null(ss$Sigma_e)) {
      Q <- ss$Sigma_e
    } else {
      warning("kalman_smoother: ss has no Sigma_e (hand-built list?); ",
              "using Q = identity, which assumes all shocks have stderr 1.",
              call. = FALSE)
      Q <- diag(n_shk)
    }
  }

  ## ---- me_extra validation ------------------------------------------------
  ## me_extra (n_obs x T) holds per-period per-observable extra ME variances
  ## (filter_tunes soft tunes: stderr^2 at tune periods, 0 elsewhere).
  if (!is.null(me_extra)) {
    if (!is.matrix(me_extra) || nrow(me_extra) != n_obs || ncol(me_extra) != TT)
      stop(sprintf(
        "kalman_smoother: me_extra must be n_obs x T (%d x %d); got %s.",
        n_obs, TT,
        if (is.matrix(me_extra)) paste0(nrow(me_extra), " x ", ncol(me_extra))
        else "non-matrix"), call. = FALSE)
  }
  has_me_extra <- !is.null(me_extra) && any(me_extra != 0)

  ## ---- Measurement error --------------------------------------------------
  ## `me_variance` is the BASE observation-noise variance, the same argument
  ## kalman_filter() takes; `me_extra` remains the optional per-period EXTRA
  ## on top of it. Until now the smoother had no base term at all, so a model
  ## filtered with me_variance > 0 could not be smoothed with the same noise
  ## model -- the two disagreed by construction.
  ##
  ## In this recursion the base term enters the innovation covariance ONLY,
  ## and that is the complete treatment: with y_t = Z s_{t-1} + D eps_t + u_t
  ## and u_t independent of the state and the shocks, u contributes to
  ## F = Z P Z' + D Q D' + H but leaves the cross-covariance
  ## Cov(s_t, y_t) = T P Z' + R Q D' untouched, and the update
  ## P_{t|t} = P_{t|t-1} - K F K' with the optimal gain already carries H
  ## through F. (kf_step() needs an explicit Joseph term because it is written
  ## in a form where the gain is not the optimal one for the noise-inclusive
  ## F; here it is.)
  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance < 0)
    stop("kalman_smoother: `me_variance` must be a non-negative finite ",
         "scalar.", call. = FALSE)
  has_me_var <- me_variance > 0

  ## ---- shock_scale validation ----------------------------------------------
  ## shock_scale (n_shk x T) holds per-period shock std scale factors.
  if (!is.null(shock_scale)) {
    if (!is.matrix(shock_scale) || nrow(shock_scale) != n_shk || ncol(shock_scale) != TT)
      stop(sprintf(
        "kalman_smoother: shock_scale must be n_shk x T (%d x %d); got %s.",
        n_shk, TT,
        if (is.matrix(shock_scale)) paste0(nrow(shock_scale), " x ", ncol(shock_scale))
        else "non-matrix"), call. = FALSE)
  }
  has_shock_scale <- !is.null(shock_scale) && !all(shock_scale == 1)

  ## Pre-compute noise covariances (baseline, used for P0 and when not scaling)
  RQR <- R_mat %*% Q %*% t(R_mat)       # state noise cov
  DQD <- D_mat %*% Q %*% t(D_mat)       # obs noise cov
  RQD <- R_mat %*% Q %*% t(D_mat)       # cross-covariance

  ## ---- Initialisation: unconditional state covariance (Lyapunov) ----
  ## Use the canonical solve_lyapunov() from stochsimul-monolith.R;
  ## fall back to large diagonal for near-unit-root / nonstationary models
  ## (solve_lyapunov returns a NaN matrix when the doubling algorithm diverges
  ## and the vec-Lyapunov system is singular).
  P_ss <- if (identical(lik_init, "kappa")) {
    ## Forced diffuse-style prior. kalman_filter() REFUSES lik_init = "auto"
    ## together with shock_scale on a nonstationary model and tells the caller
    ## to pass "kappa" or "stationary" explicitly; before this argument existed
    ## there was no way to ask the smoother for the same thing, so the two
    ## could not be made comparable even in principle.
    matrix(NA_real_, n_s, n_s)
  } else {
    solve_lyapunov(TT_mat, RQR)
  }
  if (identical(lik_init, "stationary") && anyNA(P_ss))
    stop("kalman_smoother: lik_init = \"stationary\" was requested but the ",
         "Lyapunov solve returned NaN -- TT has unit-root eigenvalues, so the ",
         "unconditional state covariance does not exist. Use lik_init = ",
         "\"kappa\" (or \"auto\") for a nonstationary model.", call. = FALSE)
  if (anyNA(P_ss)) {
    ## Unit roots detected: use a diffuse (large-diagonal) prior so the
    ## smoother does not crash.  The smoothed states will be valid but the
    ## loglik has a kappa-dependent additive offset (not suitable for
    ## cross-method comparison; use kalman_filter(lik_init="diffuse") for
    ## exact diffuse likelihood evaluation).
    if (!identical(lik_init, "kappa"))
    warning("kalman_smoother: unit root(s) detected in TT -- ",
            "solve_lyapunov() returned NaN. ",
            "Falling back to diffuse prior P0 = ", .DIFFUSE_SCALE,
            " * I(", n_s, "). ",
            "Smoothed states are valid; loglik has a kappa-dependent offset. ",
            "Use kalman_filter(lik_init=\"diffuse\") for exact diffuse loglik.",
            call. = FALSE)
    P_ss <- .DIFFUSE_SCALE * diag(n_s)
  }

  ## ---- Forward pass (Kalman filter) ----
  ## Per-period NA observations are handled by dropping the NA rows from
  ## Z_mat / DQD / RQD / the innovation vector.  An all-NA period is treated
  ## as predict-only (no update).  This mirrors the standard filter's NA
  ## branch in R/kalman-filter.R:1193-1218 so the forward quantities fed into
  ## the DK backward pass are always consistent.
  ##
  ## The DK backward pass requires: v_t (innovations), F_inv_t (inverse
  ## innovation covariance), K_t (Kalman gain), and Zt/Dt (obs/shock matrices
  ## subsetting to non-NA rows), plus the UPDATED pair (s_{t|t}, P_{t|t}) --
  ## which is DK's PREDICTION pair for alpha_{t+1} = s_t.  These are stored in
  ## lists / arrays indexed by period.
  s_filt <- matrix(0, TT, n_s)         # s_{t|t}
  ## s_pred is not consumed by the backward pass (only P_pred is, via the
  ## returned predicted_cov); kept as the mean counterpart of P_pred.
  s_pred <- matrix(0, TT, n_s)         # s_{t|t-1}
  P_filt <- array(0, dim = c(n_s, n_s, TT))
  P_pred <- array(0, dim = c(n_s, n_s, TT))
  ## Number of non-NA observations per period (for the loglik correction term).
  n_obs_t_vec <- integer(TT)

  ## Storage for DK disturbance smoother backward pass.
  dk_v     <- vector("list", TT)   # v_t (n_ok x 1)
  dk_Finv  <- vector("list", TT)   # F_t^{-1} (n_ok x n_ok)
  dk_K     <- vector("list", TT)   # K_t (n_s x n_ok)
  dk_Z     <- vector("list", TT)   # Zt  (n_ok x n_s)
  dk_D     <- vector("list", TT)   # Dt  (n_ok x n_shk)
  dk_Q     <- vector("list", TT)   # Q_t (n_shk x n_shk) -- needed when shock_scale active

  s_tt <- rep(0, n_s)
  P_tt <- P_ss
  loglik <- 0

  for (t in seq_len(TT)) {
    ## Per-period scaled shock covariance (when shock_scale is active).
    Q_t   <- if (has_shock_scale) {
      sc_t <- shock_scale[, t]
      Q * outer(sc_t, sc_t)
    } else Q
    RQR_t <- if (has_shock_scale) R_mat %*% Q_t %*% t(R_mat) else RQR
    DQD_t <- if (has_shock_scale) D_mat %*% Q_t %*% t(D_mat) else DQD
    RQD_t <- if (has_shock_scale) R_mat %*% Q_t %*% t(D_mat) else RQD

    ## ---- Predict ----
    s_tp <- as.numeric(TT_mat %*% s_tt)
    P_tp <- TT_mat %*% P_tt %*% t(TT_mat) + RQR_t

    s_pred[t, ]   <- s_tp
    P_pred[, , t] <- P_tp

    ## ---- Innovation (with per-period NA handling) ----
    y_pred <- as.numeric(Z_mat %*% s_tt)
    if (!is.null(d_obs)) y_pred <- y_pred + d_obs
    v_t    <- data[t, ] - y_pred

    obs_ok <- which(!is.na(v_t))
    n_ok   <- length(obs_ok)
    n_obs_t_vec[t] <- n_ok

    if (n_ok == 0L) {
      ## All observables missing: predict-only, no update.
      s_tt <- s_tp
      P_tt <- P_tp
      s_filt[t, ]   <- s_tt
      P_filt[, , t] <- P_tt
      ## DK: no observation this period -- v, Finv, K left NULL; Q stored.
      dk_Q[[t]] <- Q_t
      next
    }

    ## Subset to non-NA observables.
    if (n_ok < n_obs) {
      v_t   <- v_t[obs_ok]
      Zt    <- Z_mat[obs_ok, , drop = FALSE]
      Dt    <- D_mat[obs_ok, , drop = FALSE]
      DQDt  <- tcrossprod(Dt %*% Q_t, Dt)
      RQDt  <- RQD_t[, obs_ok, drop = FALSE]
    } else {
      Zt   <- Z_mat
      Dt   <- D_mat
      DQDt <- DQD_t
      RQDt <- RQD_t
    }

    ## Innovation covariance: F_t = Z_t P_{t-1|t-1} Z_t' + DQD_t [+ me_extra_t]
    F_t  <- Zt %*% P_tt %*% t(Zt) + DQDt
    if (has_me_var)   diag(F_t) <- diag(F_t) + me_variance
    if (has_me_extra) diag(F_t) <- diag(F_t) + me_extra[obs_ok, t]
    F_t  <- (F_t + t(F_t)) * 0.5

    ## base chol() THROWS on a non-PD matrix (it never returns NULL), so the
    ## singular-F branch below must catch the error to run at all.
    F_ch <- tryCatch(chol(F_t), error = function(e) NULL)

    if (is.null(F_ch)) {
      ## ---- SINGULAR / NON-PD INNOVATION COVARIANCE ------------------------
      ## This used to add JITTER on an absolute ladder (1e-8 ... 1e-2, then an
      ## UNGUARDED chol(F_t + 0.1 * I)). Two things were wrong with that.
      ##
      ## (1) The jitter was absolute while F_t is not O(1). On a unit-root
      ##     model this smoother initialises P at .DIFFUSE_SCALE * I = 1e6 * I,
      ##     and `shock_scale` multiplies Q on top, so F_t can be many orders
      ##     of magnitude larger -- a 0.1 nudge means nothing. Where it did
      ##     "work" it silently corrupted the answer instead: switching one
      ##     shock off via shock_scale (the u_k = 0 hard-tune idiom) moved the
      ##     loglik from -63.7 to -8.5e+09 while kalman_filter returned -63.7.
      ## (2) The last rung had no tryCatch, so when +0.1*I was still not PD it
      ##     THREW -- the filter succeeding where the smoother rejects.
      ##
      ## A singular F means some observation component has zero forecast
      ## variance: it is predictable exactly and carries no new information.
      ## kalman_filter already handles this the right way, by falling back to
      ## the univariate (Koopman-Durbin) filter, which processes observables
      ## one at a time and SKIPS a component whose conditional variance is
      ## below `kalman_tol` (R/kalman-filter.R, .kf_univariate_loop_R:185).
      ##
      ## The smoother now makes the same decision. A zero-variance component
      ## is dropped for this period -- which is exactly how the smoother
      ## already treats a MISSING observable -- and the multivariate update
      ## proceeds on the informative subset, whose F is positive definite by
      ## construction (its Cholesky pivots are the retained conditional
      ## variances). `n_ok` shrinks with it, so the log(2*pi) term in the
      ## likelihood adjusts automatically, matching the filter's convention
      ## that a skipped component contributes nothing.
      keep <- .smoother_informative_obs(F_t, kalman_tol)
      n_sing_periods <- n_sing_periods + 1L
      n_sing_dropped <- n_sing_dropped + sum(!keep)

      if (!any(keep)) {
        ## No component carries information: predict-only, exactly as for an
        ## all-missing period.
        s_tt <- s_tp
        P_tt <- P_tp
        s_filt[t, ]   <- s_tt
        P_filt[, , t] <- P_tt
        n_obs_t_vec[t] <- 0L
        dk_Q[[t]] <- Q_t
        next
      }

      obs_ok <- obs_ok[keep]
      n_ok   <- length(obs_ok)
      n_obs_t_vec[t] <- n_ok
      v_t    <- v_t[keep]
      Zt     <- Z_mat[obs_ok, , drop = FALSE]
      Dt     <- D_mat[obs_ok, , drop = FALSE]
      DQDt   <- tcrossprod(Dt %*% Q_t, Dt)
      RQDt   <- RQD_t[, obs_ok, drop = FALSE]
      F_t    <- Zt %*% P_tt %*% t(Zt) + DQDt
      if (has_me_var)   diag(F_t) <- diag(F_t) + me_variance
      if (has_me_extra) diag(F_t) <- diag(F_t) + me_extra[obs_ok, t]
      F_t    <- (F_t + t(F_t)) * 0.5
      F_ch   <- tryCatch(chol(F_t), error = function(e) NULL)
      if (is.null(F_ch))
        stop("kalman_smoother: the innovation covariance at period ", t,
             " is not positive definite even after dropping every ",
             "zero-variance observation component. This should not happen -- ",
             "the retained components' conditional variances are the Cholesky ",
             "pivots and were all above the tolerance. Please report it with ",
             "a reproducible model.", call. = FALSE)
    }

    F_inv     <- chol2inv(F_ch)
    log_det_F <- 2 * sum(log(diag(F_ch)))

    ## Kalman gain: K = (T P_{t-1|t-1} Z_t' + RQD_t) F_t^{-1}
    K_t <- (TT_mat %*% P_tt %*% t(Zt) + RQDt) %*% F_inv

    ## Store forward-pass quantities for DK disturbance smoother.
    dk_v[[t]]    <- v_t
    dk_Finv[[t]] <- F_inv
    dk_K[[t]]    <- K_t
    dk_Z[[t]]    <- Zt
    dk_D[[t]]    <- Dt
    dk_Q[[t]]    <- Q_t

    ## Updated state: s_{t|t} = s_{t|t-1} + K (y_t - Z_t s_{t-1|t-1})
    s_tt <- s_tp + as.numeric(K_t %*% v_t)
    P_tt <- P_tp - K_t %*% F_t %*% t(K_t)
    P_tt <- 0.5 * (P_tt + t(P_tt))

    s_filt[t, ]   <- s_tt
    P_filt[, , t] <- P_tt

    ## Log-likelihood contribution (adjusts constant for n_ok != n_obs).
    loglik <- loglik - 0.5 * (n_ok * log(2 * pi) + log_det_F +
                                as.numeric(t(v_t) %*% F_inv %*% v_t))
  }

  ## ---- Singular-F diagnostic (once, after the forward pass) -------------
  ## Loud by construction: dropping a component is the CORRECT treatment (it
  ## carries no information), but it also means the model implies some
  ## observable is predictable exactly -- usually a switched-off shock via
  ## `shock_scale`, a hard `filter_tunes` tune, or stochastic singularity --
  ## and the caller should know rather than discover it in a loglik that is
  ## not comparable with a run where nothing was dropped.
  if (n_sing_periods > 0L)
    warning(sprintf(paste0(
      "kalman_smoother: the innovation covariance was singular in %d of %d ",
      "period(s); %d zero-variance observation component(s) were dropped in ",
      "total. Those components are predictable exactly and carry no ",
      "information, so they are skipped -- the same decision kalman_filter ",
      "makes when it falls back to the univariate filter. Smoothed states ",
      "and shocks remain valid; the log-likelihood is conditioned on fewer ",
      "components and is NOT comparable with a run in which none were ",
      "dropped."), n_sing_periods, TT, n_sing_dropped), call. = FALSE)

  ## ---- DK backward pass: smoothed shocks, states and state covariances ----
  ## Uses the Durbin-Koopman (2012) adjoint backward recursion to compute
  ## eps_{t|T} directly from the backward adjoint, matching Dynare's
  ## calib_smoother convention at ALL t including t=1.
  ##
  ## State-space (lagged form, DK §4.4):
  ##   s_t   = T s_{t-1} + R eps_t       -> alpha_t = s_{t-1}; eta_t = R eps_t
  ##   y_t   = Z s_{t-1} + D eps_t       -> obs noise = D eps_t; cross-cov = RQD'
  ##
  ## DK adjoint backward recursion (initialised r_T = 0):
  ##   r_{t-1} = Zt' Finv_t v_t + (T - K_t Zt)' r_t
  ##   (for periods with no obs: r_{t-1} = T' r_t)
  ##
  ## Smoothed structural shock (DK §4.4 eqns 4.44-4.47, correlated-noise form):
  ##
  ##   eps_{t|T} = Q_t (R' r_t + Dt' u_t)
  ##
  ## where r_t is the adjoint BEFORE the backward step at period t (i.e., the
  ## value coming in from period t+1), and u_t = Finv_t v_t - K_t' r_{t-1}.
  ##
  ## NOTE on timing: in DK's alpha_t = s_{t-1} indexing, the shock eps_t
  ## drives alpha_{t+1} = T alpha_t + R eps_t and enters y_t = Z alpha_t + D eps_t.
  ## The smoothed shock is hat_eta_t = Q_eta r_t (DK 4.44) with Q_eta = RQR'.
  ## Since hat_eta_t = R eps_{t|T}, this gives eps_{t|T} = Q R' r_t (D=0 case).
  ## For D != 0, the correlated-noise correction adds Q Dt' u_t (DK §4.4).
  ##
  ## For t>=2 this is numerically identical to the previous RTS+pseudo-inverse
  ## approach (both are minimum-MSE). At t=1 the DK formula uses r_1 (the
  ## backward adjoint from t=1 forward step) to get eps_{1|T} = Q (R' r_1 + D' u_1),
  ## whereas the previous code used J_0 to compute s_{0|T} residually; these
  ## differ only in the initial-condition treatment (both internally consistent).
  ## This matches Dynare's etahat(:,1) = Q*R'*r(:,1) convention.
  ##
  ## STATES AND COVARIANCES ride the SAME recursion (this replaced an RTS pass
  ## that was not exact under the lag-1 timing -- see the roxygen block).  With
  ## alpha_t := s_{t-1} the DK prediction pair for alpha_t is this filter's
  ## UPDATED pair (s_{t-1|t-1}, P_{t-1|t-1}), so after the backward step at
  ## period t (when `r_t`/`N_t` hold r_{t-1}/N_{t-1}):
  ##   s_{t-1|T} = s_{t-1|t-1} + P_{t-1|t-1} r_{t-1}
  ##   V_{t-1|T} = P_{t-1|t-1} - P_{t-1|t-1} N_{t-1} P_{t-1|t-1}
  ## with (s_{0|0}, P_{0|0}) = (0, P_ss) at t = 1.  N_t is PSD by construction,
  ## so diag(V_{t|T}) <= diag(P_{t|t}) still holds.
  ## No observation loads s_T, hence s_{T|T} / P_{T|T} ARE the smoothed pair at
  ## t = T; the DK pass supplies t = 0 .. T-1 (row k of s_lag = s_{k-1|T}).
  eps_smooth <- matrix(0, TT, n_shk)
  s_lag      <- matrix(0, TT, n_s)                 # row k = s_{k-1|T}
  V_lag      <- array(0, dim = c(n_s, n_s, TT))    # slice k = V_{k-1|T}
  r_t <- rep(0, n_s)               # r_{T} = 0 (terminal condition)
  N_t <- matrix(0, n_s, n_s)       # N_{T} = 0 (terminal condition)

  for (t in TT:1L) {
    Q_t <- dk_Q[[t]]
    ## eps_{t|T} = Q_t (R' r_t + Dt' u_t) -- uses r_t (BEFORE backward step at t)
    if (is.null(dk_v[[t]])) {
      ## All observables missing at period t.
      ## u_t doesn't exist; eps_{t|T} = Q_t R' r_t (only state adjoint, no D term)
      eps_smooth[t, ] <- as.numeric(Q_t %*% (t(R_mat) %*% r_t))
      ## Backward step: r_{t-1} = T' r_t, N_{t-1} = T' N_t T  (L_t = T, K_t = 0)
      r_t <- as.numeric(t(TT_mat) %*% r_t)
      N_t <- crossprod(TT_mat, N_t) %*% TT_mat
    } else {
      v_t    <- dk_v[[t]]
      Finv   <- dk_Finv[[t]]
      K_t_dk <- dk_K[[t]]
      Zt     <- dk_Z[[t]]
      Dt     <- dk_D[[t]]
      ## u_t = Finv v_t - K_t' r_t  (uses r_t BEFORE the backward step)
      ## NOTE: for dynhr's lagged-state form y_t = Z s_{t-1} + D eps_t, the
      ## correct u_t uses r_t (not r_{t-1} as in the standard DK §4.3 form for
      ## y_t = Z alpha_t). See derivation: for t=T, r_T=0 gives u_T = F_T^{-1}v_T,
      ## and eps_{T|T} = Q D' F_T^{-1} v_T -- the terminal shock is identified
      ## purely from the observation (correct when D != 0).
      u_t     <- as.numeric(Finv %*% v_t - t(K_t_dk) %*% r_t)
      ## eps_{t|T} = Q_t (R' r_t + Dt' u_t)  -- both use r_t (BEFORE backward step)
      eps_smooth[t, ] <- as.numeric(Q_t %*% (t(R_mat) %*% r_t + t(Dt) %*% u_t))
      ## Backward step: r_{t-1} = Zt' Finv v_t + (T - K_t Zt)' r_t
      Lt    <- TT_mat - K_t_dk %*% Zt           # L_t = T - K_t Z_t  (n_s x n_s)
      r_t   <- as.numeric(t(Zt) %*% (Finv %*% v_t) + t(Lt) %*% r_t)
      ## N_{t-1} = Zt' Finv Zt + L_t' N_t L_t (same L_t, correlated noise and
      ## all -- correlated noise changes K_t, not the form of the recursion).
      N_t   <- crossprod(Zt, Finv %*% Zt) + crossprod(Lt, N_t) %*% Lt
    }
    N_t <- 0.5 * (N_t + t(N_t))

    ## Smoothed lagged state / covariance: r_t and N_t now hold r_{t-1}, N_{t-1}.
    ## (drop = FALSE: n_s == 1 would otherwise collapse the slice to a scalar)
    s_in <- if (t == 1L) rep(0, n_s) else s_filt[t - 1L, ]
    P_in <- if (t == 1L) P_ss        else
      matrix(P_filt[, , t - 1L], n_s, n_s)
    s_lag[t, ]   <- s_in + as.numeric(P_in %*% r_t)
    V_t          <- P_in - P_in %*% N_t %*% P_in
    V_lag[, , t] <- 0.5 * (V_t + t(V_t))
  }

  ## Re-index to the s_{t|T} convention: row t of s_smooth is s_{t|T}.
  s_smooth <- matrix(0, TT, n_s)
  V_smooth <- array(0, dim = c(n_s, n_s, TT))
  if (TT >= 2L) {
    s_smooth[seq_len(TT - 1L), ]    <- s_lag[-1L, , drop = FALSE]
    V_smooth[, , seq_len(TT - 1L)]  <- V_lag[, , -1L, drop = FALSE]
  }
  s_smooth[TT, ]   <- s_filt[TT, ]     # nothing observes s_T => filtered == smoothed
  V_smooth[, , TT] <- P_filt[, , TT]
  s0_smooth <- s_lag[1L, ]                          # s_{0|T}
  V0_smooth <- matrix(V_lag[, , 1L], n_s, n_s)      # V_{0|T}

  colnames(s_smooth)   <- ss$state_names
  colnames(s_filt)     <- ss$state_names
  colnames(eps_smooth) <- ss$shock_names
  names(s0_smooth)     <- ss$state_names
  dimnames(V0_smooth)  <- list(ss$state_names, ss$state_names)

  dimnames(P_filt)   <- list(ss$state_names, ss$state_names, NULL)
  dimnames(P_pred)   <- list(ss$state_names, ss$state_names, NULL)
  dimnames(V_smooth) <- list(ss$state_names, ss$state_names, NULL)

  list(
    smoothed_states = s_smooth,
    smoothed_shocks = eps_smooth,
    filtered_states = s_filt,
    ## Per-period state covariances (n_state x n_state x T):
    filtered_cov    = P_filt,            # P_{t|t}
    predicted_cov   = P_pred,            # P_{t|t-1}
    smoothed_cov    = V_smooth,          # P_{t|T} (DK)
    P_filt_last     = P_filt[, , TT],    # P_{T|T} (kept for back-compat)
    ## Pre-sample smoothed moments, free from the same backward recursion.
    smoothed_initial     = s0_smooth,    # s_{0|T}
    smoothed_initial_cov = V0_smooth,    # V_{0|T}
    loglik          = loglik
  )
}


# ---------------------------------------------------------------------------
#' Historical decomposition: per-shock contributions to ALL endo variables
#'
#' Uses the full ghx/ghu (not just the compact state-space) to recover
#' contributions to all n_endo variables, not just states.
#'
#' For each shock j:
#'   \deqn{s_t^{(j)} = T s_{t-1}^{(j)} + R[:,j] \epsilon_{j,t}}
#'   \deqn{y_t^{(j)} = ghx \cdot s_{t-1}^{(j)} + ghu[:,j] \epsilon_{j,t}}
#'
#' The per-shock recursions all start at \eqn{s_0^{(j)} = 0}, so on their own
#' they reconstruct the smoothed path only when the smoothed initial state is
#' the steady state. The \code{"initial"} component closes that gap: it is the
#' shock-free trajectory started at \eqn{s_{0|T}},
#'
#'   \deqn{s_t^{(0)} = T s_{t-1}^{(0)}, \quad s_0^{(0)} = s_{0|T}, \qquad
#'         y_t^{(0)} = ghx \cdot s_{t-1}^{(0)},}
#'
#' which makes the adding-up
#' \eqn{\sum_j y_t^{(j)} + y_t^{(0)} = ghx s_{t-1|T} + ghu \varepsilon_{t|T}}
#' exact at every \eqn{t} (to smoother round-off, ~1e-15), not just
#' asymptotically once the initial condition has decayed.
#'
#' @param smoothed_shocks  \code{T x n_shock} matrix of smoothed structural
#'   shocks, OR the whole \code{\link{kalman_smoother}} result list -- in
#'   which case the smoothed states and the smoothed initial state
#'   \eqn{s_{0|T}} are taken from it too (this is the form that makes the
#'   adding-up exact without any further arguments).
#' @param ss               State-space list from
#'   \code{\link{build_dsge_state_space}()}. Contributions are in
#'   DEVIATIONS from the steady state (add \code{ss$ys[v]} to read variable
#'   \code{v} back in levels); only the data going INTO the smoother is in
#'   levels.
#' @param s0               Numeric length-\code{n_state} initial state (in
#'   \code{ss$state_names} order). Default \code{NULL}: taken from the
#'   smoother result when one was passed, otherwise zero (the pre-0.9.2
#'   behaviour -- the shock columns are unaffected either way).
#' @param smoothed_states  Optional \code{T x n_state} matrix of
#'   \eqn{s_{t|T}}; when supplied (or carried by the smoother result) the
#'   returned object also gains \code{$smoothed} (the smoothed series
#'   reconstructed directly from the states) and
#'   \code{$adding_up_residual}.
#' @param shock_groups     Optional named list grouping shocks, e.g.
#'   \code{list(supply = c("e_a", "e_z"), demand = "e_g")} (Dynare
#'   \code{shock_groups} semantics). Grouped output replaces the per-shock
#'   columns by group sums; shocks in no group land in an \code{"other"}
#'   column. Defaults to \code{model$shock_groups} when \code{model} is given.
#' @param model            Optional parsed model, used only as the fallback
#'   source of \code{shock_groups} (the explicit argument wins).
#' @return List of class \code{dynhr_shock_decomposition} with
#'   \code{$contributions} (named list of \code{T x n_endo} matrices: one per
#'   shock or group, plus \code{"initial"}), \code{$total} (\code{T x n_endo},
#'   the sum of all components = the smoothed series in deviations from
#'   steady state), \code{$initial}, \code{$s0} and \code{$components}.
#' @export
# ---------------------------------------------------------------------------
historical_decomposition <- function(smoothed_shocks, ss, s0 = NULL,
                                     smoothed_states = NULL,
                                     shock_groups = NULL, model = NULL) {

  ## Accept a whole kalman_smoother() result: shocks, states and s_{0|T}.
  if (is.list(smoothed_shocks) && !is.matrix(smoothed_shocks) &&
      !is.null(smoothed_shocks$smoothed_shocks)) {
    sm <- smoothed_shocks
    if (is.null(smoothed_states)) smoothed_states <- sm$smoothed_states
    if (is.null(s0)) s0 <- smoothed_initial_state(sm, ss)
    smoothed_shocks <- sm$smoothed_shocks
  }

  TT    <- nrow(smoothed_shocks)
  n_shk <- ss$n_shock
  n_s   <- ss$n_state
  n_end <- ss$n_endo

  TT_mat <- ss$T_mat       # n_state x n_state
  R_mat  <- ss$R_mat        # n_state x n_shock
  ghx    <- ss$ghx           # n_endo x n_state
  ghu    <- ss$ghu           # n_endo x n_shock

  if (is.null(s0)) s0 <- rep(0, n_s)
  s0 <- as.numeric(s0)
  if (length(s0) != n_s)
    stop(sprintf("historical_decomposition: s0 has length %d, expected %d.",
                 length(s0), n_s), call. = FALSE)

  contributions <- setNames(
    lapply(seq_len(n_shk), function(j) matrix(0, TT, n_end)),
    ss$shock_names
  )

  for (j in seq_len(n_shk)) {
    s_j <- rep(0, n_s)        # state attributable to shock j
    r_j <- R_mat[, j]          # state impact column
    g_j <- ghu[, j]            # full endo impact column

    for (t in seq_len(TT)) {
      eps_jt <- smoothed_shocks[t, j]
      if (is.na(eps_jt)) eps_jt <- 0

      ## All endo vars at t from shock j:
      ## y_t^(j) = ghx * s_{t-1}^(j) + ghu[:,j] * eps_jt
      contributions[[j]][t, ] <- as.numeric(ghx %*% s_j + g_j * eps_jt)

      ## State transition for shock j:
      s_j <- as.numeric(TT_mat %*% s_j + r_j * eps_jt)
    }
    colnames(contributions[[j]]) <- ss$endo_names
  }

  ## ---- Initial-condition component: zero shocks, state started at s0 ------
  initial <- matrix(0, TT, n_end)
  s_i     <- s0
  for (t in seq_len(TT)) {
    initial[t, ] <- as.numeric(ghx %*% s_i)
    s_i          <- as.numeric(TT_mat %*% s_i)
  }
  colnames(initial) <- ss$endo_names

  ## ---- Optional shock groups (applied to the shock columns only) ----------
  groups <- .resolve_shock_groups(shock_groups, model, ss$shock_names)
  if (!is.null(groups)) contributions <- .group_contributions(contributions, groups)

  contributions$initial <- initial

  total <- Reduce(`+`, contributions)
  colnames(total) <- ss$endo_names

  out <- list(contributions = contributions,
              total         = total,
              initial       = initial,
              s0            = s0,
              components    = names(contributions),
              shock_groups  = groups,
              orientation   = "time_endo")

  ## ---- Adding-up check against the directly-reconstructed smoothed path ---
  if (!is.null(smoothed_states)) {
    smoothed <- matrix(0, TT, n_end)
    s_prev   <- s0
    for (t in seq_len(TT)) {
      eps_t <- smoothed_shocks[t, ]
      eps_t[is.na(eps_t)] <- 0
      smoothed[t, ] <- as.numeric(ghx %*% s_prev + ghu %*% eps_t)
      s_prev        <- smoothed_states[t, ]
    }
    colnames(smoothed)     <- ss$endo_names
    out$smoothed           <- smoothed
    out$adding_up_residual <- max(abs(total - smoothed))
  }

  structure(out, class = c("dynhr_shock_decomposition", "list"))
}
