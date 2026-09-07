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
#' @param Sigma_e  Optional \code{n_exo x n_exo} shock covariance, overriding
#'   the one implied by the model's \code{shocks;} block at \code{params}.
#'   \code{NULL} (default) derives it, which is what nearly every caller
#'   wants.
#'
#'   This is the deliberate injection point. Note that \code{dr$Sigma_e} is
#'   NOT consulted here: \code{params} stays authoritative so that one solved
#'   decision rule can be reused while the likelihood is evaluated at many
#'   parameter values -- the pattern estimation depends on. A \code{dr} whose
#'   \code{Sigma_e} disagrees with \code{params} raises a warning saying so,
#'   because \code{\link{compute_irfs}} and \code{\link{compute_moments}}
#'   DO honour that field and the asymmetry is easy to trip over.
#' @return List with T_mat, R_mat, Z_mat, D_mat, Sigma_e, ghx, ghu,
#'   indices, names, and the observation intercept \code{d =
#'   dr$ys[obs_vars]} (with the full steady state in \code{ys}). The
#'   intercept is what makes a state space self-contained: everything that
#'   consumes one takes observables in LEVELS and subtracts \code{d}.
#' @export
# ---------------------------------------------------------------------------
build_dsge_state_space <- function(m, dr, obs_vars, verbose = TRUE,
                                   params = m$param_values, Sigma_e = NULL) {
  
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
  ## `params` is authoritative (see .kf_report_sigma_e_conflict for why it is
  ## NOT dr$Sigma_e); an explicit `Sigma_e` argument overrides both, and a
  ## disagreeing dr$Sigma_e is reported rather than silently ignored.
  if (is.null(Sigma_e)) {
    Sigma_e <- .get_shock_cov(m, m$varexo_names, params)
    .kf_report_sigma_e_conflict(m, dr, m$varexo_names, params, Sigma_e)
  } else {
    Sigma_e <- as.matrix(Sigma_e)
  }

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
      ## Which shocks carry no variance at all. A per-shock `stderr 0` is
      ## legitimate -- it is what makes a deterministic known_shocks or
      ## shock_means injection meaningful -- so this is RECORDED rather than
      ## warned about, and surfaced by historical_decomposition() only when
      ## something actually fails to add up.
      zero_variance_shocks = m$varexo_names[diag(as.matrix(Sigma_e)) == 0],
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
#' @param lik_init Initialisation of \eqn{P_{0|0}}:
#'   \describe{
#'     \item{\code{"auto"}}{(default) the unconditional Lyapunov covariance
#'       when it exists, and the EXACT diffuse smoother when it does not
#'       (unit roots in \eqn{T}).}
#'     \item{\code{"stationary"}}{demand the Lyapunov solution; error if it
#'       does not exist.}
#'     \item{\code{"diffuse"}}{the exact diffuse recursion
#'       (Koopman--Durbin, sequential form on the augmented state). On a model
#'       with no unit roots this IS the stationary initialisation, and the
#'       result says so in \code{diagnostics$lik_init_used}.}
#'     \item{\code{"kappa"}}{the old large-diagonal prior
#'       \eqn{P_0 = 10^6 I}. Kept for continuity with earlier releases and
#'       for comparing against \code{kalman_filter(lik_init = "kappa")};
#'       its log-likelihood carries an arbitrary additive constant, so it is
#'       not comparable across initialisations.}
#'   }
#'   \strong{Changed in 0.9.3.2:} \code{"auto"} used to fall back to
#'   \code{"kappa"} on a unit root, with a warning saying the loglik had a
#'   kappa-dependent offset. It now runs the exact recursion, whose
#'   log-likelihood matches
#'   \code{kalman_filter(lik_init = "diffuse", method = "univariate")} to
#'   machine precision. Smoothed states move by ~1e-8 (the fallback was
#'   accurate; the likelihood was the part that was not).
#' @param kalman_tol Conditional-variance floor below which an observation
#'   component is treated as carrying no information and dropped for that
#'   period, matching \code{\link{kalman_filter}}'s univariate fallback.
#' @param a0 Initial state mean \eqn{s_{0|0}}, length \code{n_state}, in
#'   \strong{deviations from the steady state} (the convention
#'   \code{smoothed_states} is in -- \code{data} is in levels, the states are
#'   not). \code{NULL} (default) starts at the steady state. Matched by name
#'   when named. The pre-sample \code{smoothed_initial} is then
#'   \eqn{s_{0|T}} under YOUR prior rather than under the model's
#'   unconditional one.
#' @param P0 Initial state covariance \eqn{P_{0|0}}, \code{n_state x n_state}
#'   (or a scalar for a multiple of the identity), symmetric and positive
#'   semi-definite. \code{NULL} (default) uses \code{lik_init}. Supplying it
#'   replaces that whole ladder, including the unit-root fallback and its
#'   warning: the caller has said what the prior is, so there is nothing left
#'   to fall back from.
#'
#'   Together, \code{a0} and \code{P0} are how a state is carried across a
#'   sample split -- take \code{filtered_states} and \code{filtered_cov} at
#'   the last period of the first block and hand them to the second, or take
#'   \code{final_state} / \code{final_cov} from \code{\link{kalman_filter}},
#'   which are the same pair under the same convention.
#'
#'   \code{a0}/\code{P0} and \code{pre_sample} answer different questions and
#'   compose: \code{a0}/\code{P0} SET the prior at \eqn{s_0} (what you know
#'   before the sample), while \code{pre_sample} ESTIMATES periods before
#'   \eqn{s_0} from the data that follows them. Use the first to carry a
#'   state forward across a split, the second to backcast a latent history.
#'   \code{P0} and \code{lik_init = "diffuse"} are mutually exclusive: they
#'   are two different priors, and passing both is an error rather than a
#'   silent preference.
#' @param pre_sample Number of periods BEFORE the first observation to
#'   backfill (default \code{0}). The latent history is estimated from the
#'   data that follows it, and the results come back in
#'   \code{presample_states}, \code{presample_shocks} and
#'   \code{presample_cov} -- chronological, so the last row is the period
#'   immediately before \code{data}. Every other returned series stays
#'   aligned with \code{data}.
#'
#'   No new recursion is involved: an all-missing period is predict-only, so
#'   this is the ordinary backward pass run over \code{pre_sample} padded
#'   rows -- the same mechanism that has always produced the single
#'   \code{smoothed_initial} period. The log-likelihood is unchanged (missing
#'   rows contribute nothing). Exact for a stationary model, where \eqn{P_0}
#'   is the unconditional covariance, and (since 0.9.3.2) exact on a unit-root
#'   model too, where it inherits the exact diffuse initialisation rather than
#'   the kappa fallback's arbitrary constant.
#' @param known_shocks Known historical shock values: an \code{n_exo x T}
#'   matrix carrying the value where a shock is known and \code{NA} where it is
#'   not -- the \code{NA}-as-unknown convention \code{data} uses, and the
#'   \code{n_exo x T} shape \code{shock_scale} uses. Rows are matched BY NAME
#'   when the matrix has rownames. \code{NULL} (default), or a matrix that is
#'   all \code{NA}, is a no-op.
#'
#'   Use it for a shock you actually know: an announced policy change, a
#'   measured intervention, a judgemental adjustment carried over from another
#'   exercise.
#'
#'   A known shock is a deterministic part of the system, so it splits off
#'   exactly: its trajectory is subtracted from the data, the ordinary
#'   recursion runs on the remainder, and the trajectory is added back. The
#'   injected values are returned in \code{smoothed_shocks} as themselves --
#'   they are inputs, not estimates. \code{loglik} is the CONDITIONAL
#'   \eqn{\log p(y \mid \varepsilon = v)}; see \code{\link{kalman_filter}}
#'   for the joint.
#' @param shock_means Deterministic shock MEANS, an \code{n_exo x T} matrix of
#'   mean shifts (\code{NA} and \code{0} both mean "no shift here"), rows
#'   matched by name. The same argument \code{\link{kalman_filter}} takes, with
#'   the same meaning and the same mechanism -- see there for the distinction
#'   from \code{known_shocks}, which is the point of having both.
#'
#'   One thing differs, and it follows from that distinction: a known shock's
#'   REALISATION is fixed, so \code{smoothed_shocks} reports the injected value
#'   back as itself, whereas a MEAN leaves the shock random, so
#'   \code{smoothed_shocks} reports \eqn{m_t + u_{t|T}} -- the mean plus the
#'   smoothed deviation around it. A smoother that returned the mean unrevised
#'   would be ignoring the data; one that ignored the mean would be ignoring
#'   the input.
#' @param shock_timing How to read the columns of \code{shock_means}:
#'   \code{"dated"} (default) or \code{"transition_next"}. See
#'   \code{\link{kalman_filter}}, whose \emph{Matching another package's
#'   shock timing} section gives the experiment that settles which one a given
#'   reference implementation uses.
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
#'   \eqn{t} is \eqn{s_{t|T}}), \code{smoothed_shocks}, and the two forward
#'   paths under the same names \code{\link{kalman_filter}} uses --
#'   \code{updated_states} (row \eqn{t} is \eqn{s_{t|t}}) and
#'   \code{predicted_states} (row \eqn{t} is \eqn{s_{t|t-1}}), the mean
#'   counterpart of \code{predicted_cov}. \code{filtered_states} is the same
#'   matrix as \code{updated_states}. Note the orientation differs from the
#'   filter's: here rows are periods. Also
#'   \code{filtered_cov} / \code{predicted_cov} / \code{smoothed_cov},
#'   \code{smoothed_initial} (\eqn{s_{0|T}}),
#'   \code{smoothed_initial_cov} (\eqn{V_{0|T}}), \code{loglik}, and
#'   \code{diagnostics} -- the same machine-readable record
#'   \code{\link{kalman_filter}} returns (see its \code{Value} section),
#'   with \code{method_requested = NA} because the smoother takes no
#'   \code{method} argument and \code{method_used} naming the recursion that
#'   ran (\code{"durbin-koopman"} or \code{"sequential-diffuse"}).
#'
#'   During a diffuse phase the reported \code{filtered_cov} /
#'   \code{predicted_cov} are the PROPER (\eqn{P_{star}}) part; the diffuse
#'   part is unbounded by construction. \code{diagnostics$diffuse_periods}
#'   says which periods those are.
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
                            lik_init = c("auto", "stationary", "kappa",
                                         "diffuse"),
                            kalman_tol = 1e-10, a0 = NULL, P0 = NULL,
                            pre_sample = 0L, known_shocks = NULL,
                            shock_means = NULL,
                            shock_timing = c("dated", "transition_next")) {
  lik_init <- match.arg(lik_init)
  shock_timing <- match.arg(shock_timing)

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
                      lik_init = lik_init, kalman_tol = kalman_tol,
                      a0 = a0, P0 = P0, pre_sample = pre_sample,
                      known_shocks = known_shocks, shock_means = shock_means,
                      shock_timing = shock_timing)
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
                                lik_init = c("auto", "stationary", "kappa",
                                             "diffuse"),
                                kalman_tol = 1e-10, a0 = NULL, P0 = NULL,
                                pre_sample = 0L, known_shocks = NULL,
                                shock_means = NULL,
                                shock_timing = c("dated", "transition_next")) {
  lik_init <- match.arg(lik_init)
  shock_timing <- match.arg(shock_timing)
  lik_init_orig <- lik_init        # for $diagnostics (R4)

  ## ---- Pre-sample backfill ------------------------------------------------
  ## Latent states BEFORE the first observation, which is what a "backcast" of
  ## the unobserved history means. No new recursion is needed: an all-missing
  ## period is already predict-only, so prepending `pre_sample` rows of NA lets
  ## the ordinary Durbin-Koopman backward pass estimate those periods from the
  ## data that follows them. That is exactly how the single pre-sample period
  ## in `smoothed_initial` has always been produced -- this generalises it to
  ## k periods and returns them separately, rather than making the caller pad
  ## the matrix by hand and re-align every output.
  ##
  ## Exact for a stationary model, where P_0 is the unconditional covariance.
  ## On a unit-root model it inherits the kappa fallback below, so the backfill
  ## carries the same kappa-dependent offset the loglik does.
  pre_sample <- as.integer(pre_sample)
  if (length(pre_sample) != 1L || is.na(pre_sample) || pre_sample < 0L)
    stop("kalman_smoother: `pre_sample` must be a single non-negative integer.",
         call. = FALSE)
  if (pre_sample > 0L) {
    n_obs_in <- ncol(data)
    data <- rbind(matrix(NA_real_, pre_sample, n_obs_in,
                         dimnames = list(NULL, colnames(data))),
                  as.matrix(data))
    ## me_extra is n_obs x T and shock_scale is n_shk x T: pad on the LEFT with
    ## the neutral value, or the per-period columns silently shift by k.
    if (!is.null(me_extra))
      me_extra <- cbind(matrix(0, nrow(me_extra), pre_sample), me_extra)
    if (!is.null(shock_scale))
      shock_scale <- cbind(matrix(1, nrow(shock_scale), pre_sample), shock_scale)
    ## known_shocks is n_exo x T on the CALLER's sample, and the validation
    ## below runs against the padded one -- so pad it here with "unknown"
    ## rather than making the caller pad by hand (and then wonder why the
    ## injected periods moved).
    if (!is.null(known_shocks)) {
      known_shocks <- as.matrix(known_shocks)
      known_shocks <- cbind(matrix(NA_real_, nrow(known_shocks), pre_sample,
                                   dimnames = list(rownames(known_shocks), NULL)),
                            known_shocks)
    }
    ## Same for a mean path: zero is "no shift", so the padded periods are
    ## unforced -- which is what a backcast of the history before the sample
    ## means when the sample's own inputs are known.
    if (!is.null(shock_means)) {
      shock_means <- as.matrix(shock_means)
      shock_means <- cbind(matrix(0, nrow(shock_means), pre_sample,
                                  dimnames = list(rownames(shock_means), NULL)),
                           shock_means)
    }
  }

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

  ## ---- Known historical shocks -------------------------------------------
  ## A known shock is a DETERMINISTIC part of the system, so it splits off
  ## exactly and the recursion never has to know about it. With
  ## eps_{j,t} = v_t known,
  ##   s_t = T s_{t-1} + R_j v_t + R_-j eps_{-j,t}
  ##   y_t = Z s_{t-1} + D_j v_t + D_-j eps_{-j,t}
  ## so subtracting the deterministic trajectory from the data leaves an
  ## ordinary smoothing problem in the unknown shocks, and the trajectory is
  ## added back afterwards. The known shock's variance is switched off in those
  ## periods -- there is nothing left to estimate about it.
  ##
  ## This is a different mechanism from kalman_filter()'s, which observes eps
  ## directly on its augmented state. They must agree, and a test pins that.
  known_sm <- NULL
  if (!is.null(known_shocks)) {
    known_sm <- .kf_known_shocks(known_shocks, NULL, ss$shock_names, nrow(data),
                                 what = "kalman_smoother")
  }
  if (!is.null(known_sm)) {
    kv <- known_sm$values; kv[is.na(kv)] <- 0        # 0 outside the known periods
    kidx <- known_sm$idx
    s_det <- matrix(0, nrow(data) + 1L, ss$n_state)  # row t+1 holds s^det_t
    y_det <- matrix(0, nrow(data), ss$n_obs)
    for (t in seq_len(nrow(data))) {
      y_det[t, ] <- as.numeric(ss$Z_mat %*% s_det[t, ] +
                               ss$D_mat[, kidx, drop = FALSE] %*% kv[, t])
      s_det[t + 1L, ] <- as.numeric(ss$T_mat %*% s_det[t, ] +
                                    ss$R_mat[, kidx, drop = FALSE] %*% kv[, t])
    }
    data <- data - y_det
    ## Switch the known shocks off where they are known: shock_scale already
    ## carries exactly this per-period, per-shock semantics.
    sc <- if (is.null(shock_scale)) matrix(1, ss$n_shock, nrow(data))
          else shock_scale
    sc[kidx, ] <- sc[kidx, , drop = FALSE] * (is.na(known_sm$values) + 0)
    shock_scale <- sc
  }

  ## ---- Deterministic shock MEANS -----------------------------------------
  ## The same split as above, for a different statement. `known_shocks` fixes
  ## the REALISATION, so the smoother reports the injected value back as
  ## itself. `shock_means` fixes the MEAN and the shock keeps its variance, so
  ## the smoother still estimates the deviation around it and reports
  ## eps_{t|T} = m_t + u_{t|T}. Anything else would be ignoring either the
  ## input or the data. See .kf_shock_means() in R/kalman-filter.R.
  ##
  ## Both entry points run before ss_convert_timing(), which is a no-op for
  ## every state space this package builds (`timing = "lagged"`); a hand-built
  ## dsge_ss declaring another timing would need them moved after it.
  mean_path <- .kf_shock_means(shock_means, shock_timing, ss$shock_names,
                               nrow(data), what = "kalman_smoother")
  if (!is.null(mean_path) && !is.null(known_sm)) {
    clash <- !is.na(known_sm$values) &
      mean_path[known_sm$idx, , drop = FALSE] != 0
    if (any(clash))
      stop(sprintf(paste0("kalman_smoother: `shock_means` and `known_shocks` ",
                          "both specify %s. A known shock's REALISATION is ",
                          "fixed, so its mean is already determined -- pass ",
                          "one or the other for a given shock and period."),
                   paste(known_sm$names[which(apply(clash, 1L, any))],
                         collapse = ", ")), call. = FALSE)
  }
  m_det <- NULL
  if (!is.null(mean_path)) {
    m_det <- matrix(0, nrow(data) + 1L, ss$n_state)   # row t+1 holds s^det_t
    y_md  <- matrix(0, nrow(data), ss$n_obs)
    for (t in seq_len(nrow(data))) {
      y_md[t, ] <- as.numeric(ss$Z_mat %*% m_det[t, ] +
                              ss$D_mat %*% mean_path[, t])
      m_det[t + 1L, ] <- as.numeric(ss$T_mat %*% m_det[t, ] +
                                    ss$R_mat %*% mean_path[, t])
    }
    data <- data - y_md
  }

  ## ---- Structured run diagnostics (R4) -----------------------------------
  ## Same field names as kalman_filter()$diagnostics, so a parity harness can
  ## read either without a special case. `method_requested` is NA here because
  ## the smoother takes no `method` argument -- which recursion ran is reported
  ## in `method_used`.
  .smoother_diagnostics <- function(lik_init_used, d_diffuse = NA_integer_,
                                    dropped = NULL, data = NULL,
                                    method_used = "sequential-diffuse",
                                    routing = list()) {
    n_per <- if (is.null(data)) 0L else nrow(data)
    miss  <- if (is.null(data)) integer(0)
             else as.integer(rowSums(is.na(as.matrix(data))))
    if (is.null(dropped)) dropped <- integer(n_per)
    if (!identical(lik_init_orig, lik_init_used))
      routing <- c(routing, list(c(
        from = lik_init_orig, to = lik_init_used,
        reason = if (identical(lik_init_used, "diffuse"))
          "unit root(s) in T: the unconditional state covariance does not exist"
        else if (identical(lik_init_used, "stationary"))
          "no unit roots: the exact diffuse initialisation is the stationary one"
        else "requested initialisation was not available")))
    routing <- .kf_routing_df(routing)
    dd <- if (length(d_diffuse) != 1L || is.na(d_diffuse)) NA_integer_
          else as.integer(d_diffuse)
    list(method_requested   = NA_character_,
         method_used        = method_used,
         lik_init_requested = lik_init_orig,
         lik_init_used      = lik_init_used,
         routing            = routing,
         diffuse_periods    = if (is.na(dd)) integer(0) else seq_len(dd),
         missing_by_period  = miss,
         n_missing          = sum(miss),
         dropped_by_period  = as.integer(dropped),
         n_dropped          = sum(as.integer(dropped)),
         known_shocks       = if (is.null(known_sm)) NULL else
           list(names = known_sm$names,
                n_cells = sum(!is.na(known_sm$values)),
                n_applied = sum(!is.na(known_sm$values))),
         shock_means        = if (is.null(mean_path)) NULL else
           list(timing = shock_timing,
                n_cells = sum(mean_path != 0),
                names = ss$shock_names[rowSums(mean_path != 0) > 0]),
         ## The smoother conditions on an injected shock (it splits the
         ## deterministic trajectory off); kalman_filter() reports the joint.
         ## A mean path conditions nothing -- it is an input, not an event.
         loglik_type        = if (is.null(known_sm)) "marginal" else "conditional")
  }

  ## ---- Shared exit ------------------------------------------------------
  ## Both recursions (the DK backward pass below and the exact-diffuse
  ## sequential smoother) produce the same `out`, and both owe the caller the
  ## same two corrections afterwards: add the deterministic known-shock
  ## trajectory back, and split the padded pre-sample rows out so every
  ## returned series is aligned with the data that was passed in. Doing it in
  ## one place is what keeps the two paths returning the same object.
  .smoother_finish <- function(out) {
    if (!is.null(mean_path)) {
      out$smoothed_states  <- out$smoothed_states  + m_det[-1L, , drop = FALSE]
      out$filtered_states  <- out$filtered_states  + m_det[-1L, , drop = FALSE]
      out$predicted_states <- out$predicted_states + m_det[-1L, , drop = FALSE]
      out$smoothed_initial <- out$smoothed_initial + m_det[1L, ]
      ## The shock is NOT fixed by a mean, so the estimate is the mean plus
      ## the smoothed deviation around it.
      out$smoothed_shocks  <- out$smoothed_shocks + t(mean_path)
    }
    if (!is.null(known_sm)) {
      out$smoothed_states  <- out$smoothed_states + s_det[-1L, , drop = FALSE]
      out$filtered_states  <- out$filtered_states + s_det[-1L, , drop = FALSE]
      out$smoothed_initial <- out$smoothed_initial + s_det[1L, ]
      kv_na <- known_sm$values
      for (i in seq_along(known_sm$idx)) {
        hit <- which(!is.na(kv_na[i, ]))
        out$smoothed_shocks[hit, known_sm$idx[i]] <- kv_na[i, hit]
      }
    }
    if (pre_sample > 0L) {
      k  <- pre_sample
      ix <- seq_len(k)
      n_all <- nrow(out$smoothed_states)
      keep  <- (k + 1L):n_all
      ## Read from the CORRECTED series (post add-back), or a pre-sample
      ## backfill run together with known_shocks would report the padded
      ## periods without the deterministic trajectory in them.
      out$presample_states <- out$smoothed_states[ix, , drop = FALSE]
      out$presample_shocks <- out$smoothed_shocks[ix, , drop = FALSE]
      out$presample_cov    <- out$smoothed_cov[, , ix, drop = FALSE]
      out$smoothed_states  <- out$smoothed_states[keep, , drop = FALSE]
      out$smoothed_shocks  <- out$smoothed_shocks[keep, , drop = FALSE]
      out$filtered_states  <- out$filtered_states[keep, , drop = FALSE]
      out$predicted_states <- out$predicted_states[keep, , drop = FALSE]
      out$filtered_cov    <- out$filtered_cov[, , keep, drop = FALSE]
      out$predicted_cov   <- out$predicted_cov[, , keep, drop = FALSE]
      out$smoothed_cov    <- out$smoothed_cov[, , keep, drop = FALSE]
      out$pre_sample      <- k
      ## The per-period diagnostics describe the CALLER's sample too: the
      ## padded rows are all-missing by construction and counting them as
      ## missing observations would be an artefact of the padding.
      if (!is.null(out$diagnostics)) {
        d <- out$diagnostics
        d$presample_periods <- k
        d$missing_by_period <- d$missing_by_period[keep]
        d$dropped_by_period <- d$dropped_by_period[keep]
        d$n_missing <- sum(d$missing_by_period)
        d$n_dropped <- sum(d$dropped_by_period)
        out$diagnostics <- d
      }
    }
    ## The timing contract, in the names, matching kalman_filter(): row t of
    ## `updated_states` is s_{t|t} and row t of `predicted_states` is
    ## s_{t|t-1}. `filtered_states` is the same matrix as `updated_states`,
    ## under the name the rest of the package uses. Aliased LAST, after every
    ## correction above, so the two cannot drift apart.
    out$updated_states <- out$filtered_states
    out
  }

  ## Counters for the singular-F diagnostic raised after the forward pass.
  n_sing_periods <- 0L
  n_sing_dropped <- 0L
  sing_by_period <- integer(nrow(data))

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

  .kf_warn_zero_shock_cov(Q, "kalman_smoother")

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
  ## A supplied P0 short-circuits the whole initialisation ladder below --
  ## including its unit-root warning, which is about a fallback that no longer
  ## applies once the caller has said what the prior is.
  P0_user <- .kf_init_cov(P0, ss$state_names, n_s, what = "kalman_smoother")
  a0_user <- .kf_init_mean(a0, ss$state_names, n_s, what = "kalman_smoother")

  P_ss <- if (!is.null(P0_user)) {
    P0_user
  } else if (identical(lik_init, "kappa")) {
    ## Forced diffuse-style prior. kalman_filter() REFUSES lik_init = "auto"
    ## together with shock_scale on a nonstationary model and tells the caller
    ## to pass "kappa" or "stationary" explicitly; before this argument existed
    ## there was no way to ask the smoother for the same thing, so the two
    ## could not be made comparable even in principle.
    matrix(NA_real_, n_s, n_s)
  } else {
    solve_lyapunov(TT_mat, RQR)
  }
  if (!is.null(P0_user) && identical(lik_init, "diffuse"))
    stop("kalman_smoother: `P0` and lik_init = \"diffuse\" are two different ",
         "initialisations -- the exact-diffuse recursion builds its own ",
         "(P_inf, P_star) split and has nothing to do with a supplied P0. ",
         "Pass one or the other. `a0` composes with either.", call. = FALSE)
  if (is.null(P0_user) && identical(lik_init, "stationary") && anyNA(P_ss))
    stop("kalman_smoother: lik_init = \"stationary\" was requested but the ",
         "Lyapunov solve returned NaN -- TT has unit-root eigenvalues, so the ",
         "unconditional state covariance does not exist. Use lik_init = ",
         "\"kappa\" (or \"auto\") for a nonstationary model.", call. = FALSE)
  if (anyNA(P_ss)) {
    ## Unit roots detected. This used to warn and substitute a large finite
    ## prior P0 = 1e6 * I -- an approximating sequence whose smoothed states
    ## are close (measured: ~5e-8 on the local-level fixture) but whose
    ## log-likelihood carries an arbitrary kappa-dependent additive constant,
    ## so it could not be compared with anything. lik_init = "auto" now runs
    ## the EXACT diffuse smoother instead (.smoother_diffuse_seq); "kappa"
    ## still asks for the old finite prior explicitly.
    if (identical(lik_init, "kappa")) P_ss <- .DIFFUSE_SCALE * diag(n_s)
    else                              lik_init <- "diffuse"
  }

  ## ---- Exact diffuse smoothing -------------------------------------------
  ## The sequential smoother in R/smoother-diffuse.R does the whole job --
  ## forward and backward -- on the augmented state, so this branch builds its
  ## inputs and leaves through the shared exit rather than continuing into the
  ## multivariate DK pass below.
  if (identical(lik_init, "diffuse")) {
    dp <- .kf_diffuse_P0(TT_mat, RQR)
    if (dp$nunit == 0L) {
      ## Nothing is actually diffuse: the exact answer IS the stationary one.
      lik_init <- "stationary"
      P_ss     <- solve_lyapunov(TT_mat, RQR)
    } else {
      nb   <- n_s + n_shk
      s_ix <- seq_len(n_s); e_ix <- n_s + seq_len(n_shk)
      Sig_list <- lapply(seq_len(TT), function(t)
        if (has_shock_scale) { sc <- shock_scale[, t]; Q * outer(sc, sc) } else Q)
      Zb  <- cbind(Z_mat, D_mat)
      Tb  <- rbind(cbind(TT_mat, R_mat), matrix(0, n_shk, nb))
      Gm  <- cbind(TT_mat, R_mat)
      Ps1 <- matrix(0, nb, nb)
      Ps1[s_ix, s_ix] <- dp$P_star; Ps1[e_ix, e_ix] <- Sig_list[[1L]]
      Pi1 <- matrix(0, nb, nb); Pi1[s_ix, s_ix] <- dp$P_inf
      me_mat <- matrix(me_variance, n_obs, TT)
      if (has_me_extra) me_mat <- me_mat + me_extra
      Ydev <- t(as.matrix(data))
      if (!is.null(d_obs)) Ydev <- Ydev - d_obs
      ds <- .smoother_diffuse_seq(Ydev, Zb, Tb, Gm, Sig_list,
                                  c(a0_user, numeric(n_shk)), Ps1, Pi1,
                                  me_mat, s_ix, e_ix, kalman_tol = kalman_tol)
      if (isTRUE(ds$diffuse_failed))
        warning("kalman_smoother: the diffuse phase did not end within the ",
                "sample -- P_inf never decayed, so some diffuse direction is ",
                "not identified by the data. The smoothed states are the ",
                "minimum-norm answer in that direction and the reported ",
                "covariances carry the proper part only. This usually means ",
                "an unobserved unit root: check that every nonstationary ",
                "state is loaded by some observable.", call. = FALSE)
      nm  <- ss$state_names; shk <- ss$shock_names
      dn3 <- list(nm, nm, NULL)
      colnames(ds$smoothed_states)  <- nm
      colnames(ds$filtered_states)  <- nm
      colnames(ds$predicted_states) <- nm
      colnames(ds$smoothed_shocks) <- shk
      names(ds$smoothed_initial)   <- nm
      dimnames(ds$smoothed_initial_cov) <- list(nm, nm)
      dimnames(ds$filtered_cov)  <- dn3
      dimnames(ds$predicted_cov) <- dn3
      dimnames(ds$smoothed_cov)  <- dn3
      out <- list(
        smoothed_states = ds$smoothed_states,
        smoothed_shocks = ds$smoothed_shocks,
        filtered_states = ds$filtered_states,
        predicted_states = ds$predicted_states,
        filtered_cov    = ds$filtered_cov,
        predicted_cov   = ds$predicted_cov,
        smoothed_cov    = ds$smoothed_cov,
        P_filt_last     = ds$filtered_cov[, , TT],
        smoothed_initial     = ds$smoothed_initial,
        smoothed_initial_cov = ds$smoothed_initial_cov,
        loglik          = ds$loglik,
        diagnostics     = .smoother_diagnostics(
          lik_init_used = "diffuse",
          ## An unfinished diffuse phase means EVERY period is still in it.
          d_diffuse = if (isTRUE(ds$diffuse_failed)) TT else ds$d_diffuse,
          dropped = ds$n_skipped, data = data))
      return(.smoother_finish(out))
    }
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

  ## s_{0|0}. Zero is the steady state in this package's deviation convention;
  ## `a0` moves it, which is what makes a hand-off from an earlier sample (or
  ## from smoother2histval()) possible.
  s_tt <- a0_user
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
    ## A successful chol() does not mean F is safely invertible: on a
    ## stochastically singular system it can succeed with a pivot at round-off
    ## and return a badly wrong update. Measured on a 2-observable / 1-shock
    ## fixture, the drop path fired in 19 of 20 periods and the ONE period
    ## where chol() happened to succeed carried the entire error -- smoothed
    ## states exact, smoothed shocks inconsistent with them by 0.18, which
    ## surfaced downstream as a historical_decomposition() adding-up residual.
    ## The pivots are the conditional variances; test them.
    if (!is.null(F_ch) && .kf_F_singular(F_ch, F_t, kalman_tol)) F_ch <- NULL

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
      sing_by_period[t] <- sum(!keep)

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
    ## t == 1 reaches BEFORE the sample: the pair here is (s_{0|0}, P_{0|0}),
    ## i.e. the prior itself. This read zero unconditionally, which silently
    ## ignored `a0` on the backward pass while the forward pass honoured it --
    ## visible only in smoothed_initial, and only once a0 was non-zero.
    s_in <- if (t == 1L) a0_user else s_filt[t - 1L, ]
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
  colnames(s_pred)     <- ss$state_names
  colnames(eps_smooth) <- ss$shock_names
  names(s0_smooth)     <- ss$state_names
  dimnames(V0_smooth)  <- list(ss$state_names, ss$state_names)

  dimnames(P_filt)   <- list(ss$state_names, ss$state_names, NULL)
  dimnames(P_pred)   <- list(ss$state_names, ss$state_names, NULL)
  dimnames(V_smooth) <- list(ss$state_names, ss$state_names, NULL)

  out <- list(
    smoothed_states = s_smooth,
    smoothed_shocks = eps_smooth,
    filtered_states = s_filt,
    ## The mean counterpart of predicted_cov, and the same timing contract
    ## kalman_filter() reports: row t is s_{t|t-1}. It was computed all along
    ## and simply not returned.
    predicted_states = s_pred,
    ## Per-period state covariances (n_state x n_state x T):
    filtered_cov    = P_filt,            # P_{t|t}
    predicted_cov   = P_pred,            # P_{t|t-1}
    smoothed_cov    = V_smooth,          # P_{t|T} (DK)
    P_filt_last     = P_filt[, , TT],    # P_{T|T} (kept for back-compat)
    ## Pre-sample smoothed moments, free from the same backward recursion.
    smoothed_initial     = s0_smooth,    # s_{0|T}
    smoothed_initial_cov = V0_smooth,    # V_{0|T}
    loglik          = loglik,
    diagnostics     = .smoother_diagnostics(
      lik_init_used = if (!is.null(P0_user)) "user" else
                      if (identical(lik_init, "kappa")) "kappa" else "stationary",
      dropped = sing_by_period, data = data,
      method_used = "durbin-koopman")
  )

  .smoother_finish(out)
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
#' @param tol Relative tolerance for the adding-up check (default
#'   \code{1e-8}). For a linear model the residual is round-off, ~1e-15;
#'   anything above \code{tol * max(abs(path))} warns and sets
#'   \code{$adding_up_ok} to \code{FALSE}.
#' @section The integration contract:
#' Everything below is what has to line up for the components to add up. It is
#' now CHECKED rather than assumed -- names are matched where they exist and a
#' mismatch is an error -- but it is worth stating, because the failure mode
#' when it does not hold is a large \code{adding_up_residual} with every
#' dimension still correct.
#'
#' \strong{Orientation.} \code{smoothed_shocks} is \code{T x n_shock} and
#' \code{smoothed_states} is \code{T x n_state}: \strong{rows are periods}.
#' That is what \code{\link{kalman_smoother}} returns. Note
#' \code{\link{kalman_filter}}'s state matrices are the other way round
#' (\code{n_state x T}) and need \code{t()}; a transposed argument is
#' detected and refused rather than silently reinterpreted.
#'
#' \strong{Names.} Columns are matched to \code{ss$shock_names} /
#' \code{ss$state_names} by name when the matrix has dimnames, and \code{s0}
#' by name when it is named. A correctly-sized but wrongly-ORDERED input used
#' to be accepted silently; it is now reordered when it can be identified and
#' refused when it cannot.
#'
#' \strong{The state at the first contribution period.} Row \eqn{t} of every
#' contribution is the model's variables \strong{dated} \eqn{t}, namely
#' \eqn{ghx\, s_{t-1} + ghu\, \varepsilon_t}. Row 1 therefore loads the
#' PRE-SAMPLE state \eqn{s_0}, not \eqn{s_1} -- which is why the
#' \code{"initial"} column exists and why \code{s0} must be the smoother's
#' \eqn{s_{0|T}} (\code{\link{smoothed_initial_state}}), not the first
#' smoothed state. Passing the whole \code{\link{kalman_smoother}} result
#' takes it for you, and is the only call shape that cannot get this wrong.
#'
#' \strong{Pre- or post-transition.} Pre-: the state entering the period, plus
#' that period's shock. No contribution is the post-transition state.
#'
#' \strong{Deterministic inputs.} \code{shock_means} and \code{known_shocks}
#' need nothing here. Both arrive through the smoother's output -- the mean or
#' injected value is already inside \code{smoothed_shocks}, and its effect on
#' the path is already inside \code{smoothed_states} -- so a deterministic
#' input shows up in ITS OWN shock's column and the adding-up is unaffected.
#' The \code{shock_timing} choice is likewise settled upstream, in the
#' smoother call, and never re-enters here.
#'
#' \strong{The two diagnostics, and why there are two.}
#' \code{$adding_up_residual} compares the components -- propagated
#' internally by \code{T_mat}/\code{R_mat} from \code{s0} -- against a path
#' rebuilt from the smoother's OWN states. It is always present: a number when
#' \code{smoothed_states} is available, \code{NA} with
#' \code{$adding_up_note} when it is not.
#'
#' It is \strong{necessary but not sufficient}. The contemporaneous
#' \eqn{ghu\,\varepsilon_t} term appears identically on both sides and
#' cancels, so a shock error is visible only through its propagated
#' (\eqn{t+1} onward) effect -- and the LAST period is therefore not checked
#' at all. Measured on \code{nk_demo}: perturbing the final period's shock by
#' 1.0, or every shock in that period by 50\%, or the final state by 1.0,
#' leaves \code{adding_up_residual} at 5.9e-15, while the same perturbation
#' mid-sample shows as 2.6 and 1.0.
#'
#' \code{$transition_residual} has neither blind spot. It asks the direct
#' question -- does the smoother's own output satisfy its own transition,
#' \eqn{s_t = T s_{t-1} + R \varepsilon_t}, period by period? -- and catches
#' all six of those perturbations. \strong{Look at it first} when a
#' decomposition will not add up: it distinguishes "the decomposition is
#' wrong" from "its INPUTS are incoherent", and
#' \code{$transition_residual_by_period} with
#' \code{$transition_worst_period} localises the latter. That is how the
#' singular-innovation defect fixed in 0.9.3.5 was pinned to the single period
#' in which \code{chol()} had happened to succeed.
#'
#' Both compare against \code{tol} (relative to the path's scale) and warn on
#' failure. A drop of predictable observation components -- see the smoother's
#' \code{$diagnostics$dropped_by_period} -- is the usual context for a
#' transition failure.
#'
#' @section Relation to IRIS simulate(..., 'contributions', true):
#' Verified equal to \strong{2e-16} on a two-shock linear model (IRIS Toolbox
#' Release 20180308): each isolated-shock column, the initial-condition column,
#' and the total. Two of IRIS's columns have no dynhr counterpart, by
#' construction rather than omission:
#' \itemize{
#'   \item a \strong{measurement-shock} column. dynhr's \code{me_variance}
#'     is observation noise, not a structural shock, so it gets no column and
#'     the decomposition is of the MODEL variable rather than the observation.
#'     Compare against IRIS's structural columns plus its init column, not
#'     against its total, when the model has measurement shocks.
#'   \item a \strong{nonlinear} column, which is identically zero for a
#'     linear model. \code{$has_nonlinear_column} and
#'     \code{$has_residual_column} are \code{FALSE} for the same reason: the
#'     adding-up here is exact, so such a column would carry nothing.
#' }
#' IRIS's \code{Init+Const+Trends} column corresponds to \code{"initial"}
#' when the model is written in deviations (dynhr's contributions always are;
#' add \code{ss$ys[v]} to read variable \code{v} back in levels).
#'
#' @return List of class \code{dynhr_shock_decomposition} with
#'   \code{$contributions} (named list of \code{T x n_endo} matrices: one per
#'   shock or group, plus \code{"initial"}), \code{$total} (\code{T x n_endo},
#'   the sum of all components = the smoothed series in deviations from
#'   steady state), \code{$initial}, \code{$s0}, \code{$components},
#'   \code{$adding_up_residual} / \code{$adding_up_relative} /
#'   \code{$adding_up_ok} / \code{$adding_up_tol},
#'   \code{$transition_residual} / \code{$transition_residual_by_period} /
#'   \code{$transition_worst_period} / \code{$transition_ok}, \code{$timing} (the
#'   dating convention in words), and \code{$has_nonlinear_column} /
#'   \code{$has_residual_column}.
#' @export
# ---------------------------------------------------------------------------
historical_decomposition <- function(smoothed_shocks, ss, s0 = NULL,
                                     smoothed_states = NULL,
                                     shock_groups = NULL, model = NULL,
                                     tol = 1e-8) {

  ## Accept a whole kalman_smoother() result: shocks, states and s_{0|T}.
  if (is.list(smoothed_shocks) && !is.matrix(smoothed_shocks) &&
      !is.null(smoothed_shocks$smoothed_shocks)) {
    sm <- smoothed_shocks
    if (is.null(smoothed_states)) smoothed_states <- sm$smoothed_states
    if (is.null(s0)) s0 <- smoothed_initial_state(sm, ss)
    smoothed_shocks <- sm$smoothed_shocks
  }

  ## ---- The integration contract, enforced rather than documented ---------
  ## Everything here is T x k with k in the state space's own order, and the
  ## adding-up is exact only if that holds. It used to be checked by LENGTH
  ## alone, so a correctly-sized but wrongly-ordered `s0` -- or shock columns
  ## in a different order from `ss$shock_names` -- was accepted silently and
  ## came back as a large `adding_up_residual` with nothing to say why. The
  ## residual is exactly the size of the initial-condition error, so on a
  ## model with big states (or a kappa-initialised smoother, whose s_{0|T}
  ## carries the arbitrary prior) it can reach 1e9 while every dimension
  ## still checks out. Match by NAME wherever names are present, and refuse
  ## rather than answer a different question -- the same rule .kf_init_mean()
  ## applies to `a0`.
  .hd_match <- function(x, want, what, kind) {
    if (is.null(x)) return(NULL)
    nm <- if (is.matrix(x)) colnames(x) else names(x)
    if (is.null(nm) || is.null(want)) return(x)
    if (!setequal(nm, want))
      stop(sprintf(paste0("historical_decomposition: `%s` %s do not match the ",
                          "state space's %s.\n  supplied: %s\n  expected: %s"),
                   what, if (is.matrix(x)) "column names" else "names", kind,
                   paste(nm, collapse = ", "), paste(want, collapse = ", ")),
           call. = FALSE)
    if (is.matrix(x)) x[, want, drop = FALSE] else x[want]
  }
  .hd_orient <- function(x, n_col, what, kind) {
    if (is.null(x) || !is.matrix(x)) return(x)
    if (ncol(x) == n_col) return(x)
    if (nrow(x) == n_col)
      stop(sprintf(paste0("historical_decomposition: `%s` looks transposed -- ",
                          "it is %d x %d and this function takes T x %s (rows ",
                          "are periods). kalman_smoother() returns it that way ",
                          "already; kalman_filter()'s state matrices are the ",
                          "other way round (n_state x T) and need t()."),
                   what, nrow(x), ncol(x), kind), call. = FALSE)
    stop(sprintf("historical_decomposition: `%s` has %d columns, expected %d (%s).",
                 what, ncol(x), n_col, kind), call. = FALSE)
  }

  smoothed_shocks <- .hd_orient(smoothed_shocks, ss$n_shock,
                                "smoothed_shocks", "n_shock")
  smoothed_states <- .hd_orient(smoothed_states, ss$n_state,
                                "smoothed_states", "n_state")
  smoothed_shocks <- .hd_match(smoothed_shocks, ss$shock_names,
                               "smoothed_shocks", "shock names")
  smoothed_states <- .hd_match(smoothed_states, ss$state_names,
                               "smoothed_states", "state names")
  s0 <- .hd_match(s0, ss$state_names, "s0", "state names")

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
              orientation   = "time_endo",
              ## The component inventory, stated rather than inferred: for a
              ## LINEAR model these are exactly the shocks (or groups) plus
              ## "initial", and there is no nonlinear or residual column --
              ## the adding-up is exact, so such a column would be zero by
              ## construction. `adding_up_residual` is the check.
              n_components  = length(contributions),
              has_nonlinear_column = FALSE,
              has_residual_column  = FALSE,
              timing = paste(
                "contribution row t is the model's variables DATED t:",
                "ghx %*% s_{t-1} + ghu %*% eps_t. Row 1 therefore loads the",
                "PRE-SAMPLE state s_0 (the `initial` column) and the first",
                "smoothed shock. Values are pre-transition in that sense: the",
                "state entering the period, not the one leaving it."))

  ## ---- Adding-up check against the directly-reconstructed smoothed path ---
  ## ALWAYS reported. It used to appear only when `smoothed_states` was
  ## supplied, so the one call shape that can be silently wrong -- a bare
  ## shock matrix with a defaulted or mismatched s0 -- was also the one with
  ## no diagnostic at all. Absent inputs now give NA and a reason, not a
  ## missing field.
  out$adding_up_tol <- tol
  if (is.null(smoothed_states)) {
    out$adding_up_residual <- NA_real_
    out$adding_up_ok       <- NA
    out$adding_up_note     <- paste(
      "not checked: `smoothed_states` was not supplied, so there is no",
      "independent path to check the components against. Pass the whole",
      "kalman_smoother() result (or its $smoothed_states) to get the check.")
  } else {
    smoothed <- matrix(0, TT, n_end)
    s_prev   <- s0
    ## ---- Transition coherence, checked DIRECTLY ---------------------------
    ## The adding-up check compares the components (propagated internally by
    ## T_mat/R_mat from s0) against a path rebuilt from the SMOOTHER's own
    ## states -- so it sees an incoherence between those two only through the
    ## `ghx %*% state` channel, and it has two blind spots that matter:
    ##
    ##   * the contemporaneous `ghu %*% eps_t` term appears IDENTICALLY on
    ##     both sides and cancels, so a shock error is visible only through
    ##     its propagated (t+1 onward) effect, damped by T;
    ##   * consequently the LAST period is not checked at all. Measured on
    ##     nk_demo: perturbing the final period's shock by 1.0, or every shock
    ##     in that period by 50%, or the final state by 1.0, leaves
    ##     adding_up_residual at 5.9e-15 -- completely undetected, while the
    ##     same perturbation mid-sample shows up as 2.6 and 1.0.
    ##
    ## So a small adding-up residual is NECESSARY but NOT SUFFICIENT evidence
    ## that the inputs are coherent. The direct test has neither blind spot:
    ## does the smoother's own output satisfy its own transition,
    ## s_t = T s_{t-1} + R eps_t, period by period? That is the quantity that
    ## actually broke under the singular-F defect fixed in 0.9.3.5, and it is
    ## the one to look at first when a decomposition will not add up -- the
    ## per-period vector localises it, which is how that defect was pinned to
    ## the single period where chol() had succeeded.
    tres   <- numeric(TT)
    tprev  <- s0
    for (t in seq_len(TT)) {
      eps_t <- smoothed_shocks[t, ]
      eps_t[is.na(eps_t)] <- 0
      smoothed[t, ] <- as.numeric(ghx %*% s_prev + ghu %*% eps_t)
      s_prev        <- smoothed_states[t, ]
      tres[t] <- max(abs(as.numeric(smoothed_states[t, ]) -
                         as.numeric(TT_mat %*% tprev + R_mat %*% eps_t)))
      tprev   <- as.numeric(smoothed_states[t, ])
    }
    colnames(smoothed)     <- ss$endo_names
    out$smoothed           <- smoothed
    resid                  <- max(abs(total - smoothed))
    scale                  <- max(1, max(abs(smoothed)))
    out$adding_up_residual <- resid
    out$adding_up_relative <- resid / scale
    out$adding_up_ok       <- resid <= tol * scale

    sscale <- max(1, max(abs(smoothed_states)))
    out$transition_residual            <- max(tres)
    out$transition_residual_by_period  <- tres
    out$transition_worst_period        <- which.max(tres)
    out$transition_ok <- max(tres) <= tol * sscale
    if (!isTRUE(out$transition_ok))
      warning(sprintf(paste0(
        "historical_decomposition: the smoother's own states and shocks do ",
        "not satisfy the transition -- max |s_t - T s_{t-1} - R eps_t| = ",
        "%.3g (%.3g relative), worst at period %d of %d, against a tolerance ",
        "of %.3g. The decomposition is faithfully reporting an incoherence in ",
        "its INPUTS, not creating one: see $transition_residual_by_period to ",
        "localise it. A drop of predictable observation components (see the ",
        "smoother's $diagnostics$dropped_by_period) is the usual context%s."),
        max(tres), max(tres) / sscale, which.max(tres), TT, tol,
        if (length(ss$zero_variance_shocks))
          paste0("; note that these shocks carry NO variance, which is the ",
                 "usual reason and is often a missing or partial `shocks;` ",
                 "block: ",
                 paste(ss$zero_variance_shocks, collapse = ", "))
        else ""),
        call. = FALSE)

    if (!isTRUE(out$adding_up_ok))
      warning(sprintf(paste0(
        "historical_decomposition: the components do not add up -- residual ",
        "%.3g (%.3g relative to the path's scale), against a tolerance of ",
        "%.3g. For a LINEAR model this should be round-off (~1e-15). The ",
        "residual is exactly the size of the initial-condition error, so the ",
        "usual causes are: `s0` is not the smoother's s_{0|T} (pass the whole ",
        "kalman_smoother() result and it is taken for you); the shocks or ",
        "states came from a different state space than `ss`; or the model is ",
        "not linear, in which case the pruned-state-space decomposition is ",
        "the right tool."), resid, out$adding_up_relative, tol), call. = FALSE)
  }
  if (is.null(smoothed_states)) {
    out$transition_residual <- NA_real_
    out$transition_ok       <- NA
  }

  structure(out, class = c("dynhr_shock_decomposition", "list"))
}
