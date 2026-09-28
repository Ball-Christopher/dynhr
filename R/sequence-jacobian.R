## R/sequence-jacobian.R
## --------------------------------------------------------------------------
## Representative-agent sequence-space Jacobian (RA-SSJ) spike.
##
## For a model's equilibrium map H(X, eps) = 0 over horizon T (X stacks
## endogenous paths, eps is the exogenous shock path), the sequence-space
## Jacobian is:
##
##   J = dH/dX  (T*n_endo x T*n_endo, block-tridiagonal)
##
## evaluated at the deterministic steady state (all periods at ss, eps = 0).
## This is IDENTICAL to the stacked Jacobian already assembled by the
## perfect-foresight Newton solver at convergence.
##
## The general-equilibrium impulse response to a unit shock at t=1 is:
##
##   dX = -J^{-1} * (dH/deps_1) * scale
##
## where dH/deps_1 is a T*n_endo vector: in period 1 it is the partial of
## the dynamic residuals w.r.t. eps_a (extracted from jacobian_fn's exo
## columns); in all other periods it is zero (AR shock enters through the
## state, not directly as an exogenous perturbation in t > 1).
##
## Validation oracle: the SSJ IRF must match the order-1 perturbation IRF
## (ghx/ghu recursion) to ~1e-8 (exact linear identity for a linear model).
##
## Reuses WITHOUT MODIFICATION:
##   .pf_col_meta()             (R/pf-newton.R)
##   .pf_build_sparse_system()  (R/pf-newton.R)
##   .pf_make_dy_exo()          (R/perfect-foresight-solve.R)
##   pf_newton_solve()          (R/pf-newton.R)
##   compute_irfs()             (R/stochsimul-monolith.R)
## --------------------------------------------------------------------------


# =============================================================================
# Constructor: build the steady-state sequence Jacobian
# =============================================================================

#' Build a representative-agent sequence-space Jacobian (RA-SSJ) object
#'
#' Assembles the T*n_endo x T*n_endo stacked Jacobian \eqn{J = \partial H /
#' \partial X} evaluated at the deterministic steady state and the shock-impact
#' vector \eqn{\partial H / \partial \varepsilon_1} at t = 1.  Both are
#' extracted from the same compiled-derivative machinery used by
#' \code{\link{pf_newton_solve}}.
#'
#' The resulting object supports one operation, \code{\link{sj_irf}}, which
#' computes the general-equilibrium impulse response
#' \eqn{dX = -J^{-1} (\partial H/\partial\varepsilon_1) \cdot \text{scale}}.
#'
#' @param model     A parsed dynhr model (\code{\link{parse_mod}}).
#' @param compiled  A compiled model (\code{\link{compile_model}}).
#' @param ss        Named numeric vector: deterministic steady state
#'   (from \code{\link{solve_steady}}).
#' @param params    Named numeric parameter vector (default: \code{model$param_values}).
#' @param horizon   Integer: sequence horizon T (default 40).
#'
#' @return An object of class \code{ra_ssj} with components:
#'   \describe{
#'     \item{\code{J}}{Sparse dgCMatrix (T*n_endo x T*n_endo): dH/dX at ss.}
#'     \item{\code{dH_deps}}{Named list, one entry per shock: each a numeric
#'       vector of length T*n_endo giving dH/deps at t=1.}
#'     \item{\code{J_lu}}{LU factorisation of J (via \code{Matrix::lu}).}
#'     \item{\code{T}}{Integer: horizon.}
#'     \item{\code{n_endo}}{Integer: number of endogenous variables.}
#'     \item{\code{endo_names}}{Character: endogenous variable names.}
#'     \item{\code{exo_names}}{Character: exogenous variable names.}
#'   }
#'
#' @details
#' \strong{Sign/transpose/ordering conventions}
#'
#' The PF Newton system solves \eqn{H(X) = 0} where the stacked residual is
#' indexed as \eqn{H_{(t-1) n + i}} = residual of equation \eqn{i} at period
#' \eqn{t}.  The stacked state vector \eqn{X} is ordered identically.
#'
#' The steady-state Jacobian is built by evaluating
#' \code{.pf_build_sparse_system} at the all-ss path with zero shocks.  At the
#' steady state \eqn{H(X_{ss}) = 0} (by definition), so the assembled residual
#' is zero; only the Jacobian is used.
#'
#' The shock-impact vector \code{dH_deps} is extracted by calling
#' \code{jacobian_fn} at t = 1 (at the ss path) and reading the columns
#' corresponding to exogenous variables.  For periods t > 1 the direct shock
#' impact is zero (the shock enters future periods only through lagged state
#' transitions, which are already captured in \eqn{J}).
#'
#' The GE impulse response formula \eqn{dX = -J^{-1} \partial H/\partial\varepsilon}
#' follows from total differentiation of \eqn{H(X, \varepsilon) = 0}:
#' \eqn{J \, dX + (\partial H/\partial\varepsilon) \, d\varepsilon = 0}.
#'
#' @seealso \code{\link{sj_irf}}
#' @export
ra_ssj <- function(model, compiled, ss, params = NULL, horizon = 40L) {
  if (is.null(params)) params <- model$param_values

  dyn    <- compiled$dynamic
  n_endo <- length(dyn$endo_names)
  n_exo  <- length(dyn$exo_names)
  T      <- as.integer(horizon)

  # Align ss to endo_names order (same convention as pf_newton_solve)
  y_ss_num <- as.numeric(ss[dyn$endo_names])

  # Build a T x n_endo path at the steady state (all periods = ss)
  Y_ss <- matrix(rep(y_ss_num, T), nrow = T, byrow = TRUE)
  colnames(Y_ss) <- dyn$endo_names

  # Zero shock path
  eps_mat <- matrix(0, nrow = T, ncol = n_exo)
  colnames(eps_mat) <- dyn$exo_names

  # Precompute column metadata (shared with pf_newton_solve)
  meta <- .pf_col_meta(dyn)

  # Preallocate triplet buffer
  trip_buf <- .pf_preallocate_triplets(dyn, T)

  # Build the stacked Jacobian J = dH/dX at the ss path with zero shocks.
  # The named steady-state vector (for STEADY_STATE() references)
  y_ss_named <- setNames(y_ss_num, dyn$endo_names)

  sys <- .pf_build_sparse_system(
    Y        = Y_ss,
    y0_num   = y_ss_num,
    y_ss_num = y_ss_num,
    eps_mat  = eps_mat,
    meta     = meta,
    dyn      = dyn,
    params   = params,
    y_ss     = y_ss_named,
    T        = T,
    n_eq     = n_endo,
    n_endo   = n_endo,
    obc_specs    = list(),
    regime       = matrix(FALSE, nrow = 1L, ncol = T),
    spec_cur_dc  = integer(0),
    trip_buf     = trip_buf
  )

  J <- sys$J  # T*n_endo x T*n_endo sparse Jacobian

  # -- Build dH/deps for each shock: T*n_endo vector --
  # At t=1, evaluate jacobian_fn at the ss dy vector with zero shocks.
  # The exo columns of Jt[, dc] give dH_t / deps (rows = equations at t=1;
  # col = dynamic column for the exo variable).
  # For t > 1, the direct shock impact is zero.
  dy1 <- .pf_make_dy_exo(meta, Y_ss, y_ss_num, y_ss_num, eps_mat, 1L, T, n_exo)
  Jt1 <- dyn$jacobian_fn(dy1, params, y_ss_named)

  # dH_deps: for each shock, a T*n_endo vector (nonzero only at rows 1..n_endo)
  dH_deps <- vector("list", n_exo)
  names(dH_deps) <- dyn$exo_names

  for (k in seq_len(n_exo)) {
    exo_nm  <- dyn$exo_names[k]
    vec     <- numeric(T * n_endo)  # zero everywhere

    # Find all dynamic columns for this exo variable (at any lead/lag)
    # In the standard linear RBC, eps_a enters at lag 0 only.
    # We sum over all exo dynamic columns for this variable.
    for (j in seq_along(meta$kind)) {
      if (meta$kind[j] != "exo") next
      if (meta$var_idx[j] != k)  next
      dc <- meta$dyn_col[j]
      # Add partial of equations at t=1 w.r.t. this exo column
      vec[seq_len(n_endo)] <- vec[seq_len(n_endo)] + Jt1[, dc]
    }

    dH_deps[[exo_nm]] <- vec
  }

  # LU factorisation of J (reused across sj_irf calls)
  J_lu <- Matrix::lu(J)

  structure(
    list(
      J          = J,
      dH_deps    = dH_deps,
      J_lu       = J_lu,
      T          = T,
      n_endo     = n_endo,
      endo_names = dyn$endo_names,
      exo_names  = dyn$exo_names
    ),
    class = "ra_ssj"
  )
}


# =============================================================================
# sj_irf: compute the GE impulse response via SSJ inversion
# =============================================================================

#' Compute the general-equilibrium impulse response via sequence-space Jacobian
#'
#' Given an \code{\link{ra_ssj}} object, computes the linear impulse response
#' \deqn{dX = -J^{-1} \frac{\partial H}{\partial \varepsilon_1} \cdot \text{scale}}
#' and returns it as a T x n_endo matrix in the same format as
#' \code{\link{compute_irfs}} (rows = periods, cols = endogenous variables).
#'
#' @param ssj    An \code{ra_ssj} object from \code{\link{ra_ssj}}.
#' @param shock  Character: name of the shock (must be in \code{ssj$exo_names}).
#' @param scale  Numeric scalar: shock size (default 1.0).  Set to
#'   \code{model$param_values["sigma_a"]} to match the one-standard-deviation
#'   IRF convention used by \code{compute_irfs}.
#'
#' @return A T x n_endo numeric matrix (rows = periods, cols = endo vars),
#'   with column names equal to \code{ssj$endo_names}.  Compatible with
#'   the slot format of an \code{IRFCollection} object returned by
#'   \code{compute_irfs}.
#'
#' @details
#' The GE response formula follows from total differentiation of
#' \eqn{H(X, \varepsilon) = 0} at the steady state:
#' \deqn{J \, dX + \frac{\partial H}{\partial \varepsilon} \, d\varepsilon = 0
#'       \;\Rightarrow\; dX = -J^{-1} \frac{\partial H}{\partial \varepsilon} \, d\varepsilon.}
#'
#' @seealso \code{\link{ra_ssj}}, \code{\link{compute_irfs}}
#' @export
sj_irf <- function(ssj, shock, scale = 1.0) {
  stopifnot(inherits(ssj, "ra_ssj"))
  if (!shock %in% ssj$exo_names)
    stop(sprintf("sj_irf: shock '%s' not in exo_names (%s)",
                 shock, paste(ssj$exo_names, collapse = ", ")))

  dH <- ssj$dH_deps[[shock]] * scale  # T*n_endo vector

  # Solve J * dX = -dH  => dX = -J^{-1} dH
  dX_vec <- as.numeric(Matrix::solve(ssj$J_lu, -dH))

  # Reshape to T x n_endo (row t = period t deviation from ss).
  # dX_vec is indexed: dX_vec[(t-1)*n_endo + i] = deviation of var i at period t.
  # matrix(..., byrow=TRUE) fills row-major: row t = dX_vec[(t-1)*n_endo + 1..n_endo].
  T      <- ssj$T
  n_endo <- ssj$n_endo
  dX <- matrix(dX_vec, nrow = T, ncol = n_endo, byrow = TRUE)
  colnames(dX) <- ssj$endo_names
  dX
}


# =============================================================================
# print method
# =============================================================================

#' @export
print.ra_ssj <- function(x, ...) {
  cat(sprintf(
    "<ra_ssj>  horizon T=%d  n_endo=%d  n_exo=%d\n",
    x$T, x$n_endo, length(x$exo_names)
  ))
  cat("  endo:", paste(x$endo_names, collapse = ", "), "\n")
  cat("  exo: ", paste(x$exo_names,  collapse = ", "), "\n")
  invisible(x)
}


# =============================================================================
# Behavioural (non-FIRE) expectations: Lenney & Rosso (2026) forecast mapping
# =============================================================================
##
## Source: Lenney, J. and Rosso, B. (2026), "A flexible deviation from FIRE in
## the sequence space", Bank of England Staff Working Paper No. 1197.
##
## Expectations process (their eq. 1-6, the Kohlhas-Walther reduced form):
## a fraction (1 - theta) of agents updates its information set each period
## (theta = delta / (1 + delta) in their news-revision parameterisation), and
## an updating agent extrapolates from the input it observes at the update,
## with coefficient gamma (gamma < 0 = overreaction to current conditions).
## The average time-h forecast of the input path, dX^{e,h} = A_h dX, is
## (their eq. 8 / A.7; 0-based dates):
##
##   A_h[s, r] = 1                               s <= h, r = s  (observed)
##             = 1 - theta^(h+1)                 s >  h, r = s  (news, B_h)
##             = -gamma (1-theta) theta^(h-r)    s >  h, r <= h (extrapolation)
##             = 0                               otherwise
##
## theta = gamma = 0 gives A_h = I for every h (FIRE); gamma = 0 alone is the
## sticky-expectations mapping of Auclert, Rognlie & Straub (2020), with
## agents observing the current input (the identity block of A_h).
##
## Behavioural block Jacobian (their eq. 7 / A.10): with F the FIRE fake-news
## matrix (F[0,.] = J[0,.], F[.,0] = J[.,0], F[t,s] = J[t,s] - J[t-1,s-1]),
##
##   Jbar = sum_{h >= 0} P_h A_h,  P_h[t, s] = F[t-h, s-h] for t, s >= h, else 0,
##
## equivalently (their A.4, the Bardoczy-Guerreiro form)
##   Jbar = J A_0 + sum_{h >= 1} R_h (A_h - A_{h-1}),  R_h = J shifted by h.
## The transform is block-local and per input: each input column-block of a
## block's Jacobian is transformed with that input's (theta, gamma), and GE
## assembly is unchanged (their section 2.3). Constant (theta, gamma) only;
## the paper's horizon-varying generalisation (their App. A.5) and the
## extrapolation decay of their footnote 22 are not implemented.

#' Behavioural forecast matrix A_h (0-based update date h)
#' @noRd
.ssj_forecast_matrix <- function(T_h, h, theta, gamma = 0) {
  A <- diag(T_h)
  if (h + 1L < T_h) {
    fut <- (h + 2L):T_h                       # 1-based rows s = h+1 .. T-1
    A[cbind(fut, fut)] <- 1 - theta^(h + 1)
    ext <- -gamma * (1 - theta) * theta^(h - 0:h)
    A[fut, seq_len(h + 1L)] <- matrix(ext, length(fut), h + 1L, byrow = TRUE)
  }
  A
}

#' Average behavioural k-step-ahead forecast path f_t x_{t+k}, t = 0..T-1-k
#'
#' The forecast held at date t (after that date's update) of the input k
#' periods ahead, as a deviation from steady state: row t+k of A_t dx.
#' @noRd
.ssj_behavioural_forecast <- function(dx, k, theta, gamma = 0) {
  .ssj_check_expectation_pars(theta, gamma, "forecast")
  T_h <- length(dx)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 0L || k >= T_h)
    .dynhr_abort("behavioural forecast: `k` must be a single integer in ",
                 "[0, length(dx) - 1].",
                 class = "dynhr_error_behavioural_expectations")
  vapply(0:(T_h - 1L - k), function(t) {
    if (k == 0L) return(dx[t + 1L])
    (1 - theta^(t + 1)) * dx[t + k + 1L] -
      gamma * (1 - theta) * sum(theta^(t - 0:t) * dx[seq_len(t + 1L)])
  }, numeric(1))
}

#' Validate one (theta, gamma) pair
#' @noRd
.ssj_check_expectation_pars <- function(theta, gamma, what) {
  if (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta) ||
      theta < 0 || theta >= 1)
    .dynhr_abort("behavioural expectations (", what, "): `theta` must be a ",
                 "single number in [0, 1) (0 = every agent updates each ",
                 "period).", class = "dynhr_error_behavioural_expectations")
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma))
    .dynhr_abort("behavioural expectations (", what, "): `gamma` must be a ",
                 "single finite number.",
                 class = "dynhr_error_behavioural_expectations")
  invisible(TRUE)
}

#' Behavioural Jacobian of one T x T FIRE Jacobian (Lenney-Rosso eq. 7)
#'
#' theta = gamma = 0 returns J itself (bit-identical).
#' @noRd
.ssj_behavioural_jacobian <- function(J, theta, gamma = 0) {
  .ssj_check_expectation_pars(theta, gamma, "jacobian")
  if (!is.matrix(J) || nrow(J) != ncol(J))
    .dynhr_abort("behavioural expectations: the Jacobian must be a square ",
                 "T x T matrix.", class = "dynhr_error_behavioural_expectations")
  if (theta == 0 && gamma == 0) return(J)
  T_h <- nrow(J)
  ## FIRE fake-news matrix (ABRS 2021): F[t,s] = J[t,s] - J[t-1,s-1].
  Fm <- J
  if (T_h >= 2L) Fm[-1L, -1L] <- J[-1L, -1L] - J[-T_h, -T_h]
  JB <- matrix(0, T_h, T_h)
  for (h in 0:(T_h - 1L)) {
    n    <- T_h - h
    rows <- (h + 1L):T_h
    Fs   <- Fm[seq_len(n), seq_len(n), drop = FALSE]
    ## Row h of A_h is e_h: the input observed at the update date h.
    JB[rows, h + 1L] <- JB[rows, h + 1L] + Fs[, 1L]
    if (n >= 2L) {
      fut  <- Fs[, -1L, drop = FALSE]
      cols <- (h + 2L):T_h
      ## News block B_h: diagonal 1 - theta^(h+1) on the future inputs.
      JB[rows, cols] <- JB[rows, cols] + (1 - theta^(h + 1)) * fut
      ## Extrapolation block: every future forecast shifts by
      ## -gamma (1-theta) theta^(h-r) dx_r for each observed r <= h.
      if (gamma != 0) {
        ext  <- -gamma * (1 - theta) * theta^(h - 0:h)
        past <- seq_len(h + 1L)
        JB[rows, past] <- JB[rows, past] + outer(rowSums(fut), ext)
      }
    }
  }
  dimnames(JB) <- dimnames(J)
  JB
}

#' Resolve an `expectations` specification to per-input (theta, gamma)
#'
#' Returns NULL for FIRE (NULL spec, type 'fire', or every input at
#' theta = gamma = 0), otherwise list(theta =, gamma =), named numerics over
#' `inputs`.
#' @noRd
.ssj_expectations_resolve <- function(expectations, inputs, where) {
  if (is.null(expectations)) return(NULL)
  fail <- function(...) .dynhr_abort(where, ": ", ...,
                                     class = "dynhr_error_behavioural_expectations")
  if (!is.list(expectations) || is.null(names(expectations)) ||
      any(!nzchar(names(expectations))))
    fail("`expectations` must be NULL (FIRE) or a named list with elements ",
         "`theta` and optionally `gamma` / `type`.")
  bad <- setdiff(names(expectations), c("type", "theta", "gamma"))
  if (length(bad))
    fail("unknown `expectations` element(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; allowed: 'type', 'theta', 'gamma'.")
  type <- expectations$type
  if (is.null(type)) type <- "behavioural"
  if (!is.character(type) || length(type) != 1L ||
      !type %in% c("behavioural", "sticky", "fire"))
    fail("`expectations$type` must be one of 'behavioural', 'sticky', 'fire'.")
  if (type == "fire") {
    if (!is.null(expectations$theta) || !is.null(expectations$gamma))
      fail("`expectations$type = 'fire'` takes no `theta` / `gamma`.")
    return(NULL)
  }
  if (is.null(expectations$theta))
    fail("`expectations` needs `theta` (the stickiness, in [0, 1)).")
  if (type == "sticky" && !is.null(expectations$gamma) &&
      any(expectations$gamma != 0))
    fail("`expectations$type = 'sticky'` is the gamma = 0 case; use ",
         "type = 'behavioural' for extrapolation.")
  per_input <- function(v, nm) {
    out <- setNames(rep(0, length(inputs)), inputs)
    if (is.null(v)) return(out)
    if (!is.numeric(v)) fail("`expectations$", nm, "` must be numeric.")
    if (is.null(names(v))) {
      if (length(v) != 1L)
        fail("an unnamed `expectations$", nm, "` must be a scalar (applied ",
             "to every input); name the entries by input to set them per ",
             "input.")
      out[] <- as.numeric(v)
      return(out)
    }
    if (length(setdiff(names(v), inputs)) || anyDuplicated(names(v)))
      fail("`expectations$", nm, "` names must be distinct inputs of this ",
           "block (", paste0("'", inputs, "'", collapse = ", "), "); got ",
           paste0("'", names(v), "'", collapse = ", "), ".")
    out[names(v)] <- as.numeric(v)
    out
  }
  theta <- per_input(expectations$theta, "theta")
  gamma <- per_input(expectations$gamma, "gamma")
  for (i in inputs)
    .ssj_check_expectation_pars(theta[[i]], gamma[[i]],
                                paste0(where, ", input '", i, "'"))
  if (all(theta == 0 & gamma == 0)) return(NULL)
  list(theta = theta, gamma = gamma)
}

#' Apply behavioural expectations to a nested Jacobian J[[output]][[input]]
#'
#' Inputs not named in a per-input spec keep their FIRE columns.
#' @noRd
.ssj_apply_expectations <- function(J, expectations, inputs, where) {
  ex <- .ssj_expectations_resolve(expectations, inputs, where)
  if (is.null(ex)) return(J)
  for (o in names(J)) for (i in intersect(names(J[[o]]), inputs)) {
    if (ex$theta[[i]] == 0 && ex$gamma[[i]] == 0) next
    J[[o]][[i]] <- .ssj_behavioural_jacobian(J[[o]][[i]], ex$theta[[i]],
                                             ex$gamma[[i]])
  }
  J
}

#' TRUE when a block carries non-FIRE expectations
#' @noRd
.ssj_block_is_behavioural <- function(blk) {
  !is.null(blk$expectations) &&
    !is.null(.ssj_expectations_resolve(blk$expectations, blk$inputs,
                                       paste0("block '", blk$name, "'")))
}
