## R/perfect-foresight-solve.R
## --------------------------------------------------------------------------
## Perfect-foresight stacked Newton path solver for non-stationary BGP models.
##
## Provides:
##   pf_boundary_exo()          -- extract exo init/terminal levels from
##                                 a parsed model's initval/endval block
##   perfect_foresight_solve()  -- stacked Newton solver for deterministic
##                                 transition paths without steady state
##
## ALGORITHM
##
##   Solves the stacked T-period system:
##     F_i(y_{t-1}, y_t, y_{t+1}, ε_t) = 0   t = 1..T, i = 1..n_endo
##
##   Boundary conditions:
##     y_0 = y_init (given, from initval block)
##     y_{T+1} = y_terminal (given, from endval block)
##
##   The stacked Jacobian is block-tridiagonal, assembled as a sparse
##   dgCMatrix and solved via Matrix::sparseQR.
##
##   This handles the case where Dynare's perfect_foresight_solver is used
##   for non-stationary BGP models (e.g., Ramsey_Cass_Koopmans, Solow).
##   Unlike mcp_solve_path(), there are no MCP/Fischer-Burmeister
##   modifications — just plain dynamic residuals.
##
## REFERENCES
##   Adjemian & Juillard (2025) — Stochastic Extended Path
##   Dynare perfect_foresight_solver (stack_solve_algo=0/7)
## --------------------------------------------------------------------------


# =============================================================================
# pf_boundary_exo: extract exo init / terminal values from parsed model blocks
# =============================================================================

#' Extract exogenous boundary values from a parsed model's initval or endval block
#'
#' A parsed model returned by \code{\link{parse_mod}} stores the evaluated
#' \code{initval} and \code{endval} blocks as named numeric vectors over
#' \emph{all} model variables (endo + exo).  This helper filters those vectors
#' down to the exogenous variables and returns a clean named numeric vector
#' suitable for passing as \code{exo_init} or \code{exo_terminal} to
#' \code{\link{perfect_foresight_solve}}.
#'
#' @section Named-vector footgun:
#' Expression values derived from parameters (e.g. \code{p["g"]}) carry
#' stray names, so \code{c(A = p["g"])} produces \code{A.g}, not \code{A}.
#' This function guards against that by calling
#' \code{setNames(unname(values), exo_names)}, ensuring the result's names are
#' exactly \code{model$varexo_names} with no dotted suffixes.
#'
#' @param model  A parsed dynhr model, i.e. the list returned by
#'   \code{\link{parse_mod}}.
#' @param which  Character scalar: \code{"init"} (default) to extract the
#'   period-0 exogenous values from the \code{initval} block (use as
#'   \code{exo_init} in \code{perfect_foresight_solve}), or
#'   \code{"terminal"} to extract the post-terminal (t=T+1) values from the
#'   \code{endval} block (use as \code{exo_terminal}).
#'
#' @return A named numeric vector of length \code{length(model$varexo_names)}
#'   with names exactly equal to \code{model$varexo_names}, or \code{NULL}
#'   if the relevant block is absent or contains no exogenous entries.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' m  <- parse_mod("Solow_nonstationary_pp.mod", verbose = FALSE)
#' cm <- compile_model(m, verbose = FALSE)
#' T  <- 100L
#' exo_path <- cbind(A = (1 + m$param_values["g"])^(1:T),
#'                   L = (1 + m$param_values["n"])^(1:T))
#' result <- perfect_foresight_solve(
#'   cm, y0, y_terminal, exo_path, m$param_values,
#'   n_periods    = T,
#'   exo_init     = pf_boundary_exo(m, "init"),
#'   exo_terminal = pf_boundary_exo(m, "terminal")
#' )
#' }
pf_boundary_exo <- function(model, which = c("init", "terminal")) {
  which <- match.arg(which)

  exo_names <- model$varexo_names
  if (is.null(exo_names) || length(exo_names) == 0L)
    return(NULL)

  block <- if (which == "init") model$initval else model$endval

  # Block absent or empty
  if (is.null(block) || length(block) == 0L)
    return(NULL)

  # Filter to exogenous entries only
  present <- intersect(exo_names, names(block))
  if (length(present) == 0L)
    return(NULL)

  # Build a full-length vector, NA for any exo not in block
  raw_vals <- block[exo_names]        # named subset (NA where absent)

  # Guard: force names to be exactly exo_names, not dotted compound names
  # that arise when block values carry stray names from parameter look-ups
  # (e.g. `A = p["g"]` produces name "A.g", silently breaking vec[exo_names]).
  setNames(unname(raw_vals), exo_names)
}


# =============================================================================
# Extended dy builder with exogenous variable lead/lag support
# =============================================================================

#' Build named dy vector for one period, with exogenous variable lead/lag support
#'
#' The compiled residual/Jacobian functions may reference exogenous variables
#' at leads or lags (e.g., dy["A__m1"], dy["L__m1"]) even though the column map
#' only lists them at lag 0. This function extends .pf_make_dy() to detect and
#' fill those additional entries by shifting the eps_mat rows.
#'
#' @param meta       Output of .pf_col_meta()
#' @param Y          T × n_endo numeric matrix (current path)
#' @param y0_num     Length-n_endo numeric: initial state
#' @param y_term_num Length-n_endo numeric: terminal condition
#' @param eps_mat    T × n_exo numeric matrix: exogenous path
#' @param t          Period index (1-based)
#' @param T          Total horizon
#' @param n_exo      Number of exogenous variables
#' @param exo_init   Optional length-n_exo numeric: period-0 values for
#'   exogenous variables (from initval block). If NULL, uses eps_mat row 1.
#' @return Named numeric vector of length covering all columns the residuals
#'         function might reference (endo at all leads/lags + exo at lag -1..+1)
#' @noRd
.pf_make_dy_exo <- function(meta, Y, y0_num, y_term_num, eps_mat, t, T, n_exo,
                             exo_init = NULL, exo_terminal = NULL) {
  # Start with the standard dy (endo + exo at current period)
  dy <- .pf_make_dy(meta, Y, y0_num, y_term_num, eps_mat[t, ], t, T)

  # Check what exogenous variable entries the residuals function might need
  # at non-zero leads/lags. The column map only has exo at lag 0, but the
  # compiled function may reference exo at lag -1 or +1.
  # We detect this by checking if the function tries to access names we
  # haven't set, but since R returns NA for missing named elements,
  # we proactively fill all reasonable exo lead/lag combos.

  exo_names <- colnames(eps_mat)
  if (is.null(exo_names) || length(exo_names) == 0) return(dy)

  # Build the initial exo values: prefer exo_init if provided, else use eps_mat row 1
  exo_init_vec <- if (!is.null(exo_init)) {
    as.numeric(exo_init[exo_names])
  } else {
    eps_mat[1L, ]
  }
  if (anyNA(exo_init_vec)) exo_init_vec <- eps_mat[1L, ]

  # Build the post-terminal exo values A(T+1) for exo LEAD references at t=T:
  # prefer exo_terminal (e.g. the endval exogenous level), else fall back to the
  # last row of eps_mat. The fallback is WRONG for a non-constant terminal exo
  # (a growing deterministic path, e.g. Solow's A=(1+g)^t, where A(T+1) differs
  # from A(T)=eps_mat[T,]); such callers must supply exo_terminal.
  exo_term_vec <- if (!is.null(exo_terminal)) {
    as.numeric(exo_terminal[exo_names])
  } else {
    eps_mat[T, ]
  }
  if (anyNA(exo_term_vec)) exo_term_vec <- eps_mat[T, ]

  for (ei in seq_along(exo_names)) {
    enm <- exo_names[ei]
    for (ll in c(-1L, 1L)) {
      ll_key <- paste0(enm, "__", if (ll < 0) paste0("m", abs(ll)) else paste0("p", ll))
      # Check if this key already exists in dy (unlikely for exo, but be safe)
      if (ll_key %in% names(dy)) next

      # Compute the value from eps_mat at shifted period
      # lead_lag = -1 (lag): need t-1's value  → t + (-1) = t - 1
      # lead_lag = +1 (lead): need t+1's value → t + (+1) = t + 1
      shifted_t <- t + ll
      val <- if (shifted_t < 1L) {
        # Before first period: use the period-0 exogenous values
        # (from the initval block, e.g. A₀=1, L₀=1)
        exo_init_vec[ei]
      } else if (shifted_t > T) {
        # After last period: use the post-terminal exogenous level A(T+1),
        # NOT eps_mat[T,] (= A(T)) — see exo_term_vec above.
        exo_term_vec[ei]
      } else {
        eps_mat[shifted_t, ei]
      }

      # Append to dy vector
      dy[[ll_key]] <- val
    }
  }

  dy
}


# =============================================================================
# Build the stacked system (no MCP modifications)
# =============================================================================

#' Build the stacked perfect-foresight system with sparse block-tridiagonal Jacobian
#'
#' For each period t:
#'   1. Evaluate the full dynamic residual F_t and Jacobian J_t.
#'   2. Route Jacobian columns into the sparse block-tridiagonal matrix.
#'   3. Assemble into (T·n_endo) × (T·n_endo) dgCMatrix.
#'
#' @param Y          T × n_endo numeric matrix: current path
#' @param y0_num     Length-n_endo initial state at t=0 (initval)
#' @param y_term_num Length-n_endo terminal condition at t=T+1 (endval)
#' @param eps_mat    T × n_exo shock/exogenous variable matrix
#' @param meta       Column metadata from .pf_col_meta()
#' @param compiled   dynhr_compiled (for residuals_fn, jacobian_fn)
#' @param params     Named numeric parameter vector
#' @param y_ss_num   Length-n_endo numeric: steady-state / terminal values passed
#'   as the third argument to \code{residuals_fn} and \code{jacobian_fn} (the
#'   \code{ss} / \code{STEADY_STATE()} oracle). Must be the reordered terminal
#'   vector (\code{y_term_num}) rather than \code{NULL}: some compiled models
#'   reference \code{ss["varname"]} directly in the generated code (issue I11).
#' @param T, n_endo  Integer: horizon and number of endogenous variables
#' @param exo_init   Optional length-n_exo numeric: period-0 values for
#'   exogenous variables (from initval). Used for correct lagged exo values.
#' @return List with:
#'   $J  — sparse dgCMatrix (T*n_endo × T*n_endo)
#'   $R  — numeric vector (T*n_endo): stacked residuals
#' @noRd
.pf_build_stacked_system <- function(Y, y0_num, y_term_num, eps_mat,
                                      meta, compiled, params, T, n_endo,
                                      exo_init = NULL, exo_terminal = NULL,
                                      y_ss_num = NULL) {
  dyn     <- compiled$dynamic
  n_total <- T * n_endo
  n_exo   <- length(dyn$exo_names)

  # Build the named ss vector passed to residuals_fn / jacobian_fn as the
  # STEADY_STATE() oracle.  Some compiled models reference ss["varname"]
  # directly (e.g. `dy["yhat__0"] - ss["y"]` for the output gap definition).
  # Passing NULL causes ss["y"] → NULL → numeric(0) → "replacement has length
  # zero" (issue I11).  Use y_term_num if supplied, else y_term_num itself
  # (the terminal values are the natural steady-state proxy).
  y_ss_oracle <- if (!is.null(y_ss_num)) {
    setNames(y_ss_num, dyn$endo_names)
  } else {
    setNames(y_term_num, dyn$endo_names)
  }

  # Preallocate residual
  R <- numeric(n_total)

  # Preallocate sparse triplet storage using jac_triplets NNZ budget.
  # Falls back to growing vectors when jac_triplets is NULL (Landmine 3).
  jt <- dyn$jac_triplets
  if (!is.null(jt) && length(jt) > 0L) {
    nnz_est   <- T * length(jt)
    i_triplet <- integer(nnz_est)
    j_triplet <- integer(nnz_est)
    v_triplet <- numeric(nnz_est)
    pre_alloc <- TRUE
    ptr       <- 0L
  } else {
    i_triplet <- integer(0)
    j_triplet <- integer(0)
    v_triplet <- numeric(0)
    pre_alloc <- FALSE
  }

  for (t in seq_len(T)) {
    row_off <- (t - 1L) * n_endo
    col_off <- (t - 1L) * n_endo

    # Build dy vector for this period
    # Extended version of .pf_make_dy that handles exogenous variable
    # lags/leads properly. The compiled residual function may reference
    # dy entries like A__m1 or L__m1 even though the column map only
    # has A__0 and L__0. We build the dy here exhaustively.
    dy <- .pf_make_dy_exo(meta, Y, y0_num, y_term_num, eps_mat, t, T, n_exo,
                           exo_init = exo_init, exo_terminal = exo_terminal)

    # Compute full residual and Jacobian.
    # Pass y_ss_oracle (not NULL) so that STEADY_STATE() references in the
    # generated code (e.g. ss["y"]) resolve correctly (issue I11).
    Rt <- dyn$residuals_fn(dy, params, y_ss_oracle)
    Jt <- dyn$jacobian_fn(dy, params, y_ss_oracle)

    # Store residual
    R[row_off + seq_len(n_endo)] <- Rt

    # Route Jacobian columns to global sparse block-tridiagonal entries
    # Uses pf_col_meta: kind = "lag"/"cur"/"lead"/"exo"
    for (k in seq_along(meta$kind)) {
      kd <- meta$kind[k]
      if (kd == "exo") next

      dc <- meta$dyn_col[k]
      vi <- meta$var_idx[k]

      # Map to global column in the stacked system
      gcol <- switch(kd,
        lag  = if (t > 1L) (t - 2L) * n_endo + vi else NA_integer_,
        cur  = col_off + vi,
        lead = if (t < T)  t * n_endo + vi         else NA_integer_
      )
      if (is.na(gcol)) next

      col_vals <- Jt[, dc]

      if (pre_alloc) {
        for (eq in seq_len(n_endo)) {
          val <- col_vals[eq]
          if (isTRUE(val != 0) && is.finite(val)) {
            ptr <- ptr + 1L
            i_triplet[ptr] <- row_off + eq
            j_triplet[ptr] <- gcol
            v_triplet[ptr] <- val
          }
        }
      } else {
        for (eq in seq_len(n_endo)) {
          val <- col_vals[eq]
          if (isTRUE(val != 0) && is.finite(val)) {
            i_triplet <- c(i_triplet, row_off + eq)
            j_triplet <- c(j_triplet, gcol)
            v_triplet <- c(v_triplet, val)
          }
        }
      }
    }
  }

  # Trim preallocated buffers to actual NNZ
  if (pre_alloc && ptr > 0L) {
    i_triplet <- i_triplet[seq_len(ptr)]
    j_triplet <- j_triplet[seq_len(ptr)]
    v_triplet <- v_triplet[seq_len(ptr)]
  }

  # Build sparse matrix from triplets — repr = "C" replaces deprecated giveCsparse = TRUE (Landmine 5)
  if (length(i_triplet) == 0L) {
    J <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 0,
                              dims = c(n_total, n_total),
                              repr = "C")
  } else {
    J <- Matrix::sparseMatrix(
      i = i_triplet, j = j_triplet, x = v_triplet,
      dims = c(n_total, n_total),
      index1 = TRUE,
      repr = "C"
    )
  }

  list(J = J, R = R)
}


# =============================================================================
# Backtracking line search
# =============================================================================

#' Armijo backtracking line search for perfect-foresight merit function
#'
#' Finds step length α ∈ (0, 1] such that:
#'   θ(Y + α·Δ) ≤ θ(Y) + σ·α·∇θ(Y)·Δ
#' where θ(Y) = ½||R(Y)||².
#'
#' @param Y            T × n_endo matrix: current path
#' @param delta_vec    Numeric vector (length T*n_endo): Newton direction
#' @param R_stack      Numeric vector (length T*n_endo): current residual
#' @param J_stack      dgCMatrix: current Jacobian
#' @param theta_cur    Scalar: current merit value = ½||R||²
#' @param fn_merit     Function to evaluate θ at a new Y
#' @param sigma        Armijo parameter (default 1e-4)
#' @param max_ls       Maximum line search iterations (default 20)
#' @param n_endo       Integer: number of endogenous variables
#' @param T            Integer: horizon
#' @return List with:
#'   $alpha      — step length
#'   $theta_new  — merit value at Y + alpha*delta
#'   $Y_new      — updated path matrix (or NULL if no tried step reduced merit)
#'   $ls_iter    — iterations used
#'   $accepted   — logical: TRUE if the Armijo condition was met
#'
#' @details
#' On Armijo failure the search no longer returns \code{NULL} with a discarded
#' step. Instead it tracks the smallest-merit trial point actually evaluated and
#' returns it whenever that point strictly improves on the incumbent. This is
#' critical for cold-start robustness: the previous behaviour (return NULL →
#' caller applies a blind half Newton step) could send the path to ~1e9 in a
#' single iteration and singularise the next Jacobian (issue M15). Returning the
#' best vetted point guarantees monotone non-increase of the merit function.
#' @noRd
.pf_line_search <- function(Y, delta_vec, R_stack, J_stack,
                             theta_cur, fn_merit,
                             sigma = 1e-4, max_ls = 20L,
                             n_endo, T) {
  # Compute gradient: ∇θ = J' · R
  grad <- as.numeric(Matrix::crossprod(J_stack, R_stack))
  directional_deriv <- sum(grad * delta_vec)

  # If directional derivative is positive, Newton direction is not a descent
  # direction. Fall back to steepest descent: Δ = -∇θ
  if (directional_deriv >= 0) {
    delta_vec <- -grad
    directional_deriv <- -sum(grad * grad)
  }

  alpha <- 1.0

  # Track the best (lowest-merit) trial point that strictly improves on the
  # incumbent, so an Armijo failure still yields a usable (merit-decreasing)
  # step rather than discarding all work (M15).
  best_theta <- theta_cur
  best_Y     <- NULL
  best_alpha <- 0

  for (ls_iter in seq_len(max_ls)) {
    # Trial point
    Y_trial <- Y
    for (t in seq_len(T)) {
      idx_t <- (t - 1L) * n_endo + seq_len(n_endo)
      Y_trial[t, ] <- Y[t, ] + alpha * delta_vec[idx_t]
    }

    theta_new <- fn_merit(Y_trial)

    if (is.finite(theta_new) && theta_new < best_theta) {
      best_theta <- theta_new
      best_Y     <- Y_trial
      best_alpha <- alpha
    }

    # Armijo condition (isTRUE guards NaN merit values from invalid steps).
    if (isTRUE(theta_new <= theta_cur + sigma * alpha * directional_deriv)) {
      return(list(alpha = alpha, theta_new = theta_new,
                  Y_new = Y_trial, ls_iter = ls_iter, accepted = TRUE))
    }

    alpha <- alpha * 0.5
  }

  # Armijo never satisfied: fall back to the best merit-decreasing trial point
  # seen during backtracking (may be NULL if none improved — caller decides).
  list(alpha = best_alpha, theta_new = best_theta,
       Y_new = best_Y, ls_iter = max_ls, accepted = FALSE)
}


# =============================================================================
# Robust Newton-step solver (Levenberg-Marquardt fallback)
# =============================================================================

#' Solve J·δ = -R with a Levenberg-Marquardt fallback for singular Jacobians
#'
#' Attempts a plain (sparse, then dense, then QR) Newton solve. If every direct
#' solve fails or returns a non-finite step — the signature of a rank-deficient
#' / singular stacked Jacobian (issue M15) — it falls back to the regularised
#' normal-equations step
#' \deqn{(J'J + \lambda \, \mathrm{diag}(J'J)) \, \delta = -J' R,}
#' increasing \eqn{\lambda} geometrically until a finite step is obtained. The
#' LM step is always well-defined for \eqn{\lambda > 0} and reduces to the
#' Gauss-Newton/Newton direction as \eqn{\lambda \to 0}, so converged solutions
#' are unaffected — the damping only ever activates when the plain solve dies.
#'
#' @param J         dgCMatrix: stacked Jacobian
#' @param R         numeric: stacked residual (RHS is -R)
#' @param sparse_fallback Logical: allow dense/QR fallback before LM
#' @param lambda0   Initial LM damping (default 1e-8)
#' @param lambda_max Maximum LM damping before giving up (default 1e8)
#' @return List with $delta (numeric or NULL), $method (character),
#'   $lambda (numeric: damping used, 0 if a plain solve succeeded)
#' @noRd
.pf_robust_solve <- function(J, R, sparse_fallback = TRUE,
                             lambda0 = 1e-8, lambda_max = 1e8) {
  n <- length(R)

  # Coerce to a canonical column-compressed dgCMatrix. Matrix::solve(lu, .) can
  # crash on non-dgCMatrix sparse representations; production J is already a
  # dgCMatrix (sparseMatrix(repr="C")) so this is a cheap no-op guard.
  if (!methods::is(J, "dgCMatrix")) {
    J <- methods::as(J, "CsparseMatrix")
  }

  # ---- 1. Plain Newton solve (sparse LU → dense → QR) ----
  delta <- tryCatch({
    lu <- Matrix::lu(J)
    as.numeric(Matrix::solve(lu, -R))
  }, error = function(e) NULL)

  if (is.null(delta) && sparse_fallback) {
    delta <- tryCatch(as.numeric(Matrix::solve(J, -R)), error = function(e) NULL)
  }
  if (is.null(delta) && sparse_fallback) {
    delta <- tryCatch(solve(as.matrix(J), -R), error = function(e) NULL)
  }
  if (is.null(delta) && sparse_fallback) {
    delta <- tryCatch(qr.solve(as.matrix(J), -R), error = function(e) NULL)
  }

  if (!is.null(delta) && all(is.finite(delta))) {
    return(list(delta = delta, method = "newton", lambda = 0))
  }

  # ---- 2. Levenberg-Marquardt damped normal equations ----
  # (J'J + λ·diag(J'J)) δ = -J'R  — finite for any λ > 0 even if J is singular.
  JtJ  <- Matrix::crossprod(J)            # J'J  (n × n, symmetric PSD)
  Jtr  <- as.numeric(Matrix::crossprod(J, R))  # J'R
  dscale <- Matrix::diag(JtJ)
  # Marquardt scaling: damp each coordinate by its own curvature; guard zeros.
  dscale[!is.finite(dscale) | dscale <= 0] <- 1

  lambda <- lambda0
  while (lambda <= lambda_max) {
    A <- JtJ
    Matrix::diag(A) <- Matrix::diag(JtJ) + lambda * dscale
    delta <- tryCatch(
      as.numeric(Matrix::solve(A, -Jtr)),
      error = function(e) tryCatch(solve(as.matrix(A), -Jtr),
                                   error = function(e2) NULL)
    )
    if (!is.null(delta) && all(is.finite(delta))) {
      return(list(delta = delta, method = "lm", lambda = lambda))
    }
    lambda <- lambda * 10
  }

  list(delta = NULL, method = "failed", lambda = lambda)
}


# =============================================================================
# Main perfect-foresight solver
# =============================================================================

#' Perfect-foresight stacked Newton path solver for non-stationary models
#'
#' Solves the T-period deterministic transition path for DSGE models that
#' have **no stationary steady state** (balanced-growth-path models),
#' using Dynare-compatible stacked Newton on a sparse block-tridiagonal system.
#'
#' This function fills the gap left by \code{\link{pf_newton_solve}} (which
#' requires \code{n_eq == n_endo} and a steady state) and
#' \code{\link{mcp_solve_path}} (which targets MCP-constrained problems).
#' It handles the "pure" perfect-foresight case used by models like
#' \code{Ramsey_Cass_Koopmans} and \code{Solow_nonstationary}.
#'
#' **Boundary conditions:**
#' - \code{y0}: values at t=0 (from Dynare's \code{initval} block)
#' - \code{y_terminal}: values at t=T+1 (from Dynare's \code{endval} block)
#'
#' **Exogenous path:**
#' - \code{exo_path}: T × n_exo matrix giving the time path of each exogenous
#'   variable for periods 1..T. For BGP models, this captures trending
#'   technology (A) and labor (L).
#'
#' **Initial path:**
#' - If \code{Y_init} is provided, uses it as the warm start. Otherwise,
#'   initializes all periods to \code{y0} values (forward fill).
#'
#' @param compiled    dynhr_compiled (from \code{\link{compile_model}})
#' @param y0          Named numeric vector: endogenous state at t=0 (from
#'   \code{initval} block). Length must equal \code{n_endo}.
#' @param y_terminal  Named numeric vector: terminal condition at t=T+1 (from
#'   \code{endval} block). Length must equal \code{n_endo}.
#' @param exo_path    T × n_exo numeric matrix: time path of exogenous variables.
#'   Column names must match \code{varexo_names} of the compiled model. Each
#'   row t gives the exogenous state for period t.
#' @param params      Named numeric parameter vector.
#' @param n_periods   Integer: number of simulation periods T (default 100).
#' @param max_iter    Maximum Newton iterations (default 50).
#' @param tol         Convergence tolerance on max|R_stack| (default 1e-8).
#' @param step_size   Initial Newton step length (default 1.0).
#' @param line_search Logical: perform Armijo backtracking (default TRUE).
#' @param sparse_fallback Logical: use dense \code{solve()} if sparse fails
#'   (default TRUE).
#' @param homotopy    Continuation strategy for cold starts. \code{"auto"}
#'   (default) tries a direct cold solve first and, only if it fails to
#'   converge, escalates through a ladder of continuation runs (4, 16, 64
#'   stages) — so easy problems are unchanged and hard cold starts are rescued.
#'   An integer \code{>= 1} forces exactly that many continuation stages: the
#'   exogenous forcing is ramped linearly from \code{exo_init} (period-0 /
#'   baseline values) to the full \code{exo_path}, warm-starting each stage from
#'   the previous solution. \code{1} is a single direct cold solve.
#'   \code{"none"} is a single direct cold solve with no continuation fallback.
#' @param Y_init      Optional initial path (T × n_endo matrix) for warm-start.
#'   NULL → initialize each period with \code{y0}.
#' @param exo_init    Optional named numeric vector: period-0 values for
#'   exogenous variables (from the model's initval block). Used to supply
#'   the correct initial values for lagged exogenous terms like A(-1) at
#'   t=1. If NULL, defaults to the first row of \code{exo_path}.
#' @param exo_terminal Optional named numeric vector: post-terminal (t=T+1)
#'   values for exogenous variables (from the model's endval block). Used to
#'   supply the correct value for lead exogenous terms like \code{A(+1)} at the
#'   last period t=T — including those introduced as \code{AUX_EXO_LEAD_*}
#'   auxiliary variables. Required for a non-constant terminal exogenous path
#'   (e.g. a growing deterministic trend A=(1+g)^t, where A(T+1) differs from
#'   A(T)); if NULL, falls back to the last row of \code{exo_path} (= A(T)),
#'   which mis-states the terminal-period growth rates in such models.
#' @param verbose     Logical: print convergence progress (default FALSE).
#' @return List with:
#'   \item{Y}{T × n_endo solution matrix (rows = periods, cols = variables).}
#'   \item{converged}{Logical: TRUE if Newton converged.}
#'   \item{n_iter}{Integer: Newton iterations used.}
#'   \item{max_res}{Numeric: final max|R_stack|.}
#'   \item{endo_names}{Character: variable ordering of Y columns.}
#'   \item{merit_history}{Numeric vector: merit function value per iteration.}
#' @export
#'
#' @examples
#' \dontrun{
#' model <- parse_mod("Ramsey_Cass_Koopmans.mod")
#' compiled <- compile_model(model)
#' # From initval block:
#' y0 <- c(K = 12.3, C = 1.5, ...)
#' # From endval block:
#' y_term <- c(K = 45.6, C = 4.2, ...)
#' # Build exogenous path (A and L trending):
#' T <- 30
#' exo_path <- cbind(A = 1.02^(0:(T-1)), L = 1.01^(0:(T-1)))
#' result <- perfect_foresight_solve(compiled, y0, y_term,
#'                                    exo_path, model$param_values,
#'                                    n_periods = T)
#' }
perfect_foresight_solve <- function(compiled,
                                     y0,
                                     y_terminal,
                                     exo_path        = NULL,
                                     params          = NULL,
                                     n_periods       = 100L,
                                     max_iter        = 50L,
                                     tol             = 1e-8,
                                     step_size       = 1.0,
                                     line_search     = TRUE,
                                     sparse_fallback = TRUE,
                                     homotopy        = "auto",
                                     Y_init          = NULL,
                                     exo_init        = NULL,
                                     exo_terminal    = NULL,
                                     verbose         = FALSE) {

  dyn    <- compiled$dynamic
  n_endo <- length(dyn$endo_names)
  n_eq   <- dyn$n_eq

  if (is.null(params)) {
    params <- compiled$model$param_values
    if (is.null(params))
      stop("perfect_foresight_solve: params must be provided.")
  }

  if (n_eq != n_endo) {
    # When compiled$occbin is present (native OccBin model with bind/relax
    # equation pairs), delegate to occbin_solve_path() which implements the
    # GI candidate-regime-sequence algorithm.
    if (!is.null(compiled$occbin) && compiled$occbin$n_constraints > 0L) {
      message(sprintf(
        "perfect_foresight_solve: n_eq (%d) > n_endo (%d); ",
        n_eq, n_endo),
        "delegating to occbin_solve_path() for regime-switching simulation.")

      T_pf <- as.integer(n_periods)
      if (T_pf < 2L) stop("perfect_foresight_solve: n_periods must be >= 2.")

      # Build shock path matrix (T x n_exo)
      n_exo_pf <- length(dyn$exo_names)
      if (is.null(exo_path)) {
        shock_path_pf <- matrix(0, nrow = T_pf, ncol = n_exo_pf)
        colnames(shock_path_pf) <- dyn$exo_names
      } else {
        if (!is.matrix(exo_path)) {
          shock_path_pf <- matrix(exo_path, nrow = T_pf, ncol = n_exo_pf, byrow = TRUE)
        } else {
          shock_path_pf <- matrix(0, nrow = T_pf, ncol = n_exo_pf)
          colnames(shock_path_pf) <- dyn$exo_names
          if (!is.null(colnames(exo_path))) {
            for (nm in intersect(colnames(exo_path), dyn$exo_names))
              shock_path_pf[seq_len(min(nrow(exo_path), T_pf)), nm] <-
                exo_path[seq_len(min(nrow(exo_path), T_pf)), nm]
          } else {
            nc_pf <- min(ncol(exo_path), n_exo_pf)
            nr_pf <- min(nrow(exo_path), T_pf)
            shock_path_pf[seq_len(nr_pf), seq_len(nc_pf)] <-
              exo_path[seq_len(nr_pf), seq_len(nc_pf)]
          }
        }
      }

      # Extract constraints from compiled$occbin
      # compiled$occbin$parse_result$constraints is the named list from
      # occbin_parse_bind_relax(); these already carry var_name, condition,
      # bound_expr needed by .occbin_eval_bind_condition.
      pf_constraints <- compiled$occbin$parse_result$constraints

      # Use y_terminal as the relax-regime steady state (terminal condition)
      res_pf <- occbin_solve_path(
        compiled         = compiled,
        y0               = y0,
        y_ss             = y_terminal,
        shock_path       = shock_path_pf,
        params           = params,
        constraints      = pf_constraints,
        max_iter         = max_iter,
        tol              = tol,
        max_regime_iter  = 30L,
        step_size        = step_size,
        line_search      = line_search
      )

      # Return in perfect_foresight_solve output format, plus regime_path
      Y_out <- res_pf$Y
      rownames(Y_out) <- paste0("t", seq_len(T_pf))
      return(list(
        Y             = Y_out,
        converged     = res_pf$converged,
        n_iter        = res_pf$n_iter,
        max_res       = res_pf$max_res,
        endo_names    = dyn$endo_names,
        merit_history = numeric(0),
        regime_path   = res_pf$regime_path,
        outer_iter    = res_pf$outer_iter
      ))
    }
    stop(sprintf(
      "perfect_foresight_solve: n_eq (%d) != n_endo (%d). ",
      n_eq, n_endo),
      "Perfect foresight solver requires n_eq == n_endo. ",
      "For bind/relax models (n_eq > n_endo), use occbin_solve_path().")
  }

  T <- as.integer(n_periods)
  if (T < 2L) stop("perfect_foresight_solve: n_periods must be >= 2.")

  # ---- Normalize y0 and y_terminal ----
  y0_num <- as.numeric(y0[dyn$endo_names])
  y_term_num <- as.numeric(y_terminal[dyn$endo_names])

  if (anyNA(y0_num)) {
    stop("perfect_foresight_solve: y0 missing values for: ",
         paste(dyn$endo_names[is.na(y0_num)], collapse = ", "))
  }
  if (anyNA(y_term_num)) {
    stop("perfect_foresight_solve: y_terminal missing values for: ",
         paste(dyn$endo_names[is.na(y_term_num)], collapse = ", "))
  }

  # ---- Normalize exo_path ----
  n_exo <- length(dyn$exo_names)
  if (is.null(exo_path)) {
    # No exogenous path: all zeros
    eps_mat <- matrix(0, nrow = T, ncol = n_exo)
    colnames(eps_mat) <- dyn$exo_names
  } else {
    if (!is.matrix(exo_path)) {
      exo_path <- matrix(exo_path, nrow = T, ncol = n_exo, byrow = TRUE)
    }
    eps_mat <- matrix(0, nrow = T, ncol = n_exo)
    colnames(eps_mat) <- dyn$exo_names
    if (!is.null(colnames(exo_path))) {
      for (nm in intersect(colnames(exo_path), dyn$exo_names))
        eps_mat[, nm] <- exo_path[, nm]
    } else {
      nc <- min(ncol(exo_path), n_exo)
      eps_mat[, seq_len(nc)] <- exo_path[, seq_len(nc)]
    }
  }

  # ---- Precompute column metadata ----
  meta <- .pf_col_meta(dyn)

  # ---- Initialize path ----
  if (!is.null(Y_init)) {
    if (nrow(Y_init) != T || ncol(Y_init) != n_endo) {
      stop(sprintf("Y_init must be %d x %d (T x n_endo), got %d x %d",
                   T, n_endo, nrow(Y_init), ncol(Y_init)))
    }
    Y <- Y_init
    colnames(Y) <- dyn$endo_names
  } else {
    # Initialize via linear interpolation from y0 to y_terminal.
    # A flat initial path can lead to large Newton steps that produce
    # NaN (e.g., log of negative values). Linear interpolation gives
    # a smoother starting point that respects boundary conditions.
    Y <- matrix(0, nrow = T, ncol = n_endo)
    for (j in seq_len(n_endo)) {
      Y[, j] <- seq(from = y0_num[j], to = y_term_num[j], length.out = T)
    }
    colnames(Y) <- dyn$endo_names
  }

  # ==========================================================================
  # Inner Newton solve, parameterised by the exogenous path + warm start.
  # Reused by the homotopy/continuation driver below. Globalised with an
  # Armijo line search, a best-vetted-point fallback on line-search failure,
  # and a Levenberg-Marquardt step when the stacked Jacobian singularises
  # (issue M15 cold-start robustness).
  # ==========================================================================
  .pf_run_newton <- function(Y_start, eps_use) {
    Y_loc <- Y_start
    colnames(Y_loc) <- dyn$endo_names

    merit_fn <- function(Y_try) {
      ss <- .pf_build_stacked_system(
        Y_try, y0_num, y_term_num, eps_use,
        meta, compiled, params, T, n_endo,
        exo_init = exo_init, exo_terminal = exo_terminal
      )
      if (anyNA(ss$R) || any(!is.finite(ss$R))) return(1e30)  # invalid step
      0.5 * sum(ss$R * ss$R)
    }

    conv_loc      <- FALSE
    n_iter_loc    <- 0L
    max_res_loc   <- Inf
    merit_loc     <- numeric(0)

    for (iter in seq_len(max_iter)) {
      n_iter_loc <- iter

      sys <- .pf_build_stacked_system(
        Y_loc, y0_num, y_term_num, eps_use,
        meta, compiled, params, T, n_endo,
        exo_init = exo_init, exo_terminal = exo_terminal
      )

      # Non-finite residuals (e.g. log of a negative after an over-large
      # cold-start step): abort cleanly so the homotopy driver can retry with
      # finer staging. Guards against `if (NA)` in the convergence test.
      if (anyNA(sys$R) || any(!is.finite(sys$R))) {
        max_res_loc <- Inf
        break
      }

      max_res     <- max(abs(sys$R))
      max_res_loc <- max_res
      theta_cur   <- 0.5 * sum(sys$R * sys$R)
      merit_loc   <- c(merit_loc, theta_cur)

      if (verbose) {
        cat(sprintf("  PF Newton iter %3d: max|R| = %.2e, theta = %.2e\n",
                    iter, max_res, theta_cur))
      }

      if (max_res < tol) {
        conv_loc <- TRUE
        if (verbose) cat("  PF Newton converged.\n")
        break
      }

      # Robust Newton step: plain solve, then LM damping on singular Jacobian.
      sol <- .pf_robust_solve(sys$J, sys$R, sparse_fallback = sparse_fallback)
      delta_vec <- sol$delta
      if (verbose && sol$method == "lm") {
        cat(sprintf("  Singular Jacobian: Levenberg-Marquardt step (lambda=%.1e)\n",
                    sol$lambda))
      }

      if (is.null(delta_vec)) {
        if (!sparse_fallback) {
          stop("perfect_foresight_solve: singular Jacobian at iter ", iter)
        }
        warning("perfect_foresight_solve: could not compute a finite Newton ",
                "step at iter ", iter, " (Jacobian rank-deficient); aborting.")
        break
      }

      # Line search or direct step
      if (line_search) {
        ls <- .pf_line_search(
          Y_loc, delta_vec, sys$R, sys$J, theta_cur, merit_fn,
          sigma = 1e-4, max_ls = 20L, n_endo = n_endo, T = T
        )

        if (is.null(ls$Y_new)) {
          # No tried step reduced the merit. Do NOT take a blind half step
          # (the old behaviour blew the path up and singularised the next
          # Jacobian — M15). Stop this Newton run cleanly; the homotopy
          # driver can retry with continuation.
          if (verbose) cat("  Line search found no descent; stopping run.\n")
          break
        }
        Y_loc <- ls$Y_new
      } else {
        alpha <- step_size
        for (t in seq_len(T)) {
          idx_t <- (t - 1L) * n_endo + seq_len(n_endo)
          Y_loc[t, ] <- Y_loc[t, ] + alpha * delta_vec[idx_t]
        }
      }
    }

    list(Y = Y_loc, converged = conv_loc, n_iter = n_iter_loc,
         max_res = max_res_loc, merit_history = merit_loc)
  }

  # ==========================================================================
  # Homotopy / continuation driver.
  # Ramp the exogenous forcing from a baseline (exo_init, broadcast over T) to
  # the full eps_mat over `n_stage` stages, warm-starting each stage from the
  # previous solution. A single stage (n_stage == 1) is a plain cold solve.
  # ==========================================================================
  # Baseline forcing for stage 0: the period-0 exogenous values (initval), or
  # the first row of eps_mat when exo_init is absent.
  if (!is.null(exo_init) && !anyNA(as.numeric(exo_init[dyn$exo_names]))) {
    base_row <- as.numeric(exo_init[dyn$exo_names])
  } else {
    base_row <- eps_mat[1L, ]
  }
  eps_base <- matrix(rep(base_row, each = T), nrow = T, ncol = n_exo)
  colnames(eps_base) <- dyn$exo_names

  run_homotopy <- function(n_stage) {
    n_stage <- max(1L, as.integer(n_stage))
    Y_cur   <- Y
    res     <- NULL
    for (st in seq_len(n_stage)) {
      frac <- if (n_stage == 1L) 1 else st / n_stage
      eps_st <- (1 - frac) * eps_base + frac * eps_mat
      if (verbose && n_stage > 1L) {
        cat(sprintf("== Homotopy stage %d/%d (forcing fraction %.3f) ==\n",
                    st, n_stage, frac))
      }
      res   <- .pf_run_newton(Y_cur, eps_st)
      Y_cur <- res$Y
      # If an intermediate stage stalls badly, abort the ladder early.
      if (!res$converged && st < n_stage && !all(is.finite(res$Y))) break
    }
    res
  }

  # Resolve the homotopy strategy. "auto" = one direct solve then escalate on
  # failure; "none" = a single direct solve with no continuation/fallback (M27c
  # — previously errored on as.integer("none")); an integer = that many fixed
  # continuation stages.
  auto_mode <- identical(homotopy, "auto")
  none_mode <- identical(homotopy, "none")
  n_stage   <- if (auto_mode || none_mode) 1L else max(1L, as.integer(homotopy))

  result <- run_homotopy(n_stage)

  # Auto fallback: if a direct cold solve failed, escalate to continuation.
  if (auto_mode && !result$converged) {
    if (verbose) cat("== Direct solve failed; retrying with continuation ==\n")
    for (ns in c(4L, 16L, 64L)) {
      result <- run_homotopy(ns)
      if (result$converged) break
    }
  }

  Y <- result$Y
  colnames(Y) <- dyn$endo_names
  rownames(Y) <- paste0("t", seq_len(T))

  if (!result$converged && verbose) {
    cat(sprintf("  PF Newton did NOT converge. max|R| = %.2e\n", result$max_res))
  }

  list(
    Y             = Y,
    converged     = result$converged,
    n_iter        = result$n_iter,
    max_res       = result$max_res,
    endo_names    = dyn$endo_names,
    merit_history = result$merit_history
  )
}
