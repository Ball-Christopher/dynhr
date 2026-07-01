## R/pf-newton.R
## --------------------------------------------------------------------------
## Perfect-foresight Newton path solver for nonlinear DSGE models with OBCs.
##
## Provides:
##   .pf_col_meta()            -- precompute dyn-column kind/index metadata
##   .pf_make_dy()             -- build named dy vector for one period
##   .pf_preallocate_triplets() -- NNZ-budget triplet buffer (sparse path)
##   .pf_build_sparse_system() -- assemble dgCMatrix + residual (sparse path)
##   pf_newton_solve()         -- stacked Newton + active-set OBC iteration
##
## The solver stacks T*n_eq equations over the path Y = (y_1,...,y_T):
##   F_i(y_{t-1}, y_t, y_{t+1}, eps_t) = 0    t = 1..T, i = 1..n_eq
## Boundary conditions: y_0 = y0 (given), y_{T+1} = y_ss (terminal).
##
## OBCs are handled by active-set iteration:
##   1. Fix a regime_path (n_spec × T logical matrix).
##   2. Solve Newton on the modified system (binding equations replace slack).
##   3. Check complementarity; flip violations; repeat.
##
## References:
##   Adjemian & Juillard (2025) — Stochastic Extended Path
##   Dynare perfect_foresight_solver (naming/convention reference)
## --------------------------------------------------------------------------


# =============================================================================
# Internal helpers
# =============================================================================

#' Precompute dynamic-column metadata for fast dy assembly
#'
#' @param dyn  compiled$dynamic (from build_dynamic_model)
#' @return List with keys (char), kind (char), var_idx (int), n_cols (int)
#' @noRd
.pf_col_meta <- function(dyn) {
  cmap  <- dyn$dyn_col_map
  n     <- nrow(cmap)
  keys  <- character(n)
  kind  <- character(n)
  vidx  <- integer(n)

  for (k in seq_len(n)) {
    nm  <- cmap$name[k]
    ll  <- cmap$lead_lag[k]
    sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll) else paste0("__m", abs(ll))
    keys[k] <- paste0(nm, sfx)

    if (nm %in% dyn$exo_names) {
      kind[k] <- "exo"
      vidx[k] <- match(nm, dyn$exo_names)
    } else {
      vidx[k] <- match(nm, dyn$endo_names)
      kind[k]  <- if (ll < 0L) "lag" else if (ll == 0L) "cur" else "lead"
    }
  }

  list(keys = keys, kind = kind, var_idx = vidx, n_cols = dyn$total_cols,
       dyn_col = cmap$col)
}


#' Build named dy vector for one period
#'
#' @param meta      Output of .pf_col_meta()
#' @param Y         T × n_endo numeric matrix (current path)
#' @param y0_num    Length-n_endo numeric: initial state (unordered is fine)
#' @param y_ss_num  Length-n_endo numeric: terminal steady state
#' @param eps_row   Length-n_exo numeric: shocks at period t
#' @param t         Period index (1-based)
#' @param T         Total horizon
#' @return Named numeric vector of length n_cols
#' @noRd
.pf_make_dy <- function(meta, Y, y0_num, y_ss_num, eps_row, t, T) {
  dy        <- numeric(meta$n_cols)
  names(dy) <- meta$keys

  for (k in seq_along(meta$keys)) {
    vi <- meta$var_idx[k]
    dc <- meta$dyn_col[k]
    dy[dc] <- switch(meta$kind[k],
      lag  = if (t == 1L) y0_num[vi]  else Y[t - 1L, vi],
      cur  = Y[t, vi],
      lead = if (t == T)  y_ss_num[vi] else Y[t + 1L, vi],
      exo  = eps_row[vi]
    )
  }

  dy
}


#' Preallocate triplet index vectors for the sparse stacked Jacobian
#'
#' Uses dyn$jac_triplets to estimate NNZ budget. Falls back to growable
#' vectors when jac_triplets is NULL or empty (Landmine 3).
#'
#' @param dyn   compiled$dynamic
#' @param T     Integer: horizon
#' @return List(i, j, v, nnz_est, preallocated)
#' @noRd
.pf_preallocate_triplets <- function(dyn, T) {
  jt <- dyn$jac_triplets

  if (is.null(jt) || length(jt) == 0L) {
    # Fallback: growable vectors (jac_tape NULL path)
    return(list(i = integer(0), j = integer(0), v = numeric(0),
                nnz_est = NA_integer_, preallocated = FALSE))
  }

  # Upper bound: T * length(jt); boundary periods have slightly fewer entries
  # but over-allocating by ~2*n_endo and trimming is safe.
  nnz_est <- T * length(jt)

  list(
    i            = integer(nnz_est),
    j            = integer(nnz_est),
    v            = numeric(nnz_est),
    nnz_est      = nnz_est,
    preallocated = TRUE
  )
}


#' Build stacked residual vector and sparse Jacobian (dgCMatrix)
#'
#' Assembles the T*n_eq x T*n_endo block-tridiagonal Jacobian and the
#' stacked residual for the current path Y with the given OBC regime.
#' OBC row modifications are applied to Rt/Jt BEFORE the triplet scatter
#' (Landmine 2).
#'
#' @param Y          T x n_endo numeric matrix
#' @param y0_num     Length-n_endo initial state
#' @param y_ss_num   Length-n_endo terminal steady state
#' @param eps_mat    T x n_exo shock matrix
#' @param meta       Output of .pf_col_meta()
#' @param dyn        compiled$dynamic
#' @param params     Named numeric parameter vector
#' @param y_ss       Named steady-state vector (for residuals_fn / jacobian_fn)
#' @param T, n_eq, n_endo  Integers
#' @param obc_specs  List of OBC specs (may be empty)
#' @param regime     n_spec x T logical matrix
#' @param spec_cur_dc Integer vector: current-period dyn column per OBC spec
#' @param trip_buf   Pre-allocated triplet buffer from .pf_preallocate_triplets()
#' @return List($R, $J)
#' @noRd
.pf_build_sparse_system <- function(Y, y0_num, y_ss_num, eps_mat,
                                     meta, dyn, params, y_ss,
                                     T, n_eq, n_endo,
                                     obc_specs, regime, spec_cur_dc,
                                     trip_buf) {
  n_spec  <- length(obc_specs)
  n_total <- T * n_eq
  n_exo   <- length(dyn$exo_names)
  R_full  <- numeric(n_total)

  pre <- trip_buf$preallocated
  if (pre) {
    i_t <- trip_buf$i
    j_t <- trip_buf$j
    v_t <- trip_buf$v
    ptr <- 0L  # next free slot (0-based)
  } else {
    i_t <- integer(0)
    j_t <- integer(0)
    v_t <- numeric(0)
  }

  for (t in seq_len(T)) {
    row_off <- (t - 1L) * n_eq
    col_off <- (t - 1L) * n_endo

    # Use .pf_make_dy_exo to fill exogenous lead/lag columns that may be
    # referenced by auxiliary-variable equations (e.g. AUX_EXO_LEAD_x_1 = x(+1)).
    dy <- .pf_make_dy_exo(meta, Y, y0_num, y_ss_num, eps_mat, t, T, n_exo)
    Rt <- dyn$residuals_fn(dy, params, y_ss)
    Jt <- dyn$jacobian_fn(dy, params, y_ss)

    # OBC: replace binding equations in Rt/Jt BEFORE scatter (Landmine 2)
    if (n_spec > 0L) {
      for (s in seq_len(n_spec)) {
        if (!regime[s, t]) next
        spec <- obc_specs[[s]]
        ei   <- spec$eq_idx
        vi   <- spec$var_idx
        bnd  <- spec$bound
        Rt[ei]    <- Y[t, vi] - bnd
        Jt[ei, ]  <- 0
        Jt[ei, spec_cur_dc[s]] <- 1
      }
    }

    R_full[row_off + seq_len(n_eq)] <- Rt

    # Scatter Jt columns into sparse triplets
    for (k in seq_along(meta$kind)) {
      kd <- meta$kind[k]
      if (kd == "exo") next
      vi <- meta$var_idx[k]
      dc <- meta$dyn_col[k]

      gcol <- switch(kd,
        lag  = if (t > 1L) (t - 2L) * n_endo + vi else NA_integer_,
        cur  = col_off + vi,
        lead = if (t < T)  t * n_endo + vi         else NA_integer_
      )
      if (is.na(gcol)) next

      col_vals <- Jt[, dc]

      if (pre) {
        for (eq in seq_len(n_eq)) {
          val <- col_vals[eq]
          if (val != 0 && is.finite(val)) {
            ptr <- ptr + 1L
            i_t[ptr] <- row_off + eq
            j_t[ptr] <- gcol
            v_t[ptr] <- val
          }
        }
      } else {
        for (eq in seq_len(n_eq)) {
          val <- col_vals[eq]
          if (val != 0 && is.finite(val)) {
            i_t <- c(i_t, row_off + eq)
            j_t <- c(j_t, gcol)
            v_t <- c(v_t, val)
          }
        }
      }
    }
  }

  # Trim preallocated buffers to actual NNZ
  if (pre && ptr > 0L) {
    i_t <- i_t[seq_len(ptr)]
    j_t <- j_t[seq_len(ptr)]
    v_t <- v_t[seq_len(ptr)]
  }

  # Build dgCMatrix — repr = "C" replaces deprecated giveCsparse = TRUE (Landmine 5)
  if (length(i_t) == 0L) {
    J <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 0,
                              dims = c(n_total, n_total),
                              repr = "C")
  } else {
    J <- Matrix::sparseMatrix(
      i = i_t, j = j_t, x = v_t,
      dims = c(n_total, n_total),
      index1 = TRUE,
      repr = "C"
    )
  }

  list(R = R_full, J = J)
}


# =============================================================================
# Main solver
# =============================================================================

## Residual-only max-norm at a trial path Y (no Jacobian). Used by the Newton
## backtracking line search: a full Newton step can leave the model domain
## (e.g. log/sqrt of a non-positive value), giving a non-finite residual.
## Returns Inf on any non-finite residual so the caller can shrink the step.
.pf_resnorm <- function(Y, y0_num, y_ss_num, eps_mat, meta, dyn, params, y_ss,
                        T, obc_specs, regime, n_spec) {
  n_exo <- length(dyn$exo_names)
  mx <- 0
  for (t in seq_len(T)) {
    dy <- .pf_make_dy_exo(meta, Y, y0_num, y_ss_num, eps_mat, t, T, n_exo)
    Rt <- suppressWarnings(dyn$residuals_fn(dy, params, y_ss))
    if (n_spec > 0L) for (s in seq_len(n_spec)) {
      if (!regime[s, t]) next
      spec <- obc_specs[[s]]
      Rt[spec$eq_idx] <- Y[t, spec$var_idx] - spec$bound
    }
    m <- max(abs(Rt))
    if (!is.finite(m)) return(Inf)
    if (m > mx) mx <- m
  }
  mx
}

#' Perfect-foresight Newton path solver for nonlinear DSGE models with OBCs
#'
#' Solves the T-period deterministic path
#' \preformatted{
#'   F_i(y_{t-1}, y_t, y_{t+1}, eps_t) = 0   for t = 1..T, i = 1..n_eq
#' }
#' with \eqn{y_0} = y0 (given) and \eqn{y_{T+1}} = y_ss (terminal transversality),
#' subject to optional occasionally-binding constraints handled via active-set
#' (complementarity-flip) iteration over the regime path.
#'
#' Unlike the perturbation + LCP solvers (pf-lcp.R, obc-lcp.R), this operates
#' directly on the nonlinear dynamic residuals and therefore does NOT require
#' model(linear). It is the inner loop used by the SEP outer solver (obc-sep.R).
#'
#' @param compiled  dynhr_compiled (from compile_model())
#' @param y0        Named numeric vector: endogenous state at t=0 (initial
#'                  conditions). Length n_endo. Names must match var_names.
#' @param y_ss      Named numeric vector: steady state (terminal condition).
#' @param shock_path  T × n_exo numeric matrix of structural shocks. Column
#'                    names must match varexo_names; missing columns are zero.
#' @param params    Named numeric parameter vector.
#' @param obc_specs List of OBC specs from obc_parse_tags() or
#'                  obc_collect_specs() (empty list = no OBCs). Each spec must
#'                  have $eq_idx, $var_idx, $var_name, $op, $bound.
#' @param max_iter        Maximum Newton iterations per regime (default 50).
#' @param tol             Newton convergence tolerance on max|R| (default 1e-10).
#' @param max_regime_iter Maximum active-set regime switches (default 30).
#' @param step_size       Newton step-length (default 1; reduce for ill-conditioned
#'                        problems or strong nonlinearities).
#' @param method    Linear algebra method: \code{"auto"} (default — selects
#'                  sparse when T*n_endo >= 100), \code{"sparse"} (always use
#'                  dgCMatrix + Matrix::solve with LU reuse), or \code{"dense"}
#'                  (original base::solve on a dense matrix).
#' @param verbose   Logical (default FALSE): print per-iteration Newton residuals
#'                  and regime-switch progress (matches perfect_foresight_solve()).
#' @return List with:
#'   \itemize{
#'     \item \code{Y} — T x n_endo solution matrix (rows = periods, cols = variables)
#'     \item \code{regime} — n_spec x T logical matrix (TRUE = binding)
#'     \item \code{irf} — T x n_endo deviation from steady state
#'     \item \code{converged} — Logical: TRUE if Newton converged on the final regime
#'     \item \code{n_iter} — Integer: Newton iterations on the final regime
#'     \item \code{endo_names} — Character: variable ordering of Y columns
#'   }
#' @export
pf_newton_solve <- function(compiled,
                            y0,
                            y_ss,
                            shock_path,
                            params,
                            obc_specs      = list(),
                            max_iter       = 50L,
                            tol            = 1e-10,
                            max_regime_iter = 30L,
                            step_size      = 1.0,
                            method         = c("auto", "sparse", "dense"),
                            verbose        = FALSE) {

  dyn    <- compiled$dynamic
  n_endo <- length(dyn$endo_names)
  n_eq   <- dyn$n_eq
  n_exo  <- length(dyn$exo_names)
  n_spec <- length(obc_specs)

  if (n_eq != n_endo) stop(sprintf(
    "pf_newton_solve: n_eq (%d) != n_endo (%d). Check model equations.", n_eq, n_endo))

  method <- match.arg(method)

  # Normalize shock_path to T x n_exo matrix in exo_names order
  if (!is.matrix(shock_path)) shock_path <- matrix(shock_path, nrow = 1L)
  T <- nrow(shock_path)
  eps_mat <- matrix(0, nrow = T, ncol = length(dyn$exo_names))
  colnames(eps_mat) <- dyn$exo_names
  if (!is.null(colnames(shock_path))) {
    for (nm in intersect(colnames(shock_path), dyn$exo_names))
      eps_mat[, nm] <- shock_path[, nm]
  } else {
    nc <- min(ncol(shock_path), length(dyn$exo_names))
    eps_mat[, seq_len(nc)] <- shock_path[, seq_len(nc)]
  }

  # Align y0 and y_ss to endo_names order
  y0_num  <- as.numeric(y0[dyn$endo_names])
  y_ss_num <- as.numeric(y_ss[dyn$endo_names])

  # Precompute dy column metadata once
  meta <- .pf_col_meta(dyn)

  # For each OBC spec, find the dyn column index for var__0 (current period)
  cmap <- dyn$dyn_col_map
  spec_cur_dc <- if (n_spec > 0L) {
    vapply(obc_specs, function(s) {
      k <- which(cmap$name == s$var_name & cmap$lead_lag == 0L)
      if (length(k) == 0L) stop(sprintf(
        "pf_newton_solve: OBC var '%s' has no current-period column in dyn_col_map.",
        s$var_name))
      cmap$col[k[1L]]
    }, integer(1))
  } else integer(0)

  # Resolve auto method: sparse when T*n_endo >= 100 (Landmine 6)
  use_sparse <- switch(method,
    auto   = (T * n_endo >= 100L),
    sparse = TRUE,
    dense  = FALSE
  )

  # Preallocate triplet buffer once per solve (sparse path only).
  # Reused across Newton iterations within each regime.
  trip_buf <- if (use_sparse) .pf_preallocate_triplets(dyn, T) else NULL

  # Initialize path via linear interpolation from y0 to y_ss.
  # A flat-SS initial path can violate model domain constraints (e.g. log/sqrt
  # of a non-positive value) when y0 differs from y_ss, causing a non-finite
  # residual at the very first Newton iteration — before the backtracking line
  # search can engage. Linear interpolation respects both boundary conditions
  # and keeps the path on-domain for smooth nonlinear models (ramst, RBC).
  Y <- matrix(0, nrow = T, ncol = n_endo)
  for (j in seq_len(n_endo)) {
    Y[, j] <- seq(from = y0_num[j], to = y_ss_num[j], length.out = T)
  }
  colnames(Y) <- dyn$endo_names

  # Fallback: if the interpolated init is STILL off-domain (very unusual —
  # means the straight line from y0 to y_ss passes through an infeasible region),
  # fall back to the flat-SS init with a warning. The iteration will then rely on
  # the backtracking line search; if that also fails the caller sees converged=FALSE.
  {
    rn_interp <- .pf_resnorm(Y, y0_num, y_ss_num, eps_mat, meta, dyn,
                              params, y_ss, T, list(), matrix(FALSE, 1L, T), 0L)
    if (!is.finite(rn_interp)) {
      warning("pf_newton_solve: interpolated initial path is off-domain ",
              "(residual non-finite); falling back to flat-SS init. ",
              "Consider supplying a feasible Y_init or using ",
              "perfect_foresight_solve(line_search = TRUE).")
      Y <- matrix(rep(y_ss_num, T), nrow = T, byrow = TRUE)
      colnames(Y) <- dyn$endo_names
    }
  }

  regime <- matrix(FALSE, nrow = max(n_spec, 1L), ncol = T)

  converged_newton <- FALSE
  n_iter_final     <- 0L

  for (regime_iter in seq_len(max_regime_iter)) {

    if (verbose) cat(sprintf("== Regime iter %d ==\n", regime_iter))

    # ------------------------------------------------------------------
    # Newton iterations for fixed regime
    # ------------------------------------------------------------------
    converged_newton <- FALSE
    n_iter_final     <- 0L

    # LU factor cache: computed at iter==1, reused on subsequent iters,
    # invalidated on regime change (Landmine 1).
    lu_cache <- NULL
    lu_valid <- FALSE

    for (iter in seq_len(max_iter)) {

      if (use_sparse) {
        # ---- Sparse path ----
        sys <- .pf_build_sparse_system(
          Y, y0_num, y_ss_num, eps_mat,
          meta, dyn, params, y_ss,
          T, n_eq, n_endo,
          obc_specs, regime, spec_cur_dc,
          trip_buf
        )
        R_full   <- sys$R
        J_sparse <- sys$J

        res_norm <- max(abs(R_full))
        n_iter_final <- iter
        if (verbose) cat(sprintf("  Newton iter %d: max|R| = %.3e\n", iter, res_norm))
        if (!is.finite(res_norm)) {
          warning("pf_newton_solve: non-finite residual at iter ", iter,
                  " (the initial path left the model domain, e.g. log/sqrt of a ",
                  "non-positive value). Supply a feasible initial path or use ",
                  "perfect_foresight_solve(line_search = TRUE).")
          break
        }
        if (res_norm < tol) { converged_newton <- TRUE; break }

        # Compute (or reuse) LU factorisation.
        # Always recompute on iter == 1 (start of each Newton loop) — Landmine 1.
        if (iter == 1L || !lu_valid) {
          lu_cache <- Matrix::lu(J_sparse)
          lu_valid <- TRUE
        }

        delta_vec <- tryCatch(
          as.numeric(Matrix::solve(lu_cache, -R_full)),
          error = function(e) {
            # Stale factorisation: recompute and retry once
            lu_cache <<- Matrix::lu(J_sparse)
            tryCatch(
              as.numeric(Matrix::solve(lu_cache, -R_full)),
              error = function(e2) rep(NA_real_, length(R_full))
            )
          }
        )

      } else {
        # ---- Dense path (byte-for-byte identical to original) ----
        R_full <- numeric(T * n_eq)
        J_full <- matrix(0, nrow = T * n_eq, ncol = T * n_endo)

        for (t in seq_len(T)) {
          row_off <- (t - 1L) * n_eq
          col_off <- (t - 1L) * n_endo
          rows_t  <- row_off + seq_len(n_eq)

          # Use .pf_make_dy_exo to fill exo lead/lag columns for AUX equations
          dy <- .pf_make_dy_exo(meta, Y, y0_num, y_ss_num, eps_mat, t, T, n_exo)
          Rt <- dyn$residuals_fn(dy, params, y_ss)
          Jt <- dyn$jacobian_fn(dy, params, y_ss)

          # OBC: replace binding equations
          if (n_spec > 0L) {
            for (s in seq_len(n_spec)) {
              if (!regime[s, t]) next
              spec <- obc_specs[[s]]
              ei   <- spec$eq_idx
              vi   <- spec$var_idx
              bnd  <- spec$bound
              Rt[ei]    <- Y[t, vi] - bnd
              Jt[ei, ]  <- 0
              Jt[ei, spec_cur_dc[s]] <- 1
            }
          }

          R_full[rows_t] <- Rt

          # Assemble Jacobian columns into global J
          for (k in seq_along(meta$kind)) {
            kd <- meta$kind[k]
            if (kd == "exo") next
            vi <- meta$var_idx[k]
            dc <- meta$dyn_col[k]

            gcol <- switch(kd,
              lag  = if (t > 1L) (t - 2L) * n_endo + vi else NA_integer_,
              cur  = col_off + vi,
              lead = if (t < T)  t * n_endo + vi         else NA_integer_
            )
            if (is.na(gcol)) next
            J_full[rows_t, gcol] <- J_full[rows_t, gcol] + Jt[, dc]
          }
        }

        res_norm <- max(abs(R_full))
        n_iter_final <- iter
        if (verbose) cat(sprintf("  Newton iter %d: max|R| = %.3e\n", iter, res_norm))
        if (!is.finite(res_norm)) {
          warning("pf_newton_solve: non-finite residual at iter ", iter,
                  " (the initial path left the model domain, e.g. log/sqrt of a ",
                  "non-positive value). Supply a feasible initial path or use ",
                  "perfect_foresight_solve(line_search = TRUE).")
          break
        }
        if (res_norm < tol) { converged_newton <- TRUE; break }

        # Newton step: J * delta = -R
        delta_vec <- solve(J_full, -R_full)
      }

      if (anyNA(delta_vec)) {
        warning("pf_newton_solve: singular Jacobian at iter ", iter,
                "; aborting Newton.")
        break
      }

      # Update path with a DOMAIN-SAFE step. Take the full Newton step whenever
      # it yields a finite residual (full steps are what reach the correct basin
      # on hard cold starts -- damping/Armijo can stall in the wrong basin, e.g.
      # Stock_SIR). Only when the full step leaves the model domain (capital/
      # investment non-positive -> log() = NaN -> previously a cryptic NaN-norm
      # crash) do we backtrack, halving until the trial path is FINITE again.
      # This keeps full-Newton convergence on well-behaved/non-monotone problems
      # while preventing the off-domain crash; a model whose full step genuinely
      # diverges simply hits max_iter (converged = FALSE) instead of erroring.
      delta_mat <- matrix(delta_vec[seq_len(T * n_endo)], nrow = T, byrow = TRUE)
      alpha     <- step_size
      accepted  <- FALSE
      for (bt in 0:30L) {
        Y_try  <- Y + alpha * delta_mat
        rn_try <- .pf_resnorm(Y_try, y0_num, y_ss_num, eps_mat, meta, dyn,
                              params, y_ss, T, obc_specs, regime, n_spec)
        if (is.finite(rn_try)) { Y <- Y_try; accepted <- TRUE; break }
        alpha <- alpha / 2
      }
      if (!accepted) {
        warning("pf_newton_solve: no finite-residual Newton step at iter ",
                iter, " even after backtracking; aborting Newton. Try ",
                "perfect_foresight_solve(line_search = TRUE) or a feasible ",
                "initial path.")
        break
      }
    }

    if (n_spec == 0L) break  # no OBCs: done after Newton

    # ------------------------------------------------------------------
    # Active-set complementarity check
    # ------------------------------------------------------------------
    new_regime <- regime
    any_flip   <- FALSE

    for (s in seq_len(n_spec)) {
      spec  <- obc_specs[[s]]
      vi    <- spec$var_idx
      bnd   <- spec$bound
      sgn   <- if (spec$op == ">") 1 else -1  # >: lower bound; <: upper bound

      for (t in seq_len(T)) {
        slack_val <- sgn * (Y[t, vi] - bnd)   # > 0 means constraint satisfied

        if (!regime[s, t]) {
          # Slack period: flip to binding if variable violates bound
          if (slack_val < 0) { new_regime[s, t] <- TRUE; any_flip <- TRUE }
        } else {
          # Binding period: flip to slack if bound is already satisfied
          if (slack_val > 0) { new_regime[s, t] <- FALSE; any_flip <- TRUE }
        }
      }
    }

    if (!any_flip) break
    # Regime changed: invalidate LU cache (Landmine 1)
    lu_cache <- NULL
    lu_valid <- FALSE
    regime <- new_regime
  }

  if (verbose) {
    if (converged_newton)
      cat(sprintf("  pf_newton_solve converged in %d Newton iter(s).\n", n_iter_final))
    else
      cat(sprintf("  pf_newton_solve did NOT converge (max_iter=%d, max_regime_iter=%d).\n",
                  max_iter, max_regime_iter))
  }

  y_ss_mat <- matrix(rep(y_ss_num, T), nrow = T, byrow = TRUE)
  colnames(Y) <- dyn$endo_names

  list(
    Y          = Y,
    regime     = regime[seq_len(n_spec), , drop = FALSE],
    irf        = Y - y_ss_mat,
    converged  = converged_newton,
    n_iter     = n_iter_final,
    endo_names = dyn$endo_names
  )
}
