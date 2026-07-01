## R/diag-pre-d40-near-unit-root.R
## --------------------------------------------------------------------------
## D40 Near-unit-root proximity diagnostic.
##
## Diagnoses how close a DSGE model is to the determinacy boundary by
## computing the spectral radius of the stable state-transition block and
## attributing it to structural parameters via central finite differences.
## --------------------------------------------------------------------------

#' D40. Near-unit-root proximity and eigenvalue-sensitivity attribution
#'
#' Solves the model at a given parameter vector \code{theta} and reports:
#' \describe{
#'   \item{\code{spectral_radius}}{Modulus of the eigenvalue closest to the
#'     unit circle from inside (max |λ| over stable roots, i.e. the spectral
#'     radius of \code{ghx}).}
#'   \item{\code{distance_to_unit}}{1 - spectral_radius.  Near zero = near
#'     unit root = highly persistent / close to determinacy boundary.}
#'   \item{\code{critical_eigenvalue}}{The complex eigenvalue achieving the
#'     spectral radius.}
#'   \item{\code{param_sensitivity}}{data.frame with columns \code{param}
#'     and \code{d_rhomax_dparam}, sorted descending by \code{|d_rhomax_dparam|}.
#'     Shows which parameters push the critical root toward the unit circle.}
#'   \item{\code{equation_attribution}}{Best-effort list: maps sensitivity
#'     back to model equations via the structural Jacobian column norms.  May
#'     be \code{NULL} with a note if equation info is unavailable.}
#' }
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Solve the model at \code{theta} (steady state + first-order
#'     perturbation). Extract \code{ghx} — the state-to-state transition
#'     matrix.
#'   \item Compute eigenvalues of \code{ghx}.  The spectral radius is
#'     \code{rho_max = max(Mod(eigen(ghx)))}.  The critical eigenvalue is the
#'     one achieving this maximum.
#'   \item For each parameter in \code{param_names}, compute
#'     \eqn{d\rho_\max / d\theta_i} by central finite differences: re-solve
#'     the model at \code{theta ± h}, recompute \code{rho_max}, difference.
#'   \item (Best-effort) Attribute via the static Jacobian: the columns of
#'     \code{f_zero} corresponding to the critical eigenvector's state
#'     variables carry the dominant equation-level structural dependence.
#' }
#'
#' @param model        A \code{dynhr_mod} object (output of \code{parse_mod}).
#' @param compiled     A \code{dynhr_compiled} (output of
#'   \code{compile_model}).  Required for re-solves.
#' @param theta        Named numeric vector of parameter values.  If \code{NULL},
#'   uses \code{model$param_values}.
#' @param params       Alias for \code{theta} (either may be supplied).
#' @param param_names  Character vector selecting which parameters to
#'   differentiate.  Defaults to all names in \code{theta}.
#' @param h            Finite-difference step size (default \code{1e-5}).
#' @param ...          Unused; reserved for future arguments.
#'
#' @return A \code{dynhr_diagnostic} list (see \code{.make_result}) with
#'   \code{$result} containing:
#'   \itemize{
#'     \item \code{spectral_radius}  -- numeric scalar
#'     \item \code{distance_to_unit} -- numeric scalar (1 - spectral_radius)
#'     \item \code{critical_eigenvalue} -- complex scalar
#'     \item \code{param_sensitivity}  -- data.frame (param, d_rhomax_dparam)
#'     \item \code{equation_attribution} -- list or NULL
#'   }
#'
#' @export
diag_near_unit_root <- function(model,
                                compiled    = NULL,
                                theta       = NULL,
                                params      = NULL,
                                param_names = NULL,
                                h           = 1e-5,
                                ...) {

  ## ---- 0. Normalise inputs -----------------------------------------------
  ## Accept either `theta` or `params` as the parameter vector.
  if (is.null(theta) && !is.null(params)) theta <- params
  if (is.null(theta)) {
    if (!is.null(model$param_values)) {
      theta <- model$param_values
    } else {
      return(.make_result(
        result  = NULL,
        pass    = NA,
        plots   = list(),
        summary = "D40 Near-unit-root: supply theta or params."
      ))
    }
  }
  if (is.null(names(theta)) && !is.null(model$param_names)) {
    names(theta) <- model$param_names
  }
  if (is.null(param_names)) param_names <- names(theta)
  if (is.null(param_names)) param_names <- paste0("p_", seq_along(theta))

  ## ---- 1. Solve model at baseline theta ----------------------------------
  .solve_at <- function(th) {
    ## Merge perturbed theta into a full param list so that downstream
    ## solvers (extract_system_matrices, Jacobian tape) see the right values.
    params_full <- model$param_values
    for (nm in names(th)) params_full[[nm]] <- th[[nm]]

    ss <- tryCatch(
      solve_steady(compiled, params_full,
                   y0      = setNames(rep(0, length(model$var_names)),
                                      model$var_names),
                   verbose = FALSE),
      error = function(e) list(converged = FALSE)
    )
    if (!isTRUE(ss$converged)) return(NULL)

    dr <- tryCatch(
      solve_perturbation(model, compiled,
                         ss      = ss$values,
                         params  = params_full,
                         order   = 1L,
                         verbose = FALSE),
      error = function(e) NULL
    )
    dr
  }

  dr_base <- .solve_at(theta)
  if (is.null(dr_base)) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = "D40 Near-unit-root: baseline solve failed (steady state or BK)."
    ))
  }

  ## ---- 2. Spectral radius of the stable transition -----------------------
  ## ghx is n_endo x n_state.  Its COLUMNS correspond to state variables.
  ## The state-to-state transition is ghx[state_idx, ], which is n_state x n_state.
  ## That square matrix is the companion form whose spectral radius we want.
  ghx        <- as.matrix(dr_base$ghx)
  state_idx  <- dr_base$state_idx
  if (is.null(state_idx) || length(state_idx) == 0L) {
    ## Fallback: use all rows (purely static model edge case)
    state_idx <- seq_len(nrow(ghx))
  }
  T_mat <- ghx[state_idx, , drop = FALSE]   # n_state x n_state

  if (nrow(T_mat) == 0L || ncol(T_mat) == 0L) {
    ## No state variables — distance_to_unit is vacuously 1
    return(.make_result(
      result = list(
        spectral_radius      = 0,
        distance_to_unit     = 1,
        critical_eigenvalue  = complex(0),
        param_sensitivity    = data.frame(
          param           = param_names,
          d_rhomax_dparam = NA_real_,
          stringsAsFactors = FALSE
        ),
        equation_attribution = NULL
      ),
      pass    = TRUE,
      plots   = list(),
      summary = "D40 Near-unit-root: no state variables; spectral radius = 0."
    ))
  }

  .rhomax <- function(T) max(Mod(eigen(T, only.values = TRUE)$values))

  rho_base <- .rhomax(T_mat)
  ev_base  <- eigen(T_mat, only.values = TRUE)$values
  crit_ev  <- ev_base[which.max(Mod(ev_base))[1L]]

  ## ---- 3. Parameter sensitivity: d(rho_max)/d(theta_i) ------------------
  ## Central finite differences.  Re-solve at theta ± h_i, recompute rho_max.
  n_par  <- length(param_names)
  dsens  <- numeric(n_par)

  for (i in seq_len(n_par)) {
    pn <- param_names[i]
    if (is.null(theta[[pn]]) || !is.finite(theta[[pn]])) {
      dsens[i] <- NA_real_
      next
    }
    th_p <- theta;  th_p[[pn]] <- th_p[[pn]] + h
    th_m <- theta;  th_m[[pn]] <- th_m[[pn]] - h

    dr_p <- .solve_at(th_p)
    dr_m <- .solve_at(th_m)

    if (is.null(dr_p) || is.null(dr_m)) {
      dsens[i] <- NA_real_
      next
    }

    ## Rebuild T_mat at the perturbed parameters (state_idx is fixed from
    ## the baseline solve — index into the perturbed ghx the same way).
    T_p <- as.matrix(dr_p$ghx)[state_idx, , drop = FALSE]
    T_m <- as.matrix(dr_m$ghx)[state_idx, , drop = FALSE]

    rho_p <- .rhomax(T_p)
    rho_m <- .rhomax(T_m)
    dsens[i] <- (rho_p - rho_m) / (2 * h)
  }

  param_sens_df <- data.frame(
    param           = param_names,
    d_rhomax_dparam = dsens,
    stringsAsFactors = FALSE
  )
  param_sens_df <- param_sens_df[order(abs(dsens), decreasing = TRUE,
                                        na.last = TRUE), ]
  rownames(param_sens_df) <- NULL

  ## ---- 4. Equation attribution (best-effort) -----------------------------
  ## Strategy: the critical eigenvalue of T_mat has a left eigenvector v
  ## (or equivalently a right eigenvector w of T_mat').  The columns of ghx
  ## corresponding to the state variables span the structural sensitivity.
  ## We use the column norms of f_zero (the static Jacobian from
  ## extract_system_matrices) projected onto the critical eigendirection to
  ## rank equations.
  ##
  ## If extract_system_matrices is unavailable or the compiled model is not
  ## provided, skip gracefully.
  equation_attribution <- tryCatch({
    ## Get structural Jacobian at baseline.
    ## extract_system_matrices(compiled, ss, params) -- compiled first.
    if (is.null(compiled)) stop("compiled is NULL; cannot extract system matrices")
    params_full <- model$param_values
    for (nm in names(theta)) params_full[[nm]] <- theta[[nm]]
    sys <- extract_system_matrices(compiled,
                                   ss     = dr_base$ys,
                                   params = params_full)

    ## Right eigenvector of T_mat corresponding to critical eigenvalue.
    ## (Real part for a complex conjugate pair.)
    ev_full <- eigen(T_mat, only.values = FALSE)
    crit_idx <- which.max(Mod(ev_full$values))[1L]
    w_crit <- Re(ev_full$vectors[, crit_idx])  # length n_state

    ## Embed w_crit back into the full n_endo space via state_idx
    w_full <- numeric(nrow(ghx))
    w_full[state_idx] <- w_crit

    ## Project columns of f_zero onto w_full: ||f_zero_col|| * |cos(angle)|
    ## columns of f_zero correspond to endogenous variables (current-period)
    fz <- sys$f_zero
    col_proj <- apply(fz, 2L, function(col) {
      denom <- sqrt(sum(w_full^2) * sum(col^2))
      if (!is.finite(denom) || denom < .Machine$double.eps) return(0)
      abs(sum(w_full * col)) / denom
    })

    ## Row (equation) scores: max projection from any column
    row_proj <- apply(fz, 1L, function(row) {
      denom <- sqrt(sum(w_full^2) * sum(row^2))
      if (!is.finite(denom) || denom < .Machine$double.eps) return(0)
      abs(sum(w_full * row)) / denom
    })

    ## Rows of f_zero (after eq_to_decl reordering) correspond to
    ## equations in declaration-variable order; use endo_names as labels.
    eq_names <- sys$endo_names %||% paste0("eq_", seq_len(nrow(fz)))
    eq_df <- data.frame(
      equation    = eq_names,
      score       = row_proj,
      stringsAsFactors = FALSE
    )
    eq_df <- eq_df[order(eq_df$score, decreasing = TRUE), ]
    rownames(eq_df) <- NULL

    list(
      equation_scores    = eq_df,
      critical_eigenvec  = w_crit,
      note = paste0(
        "Scores are cosine similarity between the structural-Jacobian row ",
        "and the critical eigenvector. Higher = more responsible for the ",
        "near-unit-root. Best-effort: exact attribution requires model-specific ",
        "interpretation."
      )
    )
  }, error = function(e) {
    list(note = sprintf("Equation attribution unavailable: %s",
                        conditionMessage(e)))
  })

  ## ---- 5. Build summary --------------------------------------------------
  top_par <- if (nrow(param_sens_df) > 0 && !is.na(param_sens_df$d_rhomax_dparam[1L])) {
    sprintf("%s (d_rho/d_theta = %.4f)",
            param_sens_df$param[1L],
            param_sens_df$d_rhomax_dparam[1L])
  } else {
    "unknown"
  }
  top3 <- head(param_sens_df[!is.na(param_sens_df$d_rhomax_dparam), ], 3L)

  summary_text <- sprintf(
    paste0(
      "D40 Near-unit-root proximity: spectral_radius = %.6f, ",
      "distance_to_unit = %.6f, critical_eigenvalue = %s. ",
      "Top driver: %s. ",
      "(%d parameters differentiated%s)"
    ),
    rho_base,
    1 - rho_base,
    format(crit_ev),
    top_par,
    n_par,
    if (!is.null(equation_attribution$equation_scores))
      "; equation attribution available"
    else
      paste0("; equation attribution: ",
             equation_attribution$note %||% "unavailable")
  )

  llm_text <- paste(c(
    sprintf("D40 | Near-unit-root proximity | INFO"),
    sprintf("  spectral_radius=%.6f  distance_to_unit=%.6f", rho_base, 1 - rho_base),
    sprintf("  critical_eigenvalue=%s", format(crit_ev)),
    if (nrow(top3) > 0)
      sprintf("  top_drivers: %s",
              paste(sprintf("%s=%.4f", top3$param, top3$d_rhomax_dparam),
                    collapse = ", ")),
    sprintf("  action: %s",
            if (1 - rho_base < 0.05)
              "Model is near-unit-root (distance < 0.05). Check parameters with large positive d_rho/d_theta."
            else if (1 - rho_base < 0.20)
              "Moderate persistence (distance < 0.20). Monitor during estimation."
            else
              "Model is well away from unit circle. No immediate concern.")
  ), collapse = "\n")

  .make_result(
    result = list(
      spectral_radius      = rho_base,
      distance_to_unit     = 1 - rho_base,
      critical_eigenvalue  = crit_ev,
      param_sensitivity    = param_sens_df,
      equation_attribution = equation_attribution
    ),
    pass    = NA,   # informational diagnostic
    plots   = list(),
    summary = summary_text,
    llm_summary = llm_text
  )
}
