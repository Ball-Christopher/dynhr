## R/diag-pre-d0-equation-rank.R
## --------------------------------------------------------------------------
## D0: static equation-system rank pre-flight (redundant-equation check).
##
## The steady-state system F_static(y; theta) = 0 has an isolated, locally
## unique solution only if its Jacobian dF_static/dy has full column rank at the
## steady state. A rank-deficient static Jacobian means the steady-state
## equations are locally linearly dependent -- redundant equations or a
## continuum of steady states -- so the model is not locally well-posed and
## perturbation / estimation results are unreliable. This is the static analogue
## of Dynare's `model_diagnostics` rank check.
##
## Surfaced as run_all_diagnostics(...)$d0. Oracle: the numerical rank of the
## static Jacobian (svd / qr); cross-checked in test-diag-d0-equation-rank.R
## against a deliberately rank-deficient model and a well-posed one.
## --------------------------------------------------------------------------

#' D0 static equation-system rank (redundant-equation pre-flight)
#'
#' @param model    A \code{dynhr_mod}.
#' @param compiled A \code{dynhr_compiled} (provides \code{$static}).
#' @param params   Named numeric parameter vector.
#' @param ss       Steady-state values: a named numeric vector over the
#'   endogenous variables, or a list carrying \code{$values}.
#' @param tol      Relative singular-value tolerance for the rank cut
#'   (default \code{1e-8}).
#' @return A \code{dynhr_diagnostic}. \code{pass = TRUE} iff the static Jacobian
#'   is square and full column rank.
#' @noRd
d0_equation_rank <- function(model, compiled, params, ss, tol = 1e-8) {
  static <- if (!is.null(compiled)) compiled$static else NULL
  if (is.null(static) || is.null(static$jacobian_fn))
    return(.make_result(pass = NA,
      summary = "D0 equation rank: no compiled static model available."))

  endo   <- static$endo_names
  n_endo <- length(endo)

  ss_vals <- if (is.list(ss)) (ss$values %||% ss$ss %||% unlist(ss, use.names = TRUE))
             else ss
  if (is.null(names(ss_vals)) && length(ss_vals) == n_endo)
    names(ss_vals) <- endo
  y <- setNames(as.numeric(ss_vals[endo]), endo)
  if (anyNA(y))
    return(.make_result(pass = NA, errored = TRUE,
      summary = "D0 equation rank: steady-state values missing for some endogenous variables."))

  exo <- (if (!is.null(compiled$dynamic)) compiled$dynamic$exo_names else NULL) %||%
         model$varexo_names %||% character(0)
  x <- setNames(rep(0, length(exo)), exo)
  if (is.list(params)) params <- unlist(params)

  J <- tryCatch(as.matrix(static$jacobian_fn(y, x, params, y)),
                error = function(e) NULL)
  if (is.null(J) || !all(is.finite(J)))
    return(.make_result(pass = NA, errored = TRUE,
      summary = "D0 equation rank: static Jacobian evaluation failed / non-finite."))

  n_eq   <- nrow(J)
  sv_d   <- svd(J)
  sv     <- sv_d$d
  smax   <- if (length(sv)) sv[1] else 0
  thr    <- tol * smax
  rank   <- if (smax > 0) sum(sv > thr) else 0L
  deficiency <- n_endo - rank
  cond   <- if (rank > 0L) smax / min(sv[sv > thr]) else Inf
  square <- (n_eq == n_endo)

  ## Variables most involved in the (right) null space -- those not separately
  ## pinned by the steady-state equations.
  suspect <- character(0)
  if (deficiency > 0L && !is.null(sv_d$v) && ncol(sv_d$v) >= 1L) {
    null_cols <- which(sv <= thr)
    if (length(null_cols)) {
      w <- rowSums(abs(sv_d$v[, null_cols, drop = FALSE]))
      k <- min(n_endo, deficiency + 1L)
      suspect <- endo[order(w, decreasing = TRUE)[seq_len(k)]]
    }
  }

  pass <- (deficiency == 0L) && square
  summary <- if (pass) {
    sprintf("D0 equation rank: full column rank %d/%d (condition number %.2e).",
            rank, n_endo, cond)
  } else if (!square) {
    sprintf(paste0("D0 equation rank: non-square static system (%d equations, ",
                   "%d endogenous); numerical rank %d."),
            n_eq, n_endo, rank)
  } else {
    sprintf(paste0("D0 equation rank: RANK-DEFICIENT %d/%d (deficiency %d). The ",
                   "steady-state equations are locally linearly dependent; ",
                   "variables not separately pinned: %s."),
            rank, n_endo, deficiency, paste(suspect, collapse = ", "))
  }

  .make_result(
    result = list(rank = rank, n_endo = n_endo, n_eq = n_eq,
                  deficiency = deficiency, square = square,
                  condition_number = cond, singular_values = sv,
                  suspect_variables = suspect, tol = tol),
    pass    = pass,
    summary = summary)
}
