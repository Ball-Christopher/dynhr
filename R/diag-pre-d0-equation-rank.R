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
#' The rank is taken on the EQUILIBRATED Jacobian \eqn{D_r J D_c} (every row,
#' then every column, scaled to unit max-abs). Diagonal scaling leaves the
#' exact rank unchanged but removes the spurious "deficiency" a relative
#' singular-value cut reports when equations are written in very different
#' units (e.g. \code{1e9*(x - 1) = e}).
#'
#' @param model    A \code{dynhr_mod}.
#' @param compiled A \code{dynhr_compiled} (provides \code{$static}).
#' @param params   Named numeric parameter vector. May be a SUBSET (e.g. a
#'   posterior mean over the estimated parameters only); missing parameters are
#'   filled from \code{model$param_values}.
#' @param ss       Steady-state values: a named numeric vector over the
#'   endogenous variables, or a list carrying \code{$values}.
#' @param tol      Relative singular-value tolerance for the rank cut on the
#'   equilibrated Jacobian (default \code{1e-8}).
#' @param dr       Optional solved decision rule (\code{solve_perturbation()}
#'   output, or anything carrying \code{$eigenvalues}). Used ONLY to decide
#'   whether a singular static Jacobian is explained by a unit root.
#' @param unit_tol Modulus tolerance for calling an eigenvalue of the solved
#'   system a unit root (default \code{1e-6}, matching Dynare's
#'   \code{model_diagnostics.m} and \code{kalman_filter(lik_init = "auto")}).
#'
#' @section Singular Jacobian: unit root versus redundant equation:
#' A singular static Jacobian has two very different causes, and the rank test
#' alone cannot tell them apart. Dynare's \code{model_diagnostics.m} resolves
#' it by looking at the SOLVED system: if some eigenvalue has modulus within
#' \code{unit_tol} of 1, the singularity "is probably due to the presence of a
#' unit root", and "if the model is actually supposed to feature unit root
#' behaviour, such a warning is expected". D0 follows that logic exactly:
#' \describe{
#'   \item{WARN}{rank-deficient AND \code{dr} shows a unit root -- singular by
#'     construction; the steady state is a continuum along the unit-root
#'     direction, which is a modelling choice, not an error.}
#'   \item{FAIL}{rank-deficient with no unit-root explanation (or no \code{dr}
#'     supplied to check) -- a redundant equation, usually paired with a
#'     missing one.}
#' }
#' @return A \code{dynhr_diagnostic}. \code{pass = TRUE} iff the static Jacobian
#'   is square and full column rank, or the deficiency is unit-root-explained
#'   (then \code{warn = TRUE}). \code{$result} carries \code{rank},
#'   \code{deficiency}, \code{condition_number} and \code{singular_values} (of
#'   the equilibrated Jacobian), \code{singular_values_raw} (unscaled),
#'   \code{suspect_variables} / \code{suspect_equations} (support of the right /
#'   left null space), \code{redundant_equations}, \code{equation_relations}
#'   and \code{variable_relations} (see below), \code{equation_tags} (the
#'   \code{[name=...]} tag of each \code{eqN}, or NA) and
#'   \code{max_abs_residual} (static residual at the evaluation point; large
#'   means \code{ss} is not a steady state at \code{params}).
#'
#' @section Per-relation report (Dynare 7 model_diagnostics):
#' Dynare 7's \code{model_diagnostics} lists the collinear equations and
#' variables one RELATION at a time, from MATLAB's rational (row-echelon)
#' null-space basis, trimmed to the smallest support that still annihilates
#' the Jacobian. D0 does the same: each element of \code{equation_relations}
#' is \code{list(redundant, coefficients, text)}, where \code{redundant} is the
#' latest equation (in model order) of the relation, \code{coefficients} are
#' the weights (original units, \code{redundant} first with weight 1) with
#' \code{sum_i coefficients[i] * dF_i/dy = 0} at the steady state, and
#' \code{text} reads e.g. \code{"eq3 = 2*eq1 + 3*eq2"}. There is one relation
#' per dimension of the left null space; \code{redundant_equations} collects
#' their \code{redundant} labels. \code{variable_relations} is the analogue
#' for the right null space (directions along which the static residuals do
#' not change). Only the pivot choice is a convention: dropping any single
#' equation with a non-zero weight restores the relation's rank.
#' @noRd
d0_equation_rank <- function(model, compiled, params, ss, tol = 1e-8,
                             dr = NULL, unit_tol = 1e-6) {
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

  ## The generated closures index params BY NAME, so a subset vector (the
  ## posterior-mean path passes only the estimated parameters) evaluates the
  ## missing ones as NA. Overlay the supplied values on the calibration.
  if (is.list(params)) params <- unlist(params)
  base <- model$param_values
  if (is.list(base)) base <- unlist(base)
  if (!is.null(base) && !is.null(names(params))) {
    common <- intersect(names(params), names(base))
    base[common] <- params[common]
    params <- base
  }

  J <- tryCatch(as.matrix(static$jacobian_fn(y, x, params, y)),
                error = function(e) .dynhr_reraise_bug(e, NULL))
  if (is.null(J) || !all(is.finite(J)))
    return(.make_result(pass = NA, errored = TRUE,
      summary = "D0 equation rank: static Jacobian evaluation failed / non-finite."))

  max_res <- if (is.function(static$residuals_fn))
    max(abs(static$residuals_fn(y, x, params, y)), 0) else NA_real_

  n_eq   <- nrow(J)
  square <- (n_eq == n_endo)

  ## Equilibrate: rows, then columns, to unit max-abs (zero lines left as is).
  rs <- apply(abs(J), 1L, max); rs[!(rs > 0)] <- 1
  Js <- J / rs
  cs <- apply(abs(Js), 2L, max); cs[!(cs > 0)] <- 1
  Js <- sweep(Js, 2L, cs, "/")

  sv_raw <- svd(J, nu = 0L, nv = 0L)$d
  sv_d   <- svd(Js, nu = n_eq, nv = n_endo)
  sv     <- sv_d$d
  smax   <- if (length(sv)) sv[1] else 0
  thr    <- tol * smax
  rank   <- if (smax > 0) sum(sv > thr) else 0L
  deficiency <- n_endo - rank
  cond   <- if (rank > 0L) smax / sv[rank] else Inf

  ## Support of the right null space (variables not separately pinned) and of
  ## the left null space (equations that are linear combinations of others).
  ## Null-space basis vectors are unit-norm, so 1e-6 separates structural
  ## non-zeros from rounding noise.
  support <- function(B, labels) {
    if (is.null(B) || ncol(B) == 0L) return(character(0))
    w <- rowSums(abs(B))
    keep <- which(w > 1e-6)
    labels[keep[order(w[keep], decreasing = TRUE)]]
  }
  eq_labels <- paste0("eq", seq_len(n_eq))
  eq_tags <- .d0_equation_tags(model, n_eq)
  eq_display <- ifelse(is.na(eq_tags), eq_labels,
                       paste0(eq_labels, " [", eq_tags, "]"))
  names(eq_tags) <- eq_labels
  null_eq  <- if (rank < n_eq)
    sv_d$u[, seq.int(rank + 1L, n_eq), drop = FALSE] else NULL
  null_var <- if (rank < n_endo)
    sv_d$v[, seq.int(rank + 1L, n_endo), drop = FALSE] else NULL
  suspect    <- support(null_var, endo)
  suspect_eq <- support(null_eq, eq_labels)

  ## Dynare 7 model_diagnostics: one report PER collinear relation, not just
  ## the union support. `rs` / `cs` map the equilibrated relations back to the
  ## original equation / variable units.
  eq_rel  <- .d0_relations(null_eq,  Js,    1 / rs, eq_labels, eq_display)
  var_rel <- .d0_relations(null_var, t(Js), 1 / cs, endo, endo, kind = "variable")
  redundant_eq <- vapply(eq_rel, `[[`, "", "redundant")

  full_rank <- (deficiency == 0L) && square

  ## Dynare model_diagnostics.m: a singular static Jacobian whose solved system
  ## carries an eigenvalue of modulus 1 is singular BY CONSTRUCTION. We cannot
  ## distinguish that from a redundant equation from the Jacobian alone, so the
  ## disambiguation needs `dr`; without it we keep the conservative FAIL.
  ev <- if (!is.null(dr)) (dr$eigenvalues %||% dr$dr$eigenvalues) else NULL
  n_unit_root <- if (is.null(ev)) NA_integer_
                 else sum(abs(Mod(ev) - 1) < unit_tol)
  unit_root_explained <- isTRUE(n_unit_root > 0L)

  pass <- full_rank || (square && unit_root_explained)
  warn <- !full_rank && square && unit_root_explained

  summary <- if (full_rank) {
    sprintf("D0 equation rank: full column rank %d/%d (condition number %.2e, equilibrated).",
            rank, n_endo, cond)
  } else if (!square) {
    sprintf(paste0("D0 equation rank: non-square static system (%d equations, ",
                   "%d endogenous); numerical rank %d."),
            n_eq, n_endo, rank)
  } else if (warn) {
    sprintf(paste0("D0 equation rank: rank-deficient %d/%d (deficiency %d), but ",
                   "the solved system has %d eigenvalue(s) of modulus 1 ",
                   "(abs(%s - 1) < %g): the static Jacobian is SINGULAR BY ",
                   "CONSTRUCTION, as expected for a unit-root model (Dynare ",
                   "model_diagnostics reports the same). The steady state is a ",
                   "continuum along the unit-root direction(s). Variables not ",
                   "separately pinned: %s."),
            rank, n_endo, deficiency, n_unit_root, "|lambda|", unit_tol,
            paste(suspect, collapse = ", "))
  } else {
    sprintf(paste0("D0 equation rank: RANK-DEFICIENT %d/%d (deficiency %d). The ",
                   "static Jacobian is singular and %s: a redundant equation ",
                   "(often paired with a missing one) or a continuum of steady ",
                   "states. Collinear equations: %s; variables not separately ",
                   "pinned: %s."),
            rank, n_endo, deficiency,
            if (is.null(ev)) "no decision rule was supplied to test for a unit root"
            else "the solved system has NO unit root to explain it",
            paste(suspect_eq, collapse = ", "), paste(suspect, collapse = ", "))
  }
  if (length(eq_rel))
    summary <- paste0(summary, sprintf(
      " Dependent equation relation(s) at the steady state: %s.",
      paste(vapply(eq_rel, `[[`, "", "text"), collapse = "; ")))
  if (is.finite(max_res) && max_res > 1e-6)
    summary <- paste0(summary, sprintf(
      " Note: the evaluation point is not a steady state at these parameters (max |static residual| = %.2e).",
      max_res))

  .make_result(
    result = list(rank = rank, n_endo = n_endo, n_eq = n_eq,
                  deficiency = deficiency, square = square,
                  condition_number = cond, singular_values = sv,
                  singular_values_raw = sv_raw,
                  suspect_variables = suspect, suspect_equations = suspect_eq,
                  redundant_equations = redundant_eq,
                  equation_relations = eq_rel,
                  variable_relations = var_rel,
                  equation_tags = eq_tags,
                  max_abs_residual = max_res, tol = tol,
                  n_unit_root = n_unit_root, unit_tol = unit_tol,
                  unit_root_explained = unit_root_explained),
    pass    = pass,
    warn    = warn,
    summary = summary)
}


#' Equation `name` tags aligned with the static-model rows
#'
#' Returns the `[name=...]` tag of each model equation, or NA. When the number
#' of parsed equations does not match the static Jacobian's rows the alignment
#' is unknown, so every tag is NA and D0 falls back to `eqN` labels.
#' @noRd
.d0_equation_tags <- function(model, n_eq) {
  eqs <- if (is.list(model)) model$equations else NULL
  if (is.null(eqs) || length(eqs) != n_eq) return(rep(NA_character_, n_eq))
  vapply(eqs, function(e) {
    tg <- if (is.list(e)) e$tag else NULL
    if (is.character(tg) && length(tg) == 1L && !is.na(tg) && nzchar(tg)) tg
    else NA_character_
  }, character(1))
}


#' One reduced relation per null-space dimension (Dynare model_diagnostics)
#'
#' Dynare 7's `model_diagnostics.m` prints the static Jacobian's collinear
#' variables and collinear equations RELATION BY RELATION, using MATLAB's
#' rational null-space basis (`null(J', "rational")`, i.e. the reduced
#' row-echelon form) and trimming each relation to the smallest support
#' `abs(n) > 10^-j`, `j = 1..10`, that still annihilates the Jacobian to
#' `1e-6`. This does the same on the equilibrated Jacobian: Gauss-Jordan on the
#' transposed orthonormal null basis, pivoting from the LAST label backwards,
#' so each relation expresses one "redundant" entry (the pivot, latest in model
#' order) as a combination of earlier non-pivot entries. Coefficients are
#' returned in the ORIGINAL units (`back_scale` undoes the equilibration),
#' normalised so the pivot carries coefficient 1.
#'
#' @param N Orthonormal null basis (n x d), or NULL.
#' @param M Matrix with `crossprod(M[k, ], b[k]) = 0` for a relation `b`
#'   (the equilibrated Jacobian for equations, its transpose for variables).
#' @param back_scale Length-n multiplier mapping equilibrated to raw units.
#' @param labels,display Machine labels and display labels (length n).
#' @param kind "equation" (text "eq3 = 2*eq2") or "variable" (a direction).
#' @return A list of relations, each `list(redundant, coefficients, text)`;
#'   `coefficients` is named by `labels`, pivot first (coefficient 1).
#' @noRd
.d0_relations <- function(N, M, back_scale, labels, display,
                          kind = c("equation", "variable")) {
  kind <- match.arg(kind)
  if (is.null(N) || ncol(N) == 0L) return(list())
  B <- t(N)
  d <- nrow(B); n <- ncol(B)
  bmax <- max(abs(B))
  if (!(bmax > 0)) return(list())
  used  <- logical(d)
  pivot <- integer(d)
  for (col in rev(seq_len(n))) {
    if (all(used)) break
    cand <- which(!used)
    r <- cand[which.max(abs(B[cand, col]))]
    if (abs(B[r, col]) <= 1e-6 * bmax) next
    B[r, ] <- B[r, ] / B[r, col]
    for (o in setdiff(seq_len(d), r))
      B[o, ] <- B[o, ] - B[o, col] * B[r, ]
    used[r]  <- TRUE
    pivot[r] <- col
  }
  rels <- lapply(which(used), function(r) {
    p  <- pivot[r]
    bn <- B[r, ] / max(abs(B[r, ]))
    keep <- seq_len(n)
    for (j in 1:10) {
      k <- sort(union(p, which(abs(bn) > 10^-j)))
      if (max(abs(crossprod(M[k, , drop = FALSE], bn[k]))) < 1e-6) {
        keep <- k
        break
      }
    }
    cf <- B[r, keep] * back_scale[keep]
    cf <- cf / cf[keep == p]
    others <- setdiff(keep, p)
    cf <- c(cf[keep == p], cf[keep != p])
    names(cf) <- labels[c(p, others)]
    list(redundant = labels[p], coefficients = cf,
         text = if (kind == "equation")
                  .d0_relation_text(display[p], -cf[-1L], display[others])
                else sprintf("static residuals unchanged along d(%s) = (%s)",
                             paste(names(cf), collapse = ", "),
                             paste(as.character(signif(cf, 6)), collapse = ", ")))
  })
  rels[order(pivot[used])]
}


#' Render one relation: "eq3 = 2*eq1 + 3*eq2" (or a zero-row note)
#' @noRd
.d0_relation_text <- function(lhs, coef, rhs) {
  if (!length(coef))
    return(sprintf("%s has an identically zero steady-state Jacobian row", lhs))
  term <- function(a, lab) {
    mag <- abs(a)
    if (abs(mag - 1) < 1e-10) lab else paste0(format(signif(mag, 6)), "*", lab)
  }
  out <- paste0(if (coef[1] < 0) "-" else "", term(coef[1], rhs[1]))
  for (i in seq_along(coef)[-1L])
    out <- paste0(out, if (coef[i] < 0) " - " else " + ", term(coef[i], rhs[i]))
  paste0(lhs, " = ", out)
}
