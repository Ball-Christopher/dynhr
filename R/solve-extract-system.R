## R/solve-extract-system.R
## --------------------------------------------------------------------------
## System-matrix extraction from a compiled dynhr_mod: builds the Jacobian
## blocks f_minus, f_zero, f_plus, f_exo at a given steady state and
## parameter vector. cache_system_structure() pre-computes timing-column
## indices so extract_system_matrices_fast() can skip regex+classification
## inside MCMC loops.
##
## Phase-1 split from perturbation-monolith.R (no logic changes).
## --------------------------------------------------------------------------

# =====================================================================
# System matrix extraction
# =====================================================================

## Structural lag/lead masks from the lead_lag_incidence matrix.
##
## A variable is a state (predetermined/mixed) iff its lag (t-1) timing is
## STRUCTURALLY present in some equation -- i.e. it has a nonzero entry in the
## t-1 row of the lead_lag_incidence -- regardless of the numerical value of the
## evaluated Jacobian at the current steady state.  Likewise for leads (t+1).
##
## This matters when a structurally-present lag carries a Jacobian coefficient
## that vanishes at the chosen steady state.  The canonical case is a Ramsey
## augmented system: a multiplier MULT_j is zero at the deterministic SS, so any
## term MULT_j(-1) * (... x(-1) ...) contributes a zero column to f_minus, even
## though x is genuinely backward-looking.  A purely numerical classification
## (apply(abs(f_minus), 2, max) > tol) would then drop x from the state set,
## under-dimensioning the perturbation state vector relative to Dynare (which
## classifies structurally).  On NK_optimal_policy this dropped the policy
## instrument `i` -- its only lag occurrence is MULT_3(-1)*(...i(-1)...) and
## MULT_3 = 0 at SS -- leaving 6 states instead of Dynare's 7 {i,u,C,Y,pi,
## MULT_3,MULT_5}.  Models where rho = 0 makes a shock iid (e.g. ireland_2004
## `zt` with rho_z = 0) keep that variable as a (harmless, zero-root) state too,
## matching Dynare, which does not constant-fold the coefficient away.
##
## @param lli         lead_lag_incidence matrix (rows = timings, cols = endo).
## @param row_labels  integer timing per LLI row (-1, 0, +1, ...).
## @return list(has_lag, has_lead) logical vectors over the LLI columns.
.structural_lag_lead <- function(lli, row_labels) {
  lag_rows  <- which(row_labels < 0L)
  lead_rows <- which(row_labels > 0L)
  n_col <- ncol(lli)
  has_lag  <- if (length(lag_rows))
    apply(lli[lag_rows, , drop = FALSE],  2L, function(x) any(x > 0L))
  else rep(FALSE, n_col)
  has_lead <- if (length(lead_rows))
    apply(lli[lead_rows, , drop = FALSE], 2L, function(x) any(x > 0L))
  else rep(FALSE, n_col)
  list(has_lag = has_lag, has_lead = has_lead)
}

## Shared helper: walk an equation LHS AST and return the endogenous variable
## name, or NULL if none can be identified.
##
## Handles:
##   variable  -> the name directly
##   binop     -> left or right sub-node (with special case "1/q" -> right)
##   unaryop   -> the operand (e.g. "-x" -> "x")
##
## NOTE: the node type is "unaryop" (from ast_unaryop() in parse-equations.R).
## An earlier typo spelled it "uniop" in a few dead branches; those are fixed
## here and in every caller.
.lhs_endo_var <- function(ast) {
  if (is.null(ast)) return(NULL)
  if (is.list(ast) && identical(ast$type, "variable")) return(ast$name)
  if (is.list(ast) && identical(ast$type, "binop")) {
    if (ast$op == "/" && is.list(ast$right) && ast$right$type == "variable")
      return(ast$right$name)
    return(.lhs_endo_var(ast$left) %||% .lhs_endo_var(ast$right))
  }
  if (is.list(ast) && identical(ast$type, "unaryop"))
    return(.lhs_endo_var(ast$operand))
  NULL
}

## Jacobian-based fallback for compound-LHS equation mapping (H5).
##
## When the first-pass LHS heuristic leaves some equations unmapped
## (eq_to_decl[i] == 0), use the static Jacobian columns (f_zero, one column
## per endogenous variable in declaration order) to find the best-fit variable
## for each unmapped equation:  for each unmapped equation i, assign the
## still-unassigned variable j whose |f_zero[i, j]| is largest.  This is a
## greedy bipartite matching that recovers the Dynare convention for compound-
## LHS equations like `uc*(1+phi*...)*exp(g) = ...` where the first-name
## heuristic picks the wrong variable.
##
## Arguments:
##   eq_to_decl   Integer vector (length n_eq): 0 for unmapped equations.
##   f_zero       n_eq x n_endo matrix: static (t=0) Jacobian.
##   endo         Character vector of endogenous variable names (declaration order).
## Returns:
##   Updated eq_to_decl (0 entries resolved where possible).
.jacobian_match_unmapped <- function(eq_to_decl, f_zero, endo) {
  n_endo <- length(endo)
  unassigned <- setdiff(seq_len(n_endo), eq_to_decl[eq_to_decl > 0L])
  unmapped   <- which(eq_to_decl == 0L)
  if (length(unmapped) == 0L || length(unassigned) == 0L) return(eq_to_decl)

  ## Greedy maximum-absolute-value matching: process unmapped equations in
  ## order; for each, pick the unassigned variable with the largest |f_zero|
  ## entry in that equation row.
  remaining <- unassigned
  for (i in unmapped) {
    if (length(remaining) == 0L) break
    ## Restrict to unassigned columns
    col_vals <- abs(f_zero[i, remaining, drop = TRUE])
    best <- which.max(col_vals)
    if (length(best) > 0L) {
      picked <- remaining[best]
      eq_to_decl[i] <- picked
      remaining <- remaining[remaining != picked]
    }
  }
  eq_to_decl
}

## TRUE when the compiled C++ Jacobian-tape evaluator is available and the user
## has not disabled the Rcpp backend (options(dynhr.use_rcpp = FALSE)).
.HAS_RCPP_TAPE <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("eval_jac_tape_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

## TRUE when the compiled C++ triplet-tape evaluator (for the Hessian2/Hessian3
## value-vector tapes) is available and the Rcpp backend is not disabled.
.HAS_RCPP_TRIPLET_TAPE <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("eval_triplet_tape_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

## Evaluate a Hessian value-vector "tape" (hess2_tape / hess3_tape) at the
## compound steady state via the C++ stack machine, returning the dense value
## vector v[k] (one per triplet) -- the same shape the interpreted hessian*_fn
## closures return. `dy_ss` is the .build_dy_ss_o2 vector in dyn_col_map order
## (== dy_keys_tape order), so unname() is positionally correct. params/ss are
## reindexed positionally to the tape's namespaces (mirrors the jac-tape dance
## in extract_system_matrices_fast). Returns NULL if the tape is unusable.
.eval_triplet_tape <- function(dyn, tape, dy_ss, params, ss) {
  if (is.null(tape) || !.HAS_RCPP_TRIPLET_TAPE()) return(NULL)
  ## Defensive: the tape indexes dy positionally by dyn_col_map row, so the
  ## supplied dy_ss MUST have exactly one entry per column-map row.
  if (length(dy_ss) != nrow(dyn$dyn_col_map)) return(NULL)
  par_pos <- params[dyn$jac_tape_param_names]
  ss_pos  <- ss[dyn$endo_names]; ss_pos[is.na(ss_pos)] <- 0
  eval_triplet_tape_cpp(unname(dy_ss), unname(par_pos), unname(ss_pos),
                        tape$op, tape$ia, tape$da, tape$expr_len)
}

#' Evaluate the dynamic Jacobian at steady state and partition into
#' f_minus (t-1), f_zero (t), f_plus (t+1), f_exo (shocks)
#'
#' If the symbolic Jacobian produces NaN/Inf for some columns, falls
#' back to central finite differences using the residual function for
#' those columns. This handles cases where the symbolic derivative is
#' undefined at steady state (e.g., 0/0 from L'Hopital limits).
#'
#' @param compiled dynhr_compiled
#' @param ss       Named numeric steady state vector
#' @param params   Named numeric parameter vector
#' @param center   Optional named numeric vector: the linearization point.
#'   When non-NULL, `dy` is built from `center` instead of `ss`. `ss` is
#'   still passed as the STEADY_STATE() oracle for models that reference it
#'   in expressions like `ss["y"]`. Default NULL => behavior identical to
#'   the classic (center = steady state) call.
#' @return List with f_minus, f_zero, f_plus, f_exo matrices and metadata
#' @noRd
extract_system_matrices <- function(compiled, ss, params, center = NULL) {
  dyn <- compiled$dynamic
  model <- compiled$model
  endo <- model$var_names
  exo  <- model$varexo_names
  n_eq <- dyn$n_eq
  n_endo <- length(endo)
  n_exo <- length(exo)
  lli <- model$lead_lag_incidence

  # Build the dy vector: all variables at all timings = expansion center.
  # When center is supplied, use it as the linearization point; ss is kept
  # only as the STEADY_STATE() oracle argument to jacobian_fn.
  lp <- if (!is.null(center)) center else ss   # linearization point
  dy <- numeric(0)
  for (k in seq_len(nrow(dyn$dyn_col_map))) {
    nm  <- dyn$dyn_col_map$name[k]
    ll  <- dyn$dyn_col_map$lead_lag[k]
    sfx <- if (ll == 0L) "__0"
    else if (ll > 0L) paste0("__p", ll)
    else paste0("__m", abs(ll))
    key <- paste0(nm, sfx)
    val <- if (nm %in% names(lp)) lp[[nm]]
    else if (nm %in% exo) 0
    else 0
    dy[key] <- val
  }

  # Evaluate dynamic Jacobian at the expansion center.
  # ss is passed as the STEADY_STATE() oracle (third arg) unchanged.
  J <- dyn$jacobian_fn(dy, params, ss)

  # ------------------------------------------------------------------
  # OccBin relax-regime row selection.
  # For models with occbin_constraints, n_eq > n_endo: the full Jacobian
  # has rows for both bind- and relax-equations.  For the reference (all-
  # slack) perturbation we keep only the relax-regime rows (regime 0).
  # This must happen BEFORE the eq_to_decl / f_zero_raw build so that
  # the downstream code (including .build_eq_to_decl) sees exactly n_endo
  # rows.  The occbin_parse_bind_relax tag scan is used to build a
  # reduced model stub (equations list only) so .build_eq_to_decl picks
  # the correct LHS variables from the relax equations.
  model_for_eq_to_decl <- model  # used below by .build_eq_to_decl
  if (inherits(compiled, "dynhr_compiled") &&
      !is.null(compiled[["occbin"]]) && n_eq > n_endo) {
    relax_rows_J <- compiled$occbin$regime_map[[1L]]$eq_indices
    J    <- J[relax_rows_J, , drop = FALSE]
    n_eq <- n_endo  # update local count for f_zero_raw allocation etc.
    # Build a reduced model stub so .build_eq_to_decl uses only the relax
    # equations (which have distinct LHS variables, unlike the bind equations
    # that shadow the same variables).
    model_for_eq_to_decl <- model
    model_for_eq_to_decl$equations <- model$equations[relax_rows_J]
  }

  # ------------------------------------------------------------------
  # Reorder Jacobian rows from model-equation order to
  # declaration-variable order.
  #
  # The equations in the .mod file may appear in any order (e.g.
  # grouped by economic topic), but the lead_lag_incidence matrix
  # and subsequent solver steps assume that row i of the Jacobian
  # corresponds to the equation for endogenous variable i (in
  # declaration order).  We determine the equation-to-variable
  # correspondence by matching each equation's LHS variable to
  # its declaration position.  When multiple equations share the
  # same LHS variable (e.g. equations 9 and 10 both involve `q`),
  # the first equation in declaration order takes precedence.
  # ------------------------------------------------------------------
  # Equation -> declaration-variable mapping.  Build the RAW static (t=0)
  # Jacobian (rows in ORIGINAL equation order, BEFORE any reordering) so the
  # shared .build_eq_to_decl() can apply its positional-first, override-only-
  # when-provably-invalid refinement (H5 compound-LHS).  Passing the raw
  # f_zero is essential: a reordered f_zero would misalign rows i.
  row_labels_early <- vapply(rownames(lli), function(rn) {
    rn <- trimws(rn)
    if (rn == "t") return(0L)
    m2 <- regmatches(rn, regexec("^t([+-]?\\d+)$", rn))[[1]]
    if (length(m2) == 2L) as.integer(m2[2L]) else 0L
  }, integer(1L), USE.NAMES = FALSE)
  f_zero_raw <- matrix(0, nrow = n_eq, ncol = n_endo)
  colnames(f_zero_raw) <- endo
  for (j in seq_along(endo)) {
    for (ri in seq_len(nrow(lli))) {
      col_idx_j <- lli[ri, j]
      if (col_idx_j > 0L && row_labels_early[ri] == 0L) {
        f_zero_raw[, j] <- J[, col_idx_j]
        break
      }
    }
  }
  eq_to_decl <- .build_eq_to_decl(model_for_eq_to_decl, f_zero = f_zero_raw)

  if (all(eq_to_decl > 0) && !identical(eq_to_decl, seq_len(n_eq))) {
    J <- J[order(eq_to_decl), , drop = FALSE]
  }

  # ------------------------------------------------------------------
  # Numerical Jacobian fallback for non-finite columns
  # ------------------------------------------------------------------
  # If the symbolic Jacobian has NaN/Inf (typically from expressions
  # like x/x whose derivative is 0/0 at x=0), recompute those columns
  # via central finite differences using the residual function.
  if (any(!is.finite(J))) {
    bad_jac <- which(!is.finite(J), arr.ind = TRUE)
    bad_cols <- sort(unique(bad_jac[, 2]))

    # Build human-readable names for warning message
    col_names <- character(0)
    for (bc in bad_cols) {
      if (bc <= nrow(dyn$dyn_col_map)) {
        col_names <- c(col_names, paste0(dyn$dyn_col_map$name[bc],
                                         "(", dyn$dyn_col_map$lead_lag[bc], ")"))
      } else {
        eidx <- bc - dyn$n_dyn_cols
        if (eidx > 0 && eidx <= n_exo) col_names <- c(col_names, exo[eidx])
      }
    }

    # Try numerical fallback using residual function
    has_resid <- !is.null(dyn$residuals_fn)
    if (has_resid) {
      # Base residual at steady state
      f0 <- dyn$residuals_fn(dy, params, ss)
    } else {
      f0 <- NULL
    }

    if (!is.null(f0) && all(is.finite(f0))) {
      # Central finite differences for bad columns
      n_fixed <- 0L
      for (bc in bad_cols) {
        dy_p <- dy
        dy_m <- dy
        h <- max(1e-6, abs(dy[bc]) * 1e-6)
        dy_p[bc] <- dy[bc] + h
        dy_m[bc] <- dy[bc] - h
        fp <- dyn$residuals_fn(dy_p, params, ss)
        fm <- dyn$residuals_fn(dy_m, params, ss)
        if (!is.null(fp) && !is.null(fm) &&
            all(is.finite(fp)) && all(is.finite(fm))) {
          J[, bc] <- (fp - fm) / (2 * h)
          n_fixed <- n_fixed + 1L
        } else {
          # Finite difference also failed; zero out
          J[!is.finite(J[, bc]), bc] <- 0
        }
      }
      warning(sprintf(
        paste0("Non-finite values in symbolic Jacobian at steady state. ",
               "Affected columns: %s. ",
               "Repaired %d/%d columns via numerical finite differences."),
        paste(col_names, collapse = ", "), n_fixed, length(bad_cols)))
    } else if (length(bad_cols) > 2L) {
      # PERVASIVE non-finite Jacobian with no fallback: the steady state is
      # degenerate (commonly: parameters computed in an external steadystate.m
      # left at placeholder values, or solve_steady did not truly converge).
      # Zeroing would feed a corrupt Jacobian into the order-2/3 Kronecker /
      # Sylvester solves, which then either error confusingly ("system is
      # exactly singular") or run effectively unbounded (M25, the
      # Basu_Bundick_2017 order-3 "hang"). Fail loudly and early instead.
      stop(sprintf(
        paste0("Dynamic Jacobian is non-finite at the steady state in %d ",
               "columns (%s) and no residual function is available for a ",
               "numerical fallback. The steady state is degenerate -- check ",
               "the calibration and supply any parameters computed in an ",
               "external steadystate.m via inject_params() (for a superset vector) ",
               "or set_param_values() (for the exact model-param set). Cannot extract ",
               "a usable system."),
        length(bad_cols),
        paste(utils::head(col_names, 10L), collapse = ", ")))
    } else {
      # Isolated non-finite (e.g. a single 0/0 L'Hopital cell) with no
      # fallback: zero it and warn, as before.
      warning(sprintf(
        paste0("Non-finite values in dynamic Jacobian at steady state. ",
               "Affected columns: %s. ",
               "No residual function available for numerical fallback. ",
               "Replacing with 0."),
        paste(col_names, collapse = ", ")))
      J[!is.finite(J)] <- 0
    }
  }

  # Parse row labels to get timing for each column
  row_labels <- vapply(rownames(lli), function(rn) {
    rn <- trimws(rn)
    if (rn == "t") return(0L)
    m <- regmatches(rn, regexec("^t([+-]?\\d+)$", rn))[[1]]
    if (length(m) == 2) return(as.integer(m[2]))
    0L
  }, integer(1), USE.NAMES = FALSE)

  # Partition Jacobian columns by timing
  f_minus <- matrix(0, nrow = n_eq, ncol = n_endo)
  f_zero  <- matrix(0, nrow = n_eq, ncol = n_endo)
  f_plus  <- matrix(0, nrow = n_eq, ncol = n_endo)
  f_exo   <- matrix(0, nrow = n_eq, ncol = n_exo)

  colnames(f_minus) <- endo; colnames(f_zero) <- endo
  colnames(f_plus) <- endo; colnames(f_exo) <- exo

  # Endogenous columns from lead_lag_incidence
  for (j in seq_along(endo)) {
    for (i in seq_len(nrow(lli))) {
      col_idx <- lli[i, j]
      if (col_idx > 0) {
        timing <- row_labels[i]
        if (timing == -1L)      f_minus[, j] <- J[, col_idx]
        else if (timing == 0L)  f_zero[, j]  <- J[, col_idx]
        else if (timing == 1L)  f_plus[, j]  <- J[, col_idx]
      }
    }
  }

  # Exogenous columns
  for (k in seq_along(exo)) {
    col_idx <- dyn$n_dyn_cols + k
    if (col_idx <= ncol(J)) {
      f_exo[, k] <- J[, col_idx]
    }
  }

  # Classify variables STRUCTURALLY from the lead_lag_incidence (see
  # .structural_lag_lead): a variable is a state iff its lag is structurally
  # present, even when the evaluated Jacobian coefficient is zero at this SS
  # (e.g. a Ramsey multiplier that is zero at the deterministic steady state).
  # A purely numerical threshold on f_minus/f_plus would silently drop such
  # variables (dropping the instrument `i` from NK_optimal_policy's state set).
  sl <- .structural_lag_lead(lli, row_labels)
  has_lag  <- sl$has_lag
  has_lead <- sl$has_lead
  is_static <- !has_lag & !has_lead
  is_pred   <- has_lag & !has_lead
  is_fwd    <- !has_lag & has_lead
  is_mixed  <- has_lag & has_lead

  # State variables: those that appear at t-1 (predetermined + mixed)
  state_vars <- endo[has_lag]
  n_state <- length(state_vars)

  # Forward vars: those with leads (forward + mixed)
  fwd_vars <- endo[has_lead]
  n_fwd <- length(fwd_vars)

  list(
    f_minus = f_minus, f_zero = f_zero, f_plus = f_plus, f_exo = f_exo,
    jacobian = J,
    n_eq = n_eq, n_endo = n_endo, n_exo = n_exo,
    n_state = n_state, n_fwd = n_fwd,
    state_vars = state_vars, fwd_vars = fwd_vars,
    is_static = is_static, is_pred = is_pred,
    is_fwd = is_fwd, is_mixed = is_mixed,
    endo_names = endo, exo_names = exo,
    ## Compiled-equation -> declaration-variable mapping used to reorder J
    ## (line ~193).  Returned so the order-2 solver permutes the Hessian with
    ## the SAME mapping instead of recomputing it from the already-reordered
    ## f_zero -- the latter yields the identity and leaves the Hessian
    ## misaligned for compound-LHS models (H5).
    eq_to_decl = eq_to_decl
  )
}


# =====================================================================
# Cached system structure (moved from dynhr_perf.R)
# =====================================================================

#' Pre-compute structural indices (call ONCE before MCMC)
#'
#' Caches everything that depends on model structure but not on params:
#' row_labels, column mappings from lli, variable classification, dy key names.
#'
#' @param compiled dynhr_compiled
#' @return A cache list to pass to extract_system_matrices_fast
#' @noRd
## Package-private cache for system structure, keyed on the jacobian_fn
## closure identity (function objects are stable across `solve_perturbation`
## calls on the same compiled model). Avoids re-running the regex/vapply
## parsing of lli row labels on every posterior evaluation.
.get_sys_cache <- function(compiled) {
  ## jacobian_fn is built once per compile_model() call by eval(parse(...)),
  ## so its enclosing environment is unique and persistent for that compiled
  ## model. We stash the cache inside that env, which gives us
  ## per-compiled-model memoisation with no extra plumbing.
  jfn <- compiled$dynamic$jacobian_fn
  if (is.function(jfn)) {
    env <- environment(jfn)
    if (!is.null(env)) {
      hit <- env$.dynhr_sys_cache
      if (!is.null(hit)) return(hit)
      c1 <- cache_system_structure(compiled)
      env$.dynhr_sys_cache <- c1
      return(c1)
    }
  }
  cache_system_structure(compiled)
}

cache_system_structure <- function(compiled) {
  dyn    <- compiled$dynamic
  model  <- compiled$model
  endo   <- model$var_names
  exo    <- model$varexo_names
  n_eq   <- dyn$n_eq
  n_endo <- length(endo)
  n_exo  <- length(exo)
  lli    <- compiled$lead_lag_incidence %||% compiled$model$lead_lag_incidence

  row_labels <- vapply(rownames(lli), function(rn) {
    rn <- trimws(rn)
    if (rn == "t") return(0L)
    m <- regmatches(rn, regexec("^t([+-]?\\d+)$", rn))[[1]]
    if (length(m) == 2) return(as.integer(m[2]))
    0L
  }, integer(1), USE.NAMES = FALSE)

  minus_map <- list(); zero_map <- list(); plus_map <- list()
  for (j in seq_along(endo)) {
    for (i in seq_len(nrow(lli))) {
      col_idx <- lli[i, j]
      if (col_idx > 0L) {
        timing <- row_labels[i]
        if      (timing == -1L) minus_map[[length(minus_map)+1]] <- c(j, col_idx)
        else if (timing ==  0L) zero_map[[length(zero_map)+1]]   <- c(j, col_idx)
        else if (timing ==  1L) plus_map[[length(plus_map)+1]]   <- c(j, col_idx)
      }
    }
  }
  minus_idx <- if (length(minus_map)) do.call(rbind, minus_map) else matrix(integer(0), 0, 2)
  zero_idx  <- if (length(zero_map))  do.call(rbind, zero_map)  else matrix(integer(0), 0, 2)
  plus_idx  <- if (length(plus_map))  do.call(rbind, plus_map)  else matrix(integer(0), 0, 2)

  exo_jcols <- dyn$n_dyn_cols + seq_along(exo)

  has_lag  <- apply(lli, 2, function(x) any(x[row_labels == -1L] > 0))
  has_lead <- apply(lli, 2, function(x) any(x[row_labels ==  1L] > 0))
  state_vars <- endo[has_lag]
  fwd_vars   <- endo[has_lead]

  dy_keys       <- character(nrow(dyn$dyn_col_map))
  dy_is_endo    <- logical(nrow(dyn$dyn_col_map))
  dy_endo_names <- character(nrow(dyn$dyn_col_map))
  for (k in seq_len(nrow(dyn$dyn_col_map))) {
    nm  <- dyn$dyn_col_map$name[k]
    ll  <- dyn$dyn_col_map$lead_lag[k]
    sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll) else paste0("__m", abs(ll))
    dy_keys[k]       <- paste0(nm, sfx)
    dy_is_endo[k]    <- nm %in% endo
    dy_endo_names[k] <- nm
  }

  list(
    dyn = dyn, endo = endo, exo = exo,
    n_eq = n_eq, n_endo = n_endo, n_exo = n_exo,
    lli = lli, row_labels = row_labels,
    minus_idx = minus_idx, zero_idx = zero_idx, plus_idx = plus_idx,
    exo_jcols = exo_jcols,
    has_lag = has_lag, has_lead = has_lead,
    state_vars = state_vars, fwd_vars = fwd_vars,
    n_state = length(state_vars), n_fwd = length(fwd_vars),
    is_static = !has_lag & !has_lead,
    is_pred   = has_lag & !has_lead,
    is_fwd    = !has_lag & has_lead,
    is_mixed  = has_lag & has_lead,
    dy_keys = dy_keys, dy_is_endo = dy_is_endo, dy_endo_names = dy_endo_names
  )
}


#' Fast system matrix extraction using cached structure
#'
#' Only fills dy values, evaluates the Jacobian, then slices using
#' pre-built index maps. No regex, no vapply, no apply per call.
#'
#' @param cache  Output of cache_system_structure
#' @param ss     Named numeric steady state vector
#' @param params Named numeric parameter vector
#' @param center Optional named numeric vector: the linearization point.
#'   When non-NULL, `dy` is built from `center` instead of `ss`. `ss` is
#'   still passed as the STEADY_STATE() oracle. Default NULL => ss.
#' @return Same list structure as extract_system_matrices
#' @noRd
extract_system_matrices_fast <- function(cache, ss, params, center = NULL) {
  n_eq   <- cache$n_eq
  n_endo <- cache$n_endo
  n_exo  <- cache$n_exo
  endo   <- cache$endo
  exo    <- cache$exo

  dy <- numeric(length(cache$dy_keys))
  names(dy) <- cache$dy_keys
  ## Vectorized fill: endo columns whose name is in the linearization point take
  ## that value; the rest stay 0. When center is non-NULL use it as the
  ## linearization point; otherwise fall back to ss (classic behavior).
  lp <- if (!is.null(center)) center else ss
  fill <- cache$dy_is_endo & (cache$dy_endo_names %in% names(lp))
  if (any(fill)) dy[fill] <- lp[cache$dy_endo_names[fill]]

  ## Evaluate the dynamic Jacobian. Prefer the C++ stack-machine tape (positional
  ## indexing, one Rcpp call) over the interpreted jacobian_fn closure (named-
  ## vector string lookups); fall back to the closure when the tape is absent
  ## (unsupported construct) or the DLL/use_rcpp option is off. Bit-parity to
  ## machine precision is asserted by test-jac-tape-parity.R.
  ## NOTE: when center is non-NULL the tape path uses lp for dy but still passes
  ## ss as the STEADY_STATE() oracle (ss_pos). This is correct: STEADY_STATE()
  ## is a model-defined reference, not the linearization point.
  dyn <- cache$dyn
  jt  <- dyn$jac_tape
  if (!is.null(jt) && .HAS_RCPP_TAPE() && is.null(center)) {
    ## C++ tape path: only safe for center=NULL (ss linearization) because the
    ## tape packs ss positionally into dy; for an off-SS center we fall back to
    ## the R closure which uses the named dy we just built above.
    par_pos <- params[dyn$jac_tape_param_names]
    ss_pos  <- ss[dyn$endo_names]; ss_pos[is.na(ss_pos)] <- 0
    J <- eval_jac_tape_cpp(unname(dy), unname(par_pos), unname(ss_pos),
                           jt$op, jt$ia, jt$da, jt$expr_len, jt$out_row, jt$out_col,
                           dyn$n_eq, dyn$total_cols)
  } else {
    J <- dyn$jacobian_fn(dy, params, ss)
  }

  ## Mirror the slow path's non-finite Jacobian repair (NaN/Inf from
  ## expressions like x/x at x=0). Without this, downstream classification
  ## sees NA columns and breaks (e.g. nk_hs2016 'pi' column at default ss).
  if (anyNA(J) || any(!is.finite(J)))
    J <- .repair_nonfinite_jacobian(J, cache$dyn, dy, params, ss)

  f_minus <- matrix(0, nrow = n_eq, ncol = n_endo)
  f_zero  <- matrix(0, nrow = n_eq, ncol = n_endo)
  f_plus  <- matrix(0, nrow = n_eq, ncol = n_endo)
  f_exo   <- matrix(0, nrow = n_eq, ncol = n_exo)
  colnames(f_minus) <- endo; colnames(f_zero) <- endo
  colnames(f_plus)  <- endo; colnames(f_exo)  <- exo

  ## Vectorized column scatter (replaces row-by-row R loops): each *_idx is a
  ## 2-col map [target endo col, source J col]. Single assignments match the
  ## loop semantics (last write wins on any duplicate target).
  mi <- cache$minus_idx
  if (nrow(mi) > 0) f_minus[, mi[, 1]] <- J[, mi[, 2]]
  zi <- cache$zero_idx
  if (nrow(zi) > 0) f_zero[, zi[, 1]]  <- J[, zi[, 2]]
  pi <- cache$plus_idx
  if (nrow(pi) > 0) f_plus[, pi[, 1]]  <- J[, pi[, 2]]
  ej <- cache$exo_jcols
  evalid <- which(ej <= ncol(J))
  if (length(evalid)) f_exo[, evalid] <- J[, ej[evalid]]

  ## STRUCTURAL classification of has_lag / has_lead from the lead_lag_incidence
  ## (see .structural_lag_lead).  Some lli entries are structurally present but
  ## numerically zero at the current steady state / params (e.g. ireland_2004
  ## 'zt' lag with rho_z = 0, or a Ramsey multiplier MULT_j that is zero at the
  ## deterministic SS so a MULT_j(-1)*(...x(-1)...) term contributes nothing to
  ## f_minus).  A numerical threshold on f_minus/f_plus would drop such variables
  ## from the state set (e.g. the instrument `i` on NK_optimal_policy), so we
  ## classify structurally -- matching Dynare and the slow extract path.  This
  ## also keeps the fast and slow paths' state_vars identical.
  sl <- .structural_lag_lead(cache$lli, cache$row_labels)
  has_lag  <- sl$has_lag
  has_lead <- sl$has_lead
  is_static <- !has_lag & !has_lead
  is_pred   <- has_lag & !has_lead
  is_fwd    <- !has_lag & has_lead
  is_mixed  <- has_lag & has_lead
  state_vars <- endo[has_lag]
  fwd_vars   <- endo[has_lead]

  list(
    f_minus = f_minus, f_zero = f_zero, f_plus = f_plus, f_exo = f_exo,
    jacobian = J,
    n_eq = n_eq, n_endo = n_endo, n_exo = n_exo,
    n_state = length(state_vars), n_fwd = length(fwd_vars),
    state_vars = state_vars, fwd_vars = fwd_vars,
    is_static = is_static, is_pred = is_pred,
    is_fwd = is_fwd, is_mixed = is_mixed,
    endo_names = endo, exo_names = exo
  )
}

## Shared between slow + fast extractors: repair non-finite Jacobian entries
## via central finite differences against the residual function, falling back
## to zero where the residual itself is non-finite.
.repair_nonfinite_jacobian <- function(J, dyn, dy, params, ss) {
  bad_jac <- which(!is.finite(J), arr.ind = TRUE)
  if (length(bad_jac) == 0L) return(J)
  bad_cols <- sort(unique(bad_jac[, 2]))
  f0 <- if (!is.null(dyn$residuals_fn))
    dyn$residuals_fn(dy, params, ss) else NULL
  if (!is.null(f0) && all(is.finite(f0))) {
    for (bc in bad_cols) {
      dy_p <- dy; dy_m <- dy
      h <- max(1e-6, abs(dy[bc]) * 1e-6)
      dy_p[bc] <- dy[bc] + h
      dy_m[bc] <- dy[bc] - h
      fp <- dyn$residuals_fn(dy_p, params, ss)
      fm <- dyn$residuals_fn(dy_m, params, ss)
      if (!is.null(fp) && !is.null(fm) &&
          all(is.finite(fp)) && all(is.finite(fm)))
        J[, bc] <- (fp - fm) / (2 * h)
      else
        J[!is.finite(J[, bc]), bc] <- 0
    }
  } else {
    J[!is.finite(J)] <- 0
  }
  J
}
