## R/compile-dynamic.R
## --------------------------------------------------------------------------
## Build the DYNAMIC model: preserve lead/lag timing, differentiate residuals
## with respect to (y_lag, y, y_lead, x), generate residuals_fn / jacobian_fn
## closures and column-mapping metadata used by the perturbation solver.
##
## Also computes the SPARSE SYMBOLIC SECOND-ORDER HESSIAN used by the
## second-order perturbation solver (solve-perturbation-order2.R).
##
## Algorithm (FaÃ  di Bruno / Binning 2013 matrix chain rule):
##   For each equation i, and each pair of columns (c1, c2) where both
##   âˆ‚F_i/âˆ‚w_c1 and âˆ‚F_i/âˆ‚w_c2 are symbolically non-zero, compute
##   âˆ‚Â²F_i/(âˆ‚w_c1 âˆ‚w_c2) = ast_differentiate2(residual_asts[[i]], ...).
##   Store non-zero results as (eq, col1, col2, expr) triplets and generate
##   hessian2_fn: (dy, params, ss) -> numeric vector of sparse values.
##
## References:
##   Binning (2013), Norges Bank WP 2013/18 â€” matrix chain rules for 2nd order
##   Levintal (2017), JEDC â€” 5th-order perturbation using augmented state
##
## Phase-1 split from jacobian-monolith.R (no logic changes to first-order).
## --------------------------------------------------------------------------

## Wrap a generated closure so it executes with the bytecode JIT disabled.
## Large generated derivative functions (thousands of huge expressions in one
## body) are catastrophic to byte-compile (codegen-bound, minutes) but run
## fine interpreted. R's JIT (enableJIT level >=1) would otherwise compile such
## a closure on first use. We keep it interpreted by disabling the JIT for the
## duration of each call (cheap; restored on exit). Only applied to bodies
## large enough that JIT compilation is a net loss; small functions are left
## untouched so they still benefit from JIT in hot MCMC loops.
.interp_only_if_large <- function(fn, fn_text, threshold = 50000L) {
  if (nchar(fn_text) <= threshold) return(fn)
  force(fn)
  function(...) {
    old <- compiler::enableJIT(0L)
    on.exit(compiler::enableJIT(old), add = TRUE)
    fn(...)
  }
}

## Build (and JIT-guard) a generated derivative function from a list of ASTs
## using common-subexpression elimination (ast_cse_emit). Emitting each interior
## subexpression once into a temporary keeps every generated statement shallow --
## avoiding R's "contextstack overflow" when parse()ing the deeply-nested
## expressions of large higher-order models -- and computes shared subexpressions
## once (Dynare-style), shrinking and speeding the function. `lhs` is one
## left-hand side per AST (entries whose AST is NULL are skipped); `alloc` is the
## output-object allocation line; `ret` the returned object name. The result is
## wrapped with .interp_only_if_large so the (still large) function is not
## byte-compiled by the JIT.
.cse_fn <- function(asts, lhs, alloc, ret, endo, exo, params, use_timing = TRUE,
                    local_vars = list()) {
  em   <- ast_cse_emit(asts, endo, exo, params, use_timing = use_timing)
  keep <- !is.na(em$refs)
  ## Emit #-local let-bindings in topological order so LOCAL_x refs in the raw
  ## ASTs resolve at runtime without pre-substitution (I20 perf fix).
  local_lines <- if (length(local_vars) > 0L) {
      ln  <- names(local_vars)
      vapply(seq_along(local_vars), function(k) {
          ## use_timing MUST match the enclosing closure: the dynamic fn signature
          ## is function(dy, ...) with NO `y`, so a #-local that references an endo
          ## variable must emit dy["name__0"]-style refs (use_timing=TRUE), not the
          ## static y["name"] form -- else "object 'y' not found" at runtime
          ## (nz_base_nonlinear: #-locals over timed endo vars).
          rhs <- ast_to_fn_body(local_vars[[k]], endo_names = endo,
                                exo_names = exo, param_names = params,
                                use_timing = use_timing)
          paste0("    LOCAL_", ln[k], " <- ", rhs)
      }, character(1))
  } else character(0)
  txt  <- paste(c(
    "function(dy, params, ss = numeric(0)) {",
    paste0("    ", alloc),
    local_lines,
    if (length(em$decls)) paste0("    ", em$decls) else character(0),
    paste0("    ", lhs[keep], " <- ", em$refs[keep]),
    paste0("    ", ret),
    "}"
  ), collapse = "\n")
  gen <- eval(parse(text = txt))
  # The generated body references only its args (dy/params/ss) and base/stats
  # functions, so detach it from this (large) build frame: capturing `em`,
  # `asts`, etc. would bloat every compiled model in memory and on disk (slow
  # to serialize for the compile cache). Use a FRESH env (not the namespace
  # itself, which is locked) parented on the dynhr namespace so base/stats
  # generics resolve AND callers may still stash per-model scratch state in the
  # function's environment (e.g. .get_sys_cache uses jacobian_fn's env).
  environment(gen) <- new.env(parent = asNamespace("dynhr"))
  .interp_only_if_large(gen, txt)
}

#' Build dynamic model residual and Jacobian functions
#'
#' The dynamic model preserves all leads and lags. The Jacobian columns
#' follow Dynare's convention: ordered by the lead_lag_incidence matrix.
#'
#' @param model            A dynhr_mod object from parse_mod().
#' @param max_order        Maximum derivative order to compute symbolically.
#'                         1 = Jacobian, 2 = +Hessian, 3 = +3rd deriv, etc.
#'                         Default 5. Setting to 1 avoids node stack overflow in
#'                         models with many complex equations.
#' @param want_param_deriv Logical. When TRUE (the default when the model is
#'                         parameter-differentiable) the symbolic
#'                         parameter-Jacobian and second-order Hessian are
#'                         compiled so the analytic gradient is available.
#'                         When FALSE both steps are skipped; the gradient layer
#'                         falls back to finite differences.  max_order >= 2
#'                         always compiles hess2 regardless.
#' @return A list with residuals_fn, jacobian_fn, column map, and expressions.
#' @noRd
build_dynamic_model <- function(model, max_order = 1L,
                                want_param_deriv = TRUE,
                                want_param_deriv2 = FALSE) {
    endo   <- model$var_names
    exo    <- model$varexo_names
    params <- model$param_names
    n_eq   <- length(model$equations)
    n_endo <- length(endo)
    n_exo  <- length(exo)
    lli    <- model$lead_lag_incidence

    # Parse row labels: "t" -> 0, "t-1" -> -1, "t+1" -> +1
    row_labels <- vapply(rownames(lli), function(rn) {
        rn <- trimws(rn)
        if (rn == "t") return(0L)
        m <- regmatches(rn, regexec("^t([+-]?\\d+)$", rn))[[1]]
        if (length(m) == 2) return(as.integer(m[2]))
        val <- as.integer(gsub("[^0-9-]", "", rn))
        if (is.na(val)) 0L else val
    }, integer(1), USE.NAMES = FALSE)

    # Build the column map: list of (var_name, lead_lag, col_index)
    n_dyn_cols <- max(lli)
    dyn_col_map <- data.frame(
        name     = character(0),
        lead_lag = integer(0),
        col      = integer(0),
        stringsAsFactors = FALSE
    )
    for (j in seq_along(endo)) {
        for (i in seq_len(nrow(lli))) {
            if (lli[i, j] > 0) {
                dyn_col_map <- rbind(dyn_col_map, data.frame(
                    name     = endo[j],
                    lead_lag = row_labels[i],
                    col      = lli[i, j],
                    stringsAsFactors = FALSE
                ))
            }
        }
    }
    dyn_col_map <- dyn_col_map[order(dyn_col_map$col), ]

    # Add exogenous variables at the end (at time t)
    for (k in seq_along(exo)) {
        dyn_col_map <- rbind(dyn_col_map, data.frame(
            name     = exo[k],
            lead_lag = 0L,
            col      = n_dyn_cols + k,
            stringsAsFactors = FALSE
        ))
    }
    total_cols <- n_dyn_cols + n_exo

    # 1. Build residuals: lhs - rhs. NO local substitution (I20 perf fix).
    #
    # Performance: models with deep #-local chains produce enormous ASTs after
    # ast_substitute_locals (e.g. US_IN10 BB_SS -> 796 KB). Differentiating a
    # 3.7 MB residual w.r.t. 142 dynamic columns takes ~140s for one equation;
    # differentiating the raw pre-substitution AST takes <0.1s.
    #
    # Fix: keep residuals in raw form. Generated closures receive #-locals as
    # named let-bindings (via local_vars= in .cse_fn), so LOCAL_x references
    # resolve at runtime without ever expanding them symbolically.
    have_locals <- length(model$local_variables) > 0L
    local_vars  <- model$local_variables

    residual_asts <- vector("list", n_eq)
    for (i in seq_len(n_eq)) {
        eq  <- model$equations[[i]]
        res <- equation_to_residual(eq)
        res <- ast_simplify(res)
        residual_asts[[i]] <- res
    }

    # I20 correctness guard: the raw-AST fast path below treats LOCAL_x as a leaf
    # with derivative 0 w.r.t. variables. That is ONLY valid when no #-local
    # references an endo/exo variable (US_IN10's locals are param/SS chains). When
    # a local DOES depend on a variable (e.g. nz_base_nonlinear `#wedge = q - tot`),
    # the dynamic VARIABLE Jacobian must chain-rule through it — so fall back to
    # substituting locals into the residuals and disable the let-binding path.
    if (have_locals && !.locals_are_param_only(local_vars)) {
        residual_asts <- lapply(residual_asts, ast_substitute_locals, local_vars)
        residual_asts <- lapply(residual_asts, ast_simplify)
        local_vars <- list()
        have_locals <- FALSE
    }

    # 3. Jacobian: differentiate w.r.t. each (var, lead_lag) in the column map.
    #    Differentiate raw ASTs (fast path). LOCAL_x nodes in the result have
    #    derivative 0 w.r.t. endo/exo variables (correct: locals depend only on
    #    params). The non-zero result ASTs keep their LOCAL_x references; the
    #    generated closure resolves them from the let-bindings in .cse_fn.
    jac_triplets <- list()

    for (i in seq_len(n_eq)) {
        for (k in seq_len(nrow(dyn_col_map))) {
            vname <- dyn_col_map$name[k]
            vll   <- dyn_col_map$lead_lag[k]
            vcol  <- dyn_col_map$col[k]

            d <- ast_differentiate(residual_asts[[i]], vname, vll)
            if (!ast_is_zero(d)) {
                jac_triplets <- c(jac_triplets, list(list(
                    row  = i,
                    col  = vcol,
                    ast  = d
                )))
            }
        }
    }

    # -----------------------------------------------------------------------
    # 3a. Symbolic parameter-Jacobian  ∂²F_i/(∂w_c ∂θ_k)  (Tier 11 #3)
    #
    # Explicit partial of each dynamic-Jacobian entry w.r.t. each structural
    # parameter, obtained by differentiating the already-computed first-order
    # Jacobian ASTs a second time in the parameter direction. Combined at
    # solve time with the steady-state chain rule (model Hessian · dys), this
    # makes the first-order solution-derivative gradient finite-difference-free.
    #
    # Only built when every residual is smoothly parameter-differentiable
    # (no STEADY_STATE() reference, no max/min): otherwise the explicit
    # partial is incomplete and the solution-derivative layer must use FD.
    # -----------------------------------------------------------------------
    n_params <- length(params)
    local_names <- names(local_vars)
    # param_deriv_ok: model is analytically parameter-differentiable (no
    # STEADY_STATE() refs etc.) AND the caller has not opted out via
    # want_param_deriv = FALSE.  Check on raw residuals (local_variable nodes
    # are safe; STEADY_STATE/max/min appear in the pre-sub form).
    param_deriv_ok <- isTRUE(want_param_deriv) &&
        all(vapply(residual_asts, .ast_param_deriv_safe, logical(1)))

    ## Pre-compute d(LOCAL_k)/d(θ) via forward sweep once (fast).
    ## Used for chain-rule param differentiation of Jacobian ASTs that contain
    ## LOCAL_x references (jac_triplets$ast and derived higher-order ASTs).
    dL_dtheta <- if (have_locals && param_deriv_ok && n_params > 0L)
        .local_param_sweep(local_vars, params) else list()

    param_jac_triplets <- list()
    if (param_deriv_ok && n_params > 0L) {
        for (t in jac_triplets) {
            for (k in seq_len(n_params)) {
                dp <- if (have_locals)
                    .chain_rule_param_deriv(t$ast, params[k], k,
                                            local_names, dL_dtheta)
                else ast_differentiate_param(t$ast, params[k])
                if (!is.null(dp) && !ast_is_zero(dp)) {
                    param_jac_triplets <- c(param_jac_triplets, list(list(
                        row   = t$row,
                        col   = t$col,
                        param = k,
                        ast   = dp
                    )))
                }
            }
        }
    }

    # -----------------------------------------------------------------------
    # 3b. Symbolic second-order Hessian (for second-order perturbation)
    #
    # Strategy (Binning 2013 Â§3): if âˆ‚F_i/âˆ‚w_j = 0 identically, then
    # âˆ‚Â²F_i/(âˆ‚w_j âˆ‚w_k) = 0 for all k.  So we only check pairs (c1, c2)
    # where BOTH columns have non-zero first-order Jacobian in equation i.
    # This exploits sparsity and reuses the already-computed first-order ASTs.
    # -----------------------------------------------------------------------

    # Build per-equation lookup: list[[i]][[col_key]] = first-order AST
    jac_ast_by_eq <- vector("list", n_eq)
    for (i in seq_len(n_eq)) jac_ast_by_eq[[i]] <- list()
    for (t in jac_triplets) {
        jac_ast_by_eq[[t$row]][[as.character(t$col)]] <- t$ast
    }

    # Build reverse map from column index to (name, lead_lag) for O(1) lookup
    col_info_name <- character(total_cols)
    col_info_ll   <- integer(total_cols)
    for (k in seq_len(nrow(dyn_col_map))) {
        c <- dyn_col_map$col[k]
        col_info_name[c] <- dyn_col_map$name[k]
        col_info_ll[c]   <- dyn_col_map$lead_lag[k]
    }

    hess2_triplets <- list()
    hess3_triplets <- list()  # third-order: (eq, col1, col2, col3, expr) with c1<=c2<=c3

    # Build the symbolic second-order model Hessian (hess2) when EITHER the
    # perturbation order needs it (max_order >= 2) OR analytic parameter
    # derivatives are requested. The first-order solution-derivative gradient
    # (Tier 11 #3) contracts this model Hessian with dys for the steady-state
    # chain rule, so it is needed even at max_order = 1 -- without it the
    # analytic gradient would silently drop the chain term.
    # When want_param_deriv = FALSE the caller has explicitly opted out of the
    # analytic gradient, so hess2 is skipped unless max_order >= 2 needs it.
    # Third-order (hess3) and higher stay gated on the perturbation order alone.
    SKIP_HIGHER_HESSIAN <- (max_order < 4L)
    build_hess2   <- (max_order >= 2L) || param_deriv_ok
    # hess3 is needed by the order-2 analytic primitives: the model-Hessian
    # derivative dH_mat (cumulant path, max_order >= 2) AND the second total
    # Jacobian derivative d2f_ij (exact-Hessian path, requested via
    # want_param_deriv2 even at max_order = 1). Forcing it under
    # want_param_deriv2 also enables param_hess2 (its gate includes build_hess3up).
    build_hess3up <- (max_order >= 2L) || isTRUE(want_param_deriv2)
    if (build_hess2) {
    for (i in seq_len(n_eq)) {
        nz_col_keys <- names(jac_ast_by_eq[[i]])
        if (length(nz_col_keys) < 1L) next
        nz_cols <- as.integer(nz_col_keys)  # sorted by names() (as strings)

        # Upper triangle (col1 <= col2); symmetry gives the lower triangle free
        for (k1_idx in seq_along(nz_cols)) {
            c1   <- nz_cols[k1_idx]
            d1   <- jac_ast_by_eq[[i]][[nz_col_keys[k1_idx]]]
            nm2  <- col_info_name[c1]   # for self-partial: same var
            ll2  <- col_info_ll[c1]

            for (k2_idx in seq(k1_idx, length(nz_cols))) {
                c2   <- nz_cols[k2_idx]
                nm2  <- col_info_name[c2]
                ll2  <- col_info_ll[c2]

                d2 <- ast_differentiate(d1, nm2, ll2)

                if (!ast_is_zero(d2)) {
                    hess2_triplets <- c(hess2_triplets, list(list(
                        eq   = i,
                        col1 = c1,
                        col2 = c2,
                        ast  = d2   # CSE builds hessian2_fn from $ast; also used
                                    # for Tier 11 #3 param_hess2 codegen
                    )))

                    # Third-order: differentiate d2 again, with c3 >= c2 to
                    # enforce canonical ordering c1 <= c2 <= c3.  If d2 is
                    # already non-zero, only further differentiation can produce
                    # non-zero d3 -- iterate the same non-zero-Jacobian column
                    # set as a sufficient (but not necessary) sparsity filter.
                    if (build_hess3up)
                    for (k3_idx in seq(k2_idx, length(nz_cols))) {
                        c3  <- nz_cols[k3_idx]
                        nm3 <- col_info_name[c3]
                        ll3 <- col_info_ll[c3]
                        d3  <- ast_differentiate(d2, nm3, ll3)
                        if (!ast_is_zero(d3)) {
                            hess3_triplets <- c(hess3_triplets, list(list(
                                eq   = i,
                                col1 = c1,
                                col2 = c2,
                                col3 = c3,
                                ast  = d3
                            )))
                        }
                    }
                }
            }
        }
    }

    } # end if (build_hess2)

    n_hess  <- length(hess2_triplets)
    n_hess3 <- length(hess3_triplets)

    # -----------------------------------------------------------------------
    # 3a-2. Symbolic third-order parameter tensor  ∂³F/(∂w_c1 ∂w_c2 ∂θ_k)
    #       (Tier 11 #3, order-2 layer = param_hess2).
    #
    # Differentiates each second-order model-Hessian AST a third time in the
    # parameter direction.  Together with the steady-state chain term (hess3
    # contracted with dys) this gives the analytic total θ-derivative of the
    # model Hessian (dH_mat) used by the analytic order-2 solution-derivative
    # gradient, removing the per-parameter central FD of hessian2_fn.
    #
    # Built only when the explicit param-Jacobian path is viable (param_deriv_ok)
    # AND both hess2 and hess3 are actually compiled -- otherwise the consumer
    # (solution_derivatives_order2) would silently drop the chain term, so it
    # gates on param_hess2_built && hessian3_built and falls back to FD.
    # -----------------------------------------------------------------------
    build_param_hess2 <- param_deriv_ok && n_params > 0L &&
        build_hess2 && build_hess3up
    param_hess2_triplets <- list()
    if (build_param_hess2) {
        for (t in hess2_triplets) {
            for (k in seq_len(n_params)) {
                dp <- if (have_locals)
                    .chain_rule_param_deriv(t$ast, params[k], k,
                                            local_names, dL_dtheta)
                else ast_differentiate_param(t$ast, params[k])
                if (!is.null(dp) && !ast_is_zero(dp)) {
                    param_hess2_triplets <- c(param_hess2_triplets, list(list(
                        eq    = t$eq,
                        col1  = t$col1,
                        col2  = t$col2,
                        param = k,
                        ast   = dp
                    )))
                }
            }
        }
    }
    n_param_hess2 <- length(param_hess2_triplets)

    # -----------------------------------------------------------------------
    # 3a-3. Symbolic var-param-param tensor  ∂³F/(∂w_c ∂θ_a ∂θ_b)
    #       (Tier 11 #3, order-2 layer = param2_jac, used by 3b only).
    #
    # Differentiates each FIRST-order dynamic-Jacobian AST twice in the
    # parameter direction (canonical a <= b; the tensor is symmetric in
    # (a,b)).  Supplies the leading explicit term of the analytic second total
    # θ-derivative of the dynamic-Jacobian blocks (d2f_*), removing the
    # diagonal/4-corner FD stencils in solution_derivatives_2.
    #
    # Built only on explicit request (want_param_deriv2, i.e. param_deriv =
    # "second") since it is 3b-only and O(n_jac * n_params^2): the consumer
    # (solution_derivatives_2) gates on param_deriv2_ok && the static
    # second-order tensors and falls back to FD otherwise.
    # -----------------------------------------------------------------------
    build_param2_jac <- param_deriv_ok && n_params > 0L && isTRUE(want_param_deriv2)
    param2_jac_triplets <- list()
    if (build_param2_jac) {
        for (t in jac_triplets) {
            for (a in seq_len(n_params)) {
                d1 <- if (have_locals)
                    .chain_rule_param_deriv(t$ast, params[a], a,
                                            local_names, dL_dtheta)
                else ast_differentiate_param(t$ast, params[a])
                if (is.null(d1) || ast_is_zero(d1)) next
                for (b in seq(a, n_params)) {
                    d2 <- if (have_locals)
                        .chain_rule_param_deriv(d1, params[b], b,
                                                local_names, dL_dtheta)
                    else ast_differentiate_param(d1, params[b])
                    if (!is.null(d2) && !ast_is_zero(d2)) {
                        param2_jac_triplets <- c(param2_jac_triplets, list(list(
                            row  = t$row,
                            col  = t$col,
                            pa   = a,
                            pb   = b,
                            ast  = d2
                        )))
                    }
                }
            }
        }
    }
    n_param2_jac <- length(param2_jac_triplets)

    # -----------------------------------------------------------------------
    # 3b. Fourth-order Hessian triplets (c1 â‰¤ c2 â‰¤ c3 â‰¤ c4 canonical order)
    #
    # Only compute âˆ‚â´F_i/(âˆ‚w_c1 âˆ‚w_c2 âˆ‚w_c3 âˆ‚w_c4) where the 3rd-order
    # entry (c1,c2,c3) is non-zero AND âˆ‚F_i/âˆ‚w_c4 â‰  0.  Differentiate the
    # 3rd-order AST d3 â‰ˆ âˆ‚Â³F_i/(âˆ‚w_c1 âˆ‚w_c2 âˆ‚w_c3) with respect to w_c4.
    #
    # Heuristic: if n_hess3 > 500, skip Hessian4/5 because the generated
    # expression becomes too large for R's parser (context stack overflow)
    # on models with many complex equations (e.g. Caldara recursive pref.).
    # Order-4/5 solvers fall back to numerical FD when compiled Hessians are
    # unavailable.
    # -----------------------------------------------------------------------
    SKIP_HIGHER_HESSIAN <- SKIP_HIGHER_HESSIAN || (n_hess3 > 500L)
    hess4_triplets <- list()
    hess5_triplets <- list()
    if (!SKIP_HIGHER_HESSIAN) {
    for (t3 in hess3_triplets) {
      i   <- t3$eq
      c1  <- t3$col1
      c2  <- t3$col2
      c3  <- t3$col3
      d3  <- jac_ast_by_eq[[i]][[as.character(c1)]]
      # Need higher derivatives of d3, not just the Jacobian.  Re-differentiate
      # the original residual through c1,c2,c3 to get the AST for d3.
      # Actually, the pre-computed d3 is not stored in hess3_triplets (only the
      # R expression).  We need to recompute d3 by differentiating the residual.
      # Re-derive d3 from the original residual:
      d1 <- jac_ast_by_eq[[i]][[as.character(c1)]]
      nm2 <- col_info_name[c2]; ll2 <- col_info_ll[c2]
      d2  <- ast_differentiate(d1, nm2, ll2)
      nm3 <- col_info_name[c3]; ll3 <- col_info_ll[c3]
      d3  <- ast_differentiate(d2, nm3, ll3)

      nz_col_keys_i <- names(jac_ast_by_eq[[i]])
      nz_cols_i <- as.integer(nz_col_keys_i)
      c3_idx <- which(nz_cols_i == c3)
      if (length(c3_idx) != 1L) next
      for (k4_idx in seq(c3_idx, length(nz_cols_i))) {
        c4  <- nz_cols_i[k4_idx]
        nm4 <- col_info_name[c4]
        ll4 <- col_info_ll[c4]
        d4  <- ast_differentiate(d3, nm4, ll4)
        if (!ast_is_zero(d4)) {
          hess4_triplets <- c(hess4_triplets, list(list(
            eq   = i,
            col1 = c1,
            col2 = c2,
            col3 = c3,
            col4 = c4,
            ast  = d4
          )))
        }
      }
    }

    } # end if (!SKIP_HIGHER_HESSIAN)

    # -----------------------------------------------------------------------
    # 3c. Fifth-order Hessian triplets (c1 â‰¤ c2 â‰¤ c3 â‰¤ c4 â‰¤ c5 canonical)
    #
    # Differentiate the 4th-order AST with respect to w_c5 where c5 â‰¥ c4 and
    # âˆ‚F_i/âˆ‚w_c5 â‰  0.
    # -----------------------------------------------------------------------
    if (!SKIP_HIGHER_HESSIAN) {
    for (t4 in hess4_triplets) {
      i   <- t4$eq
      c1  <- t4$col1; c2 <- t4$col2; c3 <- t4$col3; c4 <- t4$col4

      # Re-derive d4 from the original residual through (c1,c2,c3,c4)
      d1 <- jac_ast_by_eq[[i]][[as.character(c1)]]
      nm2 <- col_info_name[c2]; ll2 <- col_info_ll[c2]
      d2  <- ast_differentiate(d1, nm2, ll2)
      nm3 <- col_info_name[c3]; ll3 <- col_info_ll[c3]
      d3  <- ast_differentiate(d2, nm3, ll3)
      nm4 <- col_info_name[c4]; ll4 <- col_info_ll[c4]
      d4  <- ast_differentiate(d3, nm4, ll4)

      nz_col_keys_i <- names(jac_ast_by_eq[[i]])
      nz_cols_i <- as.integer(nz_col_keys_i)
      c4_idx <- which(nz_cols_i == c4)
      if (length(c4_idx) != 1L) next
      for (k5_idx in seq(c4_idx, length(nz_cols_i))) {
        c5  <- nz_cols_i[k5_idx]
        nm5 <- col_info_name[c5]
        ll5 <- col_info_ll[c5]
        d5  <- ast_differentiate(d4, nm5, ll5)
        if (!ast_is_zero(d5)) {
          hess5_triplets <- c(hess5_triplets, list(list(
            eq   = i,
            col1 = c1,
            col2 = c2,
            col3 = c3,
            col4 = c4,
            col5 = c5,
            ast  = d5
          )))
        }
      }
    } # end for (t4 in hess4_triplets)
    } # end if (!SKIP_HIGHER_HESSIAN)

    n_hess  <- length(hess2_triplets)
    n_hess3 <- length(hess3_triplets)
    n_hess4 <- length(hess4_triplets)
    n_hess5 <- length(hess5_triplets)

    # -----------------------------------------------------------------------
    # 4. Build dynamic residuals function
    # -----------------------------------------------------------------------
    residuals_fn <- .cse_fn(
        residual_asts,
        paste0("r[", seq_len(n_eq), "]"),
        paste0("r <- numeric(", n_eq, "L)"), "r",
        endo, exo, params,
        local_vars = local_vars)

    # -----------------------------------------------------------------------
    # 5. Build dynamic Jacobian function
    # -----------------------------------------------------------------------
    jacobian_fn <- .cse_fn(
        lapply(jac_triplets, `[[`, "ast"),
        vapply(jac_triplets, function(t) paste0("J[", t$row, ",", t$col, "]"),
               character(1)),
        paste0("J <- matrix(0, nrow = ", n_eq, "L, ncol = ", total_cols, "L)"), "J",
        endo, exo, params,
        local_vars = local_vars)

    # -----------------------------------------------------------------------
    # 5a. Build the parameter-Jacobian function (Tier 11 #3)
    #
    # Returns an [n_eq x total_cols x n_params] array P with
    #   P[i, c, k] = ∂²F_i/(∂w_c ∂θ_k)   (explicit, ss held fixed).
    # NULL when parameter differentiation is unsupported for this model.
    # -----------------------------------------------------------------------
    if (param_deriv_ok && length(param_jac_triplets) > 0L) {
        param_jacobian_fn <- .cse_fn(
            lapply(param_jac_triplets, `[[`, "ast"),
            vapply(param_jac_triplets, function(t)
                paste0("P[", t$row, ",", t$col, ",", t$param, "]"), character(1)),
            paste0("P <- array(0, dim = c(", n_eq, "L, ", total_cols,
                   "L, ", n_params, "L))"), "P",
            endo, exo, params,
            local_vars = local_vars)
    } else if (param_deriv_ok) {
        # Differentiable but no entry depends on any parameter (e.g. a fully
        # log-linear model): the parameter-Jacobian is identically zero.
        param_jacobian_fn <- function(dy, params, ss = numeric(0))
            array(0, dim = c(n_eq, total_cols, n_params))
    } else {
        param_jacobian_fn <- NULL
    }

    # -----------------------------------------------------------------------
    # 5a-2. Build the param-Hessian2 function (Tier 11 #3, order-2 layer)
    #
    # Returns a numeric vector of length n_param_hess2 with the values of the
    # non-zero ∂³F_i/(∂w_c1 ∂w_c2 ∂θ_k) entries at (dy, params, ss).  Indices
    # are stored in param_hess2_triplets (eq, col1, col2, param); parallel to
    # hessian2_fn / hess2_triplets.  NULL when the path is not built.
    # -----------------------------------------------------------------------
    if (build_param_hess2 && n_param_hess2 > 0L) {
        param_hessian2_fn <- .cse_fn(
            lapply(param_hess2_triplets, `[[`, "ast"),
            paste0("v[", seq_len(n_param_hess2), "L]"),
            paste0("v <- numeric(", n_param_hess2, "L)"), "v",
            endo, exo, params,
            local_vars = local_vars)
    } else if (build_param_hess2) {
        # Hessian built but no second-order entry depends on any parameter:
        # the explicit param-Hessian2 tensor is identically zero.
        param_hessian2_fn <- function(dy, params, ss = numeric(0)) numeric(0L)
    } else {
        param_hessian2_fn <- NULL
    }

    # -----------------------------------------------------------------------
    # 5a-3. Build the param2-Jacobian function (Tier 11 #3, order-2 layer)
    #
    # Returns a numeric vector of length n_param2_jac with the values of the
    # non-zero ∂³F_i/(∂w_c ∂θ_a ∂θ_b) entries at (dy, params, ss).  Indices are
    # stored in param2_jac_triplets (row, col, pa, pb) with canonical pa <= pb;
    # parallel to param_jacobian_fn.  NULL when the path is not built.
    # -----------------------------------------------------------------------
    if (build_param2_jac && n_param2_jac > 0L) {
        param2_jacobian_fn <- .cse_fn(
            lapply(param2_jac_triplets, `[[`, "ast"),
            paste0("v[", seq_len(n_param2_jac), "L]"),
            paste0("v <- numeric(", n_param2_jac, "L)"), "v",
            endo, exo, params,
            local_vars = local_vars)
    } else if (build_param2_jac) {
        # Differentiable but no Jacobian entry has a non-zero second parameter
        # derivative (e.g. entries linear in every parameter): identically zero.
        param2_jacobian_fn <- function(dy, params, ss = numeric(0)) numeric(0L)
    } else {
        param2_jacobian_fn <- NULL
    }

    # First-order Jacobian "tape" for the C++ stack-machine evaluator
    # (src/jac_tape.cpp). Positional integer indexing replaces the interpreted
    # closure's named-vector string lookups -- the dominant per-draw cost in
    # extract_system_matrices_fast. NULL if any expression uses a construct the
    # tape VM does not support (e.g. a 2-arg funcall); extract then falls back
    # to jacobian_fn. dy_keys_tape must match the order in which
    # extract_system_matrices_fast / cache_system_structure build the dy vector
    # (dyn_col_map row order).
    dy_keys_tape <- vapply(seq_len(nrow(dyn_col_map)), function(k) {
        ll  <- dyn_col_map$lead_lag[k]
        sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll)
               else paste0("__m", abs(ll))
        paste0(dyn_col_map$name[k], sfx)
    }, character(1))
    jac_tape <- compile_jacobian_tape(jac_triplets, dy_keys_tape, params, endo)

    # Hessian "tapes" for the C++ stack-machine evaluator (eval_triplet_tape_cpp).
    # Same positional namespaces as the jac tape (dy_keys_tape / params / endo);
    # NULL if any expression uses a construct the VM does not support, in which
    # case the consumers fall back to the interpreted hessian2_fn / hessian3_fn.
    # MUST be built here, BEFORE the per-triplet $ast trees are dropped below.
    hess2_tape <- if (n_hess  > 0L)
        compile_hessian_tape(hess2_triplets, dy_keys_tape, params, endo) else NULL
    hess3_tape <- if (n_hess3 > 0L)
        compile_hessian_tape(hess3_triplets, dy_keys_tape, params, endo) else NULL

    # -----------------------------------------------------------------------
    # 6. Build symbolic second-order Hessian function
    #
    # Returns a numeric vector of length n_hess with the values of the
    # non-zero second-order derivatives at the given point (dy, params, ss).
    # The indices are stored in hess2_triplets (eq, col1, col2).
    # -----------------------------------------------------------------------
    if (n_hess > 0L) {
        hessian2_fn <- .cse_fn(
            lapply(hess2_triplets, `[[`, "ast"),
            paste0("v[", seq_len(n_hess), "L]"),
            paste0("v <- numeric(", n_hess, "L)"), "v",
            endo, exo, params,
            local_vars = local_vars)
    } else {
        # All-zero Hessian (log-linearised model): return zero-length vector
        hessian2_fn <- function(dy, params, ss = numeric(0)) numeric(0L)
    }

    # -----------------------------------------------------------------------
    # 6b. Build symbolic third-order Hessian function
    #
    # Returns a numeric vector of length n_hess3 with the values of the
    # non-zero third-order derivatives at (dy, params, ss).  Indices are
    # stored in hess3_triplets (eq, col1<=col2<=col3) with full Schwarz
    # symmetry; consumers expand the orbit (up to 6 permutations).
    # -----------------------------------------------------------------------
    if (n_hess3 > 0L) {
        hessian3_fn <- .cse_fn(
            lapply(hess3_triplets, `[[`, "ast"),
            paste0("v[", seq_len(n_hess3), "L]"),
            paste0("v <- numeric(", n_hess3, "L)"), "v",
            endo, exo, params,
            local_vars = local_vars)
    } else {
        hessian3_fn <- function(dy, params, ss = numeric(0)) numeric(0L)
    }

    # -----------------------------------------------------------------------
    # 6c. Build symbolic fourth-order Hessian function
    #
    # Returns a numeric vector of length n_hess4.  Canonical order c1â‰¤c2â‰¤c3â‰¤c4.
    # -----------------------------------------------------------------------
    if (n_hess4 > 0L) {
        hessian4_fn <- .cse_fn(
            lapply(hess4_triplets, `[[`, "ast"),
            paste0("v[", seq_len(n_hess4), "L]"),
            paste0("v <- numeric(", n_hess4, "L)"), "v",
            endo, exo, params,
            local_vars = local_vars)
    } else {
        hessian4_fn <- function(dy, params, ss = numeric(0)) numeric(0L)
    }

    # -----------------------------------------------------------------------
    # 6d. Build symbolic fifth-order Hessian function
    #
    # Returns a numeric vector of length n_hess5.  Canonical order c1â‰¤...â‰¤c5.
    # -----------------------------------------------------------------------
    if (n_hess5 > 0L) {
        hessian5_fn <- .cse_fn(
            lapply(hess5_triplets, `[[`, "ast"),
            paste0("v[", seq_len(n_hess5), "L]"),
            paste0("v <- numeric(", n_hess5, "L)"), "v",
            endo, exo, params,
            local_vars = local_vars)
    } else {
        hessian5_fn <- function(dy, params, ss = numeric(0)) numeric(0L)
    }

    # ---- Drop per-triplet ASTs from the stored object ----------------------
    # Every compile-time consumer of the triplet $ast (CSE codegen of the
    # *_fn functions above, the Jacobian tape, and the parameter-derivative
    # codegen) has already run. The $ast trees are NOT read at solve/gradient
    # time (only the index fields eq/col*/row/col and the compiled functions
    # are), yet they dominate the compiled object's footprint -- e.g. Basu
    # order-3 hess3 holds 1739 fully-expanded 3rd-derivative trees at ~470 MB.
    # Drop them so the compiled model is small in RAM and fast to (de)serialise
    # for the on-disk cache.
    .drop_field <- function(tl, field)
        lapply(tl, function(t) { t[[field]] <- NULL; t })
    jac_triplets   <- .drop_field(jac_triplets,   "ast")
    hess2_triplets <- .drop_field(hess2_triplets, "ast")
    hess3_triplets <- .drop_field(hess3_triplets, "ast")
    hess4_triplets <- .drop_field(hess4_triplets, "ast")
    hess5_triplets <- .drop_field(hess5_triplets, "ast")
    # Same for the parameter-derivative triplets ($ast used only to build
    # param_jacobian_fn / param_hessian2_fn / param2_jacobian_fn via .cse_fn
    # above; gradient consumers read only the index fields eq/col*/param/pa/pb).
    param_jac_triplets   <- .drop_field(param_jac_triplets,   "ast")
    param_hess2_triplets <- .drop_field(param_hess2_triplets, "ast")
    param2_jac_triplets  <- .drop_field(param2_jac_triplets,  "ast")

    list(
        max_order      = max_order,
        residuals_fn   = residuals_fn,
        jacobian_fn    = jacobian_fn,
        jac_tape       = jac_tape,
        hess2_tape     = hess2_tape,
        hess3_tape     = hess3_tape,
        jac_tape_param_names = params,
        hessian2_fn    = hessian2_fn,
        hess2_triplets = hess2_triplets,
        n_hess         = n_hess,
        hessian3_fn    = hessian3_fn,
        hess3_triplets = hess3_triplets,
        n_hess3        = n_hess3,
        hessian4_fn    = hessian4_fn,
        hess4_triplets = hess4_triplets,
        n_hess4        = n_hess4,
        hessian5_fn    = hessian5_fn,
        hess5_triplets = hess5_triplets,
        n_hess5        = n_hess5,
        residual_asts  = residual_asts,
        jac_triplets   = jac_triplets,
        param_jacobian_fn  = param_jacobian_fn,
        param_jac_triplets = param_jac_triplets,
        param_deriv_ok     = param_deriv_ok,
        hessian2_built     = build_hess2,
        hessian3_built     = build_hess3up,
        param_hessian2_fn  = param_hessian2_fn,
        param_hess2_triplets = param_hess2_triplets,
        param_hess2_built  = build_param_hess2,
        param2_jacobian_fn = param2_jacobian_fn,
        param2_jac_triplets = param2_jac_triplets,
        param2_jac_built   = build_param2_jac,
        param_deriv2_ok    = build_param2_jac && build_param_hess2,
        n_params           = n_params,
        dyn_col_map    = dyn_col_map,
        n_eq           = n_eq,
        n_dyn_cols     = n_dyn_cols,
        n_exo          = n_exo,
        total_cols     = total_cols,
        endo_names     = endo,
        exo_names      = exo
    )
}
