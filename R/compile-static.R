## R/compile-static.R
## --------------------------------------------------------------------------
## Build the STATIC model: collapse all timing to t, differentiate
## residuals with respect to endogenous variables at t, generate the
## residuals_fn and jacobian_fn closures used by the steady-state solver.
##
## Phase-1 split from jacobian-monolith.R (no logic changes). Phase-1d will
## extract a shared _build_residuals_asts() helper between static + dynamic
## (internal).
## --------------------------------------------------------------------------

## Build a static-model closure from a list of ASTs using CSE.
##
## Mirrors .cse_fn in compile-dynamic.R but uses the static function signature
## function(y, x, params, ss = y) and use_timing = FALSE.
##
## Using CSE (ast_cse_emit) breaks each equation into a sequence of shallow
## temporary assignments, so no single parsed statement is deeply nested.
## This prevents R's C-level parser from hitting "contextstack overflow" on
## models with long chained expressions (e.g. CMR2014 eq221 with 40 nested
## LEAD product terms whose Jacobian exceeds parse depth ~48).
##
## `asts`      -- list of AST nodes (NULL entries are skipped)
## `lhs`       -- character vector, one LHS assignment string per AST entry
##               (e.g. "r[1L]" or "J[2,3]")
## `alloc`     -- single line allocating the output object
## `ret`       -- name of the object to return
## `endo`, `exo`, `params` -- name vectors for codegen
## `local_vars`-- named list of #-local variable ASTs (optional). When supplied,
##               each local is emitted as a `LOCAL_name <- <expr>` let-binding
##               at the top of the generated function, in topological order. This
##               avoids substituting locals into the ASTs (preventing exponential
##               blowup in models with deep local chains like US_IN10).
.cse_static_fn <- function(asts, lhs, alloc, ret, endo, exo, params,
                            local_vars = list()) {
    em   <- ast_cse_emit(asts, endo, exo, params, use_timing = FALSE)
    keep <- !is.na(em$refs)
    ## NB: guard on any(keep). `paste0("    ", lhs[keep], " <- ", em$refs[keep])`
    ## with an all-FALSE/empty `keep` does NOT yield character(0) -- paste0
    ## recycles the zero-length `lhs[keep]`/`em$refs[keep]` against the constant
    ## "    "/" <- " strings to length 1, emitting a spurious "     <- " line
    ## (a syntax error). This happens whenever a generated function has no
    ## assignable entries, e.g. the param-residual Jacobian of a model whose
    ## only parameters are shock standard deviations (unit-root / local-level
    ## models): all dResidual/dparam are zero, so `asts` is empty.
    assign_lines <- if (any(keep))
        paste0("    ", lhs[keep], " <- ", em$refs[keep]) else character(0)
    ## Emit #-local let-bindings in topological order (each local's raw definition
    ## references only params and prior locals, so the order matters).
    ## ast_to_fn_body on a local's raw AST emits LOCAL_prior_name references that
    ## are resolved by earlier let-bindings -- no substitution needed.
    local_lines <- if (length(local_vars) > 0L) {
        ln  <- names(local_vars)
        vapply(seq_along(local_vars), function(k) {
            rhs <- ast_to_fn_body(local_vars[[k]], endo_names = endo,
                                  exo_names = exo, param_names = params)
            paste0("    LOCAL_", ln[k], " <- ", rhs)
        }, character(1))
    } else character(0)
    txt  <- paste(c(
        "function(y, x, params, ss = y) {",
        paste0("    ", alloc),
        local_lines,
        if (length(em$decls)) paste0("    ", em$decls) else character(0),
        assign_lines,
        paste0("    ", ret),
        "}"
    ), collapse = "\n")
    gen <- eval(parse(text = txt))
    environment(gen) <- new.env(parent = asNamespace("dynhr"))
    ## Apply the same JIT guard as .cse_fn in compile-dynamic.R: large generated
    ## static closures (huge Jacobians in models like US_IN10) take minutes to
    ## byte-compile but run fine interpreted. Keep them interpreted.
    .interp_only_if_large(gen, txt)
}


## Pre-compute the derivative of each #-local variable w.r.t. each parameter
## using a FORWARD CHAIN-RULE SWEEP over the local variable definitions.
##
## Models with deep chains of #-local definitions (e.g. US_IN10 where
## BB_SS -> llb -> llh -> ... expands to ~800 KB) produce enormous substituted
## ASTs. Differentiating a 3.7 MB merged residual w.r.t. 49 params × 76 eqs
## takes ~100s per large equation. The forward sweep computes d(LOCAL)/d(param)
## once for each (local, param) pair using only the raw (unexpanded) local
## definitions, avoiding the exponential blowup.
##
## Algorithm: for each local L_k (topological order, k = 1..m):
##   d(L_k)/d(θ) = ∂L_k/∂θ|_raw               [explicit param in L_k definition]
##               + Σ_{j<k} (∂L_k/∂LOCAL_j) * d(L_j)/d(θ)  [chain rule through prior locals]
##
## Returns a list of length n_locals, each a list of length n_params.
## dL_dtheta[[k]][[p]] = AST of d(L_k)/d(θ_p), NULL means zero.
##
.local_param_sweep <- function(local_vars, params) {
    local_names <- names(local_vars)
    n_local <- length(local_names)
    n_param  <- length(params)
    if (n_local == 0L || n_param == 0L) return(list())

    dL_dtheta <- vector("list", n_local)
    for (k in seq_len(n_local)) {
        dL_dtheta[[k]] <- vector("list", n_param)
        lk_raw <- local_vars[[k]]

        for (p in seq_len(n_param)) {
            ## Direct partial: ∂L_k/∂θ_p (treating all LOCAL_ as 0)
            direct <- ast_differentiate(lk_raw, params[p], 0L, wrt = "parameter")
            acc <- if (!ast_is_zero(direct)) direct else NULL

            ## Chain rule through prior locals
            for (j in seq_len(k - 1L)) {
                ## ∂L_k/∂LOCAL_j (treating LOCAL_j as an independent variable)
                partial_kj <- ast_differentiate(lk_raw, local_names[j], wrt = "local")
                if (ast_is_zero(partial_kj)) next
                dLj <- dL_dtheta[[j]][[p]]
                if (is.null(dLj)) next      # dL_j/dθ_p = 0
                term <- list(type = "binop", op = "*",
                             left = partial_kj, right = dLj)
                acc <- if (is.null(acc)) term else
                    list(type = "binop", op = "+", left = acc, right = term)
            }
            ## NB: do NOT assign NULL -- assigning NULL to a list slot removes
            ## it (shrinks the list), breaking later index access. Only assign
            ## when non-null; the pre-allocated NULL serves as the "zero" sentinel.
            if (!is.null(acc)) dL_dtheta[[k]][[p]] <- acc
        }
    }
    dL_dtheta
}

## Apply the forward chain-rule to compute d(expr_raw)/d(θ_p) for a raw (pre-
## substitution) expression that may reference #-local variables.
## `expr_raw`   -- pre-substitution AST (contains local_variable nodes)
## `param`      -- parameter name string
## `p_idx`      -- integer index of param in local_names order
## `local_names`-- character vector of local names (same order as dL_dtheta)
## `dL_dtheta`  -- pre-computed matrix from .local_param_sweep()
## Returns an AST (possibly NULL = zero).
.chain_rule_param_deriv <- function(expr_raw, param, p_idx, local_names, dL_dtheta) {
    ## Direct partial: d(expr)/d(θ) treating all LOCAL_ as constants (= 0)
    direct <- ast_differentiate(expr_raw, param, 0L, wrt = "parameter")
    acc <- if (!ast_is_zero(direct)) direct else NULL

    ## Chain rule: Σ_k (∂expr/∂LOCAL_k) * (dL_k/dθ)
    for (k in seq_along(local_names)) {
        dLk <- dL_dtheta[[k]][[p_idx]]
        if (is.null(dLk)) next  # dL_k/dθ = 0
        ## ∂expr/∂LOCAL_k
        partial_k <- ast_differentiate(expr_raw, local_names[k], wrt = "local")
        if (ast_is_zero(partial_k)) next
        term <- list(type = "binop", op = "*", left = partial_k, right = dLk)
        acc <- if (is.null(acc)) term else
            list(type = "binop", op = "+", left = acc, right = term)
    }
    acc  # NULL = zero
}


#' Build static model residual and Jacobian functions
#'
#' In the static model, all leads and lags are collapsed to time t.
#' This is used for steady-state computation.
#'
#' @param model            A dynhr_mod object from parse_mod().
#' @param want_param_deriv Logical. When TRUE (default) the symbolic
#'                         parameter-residual Jacobian is compiled so the
#'                         analytic gradient can compute the steady-state
#'                         sensitivity. When FALSE the step is skipped;
#'                         the gradient layer falls back to finite differences.
#' @param want_param_deriv2 Logical. When TRUE the SECOND-order static tensors
#'                         (static_hess2, static_param_jac, static_param2) are
#'                         compiled so the analytic exact-posterior-Hessian path
#'                         can compute the second-order steady-state sensitivity
#'                         d2ys. Off by default (order-2 codegen is expensive and
#'                         only the exact Hessian needs it).
#' @return A list with residuals_fn, jacobian_fn, and expression strings.
#' @noRd
build_static_model <- function(model, want_param_deriv = TRUE,
                               want_param_deriv2 = FALSE) {
    endo   <- model$var_names
    exo    <- model$varexo_names
    params <- model$param_names
    n_eq   <- length(model$equations)
    n_endo <- length(endo)

    # 1. Build residuals: lhs - rhs, collapse timing. NO local substitution.
    #
    # Performance (I20): models with deep #-local chains (e.g. US_IN10 where
    # BB_SS -> llb -> llh -> ... expands to ~800 KB) produce enormous ASTs after
    # ast_substitute_locals. Differentiating a 3.7 MB residual w.r.t. 76 endo
    # variables takes ~140s for one equation; on the raw 400-char pre-sub AST it
    # takes 0.04s. The param-deriv step has the same exponential-blowup problem.
    #
    # Fix: keep residuals in raw (pre-substitution) form throughout differentiation.
    # Generated closures receive #-locals as named let-bindings at the top of the
    # function (emitted by .cse_static_fn with local_vars=), so LOCAL_x references
    # in the raw ASTs resolve correctly at runtime without ever expanding them.
    have_locals <- length(model$local_variables) > 0L
    local_vars  <- model$local_variables          # named list (possibly empty)

    residual_asts <- vector("list", n_eq)
    for (i in seq_len(n_eq)) {
        eq  <- model$equations[[i]]
        res <- equation_to_residual(eq)
        res <- ast_collapse_timing(res)
        res <- ast_simplify(res)
        residual_asts[[i]] <- res
    }

    # I20 correctness guard (mirrors compile-dynamic.R): the raw-AST fast path
    # treats LOCAL_x as derivative 0 w.r.t. endo vars, valid ONLY when no #-local
    # references a variable. A local like `#wedge = q - tot` (nz_base_nonlinear)
    # makes the static Jacobian d(residual)/d(endo) wrong unless we chain-rule
    # through it — so substitute locals and disable the let-binding path here.
    if (have_locals && !.locals_are_param_only(local_vars)) {
        residual_asts <- lapply(residual_asts, function(r)
            ast_simplify(ast_substitute_locals(r, local_vars)))
        local_vars  <- list()
        have_locals <- FALSE
    }

    # 2. Build Jacobian ASTs: d(residual_i)/d(endo_j).
    #    Differentiate on the raw (pre-sub) ASTs. LOCAL_x nodes have derivative
    #    0 w.r.t. endogenous variables (they represent steady-state quantities,
    #    i.e. functions of params only, constant w.r.t. endo vars at time t).
    #    Non-zero Jacobian ASTs keep their LOCAL_x references; the generated
    #    closure resolves them from the let-bindings emitted by .cse_static_fn.
    jacobian_asts <- matrix(vector("list", 1), nrow = n_eq, ncol = n_endo)
    for (i in seq_len(n_eq)) {
        for (j in seq_len(n_endo)) {
            jac <- ast_differentiate(residual_asts[[i]], endo[j], 0L)
            jacobian_asts[i, j] <- list(jac)
        }
    }

    # 3. Convert to R expression strings (metadata only).
    #    Use substituted residuals for the human-readable strings (jac_exprs)
    #    so they show the full expressions. Build substituted form lazily here.
    residual_asts_sub <- if (have_locals)
        lapply(residual_asts, function(r) ast_substitute_locals(r, local_vars))
    else residual_asts

    res_exprs <- vapply(residual_asts_sub, ast_to_fn_body, character(1),
                        endo_names = endo, exo_names = exo,
                        param_names = params)

    jac_exprs <- matrix("0", nrow = n_eq, ncol = n_endo)
    for (i in seq_len(n_eq)) {
        for (j in seq_len(n_endo)) {
            jac_exprs[i, j] <- ast_to_fn_body(jacobian_asts[[i, j]],
                                               endo, exo, params)
        }
    }
    colnames(jac_exprs) <- endo
    rownames(jac_exprs) <- paste0("eq", seq_len(n_eq))

    # 4. Build residuals function via CSE.
    #    Pass local_vars so .cse_static_fn emits LOCAL_x let-bindings before the
    #    residual expressions -- raw ASTs reference LOCAL_x, which resolve at runtime.
    residuals_fn <- .cse_static_fn(
        asts   = residual_asts,
        lhs    = paste0("r[", seq_len(n_eq), "L]"),
        alloc  = paste0("r <- numeric(", n_eq, "L)"),
        ret    = "r",
        endo   = endo, exo = exo, params = params,
        local_vars = local_vars
    )

    # 5. Build Jacobian function via CSE.
    #    Collect only the non-zero (i,j) entries.
    jac_asts_flat <- vector("list", n_eq * n_endo)
    jac_lhs_flat  <- character(n_eq * n_endo)
    k <- 0L
    for (i in seq_len(n_eq)) {
        for (j in seq_len(n_endo)) {
            if (!ast_is_zero(jacobian_asts[[i, j]])) {
                k <- k + 1L
                jac_asts_flat[[k]] <- jacobian_asts[[i, j]]
                jac_lhs_flat[[k]]  <- paste0("J[", i, ",", j, "]")
            }
        }
    }
    if (k > 0L) {
        jac_asts_flat <- jac_asts_flat[seq_len(k)]
        jac_lhs_flat  <- jac_lhs_flat[seq_len(k)]
    } else {
        jac_asts_flat <- list()
        jac_lhs_flat  <- character(0)
    }
    jacobian_fn <- .cse_static_fn(
        asts   = jac_asts_flat,
        lhs    = jac_lhs_flat,
        alloc  = paste0("J <- matrix(0, nrow = ", n_eq, "L, ncol = ", n_endo, "L)"),
        ret    = "J",
        endo   = endo, exo = exo, params = params,
        local_vars = local_vars
    )

    # -----------------------------------------------------------------------
    # 6. Parameter-residual Jacobian  ∂F_static_i/∂θ_k  (Tier 11 #3)
    #
    # Explicit partial of each static residual w.r.t. each parameter. Used by
    # the solution-derivative layer to solve the steady-state sensitivity
    #   J_static · (dȳ/dθ_k) = -(∂F_static/∂θ_k)
    # analytically, replacing the two warm-started Newton re-solves per
    # parameter. NULL when not smoothly parameter-differentiable.
    #
    # Performance: if there are #-locals, we use the chain-rule forward sweep
    # (.local_param_sweep + .chain_rule_param_deriv) to avoid ever working
    # with the huge substituted ASTs. The sweep pre-computes d(LOCAL_k)/d(θ)
    # for all (k, θ) from small pre-substitution local definitions, then each
    # residual is differentiated on its raw (pre-sub) form using the chain rule.
    # -----------------------------------------------------------------------
    n_params <- length(params)
    local_names <- names(local_vars)
    # param_deriv_ok: model is analytically parameter-differentiable AND the
    # caller has not opted out via want_param_deriv = FALSE.
    # STEADY_STATE and max/min appear in the pre-sub form; local_variable nodes
    # are safe (they're treated as constants in .ast_param_deriv_safe).
    param_deriv_ok <- isTRUE(want_param_deriv) &&
        all(vapply(residual_asts, .ast_param_deriv_safe, logical(1)))

    if (param_deriv_ok && n_params > 0L) {
        ## Pre-compute d(LOCAL_k)/d(θ) via forward chain-rule sweep (once, fast).
        ## This avoids ever differentiating the huge substituted residuals.
        dL_dtheta <- if (have_locals)
            .local_param_sweep(local_vars, params) else list()

        pres_asts <- vector("list", n_eq * n_params)
        pres_lhs  <- character(n_eq * n_params)
        kk <- 0L
        for (i in seq_len(n_eq)) {
            for (k_p in seq_len(n_params)) {
                ## Chain rule on the raw (pre-sub) residual. The result may
                ## contain LOCAL_x references which .cse_static_fn resolves via
                ## its let-bindings (local_vars=); no substitution needed.
                dp <- if (have_locals)
                    .chain_rule_param_deriv(residual_asts[[i]], params[k_p],
                                            k_p, local_names, dL_dtheta)
                else
                    ast_differentiate_param(residual_asts[[i]], params[k_p])
                if (!is.null(dp) && !ast_is_zero(dp)) {
                    kk <- kk + 1L
                    pres_asts[[kk]] <- dp
                    pres_lhs[[kk]]  <- paste0("R[", i, ",", k_p, "]")
                }
            }
        }
        if (kk > 0L) {
            pres_asts <- pres_asts[seq_len(kk)]
            pres_lhs  <- pres_lhs[seq_len(kk)]
        } else {
            pres_asts <- list()
            pres_lhs  <- character(0)
        }
        param_resid_fn <- .cse_static_fn(
            asts   = pres_asts,
            lhs    = pres_lhs,
            alloc  = paste0("R <- matrix(0, nrow = ", n_eq, "L, ncol = ", n_params, "L)"),
            ret    = "R",
            endo   = endo, exo = exo, params = params,
            local_vars = local_vars
        )
    } else {
        param_resid_fn <- NULL
    }

    # -----------------------------------------------------------------------
    # 7. Second-order STATIC tensors  (Tier 11 #3, order-2 layer = 3b)
    #
    # Three sparse symbolic tensors of the static residual F_s(ȳ,θ) needed to
    # solve the SECOND-order steady-state sensitivity d2ys = d²ȳ/(dθ_a dθ_b)
    # by the implicit function theorem (see .analytic_d2ys):
    #
    #   static_hess2     H_s[i,p,q] = ∂²F_s_i/(∂y_p ∂y_q)  (canonical p<=q)
    #   static_param_jac B_s[i,p,k] = ∂²F_s_i/(∂y_p ∂θ_k)
    #   static_param2    C_s[i,a,b] = ∂²F_s_i/(∂θ_a ∂θ_b)  (canonical a<=b)
    #
    # Differentiating the static residual ASTs (all timing collapsed to ll=0):
    # H_s by differentiating the Jacobian ASTs once more in y; B_s by
    # differentiating the Jacobian ASTs in the parameter direction; C_s by
    # differentiating each residual twice in the parameter direction.
    # Built only when explicitly requested (the exact-Hessian path) AND the
    # model is smoothly parameter-differentiable; the consumer gates on
    # static_param2_built and falls back to FD otherwise.
    # -----------------------------------------------------------------------
    build_static_p2 <- isTRUE(want_param_deriv2) && param_deriv_ok &&
        n_params > 0L
    static_hess2_triplets   <- list()
    static_param_jac_triplets <- list()
    static_param2_triplets  <- list()
    if (build_static_p2) {
        for (i in seq_len(n_eq)) {
            ## H_s[i,p,q]: differentiate Jacobian AST (∂F_i/∂y_p) again in y_q,
            ## canonical p <= q (the static Hessian is symmetric in (p,q)).
            for (p in seq_len(n_endo)) {
                dp <- jacobian_asts[[i, p]]
                if (ast_is_zero(dp)) next
                for (q in seq(p, n_endo)) {
                    d2 <- ast_differentiate(dp, endo[q], 0L)
                    if (!ast_is_zero(d2)) {
                        static_hess2_triplets <- c(static_hess2_triplets,
                            list(list(eq = i, p = p, q = q,
                                      expr = ast_to_fn_body(d2, endo, exo, params),
                                      ast  = d2)))
                    }
                }
                ## B_s[i,p,k]: differentiate ∂F_i/∂y_p in the parameter direction.
                for (k in seq_len(n_params)) {
                    db <- if (have_locals)
                        .chain_rule_param_deriv(dp, params[k], k, local_names, dL_dtheta)
                    else ast_differentiate_param(dp, params[k])
                    if (!is.null(db) && !ast_is_zero(db)) {
                        static_param_jac_triplets <- c(static_param_jac_triplets,
                            list(list(eq = i, p = p, param = k,
                                      expr = ast_to_fn_body(db, endo, exo, params),
                                      ast  = db)))
                    }
                }
            }
            ## C_s[i,a,b]: differentiate the residual twice in the parameter
            ## direction, canonical a <= b (symmetric in (a,b)).
            for (a in seq_len(n_params)) {
                d1 <- if (have_locals)
                    .chain_rule_param_deriv(residual_asts[[i]], params[a], a,
                                            local_names, dL_dtheta)
                else ast_differentiate_param(residual_asts[[i]], params[a])
                if (is.null(d1) || ast_is_zero(d1)) next
                for (b in seq(a, n_params)) {
                    d2 <- if (have_locals)
                        .chain_rule_param_deriv(d1, params[b], b, local_names, dL_dtheta)
                    else ast_differentiate_param(d1, params[b])
                    if (!is.null(d2) && !ast_is_zero(d2)) {
                        static_param2_triplets <- c(static_param2_triplets,
                            list(list(eq = i, pa = a, pb = b,
                                      expr = ast_to_fn_body(d2, endo, exo, params),
                                      ast  = d2)))
                    }
                }
            }
        }
    }
    n_static_hess2     <- length(static_hess2_triplets)
    n_static_param_jac <- length(static_param_jac_triplets)
    n_static_param2    <- length(static_param2_triplets)

    ## Closure builder: sparse value vector parallel to a triplet list.
    ## Uses CSE to avoid contextstack overflow on large higher-order models.
    ## Passes local_vars so LOCAL_x references in the ASTs resolve at runtime.
    .build_static_sparse_fn <- function(triplets, n) {
        if (n == 0L)
            return(function(y, x, params, ss = y) numeric(0L))
        t_asts <- lapply(triplets, `[[`, "ast")
        t_lhs  <- paste0("v[", seq_len(n), "L]")
        .cse_static_fn(
            asts   = t_asts,
            lhs    = t_lhs,
            alloc  = paste0("v <- numeric(", n, "L)"),
            ret    = "v",
            endo   = endo, exo = exo, params = params,
            local_vars = local_vars
        )
    }
    if (build_static_p2) {
        static_hess2_fn     <- .build_static_sparse_fn(static_hess2_triplets,
                                                       n_static_hess2)
        static_param_jac_fn <- .build_static_sparse_fn(static_param_jac_triplets,
                                                       n_static_param_jac)
        static_param2_fn    <- .build_static_sparse_fn(static_param2_triplets,
                                                       n_static_param2)
    } else {
        static_hess2_fn     <- NULL
        static_param_jac_fn <- NULL
        static_param2_fn    <- NULL
    }

    list(
        residuals_fn   = residuals_fn,
        jacobian_fn    = jacobian_fn,
        param_resid_fn = param_resid_fn,
        param_deriv_ok = param_deriv_ok,
        static_hess2_fn        = static_hess2_fn,
        static_hess2_triplets  = static_hess2_triplets,
        static_param_jac_fn    = static_param_jac_fn,
        static_param_jac_triplets = static_param_jac_triplets,
        static_param2_fn       = static_param2_fn,
        static_param2_triplets = static_param2_triplets,
        static_param2_built    = build_static_p2,
        n_params       = n_params,
        residual_asts  = residual_asts,
        jacobian_asts  = jacobian_asts,
        residual_exprs = res_exprs,
        jacobian_exprs = jac_exprs,
        n_eq           = n_eq,
        n_endo         = n_endo,
        endo_names     = endo
    )
}
