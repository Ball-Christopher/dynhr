## R/gradient-primitive-deriv.R
## --------------------------------------------------------------------------
## ANALYTIC (symbolic) parameter-derivatives of the smooth model primitives
## (Tier 11 #3). Replaces the central finite differences that the
## solution-derivative layers use for
##
##   dys      = dȳ/dθ_k                         (steady-state sensitivity)
##   df_*_k   = d(f_plus,f_zero,f_minus,f_exo)/dθ_k   (dynamic-Jacobian blocks)
##
## with finite-difference-free symbolic values built from the compile-time
## parameter-Jacobian (∂²F/∂w∂θ) and the existing symbolic model Hessian
## (∂²F/∂w∂w), making the implicit-differentiation gradient fully analytic.
##
## THE MATH
## --------
## The dynamic-Jacobian primitive is J(dy(ȳ(θ)), θ, ss(ȳ(θ))) with the
## dynamic point dy[c] = ȳ[var(c)] replicating the steady state across the
## lead/lag columns. With no STEADY_STATE() reference in the equations (gated)
## the total derivative is
##
##   dJ[i,c]/dθ_k = ∂²F_i/(∂w_c ∂θ_k)                         (explicit)
##                + Σ_{c'} ∂²F_i/(∂w_c ∂w_{c'}) · ddy[c']/dθ_k  (steady-state chain)
##
## where ddy[c']/dθ_k = dys[var(c')] for endogenous columns and 0 for the
## (constant, zero) shock columns. The first term is the compiled
## param_jacobian_fn; the second contracts the compiled model Hessian
## (hessian2_fn) with dys.
##
## The steady state ȳ solves F_static(ȳ, θ) = 0, so by the implicit function
## theorem
##
##   J_static · (dȳ/dθ_k) = -(∂F_static/∂θ_k)
##
## with J_static the static Jacobian and ∂F_static/∂θ_k the compiled
## param_resid_fn — one factorization reused across all parameters.
## --------------------------------------------------------------------------


#' Can the analytic primitive-derivative path be used for this model?
#'
#' Requires that both the dynamic and static compilers produced parameter
#' derivatives (no STEADY_STATE() reference, no non-differentiable construct,
#' at least one parameter). Runtime conditions (finite base Jacobian,
#' non-singular static Jacobian) are checked when the path runs and trigger a
#' graceful fall back to finite differences.
#'
#' @param compiled dynhr_compiled.
#' @return TRUE if the symbolic parameter-derivative closures are available.
#' @noRd
.can_use_analytic_primitive_deriv <- function(compiled) {
  dyn <- compiled$dynamic
  sta <- compiled$static
  isTRUE(dyn$param_deriv_ok) && isTRUE(sta$param_deriv_ok) &&
    !is.null(dyn$param_jacobian_fn) && !is.null(sta$param_resid_fn) &&
    (dyn$n_params %||% 0L) > 0L &&
    ## The steady-state chain rule contracts the symbolic model Hessian
    ## (hessian2_fn) with dys; without it the analytic derivative silently
    ## drops that term. hessian2_built is TRUE iff the second-order Hessian
    ## was actually compiled (max_order >= 2 OR param_deriv requested).
    isTRUE(dyn$hessian2_built)
    ## NB: steady_state_model-derived-parameter models (.ssm_assigns_param) are
    ## handled at FIRST order by the augmented-Jacobian dys + computed-parameter
    ## dynamic chain (Tier 12 #2). caldara_rp still falls back to FD because its
    ## augmented static Jacobian is singular (runtime qr.solve -> NULL). The
    ## SECOND-order analytic paths self-guard on !.ssm_assigns_param (the order-2
    ## computed-parameter chain is not yet derived) and use FD for such models.
}


#' Does the steady_state_model block compute (assign to) a parameter?
#'
#' When the analytic steady-state block sets a declared parameter as a function
#' of the others, perturbing a free parameter also moves that computed one --
#' a dependence the explicit param_resid_fn (which holds all parameters fixed)
#' would miss, biasing the analytic dys.
#'
#' Tier 12 #2 closes this at FIRST order with the augmented-Jacobian chain: solve
#' \eqn{J_{aug}\,dys = -rhs_{aug}} with
#' \eqn{J_{aug} = J_{static} + \sum_{p_c}(\partial F/\partial p_c)(\partial g_{p_c}/\partial y)}
#' and \eqn{rhs_k = \partial F/\partial\theta_k + \sum_{p_c}(\partial F/\partial p_c)(\partial g_{p_c}/\partial\theta_k)}
#' in \code{.analytic_dys}, plus the matching computed-parameter dynamic-Jacobian
#' term \eqn{\sum_{p_c}(\partial^2F/\partial w\,\partial p_c)(dp_c/d\theta_k)} in
#' \code{.analytic_dprimitives} (channel 3). The SSM total-derivative building
#' blocks (\eqn{\partial g_{p_c}/\partial y}, \eqn{\partial g_{p_c}/\partial\theta},
#' chained through the assignment order) come from \code{.ssm_param_chain_derivs()}
#' (verified exact). Validated end-to-end on \code{rbc2shock_ssm} against the full
#' re-solve pipeline (re-solve steady -> re-derive p_c -> re-solve perturbation).
#'
#' This function is therefore NOT a feature gate any more; it is the predicate the
#' augmented branches key on. Two models still fall back to FD, by design:
#'  (1) caldara_rp -- its static Jacobian is *numerically singular* at gamma = 40
#'      (the \eqn{s = V^{1-\gamma}} column is annihilated, cond ~1e19); the
#'      \code{nu}-chain rank-1 update does NOT lift that singular direction, so
#'      \eqn{J_{aug}} is still singular and \code{.analytic_dys}'s qr.solve returns
#'      NULL at run time -> exact FD re-solve. (2) SECOND-order derivatives of any
#'      SSM-param model: the order-2 computed-parameter chain is not yet derived,
#'      so the order-2 paths self-guard on \code{!.ssm_assigns_param} and use FD.
#' See test-ssm-param-chain.R / test-rbc2shock-ssm.R for the codified evidence.
#'
#' @param model dynhr_mod.
#' @return TRUE if any steady_state_model assignment targets a parameter name.
#' @noRd
.ssm_assigns_param <- function(model) {
  ssm <- model$steady_state_model
  if (is.null(ssm) || length(ssm) == 0L) return(FALSE)
  lhs <- vapply(ssm, function(a) a$name %||% NA_character_, character(1))
  any(lhs %in% model$param_names)
}


#' Re-derive steady_state_model-computed parameters consistent with `params`.
#'
#' The finite-difference primitive path perturbs a FREE parameter and re-solves
#' the steady state, but \code{solve_steady()} returns only the steady VALUES --
#' it does not propagate the re-derived computed parameter \eqn{p_c} back into the
#' parameter vector. Evaluating \code{extract_system_matrices()} with the stale
#' \eqn{p_c} (its base value) silently biases the FD dynamic-Jacobian derivative
#' for steady_state_model-derived-parameter models (the same stale-\eqn{p_c} blind
#' spot the analytic chain in \code{.analytic_dprimitives} closes). This helper
#' re-evaluates the SSM block at \code{params} and overwrites every SSM-assigned
#' name, so the FD path uses a self-consistent \eqn{(\bar y, \theta, p_c)} point.
#' A no-op (returns \code{params} unchanged) for non-SSM-parameter models.
#'
#' @param model  dynhr_mod.
#' @param params Named numeric parameter vector (a free parameter already perturbed).
#' @return \code{params} with the SSM-computed parameter(s) re-derived.
#' @noRd
.ssm_consistent_params <- function(model, params) {
  if (!.ssm_assigns_param(model)) return(params)
  a <- tryCatch(solve_steady_state_analytical(model, params), error = function(e) NULL)
  if (is.null(a) || is.null(a$params)) return(params)
  common <- intersect(names(params), names(a$params))
  params[common] <- a$params[common]
  params
}


#' Collect the leaf symbol names (variables + parameters) referenced in an AST.
#'
#' @param node AST node.
#' @return character vector of names (with duplicates).
#' @noRd
.ast_leaf_syms <- function(node) {
  if (is.null(node)) return(character(0))
  switch(node$type,
    "number"         = character(0),
    "variable"       = node$name,
    "parameter"      = node$name,
    "local_variable" = node$name,
    "binop"   = c(.ast_leaf_syms(node$left), .ast_leaf_syms(node$right)),
    "unaryop" = .ast_leaf_syms(node$operand),
    "funcall" = unlist(lapply(node$args, .ast_leaf_syms), use.names = FALSE),
    character(0))
}


#' Total parameter/variable derivatives of steady_state_model-computed parameters.
#'
#' Walks the SSM assignment list IN ORDER, carrying forward, for every assigned
#' name, its total derivative w.r.t. the coordinate set \code{c(endo_vars,
#' free_params)} (free = declared parameters NOT themselves computed by the SSM).
#' Each statement \eqn{v := g(\dots)} applies the chain rule numerically at
#' \eqn{(\bar y, \theta)}:
#' \deqn{\frac{dv}{dx} = \sum_{s \in \mathrm{syms}(g)} \frac{\partial g}{\partial s}\Big|_{(\bar y,\theta)} \cdot \frac{ds}{dx}}
#' where \eqn{ds/dx} is the unit vector for a still-free coordinate \eqn{s}, or
#' the already-accumulated total derivative if \eqn{s} was assigned by an earlier
#' SSM statement. \eqn{\partial g/\partial s} is the symbolic
#' \code{ast_differentiate} of the assignment RHS (wrt a variable when \eqn{s} is
#' endogenous, wrt a parameter otherwise) evaluated at the point.
#'
#' Returns the derivative rows for the SSM-COMPUTED PARAMETERS only -- i.e. the
#' \eqn{\partial g_{p_c}/\partial y} and \eqn{\partial g_{p_c}/\partial\theta}
#' the augmented-Jacobian dys chain needs.
#'
#' @param compiled dynhr_compiled.
#' @param ys       Named steady state (declaration order).
#' @param params   Named parameter vector (incl. the SSM-computed values).
#' @return list with
#'   \code{dg_dy}    [n_computed x n_endo] (rows = computed params, cols = endo),
#'   \code{dg_dtheta}[n_computed x n_params] (cols indexed by model$param_names),
#'   \code{computed} (names of SSM-computed parameters),
#'   or NULL if the SSM assigns no parameter or contains a non-differentiable RHS.
#' @noRd
.ssm_param_chain_derivs <- function(compiled, ys, params) {
  model <- compiled$model
  ssm   <- model$steady_state_model
  if (is.null(ssm) || length(ssm) == 0L) return(NULL)
  endo  <- model$var_names
  pars  <- model$param_names

  assigned <- vapply(ssm, function(a) a$name %||% NA_character_, character(1))
  computed <- intersect(assigned, pars)
  if (length(computed) == 0L) return(NULL)

  ## Free coordinates: every endo variable + every NON-computed parameter.
  free_params <- setdiff(pars, computed)
  coords <- c(endo, free_params)
  ncoord <- length(coords)

  ## ast_eval needs variables keyed "<name>__0"; SSM RHS references lag/lead 0.
  vv <- setNames(as.numeric(ys[endo]), paste0(endo, "__0"))

  ## Accumulated total-derivative map per name (over coords). Unassigned names
  ## act as free coordinates (unit vector) or constants (zero).
  dmap <- vector("list", 0L)
  base_unit <- function(nm) {
    v <- numeric(ncoord); names(v) <- coords
    if (nm %in% coords) v[[nm]] <- 1
    v
  }
  get_d <- function(nm) {
    d <- dmap[[nm]]
    if (is.null(d)) base_unit(nm) else d
  }

  for (a in ssm) {
    ex <- a$expr
    syms <- unique(.ast_leaf_syms(ex))
    total <- numeric(ncoord); names(total) <- coords
    for (s in syms) {
      ds <- if (s %in% endo)
        ast_differentiate(ex, s, 0L, wrt = "variable")
      else
        ast_differentiate(ex, s, 0L, wrt = "parameter")
      pv <- tryCatch(ast_eval(ds, vv, params), error = function(e) NA_real_)
      if (is.na(pv) || !is.finite(pv)) return(NULL)   # non-differentiable / eval fail
      if (pv != 0) total <- total + pv * get_d(s)
    }
    dmap[[a$name]] <- total
  }

  ## Assemble the requested rows for the computed parameters.
  dg_dy <- matrix(0, length(computed), length(endo),
                  dimnames = list(computed, endo))
  dg_dtheta <- matrix(0, length(computed), length(pars),
                      dimnames = list(computed, pars))
  for (pc in computed) {
    d <- dmap[[pc]]
    dg_dy[pc, ] <- d[endo]
    ## free-param columns from the chain; computed-param columns are 0 (a
    ## computed parameter is not an independent coordinate).
    dg_dtheta[pc, free_params] <- d[free_params]
  }
  list(dg_dy = dg_dy, dg_dtheta = dg_dtheta, computed = computed)
}


#' Second total partial derivatives of steady_state_model-computed parameters.
#'
#' The order-2 companion of \code{.ssm_param_chain_derivs()}. Walks the SSM
#' assignment list IN ORDER carrying, for every assigned name, BOTH its first
#' total derivative (\code{d1[name]}, a length-\code{ncoord} vector) and its
#' second total derivative (\code{d2[name]}, a \code{ncoord x ncoord} matrix)
#' over the coordinate set \code{c(endo_vars, free_params)}. Each statement
#' \eqn{v := g(\dots)} applies the chain + product rule numerically at
#' \eqn{(\bar y, \theta)}:
#' \deqn{d^2 v = \sum_s \frac{\partial g}{\partial s} d^2 s
#'              + \sum_{s,t} \frac{\partial^2 g}{\partial s\,\partial t}\,(d s)(d t)^\top}
#' where \eqn{ds}, \eqn{d^2 s} are the accumulated first/second derivatives of an
#' earlier-assigned name (or the unit vector / zero for a still-free coordinate)
#' and the partials are symbolic \code{ast_differentiate} of the RHS evaluated at
#' the point.
#'
#' Returns, for the SSM-COMPUTED PARAMETERS only, the \code{ncoord x ncoord}
#' second-derivative block (\eqn{\partial^2 g_{p_c}} over coords) the order-2
#' augmented-Jacobian d2ys chain needs (its \eqn{yy}/\eqn{y\theta}/\eqn{\theta\theta}
#' sub-blocks). NULL if the SSM assigns no parameter or has a non-differentiable RHS.
#'
#' @param compiled dynhr_compiled.
#' @param ys       Named steady state (declaration order).
#' @param params   Named parameter vector (incl. the SSM-computed values).
#' @return list with \code{computed} (names), \code{coords}, \code{free}, and
#'   \code{d2} (named list per computed parameter; each a \code{ncoord x ncoord}
#'   symmetric matrix with \code{coords} dimnames), or NULL.
#' @noRd
.ssm_param_chain_derivs2 <- function(compiled, ys, params) {
  model <- compiled$model
  ssm   <- model$steady_state_model
  if (is.null(ssm) || length(ssm) == 0L) return(NULL)
  endo  <- model$var_names
  pars  <- model$param_names

  assigned <- vapply(ssm, function(a) a$name %||% NA_character_, character(1))
  computed <- intersect(assigned, pars)
  if (length(computed) == 0L) return(NULL)

  free   <- setdiff(pars, computed)
  coords <- c(endo, free)
  nc     <- length(coords)
  vv     <- setNames(as.numeric(ys[endo]), paste0(endo, "__0"))
  wrt    <- function(s) if (s %in% endo) "variable" else "parameter"
  unit   <- function(nm) { v <- numeric(nc); names(v) <- coords
                           if (nm %in% coords) v[[nm]] <- 1; v }

  d1 <- vector("list", 0L); d2 <- vector("list", 0L)
  g1 <- function(nm) { d <- d1[[nm]]; if (is.null(d)) unit(nm) else d }
  g2 <- function(nm) { d <- d2[[nm]]
                       if (is.null(d)) matrix(0, nc, nc, dimnames = list(coords, coords)) else d }

  for (a in ssm) {
    ex <- a$expr; syms <- unique(.ast_leaf_syms(ex))
    t1 <- setNames(numeric(nc), coords)
    t2 <- matrix(0, nc, nc, dimnames = list(coords, coords))
    ds_cache <- vector("list", 0L)
    for (s in syms) {
      ds <- ast_differentiate(ex, s, 0L, wrt = wrt(s))
      ds_cache[[s]] <- ds
      pv <- tryCatch(ast_eval(ds, vv, params), error = function(e) NA_real_)
      if (is.na(pv) || !is.finite(pv)) return(NULL)
      if (pv != 0) { t1 <- t1 + pv * g1(s); t2 <- t2 + pv * g2(s) }
    }
    for (s in syms) for (t in syms) {
      dst  <- ast_differentiate(ds_cache[[s]], t, 0L, wrt = wrt(t))
      pvst <- tryCatch(ast_eval(dst, vv, params), error = function(e) NA_real_)
      if (is.na(pvst) || !is.finite(pvst)) return(NULL)
      if (pvst != 0) t2 <- t2 + pvst * outer(g1(s), g1(t))
    }
    d1[[a$name]] <- t1; d2[[a$name]] <- t2
  }

  list(computed = computed, coords = coords, free = free, d2 = d2[computed])
}


#' Precompute the equation->declaration row permutation and column layout.
#'
#' Mirrors the row reorder + timing-column partition of
#' \code{extract_system_matrices()} so an analytic total-derivative Jacobian
#' (built in compiled/equation order) is partitioned into df_minus/df_zero/
#' df_plus/df_exo exactly as the finite-difference path would be.
#'
#' @param compiled dynhr_compiled.
#' @return list(perm, cache, col_var_idx, ...) -- structure-only, value-free.
#' @noRd
.dsys_layout <- function(compiled) {
  model <- compiled$model
  dyn   <- compiled$dynamic
  endo  <- model$var_names
  n_eq  <- dyn$n_eq

  ## eq_to_decl: same LHS-variable matching as extract_system_matrices().
  perm <- .eq_to_decl_perm(model, n_eq)

  ## Reuse the cached timing-column maps (minus/zero/plus/exo) -- identical to
  ## the partition the fast/slow extractors apply after the row reorder.
  cache <- cache_system_structure(compiled)

  ## Map each dynamic column -> endogenous index (for the dys lookup); NA for
  ## the shock columns (whose steady-state value is a parameter-free 0).
  dcm <- dyn$dyn_col_map
  col_var_idx <- match(dcm$name, endo)   # length total_cols; NA for exo cols

  list(perm = perm, cache = cache, col_var_idx = col_var_idx,
       n_eq = n_eq, n_endo = length(endo), n_exo = dyn$n_exo,
       total_cols = dyn$total_cols)
}


#' Equation->declaration-order row permutation.
#'
#' Replicates the LHS-variable matching used inside
#' \code{extract_system_matrices()} (each equation's LHS endogenous variable
#' fixes its row position; unmatched equations fill the remaining slots in
#' order). Returns \code{order(eq_to_decl)} -- the permutation that reorders
#' Jacobian rows from equation order to declaration order, or
#' \code{seq_len(n_eq)} when the orders already coincide.
#'
#' @noRd
.eq_to_decl_perm <- function(model, n_eq) {
  endo <- model$var_names
  ## Use the shared package-level LHS-variable helper (.lhs_endo_var in
  ## solve-extract-system.R) -- fixes the dead "uniop" branch and eliminates
  ## the duplicate copy of this logic.
  used <- logical(length(endo)); names(used) <- endo
  eq_to_decl <- integer(n_eq)
  for (i in seq_len(n_eq)) {
    v <- .lhs_endo_var(model$equations[[i]]$lhs)
    if (!is.null(v) && nzchar(v) && v %in% endo && !used[[v]]) {
      used[[v]] <- TRUE
      eq_to_decl[i] <- match(v, endo)
    }
  }
  unassigned <- which(!used)
  unmapped   <- which(eq_to_decl == 0)
  for (k in seq_along(unmapped))
    if (k <= length(unassigned)) eq_to_decl[unmapped[k]] <- unassigned[k]
  if (all(eq_to_decl > 0) && !identical(eq_to_decl, seq_len(n_eq)))
    order(eq_to_decl)
  else
    seq_len(n_eq)
}


#' Analytic steady-state sensitivity dȳ/dθ for all parameters.
#'
#' Solves \eqn{J_{static}\,(d\bar y/d\theta) = -\partial F_{static}/\partial\theta}
#' once (single factorization) for every parameter.
#'
#' @param compiled dynhr_compiled.
#' @param ys       Named numeric steady state (declaration order).
#' @param params   Named numeric parameter vector.
#' @return n_endo x n_params matrix of dys (columns named by parameter), or
#'   NULL if the static Jacobian is unavailable / singular.
#' @noRd
.analytic_dys <- function(compiled, ys, params) {
  sta <- compiled$static
  if (is.null(sta$param_resid_fn)) return(NULL)
  exo <- compiled$model$varexo_names
  xz  <- setNames(rep(0, length(exo)), exo)

  J_static <- sta$jacobian_fn(ys, xz, params)          # n_eq x n_endo
  dF_dtheta <- sta$param_resid_fn(ys, xz, params)      # n_eq x n_params
  if (nrow(J_static) != ncol(J_static)) return(NULL)
  if (any(!is.finite(J_static)) || any(!is.finite(dF_dtheta))) return(NULL)

  ## SSM-derived-parameter chain (Tier 12 #2): when the steady_state_model
  ## computes a parameter p_c = g(ybar, theta), perturbing a free theta also
  ## moves p_c -- a dependence param_resid_fn (which holds every parameter,
  ## p_c included, FIXED) omits. Total-differentiating F(ybar, theta, p_c)=0
  ## through the chain gives the augmented system
  ##   J_aug = J + sum_pc (dF/dp_c)(dg_pc/dy),
  ##   rhs_k = dF/dtheta_k + sum_pc (dF/dp_c)(dg_pc/dtheta_k)
  ## (dg_pc/dy, dg_pc/dtheta from .ssm_param_chain_derivs, verified exact). For
  ## caldara_rp J_aug stays singular (the chain lives off the singular
  ## s-direction) so qr.solve below returns NULL and the caller uses FD.
  if (.ssm_assigns_param(compiled$model)) {
    ch <- .ssm_param_chain_derivs(compiled, ys, params)
    if (is.null(ch)) return(NULL)
    endo <- compiled$model$var_names
    pars <- compiled$model$param_names
    free <- setdiff(pars, ch$computed)
    for (pc in ch$computed) {
      dFpc <- as.numeric(dF_dtheta[, match(pc, pars)])
      J_static <- J_static + outer(dFpc, as.numeric(ch$dg_dy[pc, endo]))
      for (k in free) {
        kc <- match(k, pars)
        dF_dtheta[, kc] <- dF_dtheta[, kc] + dFpc * ch$dg_dtheta[pc, k]
      }
    }
    if (any(!is.finite(J_static)) || any(!is.finite(dF_dtheta))) return(NULL)
  }

  out <- tryCatch(
    qr.solve(J_static, -dF_dtheta),
    error = function(e) NULL)
  if (is.null(out)) return(NULL)
  if (any(!is.finite(out))) return(NULL)
  out <- matrix(out, nrow = ncol(J_static), ncol = ncol(dF_dtheta))
  rownames(out) <- compiled$static$endo_names
  ## Columns are indexed by the compile-time parameter order (the full
  ## model$param_names), which may be a superset of names(params) -- e.g. a
  ## shock standard deviation declared as a parameter but set in the shocks
  ## block. Consumers select the parameters they need by name.
  colnames(out) <- compiled$model$param_names
  out
}


#' Analytic total parameter-derivatives of the dynamic-Jacobian primitives.
#'
#' For every parameter, returns the total first derivatives
#' \code{df_plus/df_zero/df_minus/df_exo} of the dynamic-Jacobian blocks --
#' the analytic replacement for the central-FD primitives in
#' \code{.solution_deriv_one()}.
#'
#' @param compiled dynhr_compiled.
#' @param ys       Named steady state.
#' @param params   Named parameter vector.
#' @param dys      n_endo x n_params matrix from \code{.analytic_dys()}.
#' @param layout   Output of \code{.dsys_layout()} (optional; built if NULL).
#' @return named-by-parameter list, each list(df_plus, df_zero, df_minus,
#'   df_exo); or NULL on failure (caller falls back to FD).
#' @noRd
.analytic_dprimitives <- function(compiled, ys, params, dys, layout = NULL) {
  if (is.null(dys)) return(NULL)
  dyn <- compiled$dynamic
  if (is.null(dyn$param_jacobian_fn)) return(NULL)
  if (is.null(layout)) layout <- .dsys_layout(compiled)

  n_eq       <- layout$n_eq
  n_endo     <- layout$n_endo
  n_exo      <- layout$n_exo
  total_cols <- layout$total_cols
  np         <- ncol(dys)
  perm       <- layout$perm
  cache      <- layout$cache
  col_var_idx <- layout$col_var_idx

  ## Dynamic point dy[c] = ȳ[var(c)] (shock columns 0), as extract builds it.
  dy <- numeric(total_cols)
  endo_cols <- which(!is.na(col_var_idx))
  dy[endo_cols] <- ys[compiled$model$var_names[col_var_idx[endo_cols]]]
  names(dy) <- cache$dy_keys

  ## Bail to FD if the base Jacobian is non-finite (0/0 at ss repaired by FD).
  Jbase <- dyn$jacobian_fn(dy, params, ys)
  if (any(!is.finite(Jbase))) return(NULL)

  ## Explicit channel: P[i,c,k] = ∂²F_i/(∂w_c ∂θ_k), flattened to [n_eq*tc, np].
  P <- dyn$param_jacobian_fn(dy, params, ys)
  if (any(!is.finite(P))) return(NULL)
  dJ_all <- matrix(P, nrow = n_eq * total_cols, ncol = np)

  ## SSM-computed-parameter chain (Tier 12 #2): capture the PURE param-Jacobian
  ## columns ∂²F/(∂w ∂p_c) for the SSM-derived parameters BEFORE the steady-state
  ## Hessian scatter-add mutates dJ_all (channel 3 is applied after that block).
  ssm_ch <- NULL; P_computed <- NULL
  if (.ssm_assigns_param(compiled$model)) {
    ssm_ch <- .ssm_param_chain_derivs(compiled, ys, params)
    if (is.null(ssm_ch)) return(NULL)
    P_computed <- dJ_all[, match(ssm_ch$computed, colnames(dys)), drop = FALSE]
    colnames(P_computed) <- ssm_ch$computed
  }

  ## Steady-state chain channel: contract the symbolic model Hessian with the
  ## column-broadcast steady-state sensitivity V[c,k] = dys[var(c),k] (0 for
  ## shock columns). hessian2_fn returns sparse upper-triangle (eq,c1,c2) vals.
  hv <- dyn$hessian2_fn(dy, params, ys)
  if (length(hv) > 0L && any(!is.finite(hv))) return(NULL)
  if (length(hv) > 0L) {
    trip <- dyn$hess2_triplets
    eqs  <- vapply(trip, `[[`, integer(1), "eq")
    c1s  <- vapply(trip, `[[`, integer(1), "col1")
    c2s  <- vapply(trip, `[[`, integer(1), "col2")

    V <- matrix(0, total_cols, np)
    V[endo_cols, ] <- dys[col_var_idx[endo_cols], , drop = FALSE]

    ## Contribution of (eq,c1,c2,val): val*V[c2] -> (eq,c1); val*V[c1] -> (eq,c2).
    lin1 <- eqs + (c1s - 1L) * n_eq
    A1   <- hv * V[c2s, , drop = FALSE]          # [n_hess x np]
    dJ_all <- .scatter_add(dJ_all, lin1, A1)

    off  <- which(c1s != c2s)
    if (length(off)) {
      lin2 <- eqs[off] + (c2s[off] - 1L) * n_eq
      A2   <- hv[off] * V[c1s[off], , drop = FALSE]
      dJ_all <- .scatter_add(dJ_all, lin2, A2)
    }
  }

  ## Channel 3 (Tier 12 #2): SSM-computed-parameter chain. p_c moves with theta,
  ## so each free column k gains (∂²F/∂w∂p_c)·(dp_c/dtheta_k), where the cross
  ## term is the captured pure param-Jacobian p_c column and the TOTAL
  ## steady-state derivative dp_c/dtheta_k = dg_dtheta[pc,k] + Σ_y dg_dy[pc,y]·dys[y,k].
  if (!is.null(ssm_ch)) {
    endo <- compiled$model$var_names
    pars <- compiled$model$param_names
    free <- setdiff(pars, ssm_ch$computed)
    for (pc in ssm_ch$computed) {
      pcol <- P_computed[, pc]
      for (k in free) {
        kc  <- match(k, colnames(dys))
        dpc <- ssm_ch$dg_dtheta[pc, k] + sum(ssm_ch$dg_dy[pc, endo] * dys[endo, kc])
        if (dpc != 0) dJ_all[, kc] <- dJ_all[, kc] + pcol * dpc
      }
    }
  }

  if (any(!is.finite(dJ_all))) return(NULL)

  ## Partition each parameter's total-derivative Jacobian exactly as extract.
  out <- vector("list", np)
  names(out) <- colnames(dys)
  for (k in seq_len(np))
    out[[k]] <- .partition_dJ(matrix(dJ_all[, k], n_eq, total_cols), layout, compiled)
  out
}


#' Partition a [n_eq x total_cols] dynamic-Jacobian-derivative matrix (in
#' COMPILED equation order) into the df_minus/df_zero/df_plus/df_exo blocks,
#' applying the eq->declaration row permutation and the timing-column maps --
#' exactly as \code{extract_system_matrices()} partitions the level Jacobian.
#'
#' Shared by the first-order (\code{.analytic_dprimitives}) and second-order
#' (\code{solution_derivatives_2}) layers so the reorder/partition logic lives
#' in ONE place.
#'
#' @param dJ      n_eq x total_cols, COMPILED equation order (pre-permutation).
#' @param layout  output of \code{.dsys_layout()}.
#' @param compiled dynhr_compiled (for column names).
#' @return list(df_plus, df_zero, df_minus, df_exo).
#' @noRd
.partition_dJ <- function(dJ, layout, compiled) {
  cache <- layout$cache; perm <- layout$perm
  n_eq <- layout$n_eq; n_endo <- layout$n_endo
  n_exo <- layout$n_exo; total_cols <- layout$total_cols
  mi <- cache$minus_idx; zi <- cache$zero_idx
  pri <- cache$plus_idx; ej <- cache$exo_jcols
  evalid <- which(ej <= total_cols)

  dJ <- dJ[perm, , drop = FALSE]                 # equation -> declaration order
  df_minus <- matrix(0, n_eq, n_endo)
  df_zero  <- matrix(0, n_eq, n_endo)
  df_plus  <- matrix(0, n_eq, n_endo)
  df_exo   <- matrix(0, n_eq, n_exo)
  if (nrow(mi))  df_minus[, mi[, 1]] <- dJ[, mi[, 2]]
  if (nrow(zi))  df_zero[,  zi[, 1]] <- dJ[, zi[, 2]]
  if (nrow(pri)) df_plus[,  pri[, 1]] <- dJ[, pri[, 2]]
  if (length(evalid)) df_exo[, evalid] <- dJ[, ej[evalid]]
  colnames(df_minus) <- compiled$model$var_names
  colnames(df_zero)  <- compiled$model$var_names
  colnames(df_plus)  <- compiled$model$var_names
  colnames(df_exo)   <- compiled$model$varexo_names
  list(df_plus = df_plus, df_zero = df_zero, df_minus = df_minus, df_exo = df_exo)
}


#' Scatter-add rows of A into a flat-indexed matrix, accumulating duplicates.
#'
#' Returns \code{M} with \code{A} added at rows \code{lin}, accumulating over
#' repeated \code{lin} entries (base R's matrix-index assignment keeps only
#' the last write, so duplicate targets are pre-summed with \code{rowsum}).
#'
#' @param M   numeric matrix [N x np] (linear row index space).
#' @param lin integer vector of target rows (length == nrow(A)).
#' @param A   numeric matrix [length(lin) x np] of additions.
#' @return the updated matrix M.
#' @noRd
.scatter_add <- function(M, lin, A) {
  agg  <- rowsum(A, group = lin, reorder = FALSE)
  rows <- as.integer(rownames(agg))
  M[rows, ] <- M[rows, ] + agg
  M
}


#' Densify a sparse static second-order triplet list into an [n_eq, d2, d3] array.
#' @noRd
.densify_static <- function(trip, vals, n_eq, d2, d3, ia, ib, sym) {
  A <- array(0, dim = c(n_eq, d2, d3))
  for (k in seq_along(trip)) {
    t <- trip[[k]]
    A[t$eq, t[[ia]], t[[ib]]] <- vals[k]
    if (sym && t[[ia]] != t[[ib]]) A[t$eq, t[[ib]], t[[ia]]] <- vals[k]
  }
  A
}


#' Analytic second-order steady-state sensitivity d2ys = d2(ybar)/dtheta_a dtheta_b.
#'
#' Differentiates the implicit first-order relation \eqn{J_s\,dys_a + g_a = 0}
#' (with \eqn{J_s=\partial F_s/\partial y}, \eqn{g_a=\partial F_s/\partial\theta_a})
#' a second time w.r.t. \eqn{\theta_b}:
#' \deqn{J_s\,d2ys_{ab} = -[(\partial^2F_s/\partial y\partial\theta_b + H_s\cdot dys_b)\cdot dys_a
#'        + \partial^2F_s/\partial\theta_a\partial\theta_b
#'        + (\partial^2F_s/\partial y\partial\theta_a)\cdot dys_b]}
#' reusing a single factorization of \eqn{J_s} for every (a,b) pair. Symmetric
#' in (a,b). Requires the static second-order codegen
#' (\code{compile_model(param_deriv="second")}); returns NULL otherwise so the
#' caller falls back to finite differences.
#'
#' @param compiled dynhr_compiled.
#' @param ys       Named steady state.
#' @param params   Named parameter vector.
#' @param dys      n_endo x n_params matrix from \code{.analytic_dys()}.
#' @return [n_endo x n_params x n_params] symmetric array (params indexed by
#'   \code{model$param_names}), or NULL.
#' @noRd
.analytic_d2ys <- function(compiled, ys, params, dys) {
  if (is.null(dys)) return(NULL)
  sta <- compiled$static
  if (!isTRUE(sta$static_param2_built)) return(NULL)
  if (is.null(sta$static_hess2_fn) || is.null(sta$static_param_jac_fn) ||
      is.null(sta$static_param2_fn)) return(NULL)

  exo <- compiled$model$varexo_names
  xz  <- setNames(rep(0, length(exo)), exo)
  n_eq   <- sta$n_eq
  n_endo <- sta$n_endo
  np     <- ncol(dys)

  J_static <- sta$jacobian_fn(ys, xz, params)
  if (nrow(J_static) != ncol(J_static)) return(NULL)
  if (any(!is.finite(J_static))) return(NULL)

  ## Dense static tensors (params indexed over the full model param set == np).
  Hs <- .densify_static(sta$static_hess2_triplets,    sta$static_hess2_fn(ys, xz, params),
                        n_eq, n_endo, n_endo, "p", "q",     TRUE)
  Bs <- .densify_static(sta$static_param_jac_triplets, sta$static_param_jac_fn(ys, xz, params),
                        n_eq, n_endo, np,     "p", "param", FALSE)
  Cs <- .densify_static(sta$static_param2_triplets,    sta$static_param2_fn(ys, xz, params),
                        n_eq, np,     np,     "pa", "pb",   TRUE)
  if (any(!is.finite(Hs)) || any(!is.finite(Bs)) || any(!is.finite(Cs))) return(NULL)

  ## Hdys[[b]][i,p] = sum_q Hs[i,p,q] * dys[q,b].
  ## Flatten Hs to [(i,p), q] (i fastest) so one matrix multiply does the contraction.
  Hmat <- matrix(Hs, n_eq * n_endo, n_endo)
  Hdys <- lapply(seq_len(np), function(b) matrix(Hmat %*% dys[, b], n_eq, n_endo))

  pars <- compiled$model$param_names
  out  <- array(0, dim = c(n_endo, np, np))

  if (.ssm_assigns_param(compiled$model)) {
    ## ---- SSM-computed-parameter branch (Tier 12 #2, second order) ----------
    ## The order-2 analog of .analytic_dys' augmentation + .analytic_dprimitives'
    ## channel 3. The implicit-function operator is the SAME augmented Jacobian
    ## J_aug = J + sum_pc (dF/dp_c)(dg_pc/dy); the forcing uses TOTAL parameter
    ## directions w_a = e_a + sum_pc (dp_c/dtheta_a) e_pc in the Bs/Cs
    ## contractions and adds the term (dF/dp_c)*(d2 p_c / dtheta_a dtheta_b),
    ## with the SSM second chain d2 p_c (minus its dg_dy*d2ys feedback, which the
    ## augmented operator already carries). Validated end-to-end on rbc2shock_ssm
    ## against the full re-solve (d/dtheta_b of the analytic dys, ~2.5e-9 rel).
    ch1 <- .ssm_param_chain_derivs(compiled, ys, params)
    ch2 <- .ssm_param_chain_derivs2(compiled, ys, params)
    if (is.null(ch1) || is.null(ch2)) return(NULL)
    endo <- compiled$model$var_names
    comp <- ch1$computed
    free <- setdiff(pars, comp)
    dF <- sta$param_resid_fn(ys, xz, params); colnames(dF) <- pars

    J_aug <- J_static
    for (pc in comp)
      J_aug <- J_aug + outer(as.numeric(dF[, pc]), as.numeric(ch1$dg_dy[pc, endo]))
    Jqr <- tryCatch(qr(J_aug), error = function(e) NULL)
    if (is.null(Jqr)) return(NULL)

    ## Total first derivative dp_c/dtheta_a for each computed pc, free a.
    dpmat <- matrix(0, length(comp), np, dimnames = list(comp, pars))
    for (pc in comp) for (a in free)
      dpmat[pc, a] <- ch1$dg_dtheta[pc, a] + sum(ch1$dg_dy[pc, endo] * dys[endo, a])

    ## Total parameter-direction matrix W[, a] = e_a + sum_pc dp_c_a e_pc.
    W <- matrix(0, np, np, dimnames = list(pars, pars))
    for (a in free) { W[a, a] <- 1; for (pc in comp) W[pc, a] <- W[pc, a] + dpmat[pc, a] }

    Bsmat <- matrix(Bs, n_eq * n_endo, np)        # [(i,p), k]
    Csmat <- matrix(Cs, n_eq * np, np)            # [(i,k), l]
    BsW <- function(a) matrix(Bsmat %*% W[, a], n_eq, n_endo)      # sum_k Bs[,,k] W[k,a]

    ## d2 p_c lower part (excludes the dg_dy*d2ys feedback carried by J_aug):
    ## g_thth[a,b] + g_yth[,a].dys_b + g_yth[,b].dys_a + dys_a' g_yy dys_b.
    d2pc_lower <- function(pc, a, b) {
      M <- ch2$d2[[pc]]                          # ncoord x ncoord over c(endo, free)
      M[a, b] + sum(M[endo, a] * dys[endo, b]) + sum(M[endo, b] * dys[endo, a]) +
        sum(outer(dys[endo, a], dys[endo, b]) * M[endo, endo])
    }

    for (ia in seq_along(free)) for (ib in ia:length(free)) {
      a <- free[ia]; b <- free[ib]
      ka <- match(a, pars); kb <- match(b, pars)
      CsW_b <- matrix(Csmat %*% W[, b], n_eq, np)               # sum_l Cs[,,l] W[l,b]
      rhs <- (BsW(b) + Hdys[[kb]]) %*% dys[, a] +
             as.numeric(CsW_b %*% W[, a]) +
             BsW(a) %*% dys[, b]
      for (pc in comp) rhs <- rhs + as.numeric(dF[, pc]) * d2pc_lower(pc, a, b)
      d2 <- tryCatch(qr.solve(Jqr, -rhs), error = function(e) NULL)
      if (is.null(d2) || any(!is.finite(d2))) return(NULL)
      out[, ka, kb] <- as.numeric(d2)
      if (ka != kb) out[, kb, ka] <- as.numeric(d2)
    }
    dimnames(out) <- list(sta$endo_names, pars, pars)
    return(out)
  }

  Jqr <- tryCatch(qr(J_static), error = function(e) NULL)
  if (is.null(Jqr)) return(NULL)
  for (a in seq_len(np)) {
    for (b in a:np) {
      ## RHS_i = sum_p (Bs[i,p,b] + Hdys_b[i,p]) dys[p,a] + Cs[i,a,b] + sum_p Bs[i,p,a] dys[p,b]
      rhs <- (Bs[, , b] + Hdys[[b]]) %*% dys[, a] +
             Cs[, a, b] +
             Bs[, , a] %*% dys[, b]
      d2 <- tryCatch(qr.solve(Jqr, -rhs), error = function(e) NULL)
      if (is.null(d2) || any(!is.finite(d2))) return(NULL)
      out[, a, b] <- as.numeric(d2)
      if (a != b) out[, b, a] <- as.numeric(d2)
    }
  }
  dimnames(out) <- list(sta$endo_names, compiled$model$param_names,
                        compiled$model$param_names)
  out
}
