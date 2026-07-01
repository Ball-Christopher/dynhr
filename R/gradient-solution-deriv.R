## R/gradient-solution-deriv.R
## --------------------------------------------------------------------------
## SOLUTION-DERIVATIVE layer for implicit-differentiation likelihood gradients
## (Childers, Fernandez-Villaverde, Perla, Rackauckas & Wu 2022, NBER w30573).
##
## Computes analytic derivatives of the first-order decision rule (ghx, ghu,
## ys) with respect to structural parameters via implicit differentiation of
## the perturbation fixed-point equation, instead of finite-differencing the
## full QZ solve (which is noisy and expensive per parameter).
##
## THE MATH
## --------
## Let S be the n_state x n_endo selection matrix (identity rows at
## dr$state_idx), G = ghx (n_endo x n_state), H = ghu (n_endo x n_exo). The
## first-order fixed point is
##
##   R1(G, theta) = f_plus G (S G) + f_zero G + f_minus_S = 0   (n_endo x n_state)
##
## with f_minus_S = f_minus[, state_idx], and
##
##   A H + f_u = 0,   A = f_plus G S + f_zero   =>   H = -A^{-1} f_u
##
## Implicit differentiation wrt theta_j (TOTAL derivatives of the f-matrices,
## including through the steady state):
##
##   f_plus dG (S G) + A dG = -[df_plus G (S G) + df_zero G + df_minus_S]
##
## In vec() form:
##
##   [ (S G)' (x) f_plus + I_{n_state} (x) A ] vec(dG) = -vec(RHS_j)
##
## where (x) is the Kronecker product. The coefficient matrix M does not
## depend on j, so it is built and factorized ONCE and reused for every
## parameter (the Childers et al. efficiency point).
##
##   dH = -A^{-1} ( dA H + df_u ),  dA = df_plus G S + f_plus dG S + df_zero
##
## A's factorization (the same QR used to solve A H = -f_u) is reused too.
##
## TOTAL df-matrices (df_plus, df_zero, df_minus, df_u, dys) are obtained by
## central finite differences of the SMOOTH primitive -- the dynamic Jacobian
## evaluated at the (perturbed steady state, perturbed params) -- NOT of the
## QZ solve.
## --------------------------------------------------------------------------


#' Implicit-differentiation derivatives of the first-order decision rule
#'
#' Computes dG = d(ghx)/d(theta_j), dH = d(ghu)/d(theta_j), and
#' dys = d(ys)/d(theta_j) for each requested parameter, via implicit
#' differentiation of the perturbation fixed-point equations (Childers,
#' Fernandez-Villaverde, Perla, Rackauckas & Wu 2022). Also assembles the
#' corresponding state-space derivative blocks (dTT, dRR, dZZ, dDD, dd)
#' matching the conventions of \code{kalman_filter()}.
#'
#' @param model      dynhr_mod
#' @param compiled   dynhr_compiled
#' @param dr         DecisionRules from \code{solve_perturbation} (order 1)
#' @param params     Named numeric parameter vector (the base point)
#' @param param_names Character vector of parameter names to differentiate wrt
#' @param obs_vars   Character vector of observed variable names
#' @param h_rel      Relative step size for central finite differences of the
#'   smooth f-matrix / steady-state primitives (default 1e-6)
#' @return A list with:
#'   \item{derivs}{named list (by param_names) of per-parameter results, each
#'     with dG, dH, dys, dTT, dRR, dZZ, dDD, dd, ok, and (if !ok) a message}
#'   \item{base_residual}{list with R1 and AH residual max-abs at the base
#'     solution}
#'   \item{reused_factorization}{TRUE -- M and A are each factorized once and
#'     reused across all parameters}
#' @noRd

## ---------------------------------------------------------------------------
## Dedicated k=1 generalized-Sylvester solver:  A X + fp X hx = RHS.
##
## Real-Schur (Bartels-Stewart) form, FACTORED ONCE and reused across all RHS:
## with hx = U T U' (real Schur, T quasi-triangular) and Y = X U, C = RHS U the
## equation becomes  A Y + fp Y T = C, solved column-block by column-block over
## T's 1x1 / 2x2 blocks. Each block coefficient -- (A + lambda*fp) for a 1x1
## block, a 2n x 2n system for a 2x2 complex pair -- is QR-factored ONCE in
## _factor and back-substituted per RHS in _solve. O(ns*n^3) setup + O(ns*n^2)
## per RHS, all in REAL arithmetic (no complex eigen-grid / .apply_kronk), so it
## beats the general-k .solve_kron_compact for this k=1 case. _factor self-tests
## on a fixed RHS and returns NULL (=> caller falls back to .solve_kron_compact)
## if QZ fails or the realized residual is large, so the result is never wrong.
.gen_sylvester_k1_factor <- function(A, fp, hx, tol = 1e-9) {
  n <- nrow(A); ns <- nrow(hx)
  if (ns == 0L) return(list(n = n, ns = 0L, empty = TRUE))
  if (!requireNamespace("QZ", quietly = TRUE)) return(NULL)
  sch <- tryCatch(QZ::qz.dgees(hx), error = function(e) NULL)
  if (is.null(sch) || !identical(as.integer(sch$INFO), 0L)) return(NULL)
  U <- sch$Q; Tm <- sch$T
  blocks <- list(); k <- 1L
  ok <- TRUE
  while (k <= ns) {
    is2 <- (k < ns) && (abs(Tm[k + 1L, k]) > 0)
    blk <- tryCatch({
      if (!is2) {
        list(type = 1L, k = k, fac = qr(A + Tm[k, k] * fp))
      } else {
        a <- Tm[k, k]; b <- Tm[k, k + 1L]; cc <- Tm[k + 1L, k]; d <- Tm[k + 1L, k + 1L]
        Big <- rbind(cbind(A + a * fp, cc * fp),
                     cbind(b * fp,     A + d * fp))
        list(type = 2L, k = k, fac = qr(Big))
      }
    }, error = function(e) NULL)
    if (is.null(blk)) { ok <- FALSE; break }
    blocks[[length(blocks) + 1L]] <- blk
    k <- k + if (is2) 2L else 1L
  }
  if (!ok) return(NULL)
  fac <- list(n = n, ns = ns, U = U, Tm = Tm, fp = fp, blocks = blocks,
              empty = FALSE)
  ## Deterministic self-test: solve a fixed RHS and check the Sylvester residual.
  Rtest <- matrix(1, n, ns)
  Xt <- tryCatch(.gen_sylvester_k1_solve(fac, Rtest), error = function(e) NULL)
  if (is.null(Xt)) return(NULL)
  resid <- max(abs(A %*% Xt + fp %*% Xt %*% hx - Rtest))
  if (!is.finite(resid) || resid > tol * max(1, max(abs(Rtest)))) return(NULL)
  fac
}

.gen_sylvester_k1_solve <- function(fac, RHS) {
  if (isTRUE(fac$empty)) return(matrix(0, fac$n, 0L))
  n <- fac$n; Tm <- fac$Tm; fp <- fac$fp
  C <- RHS %*% fac$U
  Y <- matrix(0, n, fac$ns)
  for (blk in fac$blocks) {
    k <- blk$k
    if (blk$type == 1L) {
      rhs <- C[, k]
      if (k > 1L)
        rhs <- rhs - fp %*% (Y[, seq_len(k - 1L), drop = FALSE] %*%
                               Tm[seq_len(k - 1L), k])
      Y[, k] <- qr.solve(blk$fac, rhs)
    } else {
      k1 <- k; k2 <- k + 1L
      rhs1 <- C[, k1]; rhs2 <- C[, k2]
      if (k1 > 1L) {
        prev <- Y[, seq_len(k1 - 1L), drop = FALSE]
        rhs1 <- rhs1 - fp %*% (prev %*% Tm[seq_len(k1 - 1L), k1])
        rhs2 <- rhs2 - fp %*% (prev %*% Tm[seq_len(k1 - 1L), k2])
      }
      sol <- qr.solve(blk$fac, c(rhs1, rhs2))
      Y[, k1] <- sol[seq_len(n)]
      Y[, k2] <- sol[n + seq_len(n)]
    }
  }
  Y %*% t(fac$U)
}

solution_derivatives <- function(model, compiled, dr, params, param_names,
                                  obs_vars, h_rel = 1e-6,
                                  use_analytic = NULL) {

  if (!all(param_names %in% names(params))) {
    missing_p <- param_names[!param_names %in% names(params)]
    stop(sprintf("solution_derivatives: parameter(s) not found in `params`: %s",
                  paste(missing_p, collapse = ", ")))
  }

  endo    <- dr$endo_names
  exo     <- dr$exo_names
  n_endo  <- length(endo)
  n_exo   <- length(exo)
  state_idx <- dr$state_idx
  n_state <- length(state_idx)

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("solution_derivatives: observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  G  <- dr$ghx                     # n_endo x n_state
  H  <- dr$ghu                     # n_endo x n_exo
  ys <- dr$ys                      # named n_endo

  ## Selection matrix S: n_state x n_endo, identity rows at state_idx.
  S <- matrix(0, nrow = n_state, ncol = n_endo)
  if (n_state > 0) S[cbind(seq_len(n_state), state_idx)] <- 1

  ## ------------------------------------------------------------------
  ## Base system matrices (f_plus, f_zero, f_minus, f_exo) at (ys, params)
  ## ------------------------------------------------------------------
  sys0 <- extract_system_matrices(compiled, ys, params)
  f_plus  <- sys0$f_plus    # n_endo x n_endo
  f_zero  <- sys0$f_zero    # n_endo x n_endo
  f_minus <- sys0$f_minus   # n_endo x n_endo
  f_u     <- sys0$f_exo     # n_endo x n_exo

  f_minus_S <- if (n_state > 0) f_minus[, state_idx, drop = FALSE]
               else matrix(0, nrow = n_endo, ncol = 0)

  SG <- if (n_state > 0) S %*% G else matrix(0, nrow = 0, ncol = n_state)  # n_state x n_state

  ## A = f_plus G S + f_zero  (n_endo x n_endo)
  GS <- if (n_state > 0) G %*% S else matrix(0, nrow = n_endo, ncol = n_endo)
  A  <- f_plus %*% GS + f_zero

  ## ------------------------------------------------------------------
  ## Base residual diagnostics
  ## ------------------------------------------------------------------
  R1_base <- if (n_state > 0) {
    f_plus %*% G %*% SG + f_zero %*% G + f_minus_S
  } else {
    matrix(0, nrow = n_endo, ncol = 0)
  }
  AH_f_u_base <- A %*% H + f_u

  max_resid_R1 <- if (length(R1_base)) max(abs(R1_base)) else 0
  max_resid_AH <- if (length(AH_f_u_base)) max(abs(AH_f_u_base)) else 0

  ## ------------------------------------------------------------------
  ## Build M = (S G)' (x) f_plus + I_{n_state} (x) A  ONCE and factorize.
  ## vec(dG) solves M vec(dG) = -vec(RHS_j) for every parameter j; only
  ## RHS_j changes across parameters, so M is built and QR-factorized a
  ## single time here and reused below (Childers et al. 2022 efficiency
  ## point: the cost of the implicit-function-theorem solve is amortized
  ## over all parameters).
  ## ------------------------------------------------------------------
  ## The dG solve no longer forms the dense Kronecker M = kron(t(SG), f_plus) +
  ## kron(I, A); each parameter's generalized Sylvester equation is solved by
  ## .solve_kron_compact (Schur/QZ) in .solution_deriv_one. M_qr kept as NULL for
  ## the helper's signature.
  M_qr <- NULL

  ## Factor the k=1 generalized-Sylvester operator (A dG + f_plus dG SG = -RHS_j)
  ## ONCE, reused for every parameter's RHS. NULL => the dedicated solver was
  ## unavailable/ill-conditioned and .solution_deriv_one falls back to
  ## .solve_kron_compact.
  sylv_fac <- if (n_state > 0) .gen_sylvester_k1_factor(A, f_plus, SG) else NULL

  ## A's factorization, reused for both the base H = -A^{-1} f_u check and
  ## every dH = -A^{-1}(...) below.
  A_qr <- qr(A)

  ## ------------------------------------------------------------------
  ## Analytic (finite-difference-free) primitive derivatives (Tier 11 #3).
  ## When the model supports symbolic parameter differentiation, build the
  ## TOTAL steady-state sensitivity dys and dynamic-Jacobian derivatives for
  ## all requested parameters ONCE (one static-Jacobian factorization + one
  ## model-Hessian evaluation), replacing the per-parameter steady-state
  ## re-solves and primitive finite differences. Falls back to FD if the
  ## option is off, the model is unsupported, or the analytic build fails.
  ## ------------------------------------------------------------------
  if (is.null(use_analytic))
    use_analytic <- isTRUE(getOption("dynhr.use_analytic_primitives", TRUE))
  analytic <- NULL
  if (use_analytic && .can_use_analytic_primitive_deriv(compiled)) {
    dys_all <- .analytic_dys(compiled, ys, params)
    dprim   <- if (!is.null(dys_all))
      .analytic_dprimitives(compiled, ys, params, dys_all) else NULL
    if (!is.null(dys_all) && !is.null(dprim))
      analytic <- list(dys = dys_all, dprim = dprim)
  }

  ## ------------------------------------------------------------------
  ## Per-parameter implicit-differentiation derivatives
  ## ------------------------------------------------------------------
  derivs <- vector("list", length(param_names))
  names(derivs) <- param_names

  for (j in seq_along(param_names)) {
    pname <- param_names[j]
    precomp <- if (!is.null(analytic)) {
      ## Preserve the endogenous-variable names: single-column matrix
      ## extraction drops them when the matrix has one row (n_endo == 1),
      ## which would later make the name-indexed dd = dys[obs_vars] return NA.
      c(list(dys = setNames(as.numeric(analytic$dys[, pname]),
                            rownames(analytic$dys))),
        analytic$dprim[[pname]])
    } else NULL
    res <- tryCatch({
      .solution_deriv_one(pname, model, compiled, ys, params,
                           f_plus, f_zero, f_minus, f_u,
                           G, H, S, SG, A, A_qr, M_qr, sylv_fac,
                           state_idx, n_state, n_endo, n_exo, h_rel,
                           precomp = precomp)
    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })

    if (!isTRUE(res$ok)) {
      warning(sprintf(
        "solution_derivatives: parameter '%s' failed (%s); returning NA derivatives.",
        pname, res$message %||% "unknown error"))
      derivs[[pname]] <- list(
        dG = matrix(NA_real_, n_endo, n_state),
        dH = matrix(NA_real_, n_endo, n_exo),
        dys = setNames(rep(NA_real_, n_endo), endo),
        dTT = matrix(NA_real_, n_state, n_state),
        dRR = matrix(NA_real_, n_state, n_exo),
        dZZ = matrix(NA_real_, length(obs_idx), n_state),
        dDD = matrix(NA_real_, length(obs_idx), n_exo),
        dd  = setNames(rep(NA_real_, length(obs_vars)), obs_vars),
        ok = FALSE,
        message = res$message %||% "unknown error"
      )
      next
    }

    dG  <- res$dG
    dH  <- res$dH
    dys <- res$dys

    derivs[[pname]] <- list(
      dG  = dG,
      dH  = dH,
      dys = dys,
      dTT = if (n_state > 0) dG[state_idx, , drop = FALSE] else matrix(0, 0, 0),
      dRR = if (n_state > 0) dH[state_idx, , drop = FALSE] else matrix(0, 0, n_exo),
      dZZ = dG[obs_idx, , drop = FALSE],
      dDD = dH[obs_idx, , drop = FALSE],
      dd  = dys[obs_vars],
      ok  = TRUE,
      message = NA_character_
    )
  }

  list(
    derivs = derivs,
    base_residual = list(R1 = max_resid_R1, AH_f_u = max_resid_AH),
    reused_factorization = TRUE,
    used_analytic = !is.null(analytic)
  )
}


#' Implicit-differentiation derivative for a single parameter.
#'
#' Computes the TOTAL df-matrices (and dys) by central finite differences of
#' the smooth dynamic-Jacobian primitive (re-evaluated at perturbed steady
#' states), then solves the linear implicit-function-theorem systems for
#' dG and dH using the pre-built/factorized M and A.
#'
#' @return list(ok = TRUE, dG, dH, dys) on success, or
#'   list(ok = FALSE, message = "...") if the perturbed steady state fails to
#'   converge.
#' @noRd
.solution_deriv_one <- function(pname, model, compiled, ys, params,
                                 f_plus, f_zero, f_minus, f_u,
                                 G, H, S, SG, A, A_qr, M_qr, sylv_fac,
                                 state_idx, n_state, n_endo, n_exo, h_rel,
                                 precomp = NULL) {

  if (!is.null(precomp)) {
    ## ANALYTIC path (Tier 11 #3): the TOTAL steady-state and dynamic-Jacobian
    ## derivatives are supplied symbolically (param_jacobian_fn + model Hessian
    ## contracted with dys), with no per-parameter steady-state re-solve or
    ## primitive finite difference. See R/gradient-primitive-deriv.R.
    dys      <- precomp$dys
    df_plus  <- precomp$df_plus
    df_zero  <- precomp$df_zero
    df_minus <- precomp$df_minus
    df_u     <- precomp$df_exo
  } else {
    theta_j <- params[[pname]]
    h <- max(h_rel * abs(theta_j), 1e-7)

    params_p <- params; params_p[[pname]] <- theta_j + h
    params_m <- params; params_m[[pname]] <- theta_j - h

    ## Perturbed steady states, warm-started from the base ss.
    ss_p <- solve_steady(compiled, params_p, y0 = ys,
                          endo_names = model$var_names,
                          exo_names = model$varexo_names, verbose = FALSE)
    if (!isTRUE(ss_p$converged))
      return(list(ok = FALSE,
                   message = sprintf("steady state did not converge at %s + h", pname)))

    ss_m <- solve_steady(compiled, params_m, y0 = ys,
                          endo_names = model$var_names,
                          exo_names = model$varexo_names, verbose = FALSE)
    if (!isTRUE(ss_m$converged))
      return(list(ok = FALSE,
                   message = sprintf("steady state did not converge at %s - h", pname)))

    ys_p <- ss_p$values
    ys_m <- ss_m$values

    ## TOTAL derivative of the steady state wrt theta_j.
    dys <- (ys_p - ys_m) / (2 * h)

    ## TOTAL derivatives of the dynamic-Jacobian blocks, evaluated via central
    ## FD of the SMOOTH primitive (the Jacobian evaluator), not of the QZ
    ## solve. Each evaluation already accounts for the steady-state shift
    ## through ss_p / ss_m and the parameter shift through params_p / params_m.
    ## Re-derive any steady_state_model-computed parameter at the perturbed point
    ## so the dynamic Jacobian is evaluated with the consistent (not stale) p_c
    ## (no-op for non-SSM-parameter models).
    sys_p <- extract_system_matrices(compiled, ys_p, .ssm_consistent_params(model, params_p))
    sys_m <- extract_system_matrices(compiled, ys_m, .ssm_consistent_params(model, params_m))

    df_plus  <- (sys_p$f_plus  - sys_m$f_plus)  / (2 * h)
    df_zero  <- (sys_p$f_zero  - sys_m$f_zero)  / (2 * h)
    df_minus <- (sys_p$f_minus - sys_m$f_minus) / (2 * h)
    df_u     <- (sys_p$f_exo   - sys_m$f_exo)   / (2 * h)
  }

  df_minus_S <- if (n_state > 0) df_minus[, state_idx, drop = FALSE]
                else matrix(0, nrow = n_endo, ncol = 0)

  if (n_state == 0) {
    ## No state variables: G is n_endo x 0, dG is trivially empty.
    dG <- matrix(0, nrow = n_endo, ncol = 0)
  } else {
    ## RHS_j = df_plus G (S G) + df_zero G + df_minus_S   (n_endo x n_state)
    RHS_j <- df_plus %*% G %*% SG + df_zero %*% G + df_minus_S

    ## Solve the generalized Sylvester equation
    ##   A dG + f_plus dG SG = -RHS_j
    ## via the compact Schur solver (complex QZ of the (A, f_plus) pencil + Schur
    ## of SG) instead of forming and factorizing the dense
    ## (n_endo*n_state) x (n_endo*n_state) Kronecker matrix
    ##   M = kron(t(SG), f_plus) + kron(I_{n_state}, A),
    ## using the identity  M vec(dG) = vec(A dG + f_plus dG SG).
    ## The dense M factorization was the dominant gradient cost (a ~2700x2700
    ## QR on NZSIM, O((n_endo*n_state)^3)); .solve_kron_compact is
    ## O(n_endo^3 + n_state^3) and falls back to the dense Schur solve when the
    ## compact route is ill-conditioned, so the result is never wrong.
    dG <- if (!is.null(sylv_fac))
      .gen_sylvester_k1_solve(sylv_fac, -RHS_j)            # dedicated real-Schur
    else
      .solve_kron_compact(A, f_plus, SG, k = 1L, RHS = -RHS_j)  # fallback
  }

  ## dA = df_plus G S + f_plus dG S + df_zero
  GdS  <- if (n_state > 0) G %*% S else matrix(0, nrow = n_endo, ncol = n_endo)
  dGS  <- if (n_state > 0) dG %*% S else matrix(0, nrow = n_endo, ncol = n_endo)
  dA <- df_plus %*% GdS + f_plus %*% dGS + df_zero

  ## dH = -A^{-1} ( dA H + df_u ), reusing A's factorization.
  rhs_H <- dA %*% H + df_u
  dH <- -qr.solve(A_qr, rhs_H)

  list(ok = TRUE, dG = dG, dH = dH, dys = dys)
}
