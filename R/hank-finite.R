## R/hank-finite.R
## --------------------------------------------------------------------------
## Promotion of the finite-HANK .mod emitter (Tier 18, milestone m5, orders
## 2/3 scope) into the package proper. Scratch derivation, GO/verification
## numbers, and the full design rationale:
##   .claude/orchestration/truncation/m5_orders23_SCOPE.md  (spec + ladder)
##   .claude/orchestration/truncation/m5_mod_emit.R          (working emitter)
##   .claude/orchestration/truncation/m5_o23_verify.R        (PF-asymmetry
##     oracle for orders 2/3)
## This file is a faithful PORT of m5_mod_emit.R's `emit_hank_mod()` and
## driver into two exported entry points; it does not redesign the emission.
##
## Family caveat (IMPORTANT for callers): the direct-Euler finite model
## emitted here is a DIFFERENT member of the coarse-grid family than the
## EGM-step model behind `hank_reiter_linearize`/`hank_reiter_statespace`
## (R/hank-truncation.R) -- the two interpolate the household problem in
## different spaces, so the EGM steady state satisfies the direct-Euler
## equations only up to an interpolation-consistency gap, O(grid^2), not
## exactly. Verified drift at n_a = 8 (`hank_ks_coarse_anchored` grid):
##   ap 5.6e-3, Dm 3.6e-4, K 1.2e-2, r 2.7e-5 (see hank_finite_solve()).
## Order-1 dr IRFs (AA, CC, K, r, w, YY) agree with the EGM-family Klein
## `hank_reiter_statespace` IRFs at the family-gap level: max rel diff
## 2.5e-4 -- 6.2e-4 (K impact exactly 0 in both -- timing is consistent).
##
## Order-2/3 cost table (n_a = 8, n_e = 2, 38 vars unless noted):
##   order-1 solve:  parse + compile + steady + dr1, well under 1s.
##   order-2 solve:  compile(max_order=2) 4.4s + solve 0.3s (ghxx 38x1024,
##                   all finite, sane magnitude). PF-asymmetry oracle
##                   (m5_o23_verify.R): pruned order-2 EVEN part vs
##                   perfect-foresight of the SAME .mod -- rel err
##                   6.4e-6..1.3e-5 at shock size 1%, 7.4e-5..9.1e-5 at 3%
##                   (residual scales ~s^2, as theory requires).
##   order-2 @ n_a=16 (70 vars): compile 8.8s + solve 81s -- fine for
##                   one-off calibration/IRF work, not per-draw MCMC.
##   order-3:        DENSE PATH WALL -- >4h CPU / 27GB RSS at 38 vars
##                   (ghxxx is 38 x 32768 with an n^3-RHS Sylvester solve).
##                   NOW SOLVED via the SPARSE complex-Schur Kronecker route
##                   (.solve_kron_compact_sparse, never materialises the
##                   ns^3 x ns^3 object): hank_finite_solve(ks, order = 3L,
##                   sparse = TRUE) runs n_a = 8 in ~20-90s, few-hundred-MB
##                   RSS. Requires an EXPLICIT sparse = TRUE opt-in; the dense
##                   path (sparse = FALSE / NULL) still fails loud.
## --------------------------------------------------------------------------

#' Emit the anchored coarse-grid finite-state HANK as a dynhr .mod
#'
#' Direct port of the Tier-18/m5 prototype emitter
#' (\code{.claude/orchestration/truncation/m5_mod_emit.R}). Writes the
#' finite-state HANK's EQUILIBRIUM CONDITIONS directly (the standard way
#' Reiter-style models are fed to Dynare-style parsers), rather than
#' symbolically differentiating the EGM backward step (which is not
#' expressible: endogenous-grid inversion has no closed form).
#'
#' Per grid node \code{i = (e, a)} (asset index fastest, matching
#' \code{.hank_mat_to_vec}): an Euler equation for unconstrained nodes with
#' the interpolation bracket FROZEN at the steady-state policy's position
#' (\code{findInterval(ap_ss(i), a_grid)}), or \code{ap_i = amin} for
#' constrained nodes (frozen constraint set). The distribution follows the
#' Young lottery with linear-in-\code{ap} weights at the SAME frozen
#' brackets, and the LAST node's mass is eliminated
#' (\code{Dm_n = 1 - sum(others)}) so the mass-conservation unit root never
#' enters the emitted system. Aggregates \code{K}, prices \code{r, w}, TFP
#' \code{z} (log AR(1)), and observables \code{AA, CC, YY} close the model.
#'
#' Frozen brackets/constraint set are the standard Reiter locality
#' assumption: they make the emitted system smooth (no kinks to
#' differentiate through) at the cost of an interpolation-consistency gap
#' between the EGM steady state and the direct-Euler steady state that grows
#' with the shock/perturbation order (see file header and
#' \code{\link{hank_finite_solve}}).
#'
#' @param ks A \code{hank_ks} object (typically from
#'   \code{\link{hank_ks_coarse_anchored}}) supplying the household block,
#'   steady-state policy/distribution, and firm-block calibration
#'   (\code{alpha}, \code{delta}, \code{Z}, \code{K}, \code{r}, \code{w}).
#' @param rho_z,sigma_z TFP log-AR(1) persistence and innovation stderr.
#' @param path Optional file path to also write the emitted .mod text to
#'   (via \code{writeLines}). When \code{NULL} (default) nothing is written
#'   and only the in-memory result is returned.
#' @param order Perturbation order passed to the emitted
#'   \code{stoch_simul(order=, irf=0)} statement (only relevant if the
#'   caller subsequently drives the .mod through \code{stoch_simul()}
#'   directly; \code{\link{hank_finite_solve}} ignores it and re-derives
#'   order via \code{compile_model(max_order=)}).
#'
#' @param instrument \code{"none"} (default; the original emission) or
#'   \code{"captax"}: adds a linear capital-income tax \code{tau} rebated
#'   lump-sum (\code{tr}) to every household. The household budget becomes
#'   \code{cc_i = (1 + (1-tau)*r)*a_i + w*e_i + tr - ap_i}, the Euler
#'   condition prices next period's return at \code{(1 + (1-tau(+1))*r(+1))},
#'   and the government runs a balanced budget \code{tr = tau*r*K}. At
#'   \code{tau = 0} the model is EXACTLY the \code{instrument = "none"}
#'   emission (the resource constraint is unchanged by the rebate). A
#'   distribution-weighted utilitarian \code{planner_objective(sum_i Dm_i
#'   u(cc_i))} is emitted alongside, so the model feeds
#'   \code{\link{ramsey_model}} directly (Aiyagari-1995-style optimal
#'   capital taxation on the finite HANK).
#'
#'   VALIDITY CAVEAT (verified 2026-07-07, n_a = 8 anchored calibration):
#'   the frozen brackets/constraint set make this emission a LOCAL model
#'   around the competitive steady state. The unrestricted utilitarian
#'   level-Ramsey optimum is NOT local: the SS-welfare gradient
#'   \code{dW/dtau < 0} at \code{tau = 0} (golden-rule under-accumulation:
#'   a capital subsidy raises K and w), and the augmented-FOC Newton follows
#'   it to \code{tau ~ -0.37}, where the frozen brackets no longer describe
#'   the true finite HANK and the augmented system fails Blanchard-Kahn.
#'   Level-Ramsey therefore requires re-solving the RAMSEY steady state in
#'   the global (EGM) family first and re-emitting around it
#'   (LeGrand-Ragot-style) -- not yet implemented. What IS valid locally:
#'   SS-PRESERVING stabilization instruments (e.g. \code{tau_rule =
#'   "tau = phi*z;"}, since \code{z_ss = 0}), whose welfare ranking is
#'   exercised in \code{tests/testthat/test-hank-finite-ramsey.R}.
#' @param tau_rule Only with \code{instrument = "captax"}. \code{NULL}
#'   (default) leaves \code{tau} WITHOUT an equation -- a free Ramsey
#'   instrument (one fewer equation than variables, the
#'   \code{\link{ramsey_augment_mod}} convention; note
#'   \code{\link{hank_finite_solve}}/\code{stoch_simul} cannot solve that
#'   non-square system). Alternatively a character equation pinning the
#'   instrument, e.g. \code{"tau = 0;"} (competitive benchmark, square
#'   model) or \code{"tau = tau_exo;"} (a \code{varexo tau_exo} is added
#'   automatically -- lets \code{\link{perfect_foresight_solve}} impose an
#'   arbitrary deterministic tax path, e.g. for direct-optimization
#'   cross-checks of the Ramsey solution).
#' @return A list with \code{text} (character vector, the .mod source,
#'   one element per line), \code{path} (the file path written to, or
#'   \code{NULL}), \code{vars} (endogenous variable names in declaration
#'   order), \code{n} (grid size \code{n_e * n_a}), \code{brk} (frozen
#'   interpolation brackets), \code{constrained} (logical, per node), and
#'   \code{lottery_check} (max absolute error of the frozen-bracket lottery
#'   reproduction of \eqn{\Lambda' D_{ss}}; should be at machine precision).
#' @export
hank_finite_mod <- function(ks, rho_z, sigma_z, path = NULL, order = 1L,
                            instrument = c("none", "captax"),
                            tau_rule = NULL) {
  if (!inherits(ks, "hank_ks"))
    stop("hank_finite_mod(): `ks` must be a hank_ks object.")
  instrument <- match.arg(instrument)
  captax <- instrument == "captax"
  if (!captax && !is.null(tau_rule))
    stop("hank_finite_mod(): `tau_rule` requires instrument = \"captax\".")

  num <- function(x) sprintf("%.17g", x)

  blk <- ks$block
  n_e <- blk$n_e; n_a <- blk$n_a; n <- n_e * n_a
  ag  <- blk$a_grid; Pi <- blk$Pi; e_g <- blk$e
  amin <- ag[1L]
  a_ss <- as.numeric(t(blk$a))    ## policy on nodes, asset fastest
  D_ss <- as.numeric(t(blk$D))
  c_ss <- as.numeric(t(blk$c))
  ei <- function(i) ((i - 1L) %/% n_a) + 1L
  ai <- function(i) ((i - 1L) %% n_a) + 1L
  idx <- function(j, k) (j - 1L) * n_a + k

  ## SS objects for prices
  s_e <- { v <- Re(eigen(t(Pi))$vectors[, 1]); v / sum(v) }
  L    <- sum(s_e * e_g)
  Zbar <- ks$Z

  ## frozen brackets: k_i s.t. a'_ss(i) in [ag[k_i], ag[k_i+1]]
  brk <- pmin(pmax(findInterval(a_ss, ag), 1L), n_a - 1L)
  constrained <- a_ss <= amin + 1e-12

  ## --- pre-check: frozen-bracket lottery == Young lottery at SS ------------
  D_next <- numeric(n)
  for (i in seq_len(n)) {
    om <- (ag[brk[i] + 1L] - a_ss[i]) / (ag[brk[i] + 1L] - ag[brk[i]])
    for (j in seq_len(n_e)) {
      p <- Pi[ei(i), j]
      if (p == 0) next
      D_next[idx(j, brk[i])]      <- D_next[idx(j, brk[i])]      + p * om * D_ss[i]
      D_next[idx(j, brk[i] + 1L)] <- D_next[idx(j, brk[i] + 1L)] + p * (1 - om) * D_ss[i]
    }
  }
  lam_push <- as.numeric(Matrix::t(blk$Lambda) %*% D_ss)
  lottery_check <- max(abs(D_next - lam_push))
  if (lottery_check > 1e-8) {
    stop(sprintf(
      "hank_finite_mod(): frozen-bracket lottery does not reproduce Lambda'D_ss (max abs err %.3e). ",
      lottery_check),
      "This indicates the `ks` steady state / grid is inconsistent with the ",
      "frozen-bracket assumption.")
  }

  ## --- symbol helpers -------------------------------------------------------
  apn <- function(i, t = 0L)
    paste0("ap_", i, c("(-1)", "", "(+1)")[t + 2L])
  ## Dm_i with the LAST node substituted out
  dmn <- function(i, t = 0L) {
    sfx <- c("(-1)", "", "(+1)")[t + 2L]
    if (i < n) paste0("Dm_", i, sfx)
    else paste0("(1", paste0(" - Dm_", seq_len(n - 1L), sfx, collapse = ""), ")")
  }
  om_expr <- function(i, t = 0L)   ## interpolation weight of ap_i on ag[brk[i]]
    sprintf("((%s - %s)/%s)", num(ag[brk[i] + 1L]), apn(i, t),
            num(ag[brk[i] + 1L] - ag[brk[i]]))
  ## after-tax gross return factor at t / t+1 (identical to the no-instrument
  ## emission when tau == 0)
  ret0 <- if (captax) "(1+(1-tau)*r)"        else "(1+r)"
  ret1 <- if (captax) "(1+(1-tau(+1))*r(+1))" else "(1+r(+1))"
  tr0  <- if (captax) " + tr"      else ""
  tr1  <- if (captax) " + tr(+1)"  else ""
  cc_expr <- function(i)           ## consumption at node i, time t
    sprintf("(%s*%s + w*%s%s - %s)", ret0, num(ag[ai(i)]), num(e_g[ei(i)]),
            tr0, apn(i))

  gam <- 1 / blk$eis

  eqs <- character(0)

  ## --- policy equations -----------------------------------------------------
  for (i in seq_len(n)) {
    if (constrained[i]) {
      eqs <- c(eqs, sprintf("ap_%d = %s;", i, num(amin)))
    } else {
      terms <- character(0)
      for (j in seq_len(n_e)) {
        p <- Pi[ei(i), j]
        if (p == 0) next
        apint <- sprintf("(%s*%s + (1-%s)*%s)",
                         om_expr(i), apn(idx(j, brk[i]), 1L),
                         om_expr(i), apn(idx(j, brk[i] + 1L), 1L))
        cint <- sprintf("(%s*%s + w(+1)*%s%s - %s)",
                        ret1, apn(i), num(e_g[j]), tr1, apint)
        terms <- c(terms, sprintf("%s*%s^(-%s)", num(p), cint, num(gam)))
      }
      eqs <- c(eqs, sprintf("%s^(-%s) = beta*%s*(%s);",
                            cc_expr(i), num(gam), ret1,
                            paste(terms, collapse = " + ")))
    }
  }

  ## --- distribution equations (targets 1..n-1; node n eliminated) ----------
  contrib <- vector("list", n)                 ## target -> character terms
  for (i in seq_len(n)) {
    for (j in seq_len(n_e)) {
      p <- Pi[ei(i), j]
      if (p == 0) next
      tlo <- idx(j, brk[i]); thi <- idx(j, brk[i] + 1L)
      contrib[[tlo]] <- c(contrib[[tlo]],
        sprintf("%s*%s*%s", num(p), om_expr(i, -1L), dmn(i, -1L)))
      contrib[[thi]] <- c(contrib[[thi]],
        sprintf("%s*(1-%s)*%s", num(p), om_expr(i, -1L), dmn(i, -1L)))
    }
  }
  for (tg in seq_len(n - 1L)) {
    rhs <- if (length(contrib[[tg]]) == 0) "0" else paste(contrib[[tg]], collapse = " + ")
    eqs <- c(eqs, sprintf("Dm_%d = %s;", tg, rhs))
  }

  ## --- aggregates + prices + shock ------------------------------------------
  Ksum <- paste(vapply(seq_len(n), function(i)
    sprintf("%s*%s", num(ag[ai(i)]), dmn(i)), ""), collapse = " + ")
  eqs <- c(eqs, sprintf("K = %s;", Ksum))
  eqs <- c(eqs, sprintf("r = alpha*%s*exp(z)*K^(alpha-1)*%s - delta;",
                        num(Zbar), num(L^(1 - ks$alpha))))
  eqs <- c(eqs, sprintf("w = (1-alpha)*%s*exp(z)*K^alpha*%s;",
                        num(Zbar), num(L^(-ks$alpha))))
  eqs <- c(eqs, "z = rho_z*z(-1) + eps_z;")
  AAsum <- paste(vapply(seq_len(n), function(i)
    sprintf("%s*ap_%d", dmn(i), i), ""), collapse = " + ")
  eqs <- c(eqs, sprintf("AA = %s;", AAsum))
  CCsum <- paste(vapply(seq_len(n), function(i)
    sprintf("%s*%s", dmn(i), cc_expr(i)), ""), collapse = " + ")
  eqs <- c(eqs, sprintf("CC = %s;", CCsum))
  eqs <- c(eqs, sprintf("YY = %s*exp(z)*K^alpha*%s;",
                        num(Zbar), num(L^(1 - ks$alpha))))

  ## --- capital-tax instrument block ------------------------------------------
  exos <- "eps_z"
  if (captax) {
    ## Balanced budget: the tax base is this period's capital income r*K
    ## (K = sum_i Dm_i a_i is beginning-of-period assets, same timing as the
    ## household's taxed return).
    eqs <- c(eqs, "tr = tau*r*K;")
    if (!is.null(tau_rule)) {
      tau_rule <- trimws(tau_rule)
      if (!grepl(";$", tau_rule)) tau_rule <- paste0(tau_rule, ";")
      eqs <- c(eqs, tau_rule)
      if (grepl("\\btau_exo\\b", tau_rule)) exos <- c(exos, "tau_exo")
    }
    ## else: tau is a FREE Ramsey instrument (n_eq = n_var - 1).
  }

  ## --- assemble .mod ---------------------------------------------------------
  vars <- c(paste0("ap_", seq_len(n)), paste0("Dm_", seq_len(n - 1L)),
            "K", "z", "r", "w", "AA", "CC", "YY",
            if (captax) c("tau", "tr"))
  ini <- c(sprintf("ap_%d = %s;", seq_len(n), vapply(a_ss, num, "")),
           sprintf("Dm_%d = %s;", seq_len(n - 1L), vapply(D_ss[seq_len(n - 1L)], num, "")),
           sprintf("K = %s;", num(ks$K)),
           "z = 0;",
           sprintf("r = %s;", num(ks$r)),
           sprintf("w = %s;", num(ks$w)),
           sprintf("AA = %s;", num(sum(D_ss * a_ss))),
           sprintf("CC = %s;", num(sum(D_ss * c_ss))),
           sprintf("YY = %s;", num(Zbar * ks$K^ks$alpha * L^(1 - ks$alpha))),
           ## initval for (tau, tr): default to the COMPETITIVE benchmark
           ## (tau = 0, tr = 0) unless `ks` itself carries a nonzero tax --
           ## e.g. a hank_ks_taxed object anchored at a RAMSEY steady state
           ## (see R/hank-finite-ramsey-ss.R), whose $tau/$tr are the values
           ## the frozen brackets are built around. Purely additive: any
           ## `ks` without a `$tau` field (every existing caller) is
           ## byte-identical to the old hardcoded "tau = 0; tr = 0;".
           if (captax) {
             tau0 <- if (!is.null(ks$tau)) ks$tau else 0
             tr0v <- if (!is.null(ks$tr))  ks$tr  else 0
             c(sprintf("tau = %s;", num(tau0)), sprintf("tr = %s;", num(tr0v)))
           })

  ## Distribution-weighted utilitarian planner objective (only meaningful
  ## with the instrument; sum_i Dm_i u(cc_i) with the last node's mass
  ## reconstructed, exactly as in the CC equation).
  planner_line <- character(0)
  if (captax) {
    u_expr <- function(i) {
      if (abs(gam - 1) < 1e-12) sprintf("ln(%s)", cc_expr(i))
      else sprintf("(%s^(1-%s))/(1-%s)", cc_expr(i), num(gam), num(gam))
    }
    obj <- paste(vapply(seq_len(n), function(i)
      sprintf("%s*%s", dmn(i), u_expr(i)), ""), collapse = " + ")
    planner_line <- c("", sprintf("planner_objective(%s);", obj))
  }

  txt <- c(
    "// AUTO-GENERATED by hank_finite_mod() -- finite-state HANK (direct-Euler,",
    "// frozen brackets/constraint set). Do not edit by hand.",
    sprintf("var %s;", paste(vars, collapse = " ")),
    sprintf("varexo %s;", paste(exos, collapse = " ")),
    "parameters beta eis alpha delta rho_z;",
    sprintf("beta = %s; eis = %s; alpha = %s; delta = %s; rho_z = %s;",
            num(blk$beta), num(blk$eis), num(ks$alpha), num(ks$delta), num(rho_z)),
    "", "model;", eqs, "end;",
    "", "initval;", ini, "end;",
    planner_line,
    if (captax && is.null(tau_rule)) character(0) else c("", "steady;"),
    "", "shocks;", sprintf("var eps_z; stderr %s;", num(sigma_z)), "end;",
    if (captax && is.null(tau_rule)) character(0)
    else c("", sprintf("stoch_simul(order=%d, irf=0);", as.integer(order))))

  if (!is.null(path)) writeLines(txt, path)

  list(text = txt, path = path, vars = vars, n = n, brk = brk,
       constrained = constrained, lottery_check = lottery_check,
       instrument = instrument, tau_rule = tau_rule)
}


#' Solve the emitted finite-state HANK through the dynhr perturbation pipeline
#'
#' Wraps \code{\link{hank_finite_mod}} + \code{parse_mod} ->
#' \code{compile_model(max_order=)} -> \code{solve_steady} (from the EGM
#' initval baked into the .mod) -> the order-appropriate
#' \code{solve_perturbation}/\code{solve_perturbation_order2}/
#' \code{solve_perturbation_order3} chain.
#'
#' \code{stoch_simul()} itself is order-1-only, so orders 2/3 MUST go through
#' \code{compile_model(max_order = k)} followed by the explicit
#' \code{solve_perturbation() -> solve_perturbation_order2(dr1=) ->
#' solve_perturbation_order3(dr2=)} chain (see
#' \code{.claude/orchestration/truncation/m5_orders23_SCOPE.md}).
#'
#' Order 3 requires the SPARSE Kronecker route (\code{sparse = TRUE}): the
#' emitted system's DENSE third-order solve (\code{ghxxx}) is a documented cost
#' wall (>4h CPU / 27GB RSS at 38 vars; see file header) because the dense
#' fallback forms an 8.6 GB ns^3 × ns^3 Kronecker matrix. The sparse route
#' (\code{.solve_kron_compact_sparse}, a complex-Schur Kronecker
#' Bartels–Stewart) never materialises it and solves n_a = 8 order-3 in ~1–2
#' min. Order 3 with the dense path (\code{sparse = FALSE}) still fails loud.
#'
#' @param ks A \code{hank_ks} object (see \code{\link{hank_finite_mod}}).
#' @param order Perturbation order: \code{1}, \code{2}, or \code{3}
#'   (integer-valued; \code{3} requires \code{sparse = TRUE}; non-integer
#'   values error).
#' @param sparse Logical or \code{NULL} (default). Passed to
#'   \code{\link{solve_perturbation_order3}} for the ghxxx solve at
#'   \code{order = 3}. Order 3 requires an EXPLICIT \code{sparse = TRUE}: the
#'   default (\code{NULL}) and \code{FALSE} both error up front so a caller
#'   never silently hits the dense >4h / 27 GB wall. With \code{TRUE} the
#'   memory-light complex-Schur Kronecker route is used (~20-90 s at
#'   \code{n_a = 8}). Ignored for orders 1 and 2.
#' @param rho_z,sigma_z TFP calibration passed to \code{\link{hank_finite_mod}}.
#' @param path Optional file path to write the emitted .mod to (passed
#'   through to \code{\link{hank_finite_mod}}); when \code{NULL} (default) a
#'   temporary file is used and removed is left to the caller's session
#'   temp-dir cleanup.
#' @param verbose Passed through to \code{parse_mod}/\code{compile_model}/
#'   \code{solve_steady}/\code{solve_perturbation*}.
#' @param instrument,tau_rule Passed through to \code{\link{hank_finite_mod}}.
#'   Because this function solves the SQUARE competitive system, the
#'   \code{"captax"} instrument requires a pinning \code{tau_rule} here
#'   (e.g. \code{"tau = 0;"}); the free-instrument emission is consumed by
#'   \code{\link{ramsey_model}} instead.
#'
#' @return A list with \code{model} (the parsed model, i.e. \code{parse_mod}
#'   output), \code{compiled}, \code{steady} (the \code{solve_steady} result,
#'   with \code{$values} the polished direct-Euler steady state), \code{dr}
#'   (the highest-order \code{DecisionRules*} object: \code{DecisionRules}
#'   for \code{order = 1}, \code{DecisionRules2} for \code{order = 2}),
#'   \code{dr1} (always the order-1 decision rules, even when
#'   \code{order = 2}, for convenience), and \code{emit} (the
#'   \code{\link{hank_finite_mod}} result: \code{vars}, \code{brk},
#'   \code{constrained}, \code{lottery_check}).
#' @export
hank_finite_solve <- function(ks, order = 1L, rho_z = 0.9, sigma_z = 0.01,
                              path = NULL, verbose = FALSE,
                              instrument = c("none", "captax"),
                              tau_rule = NULL, sparse = NULL) {
  instrument <- match.arg(instrument)
  if (instrument == "captax" && is.null(tau_rule))
    stop("hank_finite_solve(): the FREE-instrument captax emission is not a ",
         "square system -- pass a `tau_rule` (e.g. \"tau = 0;\") here, or ",
         "feed hank_finite_mod(instrument = \"captax\") to ramsey_model().")
  if (length(order) != 1L || !is.finite(order) || order != as.integer(order))
    stop("hank_finite_solve(): `order` must be a single integer value.")
  order <- as.integer(order)
  if (order == 3L && !isTRUE(sparse)) {
    # Order 3 requires an EXPLICIT sparse = TRUE opt-in. The dense path is the
    # documented >4h / 27 GB wall; the default (sparse = NULL) keeps erroring so
    # that hank_finite_solve(ks, order = 3L) never silently hits it -- callers
    # must knowingly choose the sparse route.
    stop("hank_finite_solve(): order = 3 is not supported for the emitted ",
         "finite HANK via the DENSE path -- the dense third-order solve ",
         "(ghxxx) is a documented cost wall (>4h CPU / 27GB RSS at 38 vars: ",
         "see .claude/orchestration/truncation/m5_orders23_SCOPE.md). Pass ",
         "sparse = TRUE to use the memory-light sparse Kronecker route ",
         "(deferred wall now feasible, ~20-90 s at n_a = 8).")
  }
  if (!order %in% c(1L, 2L, 3L))
    stop("hank_finite_solve(): `order` must be 1, 2, or 3 (order 3 requires ",
         "the sparse Kronecker route; see `sparse`).")

  if (is.null(path)) path <- tempfile(fileext = ".mod")
  em <- hank_finite_mod(ks, rho_z = rho_z, sigma_z = sigma_z, path = path,
                        order = 1L, instrument = instrument,
                        tau_rule = tau_rule)

  m  <- parse_mod(path, verbose = verbose)
  cm <- compile_model(m, verbose = verbose, max_order = order)
  ss <- solve_steady(cm, m$param_values, verbose = verbose)
  if (!isTRUE(ss$converged))
    stop("hank_finite_solve(): solve_steady() did not converge on the ",
         "emitted finite HANK's direct-Euler system.")

  dr1 <- solve_perturbation(m, cm, ss$values, m$param_values, verbose = verbose)

  dr <- dr1
  if (order >= 2L) {
    dr <- solve_perturbation_order2(m, cm, ss$values, m$param_values,
                                    dr1 = dr1, verbose = verbose)
  }
  if (order == 3L) {
    dr <- solve_perturbation_order3(m, cm, ss$values, m$param_values,
                                    dr2 = dr, sparse = sparse,
                                    verbose = verbose)
  }

  list(model = m, compiled = cm, steady = ss, dr = dr, dr1 = dr1, emit = em)
}
