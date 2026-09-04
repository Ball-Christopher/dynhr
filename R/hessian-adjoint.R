## R/hessian-adjoint.R
## --------------------------------------------------------------------------
## EXACT posterior Hessian of the Gaussian Kalman-filter log-likelihood
## (ROADMAP Tier 6 #2). Builds on the first-order adjoint gradient
## (.kf_loglik_adjoint) and the first/second-order solution-derivative layers
## (solution_derivatives / solution_derivatives_2).
##
## THE DECOMPOSITION
## -----------------
## The loglik L is a function of the state-space matrices X = (TT,RR,ZZ,DD,d,
## Sigma_e), which in turn depend on theta. The adjoint returns the exact
## gradient grad_j = <G_X, dX_j> where G_X = dL/dX are the filter-gradient
## matrices and dX_j = dX/dtheta_j. Differentiating once more wrt theta_i:
##
##   H[i,j] = <dG_X/dtheta_i, dX_j>  +  <G_X, d2X_ij>
##            \__ filter curvature _/    \_ solution curvature _/
##                  (T1)                        (T2)
##
## T2 is computed exactly: G_X is contracted (by .kf_loglik_adjoint) against
## the second-order solution derivatives d2X_ij from solution_derivatives_2.
##
## T1 = <(Hess_X L) dX_i, dX_j> is the directional derivative of the exact
## adjoint gradient grad_j(X) = <G_X(X), dX_j> along dX_i. It is obtained by a
## central directional difference of the (exact, analytic) adjoint gradient in
## MATRIX space along dX_i -- no model re-solve, one finite-difference level on
## an exact quantity (a Hessian-vector product). The resulting T1 is symmetric
## to machine precision (an internal correctness check); the returned Hessian
## is symmetrised.
## --------------------------------------------------------------------------


#' Exact Hessian of the Kalman-filter log-likelihood wrt the state-space
#' matrices, contracted with first- and second-order solution derivatives.
#'
#' @param Y          n_obs x n_T observation matrix (no NAs).
#' @param ss         list TT, RR, ZZ, DD, d, Sigma_e (base point X).
#' @param dX_list    length-n_par list; element i has dTT,dRR,dZZ,dDD,dd,
#'                   dSigma_e (= dX/dtheta_i; missing blocks treated as zero).
#' @param d2X_list   list keyed "i|j" (i<=j and j<=i) of second-order blocks
#'                   d2TT,d2RR,d2ZZ,d2DD,d2d,d2Sigma_e (= d2X/dtheta_i dtheta_j).
#' @param me_variance scalar measurement-error variance.
#' @param eps        directional FD step for the filter-curvature term (T1).
#' @return n_par x n_par symmetric Hessian of the loglik.
#' @noRd
#' @param t1_method how to compute the filter-curvature term T1: \code{"hvp"}
#'   (default) is a matrix-space directional finite difference of the exact
#'   adjoint gradient (an HVP, residual ~1e-8 at the default \code{eps});
#'   \code{"analytic"} uses the fully-analytic forward-over-reverse second-order
#'   adjoint (\code{.kf_loglik_dG}), which is FD-free at the filter level.
#' @param t2_method how to compute the solution-curvature term T2
#'   \code{= <G_X, d2X_ij>}. \code{"loop"} (default) is the historical path: one
#'   \code{.kf_loglik_adjoint} call per ROW i, each passing that row's
#'   \code{n_par} second-order blocks \code{d2X_ij} for \code{j = 1..np}. The
#'   base-point adjoint bar matrices \code{G_X} are IDENTICAL across those rows
#'   (they depend only on \code{Y}, \code{ss}, \code{me_variance}), so the loop
#'   re-runs the full O(T*n^2) forward+backward filter sweep \code{np} times
#'   redundantly. \code{"contract_once"} runs \code{.kf_loglik_adjoint} ONCE
#'   with \code{return_bars = TRUE} to obtain \code{G_X}, then contracts those
#'   bars against every \code{d2X_ij} block with the SAME Frobenius inner
#'   products the loop path performs internally -- numerically identical to
#'   \code{"loop"} (same bars, same \code{sum(G * d2X)} arithmetic), just
#'   without the \code{np - 1} redundant filter sweeps. Default is \code{"loop"}
#'   so pre-existing callers are byte-for-byte unchanged.
kf_loglik_hessian <- function(Y, ss, dX_list, d2X_list, me_variance = 0,
                              eps = 1e-5, t1_method = c("hvp", "analytic"),
                              t2_method = c("loop", "contract_once")) {
  t1_method <- match.arg(t1_method)
  t2_method <- match.arg(t2_method)
  np <- length(dX_list)
  nm <- names(dX_list)

  ## Normalise a derivative block to the .kf_loglik_adjoint contract.
  norm_blk <- function(b) {
    if (is.null(b)) b <- list()
    list(dTT = b$dTT %||% b$d2TT, dRR = b$dRR %||% b$d2RR,
         dZZ = b$dZZ %||% b$d2ZZ, dDD = b$dDD %||% b$d2DD,
         dd  = b$dd  %||% b$d2d,  dSigma_e = b$dSigma_e %||% b$d2Sigma_e)
  }
  dX <- lapply(dX_list, norm_blk)

  ## T1: directional derivative of grad_j(X) = <G_X, dX_j> along dX_i.
  T1 <- matrix(0, np, np)
  if (t1_method == "analytic") {
    ## Fully-analytic: dG_X along dX_i via the second-order adjoint, then
    ## contract with dX_j. FD-free at the filter level.
    for (i in seq_len(np)) {
      dGi <- .kf_loglik_dG(Y, ss, dX[[i]], me_variance = me_variance)
      for (j in seq_len(np)) {
        b <- dX[[j]]
        z <- function(M, B) if (is.null(B)) 0 else sum(M * B)
        T1[i, j] <- z(dGi$dG_TT, b$dTT) + z(dGi$dG_RR, b$dRR) +
                    z(dGi$dG_ZZ, b$dZZ) + z(dGi$dG_DD, b$dDD) +
                    z(dGi$dg_d, b$dd)   + z(dGi$dG_Sig, b$dSigma_e)
      }
    }
  } else {
    perturb <- function(sgn, i) {
      b <- dX[[i]]
      z <- function(v, like) if (is.null(v)) 0 * like else v
      list(TT = ss$TT + sgn * eps * z(b$dTT, ss$TT),
           RR = ss$RR + sgn * eps * z(b$dRR, ss$RR),
           ZZ = ss$ZZ + sgn * eps * z(b$dZZ, ss$ZZ),
           DD = ss$DD + sgn * eps * z(b$dDD, ss$DD),
           d  = as.numeric(ss$d) + sgn * eps * z(b$dd, as.numeric(ss$d)),
           Sigma_e = ss$Sigma_e + sgn * eps * z(b$dSigma_e, ss$Sigma_e))
    }
    for (i in seq_len(np)) {
      gp <- .kf_loglik_adjoint(Y, perturb(+1, i), dX, me_variance = me_variance)$grad
      gm <- .kf_loglik_adjoint(Y, perturb(-1, i), dX, me_variance = me_variance)$grad
      T1[i, ] <- (gp - gm) / (2 * eps)
    }
  }

  ## T2: <G_X, d2X_ij> -- contract the base-point adjoint gradient matrices
  ## against the second-order solution derivatives.
  T2 <- matrix(0, np, np)
  if (t2_method == "contract_once") {
    ## CONTRACT-ONCE (Consumer 2): the base-point adjoint bars G_X do NOT
    ## depend on i or j -- they are a function of (Y, ss, me_variance) alone.
    ## Run the O(T*n^2) filter sweep ONCE to get G_X, then form every
    ## T2[i,j] = <G_X, d2X_ij> by the SAME Frobenius products the "loop" path
    ## computes inside .kf_loglik_adjoint's final contraction (its grad[j] is
    ## exactly sum(G_TT*d2TT_j) + ... for the j-th block). This drops np - 1
    ## redundant filter sweeps at zero numerical cost.
    bars <- .kf_loglik_adjoint(Y, ss, list(), me_variance = me_variance,
                               return_bars = TRUE)$bars
    G_TT <- bars$G_TT; G_RR <- bars$G_RR; G_ZZ <- bars$G_ZZ
    G_DD <- bars$G_DD; g_d  <- bars$g_d;  G_Sig <- bars$G_Sig
    for (i in seq_len(np)) {
      for (j in seq_len(np)) {
        b  <- norm_blk(d2X_list[[paste(i, j, sep = "|")]])
        gj <- 0
        if (!is.null(b$dTT))      gj <- gj + sum(G_TT  * b$dTT)
        if (!is.null(b$dRR))      gj <- gj + sum(G_RR  * b$dRR)
        if (!is.null(b$dZZ))      gj <- gj + sum(G_ZZ  * b$dZZ)
        if (!is.null(b$dDD))      gj <- gj + sum(G_DD  * b$dDD)
        if (!is.null(b$dd))       gj <- gj + sum(g_d   * b$dd)
        if (!is.null(b$dSigma_e)) gj <- gj + sum(G_Sig * b$dSigma_e)
        T2[i, j] <- gj
      }
    }
  } else {
    ## LOOP (default, historical): one adjoint call per row i, passing the
    ## n_par second-order blocks d2X_ij for j = 1..np.
    for (i in seq_len(np)) {
      d2row <- lapply(seq_len(np), function(j) {
        key <- paste(i, j, sep = "|")
        norm_blk(d2X_list[[key]])
      })
      T2[i, ] <- .kf_loglik_adjoint(Y, ss, d2row, me_variance = me_variance)$grad
    }
  }

  H <- T1 + T2
  H <- 0.5 * (H + t(H))
  dimnames(H) <- list(nm, nm)
  attr(H, "t1_asymmetry") <- max(abs(T1 - t(T1)))
  H
}


## --------------------------------------------------------------------------
## posterior_hessian(): higher-level wrapper that
##   1. solves the model at `params`,
##   2. builds the full dX_list (dTT,dRR,dZZ,dDD,dd,dSigma_e per param) and
##      d2X_list (d2TT,...,d2Sigma_e per (i,j) pair),
##   3. calls kf_loglik_hessian to get the loglik Hessian,
##   4. optionally adds the analytic Hessian of the log-prior (diagonal for
##      independent priors) if include_prior = TRUE.
##
## The d2Sigma_e blocks are computed by a second FD of .get_shock_cov:
##   3-point diagonal stencil for (i,i) pairs,
##   4-corner mixed stencil for (i,j) pairs.
## This matches the total-derivative convention of the first-order layer
## (which already uses .dSigma_e_fd = central FD of .get_shock_cov).
## Parameters that also move the decision rule (non-sigma params) get BOTH
## solution-derivative blocks AND dSigma_e from FD, consistent with the
## analytic gradient's "adjoint" path (analytic-gradient.R:719-727).
## --------------------------------------------------------------------------

#' Exact posterior Hessian of the KF log-likelihood (and optionally the prior)
#'
#' @param model     parsed model (from \code{parse_mod}).
#' @param compiled  compiled model (from \code{compile_model}).
#' @param dr        first-order decision rule (from \code{solve_perturbation}).
#' @param params    named numeric vector of ALL model parameters at the evaluation
#'                  point.
#' @param param_names character vector of parameter names to differentiate
#'                    (the Hessian is \code{length(param_names) x length(param_names)}).
#' @param obs_vars  character vector of observable variable names.
#' @param data         n_obs x T observation matrix (no NAs).
#' @param me_variance scalar measurement-error variance (default 0).
#' @param include_prior logical; if TRUE, add the analytic Hessian of
#'   \code{log p(theta)} (diagonal for independent priors) to the loglik Hessian.
#' @param prior_spec prior-spec data.frame. Required when \code{include_prior =
#'   TRUE} AND required when \code{check_mode = TRUE} (the mode-criticality
#'   guard needs it to build the full posterior-gradient closure).
#' @param eps       directional FD step for the filter-curvature (T1) term
#'   (default 1e-5).
#' @param h_Sigma_e relative FD step for the dSigma_e / d2Sigma_e computation
#'   (default 1e-5).
#' @param t1_method how the filter-curvature term T1 is formed, forwarded to
#'   \code{kf_loglik_hessian}: \code{"hvp"} (default) is a matrix-space
#'   directional finite difference of the exact adjoint gradient (a
#'   Hessian-vector product); \code{"analytic"} is the fully-analytic
#'   forward-over-reverse second-order adjoint (FD-free at the filter level).
#' @param t2_method how the solution-curvature term
#'   \code{T2[i,j] = <G_X, d2X_ij>} is formed. \code{"loop"} (default) and
#'   \code{"contract_once"} are the two EXACT paths: both materialise the full
#'   second-order solution derivatives \code{d2X_ij} (via
#'   \code{solution_derivatives_2}) and contract them against the base-point
#'   adjoint bars \code{G_X} -- \code{"contract_once"} just avoids the
#'   \code{np - 1} redundant filter sweeps the loop performs (numerically
#'   identical). \code{"hvp_solution"} is the d2X-FREE path: it exploits
#'   \code{T2 = Hessian_theta phi(theta)} with \code{phi(theta) = <G_X, X(theta)>}
#'   and \code{G_X} FROZEN at the base point, so \code{grad phi = <G_X, dX/dtheta>}
#'   is EXACTLY one \code{.solution_adjoint} call (the frozen filter bars fed as
#'   its \code{bars}), and \code{T2} is the central finite-difference Jacobian of
#'   that exact analytic solution-adjoint gradient. This NEVER calls
#'   \code{solution_derivatives_2} and NEVER materialises the \code{O(np^2 . n^2)}
#'   \code{d2X_ij} blocks (a memory win), keeps \code{T1} exact, and is one FD
#'   level on an exact quantity -- the direct analogue of \code{t1_method = "hvp"}
#'   one layer down (the solution map rather than the filter). It costs
#'   \code{2 . np_struct} first-order model re-solves (steady state + order-1
#'   perturbation) and requires every perturbed point to solve (BK-satisfied);
#'   the Sigma_e second-order channel is handled exactly by the same
#'   \code{<G_Sig, d2Sigma_e>} FD stencil the other paths use. A fully-analytic
#'   (FD-free) order-2-in-parameter solution adjoint would remove that one FD
#'   level too, but is a larger reverse-over-reverse build. \code{"adjoint_solution"}
#'   IS that fully-analytic build: it differentiates \code{.solution_adjoint}'s
#'   own gradient formula ANALYTICALLY in theta_i (forward-over-reverse), via
#'   \code{.solution_adjoint_hessian} (gradient-solution-adjoint-order2-param.R).
#'   The base-point bar matrices (\code{bar_df_plus}, \code{bar_df_zero},
#'   \code{bar_df_minus_S}, \code{bar_df_u}) are themselves functions of theta
#'   through two transposed linear solves (a transposed \code{A}-solve and a
#'   transposed generalized-Sylvester solve); differentiating those solves
#'   costs ONE extra transposed solve per parameter, reusing the SAME
#'   base-point factorizations (no new factorization, no per-parameter model
#'   re-solve). The remaining theta-dependence is through the SECOND total
#'   primitive derivatives (\code{d2f_plus_ij}, ...), obtained exactly the way
#'   \code{solution_derivatives_2} obtains them (analytic 5-term contraction
#'   when the compiled model carries \code{param_deriv = "second"} codegen,
#'   else the same 3-point/4-corner FD stencil) -- so \code{"adjoint_solution"}
#'   is FD-free at the SOLUTION-ADJOINT layer even though it may still be
#'   FD-based one layer down (the smooth-primitive second derivative), exactly
#'   like \code{"loop"}/\code{"contract_once"}. It NEVER materialises the
#'   \code{O(np^2 . n^2)} d2X blocks and NEVER re-solves the model. The Sigma_e
#'   second-order channel is handled exactly by the same
#'   \code{<G_Sig, d2Sigma_e>} FD stencil the other paths use. Default
#'   \code{"loop"} keeps existing callers byte-identical.
#' @param h_t2 relative FD step for the \code{t2_method = "hvp_solution"}
#'   solution-adjoint Jacobian (default 1e-5; the per-parameter step is
#'   \code{h_t2 * max(abs(theta_k), 1e-4)}).
#' @param h_rel2 relative FD step for the \code{t2_method = "adjoint_solution"}
#'   SECOND-order smooth-primitive stencil (default 1e-4, matching
#'   \code{solution_derivatives_2}'s default; only used on the FD fallback,
#'   i.e. when the compiled model lacks the order-2 parameter codegen).
#' @param check_mode logical or \code{NULL} (default). \code{NULL} = AUTO:
#'   run the guard whenever \code{prior_spec} is supplied, silently skip when
#'   it is not (backward compatible). When \code{TRUE},
#'   evaluate the analytic posterior gradient at \code{params} (via
#'   \code{make_posterior_grad}, i.e. the SAME gradient machinery used
#'   elsewhere, not a fresh FD path) and compare its norm (restricted to
#'   \code{param_names}) against \code{mode_grad_tol}. A curvature evaluation
#'   is only meaningful AT a critical point of the posterior being
#'   differentiated: the Hessian is a second-order Taylor coefficient around
#'   \code{params}, and if the gradient there is not ~0, the leading term of
#'   the local expansion is the (unaccounted-for) linear term, not the
#'   quadratic one the Hessian describes. This guard exists because of a
#'   concrete failure: a paper computed \code{posterior_hessian} at an
#'   imported Dynare mode that was NOT a critical point of dynhr's own
#'   posterior (adjoint gradient norm ~3900 there, because dynhr's
#'   diffuse-init likelihood differs from Dynare's by a theta-dependent
#'   presample convention) -- the condition number was silently wrong for a
#'   whole draft (1.7e6 vs 6.9e5 after re-moding the point in dynhr; wrong-sign
#'   eigenvalue count 4 vs 2). One cheap adjoint-gradient evaluation would
#'   have caught it immediately.
#' @param mode_grad_tol numeric or \code{NULL} (default). The ABSOLUTE
#'   gradient-norm threshold used by the \code{check_mode} guard
#'   (\code{||grad|| > mode_grad_tol} triggers the warning/error). When
#'   \code{NULL} (the default) it is set to \code{1e-2 * sqrt(p)} where
#'   \code{p = length(param_names)}: at a polished mode the analytic gradient
#'   is typically flat to numerical noise per coordinate (~1e-6--1e-4,
#'   depending on how tightly the optimiser converged), so summed in
#'   quadrature over \code{p} coordinates that noise floor scales like
#'   \code{O(sqrt(p))} -- \code{1e-2} per coordinate is a generous margin
#'   above that floor (does not false-positive on a well-polished mode) while
#'   still being 5-6 orders of magnitude below the ~3900 gradient norm of the
#'   motivating non-critical-point failure (so it reliably flags that case;
#'   a \code{p ~ 30} problem gets a threshold of ~0.055, and ||grad|| ~ 3900
#'   trips it by more than 4 orders of magnitude, while ||grad|| ~ 1e-6 at a
#'   polished mode is ~4 orders of magnitude BELOW the threshold).
#' @param require_mode logical (default \code{FALSE}). When \code{TRUE} and
#'   the \code{check_mode} guard fires, \code{stop()} instead of
#'   \code{warning()}.
#' @return n_par x n_par symmetric matrix; the loglik Hessian
#'   (plus prior Hessian if \code{include_prior = TRUE}).
#'   Carries attribute \code{t1_asymmetry} (scalar) from \code{kf_loglik_hessian},
#'   and attribute \code{second_primitives} ("analytic" when the compiled model
#'   carries the \code{param_deriv = "second"} codegen, "fd" for the stencil
#'   fallback, NA when the T2 layer itself failed). The two are certified to
#'   agree on regular models (nk_small: identical eigen spectra to printed
#'   digits at the mode; the P2 "wrong-sign eigenvalue" was a wrong-evaluation-
#'   point artifact -- estimated shock stds not injected into \code{params} --
#'   caught by the \code{check_mode} guard, NOT a T2 accuracy problem).
#'   When \code{check_mode = TRUE} and the gradient evaluation succeeds, also
#'   carries attribute \code{grad_norm} (scalar; the posterior gradient norm
#'   at \code{params}, restricted to \code{param_names}) and attribute
#'   \code{grad_at_mode} (named numeric vector; the per-parameter gradient
#'   itself, for diagnosis).
#' @export
posterior_hessian <- function(model, compiled, dr, params, param_names,
                              obs_vars, data,
                              me_variance = 0,
                              include_prior = FALSE,
                              prior_spec = NULL,
                              eps = 1e-5,
                              h_Sigma_e = 1e-5,
                              h_t2 = 1e-5,
                              h_rel2 = 1e-4,
                              t1_method = c("hvp", "analytic"),
                              t2_method = c("loop", "contract_once",
                                            "hvp_solution", "adjoint_solution"),
                              check_mode = NULL,
                              mode_grad_tol = NULL,
                              require_mode = FALSE) {

  t1_method <- match.arg(t1_method)
  t2_method <- match.arg(t2_method)
  np  <- length(param_names)
  exo <- model$varexo_names

  ## ------------------------------------------------------------------
  ## Resolve a parameter's scalar value at `params`, with a fail-loud
  ## fallback for stderr-named entries not (yet) present in `params`:
  ## the estimated-shock-std convention (Dynare `stderr <shock>`; see
  ## .apply_theta_to_params / .get_shock_stderr Priority 0) injects the
  ## draw into `params` under the bare shock name, so the normal caller
  ## path (run_mode_finding et al.) already has it there. If a caller
  ## forgot to inject it, fall back to the model's shocks-block stderr
  ## (.get_shock_stderr); only error if the name resolves in neither.
  ## ------------------------------------------------------------------
  .resolve_param_value <- function(nm) {
    if (nm %in% names(params)) return(params[[nm]])
    if (nm %in% exo) {
      se <- .get_shock_stderr(model, exo, params)
      if (nm %in% names(se) && is.finite(se[[nm]])) return(se[[nm]])
    }
    stop(sprintf(
      paste0("posterior_hessian: parameter '%s' is neither a model parameter ",
             "(in `params`) nor a resolvable exogenous-shock stderr (in ",
             "`model$varexo_names`/shocks block). Pass its value in ",
             "`params` under this name (the estimated-shock-std convention: ",
             "inject under the bare shock name, e.g. `params[nm] <- theta[nm]`)."),
      nm), call. = FALSE)
  }

  ## ------------------------------------------------------------------
  ## Mode-critical-point guard (paper gap #1). A Hessian is only meaningful
  ## AT a critical point of the posterior it is differentiating; evaluate the
  ## SAME analytic gradient machinery (make_posterior_grad, R/analytic-
  ## gradient.R) used elsewhere in the package -- no new gradient path -- and
  ## warn (or, if require_mode, stop) when ||grad|| indicates `params` is not
  ## a critical point. See the roxygen above for the motivating incident and
  ## the tolerance derivation.
  ## ------------------------------------------------------------------
  grad_norm_val  <- NA_real_
  grad_at_mode   <- NULL
  ## check_mode = NULL (default) is AUTO: run the guard whenever prior_spec
  ## is available (the gradient closure needs it), silently skip otherwise --
  ## this keeps every pre-guard caller (no prior_spec) working unchanged.
  ## Explicit TRUE demands the guard and fails loud without prior_spec;
  ## explicit FALSE always skips.
  do_check_mode <- if (is.null(check_mode)) !is.null(prior_spec) else isTRUE(check_mode)
  if (do_check_mode) {
    if (is.null(prior_spec))
      stop("posterior_hessian: check_mode = TRUE requires `prior_spec` ",
           "(needed to build the full posterior-gradient closure via ",
           "make_posterior_grad). Pass prior_spec, or leave check_mode = NULL ",
           "to auto-skip the mode-criticality guard when prior_spec is absent.",
           call. = FALSE)

    tol <- mode_grad_tol %||% (1e-2 * sqrt(np))

    grad_ok <- tryCatch({
      full_names <- prior_spec$name
      theta_full <- vapply(full_names, .resolve_param_value, 0.0)
      names(theta_full) <- full_names

      grad_fn <- make_posterior_grad(model, data, prior_spec, obs_vars, compiled,
                                     me_variance = me_variance)
      g_full  <- grad_fn(theta_full)

      g_sub <- g_full[param_names]
      if (anyNA(g_sub))
        stop("gradient has NA/unresolved entries for param_names not in prior_spec")
      TRUE
    }, error = function(e) {
      warning("posterior_hessian: check_mode gradient evaluation failed (",
              conditionMessage(e), "); skipping the mode-criticality guard.",
              call. = FALSE)
      FALSE
    })

    if (isTRUE(grad_ok)) {
      grad_at_mode <- g_sub
      grad_norm_val <- sqrt(sum(g_sub^2))

      ## Sigma_e-mismatch diagnosis (paper P2 gap #7a). The Hessian's Sigma_e
      ## comes from `params` (.get_shock_cov below). The guard's gradient is
      ## built from `prior_spec` via make_posterior_grad, which sets Sigma_e
      ## from the shock-std entries IT sees. If a shock std is INJECTED into
      ## `params` at a non-calibration (e.g. mode) value but is ABSENT from
      ## `prior_spec`, the guard evaluates ||grad|| with Sigma_e at CALIBRATION
      ## while the curvature uses the injected Sigma_e -- an off-mode evaluation
      ## that inflates ||grad|| and cries wolf at a genuine joint mode. Detect
      ## and name it precisely rather than blaming "an unpolished mode".
      ## Convention-agnostic: compare the shock stderrs the HESSIAN uses (from
      ## `params`) against what the GUARD's gradient sees (make_posterior_grad
      ## resolves only `prior_spec$name` from `params`, leaving the rest at the
      ## model defaults). Any shock whose stderr differs between the two means
      ## a Sigma_e-affecting parameter is injected into `params` but omitted
      ## from `prior_spec` -- so ||grad|| was computed off the injected point.
      sigma_e_note <- ""
      mism <- tryCatch({
        se_hess  <- .get_shock_stderr(model, exo, params)
        pv_guard <- model$param_values
        common   <- intersect(prior_spec$name, names(params))
        pv_guard[common] <- params[common]
        se_guard <- .get_shock_stderr(model, exo, pv_guard)
        nm <- intersect(names(se_hess), names(se_guard))
        nm[is.finite(se_hess[nm]) & is.finite(se_guard[nm]) &
           abs(se_hess[nm] - se_guard[nm]) > 1e-10 * pmax(1, abs(se_hess[nm]))]
      }, error = function(e) character(0))
      if (length(mism)) {
        sigma_e_note <- sprintf(
          paste0(" NOTE: the Hessian's Sigma_e (from `params`) differs from the ",
                 "gradient guard's for shock(s) {%s} -- a Sigma_e-affecting ",
                 "parameter is injected into `params` but ABSENT from ",
                 "`prior_spec`, so this ||grad|| was evaluated with Sigma_e at ",
                 "the model default, an off-mode point. Add that parameter to ",
                 "`prior_spec` (and `param_names`) for a faithful joint-mode ",
                 "check; the reported norm is likely spuriously large."),
          paste(mism, collapse = ", "))
      }

      if (grad_norm_val > tol) {
        msg <- sprintf(
          paste0("posterior_hessian: curvature requested at a non-critical ",
                 "point (||grad|| = %.6g); polish the mode first -- an ",
                 "imported (e.g. Dynare) mode is a likelihood-LEVEL ",
                 "validation point, not a critical point of dynhr's ",
                 "posterior (presample/likelihood conventions differ in ",
                 "theta).%s"),
          grad_norm_val, sigma_e_note)
        if (isTRUE(require_mode)) {
          stop(msg, call. = FALSE)
        } else {
          warning(msg, call. = FALSE)
        }
      }
    }
  }

  ## ------------------------------------------------------------------
  ## Classify param_names into structural (move the decision rule) vs
  ## stderr/sigma-like (move ONLY Sigma_e). The naming convention an
  ## estimated shock std is given is the BARE shock name itself (Dynare
  ## `stderr <shock>`; see .apply_theta_to_params / make_log_posterior /
  ## .get_shock_stderr Priority 0) -- matched here by `nm %in% exo`. A
  ## structural parameter that ALSO only moves Sigma_e (e.g. ireland_2004's
  ## `sig_a`, declared in `parameters` and referenced only via
  ## `stderr sig_a;`) is NOT caught by this name test, but it IS present in
  ## `model$param_values` / `names(params)`, so solution_derivatives()
  ## handles it without error (its structural blocks come back exactly
  ## zero via the ordinary FD path) -- unchanged, existing behaviour.
  ## ------------------------------------------------------------------
  is_sig       <- param_names %in% exo
  sig_names    <- param_names[is_sig]
  struct_names <- param_names[!is_sig]

  ## ------------------------------------------------------------------
  ## State-space at the evaluation point.
  ## ------------------------------------------------------------------
  si  <- dr$state_idx
  oi  <- match(obs_vars, dr$endo_names)
  ss0 <- list(
    TT      = dr$ghx[si, , drop = FALSE],
    RR      = dr$ghu[si, , drop = FALSE],
    ZZ      = dr$ghx[oi, , drop = FALSE],
    DD      = dr$ghu[oi, , drop = FALSE],
    d       = dr$ys[obs_vars],
    Sigma_e = .get_shock_cov(model, exo, params)
  )

  ## ------------------------------------------------------------------
  ## First-order solution derivatives (structural blocks) for the
  ## STRUCTURAL subset only. Stderr/sigma-like params (bare shock names)
  ## are skipped here -- their dTT/dRR/dZZ/dDD/dys are exactly zero
  ## (certainty equivalence) and are left NULL in dX_list below; they are
  ## routed purely through the dSigma_e/d2Sigma_e FD channel.
  ## solution_derivatives internally distinguishes ok / !ok per param.
  ## ------------------------------------------------------------------
  sd1 <- if (length(struct_names) > 0) tryCatch(
    solution_derivatives(model, compiled, dr, params,
                         param_names = struct_names, obs_vars = obs_vars),
    error = function(e) NULL
  ) else NULL

  ## ------------------------------------------------------------------
  ## Second-order solution derivatives (structural blocks) for the
  ## STRUCTURAL subset only (same rationale as sd1 above).
  ## The "hvp_solution"/"adjoint_solution" T2 paths are d2X-FREE by
  ## construction, so they SKIP this O(np^2) second-order solve entirely
  ## (that is the point of both methods).
  ## ------------------------------------------------------------------
  d2x_free_methods <- c("hvp_solution", "adjoint_solution")
  sd2 <- if (!(t2_method %in% d2x_free_methods) && length(struct_names) > 0) tryCatch(
    solution_derivatives_2(model, compiled, dr, params,
                           param_names = struct_names, obs_vars = obs_vars),
    error = function(e) NULL
  ) else NULL

  ## ------------------------------------------------------------------
  ## dSigma_e and d2Sigma_e by FD of .get_shock_cov.
  ##
  ## Base values are resolved via .resolve_param_value() (falls back to the
  ## model's shocks-block stderr for a stderr-named param_names entry not
  ## present in `params`); the perturbed vectors are built by ASSIGNING the
  ## perturbed value under `nm` (which both updates an existing entry and
  ## APPENDS one for a not-yet-present stderr name), so .get_shock_cov's own
  ## Priority-0 "a parameter named exactly the shock" lookup picks it up
  ## exactly as it would for a normal MCMC draw.
  ## ------------------------------------------------------------------
  base_val <- vapply(param_names, .resolve_param_value, 0.0)
  names(base_val) <- param_names
  h_vec <- h_Sigma_e * pmax(abs(base_val), 1e-4)
  names(h_vec) <- param_names

  ## Single-step Sigma_e at theta +/- h_i (for dSigma_e and d2Sigma_e).
  Se_p <- Se_m <- vector("list", np); names(Se_p) <- names(Se_m) <- param_names
  for (i in seq_len(np)) {
    nm <- param_names[i]; h <- h_vec[i]
    pp <- params; pm <- params
    pp[nm] <- base_val[[nm]] + h; pm[nm] <- base_val[[nm]] - h
    Se_p[[nm]] <- .get_shock_cov(model, exo, pp)
    Se_m[[nm]] <- .get_shock_cov(model, exo, pm)
  }
  Se0 <- ss0$Sigma_e

  dSe <- lapply(param_names, function(nm) {
    h <- h_vec[nm]
    (Se_p[[nm]] - Se_m[[nm]]) / (2 * h)
  })
  names(dSe) <- param_names

  ## d2Sigma_e by a second FD:
  ##   diagonal (i=j): 3-point: (Se(+h) - 2*Se(0) + Se(-h)) / h^2
  ##   mixed   (i!=j): 4-corner: (Se(+hi,+hj) - Se(+hi,-hj) - Se(-hi,+hj) + Se(-hi,-hj)) / (4 hi hj)
  d2Se <- list()
  for (i in seq_len(np)) {
    nm_i <- param_names[i]; hi <- h_vec[nm_i]
    for (j in seq_len(np)) {
      nm_j <- param_names[j]; hj <- h_vec[nm_j]
      key <- paste(i, j, sep = "|")
      if (i == j) {
        d2Se[[key]] <- (Se_p[[nm_i]] - 2 * Se0 + Se_m[[nm_i]]) / (hi^2)
      } else {
        ## Compute Se at (theta +/- hi, +/- hj)
        pp <- params; pp[nm_i] <- base_val[[nm_i]] + hi; pp[nm_j] <- base_val[[nm_j]] + hj
        pm <- params; pm[nm_i] <- base_val[[nm_i]] + hi; pm[nm_j] <- base_val[[nm_j]] - hj
        mp <- params; mp[nm_i] <- base_val[[nm_i]] - hi; mp[nm_j] <- base_val[[nm_j]] + hj
        mm <- params; mm[nm_i] <- base_val[[nm_i]] - hi; mm[nm_j] <- base_val[[nm_j]] - hj
        Se_pp <- .get_shock_cov(model, exo, pp)
        Se_pm <- .get_shock_cov(model, exo, pm)
        Se_mp <- .get_shock_cov(model, exo, mp)
        Se_mm <- .get_shock_cov(model, exo, mm)
        d2Se[[key]] <- (Se_pp - Se_pm - Se_mp + Se_mm) / (4 * hi * hj)
      }
    }
  }

  ## ------------------------------------------------------------------
  ## Build dX_list: first-order blocks per parameter.
  ## ------------------------------------------------------------------
  dX_list <- vector("list", np); names(dX_list) <- param_names
  for (i in seq_len(np)) {
    nm <- param_names[i]
    d1 <- if (!is.null(sd1)) sd1$derivs[[nm]] else NULL
    ok <- !is.null(d1) && isTRUE(d1$ok)
    dX_list[[i]] <- list(
      dTT      = if (ok) d1$dTT      else NULL,
      dRR      = if (ok) d1$dRR      else NULL,
      dZZ      = if (ok) d1$dZZ      else NULL,
      dDD      = if (ok) d1$dDD      else NULL,
      dd       = if (ok) d1$dd       else NULL,
      dSigma_e = dSe[[nm]]
    )
  }

  ## ------------------------------------------------------------------
  ## Build d2X_list: second-order blocks per (i,j) pair.
  ## Skipped entirely for the d2X-free "hvp_solution" path (which forms T2
  ## without any materialised second-order solution block).
  ## ------------------------------------------------------------------
  d2X_list <- list()
  if (!(t2_method %in% d2x_free_methods)) {
    for (i in seq_len(np)) {
      for (j in seq_len(np)) {
        key  <- paste(i, j, sep = "|")
        sd2_blk <- if (!is.null(sd2)) sd2$d2[[key]] else NULL
        blk <- list(
          d2TT      = sd2_blk$d2TT,
          d2RR      = sd2_blk$d2RR,
          d2ZZ      = sd2_blk$d2ZZ,
          d2DD      = sd2_blk$d2DD,
          d2d       = sd2_blk$d2d,
          d2Sigma_e = d2Se[[key]]
        )
        d2X_list[[key]] <- blk
      }
    }
  }

  ## ------------------------------------------------------------------
  ## Loglik Hessian.
  ##   loop / contract_once : T1 + exact T2 = <G_X, d2X_ij> in one call.
  ##   hvp_solution         : T1 exact (zero d2X passed), then T2 formed d2X-free
  ##                          as the FD Jacobian of the exact solution-adjoint
  ##                          gradient (frozen filter bars) + the exact Sigma_e
  ##                          <G_Sig, d2Sigma_e> channel.
  ##   adjoint_solution     : T1 exact (zero d2X passed), then T2 formed d2X-free
  ##                          and FD-free (at the solution-adjoint layer) via
  ##                          the ANALYTIC Hessian-of-solution-adjoint
  ##                          (.solution_adjoint_hessian) + the exact Sigma_e
  ##                          <G_Sig, d2Sigma_e> channel.
  ## ------------------------------------------------------------------
  if (t2_method == "hvp_solution") {
    H <- .posterior_hessian_hvp_solution(
      model, compiled, params, param_names, struct_names, obs_vars, data, ss0,
      dX_list, d2Se, base_val, me_variance = me_variance, eps = eps,
      h_t2 = h_t2, t1_method = t1_method)
  } else if (t2_method == "adjoint_solution") {
    H <- .posterior_hessian_adjoint_solution(
      model, compiled, dr, params, param_names, struct_names, obs_vars, data, ss0,
      dX_list, d2Se, me_variance = me_variance, eps = eps,
      h_rel2 = h_rel2, t1_method = t1_method)
  } else {
    H <- kf_loglik_hessian(data, ss0, dX_list, d2X_list,
                           me_variance = me_variance, eps = eps,
                           t1_method = t1_method, t2_method = t2_method)
  }

  ## ------------------------------------------------------------------
  ## Optionally add the analytic prior Hessian (diagonal for independent
  ## priors; d2/dx2 log p(x) = d/dx [dlog_prior_density1]).
  ## We differentiate .dlog_prior_density1 numerically with a central
  ## difference (avoids reimplementing analytic second derivatives for
  ## every distribution).
  ## ------------------------------------------------------------------
  if (include_prior) {
    if (is.null(prior_spec))
      stop("posterior_hessian: prior_spec required when include_prior = TRUE")
    for (nm in param_names) {
      row_i <- match(nm, prior_spec$name)
      if (is.na(row_i)) next
      x0   <- base_val[[nm]]
      dist <- prior_spec$distribution[row_i]
      p1   <- prior_spec$p1[row_i]
      p2   <- prior_spec$p2[row_i]
      ## Check that this prior name maps to a param in param_names
      idx  <- match(nm, param_names)
      if (is.na(idx)) next
      ## d2 log p / dx2 by central FD of .dlog_prior_density1
      h_pr <- 1e-5 * max(abs(x0), 1e-4)
      d2pr <- (.dlog_prior_density1(x0 + h_pr, dist, p1, p2) -
               .dlog_prior_density1(x0 - h_pr, dist, p1, p2)) / (2 * h_pr)
      H[idx, idx] <- H[idx, idx] + d2pr
    }
  }

  if (isTRUE(check_mode) && !is.na(grad_norm_val)) {
    attr(H, "grad_norm")    <- grad_norm_val
    attr(H, "grad_at_mode") <- grad_at_mode
  }

  ## Which T2 second-primitive path solution_derivatives_2 used ("analytic"
  ## needs compile_model(param_deriv = "second"); "fd" is the stencil
  ## fallback; NA when sd2 itself failed and the FD-of-loglik fallback ran).
  ## "hvp_solution" is the d2X-free FD-of-solution-adjoint path (no sd2 at all).
  attr(H, "second_primitives") <-
    if (t2_method == "hvp_solution") "hvp_solution"
    else if (t2_method == "adjoint_solution") "adjoint_solution"
    else if (!is.null(sd2)) sd2$second_primitives else NA_character_
  H
}


## --------------------------------------------------------------------------
## .posterior_hessian_hvp_solution(): d2X-FREE loglik Hessian.
##
## T1 (filter curvature) is computed exactly by kf_loglik_hessian with a ZERO
## second-order block list (so its T2 term vanishes and it returns T1 alone).
##
## T2 (solution curvature) = <G_X, d2X_ij> is formed WITHOUT any materialised
## d2X block, using the identity
##     T2[i,j] = d^2/dtheta_i dtheta_j  phi(theta),   phi(theta) = <G_X, X(theta)>
## with G_X FROZEN at the base point. Then grad_theta phi(theta) = <G_X, dX/dtheta>
## is EXACTLY one .solution_adjoint call whose `bars` are the frozen filter
## gradients (G_TT,G_RR,G_ZZ,G_DD,g_d) -- the solution channel of the loglik
## gradient. Central-differencing that exact analytic gradient in theta_i gives
## row i of the solution part of T2 (struct params only; the solution map does
## not depend on shock-std params -- certainty equivalence). The Sigma_e
## curvature <G_Sig, d2Sigma_e_ij> is added exactly from the same d2Se stencil
## the loop/contract_once paths use (Sigma_e never enters .solution_adjoint).
##
## Requires every perturbed point (theta_i +/- h) to solve (SS converged, BK
## satisfied); a failed struct direction warns and contributes 0 to that row
## (documented limitation of the FD path).
## --------------------------------------------------------------------------
#' @noRd
.posterior_hessian_hvp_solution <- function(model, compiled, params,
                                            param_names, struct_names, obs_vars,
                                            Y, ss0, dX_list, d2Se, base_val,
                                            me_variance = 0, eps = 1e-5,
                                            h_t2 = 1e-5,
                                            t1_method = c("hvp", "analytic")) {
  t1_method <- match.arg(t1_method)
  np <- length(param_names)

  ## ---- T1 exact: kf_loglik_hessian with an empty d2X list (its T2 == 0) ----
  H_t1 <- kf_loglik_hessian(Y, ss0, dX_list, list(),
                            me_variance = me_variance, eps = eps,
                            t1_method = t1_method, t2_method = "contract_once")
  t1_asym <- attr(H_t1, "t1_asymmetry")
  T1 <- matrix(as.numeric(H_t1), np, np)

  ## ---- base-point frozen filter bars G_X ----------------------------------
  bars <- .kf_loglik_adjoint(Y, ss0, list(), me_variance = me_variance,
                             return_bars = TRUE)$bars
  bars_solution <- list(G_TT = bars$G_TT, G_RR = bars$G_RR,
                        G_ZZ = bars$G_ZZ, G_DD = bars$G_DD, g_d = bars$g_d)
  G_Sig <- bars$G_Sig

  ## ---- solution-adjoint gradient at a perturbed parameter point -----------
  ## Returns the (frozen-G_X) solution-channel gradient over `struct_names`, or
  ## NULL if the re-solve fails.
  sol_grad_at <- function(params_pt) {
    ss_pt <- tryCatch(
      solve_steady_state(model, compiled, params_pt, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(ss_pt) || !isTRUE(ss_pt$converged)) return(NULL)
    dr_pt <- tryCatch(
      solve_perturbation(model, compiled, ss_pt$ss, params_pt, order = 1L,
                         verbose = FALSE),
      error = function(e) NULL)
    if (is.null(dr_pt) || !isTRUE(dr_pt$bk_satisfied)) return(NULL)
    res <- tryCatch(
      .solution_adjoint(model, compiled, dr_pt, params_pt, struct_names,
                        obs_vars = obs_vars, bars = bars_solution),
      error = function(e) NULL)
    if (is.null(res)) return(NULL)
    g <- res$grad
    g[!res$ok] <- NA_real_
    g
  }

  ## ---- T2 solution part: central FD of sol_grad_at over struct params -----
  n_s <- length(struct_names)
  T2_sol <- matrix(0, n_s, n_s)
  if (n_s > 0) {
    for (i in seq_len(n_s)) {
      nm_i <- struct_names[i]
      hi   <- h_t2 * max(abs(base_val[[nm_i]]), 1e-4)
      pp <- params; pp[[nm_i]] <- base_val[[nm_i]] + hi
      pm <- params; pm[[nm_i]] <- base_val[[nm_i]] - hi
      gp <- sol_grad_at(pp)
      gm <- sol_grad_at(pm)
      if (is.null(gp) || is.null(gm)) {
        warning(sprintf(
          "posterior_hessian(hvp_solution): re-solve failed at %s +/- h; ",
          nm_i), "row set to 0 (T2 solution part).", call. = FALSE)
        next
      }
      row_i <- (gp[struct_names] - gm[struct_names]) / (2 * hi)
      row_i[!is.finite(row_i)] <- 0
      T2_sol[i, ] <- row_i
    }
  }

  ## ---- assemble full np x np T2 -------------------------------------------
  T2 <- matrix(0, np, np)
  if (n_s > 0) {
    sp <- match(struct_names, param_names)
    T2[sp, sp] <- T2_sol
  }
  ## Sigma_e second-order channel: exact <G_Sig, d2Sigma_e_ij> for ALL pairs.
  if (!is.null(G_Sig)) {
    for (i in seq_len(np)) for (j in seq_len(np)) {
      d2 <- d2Se[[paste(i, j, sep = "|")]]
      if (!is.null(d2)) T2[i, j] <- T2[i, j] + sum(G_Sig * d2)
    }
  }

  H <- T1 + T2
  H <- 0.5 * (H + t(H))
  dimnames(H) <- list(param_names, param_names)
  attr(H, "t1_asymmetry") <- t1_asym
  H
}


## --------------------------------------------------------------------------
## .posterior_hessian_adjoint_solution(): d2X-FREE loglik Hessian, FD-free at
## the solution-adjoint layer.
##
## T1 (filter curvature) is computed exactly by kf_loglik_hessian, same as the
## "hvp_solution" path (a ZERO second-order block list, so its own T2 vanishes
## and it returns T1 alone).
##
## T2 (solution curvature) = <G_X, d2X_ij> is formed in ONE call to
## \code{.solution_adjoint_hessian} (gradient-solution-adjoint-order2-param.R):
## the ANALYTIC Hessian of phi(theta) = <G_X, X(theta)> with G_X FROZEN at the
## base point -- no per-parameter model re-solve (unlike "hvp_solution", which
## needs 2 . np_struct re-solves), no materialised d2X blocks. The Sigma_e
## curvature <G_Sig, d2Sigma_e_ij> is added exactly from the same d2Se stencil
## the loop/contract_once/hvp_solution paths use (Sigma_e never enters the
## solution-adjoint chain).
## --------------------------------------------------------------------------
#' @noRd
.posterior_hessian_adjoint_solution <- function(model, compiled, dr, params,
                                                param_names, struct_names,
                                                obs_vars, Y, ss0, dX_list, d2Se,
                                                me_variance = 0, eps = 1e-5,
                                                h_rel2 = 1e-4,
                                                t1_method = c("hvp", "analytic")) {
  t1_method <- match.arg(t1_method)
  np <- length(param_names)

  ## ---- T1 exact: kf_loglik_hessian with an empty d2X list (its T2 == 0) ----
  H_t1 <- kf_loglik_hessian(Y, ss0, dX_list, list(),
                            me_variance = me_variance, eps = eps,
                            t1_method = t1_method, t2_method = "contract_once")
  t1_asym <- attr(H_t1, "t1_asymmetry")
  T1 <- matrix(as.numeric(H_t1), np, np)

  ## ---- base-point frozen filter bars G_X ----------------------------------
  bars <- .kf_loglik_adjoint(Y, ss0, list(), me_variance = me_variance,
                             return_bars = TRUE)$bars
  bars_solution <- list(G_TT = bars$G_TT, G_RR = bars$G_RR,
                        G_ZZ = bars$G_ZZ, G_DD = bars$G_DD, g_d = bars$g_d)
  G_Sig <- bars$G_Sig

  ## ---- T2 solution part: ONE analytic Hessian-of-solution-adjoint call ----
  n_s <- length(struct_names)
  T2_sol <- matrix(0, n_s, n_s)
  if (n_s > 0) {
    res <- tryCatch(
      .solution_adjoint_hessian(model, compiled, dr, params, struct_names,
                                obs_vars = obs_vars, bars = bars_solution,
                                h_rel2 = h_rel2),
      error = function(e) {
        warning("posterior_hessian(adjoint_solution): solution-Hessian ",
                "failed (", conditionMessage(e), "); T2 solution part set ",
                "to 0.", call. = FALSE)
        NULL
      })
    if (!is.null(res)) {
      T2_sol <- res$T2
      T2_sol[!is.finite(T2_sol)] <- 0
      if (any(!res$ok))
        warning(sprintf(
          "posterior_hessian(adjoint_solution): primitives failed for %s; ",
          paste(struct_names[!res$ok], collapse = ", ")),
          "corresponding rows/cols set to 0 (T2 solution part).",
          call. = FALSE)
    }
  }

  ## ---- assemble full np x np T2 -------------------------------------------
  T2 <- matrix(0, np, np)
  if (n_s > 0) {
    sp <- match(struct_names, param_names)
    T2[sp, sp] <- T2_sol
  }
  ## Sigma_e second-order channel: exact <G_Sig, d2Sigma_e_ij> for ALL pairs.
  if (!is.null(G_Sig)) {
    for (i in seq_len(np)) for (j in seq_len(np)) {
      d2 <- d2Se[[paste(i, j, sep = "|")]]
      if (!is.null(d2)) T2[i, j] <- T2[i, j] + sum(G_Sig * d2)
    }
  }

  H <- T1 + T2
  H <- 0.5 * (H + t(H))
  dimnames(H) <- list(param_names, param_names)
  attr(H, "t1_asymmetry") <- t1_asym
  H
}


## --------------------------------------------------------------------------
## posterior_hessian_fd_grad(): INDEPENDENT cross-check Hessian (paper gap #7,
## sharpened).
##
## posterior_hessian()'s own internal check (t1_asymmetry, and the
## t1_method = "hvp" vs "analytic" agreement) certifies only the T1
## filter-curvature term: T1 is a directional finite difference (or, for
## "analytic", a forward-over-reverse adjoint) of the SAME exact adjoint
## gradient that also appears in T2 = <G_X, d2X_ij>. Because T1 and T2 SHARE
## the base-point adjoint gradient machinery, an error confined to T2 (e.g. a
## wrong second-order solution-derivative block from solution_derivatives_2)
## is invisible to that internal check -- it would corrupt the Hessian's
## flattest directions (the ones T2 dominates) while T1 continues to look
## perfectly consistent with itself.
##
## The genuinely independent reference is a central finite difference OF THE
## EXACT ADJOINT GRADIENT ITSELF (make_posterior_grad), built via a completely
## different code path (repeated model re-solves + full forward likelihood
## evaluations, no solution_derivatives_2 involved at all) and therefore with
## a different error structure: it cannot share a T2 bug with
## posterior_hessian, because it never computes a T2 term.
## --------------------------------------------------------------------------

#' Independent cross-check Hessian: central FD of the exact adjoint gradient
#'
#' Builds the analytic posterior-gradient closure ONCE (via
#' \code{make_posterior_grad}) and forms the Hessian by a central finite
#' difference of that EXACT gradient (2 * length(param_names) gradient
#' evaluations, each of which re-solves the model and re-runs the filter --
#' a completely different code path from \code{posterior_hessian}'s
#' solution_derivatives_2-based T2 term).
#'
#' This is the VALIDATION reference for \code{posterior_hessian}'s T2 term,
#' not a replacement for it: \code{posterior_hessian}'s own internal
#' consistency check (\code{t1_asymmetry}, and \code{t1_method = "hvp"} vs
#' \code{"analytic"} agreement) certifies only T1, because T1 and T2 share the
#' same base-point adjoint gradient (\code{.kf_loglik_adjoint}) -- a T2-only
#' bug (e.g. in \code{solution_derivatives_2}) would pass that internal check
#' while still being wrong in the Hessian's flattest directions (typically the
#' directions T2 dominates). Central-differencing the exact gradient uses no
#' second-order solution-derivative code at all, so it cannot share a T2 bug
#' with \code{posterior_hessian}.
#'
#' Recommended use: compare eigen-decompositions of the two Hessians,
#' especially the SIGNS of the smallest-magnitude (flattest) eigenvalues --
#' agreement on the largest eigenvalues is necessary but not sufficient
#' (T2 errors tend to be small in magnitude and hence concentrated in/near the
#' null space), so a flipped sign in a flat direction is the diagnostic this
#' cross-check exists to catch.
#'
#' @param model     parsed model (from \code{parse_mod}).
#' @param data      n_obs x T observation matrix (no NAs); passed through to
#'                  \code{make_posterior_grad} as \code{data}.
#' @param prior_spec prior-spec data.frame (defines the FULL parameter vector
#'                  ordering that \code{make_posterior_grad}'s closure expects;
#'                  \code{param_names} may be any subset of \code{prior_spec$name}).
#' @param obs_vars  character vector of observable variable names.
#' @param compiled  compiled model (from \code{compile_model}).
#' @param theta     named numeric vector, the evaluation point, covering (at
#'                  least) every name in \code{prior_spec$name} (extra names
#'                  are ignored; missing names fall back to \code{prior_spec$mean}).
#' @param param_names character vector of parameter names to differentiate
#'                  (default \code{names(theta)}); the returned Hessian is
#'                  \code{length(param_names) x length(param_names)}, a
#'                  sub-block of the full \code{prior_spec}-dimensional
#'                  central-difference Hessian.
#' @param h_rel     relative FD step (default 1e-5): the step for parameter
#'                  \code{k} is \code{h_rel * max(abs(theta[k]), 1e-3)}.
#' @param me_variance scalar measurement-error variance (default 0); passed
#'                  through to \code{make_posterior_grad}.
#' @param ...       additional arguments forwarded to \code{make_posterior_grad}
#'                  (e.g. \code{grad_method}, \code{likelihood}).
#' @return \code{length(param_names) x length(param_names)} symmetric matrix:
#'   the central-difference Hessian of the exact adjoint gradient, symmetrised
#'   via \code{(H + t(H))/2} (the raw central difference is not exactly
#'   symmetric at finite \code{h_rel}; the asymmetry is itself an FD-error
#'   diagnostic and is attached as attribute \code{fd_asymmetry}).
#' @seealso \code{posterior_hessian}
#' @export
posterior_hessian_fd_grad <- function(model, data, prior_spec, obs_vars,
                                      compiled, theta,
                                      param_names = names(theta),
                                      h_rel = 1e-5,
                                      me_variance = 0,
                                      ...) {
  full_names <- prior_spec$name
  ## Build the FULL theta vector the gradient closure expects: start from
  ## prior_spec$mean, override with anything supplied in `theta`.
  theta_full <- setNames(prior_spec$mean, full_names)
  common <- intersect(names(theta), full_names)
  theta_full[common] <- theta[common]

  if (!all(param_names %in% full_names))
    stop("posterior_hessian_fd_grad: param_names must be a subset of ",
         "prior_spec$name; missing: ",
         paste(setdiff(param_names, full_names), collapse = ", "))

  grad_fn <- make_posterior_grad(model, data, prior_spec, obs_vars, compiled,
                                 me_variance = me_variance, ...)

  np <- length(param_names)
  h_vec <- h_rel * pmax(abs(theta_full[param_names]), 1e-3)
  names(h_vec) <- param_names

  H <- matrix(0, np, np, dimnames = list(param_names, param_names))
  for (k in seq_len(np)) {
    nm <- param_names[k]
    h  <- h_vec[[nm]]
    tp <- theta_full; tp[[nm]] <- theta_full[[nm]] + h
    tm <- theta_full; tm[[nm]] <- theta_full[[nm]] - h

    gp <- grad_fn(tp)[param_names]
    gm <- grad_fn(tm)[param_names]

    H[k, ] <- (gp - gm) / (2 * h)
  }

  attr(H, "fd_asymmetry") <- max(abs(H - t(H)))
  H <- 0.5 * (H + t(H))
  dimnames(H) <- list(param_names, param_names)
  H
}


#' Laplace approximation to the log marginal likelihood
#'
#' Computes the Laplace (second-order) approximation to the model log marginal
#' likelihood from a mode-finding result that carries an exact posterior Hessian
#' (\code{run_mode_finding(use_exact_hessian = TRUE)}):
#' \deqn{\log p(Y\mid M) \approx \log p(Y\mid\hat\theta) + \log p(\hat\theta\mid M)
#'        + \tfrac{d}{2}\log 2\pi - \tfrac12 \log\det(-H),}
#' where \eqn{H} is the exact Hessian of the log-posterior at the mode and the
#' first two terms are the (unnormalised) log-posterior at the mode
#' (\code{mode_result$mode$logpost}). Deterministic and exact to
#' \eqn{O(1/T)} (Tierney & Kadane, 1986); a fast, randomness-free alternative to
#' the SMC marginal-likelihood estimator for model comparison.
#'
#' @param mode_result A \code{dynhr_mode_result} from
#'   \code{\link{run_mode_finding}} with \code{hessian_exact} populated.
#' @return The Laplace log marginal likelihood (numeric scalar), or
#'   \code{NA_real_} when no exact Hessian is present or \eqn{-H} is not
#'   positive-definite at the mode (the mode is not a strict local maximum).
#' @seealso \code{posterior_hessian}, \code{\link{run_mode_finding}}
#' @export
laplace_log_marglik <- function(mode_result) {
  H <- mode_result$hessian_exact
  if (is.null(H)) {
    warning("laplace_log_marglik: mode_result has no exact Hessian; re-run ",
            "run_mode_finding(use_exact_hessian = TRUE).", call. = FALSE)
    return(NA_real_)
  }
  ev <- eigen(-H, symmetric = TRUE)$values
  if (any(ev <= 0)) {
    warning("laplace_log_marglik: -H is not positive-definite at the mode; ",
            "the Laplace marginal likelihood is undefined.", call. = FALSE)
    return(NA_real_)
  }
  d <- length(mode_result$theta_mode)
  as.numeric(mode_result$mode$logpost + 0.5 * d * log(2 * pi) - 0.5 * sum(log(ev)))
}
