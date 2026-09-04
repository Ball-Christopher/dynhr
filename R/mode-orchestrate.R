## R/mode-orchestrate.R
## --------------------------------------------------------------------------
## Phase-2 split from estimate-monolith.R.
##
## combined_optimize()  -- multi-stage mode finder (CMA-ES -> L-BFGS-B etc.)
## .run_mode_finding()  -- orchestrator: builds bounds, wraps neg-logpost,
##                         dispatches to the right optimizer by name.
## --------------------------------------------------------------------------

#' Multi-stage mode finder
#'
#' Runs a sequence of optimizers, passing the best parameter vector from each
#' stage into the next.  The iteration budget `max_iter` is split across
#' stages according to `split`.
#'
#' @param fn         Objective to MINIMISE (return scalar numeric)
#' @param par        Named starting vector
#' @param lower,upper  Box bounds (scalar or vector)
#' @param max_iter   Total iteration budget
#' @param stages     Character vector of stage names: "cmaes", "nmkb",
#'                   "jade", "lbfgsb", "nm", "newrat"
#' @param split      Fractional budget allocation (will be normalised to sum 1)
#' @param verbose    Print progress messages
#' @param progress   Show cli progress bars inside each stage
#' @param gr         Optional analytic gradient of fn (for lbfgsb and newrat stages)
#' @param newrat_H0_inv  Optional initial inverse-Hessian for the newrat stage
#'   (n x n matrix). When NULL the newrat stage uses csminwel's default 1e-4*I.
#' @param record_curvature  Opt-in (default \code{FALSE}): stash the final
#'   BFGS inverse-Hessian returned by the newrat (csminwel) stage in the
#'   result as \code{$H_bfgs} (\code{NULL} if no "newrat" stage ran). This is
#'   the H0 seed after being updated by every (s_k, y_k) curvature pair
#'   csminwel collected along its own optimisation trajectory -- free
#'   (already computed internally by csminwel), just not normally propagated
#'   upward. Off by default to keep the return shape unchanged for existing
#'   callers.
#' @return list(par, value, convergence, iterations, message), plus
#'   \code{$H_bfgs} when \code{record_curvature = TRUE}.
#' @noRd
combined_optimize <- function(fn, par, lower = -Inf, upper = Inf,
                              max_iter = 10000,
                              stages = c("cmaes", "lbfgsb"),
                              split = c(0.8, 0.2),
                              verbose = TRUE, progress = TRUE,
                              gr = NULL,
                              newrat_H0_inv = NULL,
                              record_curvature = FALSE) {
  n <- length(par)
  par_names <- names(par)
  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)

  split      <- split / sum(split)
  iter_alloc <- pmax(round(max_iter * split), 1)

  best_par <- par
  best_val <- fn(par)
  H_bfgs_out <- NULL   # populated only by a "newrat" stage when recorded

  for (i in seq_along(stages)) {
    stage <- stages[i]
    iters <- iter_alloc[i]

    if (verbose) cat(sprintf("\n--- Stage %d: %s (%d iter) ---\n",
                             i, toupper(stage), iters))

    res <- switch(stage,
      "cmaes" = cmaes_optimize(fn, best_par, lower, upper,
                               max_iter = iters, verbose = verbose,
                               progress = progress),
      "nmkb"  = nmkb_optimize(fn, best_par, lower, upper,
                               max_iter = iters, verbose = verbose),
      "jade"  = jade_optimize(fn, best_par, lower, upper,
                               max_iter = iters, verbose = verbose,
                               progress = progress),
      "lbfgsb" = {
        fn_named <- function(x) {
          names(x) <- par_names
          val <- fn(x)
          if (!is.finite(val)) 1e20 else val
        }
        ## Analytic gradient (when supplied): L-BFGS-B then needs ONE gradient
        ## eval per step instead of FD's (n+1) objective evals -- a large win on
        ## models with many parameters (P2). `gr` returns d(fn)/dx, i.e. the
        ## gradient of the MINIMISATION objective (already negated by the caller).
        ## A non-finite gradient (off-domain proposal) is replaced by zeros so
        ## L-BFGS-B backs off via the line search rather than erroring.
        gr_named <- if (!is.null(gr)) function(x) {
          names(x) <- par_names
          g <- gr(x)
          if (anyNA(g) || any(!is.finite(g))) rep(0, length(x)) else g
        } else NULL
        optres <- optim(best_par, fn_named, gr = gr_named, method = "L-BFGS-B",
                        lower = lower, upper = upper,
                        control = list(maxit = iters, factr = 1e7))
        list(par = setNames(optres$par, par_names),
             value = optres$value,
             convergence = optres$convergence,
             iterations = iters,
             message = optres$message %||% "done")
      },
      "nm" = {
        fn_named <- function(x) {
          names(x) <- par_names
          val <- fn(x)
          if (!is.finite(val)) 1e20 else val
        }
        optres <- optim(best_par, fn_named, method = "Nelder-Mead",
                        control = list(maxit = iters))
        list(par = setNames(optres$par, par_names),
             value = optres$value,
             convergence = optres$convergence,
             iterations = iters,
             message = optres$message %||% "done")
      },
      "newrat" = {
        ## csminwel (Sims's quasi-Newton) with the infeasible-backtracking line
        ## search. Maps +Inf -> +1e20 so csminit can backtrack on BK-violating
        ## neighbours; does NOT zero the gradient at infeasible points (unlike the
        ## L-BFGS-B path) -- csminwel's line search handles them via backtracking.
        fn_cs <- function(x) {
          names(x) <- par_names
          val <- fn(x)
          if (!is.finite(val)) 1e20 else val
        }
        ## When an analytic gradient is available, wrap it to handle
        ## infeasibility: return the raw gradient (finite or not) -- csminwel's
        ## get_grad wrapper marks badg=TRUE on non-finite components and triggers
        ## Hessian reset / retry, which is exactly the desired behaviour near
        ## infeasible walls. Do NOT zero the gradient here.
        gr_cs <- if (!is.null(gr)) {
          function(x) {
            names(x) <- par_names
            g <- gr(x)
            ## csminwel can handle badg gracefully; return raw gradient.
            as.numeric(g)
          }
        } else NULL
        optres <- csminwel(fn_cs, best_par, H0 = newrat_H0_inv,
                           grad = gr_cs, nit = iters, verbose = verbose)
        if (isTRUE(record_curvature)) H_bfgs_out <- optres$H
        list(par = setNames(optres$xh, par_names),
             value = optres$fh,
             convergence = optres$convergence,
             iterations = optres$itct,
             message = sprintf("csminwel retcode=%d", optres$retcode))
      }
    )

    if (res$value < best_val) {
      best_val <- res$value
      best_par <- res$par
    }

    if (verbose) cat(sprintf("  -> logpost = %.4f (%s)\n", -best_val, res$message))
  }

  out <- list(par = best_par, value = best_val,
             convergence = 0, iterations = sum(iter_alloc),
             message = sprintf("combined(%s)", paste(stages, collapse = "->")))
  if (isTRUE(record_curvature)) out$H_bfgs <- H_bfgs_out
  out
}


#' Orchestrate mode-finding for MCMC initialisation
#'
#' Builds parameter bounds from `prior_spec`, wraps the log-posterior as
#' a minimisation objective, and dispatches to the requested optimizer.
#'
#' @param log_post_fn  Function(theta) -> list(logpost, ...) or scalar
#' @param theta_init   Named starting vector
#' @param prior_spec   Prior spec data.frame (from extract_prior_spec())
#' @param nm_maxit     Total iteration budget passed to the optimizer
#' @param lbfgsb_maxit (ignored; retained for back-compat)
#' @param method       Optimizer name: "combined", "cmaes", "nmkb", "jade",
#'                     "nelder", "cmaes_nmkb", "cmaes_jade"
#' @param transform    Optional "dynhr_param_transform" object (from
#'   \code{\link{build_param_transform}}). When non-NULL (opt-in), mode
#'   finding is performed in UNCONSTRAINED eta-space:
#'   \itemize{
#'     \item `eta_init = transform$to_unconstrained(theta_init)`.
#'     \item the objective minimised is `-log_post_fn(to_constrained(eta))`,
#'       i.e. WITHOUT the change-of-variables (Jacobian) term. This preserves
#'       the same argmax as the untransformed theta-space posterior (Stan's
#'       "optimizing" convention drops the Jacobian for MAP/mode estimates),
#'       while removing the need for box constraints / boundary clamping --
#'       eta ranges over the whole real line for every parameter.
#'     \item box bounds for the bounded optimizer stages (cmaes/nmkb) are
#'       built in eta-space (see below) instead of theta-space.
#'     \item after optimisation, `theta_mode = to_constrained(eta_hat)` and
#'       `logpost` is the (Jacobian-free) theta-space log-posterior at
#'       `theta_mode`, exactly as in the untransformed path.
#'   }
#'   When NULL (default), behaviour is bit-identical to before.
#' @param verbose      Print progress messages
#' @param record_curvature  Opt-in (default \code{FALSE}); forwarded to
#'   \code{\link{combined_optimize}} -- see there. When \code{TRUE} and a
#'   "newrat" stage ran, the returned list carries \code{$H_bfgs} (the final
#'   BFGS inverse-Hessian in WORKING space, i.e. eta-space when
#'   \code{transform} is supplied).
#' @return list(theta_mode, logpost, convergence, iterations, method[, H_bfgs])
#' @noRd
.run_mode_finding <- function(log_post_fn, theta_init, prior_spec,
                              nm_maxit = 10000, lbfgsb_maxit = NULL,
                              method = "newrat", transform = NULL,
                              grad_fn = NULL,
                              hessian_fn = NULL,
                              verbose = TRUE,
                              record_curvature = FALSE) {

  n <- length(theta_init)
  par_names <- names(theta_init)

  if (is.null(transform)) {
    # ---- Untransformed (default) path: theta-space, box-constrained -------
    lower <- setNames(rep(-Inf, n), par_names)
    upper <- setNames(rep( Inf, n), par_names)
    for (i in seq_len(nrow(prior_spec))) {
      nm <- prior_spec$name[i]
      if (nm %in% par_names) {
        if (!is.na(prior_spec$lower[i])) lower[nm] <- prior_spec$lower[i]
        if (!is.na(prior_spec$upper[i])) upper[nm] <- prior_spec$upper[i]
      }
    }
    eps   <- 1e-8
    lower <- ifelse(is.finite(lower), lower + eps, lower)
    upper <- ifelse(is.finite(upper), upper - eps, upper)

    neg_lp <- function(theta) {
      names(theta) <- par_names
      res <- log_post_fn(theta)
      if (is.list(res)) -res$logpost else -res
    }

    par_init <- theta_init
  } else {
    # ---- Transformed (opt-in) path: unconstrained eta-space ----------------
    # eta_init = to_unconstrained(theta_init); the objective is the
    # Jacobian-free theta-space negative log-posterior evaluated at
    # to_constrained(eta) -- same argmax as theta-space, but unconstrained.
    eta_init <- transform$to_unconstrained(theta_init)

    neg_lp <- function(eta) {
      names(eta) <- par_names
      theta <- transform$to_constrained(eta)
      res <- log_post_fn(theta)
      if (is.list(res)) -res$logpost else -res
    }

    # Bounds for the box-constrained stages (cmaes / nmkb), built in
    # eta-space. Map each FINITE theta-bound through to_unconstrained;
    # replace +-Inf eta-bounds with eta_init_j +- 30 -- wide enough to be
    # effectively unconstrained for any optimizer step, but finite so
    # cmaes_optimize()/nmkb_optimize() (which require finite bounds for
    # their default sigma0 / simplex sizing) behave sensibly. No 1e-8
    # clamping is applied in eta-space: eta = +-Inf already maps to the open
    # interval boundary in theta-space via to_constrained(), so there is no
    # boundary to clamp away from.
    lower <- setNames(rep(NA_real_, n), par_names)
    upper <- setNames(rep(NA_real_, n), par_names)
    for (j in seq_len(n)) {
      nm <- par_names[j]
      a  <- transform$a[[nm]]
      b  <- transform$b[[nm]]
      # to_unconstrained() is vectorized over ALL par_names at once; build a
      # full theta vector (theta_init elsewhere, the bound value at `nm`)
      # and pick out the j-th eta component.
      if (is.finite(a)) {
        theta_a <- theta_init
        theta_a[[nm]] <- a
        lower[nm] <- transform$to_unconstrained(theta_a)[[nm]]
      } else {
        lower[nm] <- eta_init[[nm]] - 30
      }
      if (is.finite(b)) {
        theta_b <- theta_init
        theta_b[[nm]] <- b
        upper[nm] <- transform$to_unconstrained(theta_b)[[nm]]
      } else {
        upper[nm] <- eta_init[[nm]] + 30
      }
    }

    par_init <- eta_init
  }

  lp0 <- -neg_lp(par_init)
  if (verbose) cat(sprintf("  Initial log-posterior: %.4f\n", lp0))

  ## Negated-objective gradient for the L-BFGS-B stage (P2).
  ## grad_fn returns d(logpost)/d(theta) (theta-space, no Jacobian term).
  ## For the untransformed path the negated minimisation gradient is simply
  ## -grad_fn(theta).
  ## For the transformed (eta-space) path the chain rule gives:
  ##   d(-logpost)/d(eta_j) = -dtheta_j/deta_j * grad_fn(theta)_j
  ## where dtheta_deta(eta) is the diagonal of the Jacobian (separable
  ## transform, so off-diagonal terms are zero). No log-Jacobian correction
  ## is applied: mode-finding uses the Jacobian-free objective (pure MAP,
  ## same argmax in theta-space; the volume correction is irrelevant here).
  neg_gr <- NULL
  if (!is.null(grad_fn)) {
    if (is.null(transform)) {
      neg_gr <- function(theta) {
        names(theta) <- par_names
        -as.numeric(grad_fn(theta))
      }
    } else {
      neg_gr <- function(eta) {
        names(eta)   <- par_names
        theta        <- transform$to_constrained(eta)
        names(theta) <- par_names
        g_theta      <- as.numeric(grad_fn(theta))
        -transform$dtheta_deta(eta) * g_theta
      }
    }
  }

  ## -------------------------------------------------------------------------
  ## Wrap the theta-space analytic Hessian for the working (possibly eta) space.
  ## hessian_fn(theta) -> n x n neg-logpost Hessian (positive semi-def at mode).
  ## For the transform path we apply the diagonal Jacobian chain rule:
  ##   H_eta[i,j] = J[i] * H_theta[i,j] * J[j]
  ## where J = dtheta/deta (diagonal of the transform Jacobian).
  ## This is the exact Hessian of the Jacobian-free objective in eta-space
  ## (ignoring the log-Jacobian correction that .run_mode_finding already
  ## drops for MAP/mode finding).
  hessian_fn_working <- NULL
  if (!is.null(hessian_fn)) {
    if (is.null(transform)) {
      hessian_fn_working <- function(x) {
        names(x) <- par_names
        hessian_fn(x)
      }
    } else {
      hessian_fn_working <- function(x) {
        names(x)     <- par_names
        theta        <- transform$to_constrained(x)
        names(theta) <- par_names
        H_theta      <- hessian_fn(theta)
        J <- transform$dtheta_deta(x)   # diagonal: d(theta_j)/d(eta_j)
        ## Sandwich: H_eta = diag(J) %*% H_theta %*% diag(J)
        J * H_theta * rep(J, each = length(J))  # outer product broadcasting
      }
    }
  }

  total_iter <- nm_maxit

  ## -------------------------------------------------------------------------
  ## newrat H0 (initial inverse-Hessian for csminwel).
  ## When the method includes a newrat stage AND an analytic Hessian supplier is
  ## available (passed as `hessian_fn` -- a function(par_in_working_space) ->
  ## n x n matrix of the negative log-posterior Hessian), build H0_inv from it
  ## at the current start point with .make_pd regularisation then invert.
  ## Falls back to NULL (csminwel default 1e-4 * I) on any failure.
  ##
  ## NB: csminwel's H0 parameter is the INVERSE-Hessian (used directly as the
  ## quasi-Newton scaling matrix H in dx = -H g). Confirm: line 106 of
  ## csminwel.R: dx <- as.numeric(-H0 %*% g), where H0 is always the
  ## inverse-Hessian approximation. So we pass the INVERSE of neg-logpost
  ## Hessian, i.e. the covariance matrix at the mode.
  ## -------------------------------------------------------------------------
  newrat_uses_csminwel <- method %in% c("newrat", "cmaes_newrat")
  newrat_H0_inv <- NULL
  if (newrat_uses_csminwel && !is.null(hessian_fn_working)) {
    newrat_H0_inv <- tryCatch({
      H_neg <- hessian_fn_working(par_init)
      ## H_neg should be the neg-logpost Hessian (positive definite at mode).
      ## .make_pd needs the neg-Hessian (pos-definite); invert to get H0_inv.
      V <- .make_pd(H_neg, cond_target = 100)
      V
    }, error = function(e) {
      if (verbose)
        cat(sprintf("  (newrat H0 build failed: %s; using csminwel default 1e-4*I)\n",
                    conditionMessage(e)))
      NULL
    })
    if (!is.null(newrat_H0_inv) && verbose)
      cat("  Built initial inverse-Hessian for newrat from analytic posterior Hessian.\n")
  }

  res <- switch(method,
    "cmaes_nmkb" = combined_optimize(neg_lp, par_init, lower, upper,
                                     max_iter = total_iter,
                                     stages = c("cmaes", "nmkb"),
                                     split = c(0.8, 0.2),
                                     verbose = verbose),
    "nelder"     = combined_optimize(neg_lp, par_init, lower, upper,
                                     max_iter = total_iter,
                                     stages = c("nm", "lbfgsb"),
                                     split = c(0.8, 0.2),
                                     verbose = verbose, gr = neg_gr),
    "cmaes"      = cmaes_optimize(neg_lp, par_init, lower, upper,
                                  max_iter = total_iter, verbose = verbose),
    "nmkb"       = nmkb_optimize(neg_lp, par_init, lower, upper,
                                 max_iter = total_iter, verbose = verbose),
    "jade"       = jade_optimize(neg_lp, par_init, lower, upper,
                                 max_iter = total_iter, verbose = verbose),
    "combined"   = combined_optimize(neg_lp, par_init, lower, upper,
                                     max_iter = total_iter,
                                     stages = c("cmaes", "lbfgsb"),
                                     split = c(0.8, 0.2),
                                     verbose = verbose, gr = neg_gr),
    "cmaes_jade" = combined_optimize(neg_lp, par_init, lower, upper,
                                     max_iter = total_iter,
                                     stages = c("cmaes", "jade", "lbfgsb"),
                                     split = c(0.5, 0.3, 0.2),
                                     verbose = verbose, gr = neg_gr),
    ## newrat: csminwel-only (no global search warm-start). Suited for
    ## problems where a reasonable starting point is available (e.g. prior
    ## means for a well-scaled model, or after CMA-ES). The infeasible-
    ## backtracking line search in csminwel handles BK-violating neighbours
    ## without stalling.
    "newrat"     = combined_optimize(neg_lp, par_init, lower, upper,
                                     max_iter = total_iter,
                                     stages = c("newrat"),
                                     split = c(1.0),
                                     verbose = verbose, gr = neg_gr,
                                     newrat_H0_inv = newrat_H0_inv,
                                     record_curvature = record_curvature),
    ## cmaes_newrat: CMA-ES global search to escape poor starts, then newrat
    ## quasi-Newton polish. Best of both worlds for near-unit-root models.
    "cmaes_newrat" = combined_optimize(neg_lp, par_init, lower, upper,
                                       max_iter = total_iter,
                                       stages = c("cmaes", "newrat"),
                                       split = c(0.7, 0.3),
                                       verbose = verbose, gr = neg_gr,
                                       newrat_H0_inv = newrat_H0_inv,
                                       record_curvature = record_curvature),
    ## Unknown method: fail loud rather than silently running a different
    ## optimizer (cmaes+nmkb) than the one requested (e.g. a typo'd method).
    stop(".run_mode_finding: unknown mode-finding method \"", method, "\". ",
         "Valid: newrat, cmaes_newrat, cmaes, nmkb, jade, nelder, combined, ",
         "cmaes_nmkb, cmaes_jade. (L-BFGS-B is not a standalone method -- it ",
         "runs as the polish stage inside \"combined\" and \"cmaes_jade\".)",
         call. = FALSE)
  )

  ## Unname: csminwel (newrat) propagates the parameter-vector name onto the
  ## objective value, so -res$value can arrive named. The mode log-posterior
  ## must be a clean unnamed scalar regardless of optimizer (downstream
  ## consumers such as the Laplace marginal likelihood compare it by value).
  mode_logpost <- unname(-res$value)

  if (!is.null(transform)) {
    theta_mode <- transform$to_constrained(res$par)
    names(theta_mode) <- par_names
  } else {
    theta_mode <- res$par
  }

  if (verbose) {
    cat(sprintf("\n  Mode-finding complete (%s)\n", method))
    cat(sprintf("    Log-posterior: %.4f -> %.4f (?? = %+.4f)\n",
                lp0, mode_logpost, mode_logpost - lp0))
    cat(sprintf("    Iterations: %d  |  %s\n", res$iterations, res$message))
  }

  out <- list(
    theta_mode  = theta_mode,
    logpost     = mode_logpost,
    convergence = res$convergence,
    iterations  = res$iterations,
    method      = method
  )

  ## Opt-in: propagate the newrat/csminwel final BFGS inverse-Hessian
  ## (res$H_bfgs, working space) back to THETA-space so callers get a
  ## covariance-like matrix directly comparable to the analytic/FD Hessian.
  ## H_bfgs (working space) approximates the inverse of the working-space
  ## neg-logpost Hessian, i.e. a covariance: Cov_theta = J %*% Cov_eta %*% J'
  ## with J = diag(dtheta/deta) (separable transform -> elementwise scaling,
  ## mirrors the H_theta -> H_eta sandwich used above for hessian_fn_working).
  if (isTRUE(record_curvature) && !is.null(res$H_bfgs)) {
    H_bfgs_theta <- res$H_bfgs
    if (!is.null(transform)) {
      J <- transform$dtheta_deta(res$par)   # diagonal dtheta/deta at the mode (eta-space)
      H_bfgs_theta <- J * H_bfgs_theta * rep(J, each = length(J))
    }
    dimnames(H_bfgs_theta) <- list(par_names, par_names)
    out$H_bfgs <- H_bfgs_theta
  }
  out
}


#' RWMH proposal covariance for the one-call estimation entry points
#'
#' \code{.run_mode_finding()} deliberately returns only the mode and its
#' convergence record -- it is also the optimiser core for
#' \code{method_of_moments()} and the HANK AR runner, neither of which wants
#' to pay for a Hessian. The samplers, however, need posterior curvature.
#'
#' HISTORY (fixed 2026-09-04). \code{run_full_estimation()} and
#' \code{estimate-runner.R} both wrote
#' \code{if (!is.null(mode_res$V_mode)) ... else diag(prior_spec$std^2)}.
#' Because \code{.run_mode_finding()} never sets \code{V_mode}, that condition
#' was NEVER true: the one-call API always proposed from PRIOR variances and
#' silently discarded the posterior curvature. It survived on models whose
#' priors happen to sit at the posterior scale (fs2000: 10\% acceptance) and
#' froze outright on a well-identified one (a 9-parameter NK fixture: 0\%
#' acceptance, ZERO posterior variance -- every draw equal to the mode), where
#' the inverse-Hessian proposal samples at 23.5\%.
#'
#' This helper computes the Hessian at the mode and runs it through exactly the
#' same regularisation \code{run_mode_finding()} uses
#' (\code{.proposal_cov_from_hessian()}: \code{.make_pd} for a non-finite or
#' ill-conditioned Hessian, then eigen-basis capping at the prior scale). If
#' the Hessian cannot be computed at all it falls back to the prior-variance
#' diagonal as before -- but LOUDLY, because a silent fallback here is
#' indistinguishable from a working sampler until someone checks the
#' acceptance rate.
#'
#' @param log_post_fn Log-posterior closure.
#' @param theta_mode Named mode vector.
#' @param prior_spec Prior specification.
#' @param verbose Print progress.
#' @return An n_par x n_par proposal covariance with dimnames.
#' @noRd
.sampler_proposal_cov <- function(log_post_fn, theta_mode, prior_spec,
                                  verbose = TRUE) {
  n_par     <- length(theta_mode)
  opt_scale <- 2.38^2 / n_par
  prior_fallback <- function() {
    S <- diag(prior_spec$std^2, nrow = n_par) * opt_scale
    dimnames(S) <- list(names(theta_mode), names(theta_mode))
    S
  }

  h    <- max(1e-4, 1e-4 * max(abs(theta_mode)))
  hess <- tryCatch(num_hessian(log_post_fn, theta_mode, h = h),
                   error = function(e) NULL)
  if (is.null(hess) || !any(is.finite(hess))) {
    warning("run_full_estimation: could not evaluate the posterior Hessian at ",
            "the mode, so the RWMH proposal falls back to the PRIOR variances. ",
            "That is only a sensible proposal when the prior and posterior are ",
            "on a similar scale; on a well-identified posterior it can freeze ",
            "the chain (acceptance ~0). Check `acceptance_rate`, and supply a ",
            "proposal covariance explicitly (mcmc(..., Sigma_prop = )) if it ",
            "is low.", call. = FALSE)
    return(prior_fallback())
  }

  S <- tryCatch(
    .proposal_cov_from_hessian(hess, prior_spec, theta_mode, verbose)$Sigma,
    error = function(e) NULL)
  if (is.null(S)) return(prior_fallback())
  dimnames(S) <- list(names(theta_mode), names(theta_mode))
  S
}
