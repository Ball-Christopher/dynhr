## R/passport-run.R
## --------------------------------------------------------------------------
## run_estimation_passport() -- the RUNNER half of the passport pair.
##
## This file is the complement to R/estimation-passport.R's
## estimation_passport(), which is a thin decision layer that SURFACES
## already-computed diagnostic results and never recomputes anything.
## run_estimation_passport() is the opposite: given only (model, data), it
## COMPUTES the paper's seven pre-flight checks directly, from scratch, using
## dynhr's existing solve/estimate machinery -- one call, no assembled
## diagnostic suite required first.
##
## The two are meant to be chained: run_estimation_passport() is the cheap
## "should I even start a long estimation?" gate; estimation_passport() is the
## thorough "can I trust the posterior I already have?" verdict once a full
## diagnostic suite exists. See @seealso on both functions.
##
## THE SEVEN CHECKS (pathological-dsge-paper, DYNHR_GAPS.md gap #6, section
## "The Estimation Passport" @sec-passport; thresholds cited per-check below)
## --------------------------------------------------------------------------
##   1. eigen        -- near-unit-root top eigenvalue of the solved TT at the
##                       mode (@sec-bk-nur; threshold 0.95, paper.qmd line 851
##                       "Above ~0.95, route the filter through lik_init=auto").
##   2. hessian       -- exact posterior-Hessian condition number + spectrum at
##                       the mode (@sec-aniso; paper.qmd lines 852-857: a
##                       correctly measured 1e5-1e8 is amber/expected and calls
##                       for a full inverse-Hessian proposal; ~1e11 (the FD
##                       artifact regime) is red).
##   3. bk_share      -- BK-feasible share of the prior with a Clopper-Pearson
##                       95% CI (replication/00_diagnose_nk_small.R part (1);
##                       14_mc_error.R's `cp()` helper for the exact CI).
##   4. collision     -- proposal/parameter collision rate at the mode (BK-
##                       infeasible tuned-proposal share; 00_diagnose_nk_small.R
##                       part (4); paper.qmd tbl-pathologies: nk_small ~11-16%
##                       is the "binding" regime, sw2007-class ~0.05% is
##                       "benign").
##   5. contraction   -- prior-to-posterior contraction table (@sec-weakid-
##                       contract; posterior-sd / prior-sd per parameter;
##                       paper.qmd line 641: ratios near 1 mean "the data have
##                       not moved the parameter").
##   6. multistart    -- seeded-multistart log-posterior spread (@sec-mode;
##                       reuses d7_mode_robustness(), R/diag-mcmc-d7-mode-
##                       robustness.R, tol=3 nats "same basin" / gap_warn=50).
##   7. sbc_lite      -- an SBC smoke check (NOT a certification): a handful of
##                       replications/draws via dynhr_sbc(), reporting its own
##                       uniformity verdict (@sec-sbc). Deliberately small by
##                       design (see sbc_full_enabled()/DYNHR_SBC_FULL
##                       convention in tests/testthat/helper-tolerances.R for
##                       the "keep it fast, it's a smoke check" precedent).
##
## Every check can error or be skipped (checks=) without aborting the whole
## call: a failed/absent check is reported "not_assessed", exactly the
## never-fabricate principle R/estimation-passport.R already documents.
## --------------------------------------------------------------------------


# ---- the green/amber/red threshold mapper ----------------------------------
#
# One small, testable, pure function per check family. Never called directly
# by users; exposed (unexported, dot-prefixed) so tests can drive it with
# synthetic values without re-running any solve.
#
# check: one of "eigenvalue", "condition_number". Both are ratio-like
# quantities that are monotonically "worse" as they increase, so the mapper is
# a simple two-cutoff ladder green -> amber -> red.
.passport_grade <- function(value, check = c("eigenvalue", "condition_number"),
                            amber_at = NULL, red_at = NULL) {
  check <- match.arg(check)
  if (is.null(amber_at) || is.null(red_at)) {
    defaults <- switch(check,
      eigenvalue       = list(amber = 0.95, red = 1.0),
      condition_number = list(amber = 1e8,  red = 1e10))
    amber_at <- amber_at %||% defaults$amber
    red_at   <- red_at   %||% defaults$red
  }
  if (is.null(value) || length(value) == 0L || !is.finite(value))
    return(list(status = "not_assessed", value = NA_real_,
                amber_at = amber_at, red_at = red_at))
  status <- if (value >= red_at) "red" else if (value >= amber_at) "amber" else "green"
  list(status = status, value = value, amber_at = amber_at, red_at = red_at)
}


# ---- small shared helpers ---------------------------------------------------

## Solve the model at a theta (named subset of estimated params), returning
## BK-satisfaction and the top modulus eigenvalue of the state block of TT.
## Mirrors replication/00_diagnose_nk_small.R's solve_at(), reusing dynhr's
## own solve_steady_state / solve_perturbation (never a bespoke re-derivation).
.passport_solve_at <- function(model, compiled, pv0, theta) {
  pv <- .apply_theta_to_params(model, theta, pv0)
  ss <- tryCatch(solve_steady_state(model, compiled, pv, verbose = FALSE),
                error = function(e) NULL)
  if (is.null(ss) || !isTRUE(ss$converged))
    return(list(bk = FALSE, top_eig = NA_real_))
  sol <- tryCatch(solve_perturbation(model, compiled, ss$ss, pv, verbose = FALSE),
                  error = function(e) NULL)
  if (is.null(sol)) return(list(bk = FALSE, top_eig = NA_real_))
  bk <- isTRUE(sol$bk_satisfied)
  te <- NA_real_
  if (bk && length(sol$state_idx)) {
    TT <- sol$ghx[sol$state_idx, , drop = FALSE]
    if (nrow(TT) > 0L)
      te <- max(Mod(eigen(TT, only.values = TRUE)$values))
  }
  list(bk = bk, top_eig = te)
}

## Timed wrapper: runs `expr_fn()`, catches any error, always returns a
## runtime. Used so every check reports its own wall time and a failure never
## propagates past its own entry.
.passport_timed <- function(expr_fn) {
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch(list(ok = TRUE, value = expr_fn(), error = NULL),
                  error = function(e) list(ok = FALSE, value = NULL,
                                           error = conditionMessage(e)))
  out$runtime <- proc.time()[["elapsed"]] - t0
  out
}

## Build a single not_assessed check entry (error or opted-out).
.passport_na_entry <- function(reason, runtime = NA_real_) {
  list(status = "not_assessed", value = NULL, threshold = NULL,
       reason = reason, runtime = runtime)
}


# ---------------------------------------------------------------------------
#' Run the estimation passport's seven pre-flight diagnostics
#'
#' The paper's "estimation passport" (pathological-dsge-paper, @sec-passport)
#' bundles seven cheap-but-decisive checks that should be run BEFORE
#' committing to a long estimation. This function computes all seven directly
#' from \code{(model, data)} -- no pre-existing diagnostic suite or mode-find
#' result required -- each with a green/amber/red threshold taken from the
#' paper, and returns a single structured report.
#'
#' This is the RUNNER half of the passport pair. Its sibling,
#' \code{\link{estimation_passport}}, is a decision layer that SURFACES
#' already-computed \code{deep_parameter_passport()} / diagnostic-suite
#' results and never recomputes; this function does the opposite -- it always
#' computes from scratch. The two compose naturally: run this function first
#' as a pre-flight gate, and once a full diagnostic suite exists, use
#' \code{estimation_passport()} for the thorough per-parameter verdict.
#'
#' Every check is independently guarded: an error inside one check is caught
#' and that check is reported \code{"not_assessed"} (never a fabricated
#' grade), and the remaining checks still run. Use \code{checks=} to skip
#' checks explicitly (e.g. to omit the slow SBC-lite check during
#' development).
#'
#' @param model      Parsed model (output of \code{\link{parse_mod}}).
#' @param data       Observation matrix/data.frame (T x n_obs, or n_obs x T --
#'   passed straight through to \code{\link{run_mode_finding}}), or a path
#'   readable the same way.
#' @param prior_spec Prior-spec data.frame (output of
#'   \code{extract_prior_spec} / \code{\link{prior_spec}}). If
#'   \code{NULL}, extracted from \code{model}.
#' @param obs_vars   Character vector of observed variable names. If
#'   \code{NULL}, taken from \code{model$obs_vars} / \code{model$varobs_names}.
#' @param checks     Character vector selecting which of the seven checks to
#'   run: any of \code{"eigen"}, \code{"hessian"}, \code{"bk_share"},
#'   \code{"collision"}, \code{"contraction"}, \code{"multistart"},
#'   \code{"sbc_lite"}. Default = all seven. Checks not listed are reported
#'   \code{"not_assessed"} with reason \code{"skipped (checks=)"}.
#' @param n_prior_draws Compute budget for the BK-feasible-share check
#'   (default 500L; the paper uses 3000-10000, dial down for speed).
#' @param n_collision   Compute budget for the mode-collision check (default
#'   500L; the paper uses 4000-10000).
#' @param n_multistart  Number of dispersed starts for the multistart check
#'   (default 4L; the paper's D7-style gate needs >= 2).
#' @param sbc_n         List with \code{n_replications} and \code{n_draws} for
#'   the SBC-lite check (default \code{list(n_replications = 20L, n_draws =
#'   400L)} -- deliberately tiny; this is a smoke check, not a certification;
#'   see \code{sbc_full_enabled()} / \code{DYNHR_SBC_FULL} for the package's
#'   general "lighten SBC in tests" convention).
#' @param mode_result Optional pre-computed \code{\link{run_mode_finding}}
#'   result (class \code{dynhr_mode_result}). When supplied, checks 1/2/3/4/5
#'   reuse its mode/Hessian/Sigma_prop instead of re-running mode-finding
#'   (saves the single most expensive step when the caller already has a
#'   mode). Checks 6 (multistart) and 7 (SBC-lite) always run their own solves
#'   regardless (they are defined by re-running from dispersed starts / fresh
#'   simulated data).
#' @param n_iter      Mode-finding iteration budget passed to
#'   \code{run_mode_finding} (default 200L; dial down further for speed).
#' @param n_cores     Cores for the mode-finding / SBC steps that support
#'   parallelism (default \code{NULL} = serial).
#' @param seed        Integer seed for every stochastic sub-check (prior
#'   draws, collision draws, multistart perturbations, SBC replications).
#'   Reused per-check with a small per-check offset so checks are mutually
#'   reproducible but not identical draws.
#' @param verbose     Print progress. Default \code{FALSE}.
#' @param ... Passed through to the internal \code{run_mode_finding()} call
#'   (e.g. \code{method}, \code{me_variance}, \code{likelihood}).
#' @return An object of class \code{"dynhr_passport_run"}: a list with one
#'   named entry per check (\code{eigen}, \code{hessian}, \code{bk_share},
#'   \code{collision}, \code{contraction}, \code{multistart}, \code{sbc_lite}),
#'   each a list with \code{status} (one of \code{"green"}, \code{"amber"},
#'   \code{"red"}, \code{"not_assessed"}), \code{value}(s), \code{threshold}
#'   used, and \code{runtime} (seconds); plus \code{$meta} (model/data
#'   summary, checks requested, seed) and \code{$mode_result} (the underlying
#'   \code{run_mode_finding()} result, for reuse/inspection).
#' @seealso \code{\link{estimation_passport}} for the decision layer that
#'   surfaces an existing diagnostic suite / deep parameter passport instead
#'   of recomputing; \code{\link{run_mode_finding}}, \code{\link{dynhr_sbc}},
#'   \code{extract_prior_spec}.
#' @export
run_estimation_passport <- function(model, data,
                                    prior_spec = NULL,
                                    obs_vars   = NULL,
                                    checks = c("eigen", "hessian", "bk_share",
                                              "collision", "contraction",
                                              "multistart", "sbc_lite"),
                                    n_prior_draws = 500L,
                                    n_collision   = 500L,
                                    n_multistart  = 4L,
                                    sbc_n = list(n_replications = 20L, n_draws = 400L),
                                    mode_result = NULL,
                                    n_iter  = 200L,
                                    n_cores = NULL,
                                    seed    = 1L,
                                    verbose = FALSE,
                                    ...) {

  all_checks <- c("eigen", "hessian", "bk_share", "collision", "contraction",
                  "multistart", "sbc_lite")
  checks <- match.arg(checks, all_checks, several.ok = TRUE)
  run_check <- function(nm) nm %in% checks

  compiled <- compile_model(model, verbose = FALSE)
  if (is.null(prior_spec)) prior_spec <- extract_prior_spec(model, verbose = FALSE)
  if (is.null(obs_vars) || length(obs_vars) == 0L) {
    obs_vars <- model$obs_vars %||% model$varobs_names
    if (is.null(obs_vars) || length(obs_vars) == 0L)
      stop("run_estimation_passport: `obs_vars` not supplied and the model ",
           "declares no `varobs`.", call. = FALSE)
  }
  pv0 <- model$param_values
  pnames <- prior_spec$name

  out <- vector("list", length(all_checks))
  names(out) <- all_checks
  for (nm in all_checks) out[[nm]] <- .passport_na_entry("not run yet")

  ## ---- shared mode-find (checks 1/2/4/5 all key off it) -------------------
  need_mode <- any(run_check(c("eigen", "hessian", "collision", "contraction")))
  mf <- mode_result
  mf_timing <- NA_real_
  if (need_mode && is.null(mf)) {
    mf_run <- .passport_timed(function() {
      set.seed(seed)
      run_mode_finding(list(model = model, compiled = compiled), data, obs_vars,
                       n_iter = n_iter, verbose = verbose, n_cores = n_cores,
                       use_exact_hessian = run_check("hessian"), ...)
    })
    mf_timing <- mf_run$runtime
    if (isTRUE(mf_run$ok)) {
      mf <- mf_run$value
    } else {
      na_reason <- sprintf("mode-finding failed: %s", mf_run$error)
      for (nm in c("eigen", "hessian", "collision", "contraction"))
        if (run_check(nm)) out[[nm]] <- .passport_na_entry(na_reason, mf_timing)
    }
  }

  ## ---- check 1: near-unit-root top eigenvalue at the mode ------------------
  if (run_check("eigen")) {
    if (!is.null(mf)) {
      r <- .passport_timed(function() {
        theta_mode <- mf$theta_mode
        .passport_solve_at(model, compiled, pv0, theta_mode)
      })
      if (isTRUE(r$ok) && is.finite(r$value$top_eig)) {
        g <- .passport_grade(r$value$top_eig, "eigenvalue")
        out$eigen <- list(status = g$status, value = r$value$top_eig,
                          bk_satisfied = r$value$bk,
                          threshold = list(amber = g$amber_at, red = g$red_at),
                          runtime = mf_timing + r$runtime)
      } else if (isTRUE(r$ok)) {
        out$eigen <- .passport_na_entry(
          "mode is BK-infeasible or eigenvalue undefined", mf_timing + r$runtime)
      } else {
        out$eigen <- .passport_na_entry(r$error, mf_timing + r$runtime)
      }
    }
  }

  ## ---- check 2: exact posterior-Hessian condition number + spectrum -------
  if (run_check("hessian")) {
    if (!is.null(mf)) {
      H <- mf$hessian_exact
      if (!is.null(H) && all(is.finite(H))) {
        Hs <- (H + t(H)) / 2
        ev <- eigen(Hs, symmetric = TRUE, only.values = TRUE)$values
        cond <- max(abs(ev)) / min(abs(ev))
        g <- .passport_grade(cond, "condition_number")
        ## posterior_hessian() returns the Hessian of the LOG-POSTERIOR
        ## (unnegated); at a converged maximum it should be negative
        ## semi-definite, so a "wrong-sign" eigenvalue (paper.qmd
        ## "Limitations", 1 of 9 / 2 of 29 / 2 of 34 across the three models)
        ## is a POSITIVE eigenvalue here. Condition number is a ratio of
        ## |eigenvalues| either way, so it is unaffected by this convention.
        out$hessian <- list(status = g$status, value = cond,
                            eigenvalues = sort(abs(ev)),
                            n_wrong_sign = sum(ev > 0),
                            threshold = list(amber = g$amber_at, red = g$red_at),
                            runtime = mf_timing)
      } else {
        out$hessian <- .passport_na_entry(
          "exact Hessian unavailable at the mode (non-Gaussian/OBC likelihood, ME-extra, or singular solve)",
          mf_timing)
      }
    }
  }

  ## ---- check 3: BK-feasible share of the prior (Clopper-Pearson) ----------
  if (run_check("bk_share")) {
    r <- .passport_timed(function() {
      set.seed(seed + 1000L)
      sampler <- .smc_make_prior_sampler(prior_spec)
      n <- as.integer(n_prior_draws)
      bk <- logical(n); te <- rep(NA_real_, n)
      for (i in seq_len(n)) {
        theta <- sampler()
        names(theta) <- pnames
        s <- .passport_solve_at(model, compiled, pv0, theta)
        bk[i] <- isTRUE(s$bk); te[i] <- s$top_eig
      }
      list(bk = bk, te = te, n = n)
    })
    if (isTRUE(r$ok)) {
      v <- r$value
      x <- sum(v$bk)
      bt <- stats::binom.test(x, v$n)
      nur_of_feasible <- if (x > 0) mean(v$te[v$bk] > 0.95, na.rm = TRUE) else NA_real_
      ## bk_share itself has no red/green ladder in the paper (a large
      ## infeasible share is "benign" per tbl-pathologies as long as
      ## collisions are near zero) -- it is reported as an informational
      ## proportion + CI, status derived only from whether it was measured.
      out$bk_share <- list(
        status = "green",
        value = x / v$n,
        conf_int = as.numeric(bt$conf.int),
        n = v$n, n_feasible = x,
        near_unit_root_share_of_feasible = nur_of_feasible,
        threshold = NULL,
        runtime = r$runtime)
    } else {
      out$bk_share <- .passport_na_entry(r$error, r$runtime)
    }
  }

  ## ---- check 4: proposal/parameter collision rate at the mode -------------
  if (run_check("collision")) {
    if (!is.null(mf) && !is.null(mf$Sigma_prop)) {
      r <- .passport_timed(function() {
        set.seed(seed + 2000L)
        theta_mode <- mf$theta_mode
        Sig <- mf$Sigma_prop
        Lch <- tryCatch(chol(Sig),
                        error = function(e)
                          chol(as.matrix(Matrix::nearPD(Sig)$mat)))
        n <- as.integer(n_collision)
        coll <- logical(n)
        for (i in seq_len(n)) {
          prop <- theta_mode + as.numeric(crossprod(Lch, stats::rnorm(length(theta_mode))))
          names(prop) <- names(theta_mode)
          coll[i] <- !isTRUE(.passport_solve_at(model, compiled, pv0, prop)$bk)
        }
        list(coll = coll, n = n)
      })
      if (isTRUE(r$ok)) {
        v <- r$value
        x <- sum(v$coll)
        bt <- stats::binom.test(x, v$n)
        rate <- x / v$n
        ## Thresholds from tbl-pathologies / @sec-bk: collisions near 0
        ## ( <1%, sw2007/micro-class) are benign; the nk_small-class binding
        ## regime is 11-16%. Amber at >5%, red at >10% (inside the paper's
        ## observed binding band) -- collisions signal weak identification,
        ## not sampler brokenness, but are worth a practitioner's pause.
        status <- if (rate >= 0.10) "red" else if (rate >= 0.05) "amber" else "green"
        out$collision <- list(status = status, value = rate,
                              conf_int = as.numeric(bt$conf.int),
                              n = v$n, n_collisions = x,
                              threshold = list(amber = 0.05, red = 0.10),
                              runtime = mf_timing + r$runtime)
      } else {
        out$collision <- .passport_na_entry(r$error, mf_timing + r$runtime)
      }
    } else if (!is.null(mf)) {
      ## mf exists but carries no proposal covariance. (When mf is NULL the
      ## mode-finding-failure reason was already recorded above -- keep it.)
      out$collision <- .passport_na_entry("mode-finding did not return a Sigma_prop", mf_timing)
    }
  }

  ## ---- check 5: prior-to-posterior contraction table ----------------------
  if (run_check("contraction")) {
    if (!is.null(mf) && !is.null(mf$Sigma_prop)) {
      r <- .passport_timed(function() {
        Sig <- mf$Sigma_prop
        ## Sigma_prop is the tuned RWMH proposal covariance == (2.38^2/np) *
        ## posterior covariance (standard RWMH scaling; matches
        ## 00_diagnose_nk_small.R's use of the same object for sd_post).
        np <- length(mf$theta_mode)
        sd_post <- sqrt(diag(Sig)) * sqrt(np / 2.38^2)
        names(sd_post) <- names(mf$theta_mode)
        prior_sd <- setNames(prior_spec$std, prior_spec$name)
        common <- intersect(names(sd_post), names(prior_sd))
        ratio <- sd_post[common] / prior_sd[common]
        data.frame(param = common, prior_sd = prior_sd[common],
                  posterior_sd = sd_post[common], ratio = ratio,
                  row.names = NULL)
      })
      if (isTRUE(r$ok)) {
        tbl <- r$value
        ## Ratios near 1 flag uninformed parameters (paper.qmd line 641).
        ## amber at >= 0.7 (paper's kappa example, "barely moves"), red at
        ## >= 0.95 (essentially unmoved, e.g. gammaQ at 1.00).
        worst <- if (nrow(tbl)) tbl$param[which.max(tbl$ratio)] else NA_character_
        worst_ratio <- if (nrow(tbl)) max(tbl$ratio) else NA_real_
        status <- if (!is.finite(worst_ratio)) "not_assessed"
                  else if (worst_ratio >= 0.95) "red"
                  else if (worst_ratio >= 0.70) "amber"
                  else "green"
        out$contraction <- list(status = status, value = worst_ratio,
                                table = tbl, worst_param = worst,
                                threshold = list(amber = 0.70, red = 0.95),
                                runtime = mf_timing + r$runtime)
      } else {
        out$contraction <- .passport_na_entry(r$error, mf_timing + r$runtime)
      }
    } else if (!is.null(mf)) {
      ## mf exists but carries no proposal covariance. (When mf is NULL the
      ## mode-finding-failure reason was already recorded above -- keep it.)
      out$contraction <- .passport_na_entry("mode-finding did not return a Sigma_prop", mf_timing)
    }
  }

  ## ---- check 6: seeded-multistart log-posterior spread --------------------
  ## run_mode_finding() (the public entry point) always starts from the prior
  ## mean (extract_prior_spec()$mean, computed internally -- it takes no
  ## theta_init override), so a dispersed-start multistart cannot be built by
  ## calling it repeatedly. Instead this reuses the same internal machinery
  ## run_mode_finding() itself calls for its serial path -- make_log_posterior()
  ## + .run_mode_finding() (R/mode-orchestrate.R) -- with an explicit,
  ## per-chain dispersed theta_init. Same optimizer, same objective, just a
  ## different starting point per chain.
  if (run_check("multistart")) {
    r <- .passport_timed(function() {
      n_ms <- max(2L, as.integer(n_multistart))
      chains <- vector("list", n_ms)
      base_theta <- setNames(prior_spec$mean, pnames)
      lp_fn <- make_log_posterior(model, data, prior_spec, obs_vars, compiled)
      for (k in seq_len(n_ms)) {
        set.seed(seed + 3000L + k)
        theta_init <- base_theta
        if (k > 1L) {
          ## Disperse starts 2..n around the prior mean, scaled by prior sd
          ## (mirrors run_mode_finding's own dispersed-start convention).
          jitter <- stats::rnorm(length(theta_init), sd = 0.5 * prior_spec$std)
          theta_init <- theta_init + jitter
        }
        names(theta_init) <- pnames
        mode_k <- tryCatch(
          .run_mode_finding(lp_fn, theta_init, prior_spec,
                           nm_maxit = n_iter, method = "newrat", verbose = FALSE),
          error = function(e) NULL)
        ## A chain that never left the infeasible floor (~-1e20, the
        ## log-posterior's -Inf substitute) is a FAILED start, not a worse
        ## mode -- keeping it would inflate the spread by ~1e20 nats and
        ## drown the signal. Count it as non-converged instead.
        chains[[k]] <- if (!is.null(mode_k) && is.finite(mode_k$logpost) &&
                           mode_k$logpost > -1e19)
          list(logpost = mode_k$logpost, theta_mode = mode_k$theta_mode)
        else NULL
      }
      chains <- Filter(Negate(is.null), chains)
      chains
    })
    if (isTRUE(r$ok) && length(r$value) >= 2L) {
      d7 <- tryCatch(d7_mode_robustness(list(results = r$value)),
                    error = function(e) NULL)
      if (!is.null(d7)) {
        status <- if (isTRUE(d7$pass)) "green"
                  else if (isFALSE(d7$pass)) "amber" else "not_assessed"
        spread <- diff(range(d7$result$logposts))
        ## Red when the spread crosses the paper's own gap_warn convention
        ## (50 nats, R/diag-mcmc-d7-mode-robustness.R): "a large spread means
        ## the mode is a ridge" (paper.qmd line 865).
        if (identical(status, "amber") && is.finite(spread) && spread > 50)
          status <- "red"
        out$multistart <- list(status = status, value = spread,
                               n_chains = d7$result$n_chains,
                               n_at_best = d7$result$n_at_best,
                               n_basins = d7$result$n_basins,
                               logposts = d7$result$logposts,
                               threshold = list(tol = 3.0, gap_warn = 50),
                               runtime = r$runtime)
      } else {
        out$multistart <- .passport_na_entry("d7_mode_robustness failed on the chain set", r$runtime)
      }
    } else if (isTRUE(r$ok)) {
      out$multistart <- .passport_na_entry(
        sprintf("only %d/%d seeded start(s) converged", length(r$value), n_multistart),
        r$runtime)
    } else {
      out$multistart <- .passport_na_entry(r$error, r$runtime)
    }
  }

  ## ---- check 7: SBC-lite (smoke check, NOT a certification) ---------------
  if (run_check("sbc_lite")) {
    r <- .passport_timed(function() {
      n_rep   <- as.integer(sbc_n$n_replications %||% 20L)
      n_draws <- as.integer(sbc_n$n_draws %||% 400L)
      dynhr_sbc(model, obs_vars = obs_vars,
               T_obs = sbc_n$T_obs %||% 60L,
               n_replications = n_rep, n_draws = n_draws,
               n_burn = as.integer(n_draws / 2), thin = sbc_n$thin %||% 2L,
               sampler = sbc_n$sampler %||% "rwmh",
               lik_init = sbc_n$lik_init %||% "auto",
               seed = seed, verbose = FALSE, n_cores = n_cores)
    })
    if (isTRUE(r$ok)) {
      s <- r$value
      verdict <- s$uniformity$verdict %||% NA_character_
      succ <- s$n_replications - s$n_failed
      ## dynhr_sbc()'s own uniformity verdict (sbc_uniformity_test(), R/
      ## validate-sbc.R) is one of "calibrated" / "suspect" / "miscalibrated"
      ## / "insufficient" -- read verbatim (never re-derived), per the
      ## paper's @sec-sbc "report the package's own uniformity verdict".
      ## "insufficient" (no usable rank spread/GOF signal, e.g. every
      ## replication landed on the same rank) is a degraded-input report,
      ## not a miscalibration finding -- route it to not_assessed like the
      ## too-few-successes case rather than "red".
      status <- if (is.na(verdict)) "not_assessed"
                else if (succ < 5L) "not_assessed"  # too few reps to read anything
                else if (identical(verdict, "insufficient")) "not_assessed"
                else if (identical(verdict, "calibrated")) "green"
                else if (identical(verdict, "suspect")) "amber"
                else "red"
      out$sbc_lite <- list(status = status, value = verdict,
                           n_succeeded = succ, n_attempted = s$n_replications,
                           n_draws = s$settings$n_draws,
                           note = "SMOKE CHECK ONLY -- a handful of ranks, not a certification; rerun with a larger sbc_n for publication-grade SBC",
                           threshold = NULL,
                           runtime = r$runtime)
    } else {
      out$sbc_lite <- .passport_na_entry(r$error, r$runtime)
    }
  }

  ## ---- checks explicitly opted out ------------------------------------------
  for (nm in setdiff(all_checks, checks))
    out[[nm]] <- .passport_na_entry("skipped (checks=)")

  result <- c(out, list(
    meta = list(
      checks_requested = checks,
      seed = seed,
      n_prior_draws = n_prior_draws,
      n_collision = n_collision,
      n_multistart = n_multistart,
      sbc_n = sbc_n
    ),
    mode_result = mf
  ))
  class(result) <- "dynhr_passport_run"
  result
}


# ---------------------------------------------------------------------------
#' Print an estimation-passport run report card
#'
#' @param x A \code{dynhr_passport_run}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @method print dynhr_passport_run
#' @export
print.dynhr_passport_run <- function(x, ...) {
  bar <- strrep("=", 76L)
  badge <- function(s) switch(s,
    green = "[GREEN]", amber = "[AMBER]", red = "[RED]  ",
    not_assessed = "[N/A]  ", "[?]    ")

  labels <- c(
    eigen       = "1. Near-unit-root top eigenvalue",
    hessian     = "2. Posterior Hessian condition number",
    bk_share    = "3. BK-feasible prior share",
    collision   = "4. Proposal collision rate at mode",
    contraction = "5. Prior-to-posterior contraction (worst param)",
    multistart  = "6. Seeded-multistart logpost spread",
    sbc_lite    = "7. SBC-lite (smoke check)")

  cat(bar, "\n")
  cat("  ESTIMATION PASSPORT -- RUN REPORT\n")
  cat(bar, "\n")

  for (nm in names(labels)) {
    ck <- x[[nm]]
    st <- ck$status %||% "not_assessed"
    val_str <- if (identical(st, "not_assessed")) {
      sprintf("(not assessed: %s)", ck$reason %||% "unknown")
    } else {
      v <- ck$value
      if (is.numeric(v) && length(v) == 1L) sprintf("%.4g", v)
      else if (is.character(v)) v
      else "(see $result for detail)"
    }
    rt <- ck$runtime
    rt_str <- if (is.numeric(rt) && length(rt) == 1L && is.finite(rt))
      sprintf("  [%.1fs]", rt) else ""
    cat(sprintf("  %s %-48s %s%s\n", badge(st), labels[[nm]], val_str, rt_str))
  }
  cat(bar, "\n")
  invisible(x)
}
