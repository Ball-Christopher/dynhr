## R/posterior-closure.R
## --------------------------------------------------------------------------
## ONE closure builder for every log-posterior factory.
##
## Thirteen factories -- make_log_posterior()'s gaussian and student_t
## branches plus make_log_posterior_{pruned,pruned3,pskf,pskf_order2,tpf,
## whittle,cumulant,sv_rbpf,obc,obc_pkf,obc_ppf}() -- used to re-implement the
## SAME per-draw contract:
##
##   log_prior -> support check -> .apply_theta_to_params -> steady state
##   (warm-started, reset/persisted) -> solve -> loglik -> me-floor latch
##   -> system prior -> power tempering -> list(logpost, loglik, logprior, ...)
##
## Copy-paste is not merely repetitive here, it LOSES BUG FIXES: the `power`
## plumbing, the `ss_warm <<-` warm-start idiom and the `.me_floor_checked`
## latch each had to be re-applied to between four and seven files, and each
## round missed some. `.make_posterior_closure()` carries that skeleton ONCE;
## each factory is now a short adapter that validates its own arguments and
## supplies a `loglik_fn` (plus a `solve_fn` when its branch needs an
## order-2/3 or OBC solve instead of the shared first-order path).
##
## WHERE THE FACTORIES GENUINELY DIFFER, that difference is a HOOK ARGUMENT,
## never a second copy of the skeleton:
##   * `warm_start` / `warm_retry` -- cumulant and the three OBC factories
##     solve the steady state cold on every draw; pruned/pruned3 warm-start
##     but do NOT retry cold on a failure; the rest warm-start with a retry.
##   * `stationarity` -- "spectral" reuses the QZ eigenvalues, "eigen" always
##     re-runs eigen() on the state block (cumulant + the OBC trio did), and
##     "none" skips the guard entirely (pruned/pruned3/student_t had none).
##   * `system_prior_mode` -- "lp" folds the system-prior density INTO
##     $logprior (gaussian, tpf, whittle, cumulant); "extra" leaves $logprior
##     as the parameter prior alone and only adds the density to $logpost
##     (pruned, pruned3, pskf, pskf_order2, student_t). These are DIFFERENT
##     returned contracts, pinned by test-posterior-closure-parity.R, so they
##     are a mode flag rather than a unification.
##   * `reject_fields` -- the failure shape. Only make_log_posterior_obc_pkf()
##     carries an extra `regime_path = NULL` on its -Inf returns.
##   * `pass_dots` -- make_log_posterior_tpf()'s closure is
##     `function(theta, U_list = NULL)`; every other closure is
##     `function(theta)` and must STAY that way (the parity test pins
##     `names(formals())`).
##
## The builder is deliberately NOT exported and NOT S3: it is the internal
## spine of the factories, and the factories are the public surface.
## --------------------------------------------------------------------------


## Sentinel returned by a `solve_fn` / `loglik_fn` that wants the draw
## rejected. `NULL` is the common case (plain -Inf); the object form carries a
## non-default `logpost`/`loglik` (the infeasible_penalty branch) and any
## extra fields that reject must report ($infeasible, $violation).
#' @noRd
.posterior_reject <- function(logpost = -Inf, loglik = -Inf, extra = NULL) {
  structure(list(logpost = logpost, loglik = loglik, extra = extra),
            class = "dynhr_posterior_reject")
}

#' @noRd
.is_posterior_reject <- function(x)
  is.null(x) || inherits(x, "dynhr_posterior_reject")


## Spectral radius of the state-transition block.
##
## `reuse = TRUE`: take the stable generalized eigenvalues the QZ solve
## already produced (verified bit-equal to eigen(ghx_state) across the model
## corpus to ~1e-13), saving ~100 us per feasible draw. `reuse = FALSE`
## always re-runs eigen() -- what the cumulant and OBC factories did, kept
## because their fallback is the only path they ever took and a switch would
## move their goldens in the last ulp.
#' @noRd
.posterior_spectral_radius <- function(dr, reuse = TRUE) {
  ns <- length(dr$state_idx)
  ev <- dr$eigenvalues
  if (reuse && !is.null(ev) && length(ev) >= ns)
    max(Mod(ev[seq_len(ns)]))
  else
    max(Mod(eigen(dr$ghx[dr$state_idx, , drop = FALSE],
                  only.values = TRUE)$values))
}


## The shared first-order solve: system matrices -> decision rule -> BK check
## -> optional stationarity guard. Returns list(sys, dr) or NULL (reject).
##
## `stationarity`:
##   "none"     no guard (pruned / pruned3 / student_t)
##   "spectral" reuse the QZ eigenvalues (gaussian, pskf, tpf, whittle, sv)
##   "eigen"    always eigen() the state block (cumulant, obc, obc_pkf, ppf)
#' @noRd
.posterior_solve1 <- function(model, compiled, sys_cache, ss, params,
                              stationarity = c("spectral", "eigen", "none")) {
  stationarity <- match.arg(stationarity)
  sys <- extract_system_matrices_fast(sys_cache, ss, params)
  dr  <- .solve_from_system(sys, model, compiled, ss, params, FALSE)
  if (is.null(dr) || !isTRUE(dr$bk_satisfied)) return(NULL)
  if (stationarity != "none") {
    sr <- .posterior_spectral_radius(dr, reuse = identical(stationarity,
                                                           "spectral"))
    if (sr >= 1) return(NULL)
  }
  list(sys = sys, dr = dr)
}


## ---------------------------------------------------------------------------
## One-entry structural-parameter cache
##
## A draw that moves only covariance-side parameters (an estimated shock
## standard deviation or correlation, a measurement-error variance) leaves the
## steady state, the system matrices and the decision rule exactly as they
## were; re-solving them is most of the cost of the evaluation. The cache holds
## the last feasible draw's solve and returns it when the parameters that CAN
## reach the solve are bit-for-bit those of that draw.
##
## Which parameters those are is decided from the model, never from names:
##   * the declared parameters the dynamic equations, the planner objective
##     and the model-local definitions reference (model$equation_param_names),
##     which are the only parameters the residual and Jacobian functions read;
##   * every symbol the steady_state_model statements mention, and every name
##     they assign (the statements evaluate in an environment holding ALL
##     parameters, so they may read any name they spell out).
## Everything else -- a shock standard deviation or correlation injected under
## a shock / "corr a,b" key, a declared parameter used only by the shocks block
## -- cannot reach the steady state or the decision rule. A model whose
## structural set cannot be established (no equation_param_names, an
## unparseable steady-state statement) gets no cache.
## ---------------------------------------------------------------------------

## Names whose values the steady state / decision rule can depend on, or NULL
## when that cannot be established (caller then does not cache).
#' @noRd
.structural_param_names <- function(model) {
  eqp <- model$equation_param_names
  if (is.null(eqp)) return(NULL)
  ssm_names <- character(0)
  for (a in model$steady_state_model) {
    if (!is.list(a) || is.null(a$text) || is.null(a$name)) return(NULL)
    ex <- parse(text = a$text, keep.source = FALSE)
    ssm_names <- c(ssm_names, all.names(ex), a$name)
  }
  unique(c(eqp, ssm_names))
}

## A per-closure cache environment (NULL when no cache can be built). The
## environment is created per call, so closures never share one, and each
## mirai daemon deserialises its own copy.
#' @noRd
.structural_cache_new <- function(model) {
  nms <- .structural_param_names(model)
  if (is.null(nms)) return(NULL)
  cache <- new.env(parent = emptyenv())
  cache$names <- nms
  cache$entry <- NULL
  cache
}

## The structural slice of a parameter vector (names absent from it give NA,
## the same way on every call).
#' @noRd
.structural_cache_key <- function(cache, params) params[cache$names]

## Entry for the cached solve at `key`, or NULL. A lookup that does not hit
## DROPS the entry: it is then only ever the immediately preceding feasible
## draw, which is what makes a hit reproduce exactly what the warm-started
## re-solve would have returned.
#' @noRd
.structural_cache_lookup <- function(cache, key) {
  e <- cache$entry
  if (!is.null(e) && identical(key, e$key, num.eq = FALSE)) return(e)
  cache$entry <- NULL
  NULL
}

## Remember a feasible draw. `params_in` is the parameter vector the steady
## state was solved at; `params_out` the one after any steady_state_model
## re-derivation (equal to params_in when the model has none). Only the names
## the re-derivation changed or added are kept: they are functions of the
## structural parameters alone, so they can be written onto a later draw's own
## parameter vector.
#' @noRd
.structural_cache_store <- function(cache, key, params_in, params_out,
                                    ss_result, sol) {
  chg <- character(0)
  if (!identical(params_in, params_out)) {
    common <- intersect(names(params_out), names(params_in))
    a <- params_out[common]; b <- params_in[common]
    same <- a == b
    na_ab <- is.na(same)
    same[na_ab] <- is.na(a)[na_ab] & is.na(b)[na_ab]
    chg <- c(common[!same], setdiff(names(params_out), names(params_in)))
    ## keep params_out's own order, so appended names land where they did
    chg <- names(params_out)[names(params_out) %in% chg]
  }
  cache$entry <- list(key = key, ss_result = ss_result, sol = sol,
                      delta = chg, delta_values = params_out[chg])
  invisible(NULL)
}

## The draw's own parameter vector with the cached steady_state_model
## re-derivation applied.
#' @noRd
.structural_cache_params <- function(entry, params) {
  if (length(entry$delta)) params[entry$delta] <- entry$delta_values
  params
}

#' Build a log-posterior closure from a likelihood hook
#'
#' The shared per-draw skeleton behind every \code{make_log_posterior_*}
#' factory. See the file header for the design and for what each hook exists
#' to express.
#'
#' @param model,data,prior_spec,obs_vars,compiled As in
#'   \code{make_log_posterior}. \code{data} and \code{obs_vars} are captured
#'   for the hooks' convenience only -- the builder itself never touches the
#'   observation matrix, because each likelihood wants a different orientation
#'   and the adapter has already fixed that at factory time.
#' @param loglik_fn Hook evaluating the likelihood for one draw. Called as
#'   \code{loglik_fn(sol = , params = , ss = , theta = , me_floor_check = ,
#'   ...)} where \code{sol} is whatever \code{solve_fn} returned, \code{ss} is
#'   the converged steady-state vector, \code{params} has already absorbed any
#'   \code{steady_state_model}-computed parameter, and \code{me_floor_check}
#'   is \code{TRUE} only on the first draw that reaches the likelihood (the
#'   once-per-closure me-floor hazard guard). It returns either
#'   \code{list(loglik = , extra = list(), dr = , Sigma_e = )} or a
#'   \code{.posterior_reject()} / \code{NULL} to reject the draw.
#'   \code{$extra} becomes additional fields on the returned list, appended
#'   after \code{logprior} in order; \code{$dr} and \code{$Sigma_e} override
#'   what the system-prior hook is handed.
#' @param solve_fn Hook mapping \code{(model, compiled, sys_cache, ss,
#'   params, theta)} to whatever \code{loglik_fn} needs, or a reject.
#'   \code{NULL} (default) uses the shared first-order path
#'   \code{.posterior_solve1()} with \code{stationarity}.
#' @param sys_cache Pre-built \code{cache_system_structure(compiled)}; built
#'   here when \code{NULL}.
#' @param power Tempering exponent zeta, ALREADY resolved by the adapter (via
#'   \code{.resolve_power_posterior()} or \code{.dynhr_opt()}) so it is a
#'   fixed property of the closure.
#' @param warm_start Keep and reuse the converged steady state across draws.
#' @param warm_retry On a failed warm-started solve, retry once from the cold
#'   initval guess before declaring the draw infeasible.
#' @param stationarity Guard mode for the default \code{solve_fn}; see
#'   \code{.posterior_solve1()}.
#' @param needs_me_floor Latch the once-per-closure me-floor guard after the
#'   first likelihood evaluation.
#' @param system_prior The factory's \code{system_priors} argument (\code{NULL}
#'   = no system prior).
#' @param system_prior_fn Hook
#'   \code{function(spec, theta, sol, params, res)} returning the system-prior
#'   log-density, where \code{spec} is \code{system_prior}. \code{NULL} uses
#'   \code{.eval_system_priors()} with the standard state list; the PSKF pair
#'   supplies \code{.pskf_system_prior_sum()} instead, because its public
#'   contract is a bare list of \code{function(dr, params)} closures.
#' @param system_prior_mode \code{"lp"} folds the density into
#'   \code{$logprior} (and rejects a non-finite density early);
#'   \code{"extra"} leaves \code{$logprior} alone and adds the density to
#'   \code{$logpost} only.
#' @param reject_fields Named list of extra fields carried on every rejected
#'   draw (\code{list(regime_path = NULL)} for the PKF factory).
#' @param pass_dots Return \code{function(theta, ...)} rather than
#'   \code{function(theta)}; only the TPF closure (which takes \code{U_list})
#'   needs it.
#' @return A closure \code{function(theta)} (or \code{function(theta, ...)})
#'   returning \code{list(logpost, loglik, logprior, <extra>)}.
#' @noRd
.make_posterior_closure <- function(model, data, prior_spec, obs_vars,
                                    compiled,
                                    loglik_fn,
                                    solve_fn          = NULL,
                                    sys_cache         = NULL,
                                    power             = 1,
                                    warm_start        = TRUE,
                                    warm_retry        = TRUE,
                                    stationarity      = "spectral",
                                    needs_me_floor    = TRUE,
                                    system_prior      = NULL,
                                    system_prior_fn   = NULL,
                                    system_prior_mode = c("lp", "extra"),
                                    reject_fields     = NULL,
                                    pass_dots         = FALSE,
                                    ## TRUE only for the Gaussian Kalman
                                    ## adapter, whose kalman_filter() call
                                    ## subtracts observation_trends.
                                    obs_trends_ok     = FALSE,
                                    ## TRUE only where solve_fn is a pure
                                    ## function of (ss, params): the draw's
                                    ## solve is then reused while the
                                    ## structural parameters do not move.
                                    structural_cache  = FALSE) {
  if (!isTRUE(obs_trends_ok))
    .refuse_obs_trends(model, "this likelihood")

  ## Force everything the closure reads: without this they stay unevaluated
  ## promises pointing at the adapter's frame, and a mirai daemon shipped the
  ## closure before its first call would fail to resolve them.
  force(model); force(data); force(prior_spec); force(obs_vars)
  force(compiled); force(loglik_fn); force(power)
  force(system_prior); force(reject_fields)

  system_prior_mode <- match.arg(system_prior_mode)
  ## `sys_cache = FALSE` means "this branch never calls
  ## extract_system_matrices_fast()" (pruned3 goes straight to
  ## solve_perturbation(order = 3L)); building the cache would be factory-time
  ## work for nothing, so honour the opt-out rather than always paying it.
  if (isFALSE(sys_cache)) sys_cache <- NULL
  else if (is.null(sys_cache)) sys_cache <- cache_system_structure(compiled)
  if (is.null(solve_fn)) {
    stat_mode <- stationarity
    solve_fn  <- function(model, compiled, sys_cache, ss, params, theta)
      .posterior_solve1(model, compiled, sys_cache, ss, params, stat_mode)
  }
  if (is.null(system_prior_fn))
    system_prior_fn <- function(spec, theta, sol, params, res)
      .eval_system_priors(
        spec,
        list(theta   = theta,
             model   = model,
             dr      = res$dr %||% sol$dr,
             Sigma_e = res$Sigma_e %||%
               .get_shock_cov(model, model$varexo_names, params),
             params  = params))

  ## ---- Per-closure mutable state ----------------------------------------
  ## Steady-state warm start: consecutive MCMC proposals are close in
  ## parameter space, so the previously converged steady state is an excellent
  ## Newton seed. Per-closure (not global), so each parallel chain keeps its
  ## own. Reset to NULL whenever a solve fails, so a misleading warm guess
  ## cannot poison every subsequent draw.
  ss_warm <- NULL
  ## Structural-parameter cache (see .structural_cache_new): one entry, owned
  ## by this closure.
  ## options(dynhr.structural_cache = FALSE), read when the closure is built,
  ## turns it off (a diagnostic switch: the cached and the re-solved draws are
  ## bit-identical, which the cache's own tests compare against).
  struct_cache <- if (isTRUE(structural_cache) && warm_retry &&
                      isTRUE(getOption("dynhr.structural_cache", TRUE)))
    .structural_cache_new(model) else NULL
  ## me-floor hazard guard: warn at most once per closure, not once per draw.
  .me_floor_checked <- FALSE

  ## Build a returned list: the three mandatory fields then the extras, in
  ## order. Extras are assigned one at a time through `[<-` with a length-1
  ## list so that a NULL extra ($regime_path on a rejected PKF draw) is STORED
  ## as NULL rather than dropping the element -- `out[[nm]] <- NULL` removes.
  ##
  ## Every returned list also carries the attribute "posterior_parts" (read it
  ## with .posterior_parts()): the target's COMPONENTS, which the returned
  ## fields alone cannot always recover -- in "lp" mode $logprior already
  ## contains the system prior, and $logpost has `power` and (on an
  ## infeasible_penalty reject) a finite optimiser penalty folded in. SMC reads
  ## the parts so it can draw stage 0 from the PARAMETER prior and temper
  ##   phi(theta) = log_sysprior + power * loglik
  ##. An attribute rather than a field because the field set is
  ## a pinned contract (test-posterior-closure-parity.R).
  .emit <- function(logpost, loglik, logprior, extra,
                    lp_param = logprior, sp = 0) {
    out <- list(logpost = logpost, loglik = loglik, logprior = logprior)
    for (nm in names(extra)) out[nm] <- list(extra[[nm]])
    attr(out, "posterior_parts") <- c(logprior_param = unname(lp_param),
                                      log_sysprior   = unname(sp),
                                      loglik         = unname(loglik),
                                      power          = unname(power))
    out
  }

  .reject <- function(lp, r = NULL, sp = 0) {
    .emit(if (is.null(r)) -Inf else r$logpost,
          if (is.null(r)) -Inf else r$loglik,
          lp,
          c(reject_fields, if (is.null(r)) NULL else r$extra),
          sp = sp)
  }

  ## Theta is read BY NAME against the prior spec (.theta_by_name(),
  ## R/posterior.R): an unnamed theta is taken in prior_spec$name order, a
  ## mis-named one is a classed error.
  spec_names <- prior_spec$name
  inner <- function(theta, ...) {
    theta <- .theta_by_name(theta, spec_names)
    lp <- log_prior(theta, prior_spec)
    if (!is.finite(lp)) return(.reject(lp))

    params <- .apply_theta_to_params(model, theta)

    ## ---- Steady state + solve (reused when only covariance-side parameters
    ## moved since the previous feasible draw) -------------------------------
    hit <- NULL
    if (!is.null(struct_cache)) {
      key <- .structural_cache_key(struct_cache, params)
      hit <- .structural_cache_lookup(struct_cache, key)
    }
    if (!is.null(hit)) {
      ss_result <- hit$ss_result
      params    <- .structural_cache_params(hit, params)
      sol       <- hit$sol
    } else {
      params_in <- params
      ss_result <- solve_steady_state(model, compiled, params,
                                      y0 = if (warm_start) ss_warm else NULL,
                                      verbose = FALSE)
      if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
        ## The warm guess may have been misleading (a large proposal jump);
        ## retry once from the cold initval-based guess before giving up.
        if (warm_start && warm_retry && !is.null(ss_warm))
          ss_result <- solve_steady_state(model, compiled, params,
                                          verbose = FALSE)
        if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
          if (warm_start) ss_warm <<- NULL
          return(.reject(lp))
        }
      }
      if (warm_start) ss_warm <<- ss_result$ss

      ## Re-derive any steady_state_model-computed parameter so the
      ## linearization point uses the consistent (not stale) p_c. Without this
      ## the dynamic system is built at an invalid steady state with the wrong
      ## p_c, silently biasing the loglik.
      params <- ss_result$params %||% params

      ## ---- Solve ----------------------------------------------------------
      sol <- solve_fn(model, compiled, sys_cache, ss_result$ss, params, theta)
      if (.is_posterior_reject(sol)) return(.reject(lp, sol))
      if (!is.null(struct_cache))
        .structural_cache_store(struct_cache, key, params_in, params,
                                ss_result, sol)
    }

    ## ---- Likelihood -------------------------------------------------------
    res <- loglik_fn(sol = sol, params = params, ss = ss_result$ss,
                     theta = theta,
                     me_floor_check = !.me_floor_checked &&
                       isTRUE(getOption("dynhr.me_floor_check", TRUE)),
                     ...)
    ## Latched AFTER the call and regardless of its outcome: the guard is
    ## about having HAD the chance to warn once, not about success.
    if (needs_me_floor) .me_floor_checked <<- TRUE
    if (.is_posterior_reject(res)) return(.reject(lp, res))
    loglik <- res$loglik

    ## ---- System prior -----------------------------------------------------
    sp_lp    <- NULL
    sp_part  <- 0
    lp_param <- lp
    if (!is.null(system_prior)) {
      sp_lp   <- system_prior_fn(system_prior, theta, sol, params, res)
      sp_part <- sp_lp
      if (identical(system_prior_mode, "lp")) {
        if (!is.finite(sp_lp))
          return(.reject(lp, .posterior_reject(logpost = -Inf,
                                               loglik  = loglik),
                         sp = sp_lp))
        lp    <- lp + sp_lp
        sp_lp <- NULL
      }
    }

    ## ---- Tempering --------------------------------------------------------
    ## `power` tempers the LIKELIHOOD only: the parameter prior and the system
    ## prior are prior-side and stay untempered, and $loglik always carries the
    ## RAW likelihood so SMC tempering and marginal-likelihood estimators see
    ## the true value. The two association orders below are the ones the
    ## respective factories used and are pinned to the last ulp by
    ## test-posterior-closure-parity.R -- do not "simplify" them into one.
    logpost <- if (is.null(sp_lp)) power * loglik + lp
               else                lp + power * loglik + sp_lp

    .emit(logpost, loglik, lp, res$extra, lp_param = lp_param, sp = sp_part)
  }

  if (pass_dots) inner else function(theta) inner(theta)
}


#' Target components of one log-posterior evaluation
#'
#' Reads the \code{"posterior_parts"} attribute that
#' \code{.make_posterior_closure()} attaches to every returned list:
#' \code{logprior_param} (the PARAMETER prior only), \code{log_sysprior} (the
#' system-prior log density, 0 when there is none), \code{loglik} (the RAW
#' likelihood) and \code{power} (the tempering exponent zeta), so that
#' \code{logpost = logprior_param + log_sysprior + power * loglik} on every
#' feasible draw whatever the factory's \code{system_prior_mode}.
#'
#' A closure that does not attach the attribute (a user-written
#' \code{function(theta) list(logpost, loglik, logprior)}) is read with the
#' historical contract: \code{$logprior} is the parameter prior, the tempered
#' score is \code{$loglik}, no system prior, \code{power = 1}.
#'
#' @param r A list returned by a log-posterior closure.
#' @return Named numeric \code{c(logprior_param, log_sysprior, loglik, power)}.
#' @noRd
.posterior_parts <- function(r) {
  p <- attr(r, "posterior_parts", exact = TRUE)
  if (!is.null(p)) return(p)
  c(logprior_param = if (is.null(r$logprior)) 0 else unname(r$logprior),
    log_sysprior   = 0,
    loglik         = if (is.null(r$loglik)) -Inf else unname(r$loglik),
    power          = 1)
}
