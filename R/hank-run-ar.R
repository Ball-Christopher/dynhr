## R/hank-run-ar.R
## --------------------------------------------------------------------------
## Runner integration for the exact-AR(1) HANK likelihood: a mode-finding and
## posterior-sampling entry point over STRUCTURAL parameters (which move the
## sequence-space model) together with the per-shock (rho, sigma) block.
##
## WHY THIS IS NOT A FLAG ON run_mode_finding(). That runner is a decision-rule
## DSGE runner end to end -- it takes `solved$model` / `solved$compiled` from a
## compiled .mod, builds its posterior with make_log_posterior(), and provisions
## mirai daemons that recompile that posterior per worker. A hank_model() has
## none of those objects: it is a sequence-space DAG whose "solution" is the
## Jacobian G. So the integration is a sibling runner over the same generic
## optimizer core (.run_mode_finding, R/mode-orchestrate.R) and the same
## samplers, not a new branch inside the DSGE one.
##
## THE DESIGN REQUIREMENT (0.9.0.0024, briefs/21 sections 10-11). The exact-AR
## structural gradient is 21x an FD gradient ONLY when three things are true
## across calls, and a naive runner breaks all three:
##
##   1. ONE model rebuild per distinct STRUCTURAL theta, shared by the
##      likelihood and the gradient. The rebuild is the dominant term (0.46 s
##      at 11 shocks x 12 observables, T_h = 200, against 0.008 s for the DAG
##      propagation), so a runner that calls model_fn() once in the objective
##      and again in the gradient has already halved the gain -- and a rho- or
##      sigma-only move must not rebuild at all.
##   2. ONE `cache` environment for the whole run. The stacked gather index
##      survives a structural move (only the autocovariance slabs are
##      invalidated) and costs more to rebuild than every slab put together.
##   3. ONE kernel adjoint per gradient, shared between the structural score
##      and the (rho, sigma) scores. They contract the SAME dl/dTheta; computing
##      it twice is a second adjoint pass for nothing. This is what
##      hank_loglik_ar_structural_grad()'s `score` argument exists for.
##
## and, on the exact-derivative route, hank_dtheta_fn() must be SEEDED with the
## model already built (its own memo would otherwise rebuild it), which is why
## the closure is reconstructed on a structural move rather than once: its memo
## is keyed on theta, so a fresh seeded closure loses nothing that is not
## already stale, and gains the rebuild.
##
## REPRESENTABILITY. Every entry point here carries the `boundary` guard (see
## .hank_ar_boundary_gate in R/hank-kalman-ar.R). Mode-finding defaults to
## "reject": an optimizer following the score has no reason not to walk out
## into the terminal-boundary-contaminated region, where the likelihood is
## wrong rather than imprecise, and a mode found out there is not a mode.
## --------------------------------------------------------------------------


#' Exact-AR estimation target for a HANK model (structural + shock parameters)
#'
#' Builds the matched pair of closures a mode-finder or sampler needs --
#' \code{log_post_fn(theta)} and \code{grad_fn(theta)} -- over a parameter
#' vector that may contain STRUCTURAL parameters (which rebuild the
#' sequence-space model) alongside the per-shock \code{rho_<shock>} and
#' \code{sigma_<shock>} block. The two closures share one model memo, one
#' likelihood \code{cache}, and one kernel adjoint per gradient; see the file
#' header for why each of those is load-bearing rather than an optimization.
#'
#' With \code{model_fn = NULL} this reduces to the shock-only target and is
#' equivalent to pairing \code{\link{make_log_posterior_hank}}
#' (\code{likelihood = "exact_ar"}) with
#' \code{\link{make_posterior_grad_hank_ar}}, except that the pair shares its
#' cache.
#'
#' @section Priors: Normal on each \code{rho} (soft-truncated to
#'   \code{(-1, 1)}), half-Normal on each \code{sigma}, and Normal on each
#'   structural parameter, optionally truncated by
#'   \code{structural_lower}/\code{structural_upper}. This matches
#'   \code{\link{make_log_posterior_hank}}'s convention. For anything richer,
#'   take \code{log_post_fn}/\code{grad_fn} and add your own log-prior --
#'   \code{$loglik} is returned separately for exactly that.
#'
#' @param Y \code{T_data x n_obs} matrix of demeaned observations.
#' @param observables Character vector of observable names, in the column order
#'   of \code{Y}.
#' @param model A \code{\link{hank_model}}. Required when \code{model_fn} is
#'   \code{NULL}; otherwise optional, and used to seed the memo (so the first
#'   evaluation at \code{structural} costs no rebuild).
#' @param model_fn \code{function(theta_structural)} returning a
#'   \code{\link{hank_model}}. Supply it to estimate structural parameters.
#' @param structural Named numeric vector of structural parameters: its names
#'   and order define that block of \code{theta}, and its VALUES are the Normal
#'   prior means. \code{NULL} (default) estimates the shock block only.
#' @param structural_sd Named numeric vector or scalar: Normal prior sd for
#'   each structural parameter. Required when \code{structural} is supplied.
#' @param structural_lower,structural_upper Optional named numeric vectors (or
#'   scalars) truncating the structural prior; draws outside return
#'   \code{-Inf}. Default unbounded.
#' @param dtheta_fn Optional \code{function(theta, param, shock_specs)} giving
#'   \eqn{d\Theta/d\theta_k} exactly. \code{"auto"} (default) builds one with
#'   \code{\link{hank_dtheta_fn}} when \code{model_fn} is supplied, which is
#'   the zero-rebuild route; \code{NULL} falls back to central differences of
#'   \code{model_fn}.
#' @param verify Passed to \code{\link{hank_loglik_ar_structural_grad}} for the
#'   FIRST gradient only, then switched off: the check costs one central
#'   difference per parameter, which is worth paying once and never again
#'   inside a loop (a wrong \code{dtheta_fn} is wrong at every theta, so
#'   re-checking buys nothing).
#' @param q,me_var,me_sd As in \code{\link{make_log_posterior_hank}}.
#' @param prior_rho_mean,prior_rho_sd,prior_sigma_sd As in
#'   \code{\link{make_log_posterior_hank}}.
#' @param rho_method Passed to \code{\link{hank_loglik_ar_grad}}.
#' @param boundary,boundary_tol Representability guard, as in
#'   \code{\link{make_log_posterior_hank}}. Default \code{"warn"} here;
#'   \code{\link{hank_run_mode_finding}} defaults to \code{"reject"}.
#'
#' @return An object of class \code{hank_ar_target}: a list with
#'   \code{log_post_fn}, \code{grad_fn}, \code{theta_names},
#'   \code{structural_names}, \code{shock_names}, \code{prior_spec} (in the
#'   \code{extract_prior_spec} column layout, so the generic mode-finder and
#'   \code{build_param_transform} accept it -- and under
#'   \code{boundary = "reject"} its \code{rho} bounds are the representability
#'   cap rather than \eqn{\pm 1}, so the box-constrained optimizer stages and
#'   the eta transform respect the guard instead of discovering a
#'   \code{-Inf} cliff by walking off it), \code{cache}, and \code{stats}
#'   (an environment counting \code{rebuilds} / \code{reuses} -- the direct
#'   check that the sharing above is actually happening).
#' @seealso \code{\link{hank_run_mode_finding}}, \code{\link{hank_run_estimation}},
#'   \code{\link{hank_loglik_ar_structural_grad}}, \code{\link{hank_dtheta_fn}}
#' @export
hank_ar_target <- function(Y, observables, model = NULL, model_fn = NULL,
                           structural = NULL, structural_sd = NULL,
                           structural_lower = NULL, structural_upper = NULL,
                           dtheta_fn = "auto", verify = TRUE,
                           q = NULL, me_var = 0, me_sd = NULL,
                           prior_rho_mean = 0.5, prior_rho_sd = 0.3,
                           prior_sigma_sd = 0.05,
                           rho_method = c("fd_slab", "adjoint"),
                           boundary = c("warn", "reject", "ignore"),
                           boundary_tol = 1e-3) {
  rho_method <- match.arg(rho_method)
  boundary   <- match.arg(boundary)
  if (is.null(model) && is.null(model_fn))
    stop("hank_ar_target: supply `model`, `model_fn`, or both.")
  if (!is.null(model) && !inherits(model, "hank_model"))
    stop("hank_ar_target: `model` must be a hank_model() object.")
  if (!is.null(model_fn) && !is.function(model_fn))
    stop("hank_ar_target: `model_fn` must be a function of the structural ",
         "parameter vector returning a hank_model().")
  if (!is.null(structural)) {
    if (is.null(model_fn))
      stop("hank_ar_target: `structural` parameters need `model_fn` -- there ",
           "is no way to move them without rebuilding the model.")
    if (!is.numeric(structural) || is.null(names(structural)) ||
        any(!nzchar(names(structural))))
      stop("hank_ar_target: `structural` must be a NAMED numeric vector ",
           "(its names define that block of theta and its values are the ",
           "prior means).")
    if (is.null(structural_sd))
      stop("hank_ar_target: `structural_sd` is required when `structural` ",
           "is supplied (there is no defensible default prior width for a ",
           "structural parameter).")
  }
  Y <- as.matrix(Y)
  ## One model is needed up front for the shock names, T_h and the observable
  ## check. If the caller did not hand one over, build it ONCE at `structural`
  ## and seed the memo with it -- never build a throwaway.
  if (is.null(model)) model <- model_fn(structural)
  if (!inherits(model, "hank_model"))
    stop("hank_ar_target: `model_fn` must return a hank_model() object.")
  missing_obs <- setdiff(observables, names(model$G))
  if (length(missing_obs))
    stop("hank_ar_target: observable(s) not produced by model: ",
         paste(missing_obs, collapse = ", "))
  if (is.null(me_sd)) me_sd <- sqrt(me_var)

  rep_named <- function(x, nm, what) {
    if (is.null(x)) return(NULL)
    out <- if (is.null(names(x))) stats::setNames(rep(x, length.out = length(nm)), nm)
           else x[nm]
    if (any(is.na(out)))
      stop("hank_ar_target: `", what, "` is missing an entry for ",
           paste(nm[is.na(out)], collapse = ", "), ".")
    out
  }

  ## ---- parameter blocks and priors ---------------------------------------
  s_nm   <- names(structural)                       # may be NULL/length 0
  exo    <- model$exogenous
  rho_nm <- paste0("rho_", exo); sig_nm <- paste0("sigma_", exo)

  rho_mean <- rep_named(prior_rho_mean, exo, "prior_rho_mean")
  rho_sd   <- rep_named(prior_rho_sd, exo, "prior_rho_sd")
  sig_sd   <- rep_named(prior_sigma_sd, exo, "prior_sigma_sd")
  s_mean   <- structural
  s_sd     <- rep_named(structural_sd, s_nm, "structural_sd")
  s_lo     <- rep_named(structural_lower %||% -Inf, s_nm, "structural_lower")
  s_hi     <- rep_named(structural_upper %||%  Inf, s_nm, "structural_upper")

  theta_names <- c(s_nm, rho_nm, sig_nm)
  ## Under boundary = "reject" the representability bound IS the rho support,
  ## so put it in the prior spec rather than leaving (-1, 1) there and letting
  ## the optimizer discover a -Inf cliff by walking off it. This is what makes
  ## the box-constrained stages and the eta transform respect the guard instead
  ## of fighting it.
  rho_cap <- if (identical(boundary, "reject"))
    boundary_tol^(1 / model$T_h) else 1
  prior_spec <- data.frame(
    name = theta_names,
    distribution = c(rep("normal", length(s_nm)),
                     rep("normal", length(exo)),
                     rep("normal", length(exo))),
    p1 = c(unname(s_mean), unname(rho_mean), rep(0, length(exo))),
    p2 = c(unname(s_sd),   unname(rho_sd),   unname(sig_sd)),
    lower = c(unname(s_lo), rep(-rho_cap, length(exo)), rep(0, length(exo))),
    upper = c(unname(s_hi), rep( rho_cap, length(exo)), rep(Inf, length(exo))),
    stringsAsFactors = FALSE)

  ## ---- shared state ------------------------------------------------------
  ar_cache <- new.env(parent = emptyenv())
  boundary_state <- new.env(parent = emptyenv())
  stats_env <- new.env(parent = emptyenv())
  stats_env$rebuilds <- 0L; stats_env$reuses <- 0L
  st <- new.env(parent = emptyenv())
  ## Seeded with the model in hand at `structural`, so the first evaluation
  ## costs no rebuild.
  st$theta_s <- structural; st$model <- model; st$dt <- NULL
  st$verified <- !isTRUE(verify); st$announced <- FALSE

  ## "auto" -> a memoized exact dTheta closure SEEDED with this model (the
  ## zero-rebuild route); an explicit function -> used as given; NULL -> the
  ## finite-difference route inside hank_loglik_ar_structural_grad().
  make_dt <- function(theta_s, mod, shock_specs) {
    if (is.null(model_fn) || is.null(dtheta_fn)) return(NULL)
    if (identical(dtheta_fn, "auto"))
      hank_dtheta_fn(model_fn, shock_specs, observables,
                     model = mod, theta = theta_s)
    else dtheta_fn
  }

  ## The ONE rebuild rule: model_fn() is called only when the STRUCTURAL block
  ## actually moved. A rho/sigma-only step (every step of the shock block, and
  ## every likelihood/gradient pair at the same theta) reuses.
  base_model <- function(theta_s, shock_specs) {
    if (is.null(model_fn)) return(model)
    if (!is.null(st$model) && identical(st$theta_s, theta_s)) {
      stats_env$reuses <- stats_env$reuses + 1L
      return(st$model)
    }
    mod <- model_fn(theta_s)
    if (!inherits(mod, "hank_model"))
      stop("hank_ar_target: `model_fn` must return a hank_model() object.")
    stats_env$rebuilds <- stats_env$rebuilds + 1L
    st$theta_s <- theta_s; st$model <- mod
    ## Reconstruct the dTheta memo SEEDED with this model: its memo is keyed on
    ## theta, so nothing reusable is lost, and the seeding is what keeps the
    ## exact route at zero rebuilds per gradient.
    st$dt <- make_dt(theta_s, mod, shock_specs)
    mod
  }

  split_theta <- function(theta) {
    if (!all(theta_names %in% names(theta)))
      stop("hank_ar_target: theta must have entries ",
           paste(theta_names, collapse = ", "), ".")
    list(s = if (length(s_nm)) theta[s_nm] else NULL,
         rho = stats::setNames(theta[rho_nm], exo),
         sigma = stats::setNames(theta[sig_nm], exo))
  }

  log_prior_of <- function(p) {
    lp <- sum(stats::dnorm(p$rho, rho_mean, rho_sd, log = TRUE)) +
      sum(stats::dnorm(p$sigma, 0, sig_sd, log = TRUE) + log(2))
    if (length(s_nm))
      lp <- lp + sum(stats::dnorm(p$s, s_mean, s_sd, log = TRUE))
    lp
  }
  infeasible <- function(p)
    any(p$rho <= -1 | p$rho >= 1) || any(p$sigma <= 0) ||
      (length(s_nm) && (any(p$s < s_lo) || any(p$s > s_hi)))

  specs_of <- function(p) stats::setNames(
    lapply(exo, function(z) list(rho = p$rho[[z]], sigma = p$sigma[[z]])), exo)

  ## ---- the two closures --------------------------------------------------
  log_post_fn <- function(theta) {
    p <- split_theta(theta)
    if (infeasible(p))
      return(list(logpost = -Inf, loglik = NA_real_, logprior = -Inf))
    logprior <- log_prior_of(p)
    specs <- specs_of(p)
    mod <- base_model(p$s, specs)
    if (.hank_ar_boundary_gate(p$rho, mod$T_h, boundary, boundary_tol,
                               "hank_ar_target", boundary_state))
      return(list(logpost = -Inf, loglik = NA_real_, logprior = logprior))
    Theta_list <- .hank_theta_list(mod, specs, observables)
    loglik <- tryCatch(
      hank_loglik_ar(Y, Theta_list, rho = p$rho, sigma = p$sigma,
                     me_sd = me_sd, q = q, check_boundary = FALSE,
                     cache = ar_cache),
      error = function(e) -Inf)
    if (!is.finite(loglik))
      return(list(logpost = -Inf, loglik = loglik, logprior = logprior))
    list(logpost = loglik + logprior, loglik = loglik, logprior = logprior)
  }

  grad_full <- function(theta) {
    p <- split_theta(theta)
    if (infeasible(p))
      return(list(logpost = -Inf, loglik = NA_real_, logprior = -Inf,
                  grad = NULL))
    logprior <- log_prior_of(p)
    specs <- specs_of(p)
    mod <- base_model(p$s, specs)
    if (.hank_ar_boundary_gate(p$rho, mod$T_h, boundary, boundary_tol,
                               "hank_ar_target", boundary_state))
      return(list(logpost = -Inf, loglik = NA_real_, logprior = logprior,
                  grad = NULL))
    Theta_list <- .hank_theta_list(mod, specs, observables)

    ## The rho score needs Theta at a perturbed persistence; both routes reuse
    ## the SAME model, so neither costs a rebuild.
    irf_of <- function(z, path) {
      irf <- hank_model_irf(mod, stats::setNames(list(path), z))
      matrix(sapply(observables, function(o) irf[[o]]),
             mod$T_h, length(observables), dimnames = list(NULL, observables))
    }
    theta_fn  <- function(z, r) irf_of(z, r^(seq_len(mod$T_h) - 1L))
    dtheta_rho <- function(z, r) {
      s <- seq_len(mod$T_h) - 1L
      irf_of(z, ifelse(s == 0, 0, s * r^pmax(s - 1L, 0)))
    }

    ## ONE adjoint pass: "theta" is requested alongside the shock blocks so the
    ## structural contraction below reuses it instead of computing a second.
    wrt <- if (length(s_nm)) c("sigma", "rho", "theta") else c("sigma", "rho")
    sc <- tryCatch(
      hank_loglik_ar_grad(Y, Theta_list, rho = p$rho, sigma = p$sigma,
                          me_sd = me_sd, q = q, check_boundary = FALSE,
                          cache = ar_cache, wrt = wrt, theta_fn = theta_fn,
                          rho_method = rho_method, dtheta_fn = dtheta_rho),
      error = function(e) NULL)
    if (is.null(sc) || !is.finite(sc$loglik))
      return(list(logpost = -Inf,
                  loglik = if (is.null(sc)) NA_real_ else sc$loglik,
                  logprior = logprior, grad = NULL))

    grad <- stats::setNames(numeric(length(theta_names)), theta_names)
    grad[rho_nm] <- sc$rho[exo]   - (p$rho - rho_mean) / rho_sd^2
    grad[sig_nm] <- sc$sigma[exo] - p$sigma / sig_sd^2

    if (length(s_nm)) {
      if (is.null(st$dt)) st$dt <- make_dt(p$s, mod, specs)
      ## The memo inside st$dt is shock_specs-INDEPENDENT (see hank_dtheta_fn),
      ## so the CURRENT specs are passed per call rather than baked in: a rho
      ## move costs the cheap re-application, not the model rebuild.
      dfun <- if (is.null(st$dt)) NULL else
        function(th, k) st$dt(th, k, specs)
      call_ssg <- function() hank_loglik_ar_structural_grad(
        Y, model_fn = model_fn, theta = p$s, observables = observables,
        rho = p$rho, sigma = p$sigma, me_sd = me_sd, q = q,
        dtheta_fn = dfun, verify = !st$verified, cache = ar_cache,
        check_boundary = FALSE, boundary_tol = boundary_tol,
        Theta_list = Theta_list, score = sc)
      ## That function announces once PER CALL which parameters it is trusting
      ## unverified -- correct for a direct caller, thousands of identical
      ## messages inside a sampler. Let the first one through (it is genuinely
      ## worth seeing once) and silence the repeats.
      ssg <- if (st$announced) suppressMessages(call_ssg()) else call_ssg()
      st$verified <- TRUE; st$announced <- TRUE
      grad[s_nm] <- ssg$structural[s_nm] - (p$s - s_mean) / s_sd^2
    }

    list(logpost = sc$loglik + logprior, loglik = sc$loglik,
         logprior = logprior, grad = grad)
  }

  grad_fn <- function(theta) {
    g <- grad_full(theta)$grad
    if (is.null(g)) stats::setNames(rep(NA_real_, length(theta_names)),
                                    theta_names) else g
  }

  structure(list(log_post_fn = log_post_fn, grad_fn = grad_fn,
                 grad_full = grad_full, theta_names = theta_names,
                 structural_names = s_nm, shock_names = exo,
                 prior_spec = prior_spec, cache = ar_cache, stats = stats_env,
                 boundary = boundary, boundary_tol = boundary_tol),
            class = "hank_ar_target")
}


#' Posterior mode-finding for a HANK model under the exact-AR likelihood
#'
#' The \code{likelihood = "exact_ar"} counterpart of
#' \code{\link{run_mode_finding}} for sequence-space HANK models, which that
#' runner rejects (it is a decision-rule/\code{.mod} runner end to end; see the
#' header of \code{R/hank-run-ar.R}). It drives the SAME optimizer core and
#' returns a result object \code{\link{hank_run_estimation}} consumes.
#'
#' @param target A \code{\link{hank_ar_target}}, or the arguments to build one
#'   (passed through \code{...}).
#' @param theta_init Named starting vector over \code{target$theta_names}.
#' @param method Optimizer, passed to the generic mode-finder
#'   (\code{"newrat"} default; also \code{"combined"}, \code{"nelder"},
#'   \code{"cmaes"}, \code{"nmkb"}, \code{"jade"}).
#' @param n_iter Iteration budget.
#' @param use_grad Use \code{target$grad_fn} for the gradient-capable stages
#'   (default \code{TRUE} -- it is the whole point of the exact-AR score).
#' @param transform_params Sample/optimize in unconstrained eta-space
#'   (default \code{TRUE}, matching \code{\link{run_mode_finding}}).
#' @param boundary Overrides the target's representability guard for the
#'   optimization; default \code{"reject"}, because an optimizer following the
#'   score has no reason not to walk out into the contaminated region and a
#'   mode found out there is not a mode. Pass \code{NULL} to keep the target's
#'   own setting.
#' @param verbose Print progress.
#' @param ... Passed to \code{\link{hank_ar_target}} when \code{target} is not
#'   already one.
#'
#' @return A list of class \code{hank_ar_mode}: \code{theta_mode},
#'   \code{logpost}, \code{loglik}, \code{target}, \code{log_post_fn},
#'   \code{grad_fn}, \code{prior_spec}, \code{convergence}, \code{iterations},
#'   \code{method}, and \code{rebuilds} (model rebuilds actually paid).
#' @seealso \code{\link{hank_ar_target}}, \code{\link{hank_run_estimation}}
#' @export
hank_run_mode_finding <- function(target, theta_init, method = "newrat",
                                  n_iter = 2000L, use_grad = TRUE,
                                  transform_params = TRUE,
                                  boundary = "reject", verbose = TRUE, ...) {
  if (!inherits(target, "hank_ar_target"))
    stop("hank_run_mode_finding: `target` must be a hank_ar_target() ",
         "(build one with hank_ar_target(Y, observables, ...)).")
  if (!is.null(boundary) && !identical(boundary, target$boundary))
    stop("hank_run_mode_finding: the supplied `target` was built with ",
         "boundary = \"", target$boundary, "\" but this call asks for \"",
         boundary, "\". The guard belongs to the closure (it decides what is ",
         "in the prior support), so rebuild the target with the boundary you ",
         "want, or pass boundary = NULL to keep the target's.")
  if (!is.numeric(theta_init) || is.null(names(theta_init)))
    stop("hank_run_mode_finding: `theta_init` must be a named numeric vector.")
  miss <- setdiff(target$theta_names, names(theta_init))
  if (length(miss))
    stop("hank_run_mode_finding: `theta_init` is missing ",
         paste(miss, collapse = ", "), ".")
  theta_init <- theta_init[target$theta_names]

  lp0 <- target$log_post_fn(theta_init)
  if (!is.finite(lp0$logpost))
    stop("hank_run_mode_finding: the starting point has log-posterior ",
         format(lp0$logpost), " -- an optimizer cannot leave an infeasible ",
         "start. Check the representability guard (boundary = \"",
         target$boundary, "\", |rho|^T_h <= ", format(target$boundary_tol),
         ") and the prior bounds.")

  transform <- if (isTRUE(transform_params))
    build_param_transform(target$prior_spec, target$theta_names) else NULL

  res <- .run_mode_finding(
    log_post_fn = target$log_post_fn, theta_init = theta_init,
    prior_spec = target$prior_spec, nm_maxit = n_iter, method = method,
    transform = transform,
    grad_fn = if (isTRUE(use_grad)) target$grad_fn else NULL,
    verbose = verbose)

  theta_mode <- res$theta_mode
  names(theta_mode) <- target$theta_names
  at_mode <- target$log_post_fn(theta_mode)
  structure(list(theta_mode = theta_mode, logpost = at_mode$logpost,
                 loglik = at_mode$loglik, target = target,
                 log_post_fn = target$log_post_fn, grad_fn = target$grad_fn,
                 prior_spec = target$prior_spec,
                 convergence = res$convergence, iterations = res$iterations,
                 method = method, rebuilds = target$stats$rebuilds),
            class = "hank_ar_mode")
}


#' Posterior sampling for a HANK model under the exact-AR likelihood
#'
#' The \code{likelihood = "exact_ar"} counterpart of
#' \code{\link{run_full_estimation}} / \code{\link{run_posterior_estimation}}:
#' runs a sampler on the mode result from
#' \code{\link{hank_run_mode_finding}}, reusing that run's target -- and
#' therefore its model memo and likelihood cache -- rather than rebuilding a
#' posterior.
#'
#' @param mode_result A \code{hank_ar_mode} from
#'   \code{\link{hank_run_mode_finding}}.
#' @param method \code{"RWMH"} (default) or \code{"NUTS"}. NUTS uses the exact
#'   score (\code{target$grad_fn}); that is the case the whole exact-AR
#'   gradient stack exists for.
#' @param n_draws,n_burn Sampler lengths (NUTS reads \code{n_burn} as warmup).
#' @param Sigma_prop RWMH proposal covariance; default a diagonal built from
#'   the prior sds, scaled by \code{0.1^2}.
#' @param transform_params Sample in unconstrained eta-space (default
#'   \code{TRUE}); the change-of-variables Jacobian IS included here, unlike in
#'   mode-finding.
#' @param seed Optional RNG seed.
#' @param verbose Print progress.
#' @param ... Passed to the sampler.
#'
#' @return The sampler's own result object, with \code{target} and
#'   \code{theta_mode} attached.
#' @seealso \code{\link{hank_run_mode_finding}}, \code{\link{hank_ar_target}}
#' @export
hank_run_estimation <- function(mode_result, method = c("RWMH", "NUTS"),
                                n_draws = 2000L, n_burn = 1000L,
                                Sigma_prop = NULL, transform_params = TRUE,
                                seed = NULL, verbose = TRUE, ...) {
  method <- match.arg(method)
  if (!inherits(mode_result, "hank_ar_mode"))
    stop("hank_run_estimation: `mode_result` must come from ",
         "hank_run_mode_finding().")
  if (!is.null(seed)) set.seed(seed)
  target <- mode_result$target
  theta0 <- mode_result$theta_mode
  transform <- if (isTRUE(transform_params))
    build_param_transform(target$prior_spec, target$theta_names) else NULL

  if (method == "RWMH") {
    if (is.null(Sigma_prop)) {
      sd0 <- target$prior_spec$p2
      Sigma_prop <- diag((0.1 * sd0)^2, nrow = length(sd0))
      dimnames(Sigma_prop) <- list(target$theta_names, target$theta_names)
    }
    out <- rwmh(target$log_post_fn, theta0, Sigma_prop, n_draws = n_draws,
                n_burn = n_burn, verbose = verbose, transform = transform, ...)
  } else {
    out <- dynhr_nuts(target$log_post_fn, theta0, n_draws = n_draws,
                      n_warmup = n_burn, grad_fn = target$grad_fn,
                      verbose = verbose, transform = transform, ...)
  }
  out$target <- target
  out$theta_mode <- theta0
  out
}
