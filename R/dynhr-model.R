## R/dynhr-model.R
## --------------------------------------------------------------------------
## `dynhr_model`: a pipeline OBJECT over the functional API in R/api.R.
##
## The functional facade (`prior_spec` -> `make_posterior` -> `find_mode` ->
## `mcmc`/`nuts`/`smc`/`dime`) requires the caller to re-thread the same five
## or six pieces -- model, compiled, data, obs_vars, prior_spec, me_variance --
## through every step. `dynhr_model()` bundles them once; the `dm_*()` verbs
## then take that object as their FIRST argument (so they pipe) and return a
## NEW object with the step's result stored in a slot.
##
## Contract: this file is PURELY ADDITIVE. Every verb is a re-threading of
## stored fields into the existing exported function, with no numerical logic
## of its own, so a `dm_*()` result is bit-identical to the functional call it
## wraps (pinned by tests/testthat/test-dynhr-model-golden.R).
##
## Design notes
##   * S3, not R6: R6 is not an Import of dynhr and must not become one.
##   * Compilation is EAGER (in the constructor), not lazy: an object that has
##     been built has already paid for `compile_model()`, so no verb can
##     surprise the caller with a multi-second compile. Use `max_order=` to
##     compile deep enough for the perturbation order you intend to solve at,
##     or pass a `compiled=` you built yourself.
##   * No verb silently runs its own prerequisite. A missing step is an error
##     that names the verb to call, so a pipeline never quietly re-derives a
##     piece the caller thought it had configured.
## --------------------------------------------------------------------------


## Slots a verb writes, in pipeline order, for print()/summary().
.dm_steps <- c(solution = "dm_solve", log_post = "dm_posterior",
               mode = "dm_mode", chains = "dm_sample",
               diagnostics = "dm_diagnostics", irf = "dm_irf",
               forecast = "dm_forecast")


#' Fail with a message naming the missing pipeline step
#'
#' @param verb   Name of the verb that was called
#' @param slot   Name of the empty slot
#' @param needed Name of the verb that fills it
#' @return Never returns; signals an error
#' @noRd
.dm_require <- function(verb, slot, needed) {
  stop(sprintf("%s(): `$%s` is empty -- run %s() on this object first.",
               verb, slot, needed), call. = FALSE)
}


#' Fail when a field the constructor should have stored is absent
#'
#' @param verb  Name of the verb that was called
#' @param field Name of the missing constructor field
#' @return Never returns; signals an error
#' @noRd
.dm_require_field <- function(verb, field) {
  stop(sprintf("%s(): `$%s` is empty -- supply %s= to dynhr_model().",
               verb, field, field), call. = FALSE)
}


#' Check that an object is a dynhr_model
#'
#' @param dm   Candidate object
#' @param verb Calling verb, for the error message
#' @return `dm`, invisibly
#' @noRd
.dm_check <- function(dm, verb) {
  if (!inherits(dm, "dynhr_model"))
    stop(sprintf("%s(): `dm` must be a dynhr_model object (see dynhr_model()).",
                 verb), call. = FALSE)
  invisible(dm)
}


#' Store a verb's result on the object and return the new object
#'
#' @param dm    A `dynhr_model`
#' @param slot  Slot name to write
#' @param value Value to store
#' @return The updated `dynhr_model`
#' @noRd
.dm_set <- function(dm, slot, value) {
  dm[[slot]] <- value
  dm
}


# ============================================================================
# Constructor
# ============================================================================

#' Bundle a DSGE model, its data and its estimation settings into one object
#'
#' `dynhr_model()` packages everything a dynhr estimation pipeline needs --
#' the parsed model, its compiled derivatives, the observation matrix, the
#' observed-variable names, the prior specification, and the likelihood
#' settings -- into a single object.  The `dm_*()` verbs
#' (\code{\link{dm_solve}}, \code{\link{dm_posterior}}, \code{\link{dm_mode}},
#' \code{\link{dm_sample}}, \code{\link{dm_diagnostics}},
#' \code{\link{dm_irf}}, \code{\link{dm_forecast}}) each take that object as
#' their first argument and return a NEW object with the step's result stored
#' in a slot, so a whole pipeline reads as one chain instead of six calls that
#' re-thread the same arguments.
#'
#' The object is a convenience layer only: each verb forwards the stored
#' fields to the corresponding exported function unchanged, and the raw result
#' is stored as-is (a \code{DecisionRules*} from
#' \code{\link{solve_perturbation}}, a closure from
#' \code{\link{make_posterior}}, a \code{\link{dynhr_chains}} from
#' \code{\link{mcmc}}, ...).  Results are therefore identical to the
#' functional pipeline, bit for bit.
#'
#' @section Compilation:
#' Compilation is **eager**: unless `compiled` is supplied, the constructor
#' calls \code{\link{compile_model}} immediately with `max_order`.  Set
#' `max_order = 2L` (or higher) if you intend to call
#' \code{dm_solve(order = 2L)} or a higher-order likelihood; the solvers fail
#' loudly rather than silently consuming zero tensors when the compiled
#' derivatives are too shallow.
#'
#' @section Prerequisites:
#' No verb runs another verb for you.  Calling \code{\link{dm_mode}} before
#' \code{\link{dm_posterior}} is an error naming the missing step.  The one
#' documented default is \code{\link{dm_sample}}'s starting point: with no
#' `$mode` on the object it starts from the prior means rather than erroring
#' (SMC and DIME need no starting point at all).
#'
#' @param model       A `dynhr_mod` object from \code{\link{parse_mod}}, or a
#'   path to a `.mod` file (parsed on the spot).
#' @param data        Observation matrix (\eqn{T \times n_{\text{obs}}}) whose
#'   columns match `obs_vars` -- the same orientation contract as
#'   \code{\link{make_posterior}}.  Optional: a model object with no data
#'   still supports \code{\link{dm_solve}} and \code{\link{dm_irf}}.
#' @param obs_vars    Character vector of observed variable names.  Defaults
#'   to the model's `varobs` declaration when it has one.
#' @param compiled    A `dynhr_compiled` from \code{\link{compile_model}}.
#'   When `NULL` (default) the constructor compiles the model itself.
#' @param prior_spec  Prior spec data.frame from \code{\link{prior_spec}}.
#'   When `NULL` (default) it is extracted from the model's
#'   `estimated_params` block if there is one, and left empty otherwise.
#' @param me_variance Measurement-error variance added to the diagonal of the
#'   observation noise (default 0).
#' @param likelihood  Likelihood type, as in \code{\link{make_posterior}}
#'   (default `"gaussian"`).
#' @param params      Named numeric vector of parameter values used by
#'   \code{\link{dm_solve}} and \code{\link{dm_irf}}.  Defaults to the model's
#'   calibrated `param_values`.
#' @param max_order   Perturbation order to compile derivatives for
#'   (default `1L`); ignored when `compiled` is supplied.
#' @param verbose     Print parse/compile progress (default `FALSE`).
#' @param ...         Extra arguments for the likelihood, stored on the object
#'   and forwarded by \code{\link{dm_posterior}} to
#'   \code{\link{make_posterior}} (e.g. `order` for
#'   `likelihood = "cumulant"`).
#'
#' @return An object of class `dynhr_model`: a list with elements `model`,
#'   `compiled`, `data`, `obs_vars`, `prior_spec`, `me_variance`,
#'   `likelihood`, `lik_args`, `params`, and the initially-`NULL` result slots
#'   `steady`, `solution`, `log_post`, `mode`, `chains`, `diagnostics`, `irf`
#'   and `forecast`.
#'
#' @seealso \code{\link{dm_solve}}, \code{\link{dm_posterior}},
#'   \code{\link{dm_mode}}, \code{\link{dm_sample}},
#'   \code{\link{dm_diagnostics}}, \code{\link{dm_irf}},
#'   \code{\link{dm_forecast}}, \code{\link{make_posterior}}
#'
#' @examples
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#'
#' dm <- dynhr_model(system.file("extdata/models/nk_demo.mod",
#'                               package = "dynhr"),
#'                   data = Y, obs_vars = obs_vars)
#' dm
#'
#' ## One object carries every piece the next step needs.
#' dm <- dm_solve(dm)
#' dm <- dm_posterior(dm)
#' theta0 <- setNames(dm$prior_spec$mean, dm$prior_spec$name)
#' dm$log_post(theta0)$logpost
#' @export
dynhr_model <- function(model, data = NULL, obs_vars = NULL, compiled = NULL,
                        prior_spec = NULL, me_variance = 0,
                        likelihood = "gaussian", params = NULL,
                        max_order = 1L, verbose = FALSE, ...) {
  ## ---- model: a parsed object or a path ---------------------------------
  if (is.character(model) && length(model) == 1L) {
    if (!file.exists(model))
      stop("dynhr_model: `model` looks like a file path but does not exist: ",
           model, call. = FALSE)
    model <- parse_mod(model, verbose = verbose)
  }
  if (!inherits(model, "dynhr_mod"))
    stop("dynhr_model: `model` must be a dynhr_mod (from parse_mod()) or a ",
         "path to a .mod file.", call. = FALSE)

  ## ---- compiled derivatives (EAGER) -------------------------------------
  if (is.null(compiled)) {
    compiled <- compile_model(model, verbose = verbose,
                              max_order = as.integer(max_order))
  } else if (!inherits(compiled, "dynhr_compiled")) {
    stop("dynhr_model: `compiled` must be a dynhr_compiled object from ",
         "compile_model().", call. = FALSE)
  }

  ## ---- observed variables ------------------------------------------------
  if (is.null(obs_vars) || length(obs_vars) == 0L)
    obs_vars <- model$obs_vars %||% model$varobs_names

  ## ---- data --------------------------------------------------------------
  if (!is.null(data)) {
    if (is.data.frame(data)) data <- as.matrix(data)
    if (!is.matrix(data))
      stop("dynhr_model: `data` must be a T x n_obs matrix or data.frame ",
           "(rows = periods, columns = observables).", call. = FALSE)
    if (!is.null(obs_vars) && !is.null(colnames(data))) {
      missing_cols <- setdiff(obs_vars, colnames(data))
      if (length(missing_cols))
        stop("dynhr_model: `data` has no column for observable(s): ",
             paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    if (!is.null(obs_vars) && is.null(colnames(data)) &&
        ncol(data) != length(obs_vars))
      stop("dynhr_model: `data` has ", ncol(data), " unnamed columns but ",
           length(obs_vars), " obs_vars were supplied.", call. = FALSE)
  }

  ## ---- prior spec --------------------------------------------------------
  if (is.null(prior_spec)) {
    prior_spec <- tryCatch(extract_prior_spec(model, verbose = verbose),
                           error = function(e) NULL)
  }

  ## ---- parameter values --------------------------------------------------
  if (is.null(params)) params <- model$param_values

  obj <- list(
    model       = model,
    compiled    = compiled,
    data        = data,
    obs_vars    = obs_vars,
    prior_spec  = prior_spec,
    me_variance = me_variance,
    likelihood  = likelihood,
    lik_args    = list(...),
    params      = params,
    ## result slots, filled by the verbs
    steady      = NULL,
    solution    = NULL,
    order       = NULL,
    log_post    = NULL,
    mode        = NULL,
    chains      = NULL,
    diagnostics = NULL,
    irf         = NULL,
    forecast    = NULL
  )
  class(obj) <- "dynhr_model"
  obj
}


# ============================================================================
# Verbs
# ============================================================================

#' Solve the model's perturbation decision rules
#'
#' Computes the steady state (cached on the object as `$steady`) and calls
#' \code{\link{solve_perturbation}} at `order`, storing the resulting
#' `DecisionRules*` object in `$solution` and the order in `$order`.
#'
#' @param dm    A \code{\link{dynhr_model}}
#' @param order Perturbation order (default `1L`).  The object's compiled
#'   derivatives must be deep enough -- see `max_order` in
#'   \code{\link{dynhr_model}}.
#' @param ...   Further arguments for \code{\link{solve_perturbation}}
#'   (e.g. `Sigma_e`, `loglinear`).
#' @return The updated `dynhr_model`, with `$steady`, `$solution` and
#'   `$order` set.
#' @seealso \code{\link{dynhr_model}}, \code{\link{solve_perturbation}}
#' @examples
#' dm <- dynhr_model(system.file("extdata/models/rbc.mod", package = "dynhr"))
#' dm <- dm_solve(dm)
#' class(dm$solution)
#' @export
dm_solve <- function(dm, order = 1L, ...) {
  .dm_check(dm, "dm_solve")
  if (is.null(dm$steady)) {
    ss <- solve_steady(dm$compiled, dm$params,
                       endo_names = dm$model$var_names,
                       exo_names  = dm$model$varexo_names,
                       verbose    = FALSE)
    dm <- .dm_set(dm, "steady", ss)
  }
  dr <- solve_perturbation(dm$model, dm$compiled, dm$steady$values, dm$params,
                           order = as.integer(order), ...)
  dm <- .dm_set(dm, "solution", dr)
  .dm_set(dm, "order", as.integer(order))
}


#' Build the log-posterior closure for a dynhr_model
#'
#' Re-threads the object's `model`, `data`, `prior_spec`, `obs_vars`,
#' `compiled`, `me_variance` and `likelihood` into
#' \code{\link{make_posterior}} and stores the returned closure in
#' `$log_post`.  Arguments in `...` are merged over the likelihood arguments
#' the constructor captured.
#'
#' @param dm  A \code{\link{dynhr_model}} carrying `data` and `prior_spec`
#' @param ... Further arguments for \code{\link{make_posterior}} /
#'   `make_log_posterior()` (e.g. `likelihood`, `order`, `lik_init`).  These
#'   override the values stored on the object for this call only.
#' @return The updated `dynhr_model`, with `$log_post` set to the closure
#'   `function(theta) list(logpost, loglik, logprior)`.
#' @seealso \code{\link{dynhr_model}}, \code{\link{make_posterior}}
#' @examples
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' dm <- dynhr_model(system.file("extdata/models/nk_demo.mod", package = "dynhr"),
#'                   data = Y, obs_vars = obs_vars)
#' dm <- dm_posterior(dm)
#' theta0 <- setNames(dm$prior_spec$mean, dm$prior_spec$name)
#' dm$log_post(theta0)$logpost
#' @export
dm_posterior <- function(dm, ...) {
  .dm_check(dm, "dm_posterior")
  if (is.null(dm$data))       .dm_require_field("dm_posterior", "data")
  if (is.null(dm$prior_spec)) .dm_require_field("dm_posterior", "prior_spec")

  extra <- utils::modifyList(dm$lik_args, list(...))
  base  <- list(model       = dm$model,
                data        = dm$data,
                prior_spec  = dm$prior_spec,
                obs_vars    = dm$obs_vars,
                compiled    = dm$compiled,
                me_variance = dm$me_variance,
                likelihood  = dm$likelihood)
  ## anything in `...` (or in the constructor's stored likelihood args)
  ## overrides the stored field of the same name for this call
  .dm_set(dm, "log_post",
          do.call(make_log_posterior, utils::modifyList(base, extra)))
}


#' Find the posterior mode of a dynhr_model
#'
#' Calls \code{\link{find_mode}} on the object's `$log_post` closure and
#' stores the result in `$mode`.
#'
#' @param dm         A \code{\link{dynhr_model}} whose `$log_post` is set
#'   (see \code{\link{dm_posterior}})
#' @param theta_init Named numeric starting vector.  Defaults to the prior
#'   means from `$prior_spec`.
#' @param ...        Further arguments for \code{\link{find_mode}}
#'   (`n_iter`, `method`, `verbose`).
#' @return The updated `dynhr_model`, with `$mode` set to the list returned
#'   by \code{\link{find_mode}} (`theta_mode`, `logpost`, `convergence`, ...).
#' @seealso \code{\link{dynhr_model}}, \code{\link{find_mode}}
#' @examples
#' \donttest{
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' dm <- dynhr_model(system.file("extdata/models/nk_demo.mod", package = "dynhr"),
#'                   data = Y, obs_vars = obs_vars)
#' dm <- dm_mode(dm_posterior(dm), n_iter = 50L, verbose = FALSE)
#' dm$mode$logpost
#' }
#' @export
dm_mode <- function(dm, theta_init = NULL, ...) {
  .dm_check(dm, "dm_mode")
  if (is.null(dm$log_post)) .dm_require("dm_mode", "log_post", "dm_posterior")
  if (is.null(dm$prior_spec)) .dm_require_field("dm_mode", "prior_spec")
  if (is.null(theta_init))
    theta_init <- stats::setNames(dm$prior_spec$mean, dm$prior_spec$name)
  .dm_set(dm, "mode", find_mode(dm$log_post, theta_init, dm$prior_spec, ...))
}


#' Sample the posterior of a dynhr_model
#'
#' Dispatches to \code{\link{mcmc}} (RWMH), \code{\link{nuts}},
#' \code{\link{smc}} or \code{\link{dime}} on the object's `$log_post`
#' closure and stores the returned \code{\link{dynhr_chains}} in `$chains`.
#'
#' @section Starting point:
#' For the two point-started samplers (`"rwmh"`, `"nuts"`), `theta0` defaults
#' to `$mode$theta_mode` when \code{\link{dm_mode}} has been run and to the
#' prior means otherwise -- starting from the prior is a documented default,
#' not an error.  `Sigma_prop` (RWMH only) defaults to
#' `diag(prior_spec$std^2)`; pass the scaled inverse Hessian from the mode
#' finder explicitly for a production run.  `"smc"` and `"dime"` initialise
#' from the prior and use neither argument.
#'
#' @param dm         A \code{\link{dynhr_model}} whose `$log_post` is set
#' @param sampler    One of `"rwmh"` (default), `"nuts"`, `"smc"`, `"dime"`
#' @param theta0     Named numeric starting vector for `"rwmh"` / `"nuts"`
#'   (see Starting point)
#' @param Sigma_prop Proposal covariance for `"rwmh"` (see Starting point)
#' @param ...        Further arguments for the chosen sampler (e.g.
#'   `n_draws`, `n_warmup`, `n_particles`, `verbose`).
#' @return The updated `dynhr_model`, with `$chains` set to a
#'   \code{\link{dynhr_chains}} object.
#' @seealso \code{\link{dynhr_model}}, \code{\link{mcmc}},
#'   \code{\link{nuts}}, \code{\link{smc}}, \code{\link{dime}}
#' @examples
#' \donttest{
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' dm <- dm_posterior(dynhr_model(
#'   system.file("extdata/models/nk_demo.mod", package = "dynhr"),
#'   data = Y, obs_vars = obs_vars))
#' set.seed(1)
#' dm <- dm_sample(dm, n_draws = 50L, n_warmup = 50L, verbose = FALSE)
#' dm$chains
#' }
#' @export
dm_sample <- function(dm, sampler = c("rwmh", "nuts", "smc", "dime"),
                      theta0 = NULL, Sigma_prop = NULL, ...) {
  .dm_check(dm, "dm_sample")
  sampler <- match.arg(sampler)
  if (is.null(dm$log_post)) .dm_require("dm_sample", "log_post", "dm_posterior")
  if (is.null(dm$prior_spec)) .dm_require_field("dm_sample", "prior_spec")

  if (sampler %in% c("smc", "dime")) {
    res <- if (sampler == "smc") smc(dm$log_post, dm$prior_spec, ...)
           else                  dime(dm$log_post, dm$prior_spec, ...)
    return(.dm_set(dm, "chains", res))
  }

  if (is.null(theta0)) {
    theta0 <- if (!is.null(dm$mode)) dm$mode$theta_mode
              else stats::setNames(dm$prior_spec$mean, dm$prior_spec$name)
  }
  if (sampler == "nuts")
    return(.dm_set(dm, "chains", nuts(dm$log_post, theta0, ...)))

  if (is.null(Sigma_prop))
    Sigma_prop <- diag(dm$prior_spec$std^2, nrow = length(theta0))
  .dm_set(dm, "chains", mcmc(dm$log_post, theta0, Sigma_prop, ...))
}


#' Run the dynhr diagnostic battery on a dynhr_model
#'
#' Forwards whichever of the object's pieces are present -- model, compiled
#' derivatives, decision rules, parameters, data, observed names, prior spec
#' and posterior draws -- to \code{\link{run_diagnostics}}, and stores the
#' returned list in `$diagnostics`.  Stages whose inputs are absent are
#' skipped by the battery itself, so this works with a partially-run
#' pipeline.
#'
#' @param dm  A \code{\link{dynhr_model}}, ideally after
#'   \code{\link{dm_solve}} and \code{\link{dm_sample}}
#' @param ... Further arguments for \code{\link{run_diagnostics}}
#'   (e.g. `verbose`, `report`).  These override the forwarded fields.
#' @return The updated `dynhr_model`, with `$diagnostics` set.
#' @seealso \code{\link{dynhr_model}}, \code{\link{run_diagnostics}}
#' @examples
#' \donttest{
#' dm <- dm_solve(dynhr_model(system.file("extdata/models/rbc.mod",
#'                                        package = "dynhr")))
#' dm <- suppressWarnings(dm_diagnostics(dm, verbose = FALSE))
#' names(dm$diagnostics)
#' }
#' @export
dm_diagnostics <- function(dm, ...) {
  .dm_check(dm, "dm_diagnostics")
  base <- list(model     = dm$model,
               compiled  = dm$compiled,
               dr        = dm$solution,
               ss        = dm$steady,
               params    = dm$params,
               priors    = dm$prior_spec,
               data      = dm$data,
               obs_names = dm$obs_vars,
               draws     = if (!is.null(dm$chains)) dm$chains$chain else NULL,
               theta_mode = if (!is.null(dm$mode)) dm$mode$theta_mode else NULL)
  base <- base[!vapply(base, is.null, logical(1))]
  args <- utils::modifyList(base, list(...))
  .dm_set(dm, "diagnostics", do.call(run_diagnostics, args))
}


#' Impulse responses from a dynhr_model
#'
#' Calls the impulse-response routine matching the order the object was
#' solved at (\code{\link{compute_irfs}} at order 1,
#' \code{\link{compute_irfs_order2}} at order 2,
#' \code{\link{compute_irfs_order3}} at order 3) and stores the result in
#' `$irf`.
#'
#' @param dm        A \code{\link{dynhr_model}} whose `$solution` is set
#'   (see \code{\link{dm_solve}})
#' @param n_periods IRF horizon (default `40L`)
#' @param ...       Further arguments for the underlying IRF routine
#'   (e.g. `shock_size`, and `pruning` at order 2).
#' @return The updated `dynhr_model`, with `$irf` set.
#' @seealso \code{\link{dynhr_model}}, \code{\link{compute_irfs}}
#' @examples
#' dm <- dm_irf(dm_solve(dynhr_model(
#'   system.file("extdata/models/rbc.mod", package = "dynhr"))),
#'   n_periods = 8L)
#' names(dm$irf)
#' @export
dm_irf <- function(dm, n_periods = 40L, ...) {
  .dm_check(dm, "dm_irf")
  if (is.null(dm$solution)) .dm_require("dm_irf", "solution", "dm_solve")
  ord <- dm$order %||% 1L
  fn <- switch(as.character(ord),
               "1" = compute_irfs,
               "2" = compute_irfs_order2,
               "3" = compute_irfs_order3,
               stop("dm_irf(): no IRF routine for perturbation order ", ord,
                    "; call compute_irfs_order4()/compute_irfs_order5() on ",
                    "`dm$solution` directly.", call. = FALSE))
  .dm_set(dm, "irf",
          fn(dm$solution, dm$model, n_periods = as.integer(n_periods),
             params = dm$params, ...))
}


#' Conditional forecast from a dynhr_model
#'
#' Re-threads the object's `model`, `solution`, `data`, `obs_vars` and
#' `compiled` into \code{\link{conditional_forecast}} and stores the result
#' in `$forecast`.
#'
#' @param dm         A \code{\link{dynhr_model}} carrying `data` and a
#'   `$solution` (see \code{\link{dm_solve}})
#' @param conditions Conditioning information, as in
#'   \code{\link{conditional_forecast}} (`NULL` gives an unconditional
#'   forecast)
#' @param horizon    Forecast horizon (default `8L`)
#' @param ...        Further arguments for \code{\link{conditional_forecast}}
#'   (e.g. `type`, `method`, `n_draws`).
#' @return The updated `dynhr_model`, with `$forecast` set.
#' @seealso \code{\link{dynhr_model}}, \code{\link{conditional_forecast}}
#' @examples
#' \donttest{
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' dm <- dm_solve(dynhr_model(
#'   system.file("extdata/models/nk_demo.mod", package = "dynhr"),
#'   data = Y, obs_vars = obs_vars))
#' dm <- dm_forecast(dm, horizon = 4L)
#' dim(dm$forecast$forecast)
#' }
#' @export
dm_forecast <- function(dm, conditions = NULL, horizon = 8L, ...) {
  .dm_check(dm, "dm_forecast")
  if (is.null(dm$solution)) .dm_require("dm_forecast", "solution", "dm_solve")
  if (is.null(dm$data))     .dm_require_field("dm_forecast", "data")
  .dm_set(dm, "forecast",
          conditional_forecast(dm$model, dm$solution, dm$data,
                               conditions = conditions,
                               horizon    = as.integer(horizon),
                               obs_vars   = dm$obs_vars,
                               compiled   = dm$compiled, ...))
}


# ============================================================================
# S3 methods
# ============================================================================

#' Print a dynhr_model object
#'
#' Compact one-screen header: model size, data dimensions, observables, and
#' which pipeline steps have been run.
#'
#' @param x   A \code{\link{dynhr_model}}
#' @param ... Ignored
#' @return `x`, invisibly
#' @export
print.dynhr_model <- function(x, ...) {
  nm <- x$model$source_file %||% NA_character_
  nm <- if (length(nm) != 1L || is.na(nm)) "<in-memory model>" else basename(nm)
  cat(sprintf("\n<dynhr_model: %s>\n", nm))
  cat(sprintf("  Variables  : %d endogenous, %d exogenous, %d equations\n",
              length(x$model$var_names), length(x$model$varexo_names),
              length(x$model$equations)))
  cat(sprintf("  Compiled   : max_order %s\n",
              format(x$compiled$dynamic$max_order %||% NA)))
  if (is.null(x$data)) {
    cat("  Data       : <none>\n")
  } else {
    cat(sprintf("  Data       : %d periods x %d observables\n",
                nrow(x$data), ncol(x$data)))
  }
  cat(sprintf("  Observables: %s\n",
              if (length(x$obs_vars)) paste(x$obs_vars, collapse = ", ")
              else "<none>"))
  cat(sprintf("  Priors     : %s\n",
              if (is.null(x$prior_spec)) "<none>"
              else sprintf("%d estimated parameters", nrow(x$prior_spec))))
  cat(sprintf("  Likelihood : %s (me_variance = %s)\n",
              x$likelihood, format(x$me_variance)))

  done <- names(.dm_steps)[vapply(names(.dm_steps),
                                  function(s) !is.null(x[[s]]), logical(1))]
  cat(sprintf("  Steps run  : %s\n",
              if (length(done)) paste(unname(.dm_steps[done]), collapse = " -> ")
              else "<none>"))
  todo <- setdiff(names(.dm_steps), done)
  if (length(todo))
    cat(sprintf("  Next       : %s\n", unname(.dm_steps[todo[1]])))
  invisible(x)
}


#' Summarise a dynhr_model object
#'
#' Prints the \code{\link{print.dynhr_model}} header, then a per-step detail
#' block for every step that has been run (perturbation order and state
#' count, the log-posterior at the prior mean, the mode's log-posterior,
#' the chain dimensions, the diagnostic names).
#'
#' @param object A \code{\link{dynhr_model}}
#' @param ...    Ignored
#' @return A named list of the per-step details, invisibly
#' @export
summary.dynhr_model <- function(object, ...) {
  print(object)
  out <- list()
  cat("\n--- Step detail ---\n")
  if (!is.null(object$solution)) {
    out$order  <- object$order
    out$n_state <- length(object$solution$state_idx)
    cat(sprintf("  dm_solve       : order %d, %d states, class %s\n",
                out$order, out$n_state, class(object$solution)[1]))
  }
  if (!is.null(object$log_post) && !is.null(object$prior_spec)) {
    th <- stats::setNames(object$prior_spec$mean, object$prior_spec$name)
    lp <- tryCatch(object$log_post(th)$logpost, error = function(e) NA_real_)
    out$logpost_prior_mean <- lp
    cat(sprintf("  dm_posterior   : logpost at prior mean = %.6f\n", lp))
  }
  if (!is.null(object$mode)) {
    out$mode_logpost <- object$mode$logpost
    cat(sprintf("  dm_mode        : logpost = %.6f (%s, convergence %s)\n",
                object$mode$logpost, object$mode$method %||% "?",
                format(object$mode$convergence %||% NA)))
  }
  if (!is.null(object$chains)) {
    out$n_draws <- nrow(object$chains$chain)
    out$sampler <- object$chains$sampler
    cat(sprintf("  dm_sample      : %s, %d draws x %d parameters\n",
                toupper(out$sampler %||% "?"), out$n_draws,
                ncol(object$chains$chain)))
  }
  if (!is.null(object$diagnostics)) {
    out$diagnostics <- names(object$diagnostics)
    cat(sprintf("  dm_diagnostics : %d stages (%s)\n",
                length(out$diagnostics),
                paste(utils::head(out$diagnostics, 6L), collapse = ", ")))
  }
  if (!is.null(object$irf)) {
    out$irf <- names(object$irf)
    cat(sprintf("  dm_irf         : %d shock(s)\n", length(out$irf)))
  }
  if (!is.null(object$forecast)) {
    cat("  dm_forecast    : present\n")
    out$forecast <- TRUE
  }
  if (!length(out)) cat("  <no step has been run>\n")
  invisible(out)
}
