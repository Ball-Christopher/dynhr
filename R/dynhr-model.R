## R/dynhr-model.R
## --------------------------------------------------------------------------
## `dynhr_model`: a pipeline OBJECT holding an estimation spec (W82, brief 28
## option (b)).
##
## The object is list(spec, params, caches). `spec` is a
## `dynhr_estimation_spec` -- the SINGLE source of truth for model, compiled
## model, data, obs_vars, priors and the likelihood (me_variance, lik_init,
## heteroskedastic_shocks / filter_tunes / stochastic_volatility / OBC
## routing), and for the mode, sampler and compute settings. The `$` / `[[`
## accessors read `dm$model`, `dm$data`, `dm$obs_vars`, `dm$prior_spec`, ...
## from it; they cannot be assigned (use update(), which clears the results
## that depend on what changed).
##
## The estimation verbs are thin: dm_posterior() builds its closure with the
## spec runner's own builder (.est_build_posterior()), dm_mode() and
## dm_sample() run their stage through run_estimation() on the held spec, so
## a dm_*() result is the run_estimation() result of the same spec (pinned by
## tests/testthat/test-dynhr-model-golden.R). dm_solve() / dm_irf() /
## dm_forecast() / dm_diagnostics() call the functional solvers, at the
## calibrated `params` (default), the posterior mode or the posterior mean.
##
## Design notes
##   * S3, not R6: R6 is not an Import of dynhr and must not become one.
##   * Compilation is EAGER (in the constructor and in update(model =)): an
##     object that has been built has already paid for `compile_model()`.
##   * A verb whose own arguments edit the spec (dm_posterior(me_variance =),
##     dm_mode(n_iter =), dm_sample(n_draws =)) records the edit in the
##     returned object's spec, so the spec always describes what was run.
## --------------------------------------------------------------------------


## Slots a verb writes, in pipeline order, for print()/summary().
.dm_steps <- c(solution = "dm_solve", log_post = "dm_posterior",
               mode = "dm_mode", chains = "dm_sample",
               diagnostics = "dm_diagnostics", irf = "dm_irf",
               forecast = "dm_forecast")

## Every result (cache) slot. `*_at` record the point ("params", "mode",
## "posterior_mean") a solution / IRF / forecast was computed at;
## `solution_params` the parameter vector the solution was solved at; `fit`
## the run_estimation() result of the last dm_sample().
.dm_cache_slots <- c("steady", "solution", "order", "solution_at",
                     "solution_params", "log_post", "mode", "chains", "fit",
                     "diagnostics", "irf", "irf_at", "forecast", "forecast_at")
.dm_solution_slots <- c("steady", "solution", "order", "solution_at",
                        "solution_params")

## Fields `$` reads from the held spec (read-only).
.dm_spec_fields <- c("model", "compiled", "data", "obs_vars", "prior_spec",
                     "priors", "me_variance", "likelihood", "lik_args")

## The points a solution can be computed at.
.dm_at_choices <- c("params", "mode", "posterior_mean")


#' Fail with a message naming the missing pipeline step
#'
#' @param verb   Name of the verb that was called
#' @param slot   Name of the empty slot
#' @param needed Name of the verb that fills it
#' @param at     Optional \code{at} value the verb was called with
#' @return Never returns; signals an error
#' @noRd
.dm_require <- function(verb, slot, needed, at = NULL) {
  .dynhr_abort(sprintf("%s(%s): `$%s` is empty -- run %s() on this object first.",
                       verb, if (is.null(at)) "" else sprintf("at = \"%s\"", at),
                       slot, needed),
               class = "dynhr_error_dm_missing_step")
}


#' Fail when a field the constructor should have stored is absent
#'
#' @param verb  Name of the verb that was called
#' @param field Name of the missing constructor field
#' @return Never returns; signals an error
#' @noRd
.dm_require_field <- function(verb, field) {
  .dynhr_abort(sprintf("%s(): `$%s` is empty -- supply %s= to dynhr_model() (or update()).",
                       verb, field, field),
               class = "dynhr_error_dm_missing_field")
}


#' Check that an object is a dynhr_model
#'
#' @param dm   Candidate object
#' @param verb Calling verb, for the error message
#' @return \code{dm}, invisibly
#' @noRd
.dm_check <- function(dm, verb) {
  if (!inherits(dm, "dynhr_model"))
    .dynhr_abort(sprintf("%s(): `dm` must be a dynhr_model object (see dynhr_model()).",
                         verb), class = "dynhr_error_bad_argument")
  invisible(dm)
}


#' Store values in slots of a dynhr_model and return the new object
#'
#' @param dm     A \code{dynhr_model}
#' @param slot   Slot name to write
#' @param value  Value to store (\code{NULL} empties the slot)
#' @return The updated \code{dynhr_model}
#' @noRd
.dm_set <- function(dm, slot, value) {
  cls <- oldClass(dm)
  dm <- unclass(dm)
  dm[slot] <- list(value)
  class(dm) <- cls
  dm
}

## Several slots at once (a named list).
.dm_set_many <- function(dm, values) {
  for (nm in names(values)) dm <- .dm_set(dm, nm, values[[nm]])
  dm
}

## Empty the named slots.
.dm_drop <- function(dm, slots) {
  for (s in unique(slots)) dm <- .dm_set(dm, s, NULL)
  dm
}


#' Read a field of a dynhr_model (spec-derived or a slot)
#'
#' @param x    A \code{dynhr_model}
#' @param name Field name
#' @return The field value (\code{NULL} when absent)
#' @noRd
.dm_get <- function(x, name) {
  spec <- .subset2(x, "spec")
  switch(name,
    model       = spec$model$mod,
    compiled    = spec$model$compiled %||% spec$model$solved$compiled,
    data        = spec$data$value %||%
                    (if (!is.null(spec$data$path)) .est_load_data(spec)),
    obs_vars    = if (length(spec$obs_vars)) spec$obs_vars,
    prior_spec  = ,
    priors      = .subset2(x, ".prior_spec"),
    me_variance = spec$likelihood$me_variance,
    likelihood  = spec$likelihood$type,
    lik_args    = spec$likelihood$extra,
    .subset2(x, name))
}

## Refuse an assignment to a spec-derived field or to the spec / params.
.dm_assign <- function(x, name, value) {
  if (name %in% c(.dm_spec_fields, "spec", "params", ".prior_spec"))
    .dynhr_abort("dynhr_model: `$", name, "` is part of the object's ",
                 "estimation spec; change it with update(dm, ",
                 if (identical(name, "priors")) "prior_spec" else name,
                 " = ...), which also clears the results that depend on it.",
                 class = "dynhr_error_bad_argument")
  .dm_set(x, name, value)
}

## The prior spec in force: the spec's own, else the model's estimated_params
## (NULL for a model without one).
.dm_resolve_priors <- function(spec, verbose = FALSE) {
  spec$model$prior_spec %||%
    tryCatch(extract_prior_spec(spec$model$mod, verbose = verbose),
             error = function(e) .dynhr_reraise_bug(e, NULL))
}

## A new object around a held spec: every result slot empty.
.dm_new <- function(spec, params, verbose = FALSE) {
  obj <- list(spec = spec, params = params,
              .prior_spec = .dm_resolve_priors(spec, verbose))
  for (s in .dm_cache_slots) obj[s] <- list(NULL)
  class(obj) <- "dynhr_model"
  obj
}

## Does the spec carry data?
.dm_has_data <- function(spec) !is.null(spec$data$value) || !is.null(spec$data$path)

## The inputs every estimation verb needs, as named errors.
.dm_need_inputs <- function(dm, verb) {
  if (!.dm_has_data(dm$spec)) .dm_require_field(verb, "data")
  if (is.null(dm$prior_spec)) .dm_require_field(verb, "prior_spec")
  invisible(dm)
}


# ============================================================================
# Spec edits and cache invalidation
# ============================================================================

#' Replace the held spec and clear the results that depend on what changed
#'
#' Parts are compared by the spec's content hashes: a model, data, obs_vars,
#' likelihood, options or mode$extra change invalidates the log-posterior and
#' everything estimated from it; a mode change the mode and the draws; a
#' sampler or compute change the draws; a model or params change the
#' solution, IRFs and forecasts. A solution / IRF / forecast computed at the
#' mode or the posterior mean goes with the estimate it was computed at.
#'
#' @param dm     A \code{dynhr_model}
#' @param spec   The new spec
#' @param params New parameter vector, or \code{NULL} to keep
#' @return The updated \code{dynhr_model}
#' @noRd
.dm_respec <- function(dm, spec, params = NULL) {
  old <- dm$spec
  ch <- function(p) !identical(old$hashes[[p]], spec$hashes[[p]])
  model_ch  <- ch("model") || !identical(old$model$compiled, spec$model$compiled)
  data_ch   <- ch("data") || !identical(old$obs_vars, spec$obs_vars)
  params_ch <- !is.null(params) && !identical(params, dm$params)
  post_ch   <- model_ch || data_ch || ch("likelihood") || ch("options") ||
    !identical(old$mode$extra, spec$mode$extra)
  mode_ch   <- post_ch || ch("mode")
  samp_ch   <- mode_ch || ch("sampler") || ch("compute")
  any_ch    <- samp_ch || params_ch || ch("outputs")
  drop <- c(if (post_ch) "log_post", if (mode_ch) "mode",
            if (samp_ch) c("chains", "fit"), if (any_ch) "diagnostics")
  stale <- c(if (mode_ch) "mode", if (samp_ch) "posterior_mean")
  if (model_ch || params_ch || (dm$solution_at %||% "params") %in% stale)
    drop <- c(drop, .dm_solution_slots)
  if (model_ch || params_ch || (dm$irf_at %||% "params") %in% stale)
    drop <- c(drop, "irf", "irf_at")
  if (model_ch || params_ch || data_ch ||
      (dm$forecast_at %||% "params") %in% stale)
    drop <- c(drop, "forecast", "forecast_at")
  dm <- .dm_set(dm, "spec", spec)
  if (!is.null(params)) dm <- .dm_set(dm, "params", params)
  if (ch("model")) dm <- .dm_set(dm, ".prior_spec", .dm_resolve_priors(spec))
  .dm_drop(dm, drop)
}

## After a new mode ("mode") or new draws ("sample"): clear what was computed
## from the previous one.
.dm_drop_after <- function(dm, stage) {
  ats  <- if (identical(stage, "mode")) c("mode", "posterior_mean") else "posterior_mean"
  drop <- c(if (identical(stage, "mode")) c("chains", "fit"), "diagnostics")
  if ((dm$solution_at %||% "params") %in% ats) drop <- c(drop, .dm_solution_slots)
  if ((dm$irf_at %||% "params") %in% ats) drop <- c(drop, "irf", "irf_at")
  if ((dm$forecast_at %||% "params") %in% ats) drop <- c(drop, "forecast", "forecast_at")
  .dm_drop(dm, drop)
}

#' Edit the held spec (the engine of update() and of the verbs' own arguments)
#'
#' @param dm A \code{dynhr_model}
#' @param model,data,obs_vars,likelihood,mode,sampler,compute,outputs,options
#'   As in \code{update.dynhr_estimation_spec()}; \code{model} is compiled
#'   eagerly, \code{likelihood} may also be a type name.
#' @param params,prior_spec New parameter vector / prior spec (\code{NULL}
#'   keeps; a new model resets params to its calibration)
#' @return The updated \code{dynhr_model}
#' @noRd
.dm_edit <- function(dm, model = NULL, data = NULL, obs_vars = NULL,
                     likelihood = NULL, mode = NULL, sampler = NULL,
                     compute = NULL, outputs = NULL, options = NULL,
                     params = NULL, prior_spec = NULL, verbose = FALSE) {
  old <- dm$spec
  nz <- function(x) if (is.list(x) && !is.object(x) && !length(x)) NULL else x
  if (is.character(likelihood) && length(likelihood) == 1L)
    likelihood <- list(type = likelihood)
  mp <- NULL
  if (!is.null(model) || !is.null(prior_spec)) {
    if (is.character(model) && length(model) == 1L) {
      if (!file.exists(model))
        .dynhr_abort("update.dynhr_model: `model` looks like a file path but ",
                     "does not exist: ", model, class = "dynhr_error_bad_argument")
      model <- parse_mod(model, verbose = verbose)
    }
    if (!is.null(model) && !inherits(model, "dynhr_mod"))
      .dynhr_abort("update.dynhr_model: `model` must be a dynhr_mod (from ",
                   "parse_mod()) or a path to a .mod file.",
                   class = "dynhr_error_bad_argument")
    mod <- model %||% old$model$mod
    cm  <- if (is.null(model)) old$model$compiled %||% old$model$solved$compiled
           else compile_model(model, verbose = verbose,
                              max_order = old$model$max_order)
    ps  <- prior_spec %||% (if (is.null(model)) old$model$prior_spec)
    mp  <- .spec_model_part(mod, compiled = cm, max_order = old$model$max_order,
                            prior_spec = ps)
    if (!is.null(model) && is.null(obs_vars))
      obs_vars <- if (length(old$obs_vars)) old$obs_vars
                  else model$obs_vars %||% model$varobs_names %||% character(0)
  }
  parts <- .spec_update_parts(old, data = data, obs_vars = obs_vars,
                              likelihood = nz(likelihood), mode = nz(mode),
                              sampler = nz(sampler), compute = nz(compute),
                              outputs = nz(outputs), options = options,
                              lenient_obs = TRUE, model_part = mp)
  spec <- .spec_assemble(parts, validate = .dm_has_data(parts))
  if (is.null(params) && !is.null(model)) params <- model$param_values
  .dm_respec(dm, spec, params = params)
}

## Split a verb's `...` into the spec components that take each name.
.dm_split_args <- function(dots, comps, verb, method = NULL) {
  out <- stats::setNames(lapply(comps, function(cp) list()), comps)
  if (!length(dots)) return(out)
  nms <- names(dots)
  if (is.null(nms) || any(!nzchar(nms)))
    .dynhr_abort(verb, "(): every argument in `...` must be named.",
                 class = "dynhr_error_bad_argument")
  flds <- lapply(comps, function(cp) switch(cp,
    sampler = .spec_sampler_table[[method]],
    mode    = setdiff(names(.spec_schema$mode), "result"),
    names(.spec_schema[[cp]])))
  for (nm in nms) {
    hit <- which(vapply(flds, function(f) nm %in% f, logical(1)))
    if (!length(hit))
      .dynhr_abort(verb, "(): `", nm, "` is not a field of ",
                   paste0(ifelse(comps == "sampler",
                                 sprintf("sampler \"%s\"", method %||% ""),
                                 paste0(comps, "_spec()")), collapse = " or "),
                   ". Valid: ", paste(unlist(flds), collapse = ", "),
                   " (sampler-specific arguments go in extra = list(...)).",
                   class = "dynhr_error_spec_unknown_field")
    out[[comps[hit[1L]]]][nm] <- list(dots[[nm]])
  }
  out
}

## dm_posterior()'s `...`: likelihood fields by name or by the constructor
## argument that sets them (power -> power_posterior, likelihood -> type, ...),
## everything else merged into likelihood$extra.
.dm_likelihood_edits <- function(lik, dots) {
  flds <- setdiff(names(.spec_schema$likelihood), "extra")
  out <- list()
  extra <- lik$extra
  for (k in names(dots)) {
    if (identical(k, "extra")) {
      extra <- utils::modifyList(extra, dots[[k]])
      next
    }
    tgt <- if (k %in% flds) paste0("likelihood$", k)
           else .spec_extra_shadow(k, "likelihood", lik$type)
    if (!is.null(tgt) && startsWith(tgt, "likelihood$"))
      out[sub("^likelihood\\$", "", tgt)] <- list(dots[[k]])
    else
      extra[k] <- list(dots[[k]])
  }
  if (!identical(extra, lik$extra)) out$extra <- extra
  out
}


# ============================================================================
# Constructor
# ============================================================================

#' Bundle a DSGE model, its data and its estimation settings into one object
#'
#' \code{dynhr_model()} builds a pipeline object around a
#' \code{\link{dynhr_estimation_spec}}: the parsed and compiled model, the
#' data, the observables, the priors, the likelihood and the mode, sampler and
#' compute settings. The \code{dm_*()} verbs (\code{\link{dm_solve}},
#' \code{\link{dm_posterior}}, \code{\link{dm_mode}}, \code{\link{dm_sample}},
#' \code{\link{dm_diagnostics}}, \code{\link{dm_irf}},
#' \code{\link{dm_forecast}}) each take the object as their first argument
#' and return a NEW object with the step's result stored in a slot, so a
#' whole pipeline reads as one chain.
#'
#' @section One source of truth:
#' The object is a list of the held spec (\code{dm$spec}), the parameter
#' vector \code{dm$params} and the result slots. The fields \code{model},
#' \code{compiled}, \code{data}, \code{obs_vars}, \code{prior_spec} (alias
#' \code{priors}), \code{me_variance}, \code{likelihood} (the type) and
#' \code{lik_args} (\code{likelihood$extra}) are read from the spec by
#' \code{$} and \code{[[}; they cannot be assigned. Change them with
#' \code{\link{update.dynhr_model}}, which edits the spec and clears every
#' stored result that depends on what changed.
#'
#' The estimation verbs run through the spec runner: \code{\link{dm_posterior}}
#' builds its closure with the builder \code{\link{run_estimation}} uses (so a
#' \code{.mod} \code{heteroskedastic_shocks}, \code{filter_tunes},
#' \code{stochastic_volatility} block or an OBC constraint is applied), and
#' \code{\link{dm_mode}} / \code{\link{dm_sample}} run their stage with
#' \code{run_estimation()} on the held spec, so their results equal
#' \code{run_estimation(dynhr_estimation_spec(dm))}.
#' \code{dynhr_estimation_spec(dm)} and \code{as_estimation_spec(dm)} return
#' the held spec.
#'
#' @section Compilation:
#' Compilation is \emph{eager}: unless \code{compiled} is supplied, the constructor
#' calls \code{\link{compile_model}} immediately with \code{max_order}.  Set
#' \code{max_order = 2L} (or higher) if you intend to call
#' \code{dm_solve(order = 2L)} or a higher-order likelihood.
#'
#' @param model       A \code{dynhr_mod} object from \code{\link{parse_mod}}, a
#'   path to a \code{.mod} file (parsed on the spot), or a
#'   \code{\link{dynhr_estimation_spec}} to hold as it is (then only
#'   \code{params}, \code{max_order} and \code{verbose} may be given).
#' @param data        Observation matrix (\eqn{T \times n_{\text{obs}}}) whose
#'   columns match \code{obs_vars}.  Optional: a model object with no data
#'   still supports \code{\link{dm_solve}} and \code{\link{dm_irf}}; its spec
#'   is validated once data is added.
#' @param obs_vars    Character vector of observed variable names.  Defaults
#'   to the model's \code{varobs} declaration when it has one.
#' @param compiled    A \code{dynhr_compiled} from \code{\link{compile_model}}.
#'   When \code{NULL} (default) the constructor compiles the model itself.
#' @param prior_spec  Prior spec data.frame from \code{\link{prior_spec}}.
#'   When \code{NULL} (default) the model's \code{estimated_params} block is
#'   used (none: \code{dm$prior_spec} is \code{NULL}).
#' @param me_variance Measurement-error variance
#'   (\code{likelihood$me_variance}). \code{NULL} (default) takes the field's
#'   default, the \code{me_variance} option (0 unless set).
#' @param likelihood  The likelihood: a type name (default
#'   \code{"gaussian"}), a \code{\link{likelihood_spec}} or a named list of its
#'   fields.
#' @param params      Named numeric vector of parameter values used by
#'   \code{\link{dm_solve}} and \code{\link{dm_irf}} at \code{at = "params"}.
#'   Defaults to the model's calibrated \code{param_values}.
#' @param max_order   Perturbation order to compile derivatives for
#'   (default \code{1L}); ignored when \code{compiled} is supplied.
#' @param verbose     Print parse/compile progress (default \code{FALSE}); also
#'   the default of the spec's \code{compute$verbose}.
#' @param mode,sampler,compute,outputs The spec's \code{\link{mode_spec}},
#'   \code{\link{sampler_spec}} (\code{NULL}: none yet),
#'   \code{\link{compute_spec}} and \code{\link{outputs_spec}}, or named lists
#'   of their fields.
#' @param options     As in \code{\link{dynhr_estimation_spec}}: \code{NULL}
#'   snapshots the current package options.
#' @param ...         Further likelihood arguments. A name that sets a typed
#'   likelihood field (e.g. \code{power}, the \code{power_posterior}
#'   field) is moved to that field; the rest are stored in
#'   \code{likelihood$extra} and forwarded to the log-posterior constructor
#'   (e.g. \code{order} for \code{likelihood = "cumulant"}).
#'
#' @return An object of class \code{dynhr_model}: a list with \code{spec},
#'   \code{params} and the initially-\code{NULL} result slots \code{steady},
#'   \code{solution}, \code{order}, \code{solution_at},
#'   \code{solution_params}, \code{log_post}, \code{mode}, \code{chains},
#'   \code{fit}, \code{diagnostics}, \code{irf}, \code{irf_at},
#'   \code{forecast} and \code{forecast_at}.
#'
#' @seealso \code{\link{update.dynhr_model}}, \code{\link{dm_solve}},
#'   \code{\link{dm_posterior}}, \code{\link{dm_mode}},
#'   \code{\link{dm_sample}}, \code{\link{dm_diagnostics}},
#'   \code{\link{dm_irf}}, \code{\link{dm_forecast}},
#'   \code{\link{dynhr_estimation_spec}}, \code{\link{run_estimation}}
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
                        prior_spec = NULL, me_variance = NULL,
                        likelihood = "gaussian", params = NULL,
                        max_order = 1L, verbose = FALSE,
                        mode = mode_spec(), sampler = sampler_spec("rwmh"),
                        compute = compute_spec(verbose = verbose),
                        outputs = outputs_spec(), options = NULL, ...) {
  ## Own the message epoch for this run: repeat-suppressed warnings
  ## (`.dynhr_warn(once = TRUE)`) are keyed within it and re-arm for the
  ## next run, and the close reports what it suppressed. A nested call
  ## inherits this epoch rather than opening a second one.
  .dynhr_run_epoch <- .dynhr_epoch("dynhr_model")
  on.exit(.dynhr_close_epoch(.dynhr_run_epoch), add = TRUE)

  ## ---- a spec, held as it is -------------------------------------------
  if (inherits(model, "dynhr_estimation_spec")) {
    given <- c(data = !is.null(data), obs_vars = !is.null(obs_vars),
               compiled = !is.null(compiled), prior_spec = !is.null(prior_spec),
               me_variance = !is.null(me_variance),
               likelihood = !missing(likelihood), mode = !missing(mode),
               sampler = !missing(sampler), compute = !missing(compute),
               outputs = !missing(outputs), options = !is.null(options),
               `...` = length(list(...)) > 0L)
    if (any(given))
      .dynhr_abort("dynhr_model: `model` is an estimation spec; edit it with ",
                   "update() first instead of passing ",
                   paste(names(given)[given], collapse = ", "), ".",
                   class = "dynhr_error_bad_argument")
    spec <- model
    if (is.null(spec$model$compiled %||% spec$model$solved$compiled))
      spec$model$compiled <- compile_model(spec$model$mod, verbose = verbose,
                                           max_order = spec$model$max_order)
    if (.dm_has_data(spec)) validate_spec(spec)
    return(.dm_new(spec, params %||% spec$model$mod$param_values, verbose))
  }

  ## ---- model: a parsed object or a path ---------------------------------
  if (is.character(model) && length(model) == 1L) {
    if (!file.exists(model))
      .dynhr_abort("dynhr_model: `model` looks like a file path but does not exist: ",
                   model, class = "dynhr_error_bad_argument")
    model <- parse_mod(model, verbose = verbose)
  }
  if (!inherits(model, "dynhr_mod"))
    .dynhr_abort("dynhr_model: `model` must be a dynhr_mod (from parse_mod()), a ",
                 "path to a .mod file or a dynhr_estimation_spec.",
                 class = "dynhr_error_bad_argument")

  ## ---- compiled derivatives (EAGER) -------------------------------------
  if (is.null(compiled)) {
    compiled <- compile_model(model, verbose = verbose,
                              max_order = as.integer(max_order))
  } else if (!inherits(compiled, "dynhr_compiled")) {
    .dynhr_abort("dynhr_model: `compiled` must be a dynhr_compiled object from ",
                 "compile_model().", class = "dynhr_error_bad_argument")
  }

  ## ---- observed variables ------------------------------------------------
  if (is.null(obs_vars) || length(obs_vars) == 0L)
    obs_vars <- model$obs_vars %||% model$varobs_names

  ## ---- data --------------------------------------------------------------
  if (!is.null(data)) {
    if (is.data.frame(data)) data <- as.matrix(data)
    if (!is.matrix(data))
      .dynhr_abort("dynhr_model: `data` must be a T x n_obs matrix or data.frame ",
                   "(rows = periods, columns = observables).",
                   class = "dynhr_error_bad_argument")
    if (!is.null(obs_vars) && !is.null(colnames(data))) {
      missing_cols <- setdiff(obs_vars, colnames(data))
      if (length(missing_cols))
        .dynhr_abort("dynhr_model: `data` has no column for observable(s): ",
                     paste(missing_cols, collapse = ", "),
                     class = "dynhr_error_bad_argument")
    }
    if (!is.null(obs_vars) && is.null(colnames(data)) &&
        ncol(data) != length(obs_vars))
      .dynhr_abort("dynhr_model: `data` has ", ncol(data), " unnamed columns but ",
                   length(obs_vars), " obs_vars were supplied.",
                   class = "dynhr_error_bad_argument")
  }

  ## ---- parameter values --------------------------------------------------
  if (is.null(params)) params <- model$param_values

  ## ---- the likelihood component --------------------------------------------
  ## me_variance NULL: the field's default rule (the option), as every other
  ## entry point resolves it (brief 28 S2). `...` keys that name a typed field
  ## (power, ...) are moved to that field by .spec_route_extras() below.
  lik <- if (is.character(likelihood) && length(likelihood) == 1L)
    likelihood_spec(likelihood) else .spec_as_sub(likelihood, "likelihood")
  if (!is.null(me_variance))
    lik <- .spec_build("likelihood", list(me_variance = me_variance), base = lik)
  dots <- list(...)
  if (length(dots))
    lik <- .spec_build("likelihood", list(extra = utils::modifyList(lik$extra, dots)),
                       base = lik)

  ## ---- the held spec ----------------------------------------------------------
  has_data <- !is.null(data)
  parts <- list(
    model      = .spec_model_part(model, compiled = compiled,
                                  max_order = max_order, prior_spec = prior_spec),
    data       = .spec_data_part(data),
    obs_vars   = if (length(obs_vars) || has_data)
                   .spec_resolve_obs_vars(obs_vars, model) else character(0),
    likelihood = lik,
    mode       = .spec_as_sub(mode, "mode"),
    sampler    = .spec_as_sampler(sampler),
    compute    = .spec_as_sub(compute, "compute"),
    outputs    = .spec_as_sub(outputs, "outputs"),
    options    = .spec_options_snapshot())
  parts <- .spec_route_extras(parts)
  parts <- .spec_apply_options(parts, options)
  .dm_new(.spec_assemble(parts, validate = has_data), params, verbose)
}


#' Edit a dynhr_model's spec or parameters
#'
#' Edits the estimation spec a \code{\link{dynhr_model}} holds (as
#' \code{update()} edits a \code{\link{dynhr_estimation_spec}}) and clears
#' every stored result that depends on what changed, so no result computed
#' under the old settings survives: a model, data, observables, likelihood,
#' options or \code{mode$extra} edit clears the log-posterior and everything
#' estimated from it; a \code{mode} edit the mode and the draws; a
#' \code{sampler} or \code{compute} edit the draws; a \code{model} or
#' \code{params} edit the solution, IRFs and forecasts. A solution, IRF or
#' forecast computed at the mode or the posterior mean is cleared with that
#' estimate. \code{outputs} edits clear only the diagnostics.
#'
#' @param object A \code{\link{dynhr_model}}.
#' @param model A \code{dynhr_mod} or a \code{.mod} path (compiled at once,
#'   to the object's \code{max_order}; \code{params} reset to its calibration
#'   and the prior spec to its \code{estimated_params} unless given).
#' @param data,obs_vars New data / observables.
#' @param likelihood,mode,sampler,compute,outputs A named list of fields to
#'   change (or a complete sub-spec that replaces the component;
#'   \code{likelihood} may also be a type name; \code{sampler = FALSE} removes
#'   the sampler).
#' @param options Options to merge into the spec's snapshot.
#' @param params Named numeric vector of parameter values (\code{dm$params}).
#' @param prior_spec A prior spec data.frame.
#' @param ... Not used (an error when non-empty).
#' @return The updated \code{dynhr_model}.
#' @seealso \code{\link{dynhr_model}}, \code{\link{dynhr_estimation_spec}}
#' @examples
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' dm <- dm_posterior(dynhr_model(
#'   system.file("extdata/models/nk_demo.mod", package = "dynhr"),
#'   data = Y, obs_vars = obs_vars))
#' dm2 <- update(dm, likelihood = list(me_variance = 1e-4))
#' is.null(dm2$log_post)   # the closure of the old likelihood is gone
#' @export
update.dynhr_model <- function(object, model = NULL, data = NULL,
                               obs_vars = NULL, likelihood = NULL, mode = NULL,
                               sampler = NULL, compute = NULL, outputs = NULL,
                               options = NULL, params = NULL,
                               prior_spec = NULL, ...) {
  .dm_check(object, "update")
  if (length(list(...)))
    .dynhr_abort("update.dynhr_model: unknown argument(s) ",
                 paste(names(list(...)), collapse = ", "), ".",
                 class = "dynhr_error_spec_unknown_field")
  if (!is.null(params) &&
      (!is.numeric(params) || is.null(names(params)) || anyNA(names(params))))
    .dynhr_abort("update.dynhr_model: `params` must be a named numeric vector.",
                 class = "dynhr_error_bad_argument")
  .dm_edit(object, model = model, data = data, obs_vars = obs_vars,
           likelihood = likelihood, mode = mode, sampler = sampler,
           compute = compute, outputs = outputs, options = options,
           params = params, prior_spec = prior_spec)
}


#' @rdname dynhr_model
#' @param x A \code{dynhr_model}.
#' @param name,i Field name.
#' @param value Value to store in a result slot.
#' @export
`$.dynhr_model` <- function(x, name) .dm_get(x, name)

#' @rdname dynhr_model
#' @export
`[[.dynhr_model` <- function(x, i, ...) .dm_get(x, i)

#' @rdname dynhr_model
#' @export
`$<-.dynhr_model` <- function(x, name, value) .dm_assign(x, name, value)

#' @rdname dynhr_model
#' @export
`[[<-.dynhr_model` <- function(x, i, value) .dm_assign(x, i, value)


# ============================================================================
# Points a solution is computed at (params / mode / posterior mean)
# ============================================================================

## The estimated vector at `at` ("mode" / "posterior_mean"), or a named error.
.dm_theta_at <- function(dm, at, verb) {
  if (identical(at, "mode")) {
    th <- dm$mode$theta_mode
    if (is.null(th)) .dm_require(verb, "mode", "dm_mode", at = at)
    return(th)
  }
  ch <- dm$chains$chain
  if (is.null(ch) || !nrow(ch)) .dm_require(verb, "chains", "dm_sample", at = at)
  colMeans(ch)
}

## Solve the model at the mode / posterior mean through the same
## theta -> params -> steady state -> decision-rule pipeline the runner's
## diagnostics use. Returns list(steady, solution, params).
.dm_solve_at <- function(dm, at, order, verb, ...) {
  theta <- .dm_theta_at(dm, at, verb)
  p  <- apply_theta_to_params(dm$model, theta, params = dm$params)
  ss <- solve_steady_state(dm$model, dm$compiled, p, verbose = FALSE)
  if (is.null(ss) || !isTRUE(ss$converged))
    .dynhr_abort(verb, "(at = \"", at, "\"): the steady state does not converge ",
                 "at the ", if (identical(at, "mode")) "posterior mode"
                 else "posterior mean", ".",
                 class = "dynhr_error_dm_solve_failed")
  p <- ss$params %||% p
  dr <- solve_perturbation(dm$model, dm$compiled, ss$ss, p,
                           order = as.integer(order), ...)
  ss$values <- ss$values %||% ss$ss
  list(steady = ss, solution = dr, params = p)
}

## The solution (and its parameters) at `at` for dm_irf() / dm_forecast():
## the stored one when it was solved there, else solved afresh at the stored
## order (not stored). A missing dm_solve() is an error naming it.
.dm_point <- function(dm, at, verb) {
  if (is.null(dm$solution)) .dm_require(verb, "solution", "dm_solve")
  if (identical(dm$solution_at %||% "params", at))
    return(list(solution = dm$solution,
                params = dm$solution_params %||% dm$params))
  ord <- dm$order %||% 1L
  if (!identical(at, "params")) return(.dm_solve_at(dm, at, ord, verb))
  ss <- solve_steady(dm$compiled, dm$params, endo_names = dm$model$var_names,
                     exo_names = dm$model$varexo_names, verbose = FALSE)
  list(solution = solve_perturbation(dm$model, dm$compiled, ss$values,
                                     dm$params, order = ord),
       params = dm$params)
}


# ============================================================================
# Verbs
# ============================================================================

#' Solve the model's perturbation decision rules
#'
#' Computes the steady state (stored as \code{$steady}) and calls
#' \code{\link{solve_perturbation}} at \code{order}, storing the resulting
#' \code{DecisionRules*} object in \code{$solution}, the order in
#' \code{$order}, the point in \code{$solution_at} and the parameter vector in
#' \code{$solution_params}.
#'
#' @section Solution point:
#' \code{at = "params"} (default) solves at \code{dm$params} (the calibration
#' unless changed), exactly as \code{solve_steady()} +
#' \code{solve_perturbation()} do. \code{"mode"} and \code{"posterior_mean"}
#' solve at the posterior mode (\code{\link{dm_mode}}) or the mean of the
#' draws (\code{\link{dm_sample}}), mapped to parameters by
#' \code{\link{apply_theta_to_params}} (estimated shock standard deviations
#' and correlations included) and solved through
#' \code{solve_steady_state()}, as the estimation runner's diagnostics do;
#' they are an error naming the verb to run when that estimate is absent.
#'
#' @param dm    A \code{\link{dynhr_model}}
#' @param order Perturbation order (default \code{1L}).  The object's compiled
#'   derivatives must be deep enough -- see \code{max_order} in
#'   \code{\link{dynhr_model}}.
#' @param at    Where to solve: \code{"params"}, \code{"mode"} or
#'   \code{"posterior_mean"} (see Solution point).
#' @param ...   Further arguments for \code{\link{solve_perturbation}}
#'   (e.g. \code{Sigma_e}, \code{loglinear}).
#' @return The updated \code{dynhr_model}.
#' @seealso \code{\link{dynhr_model}}, \code{\link{solve_perturbation}}
#' @examples
#' dm <- dynhr_model(system.file("extdata/models/rbc.mod", package = "dynhr"))
#' dm <- dm_solve(dm)
#' class(dm$solution)
#' @export
dm_solve <- function(dm, order = 1L, at = c("params", "mode", "posterior_mean"),
                     ...) {
  .dm_check(dm, "dm_solve")
  at <- match.arg(at)
  order <- as.integer(order)
  if (identical(at, "params")) {
    ss <- dm$steady
    if (is.null(ss) || !identical(dm$solution_at %||% "params", "params"))
      ss <- solve_steady(dm$compiled, dm$params,
                         endo_names = dm$model$var_names,
                         exo_names  = dm$model$varexo_names,
                         verbose    = FALSE)
    pt <- list(steady = ss,
               solution = solve_perturbation(dm$model, dm$compiled, ss$values,
                                             dm$params, order = order, ...),
               params = dm$params)
  } else {
    pt <- .dm_solve_at(dm, at, order, "dm_solve", ...)
  }
  .dm_set_many(dm, list(steady = pt$steady, solution = pt$solution,
                        order = order, solution_at = at,
                        solution_params = pt$params))
}


#' Build the log-posterior closure for a dynhr_model
#'
#' Builds the log-posterior of the held spec with the builder
#' \code{\link{run_estimation}} uses, and stores the closure in
#' \code{$log_post}: the spec's likelihood fields (\code{me_variance},
#' \code{lik_init}, \code{power_posterior}, ...) and the model's
#' \code{heteroskedastic_shocks}, \code{filter_tunes} and
#' \code{stochastic_volatility} blocks and OBC constraints are applied, so the
#' closure is the runner's objective.
#'
#' Arguments in \code{...} edit the held likelihood (they are recorded in the
#' returned object's spec): a likelihood field by name, or by the constructor
#' argument that sets it (\code{power} is \code{power_posterior},
#' \code{likelihood} is \code{type}); any other name is added to
#' \code{likelihood$extra}. Such an edit clears the results computed under
#' the old likelihood (see \code{\link{update.dynhr_model}}).
#'
#' @param dm  A \code{\link{dynhr_model}} carrying \code{data} and a prior spec
#' @param ... Likelihood edits (see Details).
#' @return The updated \code{dynhr_model}, with \code{$log_post} set to the closure
#'   \code{function(theta) list(logpost, loglik, logprior)}.
#' @seealso \code{\link{dynhr_model}}, \code{\link{likelihood_spec}},
#'   \code{\link{run_estimation}}
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
  .dm_need_inputs(dm, "dm_posterior")
  dots <- list(...)
  if (length(dots)) {
    if (is.null(names(dots)) || any(!nzchar(names(dots))))
      .dynhr_abort("dm_posterior(): every argument in `...` must be named.",
                   class = "dynhr_error_bad_argument")
    dm <- .dm_edit(dm, likelihood = .dm_likelihood_edits(dm$spec$likelihood, dots))
  }
  build <- .est_build_posterior(validate_spec(dm$spec))
  .dm_set(dm, "log_post", build$log_post_fn)
}


#' Find the posterior mode of a dynhr_model
#'
#' Runs the mode stage of \code{\link{run_estimation}} on the held spec
#' (\code{outputs$form = "mode"}, no sampler) and stores the
#' \code{dynhr_mode_result} in \code{$mode}: \code{theta_mode}, the optimiser
#' result in \code{$mode$mode} (\code{logpost}, \code{convergence}, ...), the
#' proposal covariance \code{Sigma_prop} a later \code{\link{dm_sample}}
#' uses, and the run record. The mode-stage settings are the spec's
#' \code{mode} component (\code{method}, \code{n_iter},
#' \code{transform_params}, \code{theta_init} -- INITVAL / prior means when
#' \code{NULL} --, ...); the seed and progress output its \code{compute}
#' component. The log-posterior is built as \code{\link{dm_posterior}}
#' builds it; \code{$log_post} is filled when empty.
#'
#' @param dm         A \code{\link{dynhr_model}} carrying \code{data} and a
#'   prior spec
#' @param theta_init Named numeric starting vector (the \code{mode$theta_init}
#'   field); \code{NULL} keeps the spec's.
#' @param ...        Fields of \code{\link{mode_spec}} (e.g. \code{n_iter},
#'   \code{method}, \code{transform_params}, \code{options}) or of
#'   \code{\link{compute_spec}} (e.g. \code{seed}, \code{verbose}), recorded
#'   in the returned object's spec.
#' @return The updated \code{dynhr_model}, with \code{$mode} set.
#' @seealso \code{\link{dynhr_model}}, \code{\link{run_estimation}},
#'   \code{\link{mode_spec}}
#' @examples
#' \donttest{
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' dm <- dynhr_model(system.file("extdata/models/nk_demo.mod", package = "dynhr"),
#'                   data = Y, obs_vars = obs_vars)
#' dm <- dm_mode(dm, n_iter = 50L)
#' dm$mode$theta_mode
#' }
#' @export
dm_mode <- function(dm, theta_init = NULL, ...) {
  .dm_check(dm, "dm_mode")
  .dm_need_inputs(dm, "dm_mode")
  ed <- .dm_split_args(list(...), c("mode", "compute"), "dm_mode")
  if (!is.null(theta_init)) ed$mode$theta_init <- theta_init
  dm <- .dm_edit(dm, mode = ed$mode, compute = ed$compute)
  mr <- run_estimation(update(validate_spec(dm$spec), sampler = FALSE,
                              outputs = list(form = "mode")))
  dm <- .dm_drop_after(dm, "mode")
  dm <- .dm_set(dm, "mode", mr)
  if (is.null(dm$log_post)) dm <- .dm_set(dm, "log_post", mr$log_post_fn)
  dm
}


#' Sample the posterior of a dynhr_model
#'
#' Runs the sampler stage of \code{\link{run_estimation}} on the held spec
#' (\code{outputs$form = "posterior"}) and stores the draws as a
#' \code{\link{dynhr_chains}} in \code{$chains} (all chains pooled, with
#' \code{chain_list} for several) and the \code{dynhr_posterior_result} in
#' \code{$fit} (convergence / R-hat, the posterior mean, the run record).
#' Every sampler the runner runs is available, with its defaults: several
#' chains, \code{compute$seed} (\code{seed = NULL} keeps the ambient RNG
#' stream), the proposal covariance from the mode stage, the
#' unconstrained-space transform.
#'
#' @section Starting point:
#' The point-started samplers (\code{"rwmh"}, \code{"pmmh"}, \code{"nuts"},
#' \code{"hmc"}, \code{"mala"}, \code{"chees"}) start at the mode and use its
#' proposal covariance (unless \code{Sigma_prop} is given): the stored
#' \code{$mode} when \code{\link{dm_mode}} has run, else the runner's mode
#' stage runs first and its result is stored in \code{$mode}. The
#' prior-initialised samplers (\code{"smc"}, \code{"dsmh"}, \code{"dime"})
#' need no mode and run none. An SMC cloud is resampled to equal weights
#' (\code{$chains$chain}; the weighted particles stay in
#' \code{$chains$particles} / \code{$chains$smc_weights}), so the diagnostics
#' and posterior means see the weighted posterior.
#'
#' @param dm      A \code{\link{dynhr_model}} carrying \code{data} and a prior
#'   spec
#' @param sampler A method name (\code{"rwmh"}, \code{"pmmh"}, \code{"nuts"},
#'   \code{"hmc"}, \code{"mala"}, \code{"chees"}, \code{"smc"}, \code{"dsmh"},
#'   \code{"dime"}) or a \code{\link{sampler_spec}}. \code{NULL} (default)
#'   keeps the held sampler (\code{"rwmh"} unless set).
#' @param ...     Fields of that sampler (\code{n_draws}, \code{n_warmup},
#'   \code{n_chains}, \code{n_particles}, \code{Sigma_prop},
#'   \code{transform_params}, \code{extra}, ...; see
#'   \code{\link{sampler_spec}}) or of \code{\link{compute_spec}}
#'   (\code{seed}, \code{verbose}, \code{parallel}, ...), recorded in the
#'   returned object's spec. For the held sampler's method the given fields
#'   are merged into it; another method starts from that method's defaults.
#' @return The updated \code{dynhr_model}, with \code{$chains} and \code{$fit}
#'   set (and \code{$mode} when the mode stage ran).
#' @seealso \code{\link{dynhr_model}}, \code{\link{run_estimation}},
#'   \code{\link{sampler_spec}}
#' @examples
#' \donttest{
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' dm <- dynhr_model(system.file("extdata/models/nk_demo.mod", package = "dynhr"),
#'                   data = Y, obs_vars = obs_vars)
#' dm <- dm_mode(dm, n_iter = 50L)
#' dm <- dm_sample(dm, n_draws = 200L, n_warmup = 100L, n_chains = 2L, seed = 1L)
#' dm$chains
#' dm$fit$convergence$rhat
#' }
#' @export
dm_sample <- function(dm, sampler = NULL, ...) {
  .dm_check(dm, "dm_sample")
  .dm_need_inputs(dm, "dm_sample")
  held <- dm$spec$sampler
  held <- if (inherits(held, "dynhr_sampler_spec")) held
  if (inherits(sampler, "dynhr_sampler_spec")) {
    method <- sampler$method
    base   <- sampler
  } else {
    if (!is.null(sampler) && !(is.character(sampler) && length(sampler) == 1L))
      .dynhr_abort("dm_sample(): `sampler` must be a method name or a ",
                   "sampler_spec().", class = "dynhr_error_bad_argument")
    method <- tolower(sampler %||% held$method %||% "rwmh")
    if (!method %in% .spec_samplers)
      .dynhr_abort("dm_sample(): unknown sampler \"", method, "\". Valid: ",
                   paste(.spec_samplers, collapse = ", "), ".",
                   class = "dynhr_error_bad_argument")
    base <- if (identical(held$method, method)) held
  }
  ed <- .dm_split_args(list(...), c("sampler", "compute"), "dm_sample",
                       method = method)
  dm <- .dm_edit(dm, sampler = .spec_build("sampler", ed$sampler,
                                           method = method, base = base),
                 compute = ed$compute)

  ## A point-started sampler starts at the stored mode (its proposal
  ## covariance included); without one the runner runs the mode stage.
  use_mode <- !method %in% .spec_prior_init_samplers && !is.null(dm$mode)
  run <- if (use_mode)
    update(validate_spec(dm$spec), mode = list(result = dm$mode),
           outputs = list(form = "posterior"))
  else
    update(validate_spec(dm$spec), outputs = list(form = "posterior"))
  fit <- run_estimation(run)
  co  <- .est_chains_object(fit$chains[[1L]], method)
  dm  <- .dm_drop_after(dm, "sample")
  dm  <- .dm_set_many(dm, list(chains = co$chains, fit = fit))
  mr  <- fit$mode_result
  if (!use_mode && is.null(dm$mode) && !isTRUE(mr$meta$mode_skipped))
    dm <- .dm_set(dm, "mode", mr)
  if (is.null(dm$log_post)) dm <- .dm_set(dm, "log_post", mr$log_post_fn)
  dm
}


#' Run the dynhr diagnostic battery on a dynhr_model
#'
#' Forwards whichever of the object's pieces are present -- model, compiled
#' derivatives, decision rules (at the point \code{\link{dm_solve}} solved
#' them) and their parameters, data, observed names, prior spec, the
#' posterior draws (equally weighted; one matrix per chain as
#' \code{chains_list}) and the mode -- to \code{\link{run_diagnostics}}, and
#' stores the returned list in \code{$diagnostics}.  Stages whose inputs are
#' absent are skipped by the battery itself, so this works with a
#' partially-run pipeline.
#'
#' @param dm  A \code{\link{dynhr_model}}, ideally after
#'   \code{\link{dm_solve}} and \code{\link{dm_sample}}
#' @param ... Further arguments for \code{\link{run_diagnostics}}
#'   (e.g. \code{verbose}, \code{report}).  These override the forwarded fields.
#' @return The updated \code{dynhr_model}, with \code{$diagnostics} set.
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
  ch <- dm$chains
  base <- list(model       = dm$model,
               compiled    = dm$compiled,
               dr          = dm$solution,
               ss          = dm$steady,
               params      = dm$solution_params %||% dm$params,
               priors      = dm$prior_spec,
               data        = dm$data,
               obs_names   = dm$obs_vars,
               draws       = ch$chain,
               chains_list = if (!is.null(ch$chain)) {
                               if (!is.null(ch$chain_list))
                                 lapply(ch$chain_list, function(c) c$chain %||% c)
                               else list(ch$chain)
                             },
               theta_mode  = dm$mode$theta_mode)
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
#' \code{$irf} (and the point in \code{$irf_at}).
#'
#' @param dm        A \code{\link{dynhr_model}} whose \code{$solution} is set
#'   (see \code{\link{dm_solve}})
#' @param n_periods IRF horizon (default \code{40L})
#' @param at        \code{"params"} (default), \code{"mode"} or
#'   \code{"posterior_mean"}: the stored solution when \code{\link{dm_solve}}
#'   solved it there, else the model is solved at that point (at the stored
#'   order, default solver settings) for this call. See \code{\link{dm_solve}}.
#' @param ...       Further arguments for the underlying IRF routine
#'   (e.g. \code{shock_size}, and \code{pruning} at order 2).
#' @return The updated \code{dynhr_model}, with \code{$irf} set.
#' @seealso \code{\link{dynhr_model}}, \code{\link{compute_irfs}}
#' @examples
#' dm <- dm_irf(dm_solve(dynhr_model(
#'   system.file("extdata/models/rbc.mod", package = "dynhr"))),
#'   n_periods = 8L)
#' names(dm$irf)
#' @export
dm_irf <- function(dm, n_periods = 40L, at = c("params", "mode", "posterior_mean"),
                   ...) {
  .dm_check(dm, "dm_irf")
  at <- match.arg(at)
  pt <- .dm_point(dm, at, "dm_irf")
  ord <- dm$order %||% 1L
  fn <- switch(as.character(ord),
               "1" = compute_irfs,
               "2" = compute_irfs_order2,
               "3" = compute_irfs_order3,
               .dynhr_abort("dm_irf(): no IRF routine for perturbation order ", ord,
                            "; call compute_irfs_order4()/compute_irfs_order5() on ",
                            "`dm$solution` directly.",
                            class = "dynhr_error_bad_argument"))
  .dm_set_many(dm, list(
    irf = fn(pt$solution, dm$model, n_periods = as.integer(n_periods),
             params = pt$params, ...),
    irf_at = at))
}


#' Conditional forecast from a dynhr_model
#'
#' Re-threads the object's \code{model}, solution, \code{data},
#' \code{obs_vars} and \code{compiled} into \code{\link{conditional_forecast}}
#' and stores the result in \code{$forecast} (and the point in
#' \code{$forecast_at}).
#'
#' @param dm         A \code{\link{dynhr_model}} carrying \code{data} and a
#'   \code{$solution} (see \code{\link{dm_solve}})
#' @param conditions Conditioning information, as in
#'   \code{\link{conditional_forecast}} (\code{NULL} gives an unconditional
#'   forecast)
#' @param horizon    Forecast horizon (default \code{8L})
#' @param at         \code{"params"} (default), \code{"mode"} or
#'   \code{"posterior_mean"}, as in \code{\link{dm_irf}}. Away from
#'   \code{"params"} the model's parameter values (hence the shock covariance
#'   of the filter and the forecast) are those at that point.
#' @param ...        Further arguments for \code{\link{conditional_forecast}}
#'   (e.g. \code{type}, \code{method}, \code{n_draws}).
#' @return The updated \code{dynhr_model}, with \code{$forecast} set.
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
#' dim(dm$forecast$paths_point)
#' }
#' @export
dm_forecast <- function(dm, conditions = NULL, horizon = 8L,
                        at = c("params", "mode", "posterior_mean"), ...) {
  .dm_check(dm, "dm_forecast")
  at <- match.arg(at)
  if (is.null(dm$solution)) .dm_require("dm_forecast", "solution", "dm_solve")
  if (is.null(dm$data))     .dm_require_field("dm_forecast", "data")
  pt <- .dm_point(dm, at, "dm_forecast")
  model <- dm$model
  if (!identical(at, "params")) model$param_values <- pt$params
  .dm_set_many(dm, list(
    forecast = conditional_forecast(model, pt$solution, dm$data,
                                    conditions = conditions,
                                    horizon    = as.integer(horizon),
                                    obs_vars   = dm$obs_vars,
                                    compiled   = dm$compiled, ...),
    forecast_at = at))
}


# ============================================================================
# S3 methods
# ============================================================================

#' Print a dynhr_model object
#'
#' Compact one-screen header: model size, data dimensions, observables, the
#' held spec's sampler, and which pipeline steps have been run.
#'
#' @param x   A \code{\link{dynhr_model}}
#' @param ... Ignored
#' @return \code{x}, invisibly
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
  d <- x$data
  if (is.null(d)) {
    cat("  Data       : <none>\n")
  } else {
    cat(sprintf("  Data       : %d periods x %d observables\n",
                nrow(d), ncol(d)))
  }
  cat(sprintf("  Observables: %s\n",
              if (length(x$obs_vars)) paste(x$obs_vars, collapse = ", ")
              else "<none>"))
  cat(sprintf("  Priors     : %s\n",
              if (is.null(x$prior_spec)) "<none>"
              else sprintf("%d estimated parameters", nrow(x$prior_spec))))
  cat(sprintf("  Likelihood : %s (me_variance = %s)\n",
              x$likelihood, format(x$me_variance)))
  s <- x$spec$sampler
  cat(sprintf("  Sampler    : %s (seed %s)\n",
              if (inherits(s, "dynhr_sampler_spec")) s$method
              else if (is.null(s)) "<none>" else "sequence",
              format(x$spec$compute$seed %||% "NULL")))

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
#' block for every step that has been run (perturbation order, point and
#' state count, the log-posterior at the prior mean, the mode's
#' log-posterior, the draws and R-hat, the diagnostic names).
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
    out$order   <- object$order
    out$at      <- object$solution_at %||% "params"
    out$n_state <- length(object$solution$state_idx)
    cat(sprintf("  dm_solve       : order %d at %s, %d states, class %s\n",
                out$order, out$at, out$n_state, class(object$solution)[1]))
  }
  if (!is.null(object$log_post) && !is.null(object$prior_spec)) {
    th <- stats::setNames(object$prior_spec$mean, object$prior_spec$name)
    lp <- tryCatch(object$log_post(th)$logpost, error = function(e) .dynhr_reraise_bug(e, NA_real_))
    out$logpost_prior_mean <- lp
    cat(sprintf("  dm_posterior   : logpost at prior mean = %.6f\n", lp))
  }
  if (!is.null(object$mode)) {
    md <- object$mode
    out$mode_logpost <- md$mode$logpost %||% md$logpost
    cat(sprintf("  dm_mode        : logpost = %.6f (%s, convergence %s)\n",
                out$mode_logpost %||% NA_real_,
                md$meta$method %||% md$method %||% "?",
                format(md$mode$convergence %||% md$convergence %||% NA)))
  }
  if (!is.null(object$chains)) {
    out$n_draws <- nrow(object$chains$chain)
    out$sampler <- object$chains$sampler
    rh <- object$fit$convergence$rhat
    out$max_rhat <- if (length(rh)) max(rh, na.rm = TRUE)
    cat(sprintf("  dm_sample      : %s, %d draws x %d parameters%s\n",
                toupper(out$sampler %||% "?"), out$n_draws,
                ncol(object$chains$chain),
                if (length(rh)) sprintf(", max R-hat %.3f", out$max_rhat) else ""))
  }
  if (!is.null(object$diagnostics)) {
    out$diagnostics <- names(object$diagnostics)
    cat(sprintf("  dm_diagnostics : %d stages (%s)\n",
                length(out$diagnostics),
                paste(utils::head(out$diagnostics, 6L), collapse = ", ")))
  }
  if (!is.null(object$irf)) {
    out$irf <- names(object$irf)
    cat(sprintf("  dm_irf         : %d shock(s) at %s\n", length(out$irf),
                object$irf_at %||% "params"))
  }
  if (!is.null(object$forecast)) {
    cat(sprintf("  dm_forecast    : present (at %s)\n",
                object$forecast_at %||% "params"))
    out$forecast <- TRUE
  }
  if (!length(out)) cat("  <no step has been run>\n")
  invisible(out)
}
