## R/run-record.R
## ---------------------------------------------------------------------------
## Reproducible run records and dynhr_rerun().
##
## Every result of the estimation entry points run_mode_finding(),
## run_posterior_estimation() and run_full_estimation() carries
## `$run_record`, a `dynhr_run_record`: everything needed to replay the call
## -- the fully resolved arguments, the effective package options, the RNG
## state, the build that produced it -- plus hashes of the inputs.
## dynhr_rerun() replays it.
##
## The record is also the INPUT form of the planned estimation-spec object
## (E5 step C), so its layout is versioned (`schema_version`) and documented
## in ?dynhr_run_record. Bump `.dynhr_run_record_schema` on any change a
## reader must know about, and keep dynhr_rerun() able to read older schemas.
##
## Life cycle inside an entry point:
##   .rr <- .dynhr_rr_begin("<fn>", environment(), list(...))   # FIRST line,
##          # before the body touches any argument or the RNG
##   ...
##   result$run_record <- .dynhr_rr_finish(.rr, ...)             # before return
## ---------------------------------------------------------------------------

## Schema 2: the record also carries `spec`, the estimation spec the
## run executed (caches stripped), and `provenance$integrity` after a
## checkpoint resume that crossed a mismatch.
.dynhr_run_record_schema <- 2L

## The functions a record may name. dynhr_rerun() refuses anything else, so a
## record read back from an untrusted RDS cannot make it call an arbitrary
## function.
.dynhr_rr_entry_points <- c("run_mode_finding", "run_posterior_estimation",
                            "run_full_estimation", "run_estimation")

## Hash of an R object: sha256 via digest when installed (it is in Suggests,
## NOT Imports), otherwise base-R md5 of the serialised object. The algorithm
## is prefixed ("sha256:" / "md5:") so hashes from the two paths are never
## compared as if they were the same kind.
.dynhr_rr_hash <- function(x) {
  if (requireNamespace("digest", quietly = TRUE))
    return(paste0("sha256:", digest::digest(x, algo = "sha256")))
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, compress = FALSE, version = 3L)
  paste0("md5:", unname(tools::md5sum(f)))
}

## Structural model hash: the parsed model with its source path blanked, so the
## same model read from two locations (or from inline text) hashes the same.
.dynhr_rr_model_hash <- function(model) {
  if (is.null(model)) return(NA_character_)
  if (is.list(model) && !is.null(model$source_file))
    model$source_file <- NA_character_
  .dynhr_rr_hash(model)
}

## md5 of a file's bytes, NA when absent.
.dynhr_rr_file_md5 <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path))
    return(NA_character_)
  unname(tools::md5sum(path))
}

## Build identity + platform, cheap enough to take on every run (no
## sessionInfo(), no shelling out).
.dynhr_rr_provenance <- function() {
  ns   <- asNamespace("dynhr")
  path <- getNamespaceInfo(ns, "path")
  blas <- unname(extSoftVersion()["BLAS"])
  list(
    version    = unname(getNamespaceVersion(ns)),
    git_commit = .dynhr_git_stamp_at(path),
    dev_load   = exists(".__DEVTOOLS__", envir = ns, inherits = FALSE),
    r_version  = R.version.string,
    platform   = R.version$platform,
    os         = paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
    blas       = if (length(blas) == 1L && nzchar(blas)) blas else NA_character_,
    lapack     = La_library(),
    lapack_version = La_version()
  )
}

## Current default of formal `nm` of `f`, evaluated in the dynhr namespace.
## Returns list(ok = FALSE) when the default is absent or refers to other
## formals / unknown symbols (then it cannot be evaluated out of context).
.dynhr_rr_current_default <- function(f, nm) {
  fmls <- formals(f)
  if (!nm %in% names(fmls)) return(list(ok = FALSE))
  ## test the empty default BEFORE binding it: a variable holding the
  ## missing-argument marker errors when read
  if (identical(fmls[[nm]], quote(expr = ))) return(list(ok = FALSE))
  expr <- fmls[[nm]]
  ns   <- asNamespace("dynhr")
  vars <- all.vars(expr)
  if (length(intersect(vars, names(fmls))) ||
      !all(vapply(vars, exists, logical(1), envir = ns)))
    return(list(ok = FALSE))
  list(ok = TRUE, value = eval(expr, envir = ns))
}

# ---------------------------------------------------------------------------
# Capture
# ---------------------------------------------------------------------------

## Snapshot an entry point's call at entry. `env` is the entry point's frame,
## `dots` its evaluated `list(...)`. Forces every argument promise (all
## defaults of the entry points are constants) and records which formals the
## caller supplied: missing() is preserved on replay by omitting an unsupplied
## argument whose recorded value still equals the current default.
.dynhr_rr_begin <- function(fn, env, dots = list()) {
  f    <- get(fn, envir = asNamespace("dynhr"), inherits = FALSE)
  fmls <- formals(f)
  nms  <- setdiff(names(fmls), "...")
  args <- list()
  supplied <- character(0)
  for (nm in nms) {
    is_missing <- eval(call("missing", as.name(nm)), envir = env)
    if (!is_missing) supplied <- c(supplied, nm)
    ## a required formal the caller left out: nothing to record (the entry
    ## point errors on it itself)
    if (is_missing && identical(fmls[[nm]], quote(expr = ))) next
    args[nm] <- list(get(nm, envir = env, inherits = FALSE))
  }
  dot_names <- names(dots)
  if (length(dots)) {
    if (is.null(dot_names)) dot_names <- rep("", length(dots))
    args <- c(args, dots)
    supplied <- c(supplied, dot_names[nzchar(dot_names)])
  }
  base_all <- options()
  seed <- if ("seed" %in% nms) args$seed
          else if (identical(fn, "run_estimation")) args$spec$compute$seed
  list(
    schema_version = .dynhr_run_record_schema,
    created        = Sys.time(),
    fn             = fn,
    args           = args,
    supplied       = supplied,
    dot_names      = if (length(dots)) dot_names else character(0),
    inputs         = list(),
    stripped       = character(0),
    rebuild        = list(),
    options        = dynhr_get_options(effective = TRUE),
    base_options   = base_all[grep("^dynhr\\.", names(base_all))],
    rng            = list(kind = RNGkind(), seed = seed,
                          random_seed = .rng_snapshot()),
    provenance     = .dynhr_rr_provenance(),
    hashes         = list()
  )
}

## Complete a record at the end of an entry point: by-value substitutions,
## stripping of rebuildable heavy objects, hashes, parallel provenance.
## `model` / `data`: the parsed model and the data matrix the run actually used
## (run_full_estimation, which may have received a path for either).
## `spec`: the estimation spec that ran (schema 2); `integrity`: the
## checkpoint-resume verdict, when there was one.
.dynhr_rr_finish <- function(rec, model = NULL, data = NULL, spec = NULL,
                             integrity = NULL, result = NULL) {
  a <- rec$args
  switch(rec$fn,
    run_estimation = {
      ## the spec is kept once, as rec$spec (below)
      rec$args["spec"] <- list(NULL)
      rec$stripped <- c(rec$stripped, "spec (kept as $spec)")
      spec  <- spec %||% a$spec
      model <- model %||% a$spec$model$mod
      data  <- data %||% a$spec$data$value
    },
    run_mode_finding = {
      solved <- a$solved
      if (!is.null(solved$compiled)) {
        rec$rebuild$solved_compiled <-
          list(max_order = solved$compiled$max_order %||% 1L)
        solved$compiled <- NULL
        rec$args["solved"] <- list(solved)
        rec$stripped <- c(rec$stripped, "solved$compiled")
      }
      model <- solved$model
      data  <- a$data
    },
    run_posterior_estimation = {
      mr <- a$mode_result
      rec$args["mode_result"] <- list(.dynhr_rr_mode_ref(mr))
      rec$stripped <- c(rec$stripped, "mode_result")
      model <- mr$solved$model
      data  <- mr$data
    },
    run_full_estimation = {
      if (!is.null(a$compiled)) {
        rec$args["compiled"] <- list(NULL)
        rec$stripped <- c(rec$stripped, "compiled")
      }
      if (is.null(a$model) && !is.null(model)) rec$inputs$model <- model
      if (is.character(a$data) && !is.null(data)) rec$inputs$data <- data
      rec$hashes$mod_file <- .dynhr_rr_file_md5(a$mod_file)
      if (is.character(a$data)) rec$hashes$data_file <- .dynhr_rr_file_md5(a$data)
    })
  rec$hashes$data  <- if (is.null(data)) NA_character_ else .dynhr_rr_hash(data)
  rec$hashes$model <- .dynhr_rr_model_hash(model)
  est <- identical(rec$fn, "run_estimation")
  par <- isTRUE(if (est) a$spec$compute$parallel else a$parallel)
  rec$provenance$parallel <- par
  rec$provenance$n_cores  <- if (par)
    .mirai_n_cores(if (est) a$spec$compute$n_cores else a$n_cores) else 1L
  rec$provenance$mirai_version <- if (par)
    as.character(utils::packageVersion("mirai")) else NA_character_
  if (!is.null(integrity)) rec$provenance$integrity <- integrity
  if (!is.null(spec)) {
    rec$spec <- .dynhr_rr_strip_spec(spec)
    ## run_posterior_estimation(): args$mode_result already holds the mode
    ## reference (with the mode run's record); point at it, do not store it
    ## twice
    if (identical(rec$fn, "run_posterior_estimation") &&
        inherits(rec$spec$mode$result, "dynhr_mode_ref")) {
      md <- rec$spec$mode
      md["result"] <- list(structure(list(in_args = TRUE),
                                     class = "dynhr_mode_ref"))
      rec$spec["mode"] <- list(md)
    }
  }
  class(rec) <- "dynhr_run_record"
  if (!is.null(result)) rec <- .dynhr_rr_set_resolved(rec, result)
  rec
}

## Copy what the run resolved at run time into the record: the per-stage
## analytic-gradient method (`result$resolved`, filled by the sampler stage
## from each gradient closure's attr(, "grad_method")). A result with no
## gradient-based stage leaves the record untouched.
.dynhr_rr_set_resolved <- function(rec, result) {
  gm <- result$resolved$grad_method
  if (!length(gm)) return(rec)
  rec$resolved <- list(
    grad_method           = gm,
    grad_method_requested = result$resolved$grad_method_requested)
  rec
}

## "NUTS=adjoint_solution (auto)" -- one entry per gradient-based stage, the
## requested method in brackets when it was resolved to another; NULL when no
## stage built an analytic gradient.
.grad_method_record_line <- function(resolved) {
  gm <- resolved$grad_method
  if (!length(gm)) return(NULL)
  req <- resolved$grad_method_requested
  req <- if (length(req) == length(gm)) req else rep(NA_character_, length(gm))
  paste(sprintf("%s=%s%s", names(gm), gm,
                ifelse(!is.na(req) & req != gm, paste0(" (", req, ")"), "")),
        collapse = ", ")
}

## A spec as a record keeps it: the caches (compiled model, solved object)
## are dropped (rebuilt from the model) and a precomputed mode result is
## replaced by its dynhr_mode_ref (the mode run's own record + the pieces the
## sampler consumed). The hashes stay those of the spec that ran.
.dynhr_rr_strip_spec <- function(spec) {
  if (!inherits(spec, "dynhr_estimation_spec")) return(spec)
  spec$model["compiled"] <- list(NULL)
  spec$model["solved"]   <- list(NULL)
  mr <- spec$mode$result
  if (!is.null(mr) && !inherits(mr, "dynhr_mode_ref")) {
    md <- spec$mode
    md["result"] <- list(.dynhr_rr_mode_ref(mr))
    spec["mode"] <- list(md)
  }
  spec
}

## What run_posterior_estimation() records for its `mode_result`: the mode
## run's OWN record (so a replay re-runs mode finding and rebuilds the full
## object, closures and all) plus the pieces the sampler consumes, kept for
## provenance and for the drift check after a replayed mode run.
.dynhr_rr_mode_ref <- function(mode_result) {
  rr <- mode_result$run_record
  structure(list(
    run_record = if (inherits(rr, "dynhr_run_record")) rr else NULL,
    theta_mode = mode_result$theta_mode,
    Sigma_prop = mode_result$Sigma_prop,
    logpost    = mode_result$mode$logpost
  ), class = "dynhr_mode_ref")
}

# ---------------------------------------------------------------------------
# Replay
# ---------------------------------------------------------------------------

#' Re-run an estimation from its run record
#'
#' Every result of \code{\link{run_mode_finding}},
#' \code{\link{run_posterior_estimation}} and \code{\link{run_full_estimation}}
#' carries a \code{$run_record} (class \code{dynhr_run_record}) that stores
#' what is needed to replay the call. \code{dynhr_rerun()} replays it: same
#' arguments, same package options, same RNG state. With an unchanged build
#' the replay is bit-identical to the original (draws, mode, log-likelihood).
#' The record survives \code{saveRDS()} / \code{readRDS()}.
#'
#' @section The run record (schema version 2):
#' A named list of class \code{dynhr_run_record}:
#' \describe{
#'   \item{\code{schema_version}}{Integer layout version (\code{2}; version
#'     \code{1} records, which have no \code{spec}, are still read).}
#'   \item{\code{spec}}{The \code{\link{dynhr_estimation_spec}} the run
#'     executed (every entry point builds one): the full description of the
#'     run, with its caches (compiled model, solved object) dropped and a
#'     precomputed mode result replaced by a \code{dynhr_mode_ref}. Write it
#'     with \code{\link{write_spec}}; \code{\link{run_estimation}} runs it.}
#'   \item{\code{created}}{\code{POSIXct} time the run started.}
#'   \item{\code{fn}}{Entry point name.}
#'   \item{\code{args}}{The fully resolved argument list: every formal of
#'     \code{fn} with its default filled in, followed by the \code{...}
#'     arguments (their names in \code{dot_names}). Data and the parsed model
#'     are stored by value. Rebuildable heavy objects are NOT stored (listed in
#'     \code{stripped}): \code{run_mode_finding}'s \code{solved$compiled} and
#'     \code{run_full_estimation}'s \code{compiled} are rebuilt with
#'     \code{\link{compile_model}} (at the recorded \code{max_order}, default
#'     \code{param_deriv}); \code{run_posterior_estimation}'s
#'     \code{mode_result} is replaced by a \code{dynhr_mode_ref}: the mode
#'     run's own run record plus \code{theta_mode}, \code{Sigma_prop} and the
#'     mode log-posterior.}
#'   \item{\code{supplied}}{Names of the arguments the caller supplied.}
#'   \item{\code{dot_names}}{Names of the \code{...} arguments, in order.}
#'   \item{\code{inputs}}{By-value copies of inputs given by path to
#'     \code{run_full_estimation}: \code{model} (the parsed model, when
#'     \code{mod_file} was used) and \code{data} (the data matrix, when
#'     \code{data} was a CSV path).}
#'   \item{\code{stripped}}{Which arguments were not stored as passed.}
#'   \item{\code{rebuild}}{How to rebuild them (e.g. the compiled model's
#'     \code{max_order}).}
#'   \item{\code{options}}{\code{dynhr_get_options(effective = TRUE)} at entry:
#'     every registered package option with its effective value.}
#'   \item{\code{base_options}}{The base-R \code{dynhr.*} options at entry.}
#'   \item{\code{rng}}{\code{kind} (\code{RNGkind()}), \code{seed} (the
#'     \code{seed} argument, \code{NULL} where there is none) and
#'     \code{random_seed} (\code{.Random.seed} at entry, always stored, so an
#'     unseeded run is replayable too).}
#'   \item{\code{provenance}}{dynhr \code{version}, \code{git_commit} stamp,
#'     \code{dev_load}, \code{r_version}, \code{platform}, \code{os},
#'     \code{blas}, \code{lapack}, \code{lapack_version}, \code{parallel},
#'     \code{n_cores} (daemon-pool size cap used; 1 when serial),
#'     \code{mirai_version} (when parallel) and, after a checkpoint resume
#'     that crossed a spec, code or environment mismatch, \code{integrity}
#'     (see \code{\link{run_estimation}}).}
#'   \item{\code{hashes}}{Hash of the \code{data} used and of the parsed
#'     \code{model} (source path blanked), prefixed with its algorithm:
#'     \code{"sha256:"} when the suggested package \pkg{digest} is installed,
#'     else \code{"md5:"} of the serialised object; for
#'     \code{run_full_estimation} also
#'     the md5 of \code{mod_file} and, when \code{data} was a path, of that
#'     file.}
#'   \item{\code{resolved}}{What the run resolved at run time, when it
#'     differs from what a reader could read off \code{args}/\code{options}:
#'     \code{grad_method}, the analytic-gradient method each sampler stage
#'     that built one actually used (a character vector named by sampler, e.g.
#'     \code{c(NUTS = "adjoint_solution")} for the default \code{"auto"} on a
#'     Gaussian model; see \code{\link{make_posterior_grad}}), and
#'     \code{grad_method_requested}, the method each asked for. Absent when
#'     no stage built an analytic gradient.}
#' }
#'
#' @section Replay rules:
#' \itemize{
#'   \item Spec-style replay: a record of \code{\link{run_estimation}}, or
#'     overrides given as edits of a spec component (\code{likelihood},
#'     \code{mode}, \code{sampler}, \code{compute} or \code{outputs} as a named
#'     list, e.g. \code{sampler = list(n_draws = 500L)}; also \code{model},
#'     \code{data}, \code{obs_vars}), replays the record's \code{spec} through
#'     \code{update()} and \code{\link{run_estimation}}; \code{options} are
#'     routed as in \code{update()}. A recorded mode stage is replayed from its
#'     own record first. The result carries a \code{run_estimation} record.
#'     Such edits cannot be mixed with argument overrides.
#'   \item Arguments: the recorded \code{args}, with \code{...} overriding
#'     individual entries. An argument the original caller did not supply is
#'     omitted when its recorded value equals the current default, so
#'     \code{missing()} inside the entry point behaves as it did originally;
#'     if the default has since changed, the recorded value is passed.
#'   \item Options: the recorded snapshot (with \code{options} overriding
#'     entries) replaces the package option store, and the recorded
#'     \code{dynhr.*} base options are set, for the duration of the call
#'     only; both are restored on exit, also when the run errors. Parallel
#'     daemons receive the snapshot through the usual option shipping.
#'   \item Option precedence along a replay chain (a posterior record replays
#'     its mode run first): for EACH stage, an entry of \code{options} wins;
#'     otherwise that stage's OWN recorded snapshot applies -- the mode stage
#'     uses the mode record's snapshot, the sampling stage the posterior
#'     record's. So \code{options} reaches every stage (including options the
#'     mode stage bakes into the log-posterior, such as
#'     \code{power_posterior}), and a replay without \code{options} reproduces
#'     the original chain exactly even when its two stages ran under
#'     different options. When they did, in a result-changing option not
#'     covered by \code{options}, a classed warning
#'     (\code{dynhr_warning_rerun_snapshot_mismatch}) names the options.
#'   \item RNG: \code{.Random.seed} is set to the recorded entry state. The
#'     global RNG is left where the replayed run leaves it, exactly as the
#'     original call left it.
#'   \item \code{run_posterior_estimation}: unless \code{mode_result} is given
#'     in \code{...}, the mode run is replayed first from its own record.
#'     Arguments in \code{...} that belong to \code{run_mode_finding} but not
#'     to \code{run_posterior_estimation} (\code{data}, \code{me_variance},
#'     \code{likelihood}, \code{method}, ...) act through the mode stage and
#'     are routed to that replay; \code{mode_args} passes mode-stage
#'     arguments explicitly and is the only way to reach a name both
#'     functions share (\code{transform_params}, \code{parallel},
#'     \code{n_cores}, \code{verbose}), which in \code{...} goes to the
#'     sampling stage. Mode-stage overrides cannot be combined with a
#'     supplied \code{mode_result}. When nothing was overridden, a classed
#'     warning (\code{dynhr_warning_rerun_mode_drift}) fires if the replayed
#'     \code{theta_mode} differs from the recorded one. A mode result that had
#'     no run record cannot be rebuilt: then \code{mode_result} must be
#'     supplied (error class \code{dynhr_error_rerun_needs_mode_result}).
#'   \item Build and environment (see Reproducibility): when the loaded dynhr
#'     version differs from the record's and a registered result change
#'     between the two touches a component the run used, a classed warning
#'     (\code{dynhr_warning_rerun_result_change}, also of class
#'     \code{dynhr_warning_rerun_version_mismatch}) names the changes and the
#'     version to install for a bit-identical replay. Any other difference
#'     of version, \code{GIT_COMMIT} stamp or numerical environment (R,
#'     platform, OS, BLAS/LAPACK) is a message
#'     (\code{dynhr_message_rerun_build_differs}) stating what to expect.
#' }
#'
#' @section Reproducibility:
#' A record is compared with the current session in three classes:
#' \strong{A} the target and algorithm (model, data, priors, likelihood,
#' sampler and mode settings, result-changing options: the spec's content
#' hashes); \strong{B} the code (dynhr version and \code{GIT_COMMIT}, looked up
#' in the result-change registry, a table of every change that alters results
#' together with the components it touches); \strong{C} the numerical
#' environment (R, platform, BLAS/LAPACK, cores). What to expect:
#' \itemize{
#'   \item Same A, B and C: the replay is bit-identical.
#'   \item Same A; B differs without a registered change touching the
#'     components the run uses, and/or C differs: deterministic quantities
#'     (the mode, the log-posterior at a point) agree to floating-point
#'     precision, about 1e-8 relative; MCMC trajectories do not stay
#'     bit-identical (last-bit differences grow along a chain) but are
#'     equivalent in distribution, with posterior means within Monte Carlo
#'     standard error. The daemon count alone does not change seeded draws.
#'   \item B differs through a registered change touching a used component:
#'     results are expected to differ, and the warning says why.
#'   \item A differs: it is a different run.
#' }
#' \code{\link{dynhr_verify}} checks a finished result against these
#' expectations under the current build and environment.
#'
#' @param x An estimation result carrying \code{$run_record}, or a
#'   \code{dynhr_run_record} (for example read back with \code{readRDS()}).
#' @param ... Named arguments overriding the recorded ones (e.g.
#'   \code{n_draws = 500L}), or spec edits (e.g.
#'   \code{sampler = list(n_draws = 500L)}; see Replay rules).
#' @param options Optional named list overriding entries of the recorded
#'   option snapshot for this replay (e.g.
#'   \code{list(power_posterior = 0.5)}), at every stage of the replay
#'   chain; \code{NULL} entries unset the option.
#' @param mode_args Optional named list of argument overrides for the nested
#'   mode-finding replay of a \code{run_posterior_estimation} record (error
#'   for other records). Entries win over mode-stage names routed from
#'   \code{...}.
#' @return The result of re-invoking the recorded entry point; it carries its
#'   own \code{$run_record}.
#' @seealso \code{\link{dynhr_verify}}, \code{\link{run_estimation}},
#'   \code{\link{as_estimation_spec}}, \code{\link{dynhr_get_options}},
#'   \code{\link{dynhr_set_options}}
#' @examples
#' \donttest{
#' solved   <- solve_model(system.file("extdata/models/nk_demo.mod",
#'                                     package = "dynhr"), verbose = FALSE)
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' mode <- run_mode_finding(solved, Y, obs_vars = obs_vars,
#'                          n_iter = 200L, verbose = FALSE)
#' mode$run_record
#' f <- tempfile(fileext = ".rds")
#' saveRDS(mode$run_record, f)
#' again <- dynhr_rerun(readRDS(f))
#' identical(again$theta_mode, mode$theta_mode)
#' }
#' @export
dynhr_rerun <- function(x, ..., options = NULL, mode_args = NULL) {
  rec <- if (inherits(x, "dynhr_run_record")) x else if (is.list(x)) x$run_record
  if (!inherits(rec, "dynhr_run_record"))
    .dynhr_abort("dynhr_rerun: `x` is neither a dynhr_run_record nor a ",
                 "result carrying one in `$run_record`.",
                 class = "dynhr_error_no_run_record")
  overrides <- list(...)
  .dynhr_rr_check_named(overrides, "arguments in `...`")
  if (!is.null(options)) .dynhr_rr_check_named(options, "`options`")
  if (!is.null(mode_args)) .dynhr_rr_check_named(mode_args, "`mode_args`")
  .dynhr_rr_rerun(rec, overrides, options, mode_args)
}

## `x` must be a list whose elements are all named.
.dynhr_rr_check_named <- function(x, what) {
  if (!is.list(x) || (length(x) &&
        (is.null(names(x)) || any(!nzchar(names(x))))))
    .dynhr_abort("dynhr_rerun: ", what, " must be a named list.",
                 class = "dynhr_error_bad_argument")
  invisible(NULL)
}

## The replay proper. `overrides`: argument overrides for THIS record's entry
## point; `options`: option overrides applied to EVERY stage of the chain;
## `mode_args`: argument overrides for the nested mode replay of a posterior
## record.
.dynhr_rr_rerun <- function(rec, overrides, options, mode_args) {
  sv <- rec$schema_version
  if (!is.numeric(sv) || length(sv) != 1L || sv > .dynhr_run_record_schema)
    .dynhr_abort("dynhr_rerun: run record schema version ", format(sv),
                 " is not supported by this dynhr (reads up to ",
                 .dynhr_run_record_schema, ").",
                 class = "dynhr_error_run_record_schema")
  fn <- rec$fn
  if (!is.character(fn) || length(fn) != 1L || !fn %in% .dynhr_rr_entry_points)
    .dynhr_abort("dynhr_rerun: the record names an unknown entry point (",
                 paste(format(fn), collapse = " "), ").",
                 class = "dynhr_error_run_record_schema")

  ## Spec-style replay (schema 2): a run_estimation() record, or overrides
  ## given as sub-spec edits (`sampler = list(n_draws = 500L)`, ...).
  if (identical(fn, "run_estimation") || .dynhr_rr_spec_style(overrides)) {
    if (length(mode_args))
      .dynhr_abort("dynhr_rerun: `mode_args` does not combine with spec-style ",
                   "overrides; edit the mode stage with mode = list(...).",
                   class = "dynhr_error_bad_argument")
    return(.dynhr_rr_rerun_spec(rec, overrides, options))
  }

  ## Route overrides for a posterior record: an argument of run_mode_finding()
  ## that run_posterior_estimation() does not have (data, me_variance,
  ## likelihood, ...) can only act through the mode stage, so it goes there;
  ## `mode_args` does so explicitly (and is the only way to reach a name both
  ## functions share, e.g. transform_params). Anything else stays here.
  mode_over <- list()
  if (identical(fn, "run_posterior_estimation")) {
    post_f  <- names(formals(run_posterior_estimation))
    mode_f  <- setdiff(names(formals(run_mode_finding)), c(post_f, "..."))
    to_mode <- names(overrides) %in% mode_f
    mode_over <- overrides[to_mode]
    overrides <- overrides[!to_mode]
    for (nm in names(mode_args)) mode_over[nm] <- list(mode_args[[nm]])
    if (length(mode_over) && "mode_result" %in% names(overrides))
      .dynhr_abort("dynhr_rerun: mode-stage overrides (",
                   paste(names(mode_over), collapse = ", "), ") cannot be ",
                   "combined with a supplied `mode_result`.",
                   class = "dynhr_error_bad_argument")
  } else if (length(mode_args)) {
    .dynhr_abort("dynhr_rerun: `mode_args` applies only to a ",
                 "run_posterior_estimation() record (this one is ", fn, ").",
                 class = "dynhr_error_bad_argument")
  }

  .dynhr_rr_check_build(rec)

  ## Arguments first: a posterior replay may replay its mode run, which
  ## applies (and restores) its OWN option snapshot -- plus `options` -- and
  ## its own RNG state.
  args <- .dynhr_rr_replay_args(rec, overrides, options, mode_over)

  ## ---- option snapshot, for this call only ----------------------------
  snap <- .dynhr_rr_stage_options(rec, options)
  saved_opts <- as.list(.dynhr_opts)
  saved_base <- base::options(rec$base_options)
  on.exit({
    rm(list = ls(.dynhr_opts, all.names = TRUE), envir = .dynhr_opts)
    if (length(saved_opts)) list2env(saved_opts, envir = .dynhr_opts)
    if (length(saved_base)) base::options(saved_base)
  }, add = TRUE)
  rm(list = ls(.dynhr_opts, all.names = TRUE), envir = .dynhr_opts)
  if (length(snap)) list2env(snap, envir = .dynhr_opts)

  ## ---- RNG state at entry ---------------------------------------------
  if (!is.null(rec$rng$random_seed)) {
    assign(".Random.seed", rec$rng$random_seed, envir = globalenv())
  } else if (length(rec$rng$kind) == 3L) {
    RNGkind(rec$rng$kind[[1L]], rec$rng$kind[[2L]], rec$rng$kind[[3L]])
  }

  .dynhr_rr_invoke(fn, args)
}

## Spec components that dynhr_rerun() accepts as update()-style overrides.
.dynhr_rr_spec_components <- c("likelihood", "mode", "sampler", "compute",
                               "outputs")

## Are `overrides` spec edits? A list value (or sampler = FALSE) for a spec
## component name. (run_full_estimation()'s `sampler` / `likelihood` formals
## take a string, so the two readings never collide.)
.dynhr_rr_spec_style <- function(overrides) {
  nm <- intersect(names(overrides), .dynhr_rr_spec_components)
  any(vapply(nm, function(k) {
    v <- overrides[[k]]
    is.list(v) || (identical(k, "sampler") && isFALSE(v))
  }, logical(1)))
}

## The runnable spec of a record: its stored spec (schema 2), else the one
## as_estimation_spec() rebuilds from the arguments (schema 1). A recorded
## mode stage (dynhr_mode_ref) is replayed from its own record first, under
## its own snapshot and `options`.
.dynhr_rr_record_spec <- function(rec, options) {
  spec <- rec$spec %||% rec$args$spec
  if (is.null(spec)) spec <- as_estimation_spec(rec)
  ref <- spec$mode$result
  if (inherits(ref, "dynhr_mode_ref") && isTRUE(ref$in_args))
    ref <- rec$args$mode_result
  if (inherits(ref, "dynhr_mode_ref")) {
    if (!inherits(ref$run_record, "dynhr_run_record"))
      .dynhr_abort("dynhr_rerun: the recorded mode_result had no run record, ",
                   "so it cannot be rebuilt. Supply it: ",
                   "dynhr_rerun(x, mode = list(result = <run_mode_finding() result>)).",
                   class = "dynhr_error_rerun_needs_mode_result")
    mr <- .dynhr_rr_rerun(ref$run_record, list(), options, NULL)
    md <- spec$mode
    md["result"] <- list(mr)
    spec["mode"] <- list(md)
  }
  spec
}

## Spec-style replay: the record's spec, edited with update() (`overrides`
## per component, `options` routed to the typed fields / the snapshot), run
## by run_estimation() under the record's option snapshot and entry RNG state.
.dynhr_rr_rerun_spec <- function(rec, overrides, options) {
  allowed <- setdiff(names(formals(update.dynhr_estimation_spec)),
                     c("object", "options", "..."))
  bad <- setdiff(names(overrides), allowed)
  if (length(bad))
    .dynhr_abort("dynhr_rerun: spec-style overrides (",
                 paste(allowed, collapse = ", "), ") cannot be combined with ",
                 "argument overrides (", paste(bad, collapse = ", "), ").",
                 class = "dynhr_error_bad_argument")
  .dynhr_rr_check_build(rec)
  mode_given <- !is.null(overrides$mode$result)
  spec <- if (mode_given) {
    s <- rec$spec %||% rec$args$spec
    s %||% as_estimation_spec(rec)
  } else .dynhr_rr_record_spec(rec, options)
  if (inherits(spec$mode$result, "dynhr_mode_ref")) {
    md <- spec$mode
    md["result"] <- list(NULL)
    spec["mode"] <- list(md)
  }
  spec <- do.call(update, c(list(spec), overrides, list(options = options)))

  snap <- .dynhr_rr_stage_options(rec, options)
  saved_opts <- as.list(.dynhr_opts)
  saved_base <- base::options(rec$base_options)
  on.exit({
    rm(list = ls(.dynhr_opts, all.names = TRUE), envir = .dynhr_opts)
    if (length(saved_opts)) list2env(saved_opts, envir = .dynhr_opts)
    if (length(saved_base)) base::options(saved_base)
  }, add = TRUE)
  rm(list = ls(.dynhr_opts, all.names = TRUE), envir = .dynhr_opts)
  if (length(snap)) list2env(snap, envir = .dynhr_opts)
  if (!is.null(rec$rng$random_seed))
    assign(".Random.seed", rec$rng$random_seed, envir = globalenv())
  run_estimation(spec)
}

## The option store for one stage of a replay: the stage's OWN recorded
## snapshot, with the caller's `options` on top (NULL entries unset).
.dynhr_rr_stage_options <- function(rec, options) {
  snap <- rec$options
  attr(snap, "set") <- NULL
  for (nm in names(options)) snap[nm] <- list(options[[nm]])
  snap[!vapply(snap, is.null, logical(1))]
}

## Result-changing registered options on which two snapshots disagree,
## ignoring the ones the caller overrides for the whole chain.
.dynhr_rr_snapshot_diff <- function(a, b, overridden) {
  reg <- .dynhr_option_registry
  rc  <- names(reg)[vapply(reg, function(e) isTRUE(e$changes_results),
                           logical(1))]
  rc  <- setdiff(rc, overridden)
  rc[!vapply(rc, function(nm) identical(a[[nm]], b[[nm]]), logical(1))]
}

## Numerical-environment fields of the provenance (plan section 7, class C).
.dynhr_rr_env_fields <- c("r_version", "platform", "os", "blas", "lapack",
                          "lapack_version")

## Differences between a record's provenance `was` and the current one `now`:
## list(code = , environment = ) of "field: record -> loaded" lines.
.dynhr_rr_build_diffs <- function(was, now) {
  fmt <- function(v) if (is.null(v) || (length(v) == 1L && is.na(v))) "none"
                     else paste(format(v), collapse = " ")
  code <- character(0)
  if (!identical(as.character(was$version), as.character(now$version)))
    code <- c(code, sprintf("version %s (record) -> %s (loaded)",
                            fmt(was$version), fmt(now$version)))
  if (!identical(as.character(was$git_commit), as.character(now$git_commit)))
    code <- c(code, sprintf("GIT_COMMIT %s (record) -> %s (loaded)",
                            fmt(was$git_commit), fmt(now$git_commit)))
  env_d <- .dynhr_rr_env_fields[!vapply(.dynhr_rr_env_fields, function(k)
    identical(was[[k]], now[[k]]), logical(1))]
  list(code = code,
       environment = if (length(env_d))
         sprintf("%s %s -> %s", env_d,
                 vapply(env_d, function(k) fmt(was[[k]]), ""),
                 vapply(env_d, function(k) fmt(now[[k]]), ""))
       else character(0))
}

## The component tags a record's run used (all stages), from its spec; NULL
## (unknown: every registered change counts) when no spec can be rebuilt --
## a schema-1 posterior record whose mode result had no run record.
.dynhr_rr_record_tags <- function(rec) {
  spec <- rec$spec %||% rec$args$spec
  if (is.null(spec)) {
    if (identical(rec$fn, "run_posterior_estimation") &&
        !inherits(rec$args$mode_result$run_record, "dynhr_run_record"))
      return(NULL)
    spec <- suppressWarnings(suppressMessages(as_estimation_spec(rec)))
  }
  .est_component_tags(spec)
}

## Plan section 7, dynhr_rerun() row: when this build differs from the one
## that made the record, a registered result change touching a component the
## run used WARNS (naming the changes and the version to install for a
## bit-identical replay); any other code or numerical-environment difference
## is a MESSAGE stating the expectation (deterministic quantities to
## floating-point precision, draws equivalent in distribution).
.dynhr_rr_check_build <- function(rec) {
  was <- rec$provenance
  now <- .dynhr_rr_provenance()
  d <- .dynhr_rr_build_diffs(was, now)
  if (!length(d$code) && !length(d$environment)) return(invisible(NULL))
  reg <- if (!identical(as.character(was$version), as.character(now$version)))
    .dynhr_result_changes_between(was$version, now$version,
                                  .dynhr_rr_record_tags(rec))
  else character(0)
  diffs <- paste(c(d$code, d$environment), collapse = "; ")
  if (length(reg)) {
    .dynhr_warn("dynhr_rerun: this run record was made by a different dynhr ",
                "build (", diffs, "), and ", length(reg), " registered result ",
                "change(s) between the two versions touch components this run ",
                "uses, so the replay is EXPECTED to differ: ",
                paste(reg, collapse = "; "), ". Install dynhr ",
                format(was$version), " to reproduce the recorded run ",
                "bit-for-bit.",
                class = c("dynhr_warning_rerun_result_change",
                          "dynhr_warning_rerun_version_mismatch"))
  } else {
    .dynhr_inform("dynhr_rerun: the replay runs under a different ",
                  if (length(d$code)) "build" else "numerical environment",
                  " (", diffs, "). No registered result change touches this ",
                  "run: deterministic quantities (the mode, the log-posterior ",
                  "at a point) agree to floating-point precision (about 1e-8 ",
                  "relative), and MCMC draws are equivalent in distribution ",
                  "(means within MCSE), not bit-identical. dynhr_verify() ",
                  "checks both.",
                  class = "dynhr_message_rerun_build_differs")
  }
  invisible(NULL)
}

## The argument list to pass on replay.
.dynhr_rr_replay_args <- function(rec, overrides, options = NULL,
                                  mode_over = list()) {
  f    <- get(rec$fn, envir = asNamespace("dynhr"), inherits = FALSE)
  args <- rec$args
  anm  <- names(args)
  if (is.null(anm)) anm <- rep("", length(args))

  ## by-value inputs (run_full_estimation given paths)
  if (!is.null(rec$inputs$model) && !"model" %in% names(overrides))
    args["model"] <- list(rec$inputs$model)
  if (!is.null(rec$inputs$data) && !"data" %in% names(overrides))
    args["data"] <- list(rec$inputs$data)

  ## rebuild stripped objects
  if ("solved$compiled" %in% rec$stripped && !"solved" %in% names(overrides)) {
    solved <- args$solved
    solved$compiled <- compile_model(
      solved$model, verbose = FALSE,
      max_order = rec$rebuild$solved_compiled$max_order %||% 1L)
    args["solved"] <- list(solved)
  }
  if ("mode_result" %in% rec$stripped && !"mode_result" %in% names(overrides)) {
    ref <- args$mode_result
    if (!inherits(ref$run_record, "dynhr_run_record"))
      .dynhr_abort("dynhr_rerun: the recorded mode_result had no run record, ",
                   "so it cannot be rebuilt. Supply it: ",
                   "dynhr_rerun(x, mode_result = <run_mode_finding() result>).",
                   class = "dynhr_error_rerun_needs_mode_result")
    mrec <- ref$run_record
    ## Each stage replays under its OWN recorded snapshot -- that is what the
    ## original chain ran under (the mode stage bakes e.g. power_posterior
    ## into log_post_fn) -- and `options` overrides BOTH. Say so when the
    ## two stages ran under different result-changing options.
    mism <- .dynhr_rr_snapshot_diff(rec$options, mrec$options, names(options))
    if (length(mism))
      .dynhr_warn("dynhr_rerun: the mode run and the posterior run were made ",
                  "under different result-changing options (",
                  paste(vapply(mism, function(nm) sprintf(
                    "%s: mode %s, posterior %s", nm,
                    paste(format(mrec$options[[nm]]), collapse = " "),
                    paste(format(rec$options[[nm]]), collapse = " ")),
                    character(1)), collapse = "; "),
                  "). Each stage is replayed under its own recorded options; ",
                  "pass options = list(...) to force one value on both.",
                  class = "dynhr_warning_rerun_snapshot_mismatch")
    mr <- .dynhr_rr_rerun(mrec, mode_over, options, NULL)
    ## drift is only news when nothing was changed on purpose
    if (!length(mode_over) && !length(options) &&
        !identical(mr$theta_mode, ref$theta_mode))
      .dynhr_warn("dynhr_rerun: the replayed mode run reached a different ",
                  "theta_mode from the recorded one; the posterior replay ",
                  "starts from the NEW mode.",
                  class = "dynhr_warning_rerun_mode_drift")
    args["mode_result"] <- list(mr)
  }

  ## preserve missing(): drop an unsupplied formal still at its default
  fnames <- setdiff(names(formals(f)), "...")
  keep <- rep(TRUE, length(args))
  for (i in seq_along(args)) {
    nm <- anm[[i]]
    if (!nzchar(nm) || !nm %in% fnames || nm %in% rec$supplied ||
        nm %in% names(overrides))
      next
    d <- .dynhr_rr_current_default(f, nm)
    if (isTRUE(d$ok) && identical(d$value, args[[i]])) keep[i] <- FALSE
  }
  args <- args[keep]

  for (nm in names(overrides)) args[nm] <- list(overrides[[nm]])
  args
}

## Call entry point `fn` with `args`, passing each value by SYMBOL from a
## private frame so the call (sys.call(), tracebacks, error messages) shows
## `run_mode_finding(solved = solved, data = data, ...)`, not deparsed data.
.dynhr_rr_invoke <- function(fn, args) {
  e   <- new.env(parent = asNamespace("dynhr"))
  anm <- names(args)
  if (is.null(anm)) anm <- rep("", length(args))
  syms <- vector("list", length(args))
  for (i in seq_along(args)) {
    v <- if (nzchar(anm[[i]])) anm[[i]] else paste0(".rr_dot_", i)
    assign(v, args[[i]], envir = e)
    syms[[i]] <- as.name(v)
  }
  names(syms) <- anm
  do.call(fn, syms, envir = e)
}

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

#' @rdname dynhr_rerun
#' @param width Line width for \code{format()} (unused beyond truncation of
#'   long lists).
#' @export
format.dynhr_run_record <- function(x, width = 72L, ...) {
  pv  <- x$provenance
  na_or <- function(v, alt = "none")
    if (length(v) != 1L || is.na(v)) alt else as.character(v)
  short <- function(h) if (length(h) != 1L || is.na(h)) "none" else
    sub("^(\\w+):(.{12}).*$", "\\1 \\2", h)
  trunc_list <- function(v) {
    s <- paste(v, collapse = ", ")
    if (nchar(s) > width) paste0(substr(s, 1L, width - 3L), "...") else s
  }
  set_nm <- attr(x$options, "set")
  set_txt <- if (length(set_nm)) {
    trunc_list(vapply(set_nm, function(nm) {
      v <- x$options[[nm]]
      paste0(nm, "=", if (is.atomic(v) && length(v) == 1L) format(v) else
        paste0("<", class(v)[1L], ">"))
    }, character(1)))
  } else "none (all defaults)"
  mode_line <- NULL
  if (inherits(x$args$mode_result, "dynhr_mode_ref")) {
    mrr <- x$args$mode_result$run_record
    mode_line <- paste0("  mode_result : ",
                        if (inherits(mrr, "dynhr_run_record"))
                          paste0("replayable (", mrr$fn, " record of ",
                                 format(mrr$created, "%Y-%m-%d %H:%M:%S"), ")")
                        else "NOT replayable (no run record; pass mode_result=)")
  }
  grad_line <- .grad_method_record_line(x$resolved)
  c(sprintf("<dynhr_run_record> schema %s", format(x$schema_version)),
    paste0("  fn          : ", x$fn),
    paste0("  created     : ", format(x$created, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("  dynhr       : ", na_or(pv$version),
           if (!is.na(pv$git_commit %||% NA))
             paste0(" (GIT_COMMIT ", substr(pv$git_commit, 1L, 12L), ")"),
           if (isTRUE(pv$dev_load)) " [dev load]"),
    paste0("  R           : ", na_or(pv$r_version), " / ", na_or(pv$platform)),
    paste0("  parallel    : ", if (isTRUE(pv$parallel))
      paste0("yes, ", pv$n_cores, " core(s), mirai ", na_or(pv$mirai_version))
      else "no"),
    paste0("  seed        : ", if (is.null(x$rng$seed)) "none" else
      format(x$rng$seed), "  (entry .Random.seed ",
      if (is.null(x$rng$random_seed)) "not stored)" else "stored)"),
    paste0("  supplied    : ", trunc_list(x$supplied)),
    paste0("  options set : ", set_txt),
    paste0("  data        : ", short(x$hashes$data)),
    paste0("  model       : ", short(x$hashes$model)),
    mode_line,
    if (!is.null(x$spec))
      paste0("  spec        : ", short(x$spec$hashes$spec),
             " (", .run_est_form_label(x$spec), ")"),
    if (!is.null(grad_line)) paste0("  gradient    : ", grad_line),
    if (!is.null(x$provenance$integrity))
      "  integrity   : resumed across a checkpoint mismatch (see $provenance$integrity)",
    if (length(x$stripped))
      paste0("  rebuilt on rerun: ", paste(x$stripped, collapse = ", ")))
}

## "form: <form>" label of a recorded spec.
.run_est_form_label <- function(spec) {
  form <- spec$outputs$form %||% "auto"
  paste0("form ", if (identical(form, "auto")) .est_form(spec) else form)
}

#' @rdname dynhr_rerun
#' @export
print.dynhr_run_record <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname dynhr_rerun
#' @export
as.list.dynhr_run_record <- function(x, ...) unclass(x)
