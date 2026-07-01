## R/param-set.R
## --------------------------------------------------------------------------
## set_param_value() / set_param_values() — Dynare-equivalent API for
## overriding a model's calibrated parameter vector.
##
## inject_params() — lenient superset-injection wrapper (HA-2 workflow):
##   accepts named list or vector, warns on unknown names, reports NA-fills.
##
## read_mat_params() — thin .mat reader returning a named numeric vector
##   suitable for inject_params(); gated on the R.matlab Suggests package.
##
## Use-case (M17 / L11 / HA-2): parameters computed inside an external Dynare
## *_steadystate.m (e.g. via fsolve) or MATLAB steadyState.m are invisible to
## dynhr's parser.  These helpers let the caller inject those values before
## solving:
##
##   m <- parse_mod("my_model.mod")
##   m <- inject_params(m, c(chi = 0.5341, I_ss = 0.1, c1 = 0.06, A = 1))
##   solved <- solve_model(m)
##
## HA-2 two-stage workflow (Winberry-class HANK models):
##   m    <- parse_mod("my_model.mod")   # 949 params are NA after parse
##   pars <- read_mat_params("params.mat")      # load .mat params (R.matlab)
##   m    <- inject_params(m, pars)            # fills NAs; warns on extras
##   ss   <- solve_steady(compile_model(m), m$param_values, y0 = external_ss)
##
## The returned model object has its param_values updated; all downstream
## functions (solve_model, solve_steady, solve_perturbation, compute_irfs,
## compute_moments, make_log_posterior, …) read params from model$param_values
## by default, so this is the single mutation point.
## --------------------------------------------------------------------------


#' Override one parameter value in a parsed model
#'
#' Returns a copy of \code{model} with \code{model$param_values[name]} set to
#' \code{value}.  This is the R equivalent of Dynare's \code{set_param_value}
#' and is the recommended way to supply parameters that are computed inside an
#' external \code{*_steadystate.m} (e.g. via \code{fsolve}) and are therefore
#' invisible to dynhr's parser.
#'
#' @param model A \code{dynhr_mod} object returned by \code{\link{parse_mod}}.
#' @param name  Character scalar: the parameter name to set.  Must already
#'   appear in \code{model$param_names}.
#' @param value Numeric scalar: the new value.
#' @return A \code{dynhr_mod} with the updated \code{param_values}.
#'
#' @seealso \code{\link{set_param_values}} for setting multiple parameters at
#'   once; \code{\link{solve_model}} for the downstream solver.
#'
#' @examples
#' \dontrun{
#' m <- parse_mod("my_model.mod")
#' m <- set_param_value(m, "chi", 0.5341)
#' solved <- solve_model(m)
#' }
#' @export
set_param_value <- function(model, name, value) {
  if (!inherits(model, "dynhr_mod"))
    stop("'model' must be a dynhr_mod (from parse_mod()).")
  if (!is.character(name) || length(name) != 1L)
    stop("'name' must be a single character string.")
  if (!is.numeric(value) || length(value) != 1L)
    stop("'value' must be a single numeric scalar.")
  if (!(name %in% model$param_names))
    stop(sprintf("Parameter '%s' not found in model$param_names.", name))

  model$param_values[name] <- value
  model
}


#' Override multiple parameter values in a parsed model
#'
#' Returns a copy of \code{model} with \code{model$param_values} updated from
#' \code{params}.  This is the R equivalent of Dynare's \code{set_param_value}
#' applied to a whole vector and is the recommended way to supply parameters
#' computed inside an external \code{*_steadystate.m} that are invisible to
#' dynhr's parser.
#'
#' Workflow for models with an external steady-state file:
#' \enumerate{
#'   \item Call \code{\link{parse_mod}} to obtain the base model.
#'   \item Compute the externally-calibrated parameters (or read them from
#'         a Dynare export / \code{.mat} file).
#'   \item Pass them to \code{set_param_values}.
#'   \item Feed the updated model to \code{\link{solve_model}}.
#' }
#'
#' @param model  A \code{dynhr_mod} object returned by \code{\link{parse_mod}}.
#' @param params Named numeric vector of parameter overrides.  Every name must
#'   already appear in \code{model$param_names}.
#' @return A \code{dynhr_mod} with the updated \code{param_values}.
#'
#' @seealso \code{\link{set_param_value}} for a single parameter;
#'   \code{\link{solve_model}} for the downstream solver.
#'
#' @examples
#' \dontrun{
#' ## Gertler-Karadi (2011): 4 params computed by the external *_steadystate.m
#' m <- parse_mod("my_model.mod")
#' m <- set_param_values(m, c(chi = 0.5341, I_ss = 0.1817,
#'                             c1  = 0.0618, A    = 0.5040))
#' solved <- solve_model(m, order = 1)
#' }
#' @export
set_param_values <- function(model, params) {
  if (!inherits(model, "dynhr_mod"))
    stop("'model' must be a dynhr_mod (from parse_mod()).")
  if (!is.numeric(params) || is.null(names(params)))
    stop("'params' must be a named numeric vector.")
  unknown <- setdiff(names(params), model$param_names)
  if (length(unknown) > 0L)
    stop(sprintf(
      "Parameter(s) not found in model$param_names: %s",
      paste(unknown, collapse = ", ")
    ))

  model$param_values[names(params)] <- params
  model
}


#' Inject a (super)set of parameter values into a parsed model
#'
#' A lenient convenience wrapper over \code{\link{set_param_values}} designed
#' for the HA-2 two-stage workflow where an external steady-state solver
#' (MATLAB \code{steadyState.m}) produces a larger set of parameters than the
#' model may currently need.
#'
#' Differences from \code{set_param_values}:
#' \itemize{
#'   \item Accepts a named \strong{list} as well as a named numeric vector
#'         (list elements are coerced to numeric scalars via \code{as.numeric}).
#'   \item \strong{Warns} (does not error) on names that do not appear in
#'         \code{model$param_names}, and silently ignores them.  This allows
#'         passing the full output of a MATLAB parameter file even when only a
#'         subset maps to declared Dynare parameters.
#'   \item \strong{Reports} (via \code{message}) how many previously-\code{NA}
#'         entries in \code{model$param_values} were filled by this call.
#' }
#'
#' @param model  A \code{dynhr_mod} object returned by \code{\link{parse_mod}}.
#' @param params Named numeric vector \strong{or} named list of parameter values.
#'   Names not in \code{model$param_names} trigger a warning and are dropped.
#' @param .quiet Logical; if \code{TRUE} suppress the NA-fill message.
#'   Default \code{FALSE}.
#' @return A \code{dynhr_mod} with \code{param_values} updated for all
#'   names that exist in both \code{params} and \code{model$param_names}.
#'
#' @seealso \code{\link{set_param_values}} for the strict (error-on-unknown)
#'   variant; \code{\link{read_mat_params}} for loading a \code{.mat} file into
#'   a named vector; \code{\link{solve_steady}} for the \code{y0} SS-seed arg.
#'
#' @examples
#' \dontrun{
#' ## Winberry (2018) HA model: parse -> inject external params -> steady(y0)
#' m    <- parse_mod("my_model.mod")        # 949 params are NA
#' pars <- read_mat_params("params.mat")           # load from .mat
#' m    <- inject_params(m, pars)                 # fills 25 grid params; warns on extras
#' comp <- compile_model(m)
#' ss   <- solve_steady(comp, m$param_values, y0 = external_ss_vector)
#' dr   <- solve_perturbation(comp, ss, order = 1)
#' }
#' @export
inject_params <- function(model, params, .quiet = FALSE) {
  if (!inherits(model, "dynhr_mod"))
    stop("'model' must be a dynhr_mod (from parse_mod()).")

  # Coerce list to named numeric vector
  if (is.list(params)) {
    if (is.null(names(params)))
      stop("'params' list must be named.")
    params <- vapply(params, function(x) {
      v <- suppressWarnings(as.numeric(x))
      if (length(v) != 1L)
        stop("Each element of 'params' list must coerce to a single numeric scalar.")
      v
    }, numeric(1L))
  }

  if (!is.numeric(params) || is.null(names(params)))
    stop("'params' must be a named numeric vector or named list.")

  # Warn on unknown names but don't error — drop them
  unknown <- setdiff(names(params), model$param_names)
  if (length(unknown) > 0L) {
    warning(sprintf(
      "inject_params: %d name(s) not in model$param_names (ignored): %s%s",
      length(unknown),
      paste(head(unknown, 5L), collapse = ", "),
      if (length(unknown) > 5L) sprintf(", ... (%d total)", length(unknown)) else ""
    ), call. = FALSE)
    params <- params[names(params) %in% model$param_names]
  }

  if (length(params) == 0L) {
    if (!.quiet) message("inject_params: no matching parameters found; model unchanged.")
    return(model)
  }

  # Count how many were previously NA
  was_na <- sum(is.na(model$param_values[names(params)]))

  model$param_values[names(params)] <- params

  if (!.quiet && was_na > 0L)
    message(sprintf(
      "inject_params: filled %d previously-NA parameter(s) out of %d injected.",
      was_na, length(params)
    ))

  model
}


#' Read a MATLAB \code{.mat} file and return a named numeric vector of params
#'
#' Convenience reader for the HA-2 two-stage workflow: Winberry-class models
#' write their Chebyshev/quadrature parameters to \code{.mat} files from a
#' MATLAB \code{steadyState.m} script.  This function loads the file via
#' \code{R.matlab::readMat()} (listed in \code{Suggests}) and flattens all
#' numeric scalar fields into a single named numeric vector suitable for
#' \code{\link{inject_params}}.
#'
#' @section Installation:
#' \code{R.matlab} must be installed separately:
#' \preformatted{install.packages("R.matlab")}
#'
#' @section Alternative (no R.matlab):
#' If \code{R.matlab} is not available, export the \code{.mat} parameters to
#' JSON from Octave/MATLAB:
#' \preformatted{
#' ## In Octave:
#' load("params.mat");
#' fields = fieldnames(vars);
#' out = struct();
#' for i = 1:numel(fields)
#'   v = vars.(fields{i});
#'   if isscalar(v), out.(fields{i}) = v; end
#' end
#' pkg load io; jsonlite_dump(out, "grids.json");
#' ## Then in R:
#' pars <- unlist(jsonlite::fromJSON("grids.json"))
#' }
#'
#' @param path   Character scalar: path to the \code{.mat} file.
#' @param scalar_only Logical (default \code{TRUE}): if \code{TRUE}, only
#'   fields that are \eqn{1 \times 1} numeric matrices (i.e. scalar values)
#'   are included.  Set to \code{FALSE} to include all numeric fields
#'   (non-scalars are dropped with a message).
#' @return Named numeric vector with one entry per scalar field in the
#'   \code{.mat} file.
#'
#' @seealso \code{\link{inject_params}} to apply the result to a model.
#'
#' @examples
#' \dontrun{
#' pars <- read_mat_params("params.mat")
#' m    <- inject_params(m, pars)
#' }
#' @export
read_mat_params <- function(path, scalar_only = TRUE) {
  if (!requireNamespace("R.matlab", quietly = TRUE))
    stop(paste0(
      "read_mat_params() requires the 'R.matlab' package.  ",
      "Install it with: install.packages(\"R.matlab\").\n",
      "Alternatively, export .mat fields to JSON from Octave and use ",
      "jsonlite::fromJSON()."
    ))

  if (!file.exists(path))
    stop(sprintf("File not found: %s", path))

  raw <- R.matlab::readMat(path)

  out <- numeric(0L)
  for (nm in names(raw)) {
    v <- raw[[nm]]
    # readMat wraps scalars in 1x1 matrices
    if (is.numeric(v) && length(v) == 1L) {
      out <- c(out, setNames(as.numeric(v), nm))
    } else if (!scalar_only && is.numeric(v)) {
      # Flatten vector/matrix fields with indexed names
      flat <- as.numeric(v)
      idx_names <- if (length(flat) == 1L) nm else
        paste0(nm, "_", seq_along(flat))
      out <- c(out, setNames(flat, idx_names))
    }
  }

  if (length(out) == 0L)
    warning(sprintf("read_mat_params: no numeric scalar fields found in '%s'.", path),
            call. = FALSE)

  out
}
