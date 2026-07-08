## R/dynare-modefile.R
## --------------------------------------------------------------------------
## read_dynare_mode_file() -- P2 gap #7: loader for Dynare posterior-mode
## .mat files (e.g. <mod_name>_mode.mat).
##
## FOOTGUN this fixes: Dynare's posterior-mode shock standard deviations live
## in `M_.Sigma_e` / the shock-std entries of `xparam1` INSIDE the mode .mat
## file, NOT in `dynare_solution.json$params` (the steady-state/calibration
## export dynhr's replication harness normally reads). Grabbing "the params"
## from dynare_solution.json and assuming it includes the estimated mode is a
## silent-wrong-answer trap; this reader targets the actual mode artifact.
## --------------------------------------------------------------------------

#' Read a Dynare posterior-mode \code{.mat} file
#'
#' Dynare's mode-finding step (\code{mode_compute}) writes the posterior mode
#' -- INCLUDING estimated shock standard deviations/variances -- to a
#' \code{.mat} file such as \code{<mod_name>_mode.mat}. Critically, those
#' mode shock stderrs live in \code{xparam1} (the raw estimated-parameter
#' vector, in \code{estimated_params} order) and in \code{M_.Sigma_e} (the
#' assembled shock covariance matrix at the mode) -- they are NOT present in
#' \code{dynare_solution.json$params}, which only carries the calibrated/
#' steady-state parameter values dynhr's replication harness exports
#' separately. Reaching for \code{dynare_solution.json} to recover an
#' estimated mode's shock stderrs silently returns the wrong (calibration,
#' not posterior-mode) values. Use this function instead when you need the
#' actual posterior-mode estimates.
#'
#' \code{R.matlab::readMat()} returns nested-\code{struct} fields as deeply
#' nested lists (a top-level \code{M_} struct arrives named \code{"M."} with
#' its own sub-fields addressed positionally via the \code{"dimnames"}
#' attribute, and Dynare's char-matrix fields such as \code{param_names}
#' arrive as a list of 1-element character lists rather than a plain
#' character matrix). This function extracts the fields defensively rather
#' than assuming a flat structure, and returns NULL for any field that is
#' absent from a given mode file instead of erroring.
#'
#' @param path Character scalar: path to the Dynare mode \code{.mat} file.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{xparam1}}{Numeric vector of estimated parameters at the
#'       mode (in \code{estimated_params} order). Named with
#'       \code{parameter_names} when the lengths match.}
#'     \item{\code{parameter_names}}{Character vector of estimated parameter
#'       names (from a top-level \code{parameter_names} field or
#'       \code{M_.param_names}), or \code{NULL} if not present.}
#'     \item{\code{Sigma_e}}{Numeric matrix, the shock covariance at the mode
#'       (from \code{M_.Sigma_e}), or \code{NULL} if not present.}
#'     \item{\code{exo_names}}{Character vector of shock (exogenous variable)
#'       names (from \code{M_.exo_names}), or \code{NULL} if not present.}
#'     \item{\code{hh}}{Numeric matrix, the inverse Hessian at the mode, if
#'       present (Dynare's \code{hh} / \code{hessian} field); else
#'       \code{NULL}.}
#'     \item{\code{raw}}{The full object returned by \code{R.matlab::readMat()},
#'       for callers that need to dig further.}
#'   }
#'
#' @examples
#' \dontrun{
#' modefile <- read_dynare_mode_file("usmodel_mode.mat")
#' modefile$xparam1        ## posterior-mode parameter vector, named
#' modefile$Sigma_e         ## shock covariance at the mode
#' }
#' @export
read_dynare_mode_file <- function(path) {
  if (!requireNamespace("R.matlab", quietly = TRUE))
    stop(paste0(
      "read_dynare_mode_file() requires the 'R.matlab' package.  ",
      "Install it with: install.packages(\"R.matlab\")."
    ))

  if (!file.exists(path))
    stop(sprintf("read_dynare_mode_file: file not found: %s", path))

  raw <- tryCatch(
    R.matlab::readMat(path),
    error = function(e)
      stop(sprintf("read_dynare_mode_file: failed to read '%s': %s",
                    path, conditionMessage(e)))
  )

  ## ------------------------------------------------------------------
  ## Helpers for readMat's nested-struct representation.
  ## ------------------------------------------------------------------

  ## Look up a field in a (possibly struct-shaped) list by name, trying a
  ## few name variants readMat may have produced (dots preserved, dots
  ## turned into nothing, trailing "." for reserved/underscore names).
  .get_field <- function(container, name) {
    if (is.null(container)) return(NULL)
    nms <- names(container)
    if (!is.null(nms) && name %in% nms) return(container[[name]])

    ## struct fields: names live in the "dimnames" attribute, values are
    ## addressed positionally.
    dn <- attr(container, "dimnames")
    if (is.list(dn) && length(dn) >= 1L && !is.null(dn[[1]])) {
      idx <- match(name, dn[[1]])
      if (!is.na(idx)) return(container[[idx]])
    }
    NULL
  }

  ## A top-level struct variable "M_" is renamed "M." by readMat (dots are
  ## not valid R names and the trailing underscore becomes a dot); try both.
  .get_top <- function(name) {
    v <- .get_field(raw, name)
    if (!is.null(v)) return(v)
    alt <- sub("_$", ".", name)
    if (!identical(alt, name)) {
      v <- .get_field(raw, alt)
      if (!is.null(v)) return(v)
    }
    NULL
  }

  M_ <- .get_top("M_")

  ## Dynare struct field names use underscores (Sigma_e, param_names,
  ## exo_names); MATLAB struct field access through readMat may instead
  ## surface them with dots (Sigma.e, param.names, exo.names) depending on
  ## how the .mat file was produced. Try both spellings.
  .get_struct_field <- function(container, name) {
    v <- .get_field(container, name)
    if (!is.null(v)) return(v)
    alt <- gsub("_", ".", name)
    if (!identical(alt, name)) {
      v <- .get_field(container, alt)
      if (!is.null(v)) return(v)
    }
    NULL
  }

  ## Recursively unwrap readMat's char-matrix-as-list-of-1-element-lists
  ## representation into a plain character vector, trimmed.
  .unwrap_char_field <- function(v) {
    if (is.null(v)) return(NULL)
    if (is.character(v)) return(trimws(as.character(v)))
    if (is.list(v)) {
      out <- vapply(v, function(x) {
        while (is.list(x)) x <- x[[1]]
        trimws(as.character(x))
      }, character(1L))
      return(out)
    }
    NULL
  }

  ## ------------------------------------------------------------------
  ## xparam1
  ## ------------------------------------------------------------------
  xparam1_raw <- .get_top("xparam1")
  xparam1 <- if (!is.null(xparam1_raw)) as.numeric(xparam1_raw) else NULL

  ## ------------------------------------------------------------------
  ## parameter_names: top-level "parameter_names", else M_.param_names
  ## ------------------------------------------------------------------
  pn_raw <- .get_top("parameter_names")
  if (is.null(pn_raw)) pn_raw <- .get_struct_field(M_, "param_names")
  parameter_names <- .unwrap_char_field(pn_raw)

  if (!is.null(xparam1) && !is.null(parameter_names) &&
        length(xparam1) == length(parameter_names)) {
    names(xparam1) <- parameter_names
  }

  ## ------------------------------------------------------------------
  ## Sigma_e: M_.Sigma_e
  ## ------------------------------------------------------------------
  Sigma_e_raw <- .get_struct_field(M_, "Sigma_e")
  Sigma_e <- if (!is.null(Sigma_e_raw)) {
    Sm <- as.matrix(Sigma_e_raw)
    storage.mode(Sm) <- "double"
    Sm
  } else NULL

  ## ------------------------------------------------------------------
  ## exo_names: M_.exo_names
  ## ------------------------------------------------------------------
  exo_raw <- .get_struct_field(M_, "exo_names")
  exo_names <- .unwrap_char_field(exo_raw)

  if (!is.null(Sigma_e) && !is.null(exo_names) &&
        nrow(Sigma_e) == length(exo_names)) {
    dimnames(Sigma_e) <- list(exo_names, exo_names)
  }

  ## ------------------------------------------------------------------
  ## hh: inverse Hessian at the mode (top-level "hh", else "hessian")
  ## ------------------------------------------------------------------
  hh_raw <- .get_top("hh")
  if (is.null(hh_raw)) hh_raw <- .get_top("hessian")
  hh <- if (!is.null(hh_raw)) {
    Hm <- as.matrix(hh_raw)
    storage.mode(Hm) <- "double"
    Hm
  } else NULL

  list(
    xparam1         = xparam1,
    parameter_names = parameter_names,
    Sigma_e         = Sigma_e,
    exo_names       = exo_names,
    hh              = hh,
    raw             = raw
  )
}


#' Polish an imported Dynare mode with dynhr's own mode-finding
#'
#' An imported Dynare posterior mode (e.g. from \code{\link{read_dynare_mode_file}})
#' is a \strong{likelihood-LEVEL validation point}, not a critical point of
#' dynhr's own posterior: presample handling, Kalman-filter initialisation, and
#' other likelihood conventions can differ between Dynare and dynhr even when
#' both implement "the same" model. \strong{Computing curvature (a Hessian /
#' proposal covariance) directly at an imported Dynare mode is therefore
#' silently misleading} -- the gradient of dynhr's posterior need not vanish
#' there, so the curvature measured at that point does not describe dynhr's
#' posterior around its own mode. (This has already cost one paper draft.)
#'
#' This function does the safe thing instead: it validates the imported
#' parameter vector against \code{prior_spec}, then \strong{warm-starts
#' dynhr's own mode-finding} (\code{\link{find_mode}}) from it, so that any
#' downstream curvature computation happens at an actual critical point of
#' dynhr's posterior.
#'
#' @param mode  Either the list returned by \code{\link{read_dynare_mode_file}}
#'   (its \code{xparam1} element is used, and must be named) or a named
#'   numeric vector of parameter values (e.g. a hand-built Dynare
#'   \code{xparam1}).
#' @param model     Parsed model (\code{dynhr_mod} from \code{\link{parse_mod}}).
#' @param data      Observation matrix (\eqn{T \times n_{obs}}), column names
#'   matching \code{obs_vars}.
#' @param prior_spec  Prior specification data.frame (from
#'   \code{extract_prior_spec}); \code{prior_spec$name} is the
#'   authoritative list of dynhr's estimated parameters.
#' @param obs_vars  Character vector of observable variable names.
#' @param compiled  Compiled model (from \code{\link{compile_model}}).
#' @param method    Optimizer sequence forwarded to \code{\link{find_mode}}
#'   (default \code{"newrat"}).
#' @param ...       Additional arguments forwarded to \code{\link{find_mode}}
#'   (currently \code{n_iter}, \code{verbose}; an unrecognised name fails
#'   loud with "unused argument" rather than being silently dropped, since
#'   \code{find_mode} has no \code{...} of its own).
#'
#' @return A list:
#'   \describe{
#'     \item{\code{imported}}{Named numeric vector: the validated imported
#'       mode, restricted/reordered to \code{prior_spec$name}.}
#'     \item{\code{imported_logpost}}{dynhr log-posterior at the imported mode.}
#'     \item{\code{polished}}{The full \code{\link{find_mode}} result
#'       (\code{theta_mode}, \code{logpost}, \code{convergence},
#'       \code{iterations}, \code{method}), warm-started from
#'       \code{imported}.}
#'     \item{\code{polished_logpost}}{dynhr log-posterior at the polished mode
#'       (\code{polished$logpost}, duplicated here for convenience).}
#'     \item{\code{gap_nats}}{\code{polished_logpost - imported_logpost}: how
#'       far the imported point was from dynhr's own mode, in nats. Should be
#'       \code{>= 0} (up to optimizer noise); a large gap is exactly the
#'       likelihood-convention mismatch this function exists to catch.}
#'   }
#'
#' @examples
#' \dontrun{
#' dyn_mode <- read_dynare_mode_file("usmodel_mode.mat")
#' res <- polish_dynare_mode(dyn_mode, model, data, prior_spec, obs_vars, compiled)
#' res$gap_nats            ## how far Dynare's mode was from dynhr's
#' res$polished$theta_mode  ## use THIS for dynhr Hessian / MCMC, not dyn_mode
#' }
#' @seealso \code{\link{read_dynare_mode_file}}, \code{\link{find_mode}},
#'   \code{\link{run_mode_finding}}
#' @export
polish_dynare_mode <- function(mode, model, data, prior_spec, obs_vars,
                               compiled, method = "newrat", ...) {

  ## ------------------------------------------------------------------
  ## 1. Extract a named numeric vector from `mode` (read_dynare_mode_file()
  ##    result or a plain named vector), then validate/map its names against
  ##    prior_spec$name -- fail loud on anything that doesn't match, rather
  ##    than silently subsetting or recycling.
  ## ------------------------------------------------------------------
  xparam1 <- if (is.list(mode) && !is.null(mode$xparam1)) {
    mode$xparam1
  } else if (is.numeric(mode)) {
    mode
  } else {
    stop("polish_dynare_mode: `mode` must be a read_dynare_mode_file() ",
         "result (with a named $xparam1) or a named numeric vector.",
         call. = FALSE)
  }

  if (is.null(names(xparam1)) || any(!nzchar(names(xparam1))))
    stop("polish_dynare_mode: `mode` (or `mode$xparam1`) must be a NAMED ",
         "numeric vector -- cannot map unnamed values onto prior_spec's ",
         "estimated parameters.", call. = FALSE)

  est_names <- prior_spec$name
  unmatched <- setdiff(names(xparam1), est_names)
  if (length(unmatched) > 0L)
    stop(sprintf(
      "polish_dynare_mode: %d name(s) in `mode` do not match any estimated ",
      length(unmatched)),
      "parameter in `prior_spec$name`:\n",
      paste0("  - ", unmatched, collapse = "\n"),
      "\nAvailable estimated parameters: ", paste(est_names, collapse = ", "),
      call. = FALSE)

  missing_names <- setdiff(est_names, names(xparam1))
  if (length(missing_names) > 0L)
    stop(sprintf(
      "polish_dynare_mode: `mode` is missing %d estimated parameter(s) ",
      length(missing_names)),
      "present in `prior_spec$name`:\n",
      paste0("  - ", missing_names, collapse = "\n"),
      call. = FALSE)

  ## Reorder to prior_spec's canonical order.
  imported <- setNames(as.numeric(xparam1[est_names]), est_names)

  ## ------------------------------------------------------------------
  ## 2. Build dynhr's own log-posterior and evaluate it at the imported mode
  ##    -- this IS the likelihood-level validation point (NOT a critical
  ##    point of dynhr's posterior; see caveat above).
  ## ------------------------------------------------------------------
  log_post_fn <- make_log_posterior(model, data, prior_spec, obs_vars, compiled)
  imported_res <- log_post_fn(imported)
  imported_logpost <- imported_res$logpost

  ## ------------------------------------------------------------------
  ## 3. Warm-start dynhr's own mode-finding FROM the imported mode. This is
  ##    the actual fix: any curvature computed downstream should be at a
  ##    genuine critical point of dynhr's posterior, not at the (possibly
  ##    off-critical-point) imported Dynare mode.
  ## ------------------------------------------------------------------
  polished <- find_mode(log_post_fn, imported, prior_spec,
                        method = method, ...)

  if (is.null(polished) || !is.finite(polished$logpost))
    stop("polish_dynare_mode: dynhr mode-finding failed (non-finite ",
         "log-posterior) starting from the imported mode.", call. = FALSE)

  list(
    imported         = imported,
    imported_logpost = imported_logpost,
    polished         = polished,
    polished_logpost = polished$logpost,
    gap_nats         = polished$logpost - imported_logpost
  )
}
