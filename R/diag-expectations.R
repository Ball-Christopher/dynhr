## R/diag-expectations.R
## --------------------------------------------------------------------------
## diag_expectations() -- check @dynhr:expectations blocks from .mod files.
##
## Replaces model-specific Stage 4 logic in pipeline scripts. Each check
## spec is parsed from the .mod file and evaluated against observed data,
## parameter vectors, or IRF results.
##
## Supported check types:
##   data_ratio  -- mean(numerator) / mean(denominator) in [min, max]
##   data_mean   -- mean(variable) in [min, max]
##   data_sd     -- sd(variable) in [min, max]
##   param_range -- params[variable] in [min, max]
##   irf_sign    -- sign(irf[[shock]][horizon, variable]) matches expected
## --------------------------------------------------------------------------


#' Check model expectations against data and estimated parameters
#'
#' Reads the \code{@dynhr:expectations} block from \code{model$metadata} (or
#' directly from a parsed model) and evaluates each check against the
#' supplied inputs.  Returns a \code{dynhr_diagnostic} object.
#'
#' @details
#' Expectations are declared in the \code{.mod} file as a
#' \code{@dynhr:expectations} block, e.g.
#' \preformatted{
#' // @dynhr:expectations
#' // cy_ratio:  type="data_ratio",  numerator="c_obs", denominator="y_obs",
#' //            min=0.5, max=2.0, description="C/Y plausibility"
#' // rho_a_hi:  type="param_range", variable="rho_a",
#' //            min=0.5, max=0.999, description="Tech AR persistence"
#' // @dynhr:end
#' }
#'
#' @param model   dynhr_mod (from \code{parse_mod}) or a list with a
#'   \code{$metadata$expectations} slot.
#' @param data    Optional numeric matrix (\eqn{T \times n_{obs}}) with named
#'   columns.  Required for \code{data_ratio}, \code{data_mean},
#'   \code{data_sd} checks.
#' @param params  Optional named numeric vector of parameter values (estimated
#'   or calibrated).  Required for \code{param_range} checks.
#' @param irfs    Optional named list of IRF matrices (shock -> T x n_var
#'   matrix).  Required for \code{irf_sign} checks.
#' @param meta    Optional provenance descriptor from \code{diag_meta()}.
#' @return A \code{dynhr_diagnostic} object.
#' @export
diag_expectations <- function(model, data = NULL, params = NULL,
                               irfs = NULL, meta = NULL) {

  # --- Extract expectations spec ---
  specs <- NULL
  if (!is.null(model$metadata$expectations)) {
    specs <- model$metadata$expectations
  } else if (!is.null(model$expectations)) {
    specs <- model$expectations
  } else if (!is.null(model$source_file) && file.exists(model$source_file)) {
    # Try extracting metadata directly from the source .mod file
    mod_meta <- extract_mod_metadata(model$source_file)
    if (!is.null(mod_meta$expectations) && length(mod_meta$expectations) > 0) {
      specs <- mod_meta$expectations
    }
  }

  if (is.null(specs) || length(specs) == 0) {
    return(.make_result(
      pass    = NA,
      summary = "diag_expectations: no @dynhr:expectations block found in model"
    ))
  }

  # --- Evaluate each check ---
  check_results <- lapply(specs, function(s) .eval_expectation(s, data, params, irfs))

  n_total  <- length(check_results)
  n_pass   <- sum(vapply(check_results, function(r) isTRUE(r$pass), logical(1)))
  n_fail   <- sum(vapply(check_results, function(r) identical(r$pass, FALSE), logical(1)))
  n_skip   <- sum(vapply(check_results, function(r) is.na(r$pass), logical(1)))

  all_pass <- if (n_skip == n_total) NA else (n_fail == 0)

  lines <- vapply(check_results, function(r) {
    badge <- if (is.na(r$pass)) "SKIP" else if (r$pass) "PASS" else "FAIL"
    sprintf("  [%s] %s: %s", badge, r$name, r$detail)
  }, character(1))

  summary_txt <- paste0(
    sprintf("diag_expectations: %d/%d checks passed", n_pass, n_total - n_skip),
    if (n_skip > 0) sprintf(", %d skipped (inputs not provided)", n_skip) else "",
    "\n",
    paste(lines, collapse = "\n")
  )

  .make_result(
    pass    = all_pass,
    summary = summary_txt,
    result  = list(checks = check_results, n_pass = n_pass, n_fail = n_fail,
                   n_skip = n_skip),
    plots   = list()
  )
}


# ---------------------------------------------------------------------------
# Internal: evaluate a single expectation spec
# ---------------------------------------------------------------------------
.eval_expectation <- function(spec, data, params, irfs) {
  name <- spec$name %||% "unnamed"
  desc <- spec$description %||% name

  result_skip <- function(reason)
    list(name = name, pass = NA, detail = paste(desc, "--", reason))
  result_pass <- function(detail)
    list(name = name, pass = TRUE, detail = detail)
  result_fail <- function(detail)
    list(name = name, pass = FALSE, detail = detail)

  type <- spec$type %||% "unknown"

  switch(type,

      data_ratio = {
        if (is.null(data)) return(result_skip("data not supplied"))
        num_col <- spec$numerator
        den_col <- spec$denominator
        if (is.null(num_col) || is.null(den_col))
          return(result_skip("numerator/denominator not specified"))
        if (!(num_col %in% colnames(data)))
          return(result_skip(sprintf("column '%s' not in data", num_col)))
        if (!(den_col %in% colnames(data)))
          return(result_skip(sprintf("column '%s' not in data", den_col)))
        num_mean <- mean(data[, num_col], na.rm = TRUE)
        den_mean <- mean(data[, den_col], na.rm = TRUE)
        if (abs(den_mean) < .Machine$double.eps)
          return(result_skip(sprintf("denominator '%s' mean is ~0", den_col)))
        ratio <- num_mean / den_mean
        .check_range(ratio, spec, result_pass, result_fail,
                     sprintf("%s/%s = %.3f", num_col, den_col, ratio))
      },

      data_mean = {
        if (is.null(data)) return(result_skip("data not supplied"))
        col <- spec$variable
        if (is.null(col)) return(result_skip("variable not specified"))
        if (!(col %in% colnames(data)))
          return(result_skip(sprintf("column '%s' not in data", col)))
        val <- mean(data[, col], na.rm = TRUE)
        .check_range(val, spec, result_pass, result_fail,
                     sprintf("mean(%s) = %.4f", col, val))
      },

      data_sd = {
        if (is.null(data)) return(result_skip("data not supplied"))
        col <- spec$variable
        if (is.null(col)) return(result_skip("variable not specified"))
        if (!(col %in% colnames(data)))
          return(result_skip(sprintf("column '%s' not in data", col)))
        val <- sd(data[, col], na.rm = TRUE)
        .check_range(val, spec, result_pass, result_fail,
                     sprintf("sd(%s) = %.4f", col, val))
      },

      param_range = {
        if (is.null(params)) return(result_skip("params not supplied"))
        pname <- spec$variable
        if (is.null(pname)) return(result_skip("variable not specified"))
        if (!(pname %in% names(params)))
          return(result_skip(sprintf("parameter '%s' not found", pname)))
        val <- params[[pname]]
        .check_range(val, spec, result_pass, result_fail,
                     sprintf("%s = %.4f", pname, val))
      },

      irf_sign = {
        if (is.null(irfs)) return(result_skip("irfs not supplied"))
        shock    <- spec$shock
        variable <- spec$variable
        horizon  <- spec$horizon %||% 1L
        sign_exp <- spec$sign %||% "positive"
        if (is.null(shock) || is.null(variable))
          return(result_skip("shock/variable not specified"))
        if (!(shock %in% names(irfs)))
          return(result_skip(sprintf("shock '%s' not in irfs", shock)))
        irf_mat <- irfs[[shock]]
        if (is.null(colnames(irf_mat)) || !(variable %in% colnames(irf_mat)))
          return(result_skip(sprintf("variable '%s' not in irf[['%s']]", variable, shock)))
        h <- min(horizon, nrow(irf_mat))
        val <- irf_mat[h, variable]
        expected_positive <- tolower(sign_exp) %in% c("positive", "+", "pos")
        actual_positive   <- val > 0
        detail <- sprintf("irf[%s][%s][h=%d] = %.4f, expected %s",
                          shock, variable, h, val, sign_exp)
        if (expected_positive == actual_positive) result_pass(detail)
        else result_fail(detail)
      },

      result_skip(sprintf("unknown check type '%s'", type))
    )
}


.check_range <- function(val, spec, result_pass, result_fail, label) {
  lo   <- spec$min %||% -Inf
  hi   <- spec$max %||%  Inf
  ok   <- (is.na(lo) || val >= lo) && (is.na(hi) || val <= hi)
  detail <- sprintf("%s  [expected: %.3g to %.3g]", label, lo, hi)
  if (ok) result_pass(detail) else result_fail(detail)
}
