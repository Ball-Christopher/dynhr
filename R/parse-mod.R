## R/parse-mod.R
## --------------------------------------------------------------------------
## Top-level parse_mod() orchestrator, the dynhr_mod S3 constructor, S3
## print/summary methods, the global .KNOWN_FUNCTIONS recognised in
## equations, and build_variable_classification (which uses ast_collect_*
## from parse-equations.R).
##
## Phase-1 split from parser-monolith.R (no logic changes).
## --------------------------------------------------------------------------

################################################################################
# dynhr_parser.R  --  v0.3
# Phase 1 -- Pure-R parser for Dynare .mod files
#
# Part of the dynhr project: a native R implementation of core Dynare
# functionality (stoch_simul, estimation) without external dependencies
# on MATLAB, Octave, or the Dynare preprocessor.
#
# Author:  dynhr project
# Date:    2026-05-04
# License: MIT
#
# Changelog v0.3:
#   - Fix 1: extract_declaration now finds ALL matching declaration blocks
#            (not just the first) and handles optional options e.g. var(log)
#   - Fix 2: parse_calibration collapses multiline whitespace before eval(),
#            fixing crashes on expressions spanning multiple lines
#   - Fix 3: parse_initval_block applies the same multiline whitespace fix
#   - Fix 4: remove_blocks decl_kw pattern updated for optional options
#            consistency with extract_declaration
#
# Changelog v0.2:
#   - Fix 1: % comments now stripped anywhere on a line (not just at start)
#   - Fix 2: extract_declaration uses (?si) for multi-line declarations
#   - Fix 3: new_dynhr_mod uses top-level overwrite, not modifyList
################################################################################

# Known math function names recognised in model equations
.KNOWN_FUNCTIONS <- c(
  "exp", "log", "ln", "log2", "log10",
  "sqrt", "cbrt", "abs", "sign",
  "sin", "cos", "tan", "asin", "acos", "atan",
  "sinh", "cosh", "tanh",
  "max", "min",
  "normcdf", "normpdf", "erf", "erfc",
  "diff", "adl",
  "STEADY_STATE", "steady_state", "EXPECTATION"
)


#' Build the lead_lag_incidence matrix and classify variables
#'
#' @param equations  List of equation objects (each with $lhs and $rhs ASTs).
#' @param var_names  Character vector of endogenous variable names.
#' @param varexo_names Character vector of exogenous variable names.
#' @return A list with lead_lag_incidence, classification, counts, etc.
#' @noRd
build_variable_classification <- function(equations, var_names,
                                          varexo_names = character(0),
                                          predetermined_vars = character(0),
                                          local_vars = list()) {
  all_var_names <- c(var_names, varexo_names)

  all_refs <- do.call(rbind, lapply(equations, function(eq) {
    rbind(ast_collect_variables(eq$lhs, local_vars),
          ast_collect_variables(eq$rhs, local_vars))
  }))

  # `do.call(rbind, list())` returns NULL when there are no equations (or no
  # equation contributes any variable reference), and subsetting NULL yields
  # NULL whose nrow() is also NULL -- so the `nrow(endo_refs) == 0` test below
  # would error "argument is of length zero" rather than taking the intended
  # empty branch.  Normalise to a zero-row data frame so the empty case is
  # detected and classified (all declared vars -> static) instead of crashing
  # (repl H7).  This is the legitimate degenerate case: no endogenous variable
  # appears in any equation timing, e.g. an empty/exo-only model.
  if (is.null(all_refs) || nrow(all_refs) == 0) {
    all_refs <- data.frame(name = character(0), lead_lag = integer(0),
                           stringsAsFactors = FALSE)
  }

  endo_refs <- all_refs[all_refs$name %in% var_names, , drop = FALSE]

  if (nrow(endo_refs) == 0) {
    n_endo <- length(var_names)
    return(list(
      lead_lag_incidence = matrix(0L, nrow = 3, ncol = n_endo,
                                  dimnames = list(c("t-1","t","t+1"),
                                                  var_names)),
      max_lag  = 1L,
      max_lead = 1L,
      classification = list(static = var_names, predetermined = character(0),
                            forward = character(0), mixed = character(0)),
      n_static = n_endo, n_predetermined = 0L,
      n_forward = 0L, n_mixed = 0L
    ))
  }

  max_lag  <- max(1L, -min(endo_refs$lead_lag, 0L))
  max_lead <- max(1L, max(endo_refs$lead_lag, 0L))

  row_labels <- (-max_lag):max_lead
  n_rows <- length(row_labels)
  n_endo <- length(var_names)

  lli <- matrix(0L, nrow = n_rows, ncol = n_endo)
  rownames(lli) <- paste0("t", ifelse(row_labels >= 0,
                                      paste0("+", row_labels),
                                      as.character(row_labels)))
  rownames(lli)[row_labels == 0] <- "t"
  colnames(lli) <- var_names

  col_counter <- 1L

  for (i in seq_along(var_names)) {
    for (j in seq_along(row_labels)) {
      ll <- row_labels[j]
      if (any(endo_refs$name == var_names[i] &
              endo_refs$lead_lag == ll)) {
        lli[j, i] <- col_counter
        col_counter <- col_counter + 1L
      }
    }
  }

  t_minus_row <- which(row_labels == -1L)
  t_row       <- which(row_labels == 0L)
  t_plus_row  <- which(row_labels == 1L)

  # NOTE on predetermined variables: the -1 timing shift is applied to the
  # model-block text in parse_mod (step 4) BEFORE this classification runs, so by
  # the time we get here the LLI already reflects the shifted timing (a
  # predetermined stock k appears at t-1, not t).  No LLI-level adjustment is
  # needed or wanted here; doing one would double-shift.  `predetermined_vars`
  # is retained in the signature for callers/diagnostics but intentionally
  # unused for the incidence matrix.

  static_vars  <- character(0)
  pred_vars    <- character(0)
  fwd_vars     <- character(0)
  mixed_vars   <- character(0)

  for (i in seq_along(var_names)) {
    has_lag  <- length(t_minus_row) > 0 && lli[t_minus_row, i] > 0
    has_lead <- length(t_plus_row) > 0  && lli[t_plus_row, i] > 0
    has_cur  <- length(t_row) > 0       && lli[t_row, i] > 0

    if (has_lag && has_lead && has_cur) {
      # Appears at all three time periods - genuinely mixed (both backward
      # and forward looking, e.g. a jump variable appearing with all timings).
      mixed_vars <- c(mixed_vars, var_names[i])
    } else if (has_lag && has_lead) {
      # Appears at t-1 and t+1 but NOT at t. Dynare classifies these as
      # "mixed" for the forward-looking count because they have a lead.
      # They contribute to both the backward and forward variable counts
      # for the BK condition.
      mixed_vars <- c(mixed_vars, var_names[i])
    } else if (has_lag) {
      pred_vars <- c(pred_vars, var_names[i])
    } else if (has_lead) {
      fwd_vars <- c(fwd_vars, var_names[i])
    } else {
      static_vars <- c(static_vars, var_names[i])
    }
  }

  list(
    lead_lag_incidence = lli,
    max_lag            = max_lag,
    max_lead           = max_lead,
    classification     = list(
      static        = static_vars,
      predetermined = pred_vars,
      forward       = fwd_vars,
      mixed         = mixed_vars
    ),
    n_static        = length(static_vars),
    n_predetermined = length(pred_vars),
    n_forward       = length(fwd_vars),
    n_mixed         = length(mixed_vars)
  )
}


#' Create a dynhr_mod object
#'
#' @param ... Named fields to populate the model.
#' @return An object of class "dynhr_mod".
#' @noRd
new_dynhr_mod <- function(...) {
  fields <- list(...)
  defaults <- list(
    source_file         = NA_character_,
    var_names           = character(0),
    varexo_names        = character(0),
    varexo_det_names    = character(0),
    param_names         = character(0),
    predetermined_vars  = character(0),
    param_values        = numeric(0),
    equations           = list(),
    local_variables     = list(),
    model_options       = list(),
    initval             = numeric(0),
    endval              = numeric(0),
    ## histval: named list  var -> numeric vector indexed by LAG (element k is
    ## the value k periods before the first simulation period).  See
    ## parse_histval_block() for the Dynare index convention.
    histval             = list(),
    steady_state_model  = NULL,
    shocks              = list(variances = data.frame(), correlations = data.frame()),
    ## shock_groups: named list  group label -> character vector of shocks,
    ## taken from the FIRST shock_groups block (the fallback grouping used by
    ## historical_decomposition()).  shock_groups_blocks keeps every block,
    ## keyed by its `name=` option.  See parse_shock_groups().
    shock_groups        = list(),
    shock_groups_blocks = list(),
    det_shocks          = data.frame(name = character(0), period = integer(0),
                                     value = numeric(0), stringsAsFactors = FALSE),
    filter_tunes        = list(tunes = data.frame()),
    heteroskedastic_shocks = list(scales = data.frame()),
    stochastic_volatility = list(sv = data.frame()),
    estimated_params    = data.frame(),
    estimated_params_init = data.frame(),
    commands            = list(),
    occbin_constraints  = list(),
    lead_lag_incidence  = matrix(0),
    variable_classification = list(),
    n_static            = 0L,
    n_predetermined     = 0L,
    n_forward           = 0L,
    n_mixed             = 0L
    ,planner_objective  = list(text = "", ast = NULL)
    ,varobs             = character(0)
    ,varobs_names       = character(0)
    ,obs_vars           = character(0)
    ,ramsey_instruments = character(0)
    ,metadata           = list()
    ,mcp_constraints    = NULL
  )

  # Overwrite defaults with user-supplied fields (top-level only, no
  # recursive merge into nested structures like data.frames)
  for (nm in names(fields)) {
    defaults[[nm]] <- fields[[nm]]
  }
  # exo_names is an alias for varexo_names so user code like
  # m$exo_names or length(m$exo_names) returns the shock list, not NULL.
  defaults$exo_names <- defaults$varexo_names
  model <- defaults
  class(model) <- "dynhr_mod"
  model
}


#' Print a parsed dynhr model
#'
#' Prints a one-screen summary of a `dynhr_mod` object: the source file,
#' variable / shock / parameter / equation counts, the variable
#' classification (static, predetermined, forward-looking, mixed), and any
#' shock, filter-tune, estimated-parameter or command blocks that were
#' parsed.
#'
#' @param x A `dynhr_mod` object, as returned by [parse_mod()].
#' @param ... Ignored; present for S3 generic compatibility.
#'
#' @return `x`, invisibly. Called for the side effect of printing.
#'
#' @export
print.dynhr_mod <- function(x, ...) {
  cat("=== dynhr_mod ===\n")
  if (!is.na(x$source_file))
    cat("Source file:        ", x$source_file, "\n")
  cat("Endogenous vars:    ", length(x$var_names), "\n")
  cat("Exogenous shocks:   ", length(x$varexo_names), "\n")
  cat("Parameters:         ", length(x$param_names), "\n")
  cat("Equations:          ", length(x$equations), "\n")
  cat("Local (#) variables:", length(x$local_variables), "\n")
  cat("\nVariable classification:\n")
  cat("  Static:         ", x$n_static, "\n")
  cat("  Predetermined:  ", x$n_predetermined, "\n")
  cat("  Forward-looking:", x$n_forward, "\n")
  cat("  Mixed:          ", x$n_mixed, "\n")

  if (nrow(x$shocks$variances) > 0)
    cat("\nShocks defined:    ", nrow(x$shocks$variances), "\n")
  if (!is.null(x$filter_tunes$tunes) && nrow(x$filter_tunes$tunes) > 0)
    cat("Filter tunes:      ", nrow(x$filter_tunes$tunes), "\n")
  if (nrow(x$estimated_params) > 0)
    cat("Estimated params:  ", nrow(x$estimated_params), "\n")

  if (length(x$commands) > 0) {
    cmd_names <- vapply(x$commands, function(c) c$name, character(1))
    cat("Commands:          ", paste(cmd_names, collapse = ", "), "\n")
  }

  if (!is.null(x$model_options) && length(x$model_options) > 0) {
    cat("Model options:     ",
        paste(names(x$model_options), collapse = ", "), "\n")
  }

  invisible(x)
}


#' Summarise a parsed dynhr model
#'
#' Prints everything [print.dynhr_mod()] shows, then the detail behind it:
#' the variables in each timing class, calibrated parameter values, the
#' shock variance / correlation tables, the estimated-parameter block, and
#' the model equations rendered back to source form.
#'
#' @param object A `dynhr_mod` object, as returned by [parse_mod()].
#' @param ... Ignored; present for S3 generic compatibility.
#'
#' @return `object`, invisibly. Called for the side effect of printing.
#'
#' @export
summary.dynhr_mod <- function(object, ...) {
  print(object)
  cat("\n--- Variable Details ---\n")
  vc <- object$variable_classification
  if (length(vc$static) > 0)
    cat("  Static:       ", paste(vc$static, collapse = ", "), "\n")
  if (length(vc$predetermined) > 0)
    cat("  Predetermined:", paste(vc$predetermined, collapse = ", "), "\n")
  if (length(vc$forward) > 0)
    cat("  Forward:      ", paste(vc$forward, collapse = ", "), "\n")
  if (length(vc$mixed) > 0)
    cat("  Mixed:        ", paste(vc$mixed, collapse = ", "), "\n")

  if (length(object$param_values) > 0) {
    cat("\n--- Calibrated Parameters ---\n")
    for (nm in names(object$param_values)) {
      cat("  ", nm, "=", object$param_values[nm], "\n")
    }
  }

  if (nrow(object$shocks$variances) > 0) {
    cat("\n--- Shock Structure ---\n")
    print(object$shocks$variances, row.names = FALSE)
    if (nrow(object$shocks$correlations) > 0) {
      cat("  Correlations:\n")
      print(object$shocks$correlations, row.names = FALSE)
    }
  }

  if (nrow(object$estimated_params) > 0) {
    cat("\n--- Estimated Parameters ---\n")
    print(object$estimated_params, row.names = FALSE)
  }

  if (length(object$equations) > 0) {
    cat("\n--- Equations (", length(object$equations), ") ---\n")
    for (i in seq_along(object$equations)) {
      eq <- object$equations[[i]]
      tag_str <- if (!is.na(eq$tag)) paste0(" [", eq$tag, "]") else ""
      cat("  (", i, ")", tag_str, " ",
          ast_to_string(eq$lhs), " = ", ast_to_string(eq$rhs), "\n")
    }
  }

  if (!is.null(object$lead_lag_incidence) &&
      nrow(object$lead_lag_incidence) > 1) {
    cat("\n--- Lead-Lag Incidence Matrix ---\n")
    print(object$lead_lag_incidence)
  }

  invisible(object)
}


#' Locate a sibling external steady-state file for a .mod
#'
#' Looks for `<stem>_steadystate.m` / `<stem>_steadystate2.m` next to the
#' .mod. Review/harness copies are often suffixed (e.g. GK_2011_pp.mod) while
#' the external file keeps the original stem (GK_2011_steadystate.m), so the
#' exact stem is tried first, then the stem with trailing review suffixes
#' stripped.
#'
#' @param source_file Path of the .mod, or NA for inline text.
#' @return The first existing file path, or NA_character_.
#' @noRd
.external_ss_file <- function(source_file) {
  if (is.na(source_file)) return(NA_character_)
  ss_dir  <- dirname(source_file)
  ss_stem <- sub("\\.(mod|dyn)$", "", basename(source_file), ignore.case = TRUE)
  ss_stems <- unique(c(
    ss_stem,
    sub("(_pp|_export|_stochsimul|_dynhr)+$", "", ss_stem)
  ))
  ss_files <- as.vector(t(outer(
    ss_stems, c("_steadystate.m", "_steadystate2.m"),
    function(s, suf) file.path(ss_dir, paste0(s, suf))
  )))
  ss_files <- ss_files[file.exists(ss_files)]
  if (length(ss_files) == 0L) NA_character_ else ss_files[[1L]]
}


#' The `first_obs` option of the .mod's `estimation(...)` command, as text
#'
#' Read through the command-option parser, so an earlier option holding a
#' parenthesised list (`estimation(optim=('MaxIter',200), first_obs=5)`) no
#' longer hides it (the old `\\(([^)]*)\\)` capture stopped at its `)`).
#' @return Character scalar, or NULL when absent.
#' @noRd
.estimation_first_obs <- function(txt) {
  ci <- extract_command(txt, "estimation")
  if (!ci$found) return(NULL)
  fo <- parse_command_options(ci$options_str)$first_obs
  if (is.null(fo) || length(fo) != 1L) return(NULL)
  if (is.numeric(fo)) format(fo, scientific = FALSE) else as.character(fo)
}


#' Every `KEYWORD[(options)] names;` declaration statement of .mod text
#'
#' The options are read with BALANCED parentheses (`.dynhr_opts_re`), so
#' `var(deflator = A^(1/(1-alpha))) Y;` is one statement with options
#' `deflator = A^(1/(1-alpha))`.
#'
#' @param txt_decl Comment-stripped text with every paired block removed.
#' @param kw       Keyword: "var" (never matches varexo / varobs / trend_var),
#'   "trend_var" or "log_trend_var".
#' @return List of list(opts = option text ("" when none), names).
#' @noRd
.mod_decl_statements <- function(txt_decl, kw) {
  kw_pat <- if (kw == "var") "var(?!exo)" else kw
  pat <- paste0("(?si)\\b", kw_pat, "\\b\\s*", .dynhr_opts_re,
                "\\s*([^;]*?)\\s*;")
  hits <- regmatches(txt_decl, gregexpr(pat, txt_decl, perl = TRUE))[[1]]
  lapply(hits, function(h) {
    p <- regmatches(h, regexec(pat, h, perl = TRUE))[[1]]
    list(opts = trimws(p[2]), names = parse_declaration_names(p[3]))
  })
}

## `key = value` options of one declaration as a named character vector
## (a bare flag such as `log` maps to "").
.mod_decl_option_map <- function(opts) {
  if (!nzchar(opts)) return(character(0))
  parts <- trimws(.split_options_top_level(opts, ","))
  parts <- parts[nzchar(parts)]
  out <- character(0)
  for (p in parts) {
    eq <- .options_top_level_eq(p)
    if (eq > 0L) out[tolower(trimws(substr(p, 1L, eq - 1L)))] <-
      trimws(substr(p, eq + 1L, nchar(p)))
    else out[tolower(p)] <- ""
  }
  out
}


#' Names declared under `var(log)` (Dynare 6+)
#'
#' Scans every `var(<options>) names;` declaration whose options contain the
#' `log` flag (options read with balanced parentheses).
#'
#' @param txt_decl Comment-stripped text with every paired block removed.
#' @return Character vector of names, in declaration order.
#' @noRd
.mod_log_var_names <- function(txt_decl) {
  out <- character(0)
  for (st in .mod_decl_statements(txt_decl, "var")) {
    om <- .mod_decl_option_map(st$opts)
    if ("log" %in% names(om) && !nzchar(om[["log"]]))
      out <- c(out, st$names)
  }
  unique(out)
}


## ===========================================================================
## Nonstationary models: trend_var / log_trend_var / var(deflator=) /
## var(log_deflator=)  (Dynare reference manual, "Stationarizing variables")
## ===========================================================================
##
## `trend_var(growth_factor = G) A;` declares a trend variable with
## A_t = G_t * A_{t-1}; `log_trend_var(log_growth_factor = g) LA;` an additive
## one, LA_t = g_t + LA_{t-1}.  `var(deflator = D) Y;` declares Y as trending
## with the multiplicative trend D (an expression in trend variables,
## endogenous variables and parameters); `var(log_deflator = D) Y;` with the
## additive trend D.  The model block is written in the NONSTATIONARY
## variables; the preprocessor rewrites it into the stationary one, in which
## every declared name stands for its DETRENDED variable -- which is what the
## steady state, steady_state_model, initval and the decision rule refer to.
##
## The rewrite (Dynare 7.1 preprocessor, ModFile::transformPass ->
## DynamicModel::detrendEquations + removeTrendVariableFromEquations; checked
## against its `json=transform` output):
##   1. each deflated Y, LAST declared first:  Y(k) -> Y(k) * D(k)
##      (log_deflator: Y(k) + D(k)), D(k) being D with every variable's
##      timing moved by k (VariableNode::detrend);
##   2. each trend variable at a lead or lag (VariableNode::removeTrendLeadLag):
##        A(+k) -> A * G(+1) * ... * G(+k),   A(-k) -> A / (G * G(-1) * ... * G(-k+1))
##      (additive: sums, `+` for leads and `-` for lags);
##   3. every remaining A (timing 0) -> 1 (additive trend: 0)
##      (BinaryOpNode::replaceTrendVar).
## Step 3 is valid only on a balanced growth path; see
## .mod_check_balanced_growth().

#' Trend-variable and deflator declarations of .mod text
#'
#' @param txt_decl Comment-stripped text with every paired block removed.
#' @return NULL when the .mod declares no trend variable and no deflator;
#'   otherwise list(trend_vars = data.frame(name, growth_factor, log),
#'   deflated = data.frame(name, deflator, log)), each in declaration order.
#' @noRd
.mod_nonstationary_decls <- function(txt_decl) {
  bad <- function(...) .dynhr_abort("parse_mod: ", ...,
                                    class = "dynhr_error_mod_syntax")
  tv <- data.frame(name = character(0), growth_factor = character(0),
                   log = logical(0), stringsAsFactors = FALSE)
  for (kw in c("trend_var", "log_trend_var")) {
    key <- if (kw == "trend_var") "growth_factor" else "log_growth_factor"
    for (st in .mod_decl_statements(txt_decl, kw)) {
      om <- .mod_decl_option_map(st$opts)
      if (!identical(names(om), key) || !nzchar(om[[key]]))
        bad(kw, " takes exactly one option, `", key, " = EXPRESSION`; got `",
            kw, if (nzchar(st$opts)) paste0("(", st$opts, ")") else "", "`.")
      if (length(st$names) == 0L)
        bad(kw, "(", st$opts, ") declares no variable.")
      tv <- rbind(tv, data.frame(name = st$names, growth_factor = om[[key]],
                                 log = kw == "log_trend_var",
                                 stringsAsFactors = FALSE))
    }
  }
  df <- data.frame(name = character(0), deflator = character(0),
                   log = logical(0), stringsAsFactors = FALSE)
  for (st in .mod_decl_statements(txt_decl, "var")) {
    om <- .mod_decl_option_map(st$opts)
    has <- intersect(c("deflator", "log_deflator"), names(om))
    if (length(has) == 0L) next
    if (length(has) == 2L)
      bad("var(", st$opts, ") gives both `deflator` and `log_deflator`.")
    if (!nzchar(om[[has]]))
      bad("var(", st$opts, "): `", has, "` needs an expression.")
    df <- rbind(df, data.frame(name = st$names, deflator = om[[has]],
                               log = has == "log_deflator",
                               stringsAsFactors = FALSE))
  }
  if (nrow(tv) == 0L && nrow(df) == 0L) return(NULL)
  dup <- unique(c(tv$name[duplicated(tv$name)], df$name[duplicated(df$name)]))
  if (length(dup) > 0L)
    bad("variable(s) declared twice as trend variable or deflated variable: ",
        paste(dup, collapse = ", "), ".")
  rownames(tv) <- NULL
  rownames(df) <- NULL
  list(trend_vars = tv, deflated = df)
}

#' Validate trend / deflator declarations against the declared symbols
#'
#' @param ns   Result of .mod_nonstationary_decls().
#' @param endo,exo,params Declared endogenous / exogenous names, parameters.
#' @return `ns`, invisibly; aborts (dynhr_error_mod_syntax) on any problem.
#' @noRd
.mod_check_nonstationary <- function(ns, endo, exo, params) {
  bad <- function(...) .dynhr_abort("parse_mod: ", ...,
                                    class = "dynhr_error_mod_syntax")
  tv <- ns$trend_vars
  clash <- intersect(tv$name, c(endo, exo, params))
  if (length(clash) > 0L)
    bad("trend variable(s) ", paste(clash, collapse = ", "),
        " also declared as a variable or parameter.")
  ids_of <- function(e) setdiff(.osr_expr_ids(e), c(.KNOWN_FUNCTIONS, "pi"))
  for (k in seq_len(nrow(tv))) {
    ids <- ids_of(tv$growth_factor[k])
    if (length(intersect(ids, tv$name)) > 0L)
      bad("the growth factor of trend variable ", tv$name[k], " (`",
          tv$growth_factor[k], "`) uses a trend variable; trends of trends ",
          "(I(2) trends) are not supported.")
    unk <- setdiff(ids, c(endo, params))
    if (length(unk) > 0L)
      bad("the growth factor of trend variable ", tv$name[k], " (`",
          tv$growth_factor[k], "`) uses ", paste(unk, collapse = ", "),
          ", which ", if (length(unk) == 1L) "is" else "are",
          " not a declared endogenous variable or parameter.")
  }
  df <- ns$deflated
  for (k in seq_len(nrow(df))) {
    what <- if (df$log[k]) "log_deflator" else "deflator"
    if (!(df$name[k] %in% endo))
      bad("var(", what, "=...) variable ", df$name[k],
          " is not an endogenous variable.")
    unk <- setdiff(ids_of(df$deflator[k]), c(tv$name, endo, params))
    if (length(unk) > 0L)
      bad("the ", what, " of ", df$name[k], " (`", df$deflator[k], "`) uses ",
          paste(unk, collapse = ", "), ", which ",
          if (length(unk) == 1L) "is" else "are",
          " not a declared trend variable, endogenous variable or parameter.")
  }
  invisible(ns)
}

## Replace every occurrence of the NAMED variables in model-block text.
##
## `repl` is a named list of functions t -> replacement text, t being the
## occurrence's timing (0 for a bare name).  Tokenised as .mod_retime_text():
## other identifiers, equation tags, quoted strings, number literals and the
## whole argument of STEADY_STATE(...) are left untouched.
.mod_substitute_var_text <- function(text, repl) {
  if (length(repl) == 0L || !nzchar(text)) return(text)
  m <- gregexpr(.MOD_TIMING_TOKEN_RE, text, perl = TRUE)[[1]]
  if (m[1L] == -1L) return(text)
  starts <- as.integer(m)
  ends   <- starts + attr(m, "match.length") - 1L
  toks   <- substring(text, starts, ends)
  ss_fns <- c("STEADY_STATE", "steady_state")
  hits   <- which(toks %in% c(names(repl), ss_fns))
  if (length(hits) == 0L) return(text)

  n <- nchar(text)
  pieces  <- character(0)
  cursor  <- 1L
  skip_to <- 0L
  for (j in hits) {
    if (starts[j] <= skip_to) next
    after <- ends[j] + 1L
    look  <- substr(text, after, min(n, after + 200L))
    if (toks[j] %in% ss_fns) {
      op <- regexpr("^\\s*\\(", look, perl = TRUE)
      if (op > 0L) {
        close <- .mod_matching_paren(text, after + attr(op, "match.length") - 1L)
        if (!is.na(close)) skip_to <- close
      }
      next
    }
    tm <- regmatches(look, regexec("^\\s*\\(\\s*([+-]?)\\s*(\\d+)\\s*\\)",
                                   look, perl = TRUE))[[1]]
    if (length(tm) == 3L) {
      t_old    <- as.integer(tm[3L]) * (if (tm[2L] == "-") -1L else 1L)
      consumed <- nchar(tm[1L])
    } else {
      t_old    <- 0L
      consumed <- 0L
    }
    pieces <- c(pieces, substr(text, cursor, starts[j] - 1L),
                repl[[toks[j]]](t_old))
    cursor <- after + consumed
  }
  paste(c(pieces, substr(text, cursor, n)), collapse = "")
}

#' Detrend nonstationary model-block text (steps 1-2 of the rewrite above)
#'
#' @param body        Model-block text on standard timing.
#' @param ns          Validated .mod_nonstationary_decls() result.
#' @param shift_names Every name whose timing a lead/lag of an expression
#'   moves (trend, endogenous and exogenous variables).
#' @return The text with every deflated variable detrended and every trend
#'   variable at a lead/lag rewritten; trend variables at timing 0 are still
#'   present (step 3 is .mod_remove_trend_vars()).
#' @noRd
.mod_detrend_text <- function(body, ns, shift_names) {
  shift <- function(expr, t) {
    if (t == 0L) return(expr)
    .mod_retime_text(expr, stats::setNames(rep(as.integer(t),
                                               length(shift_names)),
                                           shift_names))
  }
  df <- ns$deflated
  for (k in rev(seq_len(nrow(df)))) {
    nm <- df$name[k]
    d  <- df$deflator[k]
    op <- if (df$log[k]) "+" else "*"
    fn <- function(t) paste0("(", nm, .mod_timing_suffix(t), op, "(",
                             shift(d, t), "))")
    body <- .mod_substitute_var_text(body, stats::setNames(list(fn), nm))
  }
  tv <- ns$trend_vars
  repl <- lapply(seq_len(nrow(tv)), function(k) {
    nm <- tv$name[k]
    g  <- tv$growth_factor[k]
    lg <- tv$log[k]
    function(t) {
      if (t == 0L) return(nm)
      steps <- if (t > 0L) seq_len(t) else -(seq_len(-t) - 1L)
      terms <- vapply(steps, function(s) paste0("(", shift(g, s), ")"),
                      character(1))
      if (lg) paste0("(", nm, if (t > 0L) "+" else "-", "(",
                     paste(terms, collapse = "+"), "))")
      else    paste0("(", nm, if (t > 0L) "*" else "/", "(",
                     paste(terms, collapse = "*"), "))")
    }
  })
  names(repl) <- tv$name
  .mod_substitute_var_text(body, repl)
}

## Step 3: trend variables (now all at timing 0) -> 1, additive ones -> 0.
.mod_remove_trend_vars <- function(body, ns) {
  tv <- ns$trend_vars
  repl <- lapply(tv$log, function(lg) {
    val <- if (lg) "0" else "1"
    function(t) val
  })
  names(repl) <- tv$name
  .mod_substitute_var_text(body, repl)
}

#' Balanced-growth check of a detrended model (fail loud, as Dynare does)
#'
#' Replacing the trend variables by 1 (or 0) is legitimate only if no
#' equation depends on the level of a trend variable except through a common
#' factor: F(x, A) = c(A) * F(x, 1).  Dynare tests this
#' (DynamicModel::testTrendDerivativesEqualToZero: the cross derivative
#' d2 log F / dA dx must vanish) and aborts otherwise.  Here each equation of
#' the text left by .mod_detrend_text(), with the trend variables as free
#' symbols, is evaluated at two deterministic points x1, x2 (every variable at
#' every timing in [0.6, 1.4]) and two trend values; balance requires
#'     F(x1, a) * F(x2, a0) == F(x2, a) * F(x1, a0)
#' for each trend variable moved alone (a0 = 1, a = 1.37; additive: 0, 0.37).
#' Equations whose parameters are not all valued, that use a function the
#' evaluator does not know or the EXPECTATION operator, or whose residual is
#' not finite or vanishes at the points, are not tested.
#'
#' @param body_T      Text from .mod_detrend_text().
#' @param ns          The trend declarations.
#' @param var_names   Names of variables in body_T (endogenous + exogenous).
#' @param param_names,param_values Declared parameters and their values.
#' @return Invisibly TRUE; aborts with class
#'   c("dynhr_error_unbalanced_growth", "dynhr_error_mod_syntax").
#' @noRd
.mod_check_balanced_growth <- function(body_T, ns, var_names, param_names,
                                       param_values) {
  if (grepl("\\bEXPECTATION\\b", body_T, perl = TRUE)) return(invisible(TRUE))
  tv <- ns$trend_vars
  if (nrow(tv) == 0L) return(invisible(TRUE))
  parsed <- parse_model_block(body_T, var_names, c(param_names, tv$name))
  ok_fns <- c("exp", "log", "ln", "sqrt", "abs", "sign", "sin", "cos", "tan",
              "asin", "acos", "atan", "sinh", "cosh", "tanh", "normcdf",
              "normpdf", "erf", "cbrt", "max", "min", "STEADY_STATE",
              "steady_state")
  walk <- function(node, what) {
    if (is.null(node)) return(character(0))
    here <- switch(node$type,
      "variable"  = if (what == "var") paste0(node$name, "__",
        if (node$lead_lag >= 0L) "p" else "m", abs(node$lead_lag)),
      "parameter" = if (what == "par") node$name,
      "funcall"   = if (what == "fn") node$name,
      NULL)
    kids <- switch(node$type,
      "binop"   = list(node$left, node$right),
      "unaryop" = list(node$operand),
      "funcall" = node$args,
      list())
    c(here, unlist(lapply(kids, walk, what = what), use.names = FALSE))
  }
  pv <- param_values[!is.na(param_values)]
  bad_eqs <- character(0)
  for (i in seq_along(parsed$equations)) {
    eq  <- parsed$equations[[i]]
    res <- ast_substitute_locals(ast_binop("-", eq$lhs, eq$rhs),
                                 parsed$local_vars)
    if (!all(walk(res, "fn") %in% ok_fns)) next
    pars <- setdiff(unique(walk(res, "par")), tv$name)
    if (!all(pars %in% names(pv))) next
    keys <- unique(walk(res, "var"))
    ## ast_eval() keys: x__0, x__p1, x__m1 (walk() wrote x__p0 for timing 0)
    keys_eval <- sub("__p0$", "__0", keys)
    pt <- function(offset)
      stats::setNames(0.6 + 0.8 * ((seq_along(keys_eval) + offset) *
                                     0.6180339887498949) %% 1, keys_eval)
    x1 <- pt(0L)
    x2 <- pt(7L)
    for (k in seq_len(nrow(tv))) {
      a0 <- stats::setNames(ifelse(tv$log, 0, 1), tv$name)
      a  <- a0
      a[[k]] <- a0[[k]] + 0.37
      f <- function(x, av) ast_eval(res, x, c(pv, av))
      r1a <- f(x1, a); r2a <- f(x2, a); r10 <- f(x1, a0); r20 <- f(x2, a0)
      vals <- c(r1a, r2a, r10, r20)
      if (!all(is.finite(vals)) || abs(r10) < 1e-10 || abs(r20) < 1e-10) next
      lhs <- r1a * r20
      rhs <- r2a * r10
      if (abs(lhs - rhs) > 1e-8 * max(abs(lhs), abs(rhs), 1e-300)) {
        lab <- if (!is.null(eq$tag) && !is.na(eq$tag)) sprintf("%d [%s]", i, eq$tag)
               else as.character(i)
        bad_eqs <- c(bad_eqs, sprintf("equation %s (trend variable %s)",
                                      lab, tv$name[k]))
      }
    }
  }
  if (length(bad_eqs) > 0L)
    .dynhr_abort(
      "parse_mod: the model is not on a balanced growth path: after ",
      "detrending, ", paste(bad_eqs, collapse = "; "), " still depend",
      if (length(bad_eqs) == 1L) "s" else "", " on the level of the trend, ",
      "so the trend variable cannot be removed (Dynare aborts on this too). ",
      "Check the deflator of every trending variable in these equations.",
      class = c("dynhr_error_unbalanced_growth", "dynhr_error_mod_syntax"))
  invisible(TRUE)
}


#' Parse a Dynare .mod file into a dynhr_mod object
#'
#' Reads a Dynare-format `.mod` file (or inline model text) and returns a
#' `dynhr_mod` object containing the parsed equations, declarations, parameter
#' values, shock structure, and any estimation blocks.
#'
#' @details
#' Supported .mod syntax (see \code{vignette("mod-syntax")} for the full
#' reference with verified examples):
#'
#' \emph{Declarations:} \code{var} (with \code{\%(long_name=...)}
#' annotations), \code{varexo}, \code{varexo_det},
#' \code{parameters}, \code{predetermined_variables}, \code{varobs}.
#'
#' \emph{Model block:} \code{model;} / \code{model(linear);} ... \code{end;}.
#' Full nonlinear and linearised equations.  Timing: \code{x(+1)} lead,
#' \code{x(-1)} lag, multi-period via automatic auxiliary variable expansion.
#' Model-local \code{#name = expr;} intermediates stored in
#' \code{model$local_variables} --- no need to inline them.
#' Equation \code{[name=...]} and \code{[mcp='...']} tags.
#' Math: \code{exp}, \code{log}, \code{sqrt}, \code{abs}, trig functions,
#' \code{max}/\code{min}, \code{normcdf}, \code{STEADY_STATE()},
#' standard arithmetic.
#' \code{EXPECTATION(k)(expr)} (expectation of \code{expr} given the
#' information of period \eqn{t+k}, e.g. \code{EXPECTATION(-1)(x(+1))}) is
#' replaced, as in Dynare, by an auxiliary variable
#' \code{AUX_EXPECT_LAG_<|k|>_<n>} (\code{..._LEAD_...} for \eqn{k \ge 0})
#' defined by \code{AUX = expr} shifted by \eqn{-k} periods and entering the
#' equation as \code{AUX(k)}.  A nested \code{EXPECTATION}, an argument that
#' uses a \code{#} model-local variable, and \code{EXPECTATION} outside the
#' model block are errors of class \code{dynhr_error_mod_expectation}.
#' Variables in \code{predetermined_variables} are re-timed by \eqn{-1}
#' everywhere in the model block, \code{#} local definitions included, before
#' the \code{EXPECTATION} substitution and the auxiliary expansion of
#' leads/lags beyond one period (Dynare's order).
#'
#' \emph{Steady state:} \code{steady_state_model;} ... \code{end;} ---
#' sequential assignments evaluated by \code{eval_steady_state_model()}.
#'
#' \emph{Shocks:} \code{shocks;} ... \code{end;} --- \code{stderr},
#' \code{variance}, \code{corr}, \code{covar}, deterministic period shocks,
#' \code{shocks(overwrite)}.  Multiple \code{shocks} blocks are merged.
#'
#' \emph{Shock groups:} \code{shock_groups(name = g);} \code{'supply' = e_a,
#' e_z;} ... \code{end;} --- parsed into \code{model$shock_groups} (the first
#' block, which \code{\link{historical_decomposition}} uses as its default
#' grouping) and \code{model$shock_groups_blocks} (every block, keyed by its
#' \code{name=} option).  A group that names an undeclared shock is an error.
#'
#' \emph{Macro language (Dynare 6/7):} \code{@#define} (values, and
#' functions \code{@#define f(x) = ...}), \code{@#for}/\code{@#endfor}
#' (ranges, arrays, tuple unpacking \code{@#for (i, j) in X}, \code{when}
#' filters), \code{@#if}/\code{@#elseif}/\code{@#else}/\code{@#endif},
#' \code{@#ifdef}/\code{@#ifndef}, \code{@#include}, \code{@#includepath},
#' \code{@#error} (aborts when reached in a taken branch), \code{@#echo},
#' \code{@#echomacrovars}, and \code{@\{expr\}} interpolation.  Macro
#' expressions have Dynare's reals, strings, booleans, arrays and tuples,
#' string/array concatenation with \code{+}, set operations, comprehensions,
#' casts, indexing and built-in functions; they are evaluated by dynhr's own
#' interpreter, never as R code.  Expanded before any other parsing.
#'
#' \emph{Other:} \code{initval}/\code{endval}, \code{histval} (pre-sample
#' history; \code{y(0)} is lag 1, \code{y(-1)} is lag 2, stored lag-indexed in
#' \code{model$histval}), \code{estimated_params},
#' \code{estimated_params_init}, \code{occbin_constraints},
#' \code{planner_objective}, \code{verbatim} blocks.
#' Comments: \code{//}, \code{/* ... */}, \code{\%}.
#'
#' \emph{Syntax errors:} an equation with text left over after a complete
#' expression (e.g. a stray \code{)}), an unrecognised character, or an
#' unsupported operator such as \code{!=} is an error of class
#' \code{dynhr_error_mod_syntax} naming the equation; nothing is skipped.
#'
#' \emph{var(log):} \code{var(log) c k;} declares \code{c} and \code{k} and,
#' as in Dynare 6+, also creates the endogenous variables \code{LOG_c} and
#' \code{LOG_k}. Every occurrence of \code{c} in the model block (with any
#' timing, \code{#} locals included) is replaced by \code{exp(LOG_c)} and the
#' equation \code{c = exp(LOG_c);} is appended, so the perturbation is taken
#' in \code{LOG_c} while \code{c} stays available as a (static) variable. The
#' \code{LOG_*} variables are appended after the declared ones, before any
#' \code{AUX_*} variable; \code{model$log_vars} maps each declared name to its
#' \code{LOG_} name. A \code{steady_state_model} that assigns \code{c} gets
#' \code{LOG_c = log(c);} appended (unless it assigns \code{LOG_c} itself), and
#' \code{initval}/\code{endval}/\code{histval} values of \code{c} are carried
#' to \code{LOG_c} as their logs. A \code{LOG_c} name that is already declared
#' is an error of class \code{dynhr_error_mod_syntax}.
#'
#' \emph{observation_trends:} \code{observation_trends; y (g); end;} adds the
#' deterministic linear trend \code{g * t} to the measurement equation of the
#' observed variable \code{y}, as Dynare does: the observation intercept in
#' period \code{t} is \code{ys[y] + g * t}, where \code{t = first_obs} for the
#' first row of the data and \code{first_obs} is the \code{estimation(first_obs
#' = ...)} option when present (1 otherwise). The slope may be any expression
#' in declared parameters, and is re-evaluated at the parameter vector in force
#' (so an estimated slope moves the likelihood). Stored in
#' \code{model$observation_trends} (\code{$trends}: observable -> expression
#' text; \code{$first_obs}); honoured by \code{\link{kalman_filter}},
#' \code{\link{kalman_smoother}} and everything built on them. An entry for a
#' variable that is not in \code{varobs} is an error, as in Dynare.
#'
#' \emph{Optimal simple rules:} \code{osr_params}, \code{optim_weights}
#' (variance entries \code{y w;} and covariance entries \code{y, pie w;}) and
#' \code{osr_params_bounds} (\code{param, lower, upper;}) are stored in
#' \code{model$osr} (\code{$params}; \code{$weights}, a data frame
#' \code{var1, var2, expr}; \code{$bounds}, a data frame
#' \code{name, lower, upper}; weights and bounds are parameter expressions
#' kept as text) and drive \code{\link{osr}} when its \code{free_params} /
#' \code{loss_vars} are omitted.
#'
#' \emph{model_replace / model_remove / var_remove (Dynare 6+):} applied in
#' file order before the model is parsed.  \code{model_remove('eq1', key =
#' 'value');} removes every equation carrying a listed tag (a bare string is
#' the \code{name} tag; a tag that matches nothing is an error).  Each removed
#' equation must have an \code{endogenous = 'v'} tag or a single endogenous
#' variable on its left-hand side; that variable becomes exogenous if the
#' remaining equations still use it (placed among the exogenous variables in
#' declaration order) and is dropped otherwise.
#' \code{model_replace('eq1'); EQUATIONS end;} removes the tagged equations
#' and appends \code{EQUATIONS} at the end of the model, with no change of
#' variable type.  \code{var_remove a b;} removes declared symbols
#' (variables or parameters); removing one the model still uses is an error.
#' Errors are of class \code{dynhr_error_mod_syntax}.
#'
#' \emph{diff():} \code{diff(EXPR)} is \code{EXPR - EXPR(-1)}, the lag moving
#' every variable and \code{#} local in \code{EXPR} (nested \code{diff} and
#' exogenous arguments allowed).  Dynare's preprocessor introduces
#' \code{AUX_DIFF_*} variables for it; dynhr expands it inline, which gives
#' the same equations in the model's own variables, with the usual
#' \code{AUX_ENDO_LAG}/\code{AUX_EXO_LAG} auxiliaries for the lags.
#' \code{adl()} is not an operator of Dynare 7 and is an error.
#'
#' \emph{ramsey_constraints:} \code{ramsey_constraints; i > 0; end;} bounds
#' a variable of the Ramsey problem by a parameter expression (\code{>} /
#' \code{>=} lower, \code{<} / \code{<=} upper).  Stored in
#' \code{model$ramsey_constraints} (a data frame \code{var, op, bound}, the
#' bound kept as text) and honoured by \code{\link{ramsey_obc_pf}}: each
#' constraint is complementary to the Ramsey first-order condition with
#' respect to its variable, as in Dynare 7, which uses the block only in its
#' mixed-complementarity perfect-foresight solver
#' (\code{perfect_foresight_solver(lmmcp)}); stochastic Ramsey solutions
#' ignore it there too.  A variable that is not declared, or a bound that uses
#' anything but parameters, is an error of class \code{dynhr_error_mod_syntax}.
#'
#' \emph{Command options} are split on top-level commas only, so list values
#' such as \code{instruments=(i,tau)}, \code{irf_shocks=(e,u)} or
#' \code{bandpass_filter=[6 32]} are read whole: a list of names becomes a
#' character vector, a list of numbers (or integer ranges \code{m:n}) a
#' numeric vector.  \code{\link{discretionary_policy}} takes its instruments
#' and discount from the \code{discretionary_policy} command and
#' \code{\link{osr}} its planner discount from the \code{osr} command.
#'
#' \emph{Nonstationary models (trend_var, var(deflator=...)):}
#' \code{trend_var(growth_factor = G) A;} declares a trend variable
#' (\code{A_t = G_t A_{t-1}}), \code{log_trend_var(log_growth_factor = g)
#' LA;} an additive one (\code{LA_t = g_t + LA_{t-1}}); \code{var(deflator =
#' D) Y;} declares \code{Y} as trending with the multiplicative trend
#' \code{D} (an expression in trend variables, endogenous variables and
#' parameters) and \code{var(log_deflator = D) Y;} with the additive trend
#' \code{D}.  The model block is written in the nonstationary variables and,
#' as Dynare's preprocessor does, rewritten into the stationary model before
#' anything else (after the \code{predetermined_variables} re-timing): every
#' \code{Y(k)} becomes \code{Y(k)*D(k)} (\code{Y(k)+D(k)}), a trend variable
#' at a lead or lag becomes \code{A*G(+1)*...*G(+k)} or
#' \code{A/(G*G(-1)*...*G(-k+1))} (sums for an additive trend), and then
#' every trend variable becomes 1 (additive: 0).  So every declared name
#' stands for its DETRENDED variable, which is what the steady state,
#' \code{steady_state_model}, \code{initval} and the decision rule refer to.
#' The rewrite is valid only on a balanced growth path: an equation that
#' still depends on the level of a trend after detrending is an error of
#' class \code{dynhr_error_unbalanced_growth} (Dynare aborts too).  A growth
#' factor that uses a trend variable (an I(2) trend) is an error.  The
#' declarations are stored in \code{model$nonstationary}
#' (\code{$trend_vars}: \code{name, growth_factor, log};
#' \code{$deflated}: \code{name, deflator, log}).
#'
#' \emph{deterministic_trends:} \code{deterministic_trends; y (g); end;} has
#' the \code{observation_trends} entry syntax but names any declared
#' endogenous variable; it is stored in \code{model$deterministic_trends}, and
#' each entry for an observed variable (every entry when there is no
#' \code{varobs}) adds the linear trend \code{g * t} to that observable's
#' measurement equation exactly as \code{observation_trends} does (it is
#' merged into \code{model$observation_trends}).  A variable in both blocks is
#' an error.  (Dynare 7.1 stores the block in
#' \code{options_.deterministic_trend_coeffs}, which none of its routines
#' reads, and its generated driver fails on the statement.)
#'
#' @param file_or_text Either a file path (ending in \code{.mod} or
#'   \code{.dyn}) or a character string containing \code{.mod} source text.
#'   Paths that look like filenames but do not exist produce a clear
#'   "file not found" error.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#' @return A \code{dynhr_mod} object.  Key fields: \code{$var_names},
#'   \code{$varexo_names}, \code{$param_names}, \code{$param_values},
#'   \code{$equations}, \code{$local_variables}, \code{$shocks},
#'   \code{$steady_state_model}, \code{$initval}, \code{$estimated_params},
#'   \code{$lead_lag_incidence}, \code{$variable_classification}.
#' @seealso \code{vignette("mod-syntax")} for the complete syntax reference
#'   with verified examples; \code{vignette("mod-conversion")} for converting
#'   Dynare models and using dynhr-specific annotations.
#' @export
parse_mod <- function(file_or_text, verbose = FALSE) {

  # ---- 1. Read input ---------------------------------------------------
  if (length(file_or_text) != 1L || !is.character(file_or_text) ||
      is.na(file_or_text)) {
    stop("parse_mod(): `file_or_text` must be a single non-NA character ",
         "string (a .mod path or model source text).", call. = FALSE)
  }
  if (file.exists(file_or_text)) {
    source_file <- file_or_text
    # Read raw bytes (no encoding conversion), then explicitly convert
    # from Latin-1 to UTF-8.  This handles .mod files with non-ASCII
    # characters (e.g. accented author names) that would otherwise cause
    # gsub(perl=TRUE) to fail with "invalid multibyte string".
    txt <- paste(readLines(file_or_text, warn = FALSE), collapse = "\n")
    ## A file that is already valid UTF-8 (pure ASCII included) is kept as
    ## is: re-reading it as Latin-1 would mangle Dynare 7's `U+27C2`
    ## complementarity separator (U+27C2, three UTF-8 bytes) into three
    ## Latin-1 characters.  Anything else is taken to be Latin-1, as before.
    if (isTRUE(validUTF8(txt))) Encoding(txt) <- "UTF-8"
    else txt <- iconv(txt, from = "LATIN1", to = "UTF-8", sub = "?")
    # Remove BOM if present
    txt <- gsub("^\uFEFF", "", txt)
    if (verbose) .dynhr_cat("Read", nchar(txt), "characters from", source_file, "\n")
  } else {
    # An input that LOOKS like a file path (no newlines, ends in a model
    # extension) but does not exist is almost certainly a mistyped path, not
    # inline source text.  Fail loudly here -- otherwise it falls through and
    # is parsed as "text" with zero declarations/equations, producing the
    # cryptic "argument is of length zero" crash downstream (repl H7).
    looks_like_path <- !grepl("\n", file_or_text, fixed = TRUE) &&
      grepl("\\.(mod|dyn)$", file_or_text, ignore.case = TRUE)
    if (looks_like_path) {
      stop("parse_mod(): file not found: ", file_or_text,
           "\n(Pass an existing .mod/.dyn path, or inline model source text.)",
           call. = FALSE)
    }
    source_file <- NA_character_
    txt <- file_or_text
  }

  # ---- 2. Expand Dynare macro directives -------------------------
  # Evaluate the supported @#for / @#if / @#ifdef / @#ifndef / @#define macro
  # subset and splice in @{...} interpolations BEFORE the lexer and any block
  # extraction run.  This is an exact byte-for-byte no-op for macro-free
  # source, so models without @# directives are entirely unaffected.
  # Unsupported directives FAIL LOUD (see R/parse-macro.R) rather than being
  # silently dropped, which would yield an under-specified model.
  txt <- expand_macros(txt,
                       mod_dir = if (!is.na(source_file)) dirname(source_file) else NULL)
  if (verbose) .dynhr_cat("Expanded Dynare macro directives\n")

  # ---- 3. Strip comments and macros ------------------------------------
  txt <- strip_comments_and_macros(txt)
  if (verbose) .dynhr_cat("Stripped comments and macros\n")

  # ---- 3. Extract declarations -----------------------------------------
  # Strip block constructs to prevent internal 'var' keywords
  # (e.g. "var eps_a = sig_a^2;" inside shocks block) from being
  # mistaken for variable declarations.
  txt_decl <- .strip_blocks(.strip_mod_comments(txt))

  var_text       <- extract_declaration(txt_decl, "var")
  varexo_text     <- extract_declaration(txt_decl, "varexo")
  varexo_det_text <- extract_declaration(txt_decl, "varexo_det")
  param_text      <- extract_declaration(txt_decl, "parameters")
  predet_text     <- extract_declaration(txt_decl, "predetermined_variables")
  varobs_text     <- extract_declaration(txt_decl, "varobs")

  var_names       <- parse_declaration_names(var_text)
  varexo_names    <- parse_declaration_names(varexo_text)
  varexo_det_names <- parse_declaration_names(varexo_det_text)
  param_names     <- parse_declaration_names(param_text)
  predet_vars     <- parse_declaration_names(predet_text)
  varobs_names    <- parse_declaration_names(varobs_text)

  # ---- Deduplicate: exo/exo_det names must not appear in endo list ----
  exo_overlap     <- intersect(var_names, varexo_names)
  exo_det_overlap <- intersect(var_names, varexo_det_names)
  if (length(exo_overlap) > 0) {
    .dynhr_inform(sprintf("Parser: removing %d varexo names from var: %s",
                    length(exo_overlap), paste(exo_overlap, collapse = ", ")))
    var_names <- setdiff(var_names, varexo_names)
  }
  if (length(exo_det_overlap) > 0) {
    .dynhr_inform(sprintf("Parser: removing %d varexo_det names from var: %s",
                    length(exo_det_overlap), paste(exo_det_overlap, collapse = ", ")))
    var_names <- setdiff(var_names, varexo_det_names)
  }
  # Also protect varexo from varexo_det contamination
  exo_cross <- intersect(varexo_names, varexo_det_names)
  if (length(exo_cross) > 0) {
    .dynhr_inform(sprintf("Parser: removing %d varexo_det names from varexo: %s",
                    length(exo_cross), paste(exo_cross, collapse = ", ")))
    varexo_names <- setdiff(varexo_names, varexo_det_names)
  }

  all_var_names <- c(var_names, varexo_names, varexo_det_names)

  # ---- 3b. model_replace / model_remove / var_remove (Dynare 6+) -------
  # Applied to the model-block TEXT, in file order, before anything else
  # touches the model: they choose which equations exist and may turn an
  # endogenous variable exogenous or remove symbols.
  model_block <- extract_paired_block(txt, "model")
  .edits <- .mod_apply_model_edits(
    txt, if (model_block$found) model_block$body else "",
    list(endo = var_names, exo = varexo_names, exo_det = varexo_det_names,
         params = param_names), txt_decl)
  if (.edits$changed) {
    model_block$body <- .edits$body
    var_names        <- .edits$endo
    varexo_names     <- .edits$exo
    varexo_det_names <- .edits$exo_det
    param_names      <- .edits$params
    all_var_names    <- c(var_names, varexo_names, varexo_det_names)
    if (verbose) .dynhr_cat("Applied model_replace/model_remove/var_remove\n")
  }

  ## Nonstationary model: trend_var / log_trend_var and var(deflator=) /
  ## var(log_deflator=) declarations (NULL when there are none). The model
  ## block is detrended in step 4; the balanced-growth check runs after the
  ## calibration (step 5b), which it needs.
  nonstationary <- .mod_nonstationary_decls(txt_decl)
  if (!is.null(nonstationary))
    .mod_check_nonstationary(nonstationary, var_names,
                             c(varexo_names, varexo_det_names), param_names)
  .ns_body_T <- NULL

  ## var(log) declarations: declared name -> its LOG_ variable (Dynare's
  ## naming). Filled in step 4 only when there is a model block to rewrite.
  log_map <- character(0)
  .log_decl <- intersect(.mod_log_var_names(txt_decl), var_names)

  if (verbose) {
    .dynhr_cat("Variables:  ", length(var_names), "endogenous,",
        length(varexo_names), "exogenous,",
        length(varexo_det_names), "exo deterministic\n")
    .dynhr_cat("Parameters:", length(param_names), "\n")
  }

  # ---- 4. Parse the model block (extracted, and edited, in step 3b) -----
  model_opts  <- parse_command_options(model_block$options_str)

  equations   <- list()
  local_vars  <- list()

  if (model_block$found) {

    # ---- ONE ordered timing pass over the model-block TEXT ----
    # The model body -- equations AND `#` model-local definitions, which are
    # part of the same text -- goes through these steps IN THIS ORDER, before
    # anything is parsed:
    #   (1) predetermined re-timing: every occurrence of a variable declared in
    #       `predetermined_variables` is shifted by -1 (k -> k(-1),
    #       k(+1) -> k, k(-1) -> k(-2)).  Locals are re-timed too, because
    #       they are text in the same body.
    #   (2) EXPECTATION(k)(expr) -> auxiliary endogenous variable.  It
    #       runs on standard timing, so a predetermined variable inside expr
    #       is already re-timed.
    #   (3) auxiliary expansion of |lead|/|lag| > 1 (and exogenous +-1) into
    #       AUX chains.  It sees standard timing, so a predetermined `k(-1)`
    #       (standard k(-2)) gets its AUX_LAG chain and a predetermined
    #       `k(+2)` (standard k(+1)) needs none.
    # This is Dynare's order: ModFile::transformPass re-times predetermined
    # variables (DynamicModel::transformPredeterminedVariables, which also
    # re-times the model-local variables), then substitutes the expectation
    # operator (DynamicModel::substituteExpectation), then the endo/exo
    # leads and lags (substituteEndoLeadGreaterThanTwo, ...LagGreaterThanTwo,
    # substituteExoLead, substituteExoLag).  Verified against the Dynare 7.1
    # preprocessor's `json=transform` output.  The equation ASTs built below
    # are therefore already in standard timing; nothing is re-timed after
    # parsing.
    # (0a) Dynare 7 complementarity conditions `EQ U+27C2 COND;` / `EQ _|_ COND;`
    # are cut off the equations FIRST, while the body is still the text the
    # user wrote: the condition names the constrained variable, which none of
    # the timing passes below may touch (a predetermined `i` must not become
    # `i(-1)` in `U+27C2 i > 0`).  They are re-attached to their equations after
    # parsing.
    .cc <- .mod_extract_complementarity(
      model_block$body, endo_names = var_names,
      exo_names = c(varexo_names, varexo_det_names), param_names = param_names)
    model_block$body <- .cc$body
    # (0a') diff(EXPR) -> ((EXPR) - (EXPR shifted by -1)).  A `#` local inside
    # EXPR is shifted as `local(-1)`, which (0b) then resolves; every later
    # pass is a uniform re-timing or a substitution, which commutes with it.
    .loc_d <- names(.mod_local_defs(.mod_model_units(model_block$body)))
    model_block$body <- .mod_expand_diff_text(
      model_block$body, c(var_names, varexo_names, varexo_det_names, .loc_d),
      declared = c(all_var_names, param_names))
    # (0b) Dynare 7.0: a `#` model-local used with a lead or lag, `r(-1)`, is
    # replaced by its definition shifted by that lead/lag.  Done before the
    # other passes, on the user's timing; every later pass is a uniform
    # re-timing or a substitution that commutes with this one.
    model_block$body <- .mod_substitute_timed_locals(
      model_block$body, c(var_names, varexo_names, varexo_det_names))

    .bad_predet <- setdiff(predet_vars, var_names)
    if (length(.bad_predet) > 0L)
      .dynhr_abort(
        "parse_mod: `predetermined_variables` must list declared endogenous ",
        "(`var`) variables; not endogenous: ",
        paste(.bad_predet, collapse = ", "), ".",
        class = "dynhr_error_mod_syntax")
    # (1) predetermined re-timing
    if (length(predet_vars) > 0L) {
      model_block$body <- .mod_retime_text(
        model_block$body, stats::setNames(rep(-1L, length(predet_vars)),
                                          predet_vars))
      if (verbose) .dynhr_cat("  Re-timed", length(predet_vars),
                              "predetermined var(s) by -1\n")
    }
    # (1a) nonstationary model -> stationary (detrended) model, as Dynare's
    # preprocessor does right after the predetermined re-timing
    # (DynamicModel::detrendEquations, removeTrendVariableFromEquations);
    # see the comment above .mod_nonstationary_decls().  Before var(log),
    # EXPECTATION and the aux expansion, which all see the detrended model.
    if (!is.null(nonstationary)) {
      .ns_body_T <- .mod_detrend_text(
        model_block$body, nonstationary,
        shift_names = c(nonstationary$trend_vars$name, var_names,
                        varexo_names, varexo_det_names))
      model_block$body <- .mod_remove_trend_vars(.ns_body_T, nonstationary)
      if (verbose) .dynhr_cat("  Detrended",
                              nrow(nonstationary$deflated), "variable(s),",
                              nrow(nonstationary$trend_vars),
                              "trend variable(s)\n")
    }
    # (1b) var(log): x -> exp(LOG_x) everywhere, plus `x = exp(LOG_x);`.
    # After the predetermined re-timing, so the appended definitional
    # equation is never re-timed; before EXPECTATION and the aux expansion, so
    # a LOG_x at |timing| > 1 gets its AUX chain like any other variable.
    if (length(.log_decl) > 0L) {
      log_map <- stats::setNames(paste0("LOG_", .log_decl), .log_decl)
      .loc_m <- regmatches(model_block$body, gregexpr(
        "(?:^|;)\\s*(?:\\[[^\\]]*\\]\\s*)?#\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=",
        model_block$body, perl = TRUE))[[1]]
      .loc_names <- sub("^.*#\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=$", "\\1",
                        .loc_m, perl = TRUE)
      .clash <- intersect(unname(log_map),
                          c(all_var_names, param_names, .loc_names))
      if (length(.clash) > 0L)
        .dynhr_abort(
          "parse_mod: var(log) creates the variable(s) ",
          paste(.clash, collapse = ", "), ", but that name is already ",
          "declared or used. Rename it; Dynare reserves LOG_<name> for the ",
          "log of a var(log) variable.", class = "dynhr_error_mod_syntax")
      model_block$body <- paste0(
        .mod_log_substitute_text(model_block$body, log_map), "\n",
        paste0(names(log_map), " = exp(", log_map, ");", collapse = "\n"),
        "\n")
      var_names     <- c(var_names, unname(log_map))
      all_var_names <- c(var_names, varexo_names, varexo_det_names)
      if (verbose) .dynhr_cat("  var(log):", length(log_map),
                              "LOG_ variable(s)\n")
    }
    # (2) EXPECTATION operator
    exp_result <- .substitute_expectation_text(
      model_block$body, var_names,
      exo_names = c(varexo_names, varexo_det_names),
      reserved_names = c(all_var_names, param_names))
    if (length(exp_result$aux_var_names) > 0L) {
      model_block$body <- exp_result$model_body
      var_names     <- c(var_names, exp_result$aux_var_names)
      all_var_names <- c(var_names, varexo_names, varexo_det_names)
      if (verbose) .dynhr_cat("  EXPECTATION operator:",
                              length(exp_result$aux_var_names),
                              "auxiliary variable(s)\n")
    }

    # ---- (3) Auxiliary variable expansion for leads/lags > 1 ----
    # Endogenous names need aux vars only for |timing| > 1; EXOGENOUS names
    # (passed via exo_names) need an AUX_EXO_LEAD/LAG endo chain for ANY nonzero
    # timing (incl. +-1), since the perturbation state space cannot carry an exo
    # at a nonzero timing (H8; e.g. eps_z_news(-8), ed(+1)).
    exo_search_names <- c(varexo_names, varexo_det_names)
    aux_result <- .expand_aux_timing(model_block$body, var_names,
                                     exo_names = exo_search_names,
                                     verbose = verbose)
    if (length(aux_result$aux_var_names) > 0) {
      model_block$body <- aux_result$model_body
      # AUX variables become endogenous vars regardless of whether they
      # originated from an endo or exo name
      var_names     <- c(var_names, aux_result$aux_var_names)
      all_var_names <- c(var_names, varexo_names, varexo_det_names)
    }
    # ---- End auxiliary expansion ----

    parsed <- parse_model_block(model_block$body,
                                all_var_names, param_names)
    equations  <- parsed$equations
    local_vars <- parsed$local_vars
    # Re-attach the complementarity conditions cut off in (0a).  Every pass in
    # between only APPENDS equations (var(log) definitions, EXPECTATION and
    # lead/lag auxiliaries), so the k-th original equation is still
    # equations[[k]].
    for (.k in names(.cc$conditions))
      equations[[as.integer(.k)]]$complementarity <- .cc$conditions[[.k]]
    if (verbose) .dynhr_cat("Parsed", length(equations), "equations\n")
  }

  # ---- 5. Calibration (top-level param assignments) --------------------
  remaining_txt <- remove_blocks(txt)
  # Pre-evaluate verbatim; blocks (which remove_blocks strips) so derived
  # scalar params and their matrix intermediates are available.  We do a quick
  # first calibration pass to seed the verbatim env with base param values
  # (verbatim RHS may reference declared params like z_bar, rho_zz, ...),
  # evaluate the verbatim blocks, then re-run calibration with that env so
  # both verbatim-internal scalar params AND top-level assignments that
  # reference verbatim intermediates resolve.
  if (grepl("(?si)\\bverbatim\\s*;", txt, perl = TRUE)) {
    base_vals    <- parse_calibration(remaining_txt, param_names, quiet = TRUE)
    verbatim_env <- eval_verbatim_blocks(txt, base_vals)
    param_values <- parse_calibration(remaining_txt, param_names,
                                      seed_env = verbatim_env, quiet = TRUE)
    # Also harvest scalar params that were assigned ONLY inside the verbatim
    # block (e.g. P0_z_bar0) and never re-stated at top level.
    for (pn in setdiff(param_names, names(param_values))) {
      if (exists(pn, envir = verbatim_env, inherits = FALSE)) {
        v <- get(pn, envir = verbatim_env)
        if (is.numeric(v) && length(v) == 1L && is.finite(v))
          param_values[pn] <- v
      }
    }
  } else {
    param_values <- parse_calibration(remaining_txt, param_names)
  }
  if (verbose && length(param_values) > 0)
    .dynhr_cat("Calibrated", length(param_values), "parameters\n")

  # ---- 5b. balanced growth of a detrended model ------------------------
  # Needs the calibrated parameters; aborts when an equation still depends on
  # the level of a trend variable (Dynare aborts there too).
  if (!is.null(.ns_body_T))
    .mod_check_balanced_growth(
      .ns_body_T, nonstationary,
      var_names = all_var_names,
      param_names = param_names, param_values = param_values)

  # ---- fail-loud on declared-but-unvalued parameters --------------------
  # Parameters that are DECLARED in the `parameters` block but never assigned a
  # value in the .mod (no top-level calibration, no verbatim; assignment) are
  # absent from / NA in param_values. This is the single most pervasive
  # replication pitfall: they are typically computed in an external
  # `*_steadystate.m` (which dynhr cannot run), so the solution is silently
  # wrong. Warn and name them, pointing at set_param_values(). We stay silent
  # when the .mod carries a `steady_state_model` block (those params are adopted
  # at solve time via ss$params, so the absence is expected, not an error).
  .unset_params <- union(setdiff(param_names, names(param_values)),
                         names(param_values)[is.na(param_values)])
  if (length(.unset_params) > 0L &&
      !grepl("(?si)\\bsteady_state_model\\b", txt, perl = TRUE)) {
    # Only blame an external steady-state file when one actually sits next to
    # the .mod; otherwise the parameter is simply never assigned (or its
    # assignment could not be evaluated -- parse_calibration() says which).
    .ext_ss <- .external_ss_file(source_file)
    .why <- if (!is.na(.ext_ss))
      sprintf("they may be computed in the external steady-state file '%s', which dynhr cannot run",
              basename(.ext_ss))
    else
      "no top-level assignment of them could be found and evaluated"
    .dynhr_warn(sprintf(
      paste0("parse_mod: %d declared parameter(s) have no value in the .mod ",
             "(%s): %s. ",
             "Use inject_params(model, named_vec) to fill them from an ",
             "external (possibly superset) vector without errors; or use ",
             "set_param_values() if you have the exact model-parameter set. ",
             "Unset parameters will produce silently wrong decision rules."),
      length(.unset_params), .why,
      paste(.unset_params, collapse = ", ")),
      call. = FALSE)
  }

  # ---- 6. Initval / Endval ---------------------------------------------
  # Create an evaluation environment with known parameter values so that
  # initval/endval expressions can reference parameter names directly.
  ## Allowlist sandbox, not a child of baseenv() -- see
  ## `.dynhr_sandbox_env()` in parse-blocks.R.
  initval_env <- .dynhr_sandbox_env(param_values, .dynhr_safe_matrix_fn_names)
  initval_block <- extract_paired_block(txt, "initval")
  initval <- if (initval_block$found)
    parse_initval_block(initval_block$body, env = initval_env)
  else numeric(0)

  endval_block <- extract_paired_block(txt, "endval")
  endval <- if (endval_block$found)
    parse_initval_block(endval_block$body, env = initval_env)
  else numeric(0)
  ## var(log): a starting value for x is a starting value for LOG_x = log(x).
  initval <- .mod_log_carry(initval, log_map)
  endval  <- .mod_log_carry(endval,  log_map)

  # ---- 6b. Histval (pre-sample history) --------------------------------
  # A separate env so that histval assignments do not leak into initval/endval
  # resolution (histval names ARE variables, not parameters).
  histval_env <- .dynhr_sandbox_env(param_values, .dynhr_safe_matrix_fn_names)
  histval_block <- extract_paired_block(txt, "histval")
  histval <- if (histval_block$found)
    parse_histval_block(histval_block$body, env = histval_env)
  else list()
  for (.lv in names(log_map)) {
    if (!is.null(histval[[.lv]]) && is.null(histval[[log_map[[.lv]]]]) &&
        all(is.na(histval[[.lv]]) | histval[[.lv]] > 0))
      histval[[log_map[[.lv]]]] <- log(histval[[.lv]])
  }

  # ---- 7. Steady state model -------------------------------------------
  ssm_block <- extract_paired_block(txt, "steady_state_model")
  ## var(log): LOG_x = log(x) after the block's own assignment of x (Dynare
  ## fills the LOG_ auxiliary's steady state the same way), unless the block
  ## already assigns LOG_x itself.
  if (ssm_block$found && length(log_map) > 0L) {
    .ssm_lhs_txt <- regmatches(ssm_block$body, gregexpr(
      "(?:^|;)\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=(?!=)", ssm_block$body,
      perl = TRUE))[[1]]
    .ssm_lhs_txt <- sub("^;?\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=$", "\\1",
                        .ssm_lhs_txt, perl = TRUE)
    .add <- names(log_map)[names(log_map) %in% .ssm_lhs_txt &
                             !(log_map %in% .ssm_lhs_txt)]
    if (length(.add) > 0L)
      ssm_block$body <- paste0(ssm_block$body, "\n",
                               paste0(log_map[.add], " = log(", .add, ");",
                                      collapse = "\n"))
  }
  ssm <- if (ssm_block$found)
    parse_steady_state_model(ssm_block$body, all_var_names,
                             param_names)
  else NULL

  # ---- M17b: warn on placeholder params overridden by steady_state_model -----
  # A parameter assigned BOTH a top-level .mod calibration value AND a value in
  # the steady_state_model block is a placeholder: the SSM value wins at solve
  # time (eval_steady_state_model overwrites it), so the .mod value is dead and
  # silently misleading anyone reading model$param_values. Name them so the
  # override is visible (the wrong-placeholder pitfall: Gertler_Karadi etc.).
  if (!is.null(ssm) && length(ssm) > 0L) {
    .ssm_lhs <- vapply(ssm, function(a) a$name %||% NA_character_, character(1))
    .ssm_params <- intersect(.ssm_lhs[!is.na(.ssm_lhs)], param_names)
    .ssm_placeholders <- intersect(
      .ssm_params,
      names(param_values)[!is.na(param_values)]
    )
    if (length(.ssm_placeholders) > 0L) {
      .dynhr_warn(sprintf(
        paste0("parse_mod: %d parameter(s) have a .mod calibration value that ",
               "is overridden by the steady_state_model block at solve time ",
               "(the .mod value is a placeholder): %s. Read the solved value ",
               "from the steady state, not model$param_values."),
        length(.ssm_placeholders),
        paste(.ssm_placeholders, collapse = ", ")),
        call. = FALSE)
    }
  }

  # ---- M17c: warn on placeholder params calibrated in external *_steadystate.m
  # Parameters that DO have a .mod calibration value (typically a round-number
  # placeholder like `chi=1;`) but are ALSO assigned inside a sibling
  # `<Model>_steadystate.m` / `<Model>_steadystate2.m` file are silently wrong:
  # Dynare overwrites them at runtime but dynhr cannot run the .m file.
  # We detect the external file by convention (same directory, stem + _steadystate),
  # then scan it for lines that assign any declared parameter name.
  # Heuristic: a line is considered an assignment of param P when it matches
  #   \bP\s*= (bare LHS) OR M_.params(...) assignment containing 'P' in a
  #   strmatch/strcmp lookup.  False-positive risk: a commented-out line or a RHS
  #   reference to P could match; this is acceptable (over-warn, not under-warn).
  if (!is.na(source_file)) {
    .ss_file <- .external_ss_file(source_file)
    if (!is.na(.ss_file)) {
      .ss_txt  <- paste(readLines(.ss_file, warn = FALSE), collapse = "\n")
      # Params that have a .mod calibration value (the placeholder set)
      .calibrated_params <- names(param_values)[!is.na(param_values)]
      # Detect which of those appear to be assigned in the .m file
      .ext_assigned <- Filter(function(p) {
        # bare LHS assignment: \bP\s*=  (not == or ~=)
        bare_pat  <- paste0("\\b", p, "\\s*=[^=~]")
        # M_.params(...) pattern with 'P' in a string lookup
        mp_pat    <- paste0("M_\\.params\\([^)]*['\"]", p, "['\"]")
        grepl(bare_pat, .ss_txt, perl = TRUE) ||
          grepl(mp_pat,   .ss_txt, perl = TRUE)
      }, .calibrated_params)
      if (length(.ext_assigned) > 0L) {
        .dynhr_warn(sprintf(
          paste0("parse_mod: %d parameter(s) have a .mod calibration value that ",
                 "is overridden by the external steady-state file '%s' at solve ",
                 "time (the .mod value is a placeholder): %s. ",
                 "Use inject_params(model, named_vec) to overwrite these from an ",
                 "external (possibly superset) vector without errors; or use ",
                 "set_param_values() if you have the exact model-parameter set. ",
                 "Placeholder values will produce silently wrong decision rules."),
          length(.ext_assigned),
          basename(.ss_file),
          paste(.ext_assigned, collapse = ", ")),
          call. = FALSE)
      }
    }
  }

  # ---- 8. Shocks -------------------------------------------------------
  # Build an augmented eval environment that includes BOTH declared parameter
  # values AND any model-local numeric constants (e.g. `phi = 0.1;` not in
  # `parameters`).  This lets shocks-block expressions like
  # `var e_a, e_b = phi*sig_a*sig_b;` resolve correctly.
  local_const_env <- .dynhr_sandbox_env(param_values, .dynhr_safe_matrix_fn_names)
  .lc_pat <- "([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*([0-9][0-9.eE+\\-]*)\\s*;"
  .lc_ms  <- regmatches(remaining_txt,
                         gregexpr(.lc_pat, remaining_txt, perl = TRUE))[[1]]
  for (.lc_m in .lc_ms) {
    .lc_parts <- regmatches(.lc_m, regexec(.lc_pat, .lc_m, perl = TRUE))[[1]]
    .lc_val   <- suppressWarnings(as.numeric(.lc_parts[3]))
    if (!is.na(.lc_val) && !exists(.lc_parts[2], envir = local_const_env,
                                    inherits = FALSE))
      assign(.lc_parts[2], .lc_val, envir = local_const_env)
  }

  # Parse ALL shocks blocks and merge them cumulatively (M5 fix).
  # `shocks(overwrite)` replaces matching vars; plain `shocks` appends /
  # last-declaration-wins for variances.
  .empty_variances    <- data.frame(name = character(0),
                                    stderr = numeric(0),
                                    variance = numeric(0),
                                    stderr_expr = character(0),
                                    variance_expr = character(0),
                                    skew = numeric(0),
                                    skew_expr = character(0),
                                    stringsAsFactors = FALSE)
  .empty_correlations <- data.frame(var1 = character(0),
                                    var2 = character(0),
                                    corr = numeric(0),
                                    corr_expr = character(0),
                                    cov  = numeric(0),
                                    cov_expr  = character(0),
                                    stringsAsFactors = FALSE)

  .merge_variances <- function(base_df, new_df, overwrite = FALSE) {
    if (nrow(new_df) == 0) return(base_df)
    if (nrow(base_df) == 0) return(new_df)
    # For each name in new_df, replace or append
    for (nm in new_df$name) {
      idx <- match(nm, base_df$name)
      if (!is.na(idx)) {
        nr <- new_df[new_df$name == nm, ][1, ]
        # A later shocks block silently overwriting an earlier variance for
        # the same shock loses the earlier value (e.g. Ascari keeps only the last
        # block's stderrs). Warn on a genuine conflict, mirroring the correlation
        # dedup warning. Re-stating the same value is silent (byte-compatible).
        old_v <- base_df$variance[idx]
        new_v <- nr$variance
        differs <- (is.na(old_v) != is.na(new_v)) ||
          (!is.na(old_v) && !is.na(new_v) &&
             !isTRUE(all.equal(old_v, new_v)))
        if (differs) {
          .dynhr_warn(sprintf(
            paste0("parse_mod: shock variance for '%s' is defined in more than ",
                   "one shocks block; using the last. Scope each shocks block to ",
                   "its stoch_simul if they are meant to differ."),
            nm), call. = FALSE)
        }
        base_df[idx, ] <- nr
      } else {
        base_df <- rbind(base_df, new_df[new_df$name == nm, ])
      }
    }
    base_df
  }

  .merge_correlations <- function(base_df, new_df) {
    if (nrow(new_df) == 0) return(base_df)
    if (nrow(base_df) == 0) return(new_df)
    # A second `shocks` block must not silently append a duplicate row for a
    # pair already defined (plain rbind let last-write-wins apply downstream with
    # no signal -- e.g. BKK's corr=0 IRF block vs corr=0.258 sim block). Replace
    # per (unordered) pair, keeping the last definition, and warn on a conflict.
    for (k in seq_len(nrow(new_df))) {
      nr  <- new_df[k, , drop = FALSE]
      idx <- which((base_df$var1 == nr$var1 & base_df$var2 == nr$var2) |
                   (base_df$var1 == nr$var2 & base_df$var2 == nr$var1))
      if (length(idx) > 0) {
        .dynhr_warn(sprintf(
          paste0("parse_mod: shock correlation for (%s, %s) is defined in more ",
                 "than one shocks block; using the last. Scope each shocks block ",
                 "to its stoch_simul if they are meant to differ."),
          nr$var1, nr$var2), call. = FALSE)
        base_df[idx[1], ] <- nr
        if (length(idx) > 1L) base_df <- base_df[-idx[-1L], , drop = FALSE]
      } else {
        base_df <- rbind(base_df, nr)
      }
    }
    base_df
  }

  shocks_blocks <- extract_all_paired_blocks(txt, "shocks")
  .empty_det <- data.frame(name = character(0), period = integer(0),
                            value = numeric(0), stringsAsFactors = FALSE)
  shocks <- list(variances = .empty_variances,
                 correlations = .empty_correlations)
  det_shocks <- .empty_det
  # Retain a per-block breakdown so per-`stoch_simul` scoping can be added
  # later without another parser rewrite. Each entry carries the block's own
  # parsed contents plus its `overwrite` flag. The merged `shocks` remains the
  # primary field for backward compatibility.
  shocks_blocks_parsed <- vector("list", length(shocks_blocks))
  for (.bi in seq_along(shocks_blocks)) {
    .blk    <- shocks_blocks[[.bi]]
    .ovr    <- grepl("overwrite", .blk$options_str, ignore.case = TRUE)
    .parsed <- parse_shocks_block(.blk$body, param_env = local_const_env)
    if (.ovr) {
      # Dynare `shocks(overwrite)` discards the entire prior shock
      # specification, then applies its own entries fresh.
      shocks$variances    <- .empty_variances
      shocks$correlations <- .empty_correlations
      det_shocks          <- .empty_det
    }
    shocks$variances    <- .merge_variances(shocks$variances,
                                            .parsed$variances)
    shocks$correlations <- .merge_correlations(shocks$correlations,
                                               .parsed$correlations)
    if (!is.null(.parsed$deterministic) && nrow(.parsed$deterministic) > 0)
      det_shocks <- rbind(det_shocks, .parsed$deterministic)
    shocks_blocks_parsed[[.bi]] <- list(
      overwrite     = .ovr,
      options       = .blk$options_str,
      variances     = .parsed$variances,
      correlations  = .parsed$correlations,
      deterministic = .parsed$deterministic
    )
  }
  if (length(shocks_blocks) > 1L) {
    # M5 (fail-loud): the union is merged with later-block-wins / overwrite
    # semantics, but dynhr does not yet scope shocks to individual
    # `stoch_simul` calls. Surface the approximation and point at the breakdown.
    .dynhr_warn(sprintf(
      paste0("parse_mod: %d shocks blocks merged into m$shocks (later blocks ",
             "win per shock; shocks(overwrite) resets first). Per-stoch_simul ",
             "scoping is not modelled - inspect m$shocks_blocks for the ",
             "per-block breakdown if separate experiments use different shocks."),
      length(shocks_blocks)), call. = FALSE)
  }

  # ---- 8b. filter_tunes -------------------------------------------------
  # Sample-start date for resolving date-literal periods: taken from
  # estimation(first_obs=...) if present (e.g. "1990q1"). Purely a string
  # lookup -- the `commands` list (step 10) is built after this point.
  ft_block <- extract_paired_block(txt, "filter_tunes")
  filter_tunes <- if (ft_block$found) {
    sample_start <- NULL
    fo_ft <- .estimation_first_obs(txt)
    if (!is.null(fo_ft) && !is.null(.parse_dynare_date(fo_ft)))
      sample_start <- fo_ft
    parse_filter_tunes_block(ft_block$body, sample_start = sample_start)
  } else {
    empty <- data.frame(var = character(0), stringsAsFactors = FALSE)
    empty$periods <- list()
    empty$values  <- list()
    empty$stderr  <- list()
    list(tunes = empty)
  }

  # ---- 8c. heteroskedastic_shocks block -----------------------------------
  hs_block <- extract_paired_block(txt, "heteroskedastic_shocks")
  heteroskedastic_shocks_parsed <- if (hs_block$found) {
    # Dynare 7 (dynare_estimation_init.m): an integer period p is sample
    # period p - first_obs + 1; a date d is d - date(first_obs) + 1.  So an
    # integer first_obs offsets integer periods, and a date first_obs is the
    # date index that resolves dated periods.
    sample_start_hs <- NULL
    first_obs_hs    <- 1L
    fo_hs <- .estimation_first_obs(txt)
    if (!is.null(fo_hs)) {
      if (!is.null(.parse_dynare_date(fo_hs)) && !grepl("^[0-9]+$", fo_hs)) {
        sample_start_hs <- fo_hs
        first_obs_hs    <- NA_integer_
      } else if (grepl("^[0-9]+$", fo_hs) && as.integer(fo_hs) >= 1L) {
        first_obs_hs <- as.integer(fo_hs)
      } else {
        .dynhr_abort(
          "parse_mod: heteroskedastic_shocks needs estimation(first_obs=) ",
          "as a positive integer or a date; got `", fo_hs, "`.",
          class = "dynhr_error_mod_syntax")
      }
    }
    parse_heteroskedastic_shocks_block(hs_block$body,
                                       sample_start = sample_start_hs,
                                       first_obs = first_obs_hs)
  } else {
    empty_hs <- data.frame(var = character(0), stringsAsFactors = FALSE)
    empty_hs$periods <- list()
    empty_hs$scales  <- list()
    list(scales = empty_hs)
  }

  # ---- 8d. stochastic_volatility block ------------------------------------
  sv_block <- extract_paired_block(txt, "stochastic_volatility")
  stochastic_volatility_parsed <- if (sv_block$found) {
    out_sv <- parse_stochastic_volatility_block(sv_block$body)
    class(out_sv) <- "sv_spec"
    out_sv
  } else {
    empty_sv <- data.frame(shock = character(0), stringsAsFactors = FALSE)
    empty_sv$mu <- list(); empty_sv$rho <- list(); empty_sv$sigma_eta <- list()
    out_sv <- list(sv = empty_sv)
    class(out_sv) <- "sv_spec"
    out_sv
  }

  # ---- 8e. shock_groups block(s) -------------------------------------------
  # Dynare allows several named shock_groups blocks; keep them all keyed by
  # their name= option and expose the first as m$shock_groups, which is the
  # fallback historical_decomposition() reads when no shock_groups argument
  # is supplied.  Membership is validated against the declared shocks.
  shock_groups_blocks <- parse_shock_groups(
    txt, shock_names = c(varexo_names, varexo_det_names))
  shock_groups <- if (length(shock_groups_blocks) > 0L)
    shock_groups_blocks[[1L]] else list()

  # ---- 8f. observation_trends block ----------------------------------------
  # Dynare (dsge_likelihood.m / compute_trend_coefficients.m): the measurement
  # constant of observable i in sample period t is
  #     ys_i + trend_coeff_i * (first_obs + t - 1),   t = 1..T,
  # i.e. the trend index of the FIRST data row is estimation's first_obs
  # (default 1). The slopes are parameter expressions, kept as text.
  #
  # deterministic_trends (the nonstationary-model counterpart, keyed by
  # endogenous variable): kept in model$deterministic_trends, and every entry
  # for an observed variable (every entry when there is no varobs) is a trend
  # of that observable exactly as if it were in observation_trends.  Dynare
  # 7.1's preprocessor stores the block in
  # options_.deterministic_trend_coeffs, which no Dynare routine reads (its
  # driver even fails on the statement, `M_.endogenous_names`); dynhr gives
  # it the observation_trends meaning.  A variable in both blocks is an error.
  ot_block <- extract_paired_block(txt, "observation_trends")
  dt_block <- extract_paired_block(txt, "deterministic_trends")
  observation_trends <- NULL
  deterministic_trends <- NULL
  if (ot_block$found || dt_block$found) {
    .ot <- if (ot_block$found)
      parse_observation_trends_block(ot_block$body, var_names,
                                     param_names, varobs_names)
    else character(0)
    if (dt_block$found) {
      .dt <- parse_observation_trends_block(dt_block$body, var_names,
                                            param_names, character(0),
                                            block = "deterministic_trends")
      .both <- intersect(names(.dt), names(.ot))
      if (length(.both) > 0L)
        .dynhr_abort(
          "parse_mod: ", paste(.both, collapse = ", "), " has a trend in ",
          "both observation_trends and deterministic_trends; give it in one.",
          class = "dynhr_error_mod_syntax")
      if (length(.dt) > 0L) deterministic_trends <- .dt
      .dt_obs <- if (length(varobs_names) > 0L)
        .dt[names(.dt) %in% varobs_names] else .dt
      .ot <- c(.ot, .dt_obs)
    }
    .fo <- 1L
    .fo_ot <- .estimation_first_obs(txt)
    if (!is.null(.fo_ot)) {
      .fo_v <- suppressWarnings(as.numeric(.fo_ot))
      if (is.na(.fo_v) || .fo_v < 1 || .fo_v != round(.fo_v))
        .dynhr_abort(
          "parse_mod: ", if (ot_block$found) "observation_trends"
          else "deterministic_trends", " needs the trend index of the first ",
          "observation, which Dynare takes from estimation(first_obs=...); ",
          "`first_obs = ", .fo_ot, "` is not a positive integer period ",
          "index (a date cannot be resolved without the data file). Write ",
          "first_obs as an integer.", class = "dynhr_error_mod_syntax")
      .fo <- as.integer(.fo_v)
    }
    if (length(.ot) > 0L)
      observation_trends <- list(trends = .ot, first_obs = .fo)
  }

  # ---- 8g. optimal simple rules: osr_params / optim_weights / bounds ------
  # Dynare's OSR inputs (M_.osr.param_names, M_.osr.variable_weights,
  # M_.osr.param_bounds). Weights and bounds are parameter expressions, kept
  # as text and evaluated by osr() at the parameter vector in force.
  .osr_par <- parse_declaration_names(extract_declaration(txt_decl, "osr_params"))
  .osr_bad <- setdiff(.osr_par, param_names)
  if (length(.osr_bad) > 0L)
    .dynhr_abort("parse_mod: osr_params: ", paste(.osr_bad, collapse = ", "),
                 if (length(.osr_bad) == 1L) " is not a declared parameter."
                 else " are not declared parameters.",
                 class = "dynhr_error_mod_syntax")
  .ow_block <- extract_paired_block(txt, "optim_weights")
  .ob_block <- extract_paired_block(txt, "osr_params_bounds")
  osr_spec <- NULL
  if (length(.osr_par) > 0L || .ow_block$found || .ob_block$found) {
    osr_spec <- list(
      params  = .osr_par,
      weights = parse_optim_weights_block(
        if (.ow_block$found) .ow_block$body else "", var_names, param_names),
      bounds  = parse_osr_params_bounds_block(
        if (.ob_block$found) .ob_block$body else "", param_names, .osr_par))
    if (verbose)
      .dynhr_cat("OSR:", length(.osr_par), "osr_params,",
                 nrow(osr_spec$weights), "optim_weights entries,",
                 nrow(osr_spec$bounds), "bounds\n")
  }

  # ---- 9. Estimated params ---------------------------------------------
  ep_block <- extract_paired_block(txt, "estimated_params")
  estimated_params <- if (ep_block$found)
    parse_estimated_params_block(ep_block$body)
  else data.frame(type = character(0),
                  name = character(0),
                  name2 = character(0),
                  prior = character(0),
                  p1 = numeric(0), p2 = numeric(0),
                  p3 = numeric(0), p4 = numeric(0),
                  init = numeric(0), lb = numeric(0), ub = numeric(0),
                  stringsAsFactors = FALSE)

  # Backfill calibration for parameters that are declared and estimated but
  # never assigned a value in the calibration section (e.g. constepinf,
  # constebeta, ctrend in Smets-Wouters 2007). Dynare uses the estimated_params
  # INITVAL as the starting parameter value in this case.
  if (nrow(estimated_params) > 0 && "init" %in% names(estimated_params)) {
    ep_par <- estimated_params[estimated_params$type == "parameter", , drop = FALSE]
    for (i in seq_len(nrow(ep_par))) {
      nm <- ep_par$name[i]
      iv <- ep_par$init[i]
      if (nm %in% param_names && !(nm %in% names(param_values)) &&
          length(iv) == 1L && is.finite(iv)) {
        param_values[nm] <- iv
      }
    }
  }

  ep_init_block <- extract_paired_block(txt, "estimated_params_init")
  ep_init <- if (ep_init_block$found) {
    use_cal <- grepl("use_calibration", ep_init_block$options_str, ignore.case = TRUE)
    df <- parse_estimated_params_block(ep_init_block$body)
    # When use_calibration is set, seed any estimated param that has a
    # calibrated value but no explicit init row.  This ensures the optimizer
    # starts from the calibration rather than from prior means.
    if (use_cal && nrow(estimated_params) > 0 && length(param_values) > 0) {
      ep_par_names <- estimated_params$name[estimated_params$type == "parameter"]
      for (pn in ep_par_names) {
        already_init <- nrow(df) > 0 && "name" %in% names(df) && pn %in% df$name
        if (!already_init && pn %in% names(param_values)) {
          new_row <- data.frame(
            type = "parameter", name = pn, name2 = NA_character_,
            prior = NA_character_,
            p1 = NA_real_, p2 = NA_real_, p3 = NA_real_, p4 = NA_real_,
            init = unname(param_values[pn]),
            lb = NA_real_, ub = NA_real_,
            stringsAsFactors = FALSE
          )
          df <- if (nrow(df) == 0) new_row else rbind(df, new_row)
        }
      }
      # Also backfill stderr estimated params from calibrated shock stds
      ep_std_names <- estimated_params$name[estimated_params$type == "stderr"]
      for (sn in ep_std_names) {
        already_init <- nrow(df) > 0 && "name" %in% names(df) && sn %in% df$name
        if (!already_init && sn %in% names(param_values)) {
          new_row <- data.frame(
            type = "stderr", name = sn, name2 = NA_character_,
            prior = NA_character_,
            p1 = NA_real_, p2 = NA_real_, p3 = NA_real_, p4 = NA_real_,
            init = unname(param_values[sn]),
            lb = NA_real_, ub = NA_real_,
            stringsAsFactors = FALSE
          )
          df <- if (nrow(df) == 0) new_row else rbind(df, new_row)
        }
      }
    }
    attr(df, "use_calibration") <- use_cal
    df
  } else data.frame()

  # ---- 9b. occbin_constraints block ------------------------------------
  # Use nested extractor: the block may contain equations;...end; sub-blocks
  # that would prematurely close a non-greedy match in extract_paired_block.
  occbin_block <- extract_paired_block_nested(txt, "occbin_constraints", "equations")
  occbin_constraints <- if (occbin_block$found)
    parse_occbin_constraints_block(occbin_block$body)
  else list()
  if (verbose && length(occbin_constraints) > 0)
    .dynhr_cat("occbin_constraints:", length(occbin_constraints), "named constraints\n")

  # ---- 10. Commands ----------------------------------------------------
  commands <- list()
  for (cmd in c("stoch_simul", "estimation", "steady", "check",
                "model_diagnostics", "model_info", "forecast",
                "calib_smoother", "shock_decomposition",
                "conditional_forecast", "osr", "ramsey_model",
                "ramsey_policy", "planner_objective", "discretionary_policy",
                "perfect_foresight_setup", "perfect_foresight_solver",
                "extended_path", "simul")) {
    cmd_info <- extract_command(txt, cmd)
    if (cmd_info$found) {
      cmd_opts <- parse_command_options(cmd_info$options_str)
      cmd_vars <- if (nchar(cmd_info$var_list) > 0)
        parse_declaration_names(cmd_info$var_list)
      else character(0)
      commands <- c(commands, list(list(
        name     = cmd,
        options  = cmd_opts,
        var_list = cmd_vars
      )))
      if (verbose) .dynhr_cat("Found command:", cmd, "\n")
    }
  }

  # ---- 10a. shock_paths (Dynare 7) --------------------------------------
  # Resolved into the SAME fields as the `shocks` block's periods/values path
  # (det_shocks) and the endval block (terminal exogenous values), see
  # parse_shock_paths_blocks().  Dynare's manual: shock_paths "cannot be
  # used in conjunction with" shocks periods / mshocks / endval.
  shock_paths <- NULL
  controlled_paths <- NULL
  .ctl <- .dynhr_empty_controlled()
  .sp_blocks <- .extract_shock_paths_blocks(txt)
  .cp_blocks <- .extract_controlled_paths_blocks(txt)
  .pf <- list()
  if (length(.sp_blocks) > 0L || length(.cp_blocks) > 0L) {
    for (.cn in c("perfect_foresight_setup",
                  "perfect_foresight_with_expectation_errors_setup")) {
      .ci <- extract_command(txt, .cn)
      if (.ci$found) .pf <- parse_command_options(.ci$options_str)
    }
  }
  if (length(.sp_blocks) > 0L) {
    if (endval_block$found || nrow(det_shocks) > 0L ||
        extract_paired_block(txt, "mshocks")$found)
      .dynhr_abort(
        "parse_mod: a shock_paths block cannot be combined with an endval ",
        "block or with deterministic shocks in a shocks/mshocks block ",
        "(Dynare 7 reference manual, shock_paths); describe the whole ",
        "scenario in shock_paths.", class = "dynhr_error_mod_syntax")
    .first_sim <- if (!is.null(.pf$first_simulation_period))
      as.character(.pf$first_simulation_period) else NULL
    .n_per <- if (is.numeric(.pf$periods)) as.integer(.pf$periods)
      else if (!is.null(.first_sim) && !is.null(.pf$last_simulation_period))
        .date_to_period(as.character(.pf$last_simulation_period), .first_sim,
                        context = "shock_paths")
      else NA_integer_
    .sp <- parse_shock_paths_blocks(
      .sp_blocks, exo_names = c(varexo_names, varexo_det_names),
      endo_names = var_names, param_values = param_values,
      initval = initval,
      has_steady = any(vapply(commands, function(cm) identical(cm$name, "steady"),
                              logical(1))),
      n_periods = .n_per, first_sim = .first_sim)
    det_shocks  <- .sp$det
    endval      <- .sp$terminal
    .ctl        <- .sp$controlled
    shock_paths <- list(n_periods = .n_per,
                        first_simulation_period = .first_sim %||% NA_character_)
  }

  # ---- 10a'. perfect_foresight_controlled_paths (Dynare 7) ----------------
  # Together with the exogenize/endogenize stanzas of shock_paths: one row
  # per (period, exogenized endogenous, endogenized shock), consumed by
  # perfect_foresight_solve(controlled_paths =).  NULL when absent.
  if (length(.cp_blocks) > 0L) {
    .ctl <- rbind(.ctl, parse_controlled_paths_blocks(
      .cp_blocks, endo_names = var_names,
      exo_names = c(varexo_names, varexo_det_names),
      param_values = param_values,
      first_sim = if (!is.null(.pf$first_simulation_period))
        as.character(.pf$first_simulation_period) else NULL))
  }
  if (nrow(.ctl) > 0L) {
    .n_per_c <- if (is.numeric(.pf$periods)) as.integer(.pf$periods)
                else NA_integer_
    controlled_paths <- .dynhr_check_controlled(.ctl, .n_per_c, det_shocks)
  }

  # ---- 10b. planner_objective -----------------------------------------
  po <- extract_planner_objective(txt)
  planner_objective <- list(text = "", ast = NULL)
  if (isTRUE(po$found) && nzchar(po$text)) {
    planner_objective$text <- po$text
    planner_objective$ast <- parse_expression(
      po$text,
      var_names = all_var_names,
      param_names = param_names
    )
    if (verbose) .dynhr_cat("Found command: planner_objective\n")
  }

  # ---- 10c. ramsey_policy instruments= --------------------------------
  ramsey_instruments <- extract_ramsey_instruments(txt)

  # ---- 10c'. ramsey_constraints block -----------------------------------
  # Bounds on Ramsey variables (`i > 0;`), complementary to the Ramsey FOC
  # with respect to the variable; honoured by ramsey_obc_pf(), as Dynare
  # honours them only in its mixed-complementarity perfect-foresight solver.
  .rc_block <- extract_paired_block(txt, "ramsey_constraints")
  ramsey_constraints <- if (.rc_block$found)
    parse_ramsey_constraints_block(.rc_block$body, var_names, param_names)
  else NULL

  # ---- 10d. metadata (@dynhr: blocks) ---------------------------------
  ## No tryCatch and no temp file. `extract_mod_metadata()` takes a path
  ## OR lines directly now (`.mod_lines()`), so the text branch no longer needs
  ## a round trip through the filesystem -- and a genuine metadata parse error
  ## SURFACES instead of silently yielding `list()`, which is what hid
  ## `.extract_narratives()`'s readLines() failure for the whole of 0.9.3.
  metadata <- if (!is.na(source_file)) extract_mod_metadata(source_file)
              else extract_mod_metadata(strsplit(file_or_text, "\n")[[1]])

  # ---- 11. (predetermined re-timing) ------------------------------------
  # Dynare's `predetermined_variables` convention (a declared variable is a
  # BEGINNING-of-period stock: `k` means standard `k(-1)`, `k(+1)` means `k`,
  # ...) is applied in step 4's single ordered timing pass, on the model-block
  # text BEFORE auxiliary expansion and parsing -- see the D3 comment there.
  # It used to be applied here, on the parsed equation ASTs, which (a) missed
  # the `#` model-local definitions and (b) ran after the aux expansion, so a
  # predetermined `k(-1)` (standard k(-2)) never got its AUX_LAG chain.

  # ---- 12. Variable classification -------------------------------------
  vc <- build_variable_classification(equations, var_names, varexo_names,
                                      predetermined_vars = predet_vars,
                                      local_vars = local_vars)

  # ---- 12. Assemble and return -----------------------------------------
  model <- new_dynhr_mod(
    source_file          = source_file,
    var_names            = var_names,
    varexo_names         = varexo_names,
    varexo_det_names     = varexo_det_names,
    param_names          = param_names,
    predetermined_vars   = predet_vars,
    param_values         = param_values,
    equations            = equations,
    local_variables      = local_vars,
    model_options        = model_opts,
    initval              = initval,
    endval               = endval,
    histval              = histval,
    steady_state_model   = ssm,
    shocks               = shocks,
    shocks_blocks        = shocks_blocks_parsed,
    shock_groups         = shock_groups,
    shock_groups_blocks  = shock_groups_blocks,
    det_shocks           = det_shocks,
    filter_tunes         = filter_tunes,
    heteroskedastic_shocks = heteroskedastic_shocks_parsed,
    stochastic_volatility = stochastic_volatility_parsed,
    estimated_params     = estimated_params,
    estimated_params_init = ep_init,
    commands             = commands,
    planner_objective    = planner_objective,
    occbin_constraints   = occbin_constraints,
    lead_lag_incidence   = vc$lead_lag_incidence,
    variable_classification = vc$classification,
    n_static             = vc$n_static,
    n_predetermined      = vc$n_predetermined,
    n_forward            = vc$n_forward,
    n_mixed              = vc$n_mixed,
    varobs               = varobs_names,
    varobs_names         = varobs_names,
    ## alias: estimation entry points take an `obs_vars` argument -- populate
    ## the same-named model field from the .mod's varobs line so callers can
    ## omit it
    obs_vars             = varobs_names,
    ramsey_instruments   = ramsey_instruments,
    ramsey_constraints   = ramsey_constraints,
    metadata             = metadata,
    ## Present only when the .mod uses the construct (NULL leaves the field
    ## out), so models without them keep exactly their old field set.
    observation_trends   = observation_trends,
    deterministic_trends = deterministic_trends,
    nonstationary        = nonstationary,
    osr                  = osr_spec,
    log_vars             = if (length(log_map) > 0L) log_map else NULL,
    shock_paths          = shock_paths,
    controlled_paths     = controlled_paths
  )

  # ---- 13. MCP constraint tags ----------------------------------
  # mcp_parse_tags() reads equation $tag_raw fields and builds MCP spec
  # objects.  Wiring it here mirrors the M6/M13 metadata-wiring pattern
  # (same parse-wiring family).  Errors are caught so a tag-syntax problem
  # does not crash parse_mod() for the whole model.
  model$mcp_constraints <- tryCatch(
    mcp_parse_tags(model, verbose = FALSE),
    error = function(e) {
      .dynhr_warn("parse_mod(): mcp_parse_tags() failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
  if (verbose && length(model$mcp_constraints) > 0)
    .dynhr_cat("MCP constraints:", length(model$mcp_constraints), "\n")

  # ---- 14. Undeclared parameter detection (M24b) ----------------------
  # Collect every parameter-typed AST node name from equations and the
  # planner objective.  Any that are NOT in the declared param_names (and
  # are not known math functions or model-local #-defines) are parameters
  # computed externally (e.g. in *_steadystate.m) and absent from
  # m$param_values.  Warn clearly; do NOT error (legitimate pattern).
  .ast_collect_param_names <- function(node) {
    if (is.null(node)) return(character(0))
    switch(node$type,
      "parameter"      = node$name,
      "variable"       = character(0),
      "number"         = character(0),
      "local_variable" = character(0),
      "binop"          = c(.ast_collect_param_names(node$left),
                           .ast_collect_param_names(node$right)),
      "unaryop"        = .ast_collect_param_names(node$operand),
      "funcall"        = {
        # Function name itself is not a parameter; recurse into args only
        unlist(lapply(node$args, .ast_collect_param_names), use.names = FALSE)
      },
      character(0)
    )
  }

  .eq_param_names <- unique(unlist(lapply(model$equations, function(eq) {
    c(.ast_collect_param_names(eq$lhs), .ast_collect_param_names(eq$rhs))
  }), use.names = FALSE))

  .po_param_names <- if (!is.null(model$planner_objective$ast))
    .ast_collect_param_names(model$planner_objective$ast)
  else character(0)

  .all_eq_params <- unique(c(.eq_param_names, .po_param_names))

  # Exclude: declared params, known math functions, model-local #-defines
  .local_define_names <- names(model$local_variables)
  .undeclared <- setdiff(
    .all_eq_params,
    c(param_names, .KNOWN_FUNCTIONS, .local_define_names)
  )

  if (length(.undeclared) > 0L) {
    .dynhr_warn(
      "Undeclared parameters referenced in equations: ",
      paste(.undeclared, collapse = ", "),
      ".\nDeclare them in `parameters`, or supply via inject_params() (for an ",
      "external/superset vector) or set_param_values() (for the exact model-param set).",
      call. = FALSE
    )
  }

  # Declared parameters that actually feed the dynamic Jacobian (referenced in
  # an equation, the planner objective, or a model-local #-define that is later
  # inlined).  Used by solve_perturbation's I2-extended fail-loud guard so it
  # only fires on params whose absence would corrupt the decision rule -- NOT on
  # declared-but-equation-absent params such as shock standard deviations
  # (`stderr <p>` in the shocks block), which never enter the Jacobian.
  .localdef_param_names <- unique(unlist(
    lapply(model$local_variables, .ast_collect_param_names), use.names = FALSE))
  model$equation_param_names <- intersect(
    unique(c(.all_eq_params, .localdef_param_names)), param_names)

  if (verbose) .dynhr_cat("Done. dynhr_mod object created.\n")
  model
}


# ===========================================================================
# Model-block timing helpers
#
# They rewrite the model-block TEXT, so equations and `#` model-local
# definitions are treated alike (the locals are part of the same text).  The
# order in which parse_mod() applies them is documented at its step 4.
# ===========================================================================

## One token of model-block text.  ORDER matters: an equation tag `[...]` and a
## quoted string are consumed whole and never rewritten; a number literal is
## consumed whole so the `e3` of `1e3` is never read as an identifier; only
## then identifiers.
.MOD_TIMING_TOKEN_RE <- paste0(
  "\\[[^\\]]*\\]",
  "|'[^']*'|\"[^\"]*\"",
  "|(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?",
  "|[A-Za-z_][A-Za-z0-9_]*")

## `x(t)` timing suffix in Dynare syntax: "" for t = 0.
.mod_timing_suffix <- function(t) {
  if (t == 0L) "" else if (t > 0L) sprintf("(+%d)", t) else sprintf("(%d)", t)
}

## Position of the `)` that closes the `(` at `open_pos` of `text`, or NA.
.mod_matching_paren <- function(text, open_pos) {
  ch <- strsplit(substr(text, open_pos, nchar(text)), "", fixed = TRUE)[[1]]
  depth <- cumsum((ch == "(") - (ch == ")"))
  hit <- which(depth == 0L)[1L]
  if (is.na(hit)) NA_integer_ else open_pos + hit - 1L
}

## Shift the timing of NAMED variables in model-block text.
##
## `shifts` is a named integer vector (variable -> offset).  Every occurrence
## of such a variable -- bare `x` (timing 0) or with an explicit integer
## timing `x(+2)` / `x( - 1 )` -- is rewritten to its new timing in canonical
## form (`x`, `x(+t)`, `x(-t)`).  Left untouched: other identifiers,
## equation tags, quoted strings, number literals, and the whole argument of
## `STEADY_STATE(...)` / `steady_state(...)` (a steady-state value has no
## timing).
.mod_retime_text <- function(text, shifts) {
  shifts <- shifts[shifts != 0L]
  if (length(shifts) == 0L || !nzchar(text)) return(text)
  m <- gregexpr(.MOD_TIMING_TOKEN_RE, text, perl = TRUE)[[1]]
  if (m[1L] == -1L) return(text)
  starts <- as.integer(m)
  ends   <- starts + attr(m, "match.length") - 1L
  toks   <- substring(text, starts, ends)
  ss_fns <- c("STEADY_STATE", "steady_state")
  hits   <- which(toks %in% c(names(shifts), ss_fns))
  if (length(hits) == 0L) return(text)

  n <- nchar(text)
  pieces  <- character(0)
  cursor  <- 1L
  skip_to <- 0L
  for (j in hits) {
    if (starts[j] <= skip_to) next
    after <- ends[j] + 1L
    look  <- substr(text, after, min(n, after + 200L))
    if (toks[j] %in% ss_fns) {
      op <- regexpr("^\\s*\\(", look, perl = TRUE)
      if (op > 0L) {
        close <- .mod_matching_paren(text, after + attr(op, "match.length") - 1L)
        if (!is.na(close)) skip_to <- close
      }
      next
    }
    tm <- regmatches(look, regexec("^\\s*\\(\\s*([+-]?)\\s*(\\d+)\\s*\\)",
                                   look, perl = TRUE))[[1]]
    if (length(tm) == 3L) {
      old      <- as.integer(tm[3L]) * (if (tm[2L] == "-") -1L else 1L)
      consumed <- nchar(tm[1L])
    } else {
      old      <- 0L
      consumed <- 0L
    }
    pieces <- c(pieces, substr(text, cursor, starts[j] - 1L), toks[j],
                .mod_timing_suffix(old + shifts[[toks[j]]]))
    cursor <- after + consumed
  }
  paste(c(pieces, substr(text, cursor, n)), collapse = "")
}

## Dynare's EXPECTATION operator on model-block text.
##
## `EXPECTATION(k)(expr)` is the expectation of `expr` conditional on the
## information of period t+k (k is a signed integer; k = -1 is "as of t-1").
## As in the Dynare preprocessor (UnaryOpNode::substituteExpectation, reached
## from DynamicModel::substituteExpectation), each distinct operator becomes
## an auxiliary ENDOGENOUS variable AUX defined by the new equation
##     AUX = expr shifted by -k   (every variable's timing moved by -k)
## and the operator is replaced by AUX(k).  E.g. EXPECTATION(-1)(y(+1)) ->
## AUX(-1) with AUX = y(+2): AUX_t = E_t y_{t+2}, so AUX_{t-1} = E_{t-1}
## y_{t+1}.  Names follow Dynare's SymbolTable::addExpectationAuxiliaryVar,
## AUX_EXPECT_{LAG|LEAD}_<|k|>_<index> (LAG for k < 0).  Identical operators
## share one auxiliary variable.  Checked against the Dynare 7.1
## preprocessor's `json=transform` output for k = -2, -1, 0, +1.
##
## Not implemented, each an error of class `dynhr_error_mod_expectation`
## (never an approximation): a nested EXPECTATION, and an argument that uses
## a `#` model-local variable (Dynare inlines the local first).  A malformed
## operator is a `dynhr_error_mod_syntax` error.
##
## Returns list(model_body, aux_var_names).
.substitute_expectation_text <- function(body, endo_names,
                                         exo_names = character(0),
                                         reserved_names = character(0)) {
  none <- list(model_body = body, aux_var_names = character(0))
  if (!grepl("\\bEXPECTATION\\b", body, perl = TRUE)) return(none)

  unsupported <- function(...)
    .dynhr_abort("parse_mod: ", ..., class = c("dynhr_error_mod_expectation",
                                               "dynhr_error_mod_syntax"))
  # `#` model-local names (statement start, optionally after a tag).
  loc_m <- regmatches(body, gregexpr(
    "(?:^|;)\\s*(?:\\[[^\\]]*\\]\\s*)?#\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=",
    body, perl = TRUE))[[1]]
  local_names <- sub("^.*#\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=$", "\\1", loc_m,
                     perl = TRUE)
  var_all <- unique(c(endo_names, exo_names))

  aux_names <- character(0)
  aux_eqs   <- character(0)
  aux_keys  <- character(0)
  used      <- c(reserved_names, endo_names, exo_names)
  index     <- 0L
  repeat {
    pos <- regexpr("\\bEXPECTATION\\b", body, perl = TRUE)
    if (pos < 0L) break
    rest <- substr(body, pos, nchar(body))
    hd <- regmatches(rest, regexec(
      "^EXPECTATION\\s*\\(\\s*([+-]?)\\s*(\\d+)\\s*\\)\\s*\\(", rest,
      perl = TRUE))[[1]]
    if (length(hd) != 3L)
      .dynhr_abort(
        "parse_mod: malformed EXPECTATION operator near `",
        substr(gsub("\\s+", " ", rest), 1L, 60L), "`: the form is ",
        "`EXPECTATION(k)(expression)` with k an integer, e.g. ",
        "EXPECTATION(-1)(x(+1)).", class = "dynhr_error_mod_syntax")
    k <- as.integer(hd[3L]) * (if (hd[2L] == "-") -1L else 1L)
    open_pos  <- pos + nchar(hd[1L]) - 1L
    close_pos <- .mod_matching_paren(body, open_pos)
    if (is.na(close_pos))
      .dynhr_abort("parse_mod: unbalanced parentheses in the argument of ",
                   "EXPECTATION(", k, ").", class = "dynhr_error_mod_syntax")
    arg <- substr(body, open_pos + 1L, close_pos - 1L)
    if (grepl("\\bEXPECTATION\\b", arg, perl = TRUE))
      unsupported("nested EXPECTATION operators are not implemented (in `",
                  "EXPECTATION(", k, ")(", trimws(arg), ")`).")
    arg_ids <- regmatches(arg, gregexpr(.MOD_TIMING_TOKEN_RE, arg,
                                        perl = TRUE))[[1]]
    arg_loc <- intersect(arg_ids, local_names)
    if (length(arg_loc) > 0L)
      unsupported("EXPECTATION(", k, ")(", trimws(arg), ") uses the model-",
                  "local variable(s) ", paste(arg_loc, collapse = ", "),
                  "; an EXPECTATION argument that uses a `#` local is not ",
                  "implemented. Write the local's expression out in the ",
                  "argument.")

    key <- paste0(k, "|", gsub("\\s+", "", arg))
    hit <- match(key, aux_keys)
    if (is.na(hit)) {
      repeat {
        index <- index + 1L
        nm <- sprintf("AUX_EXPECT_%s_%d_%d", if (k < 0L) "LAG" else "LEAD",
                      abs(k), index)
        if (!(nm %in% used)) break
      }
      used      <- c(used, nm)
      aux_names <- c(aux_names, nm)
      aux_keys  <- c(aux_keys, key)
      shifted   <- .mod_retime_text(
        arg, stats::setNames(rep(-k, length(var_all)), var_all))
      aux_eqs   <- c(aux_eqs, paste0(nm, " = ", trimws(shifted), ";"))
    } else {
      nm <- aux_names[hit]
    }
    body <- paste0(substr(body, 1L, pos - 1L), nm, .mod_timing_suffix(k),
                   substr(body, close_pos + 1L, nchar(body)))
  }
  list(model_body = paste0(body, "\n", paste(aux_eqs, collapse = "\n"), "\n"),
       aux_var_names = aux_names)
}

## Dynare's var(log) substitution on model-block text (Dynare 6+).
##
## `log_map` maps a declared name x to its log variable LOG_x. Every
## occurrence of x -- bare or with an explicit integer timing, in equations and
## `#` local definitions alike -- becomes `exp(LOG_x)` at the same timing:
## x(+1) -> exp(LOG_x(+1)). Left untouched, as in .mod_retime_text(): other
## identifiers, equation tags, quoted strings, number literals, and the whole
## argument of STEADY_STATE(...) (x itself stays a declared variable, pinned
## by the appended `x = exp(LOG_x)`, so its steady state is unchanged).
## Checked against the Dynare 7.1 preprocessor's `json=transform` output,
## which rewrites `c(+1)` to `exp(LOG_c(1))` and appends `c = exp(LOG_c)`.
.mod_log_substitute_text <- function(text, log_map) {
  repl <- lapply(unname(log_map), function(lv) {
    force(lv)
    function(t) paste0("exp(", lv, .mod_timing_suffix(t), ")")
  })
  names(repl) <- names(log_map)
  .mod_substitute_var_text(text, repl)
}

## var(log): carry a named starting value of x (initval / endval) over to
## LOG_x = log(x), unless LOG_x is set explicitly. A non-positive x has no log
## and is left to the steady-state solver's default guess.
.mod_log_carry <- function(values, log_map) {
  for (v in names(log_map)) {
    lv <- log_map[[v]]
    if (v %in% names(values) && !(lv %in% names(values)) &&
        is.finite(values[[v]]) && values[[v]] > 0)
      values[lv] <- log(values[[v]])
  }
  values
}

## ---------------------------------------------------------------------------
## Dynare 7 complementarity conditions (`EQ U+27C2 COND;`)
## ---------------------------------------------------------------------------
## Dynare 7.2 reference manual, "The model file" > perfect_foresight_solver >
## option `lmmcp`: a complementarity condition follows its equation, separated
## by the perpendicular symbol U+27C2 or the ASCII `_|_`.  Both bounds
## may be given at once, and the bounds may be arbitrary functions of
## parameters:  `mu = 0 U+27C2 0 < i < 1+2*alpha;`.  The old `[mcp = 'i > 0']`
## tag is still accepted (deprecated in Dynare 7).
##
## Forms accepted (checked against the Dynare 7.1 preprocessor's JSON
## `complementarity_condition` = {variable, lower_bound, upper_bound}):
##   x > B   x >= B   B < x   B <= x           -> lower bound B
##   x < B   x <= B   B > x   B >= x           -> upper bound B
##   L < x < U  (or <=)   U > x > L  (or >=)   -> both
## `>` and `>=` mean the same thing (the JSON does not distinguish them).
## Rejected, as by the preprocessor (dynhr_error_mod_syntax): a bound that
## contains an endogenous or exogenous variable (including
## `steady_state(y)`), a constrained side that is not a bare endogenous
## variable (`x(-1) > 0`, `2*x > 0`), and mixed directions (`0 < x > 1`).

## Split a complementarity condition at its top-level relational operators.
## Returns list(parts, ops).
.mod_cc_split <- function(cond) {
  ch <- strsplit(cond, "", fixed = TRUE)[[1]]
  n  <- length(ch)
  depth <- 0L
  parts <- character(0); ops <- character(0)
  cur <- ""
  k <- 1L
  while (k <= n) {
    c1 <- ch[k]
    if (c1 == "(") depth <- depth + 1L
    if (c1 == ")") depth <- depth - 1L
    if (depth == 0L && c1 %in% c("<", ">")) {
      op <- c1
      if (k < n && ch[k + 1L] == "=") { op <- paste0(c1, "="); k <- k + 1L }
      parts <- c(parts, trimws(cur)); ops <- c(ops, op); cur <- ""
    } else {
      cur <- paste0(cur, c1)
    }
    k <- k + 1L
  }
  list(parts = c(parts, trimws(cur)), ops = ops)
}

## Parse one condition into list(variable, lower, upper) (bounds as text, NA
## when absent).  `label` names the equation for error messages.
.mod_parse_complementarity <- function(cond, label, endo_names, exo_names,
                                       param_names) {
  bad <- function(why)
    .dynhr_abort("parse_mod: complementarity condition `", cond, "` on ",
                 label, " has an incorrect form", why, ".",
                 class = "dynhr_error_mod_syntax")
  sp <- .mod_cc_split(cond)
  parts <- sp$parts
  dirs  <- substr(sp$ops, 1L, 1L)
  is_var <- function(s) grepl("^[A-Za-z_][A-Za-z0-9_]*$", s) &&
    s %in% endo_names
  lower <- NA_character_
  upper <- NA_character_
  if (length(parts) == 2L && all(nzchar(parts))) {
    if (is_var(parts[1L])) {
      v <- parts[1L]; b <- parts[2L]
      if (dirs == ">") lower <- b else upper <- b
    } else if (is_var(parts[2L])) {
      v <- parts[2L]; b <- parts[1L]
      if (dirs == "<") lower <- b else upper <- b
    } else {
      bad(": one side must be a bare endogenous variable")
    }
  } else if (length(parts) == 3L && all(nzchar(parts))) {
    if (!is_var(parts[2L]))
      bad(": the middle term must be a bare endogenous variable")
    if (dirs[1L] != dirs[2L])
      bad(": both inequalities must point the same way")
    v <- parts[2L]
    if (dirs[1L] == "<") { lower <- parts[1L]; upper <- parts[3L] }
    else                 { upper <- parts[1L]; lower <- parts[3L] }
  } else {
    bad(": expected `x > L`, `x < U` or `L < x < U`")
  }
  for (b in c(lower, upper)) {
    if (is.na(b)) next
    ast  <- parse_expression(b, c(endo_names, exo_names), param_names,
                             context = paste0("the complementarity bound `",
                                              b, "` on ", label))
    used <- ast_collect_variables(ast)
    if (!is.null(used) && nrow(used) > 0L)
      bad(paste0(": bounds must not contain any endogenous or exogenous ",
                 "variable (`", b, "` uses ",
                 paste(unique(used$name), collapse = ", "), ")"))
  }
  list(variable = v, lower = lower, upper = upper)
}

## Cut every `U+27C2 COND` / `_|_ COND` off the model-block body.
##
## Returns list(body, conditions): `body` is the input with each condition
## (separator included) removed and everything else untouched; `conditions`
## is a list named by the equation's ordinal among the block's equations --
## counted exactly as parse_model_block() counts them (one tag stripped; `#`
## locals, empty statements and `end` skipped).
.mod_extract_complementarity <- function(body, endo_names, exo_names,
                                         param_names) {
  sep_re <- "\u27c2|_\\|_"
  if (!grepl(sep_re, body, perl = TRUE))
    return(list(body = body, conditions = list()))
  pieces <- strsplit(body, ";", fixed = TRUE)[[1]]
  conds  <- list()
  n_eq   <- 0L
  for (p in seq_along(pieces)) {
    raw  <- pieces[p]
    stmt <- trimws(raw)
    if (!nzchar(stmt) || grepl("^end$", stmt, ignore.case = TRUE)) next
    ## blank the first equation tag (parse_model_block() strips exactly one)
    tg <- regexpr("\\[[^\\]]+\\]", raw, perl = TRUE)
    masked <- raw
    if (tg > 0L)
      substr(masked, tg, tg + attr(tg, "match.length") - 1L) <-
        strrep(" ", attr(tg, "match.length"))
    rest <- trimws(if (tg > 0L) sub("\\[[^\\]]+\\]", " ", raw, perl = TRUE)
                   else raw)
    if (!nzchar(rest)) next
    is_local <- startsWith(rest, "#")
    if (!is_local) n_eq <- n_eq + 1L
    hits <- gregexpr(sep_re, masked, perl = TRUE)[[1]]
    if (hits[1L] == -1L) next
    label <- sprintf("model equation %d (`%s`)", n_eq,
                     gsub("\\s+", " ", rest))
    if (is_local)
      .dynhr_abort("parse_mod: a complementarity condition cannot follow a ",
                   "`#` model-local definition (`", gsub("\\s+", " ", rest),
                   "`).", class = "dynhr_error_mod_syntax")
    if (length(hits) > 1L)
      .dynhr_abort("parse_mod: ", label, " has more than one complementarity ",
                   "separator (`\u27c2` / `_|_`).",
                   class = "dynhr_error_mod_syntax")
    if (tg > 0L && grepl("\\bmcp\\s*=", regmatches(raw, tg), perl = TRUE))
      .dynhr_abort("parse_mod: ", label, " has both an `mcp` tag and a ",
                   "`\u27c2` complementarity condition; give only one.",
                   class = "dynhr_error_mod_syntax")
    cond <- trimws(substr(raw, hits[1L] + attr(hits, "match.length")[1L],
                          nchar(raw)))
    conds[[as.character(n_eq)]] <- .mod_parse_complementarity(
      cond, label, endo_names, exo_names, param_names)
    pieces[p] <- substr(raw, 1L, hits[1L] - 1L)
  }
  ## strsplit() drops a trailing empty piece, so restore a final `;`.
  out <- paste(pieces, collapse = ";")
  if (endsWith(body, ";")) out <- paste0(out, ";")
  list(body = out, conditions = conds)
}

## ---------------------------------------------------------------------------
## Model-local variables with a lead or lag (Dynare 7.0)
## ---------------------------------------------------------------------------
## Dynare 7.2 reference manual, "The model file" > model block: "if the
## model-local variable appears with a lead or a lag attached to it between
## parenthesis, the substitution will be done by shifting the expression
## accordingly".  So with `# r_real = i - pi(+1) + u;`, `r_real(-1)` is
## `(i(-1) - pi + u(-1))`.  Checked against the Dynare 7.1 preprocessor's
## `json=transform` output, including a local defined through another local
## (`# gap = r_real(-1) - y(-1)*beta;` used as `gap(+1)` and `gap(-2)`), an
## exogenous variable inside the definition, and a predetermined variable.
##
## Every `L(k)` (k != 0) of a local L is replaced by `(DEF_L shifted by k)`,
## where the shift moves every endogenous/exogenous variable AND every other
## local by k (`.mod_retime_text`); a shifted reference to another local is
## resolved on the next pass.  `L(0)` is `L`.  A bare `L` is left alone (it
## stays a model-local).  Locals are acyclic in Dynare (defined before use),
## so the passes terminate; a cycle aborts.
##
## @param body        Model-block text.
## @param shift_names Names whose timing moves with the local (declared
##   endogenous + exogenous names).
.mod_substitute_timed_locals <- function(body, shift_names) {
  loc_re <- paste0("(?s)^\\s*(?:\\[[^\\]]*\\]\\s*)?#\\s*",
                   "([A-Za-z_][A-Za-z0-9_]*)\\s*=(?!=)(.*)$")
  defs_of <- function(txt) {
    pcs <- strsplit(txt, ";", fixed = TRUE)[[1]]
    hit <- regmatches(pcs, regexec(loc_re, pcs, perl = TRUE))
    hit <- hit[lengths(hit) == 3L]
    stats::setNames(vapply(hit, `[`, "", 3L), vapply(hit, `[`, "", 2L))
  }
  defs <- defs_of(body)
  if (length(defs) == 0L) return(body)
  lnames <- names(defs)
  all_shift <- c(shift_names, lnames)
  for (pass in seq_len(length(lnames) + 2L)) {
    m <- gregexpr(.MOD_TIMING_TOKEN_RE, body, perl = TRUE)[[1]]
    starts <- as.integer(m)
    ends   <- starts + attr(m, "match.length") - 1L
    toks   <- substring(body, starts, ends)
    hits   <- which(toks %in% lnames)
    n <- nchar(body)
    pieces  <- character(0)
    cursor  <- 1L
    changed <- FALSE
    for (j in hits) {
      after <- ends[j] + 1L
      look  <- substr(body, after, min(n, after + 200L))
      tm <- regmatches(look, regexec("^\\s*\\(\\s*([+-]?)\\s*(\\d+)\\s*\\)",
                                     look, perl = TRUE))[[1]]
      if (length(tm) != 3L) next
      k <- as.integer(tm[3L]) * (if (tm[2L] == "-") -1L else 1L)
      rep_txt <- if (k == 0L) toks[j] else
        paste0("(", trimws(.mod_retime_text(
          defs[[toks[j]]],
          stats::setNames(rep(k, length(all_shift)), all_shift))), ")")
      pieces  <- c(pieces, substr(body, cursor, starts[j] - 1L), rep_txt)
      cursor  <- after + nchar(tm[1L])
      changed <- TRUE
    }
    if (!changed) return(body)
    body <- paste(c(pieces, substr(body, cursor, n)), collapse = "")
    defs <- defs_of(body)
  }
  .dynhr_abort("parse_mod: the `#` model-local variables ",
               paste(lnames, collapse = ", "), " refer to each other in a ",
               "cycle through their leads/lags; Dynare requires a local to be ",
               "defined from earlier locals only.",
               class = "dynhr_error_mod_syntax")
}


# ===========================================================================
# model_replace / model_remove / var_remove (Dynare 6+)
# ===========================================================================

## Parse a comma-separated tag list: `key='value'` items, bare `key` items
## (equation tags such as `[static]`), and -- when `bare_name` -- bare quoted
## strings `'value'`, which designate the `name` tag (model_remove /
## model_replace syntax).  Returns a named character vector (names = keys).
.mod_parse_tag_list <- function(txt, what, bare_name = FALSE) {
  items <- regmatches(txt, gregexpr("(?:'[^']*'|\"[^\"]*\"|[^,'\"])+", txt,
                                    perl = TRUE))[[1]]
  items <- trimws(items)
  items <- items[nzchar(items)]
  out <- character(0)
  for (it in items) {
    kv <- regmatches(it, regexec(
      "^([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(?:'([^']*)'|\"([^\"]*)\")$", it,
      perl = TRUE))[[1]]
    if (length(kv) == 4L) {
      out <- c(out, stats::setNames(paste0(kv[3L], kv[4L]), kv[2L]))
    } else if (bare_name && grepl("^(?:'[^']*'|\"[^\"]*\")$", it, perl = TRUE)) {
      out <- c(out, c(name = substr(it, 2L, nchar(it) - 1L)))
    } else if (!bare_name) {
      ## Equation tags are read leniently (parse_model_block() accepts any
      ## tag text): a bare `key`, or `key = value` with an unquoted value.
      ## Anything else cannot be selected by a tag list and is skipped.
      kv <- regmatches(it, regexec(
        "^([A-Za-z_][A-Za-z0-9_]*)\\s*(?:=\\s*(.*))?$", it, perl = TRUE))[[1]]
      if (length(kv) == 3L)
        out <- c(out, stats::setNames(trimws(kv[3L]), kv[2L]))
    } else {
      .dynhr_abort("parse_mod: malformed item `", it, "` in ", what, ".",
                   class = "dynhr_error_mod_syntax")
    }
  }
  out
}

## Model-block text as `;`-terminated statements: list(text, tags, local,
## rest) with `tags` the statement's leading `[...]` equation tags, `rest` the
## text after them and `local` TRUE for a `#` model-local definition.
.mod_model_units <- function(body) {
  pcs <- trimws(strsplit(body, ";", fixed = TRUE)[[1]])
  pcs <- pcs[nzchar(pcs)]
  lapply(pcs, function(s) {
    tg <- regmatches(s, regexec("^\\[([^\\]]*)\\]", s, perl = TRUE))[[1]]
    has_tag <- length(tg) == 2L
    rest <- if (has_tag) trimws(substring(s, nchar(tg[1L]) + 1L)) else s
    list(text = s,
         tags = if (has_tag) .mod_parse_tag_list(tg[2L], "an equation tag")
                else character(0),
         local = startsWith(rest, "#"), rest = rest)
  })
}

## Identifier tokens of model-block text (tags, strings and numbers skipped).
.mod_text_idents <- function(text) {
  toks <- regmatches(text, gregexpr(.MOD_TIMING_TOKEN_RE, text,
                                    perl = TRUE))[[1]]
  toks[grepl("^[A-Za-z_]", toks)]
}

## Identifiers of `text` with every `#` model-local it uses expanded (Dynare's
## collectVariables looks through model-local variables).
.mod_idents_thru_locals <- function(text, local_defs) {
  seen <- character(0)
  ids  <- .mod_text_idents(text)
  todo <- intersect(ids, names(local_defs))
  while (length(todo) > 0L) {
    l <- todo[[1L]]
    todo <- todo[-1L]
    if (l %in% seen) next
    seen <- c(seen, l)
    more <- .mod_text_idents(local_defs[[l]])
    ids  <- c(ids, more)
    todo <- c(todo, setdiff(intersect(more, names(local_defs)), seen))
  }
  unique(ids)
}

## `#` model-local definitions of a unit list: name -> definition text.
.mod_local_defs <- function(units) {
  defs <- character(0)
  for (u in units) if (u$local) {
    m <- regmatches(u$rest, regexec(
      "(?s)^#\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=(?!=)(.*)$", u$rest,
      perl = TRUE))[[1]]
    if (length(m) == 3L) defs[m[2L]] <- m[3L]
  }
  defs
}

## Every declared symbol, in declaration order (Dynare's symbol-table order,
## which fixes where an endogenous variable turned exogenous lands).
.mod_decl_order <- function(txt_decl) {
  ## options in balanced parentheses (group 1), names in group 2
  pat <- paste0("(?si)\\b(?:var(?!exo)|varexo(?!_det)|varexo_det|",
                "parameters)\\b(?:\\s*(\\((?:[^()]++|(?1))*\\)))?\\s+(.*?)\\s*;")
  hits <- regmatches(txt_decl, gregexpr(pat, txt_decl, perl = TRUE))[[1]]
  unlist(lapply(hits, function(h)
    parse_declaration_names(sub(pat, "\\2", h, perl = TRUE))))
}

## Dynare's model_replace / model_remove / var_remove statements (Dynare 6+),
## applied to the model-block text in the order they appear in the file.
##
## * `model_remove(TAGS);` -- TAGS is a comma-separated list of `'value'`
##   (the `name` tag) and `key = 'value'` items; every equation carrying a
##   listed tag is removed, and a listed tag that matches no equation is an
##   error.  Each removed equation must have an `endogenous = 'v'` tag or a
##   single endogenous variable on its left-hand side; that variable becomes
##   EXOGENOUS if it still appears in the remaining equations, otherwise it is
##   removed from the model (DynamicModel::removeEquations with
##   excluded_vars_change_type).  An endogenous variable turned exogenous takes
##   its declaration-order place among the exogenous variables.
## * `model_replace(TAGS); EQUATIONS end;` -- removes the tagged equations
##   (no variable changes type) and appends EQUATIONS at the end of the model.
## * `var_remove NAMES;` -- removes declared symbols (endogenous, exogenous or
##   parameters).  A removed symbol that is still used in the model is an
##   error (the Dynare preprocessor would emit an inconsistent model).
## Checked against the Dynare 7.1 preprocessor's `json=transform` output.
##
## @param txt      Comment-stripped .mod text (blocks present).
## @param body     The model-block body ("" when there is no model block).
## @param decl     list(endo, exo, exo_det, params) of declared names.
## @param txt_decl Declaration text (paired blocks stripped).
## @return list(body, changed, endo, exo, exo_det, params).
.mod_apply_model_edits <- function(txt, body, decl, txt_decl) {
  out <- c(list(body = body, changed = FALSE), decl)
  q <- "(?:'[^']*'|\"[^\"]*\"|[^)'\"])*"
  pats <- c(
    model_replace = paste0("(?si)\\bmodel_replace\\b\\s*\\((", q,
                           ")\\)\\s*;(.*?)\\bend\\s*;"),
    model_remove  = paste0("(?si)\\bmodel_remove\\b\\s*\\((", q, ")\\)\\s*;"),
    var_remove    = "(?si)\\bvar_remove\\b\\s+([^;]*);")
  edits <- list()
  for (ty in names(pats)) {
    n_kw <- length(regmatches(txt, gregexpr(paste0("\\b", ty, "\\b"), txt,
                                            perl = TRUE))[[1]])
    m <- gregexpr(pats[[ty]], txt, perl = TRUE)[[1]]
    n_ok <- if (m[1L] == -1L) 0L else length(m)
    if (n_kw != n_ok)
      .dynhr_abort("parse_mod: malformed `", ty, "` statement; the forms are ",
                   "`model_remove('name', key = 'value', ...);`, ",
                   "`model_replace('name', ...); EQUATIONS end;` and ",
                   "`var_remove NAMES;`.", class = "dynhr_error_mod_syntax")
    if (n_ok == 0L) next
    hits <- regmatches(txt, list(m))[[1]]
    for (k in seq_along(hits)) {
      g <- regmatches(hits[k], regexec(pats[[ty]], hits[k], perl = TRUE))[[1]]
      edits[[length(edits) + 1L]] <- list(
        type = ty, at = as.integer(m[k]), args = g[2L],
        body = if (length(g) >= 3L) g[3L] else "")
    }
  }
  if (length(edits) == 0L) return(out)
  edits <- edits[order(vapply(edits, `[[`, 0L, "at"))]

  units <- .mod_model_units(body)
  decl_order <- .mod_decl_order(txt_decl)
  excluded <- character(0)
  removed  <- character(0)
  for (ed in edits) {
    if (ed$type == "var_remove") {
      for (v in parse_declaration_names(ed$args)) {
        slot <- c("endo", "exo", "exo_det", "params")[
          c(v %in% out$endo, v %in% out$exo, v %in% out$exo_det,
            v %in% out$params)]
        if (length(slot) == 0L && !(v %in% excluded))
          .dynhr_abort("parse_mod: `var_remove` names `", v, "`, which is ",
                       "not a declared symbol.",
                       class = "dynhr_error_mod_syntax")
        for (s in slot) out[[s]] <- setdiff(out[[s]], v)
        removed <- c(removed, v)
      }
      next
    }
    sel <- .mod_parse_tag_list(ed$args, paste0("the ", ed$type, "(...) list"),
                               bare_name = TRUE)
    hit <- logical(length(units))
    missing <- character(0)
    for (k in seq_along(sel)) {
      mk <- vapply(units, function(u) !u$local &&
                     any(names(u$tags) == names(sel)[k] & u$tags == sel[[k]]),
                   logical(1))
      if (!any(mk))
        missing <- c(missing, paste0(names(sel)[k], "='", sel[[k]], "'"))
      hit <- hit | mk
    }
    if (length(missing) > 0L)
      .dynhr_abort("parse_mod: ", ed$type, ": no equation of the model ",
                   "carries the tag(s) ", paste(missing, collapse = ", "), ".",
                   class = "dynhr_error_mod_syntax")
    if (ed$type == "model_remove") {
      defs  <- .mod_local_defs(units)
      eq_no <- cumsum(!vapply(units, `[[`, logical(1), "local"))
      new_excl <- character(0)
      for (j in which(hit)) {
        u <- units[[j]]
        v <- if ("endogenous" %in% names(u$tags)) {
          u$tags[names(u$tags) == "endogenous"][[1L]]
        } else {
          eqp <- regexpr("(?<![<>!=])=(?!=)", u$rest, perl = TRUE)
          lhs <- if (eqp > 0L) substr(u$rest, 1L, eqp - 1L) else u$rest
          intersect(.mod_idents_thru_locals(lhs, defs), out$endo)
        }
        if (length(v) != 1L || !(v %in% out$endo))
          .dynhr_abort("parse_mod: model_remove: equation ", eq_no[j], " (`",
                       gsub("\\s+", " ", u$rest), "`) has been excluded but ",
                       "it does not have a single endogenous variable on its ",
                       "left-hand side or an `endogenous` tag naming one.",
                       class = "dynhr_error_mod_syntax")
        new_excl <- c(new_excl, v)
      }
      units <- units[!hit]
      defs  <- .mod_local_defs(units)
      used  <- unique(unlist(lapply(units, function(u)
        if (u$local) character(0) else .mod_idents_thru_locals(u$rest, defs))))
      for (v in unique(new_excl)) {
        out$endo <- setdiff(out$endo, v)
        if (v %in% used) {
          ex <- c(out$exo, v)
          out$exo <- ex[order(match(ex, decl_order), na.last = TRUE)]
        } else {
          excluded <- c(excluded, v)
        }
      }
    } else {
      units <- c(units[!hit], .mod_model_units(ed$body))
    }
  }
  if (length(removed) > 0L) {
    still <- intersect(removed, unique(unlist(lapply(units, function(u)
      .mod_text_idents(u$rest)))))
    if (length(still) > 0L)
      .dynhr_abort("parse_mod: `var_remove` removes ",
                   paste(still, collapse = ", "), ", which the model still ",
                   "uses.", class = "dynhr_error_mod_syntax")
  }
  out$body <- paste0(vapply(units, `[[`, "", "text"), ";", collapse = "\n")
  out$changed <- TRUE
  out
}


# ===========================================================================
# The `diff` operator
# ===========================================================================

## Dynare's `diff(EXPR)` = EXPR - EXPR(-1), where EXPR(-1) moves the timing of
## every variable (and `#` model-local) in EXPR by -1.  The Dynare preprocessor
## (DynamicModel::substituteDiff) introduces AUX_DIFF_<n> auxiliary variables
## for a lead-free argument and inlines an argument with a lead; either way the
## equations in the model's own variables are the same, and they are what this
## inline expansion produces (the `|lag| > 1` and exogenous-lag auxiliaries are
## then created by the usual expansion).  Checked against the Dynare 7.1
## preprocessor's `json=transform` output, nested diff() and diff() of an
## exogenous variable included.  A declared symbol named `diff` is left alone.
##
## `adl()` is NOT an operator of Dynare 7 -- the 7.1 preprocessor parses
## `adl(x, 'a', 2)` as an undeclared external function and rejects it -- so it
## is an error here too, with a pointer to the explicit form.
##
## @param body        Model-block text.
## @param shift_names Names whose timing the lag moves (endogenous, exogenous
##   and model-local names).
## @param declared    Every declared symbol name.
.mod_expand_diff_text <- function(body, shift_names, declared) {
  if (!("adl" %in% declared) &&
      grepl("\\badl\\s*\\(", gsub("\\[[^\\]]*\\]", " ", body, perl = TRUE),
            perl = TRUE))
    .dynhr_abort("parse_mod: `adl()` is not a Dynare 7 operator (the Dynare ",
                 "7.1 preprocessor rejects it as an undeclared external ",
                 "function). Write the distributed lag out, e.g. ",
                 "`a_lag_1*x(-1) + a_lag_2*x(-2)` with declared parameters.",
                 class = "dynhr_error_mod_syntax")
  if ("diff" %in% declared) return(body)
  shifts <- stats::setNames(rep(-1L, length(shift_names)), shift_names)
  repeat {
    m <- gregexpr(.MOD_TIMING_TOKEN_RE, body, perl = TRUE)[[1]]
    if (m[1L] == -1L) break
    starts <- as.integer(m)
    toks <- substring(body, starts, starts + attr(m, "match.length") - 1L)
    open <- NA_integer_
    for (j in which(toks == "diff")) {
      op <- regexpr("^\\s*\\(", substr(body, starts[j] + 4L, nchar(body)),
                    perl = TRUE)
      if (op > 0L) {
        open  <- starts[j] + 3L + attr(op, "match.length")
        start <- starts[j]
        break
      }
    }
    if (is.na(open)) break
    close <- .mod_matching_paren(body, open)
    if (is.na(close))
      .dynhr_abort("parse_mod: unbalanced parentheses in the argument of ",
                   "diff().", class = "dynhr_error_mod_syntax")
    arg <- trimws(substr(body, open + 1L, close - 1L))
    flat <- arg
    repeat {
      nxt <- gsub("\\([^()]*\\)", "", flat)
      if (identical(nxt, flat)) break
      flat <- nxt
    }
    if (!nzchar(arg) || grepl(",", flat, fixed = TRUE))
      .dynhr_abort("parse_mod: diff() takes exactly one argument (got `diff(",
                   arg, ")`).", class = "dynhr_error_mod_syntax")
    if (grepl("\\bEXPECTATION\\b", arg, perl = TRUE))
      .dynhr_abort("parse_mod: diff() of an EXPECTATION operator is not ",
                   "implemented (in `diff(", arg, ")`).",
                   class = "dynhr_error_mod_syntax")
    body <- paste0(substr(body, 1L, start - 1L),
                   "((", arg, ") - (",
                   trimws(.mod_retime_text(arg, shifts)), "))",
                   substr(body, close + 1L, nchar(body)))
  }
  body
}
