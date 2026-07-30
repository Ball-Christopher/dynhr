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
  # equation ASTs in parse_mod (step 11) BEFORE this classification runs, so by
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
    steady_state_model  = NULL,
    shocks              = list(variances = data.frame(), correlations = data.frame()),
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


#' Print method for dynhr_mod
#' @noRd
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


#' Summary method for dynhr_mod
#' @noRd
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
#' \emph{Declarations:} \code{var} (with \code{var(log)} and
#' \code{\%(long_name=...)} annotations), \code{varexo}, \code{varexo_det},
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
#' \code{EXPECTATION()}, standard arithmetic.
#'
#' \emph{Steady state:} \code{steady_state_model;} ... \code{end;} ---
#' sequential assignments evaluated by \code{eval_steady_state_model()}.
#'
#' \emph{Shocks:} \code{shocks;} ... \code{end;} --- \code{stderr},
#' \code{variance}, \code{corr}, \code{covar}, deterministic period shocks,
#' \code{shocks(overwrite)}.  Multiple \code{shocks} blocks are merged.
#'
#' \emph{Macro language:} \code{@#define}, \code{@#for}/\code{@#endfor},
#' \code{@#if}/\code{@#else}/\code{@#endif},
#' \code{@#ifdef}/\code{@#ifndef}, \code{@#include}, and
#' \code{@\{expr\}} interpolation.  Expanded before any other parsing.
#'
#' \emph{Other:} \code{initval}/\code{endval}, \code{estimated_params},
#' \code{estimated_params_init}, \code{occbin_constraints},
#' \code{planner_objective}, \code{verbatim} blocks.
#' Comments: \code{//}, \code{/* ... */}, \code{\%}.
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
    txt <- iconv(txt, from = "LATIN1", to = "UTF-8", sub = "?")
    # Remove BOM if present
    txt <- gsub("^\uFEFF", "", txt)
    if (verbose) cat("Read", nchar(txt), "characters from", source_file, "\n")
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

  # ---- 2. Expand Dynare macro directives (M18) -------------------------
  # Evaluate the supported @#for / @#if / @#ifdef / @#ifndef / @#define macro
  # subset and splice in @{...} interpolations BEFORE the lexer and any block
  # extraction run.  This is an exact byte-for-byte no-op for macro-free
  # source, so models without @# directives are entirely unaffected.
  # Unsupported directives FAIL LOUD (see R/parse-macro.R) rather than being
  # silently dropped, which would yield an under-specified model.
  txt <- expand_macros(txt,
                       mod_dir = if (!is.na(source_file)) dirname(source_file) else NULL)
  if (verbose) cat("Expanded Dynare macro directives\n")

  # ---- 3. Strip comments and macros ------------------------------------
  txt <- strip_comments_and_macros(txt)
  if (verbose) cat("Stripped comments and macros\n")

  # ---- 3. Extract declarations -----------------------------------------
  # Strip block constructs to prevent internal 'var' keywords
  # (e.g. "var eps_a = sig_a^2;" inside shocks block) from being
  # mistaken for variable declarations.
  txt_decl <- .strip_blocks(.strip_mod_comments(txt))

  var_text        <- extract_declaration(txt_decl, "var")
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
    message(sprintf("Parser: removing %d varexo names from var: %s",
                    length(exo_overlap), paste(exo_overlap, collapse = ", ")))
    var_names <- setdiff(var_names, varexo_names)
  }
  if (length(exo_det_overlap) > 0) {
    message(sprintf("Parser: removing %d varexo_det names from var: %s",
                    length(exo_det_overlap), paste(exo_det_overlap, collapse = ", ")))
    var_names <- setdiff(var_names, varexo_det_names)
  }
  # Also protect varexo from varexo_det contamination
  exo_cross <- intersect(varexo_names, varexo_det_names)
  if (length(exo_cross) > 0) {
    message(sprintf("Parser: removing %d varexo_det names from varexo: %s",
                    length(exo_cross), paste(exo_cross, collapse = ", ")))
    varexo_names <- setdiff(varexo_names, varexo_det_names)
  }

  all_var_names <- c(var_names, varexo_names, varexo_det_names)

  if (verbose) {
    cat("Variables:  ", length(var_names), "endogenous,",
        length(varexo_names), "exogenous,",
        length(varexo_det_names), "exo deterministic\n")
    cat("Parameters:", length(param_names), "\n")
  }

  # ---- 4. Extract and parse model block --------------------------------
  model_block <- extract_paired_block(txt, "model")
  model_opts  <- parse_command_options(model_block$options_str)

  equations   <- list()
  local_vars  <- list()

  if (model_block$found) {

    # ---- Auxiliary variable expansion for leads/lags > 1 ----
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
    if (verbose) cat("Parsed", length(equations), "equations\n")
  }

  # ---- 5. Calibration (top-level param assignments) --------------------
  remaining_txt <- remove_blocks(txt)
  # M19: pre-evaluate verbatim; blocks (which remove_blocks strips) so derived
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
    cat("Calibrated", length(param_values), "parameters\n")

  # ---- M17: fail-loud on declared-but-unvalued parameters --------------------
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
    warning(sprintf(
      paste0("parse_mod: %d declared parameter(s) have no value in the .mod ",
             "(likely computed in an external *_steadystate.m): %s. ",
             "Use inject_params(model, named_vec) to fill them from an ",
             "external (possibly superset) vector without errors; or use ",
             "set_param_values() if you have the exact model-parameter set. ",
             "Unset parameters will produce silently wrong decision rules."),
      length(.unset_params),
      paste(.unset_params, collapse = ", ")),
      call. = FALSE)
  }

  # ---- 6. Initval / Endval ---------------------------------------------
  # Create an evaluation environment with known parameter values so that
  # initval/endval expressions can reference parameter names directly.
  initval_env <- list2env(as.list(param_values), parent = baseenv())
  initval_block <- extract_paired_block(txt, "initval")
  initval <- if (initval_block$found)
    parse_initval_block(initval_block$body, env = initval_env)
  else numeric(0)

  endval_block <- extract_paired_block(txt, "endval")
  endval <- if (endval_block$found)
    parse_initval_block(endval_block$body, env = initval_env)
  else numeric(0)

  # ---- 7. Steady state model -------------------------------------------
  ssm_block <- extract_paired_block(txt, "steady_state_model")
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
      warning(sprintf(
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
    .ss_dir  <- dirname(source_file)
    .ss_stem <- sub("\\.(mod|dyn)$", "", basename(source_file), ignore.case = TRUE)
    # Review/harness copies are often suffixed (e.g. GK_2011_pp.mod) while the
    # external SS file keeps the original stem (GK_2011_steadystate.m). Try the
    # exact stem first, then the stem with trailing review suffixes stripped.
    .ss_stems <- unique(c(
      .ss_stem,
      sub("(_pp|_export|_stochsimul|_dynhr)+$", "", .ss_stem)
    ))
    .ss_files <- as.vector(t(outer(
      .ss_stems, c("_steadystate.m", "_steadystate2.m"),
      function(s, suf) file.path(.ss_dir, paste0(s, suf))
    )))
    .ss_file <- .ss_files[file.exists(.ss_files)]
    if (length(.ss_file) > 0L) {
      .ss_file <- .ss_file[[1L]]
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
        warning(sprintf(
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
  local_const_env <- list2env(as.list(param_values), parent = baseenv())
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
        # M5: a later shocks block silently overwriting an earlier variance for
        # the same shock loses the earlier value (e.g. Ascari keeps only the last
        # block's stderrs). Warn on a genuine conflict, mirroring the correlation
        # dedup warning. Re-stating the same value is silent (byte-compatible).
        old_v <- base_df$variance[idx]
        new_v <- nr$variance
        differs <- (is.na(old_v) != is.na(new_v)) ||
          (!is.na(old_v) && !is.na(new_v) &&
             !isTRUE(all.equal(old_v, new_v)))
        if (differs) {
          warning(sprintf(
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
    # M5: a second `shocks` block must not silently append a duplicate row for a
    # pair already defined (plain rbind let last-write-wins apply downstream with
    # no signal -- e.g. BKK's corr=0 IRF block vs corr=0.258 sim block). Replace
    # per (unordered) pair, keeping the last definition, and warn on a conflict.
    for (k in seq_len(nrow(new_df))) {
      nr  <- new_df[k, , drop = FALSE]
      idx <- which((base_df$var1 == nr$var1 & base_df$var2 == nr$var2) |
                   (base_df$var1 == nr$var2 & base_df$var2 == nr$var1))
      if (length(idx) > 0) {
        warning(sprintf(
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
  # M5: retain a per-block breakdown so per-`stoch_simul` scoping can be added
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
    warning(sprintf(
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
    est_m <- regmatches(txt, regexec(
      "(?si)\\bestimation\\s*\\(([^)]*)\\)", txt, perl = TRUE))[[1]]
    if (length(est_m) > 0 && nchar(est_m[1]) > 0) {
      fo_m <- regmatches(est_m[2], regexec(
        "(?i)\\bfirst_obs\\s*=\\s*([A-Za-z0-9]+)", est_m[2], perl = TRUE))[[1]]
      if (length(fo_m) > 0 && nchar(fo_m[1]) > 0 &&
          !is.null(.parse_dynare_date(fo_m[2])))
        sample_start <- fo_m[2]
    }
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
    sample_start_hs <- NULL
    est_m_hs <- regmatches(txt, regexec(
      "(?si)\\bestimation\\s*\\(([^)]*)\\)", txt, perl = TRUE))[[1]]
    if (length(est_m_hs) > 0 && nchar(est_m_hs[1]) > 0) {
      fo_m_hs <- regmatches(est_m_hs[2], regexec(
        "(?i)\\bfirst_obs\\s*=\\s*([A-Za-z0-9]+)", est_m_hs[2], perl = TRUE))[[1]]
      if (length(fo_m_hs) > 0 && nchar(fo_m_hs[1]) > 0 &&
          !is.null(.parse_dynare_date(fo_m_hs[2])))
        sample_start_hs <- fo_m_hs[2]
    }
    parse_heteroskedastic_shocks_block(hs_block$body, sample_start = sample_start_hs)
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
    cat("occbin_constraints:", length(occbin_constraints), "named constraints\n")

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
      if (verbose) cat("Found command:", cmd, "\n")
    }
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
    if (verbose) cat("Found command: planner_objective\n")
  }

  # ---- 10c. ramsey_policy instruments= --------------------------------
  ramsey_instruments <- extract_ramsey_instruments(txt)

  # ---- 10d. metadata (@dynhr: blocks) ---------------------------------
  metadata <- if (!is.na(source_file)) {
    tryCatch(extract_mod_metadata(source_file), error = function(e) list())
  } else {
    # Text input: write to a temp file so extract_mod_metadata can read lines
    tryCatch({
      tmp_meta <- tempfile(fileext = ".mod")
      on.exit(unlink(tmp_meta), add = TRUE)
      writeLines(strsplit(file_or_text, "\n")[[1]], tmp_meta)
      extract_mod_metadata(tmp_meta)
    }, error = function(e) list())
  }

  # ---- 11. Adjust predetermined variable timings in equation ASTs -------
  # Dynare's `predetermined_variables` convention: a variable declared
  # predetermined is a BEGINNING-of-period stock, so writing `k` in the model
  # means standard-timing `k(-1)`, `k(+1)` means standard `k`, `k(+2)` means
  # `k(+1)`, and so on.  EVERY occurrence of a predetermined variable must
  # therefore be re-timed by -1 (not just the lead_lag==0 occurrence).  Doing
  # this on the equation ASTs keeps the lead_lag_incidence, dynamic Jacobian,
  # variable classification, and QZ solver all consistent -- so a predetermined
  # stock is correctly seen as a plain lagged STATE rather than (when only the
  # t entry was shifted) a spurious "mixed" variable appearing at both t-1 and
  # t+1.  Models with no predetermined_variables are byte-for-byte unaffected
  # (the loop body never runs).
  #
  # Aux-var interaction: for leads > 1 on a predetermined var (e.g. k(+2) in
  # Aguiar-Gopinath), text-level auxiliary expansion (step 4) already ran and
  # turned k(+2) into AUX_LEAD_k_1(+1) with AUX_LEAD_k_1 = k(+1).  The selective
  # -1 shift below then re-times the k(+1) inside that aux definition to k(0),
  # i.e. AUX_LEAD_k_1 = k -- the correct standard timing -- since AUX_* names
  # are NOT predetermined and are left untouched.
  if (length(predet_vars) > 0) {
    for (i in seq_along(equations)) {
      equations[[i]]$lhs <- ast_shift_named_timing(equations[[i]]$lhs,
                                                   predet_vars, shift = -1L)
      equations[[i]]$rhs <- ast_shift_named_timing(equations[[i]]$rhs,
                                                   predet_vars, shift = -1L)
    }
    if (verbose) cat("  Re-timed", length(predet_vars),
                     "predetermined var(s) by -1 in equation ASTs\n")
  }

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
    steady_state_model   = ssm,
    shocks               = shocks,
    shocks_blocks        = shocks_blocks_parsed,
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
    ## alias: estimation entry points take an `obs_vars` argument — populate
    ## the same-named model field from the .mod's varobs line so callers can
    ## omit it (pathological-DSGE paper gap #4, 2026-07)
    obs_vars             = varobs_names,
    ramsey_instruments   = ramsey_instruments,
    metadata             = metadata
  )

  # ---- 13. MCP constraint tags (M16) ----------------------------------
  # mcp_parse_tags() reads equation $tag_raw fields and builds MCP spec
  # objects.  Wiring it here mirrors the M6/M13 metadata-wiring pattern
  # (same parse-wiring family).  Errors are caught so a tag-syntax problem
  # does not crash parse_mod() for the whole model.
  model$mcp_constraints <- tryCatch(
    mcp_parse_tags(model, verbose = FALSE),
    error = function(e) {
      warning("parse_mod(): mcp_parse_tags() failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
  if (verbose && length(model$mcp_constraints) > 0)
    cat("MCP constraints:", length(model$mcp_constraints), "\n")

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
    warning(
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

  if (verbose) cat("Done. dynhr_mod object created.\n")
  model
}
