## R/parse-blocks.R
## --------------------------------------------------------------------------
## Declaration and block extraction for .mod files: var / varexo /
## varexo_det / parameters declarations, shocks block, estimated_params,
## initval, endval, steady_state_model, and command-line options.
##
## Phase-1 split from parser-monolith.R (no logic changes).
## --------------------------------------------------------------------------

#' Extract the content of paired blocks (keyword; ... end;)
#'
#' @param txt       Cleaned .mod text (comments already stripped).
#' @param keyword   Block keyword (e.g. "model", "shocks", "initval").
#' @return A list with:
#'   - options_str: text inside optional parentheses after keyword
#'   - body:       text between keyword; and end;
#'   - found:      logical
#' @noRd
extract_paired_block <- function(txt, keyword) {
  pat <- paste0(
    "(?si)\\b", keyword,
    "\\s*(?:\\(([^)]*)\\))?\\s*;",  # optional (options);
    "(.*?)",                         # body (non-greedy)
    "\\bend\\s*;"                    # end;
  )
  m <- regmatches(txt, regexec(pat, txt, perl = TRUE))[[1]]
  if (length(m) == 0)
    return(list(options_str = "", body = "", found = FALSE))
  list(
    options_str = trimws(m[2] %||% ""),
    body        = trimws(m[3]),
    found       = TRUE
  )
}


#' Extract ALL occurrences of a paired block (keyword; ... end;)
#'
#' Unlike \code{extract_paired_block()} which returns only the first match,
#' this function uses \code{gregexpr} to find every occurrence.  Used for
#' \code{shocks} blocks, which may appear more than once in a .mod file.
#'
#' @param txt     Cleaned .mod text (comments already stripped).
#' @param keyword Block keyword (e.g. "shocks").
#' @return A (possibly empty) list of block objects, each with slots
#'   \code{options_str}, \code{body}, and \code{found = TRUE}.
#' @noRd
extract_all_paired_blocks <- function(txt, keyword) {
  pat <- paste0(
    "(?si)\\b", keyword,
    "\\s*(?:\\(([^)]*)\\))?\\s*;",   # optional (options);
    "(.*?)",                           # body (non-greedy)
    "\\bend\\s*;"                      # end;
  )
  ms   <- gregexpr(pat, txt, perl = TRUE)[[1]]
  if (ms[1] == -1L) return(list())
  lens <- attr(ms, "match.length")
  lapply(seq_along(ms), function(k) {
    chunk <- substr(txt, ms[k], ms[k] + lens[k] - 1L)
    parts <- regmatches(chunk, regexec(pat, chunk, perl = TRUE))[[1]]
    list(options_str = trimws(if (length(parts) >= 2) parts[2] else ""),
         body        = trimws(if (length(parts) >= 3) parts[3] else ""),
         found       = TRUE)
  })
}


#' Extract a paired block that may contain nested keyword; ... end; sub-blocks
#'
#' Like extract_paired_block, but matches the depth-0 end; rather than the
#' first end;. Needed for occbin_constraints which contains equations; ... end;
#' sub-blocks that would otherwise absorb the outer end;.
#'
#' @param txt            Cleaned .mod text.
#' @param keyword        Outer block keyword (e.g. "occbin_constraints").
#' @param nested_keyword Keyword of inner sub-blocks (e.g. "equations"). Only
#'                       one level of nesting is supported.
#' @return List with options_str, body, found (same shape as extract_paired_block).
#' @noRd
extract_paired_block_nested <- function(txt, keyword, nested_keyword = "equations") {
  # Find the outer opening tag: keyword [( options )] ;
  open_pat <- paste0("(?si)\\b", keyword, "\\s*(?:\\(([^)]*)\\))?\\s*;")
  open_m   <- regexpr(open_pat, txt, perl = TRUE)
  if (open_m == -1L)
    return(list(options_str = "", body = "", found = FALSE))

  opts_capture <- regmatches(txt, regexec(open_pat, txt, perl = TRUE))[[1]]
  options_str  <- trimws(opts_capture[2] %||% "")

  body_start <- open_m + attr(open_m, "match.length")   # char after opening ;

  # Find all nested-open and end tokens after body_start.
  # The nested keyword may or may not be followed by ; (Dynare uses
  # "equations\n" without semicolon inside occbin_constraints).
  nested_pat <- paste0("(?si)\\b", nested_keyword, "\\b")
  end_pat    <- "(?i)\\bend\\s*;"

  nested_pos <- gregexpr(nested_pat, txt, perl = TRUE)[[1]]
  end_pos    <- gregexpr(end_pat,    txt, perl = TRUE)[[1]]

  # Filter to positions >= body_start
  if (nested_pos[1] != -1L) nested_pos <- nested_pos[nested_pos >= body_start]
  else                       nested_pos <- integer(0)
  if (end_pos[1] != -1L)    end_pos    <- end_pos[end_pos >= body_start]
  else                       end_pos    <- integer(0)

  # Build sorted event sequence: type +1 (open) or -1 (close), position
  opens  <- if (length(nested_pos) > 0) data.frame(pos = nested_pos, delta = +1L) else data.frame(pos = integer(0), delta = integer(0))
  closes <- if (length(end_pos)    > 0) data.frame(pos = end_pos,    delta = -1L) else data.frame(pos = integer(0), delta = integer(0))
  events <- rbind(opens, closes)
  events <- events[order(events$pos), , drop = FALSE]

  depth    <- 1L   # we are inside the outer block
  body_end <- NA_integer_

  for (k in seq_len(nrow(events))) {
    depth <- depth + events$delta[k]
    if (depth == 0L) {
      body_end <- events$pos[k] - 1L   # just before this end;
      break
    }
  }

  if (is.na(body_end))
    return(list(options_str = "", body = "", found = FALSE))

  body <- trimws(substr(txt, body_start, body_end))
  list(options_str = options_str, body = body, found = TRUE)
}


#' Extract a declaration block (keyword names... ;)
#'
#' Uses (?si) so the match can span multiple lines (NZSim has 71-variable
#' declarations spread across many lines).
#'
#' v0.3: Now finds ALL matching declarations (not just the first) and
#' handles optional parenthesised options after keyword, e.g. var(log).
#'
#' @param txt     Cleaned .mod text.
#' @param keyword Declaration keyword (e.g. "var", "varexo", "parameters").
#' @return Character string of the declaration body, or "" if not found.
#' @noRd
extract_declaration <- function(txt, keyword, debug = FALSE) {
  if (keyword == "var") {
    kw_pat <- "\\bvar(?!exo)\\b"
  } else if (keyword == "varexo") {
    kw_pat <- "\\bvarexo(?!_det)\\b"
  } else {
    kw_pat <- paste0("\\b", keyword, "\\b")
  }

  pat <- paste0("(?si)", kw_pat, "(?:\\s*\\([^)]*\\))?\\s+(.*?)\\s*;")

  all_positions <- gregexpr(pat, txt, perl = TRUE)
  all_matches   <- regmatches(txt, all_positions)[[1]]

  # -- DIAGNOSTIC --
  if (debug){
    cat(sprintf("  [extract_declaration] keyword='%s'  matches=%d\n", keyword, length(all_matches)))
    for (j in seq_along(all_matches)) {
      cat(sprintf("    match %d (nchar=%d): %.120s\n", j, nchar(all_matches[j]),
                  gsub("\\s+", " ", all_matches[j])))
    }
  }
  # ----------------

  if (length(all_matches) == 0) return("")

  bodies <- vapply(all_matches, function(m) {
    parts <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
    if (length(parts) >= 2) trimws(parts[2]) else ""
  }, character(1))

  result <- paste(bodies[nchar(bodies) > 0], collapse = " ")

  # -- DIAGNOSTIC --
  if (debug) cat(sprintf("    -> extracted names text: %.200s\n", gsub("\\s+", " ", result)))
  # ----------------

  result
}

#' Extract a command invocation (e.g. stoch_simul(options) var_list ;)
#'
#' @param txt     Cleaned .mod text.
#' @param command Command name (e.g. "stoch_simul", "estimation", "steady").
#' @return A list with options_str, var_list, found.
#' @noRd
extract_command <- function(txt, command) {
  pat <- paste0(
    "(?i)\\b", command,
    "\\s*(?:\\(([^)]*)\\))?",  # optional (options)
    "\\s*([^;]*?)\\s*;"        # optional variable list, then ;
  )
  m <- regmatches(txt, regexec(pat, txt, perl = TRUE))[[1]]
  if (length(m) == 0)
    return(list(options_str = "", var_list = "", found = FALSE))
  list(
    options_str = trimws(m[2] %||% ""),
    var_list    = trimws(m[3] %||% ""),
    found       = TRUE
  )
}

#' Extract planner_objective expression
#'
#' Dynare's syntax is a bare command, NOT a function call:
#'   planner_objective EXPRESSION;
#' (e.g. \code{planner_objective log(C) - L^(1+epsilon)/(1+epsilon);}). The
#' expression itself may contain balanced parentheses (\code{log(C)}); only the
#' terminating \code{;} delimits it. We therefore capture everything between the
#' keyword and the first \code{;}. A parenthesised form \code{planner_objective(EXPR);}
#' is also tolerated — the outer parens become part of EXPR, which the
#' expression parser handles transparently.
#'
#' @param txt Cleaned .mod text.
#' @return List with \code{text} and \code{found}.
#' @noRd
extract_planner_objective <- function(txt) {
  m <- regmatches(
    txt,
    regexec("(?si)\\bplanner_objective\\b\\s*(.+?)\\s*;", txt, perl = TRUE)
  )[[1]]
  if (length(m) < 2) return(list(text = "", found = FALSE))
  list(text = trimws(m[2]), found = TRUE)
}


#' Extract instrument names from ramsey_policy(instruments=(...))
#'
#' Parses the \code{instruments=(i)} option from a \code{ramsey_policy} or
#' \code{ramsey_model} command line into a character vector of variable names.
#'
#' @param txt Cleaned .mod text.
#' @return Character vector of instrument names, or \code{character(0)} if none.
#' @noRd
extract_ramsey_instruments <- function(txt) {
  # Match ramsey_policy(...instruments=(a b c)...) or instruments=(a,b,c)
  m <- regmatches(
    txt,
    regexec(
      "(?si)\\bramsey_(?:policy|model)\\s*\\([^)]*\\binstruments\\s*=\\s*\\(([^)]*)\\)",
      txt, perl = TRUE
    )
  )[[1]]
  if (length(m) < 2 || !nzchar(trimws(m[2])))
    return(character(0))
  parse_declaration_names(m[2])
}


#' Remove all recognised blocks from .mod text (to expose top-level calibration)
#'
#' @param txt Cleaned .mod text.
#' @return Text with all blocks replaced by whitespace.
#' @noRd
remove_blocks <- function(txt) {
  paired_kw <- c("model", "initval", "endval", "steady_state_model",
                 "shocks", "mshocks", "estimated_params_init",
                 "estimated_params_bounds", "estimated_params",
                 "observation_trends", "optim_weights",
                 "osr_params_bounds", "ramsey_constraints",
                 "moment_calibration", "irf_calibration",
                 "matched_moments", "occbin_constraints",
                 "epilogue", "pac_model", "var_model",
                 "trend_component_model", "filter_initial_state",
                 "homotopy_setup", "histval", "shock_groups",
                 "conditional_forecast_paths", "deterministic_trends",
                 "svar_identification", "identification",
                 "verbatim")
  for (kw in paired_kw) {
    pat <- paste0("(?si)\\b", kw, "\\s*(?:\\([^)]*\\))?\\s*;.*?\\bend\\s*;")
    txt <- gsub(pat, " ", txt, perl = TRUE)
  }
  decl_kw <- c("var", "varexo_det", "varexo", "varobs", "parameters",
               "predetermined_variables", "trend_var",
               "model_local_variable", "var_expectation",
               "var_remove", "model_replace", "model_remove")
  # v0.3: handle optional parenthesised options, e.g. var(log) x y;
  for (kw in decl_kw) {
    pat <- paste0("(?si)\\b", kw, "\\b(?:\\s*\\([^)]*\\))?\\s+.*?;")
    txt <- gsub(pat, " ", txt, perl = TRUE)
  }
  cmd_kw <- c("stoch_simul", "estimation", "steady", "check",
              "model_diagnostics", "model_info", "simul",
              "perfect_foresight_setup", "perfect_foresight_solver",
              "forecast", "conditional_forecast",
              "smoother2histval", "shock_decomposition",
              "calib_smoother", "extended_path",
              "osr", "ramsey_model", "ramsey_policy",
              "discretionary_policy", "planner_objective",
              "dynare_sensitivity", "bvar_density",
              "bvar_forecast", "dsample", "Sigma_e")
  for (kw in cmd_kw) {
    pat <- paste0("(?i)\\b", kw, "\\b\\s*(?:\\([^)]*\\))?\\s*[^;]*;")
    txt <- gsub(pat, " ", txt, perl = TRUE)
  }
  txt
}


#' Parse a variable/parameter declaration body into a character vector of names
#'
#' Handles:
#'   - Space-separated or comma-separated names
#'   - LaTeX names like ${y}$ or (long_name='output') annotations
#'   - Parenthesised options after variable names
#'
#' @param decl_text Body text from extract_declaration().
#' @return Character vector of names.
#' @noRd
parse_declaration_names <- function(decl_text) {
  if (nchar(decl_text) == 0) return(character(0))
  # Remove LaTeX names: $...$  or ${...}$
  txt <- gsub("\\$[^$]*\\$", " ", decl_text, perl = TRUE)
  # Remove parenthesised annotations: (long_name='...' , ...)
  # Handle one level of nested parens: (AR(1)) inside (long_name='AR(1)...').
  txt <- gsub("\\([^()]*(?:\\([^()]*\\))?[^()]*\\)", " ", txt, perl = TRUE)
  # Replace commas with spaces
  txt <- gsub(",", " ", txt)
  # Split on whitespace
  parts <- strsplit(trimws(txt), "\\s+")[[1]]
  parts <- parts[nchar(parts) > 0]
  # Filter: names must be valid identifiers
  parts[grepl("^[A-Za-z_][A-Za-z0-9_]*$", parts)]
}


#' Parse top-level parameter assignments outside any block
#'
#' Matches lines like:  alpha = 0.33;  beta = 1/(1+0.02/4);
#' v0.3: Now handles multiline expressions by collapsing whitespace.
#'
#' @param txt         Cleaned .mod text with blocks removed.
#' @param param_names Declared parameter names.
#' @return Named numeric vector of parameter values.
#' @noRd
parse_calibration <- function(txt, param_names, seed_env = NULL, quiet = FALSE) {
  values <- numeric(0)

  # Find statements of the form: IDENT = expression ;
  # Use [^;\\n] to prevent matching across newlines — this stops a
  # statement like  title_string='foo' (no trailing ;) from greedily
  # consuming  beta = 0.99;  on a later line (both lack ;, so the
  # old [^;]+ would merge them into one giant match).
  pat <- "([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*([^;\\n]+)\\s*;"
  matches <- gregexpr(pat, txt, perl = TRUE)
  all_matches <- regmatches(txt, matches)[[1]]

  # `seed_env` (optional) carries intermediate values computed in verbatim;
  # blocks (matrices like V, Correlation_matrix, and scalars like P0_z_bar0)
  # so that top-level assignments referencing them -- e.g.
  # `sigma_z = sqrt(V(1,1));` -- resolve instead of erroring to NA (repl M19).
  # MATLAB-style indexing in such RHS expressions is translated to R below.
  .pcal_env <- if (is.null(seed_env)) new.env(parent = baseenv())
               else list2env(as.list(seed_env), parent = baseenv())

  for (m in all_matches) {
    parts <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
    name <- parts[2]
    expr_text <- parts[3]

    # v0.3: collapse multiline expressions to a single line so that
    # R's parse() sees e.g. "(1 - theta) * (1 - beta*theta) / theta"
    # instead of "(1 - theta)\n      * (1 - beta*theta)\n      / theta"
    expr_text <- gsub("\\s+", " ", trimws(expr_text))

    # Evaluate expression in the calibration environment.
    # I2: evaluate ALL assignments in source order (not just declared params)
    # so that undeclared intermediate names (e.g. `alpha`, `theta`, `omega` in
    # Adam-Billi 2006) accumulate in .pcal_env before the derived declared
    # param expressions that reference them are reached.
    val <- tryCatch(eval(parse(text = expr_text), envir = .pcal_env),
                    error = function(e) NULL)

    # M19: when a verbatim seed_env is in play, the RHS may use MATLAB
    # indexing (e.g. `sqrt(V(1,1))`).  Retry with a MATLAB->R translation.
    if (is.null(val) && !is.null(seed_env)) {
      r_expr <- .matlab_stmt_to_r(expr_text)
      if (!is.null(r_expr))
        val <- tryCatch(eval(parse(text = r_expr), envir = .pcal_env),
                        error = function(e) NULL)
    }

    if (is.numeric(val) && length(val) == 1) {
      # Always accumulate in the eval env (intermediate or declared param).
      assign(name, val, envir = .pcal_env)
      # Only record in the output vector when it is a declared parameter.
      if (name %in% param_names)
        values[name] <- val
    }
  }

  # I2 fail-loud: declared params that are still NA (or absent) after
  # evaluation were referenced in an expression whose RHS could not be
  # resolved.  Name the unresolved params so the failure is not silent.
  .missing <- setdiff(param_names, names(values))
  .still_na <- names(values)[is.na(values)]
  .unresolved <- union(.missing, .still_na)
  if (length(.unresolved) > 0L && !isTRUE(quiet)) {
    # Identify which names in each failing RHS were undefined in .pcal_env
    .diagnose_unresolved <- function(pn) {
      # Re-scan the assignment list for this param's last RHS
      for (m in rev(all_matches)) {
        pts <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
        if (pts[2] != pn) next
        rhs <- gsub("\\s+", " ", trimws(pts[3]))
        # Collect identifier tokens in the RHS
        toks <- regmatches(rhs, gregexpr("[A-Za-z_][A-Za-z0-9_]*", rhs))[[1]]
        undef <- Filter(function(tok) {
          !exists(tok, envir = .pcal_env, inherits = FALSE) &&
            !(tok %in% c("TRUE","FALSE","Inf","NaN","NA","NULL")) &&
            is.null(tryCatch(get(tok, envir = baseenv()), error = function(e) NULL))
        }, toks)
        if (length(undef) > 0L) return(paste(unique(undef), collapse = ", "))
        return(NULL)
      }
      NULL
    }
    msgs <- vapply(.unresolved, function(pn) {
      diag <- .diagnose_unresolved(pn)
      if (!is.null(diag)) paste0(pn, " (undefined: ", diag, ")")
      else pn
    }, character(1))
    warning(sprintf(
      paste0("parse_calibration: %d declared parameter(s) could not be ",
             "resolved -- their values will be NA: %s"),
      length(.unresolved), paste(msgs, collapse = "; ")),
      call. = FALSE)
  }

  values
}


#' Parse an initval or endval block into a named numeric vector
#'
#' Expects lines of the form:  var_name = value ;
#' v0.3: Now handles multiline expressions by collapsing whitespace.
#'
#' @param body Body text of the initval/endval block.
#' @return Named numeric vector.
#' @noRd
parse_initval_block <- function(body, env = parent.frame()) {
  values <- numeric(0)
  # Use [^;\\n] to prevent cross-line greedy consumption (same fix as
  # parse_calibration).  Assign each result back to env so that subsequent
  # statements can reference earlier ones (e.g. Y=1.05; NX=exp(Y)-...).
  pat <- "([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*([^;\\n]+)\\s*;"
  matches <- gregexpr(pat, body, perl = TRUE)
  all_matches <- regmatches(body, matches)[[1]]

  for (m in all_matches) {
    parts <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
    name <- parts[2]
    expr_text <- parts[3]
    # v0.3: collapse multiline expressions to a single line
    expr_text <- gsub("\\s+", " ", trimws(expr_text))
    val <- tryCatch(eval(parse(text = expr_text), envir = env),
                    error = function(e) NULL)
    if (is.numeric(val) && length(val) == 1) {
      values[name] <- val
      assign(name, val, envir = env)
    }
  }
  values
}


#' Parse the steady_state_model block into a list of assignment ASTs
#'
#' @param body        Body text of the steady_state_model block.
#' @param var_names   Declared variable names.
#' @param param_names Declared parameter names.
#' @return List of assignment objects.
#' @noRd
parse_steady_state_model <- function(body, var_names, param_names) {
  assignments <- list()

  # Join Dynare continuation lines: a line starting with "^" continues
  # the previous line AND the "^" is itself the power operator.
  body_lines <- strsplit(body, "\n")[[1]]
  joined <- character(0)
  for (ln in body_lines) {
    if (grepl("^\\s*\\^", ln) && length(joined) > 0) {
      # Continuation: the leading ^ is the power operator.  Append it
      # followed by the rest of the line (without leading whitespace).
      rest <- sub("^\\s*\\^\\s*", "", ln)
      joined[length(joined)] <- paste0(joined[length(joined)], "^", rest)
    } else {
      joined <- c(joined, ln)
    }
  }
  body <- paste(joined, collapse = "\n")

  stmts <- strsplit(body, ";")[[1]]
  stmts <- trimws(stmts)
  stmts <- stmts[nchar(stmts) > 0]

  for (stmt in stmts) {
    if (grepl("^\\s*end\\s*$", stmt, ignore.case = TRUE)) next
    if (!grepl("=", stmt)) next

    # Clean the statement text for eval: remove comments and collapse
    # multi-line expressions (Dynare continuation via ^).
    stmt_clean <- gsub("%[^\n]*", "", stmt)          # remove % comments
    stmt_clean <- gsub("\\s+", " ", stmt_clean)       # collapse whitespace

    parts <- strsplit(stmt_clean, "\\s*=\\s*", perl = TRUE)[[1]]
    name <- trimws(parts[1])
    expr_text <- trimws(paste(parts[-1], collapse = "="))

    expr_ast <- parse_expression(expr_text, var_names, param_names)
    assignments <- c(assignments, list(list(
      name = name,
      expr = expr_ast,
      text = stmt_clean
    )))
  }
  assignments
}


#' Parse the shocks block into a structured list
#'
#' Handles:
#'   - var e;  stderr value;
#'   - var e = value;
#'   - corr e1, e2 = value;
#'
#' @param body Body text of the shocks block.
#' @return A list with:
#'   - variances:    data.frame (name, stderr, variance, stderr_expr,
#'                   variance_expr). The numeric columns are snapshots
#'                   evaluated against the calibrated param_env; the *_expr
#'                   columns keep the raw expression text so
#'                   .get_shock_stderr() can re-evaluate it against the
#'                   current parameter vector (estimated shock stds must
#'                   track theta, not the calibration).
#'   - correlations: data.frame (var1, var2, corr, corr_expr, cov, cov_expr).
#'                   For \code{corr a,b=expr} rows: corr / corr_expr hold the
#'                   ratio snapshot and expression; cov / cov_expr are NA.
#'                   For \code{var a,b=expr} rows (absolute covariance):
#'                   cov / cov_expr hold the snapshot and expression; corr /
#'                   corr_expr are NA.  .get_shock_cov() distinguishes the two
#'                   types by checking which columns are non-NA.
#' @noRd
parse_shocks_block <- function(body, param_env = NULL) {
  # If parameters available, eval expressions in that environment so
  # `var eps_a = sig_a^2;` works. Otherwise default to base env (numeric only).
  eval_env <- if (is.environment(param_env)) param_env
              else if (is.list(param_env)) list2env(param_env, parent = baseenv())
              else if (is.numeric(param_env) && !is.null(names(param_env)))
                list2env(as.list(param_env), parent = baseenv())
              else baseenv()
  safe_eval <- function(text) {
    tryCatch(eval(parse(text = text), envir = eval_env),
             error = function(e) NA_real_)
  }
  variances    <- data.frame(name = character(0),
                             stderr = numeric(0),
                             variance = numeric(0),
                             stderr_expr = character(0),
                             variance_expr = character(0),
                             skew = numeric(0),
                             skew_expr = character(0),
                             stringsAsFactors = FALSE)
  correlations <- data.frame(var1 = character(0),
                             var2 = character(0),
                             corr = numeric(0),
                             corr_expr = character(0),
                             cov  = numeric(0),
                             cov_expr  = character(0),
                             stringsAsFactors = FALSE)

  stmts <- strsplit(body, ";")[[1]]
  stmts <- trimws(stmts)
  stmts <- stmts[nchar(stmts) > 0]

  # Helper: expand period tokens like "1:4 2 5" -> integer vector
  .parse_periods <- function(text) {
    tokens <- strsplit(trimws(text), "[[:space:],]+")[[1]]
    tokens <- tokens[nchar(tokens) > 0]
    result <- integer(0)
    for (tok in tokens) {
      rng <- regmatches(tok, regexec("^([0-9]+):([0-9]+)$", tok))[[1]]
      if (length(rng) == 3) {
        result <- c(result, seq(as.integer(rng[2]), as.integer(rng[3])))
      } else {
        result <- c(result, as.integer(tok))
      }
    }
    result
  }

  deterministic <- data.frame(name   = character(0),
                               period = integer(0),
                               value  = numeric(0),
                               stringsAsFactors = FALSE)

  i <- 1L
  while (i <= length(stmts)) {
    s <- stmts[i]

    if (grepl("^\\s*end\\s*$", s, ignore.case = TRUE)) { i <- i+1; next }

    # Pattern: var IDENT = value  (variance directly)
    m1 <- regmatches(s, regexec(
      "^\\s*var\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.+)$",
      s, perl = TRUE))[[1]]
    if (length(m1) > 0 && nchar(m1[1]) > 0) {
      vname <- m1[2]
      val <- safe_eval(m1[3])
      variances <- rbind(variances,
                         data.frame(name = vname, stderr = sqrt(abs(val)),
                                    variance = val,
                                    stderr_expr = NA_character_,
                                    variance_expr = trimws(m1[3]),
                                    skew = 0,
                                    skew_expr = NA_character_,
                                    stringsAsFactors = FALSE))
      i <- i + 1; next
    }

    # Pattern: var IDENT, IDENT = expr  (off-diagonal covariance; Dynare syntax)
    # Must be checked BEFORE the bare m2 pattern so it does not silently fall
    # through to m2's match on "var IDENT" (which would drop ", IDENT = expr").
    m1b <- regmatches(s, regexec(
      paste0("^\\s*var\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*,\\s*",
             "([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.+)$"),
      s, perl = TRUE))[[1]]
    if (length(m1b) > 0 && nchar(m1b[1]) > 0) {
      cov_val <- safe_eval(m1b[4])
      correlations <- rbind(correlations,
                            data.frame(var1 = m1b[2], var2 = m1b[3],
                                       corr = NA_real_,
                                       corr_expr = NA_character_,
                                       cov  = cov_val,
                                       cov_expr  = trimws(m1b[4]),
                                       stringsAsFactors = FALSE))
      i <- i + 1; next
    }

    # Pattern: var IDENT  -- followed by stderr (stochastic) OR
    #                         periods + values (deterministic PF shocks)
    m2 <- regmatches(s, regexec(
      "^\\s*var\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*$",
      s, perl = TRUE))[[1]]
    if (length(m2) > 0 && nchar(m2[1]) > 0) {
      vname <- m2[2]
      if (i + 1 <= length(stmts)) {
        s2 <- stmts[i + 1]
        # --- stochastic: var IDENT; stderr VALUE; ---
        m3 <- regmatches(s2, regexec(
          "^\\s*stderr\\s+(.+)$", s2, perl = TRUE))[[1]]
        if (length(m3) > 0 && nchar(m3[1]) > 0) {
          val <- safe_eval(m3[2])
          variances <- rbind(variances,
                             data.frame(name = vname, stderr = val,
                                        variance = val^2,
                                        stderr_expr = trimws(m3[2]),
                                        variance_expr = NA_character_,
                                        skew = 0,
                                        skew_expr = NA_character_,
                                        stringsAsFactors = FALSE))
          i <- i + 2; next
        }
        # --- deterministic: var IDENT; periods ...; values ...; ---
        mp <- regmatches(s2, regexec(
          "^\\s*periods\\s+(.+)$", s2, perl = TRUE))[[1]]
        if (length(mp) > 0 && nchar(mp[1]) > 0 && i + 2 <= length(stmts)) {
          s3 <- stmts[i + 2]
          mv <- regmatches(s3, regexec(
            "^\\s*values\\s+(.+)$", s3, perl = TRUE))[[1]]
          if (length(mv) > 0 && nchar(mv[1]) > 0) {
            ## Parse periods into per-ENTRY groups: each comma/space-separated
            ## token is one entry; a range `a:b` is a single entry spanning
            ## several periods. This preserves the entry<->value pairing Dynare
            ## uses for the grouped syntax `periods 1:12, 13:24; values 2.1, 0.66`
            ## (one value per period GROUP), which a flat period vector loses.
            p_tokens <- strsplit(trimws(mp[2]), "[[:space:],]+")[[1]]
            p_tokens <- p_tokens[nchar(p_tokens) > 0]
            token_periods <- lapply(p_tokens, function(tok) {
              rng <- regmatches(tok, regexec("^([0-9]+):([0-9]+)$", tok))[[1]]
              if (length(rng) == 3)
                seq(as.integer(rng[2]), as.integer(rng[3]))
              else as.integer(tok)
            })
            periods_vec <- unlist(token_periods)
            val_tokens  <- strsplit(trimws(mv[2]), "[[:space:],]+")[[1]]
            val_tokens  <- val_tokens[nchar(val_tokens) > 0]
            values_vec  <- vapply(val_tokens, safe_eval, numeric(1))
            ## Map values to periods (Dynare semantics):
            ##   1 value           -> recycled to all periods;
            ##   #values==#periods -> one value per period (a range with several
            ##                        values is per-period within the range);
            ##   #values==#entries -> one value per period GROUP (expand each
            ##                        value across its entry's periods).
            if (length(values_vec) == 1L) {
              values_vec <- rep(values_vec, length(periods_vec))
            } else if (length(values_vec) == length(periods_vec)) {
              # per-period: already aligned
            } else if (length(values_vec) == length(token_periods)) {
              values_vec <- unlist(Map(function(p, v) rep(v, length(p)),
                                       token_periods, values_vec))
            } else {
              stop(sprintf(paste0("parse_shocks_block: deterministic shock '%s' ",
                "has %d value(s) for %d period(s) in %d group(s); expected 1, ",
                "one per period, or one per group."), vname, length(values_vec),
                length(periods_vec), length(token_periods)), call. = FALSE)
            }
            deterministic <- rbind(deterministic,
                                   data.frame(name   = vname,
                                              period = periods_vec,
                                              value  = values_vec,
                                              stringsAsFactors = FALSE))
            i <- i + 3; next
          }
        }
      }
      i <- i + 1; next
    }

    # Pattern: corr IDENT, IDENT = expr
    # Store both the parse-time snapshot (corr) and the raw expression text
    # (corr_expr) so .get_shock_cov() can re-evaluate it against the current
    # parameter vector during MCMC (mirrors stderr_expr / variance_expr).
    m4 <- regmatches(s, regexec(
      paste0("^\\s*corr\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*,\\s*",
             "([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.+)$"),
      s, perl = TRUE))[[1]]
    if (length(m4) > 0 && nchar(m4[1]) > 0) {
      val <- safe_eval(m4[4])
      correlations <- rbind(correlations,
                            data.frame(var1 = m4[2], var2 = m4[3], corr = val,
                                       corr_expr = trimws(m4[4]),
                                       cov  = NA_real_,
                                       cov_expr  = NA_character_,
                                       stringsAsFactors = FALSE))
      i <- i + 1; next
    }

    # Pattern: skew IDENT = expr  (skewness shape parameter alpha for shock IDENT)
    # Both a numeric snapshot and raw expression text (skew_expr) are stored so
    # .get_shock_skewness() can re-evaluate against the current parameter vector
    # during MCMC (mirrors the stderr_expr pattern at lines 469-474).
    # alpha can be any real number (negative allowed -- see brief Landmine 7).
    m5 <- regmatches(s, regexec(
      "^\\s*skew\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.+)$",
      s, perl = TRUE))[[1]]
    if (length(m5) > 0 && nchar(m5[1]) > 0) {
      vname <- m5[2]
      val   <- safe_eval(m5[3])
      # Find existing row for this shock and update its skew columns;
      # if no row exists yet (e.g. skew appears before var), create a stub.
      idx <- which(variances$name == vname)
      if (length(idx) > 0) {
        variances$skew[idx[1]]      <- val
        variances$skew_expr[idx[1]] <- trimws(m5[3])
      } else {
        variances <- rbind(variances,
                           data.frame(name = vname, stderr = 0,
                                      variance = 0,
                                      stderr_expr = NA_character_,
                                      variance_expr = NA_character_,
                                      skew = val,
                                      skew_expr = trimws(m5[3]),
                                      stringsAsFactors = FALSE))
      }
      i <- i + 1; next
    }

    i <- i + 1
  }

  list(variances = variances, correlations = correlations,
       deterministic = deterministic)
}


#' Parse the estimated_params block
#'
#' Dynare allows several syntaxes per line. The general form is
#'   NAME [, INITVAL [, LB, UB]] , PRIOR_SHAPE, P1, P2 [, P3, P4 [, JSCALE]] ;
#' where the leading INITVAL/LB/UB and the trailing P3/P4/JSCALE are optional.
#' Both the "short" form (no INITVAL block, as in the bundled models)
#'   alp, beta_pdf, 0.356, 0.02 ;
#' and the "long" form (with INITVAL, LB, UB, as in Smets-Wouters 2007)
#'   stderr ea, 0.4618, 0.01, 3, INV_GAMMA_PDF, 0.1, 2 ;
#' are supported. Pure-calibration forms with no prior shape
#'   NAME, INITVAL [, LB, UB] ;
#' are also accepted (prior left NA).
#'
#' The prior shape is located as the first non-numeric token after the
#' name(s); fields before it are INITVAL[, LB, UB] and fields after it are
#' the prior parameters P1, P2, P3, P4 (any trailing JSCALE is ignored).
#'
#' @param body Body text of the estimated_params block.
#' @return data.frame with columns: type, name, name2, prior, p1, p2, p3, p4,
#'   init, lb, ub.
#' @noRd
parse_estimated_params_block <- function(body) {
  empty <- data.frame(
    type  = character(0),
    name  = character(0),
    name2 = character(0),
    prior = character(0),
    p1 = numeric(0), p2 = numeric(0),
    p3 = numeric(0), p4 = numeric(0),
    init = numeric(0), lb = numeric(0), ub = numeric(0),
    stringsAsFactors = FALSE
  )

  stmts <- strsplit(body, ";")[[1]]
  stmts <- trimws(stmts)
  stmts <- stmts[nchar(stmts) > 0]

  rows <- list()
  num <- function(x) suppressWarnings(as.numeric(x))

  for (s in stmts) {
    if (grepl("^\\s*end\\s*$", s, ignore.case = TRUE)) next

    parts <- trimws(strsplit(s, ",")[[1]])

    etype  <- "parameter"
    ename  <- parts[1]

    # ---- M11: bare entry (name only, no prior / bounds) ----------------
    # e.g. "omega;" — treated as estimated parameter with unspecified prior.
    if (length(parts) == 1 && grepl("^[A-Za-z_][A-Za-z0-9_]*$", ename)) {
      rows[[length(rows) + 1L]] <- data.frame(
        type = etype, name = ename, name2 = NA_character_,
        prior = NA_character_,
        p1 = NA_real_, p2 = NA_real_, p3 = NA_real_, p4 = NA_real_,
        init = NA_real_, lb = NA_real_, ub = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    if (length(parts) < 2) next
    ename2 <- NA_character_
    value_parts <- parts[-1]

    if (grepl("^stderr\\s+", ename)) {
      etype <- "stderr"
      ename <- trimws(sub("^stderr\\s+", "", ename))
    } else if (grepl("^skew\\s+", ename)) {
      ## Two-token treatment: "skew SHOCKNAME" sets the skewness shape parameter
      ## alpha for shock SHOCKNAME (Dynare 7 convention; brief Landmine 7:
      ## alpha can be negative, so priors on (-Inf, Inf) are natural).
      etype <- "skew"
      ename <- trimws(sub("^skew\\s+", "", ename))
    } else if (grepl("^corr\\s+", ename)) {
      etype <- "corr"
      ename <- trimws(sub("^corr\\s+", "", ename))
      # corr takes two names: "corr a, b, prior, ..."  ->  name2 = parts[2]
      ename2 <- trimws(parts[2])
      value_parts <- parts[-(1:2)]
    }

    # Locate the prior shape: the first non-numeric value token.
    is_num <- !is.na(num(value_parts))
    shape_idx <- which(!is_num)[1]

    init <- NA_real_; lb <- NA_real_; ub <- NA_real_
    prior_name <- NA_character_
    p1 <- NA_real_; p2 <- NA_real_; p3 <- NA_real_; p4 <- NA_real_

    if (is.na(shape_idx)) {
      # No prior shape: pure calibration form NAME, INITVAL[, LB, UB].
      if (length(value_parts) >= 1) init <- num(value_parts[1])
      if (length(value_parts) >= 3) {
        lb <- num(value_parts[2]); ub <- num(value_parts[3])
      }
    } else {
      prior_name <- value_parts[shape_idx]
      before <- value_parts[seq_len(shape_idx - 1L)]
      after  <- if (shape_idx < length(value_parts))
        value_parts[(shape_idx + 1L):length(value_parts)] else character(0)

      # Fields before the shape: INITVAL [, LB, UB].
      if (length(before) >= 1) init <- num(before[1])
      if (length(before) >= 3) { lb <- num(before[2]); ub <- num(before[3]) }

      # Fields after the shape: P1, P2, P3, P4 (trailing JSCALE ignored).
      if (length(after) >= 1) p1 <- num(after[1])
      if (length(after) >= 2) p2 <- num(after[2])
      if (length(after) >= 3) p3 <- num(after[3])
      if (length(after) >= 4) p4 <- num(after[4])
    }

    rows[[length(rows) + 1L]] <- data.frame(
      type = etype, name = ename, name2 = ename2,
      prior = prior_name,
      p1 = p1, p2 = p2, p3 = p3, p4 = p4,
      init = init, lb = lb, ub = ub,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}


#' Parse the body of an occbin_constraints block
#'
#' Handles the Dynare-6 occbin_constraints block syntax:
#'
#'   occbin_constraints;
#'     name 'zlb'; bind r <= -0.004; [relax r > -0.004;]
#'       [equations  r = -0.004;  end;]
#'     name 'cap'; bind y >= 0.01;
#'   end;
#'
#' Each named constraint produces one spec in the returned list. var_idx and
#' eq_idx are NOT resolved here (they require a compiled model); resolution
#' happens in obc_resolve_block_specs() inside obc-spec.R.
#'
#' @param body Body text of the occbin_constraints block (between ; and end;).
#' @return List of raw constraint specs, each containing:
#'   $name       -- constraint name string (from name '...')
#'   $var_name   -- constrained variable name (from bind clause)
#'   $op         -- ">" (lower bound) or "<" (upper bound), normalised
#'   $bound      -- numeric bound value; NA if bound is a parameter reference
#'   $bound_expr -- raw bound expression string (parameter name or NA if numeric)
#'   $relax_str  -- relax clause text, or "" if absent
#'   $bind_eqs   -- equations sub-block body text, or "" if absent
#' @noRd
parse_occbin_constraints_block <- function(body) {
  specs <- list()

  # Split body into per-constraint chunks using `name` as the delimiter.
  # Each chunk spans from one `name '...' ;` to the next (or end of body).
  # Strategy: find all `name '...' ;` positions, then slice between them.
  name_pat <- "(?si)\\bname\\s+['\"]([^'\"]+)['\"]\\s*;"
  name_pos  <- gregexpr(name_pat, body, perl = TRUE)[[1]]
  if (name_pos[1] == -1L) return(specs)   # no constraints found

  name_starts  <- as.integer(name_pos)
  name_lengths <- attr(name_pos, "match.length")

  # Build chunk boundaries: each chunk runs from end-of-name-tag to start of
  # the next name tag (or end of body).
  n_constraints <- length(name_starts)
  chunk_ends    <- c(name_starts[-1L] - 1L, nchar(body))

  for (i in seq_len(n_constraints)) {
    # Extract the name
    name_match <- regmatches(
      substr(body, name_starts[i],
             name_starts[i] + name_lengths[i] - 1L),
      regexec(name_pat, substr(body, name_starts[i],
                               name_starts[i] + name_lengths[i] - 1L),
              perl = TRUE)
    )[[1]]
    constraint_name <- trimws(name_match[2])

    # Chunk of text that belongs to this constraint (after the name tag)
    chunk_start <- name_starts[i] + name_lengths[i]
    chunk <- substr(body, chunk_start, chunk_ends[i])

    # ---- bind clause: bind VAR OP EXPR ; -----------------------------------
    # EXPR can be: numeric literal, parameter name, or an arbitrary expression
    # (e.g. PHI*steady_state(iv)).  We match VAR and OP as structured fields
    # and capture EXPR as "everything between OP and ;" (trimmed).
    bind_m <- regmatches(chunk, regexec(
      "(?i)\\bbind\\s+(\\w+)\\s*([<>]=?)\\s*([^;]+);",
      chunk, perl = TRUE
    ))[[1]]

    if (length(bind_m) < 4) {
      warning(sprintf(
        "parse_occbin_constraints_block: constraint '%s' has no valid bind clause; skipping.",
        constraint_name))
      next
    }

    var_name  <- trimws(bind_m[2])
    op_raw    <- bind_m[3]
    bound_raw <- trimws(bind_m[4])
    # Try to parse as a pure numeric literal; if not, store the raw expression
    bound_val <- suppressWarnings(as.numeric(bound_raw))
    if (is.na(bound_val)) {
      # Expression or parameter name — store NA and keep the expression for
      # later resolution (Phase 2 bind-condition evaluation)
      bound_val <- NA_real_
      bound_expr <- bound_raw
    } else {
      bound_expr <- NA_character_
    }
    # bind VAR OP means "constraint activates when VAR OP bound" â€" the
    # constraint direction is the OPPOSITE: bind <= â†’ constraint >, etc.
    op       <- if (startsWith(op_raw, "<")) ">" else "<"

    # ---- relax clause (optional): relax VAR OP VALUE ; ----------------------
    # The relax clause mirrors the bind clause and may name a DIFFERENT variable
    # (e.g. RBC: bind iv, relax lam).  We first try the structured form
    # "VAR OP VALUE" where VALUE is numeric or a parameter name / steady_state()
    # expression.  If that fails we fall back to storing the raw string only.
    relax_m <- regmatches(chunk, regexec(
      "(?i)\\brelax\\s+([^;]+);", chunk, perl = TRUE
    ))[[1]]
    relax_str <- if (length(relax_m) >= 2) trimws(relax_m[2]) else ""

    # Parse the structured relax clause: VAR OP VALUE
    # VALUE can be: numeric literal, parameter name, or steady_state(EXPR)
    relax_var_name   <- NA_character_
    relax_op         <- NA_character_
    relax_bound_expr <- NA_character_
    if (nzchar(relax_str)) {
      relax_struct <- regmatches(relax_str, regexec(
        "^(\\w+)\\s*([<>]=?)\\s*(.+)$", relax_str, perl = TRUE
      ))[[1]]
      if (length(relax_struct) == 4L) {
        relax_var_name   <- trimws(relax_struct[2L])
        relax_op_raw     <- relax_struct[3L]
        relax_bound_expr <- trimws(relax_struct[4L])
        # Store the relax op as-is (not flipped — it is "constraint inactive when")
        relax_op         <- relax_op_raw
      }
    }

    # ---- equations sub-block (optional) -------------------------------------
    # Dynare uses "equations\n ... end;" without ; after the keyword,
    # so we use a custom regex rather than extract_paired_block.
    eq_m <- regmatches(chunk, regexec(
      "(?si)\\bequations\\b\\s*;?(.*?)\\bend\\s*;", chunk, perl = TRUE
    ))[[1]]
    bind_eqs <- if (length(eq_m) >= 2) trimws(eq_m[2]) else ""

    specs <- c(specs, list(list(
      name             = constraint_name,
      var_name         = var_name,
      op               = op,
      bound            = bound_val,
      bound_expr       = bound_expr,
      relax_str        = relax_str,
      relax_var_name   = relax_var_name,
      relax_op         = relax_op,
      relax_bound_expr = relax_bound_expr,
      bind_eqs         = bind_eqs
    )))
  }

  specs
}


#' Parse a command option string like "order=1, irf=40, nograph"
#'
#' @param options_str The text inside the parentheses.
#' @return Named list. Bare flags have value TRUE; key=value pairs are parsed.
#' @noRd
parse_command_options <- function(options_str) {
  if (is.null(options_str) || nchar(trimws(options_str)) == 0)
    return(list())

  opts <- list()
  parts <- trimws(strsplit(options_str, ",")[[1]])

  for (p in parts) {
    if (nchar(p) == 0) next
    if (grepl("=", p)) {
      kv <- strsplit(p, "\\s*=\\s*", perl = TRUE)[[1]]
      key <- trimws(kv[1])
      val_str <- trimws(kv[2])
      val <- suppressWarnings(as.numeric(val_str))
      if (is.na(val)) {
        if (tolower(val_str) == "true")  val <- TRUE
        else if (tolower(val_str) == "false") val <- FALSE
        else val <- val_str
      }
      opts[[key]] <- val
    } else {
      opts[[p]] <- TRUE
    }
  }
  opts
}


#' Translate a single MATLAB statement (from a verbatim; block) to R
#'
#' Conservative: handles the common scalar/matrix idioms that appear in the
#' parameter-setup verbatim blocks of replication models (matrix literals,
#' transpose, indexing, and a handful of element/matrix functions).  Returns
#' \code{NULL} when the statement uses syntax we do not understand (e.g.
#' \code{strmatch}, \code{M_.params}, control flow) so the caller can skip it.
#'
#' @param stmt A single MATLAB statement (no trailing semicolon).
#' @return R source text, or \code{NULL} if untranslatable.
#' @noRd
.matlab_stmt_to_r <- function(stmt) {
  s <- trimws(stmt)
  if (nchar(s) == 0) return(NULL)

  # Bail out on anything we deliberately do not support: Dynare globals,
  # string functions, control flow, ranges, anonymous functions, logical ops.
  # (Note: a bare "'" is NOT excluded here -- it is MATLAB transpose, handled
  # below.  String literals are caught after transpose conversion.)
  if (grepl("M_\\.|oo_\\.|options_\\.|strmatch|deblank|\\bfor\\b|\\bif\\b|\\bwhile\\b|\\bend\\b|@|\\.\\.\\.|:|&&|\\|\\||~", s))
    return(NULL)

  # Matrix literal:  [a b; c d]  ->  rbind(c(a,b), c(c,d))
  # Rows are split on ';' or newlines; entries on whitespace and/or commas.
  conv_matrix <- function(inner) {
    inner <- trimws(inner)
    rows <- strsplit(inner, "[;\n]")[[1]]
    rows <- trimws(rows)
    rows <- rows[nchar(rows) > 0]
    if (length(rows) == 0) return("numeric(0)")
    row_vecs <- vapply(rows, function(r) {
      # split on commas or runs of whitespace, but NOT inside nested parens
      toks <- strsplit(r, "(?<![,([{])[[:space:]]+|\\s*,\\s*", perl = TRUE)[[1]]
      toks <- toks[nchar(trimws(toks)) > 0]
      paste0("c(", paste(toks, collapse = ", "), ")")
    }, character(1))
    if (length(row_vecs) == 1L) return(row_vecs[[1]])
    paste0("rbind(", paste(row_vecs, collapse = ", "), ")")
  }
  # Replace every top-level [...] literal.  We scan left to right; nested
  # brackets are not used in these blocks.
  while (grepl("\\[", s)) {
    op <- regexpr("\\[", s)[1]
    # find matching close bracket (no nesting expected)
    cl <- regexpr("\\]", substr(s, op, nchar(s)))[1]
    if (cl == -1L) return(NULL)
    cl <- op + cl - 1L
    inner <- substr(s, op + 1L, cl - 1L)
    repl  <- conv_matrix(inner)
    s <- paste0(substr(s, 1L, op - 1L), repl, substr(s, cl + 1L, nchar(s)))
  }

  # MATLAB transpose  X'  ->  t(X)   (identifier, indexed name, or paren group
  # immediately followed by a single quote).  Iterate so chained forms resolve.
  repeat {
    s2 <- gsub("([A-Za-z_][A-Za-z0-9_]*(?:\\[[^]]*\\])?|\\))'",
               "t(\\1)", s, perl = TRUE)
    if (identical(s2, s)) break
    s <- s2
  }
  # Any quote that survives transpose conversion is a string literal we do not
  # support -- bail out rather than produce broken R.
  if (grepl("'", s, fixed = TRUE)) return(NULL)

  # eye(n)  / eye(n,n)  ->  diag(n)
  s <- gsub("eye\\(\\s*([0-9]+)\\s*(?:,\\s*[0-9]+\\s*)?\\)", "diag(\\1)",
            s, perl = TRUE)

  # Element/matrix functions map 1:1 to R: sqrt, diag, chol.
  # MATLAB chol() is UPPER-triangular like R's chol(); leave as-is.

  # Indexing:  X(i)  and  X(i,j)  for already-defined matrices/vectors.
  # We cannot statically tell a function call from indexing, so only rewrite
  # parens that follow a name which is NOT a known function.  Easiest robust
  # approach: rewrite ALL `name(args)` where name is not in a small builtin
  # whitelist into `name[args]`.  Builtins keep round parens.
  builtins <- c("sqrt", "diag", "chol", "exp", "log", "abs", "sum", "prod",
                "c", "rbind", "cbind", "t", "solve", "max", "min")
  # Rewrite MATLAB-style indexing name(args) -> name[args] (skip builtins),
  # via a manual paren-matching scan (gsub cannot match balanced parens).
  s <- .matlab_index_rewrite(s, builtins)

  # Operators.  MATLAB distinguishes matrix ops (* / ^) from element-wise ops
  # (.* ./ .^).  R's *, /, ^ are element-wise and %*% is matrix mult.  We map:
  #   .*  ./  .^   -> *  /  ^        (element-wise, direct)
  #   X^(-1)       -> solve(X)       (matrix inverse; operand balanced-matched)
  #   *            -> .mtimes(a, b)  (matrix-or-scalar-aware multiply helper)
  # The .mtimes helper (defined in the eval env) does %*% when both operands
  # are non-1x1 matrices, else element-wise *, so scalar*matrix still works.
  # First protect element-wise dot-operators with placeholders.
  s <- gsub("\\.\\*", "", s, perl = TRUE)
  s <- gsub("\\./",  "", s, perl = TRUE)
  s <- gsub("\\.\\^", "", s, perl = TRUE)

  # X^(-1) inverse: find each "^(-1)" and wrap the preceding balanced operand.
  s <- .matlab_inv_rewrite(s)

  # Remaining bare '*' is MATLAB matrix multiply.  Replace the binary operator
  # by a placeholder, then fold into .mtimes() pairwise (left-associative).
  s <- .matlab_mtimes_rewrite(s)

  # Restore element-wise operators.
  s <- gsub("", "*", s, fixed = TRUE)
  s <- gsub("", "/", s, fixed = TRUE)
  s <- gsub("", "^", s, fixed = TRUE)

  s
}


#' Rewrite MATLAB X^(-1) (matrix inverse) -> solve(X), balanced operand
#' @noRd
.matlab_inv_rewrite <- function(s) {
  repeat {
    m <- regexpr("\\^\\s*\\(\\s*-\\s*1\\s*\\)", s, perl = TRUE)
    if (m == -1L) break
    end_op <- m - 1L                         # last char of the operand
    # walk back over whitespace
    j <- end_op
    while (j >= 1 && substr(s, j, j) == " ") j <- j - 1L
    end_op <- j
    last <- substr(s, end_op, end_op)
    if (last == ")" || last == "]") {
      # balanced group: walk back to its opener
      openc <- if (last == ")") "(" else "["
      depth <- 0L; k <- end_op
      while (k >= 1) {
        ch <- substr(s, k, k)
        if (ch == last) depth <- depth + 1L
        else if (ch == openc) { depth <- depth - 1L; if (depth == 0L) break }
        k <- k - 1L
      }
      # include a leading function/index name if present (e.g. diag(...))
      start_op <- k
      p <- k - 1L
      while (p >= 1 && grepl("[A-Za-z0-9_.]", substr(s, p, p))) p <- p - 1L
      start_op <- p + 1L
    } else {
      # bare identifier/number operand
      p <- end_op
      while (p >= 1 && grepl("[A-Za-z0-9_.]", substr(s, p, p))) p <- p - 1L
      start_op <- p + 1L
    }
    operand <- substr(s, start_op, end_op)
    inv_end <- m + attr(m, "match.length") - 1L
    s <- paste0(substr(s, 1L, start_op - 1L),
                "solve(", operand, ")",
                substr(s, inv_end + 1L, nchar(s)))
  }
  s
}


#' Rewrite MATLAB matrix-multiply a*b*c -> .mtimes(.mtimes(a,b),c)
#' Splits the expression on top-level '*' (depth 0 w.r.t. ()/[]) and folds.
#' @noRd
.matlab_mtimes_rewrite <- function(s) {
  if (!grepl("\\*", s)) return(s)
  parts <- .split_top_level_multi(s, "*")
  if (length(parts) <= 1L) return(s)
  parts <- vapply(parts, .matlab_mtimes_rewrite, character(1))  # recurse
  Reduce(function(a, b) paste0(".mtimes(", a, ", ", b, ")"), parts)
}


#' Split on a single-char operator at bracket depth 0 (handles () and [])
#' @noRd
.split_top_level_multi <- function(s, op) {
  out <- character(0); cur <- ""; depth <- 0L
  for (k in seq_len(nchar(s))) {
    ch <- substr(s, k, k)
    if (ch == "(" || ch == "[") depth <- depth + 1L
    else if (ch == ")" || ch == "]") depth <- max(0L, depth - 1L)
    if (ch == op && depth == 0L) { out <- c(out, cur); cur <- "" }
    else cur <- paste0(cur, ch)
  }
  c(out, cur)
}


#' Rewrite MATLAB-style indexing name(args) -> name[args] (skip builtins)
#' @noRd
.matlab_index_rewrite <- function(s, builtins) {
  out <- ""
  i <- 1L
  n <- nchar(s)
  while (i <= n) {
    # match an identifier followed by '('
    rest <- substr(s, i, n)
    m <- regexec("^([A-Za-z_][A-Za-z0-9_]*)\\(", rest, perl = TRUE)
    cap <- regmatches(rest, m)[[1]]
    if (length(cap) >= 2) {
      name <- cap[2]
      # locate matching close paren
      open <- i + nchar(name)              # position of '('
      depth <- 0L
      j <- open
      close <- NA_integer_
      while (j <= n) {
        ch <- substr(s, j, j)
        if (ch == "(") depth <- depth + 1L
        else if (ch == ")") { depth <- depth - 1L; if (depth == 0L) { close <- j; break } }
        j <- j + 1L
      }
      if (is.na(close)) { out <- paste0(out, substr(s, i, n)); break }
      args <- substr(s, open + 1L, close - 1L)
      if (name %in% builtins) {
        # keep round parens; recurse into args for nested indexing
        out <- paste0(out, name, "(", .matlab_index_rewrite(args, builtins), ")")
      } else {
        # treat as indexing: name[args]
        out <- paste0(out, name, "[", .matlab_index_rewrite(args, builtins), "]")
      }
      i <- close + 1L
    } else {
      out <- paste0(out, substr(s, i, i))
      i <- i + 1L
    }
  }
  out
}


#' Split a string on a separator char, ignoring separators nested in brackets
#' @noRd
.split_top_level <- function(s, sep = ";", open = "[", close = "]") {
  out   <- character(0)
  cur   <- ""
  depth <- 0L
  for (k in seq_len(nchar(s))) {
    ch <- substr(s, k, k)
    if (ch == open)  depth <- depth + 1L
    else if (ch == close) depth <- max(0L, depth - 1L)
    if (ch == sep && depth == 0L) {
      out <- c(out, cur); cur <- ""
    } else {
      cur <- paste0(cur, ch)
    }
  }
  c(out, cur)
}


#' Evaluate verbatim; blocks for scalar parameter values
#'
#' Dynare \code{verbatim;} blocks contain raw MATLAB that the preprocessor
#' passes through untouched; replication models routinely use them to compute
#' derived parameter values (e.g. a VAR-implied mean \code{P0}, or shock
#' std-devs/correlations from a covariance matrix) that the model then
#' references.  dynhr otherwise strips these blocks entirely, leaving those
#' parameters \code{NA} (repl M19 -> e.g. \code{Sigma_e = 0}, BK violations).
#'
#' This evaluates each verbatim block's statements in a shared R environment
#' (seeded with the already-parsed scalar \code{param_values}), translating a
#' conservative subset of MATLAB (matrix literals, transpose, indexing,
#' \code{eye/sqrt/diag/chol}, scalar power, \code{M^(-1)} inverse).  Statements
#' it cannot translate or that error are skipped silently -- it never throws.
#' The returned environment also lets downstream top-level calibration
#' assignments (which may live OUTSIDE the verbatim block but reference its
#' intermediates, e.g. \code{sigma_z = sqrt(V(1,1))}) resolve.
#'
#' @param txt          Cleaned .mod text (comments stripped).
#' @param param_values Named numeric vector of already-known scalar params.
#' @return An environment containing every successfully evaluated verbatim
#'   binding (scalars and matrices), parented at \code{baseenv()}.
#' @noRd
# Matrix-or-scalar-aware multiply used by translated verbatim expressions.
# MATLAB '*' is matrix multiply; R '*' is element-wise.  Use %*% when both
# operands are genuine (non-1x1) matrices/vectors, else fall back to '*' so
# scalar*matrix and scalar*scalar still behave.
.mtimes <- function(a, b) {
  is_mat <- function(x) (is.matrix(x) && length(x) > 1L) ||
                        (is.numeric(x) && length(x) > 1L && is.null(dim(x)))
  if (is_mat(a) && is_mat(b)) {
    am <- if (is.matrix(a)) a else matrix(a, nrow = 1L)
    bm <- if (is.matrix(b)) b else matrix(b, ncol = 1L)
    if (ncol(am) == nrow(bm)) return(am %*% bm)
    # vector*vector with mismatched orientation: treat a as row, b as column
    return(as.numeric(a) %*% as.numeric(b))
  }
  a * b
}

eval_verbatim_blocks <- function(txt, param_values = numeric(0)) {
  env <- list2env(as.list(param_values), parent = baseenv())
  assign(".mtimes", .mtimes, envir = env)
  blocks <- extract_all_paired_blocks(txt, "verbatim")
  if (length(blocks) == 0) return(env)

  for (blk in blocks) {
    body <- gsub("%[^\n]*", "", blk$body)   # strip MATLAB % comments
    # Split into statements on ';' -- but NOT a ';' inside a [...] matrix
    # literal, where ';' is MATLAB's row separator (e.g. [a;b;c]).  Scan and
    # cut only at bracket-depth 0.
    stmts <- .split_top_level(body, sep = ";", open = "[", close = "]")
    for (stmt in stmts) {
      stmt <- trimws(stmt)
      if (nchar(stmt) == 0) next
      if (!grepl("=", stmt)) next            # only assignments
      # split LHS = RHS at the first top-level '=' that is not '=='/'<='/'>='
      eqpos <- regexpr("(?<![<>=~!])=(?!=)", stmt, perl = TRUE)
      if (eqpos == -1L) next
      lhs <- trimws(substr(stmt, 1L, eqpos - 1L))
      rhs <- trimws(substr(stmt, eqpos + 1L, nchar(stmt)))
      if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", lhs)) next   # scalar LHS name only
      r_rhs <- .matlab_stmt_to_r(rhs)
      if (is.null(r_rhs)) next
      val <- tryCatch(eval(parse(text = r_rhs), envir = env),
                      error = function(e) NULL)
      if (!is.null(val) && is.numeric(val))
        assign(lhs, val, envir = env)
    }
  }
  env
}
