## R/parse-blocks.R
## --------------------------------------------------------------------------
## Declaration and block extraction for .mod files: var / varexo /
## varexo_det / parameters declarations, shocks block, estimated_params,
## initval, endval, steady_state_model, and command-line options.
##
## Phase-1 split from parser-monolith.R (no logic changes).
## --------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# .mod EXPRESSION SANDBOX
# ---------------------------------------------------------------------------
# Parsing a .mod file EVALUATES R code that came out of that file: calibration
# assignments, initval/endval/histval right-hand sides, shocks-block variances,
# verbatim blocks, and macro directives are all `eval(parse(text = ...))`.
#
# Until 0.9.4 every one of those evaluation environments was a child of
# `baseenv()`, which means a .mod file could run ARBITRARY CODE simply by
# being parsed -- `alpha = system("curl ...");` or `beta = unlink("~", TRUE)`
# would execute inside `parse_mod()`.  Model files are routinely downloaded
# from replication archives, so "reading a model" must not be a code-execution
# primitive.
#
# The fix is one shared sandbox constructor.  Its environments are parented at
# `emptyenv()` and pre-populated with an EXPLICIT ALLOWLIST of functions, so an
# unlisted name is not reachable at all; on top of that every expression is
# AST-checked before it is evaluated, so a disallowed call fails LOUDLY with a
# classed condition naming the symbol instead of failing obscurely at eval time
# (or, worse, resolving to something harmless-looking).
#
# Three allowlists, because the three dialects need different vocabularies:
#   * `.dynhr_safe_fn_names`        arithmetic + elementary maths.  Calibration,
#                                   initval/endval, histval, shocks, and D33's
#                                   reduced-form re-evaluation.
#   * `.dynhr_safe_matrix_fn_names` + matrix constructors/algebra, for the
#                                   MATLAB-translated `verbatim` blocks.
#   * `.dynhr_safe_macro_fn_names`  + comparison/logical/membership operators.
#                                   The `@#` macro language no longer uses it:
#                                   since W30 (2026-09-25) macro expressions
#                                   are run by dynhr's own interpreter
#                                   (R/parse-macro.R), which evaluates no R
#                                   code at all.  Kept as a vetted allowlist.
# None of them contains a function that touches the filesystem, the network,
# the process table, or any environment outside the sandbox.

.dynhr_safe_fn_names <- c(
  "(", "+", "-", "*", "/", "^", "%%", "%/%",
  "exp", "log", "log10", "log2", "log1p", "expm1", "sqrt", "abs", "sign",
  "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
  "sinh", "cosh", "tanh", "min", "max", "sum", "prod", "mean",
  "round", "floor", "ceiling", "trunc",
  "gamma", "lgamma", "beta", "lbeta", "digamma", "factorial", "choose",
  "pnorm", "qnorm", "dnorm"
)

.dynhr_safe_matrix_fn_names <- c(
  .dynhr_safe_fn_names,
  "c", "[", "[[", "matrix", "rbind", "cbind", "t", "solve", "diag", "chol",
  "crossprod", "tcrossprod", "kronecker", "%*%", "%o%",
  "nrow", "ncol", "dim", "length", "numeric", "rep", "seq_len", "as.numeric",
  "det", "outer", ".mtimes"
)

.dynhr_safe_macro_fn_names <- c(
  .dynhr_safe_fn_names,
  "c", "[", "length", "%in%", "paste0", "paste", "nchar",
  "==", "!=", "<", ">", "<=", ">=", "&", "|", "&&", "||", "!", "xor",
  "isTRUE", "isFALSE", "any", "all", ":"
)

#' Build an allowlisted sandbox environment for evaluating .mod expressions
#'
#' @param seed Optional named list / named numeric vector / environment whose
#'   bindings seed the sandbox (already-known parameter values).
#' @param fns  Character vector of allowed function names.
#' @return An environment whose parent is `emptyenv()`.
#' @noRd
.dynhr_sandbox_env <- function(seed = NULL, fns = .dynhr_safe_fn_names) {
  env <- new.env(parent = emptyenv())
  for (fn in fns) {
    if (identical(fn, ".mtimes")) {
      assign(".mtimes", .mtimes, envir = env)
      next
    }
    ## A handful of the allowlisted maths functions live in stats, not base.
    src <- if (exists(fn, envir = baseenv(), inherits = FALSE)) baseenv()
           else asNamespace("stats")
    assign(fn, get(fn, envir = src), envir = env)
  }
  assign("pi", base::pi, envir = env)
  if (is.environment(seed)) seed <- as.list(seed)
  if (is.numeric(seed) || is.character(seed)) seed <- as.list(seed)
  if (is.list(seed) && length(seed) > 0L && !is.null(names(seed)))
    for (nm in names(seed)) assign(nm, seed[[nm]], envir = env)
  env
}

#' Function names called by an expression (AST, function positions only)
#'
#' Walks the AST and collects the head of every call. The earlier
#' `setdiff(all.names(e, functions = TRUE), all.vars(e))` dropped a function
#' name that ALSO appeared as a variable (`system("x") + system`), so the
#' allowlist check never saw it. A non-symbol call head (`f()()`,
#' `(function(x) x)(1)`) is reported as "<call>" so it can never pass.
#' @noRd
.dynhr_expr_calls <- function(e) {
  if (is.call(e)) {
    h <- e[[1L]]
    head <- if (is.symbol(h)) as.character(h) else "<call>"
    return(unique(c(head, unlist(lapply(as.list(e), .dynhr_expr_calls)))))
  }
  if (is.expression(e) || is.list(e))
    return(unique(unlist(lapply(as.list(e), .dynhr_expr_calls))))
  character(0)
}

#' Evaluate one .mod expression inside the sandbox
#'
#' Security contract: the expression is parsed, its call names are checked
#' against `allowed`, and only then is it evaluated.  A call to anything
#' outside the allowlist ABORTS with a `dynhr_error_unsafe_mod_expression`
#' condition naming the offending symbol -- it is never evaluated.
#'
#' Resolution contract (unchanged from the pre-0.9.4 behaviour): text that is
#' not valid R, or that references a symbol not yet bound in `envir`, yields
#' `NULL`.  Both are legitimate during parsing (MATLAB-flavoured right-hand
#' sides, forward references), so they are reported by the caller's own
#' fail-loud pass rather than aborting the parse here.
#'
#' @param text     Expression text from the .mod file.
#' @param envir    Sandbox environment (see `.dynhr_sandbox_env()`).
#' @param allowed  Allowlisted function names.
#' @param context  Short label naming the block, for the error message.
#' @return The evaluated value, or `NULL` when it could not be resolved.
#' @noRd
.dynhr_sandbox_eval <- function(text, envir, allowed = .dynhr_safe_fn_names,
                                context = ".mod expression") {
  ## parse() is the one step that can legitimately fail on text that never was
  ## R (MATLAB matrix literals, Dynare string options). There is no
  ## non-signalling parser in base R, so this single handler stays; everything
  ## after it is a plain conditional.
  parsed <- tryCatch(parse(text = text, keep.source = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed) || length(parsed) != 1L) return(NULL)
  e <- parsed[[1L]]

  bad <- setdiff(.dynhr_expr_calls(e), allowed)
  if (length(bad) > 0L)
    .dynhr_abort(
      "dynhr refuses to evaluate ", context, " `", text, "`: it calls `",
      bad[1L], "`, which is not on the .mod expression allowlist. Parsing a ",
      ".mod file must never execute arbitrary code; move any computation that ",
      "needs `", bad[1L], "` out of the model file.",
      class = "dynhr_error_unsafe_mod_expression")

  syms <- setdiff(all.vars(e), "pi")
  for (s in syms) if (!exists(s, envir = envir)) return(NULL)

  eval(e, envir = envir)
}

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
    "\\s*", .dynhr_opts_re, "\\s*;",  # optional (options), bracket-aware;
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
    "\\s*", .dynhr_opts_re, "\\s*;",   # optional (options), bracket-aware;
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
  open_pat <- paste0("(?si)\\b", keyword, "\\s*", .dynhr_opts_re, "\\s*;")
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

  ## Options are BALANCED parentheses (group 1, PCRE recursion): a
  ## `var(deflator = A^(1/(1-alpha))) Y;` used to stop at the first `)`,
  ## fail the match and silently drop the whole declaration.
  pat <- paste0("(?si)", kw_pat,
                "(?:\\s*(\\((?:[^()]++|(?1))*\\)))?\\s+(.*?)\\s*;")

  all_positions <- gregexpr(pat, txt, perl = TRUE)
  all_matches   <- regmatches(txt, all_positions)[[1]]

  # -- DIAGNOSTIC --
  if (debug){
    .dynhr_cat(sprintf("  [extract_declaration] keyword='%s'  matches=%d\n", keyword, length(all_matches)))
    for (j in seq_along(all_matches)) {
      .dynhr_cat(sprintf("    match %d (nchar=%d): %.120s\n", j, nchar(all_matches[j]),
                  gsub("\\s+", " ", all_matches[j])))
    }
  }
  # ----------------

  if (length(all_matches) == 0) return("")

  bodies <- vapply(all_matches, function(m) {
    parts <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
    if (length(parts) >= 3) trimws(parts[3]) else ""
  }, character(1))

  result <- paste(bodies[nchar(bodies) > 0], collapse = " ")

  # -- DIAGNOSTIC --
  if (debug) .dynhr_cat(sprintf("    -> extracted names text: %.200s\n", gsub("\\s+", " ", result)))
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
  ## `\b` on BOTH sides: without the trailing one, `steady` matched the start
  ## of `steady_state_model;` (and `steady_state(x)` in an equation), giving a
  ## bogus `steady` command whose var list was `_state_model`, which
  ## write_mod() then emitted as `steady _state_model;` (review 2026-09-25).
  ## The option list may itself hold parenthesised lists
  ## (`instruments=(i,tau)`, `irf_shocks=(e,u)`): match balanced parentheses
  ## (PCRE recursion into group 1), not "up to the first `)`", which cut the
  ## options at the first inner `)` and pushed the rest into the var list.
  pat <- paste0(
    "(?i)\\b", command, "\\b",
    "\\s*", .dynhr_opts_re,                    # optional (options)
    "\\s*([^;]*?)\\s*;"                          # optional var list, then ;
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
  # ramsey_policy(...instruments=(a b c)...) or instruments=(a,b,c), read
  # through the command-option parser: a regex stopping at the first `)`
  # missed `instruments` whenever an earlier option held a list
  # (`ramsey_policy(irf_shocks=(e,u), instruments=(i))`).
  for (cmd in c("ramsey_policy", "ramsey_model")) {
    ci <- extract_command(txt, cmd)
    if (!ci$found) next
    ins <- parse_command_options(ci$options_str)$instruments
    if (!is.null(ins)) {
      ins <- as.character(ins)
      return(ins[nzchar(ins)])
    }
  }
  character(0)
}


#' Remove all recognised blocks from .mod text (to expose top-level calibration)
#'
#' @param txt Cleaned .mod text.
#' @return Text with all blocks replaced by whitespace.
#' @noRd
remove_blocks <- function(txt) {
  # Paired blocks: the SAME list the declaration scan strips (parse-lexer.R).
  txt <- .dynhr_strip_paired_blocks(txt)
  decl_kw <- c("var", "varexo_det", "varexo", "varobs", "parameters",
               "predetermined_variables", "trend_var", "log_trend_var",
               "model_local_variable", "var_expectation",
               "var_remove", "model_remove")
  # v0.3: handle optional parenthesised options, e.g. var(log) x y;
  for (kw in decl_kw) {
    pat <- paste0("(?si)\\b", kw, "\\b\\s*", .dynhr_opts_re, "\\s+.*?;")
    txt <- gsub(pat, " ", txt, perl = TRUE)
  }
  ## `model_remove('eq');` has no space before its `;`, so the pattern above
  ## misses it; its quoted tags may contain `)`.
  txt <- gsub("(?si)\\bmodel_remove\\b\\s*\\((?:'[^']*'|\"[^\"]*\"|[^)'\"])*\\)\\s*;",
              " ", txt, perl = TRUE)
  cmd_kw <- c("stoch_simul", "estimation", "steady", "check",
              "model_diagnostics", "model_info", "simul",
              "perfect_foresight_setup", "perfect_foresight_solver",
              "forecast", "conditional_forecast",
              "smoother2histval", "shock_decomposition",
              "calib_smoother", "extended_path",
              "osr", "ramsey_model", "ramsey_policy",
              "discretionary_policy", "planner_objective",
              "dynare_sensitivity", "bvar_density",
              "bvar_forecast", "dsample")
  for (kw in cmd_kw) {
    pat <- paste0("(?i)\\b", kw, "\\b\\s*", .dynhr_opts_re, "\\s*[^;]*;")
    txt <- gsub(pat, " ", txt, perl = TRUE)
  }
  # Option-only commands (no variable list) that were formerly mis-listed as
  # `...; end;` blocks -- which deleted everything up to the next `end;`.
  for (kw in c("identification", "pac_model", "var_model",
               "trend_component_model")) {
    pat <- paste0("(?i)\\b", kw, "\\b\\s*", .dynhr_opts_re, "\\s*;")
    txt <- gsub(pat, " ", txt, perl = TRUE)
  }

  ## `Sigma_e` is CASE-SENSITIVE, unlike every keyword above.
  ##
  ## Dynare's shock-covariance command is spelled `Sigma_e` (capital S) and
  ## Dynare identifiers are case-sensitive, so `sigma_e` is NOT that command --
  ## it is an ordinary, and extremely natural, parameter name. Matching it
  ## case-insensitively with the rest of `cmd_kw` deleted the user's
  ## `sigma_e = 1;` calibration line before parse_calibration() ever saw it,
  ## leaving the parameter declared but NA. Reported against 0.9.1: the model
  ## parsed cleanly after renaming the parameter to `sigma_obs_e`, which is
  ## the signature of a name collision rather than a syntax error.
  ##
  ## The other keywords stay case-insensitive because they are genuine Dynare
  ## reserved words -- a user cannot name a parameter `stoch_simul` in Dynare
  ## either, so being permissive there costs nothing. `Sigma_e` is the only
  ## entry whose lower-case form is a legal user identifier.
  txt <- gsub("\\bSigma_e\\b\\s*(?:\\([^)]*\\))?\\s*[^;]*;", " ", txt,
              perl = TRUE)
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


#' Split top-level .mod text into `NAME = RHS` calibration assignments
#'
#' Statements are delimited by `;`, NOT by newlines, so an expression split
#' over lines (`alpha = 0.3 +\n 0.03;`) is one assignment. The old
#' single-line regex silently dropped such a statement and the resulting
#' "no value" warning then blamed an external *_steadystate.m (review
#' 2026-09-25 B13).
#'
#' Within one `;`-terminated statement:
#' 1. if its LAST line contains `IDENT = rhs`, that is the assignment (the old
#'    single-line behaviour -- this also skips an unterminated junk line such
#'    as `title_string='foo'` sitting above a real `beta = 0.99;`);
#' 2. otherwise the assignment starts at the last line that BEGINS with
#'    `IDENT =` and its RHS runs to the `;`, newlines included.
#' Text after the final `;` is unterminated and ignored, as before.
#'
#' @param txt Top-level .mod text (blocks, declarations, commands removed).
#' @return List of length-2 character vectors `c(name, rhs)`, in source order.
#' @noRd
.calibration_assignments <- function(txt) {
  pieces <- strsplit(paste0(txt, "\n"), ";", fixed = TRUE)[[1]]
  pieces <- pieces[-length(pieces)]          # unterminated remainder
  out <- list()
  for (st in pieces) {
    lines <- strsplit(st, "\n", fixed = TRUE)[[1]]
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) == 0L) next
    last <- lines[length(lines)]
    m <- regmatches(last, regexec("([A-Za-z_][A-Za-z0-9_]*)\\s*=(?!=)(.*)$",
                                  last, perl = TRUE))[[1]]
    if (length(m) == 3L && nzchar(trimws(m[3]))) {
      out[[length(out) + 1L]] <- c(m[2], trimws(m[3]))
      next
    }
    starts <- gregexpr("(?m)^[ \\t]*[A-Za-z_][A-Za-z0-9_]*\\s*=(?!=)", st,
                       perl = TRUE)[[1]]
    if (starts[1L] < 0L) next
    frag <- substring(st, starts[length(starts)])
    m <- regmatches(frag, regexec("(?s)^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=(.*)$",
                                  frag, perl = TRUE))[[1]]
    rhs <- trimws(m[3])
    if (nzchar(rhs)) out[[length(out) + 1L]] <- c(m[2], rhs)
  }
  out
}


#' Parse top-level parameter assignments outside any block
#'
#' Matches statements like:  alpha = 0.33;  beta = 1/(1+0.02/4);
#' Statements are split on `;`, so an expression may span several lines.
#'
#' @param txt         Cleaned .mod text with blocks removed.
#' @param param_names Declared parameter names.
#' @return Named numeric vector of parameter values.
#' @noRd
parse_calibration <- function(txt, param_names, seed_env = NULL, quiet = FALSE) {
  values <- numeric(0)

  # Statements of the form  IDENT = expression ;  -- split on `;`, so an
  # expression may span lines (see .calibration_assignments()).
  all_assign <- .calibration_assignments(txt)

  # `seed_env` (optional) carries intermediate values computed in verbatim;
  # blocks (matrices like V, Correlation_matrix, and scalars like P0_z_bar0)
  # so that top-level assignments referencing them -- e.g.
  # `sigma_z = sqrt(V(1,1));` -- resolve instead of erroring to NA (repl M19).
  # MATLAB-style indexing in such RHS expressions is translated to R below.
  ## A2: the calibration environment is an ALLOWLIST SANDBOX parented at
  ## emptyenv(), not a child of baseenv() -- a .mod file must not be able to
  ## call system()/unlink()/file.remove() merely by being parsed.
  .pcal_env <- .dynhr_sandbox_env(seed_env, .dynhr_safe_matrix_fn_names)

  for (a in all_assign) {
    name <- a[[1L]]
    expr_text <- a[[2L]]

    # v0.3: collapse multiline expressions to a single line so that
    # R's parse() sees e.g. "(1 - theta) * (1 - beta*theta) / theta"
    # instead of "(1 - theta)\n      * (1 - beta*theta)\n      / theta"
    expr_text <- gsub("\\s+", " ", trimws(expr_text))

    # Evaluate expression in the calibration environment.
    # I2: evaluate ALL assignments in source order (not just declared params)
    # so that undeclared intermediate names (e.g. `alpha`, `theta`, `omega` in
    # Adam-Billi 2006) accumulate in .pcal_env before the derived declared
    # param expressions that reference them are reached.
    val <- .dynhr_sandbox_eval(expr_text, .pcal_env,
                               .dynhr_safe_matrix_fn_names,
                               context = "the calibration assignment")

    # M19: when a verbatim seed_env is in play, the RHS may use MATLAB
    # indexing (e.g. `sqrt(V(1,1))`).  Retry with a MATLAB->R translation.
    if (is.null(val) && !is.null(seed_env)) {
      r_expr <- .matlab_stmt_to_r(expr_text)
      if (!is.null(r_expr))
        val <- .dynhr_sandbox_eval(r_expr, .pcal_env,
                                   .dynhr_safe_matrix_fn_names,
                                   context = "the calibration assignment")
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
      for (a in rev(all_assign)) {
        if (a[[1L]] != pn) next
        rhs <- gsub("\\s+", " ", trimws(a[[2L]]))
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
    .dynhr_warn(sprintf(
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
    val <- .dynhr_sandbox_eval(expr_text, env, .dynhr_safe_matrix_fn_names,
                               context = "the initval/endval assignment")
    if (is.numeric(val) && length(val) == 1) {
      values[name] <- val
      assign(name, val, envir = env)
    }
  }
  values
}


#' Parse a histval block into a per-variable history over lags
#'
#' Dynare's \code{histval} block sets the pre-sample history of endogenous
#' (and exogenous) variables.  The index in \code{name(index)} is the period
#' RELATIVE TO THE FIRST SIMULATION PERIOD, and it is \code{<= 0}:
#' \code{y(0)} is the last pre-sample period (i.e. lag 1 as seen from period
#' 1), \code{y(-1)} is the one before that (lag 2), and so on -- exactly the
#' convention of the Dynare manual's histval example.  A bare
#' \code{y = value;} is accepted as shorthand for \code{y(0) = value;}.
#'
#' The parsed result flips that into a LAG index, which is what every
#' downstream consumer wants: element \code{k} of \code{histval[[v]]} is the
#' value of \code{v} \code{k} periods before the first simulation period, so
#' \code{histval[[v]][1]} is the initial condition \code{y0} of a one-lag
#' model.  Missing intermediate lags are \code{NA}.
#'
#' @param body Body text of the histval block.
#' @param env  Environment used to evaluate the right-hand sides (usually
#'   seeded with the calibrated parameter values).
#' @return Named list: variable name -> numeric vector indexed by lag.
#' @noRd
parse_histval_block <- function(body, env = parent.frame()) {
  out <- list()
  pat <- paste0("([A-Za-z_][A-Za-z0-9_]*)\\s*",       # name
                "(?:\\(\\s*([+-]?[0-9]+)\\s*\\))?",   # optional (index)
                "\\s*=\\s*([^;\\n]+)\\s*;")
  all_matches <- regmatches(body, gregexpr(pat, body, perl = TRUE))[[1]]

  for (m in all_matches) {
    parts <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
    name  <- parts[2]
    idx_s <- parts[3]
    idx   <- if (is.na(idx_s) || !nzchar(idx_s)) 0L else as.integer(idx_s)
    if (idx > 0L)
      stop(sprintf(paste0("histval: '%s(%+d)' uses a positive period index. ",
                          "histval indexes periods at or before the first ",
                          "simulation period, so the index must be <= 0 ",
                          "(y(0) is lag 1, y(-1) is lag 2)."),
                   name, idx), call. = FALSE)
    lag <- 1L - idx
    expr_text <- gsub("\\s+", " ", trimws(parts[4]))
    val <- .dynhr_sandbox_eval(expr_text, env, .dynhr_safe_matrix_fn_names,
                               context = "the histval assignment")
    if (!is.numeric(val) || length(val) != 1L)
      stop(sprintf(paste0("histval: could not evaluate the value for '%s' ",
                          "('%s') to a single number."), name, expr_text),
           call. = FALSE)
    v <- out[[name]]
    if (is.null(v)) v <- numeric(0)
    if (length(v) < lag) v <- c(v, rep(NA_real_, lag - length(v)))
    v[lag] <- val
    out[[name]] <- v
  }
  out
}


#' Parse a shock_groups block body into a named list of shock memberships
#'
#' Dynare syntax (each entry ends with a semicolon):
#' \preformatted{
#'   shock_groups(name = groupname);
#'   'supply' = e_a, e_z;
#'   'demand' = e_g;
#'   end;
#' }
#' The group label is normally single-quoted; double quotes and a bare
#' identifier are accepted too.  Members may be separated by commas or
#' whitespace.  The result is the named list of character vectors that
#' \code{\link{historical_decomposition}} takes as \code{shock_groups}.
#'
#' Unknown shock names are a hard error naming the offender: a mistyped shock
#' would otherwise silently drop out of the decomposition and every remaining
#' contribution would still add up, so nothing downstream could catch it.
#'
#' @param body        Body text of the shock_groups block.
#' @param shock_names Declared exogenous names used to validate membership.
#'   \code{NULL} skips validation.
#' @param block_name  Name of the enclosing block (for error messages).
#' @return Named list: group label -> character vector of shock names.
#' @noRd
parse_shock_groups_block <- function(body, shock_names = NULL,
                                     block_name = NULL) {
  where <- if (is.null(block_name)) "shock_groups"
           else sprintf("shock_groups(name = %s)", block_name)

  out <- list()
  ## Group label: 'quoted', "quoted" or a bare identifier; then = members ;
  pat <- paste0("(?:'([^']*)'|\"([^\"]*)\"|([A-Za-z_][A-Za-z0-9_]*))",
                "\\s*=\\s*([^;]*);")
  all_matches <- regmatches(body, gregexpr(pat, body, perl = TRUE))[[1]]

  for (m in all_matches) {
    parts <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
    label <- parts[2]
    if (!nzchar(label)) label <- parts[3]
    if (!nzchar(label)) label <- parts[4]
    label <- trimws(label)
    if (!nzchar(label))
      stop(where, ": a group entry has an empty name.", call. = FALSE)

    members <- strsplit(trimws(gsub(",", " ", parts[5])), "\\s+")[[1]]
    members <- members[nzchar(members)]
    if (length(members) == 0L)
      stop(sprintf("%s: group '%s' lists no shocks.", where, label),
           call. = FALSE)

    bad <- members[!grepl("^[A-Za-z_][A-Za-z0-9_]*$", members)]
    if (length(bad))
      stop(sprintf("%s: group '%s' has non-identifier member(s): %s.",
                   where, label, paste(bad, collapse = ", ")), call. = FALSE)

    if (!is.null(shock_names)) {
      unknown <- setdiff(members, shock_names)
      if (length(unknown))
        stop(sprintf(paste0("%s: group '%s' references unknown shock(s): %s. ",
                            "Declared shocks: %s."),
                     where, label, paste(unknown, collapse = ", "),
                     paste(shock_names, collapse = ", ")), call. = FALSE)
    }

    if (!is.null(out[[label]]))
      stop(sprintf("%s: group '%s' is defined twice.", where, label),
           call. = FALSE)
    out[[label]] <- members
  }

  dup <- unique(unlist(out, use.names = FALSE)[
    duplicated(unlist(out, use.names = FALSE))])
  if (length(dup))
    stop(sprintf(paste0("%s: shock(s) %s appear in more than one group; ",
                        "groups must partition the shocks."),
                 where, paste(dup, collapse = ", ")), call. = FALSE)

  out
}


#' Extract and parse every shock_groups block in a .mod file
#'
#' Dynare allows several \code{shock_groups} blocks, each labelled by its
#' \code{name=} option, so a file can carry alternative groupings and pick one
#' per \code{shock_decomposition} call.  A block with no \code{name=} option is
#' keyed \code{"default"}, matching Dynare's own default group name.
#'
#' @param txt         Cleaned .mod text.
#' @param shock_names Declared exogenous names (validation; \code{NULL} skips).
#' @return Named list: block name -> named list of group -> shock names.
#'   Empty list when the file has no \code{shock_groups} block.
#' @noRd
parse_shock_groups <- function(txt, shock_names = NULL) {
  blocks <- extract_all_paired_blocks(txt, "shock_groups")
  if (length(blocks) == 0L) return(list())

  out <- list()
  for (blk in blocks) {
    nm <- regmatches(
      blk$options_str,
      regexec("(?i)\\bname\\s*=\\s*['\"]?([A-Za-z_][A-Za-z0-9_]*)['\"]?",
              blk$options_str, perl = TRUE))[[1]]
    key <- if (length(nm) >= 2L) nm[2] else "default"
    if (!is.null(out[[key]]))
      stop(sprintf(paste0("parse_mod: two shock_groups blocks share the name ",
                          "'%s'; give each block a distinct name= option."),
                   key), call. = FALSE)
    out[[key]] <- parse_shock_groups_block(blk$body, shock_names, key)
  }
  out
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


#' Parse an observation_trends block into per-observable slope expressions
#'
#' Dynare syntax (one entry per observed variable, each ending with `;`):
#' \preformatted{
#'   observation_trends;
#'   dy (ctrend/100);
#'   pinfobs (0);
#'   end;
#' }
#' The expression is the SLOPE of a deterministic linear trend added to the
#' observed variable's measurement equation; it may use declared parameters
#' (and so be estimated). It is kept as TEXT and re-evaluated against the
#' parameter vector in force, exactly as the shocks-block expressions are.
#'
#' Dynare's preprocessor rejects an entry for a variable that is not observed
#' ("variable ... in observation_trends block is not an observed variable");
#' so does this, whenever a `varobs` list exists. Without one, the variable
#' must at least be a declared endogenous variable.
#'
#' The \code{deterministic_trends} block has the same entry syntax
#' (\code{block = "deterministic_trends"}); its entries may name any declared
#' endogenous variable, so it is called with \code{varobs_names =
#' character(0)}.
#'
#' @param body          Body text of the block.
#' @param var_names     Declared endogenous variable names.
#' @param param_names   Declared parameter names.
#' @param varobs_names  Declared observables (`character(0)` if none).
#' @param block         Block name used in error messages.
#' @return Named character vector: observable -> slope expression text.
#' @noRd
parse_observation_trends_block <- function(body, var_names, param_names,
                                           varobs_names = character(0),
                                           block = "observation_trends") {
  bad_syntax <- function(...)
    .dynhr_abort("parse_mod: ", block, ": ", ...,
                 class = "dynhr_error_mod_syntax")
  stmts <- trimws(strsplit(body, ";", fixed = TRUE)[[1]])
  stmts <- stmts[nzchar(stmts)]
  out <- character(0)
  for (st in stmts) {
    st <- gsub("\\s+", " ", st)
    m <- regmatches(st, regexec("^([A-Za-z_][A-Za-z0-9_]*)\\s*\\((.*)\\)$",
                                st, perl = TRUE))[[1]]
    if (length(m) != 3L || !nzchar(trimws(m[3])))
      bad_syntax("cannot read the entry `", st, "`; the form is ",
                 "`VARIABLE (EXPRESSION);`, e.g. `dy (ctrend/100);`.")
    nm   <- m[2]
    expr <- trimws(m[3])
    if (nm %in% names(out))
      bad_syntax("variable ", nm, " has more than one entry.")
    if (length(varobs_names) > 0L && !(nm %in% varobs_names))
      bad_syntax("variable ", nm, " in observation_trends block is not an ",
                 "observed variable (declared varobs: ",
                 paste(varobs_names, collapse = ", "), ").")
    if (!(nm %in% var_names))
      bad_syntax("variable ", nm, " is not a declared endogenous variable.")
    ## A number like 1e-3 yields the identifier-looking fragment "e"; strip
    ## numeric literals before collecting identifiers.
    expr_nonum <- gsub("(?<![A-Za-z0-9_])(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?", " ",
                       expr, perl = TRUE)
    ids <- regmatches(expr_nonum, gregexpr("[A-Za-z_][A-Za-z0-9_]*",
                                           expr_nonum, perl = TRUE))[[1]]
    unknown <- setdiff(unique(ids), c(param_names, .KNOWN_FUNCTIONS, "pi"))
    if (length(unknown) > 0L)
      bad_syntax("the trend of ", nm, " (`", expr, "`) uses ",
                 paste(unknown, collapse = ", "), ", which ",
                 if (length(unknown) == 1L) "is not a" else "are not",
                 " declared parameter", if (length(unknown) == 1L) "" else "s",
                 ". A trend slope may only depend on parameters.")
    out[nm] <- expr
  }
  out
}


#' Identifiers used by a .mod expression (numeric literals stripped first)
#' @noRd
.osr_expr_ids <- function(expr) {
  expr_nonum <- gsub("(?<![A-Za-z0-9_])(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?",
                     " ", expr, perl = TRUE)
  unique(regmatches(expr_nonum, gregexpr("[A-Za-z_][A-Za-z0-9_]*",
                                         expr_nonum, perl = TRUE))[[1]])
}


#' Parse Dynare's optim_weights block (OSR loss weights)
#'
#' Dynare syntax, one entry per statement:
#' \preformatted{
#'   optim_weights;
#'     pie 1;          // weight on Var(pie)
#'     y 0.5;          // weight on Var(y)
#'     y, pie 0.2;     // weight on Cov(y, pie)
#'   end;
#' }
#' The Dynare 7.1 preprocessor writes each entry into ONE cell of
#' \code{M_.osr.variable_weights} (a cross entry \code{y, pie w} sets cell
#' (y, pie) only) and the loss is \code{sum(W(:) .* vx(:))}, so a cross entry
#' adds \code{w * Cov(y, pie)} to the loss exactly once.  The weights are
#' parameter expressions, kept as text and evaluated when \code{osr()} runs.
#'
#' @param body Block body text.
#' @param var_names Declared endogenous variables.
#' @param param_names Declared parameters.
#' @return data.frame(var1, var2, expr); \code{var2 == var1} for a variance
#'   weight.
#' @noRd
parse_optim_weights_block <- function(body, var_names, param_names) {
  bad_syntax <- function(...)
    .dynhr_abort("parse_mod: optim_weights: ", ...,
                 class = "dynhr_error_mod_syntax")
  stmts <- trimws(strsplit(body, ";", fixed = TRUE)[[1]])
  stmts <- stmts[nzchar(stmts)]
  v1 <- character(0); v2 <- character(0); ex <- character(0)
  for (st in stmts) {
    st <- gsub("\\s+", " ", st)
    m <- regmatches(st, regexec(
      "^([A-Za-z_][A-Za-z0-9_]*)\\s*(?:,\\s*([A-Za-z_][A-Za-z0-9_]*))?\\s+(\\S.*)$",
      st, perl = TRUE))[[1]]
    if (length(m) != 4L)
      bad_syntax("cannot read the entry `", st, "`; the forms are ",
                 "`VARIABLE EXPRESSION;` and `VARIABLE, VARIABLE EXPRESSION;`.")
    a <- m[2]; b <- if (nzchar(m[3])) m[3] else m[2]; expr <- trimws(m[4])
    for (nm in unique(c(a, b)))
      if (!(nm %in% var_names))
        bad_syntax("`", nm, "` is not a declared endogenous variable.")
    unknown <- setdiff(.osr_expr_ids(expr),
                       c(param_names, .KNOWN_FUNCTIONS, "pi"))
    if (length(unknown) > 0L)
      bad_syntax("the weight `", expr, "` of ", a,
                 if (b != a) paste0(", ", b) else "", " uses ",
                 paste(unknown, collapse = ", "),
                 ", which is not a declared parameter. A weight may only ",
                 "depend on parameters.")
    if (any(v1 == a & v2 == b))
      bad_syntax("the entry for ", a, if (b != a) paste0(", ", b) else "",
                 " appears more than once.")
    v1 <- c(v1, a); v2 <- c(v2, b); ex <- c(ex, expr)
  }
  data.frame(var1 = v1, var2 = v2, expr = ex, stringsAsFactors = FALSE)
}


#' Parse Dynare's ramsey_constraints block
#'
#' Dynare syntax, one constraint per statement:
#' \preformatted{
#'   ramsey_constraints;
#'     i > 0;
#'     tau < tau_max;
#'   end;
#' }
#' A constraint bounds an endogenous variable of the Ramsey problem by an
#' expression in the parameters.  Dynare's preprocessor turns each into a
#' bound on that variable, complementary to the Ramsey FOC with respect to it
#' (`dynamic_complementarity_conditions`), which only its mixed-complementarity
#' perfect-foresight solver uses.  `>=` / `<=` are read as `>` / `<` (the
#' complementarity problem does not distinguish them).
#'
#' @param body Block body text.
#' @param var_names Declared endogenous variables.
#' @param param_names Declared parameters.
#' @return data.frame(var, op, bound): op is ">" (lower bound) or "<" (upper
#'   bound); bound is expression text.
#' @noRd
parse_ramsey_constraints_block <- function(body, var_names, param_names) {
  bad_syntax <- function(...)
    .dynhr_abort("parse_mod: ramsey_constraints: ", ...,
                 class = "dynhr_error_mod_syntax")
  stmts <- trimws(strsplit(body, ";", fixed = TRUE)[[1]])
  stmts <- stmts[nzchar(stmts)]
  v_out <- character(0); op_out <- character(0); b_out <- character(0)
  for (st in stmts) {
    st <- gsub("\\s+", " ", st)
    m <- regmatches(st, regexec(
      "^([A-Za-z_][A-Za-z0-9_]*)\\s*(>=|<=|>|<)\\s*(\\S.*)$", st,
      perl = TRUE))[[1]]
    if (length(m) != 4L)
      bad_syntax("cannot read the constraint `", st, "`; the form is ",
                 "`VARIABLE > EXPRESSION;` or `VARIABLE < EXPRESSION;`.")
    v <- m[2]; op <- substr(m[3], 1L, 1L); bx <- trimws(m[4])
    if (!(v %in% var_names))
      bad_syntax("`", v, "` is not a declared endogenous variable.")
    ## `pi` is the constant only when no variable/parameter is called pi
    unknown <- setdiff(.osr_expr_ids(bx),
                       c(param_names, .KNOWN_FUNCTIONS, "Inf", "inf",
                         setdiff("pi", var_names)))
    if (length(unknown) > 0L)
      bad_syntax("the bound `", bx, "` of ", v, " uses ",
                 paste(unknown, collapse = ", "),
                 ", which is not a declared parameter. A bound may only ",
                 "depend on parameters.")
    if (any(v_out == v & op_out == op))
      bad_syntax(v, " has more than one ", if (op == ">") "lower" else "upper",
                 " bound.")
    v_out <- c(v_out, v); op_out <- c(op_out, op); b_out <- c(b_out, bx)
  }
  data.frame(var = v_out, op = op_out, bound = b_out, stringsAsFactors = FALSE)
}


#' Parse Dynare's osr_params_bounds block
#'
#' Dynare syntax: \code{PARAMETER, LOWER, UPPER;} per statement.  Bounds are
#' parameter expressions (\code{Inf} / \code{-Inf} allowed), kept as text.
#' A parameter of \code{osr_params} without an entry is unbounded, as in
#' Dynare (\code{M_.osr.param_bounds} defaults to \code{[-Inf, Inf]}).
#'
#' @param body Block body text.
#' @param param_names Declared parameters.
#' @param osr_params Parameters named by \code{osr_params} (may be empty).
#' @return data.frame(name, lower, upper) holding expression text.
#' @noRd
parse_osr_params_bounds_block <- function(body, param_names,
                                          osr_params = character(0)) {
  bad_syntax <- function(...)
    .dynhr_abort("parse_mod: osr_params_bounds: ", ...,
                 class = "dynhr_error_mod_syntax")
  stmts <- trimws(strsplit(body, ";", fixed = TRUE)[[1]])
  stmts <- stmts[nzchar(stmts)]
  nm_out <- character(0); lo_out <- character(0); hi_out <- character(0)
  for (st in stmts) {
    st <- gsub("\\s+", " ", st)
    parts <- trimws(.split_top_level(st, sep = ",", open = "(", close = ")"))
    if (length(parts) != 3L || !all(nzchar(parts)) ||
        !grepl("^[A-Za-z_][A-Za-z0-9_]*$", parts[1]))
      bad_syntax("cannot read the entry `", st, "`; the form is ",
                 "`PARAMETER, LOWER_BOUND, UPPER_BOUND;`.")
    nm <- parts[1]
    if (!(nm %in% param_names))
      bad_syntax("`", nm, "` is not a declared parameter.")
    if (length(osr_params) > 0L && !(nm %in% osr_params))
      bad_syntax("`", nm, "` is not listed in osr_params (",
                 paste(osr_params, collapse = ", "), ").")
    if (nm %in% nm_out)
      bad_syntax("parameter ", nm, " has more than one entry.")
    for (bx in parts[2:3]) {
      unknown <- setdiff(.osr_expr_ids(bx),
                         c(param_names, .KNOWN_FUNCTIONS, "pi", "Inf", "inf"))
      if (length(unknown) > 0L)
        bad_syntax("the bound `", bx, "` of ", nm, " uses ",
                   paste(unknown, collapse = ", "),
                   ", which is not a declared parameter.")
    }
    nm_out <- c(nm_out, nm); lo_out <- c(lo_out, parts[2])
    hi_out <- c(hi_out, parts[3])
  }
  data.frame(name = nm_out, lower = lo_out, upper = hi_out,
             stringsAsFactors = FALSE)
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
  ## A2: every one of these environments is an allowlist sandbox (parented at
  ## emptyenv()), so a shocks-block RHS cannot reach system() and friends.
  eval_env <- if (is.environment(param_env)) param_env
              else .dynhr_sandbox_env(param_env, .dynhr_safe_matrix_fn_names)
  safe_eval <- function(text) {
    val <- .dynhr_sandbox_eval(text, eval_env, .dynhr_safe_matrix_fn_names,
                               context = "the shocks-block expression")
    if (is.null(val)) NA_real_ else val
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


#' Parse Dynare 7 `shock_paths` blocks into deterministic shock paths
#'
#' Dynare 7.2 reference manual, "The model file" > `shock_paths` block (with
#' the evaluation order of `make_ex_.m` and the preprocessor-generated
#' `+MODEL/shock_paths_N.m`, checked by running Dynare 7.1).  A stanza
#' \preformatted{  var X;
#'   periods 1, 2:4, 6:end;
#'   values  initval.X*1.05, self.X(-1)*1.05 + a/100, self.X(-1);}
#' gives one value EXPRESSION per period entry.  Dynare fills a matrix of exo
#' paths over the simulation periods 1..T plus, when `end` is used, a
#' terminal column T+1 (the terminal steady state, i.e. `endval`), starting
#' from the exogenous steady state (the `initval` value, else 0).  For each
#' period p = 1, 2, ... it evaluates every stanza in order, so an expression
#' sees earlier periods and, for earlier stanzas, the current one:
#' \itemize{
#'   \item `self.Y` / `self.Y(-k)`: the path of exogenous Y at p / p-k;
#'   \item `initval.Y` / `init.Y`: Y's `initval` value (0 if unset; an
#'     endogenous Y is its steady state after a `steady` command, which the
#'     parser cannot compute, so that case aborts);
#'   \item parameters by name.
#' }
#' The result is expressed in the same structure as the `shocks` block's
#' `periods`/`values` path (`model$det_shocks`: one row per set period
#' 1..T) plus the terminal values (`model$endval`), so every consumer of those
#' fields sees a `shock_paths` scenario exactly like the equivalent
#' `shocks` + `endval` one.
#'
#' Needs T (from `perfect_foresight_setup(periods=)` or `first/last_
#' simulation_period`) only for a range ending at `end` (`6:end`) or an `end`
#' value that looks back with `self.`; dates need `first_simulation_period`.
#' Not implemented, each a `dynhr_error_mod_syntax` abort: database scopes
#' (`db.foo`), `prev.` / `learnt_in(...).` scopes, `learnt_in` other than 1,
#' and a `shock_paths` block combined with an `endval` block or with
#' deterministic `shocks` periods (Dynare: "cannot be used in conjunction").
#' Controlled stanzas (`exogenize Y; periods ...; values ...; endogenize E;`,
#' see `.dynhr_parse_controlled_stanza()`) are returned in `controlled`, the
#' same table a `perfect_foresight_controlled_paths` block gives; their values
#' may use parameters and `initval.` but not `self.` (Dynare rejects it).
#'
#' @param blocks      List of blocks from \code{.extract_shock_paths_blocks()}.
#' @param exo_names   Declared exogenous (varexo + varexo_det) names.
#' @param endo_names  Declared endogenous names.
#' @param param_values Named numeric parameter values.
#' @param initval     Named numeric initval values (endogenous and exogenous).
#' @param has_steady  Logical: the .mod has a `steady` command.
#' @param n_periods   Simulation length T, or NA.
#' @param first_sim   First simulation period as a date literal, or NULL.
#' @return list(det = data.frame(name, period, value), terminal = named
#'   numeric of the exogenous values set at `end`, controlled =
#'   data.frame(exogenize, endogenize, period, value)).
#' @noRd
parse_shock_paths_blocks <- function(blocks, exo_names, endo_names,
                                     param_values, initval, has_steady,
                                     n_periods = NA_integer_,
                                     first_sim = NULL) {
  bad <- function(...) .dynhr_abort("parse_mod: shock_paths: ", ...,
                                    class = "dynhr_error_mod_syntax")
  env <- .dynhr_sandbox_env(param_values, .dynhr_safe_fn_names)
  base <- stats::setNames(rep(0, length(exo_names)), exo_names)
  iv_exo <- intersect(names(initval), exo_names)
  base[iv_exo] <- initval[iv_exo]

  ## -- parse the stanzas of every block ---------------------------------
  period_entry <- function(tok, vname) {
    tok <- trimws(tok)
    parts <- trimws(strsplit(tok, ":", fixed = TRUE)[[1]])
    if (length(parts) == 1L) parts <- c(parts, parts)
    if (length(parts) != 2L || !all(nzchar(parts)))
      bad("malformed period '", tok, "' for ", vname, ".")
    one <- function(s) {
      if (grepl("^[0-9]+$", s)) return(as.integer(s))
      if (tolower(s) == "end") return(Inf)
      if (is.null(.parse_dynare_date(s)))
        bad("period '", s, "' for ", vname, " is neither an integer, a date ",
            "nor `end`.")
      if (is.null(first_sim))
        .dynhr_abort(
          "parse_mod: shock_paths: the date '", s, "' (", vname, ") needs ",
          "perfect_foresight_setup(first_simulation_period = DATE), which ",
          "gives the date of simulation period 1.",
          class = c("dynhr_error_mod_date_unresolved", "dynhr_error_mod_syntax"))
      .date_to_period(s, first_sim, context = "shock_paths")
    }
    lo <- one(parts[1L]); hi <- one(parts[2L])
    if (is.infinite(lo) && !is.infinite(hi))
      bad("period range '", tok, "' for ", vname, " starts at `end`.")
    if (lo < 1L) bad("period '", tok, "' for ", vname, " is before period 1.")
    if (hi < lo) bad("empty period range '", tok, "' for ", vname, ".")
    c(lo, hi)
  }
  split_list <- function(s) {
    out <- trimws(.split_top_level(s, sep = ",", open = "(", close = ")"))
    out[nzchar(out)]
  }

  stanzas_all <- list()
  controlled <- .dynhr_empty_controlled(expr = TRUE)
  for (b in blocks) {
    opts <- parse_command_options(b$options_str)
    li <- opts$learnt_in
    if (!is.null(li) && !identical(as.character(li), "1"))
      bad("learnt_in = ", li, " (a simulation with expectation errors) is ",
          "not implemented.")
    stmts <- trimws(strsplit(b$body, ";", fixed = TRUE)[[1]])
    stmts <- stmts[nzchar(stmts)]
    stz <- list()
    i <- 1L
    while (i <= length(stmts)) {
      s <- stmts[i]
      mv <- regmatches(s, regexec("^var\\s+([A-Za-z_][A-Za-z0-9_]*)$", s,
                                  perl = TRUE))[[1]]
      mx <- regmatches(s, regexec("^exogenize\\s+", s, perl = TRUE))[[1]]
      if (length(mv) == 2L) {
        vname <- mv[2L]
        if (!(vname %in% exo_names))
          bad("`var ", vname, "` is not a declared exogenous variable.")
        if (i + 2L > length(stmts) ||
            !grepl("^periods\\s", stmts[i + 1L]) ||
            !grepl("^values\\s", stmts[i + 2L]))
          bad("`var ", vname, ";` must be followed by `periods ...;` and ",
              "`values ...;`.")
        ptoks <- split_list(sub("^periods\\s+", "", stmts[i + 1L]))
        vtoks <- split_list(sub("^values\\s+", "", stmts[i + 2L]))
        if (length(ptoks) != length(vtoks))
          bad(vname, " has ", length(ptoks), " period entries but ",
              length(vtoks), " values (they must match one to one).")
        rng <- lapply(ptoks, period_entry, vname = vname)
        stz[[length(stz) + 1L]] <- list(var = vname, lo = vapply(rng, `[`, 0, 1L),
                                        hi = vapply(rng, `[`, 0, 2L),
                                        values = vtoks)
        i <- i + 3L
      } else if (length(mx) == 1L) {
        ## Controlled stanza: values are evaluated below, once eval_value()
        ## exists; Dynare rejects `self.` here.
        cs <- .dynhr_parse_controlled_stanza(
          stmts, i, endo_names = endo_names, exo_names = exo_names,
          first_sim = first_sim, bad = bad)
        if (any(grepl("\\bself\\s*\\.", cs$rows$expr, perl = TRUE)))
          bad("`self.` is not accepted in an exogenize/endogenize stanza ",
              "(exogenize ", cs$rows$exogenize[1L], "), as in Dynare.")
        controlled <- rbind(controlled, cs$rows)
        i <- cs$i
      } else {
        bad("unrecognised statement `", s, "`.")
      }
    }
    stanzas_all[[length(stanzas_all) + 1L]] <-
      list(stanzas = stz, overwrite = isTRUE(opts$overwrite))
  }
  ## -- horizon ------------------------------------------------------------
  all_stz <- unlist(lapply(stanzas_all, `[[`, "stanzas"), recursive = FALSE)
  if (length(all_stz) == 0L && nrow(controlled) == 0L)
    return(list(det = data.frame(name = character(0), period = integer(0),
                                 value = numeric(0), stringsAsFactors = FALSE),
                terminal = numeric(0),
                controlled = .dynhr_empty_controlled()))
  his <- unlist(lapply(all_stz, `[[`, "hi"))
  los <- unlist(lapply(all_stz, `[[`, "lo"))
  has_end <- any(is.infinite(his))
  open_range <- any(is.infinite(his) & is.finite(los))
  T_known <- !is.na(n_periods)
  if (open_range && !T_known)
    bad("a period range ending at `end` needs the simulation length: add ",
        "perfect_foresight_setup(periods = T).")
  T_sim <- if (T_known) as.integer(n_periods) else as.integer(max(c(0, his[is.finite(his)])))
  if (T_known && any(his[is.finite(his)] > T_sim))
    bad("a period is beyond the simulation length T = ", T_sim, ".")
  n_col <- T_sim + as.integer(has_end)
  end_col <- if (has_end) n_col else NA_integer_

  ## -- evaluation ---------------------------------------------------------
  P   <- matrix(base, nrow = length(exo_names), ncol = max(n_col, 1L),
                dimnames = list(exo_names, NULL))
  SET <- matrix(FALSE, nrow = length(exo_names), ncol = max(n_col, 1L),
                dimnames = list(exo_names, NULL))
  scope_value <- function(nm, expr) {
    if (nm %in% exo_names) return(if (nm %in% names(initval)) initval[[nm]] else 0)
    if (nm %in% endo_names) {
      if (has_steady)
        bad("`initval.", nm, "` in `", expr, "` is the steady state of ",
            "endogenous ", nm, " (a `steady` command follows initval), which ",
            "the parser cannot compute; use a number or a parameter.")
      return(if (nm %in% names(initval)) initval[[nm]] else 0)
    }
    bad("`initval.", nm, "` in `", expr, "` is not a declared variable.")
  }
  eval_value <- function(expr, p, vname) {
    txt <- expr
    if (grepl("\\b(?:prev|learnt_in)\\b\\s*[.(]", txt, perl = TRUE))
      bad("`", expr, "` uses a `prev.` / `learnt_in(...)` scope (expectation ",
          "errors), which is not implemented.")
    ## initval.X / init.X
    txt <- .dynhr_gsub_fn(txt, "\\b(?:initval|init)\\.([A-Za-z_][A-Za-z0-9_]*)",
                   function(m) sprintf("(%.17g)", scope_value(m[2L], expr)))
    ## self.X(k) / self.X
    txt <- .dynhr_gsub_fn(txt, paste0("\\bself\\.([A-Za-z_][A-Za-z0-9_]*)",
                               "(?:\\s*\\(\\s*([+-]?)\\s*(\\d+)\\s*\\))?"),
                   function(m) {
      nm <- m[2L]
      if (!(nm %in% exo_names))
        bad("`self.", nm, "` in `", expr, "`: not a declared exogenous variable.")
      k <- if (nzchar(m[4L])) as.integer(m[4L]) * (if (m[3L] == "-") -1L else 1L)
           else 0L
      if (k > 0L)
        bad("`self.", nm, "(", k, ")` in `", expr, "` looks ahead; only the ",
            "current or a previous period can be referenced.")
      if (nm == vname && k == 0L)
        bad("`self.", vname, "` in its own values must reference a previous ",
            "period (`self.", vname, "(-1)`).")
      if (is.na(p))
        bad("`", expr, "` at `end` looks back with `self.` and needs the ",
            "simulation length: add perfect_foresight_setup(periods = T).")
      q <- p + k
      if (q < 1L)
        bad("`self.", nm, "(", k, ")` in `", expr, "` at period ", p,
            " refers to a period before 1.")
      sprintf("(%.17g)", P[nm, q])
    })
    if (grepl("\\b[A-Za-z_][A-Za-z0-9_]*\\s*\\.\\s*[A-Za-z_]", txt, perl = TRUE))
      bad("`", expr, "` refers to a database (`DB.VAR`), which is not ",
          "implemented; give the values explicitly.")
    v <- .dynhr_sandbox_eval(txt, env, .dynhr_safe_fn_names,
                             context = "the shock_paths value")
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
      bad("value `", expr, "` for ", vname, " does not evaluate to a finite ",
          "number (unknown name?).")
    as.numeric(v)
  }

  for (blk in stanzas_all) {
    if (blk$overwrite) {
      P[] <- base
      SET[] <- FALSE
    }
    cols <- seq_len(n_col)
    for (p in cols) {
      is_end <- has_end && p == end_col
      for (st in blk$stanzas) {
        k <- if (is_end) which(is.infinite(st$hi))
             else which(st$lo <= p & (st$hi >= p & (is.finite(st$lo))))
        if (length(k) == 0L) next
        k <- k[length(k)]
        p_eval <- if (is_end && !T_known) NA_integer_ else p
        P[st$var, p]   <- eval_value(st$values[[k]], p_eval, st$var)
        SET[st$var, p] <- TRUE
      }
    }
  }

  ## -- express as det_shocks rows + terminal values ------------------------
  vorder <- unique(vapply(all_stz, `[[`, "", "var"))
  det <- data.frame(name = character(0), period = integer(0),
                    value = numeric(0), stringsAsFactors = FALSE)
  terminal <- numeric(0)
  sim_cols <- seq_len(T_sim)
  for (v in vorder) {
    ps <- sim_cols[SET[v, sim_cols]]
    if (length(ps) > 0L)
      det <- rbind(det, data.frame(name = v, period = ps, value = P[v, ps],
                                   stringsAsFactors = FALSE))
    if (has_end && SET[v, end_col]) terminal[v] <- P[v, end_col]
  }
  rownames(det) <- NULL

  ## -- controlled (exogenize / endogenize) stanzas --------------------------
  ctl <- .dynhr_empty_controlled()
  if (nrow(controlled) > 0L)
    ctl <- data.frame(
      exogenize  = controlled$exogenize,
      endogenize = controlled$endogenize,
      period     = controlled$period,
      value      = vapply(seq_len(nrow(controlled)), function(r)
        eval_value(controlled$expr[r], NA_integer_, controlled$exogenize[r]),
        numeric(1)),
      stringsAsFactors = FALSE)
  list(det = det, terminal = terminal, controlled = ctl)
}


## Empty controlled-paths table (optionally with the unevaluated `expr`).
.dynhr_empty_controlled <- function(expr = FALSE) {
  out <- data.frame(exogenize = character(0), endogenize = character(0),
                    period = integer(0), value = numeric(0),
                    stringsAsFactors = FALSE)
  if (expr) {
    out$value <- NULL
    out$expr <- character(0)
  }
  out
}

#' Parse one controlled-paths stanza (Dynare 7)
#'
#' \preformatted{  exogenize Y;
#'   periods 1:2, 3;
#'   values 0.5, (a + 1);
#'   endogenize E;}
#' as it appears in a `perfect_foresight_controlled_paths` block or inside a
#' `shock_paths` block (Dynare 7.2 reference manual; statement order checked
#' against the Dynare 7.1 preprocessor, which stores one
#' `struct('exogenize_id', 'periods', 'value', 'endogenize_id')` per period
#' entry).  Periods are integers or ranges `a:b` (or dates, which need
#' `first_sim`); `end` is rejected, as by Dynare's grammar.  Each value
#' applies to every period of its entry.
#'
#' @param stmts Statements of the block (split at `;`, trimmed).
#' @param i     Index of the `exogenize` statement.
#' @param bad   Abort function of the caller (prefixes the message).
#' @return list(rows = data.frame(exogenize, endogenize, period, expr), one
#'   row per period; i = index of the statement after `endogenize`).
#' @noRd
.dynhr_parse_controlled_stanza <- function(stmts, i, endo_names, exo_names,
                                           first_sim, bad) {
  mx <- regmatches(stmts[i], regexec(
    "^exogenize\\s+([A-Za-z_][A-Za-z0-9_]*)$", stmts[i], perl = TRUE))[[1]]
  if (length(mx) != 2L)
    bad("malformed statement `", stmts[i], "` (expected `exogenize VAR;`, ",
        "one endogenous variable).")
  vname <- mx[2L]
  if (!(vname %in% endo_names))
    bad("`exogenize ", vname, "`: not a declared endogenous variable.")
  if (i + 3L > length(stmts) ||
      !grepl("^periods\\s", stmts[i + 1L]) ||
      !grepl("^values\\s", stmts[i + 2L]))
    bad("`exogenize ", vname, ";` must be followed by `periods ...;`, ",
        "`values ...;` and `endogenize SHOCK;`.")
  me <- regmatches(stmts[i + 3L], regexec(
    "^endogenize\\s+([A-Za-z_][A-Za-z0-9_]*)$", stmts[i + 3L], perl = TRUE))[[1]]
  if (length(me) != 2L)
    bad("`exogenize ", vname, "`: expected `endogenize SHOCK;` after its ",
        "values, got `", stmts[i + 3L], "`.")
  ename <- me[2L]
  if (!(ename %in% exo_names))
    bad("`endogenize ", ename, "`: not a declared exogenous variable.")
  split_list <- function(s) {
    out <- trimws(.split_top_level(s, sep = ",", open = "(", close = ")"))
    out[nzchar(out)]
  }
  ptoks <- split_list(sub("^periods\\s+", "", stmts[i + 1L]))
  vtoks <- split_list(sub("^values\\s+", "", stmts[i + 2L]))
  if (length(ptoks) == 0L || length(ptoks) != length(vtoks))
    bad("exogenize ", vname, " has ", length(ptoks), " period entries but ",
        length(vtoks), " values (they must match one to one).")
  one <- function(s) {
    if (grepl("^[0-9]+$", s)) return(as.integer(s))
    if (tolower(s) == "end")
      bad("`end` is not a valid period for exogenize ", vname, " (Dynare ",
          "accepts only integers, ranges and dates here).")
    if (is.null(.parse_dynare_date(s)))
      bad("period '", s, "' for exogenize ", vname, " is neither an integer ",
          "nor a date.")
    if (is.null(first_sim))
      .dynhr_abort(
        "parse_mod: controlled paths: the date '", s, "' (exogenize ", vname,
        ") needs perfect_foresight_setup(first_simulation_period = DATE), ",
        "which gives the date of simulation period 1.",
        class = c("dynhr_error_mod_date_unresolved", "dynhr_error_mod_syntax"))
    as.integer(.date_to_period(s, first_sim, context = "controlled paths"))
  }
  per <- integer(0); ex <- character(0)
  for (k in seq_along(ptoks)) {
    parts <- trimws(strsplit(ptoks[k], ":", fixed = TRUE)[[1]])
    if (length(parts) == 1L) parts <- c(parts, parts)
    if (length(parts) != 2L || !all(nzchar(parts)))
      bad("malformed period '", ptoks[k], "' for exogenize ", vname, ".")
    lo <- one(parts[1L]); hi <- one(parts[2L])
    if (lo < 1L)
      bad("period '", ptoks[k], "' for exogenize ", vname, " is before period 1.")
    if (hi < lo)
      bad("empty period range '", ptoks[k], "' for exogenize ", vname, ".")
    per <- c(per, seq.int(lo, hi))
    ex  <- c(ex, rep(vtoks[k], hi - lo + 1L))
  }
  list(rows = data.frame(exogenize = rep(vname, length(per)),
                         endogenize = rep(ename, length(per)),
                         period = per, expr = ex, stringsAsFactors = FALSE),
       i = i + 4L)
}

#' Parse Dynare 7 `perfect_foresight_controlled_paths` blocks
#'
#' Dynare 7.2 reference manual, `perfect_foresight_controlled_paths`: a run of
#' `exogenize` / `periods` / `values` / `endogenize` stanzas (see
#' `.dynhr_parse_controlled_stanza()`); values are numbers or expressions in
#' the parameters.  Only `learnt_in = 1` (no expectation errors) is
#' implemented.
#'
#' @param blocks List(options_str, body) from `extract_paired_block`-style
#'   extraction (`.extract_controlled_paths_blocks()`).
#' @return data.frame(exogenize, endogenize, period, value).
#' @noRd
parse_controlled_paths_blocks <- function(blocks, endo_names, exo_names,
                                          param_values, first_sim = NULL) {
  bad <- function(...) .dynhr_abort("parse_mod: perfect_foresight_controlled_paths: ",
                                    ..., class = "dynhr_error_mod_syntax")
  env <- .dynhr_sandbox_env(param_values, .dynhr_safe_fn_names)
  out <- .dynhr_empty_controlled()
  for (b in blocks) {
    opts <- parse_command_options(b$options_str)
    li <- opts$learnt_in
    if (!is.null(li) && !identical(as.character(li), "1"))
      bad("learnt_in = ", li, " (a simulation with expectation errors) is ",
          "not implemented.")
    stmts <- trimws(strsplit(b$body, ";", fixed = TRUE)[[1]])
    stmts <- stmts[nzchar(stmts)]
    i <- 1L
    while (i <= length(stmts)) {
      if (!grepl("^exogenize\\s", stmts[i]))
        bad("unrecognised statement `", stmts[i], "` (expected `exogenize`).")
      cs <- .dynhr_parse_controlled_stanza(stmts, i, endo_names, exo_names,
                                           first_sim, bad)
      rows <- cs$rows
      vals <- vapply(rows$expr, function(e) {
        v <- .dynhr_sandbox_eval(e, env, .dynhr_safe_fn_names,
                                 context = "the perfect_foresight_controlled_paths value")
        if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
          bad("value `", e, "` for exogenize ", rows$exogenize[1L], " does not ",
              "evaluate to a finite number (unknown name?).")
        as.numeric(v)
      }, numeric(1), USE.NAMES = FALSE)
      out <- rbind(out, data.frame(exogenize = rows$exogenize,
                                   endogenize = rows$endogenize,
                                   period = rows$period, value = vals,
                                   stringsAsFactors = FALSE))
      i <- cs$i
    }
  }
  out
}

#' Every `perfect_foresight_controlled_paths[(options)]; ... end;` block
#'
#' Closes at the first statement that is exactly `end` (as for shock_paths),
#' so a stray `periods 2:end;` reaches the stanza parser's own error.
#' @noRd
.extract_controlled_paths_blocks <- function(txt) {
  pat <- paste0("(?si)\\bperfect_foresight_controlled_paths\\b\\s*",
                .dynhr_opts_re, "\\s*;((?:[^;]*;)*?)\\s*\\bend\\s*;")
  ms  <- gregexpr(pat, txt, perl = TRUE)[[1]]
  if (ms[1L] == -1L) return(list())
  lens <- attr(ms, "match.length")
  lapply(seq_along(ms), function(k) {
    chunk <- substr(txt, ms[k], ms[k] + lens[k] - 1L)
    parts <- regmatches(chunk, regexec(pat, chunk, perl = TRUE))[[1]]
    list(options_str = trimws(parts[2L]), body = trimws(parts[3L]))
  })
}

#' Canonicalise and check the combined controlled-paths table
#'
#' Dynare 7.1 `controlled_paths_by_period.m`: a variable exogenized twice, or
#' a shock endogenized twice, in one period is an error, as is a shock that is
#' both endogenized and given a deterministic value in the same period.
#' Rows are sorted by (period, exogenize, endogenize) so the table does not
#' depend on stanza order (which Dynare ignores) and write_mod round-trips.
#'
#' @param cp data.frame(exogenize, endogenize, period, value).
#' @param n_periods Simulation length or NA.
#' @param det det_shocks data.frame(name, period, value).
#' @return The sorted data.frame (row names reset).
#' @noRd
.dynhr_check_controlled <- function(cp, n_periods, det) {
  bad <- function(...) .dynhr_abort("parse_mod: controlled paths: ", ...,
                                    class = "dynhr_error_mod_syntax")
  cp <- cp[order(cp$period, cp$exogenize, cp$endogenize), , drop = FALSE]
  rownames(cp) <- NULL
  if (!is.na(n_periods) && any(cp$period > n_periods))
    bad("period ", max(cp$period), " is beyond the simulation length T = ",
        n_periods, ".")
  d <- duplicated(cp[, c("period", "exogenize")])
  if (any(d))
    bad("variable ", cp$exogenize[d][1L], " is exogenized two times in ",
        "period ", cp$period[d][1L], ".")
  d <- duplicated(cp[, c("period", "endogenize")])
  if (any(d))
    bad("shock ", cp$endogenize[d][1L], " is endogenized two times in ",
        "period ", cp$period[d][1L], ".")
  if (is.data.frame(det) && nrow(det) > 0L) {
    hit <- paste(cp$endogenize, cp$period) %in% paste(det$name, det$period)
    if (any(hit))
      bad("shock ", cp$endogenize[hit][1L], " is both given a deterministic ",
          "value and endogenized in period ", cp$period[hit][1L], ".")
  }
  cp
}

## gsub() with a function of the match's capture groups (m[1] = whole match).
.dynhr_gsub_fn <- function(txt, pattern, fn) {
  mm <- gregexpr(pattern, txt, perl = TRUE)[[1]]
  if (mm[1L] == -1L) return(txt)
  starts <- as.integer(mm)
  lens   <- attr(mm, "match.length")
  out <- character(0)
  cursor <- 1L
  for (j in seq_along(starts)) {
    whole <- substr(txt, starts[j], starts[j] + lens[j] - 1L)
    m <- regmatches(whole, regexec(pattern, whole, perl = TRUE))[[1]]
    out <- c(out, substr(txt, cursor, starts[j] - 1L), fn(m))
    cursor <- starts[j] + lens[j]
  }
  paste(c(out, substr(txt, cursor, nchar(txt))), collapse = "")
}

#' Every `shock_paths[(options)]; ... end;` block, in file order
#'
#' The block closes at the first statement that is exactly `end` (a period
#' list may say `6:end` or `end`); see `.dynhr_paired_block_re()`.
#' @noRd
.extract_shock_paths_blocks <- function(txt) {
  pat <- .dynhr_paired_block_re("shock_paths")
  ms  <- gregexpr(pat, txt, perl = TRUE)[[1]]
  if (ms[1L] == -1L) return(list())
  lens <- attr(ms, "match.length")
  lapply(seq_along(ms), function(k) {
    chunk <- substr(txt, ms[k], ms[k] + lens[k] - 1L)
    parts <- regmatches(chunk, regexec(pat, chunk, perl = TRUE))[[1]]
    list(options_str = trimws(parts[2L]), body = trimws(parts[3L]))
  })
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
      .dynhr_warn(sprintf(
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
#' Options are separated by TOP-LEVEL commas only: a comma nested in
#' \code{()}, \code{[]}, \code{\{\}} or inside a quoted string belongs to the
#' option value (review 2026-09-25: `discretionary_policy(instruments=(i,tau))`
#' used to split into `instruments=(i` and a bogus flag `tau)`).
#'
#' Values:
#' \itemize{
#'   \item a number -> numeric; \code{true}/\code{false} -> logical;
#'   \item a list \code{(a, b)} / \code{(a b)} / \code{[1 2 3]} (commas or
#'     blanks) -> a numeric vector when every element is a number or an
#'     integer range \code{m:n} (expanded), otherwise a character vector of
#'     the names.  The vector carries attribute \code{"dynare_list"} (the
#'     opening bracket) so write_mod() re-emits the list form even for one
#'     element (\code{instruments=(i)});
#'   \item a list containing a quoted string (\code{optim=('MaxIter',200)},
#'     a MATLAB cell of name/value pairs) and anything else -> the value text,
#'     verbatim.
#' }
#'
#' @param options_str The text inside the parentheses.
#' @return Named list. Bare flags have value TRUE; key=value pairs are parsed.
#' @noRd
parse_command_options <- function(options_str) {
  if (is.null(options_str) || nchar(trimws(options_str)) == 0)
    return(list())

  opts <- list()
  parts <- trimws(.split_options_top_level(options_str, ","))

  for (p in parts) {
    if (nchar(p) == 0) next
    eq <- .options_top_level_eq(p)
    if (eq > 0L) {
      key <- trimws(substr(p, 1L, eq - 1L))
      val_str <- trimws(substr(p, eq + 1L, nchar(p)))
      opts[[key]] <- .parse_option_value(val_str)
    } else {
      opts[[p]] <- TRUE
    }
  }
  opts
}

#' Split on a separator outside (), [], {} and quotes
#' @noRd
.split_options_top_level <- function(s, sep = ",") {
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  out <- character(0)
  cur <- character(0)
  depth <- 0L
  quote <- ""
  for (ch in chars) {
    if (nzchar(quote)) {
      if (ch == quote) quote <- ""
    } else if (ch == "'" || ch == "\"") {
      quote <- ch
    } else if (ch %in% c("(", "[", "{")) {
      depth <- depth + 1L
    } else if (ch %in% c(")", "]", "}")) {
      depth <- max(0L, depth - 1L)
    } else if (ch == sep && depth == 0L) {
      out <- c(out, paste(cur, collapse = ""))
      cur <- character(0)
      next
    }
    cur <- c(cur, ch)
  }
  c(out, paste(cur, collapse = ""))
}

#' Position of the first top-level `=` of one option (0 when none)
#' @noRd
.options_top_level_eq <- function(p) {
  chars <- strsplit(p, "", fixed = TRUE)[[1]]
  depth <- 0L
  quote <- ""
  for (k in seq_along(chars)) {
    ch <- chars[k]
    if (nzchar(quote)) {
      if (ch == quote) quote <- ""
    } else if (ch == "'" || ch == "\"") {
      quote <- ch
    } else if (ch %in% c("(", "[", "{")) {
      depth <- depth + 1L
    } else if (ch %in% c(")", "]", "}")) {
      depth <- max(0L, depth - 1L)
    } else if (ch == "=" && depth == 0L) {
      return(k)
    }
  }
  0L
}

#' Convert one option value text to its R value (see parse_command_options)
#' @noRd
.parse_option_value <- function(val_str) {
  num <- suppressWarnings(as.numeric(val_str))
  if (!is.na(num)) return(num)
  if (tolower(val_str) == "true")  return(TRUE)
  if (tolower(val_str) == "false") return(FALSE)

  n <- nchar(val_str)
  open <- substr(val_str, 1L, 1L)
  close <- substr(val_str, n, n)
  is_list <- n >= 2L && ((open == "(" && close == ")") ||
                         (open == "[" && close == "]"))
  ## A list whose brackets close before the end, e.g. `(a)*(b)`, is an
  ## expression, not a list.
  if (is_list) {
    ch <- strsplit(substr(val_str, 2L, n - 1L), "", fixed = TRUE)[[1]]
    step <- (ch %in% c("(", "[", "{")) - (ch %in% c(")", "]", "}"))
    if (any(cumsum(step) < 0L)) is_list <- FALSE
  }
  if (!is_list || grepl("['\"]", val_str)) return(val_str)

  inner <- trimws(substr(val_str, 2L, n - 1L))
  toks <- unlist(lapply(.split_options_top_level(inner, ","), function(piece) {
    piece <- trimws(piece)
    ## blanks separate elements too (`(a b)`, `[6 32]`), except inside
    ## nested brackets
    if (!nzchar(piece)) return(character(0))
    if (grepl("[([{]", piece)) return(piece)
    strsplit(piece, "\\s+", perl = TRUE)[[1]]
  }))
  toks <- toks[nzchar(toks)]
  if (length(toks) == 0L) {
    out <- character(0)
  } else {
    rng <- regmatches(toks, regexec("^([+-]?[0-9]+):([+-]?[0-9]+)$", toks))
    is_rng <- lengths(rng) == 3L
    nums <- suppressWarnings(as.numeric(toks))
    if (all(!is.na(nums) | is_rng)) {
      out <- as.numeric(unlist(lapply(seq_along(toks), function(k) {
        if (is_rng[k]) seq(as.numeric(rng[[k]][2]), as.numeric(rng[[k]][3]))
        else nums[k]
      })))
    } else {
      out <- toks
    }
  }
  attr(out, "dynare_list") <- open
  out
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
  env <- .dynhr_sandbox_env(param_values, .dynhr_safe_matrix_fn_names)
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
      ## Matrix algebra CAN fail at eval time on legitimate input (singular
      ## solve(), non-conformable %*%), and the documented contract of this
      ## function is that it never throws, so the handler here is a genuine
      ## need rather than a swallowed parse error. The SECURITY check happened
      ## before eval, inside .dynhr_sandbox_eval().
      val <- tryCatch(.dynhr_sandbox_eval(r_rhs, env,
                                          .dynhr_safe_matrix_fn_names,
                                          context = "the verbatim statement"),
                      error = function(e) if (inherits(
                        e, "dynhr_error_unsafe_mod_expression")) stop(e) else NULL)
      if (!is.null(val) && is.numeric(val))
        assign(lhs, val, envir = env)
    }
  }
  env
}
