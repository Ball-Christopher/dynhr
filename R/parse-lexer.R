## R/parse-lexer.R
## --------------------------------------------------------------------------
## Comment and macro stripping for .mod file source text. Pure regex
## utilities consumed by parse_mod() before any block extraction.
##
## Phase-1 split from parser-monolith.R (no logic changes).
## --------------------------------------------------------------------------

#' Remove C-style and line comments, and macro directives from .mod text
#'
#' Handles:
#'   - Single-line comments:  // ...
#'   - Block comments:        /* ... */
#'   - Dynare macro lines:    @#...
#'   - MATLAB-style comments: % ... (anywhere on a line)
#'
#' @param txt Character string (full .mod file content).
#' @return Cleaned character string.
#' @noRd
strip_comments_and_macros <- function(txt) {
  # Block comments (non-greedy, DOTALL via (?s))
  txt <- gsub("(?s)/\\*.*?\\*/", " ", txt, perl = TRUE)
  # Single-line comments: // ...
  txt <- gsub("//[^\n]*", " ", txt, perl = TRUE)
  # MATLAB-style comments: % ... (anywhere on line)
  txt <- gsub("%[^\n]*", " ", txt, perl = TRUE)
  # Dynare macro directives: @#...
  txt <- gsub("@#[^\n]*", " ", txt, perl = TRUE)
  txt
}

## The ONE list of `keyword[(options)]; ... end;` blocks known to the parser.
##
## Both the declaration scan (`.strip_blocks()`, so a `var y;` INSIDE a block
## is not read as a declaration) and the calibration scan (`remove_blocks()`,
## so `x = ...;` inside a block is not read as a parameter assignment) strip
## exactly these. They used to keep two different lists, and a block on only
## one of them leaked: `conditional_forecast_paths; var y; ...` declared `y` a
## second time (review 2026-09-25 B12).
##
## Only genuine blocks (closed by `end;`) belong here. A COMMAND such as
## `identification(...)`, `var_model(...)` or `pac_model(...)` ends at its own
## `;`, and treating it as a block would delete everything up to the next
## `end;` (it is handled in remove_blocks()'s command list instead).
## Longest keyword first; every pattern is also `\b`-delimited on both sides.
.dynhr_paired_block_kws <- c(
  "steady_state_model",
  "estimated_params_init", "estimated_params_bounds",
  "estimated_params_remove", "estimated_params",
  "observation_trends", "deterministic_trends",
  "optim_weights", "osr_params_bounds", "ramsey_constraints",
  "heteroskedastic_shocks", "stochastic_volatility", "mshocks", "shocks",
  "filter_tunes", "filter_initial_state",
  "model_replace", "model",
  "initval", "endval", "histval", "homotopy_setup",
  "moment_calibration", "irf_calibration", "matched_moments",
  "occbin_constraints", "shock_groups", "conditional_forecast_paths",
  "perfect_foresight_controlled_paths",
  "svar_identification", "generate_irfs", "epilogue", "verbatim"
)

## Blocks whose statements may themselves contain the word `end` before the
## closing `end;` -- Dynare 7's `shock_paths` (`periods 6:end;`,
## `periods end;`).  These close at the first statement that is exactly
## `end`, never at an `end;` inside a statement.
.dynhr_stmt_end_block_kws <- c("shock_paths")

## Optional `(options)` of a block or command head, BRACKET-AWARE: balanced
## parentheses via PCRE recursion into capture group 1, which holds the text
## inside the outer parentheses (read by parse_command_options(), whose
## top-level splitter then separates the options).  The old `\(([^)]*)\)`
## stopped at the first `)`, so `model(differentiate_forward_vars=(x y));`
## was not recognised as the model block at all.  It MUST be the first
## capture group of any pattern it is pasted into (`(?1)` recurses into it).
## Shared by the block finders here and in parse-blocks.R and by
## extract_command().
.dynhr_opts_re <- "(?:\\(((?:[^()]++|\\((?1)\\))*)\\))?"

## Regex for one `keyword[(options)]; ... end;` block.  Groups: 1 = options,
## 2 = body.  For .dynhr_stmt_end_block_kws the body is a run of whole
## statements and the block ends at a statement that is exactly `end`.
.dynhr_paired_block_re <- function(kw) {
  head <- paste0("(?si)\\b", kw, "\\b\\s*", .dynhr_opts_re, "\\s*;")
  if (kw %in% .dynhr_stmt_end_block_kws)
    paste0(head, "((?:[^;]*;)*?)\\s*\\bend\\s*;")
  else
    paste0(head, "(.*?)\\bend\\s*;")
}

#' Replace every known `keyword; ... end;` block with a space
#' @noRd
.dynhr_strip_paired_blocks <- function(txt) {
  for (kw in c(.dynhr_stmt_end_block_kws, .dynhr_paired_block_kws)) {
    txt <- gsub(.dynhr_paired_block_re(kw), " ", txt, perl = TRUE)
  }
  txt
}

.strip_blocks <- function(txt) {
  # Remove all  keyword; ... end;  block constructs that may contain 'var' or
  # other declaration keywords internally (shared list, see above).
  .dynhr_strip_paired_blocks(txt)
}

.strip_mod_comments <- function(txt) {
  # Remove block comments /* ... */ (non-greedy, handles multiline)
  txt <- gsub("/\\*.*?\\*/", " ", txt, perl = TRUE)
  # Remove line comments // ... (but not inside strings)
  txt <- gsub("//[^\n]*", " ", txt, perl = TRUE)
  # Remove line comments % ... (Dynare also supports %)
  txt <- gsub("%[^\n]*", " ", txt, perl = TRUE)
  txt
}
