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

.strip_blocks <- function(txt) {
  # Remove all  keyword; ... end;  block constructs that may

  # contain 'var' or other declaration keywords internally.
  # Order matters: longest keyword first to avoid partial matches.
  block_kws <- c(
    "steady_state_model",
    "estimated_params_init",
    "estimated_params",
    "observation_trends",
    "optim_weights",
    "osr_params_bounds",
    "heteroskedastic_shocks",
    "shocks",
    "mshocks",
    "filter_tunes",
    "model",
    "initval",
    "endval"
  )
  for (kw in block_kws) {
    pat <- paste0("(?si)\\b", kw, "\\b\\s*(?:\\([^)]*\\))?\\s*;.*?\\bend\\s*;")
    txt <- gsub(pat, " ", txt, perl = TRUE)
  }
  txt
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
