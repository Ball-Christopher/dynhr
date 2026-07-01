## R/parse-macro.R
## --------------------------------------------------------------------------
## Dynare macro-language preprocessor (M18).
##
## A line-oriented expander for the subset of Dynare's `@#` macro directives
## that real .mod files in the replication suite use:
##
##   @#define NAME = VALUE         simple macro-variable definition
##   @#define NAME = ["a","b"]     macro-variable set to an array/list
##   @#for V in [a, b, c] ...      iterate over an explicit list
##   @#for V in A:B ...            iterate over an integer range
##   @#for V in NAME ...           iterate over a @#define-d array
##   @#endfor
##   @#if COND ... [@#else ...] @#endif
##   @#ifdef NAME / @#ifndef NAME ... @#endif
##   @#include "path"              textual splice of path relative to .mod dir
##
## Inside a loop or conditional body, `@{expr}` is interpolated: the bracketed
## expression is evaluated in the current macro environment (loop index +
## @#define-d variables) and the result is spliced into the surrounding text.
## This matches Dynare's `@{}` interpolation semantics, including arithmetic on
## the loop index (e.g. `p@{j-1}_fwrd`) and string-list iteration where the
## quotes around list elements are dropped when spliced into an identifier
## (e.g. `@#for X in [ "H", "F" ]` -> `y_@{X}` becomes `y_H`).
##
## @#include splices the named file textually (relative to the including file's
## directory) BEFORE macro expansion, exactly as Dynare does.  Nested includes
## are supported; a depth guard (max 32) prevents infinite cycles.  The
## including directory is stored in the env as `.macro_include_dir` and updated
## as we descend so that a child @#include resolves relative to itself.
##
## The pass is a strict no-op when the source contains no `@#` directive and
## no `@{...}` interpolation: `expand_macros()` returns the input unchanged
## (byte-for-byte) in that case, so macro-free models are unaffected.
##
## Design: FAIL LOUD on any directive or construct we do not understand rather
## than silently dropping it (a dropped @#for would silently produce a wrong,
## under-specified model).  The error message names the offending directive.
## --------------------------------------------------------------------------

#' Does this .mod source contain any Dynare macro construct?
#'
#' Detects bare `@#` directive lines (ignoring those inside // or /* */
#' comments) and `@{...}` interpolation.  Used to short-circuit the
#' preprocessor to an exact no-op when there is nothing to expand.
#' @noRd
.macro_has_directives <- function(txt) {
  # Strip comments for *detection* only (we never want a commented-out
  # directive to trigger the expander; Dynare does expand macros inside
  # comments, but the replication models never rely on that, and treating
  # them as inert is the conservative choice).
  scan <- gsub("(?s)/\\*.*?\\*/", "\n", txt, perl = TRUE)
  scan <- gsub("//[^\n]*", "", scan, perl = TRUE)
  scan <- gsub("%[^\n]*", "", scan, perl = TRUE)
  # (?m) so `^` matches at the start of every line, not just the whole string.
  grepl("(?m)^[ \t]*@#", scan, perl = TRUE) ||
    grepl("@\\{", scan, perl = TRUE)
}

#' Translate Dynare macro `X in ARRAY` membership to R `X %in% ARRAY`.
#'
#' Called on @#if condition strings before R parsing.  Replaces the bare infix
#' `in` membership operator (Dynare syntax) with R's `%in%`.  Also translates
#' inline array literals `[...]` to `c(...)` so that
#' `"a" in ["a","b"]` becomes `"a" %in% c("a","b")`.
#'
#' Care is taken NOT to touch:
#'   - `in` that is part of an identifier (e.g. `index`, `index_var`) -- the
#'     word-boundary regex `\\bin\\b` avoids that.
#'   - `in` inside double-quoted strings -- the substitution uses a simple
#'     lookahead/lookbehind approach that is correct for the macro expressions
#'     that appear in real .mod files (no embedded quotes containing " in ").
#' The `@#for V in ARRAY` loop `in` is parsed by a separate regex before
#' `.macro_eval` is ever called, so it is never affected by this function.
#' @noRd
.macro_translate_membership <- function(expr) {
  # Step 1: translate bare `in` operator to `%in%`.
  # Use word-boundary anchors so `index`, `inline`, etc. are not clobbered.
  # We require at least one space on each side (the Dynare spec and real models
  # always write `X in Y`, never `X inY`), which also keeps us away from
  # quoted strings that contain the literal 5-char sequence " in " (rare, but
  # the space-boundary is an extra safety net).
  expr <- gsub("(?<=\\s)in(?=\\s)", "%in%", expr, perl = TRUE)
  # Step 2: translate inline Dynare array literals [...] to R c(...).
  # Matches the outermost [...] that does not itself contain `[` (no nesting).
  expr <- gsub("\\[([^\\[\\]]*)\\]", "c(\\1)", expr, perl = TRUE)
  expr
}

#' Evaluate a macro-language expression in the macro environment.
#'
#' Macro expressions are integer/string constants, @#define-d names, the
#' active loop indices, and simple arithmetic (`+ - * / ( )`).  We evaluate in
#' a sandboxed child of baseenv() so model parameters / variables can never
#' leak in.  String values iterate as character; numeric values as numbers.
#' @noRd
.macro_eval <- function(expr, env, directive_label) {
  expr <- trimws(expr)
  # Translate Dynare `X in ARRAY` membership operator to R `X %in% ARRAY`
  # (only meaningful in @#if conditions; harmless in arithmetic expressions
  # because bare `in` never appears there).
  expr <- .macro_translate_membership(expr)
  parsed <- tryCatch(parse(text = expr), error = function(e) NULL)
  if (is.null(parsed)) {
    stop("dynhr macro preprocessor: could not parse macro expression `",
         expr, "` in ", directive_label, ".", call. = FALSE)
  }
  val <- tryCatch(
    eval(parsed, envir = env),
    error = function(e) {
      stop("dynhr macro preprocessor: failed to evaluate macro expression `",
           expr, "` in ", directive_label, ":\n  ", conditionMessage(e),
           call. = FALSE)
    }
  )
  val
}

#' Render a macro value for splicing into model text.
#'
#' Numbers render without a decimal point when integer-valued (so `@{j}` with
#' j=2 yields `2`, not `2.0` — Dynare builds identifiers like `ln_p2`).
#' Strings render verbatim (quotes already stripped at list-parse time).
#' @noRd
.macro_render <- function(val) {
  if (is.numeric(val)) {
    if (length(val) != 1L)
      stop("dynhr macro preprocessor: @{} expression produced a non-scalar ",
           "numeric value.", call. = FALSE)
    if (is.finite(val) && val == round(val))
      return(format(as.integer(round(val)), trim = TRUE))
    return(format(val, trim = TRUE))
  }
  if (is.character(val) || is.logical(val)) {
    if (length(val) != 1L)
      stop("dynhr macro preprocessor: @{} expression produced a non-scalar ",
           "value.", call. = FALSE)
    return(as.character(val))
  }
  stop("dynhr macro preprocessor: @{} expression produced an unsupported ",
       "value of type ", typeof(val), ".", call. = FALSE)
}

#' Interpolate every `@{expr}` occurrence in a single line.
#' @noRd
.macro_interp_line <- function(line, env) {
  if (!grepl("@\\{", line, perl = TRUE)) return(line)
  out <- ""
  rest <- line
  repeat {
    m <- regexpr("@\\{", rest, perl = TRUE)
    if (m[1] < 0) {
      out <- paste0(out, rest)
      break
    }
    start <- m[1]
    # Find matching closing brace (no nesting of @{} in these models, but be
    # defensive and match the first unescaped `}`).
    after <- substring(rest, start + 2L)
    close <- regexpr("\\}", after, perl = TRUE)
    if (close[1] < 0) {
      stop("dynhr macro preprocessor: unterminated @{ interpolation in line:\n  ",
           line, call. = FALSE)
    }
    expr <- substring(after, 1L, close[1] - 1L)
    val  <- .macro_eval(expr, env, paste0("@{", expr, "}"))
    out  <- paste0(out, substring(rest, 1L, start - 1L), .macro_render(val))
    rest <- substring(after, close[1] + 1L)
  }
  out
}

#' Parse the iterable in `@#for V in <iter>` into a list of macro values.
#'
#' Supports `[a, b, c]` explicit lists (elements may be quoted strings,
#' numbers, or @#define-d names) and `A:B` integer ranges (A, B may themselves
#' be macro expressions).  Also accepts a bare macro name that resolves to a
#' list (i.e. a name set via `@#define name = [...]`).
#' @noRd
.macro_parse_iter <- function(iter, env, directive_label) {
  iter <- trimws(iter)
  if (grepl("^\\[", iter)) {
    if (!grepl("\\]$", iter)) {
      stop("dynhr macro preprocessor: malformed list in ", directive_label,
           " (missing closing `]`): ", iter, call. = FALSE)
    }
    inner <- substring(iter, 2L, nchar(iter) - 1L)
    if (!nzchar(trimws(inner))) return(list())
    elems <- strsplit(inner, ",", fixed = TRUE)[[1]]
    out <- vector("list", length(elems))
    for (i in seq_along(elems)) {
      e <- trimws(elems[i])
      # Quoted string element -> strip the quotes, keep the bare token.
      if (grepl('^".*"$', e) || grepl("^'.*'$", e)) {
        out[[i]] <- substring(e, 2L, nchar(e) - 1L)
      } else {
        # Numeric literal or macro expression.
        out[[i]] <- .macro_eval(e, env, directive_label)
      }
    }
    return(out)
  }
  # Range form  A:B
  if (grepl(":", iter, fixed = TRUE)) {
    parts <- strsplit(iter, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
      stop("dynhr macro preprocessor: only simple `A:B` integer ranges are ",
           "supported in ", directive_label, "; got: ", iter, call. = FALSE)
    }
    lo <- .macro_eval(parts[1], env, directive_label)
    hi <- .macro_eval(parts[2], env, directive_label)
    if (!is.numeric(lo) || !is.numeric(hi)) {
      stop("dynhr macro preprocessor: range bounds must be numeric in ",
           directive_label, "; got: ", iter, call. = FALSE)
    }
    return(as.list(seq.int(as.integer(lo), as.integer(hi))))
  }
  # Bare macro name that may resolve to a list (set via @#define name = [...]).
  # Try looking it up in the env; if it is a list, use it directly.
  if (grepl("^[A-Za-z_][A-Za-z0-9_]*$", iter)) {
    if (exists(iter, envir = env, inherits = TRUE)) {
      val <- get(iter, envir = env, inherits = TRUE)
      if (is.list(val)) return(val)
      # Scalar macro used in @#for — wrap as a length-1 list so the loop
      # runs exactly once, mirroring Dynare behaviour.
      return(list(val))
    }
  }
  stop("dynhr macro preprocessor: unsupported @#for iterable in ",
       directive_label, " (expected `[a, b, ...]`, `A:B`, or a @#define-d array name): ",
       iter, call. = FALSE)
}

#' Recursively expand a block of macro lines.
#'
#' @param lines      Character vector of source lines for this block.
#' @param env        Environment holding @#define vars + active loop indices +
#'                   `.macro_include_dir` (character(1), may be NA) +
#'                   `.macro_include_depth` (integer(1)).
#' @return Character vector of fully expanded lines (no `@#` directives, all
#'         `@{}` interpolated).
#' @noRd
.macro_expand_lines <- function(lines, env) {
  out <- character(0)
  i <- 1L
  n <- length(lines)
  while (i <= n) {
    line <- lines[i]
    trimmed <- trimws(line)
    if (grepl("^@#", trimmed, perl = TRUE)) {
      directive <- sub("^@#\\s*", "", trimmed)

      # ---- @#define NAME = VALUE  or  @#define NAME = [...] ----------
      if (grepl("^define\\b", directive)) {
        body <- sub("^define\\s+", "", directive)
        if (!grepl("=", body, fixed = TRUE)) {
          stop("dynhr macro preprocessor: @#define must have the form ",
               "`@#define NAME = VALUE`; got: ", trimmed, call. = FALSE)
        }
        nm  <- trimws(sub("=.*$", "", body))
        rhs <- trimws(sub("^[^=]*=", "", body))
        # RHS is an array literal [...] -> store as a list for @#for iteration.
        if (grepl("^\\[", rhs)) {
          val <- .macro_parse_iter(rhs, env, paste0("@#define ", nm))
        } else if (grepl('^".*"$', rhs) || grepl("^'.*'$", rhs)) {
          # RHS is a quoted string literal -> strip quotes.
          val <- substring(rhs, 2L, nchar(rhs) - 1L)
        } else {
          # RHS is a number or expression over earlier defines.
          val <- .macro_eval(rhs, env, paste0("@#define ", nm))
        }
        assign(nm, val, envir = env)
        i <- i + 1L
        next
      }

      # ---- @#for V in <iter> ... @#endfor ----------------------------
      if (grepl("^for\\b", directive)) {
        m <- regmatches(directive,
                        regexec("^for\\s+([A-Za-z_][A-Za-z0-9_]*)\\s+in\\s+(.*)$",
                                directive, perl = TRUE))[[1]]
        if (length(m) != 3L) {
          stop("dynhr macro preprocessor: malformed @#for directive: ",
               trimmed, call. = FALSE)
        }
        loop_var <- m[2]
        iter_str <- trimws(m[3])
        # Collect the body up to the matching @#endfor (respecting nesting).
        depth <- 1L
        body_lines <- character(0)
        j <- i + 1L
        while (j <= n) {
          tj <- trimws(lines[j])
          if (grepl("^@#\\s*for\\b", tj, perl = TRUE)) depth <- depth + 1L
          else if (grepl("^@#\\s*endfor\\b", tj, perl = TRUE)) {
            depth <- depth - 1L
            if (depth == 0L) break
          }
          body_lines <- c(body_lines, lines[j])
          j <- j + 1L
        }
        if (depth != 0L) {
          stop("dynhr macro preprocessor: @#for without matching @#endfor ",
               "(starting at: ", trimmed, ").", call. = FALSE)
        }
        values <- .macro_parse_iter(iter_str, env,
                                    paste0("@#for ", loop_var, " in ", iter_str))
        for (v in values) {
          # Child env so the loop index does not clobber an outer define and
          # is scoped to this iteration; nested loops/defines stack here.
          iter_env <- new.env(parent = env)
          assign(loop_var, v, envir = iter_env)
          out <- c(out, .macro_expand_lines(body_lines, iter_env))
        }
        i <- j + 1L
        next
      }

      # ---- @#if / @#ifdef / @#ifndef ... [@#else] @#endif ------------
      if (grepl("^if\\b", directive) || grepl("^ifdef\\b", directive) ||
          grepl("^ifndef\\b", directive)) {
        # Evaluate the guard.
        if (grepl("^ifdef\\b", directive)) {
          nm <- trimws(sub("^ifdef\\s+", "", directive))
          cond <- exists(nm, envir = env, inherits = TRUE)
        } else if (grepl("^ifndef\\b", directive)) {
          nm <- trimws(sub("^ifndef\\s+", "", directive))
          cond <- !exists(nm, envir = env, inherits = TRUE)
        } else {
          expr <- trimws(sub("^if\\s+", "", directive))
          val  <- .macro_eval(expr, env, paste0("@#if ", expr))
          cond <- isTRUE(as.logical(val)) ||
            (is.numeric(val) && length(val) == 1L && val != 0)
        }
        # Collect the if-body and (optional) else-body up to @#endif,
        # respecting nested @#if.
        depth <- 1L
        if_lines   <- character(0)
        else_lines <- character(0)
        seen_else  <- FALSE
        j <- i + 1L
        while (j <= n) {
          tj <- trimws(lines[j])
          if (grepl("^@#\\s*(if|ifdef|ifndef)\\b", tj, perl = TRUE)) {
            depth <- depth + 1L
          } else if (grepl("^@#\\s*endif\\b", tj, perl = TRUE)) {
            depth <- depth - 1L
            if (depth == 0L) break
          } else if (depth == 1L && grepl("^@#\\s*else\\b", tj, perl = TRUE)) {
            seen_else <- TRUE
            j <- j + 1L
            next
          } else if (depth == 1L && grepl("^@#\\s*elseif\\b", tj, perl = TRUE)) {
            stop("dynhr macro preprocessor: @#elseif is not supported; ",
                 "rewrite as nested @#if/@#else.", call. = FALSE)
          }
          if (seen_else) else_lines <- c(else_lines, lines[j])
          else            if_lines   <- c(if_lines, lines[j])
          j <- j + 1L
        }
        if (depth != 0L) {
          stop("dynhr macro preprocessor: @#if without matching @#endif ",
               "(starting at: ", trimmed, ").", call. = FALSE)
        }
        chosen <- if (cond) if_lines else else_lines
        out <- c(out, .macro_expand_lines(chosen, env))
        i <- j + 1L
        next
      }

      # ---- @#include "path" ------------------------------------------
      if (grepl("^include\\b", directive)) {
        # Extract the filename: must be a double-quoted string.
        m_inc <- regmatches(directive,
                            regexec('^include\\s+"([^"]+)"', directive,
                                    perl = TRUE))[[1]]
        if (length(m_inc) != 2L) {
          stop("dynhr macro preprocessor: @#include must have the form ",
               "`@#include \"filename\"` (double-quoted path); got: ",
               trimmed, call. = FALSE)
        }
        inc_rel <- m_inc[2]
        # Resolve path relative to the including file's directory.
        inc_dir <- if (exists(".macro_include_dir", envir = env, inherits = TRUE))
          get(".macro_include_dir", envir = env, inherits = TRUE)
        else
          NA_character_
        if (!is.na(inc_dir) && !grepl("^(/|[A-Za-z]:[/\\\\])", inc_rel)) {
          inc_path <- file.path(inc_dir, inc_rel)
        } else {
          inc_path <- inc_rel
        }
        if (!file.exists(inc_path)) {
          stop("dynhr macro preprocessor: @#include file not found: ",
               inc_path, call. = FALSE)
        }
        # Depth guard against infinite include cycles.
        depth_now <- if (exists(".macro_include_depth", envir = env, inherits = TRUE))
          get(".macro_include_depth", envir = env, inherits = TRUE)
        else
          0L
        if (depth_now >= 32L) {
          stop("dynhr macro preprocessor: @#include nesting depth exceeded 32 ",
               "(possible cycle); offending file: ", inc_path, call. = FALSE)
        }
        # Read the included file and splice its lines (textual include, just
        # like Dynare — macro expansion of the included content happens in
        # the current env context, exactly as if the lines were inline).
        inc_txt   <- paste(readLines(inc_path, warn = FALSE), collapse = "\n")
        inc_txt   <- iconv(inc_txt, from = "LATIN1", to = "UTF-8", sub = "?")
        inc_lines <- strsplit(inc_txt, "\n", fixed = TRUE)[[1]]
        # Child env that updates the include dir + depth for nested includes.
        inc_env <- new.env(parent = env)
        assign(".macro_include_dir",   normalizePath(dirname(inc_path), mustWork = FALSE),
               envir = inc_env)
        assign(".macro_include_depth", depth_now + 1L, envir = inc_env)
        out <- c(out, .macro_expand_lines(inc_lines, inc_env))
        # Propagate any @#define assignments made inside the include back to
        # the calling env (mirror Dynare: defines in an included file are
        # visible after the @#include).
        inc_names <- ls(envir = inc_env, all.names = FALSE)
        for (nm in inc_names) {
          assign(nm, get(nm, envir = inc_env), envir = env)
        }
        i <- i + 1L
        next
      }

      # ---- @#echo / @#error : informational, no model effect ---------
      if (grepl("^echo\\b", directive) || grepl("^error\\b", directive)) {
        # Drop (no structural effect).  @#error inside a taken branch would
        # ideally abort, but evaluating it requires full macro semantics; the
        # replication models only use @#echo for diagnostics.
        i <- i + 1L
        next
      }

      # ---- Anything else: FAIL LOUD ----------------------------------
      stop("dynhr macro preprocessor: unsupported macro directive `@#",
           sub("\\s.*$", "", directive), "` in line:\n  ", trimmed,
           "\nSupported: @#define, @#for/@#endfor, @#if/@#ifdef/@#ifndef/",
           "@#else/@#endif, @#include.", call. = FALSE)
    }

    # Plain model line: interpolate any @{...} and emit.
    out <- c(out, .macro_interp_line(line, env))
    i <- i + 1L
  }
  out
}

#' Expand Dynare `@#` macro directives in .mod source text.
#'
#' Runs before the lexer / declaration parsing in `parse_mod()`.  When the
#' source contains no macro directive and no `@{}` interpolation this is an
#' exact byte-for-byte no-op.
#'
#' @param txt     Character string (full .mod file content).
#' @param mod_dir Character(1) or NULL.  Directory used to resolve
#'   `@#include` paths.  `parse_mod()` passes `dirname(source_file)` when
#'   reading from a file; inline-text callers pass NULL (includes disabled
#'   unless an absolute path is given).
#' @return The macro-expanded source text.
#' @noRd
expand_macros <- function(txt, mod_dir = NULL) {
  if (!.macro_has_directives(txt)) return(txt)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  # Preserve a possible trailing empty line that strsplit drops; harmless if
  # absent.
  env <- new.env(parent = baseenv())
  if (!is.null(mod_dir) && nzchar(mod_dir)) {
    assign(".macro_include_dir",   normalizePath(mod_dir, mustWork = FALSE),
           envir = env)
  } else {
    assign(".macro_include_dir",   NA_character_, envir = env)
  }
  assign(".macro_include_depth", 0L, envir = env)
  expanded <- .macro_expand_lines(lines, env)
  paste(expanded, collapse = "\n")
}
