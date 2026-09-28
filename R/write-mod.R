## R/write-mod.R
## --------------------------------------------------------------------------
## write_mod(): render a parsed `dynhr_mod` object back to Dynare .mod source.
##
## Design contract (E1-C):
##   * write_mod() is faithful to the PARSED model, not to the original file
##     text.  parse_mod(write_mod(parse_mod(f))) must equal parse_mod(f) for
##     every field the parser populates.
##   * Anything the parser stored that this exporter cannot render must
##     stop() naming the field.  A silently lossy exporter is worse than none,
##     so there is an explicit whitelist of known fields (.WM_KNOWN_FIELDS)
##     and an unknown field is a hard error.
##
## Also here: smoother2histval(), which turns a smoother result into the
## `histval` list structure that parse_mod() produces and write_mod() renders.
## --------------------------------------------------------------------------


## Every field new_dynhr_mod() can populate.  Fields are either RENDERED
## (written out as .mod syntax) or DERIVED (recomputed by parse_mod from the
## rendered content, so nothing needs to be written for them).  A field that
## is neither is a bug -- write_mod() stops on it.
.WM_RENDERED_FIELDS <- c(
  "var_names", "varexo_names", "varexo_det_names", "param_names",
  "predetermined_vars", "param_values", "equations", "local_variables",
  "model_options", "initval", "endval", "histval", "steady_state_model",
  "shocks", "det_shocks", "shock_groups_blocks", "filter_tunes",
  "heteroskedastic_shocks",
  "stochastic_volatility", "estimated_params", "estimated_params_init",
  "commands", "occbin_constraints", "planner_objective", "varobs",
  ## present only when the .mod uses the construct
  "observation_trends", "log_vars", "deterministic_trends",
  ## trend_var / var(deflator=) declarations: rendered as the DETRENDED model
  ## they produced, and not re-declared (that would detrend a second time);
  ## like predetermined_vars, the field is absent after the round trip
  "nonstationary",
  ## osr_params / optim_weights / osr_params_bounds (present only when used)
  "osr",
  ## ramsey_constraints block (present only when used)
  "ramsey_constraints",
  ## marks det_shocks / endval as coming from a Dynare 7 shock_paths block
  "shock_paths",
  ## Dynare 7 perfect_foresight_controlled_paths / shock_paths exogenize
  "controlled_paths"
)

.WM_DERIVED_FIELDS <- c(
  ## purely bookkeeping / not part of the model text
  "source_file",
  ## aliases of rendered fields
  "exo_names", "varobs_names", "obs_vars",
  ## recomputed by parse_mod from the rendered content
  "lead_lag_incidence", "variable_classification",
  "n_static", "n_predetermined", "n_forward", "n_mixed",
  "equation_param_names", "mcp_constraints", "shocks_blocks",
  ## shock_groups is the FIRST entry of shock_groups_blocks, which is rendered
  "shock_groups",
  "ramsey_instruments",
  ## dynhr-specific `@dynhr:` comment annotations; carried over verbatim from
  ## the source file when there are any (see .wm_metadata_block)
  "metadata"
)

.WM_KNOWN_FIELDS <- c(.WM_RENDERED_FIELDS, .WM_DERIVED_FIELDS)


# ---------------------------------------------------------------------------
# Numeric / expression rendering
# ---------------------------------------------------------------------------

#' Render a scalar double as .mod source text that reads back bit-identically
#'
#' Uses the shortest decimal representation that `as.numeric()` maps back to
#' the exact same double (1..17 significant digits), so
#' `parse_mod(write_mod(m))` reproduces every literal exactly.
#'
#' @param x Length-1 finite numeric.
#' @param what Label used in error messages.
#' @return Character string.
#' @noRd
.wm_num <- function(x, what = "value") {
  if (length(x) != 1L || !is.numeric(x))
    stop("write_mod(): ", what, " must be a length-1 numeric.", call. = FALSE)
  if (is.na(x))
    stop("write_mod(): ", what, " is NA and cannot be written to .mod source.",
         call. = FALSE)
  if (!is.finite(x))
    stop("write_mod(): ", what, " is non-finite (", x,
         ") and cannot be written to .mod source.", call. = FALSE)
  if (x == 0) return("0")
  sci <- abs(x) < 1e-4 || abs(x) >= 1e15
  for (d in 1:17) {
    s <- format(x, digits = d, scientific = sci)
    if (isTRUE(suppressWarnings(as.numeric(s)) == x)) return(s)
  }
  ## R's string -> double conversion (R_strtod, which is what parse_mod()'s
  ## tokenizer ultimately calls through as.numeric()) is NOT correctly rounded
  ## beyond ~15 significant digits: `as.numeric(sprintf("%.30e", x))` can miss
  ## x by an ulp.  So `format()`'s shortest-round-trip ladder is not the whole
  ## search space -- sweep explicit scientific widths, and nudge the decimal
  ## by a few ulps so R's own (slightly off) accumulation lands ON x.  This
  ## recovers all but a handful of computed doubles.
  ulp <- 2^(floor(log2(abs(x))) - 52L)
  for (k in 15:25) {
    for (j in c(0L, 1L, -1L, 2L, -2L, 3L, -3L)) {
      s <- sprintf(paste0("%.", k, "e"), x + j * ulp)
      if (isTRUE(suppressWarnings(as.numeric(s)) == x)) return(s)
    }
  }
  ## Genuinely not expressible: write the closest text and say so.  Never
  ## silently drop the last bits.
  s <- format(x, digits = 17, scientific = sci)
  .dynhr_warn(sprintf(paste0(
    "write_mod(): %s = %.17g cannot be written as .mod text that R reads ",
    "back bit-identically (R's as.numeric() is not correctly rounded beyond ",
    "~15 significant digits). Wrote '%s', which reads back as %.17g ",
    "(relative error %.2g)."),
    what, x, s, as.numeric(s), abs(as.numeric(s) - x) / abs(x)),
    call. = FALSE)
  s
}


## Operator precedence used to decide where parentheses are required.
## relational < +,- < *,/ < unary < ^ < atoms.
.wm_prec <- function(node) {
  switch(node$type,
         "binop" = switch(node$op,
                          "<" = 1L, ">" = 1L, "<=" = 1L, ">=" = 1L,
                          "+" = 2L, "-" = 2L,
                          "*" = 3L, "/" = 3L,
                          "^" = 5L,
                          1L),
         "unaryop" = 4L,
         100L)
}


#' Render an AST node as Dynare .mod expression text
#'
#' Parenthesisation is chosen so that re-parsing the string reproduces the
#' SAME tree (not merely a mathematically equivalent one): the right operand
#' of a same-precedence `+ - * /` is parenthesised, and the left operand of
#' `^` is parenthesised, because those are the cases where re-parsing would
#' otherwise re-associate and change floating-point evaluation order.
#'
#' @param node AST node.
#' @return Character string.
#' @noRd
.wm_expr <- function(node) {
  if (is.null(node)) return("0")
  if (!is.list(node) || is.null(node$type))
    stop("write_mod(): malformed AST node (no $type).", call. = FALSE)
  switch(node$type,
         "number" = {
           s <- .wm_num(node$value, "numeric literal")
           if (node$value < 0) paste0("(", s, ")") else s
         },
         "variable" = {
           if (node$lead_lag == 0L) node$name
           else if (node$lead_lag > 0L)
             sprintf("%s(+%d)", node$name, as.integer(node$lead_lag))
           else sprintf("%s(%d)", node$name, as.integer(node$lead_lag))
         },
         "parameter"      = node$name,
         ## A model-local variable is REFERENCED by its bare name; the "#" only
         ## introduces its definition.  (ast_to_string() prints "#name", which
         ## is a comment in .mod source -- never use it for emission.)
         "local_variable" = node$name,
         "unaryop" = {
           o <- .wm_expr(node$operand)
           if (.wm_prec(node$operand) <= 4L) o <- paste0("(", o, ")")
           paste0(node$op, o)
         },
         "binop" = {
           p  <- .wm_prec(node)
           l  <- .wm_expr(node$left)
           r  <- .wm_expr(node$right)
           lp <- .wm_prec(node$left)
           rp <- .wm_prec(node$right)
           if (lp < p || (lp == p &&
                          node$op %in% c("^", "<", ">", "<=", ">=")))
             l <- paste0("(", l, ")")
           if (rp < p || (rp == p && node$op != "^"))
             r <- paste0("(", r, ")")
           paste0(l, " ", node$op, " ", r)
         },
         "funcall" = {
           ## EXPECTATION(k)(expr) has two argument groups in .mod syntax but
           ## is stored as a 2-argument funcall.
           if (identical(node$name, "EXPECTATION") && length(node$args) == 2L &&
               identical(node$args[[1]]$type, "number")) {
             return(paste0("EXPECTATION(",
                           .wm_num(node$args[[1]]$value, "EXPECTATION horizon"),
                           ")(", .wm_expr(node$args[[2]]), ")"))
           }
           paste0(node$name, "(",
                  paste(vapply(node$args, .wm_expr, character(1)),
                        collapse = ", "), ")")
         },
         stop("write_mod(): cannot render AST node of type '", node$type,
              "'.", call. = FALSE)
  )
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

.wm_nonempty_df <- function(x) is.data.frame(x) && nrow(x) > 0L

.wm_na <- function(x) length(x) != 1L || is.na(x)

## A shocks-block expression (`stderr sig_a;`, `corr a,b = rho*s1*s2;`) is
## stored as RAW TEXT and re-evaluated later against the current parameter
## vector.  It therefore only survives a round trip when every identifier it
## names is a DECLARED parameter: .mod-level "model-local constants" (a bare
## `phi = 0.1;` that is not in the `parameters` list) are consumed by the
## parser and never stored, so an expression referencing one cannot be
## written back.  Fail loud rather than silently substituting the numeric
## snapshot -- that would freeze a value that may depend on estimated
## parameters.
.wm_check_expr_ids <- function(expr_text, allowed, what) {
  if (.wm_na(expr_text)) return(invisible(NULL))
  ids <- regmatches(expr_text,
                    gregexpr("[A-Za-z_][A-Za-z0-9_]*", expr_text,
                             perl = TRUE))[[1]]
  bad <- setdiff(unique(ids), c(allowed, .KNOWN_FUNCTIONS))
  if (length(bad) > 0L)
    stop("write_mod(): the ", what, " expression '", expr_text,
         "' references identifier(s) that are not declared parameters: ",
         paste(bad, collapse = ", "),
         ". These are .mod-level local constants that parse_mod() consumes ",
         "without storing, so the expression cannot be written back. ",
         "Declare them in the `parameters` block.", call. = FALSE)
  invisible(NULL)
}


## ---- exogenous auxiliary-variable un-expansion ---------------------------
## parse_mod() rewrites every non-zero timing on an EXOGENOUS variable into an
## AUX_EXO_LEAD_*/AUX_EXO_LAG_* endogenous chain -- and that rewrite is not
## idempotent: feeding it its own output would rewrite `x(+1)` inside the
## chain's own defining equation and emit a duplicate (singular) definition.
## write_mod() therefore folds the exogenous chains back to plain `x(+k)` /
## `x(-k)` references and drops their defining equations; parse_mod()
## regenerates both, in the same order, on the next parse.
##
## Endogenous AUX_LEAD_*/AUX_LAG_* chains need no such treatment: they only
## reference timings within [-1, +1], which the expander leaves alone.

.wm_exo_aux_map <- function(model) {
  exo <- c(model$varexo_names, model$varexo_det_names)
  map <- list()
  for (nm in model$var_names) {
    m <- regmatches(nm, regexec("^AUX_EXO_(LEAD|LAG)_(.+)_([0-9]+)$", nm))[[1]]
    if (length(m) != 4L) next
    base <- m[3]
    if (!(base %in% exo)) next
    k <- as.integer(m[4])
    ## AUX_EXO_LEAD_x_k stands for x at timing +k;
    ## AUX_EXO_LAG_x_k  stands for x at timing -(k - 1) (its entry node is at t).
    map[[nm]] <- list(base = base,
                      offset = if (m[2] == "LEAD") k else -(k - 1L))
  }
  map
}

.wm_sub_exo_aux <- function(node, map) {
  if (is.null(node)) return(node)
  switch(node$type,
         "variable" = {
           e <- map[[node$name]]
           if (is.null(e)) return(node)
           ast_variable(e$base, node$lead_lag + e$offset)
         },
         "binop" = {
           node$left  <- .wm_sub_exo_aux(node$left,  map)
           node$right <- .wm_sub_exo_aux(node$right, map)
           node
         },
         "unaryop" = {
           node$operand <- .wm_sub_exo_aux(node$operand, map)
           node
         },
         "funcall" = {
           node$args <- lapply(node$args, .wm_sub_exo_aux, map = map)
           node
         },
         node)
}

.wm_unexpand_exo_aux <- function(model) {
  map <- .wm_exo_aux_map(model)
  if (length(map) == 0L) return(model)
  aux <- names(map)

  keep <- vapply(model$equations, function(eq) {
    ## Drop the chain's own defining equation: `AUX_EXO_* = <its own meaning>`,
    ## which folds to the identity `x(+k) = x(+k)`.
    if (!identical(eq$lhs$type, "variable") ||
        !(eq$lhs$name %in% aux)) return(TRUE)
    l <- .wm_sub_exo_aux(eq$lhs, map)
    r <- .wm_sub_exo_aux(eq$rhs, map)
    !identical(l, r)
  }, logical(1))

  model$equations <- lapply(model$equations[keep], function(eq) {
    eq$lhs <- .wm_sub_exo_aux(eq$lhs, map)
    eq$rhs <- .wm_sub_exo_aux(eq$rhs, map)
    eq
  })
  model$local_variables <- lapply(model$local_variables, .wm_sub_exo_aux,
                                  map = map)
  model$var_names <- setdiff(model$var_names, aux)
  model
}

## ---- var(log) un-expansion ------------------------------------------------
## parse_mod() rewrites every `x` of a `var(log) x;` variable into
## `exp(LOG_x)`, declares LOG_x after the declared variables and appends
## `x = exp(LOG_x);`. Like the exogenous chains above, that rewrite is not
## idempotent (a second pass would find no `x` to rewrite but would append the
## definition again), so write_mod() folds it back: `exp(LOG_x(t))` -> `x(t)`,
## the definitional equation is dropped, LOG_x leaves the plain `var` line and
## `x` is declared under `var(log)`. parse_mod() rebuilds all of it, in the
## same order.

.wm_sub_log <- function(node, inv) {
  if (is.null(node)) return(node)
  switch(node$type,
         "funcall" = {
           if (identical(node$name, "exp") && length(node$args) == 1L &&
               identical(node$args[[1L]]$type, "variable") &&
               node$args[[1L]]$name %in% names(inv))
             return(ast_variable(inv[[node$args[[1L]]$name]],
                                 node$args[[1L]]$lead_lag))
           node$args <- lapply(node$args, .wm_sub_log, inv = inv)
           node
         },
         "binop" = {
           node$left  <- .wm_sub_log(node$left,  inv)
           node$right <- .wm_sub_log(node$right, inv)
           node
         },
         "unaryop" = {
           node$operand <- .wm_sub_log(node$operand, inv)
           node
         },
         node)
}

## Variable names an AST still references (to catch a LOG_x the fold missed).
.wm_ast_var_names <- function(node) {
  if (is.null(node)) return(character(0))
  switch(node$type,
         "variable" = node$name,
         "binop"    = c(.wm_ast_var_names(node$left),
                        .wm_ast_var_names(node$right)),
         "unaryop"  = .wm_ast_var_names(node$operand),
         "funcall"  = unlist(lapply(node$args, .wm_ast_var_names),
                             use.names = FALSE),
         character(0))
}

.wm_unexpand_log_vars <- function(model) {
  lm <- model$log_vars
  if (is.null(lm) || length(lm) == 0L) return(model)
  inv <- stats::setNames(names(lm), unname(lm))      # LOG_x -> x
  n_v <- length(model$var_names)
  ## parse_mod appends the LOG_ variables and their definitional equations
  ## right after the declared ones; anything it appends LATER (EXPECTATION or
  ## |lead/lag| > 1 AUX variables) is written out explicitly and would then
  ## land BEFORE the regenerated LOG_ block on the next parse.
  if (!identical(utils::tail(model$var_names, length(lm)), unname(lm)))
    stop("write_mod(): the var(log) variables (", paste(lm, collapse = ", "),
         ") are followed by auxiliary variables in the parsed model, so the ",
         "model cannot be written back with var(log) in the same variable ",
         "order.", call. = FALSE)

  keep <- vapply(model$equations, function(eq) {
    if (!identical(eq$lhs$type, "variable") ||
        !(eq$lhs$name %in% names(lm)) || !.wm_na(eq$tag_raw %||% NA))
      return(TRUE)
    !identical(.wm_sub_log(eq$lhs, inv), .wm_sub_log(eq$rhs, inv))
  }, logical(1))
  dropped <- vapply(model$equations[!keep], function(eq) eq$lhs$name,
                    character(1))
  if (!setequal(dropped, names(lm)) || length(dropped) != length(lm))
    stop("write_mod(): could not find the definitional equation ",
         "`x = exp(LOG_x)` of every var(log) variable.", call. = FALSE)
  model$equations <- lapply(model$equations[keep], function(eq) {
    eq$lhs <- .wm_sub_log(eq$lhs, inv)
    eq$rhs <- .wm_sub_log(eq$rhs, inv)
    eq
  })
  model$local_variables <- lapply(model$local_variables, .wm_sub_log,
                                  inv = inv)
  left <- intersect(unique(c(
    unlist(lapply(model$equations, function(eq)
      c(.wm_ast_var_names(eq$lhs), .wm_ast_var_names(eq$rhs))),
      use.names = FALSE),
    unlist(lapply(model$local_variables, .wm_ast_var_names),
           use.names = FALSE))), names(inv))
  if (length(left) > 0L)
    stop("write_mod(): the var(log) variable(s) ", paste(left, collapse = ", "),
         " are referenced outside exp(...) after parsing (e.g. through an ",
         "auxiliary lead/lag chain), which cannot be written back as ",
         "var(log).", call. = FALSE)
  model$var_names <- model$var_names[seq_len(n_v - length(lm))]
  model
}

## `var` declaration lines, with the var(log) variables under `var(log)`.
## Consecutive runs keep the declaration order, which is what parse_mod()
## rebuilds var_names from.
.wm_var_lines <- function(var_names, log_names, annotations) {
  if (length(var_names) == 0L) return(character(0))
  is_log <- var_names %in% log_names
  if (!any(is_log)) return(.wm_decl_line("var", var_names, annotations))
  runs <- rle(is_log)
  ends <- cumsum(runs$lengths)
  out <- character(0)
  for (k in seq_along(ends)) {
    idx <- (ends[k] - runs$lengths[k] + 1L):ends[k]
    out <- c(out, .wm_decl_line(if (runs$values[k]) "var(log)" else "var",
                                var_names[idx], annotations))
  }
  out
}

## observation_trends block: `name (expression);` per observable.  Entries
## that parse_mod() took from the deterministic_trends block `dt` (a variable
## may not be in both) are written back there instead.
.wm_obs_trends_block <- function(ot, dt = NULL) {
  out <- character(0)
  tr <- if (is.null(ot)) character(0) else ot$trends
  tr <- tr[!(names(tr) %in% names(dt))]
  if (length(tr) > 0L)
    out <- c(out, "observation_trends;",
             paste0(names(tr), " (", unname(tr), ");"), "end;", "")
  if (length(dt) > 0L)
    out <- c(out, "deterministic_trends;",
             paste0(names(dt), " (", unname(dt), ");"), "end;", "")
  out
}

## Render one `key = value` / bare-flag option list back to the text that
## parse_command_options() consumes.  A list value (attribute "dynare_list",
## e.g. `instruments=(i,tau)` or `bandpass_filter=[6 32]`) is re-emitted in
## its own bracket form, so a one-element list stays `instruments=(i)`.
.wm_options_str <- function(opts) {
  if (is.null(opts) || length(opts) == 0L) return("")
  parts <- character(0)
  for (k in names(opts)) {
    v <- opts[[k]]
    br <- attr(v, "dynare_list")
    if (!is.null(br)) {
      closer <- if (identical(br, "[")) "]" else ")"
      els <- if (is.numeric(v)) vapply(as.numeric(v), .wm_num, "", what = k)
             else as.character(v)
      sepr <- if (identical(br, "[")) " " else ", "
      parts <- c(parts, paste0(k, "=", br, paste(els, collapse = sepr), closer))
      next
    }
    if (isTRUE(v)) { parts <- c(parts, k); next }
    if (isFALSE(v)) { parts <- c(parts, paste0(k, "=false")); next }
    if (is.numeric(v) && length(v) == 1L) {
      parts <- c(parts, paste0(k, "=", .wm_num(v, k))); next
    }
    if (is.character(v) && length(v) == 1L) {
      parts <- c(parts, paste0(k, "=", v)); next
    }
    stop("write_mod(): cannot render command option '", k,
         "' of class ", class(v)[1], ".", call. = FALSE)
  }
  paste(parts, collapse = ", ")
}


## Verbatim `// @dynhr:...` metadata blocks, harvested from the source file.
## The parsed `metadata` field is a lossy digest of comment text, so the only
## faithful way to carry it across a round trip is to copy the comment block.
.wm_metadata_block <- function(model) {
  md <- model$metadata
  if (is.null(md)) return(character(0))
  has <- any(vapply(md, function(x) length(x) > 0L, logical(1)))
  if (!has) return(character(0))
  sf <- model$source_file
  if (.wm_na(sf) || !file.exists(sf)) {
    stop("write_mod(): the model carries `@dynhr:` metadata (",
         paste(names(md)[vapply(md, function(x) length(x) > 0L, logical(1))],
               collapse = ", "),
         ") but its `source_file` is not available, so the metadata comment ",
         "blocks cannot be reproduced. Re-parse from the .mod file, or drop ",
         "`model$metadata` if the annotations are not needed.", call. = FALSE)
  }
  lines <- readLines(sf, warn = FALSE)
  lines <- iconv(lines, from = "LATIN1", to = "UTF-8", sub = "?")
  out <- character(0)
  keep <- FALSE
  for (ln in lines) {
    st <- trimws(ln)
    tag <- regmatches(st, regexec("^//\\s*@dynhr:(\\w+)", st))[[1]]
    if (length(tag) == 2L) {
      if (tag[2] == "end") { if (keep) out <- c(out, st); keep <- FALSE; next }
      keep <- TRUE
      out  <- c(out, st)
      next
    }
    if (keep && grepl("^//", st)) out <- c(out, st)
  }
  ## `%(key='value')` annotations sit on the var/varexo declaration lines and
  ## are re-emitted there (see .wm_decl_line), not here.
  if (length(out) == 0L && length(md$annotations %||% list()) == 0L) {
    stop("write_mod(): `model$metadata` is populated but no `// @dynhr:` ",
         "blocks were found in '", sf, "'.", call. = FALSE)
  }
  c(out, "")
}


## `var`/`varexo`/... declaration line, re-attaching any %(...) annotations.
.wm_decl_line <- function(keyword, names_vec, annotations = list()) {
  if (length(names_vec) == 0L) return(character(0))
  toks <- vapply(names_vec, function(nm) {
    a <- annotations[[nm]]
    if (is.null(a) || length(a) == 0L) return(nm)
    kv <- paste(vapply(names(a),
                       function(k) sprintf("%s='%s'", k, a[[k]]),
                       character(1)), collapse = ", ")
    sprintf("%s %%(%s)", nm, kv)
  }, character(1), USE.NAMES = FALSE)
  paste0(keyword, " ", paste(toks, collapse = " "), ";")
}


# ---------------------------------------------------------------------------
# Block renderers
# ---------------------------------------------------------------------------

.wm_model_block <- function(model) {
  opts <- .wm_options_str(model$model_options)
  head <- if (nzchar(opts)) paste0("model(", opts, ");") else "model;"
  out  <- head

  for (nm in names(model$local_variables)) {
    out <- c(out, paste0("#", nm, " = ",
                         .wm_expr(model$local_variables[[nm]]), ";"))
  }
  if (length(model$local_variables) > 0L) out <- c(out, "")

  for (eq in model$equations) {
    tag <- eq$tag_raw
    if (!.wm_na(tag) && grepl("]", tag, fixed = TRUE))
      stop("write_mod(): equation tag '", tag, "' contains ']', which cannot ",
           "be rendered as a .mod equation tag.", call. = FALSE)
    prefix <- if (.wm_na(tag)) "" else paste0("[", tag, "] ")
    out <- c(out, paste0(prefix, .wm_expr(eq$lhs), " = ",
                         .wm_expr(eq$rhs), .wm_complementarity(eq), ";"))
  }
  c(out, "end;", "")
}


## A Dynare 7 complementarity condition (eq$complementarity, see
## .mod_extract_complementarity()) as ` _|_ L < x < U` (ASCII separator, so
## the written file does not depend on the reader's encoding).  The bounds
## are the text the parser stored.
.wm_complementarity <- function(eq) {
  cc <- eq$complementarity
  if (is.null(cc)) return("")
  lo <- cc$lower
  up <- cc$upper
  body <- if (!.wm_na(lo) && !.wm_na(up))
    paste0(lo, " < ", cc$variable, " < ", up)
  else if (!.wm_na(lo)) paste0(cc$variable, " > ", lo)
  else if (!.wm_na(up)) paste0(cc$variable, " < ", up)
  else stop("write_mod(): the complementarity condition on variable '",
            cc$variable, "' has no bound.", call. = FALSE)
  paste0(" _|_ ", body)
}


.wm_value_block <- function(keyword, values) {
  if (length(values) == 0L) return(character(0))
  if (is.null(names(values)))
    stop("write_mod(): `", keyword, "` must be a named numeric vector.",
         call. = FALSE)
  body <- vapply(names(values), function(nm)
    paste0(nm, " = ", .wm_num(unname(values[nm]), paste0(keyword, " ", nm)), ";"),
    character(1), USE.NAMES = FALSE)
  c(paste0(keyword, ";"), body, "end;", "")
}


.wm_histval_block <- function(histval) {
  if (is.null(histval) || length(histval) == 0L) return(character(0))
  if (is.null(names(histval)))
    stop("write_mod(): `histval` must be a named list.", call. = FALSE)
  body <- character(0)
  for (nm in names(histval)) {
    v <- histval[[nm]]
    if (!is.numeric(v))
      stop("write_mod(): histval entry '", nm, "' is not numeric.",
           call. = FALSE)
    for (k in seq_along(v)) {
      if (is.na(v[k])) next
      ## lag k (k periods before the first simulation period) is written
      ## `name(1 - k)`: lag 1 -> name(0), lag 2 -> name(-1).  See parse_histval_block().
      body <- c(body, sprintf("%s(%d) = %s;", nm, 1L - k,
                              .wm_num(v[k], paste0("histval ", nm))))
    }
  }
  if (length(body) == 0L) return(character(0))
  c("histval;", body, "end;", "")
}


.wm_ssm_block <- function(ssm) {
  if (is.null(ssm) || length(ssm) == 0L) return(character(0))
  body <- vapply(ssm, function(a) {
    if (is.null(a$name) || is.null(a$expr))
      stop("write_mod(): malformed steady_state_model assignment.",
           call. = FALSE)
    paste0(a$name, " = ", .wm_expr(a$expr), ";")
  }, character(1), USE.NAMES = FALSE)
  c("steady_state_model;", body, "end;", "")
}


## Render every shock_groups block, one per entry of model$shock_groups_blocks
## (keyed by its `name=` option).  model$shock_groups is the first block and is
## therefore DERIVED -- rendering it too would emit the same block twice.
.wm_shock_groups_block <- function(blocks) {
  if (is.null(blocks) || length(blocks) == 0L) return(character(0))
  if (is.null(names(blocks)) || any(!nzchar(names(blocks))))
    stop("write_mod(): `shock_groups_blocks` must be a NAMED list (one entry ",
         "per shock_groups(name=...) block).", call. = FALSE)

  L <- character(0)
  for (nm in names(blocks)) {
    grp <- blocks[[nm]]
    if (length(grp) == 0L) next
    if (is.null(names(grp)) || any(!nzchar(names(grp))))
      stop("write_mod(): shock_groups block '", nm, "' must be a NAMED list ",
           "of shock memberships.", call. = FALSE)
    if (any(grepl("'", names(grp), fixed = TRUE)))
      stop("write_mod(): shock group label(s) in block '", nm, "' contain a ",
           "single quote, which cannot be written as .mod source: ",
           paste(grep("'", names(grp), fixed = TRUE, value = TRUE),
                 collapse = ", "), call. = FALSE)
    L <- c(L, paste0("shock_groups(name = ", nm, ");"))
    for (g in names(grp))
      L <- c(L, paste0("'", g, "' = ",
                       paste(as.character(grp[[g]]), collapse = ", "), ";"))
    L <- c(L, "end;", "")
  }
  L
}


.wm_shocks_block <- function(model) {
  v  <- model$shocks$variances
  cr <- model$shocks$correlations
  dt <- model$det_shocks
  if (!.wm_nonempty_df(v) && !.wm_nonempty_df(cr) && !.wm_nonempty_df(dt))
    return(character(0))

  allowed <- model$param_names
  if (.wm_nonempty_df(v)) {
    for (i in seq_len(nrow(v))) {
      .wm_check_expr_ids(v$stderr_expr[i],   allowed, paste0("stderr ", v$name[i]))
      .wm_check_expr_ids(v$variance_expr[i], allowed, paste0("var ", v$name[i]))
      .wm_check_expr_ids(v$skew_expr[i],     allowed, paste0("skew ", v$name[i]))
    }
  }
  if (.wm_nonempty_df(cr)) {
    for (i in seq_len(nrow(cr))) {
      lbl <- paste0("(", cr$var1[i], ", ", cr$var2[i], ")")
      .wm_check_expr_ids(cr$corr_expr[i], allowed, paste0("corr ", lbl))
      .wm_check_expr_ids(cr$cov_expr[i],  allowed, paste0("covar ", lbl))
    }
  }

  body <- character(0)

  if (.wm_nonempty_df(v)) {
    for (i in seq_len(nrow(v))) {
      nm <- v$name[i]
      has_sd  <- !.wm_na(v$stderr_expr[i])
      has_var <- !.wm_na(v$variance_expr[i])
      if (has_sd && has_var)
        stop("write_mod(): shock '", nm, "' carries both a stderr and a ",
             "variance expression; cannot render unambiguously.", call. = FALSE)
      if (has_sd) {
        body <- c(body, paste0("var ", nm, ";"),
                  paste0("stderr ", v$stderr_expr[i], ";"))
      } else if (has_var) {
        body <- c(body, paste0("var ", nm, " = ", v$variance_expr[i], ";"))
      } else if (!.wm_na(v$stderr[i])) {
        body <- c(body, paste0("var ", nm, ";"),
                  paste0("stderr ", .wm_num(v$stderr[i], paste0("stderr ", nm)),
                         ";"))
      } else if (!.wm_na(v$variance[i])) {
        body <- c(body, paste0("var ", nm, " = ",
                               .wm_num(v$variance[i], paste0("var ", nm)), ";"))
      } else if (!(!.wm_na(v$skew[i]) && v$skew[i] != 0)) {
        stop("write_mod(): shock '", nm, "' has neither a stderr nor a ",
             "variance to render.", call. = FALSE)
      }
      if (!.wm_na(v$skew[i]) && v$skew[i] != 0) {
        sk <- if (!.wm_na(v$skew_expr[i])) v$skew_expr[i]
              else .wm_num(v$skew[i], paste0("skew ", nm))
        body <- c(body, paste0("skew ", nm, " = ", sk, ";"))
      }
    }
  }

  if (.wm_nonempty_df(cr)) {
    for (i in seq_len(nrow(cr))) {
      is_corr <- !.wm_na(cr$corr[i]) || !.wm_na(cr$corr_expr[i])
      is_cov  <- !.wm_na(cr$cov[i])  || !.wm_na(cr$cov_expr[i])
      if (is_corr && is_cov)
        stop("write_mod(): shock pair (", cr$var1[i], ", ", cr$var2[i],
             ") carries both a correlation and a covariance; cannot render ",
             "unambiguously.", call. = FALSE)
      if (is_corr) {
        rhs <- if (!.wm_na(cr$corr_expr[i])) cr$corr_expr[i]
               else .wm_num(cr$corr[i], "corr")
        body <- c(body, paste0("corr ", cr$var1[i], ", ", cr$var2[i],
                               " = ", rhs, ";"))
      } else if (is_cov) {
        rhs <- if (!.wm_na(cr$cov_expr[i])) cr$cov_expr[i]
               else .wm_num(cr$cov[i], "covar")
        body <- c(body, paste0("var ", cr$var1[i], ", ", cr$var2[i],
                               " = ", rhs, ";"))
      } else {
        stop("write_mod(): shock pair (", cr$var1[i], ", ", cr$var2[i],
             ") has no correlation or covariance value to render.",
             call. = FALSE)
      }
    }
  }

  if (.wm_nonempty_df(dt)) {
    for (nm in unique(dt$name)) {
      sub <- dt[dt$name == nm, , drop = FALSE]
      body <- c(body,
                paste0("var ", nm, ";"),
                paste0("periods ", paste(as.integer(sub$period),
                                         collapse = ", "), ";"),
                paste0("values ",
                       paste(vapply(sub$value, .wm_num, character(1),
                                    what = paste0("det shock ", nm)),
                             collapse = ", "), ";"))
    }
  }

  c("shocks;", body, "end;", "")
}


## Dynare 7 `shock_paths` block from det_shocks (explicit periods) and the
## terminal values in endval (period `end`).  Values are the resolved numbers,
## so the block re-parses to the same det_shocks / endval whatever
## expressions, `self.` / `initval.` scopes or `a:end` ranges produced them.
.wm_shock_paths_block <- function(model) {
  if (is.null(model$shock_paths)) return(character(0))
  dt <- model$det_shocks
  tv <- model$endval
  nms <- unique(c(if (.wm_nonempty_df(dt)) dt$name, names(tv)))
  ctl <- .wm_controlled_stanzas(model$controlled_paths)
  if (length(nms) == 0L && length(ctl) == 0L) return(character(0))
  body <- character(0)
  for (nm in nms) {
    sub <- if (.wm_nonempty_df(dt)) dt[dt$name == nm, , drop = FALSE] else dt
    per <- as.character(as.integer(sub$period))
    val <- vapply(sub$value, .wm_num, character(1),
                  what = paste0("shock_paths ", nm))
    if (nm %in% names(tv)) {
      per <- c(per, "end")
      val <- c(val, .wm_num(unname(tv[[nm]]), paste0("shock_paths ", nm)))
    }
    body <- c(body, paste0("var ", nm, ";"),
              paste0("periods ", paste(per, collapse = ", "), ";"),
              paste0("values ", paste(val, collapse = ", "), ";"))
  }
  c("shock_paths;", body, ctl, "end;", "")
}


## Controlled-paths stanzas (`exogenize Y; periods ...; values ...;
## endogenize E;`), one per (exogenize, endogenize) pair with explicit
## periods and the resolved values, so the table re-parses identically.
.wm_controlled_stanzas <- function(cp) {
  if (!.wm_nonempty_df(cp)) return(character(0))
  key <- paste(cp$exogenize, cp$endogenize, sep = "\r")
  body <- character(0)
  for (k in unique(key)) {
    sub <- cp[key == k, , drop = FALSE]
    val <- vapply(sub$value, function(v) {
      s <- .wm_num(v, paste0("controlled path ", sub$exogenize[1L]))
      if (v < 0) paste0("(", s, ")") else s
    }, character(1))
    body <- c(body, paste0("exogenize ", sub$exogenize[1L], ";"),
              paste0("periods ", paste(as.integer(sub$period),
                                       collapse = ", "), ";"),
              paste0("values ", paste(val, collapse = ", "), ";"),
              paste0("endogenize ", sub$endogenize[1L], ";"))
  }
  body
}

## A standalone Dynare 7 perfect_foresight_controlled_paths block (a model
## with a shock_paths block carries its controlled stanzas there instead).
.wm_controlled_paths_block <- function(model) {
  if (!is.null(model$shock_paths)) return(character(0))
  body <- .wm_controlled_stanzas(model$controlled_paths)
  if (length(body) == 0L) return(character(0))
  c("perfect_foresight_controlled_paths;", body, "end;", "")
}


## One `estimated_params` / `estimated_params_init` row -> a .mod statement.
.wm_ep_row <- function(r, i) {
  nm <- switch(as.character(r$type),
               "parameter" = r$name,
               "stderr"    = paste0("stderr ", r$name),
               "skew"      = paste0("skew ", r$name),
               "corr"      = {
                 if (.wm_na(r$name2))
                   stop("write_mod(): estimated_params row ", i,
                        " is type 'corr' but has no second name.",
                        call. = FALSE)
                 paste0("corr ", r$name, ", ", r$name2)
               },
               stop("write_mod(): unknown estimated_params type '",
                    r$type, "' in row ", i, ".", call. = FALSE))

  has_init <- !.wm_na(r$init)
  has_lb   <- !.wm_na(r$lb)
  has_ub   <- !.wm_na(r$ub)
  if (has_lb != has_ub)
    stop("write_mod(): estimated_params row ", i,
         " has only one of (lb, ub); Dynare syntax requires both or neither.",
         call. = FALSE)
  if (!has_init && has_lb)
    stop("write_mod(): estimated_params row ", i,
         " has bounds but no init value; the `NAME, INIT, LB, UB` syntax ",
         "cannot express that.", call. = FALSE)

  before <- character(0)
  if (has_init) before <- .wm_num(r$init, "estimated_params init")
  if (has_lb)
    before <- c(before, .wm_num(r$lb, "estimated_params lb"),
                .wm_num(r$ub, "estimated_params ub"))

  ps <- c(r$p1, r$p2, r$p3, r$p4)
  keep <- which(!is.na(ps))
  n_after <- if (length(keep) == 0L) 0L else max(keep)
  ## Dynare allows EMPTY positional fields ("uniform_pdf, , , 0.0005, 0.5"
  ## puts the bounds in p3/p4 with p1/p2 blank), and parse_mod() stores
  ## those as NA. Render a leading/interior NA as an empty field so the
  ## round trip reproduces the original positions exactly.
  after <- if (n_after == 0L) character(0)
           else vapply(ps[seq_len(n_after)], function(v)
             if (is.na(v)) "" else .wm_num(v, "estimated_params prior parameter"),
             character(1))

  if (.wm_na(r$prior)) {
    if (n_after > 0L)
      stop("write_mod(): estimated_params row ", i,
           " has prior parameters but no prior shape.", call. = FALSE)
    if (length(before) == 0L) return(paste0(nm, ";"))
    return(paste0(paste(c(nm, before), collapse = ", "), ";"))
  }
  paste0(paste(c(nm, before, r$prior, after), collapse = ", "), ";")
}


.wm_ep_block <- function(df, keyword) {
  if (!.wm_nonempty_df(df)) return(character(0))
  missing_cols <- setdiff(c("type", "name", "name2", "prior",
                            "p1", "p2", "p3", "p4", "init", "lb", "ub"),
                          names(df))
  if (length(missing_cols) > 0L)
    stop("write_mod(): `", keyword, "` is missing column(s): ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  body <- vapply(seq_len(nrow(df)),
                 function(i) .wm_ep_row(df[i, , drop = FALSE], i),
                 character(1))
  head <- if (isTRUE(attr(df, "use_calibration")))
    paste0(keyword, "(use_calibration);") else paste0(keyword, ";")
  c(head, body, "end;", "")
}


## periods/values-style blocks (filter_tunes, heteroskedastic_shocks).
.wm_pv_block <- function(df, keyword, value_kw, extra_kw = NULL) {
  if (!.wm_nonempty_df(df)) return(character(0))
  body <- character(0)
  for (i in seq_len(nrow(df))) {
    pr <- df$periods[[i]]
    body <- c(body, paste0("var ", df$var[i], ";"),
              paste0("periods ",
                     paste(if (is.character(pr)) pr else as.integer(pr),
                           collapse = ", "), ";"),
              paste0(value_kw, " ",
                     paste(vapply(df[[value_kw]][[i]], .wm_num, character(1),
                                  what = paste0(keyword, " ", value_kw)),
                           collapse = ", "), ";"))
    if (!is.null(extra_kw)) {
      ex <- df[[extra_kw]][[i]]
      if (!is.null(ex) && length(ex) > 0L && !all(is.na(ex)))
        body <- c(body, paste0(extra_kw, " ",
                               paste(vapply(ex, .wm_num, character(1),
                                            what = extra_kw),
                                     collapse = ", "), ";"))
    }
  }
  c(paste0(keyword, ";"), body, "end;", "")
}


## heteroskedastic_shocks periods are stored as SAMPLE periods (1 = first data
## row); in the .mod they index the original dataset (integer p = sample
## period + first_obs - 1) or are dates when estimation(first_obs=) is a date
## (see parse_heteroskedastic_shocks_block()).  Map them back so the written
## file re-parses to the same sample periods.
.wm_hs_source_periods <- function(model) {
  df <- model$heteroskedastic_shocks$scales
  if (!.wm_nonempty_df(df)) return(df)
  fo <- NULL
  for (cmd in model$commands)
    if (identical(cmd$name, "estimation") && !is.null(cmd$options$first_obs))
      fo <- cmd$options$first_obs
  if (is.null(fo)) return(df)
  if (is.numeric(fo)) {
    df$periods <- lapply(df$periods, function(p) as.integer(p) + as.integer(fo) - 1L)
  } else {
    df$periods <- lapply(df$periods, .period_to_date, sample_start = as.character(fo))
  }
  df
}


.wm_sv_block <- function(sv_spec) {
  df <- sv_spec$sv
  if (!.wm_nonempty_df(df)) return(character(0))
  body <- character(0)
  for (i in seq_len(nrow(df))) {
    body <- c(body, paste0("var ", df$shock[i], ";"))
    for (k in c("mu", "rho", "sigma_eta")) {
      v <- df[[k]][[i]]
      rhs <- if (is.character(v)) v else .wm_num(v, paste0("sv ", k))
      body <- c(body, paste0(k, " = ", rhs, ";"))
    }
  }
  c("stochastic_volatility;", body, "end;", "")
}


.wm_occbin_block <- function(specs) {
  if (is.null(specs) || length(specs) == 0L) return(character(0))
  body <- character(0)
  for (s in specs) {
    body <- c(body, sprintf("name '%s';", s$name))
    ## `op` is stored FLIPPED relative to the source bind clause
    ## (parse_occbin_constraints_block: "bind <= " -> op ">"), so invert it back.
    bind_op <- if (identical(s$op, ">")) "<" else ">"
    bound <- if (!.wm_na(s$bound)) .wm_num(s$bound, "occbin bound")
             else if (!.wm_na(s$bound_expr)) s$bound_expr
             else stop("write_mod(): occbin constraint '", s$name,
                       "' has no bound to render.", call. = FALSE)
    body <- c(body, sprintf("bind %s %s %s;", s$var_name, bind_op, bound))
    if (!is.null(s$relax_str) && nzchar(s$relax_str))
      body <- c(body, sprintf("relax %s;", s$relax_str))
    if (!is.null(s$bind_eqs) && nzchar(s$bind_eqs))
      body <- c(body, "equations;", s$bind_eqs, "end;")
  }
  c("occbin_constraints;", body, "end;", "")
}


## OSR inputs: `optim_weights` (variance `y w;` / covariance `y, pie w;`),
## `osr_params` and `osr_params_bounds` (`param, lower, upper;`). Weights and
## bounds are expression text, written back verbatim.
.wm_osr_block <- function(spec) {
  if (is.null(spec)) return(character(0))
  out <- character(0)
  w <- spec$weights
  if (.wm_nonempty_df(w)) {
    body <- ifelse(w$var1 == w$var2,
                   sprintf("%s %s;", w$var1, w$expr),
                   sprintf("%s, %s %s;", w$var1, w$var2, w$expr))
    out <- c(out, "optim_weights;", body, "end;", "")
  }
  if (length(spec$params) > 0L)
    out <- c(out, paste0("osr_params ", paste(spec$params, collapse = " "), ";"),
             "")
  b <- spec$bounds
  if (.wm_nonempty_df(b))
    out <- c(out, "osr_params_bounds;",
             sprintf("%s, %s, %s;", b$name, b$lower, b$upper), "end;", "")
  out
}


## ramsey_constraints block: `VARIABLE > BOUND;` / `VARIABLE < BOUND;`, the
## bound written back verbatim.
.wm_ramsey_constraints_block <- function(rc) {
  if (!.wm_nonempty_df(rc)) return(character(0))
  c("ramsey_constraints;", sprintf("%s %s %s;", rc$var, rc$op, rc$bound),
    "end;", "")
}


## Dynare EXECUTES command lines in file order, so they are emitted in a
## sane execution order (policy setup -> steady/check -> simulation), not in
## the order parse_mod() happened to discover them.  This is free: the
## parser's own command loop has a FIXED keyword order, so `model$commands`
## comes back in the same order however the file is laid out.
.WM_CMD_ORDER <- c(
  "ramsey_model", "ramsey_policy", "discretionary_policy", "osr",
  "steady", "check", "model_info", "model_diagnostics",
  "perfect_foresight_setup", "perfect_foresight_solver", "simul",
  "extended_path",
  "stoch_simul", "estimation", "calib_smoother", "shock_decomposition",
  "forecast", "conditional_forecast"
)

.wm_commands <- function(model) {
  out <- character(0)
  po <- model$planner_objective
  if (!is.null(po) && nzchar(po$text %||% ""))
    out <- c(out, paste0("planner_objective ", po$text, ";"))

  cmds <- model$commands
  if (length(cmds) > 0L) {
    nms  <- vapply(cmds, function(c) c$name, character(1))
    rank <- match(nms, .WM_CMD_ORDER, nomatch = length(.WM_CMD_ORDER) + 1L)
    cmds <- cmds[order(rank, seq_along(rank))]
  }

  for (cmd in cmds) {
    ## planner_objective is emitted from model$planner_objective above; the
    ## commands list only carries the (name-stripped) echo of the same line.
    if (identical(cmd$name, "planner_objective")) next
    opts <- .wm_options_str(cmd$options)
    vars <- if (length(cmd$var_list) > 0L)
      paste0(" ", paste(cmd$var_list, collapse = " ")) else ""
    out <- c(out, paste0(cmd$name,
                         if (nzchar(opts)) paste0("(", opts, ")") else "",
                         vars, ";"))
  }
  if (length(out) > 0L) out <- c(out, "")
  out
}


# ---------------------------------------------------------------------------
# write_mod()
# ---------------------------------------------------------------------------

#' Write a parsed dynhr model back to Dynare .mod source
#'
#' Renders a \code{dynhr_mod} object (as returned by \code{\link{parse_mod}})
#' as Dynare \code{.mod} source text, in Dynare's canonical block order:
#' declarations, parameter values, \code{model}, \code{initval} /
#' \code{endval} / \code{histval} / \code{steady_state_model}, \code{shocks},
#' \code{estimated_params}, \code{filter_tunes},
#' \code{heteroskedastic_shocks}, \code{stochastic_volatility},
#' \code{occbin_constraints}, \code{observation_trends},
#' \code{deterministic_trends}, the OSR inputs (\code{optim_weights},
#' \code{osr_params}, \code{osr_params_bounds}), and the command lines.
#'
#' @details
#' \strong{What round-trips.} The exporter is faithful to the \emph{parsed
#' model}, not to the original file text: the guarantee it provides is
#' \code{parse_mod(write_mod(parse_mod(f)))} equals \code{parse_mod(f)} for
#' every field the parser populates.  Comments, whitespace and the original
#' operator spelling are not preserved.
#'
#' \strong{Macro expansion.}  \code{parse_mod()} expands \code{@#include},
#' \code{@#for}, \code{@#if} and \code{@#define} before anything else, so a
#' \code{dynhr_mod} never carries the unexpanded directives.  \code{write_mod()}
#' therefore emits the \emph{expanded} model and marks the result with
#' \code{attr(txt, "macro_expanded") = TRUE}.
#'
#' \strong{Two deliberate normalisations}, both semantics-preserving, both
#' announced in a comment in the emitted file:
#' \itemize{
#'   \item \emph{Auxiliary variables.}  Leads/lags beyond one period (and any
#'     non-zero timing on an exogenous variable) have already been rewritten
#'     into \code{AUX_*} variables by the parser; those are declared and
#'     emitted as ordinary endogenous variables.
#'   \item \emph{predetermined_variables.}  The \code{-1} re-timing that the
#'     declaration implies has already been applied to the equation ASTs, so
#'     the emitted model states the shifted timing explicitly and does
#'     \strong{not} re-emit the \code{predetermined_variables} line (which
#'     would shift a second time).
#'   \item \emph{Nonstationary models.}  A model declared with
#'     \code{trend_var} / \code{log_trend_var} and \code{var(deflator=...)} /
#'     \code{var(log_deflator=...)} has already been detrended by the parser;
#'     the emitted model is the stationary (detrended) one, the trend
#'     declarations are \strong{not} re-emitted (which would detrend a second
#'     time) and \code{model$nonstationary} is absent after the round trip.
#'     A \code{deterministic_trends} block is written back as such.
#'   \item \emph{Exogenous auxiliary chains.}  A non-zero timing on an
#'     exogenous variable is rewritten by the parser into an
#'     \code{AUX_EXO_LEAD_*}/\code{AUX_EXO_LAG_*} chain, and that rewrite is
#'     not idempotent, so \code{write_mod()} folds the chain back to a plain
#'     \code{x(+k)} / \code{x(-k)} reference and lets \code{parse_mod()}
#'     rebuild it.  The rebuilt chain is identical, in the same order.
#'   \item \emph{var(log).}  The parser rewrites every \code{x} of a
#'     \code{var(log) x;} variable into \code{exp(LOG_x)} and appends
#'     \code{x = exp(LOG_x);}; \code{write_mod()} folds that back and
#'     re-declares \code{x} under \code{var(log)}.  The \code{LOG_x = log(x)}
#'     line the parser adds to a \code{steady_state_model} (and the
#'     \code{LOG_x} starting values it adds to \code{initval}) are written
#'     explicitly and read back unchanged.  A \code{var(log)} variable that
#'     ends up inside an auxiliary lead/lag chain cannot be folded back and is
#'     an error.
#'   \item \emph{Multiple shocks blocks.}  Several \code{shocks} blocks are
#'     merged by the parser into one \code{model$shocks}; that merged block is
#'     what gets written, so the per-block \code{model$shocks_blocks}
#'     breakdown collapses to a single entry on the next parse.  The effective
#'     shock structure is unchanged.
#' }
#'
#' \strong{Command order.}  Dynare executes command lines in file order, so
#' they are emitted in execution order (policy setup, then
#' \code{steady}/\code{check}, then simulation/estimation) rather than in the
#' order \code{parse_mod()} discovered them.  This costs nothing: the parser
#' scans a fixed keyword list, so \code{model$commands} comes back in the same
#' order however the file is laid out.
#'
#' \strong{Fail-loud.}  Any field present on the object that this exporter
#' does not know how to render (or account for as derived) raises an error
#' naming the field, as does any block content that cannot be expressed in
#' \code{.mod} syntax.  A silently lossy exporter is worse than none.
#'
#' \strong{Not stored by the parser} (and hence not emitted):
#' \code{estimated_params_bounds} and the other blocks that
#' \code{parse_mod()} strips without parsing (\code{observation_trends} IS
#' stored and written back).  \code{estimation(...)} options
#' that contain a parenthesis (e.g. \code{optim=('MaxIter',200)}) are captured
#' only partially by the parser; \code{write_mod()} reproduces exactly what
#' was captured, so the round trip is stable but the emitted command line can
#' be shorter than the original.
#'
#' @param model A \code{dynhr_mod} object.
#' @param file Optional path to write to.  When \code{NULL} (default) nothing
#'   is written and the text is returned.
#' @param header Logical; prepend a provenance comment header (default
#'   \code{TRUE}).
#'
#' @return A length-1 character string holding the complete \code{.mod} source,
#'   with attribute \code{macro_expanded = TRUE}.  Returned invisibly when
#'   \code{file} is given.
#'
#' @seealso \code{\link{parse_mod}}, \code{\link{smoother2histval}}
#'
#' @examples
#' m   <- parse_mod(system.file("extdata", "models", "rbc.mod",
#'                              package = "dynhr"))
#' txt <- write_mod(m)
#' m2  <- parse_mod(txt)
#' identical(m$var_names, m2$var_names)
#'
#' @export
write_mod <- function(model, file = NULL, header = TRUE) {
  if (!inherits(model, "dynhr_mod"))
    stop("write_mod(): `model` must be a dynhr_mod object (see parse_mod()).",
         call. = FALSE)

  unknown <- setdiff(names(model), .WM_KNOWN_FIELDS)
  if (length(unknown) > 0L)
    stop("write_mod(): the model carries field(s) this exporter cannot ",
         "render: ", paste(unknown, collapse = ", "),
         ". Refusing to write a silently lossy .mod file.", call. = FALSE)

  model <- .wm_unexpand_exo_aux(model)
  model <- .wm_unexpand_log_vars(model)

  L <- character(0)

  if (isTRUE(header)) {
    L <- c(L,
           "// Generated by dynhr::write_mod() -- do not hand-edit.",
           "// Dynare macro directives (@#include / @#for / @#define) are",
           "// already expanded; auxiliary AUX_* variables for |lead/lag| > 1",
           "// are declared explicitly.",
           if (length(model$predetermined_vars) > 0L)
             paste0("// predetermined_variables (",
                    paste(model$predetermined_vars, collapse = " "),
                    ") were re-timed by -1 into the equations below and are",
                    " NOT re-declared."),
           if (!is.null(model$nonstationary))
             paste0("// The model is written DETRENDED: the trend_var / ",
                    "log_trend_var (",
                    paste(model$nonstationary$trend_vars$name, collapse = " "),
                    ") and var(deflator=) / var(log_deflator=) (",
                    paste(model$nonstationary$deflated$name, collapse = " "),
                    ") declarations are NOT re-emitted."),
           "")
  }

  L <- c(L, .wm_metadata_block(model))

  ann <- model$metadata$annotations %||% list()
  L <- c(L,
         .wm_var_lines(model$var_names, names(model$log_vars), ann),
         .wm_decl_line("varexo",     model$varexo_names,     ann),
         .wm_decl_line("varexo_det", model$varexo_det_names, ann),
         .wm_decl_line("parameters", model$param_names,      ann),
         .wm_decl_line("varobs",     model$varobs))
  L <- c(L, "")

  pv <- model$param_values
  if (length(pv) > 0L) {
    if (is.null(names(pv)))
      stop("write_mod(): `param_values` must be a named numeric vector.",
           call. = FALSE)
    for (nm in names(pv)) {
      if (is.na(pv[[nm]])) next          # declared but never valued
      L <- c(L, paste0(nm, " = ", .wm_num(unname(pv[nm]), nm), ";"))
    }
    L <- c(L, "")
  }

  ## A shock_paths model keeps its scenario in det_shocks + endval (see
  ## parse_shock_paths_blocks()); it is written back as a shock_paths block
  ## with explicit periods, and neither as shocks periods nor as endval.
  sp_model <- model
  if (!is.null(model$shock_paths)) {
    sp_model$det_shocks <- model$det_shocks[0L, , drop = FALSE]
    sp_model$endval     <- numeric(0)
  }

  L <- c(L, .wm_model_block(model))
  L <- c(L, .wm_value_block("initval", model$initval))
  L <- c(L, .wm_value_block("endval",  sp_model$endval))
  L <- c(L, .wm_histval_block(model$histval))
  L <- c(L, .wm_ssm_block(model$steady_state_model))
  L <- c(L, .wm_shocks_block(sp_model))
  L <- c(L, .wm_shock_paths_block(model))
  L <- c(L, .wm_controlled_paths_block(model))
  L <- c(L, .wm_shock_groups_block(model$shock_groups_blocks))
  L <- c(L, .wm_ep_block(model$estimated_params, "estimated_params"))
  L <- c(L, .wm_ep_block(model$estimated_params_init, "estimated_params_init"))
  L <- c(L, .wm_pv_block(model$filter_tunes$tunes, "filter_tunes",
                         "values", extra_kw = "stderr"))
  L <- c(L, .wm_pv_block(.wm_hs_source_periods(model),
                         "heteroskedastic_shocks", "scales"))
  L <- c(L, .wm_sv_block(model$stochastic_volatility))
  L <- c(L, .wm_occbin_block(model$occbin_constraints))
  L <- c(L, .wm_obs_trends_block(model$observation_trends,
                                 model$deterministic_trends))
  L <- c(L, .wm_osr_block(model$osr))
  L <- c(L, .wm_commands(model))
  ## after the commands: Dynare reads the block as part of the Ramsey
  ## problem that ramsey_model / ramsey_policy declares
  L <- c(L, .wm_ramsey_constraints_block(model$ramsey_constraints))

  txt <- paste0(paste(L, collapse = "\n"), "\n")
  attr(txt, "macro_expanded") <- TRUE

  if (!is.null(file)) {
    writeLines(L, con = file)
    return(invisible(txt))
  }
  txt
}


# ---------------------------------------------------------------------------
# smoother2histval()
# ---------------------------------------------------------------------------

#' Turn a smoother result into a histval history
#'
#' Dynare's \code{smoother2histval} command writes the smoothed state at a
#' chosen period into a \code{histval} block so a perfect-foresight or
#' forecast run can start from the estimated history.  This is the R
#' equivalent: it returns the same named-list structure that
#' \code{parse_mod()} stores in \code{model$histval} and that
#' \code{\link{write_mod}} renders.
#'
#' @details
#' The returned list maps a variable name to a numeric vector indexed by
#' \emph{lag}: element \code{k} is the value \code{k} periods before the first
#' simulation period.  With \code{period = p} and \code{lags = 2}, element 1 is
#' the smoothed value at \code{p} and element 2 the value at \code{p - 1}, so
#' a simulation started at \code{p + 1} sees the estimated history.
#'
#' @param smoother_result Either a list carrying a \code{smoothed_states}
#'   matrix (as returned by \code{kalman_smoother()} or
#'   \code{run_full_estimation()$smoother}) or the matrix itself.  Both
#'   orientations are accepted: \code{T x n_state} (variables in columns, the
#'   \code{kalman_smoother()} convention) and \code{n_state x T} (variables in
#'   rows).
#' @param period Integer period whose smoothed state becomes lag 1.  Defaults
#'   to the last period in the smoother output.
#' @param vars Optional character vector selecting (and ordering) the
#'   variables to keep.  Defaults to all of them.
#' @param lags Number of lags to record (default 1).  \code{lags = 2} also
#'   records \code{period - 1}, and so on.
#'
#' @return A named list of numeric vectors of length \code{lags}, suitable for
#'   \code{model$histval} and for \code{\link{write_mod}}.
#'
#' @seealso \code{\link{write_mod}}, \code{\link{parse_mod}}
#' @export
smoother2histval <- function(smoother_result, period = NULL, vars = NULL,
                             lags = 1L) {
  S <- if (is.matrix(smoother_result)) smoother_result
       else if (is.list(smoother_result) &&
                !is.null(smoother_result$smoothed_states))
         smoother_result$smoothed_states
       else if (is.list(smoother_result) &&
                is.list(smoother_result$smoother) &&
                !is.null(smoother_result$smoother$smoothed_states))
         smoother_result$smoother$smoothed_states
       else NULL
  if (!is.matrix(S))
    stop("smoother2histval(): could not find a `smoothed_states` matrix in ",
         "`smoother_result`.", call. = FALSE)

  ## Orientation: kalman_smoother() names the COLUMNS (T x n_state); the OBC
  ## smoother names the ROWS (n_state x T).
  by_col <- !is.null(colnames(S))
  if (by_col && !is.null(rownames(S)) && is.null(colnames(S))) by_col <- FALSE
  if (!by_col && is.null(rownames(S)))
    stop("smoother2histval(): `smoothed_states` has no dimnames, so the ",
         "variable names cannot be recovered.", call. = FALSE)
  if (!by_col) S <- t(S)          # normalise to T x n_state

  nT <- nrow(S)
  if (is.null(period)) period <- nT
  period <- as.integer(period)
  lags   <- as.integer(lags)
  if (length(period) != 1L || is.na(period) || period < 1L || period > nT)
    stop("smoother2histval(): `period` must be in 1..", nT, ".", call. = FALSE)
  if (length(lags) != 1L || is.na(lags) || lags < 1L)
    stop("smoother2histval(): `lags` must be a positive integer.",
         call. = FALSE)
  if (period - lags + 1L < 1L)
    stop("smoother2histval(): `lags` = ", lags, " reaches before period 1 ",
         "for `period` = ", period, ".", call. = FALSE)

  nms <- colnames(S)
  if (!is.null(vars)) {
    miss <- setdiff(vars, nms)
    if (length(miss) > 0L)
      stop("smoother2histval(): variable(s) not in the smoother output: ",
           paste(miss, collapse = ", "), ".", call. = FALSE)
    nms <- vars
  }

  out <- lapply(nms, function(nm) {
    v <- vapply(seq_len(lags), function(k) S[period - k + 1L, nm], numeric(1))
    unname(v)
  })
  names(out) <- nms
  out
}
