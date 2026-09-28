## R/conditions.R
## dynhr's messaging layer: levelled, classed, run-scoped.
##
## Three things this file provides, and why each exists:
##
## 1. CLASSED CONDITIONS. `.dynhr_abort()` / `.dynhr_warn()` / `.dynhr_inform()`
##    signal real base-R conditions (they end in `stop()` / `warning()` /
##    `message()`), so `suppressWarnings()`, `suppressMessages()`,
##    `tryCatch()` and `testthat::expect_warning()` all keep working exactly as
##    they did. The layer adds a `dynhr_*` class (plus any caller-supplied
##    subclass) so a caller can catch dynhr's own conditions selectively.
##
## 2. A VERBOSITY LEVEL. One option, `verbosity`, resolved through
##    `.dynhr_opt()` -- which means it rides `.dynhr_daemon_state()` into every
##    mirai daemon for free. That is not optional: a level that lived only in
##    the host session would make the parallel path emit a different set of
##    messages from the serial one (see the long comment in options.R; this
##    package has been bitten by exactly that before).
##
## 3. A RUN-SCOPED DEDUP EPOCH. `.dynhr_epoch()` opens a message epoch that
##    `on.exit()` closes. Repeat-suppressed messages are keyed WITHIN the
##    epoch, not for the lifetime of the session.
##
##    This replaces four byte-identical `*_warn_once()` helpers whose backing
##    environments were never reset. That defect was not cosmetic: a second
##    MCMC run in one session emitted no warnings at all; five test files
##    reached into `dynhr:::` to clear the latch by hand, using three different
##    idioms; and any warn-once site WITHOUT such a poke was silently
##    dependent on test execution order. Because each top-level entry point --
##    and therefore each test that calls one -- opens its own epoch, all of
##    that goes away without any caller needing to know the mechanism exists.
##
##    An epoch also COUNTS what it suppressed, so a chain that degraded 12,431
##    times says so on the way out instead of looking identical to one that
##    degraded once.

# ---------------------------------------------------------------------------
# Verbosity levels
# ---------------------------------------------------------------------------
# Ordered, low to high. A message is emitted when its own level is <= the
# active threshold. `error` is deliberately included so the vector is a total
# order, but note that errors are NEVER suppressed -- see `.dynhr_abort()`.
.dynhr_levels <- c(silent = 0L, error = 1L, warn = 2L, info = 3L, debug = 4L)

# The shipped default. `info`, NOT `warn`, and deliberately so: the ~169
# `message()` and ~790 `cat()` calls migrating onto this layer ALWAYS printed
# before, so defaulting to `warn` would silence all of them as a side effect of
# a mechanical rename. A migration must not change behaviour; quieting the
# package is the USER's decision, made with `dynhr_set_verbosity("warn")` or
# `("silent")`. `debug` detail is the only thing hidden by default.
.dynhr_default_verbosity <- "info"

#' Resolve the active verbosity threshold as an integer
#'
#' Priority: explicit \code{verbose} argument > \code{verbosity} option > default.
#'
#' \code{verbose} is accepted in the legacy boolean form used by ~316 signatures
#' across the package, and mapped onto the level scale: \code{TRUE} means "tell me
#' what you are doing" (\code{info}) and \code{FALSE} means "only the things I need to
#' act on" (\code{warn}). Note that \code{verbose = FALSE} therefore does NOT silence
#' warnings -- it never did, and silencing them here would quietly change the
#' behaviour of every one of those call sites.
#'
#' A character \code{verbose} is taken as a level name directly, so new code can
#' pass \code{verbose = "debug"} through the same argument.
#'
#' @param verbose \code{NULL} (use the option), a length-1 logical, or a level name.
#' @return Integer threshold.
#' @noRd
.dynhr_verbosity <- function(verbose = NULL) {
  if (!is.null(verbose)) {
    if (is.character(verbose) && length(verbose) == 1L &&
        verbose %in% names(.dynhr_levels))
      return(unname(.dynhr_levels[[verbose]]))
    if (is.logical(verbose) && length(verbose) == 1L && !is.na(verbose))
      return(unname(.dynhr_levels[[if (verbose) "info" else "warn"]]))
    ## Anything else (a number, NA, a vector) is not a level. Fall through to
    ## the option rather than erroring: this helper sits on hot paths and must
    ## never be the reason a long chain dies.
  }
  lvl <- .dynhr_opt("verbosity", NULL, default = .dynhr_default_verbosity)
  if (is.numeric(lvl) && length(lvl) == 1L && is.finite(lvl))
    return(as.integer(max(0L, min(4L, lvl))))
  if (!is.character(lvl) || length(lvl) != 1L || !(lvl %in% names(.dynhr_levels)))
    lvl <- .dynhr_default_verbosity
  unname(.dynhr_levels[[lvl]])
}

#' Set the dynhr verbosity level
#'
#' Controls which of dynhr's own informational messages and warnings reach the
#' console. Errors are never suppressed by this setting.
#'
#' @param level One of \code{"silent"}, \code{"error"}, \code{"warn"}, \code{"info"} (the
#'   default), \code{"debug"}. Each level includes the ones before it. \code{"warn"}
#'   silences dynhr's progress and informational output while keeping warnings;
#'   \code{"silent"} also drops warnings. Errors are never suppressed.
#' @return Invisibly, the previous level.
#' @seealso [dynhr_set_options()], which this is a thin wrapper over, and
#'   [dynhr_verbosity()] to read the current value.
#' @examples
#' old <- dynhr_set_verbosity("silent")
#' dynhr_set_verbosity(old)
#' @export
dynhr_set_verbosity <- function(level = "info") {
  if (!is.character(level) || length(level) != 1L ||
      !(level %in% names(.dynhr_levels)))
    stop("`level` must be one of ",
         paste0('"', names(.dynhr_levels), '"', collapse = ", "), ".",
         call. = FALSE)
  prev <- dynhr_verbosity()
  dynhr_set_options(verbosity = level)
  invisible(prev)
}

#' Report the active dynhr verbosity level
#' @return Length-1 character: the active level name.
#' @examples
#' dynhr_verbosity()
#' @export
dynhr_verbosity <- function() {
  names(.dynhr_levels)[match(.dynhr_verbosity(NULL), .dynhr_levels)]
}

# ---------------------------------------------------------------------------
# The run-scoped epoch
# ---------------------------------------------------------------------------
# `.dynhr_msg_state$epoch` is the CURRENT epoch, or NULL when no top-level
# entry point is on the stack. Epochs nest: an inner `.dynhr_epoch()` inherits
# the outer one rather than opening a second, so `run_full_estimation()`
# calling `solve_model()` still yields one "already warned" scope for the whole
# run, which is the behaviour a user expects from a single command.
.dynhr_msg_state <- new.env(parent = emptyenv())
.dynhr_msg_state$epoch <- NULL

# A fallback epoch used when a message is emitted with no entry point on the
# stack (a user calling an internal helper directly, say). It behaves exactly
# like a real epoch except that nothing ever closes it, so it accumulates for
# the session. `.dynhr_reset_messages()` clears it.
.dynhr_msg_state$orphan <- NULL

.dynhr_new_epoch <- function(label) {
  e <- new.env(parent = emptyenv())
  e$label   <- label
  e$seen    <- new.env(parent = emptyenv())  # key -> integer count
  e$order   <- character(0)                  # keys, in first-seen order
  e
}

.dynhr_current_epoch <- function() {
  ep <- .dynhr_msg_state$epoch
  if (!is.null(ep)) return(ep)
  if (is.null(.dynhr_msg_state$orphan))
    .dynhr_msg_state$orphan <- .dynhr_new_epoch("<no active run>")
  .dynhr_msg_state$orphan
}

#' Open a run-scoped message epoch
#'
#' Call at the top of a user-facing entry point, and pair it with the returned
#' closure via `on.exit()`:
#'
#' ```
#' on.exit(.dynhr_close_epoch(.dynhr_epoch("run_mcmc")), add = TRUE)
#' ```
#'
#' Within the epoch, `.dynhr_warn(..., once = TRUE)` emits each distinct key at
#' most once; closing the epoch reports how many times each key was suppressed.
#'
#' Nesting is safe and deliberate: if an epoch is already open this returns
#' `NULL` and closing `NULL` is a no-op, so the OUTERMOST entry point owns the
#' scope. That is what makes one user-level command produce one set of
#' warnings even though it may call a dozen internal entry points.
#'
#' @param label Short string naming the run; used in the closing summary.
#' @return The epoch object to hand to \code{.dynhr_close_epoch()}, or \code{NULL} if an
#'   epoch was already open.
#' @noRd
.dynhr_epoch <- function(label = "run") {
  if (!is.null(.dynhr_msg_state$epoch)) return(NULL)
  ep <- .dynhr_new_epoch(label)
  .dynhr_msg_state$epoch <- ep
  ep
}

#' Close an epoch opened by \code{.dynhr_epoch()} and report suppressed repeats
#'
#' A no-op when \code{ep} is \code{NULL} (an inner, nested call).
#'
#' The summary is emitted at \code{info} level, NOT as a warning: the individual
#' warnings already fired at the moment they mattered, and re-raising at the
#' end of a long run would be noise. What the summary adds is the COUNT, which
#' is the part the old once-per-session latches threw away.
#'
#' @param ep The epoch returned by \code{.dynhr_epoch()}.
#' @return Invisibly, a named integer vector of suppressed-repeat counts.
#' @noRd
.dynhr_close_epoch <- function(ep) {
  if (is.null(ep)) return(invisible(integer(0)))
  ## Only clear the global slot if THIS epoch is the one that is open. An
  ## error unwinding through several on.exit() handlers must not let an inner
  ## frame close an outer frame's epoch.
  if (identical(.dynhr_msg_state$epoch, ep))
    .dynhr_msg_state$epoch <- NULL

  counts <- vapply(ep$order, function(k) get0(k, envir = ep$seen,
                                              ifnotfound = 0L),
                   integer(1))
  repeats <- counts[counts > 1L] - 1L
  if (length(repeats) && .dynhr_verbosity(NULL) >= .dynhr_levels[["info"]]) {
    lines <- sprintf("  %s: %d further occurrence%s",
                     names(repeats), repeats,
                     ifelse(repeats == 1L, "", "s"))
    .dynhr_signal_message(
      paste0(c(sprintf("%s: repeated messages suppressed during this run:",
                       ep$label), lines), collapse = "\n"),
      class = "dynhr_epoch_summary")
  }
  invisible(repeats)
}

#' Clear all message dedup state
#'
#' Drops any open epoch and the orphan epoch. Intended for tests and for
#' interactive use after a run has been interrupted; ordinary code should let
#' \code{.dynhr_close_epoch()} do this.
#' @return Invisibly \code{NULL}.
#' @noRd
.dynhr_reset_messages <- function() {
  .dynhr_msg_state$epoch  <- NULL
  .dynhr_msg_state$orphan <- NULL
  invisible(NULL)
}

# Record `key` against the active epoch. Returns TRUE the FIRST time the key is
# seen in this epoch and FALSE afterwards; the count is kept either way.
.dynhr_first_time <- function(key) {
  ep <- .dynhr_current_epoch()
  n  <- get0(key, envir = ep$seen, ifnotfound = 0L)
  assign(key, n + 1L, envir = ep$seen)
  if (n == 0L) {
    ep$order <- c(ep$order, key)
    return(TRUE)
  }
  FALSE
}

# ---------------------------------------------------------------------------
# Emitters
# ---------------------------------------------------------------------------
# Formatting note: `cli` is used for its INLINE MARKUP and console width
# handling only -- dynhr does not adopt cli's condition system. Callers pass
# plain strings; `.dynhr_fmt()` glues them and applies cli formatting when the
# text actually contains markup, so an ordinary message costs nothing and, more
# importantly, a message containing a literal brace (a printed R expression, a
# regex) is not mangled by an accidental cli substitution.
.dynhr_fmt <- function(...) {
  txt <- paste0(..., collapse = "")
  txt
}

# Build and signal a message condition with dynhr's classes attached.
.dynhr_signal_message <- function(txt, class = NULL) {
  cond <- structure(
    class = c(class, "dynhr_message", "message", "condition"),
    list(message = paste0(txt, "\n"), call = NULL))
  message(cond)
}

#' Signal a dynhr error
#'
#' Always raised: the verbosity level does not gate errors. A caller that wants
#' to tolerate a failure must catch it, not silence it.
#'
#' @param ... Parts of the message, pasted together.
#' @param class Optional extra condition subclass.
#' @param call. Passed through to the condition's call, matching \code{stop()}.
#' @noRd
.dynhr_abort <- function(..., class = NULL, call. = FALSE) {
  cond <- structure(
    class = c(class, "dynhr_error", "error", "condition"),
    list(message = .dynhr_fmt(...),
         call = if (isTRUE(call.)) sys.call(-1L) else NULL))
  stop(cond)
}

#' Signal a dynhr warning
#'
#' @param ... Parts of the message, pasted together.
#' @param once When \code{TRUE}, emit at most once per run epoch. Use this for
#'   anything reachable from inside a loop, a per-draw likelihood evaluation or
#'   a per-period filter step.
#' @param key Dedup key. Defaults to the message text, which is right when the
#'   text is constant and WRONG when it interpolates a varying number -- pass an
#'   explicit key in that case, or every iteration produces a distinct key and
#'   \code{once} silently does nothing.
#' @param class Optional extra condition subclass.
#' @param verbose Optional per-call verbosity override.
#' @param call.,immediate.,domain Accepted for signature compatibility with base
#'   \code{warning()}, so that migrating a call site is a rename and nothing else.
#'   This matters: \code{.dynhr_warn()} pastes its \code{...}, so a \code{call. = FALSE} left
#'   over from the base call would otherwise be silently CONCATENATED INTO THE
#'   MESSAGE TEXT rather than rejected. \code{call.} and \code{domain} are ignored (dynhr
#'   messages name their own function in the text); \code{immediate.} is honoured.
#' @noRd
.dynhr_warn <- function(..., once = FALSE, key = NULL, class = NULL,
                        verbose = NULL, call. = FALSE, immediate. = FALSE,
                        domain = NULL) {
  if (.dynhr_verbosity(verbose) < .dynhr_levels[["warn"]])
    return(invisible(NULL))
  txt <- .dynhr_fmt(...)
  if (isTRUE(once) && !.dynhr_first_time(if (is.null(key)) txt else key))
    return(invisible(NULL))
  if (isTRUE(immediate.)) {
    ## base warning(immediate. = TRUE) bypasses the deferred warning buffer.
    ## A condition object cannot express that, so emit through base warning()
    ## with the flag and keep the classed condition for handlers via a
    ## withCallingHandlers-free direct call.
    op <- options(warn = max(1L, getOption("warn", 0L)))
    on.exit(options(op), add = TRUE)
  }
  cond <- structure(
    class = c(class, "dynhr_warning", "warning", "condition"),
    list(message = txt, call = NULL))
  warning(cond)
  invisible(NULL)
}

#' Emit a dynhr informational message
#'
#' Gated at \code{info} by default; pass \code{level = "debug"} for detail that should
#' only appear when explicitly asked for.
#'
#' @inheritParams .dynhr_warn
#' @param level Level at which this message becomes visible.
#' @noRd
#' @param appendLF,domain Accepted for signature compatibility with base
#'   \code{message()} so a migrated call site is a rename and nothing else; see the
#'   note on \code{.dynhr_warn()}'s \code{call.}. Both are ignored.
.dynhr_inform <- function(..., level = "info", once = FALSE, key = NULL,
                          class = NULL, verbose = NULL, appendLF = TRUE,
                          domain = NULL) {
  want <- .dynhr_levels[[match.arg(level, names(.dynhr_levels))]]
  if (.dynhr_verbosity(verbose) < want) return(invisible(NULL))
  txt <- .dynhr_fmt(...)
  if (isTRUE(once) && !.dynhr_first_time(if (is.null(key)) txt else key))
    return(invisible(NULL))
  .dynhr_signal_message(txt, class = class)
  invisible(NULL)
}

#' Emit pre-formatted console output at \code{info} level
#'
#' The migration target for the package's ad-hoc \code{cat()} calls. \code{cat()} writes
#' to stdout and cannot be silenced, levelled or captured as a condition; this
#' routes the same text through the message stream instead.
#'
#' Deliberately NOT for \code{print.*} / \code{summary.*} / \code{format.*} methods: an object
#' rendering itself is output, not messaging, and must keep going to stdout so
#' that \code{capture.output()} and sinking behave as users expect.
#'
#' \code{sep} and the trailing newline follow \code{cat()}'s conventions so a call site
#' can be converted by swapping the function name.
#'
#' @param ... Objects to paste, as \code{cat()}.
#' @param sep Separator between elements.
#' @inheritParams .dynhr_inform
#' @noRd
.dynhr_cat <- function(..., sep = " ", level = "info", verbose = NULL,
                       class = NULL) {
  want <- .dynhr_levels[[match.arg(level, names(.dynhr_levels))]]
  if (.dynhr_verbosity(verbose) < want) return(invisible(NULL))
  txt <- paste(unlist(lapply(list(...), as.character)), collapse = sep)
  ## cat() callers embed their own newlines; strip one trailing newline so the
  ## condition does not end up double-spaced once message() adds its own.
  txt <- sub("\n$", "", txt)
  .dynhr_signal_message(txt, class = class)
  invisible(NULL)
}

# ============================================================================
# Numerical-fallback handlers that do NOT swallow programming errors
# ============================================================================
## (Package-wide helper; its main consumers are the likelihood/gradient
## fallbacks in R/analytic-gradient.R, R/whittle-likelihood.R, R/posterior.R
## and the gradient-* files.)
##
## Many likelihood / gradient / sampler call sites map an error from a deep
## numerical computation to a fallback value (NULL -> FD gradient, -Inf ->
## rejected draw). That is right for a NUMERICAL failure at an extreme draw,
## but the same catch-all handler also converted genuine code bugs into
## silent fallbacks: the 0.9.3.84 debiased-Whittle and 0.9.3.85 adjoint-uni
## "[[<- NULL" index-shift bugs raised "subscript out of bounds" on every
## draw, were caught here, and surfaced only as a slower (or, for the debias
## block, a silently non-debiased) gradient.
##
## Base-R messages that can only come from a programming error. They are
## matched through gettext(domain = "R") as well, so the classification also
## works under a translated locale.
.dynhr_bug_error_templates <- c(
  "subscript out of bounds",
  "object '%s' not found",
  "could not find function \"%s\"",
  "unused argument %s",
  "unused arguments %s",
  "argument \"%s\" is missing, with no default",
  "non-conformable arguments",
  "non-conformable arrays",
  "$ operator is invalid for atomic vectors",
  "attempt to apply non-function",
  "incorrect number of dimensions",
  "incorrect number of subscripts",
  "incorrect number of subscripts on matrix",
  "non-numeric argument to binary operator",
  "invalid argument to unary operator",
  "argument is of length zero",
  "attempt to select less than one element in %s",
  "attempt to select more than one element in %s")

#' TRUE when an error condition is a programming error (not a numerical one)
#'
#' Deliberately a DENY-list: an error is treated as a bug only when its class
#' or base-R message template can never result from valid numerical input
#' (a subscript/argument/object/dimension error). Every other error -- a
#' failed Cholesky, a singular solve, a non-converged steady state, a
#' package-signalled infeasibility -- keeps its documented numerical fallback,
#' so results on valid inputs are unchanged. An unsafe-.mod-expression
#' refusal is also re-raised (it is a security refusal, never a fallback).
#' @param e A condition object.
#' @return Length-1 logical.
#' @noRd
.dynhr_is_programming_error <- function(e) {
  ## dynhr_error_theta_names: a caller passed a theta the objective cannot
  ## map (wrong length / names) -- a programming error at the call site, so
  ## it must propagate, never be swallowed into a -Inf draw (0.9.4).
  if (inherits(e, c("subscriptOutOfBoundsError", "missingArgError",
                    "dynhr_error_unsafe_mod_expression",
                    "dynhr_error_theta_names")))
    return(TRUE)
  if (inherits(e, "dynhr_error")) return(FALSE)
  msg <- conditionMessage(e)
  if (!is.character(msg) || length(msg) != 1L || is.na(msg)) return(FALSE)
  tpls <- unique(c(.dynhr_bug_error_templates,
                   gettext(.dynhr_bug_error_templates, domain = "R")))
  ## A template matches when every literal piece around its "%s" slots
  ## occurs in the message (fixed-string matching: no regex escaping).
  any(vapply(strsplit(tpls, "%s", fixed = TRUE), function(pieces)
    all(vapply(pieces, grepl, logical(1), x = msg, fixed = TRUE)),
    logical(1)))
}

#' Error-handler body for numerical fallbacks: re-raise a programming error,
#' return \code{value} otherwise.
#'
#' Use INSIDE an existing \code{error =} handler, as
#' \code{error = function(e) .dynhr_reraise_bug(e, NULL)}; mirrors
#' \code{.dynhr_reraise_unsafe()} (R/steady-monolith.R).
#' @noRd
.dynhr_reraise_bug <- function(e, value) {
  if (.dynhr_is_programming_error(e)) stop(e)
  value
}
