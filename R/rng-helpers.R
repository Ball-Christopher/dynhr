## R/rng-helpers.R
## --------------------------------------------------------------------------
## RNG hygiene: seed locally, leave the caller's global stream untouched.
##
## Many exported functions take a `seed` (often with a constant default, which
## is the reproducibility contract: the same call gives the same answer). A bare
## `set.seed(seed)` inside them also RESETS the user's global stream, so the
## user's next `runif()` after the call depends on dynhr's seed rather than on
## the user's own history. `.local_seed()` keeps the reproducibility contract
## (the function body sees exactly the stream `set.seed(seed)` gives) and
## restores the caller's `.Random.seed` -- or its ABSENCE -- when the calling
## function exits.
##
## The expression-scoped variant, `.with_local_seed(seed, expr)`, lives in
## R/tpf-likelihood.R.
## --------------------------------------------------------------------------

#' Seed the RNG for the rest of the calling function, restoring on its exit
#'
#' Saves the global `.Random.seed` (or records that there is none), calls
#' `set.seed(seed)`, and registers an exit handler on `envir` (by default the
#' caller's frame) that puts the saved state back -- or removes `.Random.seed`
#' again when it did not exist. `seed = NULL` is a no-op (no seeding, no
#' restore), matching the `if (!is.null(seed)) set.seed(seed)` idiom it
#' replaces.
#'
#' The handler is registered with `after = FALSE`, so when a frame calls
#' `.local_seed()` more than once (e.g. once per stage) the FIRST call's saved
#' state -- the caller's original stream -- is the one restored last, i.e. the
#' one that survives.
#'
#' Scope matters: the stream reverts when `envir`'s function returns, so call
#' it in the frame that owns ALL the seeded draws (normally the exported
#' function), not in a helper whose caller keeps drawing afterwards.
#'
#' @param seed Seed passed to `set.seed()`, or `NULL`.
#' @param envir Frame whose exit restores the stream.
#' @return `invisible(NULL)`.
#' @noRd
.local_seed <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(NULL))
  ge  <- globalenv()
  had <- exists(".Random.seed", envir = ge, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = ge, inherits = FALSE) else NULL
  restore <- function() {
    if (had) {
      assign(".Random.seed", old, envir = ge)
    } else if (exists(".Random.seed", envir = ge, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = ge)
    }
  }
  ## A call whose function slot is the closure itself, so it evaluates in
  ## `envir` without `restore` having to be visible there.
  do.call(base::on.exit, list(as.call(list(restore)), add = TRUE, after = FALSE),
          envir = envir)
  set.seed(seed)
  invisible(NULL)
}
