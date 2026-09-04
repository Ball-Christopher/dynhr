## R/chain-state.R
## ---------------------------------------------------------------------------
## Sampler-agnostic resumable-chain state packs.
##
## The internal checkpoint machinery (R/checkpoint.R) already gives every
## built-in sampler an exact, bit-identical restart. This file exposes the
## same idea as a small PUBLIC API, so bespoke user samplers (e.g. a
## two-block Metropolis-within-Gibbs written outside the package) can save
## and resume chains under the same contract:
##
##   Restoring a state saved after sweep n and continuing for m more sweeps
##   must be bit-identical to having run n + m sweeps in one process,
##   provided adaptation had ended by sweep n. The two ingredients a naive
##   "restart from the last draw" loses are captured here explicitly:
##   `.Random.seed` (the exact RNG stream position) and the adaptation
##   state (proposal scales / covariance accumulators / frozen flag).
##
## On disk a pack is a single .rds holding the SERIALIZED state as a raw
## payload plus its md5 checksum, so a corrupted or tampered file (including
## a tampered RNG state) is refused at restore rather than silently
## producing a wrong-but-plausible chain.
## ---------------------------------------------------------------------------

## Current on-disk pack format version.
.chain_state_version <- 1L

## Snapshot the global RNG state, initialising the generator first if this
## session has never drawn (a fresh R process has no .Random.seed yet).
.rng_snapshot <- function() {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    stats::runif(1L)
  get(".Random.seed", envir = .GlobalEnv)
}

## md5 of a raw vector (base R only: tools::md5sum works on files).
.raw_md5 <- function(raw_vec) {
  tmp <- tempfile("dynhr_md5_")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(raw_vec, tmp)
  unname(tools::md5sum(tmp))
}


#' Capture a resumable MCMC chain state
#'
#' Builds a sampler-agnostic state pack: everything a Metropolis-type
#' sampler needs to continue a chain \emph{exactly} where it stopped. Use
#' with \code{\link{mcmc_chain_save}} / \code{\link{mcmc_chain_restore}}
#' from inside your own sampler loop; the built-in samplers get the same
#' behaviour through \code{\link{mcmc}}'s \code{checkpoint_dir} argument.
#'
#' @section The resume contract:
#' Restoring a state saved after sweep \eqn{n} and running \eqn{m} more
#' sweeps must be bit-identical to running \eqn{n + m} sweeps in one
#' process. That holds iff (a) the RNG state is restored (done for you by
#' \code{\link{mcmc_chain_restore}}), (b) the sampler resumes from
#' \code{position}/\code{lp} without re-evaluating any random step, and
#' (c) adaptation had already ended at sweep \eqn{n} \emph{or} the full
#' adapter state is carried in \code{adapt_state} and the sampler resumes
#' it exactly. If adaptation was still running and you did not capture its
#' full state, set \code{adapt_frozen = FALSE}: restore will then refuse
#' by default rather than let a subtly different chain masquerade as a
#' continuation.
#'
#' @param position    Named numeric vector: the current chain position
#'   (the last accepted draw).
#' @param lp          Log-posterior at \code{position} (avoids one
#'   re-evaluation on resume; also a cheap consistency check).
#' @param sweep       Integer: number of sweeps completed when the state
#'   was captured.
#' @param scales      Proposal scale(s) -- scalar, vector, or list; stored
#'   as-is.
#' @param adapt_state Any further adaptation state (covariance
#'   accumulators, dual-averaging accumulators, ...); stored as-is.
#' @param adapt_frozen Logical: has adaptation ended (\code{TRUE},
#'   default), or is \code{adapt_state} a mid-adaptation snapshot that the
#'   resuming sampler must continue exactly (\code{FALSE})?
#' @param rng         RNG state to store. Default: the current
#'   \code{.Random.seed} (initialising the generator if needed). Capture
#'   this \emph{after} the last proposal/accept draw of sweep \code{sweep}.
#' @param meta        Free-form metadata (chain id, arm tag, model
#'   fingerprint, ...); compared by the caller, not by the package.
#' @return An object of class \code{"dynhr_chain_state"}.
#' @seealso \code{\link{mcmc_chain_save}}, \code{\link{mcmc_chain_restore}},
#'   \code{\link{mcmc_chain_extend}}, \code{\link{mcmc}}
#' @export
mcmc_chain_state <- function(position, lp = NA_real_, sweep = 0L,
                             scales = NULL, adapt_state = NULL,
                             adapt_frozen = TRUE,
                             rng = NULL, meta = NULL) {
  if (!is.numeric(position))
    stop("mcmc_chain_state: 'position' must be a numeric vector.", call. = FALSE)
  if (!is.logical(adapt_frozen) || length(adapt_frozen) != 1L || is.na(adapt_frozen))
    stop("mcmc_chain_state: 'adapt_frozen' must be TRUE or FALSE.", call. = FALSE)
  if (is.null(rng)) rng <- .rng_snapshot()
  st <- list(position = position, lp = lp, sweep = as.integer(sweep),
             scales = scales, adapt_state = adapt_state,
             adapt_frozen = adapt_frozen, rng = rng, meta = meta)
  class(st) <- "dynhr_chain_state"
  st
}


#' Save a chain state pack to disk
#'
#' Writes atomically (temp file + rename), so a crash mid-write can never
#' leave a truncated state file behind. The state is stored as a
#' serialized raw payload together with its md5 checksum;
#' \code{\link{mcmc_chain_restore}} refuses the file if the two disagree.
#'
#' @param state A \code{\link{mcmc_chain_state}} object.
#' @param file  Path to write (conventionally \code{*.state.rds}).
#' @return \code{file}, invisibly.
#' @seealso \code{\link{mcmc_chain_state}}, \code{\link{mcmc_chain_restore}},
#'   \code{\link{mcmc_chain_extend}}, \code{\link{mcmc}} (whose
#'   \code{checkpoint_dir} / \code{resume} arguments use this contract)
#' @examples
#' state <- mcmc_chain_state(position = c(a = 0.1, b = -0.4),
#'                           lp = -12.3, sweep = 100L)
#' f <- file.path(tempdir(), "chain.state.rds")
#' mcmc_chain_save(state, f)
#'
#' back <- mcmc_chain_restore(f)
#' identical(back$position, state$position)
#' back$sweep
#' @export
mcmc_chain_save <- function(state, file) {
  if (!inherits(state, "dynhr_chain_state"))
    stop("mcmc_chain_save: 'state' must come from mcmc_chain_state().",
         call. = FALSE)
  payload <- serialize(state, NULL)
  pack <- list(format = "dynhr_chain_state",
               version = .chain_state_version,
               md5 = .raw_md5(payload),
               payload = payload)
  tmp <- paste0(file, ".tmp")
  saveRDS(pack, tmp)
  file.rename(tmp, file)
  invisible(file)
}


#' Restore a chain state pack
#'
#' Reads a pack written by \code{\link{mcmc_chain_save}}, verifies its
#' checksum (a corrupted or tampered file -- including a tampered RNG
#' state -- is an error, never a silent misresume), and by default
#' restores \code{.Random.seed}, so the very next \code{rnorm()} /
#' \code{runif()} call in the resuming sampler continues the exact random
#' stream of the saved run.
#'
#' @param file  Path written by \code{\link{mcmc_chain_save}}.
#' @param allow_mid_adaptation Logical (default \code{FALSE}): a state
#'   saved with \code{adapt_frozen = FALSE} is refused unless this is
#'   \code{TRUE} \emph{and} your sampler genuinely resumes the stored
#'   adapter state exactly (see the resume contract in
#'   \code{\link{mcmc_chain_state}}).
#' @param restore_rng Logical (default \code{TRUE}): assign the stored RNG
#'   state to \code{.Random.seed}. Set \code{FALSE} only for inspection.
#' @return The \code{\link{mcmc_chain_state}} object.
#' @seealso \code{\link{mcmc_chain_state}}, \code{\link{mcmc_chain_save}},
#'   \code{\link{mcmc_chain_extend}}
#' @examples
#' set.seed(42)
#' state <- mcmc_chain_state(position = c(a = 0.1, b = -0.4), sweep = 100L)
#' f <- file.path(tempdir(), "chain.state.rds")
#' mcmc_chain_save(state, f)
#'
#' ## restore_rng = TRUE (the default) rewinds .Random.seed to the saved
#' ## stream, so the resuming sampler draws exactly what the original would.
#' x1 <- rnorm(3)
#' back <- mcmc_chain_restore(f)
#' identical(rnorm(3), x1)
#'
#' ## restore_rng = FALSE inspects the pack without touching the RNG
#' peek <- mcmc_chain_restore(f, restore_rng = FALSE)
#' peek$position
#' @export
mcmc_chain_restore <- function(file, allow_mid_adaptation = FALSE,
                               restore_rng = TRUE) {
  if (!file.exists(file))
    stop("mcmc_chain_restore: file '", file, "' not found.", call. = FALSE)
  pack <- readRDS(file)
  if (!is.list(pack) || !identical(pack$format, "dynhr_chain_state"))
    stop("mcmc_chain_restore: '", file,
         "' is not a dynhr chain state pack.", call. = FALSE)
  if (!identical(pack$version, .chain_state_version))
    stop("mcmc_chain_restore: pack version ", pack$version,
         " != supported version ", .chain_state_version, ".", call. = FALSE)
  if (!identical(.raw_md5(pack$payload), pack$md5))
    stop("mcmc_chain_restore: checksum mismatch -- the state file is ",
         "corrupted or has been tampered with. Refusing to resume.",
         call. = FALSE)
  st <- unserialize(pack$payload)
  if (!inherits(st, "dynhr_chain_state"))
    stop("mcmc_chain_restore: payload is not a dynhr_chain_state object.",
         call. = FALSE)
  if (!isTRUE(st$adapt_frozen) && !isTRUE(allow_mid_adaptation))
    stop("mcmc_chain_restore: this state was saved MID-ADAPTATION ",
         "(adapt_frozen = FALSE). Resuming it does not satisfy the ",
         "bit-identical contract unless the sampler continues the stored ",
         "adapter state exactly; pass allow_mid_adaptation = TRUE only if ",
         "yours does.", call. = FALSE)
  if (isTRUE(restore_rng) && !is.null(st$rng))
    assign(".Random.seed", st$rng, envir = .GlobalEnv)
  st
}


#' Append new draws to a saved chain .rds
#'
#' Concatenation helper for downstream diagnostics that glob per-chain
#' \code{*.rds} files: appends freshly drawn sweeps to the draw matrix in
#' an existing checkpoint list, re-stamps its sweep counter / done flag,
#' and rewrites the file atomically. The file must hold a plain
#' \code{list()} with the draw matrix under \code{draws_field}; all other
#' fields are preserved untouched.
#'
#' @param file        Path to the existing chain .rds.
#' @param new_draws   Matrix (or vector for one sweep) of new draws;
#'   column count must match the stored draw matrix.
#' @param new_lp      Optional numeric vector of log-posterior values,
#'   one per new sweep, appended to \code{lp_field}.
#' @param draws_field Name of the draw-matrix field (default \code{"z"}).
#' @param lp_field    Name of the log-posterior field (default \code{"lp"}).
#' @param sweep_field Name of the sweep-counter field (default
#'   \code{"sweep"}); re-stamped to the new total row count when present
#'   in the list.
#' @param done        Value to stamp into the list's \code{done} field
#'   (default \code{TRUE}); set \code{NULL} to leave it untouched.
#' @return The updated list, invisibly (the file is rewritten in place).
#' @export
mcmc_chain_extend <- function(file, new_draws, new_lp = NULL,
                              draws_field = "z", lp_field = "lp",
                              sweep_field = "sweep", done = TRUE) {
  if (!file.exists(file))
    stop("mcmc_chain_extend: file '", file, "' not found.", call. = FALSE)
  obj <- readRDS(file)
  if (!is.list(obj) || is.null(obj[[draws_field]]))
    stop("mcmc_chain_extend: '", file, "' has no '", draws_field,
         "' field.", call. = FALSE)
  old <- obj[[draws_field]]
  if (!is.matrix(old)) old <- matrix(old, nrow = NROW(old))
  if (!is.matrix(new_draws)) new_draws <- matrix(new_draws, ncol = ncol(old))
  if (ncol(new_draws) != ncol(old))
    stop("mcmc_chain_extend: new draws have ", ncol(new_draws),
         " columns; the stored chain has ", ncol(old), ".", call. = FALSE)
  obj[[draws_field]] <- rbind(old, new_draws)
  if (!is.null(new_lp)) {
    if (length(new_lp) != nrow(new_draws))
      stop("mcmc_chain_extend: length(new_lp) must equal nrow(new_draws).",
           call. = FALSE)
    obj[[lp_field]] <- c(obj[[lp_field]], as.numeric(new_lp))
  }
  if (sweep_field %in% names(obj))
    obj[[sweep_field]] <- nrow(obj[[draws_field]])
  if (!is.null(done)) obj$done <- done
  tmp <- paste0(file, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, file)
  invisible(obj)
}
