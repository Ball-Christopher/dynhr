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

## --------------------------------------------------------------------------
## Child seeds: one base seed -> many statistically independent streams.
##
## `set.seed(base + k)` is NOT a family of independent streams under R's
## Mersenne-Twister: R initialises the 624-word state from a linear
## congruential recurrence in the seed, so streams from CONSECUTIVE (or evenly
## spaced) integer seeds are dependent in their first draws. Over 2e5
## consecutive seeds the across-seed lag-1 correlation of the FIRST rnorm() is
## -0.05 (z = -22.8), of the first runif() z = -26; a stride of 1000 still
## gives z = -5.5. Draws WITHIN one stream are fine. Every per-chain /
## per-replication / per-task seed in dynhr is therefore derived here instead.
##
## The derivation is counter-based: child k of (seed, stream) is a keyed hash
## of k, so it is O(1) per child, PREFIX-STABLE (child k is the same whatever
## n is -- adding chains, resuming, or growing an SBC run keeps the existing
## children), never touches the RNG (no global side effect, independent of
## RNGkind()), and identical on every platform (exact integer arithmetic in
## doubles). Children of one (seed, stream) are pairwise DISTINCT for
## k < 2^32 (a bijection of k); `stream` separates call sites that share a
## base seed. A child is itself a valid base seed, so (stage, task) indices
## nest: .child_seeds(.child_seed_at(seed, stage, "a"), n, "b").
## --------------------------------------------------------------------------

## Version of the seed derivation, recorded in checkpoint metadata. 1 = the
## old `seed + k` arithmetic (absent field), 2 = .child_seed_at().
.child_seed_scheme <- 2L

.U32 <- 4294967296   # 2^32

## Unsigned 32-bit XOR / multiply on doubles in [0, 2^32): split into 16-bit
## halves so every intermediate stays an exact double (< 2^49) and bitwXor()
## only ever sees non-negative integers < 2^16.
.u32_xor <- function(a, b) {
  bitwXor(a %/% 65536, b %/% 65536) * 65536 + bitwXor(a %% 65536, b %% 65536)
}
.u32_mul <- function(a, c) {
  ((a %/% 65536 * c) %% 65536 * 65536 + (a %% 65536) * c) %% .U32
}
## MurmurHash3 32-bit finaliser: a bijection of [0, 2^32) with full avalanche.
.u32_fmix <- function(h) {
  h <- .u32_xor(h, h %/% 65536)
  h <- .u32_mul(h, 2246822507)            # 0x85ebca6b
  h <- .u32_xor(h, h %/% 8192)
  h <- .u32_mul(h, 3266489909)            # 0xc2b2ae35
  .u32_xor(h, h %/% 65536)
}
## FNV-1a hash of a stream label (UTF-8 bytes).
.u32_str_hash <- function(s) {
  h <- 2166136261
  for (b in as.integer(charToRaw(enc2utf8(s))))
    h <- .u32_mul(.u32_xor(h, b), 16777619)
  h
}

#' Child seed(s) of a base seed at given indices
#'
#' @param seed Base seed: a finite integer-valued scalar in the
#'   `set.seed()` range, or `NULL`. `NULL` draws ONE base from the caller's
#'   current RNG stream (`sample.int(.Machine$integer.max, 1L)`, advancing it as
#'   any unseeded stochastic call does), so the children are reproducible
#'   under an outer `set.seed()`; it is the only case that touches the RNG.
#' @param k Non-negative integer-valued index vector (any length, < 2^53).
#' @param stream Single string naming the call site; different streams give
#'   unrelated children from the same base seed.
#' @return Integer vector, `length(k)`, valid `set.seed()` seeds.
#' @noRd
.child_seed_at <- function(seed, k, stream) {
  if (!is.character(stream) || length(stream) != 1L || is.na(stream) ||
      !nzchar(stream))
    .dynhr_abort("`stream` must be a single non-empty string.",
                 class = "dynhr_error_invalid_argument")
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != trunc(seed) || abs(seed) > .Machine$integer.max)
    .dynhr_abort("`seed` must be NULL or a finite integer-valued scalar with ",
                 "|seed| <= .Machine$integer.max (the set.seed() range).",
                 class = "dynhr_error_invalid_argument")
  if (!is.numeric(k) || anyNA(k) || any(k < 0) || any(k != trunc(k)) ||
      any(k >= 2^53))
    .dynhr_abort("`k` must be non-negative integer-valued indices.",
                 class = "dynhr_error_invalid_argument")
  base <- as.numeric(seed) %% .U32                     # two's complement
  key1 <- .u32_fmix(.u32_xor(.u32_fmix(.u32_xor(base, 2654435769)),
                             .u32_str_hash(stream)))
  key2 <- .u32_fmix(.u32_xor(key1, 1779033703))
  ## Weyl sequence in k (odd multiplier => a bijection of k mod 2^32), then two
  ## keyed finaliser rounds.
  x   <- (key1 + .u32_mul(as.numeric(k) %% .U32, 2654435769)) %% .U32
  out <- .u32_fmix(.u32_xor(.u32_fmix(x), key2))
  ## To a signed 32-bit seed (set.seed() uses all 32 bits); -2^31 is NA_integer_
  ## in R, so its single preimage is folded onto 0.
  out <- ifelse(out >= 2147483648, out - .U32, out)
  out[out == -2147483648] <- 0
  as.integer(out)
}

#' The first n child seeds of a base seed (prefix-stable)
#'
#' `.child_seeds(seed, n, stream)[k]` is the same for every `n >= k`.
#' @inheritParams .child_seed_at
#' @param n Number of children (non-negative integer).
#' @return Integer vector of length `n`.
#' @noRd
.child_seeds <- function(seed, n, stream) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 || n != trunc(n))
    .dynhr_abort("`n` must be a non-negative integer scalar.",
                 class = "dynhr_error_invalid_argument")
  .child_seed_at(seed, seq_len(n), stream)
}
