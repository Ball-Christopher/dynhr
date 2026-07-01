## R/checkpoint.R
## ---------------------------------------------------------------------------
## Streaming + restart support for the MCMC samplers.
##
## A checkpointed chain owns its own files in the checkpoint directory, so
## parallel chains never write to the same file (parallel-safe by construction):
##
##   chain_<id>.draws      raw little-endian float64, row-major: one P-vector per
##                         stored draw. Appended in flush-sized chunks, so the
##                         in-RAM draw footprint is bounded by flush_every * P
##                         instead of n_draws * P.
##   chain_<id>.lp         raw float64: the theta-space log-posterior trace,
##                         streamed alongside the draws.
##   chain_<id>.state.rds  the restart state (sampler position, RNG seed, scale,
##                         proposal covariance, draws-done count, accept count).
##                         Rewritten after each flush, so an interrupted run
##                         loses at most one flush window.
##   meta.rds              shared, written once at run start: sampler name and a
##                         config fingerprint. A resume that does not match it is
##                         refused (same parameters are forced).
##
## All helpers are internal (.ckpt_*); the user-facing entry points are the
## `checkpoint = ` / `resume = ` arguments of run_posterior_estimation() and the
## standalone resume_estimation() wrapper.
## ---------------------------------------------------------------------------

## File paths for one chain. `chain_id` may be NULL (single-chain runs -> "1").
.ckpt_paths <- function(dir, chain_id = NULL) {
  id <- if (is.null(chain_id)) "1" else as.character(chain_id)
  list(
    draws = file.path(dir, sprintf("chain_%s.draws", id)),
    lp    = file.path(dir, sprintf("chain_%s.lp",    id)),
    state = file.path(dir, sprintf("chain_%s.state.rds", id)),
    meta  = file.path(dir, "meta.rds"))
}

## Append an (n x P) block of draws as row-major float64.
.ckpt_append_draws <- function(path, mat) {
  con <- file(path, open = "ab")
  on.exit(close(con))
  writeBin(as.double(t(mat)), con, size = 8L, endian = "little")
  invisible(NULL)
}

## Append a length-n block of log-posterior values.
.ckpt_append_lp <- function(path, v) {
  con <- file(path, open = "ab")
  on.exit(close(con))
  writeBin(as.double(v), con, size = 8L, endian = "little")
  invisible(NULL)
}

## Read all streamed draws back as an (N x P) matrix (N inferred from file size).
.ckpt_read_draws <- function(path, n_par, par_names = NULL) {
  if (!file.exists(path) || file.size(path) == 0)
    return(matrix(numeric(0), 0L, n_par, dimnames = list(NULL, par_names)))
  n <- as.integer(file.size(path) %/% 8L)
  v <- readBin(path, "double", n = n, size = 8L, endian = "little")
  matrix(v, ncol = n_par, byrow = TRUE, dimnames = list(NULL, par_names))
}

.ckpt_read_lp <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(numeric(0))
  readBin(path, "double", n = as.integer(file.size(path) %/% 8L),
          size = 8L, endian = "little")
}

## Truncate the streamed draw/lp files back to exactly `n_done` rows. Used on
## resume to discard any partial flush written after the last saved state (so the
## draw count and the state always agree).
.ckpt_truncate <- function(paths, n_done, n_par) {
  want_draws <- as.numeric(n_done) * n_par * 8
  want_lp    <- as.numeric(n_done) * 8
  if (file.exists(paths$draws) && file.size(paths$draws) > want_draws) {
    v <- readBin(paths$draws, "raw", n = want_draws)
    writeBin(v, paths$draws)
  }
  if (file.exists(paths$lp) && file.size(paths$lp) > want_lp) {
    v <- readBin(paths$lp, "raw", n = want_lp)
    writeBin(v, paths$lp)
  }
  invisible(NULL)
}

## Config fingerprint: the invariants a resume must match.
.ckpt_fingerprint <- function(par_names, prior_spec = NULL, obs_vars = NULL,
                              extra = NULL) {
  list(par_names = par_names, prior_spec = prior_spec,
       obs_vars = obs_vars, extra = extra)
}

.ckpt_meta_write <- function(path, sampler, fingerprint) {
  saveRDS(list(sampler = sampler, fingerprint = fingerprint,
               dynhr_ckpt_version = 1L), path)
  invisible(NULL)
}

## Verify the saved meta matches the current call; stop() on any mismatch so a
## resume can never silently mix incompatible configurations.
.ckpt_meta_verify <- function(path, sampler, fingerprint) {
  if (!file.exists(path))
    stop("checkpoint resume: 'meta.rds' not found in the checkpoint directory -- ",
         "there is no run to resume.", call. = FALSE)
  m <- readRDS(path)
  if (!identical(m$sampler, sampler))
    stop("checkpoint resume: sampler mismatch (saved '", m$sampler,
         "', requested '", sampler, "'). Resume requires the same sampler.",
         call. = FALSE)
  if (!isTRUE(all.equal(m$fingerprint, fingerprint)))
    stop("checkpoint resume: the model / prior / parameter configuration differs ",
         "from the saved run. Resume forces identical parameters; start a fresh ",
         "run (new checkpoint directory) for a different configuration.",
         call. = FALSE)
  invisible(m)
}

## Save the per-chain restart state. Called after each flush; writes to a temp
## file then renames, so a crash mid-write cannot corrupt the state file.
.ckpt_save_state <- function(path, state) {
  tmp <- paste0(path, ".tmp")
  saveRDS(state, tmp)
  file.rename(tmp, path)
  invisible(NULL)
}

.ckpt_load_state <- function(path) {
  if (!file.exists(path))
    stop("checkpoint resume: per-chain state file '", basename(path),
         "' not found.", call. = FALSE)
  readRDS(path)
}
