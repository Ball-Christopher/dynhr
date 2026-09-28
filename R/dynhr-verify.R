## R/dynhr-verify.R
## ---------------------------------------------------------------------------
## dynhr_verify(): is a finished estimation result reproduced by the CURRENT
## build and numerical environment? (E5 C3, plan section 7.)
##
## The recorded spec (the result's run record) is re-evaluated here:
##   * deterministic quantities -- the log-posterior at the recorded mode and
##     at a sample of recorded draws, and the mode itself through a short
##     re-optimisation started AT the recorded mode -- to `tol` relative;
##   * a short distributional replay of the sampler(s), whose posterior means
##     must agree with the recorded ones within Monte Carlo error.
## The verdict is reported per class: A target/algorithm (the spec's content
## hashes recomputed), B code (version / commit, and the registered result
## changes between the two builds that touch the run), C numerical
## environment (R, platform, BLAS/LAPACK).
## ---------------------------------------------------------------------------

#' Verify an estimation result under the current build and environment
#'
#' Re-evaluates what a finished estimation result claims, using the current
#' dynhr build and numerical environment, and reports PASS or FAIL for each
#' reproducibility class of its run record. It answers "I am on a different
#' OS / BLAS / dynhr build: are my results the same?".
#'
#' @section What is checked:
#' The run is rebuilt from the result's run record (its estimation spec;
#' schema-1 records are converted with \code{\link{as_estimation_spec}}).
#' \describe{
#'   \item{\code{logpost_at_mode}}{The log-posterior, rebuilt now, at the
#'     recorded mode against the recorded mode log-posterior.}
#'   \item{\code{logpost_at_draws}}{The log-posterior at \code{n_points}
#'     recorded draws (evenly spaced over the pooled draws) against the
#'     values the sampler recorded with them. Skipped for noisy particle
#'     likelihoods and for samplers that keep no per-draw log-posterior.}
#'   \item{\code{mode_reoptimisation}}{A short re-optimisation
#'     (\code{mode_iter} iterations of the recorded optimiser, serial) started
#'     at the recorded mode: the recorded mode must still be a mode, i.e. the
#'     re-optimised log-posterior may exceed the recorded one by at most
#'     \code{tol} (relative). The largest parameter move is reported.}
#'   \item{\code{posterior_means}}{(with \code{replay = TRUE} and recorded
#'     draws) the sampler stage is replayed from the re-optimised mode with
#'     \code{replay_draws} draws per sampler (samplers without a draw count,
#'     such as SMC, keep theirs), and each posterior mean must agree with the
#'     recorded one within Monte Carlo error: \eqn{|z| \le} the two-sided
#'     Bonferroni bound \code{qnorm(1 - alpha / (2 p))} over the \eqn{p}
#'     parameters, with \eqn{z} the mean difference over
#'     \eqn{\sqrt{MCSE_1^2 + MCSE_2^2}} and each MCSE the standard deviation
#'     over the square root of the split-chain effective sample size.}
#' }
#' The deterministic checks use the relative difference
#' \eqn{|a - b| / \max(|b|, 1)}, so a log-posterior near zero is compared
#' absolutely.
#'
#' @section Classes and verdict:
#' \describe{
#'   \item{A (target/algorithm)}{The spec's content hashes of the model, data
#'     (a data file is re-read), likelihood, sampler and option snapshot are
#'     recomputed: FAIL when any differs from the recorded hash (the inputs
#'     changed since the run).}
#'   \item{B (code)}{The recorded dynhr version and \code{GIT_COMMIT} against
#'     the loaded ones. FAIL when a registered result change between the two
#'     versions touches a component the run uses (the results are then
#'     expected to differ; the changes are listed), or when the build differs
#'     and a check fails.}
#'   \item{C (numerical environment)}{R version, platform, OS, BLAS and
#'     LAPACK. FAIL when the environment differs and a check fails.}
#' }
#' The overall verdict is PASS when class A passes, no registered change
#' touches the run and every check that ran passed. With the same A, B and C
#' the deterministic checks agree to the last bit; across builds without a
#' registered change, or across environments, they agree to about 1e-8
#' relative and the draws are equivalent in distribution only (see
#' \code{\link{dynhr_rerun}}, section Reproducibility).
#'
#' @param x An estimation result carrying \code{$run_record} (from
#'   \code{\link{run_estimation}}, \code{\link{run_full_estimation}},
#'   \code{\link{run_posterior_estimation}} or \code{\link{run_mode_finding}}).
#' @param n_points Number of recorded draws at which the log-posterior is
#'   re-evaluated (default 20).
#' @param replay Run the distributional replay (default \code{TRUE}).
#' @param replay_draws Post-warmup draws per sampler in the replay (default
#'   500; never more than recorded).
#' @param mode_iter Iteration budget of the re-optimisation from the recorded
#'   mode (default 50).
#' @param tol Relative tolerance of the deterministic checks (default 1e-8).
#' @param alpha Family-wise level of the posterior-mean check (default 0.01).
#' @return A \code{dynhr_verify} object: \code{verdict} (\code{"PASS"} or
#'   \code{"FAIL"}), \code{classes} (a data frame: class, what, status
#'   \code{"same"} / \code{"differs"} / \code{"not recorded"}, verdict,
#'   details), \code{checks} (a data frame: check, value, threshold, pass,
#'   note), \code{expectation}, \code{registered_changes},
#'   \code{differences} (\code{code}, \code{environment}), \code{mode} (the
#'   recorded and re-optimised modes and log-posteriors) and \code{replay}
#'   (the posterior means, MCSEs and z-scores, or \code{NULL}).
#' @seealso \code{\link{dynhr_rerun}}, \code{\link{run_estimation}}
#' @examples
#' \donttest{
#' solved   <- solve_model(system.file("extdata/models/nk_demo.mod",
#'                                     package = "dynhr"), verbose = FALSE)
#' obs_vars <- c("ygr", "infl", "intr")
#' Y <- as.matrix(read.csv(system.file("extdata/models/nk_demo_data.csv",
#'                                     package = "dynhr"))[, obs_vars])
#' mode <- run_mode_finding(solved, Y, obs_vars = obs_vars,
#'                          n_iter = 200L, verbose = FALSE)
#' dynhr_verify(mode)
#' }
#' @export
dynhr_verify <- function(x, n_points = 20L, replay = TRUE, replay_draws = 500L,
                         mode_iter = 50L, tol = 1e-8, alpha = 0.01) {
  rec <- if (is.list(x)) x$run_record
  if (!inherits(rec, "dynhr_run_record"))
    .dynhr_abort("dynhr_verify: `x` must be an estimation result carrying ",
                 "$run_record (from run_estimation(), run_full_estimation(), ",
                 "run_posterior_estimation() or run_mode_finding()).",
                 class = "dynhr_error_no_run_record")
  counts <- list(n_points = n_points, replay_draws = replay_draws,
                 mode_iter = mode_iter)
  for (nm in names(counts)) {
    v <- counts[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1)
      .dynhr_abort("dynhr_verify: `", nm, "` must be a positive number.",
                   class = "dynhr_error_bad_argument")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0)
    .dynhr_abort("dynhr_verify: `tol` must be a positive number.",
                 class = "dynhr_error_bad_argument")
  if (!is.numeric(alpha) || length(alpha) != 1L || !(alpha > 0 && alpha < 1))
    .dynhr_abort("dynhr_verify: `alpha` must be in (0, 1).",
                 class = "dynhr_error_bad_argument")
  ## the caller's RNG stream is restored on exit (the stages re-seed)
  .local_seed(1L)

  spec <- rec$spec %||% rec$args$spec
  schema1 <- is.null(spec)
  if (schema1) spec <- suppressWarnings(suppressMessages(as_estimation_spec(rec)))
  recd <- .verify_recorded(x)

  ## ---- classes A, B, C ----------------------------------------------------
  hash_diff <- if (schema1) NULL else .verify_hash_diffs(spec)
  was <- rec$provenance
  now <- .dynhr_rr_provenance()
  d   <- .dynhr_rr_build_diffs(was, now)
  reg <- if (!identical(as.character(was$version), as.character(now$version)))
    .dynhr_result_changes_between(was$version, now$version,
                                  .est_component_tags(spec))
  else character(0)

  ## ---- deterministic: rebuild the posterior, re-optimise from the mode ----
  base <- spec
  md <- base$mode
  md["result"] <- list(NULL)
  base["mode"] <- list(md)
  quiet_out <- list(form = "mode", save = FALSE, diagnostics = FALSE,
                    stoch_simul = FALSE, smoother = FALSE,
                    obc_ppf_reweight = FALSE, ramsey = FALSE)
  quiet_cmp <- list(verbose = FALSE, checkpoint_dir = NULL, resume = FALSE)
  vs <- suppressWarnings(update(base,
    mode = list(theta_init = recd$theta_mode, n_iter = as.integer(mode_iter)),
    sampler = FALSE, compute = c(quiet_cmp, list(parallel = FALSE)),
    outputs = quiet_out))
  mv <- suppressMessages(.run_estimation_impl(vs, rr = NULL))
  lp_fn <- mv$log_post_fn
  lp_at <- function(th) {
    v <- lp_fn(th)
    if (is.list(v)) v$logpost else as.numeric(v)
  }
  rel <- function(a, b) abs(a - b) / max(abs(b), 1)

  lp_now_mode <- lp_at(recd$theta_mode)
  r_mode <- rel(lp_now_mode, recd$logpost_mode)
  checks <- .verify_row("logpost_at_mode", r_mode, tol, isTRUE(r_mode <= tol),
    sprintf("recorded %.12g, now %.12g", recd$logpost_mode, lp_now_mode))

  eff <- if (.spec_obc_active(spec)) spec$likelihood$obc_filter else spec$likelihood$type
  noisy <- eff %in% c("tpf", "ppf", "copf", "sv_rbpf", "global_pf")
  row_d <- if (noisy) {
    .verify_row("logpost_at_draws", NA_real_, tol, NA,
      paste0("skipped: the \"", eff, "\" likelihood is a noisy particle estimate"))
  } else if (is.null(recd$draws_lp)) {
    .verify_row("logpost_at_draws", NA_real_, tol, NA,
      if (is.null(recd$draws)) "skipped: the result has no draws"
      else "skipped: the sampler kept no per-draw log-posterior")
  } else {
    n <- nrow(recd$draws_lp$theta)
    idx <- unique(round(seq(1, n, length.out = min(as.integer(n_points), n))))
    r_d <- vapply(idx, function(i)
      rel(lp_at(recd$draws_lp$theta[i, ]), recd$draws_lp$lp[[i]]), numeric(1))
    .verify_row("logpost_at_draws", max(r_d), tol, isTRUE(all(r_d <= tol)),
                sprintf("%d recorded draws", length(idx)))
  }
  checks <- rbind(checks, row_d)

  lp_reopt <- mv$mode$logpost %||% lp_at(mv$theta_mode)
  gain <- (lp_reopt - recd$logpost_mode) / max(abs(recd$logpost_mode), 1)
  th_move <- max(abs(mv$theta_mode[names(recd$theta_mode)] - recd$theta_mode) /
                   pmax(abs(recd$theta_mode), 1))
  checks <- rbind(checks, .verify_row("mode_reoptimisation", gain, tol,
    isTRUE(gain <= tol),
    sprintf("%d-iteration %s restart at the recorded mode; largest parameter move %.3g (relative)",
            as.integer(mode_iter), vs$mode$method, th_move)))

  ## ---- distributional replay ------------------------------------------------
  rp <- NULL
  row_r <- if (!isTRUE(replay)) {
    .verify_row("posterior_means", NA_real_, NA_real_, NA, "skipped: replay = FALSE")
  } else if (is.null(recd$draws) || is.null(spec$sampler)) {
    .verify_row("posterior_means", NA_real_, NA_real_, NA,
                "skipped: the result has no draws")
  } else {
    rp <- .verify_replay(base, mv, recd$draws, replay_draws, quiet_cmp,
                         quiet_out, alpha)
    .verify_row("posterior_means", rp$max_abs_z, rp$z_crit,
                isTRUE(rp$max_abs_z <= rp$z_crit),
                sprintf("%d recorded vs %d replayed draws; max |z| at %s",
                        nrow(recd$draws), rp$n_replay, rp$worst))
  }
  checks <- rbind(checks, row_r)
  rownames(checks) <- NULL
  ran <- !is.na(checks$pass)
  checks_ok <- all(checks$pass[ran])

  ## ---- verdicts -------------------------------------------------------------
  a_status <- if (schema1) "not recorded" else if (length(hash_diff)) "differs" else "same"
  a_ok <- !length(hash_diff)
  b_status <- if (length(d$code)) "differs" else "same"
  c_status <- if (length(d$environment)) "differs" else "same"
  b_ok <- !length(reg) && (identical(b_status, "same") || checks_ok)
  c_ok <- identical(c_status, "same") || checks_ok
  classes <- data.frame(
    class   = c("A", "B", "C"),
    what    = c("target/algorithm", "code", "numerical environment"),
    status  = c(a_status, b_status, c_status),
    verdict = ifelse(c(a_ok, b_ok, c_ok), "PASS", "FAIL"),
    details = c(
      if (schema1) "schema-1 record: no stored spec hashes (spec rebuilt from the arguments)"
      else if (length(hash_diff)) paste0("hash differs: ", paste(hash_diff, collapse = ", "))
      else "spec hashes reproduce",
      paste(c(d$code, if (length(reg)) paste0(length(reg),
              " registered result change(s) touch this run")), collapse = "; "),
      paste(d$environment, collapse = "; ")),
    stringsAsFactors = FALSE)
  verdict <- if (a_ok && !length(reg) && checks_ok) "PASS" else "FAIL"
  expectation <- if (length(reg))
    paste0("results are expected to differ: registered result changes touch ",
           "this run; install dynhr ", format(was$version), " to reproduce it")
  else if (length(d$code) || length(d$environment))
    "deterministic quantities agree to ~1e-8 relative; draws are equivalent in distribution"
  else "bit-identical"
  if (identical(verdict, "FAIL") && a_ok && !length(reg) &&
      !length(d$code) && !length(d$environment))
    expectation <- paste0(expectation, " (but a check failed, and no recorded ",
                          "difference explains it)")

  structure(list(
    verdict = verdict, classes = classes, checks = checks,
    expectation = expectation, registered_changes = reg,
    differences = d,
    mode = list(theta_recorded = recd$theta_mode,
                theta_reoptimised = mv$theta_mode,
                logpost_recorded = recd$logpost_mode,
                logpost_now = lp_now_mode,
                logpost_reoptimised = lp_reopt),
    replay = rp), class = "dynhr_verify")
}

## One row of the checks table.
.verify_row <- function(check, value, threshold, pass, note)
  data.frame(check = check, value = value, threshold = threshold, pass = pass,
             note = note, stringsAsFactors = FALSE)

## What a result recorded: the mode (theta, log-posterior), the pooled draws,
## and the draws that carry their own log-posterior (list(theta, lp)).
.verify_recorded <- function(x) {
  mr <- if (inherits(x, "dynhr_posterior_result")) x$mode_result else x
  if (inherits(x, "dynhr_estimation_result")) {
    theta <- x$mode$theta_mode
    lp    <- x$mode$logpost
  } else {
    theta <- mr$theta_mode
    lp    <- mr$mode$logpost %||% mr$mode$value
  }
  if (is.null(theta) || !is.numeric(lp) || length(lp) != 1L)
    .dynhr_abort("dynhr_verify: the result carries no mode (theta_mode and ",
                 "its log-posterior).", class = "dynhr_error_bad_argument")
  runs <- list()
  draws <- NULL
  if (inherits(x, "dynhr_posterior_result")) {
    draws <- x$pooled_draws
    for (m in names(x$chains)) runs <- c(runs, x$chains[[m]]$chains)
  } else if (inherits(x, "dynhr_estimation_result") && !is.null(x$chains)) {
    draws <- x$chains$chain
    runs <- if (!is.null(x$chains$chain_list)) x$chains$chain_list else list(x$chains)
  }
  th <- list(); lps <- list()
  for (r in runs) {
    if (is.list(r) && is.matrix(r$chain) && is.numeric(r$post_logpost) &&
        length(r$post_logpost) == nrow(r$chain)) {
      th  <- c(th, list(r$chain))
      lps <- c(lps, list(r$post_logpost))
    }
  }
  draws_lp <- if (length(th)) {
    ok <- is.finite(unlist(lps))
    tm <- do.call(rbind, th)
    if (any(ok)) list(theta = tm[ok, , drop = FALSE], lp = unlist(lps)[ok])
  }
  list(theta_mode = theta, logpost_mode = lp,
       draws = if (is.matrix(draws) && nrow(draws)) draws, draws_lp = draws_lp)
}

## Target parts whose recomputed content hash differs from the recorded one.
.verify_hash_diffs <- function(spec) {
  h <- .spec_hashes(.spec_parts(spec))
  parts <- c("model", "data", "likelihood", "sampler", "options")
  if (is.null(spec$mode$result)) parts <- c(parts, "mode")
  parts[!vapply(parts, function(p) identical(h[[p]], spec$hashes[[p]]), logical(1))]
}

## Posterior mean and MCSE (sd / sqrt(split-chain ESS)) of each column.
.verify_mean_mcse <- function(draws) {
  m <- colMeans(draws)
  se <- vapply(seq_len(ncol(draws)), function(j) {
    x <- draws[, j]
    ess <- if (length(x) >= 4L) .d5_ess_basic(.d5_split(matrix(x, ncol = 1L)))
           else NA_real_
    if (!is.finite(ess) || ess <= 0) ess <- length(x)
    stats::sd(x) / sqrt(ess)
  }, numeric(1))
  list(mean = m, mcse = stats::setNames(se, colnames(draws)))
}

## The distributional replay: the sampler stage of `base` from the fresh mode
## result `mv`, with at most `replay_draws` draws per sampler.
.verify_replay <- function(base, mv, rec_draws, replay_draws, quiet_cmp,
                           quiet_out, alpha) {
  sl <- .spec_sampler_list(base$sampler)
  sl <- lapply(sl, function(s) {
    if ("n_draws" %in% names(s) && s$n_draws > replay_draws)
      s <- .spec_build("sampler", list(n_draws = as.integer(replay_draws)),
                       method = s$method, base = s)
    s
  })
  md <- base$mode
  md["result"] <- list(mv)
  base["mode"] <- list(md)
  out <- quiet_out
  out$form <- "posterior"
  rs <- suppressWarnings(update(base,
    sampler = if (length(sl) == 1L) sl[[1L]] else
      structure(sl, class = "dynhr_sampler_sequence"),
    compute = quiet_cmp, outputs = out))
  pr <- suppressMessages(suppressWarnings(.run_estimation_impl(rs, rr = NULL)))
  rep_draws <- pr$pooled_draws
  nm <- intersect(colnames(rec_draws), colnames(rep_draws))
  a <- .verify_mean_mcse(rec_draws[, nm, drop = FALSE])
  b <- .verify_mean_mcse(rep_draws[, nm, drop = FALSE])
  z <- (a$mean - b$mean) / sqrt(a$mcse^2 + b$mcse^2)
  z[!is.finite(z)] <- ifelse(a$mean[!is.finite(z)] == b$mean[!is.finite(z)], 0, Inf)
  list(mean_recorded = a$mean, mean_replay = b$mean,
       mcse_recorded = a$mcse, mcse_replay = b$mcse, z = z,
       max_abs_z = max(abs(z)), worst = nm[which.max(abs(z))],
       z_crit = stats::qnorm(1 - alpha / (2 * length(nm))),
       n_replay = nrow(rep_draws))
}

#' @rdname dynhr_verify
#' @param ... Unused.
#' @export
print.dynhr_verify <- function(x, ...) {
  cat(sprintf("<dynhr_verify> %s\n", x$verdict))
  cat(sprintf("  expectation : %s\n", x$expectation))
  for (i in seq_len(nrow(x$classes))) {
    r <- x$classes[i, ]
    cat(sprintf("  class %s %-25s %-12s %s%s\n", r$class,
                paste0("(", r$what, ")"), r$status, r$verdict,
                if (nzchar(r$details)) paste0("  -- ", r$details) else ""))
  }
  for (i in seq_len(nrow(x$checks))) {
    r <- x$checks[i, ]
    res <- if (is.na(r$pass)) "skip" else if (r$pass) "ok  " else "FAIL"
    cat(sprintf("  [%s] %-20s %s  (%s)\n", res, r$check,
                if (is.na(r$value)) "" else
                  sprintf("%.3g <= %.3g", r$value, r$threshold),
                r$note))
  }
  if (length(x$registered_changes))
    cat("  registered result changes:\n",
        paste0("    ", x$registered_changes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
