## R/sampler-diagnostics.R
## --------------------------------------------------------------------------
## Unified sampler diagnostics summary object.
##
## sampler_diagnostics(fit) accepts the return list of any dynhr sampler
## (dynhr_nuts, dynhr_hmc, dynhr_mala, dynhr_chees, dynhr_smc, dynhr_dime,
## or the RWMH path of run_posterior_estimation) and returns a uniform S3
## object of class "dynhr_sampler_diagnostics".
##
## Fields that the sampler does not record are returned as NA rather than
## raising an error.
## --------------------------------------------------------------------------


# ============================================================================
# Constructor
# ============================================================================

#' Compute a unified diagnostics summary for any dynhr sampler result
#'
#' Accepts the list returned by \code{dynhr_nuts()}, \code{dynhr_hmc()},
#' \code{dynhr_mala()}, \code{dynhr_chees()}, \code{dynhr_smc()},
#' \code{dynhr_dime()}, or the RWMH path of
#' \code{run_posterior_estimation()}.  Fields that the sampler does not
#' record are returned as \code{NA} rather than raising an error.
#'
#' @param fit List returned by a dynhr sampler.
#'
#' @return An S3 object of class \code{"dynhr_sampler_diagnostics"} with
#'   the following fields (all scalars unless noted):
#'   \describe{
#'     \item{sampler}{Character; sampler name, e.g. \code{"nuts"}.
#'       \code{NA} if not recorded.}
#'     \item{n_draws}{Integer; post-warmup draw count.}
#'     \item{n_burn}{Integer; warmup / burn-in draw count.}
#'     \item{acceptance_rate}{Numeric; overall MH acceptance rate
#'       (\code{NA} when not applicable, e.g. SMC).}
#'     \item{ess_bulk}{Named numeric vector; bulk ESS per parameter.
#'       \code{NA} when chains are not available.}
#'     \item{ess_tail}{Named numeric vector; tail ESS per parameter.
#'       \code{NA} when chains are not available.}
#'     \item{ess_bulk_min}{Numeric; minimum bulk ESS across parameters.}
#'     \item{ess_tail_min}{Numeric; minimum tail ESS across parameters.}
#'     \item{rhat}{Named numeric vector; rank-normalised split R-hat.
#'       \code{NA} when fewer than 2 chains are provided.}
#'     \item{rhat_max}{Numeric; maximum R-hat across parameters.}
#'     \item{elapsed_secs}{Numeric; wall time in seconds.
#'       \code{NA} if not recorded.}
#'     \item{ess_per_sec}{Numeric; minimum bulk ESS / elapsed_secs.
#'       \code{NA} if either is unavailable.}
#'     \item{n_gradient_evals}{Integer; total gradient evaluations.
#'       \code{NA} for samplers that do not count them (RWMH, MALA, SMC,
#'       DIME).}
#'     \item{n_divergent}{Integer; number of divergent transitions.
#'       \code{NA} for samplers without divergence tracking (HMC, MALA,
#'       RWMH, SMC, DIME).}
#'     \item{mean_treedepth}{Numeric; mean NUTS tree depth.
#'       \code{NA} for non-NUTS samplers.}
#'     \item{max_treedepth}{Integer; maximum NUTS tree depth observed.
#'       \code{NA} for non-NUTS samplers.}
#'     \item{bfmi}{Numeric; Bayesian Fraction of Missing Information.
#'       \code{NA} when energy trace is unavailable (all samplers except
#'       NUTS).}
#'     \item{flags}{Character vector; human-readable warning strings for
#'       concerning values (low BFMI, divergences, high R-hat).  Empty
#'       when all diagnostics are within acceptable range.}
#'   }
#'
#' @param chains_list Optional list of per-chain draw matrices
#'   (\code{n_draws x n_params}), each named by parameter.  When provided,
#'   multi-chain R-hat is computed.  When omitted, only single-chain bulk/
#'   tail ESS is reported and R-hat is \code{NA}.
#'
#' @examples
#' \dontrun{
#'   log_post <- function(theta) -0.5 * sum(theta^2)
#'   fit <- dynhr_hmc(log_post, theta_init = c(mu = 0, sigma = 1),
#'                    n_draws = 200, n_warmup = 100, verbose = FALSE)
#'   diag <- sampler_diagnostics(fit)
#'   print(diag)
#' }
#'
#' @export
sampler_diagnostics <- function(fit, chains_list = NULL) {
  stopifnot(is.list(fit))

  # ------------------------------------------------------------------
  # 1. Identify sampler
  # ------------------------------------------------------------------
  sampler <- fit$sampler %||% NA_character_

  # ------------------------------------------------------------------
  # 2. Basic counts
  # ------------------------------------------------------------------
  n_draws <- as.integer(fit$n_draws %||% NA_integer_)
  n_burn  <- as.integer(fit$n_burn  %||% NA_integer_)

  # ------------------------------------------------------------------
  # 3. Acceptance rate
  # ------------------------------------------------------------------
  accept_rate <- fit$acceptance_rate %||% NA_real_

  # ------------------------------------------------------------------
  # 4. ESS and R-hat from post-warmup chain
  # ------------------------------------------------------------------
  draws <- fit$chain   # NULL for NUTS/MALA with checkpoint$return_chain=FALSE

  ess_bulk_vec <- NA_real_
  ess_tail_vec <- NA_real_
  rhat_vec     <- NA_real_

  if (!is.null(draws) && is.matrix(draws) && nrow(draws) >= 4L) {
    draws_m <- as.matrix(draws)
    n_par   <- ncol(draws_m)
    p_names <- colnames(draws_m)
    if (is.null(p_names)) p_names <- paste0("theta_", seq_len(n_par))
    colnames(draws_m) <- p_names

    # Bulk ESS (rank-normalised)
    ess_bulk_vec <- vapply(seq_len(n_par), function(j) {
      z <- .rank_normalise(list(draws_m[, j]))[[1]]
      .effective_sample_size(z)
    }, numeric(1L))
    names(ess_bulk_vec) <- p_names

    # Tail ESS (min of 5th / 95th quantile indicators)
    ess_tail_vec <- vapply(seq_len(n_par), function(j) {
      x   <- draws_m[, j]
      q05 <- quantile(x, 0.05); q95 <- quantile(x, 0.95)
      min(.effective_sample_size(as.numeric(x <= q05)),
          .effective_sample_size(as.numeric(x <= q95)))
    }, numeric(1L))
    names(ess_tail_vec) <- p_names

    # Multi-chain R-hat (requires chains_list)
    if (!is.null(chains_list) && length(chains_list) >= 2L) {
      chains_list_m <- lapply(chains_list, function(m) {
        m2 <- as.matrix(m)
        colnames(m2) <- p_names
        m2
      })
      cs      <- .convergence_summary(chains_list_m)
      rhat_vec <- cs$rhat
      names(rhat_vec) <- cs$param
      # Also pool all chains for better ESS
      pooled   <- do.call(rbind, chains_list_m)
      colnames(pooled) <- p_names
      ess_bulk_vec <- vapply(seq_len(n_par), function(j) {
        z <- .rank_normalise(lapply(chains_list_m, function(m) m[, j]))[[1]]
        .effective_sample_size(z)
      }, numeric(1L))
      # The pooled-chain ESS from convergence_summary is more standard
      ess_bulk_vec <- cs$ess_bulk
      ess_tail_vec <- cs$ess_tail
      names(ess_bulk_vec) <- cs$param
      names(ess_tail_vec) <- cs$param
    }
  }

  ess_bulk_min <- if (all(is.na(ess_bulk_vec))) NA_real_ else
    min(ess_bulk_vec, na.rm = TRUE)
  ess_tail_min <- if (all(is.na(ess_tail_vec))) NA_real_ else
    min(ess_tail_vec, na.rm = TRUE)
  rhat_max     <- if (all(is.na(rhat_vec))) NA_real_ else
    max(rhat_vec, na.rm = TRUE)

  # ------------------------------------------------------------------
  # 5. Timing and ESS/sec
  # ------------------------------------------------------------------
  elapsed <- fit$elapsed_secs %||% NA_real_

  ess_per_sec <- if (!is.na(ess_bulk_min) && !is.na(elapsed) && elapsed > 0)
    ess_bulk_min / elapsed
  else NA_real_

  # ------------------------------------------------------------------
  # 6. Gradient evaluations
  # ------------------------------------------------------------------
  n_gradient_evals <- as.integer(
    fit$n_grad_evals %||% fit$n_eval %||% NA_integer_
  )

  # ------------------------------------------------------------------
  # 7. HMC/NUTS-specific: divergences, treedepth, BFMI
  # ------------------------------------------------------------------
  n_divergent <- as.integer(fit$n_divergent %||% NA_integer_)

  # NUTS has treedepths vector; HMC/ChEES do not
  treedepths <- fit$treedepths   # NULL for non-NUTS
  mean_treedepth <- if (!is.null(treedepths) && length(treedepths) > 0L)
    mean(treedepths, na.rm = TRUE)
  else NA_real_
  max_treedepth <- if (!is.null(treedepths) && length(treedepths) > 0L)
    max(treedepths, na.rm = TRUE)
  else NA_integer_

  # BFMI from energy trace (NUTS only)
  energy_trace <- fit$energy_trace
  bfmi <- if (!is.null(energy_trace) && length(energy_trace) >= 3L)
    .bfmi(energy_trace)
  else NA_real_

  # ------------------------------------------------------------------
  # 8. Flags
  # ------------------------------------------------------------------
  flags <- character(0L)

  if (!is.na(bfmi) && bfmi < 0.3)
    flags <- c(flags, sprintf("Low BFMI = %.3f (< 0.3): sampler may not explore tails adequately", bfmi))

  if (!is.na(n_divergent) && n_divergent > 0L)
    flags <- c(flags, sprintf("%d divergent transition(s): posterior geometry may be problematic", n_divergent))

  if (!is.na(rhat_max) && rhat_max > 1.01)
    flags <- c(flags, sprintf("Max R-hat = %.3f (> 1.01): chains may not have converged", rhat_max))

  if (!is.na(ess_bulk_min) && ess_bulk_min < 100)
    flags <- c(flags, sprintf("Min bulk ESS = %.0f (< 100): estimates may be unreliable", ess_bulk_min))

  # ------------------------------------------------------------------
  # 9. Assemble object
  # ------------------------------------------------------------------
  structure(
    list(
      sampler          = sampler,
      n_draws          = n_draws,
      n_burn           = n_burn,
      acceptance_rate  = accept_rate,
      ess_bulk         = ess_bulk_vec,
      ess_tail         = ess_tail_vec,
      ess_bulk_min     = ess_bulk_min,
      ess_tail_min     = ess_tail_min,
      rhat             = rhat_vec,
      rhat_max         = rhat_max,
      elapsed_secs     = elapsed,
      ess_per_sec      = ess_per_sec,
      n_gradient_evals = n_gradient_evals,
      n_divergent      = n_divergent,
      mean_treedepth   = mean_treedepth,
      max_treedepth    = max_treedepth,
      bfmi             = bfmi,
      flags            = flags
    ),
    class = "dynhr_sampler_diagnostics"
  )
}


# ============================================================================
# print method
# ============================================================================

#' @export
print.dynhr_sampler_diagnostics <- function(x, ...) {
  sampler_str <- if (is.na(x$sampler)) "unknown" else toupper(x$sampler)
  cat(sprintf("=== dynhr sampler diagnostics: %s ===\n", sampler_str))

  # Basic counts
  if (!is.na(x$n_draws)) {
    burn_str <- if (!is.na(x$n_burn)) sprintf(" (%d warmup)", x$n_burn) else ""
    cat(sprintf("  Draws        : %d post-warmup%s\n", x$n_draws, burn_str))
  }

  # Timing
  if (!is.na(x$elapsed_secs)) {
    cat(sprintf("  Wall time    : %.2f sec\n", x$elapsed_secs))
  }

  # Acceptance
  if (!is.na(x$acceptance_rate)) {
    cat(sprintf("  Accept rate  : %.1f%%\n", x$acceptance_rate * 100))
  }

  # Gradient evals
  if (!is.na(x$n_gradient_evals)) {
    cat(sprintf("  Grad evals   : %d\n", x$n_gradient_evals))
  }

  # ESS
  if (!is.na(x$ess_bulk_min)) {
    ess_sec_str <- if (!is.na(x$ess_per_sec))
      sprintf("  (%.1f ESS/sec)", x$ess_per_sec) else ""
    cat(sprintf("  ESS bulk min : %.0f%s\n", x$ess_bulk_min, ess_sec_str))
  }
  if (!is.na(x$ess_tail_min)) {
    cat(sprintf("  ESS tail min : %.0f\n", x$ess_tail_min))
  }

  # R-hat
  if (!is.na(x$rhat_max)) {
    cat(sprintf("  R-hat max    : %.4f\n", x$rhat_max))
  } else {
    cat("  R-hat        : NA (single chain)\n")
  }

  # HMC/NUTS-specific
  if (!is.na(x$bfmi)) {
    cat(sprintf("  BFMI         : %.3f%s\n", x$bfmi,
                if (x$bfmi < 0.3) " [LOW]" else ""))
  }
  if (!is.na(x$n_divergent)) {
    cat(sprintf("  Divergences  : %d%s\n", x$n_divergent,
                if (x$n_divergent > 0L) " [WARN]" else ""))
  }
  if (!is.na(x$mean_treedepth)) {
    cat(sprintf("  Tree depth   : mean=%.1f  max=%d\n",
                x$mean_treedepth, x$max_treedepth))
  }

  # Flags
  if (length(x$flags) > 0L) {
    cat("\n  FLAGS:\n")
    for (f in x$flags) cat(sprintf("    * %s\n", f))
  } else {
    cat("\n  No diagnostic flags.\n")
  }

  invisible(x)
}
