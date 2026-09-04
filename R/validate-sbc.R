## R/validate-sbc.R
## --------------------------------------------------------------------------
## Simulation-Based Calibration (SBC) for the full
## prior -> simulate -> estimate pipeline.
##
## Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
## Validating Bayesian Inference Algorithms with Simulation-Based Calibration.
## arXiv:1804.06788.
##
## LAYER 1 (pure statistics, individually testable):
##   .sbc_ranks()        -- rank statistic for one replication (Talts eq. 1)
##   sbc_uniformity_test() -- chi-squared GOF test of rank uniformity
##   .sbc_rank_hist_plot() -- rank-histogram facet plot with the 99% band
##
## LAYER 2 (pipeline driver):
##   dynhr_sbc()          -- exported entry point
##   .sbc_one_replication() -- single replication (prior -> simulate -> RWMH)
##
## LAYER 3 (sampler adapters):
##   sbc_draws_from_smc()  -- extract L_target i.i.d.-equivalent draws from
##                            a weighted SMC result via systematic resampling
##   sbc_draws_from_dime() -- extract L_target draws from a DIME ensemble by
##                            per-walker thinning + concatenation
##
## LAYER 4 (matrix publication):
##   sbc_matrix_result()  -- build the sampler x likelihood status data.frame
##                           and write inst/extdata/sbc_matrix.rds
## --------------------------------------------------------------------------


## ===========================================================================
## LAYER 3: sampler adapters  (defined before LAYER 1 so they are available
## to .sbc_one_replication without forward-reference issues)
## ===========================================================================

#' Extract SBC draws from a weighted SMC result
#'
#' SMC returns a weighted particle cloud, not a Markov chain. To compute SBC
#' ranks with \code{\link{.sbc_ranks}}, the weighted particles must first be
#' resampled to yield approximately i.i.d. draws. This function applies
#' systematic resampling (\code{.smc_systematic_resample}) to the
#' \code{smc_weights} returned by \code{\link{dynhr_smc}} and returns the
#' resampled particle matrix.
#'
#' An ESS floor of \code{0.3 * n_particles} is enforced before resampling.
#' If \code{ESS < 0.3 * n_particles} the replication should be marked failed
#' (the returned list has \code{$ok = FALSE, $reason = "SMC ESS too low"}).
#'
#' @param smc_result  List returned by \code{\link{dynhr_smc}}. Must contain
#'   \code{$chain} (n_particles x d matrix) and \code{$smc_weights}
#'   (normalized weight vector summing to 1).
#' @param L_target    Positive integer: number of resampled rows to return.
#'   Should be \code{<= n_particles} for efficiency, but any positive integer
#'   is accepted (systematic resampling handles it with replacement).
#' @return A list with either:
#'   \describe{
#'     \item{\code{ok = TRUE, draws = <L_target x d matrix>}}{Resampled
#'       draw matrix, suitable for passing directly to
#'       \code{.sbc_ranks(..., thin = 1L)}.}
#'     \item{\code{ok = FALSE, reason = "SMC ESS too low"}}{ESS floor
#'       not met; caller should count the replication as failed.}
#'   }
#' @noRd
sbc_draws_from_smc <- function(smc_result, L_target) {
  stopifnot(is.list(smc_result))
  stopifnot(!is.null(smc_result$chain), !is.null(smc_result$smc_weights))

  w           <- smc_result$smc_weights
  n_particles <- length(w)
  L_target    <- as.integer(L_target)
  if (L_target < 1L) stop("sbc_draws_from_smc: L_target must be >= 1.")

  ## ESS floor: require ESS > 0.3 * n_particles (Landmines section of brief)
  ess <- 1 / sum(w^2)
  if (ess < 0.3 * n_particles)
    return(list(ok = FALSE, reason = "SMC ESS too low"))

  idx   <- .smc_systematic_resample(w, L_target)
  draws <- smc_result$chain[idx, , drop = FALSE]
  list(ok = TRUE, draws = draws)
}


#' Extract SBC draws from a DIME ensemble result
#'
#' DIME returns a flattened \code{[n_iter * n_chain, d]} post-burn matrix
#' stored row-major by iteration then walker (i.e. rows 1:n_chain are all
#' walkers at iteration 1, rows (n_chain+1):(2*n_chain) are iteration 2,
#' etc.). Within-walker autocorrelation and within-iteration cross-walker
#' correlation both need to be handled.
#'
#' This function thin each walker's sub-chain independently by
#' \code{thin_k = ceiling(n_iter / ceil(L_target / n_chain))} and
#' concatenates, yielding approximately \code{L_target} approximately
#' independent draws.
#'
#' @param dime_result  List returned by \code{\link{run_dime}}. Must contain
#'   \code{$chain} (n_iter * n_chain x d matrix) and \code{$n_walkers}
#'   (integer = n_chain).
#' @param L_target     Positive integer: approximate number of draws to
#'   return.
#' @return An approximately \code{L_target x d} draw matrix (may be slightly
#'   more or fewer due to integer rounding), suitable for
#'   \code{.sbc_ranks(..., thin = 1L)}.
#' @noRd
sbc_draws_from_dime <- function(dime_result, L_target) {
  stopifnot(is.list(dime_result))
  stopifnot(!is.null(dime_result$chain), !is.null(dime_result$n_walkers))

  post_chain <- dime_result$chain
  n_walkers  <- as.integer(dime_result$n_walkers)
  L_target   <- as.integer(L_target)
  n_total    <- nrow(post_chain)
  n_iter     <- n_total %/% n_walkers

  if (n_iter < 1L)
    stop("sbc_draws_from_dime: post_chain has fewer rows than n_walkers.")
  if (L_target < 1L)
    stop("sbc_draws_from_dime: L_target must be >= 1.")

  ## Per-walker thinning: each walker contributes ceil(L_target / n_walkers)
  ## rows. Thin by stride = ceiling(n_iter / ceil(L_target / n_walkers)).
  per_walker <- ceiling(L_target / n_walkers)
  thin_k     <- max(1L, ceiling(n_iter / per_walker))

  ## Reshape [n_iter, n_walkers, d] and thin each walker independently.
  d         <- ncol(post_chain)
  par_names <- colnames(post_chain)
  kept_rows <- list()
  for (j in seq_len(n_walkers)) {
    ## rows for walker j (0-indexed walker): all rows t where (row-1) %% n_walkers == j-1
    walker_rows <- seq(j, n_total, by = n_walkers)  ## rows 1,1+n_walkers,... for walker 1
    thinned     <- walker_rows[seq(1L, length(walker_rows), by = thin_k)]
    kept_rows[[j]] <- thinned
  }

  all_rows  <- unlist(kept_rows, use.names = FALSE)
  draws     <- post_chain[all_rows, , drop = FALSE]
  rownames(draws) <- NULL
  draws
}


## ===========================================================================
## LAYER 1: rank statistic
## ===========================================================================

#' Compute the SBC rank statistic for a single replication
#'
#' For each parameter j, the rank is the number of posterior draws (after
#' thinning) that fall below the "true" (prior-drawn) value
#' \code{theta_tilde[j]} (Talts et al. 2018, eq. 1):
#' \deqn{r_j = \sum_{l=1}^{L'} \mathbb{1}[\theta^{(l)}_j < \tilde\theta_j]}
#' where \eqn{\theta^{(l)}} are the (thinned) posterior draws. For a
#' well-calibrated, continuous posterior, \eqn{r_j} is uniformly distributed
#' on \eqn{\{0, 1, \ldots, L'\}}. Ties (a draw exactly equal to
#' \code{theta_tilde[j]}) are measure-zero for continuous posteriors and are
#' simply not counted (i.e. treated as \code{>= theta_tilde[j]}).
#'
#' @param theta_tilde Named numeric vector of length d: the "true" parameter
#'   draw from the prior used to simulate the data for this replication.
#' @param draws_mat   L x d numeric matrix of posterior draws, columns
#'   matching (a superset of, in any order) \code{names(theta_tilde)}.
#' @param thin        Integer thinning interval (default 1L, i.e. no
#'   thinning). Draws are thinned by taking every \code{thin}-th row,
#'   \code{draws_mat[seq(1, L, by = thin), , drop = FALSE]}.
#' @return Named integer vector of length d (names = \code{names(theta_tilde)})
#'   with each rank in \code{0:L'}, where \code{L'} is the number of
#'   thinned draws.
#' @noRd
.sbc_ranks <- function(theta_tilde, draws_mat, thin = 1L) {
  if (is.null(names(theta_tilde)) || any(names(theta_tilde) == ""))
    stop(".sbc_ranks: theta_tilde must be a fully named vector.")
  if (!all(names(theta_tilde) %in% colnames(draws_mat)))
    stop(".sbc_ranks: draws_mat is missing column(s) for: ",
         paste(setdiff(names(theta_tilde), colnames(draws_mat)), collapse = ", "))

  thin <- as.integer(thin)
  if (is.na(thin) || thin < 1L) stop(".sbc_ranks: thin must be a positive integer.")

  L <- nrow(draws_mat)
  keep <- seq.int(1L, L, by = thin)
  draws_thin <- draws_mat[keep, , drop = FALSE]

  ranks <- vapply(names(theta_tilde), function(nm) {
    sum(draws_thin[, nm] < theta_tilde[[nm]])
  }, integer(1))
  names(ranks) <- names(theta_tilde)
  ranks
}


## ===========================================================================
## LAYER 1: uniformity test
## ===========================================================================

#' Chi-squared goodness-of-fit test for SBC rank uniformity, plus
#' shift/tail-asymmetry/saturation diagnostics
#'
#' For each parameter (column of \code{ranks_mat}), bins the ranks into
#' \code{n_bins} bins via \code{floor(rank * n_bins / (L' + 1))} (so each
#' rank in \code{0:L'} maps to a bin in \code{0:(n_bins - 1)}), and runs a
#' chi-squared goodness-of-fit test against the (exact, generally unequal --
#' see "Bin-count exactness" below) expected counts under uniformity.
#'
#' \strong{Motivation (2026-07-02 P2c incident;
#' \code{ORDER3_PRUNED_SS_FOLLOWUP.md} "FINAL VERDICT").} The bare chi-squared
#' verdict both FALSE-ALARMED (a batch with p = 0.011 driven partly by rank
#' noise) and UNDER-DETECTED (a "calibrated" batch with a one-sided top-bin
#' excess that replicated across reruns). The decisive signal in both cases
#' was structural -- an elevated MEAN rank (shifted posterior) and a
#' ONE-SIDED tail excess -- neither of which the omnibus chi-squared
#' statistic is powerful against (it sees only squared deviations, blind to
#' sign and to a whole-distribution shift's location vs. spread). This
#' function therefore adds two targeted diagnostics (\code{mean_rank_z},
#' \code{tail_asym_z}) computed per parameter alongside the existing
#' chi-squared table, and folds them into the verdict logic (see below).
#'
#' \strong{Caveat} (Talts et al. 2018, sec 4.1; Sailynoja et al. 2022): the
#' chi-squared test on SBC rank bins is an approximate, omnibus check. It can
#' have low power against specific, structured miscalibrations (e.g. ranks
#' that are uniform marginally per parameter but jointly non-uniform), and a
#' single non-significant p-value does not certify calibration. The rank
#' histograms (see \code{\link{dynhr_sbc}}'s \code{plot} element) should
#' always be visually inspected for systematic departures from uniformity
#' (U-shapes indicating posteriors too narrow, hump shapes indicating
#' posteriors too wide, or trends indicating bias) -- this is precisely what
#' the chi-squared test, applied per-parameter, can miss when looking only at
#' p-values.
#'
#' \strong{Bin-count exactness.} Ranks live on \code{0:L'} (\code{L' + 1}
#' distinct integer values); when \code{L' + 1} is not divisible by
#' \code{n_bins}, the \code{floor(rank * n_bins / (L' + 1))} binning gives
#' bins unequal *integer* coverage (e.g. \code{L' = 160}, \code{n_bins = 9}
#' gives 8 bins of width 18 and one of width 17). The expected per-bin count
#' under uniformity is computed EXACTLY from this integer partition
#' (\code{n_repl * (bin width) / (L' + 1)}), not the naive \code{n_repl /
#' n_bins}.
#'
#' \strong{Diagnostics (per parameter).}
#' \describe{
#'   \item{\code{mean_rank_z}}{\code{z = (mean(rank)/L' - 0.5) *
#'     sqrt(12 * n_repl)}. Under uniformity, \code{rank/L'} has mean 0.5 and
#'     variance \code{1/12}, so this is (asymptotically) a standard normal
#'     statistic sensitive to a SHIFTED posterior (elevated or depressed mean
#'     rank) -- a signature the chi-squared statistic, which only sees
#'     squared per-bin deviations, is not targeted at.}
#'   \item{\code{tail_asym_z}}{With \code{n_top} / \code{n_bottom} the counts
#'     in the last/first bin, \code{z = (n_top - n_bottom) /
#'     sqrt(n_top + n_bottom)} (0 if \code{n_top + n_bottom == 0}). Large
#'     \code{|z|} indicates a ONE-SIDED tail excess (consistent with a shift);
#'     a symmetric U-shape (both tails inflated) gives \code{z ~ 0} even
#'     though the chi-squared statistic flags it.}
#'   \item{\code{extreme_frac}}{Fraction of ranks exactly equal to 0 or
#'     \code{L'} -- saturated ranks, a symptom of chain ESS much smaller than
#'     \code{L'} (ESS << L collapses the rank distribution onto the
#'     endpoints).}
#' }
#'
#' \strong{Verdict logic (composite).} \code{"miscalibrated"} if EITHER the
#' chi-squared p-value is Bonferroni-significant for any parameter, OR
#' \code{|mean_rank_z| > 3}, OR \code{|tail_asym_z| > 3} for any parameter.
#' \code{"suspect"} (new) if none of those trip but any of
#' \code{|mean_rank_z|}, \code{|tail_asym_z|} lies in \code{[2, 3]} for some
#' parameter. Else \code{"calibrated"}. This composite catches shift/tail
#' patterns the chi-squared alone can miss at moderate \code{n_repl} (see the
#' P2c motivation above), while still flagging anything the chi-squared alone
#' would have caught. \code{"insufficient"} (new) if the input carries no
#' usable rank-uniformity signal: either every observed rank is identical
#' (\code{L' == 0}, e.g. a single posterior draw or too few replications to
#' see any spread) or every parameter's chi-squared test degenerated to
#' \code{NA} (e.g. \code{n_bins} finer than the number of distinct integer
#' ranks). This never crashes -- it degrades gracefully instead of leaking
#' an \code{NA} into the verdict comparison. Zero replications (\code{nrow(
#' ranks_mat) == 0}) is instead treated as caller error and raises directly.
#'
#' @param ranks_mat \code{n_repl x d} integer matrix of SBC ranks (one row
#'   per replication, one column per parameter), each entry in
#'   \code{0:L'}. Column names are taken as parameter names.
#' @param L Integer: the TRUE rank support \code{L'} (ranks live in
#'   \code{0:L}), e.g. the harness's actual kept-draw count
#'   (\code{thin_L}/\code{L_effective}). \code{NULL} (default) falls back to
#'   \code{max(ranks_mat)} with a warning -- inferring the support from the
#'   OBSERVED ranks is WRONG whenever the top rank never appears by chance
#'   (e.g. no replication's kept draws all landed above theta*), which
#'   silently narrows the bins and can hide a real miscalibration signal.
#'   Callers that know their harness's rank support should always pass it
#'   explicitly.
#' @param n_bins    Number of histogram bins. Default: a divisor-friendly
#'   choice, \code{min(20, floor(n_repl / 5))}, with a floor of 2.
#' @param draws     Optional: a list of per-replication chain matrices (or a
#'   single representative \code{n_draws x d} chain matrix) used ONLY to
#'   estimate per-parameter chain ESS via \code{.effective_sample_size()}
#'   (initial-positive-sequence estimator). When supplied, \code{$ess} is
#'   populated and the print method warns if \code{ess < 5 * L'} ("ranks are
#'   noise-dominated"). Default \code{NULL} (no ESS diagnostic).
#' @return A list with class \code{"dynhr_sbc_uniformity"}:
#'   \describe{
#'     \item{table}{\code{data.frame(parameter, chisq, df, p_value,
#'       mean_rank_z, tail_asym_z, extreme_frac)}.}
#'     \item{n_bins}{Number of bins used.}
#'     \item{alpha}{Nominal family-wise significance level (0.05).}
#'     \item{alpha_bonferroni}{Bonferroni-adjusted per-test threshold,
#'       \code{alpha / d}.}
#'     \item{alpha_sidak}{Sidak-adjusted per-test threshold,
#'       \code{1 - (1 - alpha)^(1/d)}.}
#'     \item{verdict}{\code{"calibrated"}, \code{"suspect"},
#'       \code{"miscalibrated"}, or \code{"insufficient"} -- see "Verdict
#'       logic" above.}
#'     \item{ess}{Named numeric vector of per-parameter chain ESS, or
#'       \code{NULL} if \code{draws} was not supplied.}
#'   }
#' @noRd
sbc_uniformity_test <- function(ranks_mat, L = NULL, n_bins = NULL, draws = NULL) {
  if (is.null(dim(ranks_mat))) ranks_mat <- matrix(ranks_mat, ncol = 1)
  n_repl <- nrow(ranks_mat)
  d      <- ncol(ranks_mat)
  par_names <- colnames(ranks_mat)
  if (is.null(par_names)) par_names <- paste0("theta", seq_len(d))

  ## Fail loud on truly unusable input (garbage in, not a graceful verdict):
  ## no replications means there is nothing to test.
  if (n_repl < 1L) {
    stop("sbc_uniformity_test: ranks_mat has 0 rows -- need at least 1 ",
         "SBC replication to test.")
  }

  if (is.null(n_bins)) {
    n_bins <- max(2L, min(20L, floor(n_repl / 5)))
  }
  n_bins <- as.integer(n_bins)
  if (n_bins < 2L) stop("sbc_uniformity_test: n_bins must be >= 2.")

  ## True rank support L' (#7 one-liner, adversarial review): inferring it
  ## from max(ranks_mat) is wrong whenever the top rank never appears by
  ## chance -- silently narrowing the bins and potentially hiding a real
  ## miscalibration signal. Callers that know their harness's rank support
  ## (thin_L / L_effective / the kept-draw count) should pass it via `L`.
  if (is.null(L)) {
    warning("sbc_uniformity_test: `L` (true rank support) not supplied -- ",
            "inferring L = max(ranks_mat) from the OBSERVED ranks. This is ",
            "WRONG whenever the top rank never appears by chance (no ",
            "replication's kept draws happened to land above theta*), which ",
            "silently narrows the bins and can hide a real miscalibration ",
            "signal. Pass the harness's actual rank support explicitly via ",
            "the `L` argument.", call. = FALSE)
    L_eff <- max(ranks_mat, na.rm = TRUE)         ## L' (observed fallback)
  } else {
    if (!is.numeric(L) || length(L) != 1L || !is.finite(L) || L < 0 ||
        L != as.integer(L)) {
      stop("sbc_uniformity_test: `L` must be a non-negative integer scalar ",
           "(the true rank support; ranks live in 0:L).")
    }
    L_eff <- as.integer(L)
    if (max(ranks_mat, na.rm = TRUE) > L_eff) {
      stop("sbc_uniformity_test: an observed rank (", max(ranks_mat, na.rm = TRUE),
           ") exceeds the supplied `L` (", L_eff, ") -- check that `L` matches ",
           "the harness's actual rank support.")
    }
  }
  L_plus1 <- L_eff + 1L                         ## L' + 1

  ## Degenerate case: every observed rank is identical (e.g. L' == 0, or a
  ## tiny n_repl that happened to land on one bucket) -- the exact-binning
  ## scheme below then produces expected_exact == 0 for some bin with
  ## counts == 0 there too, i.e. a 0/0 chisq contribution (NaN), which would
  ## otherwise leak an NA into the verdict if()-chain. There is no usable
  ## rank *distribution* to test against uniformity in this case, so report
  ## an "insufficient" verdict rather than crash or silently mis-flag.
  if (L_eff < 1L) {
    tab <- data.frame(parameter = par_names, chisq = NA_real_, df = NA_integer_,
                       p_value = NA_real_, mean_rank_z = NA_real_,
                       tail_asym_z = NA_real_, extreme_frac = NA_real_,
                       stringsAsFactors = FALSE)
    return(structure(
      list(
        table            = tab,
        n_bins           = n_bins,
        alpha            = 0.05,
        alpha_bonferroni = NA_real_,
        alpha_sidak      = NA_real_,
        verdict          = "insufficient",
        ess              = NULL,
        L_effective      = L_eff
      ),
      class = "dynhr_sbc_uniformity"
    ))
  }

  ## Exact expected bin counts: partition the L'+1 integer ranks 0:L' by the
  ## SAME binning rule used for the data, then expected count for bin b is
  ## n_repl * (# integer ranks mapping to b) / (L' + 1). This differs from
  ## n_repl / n_bins whenever (L' + 1) is not divisible by n_bins.
  all_ranks_bins <- floor((0:L_eff) * n_bins / L_plus1)
  all_ranks_bins <- pmin(all_ranks_bins, n_bins - 1L)
  bin_widths     <- tabulate(all_ranks_bins + 1L, nbins = n_bins)  ## integer counts per bin, sums to L'+1
  expected_exact <- n_repl * bin_widths / L_plus1

  rows <- lapply(seq_len(d), function(j) {
    ranks_j <- ranks_mat[, j]
    bins <- floor(ranks_j * n_bins / L_plus1)
    bins <- pmin(bins, n_bins - 1L)  ## guard rank == L' edge case
    counts <- tabulate(bins + 1L, nbins = n_bins)
    ## Guard 0/0: a bin with expected_exact == 0 (n_bins finer than the
    ## number of distinct integer ranks, L' + 1) necessarily also has
    ## counts == 0 there (no integer rank maps to it), so its GOF
    ## contribution is a no-information 0, not NaN -- without this guard,
    ## an NaN chisq_val propagates to p_value = NA and then leaks into the
    ## verdict if()-chain below (missing value where TRUE/FALSE needed).
    zero_exp <- expected_exact == 0
    chisq_val <- sum(((counts - expected_exact)^2 / expected_exact)[!zero_exp])
    df <- sum(!zero_exp) - 1L
    ## df < 1 (<= 1 bin with any expected mass) leaves no GOF test to run --
    ## NA out the chisq/p_value rather than call pchisq() with a degenerate
    ## df (which itself can return NaN and leak into the verdict below).
    p_val <- if (df < 1L) NA_real_ else
      stats::pchisq(chisq_val, df = df, lower.tail = FALSE)

    ## mean_rank_z: detects a SHIFTED posterior (elevated/depressed mean rank).
    ## rank/L' ~ Uniform[0,1] under calibration => mean 0.5, sd 1/sqrt(12).
    mean_rank_z <- (mean(ranks_j) / L_eff - 0.5) * sqrt(12 * n_repl)

    ## tail_asym_z: one-sided tail excess (large |z|) vs. a symmetric
    ## U-shape (z ~ 0, caught by chisq instead). The first/last bins need
    ## NOT hold the same number of integer ranks (floor-binning of L'+1
    ## ranks into n_bins), so the raw difference n_top - n_bottom has a
    ## NONZERO null mean whenever the edge bin widths differ -- with
    ## thin_L = 25 and n_bins = 20 the null mean is approx -1.6 z-units,
    ## enough to flip perfectly uniform ranks to "suspect"/"miscalibrated"
    ## (caught 2026-08-04 by an exact-likelihood SBC design oracle; the
    ## chisq path was always width-exact via expected_exact). Center by the
    ## exact expected counts and scale by the binomial-difference sd.
    n_bottom <- counts[1L]
    n_top    <- counts[n_bins]
    e_bottom <- expected_exact[1L]
    e_top    <- expected_exact[n_bins]
    p_bottom <- bin_widths[1L] / L_plus1
    p_top    <- bin_widths[n_bins] / L_plus1
    var_diff <- n_repl * (p_top * (1 - p_top) + p_bottom * (1 - p_bottom) +
                            2 * p_top * p_bottom)
    tail_asym_z <- if (var_diff <= 0) 0 else
      ((n_top - e_top) - (n_bottom - e_bottom)) / sqrt(var_diff)

    ## extreme_frac: saturated ranks (rank == 0 or L'), symptom of ESS << L'.
    extreme_frac <- mean(ranks_j == 0L | ranks_j == L_eff)

    data.frame(parameter = par_names[j], chisq = chisq_val, df = df,
               p_value = p_val, mean_rank_z = mean_rank_z,
               tail_asym_z = tail_asym_z, extreme_frac = extreme_frac,
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)

  alpha <- 0.05
  alpha_bonf  <- alpha / d
  alpha_sidak <- 1 - (1 - alpha)^(1 / d)

  ## NA-safe backstop: a per-parameter p_value can still be NA here (e.g.
  ## df < 1 for a pathologically bin-starved parameter column, guarded
  ## above) even when the overall L_eff >= 1 gate passed. Treat NA flags as
  ## FALSE (no evidence, not a miscalibration signal) rather than let them
  ## propagate into any(...) -- any(NA) is NA and would crash the if()
  ## below (the original bug report's exact symptom).
  chisq_flag <- isTRUE_vec(tab$p_value <= alpha_bonf)
  shift_flag <- isTRUE_vec(abs(tab$mean_rank_z) > 3 | abs(tab$tail_asym_z) > 3)
  suspect_flag <- !(chisq_flag | shift_flag) &
    isTRUE_vec(abs(tab$mean_rank_z) >= 2 | abs(tab$tail_asym_z) >= 2)

  ## If EVERY parameter's chisq test was unusable (all p_value NA), there is
  ## no uniformity evidence at all -- report "insufficient" rather than the
  ## misleadingly confident "calibrated" a naive all-FALSE flag vector would
  ## otherwise produce.
  verdict <- if (all(is.na(tab$p_value))) {
    "insufficient"
  } else if (any(chisq_flag | shift_flag)) {
    "miscalibrated"
  } else if (any(suspect_flag)) {
    "suspect"
  } else {
    "calibrated"
  }

  ## Optional ESS diagnostic. `draws` can be a single n_draws x d chain
  ## matrix (representative chain) or a list of such matrices (one per
  ## replication, e.g. for a spot-check subset) -- in the list case, ESS is
  ## computed per-chain per-parameter and averaged.
  ess <- NULL
  if (!is.null(draws)) {
    draws_list <- if (is.list(draws) && !is.matrix(draws)) draws else list(draws)
    ess_mat <- vapply(draws_list, function(ch) {
      ch <- as.matrix(ch)
      cn <- colnames(ch)
      vapply(seq_len(d), function(j) {
        col <- if (!is.null(cn) && par_names[j] %in% cn) ch[, par_names[j]] else ch[, j]
        .effective_sample_size(as.numeric(col))
      }, numeric(1))
    }, numeric(d))
    ess <- if (is.matrix(ess_mat)) rowMeans(ess_mat) else mean(ess_mat)
    names(ess) <- par_names
  }

  structure(
    list(
      table            = tab,
      n_bins           = n_bins,
      alpha            = alpha,
      alpha_bonferroni = alpha_bonf,
      alpha_sidak      = alpha_sidak,
      verdict          = verdict,
      ess              = ess,
      L_effective      = L_eff
    ),
    class = "dynhr_sbc_uniformity"
  )
}


#' @exportS3Method
#' @noRd
print.dynhr_sbc_uniformity <- function(x, ...) {
  cat("SBC rank uniformity test\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("Bins: %d\n\n", x$n_bins))

  tab <- x$table
  tab$chisq       <- round(tab$chisq, 3)
  tab$p_value     <- signif(tab$p_value, 4)
  tab$mean_rank_z <- round(tab$mean_rank_z, 3)
  tab$tail_asym_z <- round(tab$tail_asym_z, 3)
  tab$extreme_frac <- signif(tab$extreme_frac, 3)
  print(tab, row.names = FALSE)

  cat(sprintf("\nBonferroni-adjusted alpha: %.5f (Sidak: %.5f)\n",
              x$alpha_bonferroni, x$alpha_sidak))
  cat(sprintf("Verdict: %s\n", x$verdict))
  cat("(miscalibrated: Bonferroni-significant chisq OR |mean_rank_z| > 3 OR\n")
  cat(" |tail_asym_z| > 3; suspect: none of those, but one lies in [2, 3];\n")
  cat(" insufficient: no usable rank spread/GOF signal -- more replications\n")
  cat(" needed, not a calibration finding)\n")

  if (!is.null(x$ess)) {
    l_floor <- 5 * x$L_effective
    low_ess <- x$ess[x$ess < l_floor]
    if (length(low_ess) > 0) {
      cat(sprintf(
        "\nWARNING: chain ESS < 5 * L' (%.0f) for: %s -- ranks are\n",
        l_floor, paste(names(low_ess), sprintf("(ESS=%.0f)", low_ess), collapse = ", ")))
      cat("noise-dominated; lengthen/thin chains or use exact-quadrature PITs.\n")
    }
  }

  cat("\nNote: chi-squared p-values on SBC ranks are an approximate, omnibus\n")
  cat("check. Always inspect the rank histograms for systematic departures\n")
  cat("from uniformity.\n")

  invisible(x)
}


## ===========================================================================
## LAYER 1: rank-histogram plot
## ===========================================================================

#' Build SBC rank-histogram facet plot with 99% expected-count band
#'
#' One facet per parameter. Each facet shows a histogram of that parameter's
#' SBC ranks (binned as in \code{sbc_uniformity_test}), overlaid with
#' a shaded band giving the 99% expected-count interval for a uniform
#' histogram, \code{qbinom(c(.005, .995), n_repl, 1 / n_bins)}.
#'
#' @param ranks_mat \code{n_repl x d} integer matrix of SBC ranks.
#' @param n_bins    Number of bins (see \code{sbc_uniformity_test}).
#' @return A ggplot object, or \code{NULL} if the \pkg{ggplot2} package is
#'   not available.
#' @noRd
.sbc_rank_hist_plot <- function(ranks_mat, n_bins = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

  if (is.null(dim(ranks_mat))) ranks_mat <- matrix(ranks_mat, ncol = 1)
  n_repl <- nrow(ranks_mat)
  d      <- ncol(ranks_mat)
  par_names <- colnames(ranks_mat)
  if (is.null(par_names)) par_names <- paste0("theta", seq_len(d))

  if (is.null(n_bins)) n_bins <- max(2L, min(20L, floor(n_repl / 5)))
  n_bins <- as.integer(n_bins)

  L_plus1 <- max(ranks_mat, na.rm = TRUE) + 1L

  hist_df <- do.call(rbind, lapply(seq_len(d), function(j) {
    ranks_j <- ranks_mat[, j]
    bins <- floor(ranks_j * n_bins / L_plus1)
    bins <- pmin(bins, n_bins - 1L)
    counts <- tabulate(bins + 1L, nbins = n_bins)
    data.frame(parameter = par_names[j], bin = seq_len(n_bins) - 1L,
               count = counts, stringsAsFactors = FALSE)
  }))
  hist_df$parameter <- factor(hist_df$parameter, levels = par_names)

  expected <- n_repl / n_bins
  band <- stats::qbinom(c(0.005, 0.995), size = n_repl, prob = 1 / n_bins)

  band_df <- data.frame(ymin = band[1], ymax = band[2])

  col_band <- "#CCCCCC"
  col_bar  <- if (exists("dynhr_colours", mode = "list")) {
    tryCatch(dynhr_colours$mid_blue %||% "#1B7CB6", error = function(e) "#1B7CB6")
  } else {
    "#1B7CB6"
  }

  p <- ggplot2::ggplot(hist_df, ggplot2::aes(x = bin, y = count)) +
    ggplot2::geom_rect(
      data = band_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = col_band, alpha = 0.4
    ) +
    ggplot2::geom_hline(yintercept = expected, linetype = "dashed",
                        colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_col(fill = col_bar, width = 0.9) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(
      title = "SBC rank histograms",
      subtitle = sprintf(
        "%d replications, %d bins; shaded band = 99%% expected range under uniformity",
        n_repl, n_bins),
      x = "Rank bin", y = "Count")

  if (exists("theme_dynhr", mode = "function")) {
    p <- p + theme_dynhr()
  }
  p
}


## ===========================================================================
## LAYER 2: pipeline driver
## ===========================================================================

#' Single SBC replication: prior draw -> simulate -> posterior -> ranks
#'
#' Structured as a single function so a future parallel driver can dispatch
#' replications via \code{mirai} with a one-line change (each call is
#' self-contained given \code{model}/\code{compiled}/\code{prior_spec}).
#'
#' @param i             Replication index (used only for messages/seeding by
#'   the caller).
#' @param model         Parsed dynhr_mod.
#' @param compiled      dynhr_compiled (from \code{compile_model()}).
#' @param prior_spec    Prior specification data.frame.
#' @param prior_sampler Function() -> named numeric vector, drawn from the
#'   prior (from \code{.smc_make_prior_sampler()}).
#' @param obs_vars      Character vector of observed variable names.
#' @param T_obs         Number of observation periods to simulate.
#' @param presample     Number of pre-sample burn-in periods for the latent
#'   state (default 50L); skipped (forced to 0) for unit-root models, but
#'   unit-root models are rejected upstream in \code{dynhr_sbc()}.
#' @param n_draws,n_burn,thin Sampler / SBC settings, see \code{dynhr_sbc()}.
#' @param sampler       Currently only \code{"rwmh"} is supported.
#' @param me_variance,lik_init,transform_params,adapt_cov Forwarded to
#'   \code{make_log_posterior()} / \code{rwmh()} / \code{build_param_transform()}.
#' @param likelihood    Likelihood type: \code{"gaussian"} (default),
#'   \code{"cumulant"}, \code{"whittle"}, or \code{"tpf"}.
#' @param order         Perturbation order for simulating the DGP. For
#'   \code{likelihood = "tpf"} the DGP should be order 2 (pruned) to match
#'   the likelihood order; default 1 matches the Kalman filter.
#' @param seed          RNG seed for this replication (caller passes
#'   \code{seed_base + i}).
#' @return A list with either \code{$ok = TRUE, $ranks = <named integer
#'   vector>, $L_effective = <int>} or \code{$ok = FALSE, $reason =
#'   <character>}.
#' @noRd
.sbc_one_replication <- function(i, model, compiled, prior_spec, prior_sampler,
                                  obs_vars, T_obs, presample,
                                  n_draws, n_burn, thin, sampler,
                                  me_variance, lik_init,
                                  transform_params, adapt_cov,
                                  likelihood = "gaussian",
                                  order = 1L,
                                  seed, verbose,
                                  innovation_check = FALSE,
                                  ...) {

  set.seed(seed)

  ## ---- Step 1: prior draw -----------------------------------------------
  theta_tilde <- prior_sampler()

  ## Inject estimated shock stds (Dynare `stderr <shock>`) as well as structural
  ## params, so the SBC data-generating process uses the DRAWN std -- not the
  ## frozen calibration. Without this, SBC for shock-std params is invalid (the
  ## synthetic data ignores theta_tilde's stds; see the P0 likelihood bug).
  params <- .apply_theta_to_params(model, theta_tilde)

  ## ---- Step 2: solve the model at theta_tilde ---------------------------
  ss <- tryCatch(
    solve_steady(compiled, params, endo_names = model$var_names,
                 exo_names = model$varexo_names, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(ss) || !isTRUE(ss$converged))
    return(list(ok = FALSE, reason = "steady state did not converge"))

  dr <- tryCatch(
    solve_perturbation(model, compiled, ss$values, params, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(dr) || !isTRUE(dr$bk_satisfied))
    return(list(ok = FALSE, reason = "perturbation solve failed or BK violated"))

  ## ---- Step 3: simulate Y from the state space at theta_tilde -----------
  ## Mirror kalman_filter()'s state-space convention exactly (R/kalman-filter.R,
  ## ~lines 460-480): TT/RR from the state rows of ghx/ghu, ZZ/DD from the
  ## observed-variable rows, d = steady state of the observed variables.
  state_idx <- dr$state_idx
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  n_state   <- length(state_idx)
  n_exo     <- length(exo)

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("Observed variables not found in model: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  ghx <- dr$ghx; ghu <- dr$ghu
  TT  <- ghx[state_idx, , drop = FALSE]
  RR  <- ghu[state_idx, , drop = FALSE]
  ZZ  <- ghx[obs_idx,   , drop = FALSE]
  DD  <- ghu[obs_idx,   , drop = FALSE]
  d_obs <- dr$ys[obs_vars]

  Sigma_e <- .get_shock_cov(model, exo, params)

  ## Unit-root check, mirroring kalman_filter's lik_init = "auto" resolution.
  tt_evals <- eigen(TT, only.values = TRUE)$values
  is_unit_root <- any(Mod(tt_evals) > 1 - 1e-6)
  if (is_unit_root) {
    stop("dynhr_sbc: model has a unit root (max |eigenvalue(TT)| = ",
         signif(max(Mod(tt_evals)), 6), " >= 1 - 1e-6). SBC requires a ",
         "proper (stationary) data-generating process to draw s_0 from; ",
         "the diffuse prior used by the Kalman filter for unit-root models ",
         "does not correspond to any proper simulation distribution. ",
         "SBC is not applicable to this model as specified.")
  }

  n_total <- presample + T_obs
  n_obs   <- length(obs_vars)

  ## Cholesky factor of Sigma_e for drawing e_t ~ N(0, Sigma_e). Sigma_e may
  ## be singular if some shocks have zero variance; use a robust factor.
  L_e <- tryCatch(t(chol(Sigma_e)), error = function(e) {
    eg <- eigen(Sigma_e, symmetric = TRUE)
    vals <- pmax(eg$values, 0)
    eg$vectors %*% diag(sqrt(vals), nrow = length(vals))
  })

  ## ---- Simulate DGP at the requested perturbation order ------------------
  ## order = 1 (default): linear state-space simulation matching the KF.
  ## order = 2: pruned second-order simulation matching the TPF likelihood.
  ## For SBC with likelihood = "tpf", pass order = 2L so the DGP and
  ## likelihood are at the same approximation order (brief section 6).
  Y <- matrix(NA_real_, nrow = n_obs, ncol = T_obs)

  if (as.integer(order) >= 2L) {
    ## Order-2 DGP: solve second-order and simulate via simulate_model_order2()
    dr2 <- tryCatch(
      solve_perturbation_order2(model, compiled, ss$values, params,
                                 dr, Sigma_e = Sigma_e, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(dr2))
      return(list(ok = FALSE, reason = "order-2 perturbation solve failed"))

    ## simulate_model_order2 returns a (n_total - burn_in) x n_endo matrix
    ## of deviations from SS (levels stored in attr(., "levels")).
    sim_full <- tryCatch(
      simulate_model_order2(dr2, n_periods = T_obs, model = model,
                             burn_in = presample),
      error = function(e) NULL
    )
    if (is.null(sim_full))
      return(list(ok = FALSE, reason = "order-2 simulation failed"))

    ## Extract observed variables (deviations + SS = levels already in attr)
    sim_levels <- attr(sim_full, "levels")
    Y <- t(sim_levels[, endo[obs_idx], drop = FALSE])  # n_obs x T_obs

  } else {
    ## Order-1 DGP: linear state-space simulation (original SBC code)
    s <- numeric(n_state)

    if (likelihood == "pskf") {
      ## PSKF DGP: draw CSN shocks consistent with the drawn alpha.
      ## DGP must match likelihood: if alpha != 0 the filter expects skewed
      ## shocks; Gaussian shocks would cause a DGP-likelihood mismatch and
      ## the SBC rank distribution would be degenerate.
      ##
      ## JOINT CSN draw (Tier 10 item 3): e ~ CSN(mu_e, Sigma_e, Gamma_e, 0, I)
      ## with Gamma_e = diag(alpha_i / sigma_i), seed cov = FULL Sigma_e.  This
      ## is the SAME law the likelihood builds in .get_csn_shock_params (off-
      ## diagonals carried entirely by the seed Sigma_e), so the DGP and the
      ## likelihood are consistent for correlated + skewed shocks (no refusal).
      ##
      ## Stochastic representation: draw [e; U] jointly Gaussian, keep draws with
      ## U >= 0 componentwise.  Joint covariance:
      ##   Cov = [[Sigma_e,         Sigma_e Gamma_e'],
      ##          [Gamma_e Sigma_e, I + Gamma_e Sigma_e Gamma_e']].
      ## The realized marginal mean is subtracted so E[e] = 0 (preserves SS).
      ## When Sigma_e is diagonal this reduces to the old per-shock independent
      ## draw e_i = sigma_i(delta_i|z_i| + sqrt(1-delta_i^2)w_i).
      sigma_e <- sqrt(diag(Sigma_e))          # length n_exo
      alpha_e <- .get_shock_skewness(model, exo, params)  # length n_exo (named)
      Gamma_e <- diag(alpha_e / sigma_e, nrow = n_exo)    # n_exo x n_exo
      SG_t      <- Sigma_e %*% t(Gamma_e)                 # Sigma_e Gamma_e'
      M_U       <- diag(n_exo) + Gamma_e %*% Sigma_e %*% t(Gamma_e)  # Cov(U)
      Joint_cov <- rbind(cbind(Sigma_e, SG_t),
                         cbind(t(SG_t), M_U))             # 2 n_exo x 2 n_exo
      Joint_cov <- 0.5 * (Joint_cov + t(Joint_cov))
      L_joint   <- chol(Joint_cov + diag(1e-12, 2 * n_exo))  # upper-tri factor
      .draw_joint_csn <- function(n_draws) {
        Zc   <- matrix(stats::rnorm(n_draws * 2 * n_exo), n_draws, 2 * n_exo) %*%
                  L_joint
        e_d  <- Zc[, seq_len(n_exo),         drop = FALSE]
        U_d  <- Zc[, n_exo + seq_len(n_exo), drop = FALSE]
        keep <- apply(U_d, 1L, function(u) all(u >= 0))
        e_d[keep, , drop = FALSE]
      }
      ## Mean shift: use the SAME closed-form correction the likelihood applies
      ## in .get_csn_shock_params (mu_eta = -RR (sigma_i delta_i sqrt(2/pi))).
      ## This keeps the DGP data level EXACTLY consistent with where the filter
      ## places its steady state, regardless of rho.  (The realized joint-CSN
      ## marginal mean differs slightly from this closed form at rho != 0, but
      ## both DGP and likelihood use the identical closed-form shift, so the
      ## level is consistent -- which is what the SBC rank distribution needs.)
      delta_e <- alpha_e / sqrt(1 + alpha_e^2)
      mu_e    <- sigma_e * delta_e * sqrt(2 / pi)   # length n_exo

      for (t in seq_len(n_total)) {
        e_one <- .draw_joint_csn(64L)
        while (nrow(e_one) < 1L) e_one <- .draw_joint_csn(64L)
        e_t   <- as.numeric(e_one[1L, ]) - mu_e   # zero-mean joint CSN draw

        if (t > presample) {
          y_t <- as.numeric(ZZ %*% s) + as.numeric(DD %*% e_t) + d_obs
          Y[, t - presample] <- y_t
        }
        s <- as.numeric(TT %*% s) + as.numeric(RR %*% e_t)
      }

    } else {
      for (t in seq_len(n_total)) {
        e_t <- as.numeric(L_e %*% rnorm(n_exo))
        if (t > presample) {
          ## y_t = Z s_{t-1} + D e_t + ybar
          y_t <- as.numeric(ZZ %*% s) + as.numeric(DD %*% e_t) + d_obs
          Y[, t - presample] <- y_t
        }
        ## s_t = T s_{t-1} + R e_t
        s <- as.numeric(TT %*% s) + as.numeric(RR %*% e_t)
      }
    }
  }

  ## ---- Step 3a: measurement error -----------------------------------------
  ## Every likelihood in the package treats `me_variance` as GENUINE i.i.d.
  ## observation noise (the multivariate Kalman filter since F3-D, 2026-09-03;
  ## the univariate KF, pruned, TPF, PSKF, OBC and MS filters before that), so
  ## the DGP must carry the same noise or the SBC tests model mismatch rather
  ## than calibration. Before this step the harness simulated noise-free data
  ## and scored it with a noisy likelihood: at me_variance = 1e-3 on the
  ## Ireland (2004) NK model the shock-sd ranks had mean-rank z = +6
  ## (posterior below the truth, the signature of variance attributed to
  ## measurement error that the data never contained).
  Y <- .sbc_add_me(Y, me_variance)

  ## ---- Step 3b (optional): per-draw KF innovation whiteness fast-fail ----
  ## Y was simulated from exactly this (dr, params) state-space -- this is
  ## the correctly-specified case for kf_innovation_diagnostics(), and the
  ## unit-root guard above (is_unit_root) already establishes lik_init =
  ## "stationary" is valid, so no further routing check is needed here.
  ## Failures caught internally (never propagated as a replication failure):
  ## this is a diagnostic add-on, not a gate on SBC itself.
  innovation_diagnostics <- NULL
  if (isTRUE(innovation_check)) {
    innovation_diagnostics <- tryCatch(
      kf_innovation_diagnostics(Y, dr = dr, model = model, params = params,
                                obs_vars = obs_vars, lik_init = "stationary",
                                me_variance = me_variance),
      error = function(e) NULL
    )
  }

  ## ---- Step 4: posterior --------------------------------------------------
  ## Additional args (e.g. n_particles for TPF) are forwarded via ...
  log_post_fn <- make_log_posterior(model, Y, prior_spec, obs_vars,
                                     compiled, me_variance = me_variance,
                                     likelihood = likelihood,
                                     lik_init = lik_init,
                                     ...)

  theta0 <- theta_tilde[prior_spec$name]
  names(theta0) <- prior_spec$name
  d      <- length(theta0)  ## parameter dimension (used by DIME walker count)

  lp0 <- log_post_fn(theta0)
  if (!is.finite(lp0$logpost))
    return(list(ok = FALSE, reason = "log-posterior at theta_tilde is non-finite"))

  prior_sd <- prior_spec$std
  names(prior_sd) <- prior_spec$name
  prior_sd <- prior_sd[names(theta0)]
  Sigma_prop <- diag(prior_sd^2, nrow = length(theta0))
  rownames(Sigma_prop) <- colnames(Sigma_prop) <- names(theta0)

  transform <- if (isTRUE(transform_params))
    build_param_transform(prior_spec, names(theta0))
  else
    NULL

  if (!is.null(transform)) {
    ## Sigma_prop must be an ETA-SPACE covariance for the transformed sampler.
    Sigma_prop <- .cov_theta_to_eta(Sigma_prop, transform, theta0)
  }

  ## Explicit dispatch on sampler: adding "nuts" to match.arg without the
  ## branch here would be a silent no-op (Landmine 5).
  if (sampler == "rwmh") {
    mcmc <- rwmh(log_post_fn, theta0, Sigma_prop,
                 n_draws = n_draws, n_burn = n_burn,
                 verbose = FALSE, transform = transform,
                 adapt_cov = isTRUE(adapt_cov))

    draws <- mcmc$chain  ## post-burn-in rows
    ranks <- .sbc_ranks(theta_tilde[names(theta0)], draws, thin = thin)
    L_effective <- floor(nrow(draws) / thin)

  } else if (sampler == "nuts") {
    ## Use prior-scale diagonal mass (1/var) as a reasonable starting point.
    ## Windowed warmup will overwrite this during adaptation.
    mass_diag_init <- 1 / pmax(diag(Sigma_prop), 1e-12)

    ## For NUTS + Whittle: build the analytic Whittle gradient inside the
    ## replication so parallel daemons need no API-surface changes.  The
    ## grad_fn captures (model, Y, prior_spec, obs_vars, compiled) which are
    ## all present in this replication's scope.  Option (b) from sbc-refresh
    ## brief S2: no change to dynhr_sbc() signature needed.
    ## NOTE: any future multi-obs Whittle battery cell must use obs_vars with
    ## >= 2 entries AND correlated shocks to exercise the complex Hermitian path.
    nuts_grad_fn <- if (identical(likelihood, "whittle")) {
      make_posterior_grad(
        model       = model,
        data        = Y,
        prior_spec  = prior_spec,
        obs_vars    = obs_vars,
        compiled    = compiled,
        me_variance = me_variance,
        likelihood  = "whittle",
        freq_band   = list(...)$freq_band %||% c(0, pi)
      )
    } else NULL

    mcmc <- dynhr_nuts(
      log_post_fn, theta0,
      n_draws  = n_draws,
      n_warmup = n_burn,
      adapt_mass = TRUE,
      mass_diag  = mass_diag_init,
      transform  = transform,
      grad_fn    = nuts_grad_fn,
      verbose    = FALSE
    )

    draws <- mcmc$chain  ## post-warmup rows
    ranks <- .sbc_ranks(theta_tilde[names(theta0)], draws, thin = thin)
    L_effective <- floor(nrow(draws) / thin)

  } else if (sampler == "smc") {
    ## NOTE: for SMC, n_draws is repurposed as n_particles. n_burn is ignored
    ## (SMC has no burn-in phase; it tempers from the prior). thin is used only
    ## to compute L_target = floor(n_draws / thin).
    ## L_target approximate i.i.d. draws are obtained via systematic
    ## resampling of the weighted particle cloud.
    L_target <- max(1L, floor(n_draws / thin))

    mcmc <- dynhr_smc(
      log_post_fn,
      prior_sampler = prior_sampler,
      n_particles   = as.integer(n_draws),
      verbose       = FALSE
    )

    adapter_result <- sbc_draws_from_smc(mcmc, L_target)
    if (!isTRUE(adapter_result$ok))
      return(list(ok = FALSE, reason = adapter_result$reason))

    draws       <- adapter_result$draws
    ranks       <- .sbc_ranks(theta_tilde[names(theta0)], draws, thin = 1L)
    L_effective <- nrow(draws)

  } else if (sampler == "dime") {
    ## For DIME, n_draws is used as n_iter (total post-burn iterations per
    ## walker) and n_burn is passed directly. The number of walkers defaults
    ## to max(5*d, 20). L_target approximate i.i.d. draws are obtained by
    ## per-walker thinning + concatenation.
    L_target <- max(1L, floor((n_draws - n_burn) / thin))

    mcmc <- run_dime(
      log_post_fn,
      prior_spec    = prior_spec,
      prior_sampler = prior_sampler,
      n_chain       = max(5L * d, 20L),
      n_iter        = as.integer(n_draws),
      n_burn        = as.integer(n_burn),
      verbose       = FALSE
    )

    draws       <- sbc_draws_from_dime(mcmc, L_target)
    ranks       <- .sbc_ranks(theta_tilde[names(theta0)], draws, thin = 1L)
    L_effective <- nrow(draws)

  } else if (sampler %in% c("hmc", "mala", "chees")) {
    ## Gradient-based samplers (HMC / MALA / ChEES-HMC). Build the analytic
    ## gradient for the likelihoods that expose one (gaussian / whittle /
    ## cumulant -- see .ctx_allows_analytic_gradient); the dispatch guard in
    ## dynhr_sbc() rejects tpf / pskf / student_t for these samplers, so a usable
    ## grad_fn is always available on the reachable paths. n_draws / n_burn map
    ## to post-warmup draws / warmup, as for nuts.
    grad_fn <- if (likelihood %in% c("gaussian", "whittle", "cumulant")) {
      make_posterior_grad(
        model      = model, data = Y, prior_spec = prior_spec,
        obs_vars   = obs_vars, compiled = compiled, me_variance = me_variance,
        likelihood = likelihood,
        freq_band  = list(...)$freq_band %||% c(0, pi))
    } else NULL

    mcmc <- if (sampler == "hmc") {
      dynhr_hmc(log_post_fn, theta0, n_draws = n_draws, n_warmup = n_burn,
                grad_fn = grad_fn, transform = transform, verbose = FALSE)
    } else if (sampler == "mala") {
      dynhr_mala(log_post_fn, theta0, n_draws = n_draws, n_warmup = n_burn,
                 grad_fn = grad_fn, transform = transform, verbose = FALSE)
    } else {
      dynhr_chees(log_post_fn, theta0, n_draws = n_draws, n_warmup = n_burn,
                  mass_diag = 1 / pmax(diag(Sigma_prop), 1e-12),
                  grad_fn = grad_fn, transform = transform, verbose = FALSE)
    }

    draws       <- mcmc$chain  ## post-warmup rows
    ranks       <- .sbc_ranks(theta_tilde[names(theta0)], draws, thin = thin)
    L_effective <- floor(nrow(draws) / thin)

  } else {
    stop(".sbc_one_replication: unknown sampler '", sampler, "'")
  }

  list(ok = TRUE, ranks = ranks, L_effective = L_effective,
       innovation_diagnostics = innovation_diagnostics)
}


#' Simulation-Based Calibration for the dynhr estimation pipeline
#'
#' Runs Simulation-Based Calibration (SBC; Talts, Betancourt, Simpson,
#' Vehtari & Gelman, 2018) on the full prior -> simulate -> estimate
#' pipeline for \code{model}: for each of \code{n_replications} independent
#' replications, a parameter vector \eqn{\tilde\theta} is drawn from the
#' prior, data \eqn{Y} is simulated from the model's linearised state space at
#' \eqn{\tilde\theta}, and \code{rwmh()} is run on the resulting posterior.
#' The rank of \eqn{\tilde\theta_j} among the (thinned) posterior draws of
#' \eqn{\theta_j} is recorded for each parameter \eqn{j}; under correct
#' calibration these ranks are uniformly distributed.
#'
#' \strong{Conditioning on the prior solving region.} Replications whose
#' \eqn{\tilde\theta} draw fails to produce a converged steady state or a
#' Blanchard-Kahn-satisfying perturbation solution are skipped (counted in
#' \code{n_failed}). The resulting calibration check is therefore
#' conditional on the region of the prior where the model solves -- this is
#' the practical scope of SBC for nonlinear DSGE models, and is reported
#' explicitly via \code{n_failed} / \code{n_replications}.
#'
#' \strong{Unit-root models.} If the model's state-transition matrix has an
#' eigenvalue with modulus \code{>= 1 - 1e-6} (i.e.
#' \code{make_log_posterior(..., lik_init = "auto")} would resolve to the
#' exact diffuse Kalman filter), \code{dynhr_sbc()} stops with an error: the
#' diffuse prior used by the filter for unit-root models does not correspond
#' to any proper data-generating distribution for \eqn{s_0}, so there is no
#' well-defined simulation procedure to validate against. SBC is therefore
#' restricted to stationary models.
#'
#' \strong{Pre-sample burn-in.} Because the Kalman filter initialises
#' \eqn{P_0} from the model's stationary distribution
#' (\code{lik_init = "stationary"}), the simulator draws \code{presample}
#' extra periods (default 50) before \eqn{t = 1} so that \eqn{s_0} (the
#' state entering period 1) is itself approximately a draw from the
#' stationary distribution, matching the filter's assumption. Starting from
#' \eqn{s_0 = 0} with no burn-in would create a simulator/filter mismatch
#' that SBC would (correctly) flag as miscalibration.
#'
#' \strong{Serial execution.} Replications are run serially via the internal
#' helper \code{.sbc_one_replication()}. That helper is self-contained given
#' \code{model}/\code{compiled}/\code{prior_spec}, so dispatching
#' replications over \code{mirai} workers is a future one-line change (no
#' refactor of the per-replication logic required).
#'
#' @param model            A parsed \code{dynhr_mod} (from
#'   \code{\link{parse_mod}}) carrying an \code{estimated_params} block.
#' @param obs_vars         Character vector of observed variable names.
#' @param T_obs            Number of simulated observation periods per
#'   replication (default 100).
#' @param n_replications   Number of SBC replications (default 50).
#' @param n_draws          Total RWMH draws per replication, including
#'   burn-in (default 3000).
#' @param n_burn           RWMH burn-in draws to discard (default 1500).
#' @param thin             Thinning interval applied to the retained draws
#'   before computing ranks (default 10).
#' @param sampler          Sampler to use. Gradient-free (work with any
#'   likelihood): \code{"rwmh"} (default), \code{"smc"}, \code{"dime"}.
#'   Gradient-based (require an analytic gradient, available only for the
#'   \code{"gaussian"}, \code{"whittle"} and \code{"cumulant"} likelihoods):
#'   \code{"nuts"}, \code{"hmc"}, \code{"mala"}, \code{"chees"}.
#'   The gradient-based samplers use windowed mass-matrix adaptation with
#'   \code{n_burn} warmup iterations. For \code{"smc"}, \code{n_draws} is
#'   repurposed as \code{n_particles} and \code{n_burn} / \code{thin} are used
#'   only to derive \code{L_target = floor(n_draws / thin)}; the ESS floor
#'   \code{>= 0.3 * n_particles} is enforced. For \code{"dime"}, \code{n_draws}
#'   is used as \code{n_iter} (post-burn iterations per walker). A gradient-based
#'   sampler with a non-differentiable likelihood (\code{"tpf"}, \code{"pskf"},
#'   \code{"student_t"}) is unsupported and raises an error.
#' @param me_variance      Measurement-error variance, forwarded to
#'   \code{make_log_posterior} (default 0).
#' @param lik_init         Kalman filter \code{P0} initialisation, forwarded
#'   to \code{make_log_posterior} (default \code{"auto"}).
#' @param transform_params Logical (default \code{FALSE}): if \code{TRUE},
#'   run RWMH in unconstrained eta-space via
#'   \code{build_param_transform}.
#' @param adapt_cov        Logical (default \code{FALSE}): if \code{TRUE},
#'   use Haario et al. (2001) adaptive proposal covariance in \code{rwmh()}.
#' @param seed             Base RNG seed (default 1L). Replication \code{i}
#'   uses seed \code{seed + i} for reproducibility.
#' @param n_cores          Number of parallel workers when running replications
#'   in parallel (\code{NULL} = auto-detect).  Currently reserved for a future
#'   mirai-based parallel path; replications are run serially in this version.
#' @param likelihood       Likelihood type forwarded to
#'   \code{make_log_posterior}: \code{"gaussian"} (default), \code{"cumulant"},
#'   \code{"whittle"}, \code{"tpf"}, \code{"pskf"}, or \code{"student_t"}
#'   (\code{"student_t"} also needs \code{student_df} supplied via \code{...}).
#' @param order            Perturbation order used when simulating the data
#'   generating process (default \code{1L}).  Pass \code{2L} when
#'   \code{likelihood = "tpf"} so the DGP and likelihood are at the same
#'   approximation order.
#' @param ctx              Optional \code{dynhr_estimation_context} object
#'   (from \code{\link{estimation_context}}).  When supplied it overrides the
#'   individual \code{likelihood}, \code{me_variance}, and \code{lik_init}
#'   arguments and merges any \code{tpf_options} with \code{...}.
#' @param innovation_check Logical (default \code{FALSE}): if \code{TRUE},
#'   run \code{\link{kf_innovation_diagnostics}} on each replication's
#'   simulated data \eqn{Y} at the DRAWN \eqn{\tilde\theta} (the replication's
#'   own DGP parameters -- the correctly-specified case, since \eqn{Y} was
#'   simulated from exactly this \code{dr}/\code{params}) as a cheap per-draw
#'   fast-fail sanity oracle: a whiteness failure here would mean the
#'   simulator and the thin diagnostic filter disagree about the model's own
#'   state-space form, which is a bug independent of anything SBC itself is
#'   testing. Only applies on the serial path (\code{n_cores = NULL}); ignored
#'   with a message when \code{n_cores} is set, since \code{run_sbc_mirai}
#'   is unaffected by this argument. Purely additive: when \code{FALSE}
#'   (default), \code{.sbc_one_replication()} takes an identical code path to
#'   before this argument existed, so \code{ranks}/\code{uniformity}/\code{plot}
#'   are bit-identical. When \code{TRUE}, adds
#'   \code{$innovation_diagnostics} (list of per-replication results, one per
#'   successful replication) and \code{$innovation_summary} (a data.frame with
#'   \code{n_flagged_replications} / \code{n_replications_checked}) to the
#'   returned \code{dynhr_sbc} object.
#' @param verbose          Print progress messages (default \code{TRUE}).
#' @param ...              Additional arguments forwarded to
#'   \code{make_log_posterior} and the sampler
#'   (e.g. \code{n_particles} for \code{likelihood = "tpf"}).
#' @return An object of class \code{"dynhr_sbc"}: a list with elements
#'   \describe{
#'     \item{ranks}{\code{n_ok x d} integer matrix of SBC ranks (one row per
#'       successful replication).}
#'     \item{n_failed}{Number of replications skipped because the model did
#'       not solve at the prior draw.}
#'     \item{n_replications}{The requested number of replications.}
#'     \item{L_effective}{Number of thinned posterior draws per
#'       replication.}
#'     \item{uniformity}{Result of \code{sbc_uniformity_test(ranks)}.}
#'     \item{settings}{List of the call's settings, for provenance.}
#'     \item{plot}{A ggplot rank-histogram (see
#'       \code{sbc_uniformity_test}), or \code{NULL} if \pkg{ggplot2}
#'       is unavailable or there are no successful replications.}
#'     \item{innovation_diagnostics}{Only present when
#'       \code{innovation_check = TRUE}: list of per-replication
#'       \code{kf_innovation_diagnostics} objects (successful replications
#'       only, in replication order).}
#'     \item{innovation_summary}{Only present when
#'       \code{innovation_check = TRUE}: single-row data.frame with
#'       \code{n_flagged_replications} (replications with at least one
#'       \code{|z| > 4}) and \code{n_replications_checked}.}
#'   }
#' @references
#' Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
#'   Validating Bayesian Inference Algorithms with Simulation-Based
#'   Calibration. \emph{arXiv preprint arXiv:1804.06788}.
#'
#' Geweke, J. (2004). Getting It Right: Joint Distribution Tests of Posterior
#'   Simulators. \emph{Journal of the American Statistical Association},
#'   99(467), 799-804. (An alternative joint-distribution calibration check.)
#' @export
dynhr_sbc <- function(model, obs_vars, T_obs = 100L, n_replications = 50L,
                      n_draws = 3000L, n_burn = 1500L, thin = 10L,
                      sampler = c("rwmh", "nuts", "smc", "dime",
                                  "hmc", "mala", "chees"),
                      me_variance = 0,
                      lik_init = "auto", transform_params = FALSE,
                      adapt_cov = FALSE, seed = 1L, verbose = TRUE,
                      n_cores = NULL,
                      likelihood = c("gaussian", "cumulant", "whittle", "tpf",
                                     "pskf", "student_t"),
                      order = 1L,
                      ctx = NULL,
                      innovation_check = FALSE,
                      ...) {
  ## Unpack ctx when provided -- overrides individual args.
  ## tpf_options from ctx are merged into ... for the .sbc_one_replication call.
  if (!is.null(ctx) && inherits(ctx, "dynhr_estimation_context")) {
    me_variance <- ctx$me_variance
    lik_init    <- ctx$lik_init
    likelihood  <- ctx$likelihood
    ## tpf_options: merge into dots (ctx entries first, ... can override)
    if (length(ctx$tpf_options) > 0L) {
      existing_dots <- list(...)
      merged_tpf <- modifyList(ctx$tpf_options, existing_dots)
      ## We can't modify ... directly; store for later use
      .ctx_tpf_extra <- merged_tpf
    }
  }
  likelihood <- match.arg(likelihood)

  sampler <- match.arg(sampler)

  ## Guard unsupported sampler x likelihood combinations. Gradient-based samplers
  ## (nuts/hmc/mala/chees) require an ANALYTIC gradient; only gaussian, whittle and
  ## cumulant expose one (see .ctx_allows_analytic_gradient). tpf is stochastic and
  ## non-differentiable, pskf has fragile CDF-based finite differences, and student_t
  ## has no adjoint path -- finite-difference gradients of these are noisy/unreliable
  ## (O(sigma_noise / h) error), so the sampler's trajectories diverge. Use a
  ## gradient-free sampler (rwmh, smc, or dime) for those likelihoods.
  if (sampler %in% c("nuts", "hmc", "mala", "chees") &&
      likelihood %in% c("tpf", "pskf", "student_t"))
    stop("dynhr_sbc: sampler = '", sampler, "' is unsupported with likelihood = '",
         likelihood, "'. Gradient-based samplers (nuts/hmc/mala/chees) require an ",
         "analytic gradient; '", likelihood, "' exposes none (noisy/fragile finite ",
         "differences). Use a gradient-free sampler (rwmh, smc, or dime).")

  if (is.null(model$estimated_params) || nrow(model$estimated_params) == 0)
    stop("dynhr_sbc: 'model' has no estimated_params block. ",
         "SBC requires a prior over the parameters being validated -- ",
         "add an estimated_params block to the .mod file.")

  prior_spec    <- extract_prior_spec(model, verbose = FALSE)
  prior_sampler <- .smc_make_prior_sampler(prior_spec)
  par_names     <- prior_spec$name
  d             <- length(par_names)

  presample <- 50L

  n_replications <- as.integer(n_replications)

  if (verbose)
    cat(sprintf("Running SBC: %d replications, T_obs = %d, %d parameter(s)\n",
                n_replications, T_obs, d))

  ## -- parallel path -------------------------------------------------------
  if (!is.null(n_cores)) {
    if (isTRUE(innovation_check) && verbose)
      cat("Note: innovation_check is only supported on the serial path ",
          "(n_cores = NULL); ignored for this parallel run.\n", sep = "")
    raw_list <- run_sbc_mirai(
      model            = model,
      prior_spec       = prior_spec,
      prior_sampler    = prior_sampler,
      obs_vars         = obs_vars,
      T_obs            = T_obs,
      presample        = presample,
      n_replications   = n_replications,
      n_draws          = n_draws,
      n_burn           = n_burn,
      thin             = thin,
      sampler          = sampler,
      me_variance      = me_variance,
      lik_init         = lik_init,
      transform_params = transform_params,
      adapt_cov        = adapt_cov,
      seed_base        = seed,
      n_cores          = as.integer(n_cores),
      verbose          = verbose,
      likelihood       = likelihood,
      order            = order
    )
    ranks_list  <- vector("list", n_replications)
    n_failed    <- 0L
    L_effective <- NULL
    for (i in seq_len(n_replications)) {
      res <- raw_list[[i]]
      if (isTRUE(res$ok)) {
        ranks_list[[i]] <- res$ranks
        if (is.null(L_effective)) L_effective <- res$L_effective
        if (verbose) cat(sprintf("  [%d/%d] ok\n", i, n_replications))
      } else {
        n_failed <- n_failed + 1L
        if (verbose)
          cat(sprintf("  [%d/%d] failed: %s\n", i, n_replications, res$reason))
      }
    }
    ## jump to assembly below
  } else {
    ## -- serial path (n_cores = NULL) ----------------------------------------
    compiled <- compile_model(model, verbose = FALSE)

    ranks_list  <- vector("list", n_replications)
    n_failed    <- 0L
    L_effective <- NULL
    innov_list  <- if (isTRUE(innovation_check)) vector("list", n_replications) else NULL

    ## Merge ctx$tpf_options (if set) into the dots for .sbc_one_replication.
    .sbc_extra_args <- if (exists(".ctx_tpf_extra", inherits = FALSE)) {
      modifyList(.ctx_tpf_extra, list(...))
    } else {
      list(...)
    }

    for (i in seq_len(n_replications)) {
      rep_seed <- seed + i
      res <- do.call(.sbc_one_replication, c(
        list(
          i = i, model = model, compiled = compiled, prior_spec = prior_spec,
          prior_sampler = prior_sampler, obs_vars = obs_vars, T_obs = T_obs,
          presample = presample, n_draws = n_draws, n_burn = n_burn, thin = thin,
          sampler = sampler, me_variance = me_variance, lik_init = lik_init,
          transform_params = transform_params, adapt_cov = adapt_cov,
          likelihood = likelihood, order = order,
          seed = rep_seed, verbose = verbose,
          innovation_check = isTRUE(innovation_check)
        ),
        .sbc_extra_args
      ))

      if (isTRUE(res$ok)) {
        ranks_list[[i]] <- res$ranks
        if (is.null(L_effective)) L_effective <- res$L_effective
        if (isTRUE(innovation_check)) innov_list[[i]] <- res$innovation_diagnostics
        if (verbose)
          cat(sprintf("  [%d/%d] ok\n", i, n_replications))
      } else {
        n_failed <- n_failed + 1L
        if (verbose)
          cat(sprintf("  [%d/%d] failed: %s\n", i, n_replications, res$reason))
      }
    }
  }

  ok <- !vapply(ranks_list, is.null, logical(1))
  ranks <- if (any(ok)) {
    do.call(rbind, ranks_list[ok])
  } else {
    matrix(integer(0), nrow = 0, ncol = d, dimnames = list(NULL, par_names))
  }
  rownames(ranks) <- NULL

  ## L_effective was already captured above from the first successful
  ## replication's res$L_effective (the harness's own known rank support) --
  ## pass it through explicitly rather than letting sbc_uniformity_test fall
  ## back to inferring it from the observed ranks (#7 one-liner, adversarial
  ## review).
  uniformity <- if (nrow(ranks) > 0) {
    sbc_uniformity_test(ranks, L = L_effective)
  } else {
    NULL
  }

  plot_obj <- if (nrow(ranks) > 0) .sbc_rank_hist_plot(ranks) else NULL

  settings <- list(
    obs_vars = obs_vars, T_obs = T_obs, n_replications = n_replications,
    n_draws = n_draws, n_burn = n_burn, thin = thin, sampler = sampler,
    me_variance = me_variance, lik_init = lik_init,
    transform_params = transform_params, adapt_cov = adapt_cov,
    seed = seed, presample = presample, innovation_check = isTRUE(innovation_check)
  )

  out <- list(
    ranks          = ranks,
    n_failed       = n_failed,
    n_replications = n_replications,
    L_effective    = L_effective %||% NA_integer_,
    uniformity     = uniformity,
    settings       = settings,
    plot           = plot_obj
  )

  ## Additive only: these two fields exist ONLY when innovation_check = TRUE,
  ## so a default (FALSE) call's result list has the exact same names/values
  ## as before this argument was added.
  if (isTRUE(innovation_check) && exists("innov_list", inherits = FALSE)) {
    innov_ok  <- Filter(Negate(is.null), innov_list)
    n_flagged <- sum(vapply(innov_ok, function(d) isTRUE(d$joint$n_flagged > 0),
                            logical(1)))
    out$innovation_diagnostics <- innov_ok
    out$innovation_summary <- data.frame(
      n_flagged_replications  = n_flagged,
      n_replications_checked  = length(innov_ok),
      stringsAsFactors = FALSE
    )
  }

  structure(out, class = "dynhr_sbc")
}


## ===========================================================================
## LAYER 4: SBC coverage matrix
## ===========================================================================

#' Build the SBC sampler x likelihood coverage matrix and write it to disk
#'
#' Constructs a \code{data.frame} summarising the certification status of
#' every (sampler, likelihood) cell and writes it to
#' \code{inst/extdata/sbc_matrix.rds} so that the companion vignette
#' \code{vignettes/sbc-matrix.Rmd} can render a published table.
#'
#' \strong{Status vocabulary:}
#' \describe{
#'   \item{\code{"certified"}}{Passed the battery uniformity gate
#'     (all per-parameter Bonferroni-corrected p-values \eqn{> 0.05}) for
#'     at least one model.}
#'   \item{\code{"to-certify"}}{Sampler is implemented and the dispatch
#'     wired up, but the battery cell has not yet run.}
#'   \item{\code{"characterized"}}{Runs, but the likelihood is approximate
#'     (e.g. Whittle at small T); p-values stored as provenance but not
#'     gated.}
#'   \item{\code{"incompatible"}}{Mathematically impossible: a gradient-based
#'     sampler on a STOCHASTIC, non-differentiable likelihood (tpf) -- no usable
#'     gradient exists.}
#'   \item{\code{"needs-gradient"}}{A gradient-based sampler on a deterministic
#'     likelihood that is differentiable in principle but has no analytic-gradient
#'     path implemented yet (pskf, student_t) -- possible, not yet done.}
#'   \item{\code{"pending"}}{Likelihood not yet implemented in dynhr.}
#' }
#'
#' @param sbc_results Named list of \code{dynhr_sbc} results, keyed by
#'   \code{"<sampler>_<likelihood>"} (e.g. \code{"rwmh_gaussian"}).
#'   Pass \code{NULL} to build the static status table only (no provenance
#'   from actual battery runs).
#' @param out_path   Path to write the RDS.  Defaults to
#'   \code{inst/extdata/sbc_matrix.rds} relative to the working directory.
#'   Set to \code{NULL} to suppress writing.
#' @return Invisibly, a \code{data.frame} with columns
#'   \code{sampler, likelihood, status, notes, verdict, p_min}.
#' @noRd
sbc_matrix_result <- function(sbc_results = NULL, out_path = "inst/extdata/sbc_matrix.rds") {

  ## ---- Static status table: 7 samplers x 6 likelihoods --------------------
  ## Built from the support RULE -- gradient-based samplers need an analytic
  ## gradient, which only gaussian/whittle/cumulant expose (see
  ## .ctx_allows_analytic_gradient) -- plus bespoke overrides carrying the
  ## provenance from the 2026-06-13 local battery runs.
  samplers <- c("rwmh", "smc", "dime",            # gradient-free: any likelihood
                "nuts", "hmc", "mala", "chees")   # gradient-based
  liks      <- c("gaussian_kf", "tpf", "whittle", "pskf", "cumulant", "student_t")
  grad_samplers <- c("nuts", "hmc", "mala", "chees")        # require analytic gradient
  grad_liks     <- c("gaussian_kf", "whittle", "cumulant")  # expose an analytic gradient

  ## Grouped by likelihood (7 sampler rows each), matching the wide vignette table.
  static <- expand.grid(sampler = samplers, likelihood = liks,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  ## Default: a gradient sampler on a likelihood with no analytic gradient is
  ## mathematically unsupported (noisy/fragile finite differences); every other
  ## (implemented) cell is wired and awaiting a battery run.
  is_no_grad <- static$sampler %in% grad_samplers &
                !(static$likelihood %in% grad_liks)
  ## Split the "no analytic gradient" cells by WHY: tpf is mathematically
  ## INCOMPATIBLE (a stochastic, non-differentiable particle estimate -- no usable
  ## gradient exists), whereas pskf and student_t are deterministic + differentiable
  ## in principle but simply have no analytic-gradient path implemented yet.
  static$status <- ifelse(
    !is_no_grad, "to-certify",
    ifelse(static$likelihood == "tpf", "incompatible", "needs-gradient"))
  static$notes <- "dispatch wired; battery cell not yet run"
  static$notes[is_no_grad & static$likelihood == "tpf"] <-
    "gradient sampler x a STOCHASTIC, non-differentiable particle likelihood -- no usable gradient exists (mathematically incompatible)"
  static$notes[is_no_grad & static$likelihood == "pskf"] <-
    "no analytic gradient implemented for pskf (CDF/orthant-based, fragile FD) -- derivable in principle, not yet coded"
  static$notes[is_no_grad & static$likelihood == "student_t"] <-
    "no analytic gradient implemented for student_t (= the gaussian KF adjoint + t-density derivative) -- straightforward to add"

  ## Bespoke overrides (provenance from actual runs / inherent approximation).
  ## Functional setter: returns the modified frame (no superassignment).
  .set <- function(df, s, l, status, notes) {
    i <- which(df$sampler == s & df$likelihood == l)
    df$status[i] <- status
    df$notes[i]  <- notes
    df
  }
  static <- .set(static, "rwmh", "gaussian_kf", "characterized",
       "3-model battery; 2026-06-13 local run: AR1 calibrated, Ireland + FS2000 miscalibrated (pre-existing; investigation pending)")
  static <- .set(static, "nuts", "gaussian_kf", "certified", "AR1 battery; Tier A weekly CI")
  ## gaussian_kf certifications -- AR(1) battery 2026-06-27 (T=50, 50 reps,
  ## seed=42, n_draws=1200 / n_burn=600 / thin=4). "calibrated" = every
  ## per-parameter Bonferroni rank p-value > 0.025, with 0/50 replications
  ## failing to solve at the prior draw.
  static <- .set(static, "smc",  "gaussian_kf", "certified",
       "AR(1) battery 2026-06-27: calibrated, min p=0.494, 0/50 failed (50 reps, seed 42)")
  static <- .set(static, "dime", "gaussian_kf", "certified",
       "AR(1) battery 2026-06-27: calibrated, min p=0.262, 0/50 failed (50 reps, seed 42)")
  static <- .set(static, "mala", "gaussian_kf", "certified",
       "AR(1) battery 2026-06-27: calibrated, min p=0.616, 0/50 failed (50 reps, seed 42)")
  static <- .set(static, "hmc",  "gaussian_kf", "certified",
       "AR(1) battery 2026-06-27 (POST dual-averaging fix, commit 2acbd67): calibrated, min p=0.290, 0/50 failed -- was 8.5e-48 (frozen chain) pre-fix (50 reps, seed 42)")
  ## chees stays to-certify: the AR(1) cert was attempted but abandoned (~3h),
  ## and controlled experiments show a slow-mixing weakness on the ill-conditioned
  ## (near-unit-root) replications that this fixture stresses. See the GAP note.
  static <- .set(static, "chees", "gaussian_kf", "to-certify",
       "AR(1) 50-rep cert attempted 2026-06-27 but KILLED at ~3h -- the cost is the per-eval Kalman/diffuse-filter on near-unit-root prior draws (~= hmc 2.6h), NOT the sampler. Controlled synthetic Gaussians instead show chees is SLOW-MIXING (converges to truth, NOT biased) on ill-conditioned correlated posteriors: its ChEES trajectory adaptation under-shoots (L~5 vs optimal ~sqrt(kappa)=32 at kappa=1e3; big-axis variance 0.58/0.69/0.77 at n=1.5k/6k/24k vs truth 1.0; NUTS reaches 2% in 6k). GAP: the AR(1) battery is both too slow AND stresses exactly chees's weak region, so a clean cert needs a cheaper / better-conditioned fixture or tuned chees adaptation; prefer NUTS for near-unit-root/stiff posteriors. Evidence: .claude/orchestration/chees-slow-mixing-2026-06-27.md")
  static <- .set(static, "rwmh", "tpf", "characterized",
       "AR1 order=2 battery; 2026-06-13 local run: sig p = 1e-4 (rho calibrated); pre-existing; N-sweep DONE (bias -1.50->-0.24 nat over N=250..5000; N=5000 intractable ~10h; characterized as finite-N variance bias, not an implementation bug)")
  ## whittle is asymptotically approximate -> never gated; the four originally-run
  ## cells are characterized, the newly-wired gradient cells stay to-certify.
  .wn <- "Whittle is asymptotically approximate; p-values stored as provenance only; 1-obs AR1 cell uses exact complex path"
  static <- .set(static, "rwmh", "whittle", "characterized", .wn)
  static <- .set(static, "smc",  "whittle", "characterized", .wn)
  static <- .set(static, "dime", "whittle", "characterized", .wn)
  static <- .set(static, "nuts", "whittle", "characterized",
       "AR1 NUTS + analytic whittle gradient; T=200; seed=801; debiased Whittle still miscalibrated (min p < 1e-4) -- posterior width, not location; characterized like the other whittle rows")
  for (s in c("hmc", "mala", "chees"))
    static <- .set(static, s, "whittle", "to-certify",
         "analytic whittle gradient wired; battery cell not yet run (whittle approximate -> characterized, not gated, once run)")
  static <- .set(static, "rwmh", "pskf", "to-certify",
       "AR1 estimated-alpha battery cell; 2026-06-13 local run miscalibrated (pre-existing); investigation pending")
  ## cumulant exposes an analytic gradient (cumulant_loglik_grad), so the
  ## gradient-based samplers are SUPPORTED here (not unsupported like tpf/pskf).
  for (s in grad_samplers)
    static <- .set(static, s, "cumulant", "to-certify",
         "cumulant has an analytic gradient (cumulant_loglik_grad); dispatch wired; battery cell not yet run")

  static$verdict <- NA_character_
  static$p_min   <- NA_real_

  ## Overlay verdicts from any supplied battery results
  if (!is.null(sbc_results)) {
    for (key in names(sbc_results)) {
      res <- sbc_results[[key]]
      parts <- strsplit(key, "_", fixed = TRUE)[[1L]]
      if (length(parts) < 2L) next
      samp <- parts[1L]
      lik  <- paste(parts[-1L], collapse = "_")
      idx  <- which(static$sampler == samp & static$likelihood == lik)
      if (length(idx) == 0L) next
      if (!is.null(res$uniformity)) {
        static$verdict[idx] <- res$uniformity$verdict
        static$p_min[idx]   <- min(res$uniformity$table$p_value, na.rm = TRUE)
      }
    }
  }

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(static, out_path)
    message("sbc_matrix_result: wrote ", out_path)
  }

  invisible(static)
}


#' @exportS3Method
#' @noRd
print.dynhr_sbc <- function(x, ...) {
  cat("Simulation-Based Calibration (Talts et al. 2018)\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("Replications: %d (failed: %d)\n", x$n_replications, x$n_failed))
  cat(sprintf("Posterior draws per replication (thinned): %s\n",
              if (is.na(x$L_effective)) "NA" else as.character(x$L_effective)))

  if (is.null(x$uniformity)) {
    cat("\nNo successful replications -- no ranks to report.\n")
    return(invisible(x))
  }

  print(x$uniformity)
  cat("(see x$plot for the rank histograms)\n")

  invisible(x)
}


## Add i.i.d. N(0, me_variance) measurement error to an n_obs x T panel.
## Identity when me_variance == 0 (no RNG draw, so seeds are unchanged).
#' @noRd
.sbc_add_me <- function(Y, me_variance) {
  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance < 0)
    stop(".sbc_add_me: `me_variance` must be a finite non-negative scalar.",
         call. = FALSE)
  if (me_variance == 0) return(Y)
  Y + matrix(stats::rnorm(length(Y), sd = sqrt(me_variance)),
             nrow = nrow(Y), ncol = ncol(Y))
}
