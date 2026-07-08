## R/hank-welfare-posterior.R
## --------------------------------------------------------------------------
## Mixture-level and posterior-level welfare STATEMENTS, built on top of the
## per-block welfare machinery in R/hank-welfare.R (hank_welfare_decompose(),
## hank_welfare_response(), hank_cev()) and the discount-factor mixture
## steady state in R/hank-ge.R (hank_mixture_ks_steady(), class
## "hank_mixture_ks": fields omega, blocks, betas).
##
## hank_welfare_decompose() already implements the Dávila-Schaab efficiency/
## redistribution split for a SINGLE population measure D. A mixture economy
## has one block (and one D) PER TYPE k, with population share omega_k; the
## survey (references/HANK_WELFARE_SURVEY.md, section 5.2) flags the missing
## piece as assembling the POOLED population measure and welfare vector
## across types before handing them to hank_welfare_decompose() -- that is
## hank_mixture_welfare_pool() below.
##
## hank_welfare_posterior() then answers a different question: given an
## ESTIMATED posterior over the model's structural parameters (a matrix of
## draws), what is the posterior distribution of a welfare statement (e.g.
## the pooled efficiency/redistribution split, or a per-type CEV), rather
## than the welfare statement's value at a single point calibration?
## --------------------------------------------------------------------------


#' Pool per-type welfare vectors of a discount-factor mixture into a single
#' population-level welfare decomposition
#'
#' Assembles the POOLED population measure and per-cell welfare-gain vector
#' across the types of a \code{\link{hank_mixture_ks_steady}} mixture economy,
#' and hands them to \code{\link{hank_welfare_decompose}} for the pooled
#' Dávila-Schaab efficiency/redistribution split. Each type \code{k} occupies
#' population share \code{mks$omega[k]} and has its own steady-state
#' population measure \code{mks$blocks[[k]]$D} (mass 1 WITHIN the type); the
#' pooled economy-wide measure re-weights each type's mass by its population
#' share, \eqn{D^{pool} = \mathrm{concat}_k(\omega_k\, D_k)}, which sums to 1
#' overall (since each \code{D_k} sums to 1 and \code{omega} sums to 1). The
#' per-cell welfare gains simply concatenate, \eqn{\lambda^{pool} =
#' \mathrm{concat}_k(\lambda_k)}, since CEV/welfare units do not depend on
#' population mass.
#'
#' This is exactly the "missing piece" flagged in
#' \code{references/HANK_WELFARE_SURVEY.md} section 5.2 for computing
#' \eqn{\Delta W = \sum_i D_i \lambda_i} at the pooled (mixture) level rather
#' than within a single type's block.
#'
#' @param mks A \code{\link{hank_mixture_ks_steady}} result (class
#'   \code{"hank_mixture_ks"}: fields \code{omega} (per-type population
#'   weights) and \code{blocks} (list of \code{hank_het_block}s, one per
#'   type, each carrying its own \code{$D})), OR a
#'   \code{\link{hank_nk_hank_mixture}} result (class
#'   \code{c("hank_nk_mixture", "hank_nk")}), which wraps the SAME
#'   \code{omega}/\code{blocks} fields directly (plus the assembled sticky-
#'   price \code{model} and other NK calibration fields, which are ignored
#'   here). Either shape is accepted transparently -- both simply need
#'   \code{$blocks} (list of \code{K} \code{hank_het_block}s) and
#'   \code{$omega} (length-\code{K} population weights).
#' @param lambda_list List of per-type per-cell CEV/welfare-gain vectors, one
#'   per \code{mks$blocks[[k]]}, each of length \code{n_e*n_a} for that
#'   block (e.g. \code{hank_welfare_response(block)$lambda},
#'   \code{hank_welfare_channels(block)$lambda}, or a per-cell
#'   \code{\link{hank_cev}}), in the SAME cell order as \code{block$D}.
#' @param weights Optional list of per-type per-cell planner weights, same
#'   layout as \code{lambda_list} (one vector per type, matching that type's
#'   cell count); concatenated in the same pooled order and passed on to
#'   \code{\link{hank_welfare_decompose}}. \code{NULL} (default) = uniform
#'   planner weights (redistribution component 0).
#'
#' @return A list with:
#'   \item{lambda_pool}{Numeric: the concatenated per-cell welfare gains
#'     across all types, pooled cell order (type-major, i.e. type 1's cells
#'     first, in \code{mks$blocks[[1]]}'s own cell order, then type 2's, ...).}
#'   \item{D_pool}{Numeric, same length as \code{lambda_pool}: the pooled
#'     population measure \eqn{\omega_k\, D_k}, summing to 1 overall.}
#'   \item{type_mean}{Numeric length-\code{K} vector (\code{K =
#'     length(mks$blocks)}): the block-\code{D}-mass-weighted mean CEV
#'     WITHIN each type, \eqn{\sum_i D_{k,i} \lambda_{k,i} / \sum_i D_{k,i}}.}
#'   \item{type_share}{Numeric length-\code{K} vector: \code{mks$omega}
#'     renormalised to sum 1 (a no-op copy when \code{mks$omega} already
#'     sums to 1, as \code{\link{hank_mixture_ks_steady}} enforces).}
#'   \item{decompose}{The \code{\link{hank_welfare_decompose}} result
#'     (\code{efficiency}, \code{redistribution}, \code{total},
#'     \code{dispersion}) evaluated at the POOLED level, i.e. across the
#'     WHOLE mixture population rather than within a single type.}
#' @seealso \code{\link{hank_welfare_decompose}},
#'   \code{\link{hank_mixture_ks_steady}}, \code{\link{hank_welfare_response}}
#' @export
hank_mixture_welfare_pool <- function(mks, lambda_list, weights = NULL) {
  ## Input normalization: accept either a hank_mixture_ks_steady() result or
  ## a hank_nk_hank_mixture() result directly -- both expose $blocks (list of
  ## K hank_het_blocks) and $omega (length-K weights) at the top level, so
  ## the only real work here is a clear error when neither shape matches
  ## (e.g. a plain hank_ks/hank_nk single-type object, which has neither
  ## field) rather than an opaque NULL-length downstream failure.
  if (is.null(mks$blocks) || is.null(mks$omega))
    stop(paste0(
      "hank_mixture_welfare_pool(): 'mks' must be a hank_mixture_ks_steady() ",
      "result or a hank_nk_hank_mixture() result (or any object exposing ",
      "$blocks, a list of per-type hank_het_blocks, and $omega, the per-type ",
      "population weights) -- got an object with class(es) '",
      paste(class(mks), collapse = "/"), "' and no $blocks/$omega fields."))

  blocks <- mks$blocks
  K <- length(blocks)
  if (length(lambda_list) != K)
    stop(sprintf(
      "hank_mixture_welfare_pool(): length(lambda_list) (%d) must equal length(mks$blocks) (%d).",
      length(lambda_list), K))
  for (k in seq_len(K)) {
    if (length(lambda_list[[k]]) != length(blocks[[k]]$D))
      stop(sprintf(paste0(
        "hank_mixture_welfare_pool(): lambda_list[[%d]] has length %d but ",
        "mks$blocks[[%d]]$D has length %d -- they must match (same cell ",
        "order as that type's het block)."),
        k, length(lambda_list[[k]]), k, length(blocks[[k]]$D)))
  }
  if (!is.null(weights) && length(weights) != K)
    stop(sprintf(
      "hank_mixture_welfare_pool(): length(weights) (%d) must equal length(mks$blocks) (%d) when supplied.",
      length(weights), K))

  omega <- mks$omega
  type_share <- omega / sum(omega)

  lambda_pool <- unlist(lambda_list, use.names = FALSE)
  D_pool <- unlist(lapply(seq_len(K), function(k) omega[k] * blocks[[k]]$D),
                    use.names = FALSE)

  type_mean <- vapply(seq_len(K), function(k) {
    Dk <- blocks[[k]]$D
    sum(Dk * lambda_list[[k]]) / sum(Dk)
  }, numeric(1))

  weights_pool <- if (is.null(weights)) {
    NULL
  } else {
    for (k in seq_len(K)) {
      if (length(weights[[k]]) != length(blocks[[k]]$D))
        stop(sprintf(paste0(
          "hank_mixture_welfare_pool(): weights[[%d]] has length %d but ",
          "mks$blocks[[%d]]$D has length %d -- they must match."),
          k, length(weights[[k]]), k, length(blocks[[k]]$D)))
    }
    unlist(weights, use.names = FALSE)
  }

  decompose <- hank_welfare_decompose(lambda_pool, D_pool, weights = weights_pool)

  list(lambda_pool = lambda_pool, D_pool = D_pool,
       type_mean = type_mean, type_share = type_share,
       decompose = decompose)
}


#' Posterior-weighted welfare incidence statements
#'
#' Evaluates a user-supplied welfare functional over each draw of an
#' ESTIMATED posterior sample of structural parameters, and summarizes the
#' resulting POSTERIOR DISTRIBUTION of the welfare statement(s) -- as opposed
#' to a welfare statement computed once, at a single point calibration.
#' \code{welfare_fn} is expected to rebuild whatever mixture/het-block economy
#' is needed at the draw's parameter values (e.g. via
#' \code{\link{hank_mixture_ks_steady}}), compute its welfare object(s) (e.g.
#' via \code{\link{hank_welfare_response}} and
#' \code{\link{hank_mixture_welfare_pool}}), and return a NAMED numeric
#' vector of scalar welfare statistics for that draw (e.g.
#' \code{c(efficiency = ..., dispersion = ..., cev_patient = ...)}).
#'
#' Evaluation is FAIL-LOUD: if \code{welfare_fn} errors on any draw, this
#' stops immediately, reporting the offending draw's row index and the
#' original error message, rather than silently coercing that draw to
#' \code{NA} and continuing (matching the package's fail-loud convention for
#' estimation-adjacent code; see \code{p0-estimated-shock-std-likelihood}).
#' All draws must return a welfare vector with the SAME set of names (in any
#' order); a draw returning a different name set also stops loudly, since a
#' silently-changing statistic set would otherwise corrupt the summary
#' data.frame's columns.
#'
#' @param draws Matrix or data.frame of posterior draws: one row per draw,
#'   one named column per structural parameter. Row \code{i} is passed to
#'   \code{welfare_fn} as a named numeric vector.
#' @param welfare_fn A function of one argument (a named numeric vector: one
#'   posterior draw), returning a named numeric vector of welfare statistics.
#' @param probs Numeric vector of posterior quantile levels to summarize
#'   (default \code{c(0.05, 0.25, 0.5, 0.75, 0.95)}).
#' @param cores Integer: number of cores. \code{1L} (default) evaluates
#'   serially; \code{> 1} dispatches via \code{parallel::mclapply()} (not
#'   available on Windows).
#' @param keep_draws Logical: if \code{TRUE} (default), also return the full
#'   per-draw welfare matrix.
#' @param weights Optional numeric vector of length \code{nrow(draws)}: draw
#'   weights, e.g. importance-sampling weights or a normalized posterior
#'   evaluated on a parameter GRID (as opposed to an equal-weight MCMC
#'   sample). \code{NULL} (default) reproduces today's behavior EXACTLY --
#'   equal weights, computed with the unweighted \code{mean}/\code{sd}/
#'   \code{\link[stats]{quantile}} (this is the byte-identical regression
#'   pin; do not change the unweighted code path). When supplied, weights
#'   are validated fail-loud (see below) and then normalized to sum to 1;
#'   the summary's \code{mean}/\code{sd}/quantile columns are then the
#'   WEIGHTED mean, WEIGHTED sd, and WEIGHTED quantiles of each statistic
#'   across draws, using the normalized weights as probability masses.
#'
#'   Weighted quantiles use the standard inverse-CDF of the weighted
#'   empirical distribution: draws are sorted, weights are cumulatively
#'   summed alongside them into a step CDF \code{F(x) = sum_{i: x_i <= x}
#'   w_i}, and the quantile at level \code{p} is the smallest \code{x_i}
#'   with \code{F(x_i) >= p} (a right-continuous step-function inverse --
#'   the weighted analog of \code{\link[stats]{quantile}}'s \code{type = 1}).
#'   With equal weights \code{1/n}, this convention need not coincide
#'   exactly with \code{quantile(..., type = 7)}'s interpolated default (a
#'   deliberate consequence of using a genuinely different, weight-aware
#'   estimator); this is why \code{weights = NULL} is routed through the
#'   original unweighted \code{quantile()} call instead, to guarantee the
#'   byte-identical regression pin.
#'
#'   Validation (fail-loud): \code{length(weights)} must equal
#'   \code{nrow(draws)} (error names both lengths); all weights must be
#'   finite and non-negative; not all weights may be zero (their sum must
#'   be \code{> 0}, since a zero-sum vector cannot be normalized to a
#'   probability distribution).
#'
#' @return A list with:
#'   \item{summary}{A data.frame, one row per welfare statistic, with columns
#'     \code{stat}, \code{mean}, \code{sd}, and one column per entry of
#'     \code{probs} (named \code{"q<pct>"}, e.g. \code{"q50"}). When
#'     \code{weights} is supplied, these are the weighted versions (see
#'     \code{weights} above); otherwise they are the plain unweighted
#'     \code{mean}/\code{sd}/\code{quantile}, byte-identical to before the
#'     \code{weights} argument existed.}
#'   \item{draws_out}{If \code{keep_draws}: numeric \code{n_draws x n_stats}
#'     matrix of the raw per-draw welfare evaluations (column names = the
#'     welfare statistic names); omitted (\code{NULL}) otherwise.}
#'   \item{n_draws}{Integer: number of posterior draws evaluated.}
#'   \item{weights}{Only present when \code{keep_draws = TRUE}: the
#'     NORMALIZED (sum-to-1) weight vector used, length \code{n_draws} --
#'     equal weights \code{rep(1/n_draws, n_draws)} when \code{weights =
#'     NULL} was passed in.}
#' @seealso \code{\link{hank_mixture_welfare_pool}},
#'   \code{\link{hank_welfare_response}}
#' @export
hank_welfare_posterior <- function(draws, welfare_fn,
                                   probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                                   cores = 1L, keep_draws = TRUE,
                                   weights = NULL) {
  draws_mat <- as.matrix(draws)
  n_draws <- nrow(draws_mat)
  col_names <- colnames(draws_mat)

  if (!is.null(weights)) {
    if (length(weights) != n_draws)
      stop(sprintf(
        "hank_welfare_posterior(): length(weights) (%d) must equal nrow(draws) (%d).",
        length(weights), n_draws))
    if (any(!is.finite(weights)) || any(weights < 0))
      stop("hank_welfare_posterior(): 'weights' must be finite and non-negative.")
    if (sum(weights) <= 0)
      stop("hank_welfare_posterior(): 'weights' must not all be zero (sum(weights) must be > 0).")
    weights <- weights / sum(weights)
  }

  eval_one <- function(i) {
    theta <- draws_mat[i, ]
    names(theta) <- col_names
    tryCatch(
      welfare_fn(theta),
      error = function(e)
        stop(sprintf("hank_welfare_posterior(): welfare_fn() failed on draw %d: %s",
                     i, conditionMessage(e)), call. = FALSE)
    )
  }

  results <- if (cores > 1L) {
    parallel::mclapply(seq_len(n_draws), eval_one, mc.cores = cores)
  } else {
    lapply(seq_len(n_draws), eval_one)
  }

  ## Fail loud on any draw whose result isn't a plain named numeric vector,
  ## and on any name-set mismatch across draws (both would otherwise corrupt
  ## the summary data.frame silently).
  stat_names <- names(results[[1L]])
  if (is.null(stat_names) || any(stat_names == ""))
    stop("hank_welfare_posterior(): welfare_fn() must return a NAMED numeric vector ",
         "(draw 1 returned unnamed or partially-named entries).")
  for (i in seq_len(n_draws)) {
    nm_i <- names(results[[i]])
    if (is.null(nm_i) || !setequal(nm_i, stat_names))
      stop(sprintf(paste0(
        "hank_welfare_posterior(): welfare_fn() returned inconsistent ",
        "statistic names on draw %d (got: %s; expected the same name set ",
        "as draw 1: %s)."),
        i, paste(nm_i, collapse = ", "), paste(stat_names, collapse = ", ")))
  }

  draws_out <- do.call(rbind, lapply(results, function(r) r[stat_names]))
  rownames(draws_out) <- NULL
  colnames(draws_out) <- stat_names

  q_labels <- paste0("q", probs * 100)

  if (is.null(weights)) {
    ## Unweighted path: UNCHANGED from before the 'weights' argument existed
    ## -- this is the byte-identical regression pin.
    q_mat <- apply(draws_out, 2L, stats::quantile, probs = probs, names = FALSE)
    ## apply() returns length(probs) x n_stats when probs has length > 1, but
    ## collapses to a plain vector when length(probs) == 1; normalize to a
    ## matrix either way so the summary assembly below is uniform.
    q_mat <- matrix(q_mat, nrow = length(probs), ncol = ncol(draws_out),
                    dimnames = list(q_labels, stat_names))

    summary_df <- data.frame(
      stat = stat_names,
      mean = colMeans(draws_out),
      sd   = apply(draws_out, 2L, stats::sd),
      stringsAsFactors = FALSE
    )
    for (j in seq_along(q_labels))
      summary_df[[q_labels[j]]] <- q_mat[j, ]
    rownames(summary_df) <- NULL
  } else {
    ## Weighted path: weighted mean, weighted sd (frequency-weights formula,
    ## normalized weights so this is just sum(w*(x-mean)^2) since sum(w)=1),
    ## and weighted quantiles (inverse-CDF of the weighted empirical dist.,
    ## see .hank_weighted_quantile() below).
    w_mean <- as.numeric(weights %*% draws_out)
    w_var  <- vapply(seq_len(ncol(draws_out)), function(j) {
      sum(weights * (draws_out[, j] - w_mean[j])^2)
    }, numeric(1))
    w_sd <- sqrt(w_var)

    q_mat <- vapply(seq_len(ncol(draws_out)), function(j) {
      .hank_weighted_quantile(draws_out[, j], weights, probs)
    }, numeric(length(probs)))
    q_mat <- matrix(q_mat, nrow = length(probs), ncol = ncol(draws_out),
                    dimnames = list(q_labels, stat_names))

    summary_df <- data.frame(
      stat = stat_names,
      mean = w_mean,
      sd   = w_sd,
      stringsAsFactors = FALSE
    )
    for (j in seq_along(q_labels))
      summary_df[[q_labels[j]]] <- q_mat[j, ]
    rownames(summary_df) <- NULL
  }

  out <- list(summary = summary_df,
              draws_out = if (keep_draws) draws_out else NULL,
              n_draws = n_draws)
  if (keep_draws)
    out$weights <- if (is.null(weights)) rep(1 / n_draws, n_draws) else weights
  out
}


#' Weighted quantiles via the inverse-CDF of a weighted empirical distribution
#'
#' The weighted analog of \code{\link[stats]{quantile}}'s \code{type = 1}:
#' sorts \code{x}, cumulatively sums the (already-normalized, sum-to-1)
#' weights alongside it into a right-continuous step CDF \code{F(x_(i)) =
#' sum_{j <= i} w_(j)}, and returns, for each level \code{p} in \code{probs},
#' the smallest order statistic \code{x_(i)} with \code{F(x_(i)) >= p}.
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of NORMALIZED (sum-to-1) weights, same length as
#'   \code{x}.
#' @param probs Numeric vector of quantile levels in \code{[0, 1]}.
#' @return Numeric vector, same length as \code{probs}.
#' @keywords internal
.hank_weighted_quantile <- function(x, w, probs) {
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cdf <- cumsum(w_sorted)
  ## Guard the top against floating-point sum(w) slightly < 1: clamp so a
  ## probs == 1 (or a p landing just past the last cumulative-sum tick due
  ## to rounding) still resolves to the last order statistic.
  cdf[length(cdf)] <- max(cdf[length(cdf)], 1)
  vapply(probs, function(p) x_sorted[which(cdf >= p)[1L]], numeric(1))
}
