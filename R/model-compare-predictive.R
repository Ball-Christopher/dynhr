## R/model-compare-predictive.R
## --------------------------------------------------------------------------
## Predictive model comparison:
##
##   posterior_log_scores() -- per-period one-step-ahead log predictive
##                             densities averaged over posterior draws
##                             (log pointwise predictive density, WAIC terms)
##   psis_lfo()             -- approximate leave-FUTURE-out CV (Buerkner,
##                             Gabry & Vehtari 2020) with a local PSIS
##   marginal_likelihoods() -- Laplace / modified harmonic mean / THAMES /
##                             SMC log marginal likelihoods side by side
##
## The per-period contributions come from make_loglik_contrib() (the SAME
## solve + Kalman filter make_log_posterior() uses), or from a caller-supplied
## theta -> contributions closure. The modified harmonic mean estimator lives
## next to THAMES in R/mdd-thames.R.
## --------------------------------------------------------------------------


## ---- internal: what a posterior object carries ------------------------------

## Collect draws (+ their log-posterior values, when paired) and the model /
## data / likelihood configuration from a posterior result. Plain matrices
## carry draws only; everything else must then come through `...`.
.pc_source <- function(result, where) {
  src <- list(draws = NULL, logpost = NULL, support_lower = NULL,
              support_upper = NULL, model = NULL, compiled = NULL,
              data = NULL, obs_vars = NULL, ctx = NULL,
              laplace = NA_real_, smc = NA_real_)

  if (is.matrix(result) || is.data.frame(result)) {
    src$draws <- as.matrix(result)
    return(src)
  }

  if (inherits(result, "dynhr_posterior_result")) {
    parts <- list()
    lps   <- list()
    lp_ok <- TRUE
    for (m in names(result$chains)) {
      for (ch in result$chains[[m]]$chains) {
        e <- .mdd_chain_draws_logpost(ch)
        if (is.null(e$draws)) next
        parts <- c(parts, list(e$draws))
        if (is.null(e$logpost)) lp_ok <- FALSE else lps <- c(lps, list(e$logpost))
        if (!is.finite(src$smc) && is.numeric(ch$log_marginal_lik) &&
            length(ch$log_marginal_lik) == 1L)
          src$smc <- ch$log_marginal_lik
      }
    }
    if (length(parts) > 0L) src$draws <- do.call(rbind, parts)
    if (lp_ok && length(lps) > 0L) src$logpost <- unlist(lps)
    mr <- result$mode_result
    src$model    <- mr$solved$model
    src$compiled <- mr$solved$compiled
    src$data     <- mr$data
    src$obs_vars <- mr$obs_vars
    src$ctx      <- mr$ctx
    if (is.numeric(mr$log_marglik_laplace)) src$laplace <- mr$log_marglik_laplace
    return(src)
  }

  chains <- if (inherits(result, "dynhr_estimation_result")) result$chains
            else result
  if (!is.list(chains) || !is.matrix(chains$chain))
    .dynhr_abort(where, ": `result` must be a dynhr_posterior_result, a ",
                 "dynhr_estimation_result, a dynhr_chains object or a draws ",
                 "matrix (rows = draws, named columns = parameters).",
                 class = "dynhr_error_predictive_input")
  e <- .mdd_chain_draws_logpost(chains)
  src$draws   <- e$draws
  src$logpost <- e$logpost
  src$support_lower <- e$chains$support_lower
  src$support_upper <- e$chains$support_upper
  if (is.numeric(chains$log_marginal_lik) && length(chains$log_marginal_lik) == 1L)
    src$smc <- chains$log_marginal_lik
  if (inherits(result, "dynhr_estimation_result")) {
    src$model    <- result$model
    src$compiled <- result$compiled
    src$data     <- result$data
    src$obs_vars <- result$meta$obs_vars
    src$ctx      <- result$ctx
    if (is.numeric(result$mode$log_marglik_laplace))
      src$laplace <- result$mode$log_marglik_laplace
  }
  src
}

## The theta -> per-period contribution closure: the caller's, or
## make_loglik_contrib() built from the result's configuration with any
## `...` overrides (model, data, obs_vars, compiled, me_variance, lik_init,
## me_extra, shock_scale, likelihood).
.pc_contrib_fn <- function(src, dots, where) {
  if (length(dots) > 0L &&
      (is.null(names(dots)) || any(!nzchar(names(dots)))))
    .dynhr_abort(where, ": arguments in `...` must be named (they override ",
                 "make_loglik_contrib() arguments).",
                 class = "dynhr_error_predictive_input")
  ctx  <- src$ctx %||% list()
  args <- list(model       = src$model,
               data        = src$data,
               obs_vars    = src$obs_vars,
               compiled    = src$compiled,
               me_variance = ctx$me_variance %||% 0,
               lik_init    = ctx$lik_init %||% "auto",
               me_extra    = ctx$me_extra,
               shock_scale = ctx$shock_scale,
               likelihood  = ctx$likelihood %||% "gaussian")
  for (nm in names(dots)) args[nm] <- list(dots[[nm]])
  need <- c("model", "data", "compiled")
  miss <- need[vapply(args[need], is.null, logical(1))]
  if (length(miss) > 0L)
    .dynhr_abort(where, ": no ", paste(miss, collapse = "/"), " to build the ",
                 "per-period likelihood from. Pass a posterior result that ",
                 "carries them, supply them through `...`, or pass ",
                 "`loglik_contrib_fn`.", class = "dynhr_error_predictive_input")
  ## Estimation results store data T x n_obs; accept the transpose as the
  ## estimation drivers do.
  obs <- args$obs_vars %||% args$model$obs_vars %||% args$model$varobs_names
  dat <- as.matrix(args$data)
  if (length(obs) > 0L && ncol(dat) != length(obs) && nrow(dat) == length(obs))
    args$data <- t(dat)
  do.call(make_loglik_contrib, args)
}

## Evenly spaced (in chain order) subset of at most n_draws rows.
.pc_thin_index <- function(S, n_draws, where) {
  if (length(n_draws) != 1L || is.na(n_draws) || n_draws < 2)
    .dynhr_abort(where, ": `n_draws` must be a single number >= 2.",
                 class = "dynhr_error_predictive_input")
  if (n_draws >= S) return(seq_len(S))
  unique(round(seq(1, S, length.out = n_draws)))
}

## Evaluate the contribution closure at every draw. Returns the S x T matrix
## of per-period log contributions and, when the closure returns a T x J
## matrix (one column per observable), the S x T x J array of cells.
.pc_loglik_eval <- function(draws, fn, where) {
  pn <- colnames(draws)
  if (is.null(pn))
    .dynhr_abort(where, ": the draws need parameter names (column names).",
                 class = "dynhr_error_predictive_input")
  vals <- lapply(seq_len(nrow(draws)), function(s)
    fn(stats::setNames(as.numeric(draws[s, ]), pn)))
  first <- vals[[1L]]
  cells <- NULL
  if (is.matrix(first)) {
    dims <- dim(first)
    ok <- vapply(vals, function(v) is.matrix(v) && identical(dim(v), dims),
                 logical(1))
    if (!all(ok))
      .dynhr_abort(where, ": `loglik_contrib_fn` returned matrices of ",
                   "different shapes across draws.",
                   class = "dynhr_error_predictive_input")
    cells <- array(unlist(lapply(vals, as.numeric)),
                   dim = c(dims, length(vals)))
    cells <- aperm(cells, c(3L, 1L, 2L))           # S x T x J
    dimnames(cells) <- list(NULL, rownames(first), colnames(first))
    ll <- apply(cells, c(1L, 2L), sum)
  } else {
    n_T <- length(first)
    ok  <- vapply(vals, function(v) length(v) == n_T, logical(1))
    if (!all(ok))
      .dynhr_abort(where, ": `loglik_contrib_fn` returned vectors of ",
                   "different lengths across draws.",
                   class = "dynhr_error_predictive_input")
    ll <- matrix(unlist(vals), nrow = length(vals), ncol = n_T, byrow = TRUE)
  }
  if (anyNA(ll))
    .dynhr_abort(where, ": `loglik_contrib_fn` returned NA/NaN contributions.",
                 class = "dynhr_error_predictive_input")
  if (length(ll) == 0L || ncol(ll) < 1L)
    .dynhr_abort(where, ": the contribution closure returned no periods.",
                 class = "dynhr_error_predictive_input")
  list(ll = ll, cells = cells)
}

## Column-wise log-mean-exp over draws, with the relative influence values
## p_st / mean_s(p_st) - 1 that give delta-method Monte Carlo errors.
.pc_lme_cols <- function(ll) {
  S   <- nrow(ll)
  m   <- apply(ll, 2L, max)
  fin <- is.finite(m)
  p   <- exp(sweep(ll, 2L, ifelse(fin, m, 0), "-"))
  pbar <- colMeans(p)
  lpd <- ifelse(fin, m + log(pbar), m)
  rel <- sweep(p, 2L, ifelse(fin, pbar, 1), "/") - 1
  rel[, !fin] <- 0
  list(lpd = lpd, rel = rel,
       se = if (S > 1L) apply(rel, 2L, stats::sd) / sqrt(S) else rep(NA_real_, ncol(ll)))
}


## ---- posterior_log_scores ---------------------------------------------------

#' Per-period one-step-ahead posterior log predictive scores
#'
#' For every period \eqn{t}, the log of the posterior-averaged one-step-ahead
#' predictive density
#' \deqn{\mathrm{lpd}_t = \log \frac{1}{S}\sum_{s=1}^S p(y_t \mid y_{1:t-1},
#' \theta^{(s)}),}
#' computed with log-sum-exp over a thinned subset of posterior draws, where
#' \eqn{\log p(y_t \mid y_{1:t-1}, \theta)} is the Kalman filter's
#' prediction-error decomposition from \code{\link{make_loglik_contrib}}
#' (so \eqn{\sum_t} of the per-draw contributions is exactly the likelihood
#' \code{make_log_posterior()} scores).
#'
#' @details
#' \strong{What the total measures.} \eqn{\sum_t \mathrm{lpd}_t} is the
#' in-sample log pointwise predictive density (Vehtari, Gelman & Gabry 2017):
#' the posterior uses all \eqn{T} periods, so it is optimistic as an estimate
#' of out-of-sample predictive accuracy. \code{elpd_waic} subtracts the WAIC
#' effective number of parameters \eqn{p_{\mathrm{waic}} = \sum_t
#' \mathrm{Var}_s[\log p(y_t\mid y_{1:t-1},\theta^{(s)})]}; for a genuinely
#' out-of-sample (leave-future-out) score use \code{\link{psis_lfo}}.
#' Both are on the natural-log scale and comparable across models fitted to
#' the SAME data.
#'
#' \strong{Per observable.} The Kalman filter exposes one joint contribution
#' per period, so the default closure gives per-period values only
#' (\code{per_observable = NULL}). A \code{loglik_contrib_fn} returning a
#' \eqn{T \times J} matrix (e.g. a sequential, observable-by-observable
#' decomposition of each period's density) is accepted: \code{per_observable}
#' then holds the \eqn{T \times J} cell scores, while the per-period values
#' and the total are computed from the row sums (the joint density of each
#' period), not by adding the cell scores.
#'
#' \strong{Monte Carlo error.} Delta method on the draws, treating the thinned
#' subset as independent: \code{pointwise$mc_se} per period, and \code{mc_se}
#' for the total, which accounts for the correlation between periods (every
#' period is averaged over the same draws).
#'
#' @param result A \code{dynhr_posterior_result}
#'   (\code{\link{run_posterior_estimation}}), a \code{dynhr_estimation_result}
#'   (\code{\link{run_full_estimation}}), a \code{dynhr_chains} object, or a
#'   numeric draws matrix with named columns. The first two carry the model,
#'   data and likelihood configuration; for the others supply them through
#'   \code{...} or pass \code{loglik_contrib_fn}.
#' @param n_draws Number of posterior draws to average over (an evenly spaced
#'   subset in chain order when more are available). Default 1000.
#' @param loglik_contrib_fn Optional \code{function(theta)} returning the
#'   per-period log contributions \eqn{\log p(y_t \mid y_{1:t-1}, \theta)}
#'   (length \eqn{T}, or a \eqn{T \times J} matrix). Default: built with
#'   \code{\link{make_loglik_contrib}} from \code{result}.
#' @param ... Named overrides of the \code{\link{make_loglik_contrib}}
#'   arguments (\code{model}, \code{data}, \code{obs_vars}, \code{compiled},
#'   \code{me_variance}, \code{lik_init}, \code{me_extra}, \code{shock_scale},
#'   \code{likelihood}); ignored when \code{loglik_contrib_fn} is given.
#' @return A list with \code{lpd} (total), \code{mc_se}, \code{p_waic},
#'   \code{elpd_waic}, \code{pointwise} (data frame: \code{period},
#'   \code{lpd}, \code{mc_se}, \code{p_waic}; its \code{lpd} column sums to
#'   \code{lpd}), \code{per_observable} (\code{NULL} or the \eqn{T \times J}
#'   matrix described above), \code{loglik} (the \eqn{S \times T} matrix of
#'   contributions), \code{draw_index} (rows of the draws used) and
#'   \code{n_draws}.
#' @references Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian
#'   model evaluation using leave-one-out cross-validation and WAIC.
#'   \emph{Statistics and Computing}, 27, 1413-1432.
#' @seealso \code{\link{psis_lfo}}, \code{\link{marginal_likelihoods}},
#'   \code{\link{make_loglik_contrib}}, \code{\link{score_forecast}}
#' @examples
#' ## Conjugate AR(1) with known noise: y_t = rho * y_{t-1} + e_t.
#' set.seed(1)
#' y <- as.numeric(stats::arima.sim(list(ar = 0.6), 60))
#' contrib <- function(theta)
#'   stats::dnorm(y[-1], theta[["rho"]] * y[-length(y)], 1, log = TRUE)
#' draws <- cbind(rho = stats::rnorm(500, 0.6, 0.1))
#' posterior_log_scores(draws, loglik_contrib_fn = contrib)$lpd
#' @export
posterior_log_scores <- function(result, n_draws = 1000L,
                                 loglik_contrib_fn = NULL, ...) {
  where <- "posterior_log_scores"
  src <- .pc_source(result, where)
  if (is.null(src$draws) || nrow(src$draws) < 1L)
    .dynhr_abort(where, ": `result` holds no posterior draws.",
                 class = "dynhr_error_predictive_input")
  fn <- if (is.null(loglik_contrib_fn)) .pc_contrib_fn(src, list(...), where)
        else loglik_contrib_fn
  if (!is.function(fn))
    .dynhr_abort(where, ": `loglik_contrib_fn` must be a function(theta).",
                 class = "dynhr_error_predictive_input")

  idx   <- .pc_thin_index(nrow(src$draws), n_draws, where)
  draws <- src$draws[idx, , drop = FALSE]
  ev    <- .pc_loglik_eval(draws, fn, where)
  ll    <- ev$ll
  n_bad <- sum(apply(ll, 1L, function(r) all(r == -Inf)))
  if (n_bad > 0L)
    .dynhr_warn(sprintf(paste0(
      "%s: %d of %d posterior draws have zero likelihood (the model does not ",
      "solve there); they count as zero predictive density. A posterior draw ",
      "should never do this -- check that the likelihood configuration ",
      "matches the one the posterior was sampled with."),
      where, n_bad, nrow(ll)), class = "dynhr_warning_predictive_infeasible")

  lme    <- .pc_lme_cols(ll)
  p_waic <- apply(ll, 2L, function(v) if (all(is.finite(v)) && length(v) > 1L)
                                        stats::var(v) else NA_real_)
  per_obs <- NULL
  if (!is.null(ev$cells)) {
    J <- dim(ev$cells)[3L]
    per_obs <- vapply(seq_len(J), function(j) .pc_lme_cols(ev$cells[, , j])$lpd,
                      numeric(ncol(ll)))
    per_obs <- matrix(per_obs, ncol = J,
                      dimnames = list(NULL, dimnames(ev$cells)[[3L]]))
  }
  total <- sum(lme$lpd)
  list(lpd            = total,
       mc_se          = .mdd_mean_se(rowSums(lme$rel), "iid"),
       p_waic         = sum(p_waic),
       elpd_waic      = total - sum(p_waic),
       pointwise      = data.frame(period = seq_len(ncol(ll)), lpd = lme$lpd,
                                   mc_se = lme$se, p_waic = p_waic),
       per_observable = per_obs,
       loglik         = ll,
       draw_index     = idx,
       n_draws        = nrow(ll))
}


## ---- PSIS ---------------------------------------------------------------------

## Generalised Pareto fit of the exceedances `x` (sorted ascending, >= 0) by
## the Zhang & Stephens (2009) empirical-Bayes profile, with the weakly
## informative shrinkage of k towards 0.5 (Vehtari et al. 2024), exactly as
## loo:::gpdfit(). Returns the shape k (> 0 = heavy tail) and scale sigma.
.psis_gpdfit <- function(x) {
  N  <- length(x)
  M  <- 30 + floor(sqrt(N))
  jj <- seq_len(M)
  xstar <- x[max(1L, floor(N / 4 + 0.5))]
  theta <- 1 / x[N] + (1 - sqrt(M / (jj - 0.5))) / 3 / xstar
  l_theta <- N * vapply(theta, function(th) {
    kk <- mean(log1p(-th * x))
    log(-th / kk) - kk - 1
  }, numeric(1))
  l_theta[!is.finite(l_theta)] <- -Inf
  if (!any(is.finite(l_theta))) return(list(k = Inf, sigma = NA_real_))
  w <- exp(l_theta - .logsumexp(l_theta))
  theta_hat <- sum(theta * w)
  k     <- mean(log1p(-theta_hat * x))
  sigma <- -k / theta_hat
  k <- k * N / (N + 10) + 10 * 0.5 / (N + 10)
  if (is.nan(k)) k <- Inf
  list(k = k, sigma = sigma)
}

## Pareto-smoothed importance sampling (Vehtari, Simpson, Gelman, Yao & Gabry
## 2024) of a vector of log importance ratios, following loo::psis() with
## r_eff = 1: the M = ceiling(min(0.2 S, 3 sqrt(S))) largest ratios are
## replaced by the expected order statistics of the fitted generalised Pareto
## tail and every weight is truncated at the largest raw ratio. Returns the
## smoothed (unnormalised) log weights and the Pareto k-hat: < 0.5 fine,
## 0.5-0.7 usable, > 0.7 the importance-sampling estimate is unreliable.
## Bounded ratios (constant tail) give k = -Inf.
.psis_smooth <- function(log_ratios) {
  S  <- length(log_ratios)
  mx <- max(log_ratios)
  lw <- log_ratios - mx
  M  <- ceiling(min(0.2 * S, 3 * sqrt(S)))
  k  <- Inf
  if (M >= 5L && S > M) {
    ord      <- order(lw)
    tail_ids <- ord[(S - M + 1L):S]
    lw_tail  <- lw[tail_ids]
    if (max(lw_tail) - min(lw_tail) <= .Machine$double.eps / 100) {
      k <- -Inf
    } else {
      cutoff <- lw[ord[S - M]]
      fit <- .psis_gpdfit(exp(lw_tail) - exp(cutoff))
      k   <- fit$k
      if (is.finite(k) && is.finite(fit$sigma) && fit$sigma > 0) {
        p  <- (seq_len(M) - 0.5) / M
        qq <- if (abs(k) < 1e-12) -fit$sigma * log1p(-p)
              else fit$sigma * expm1(-k * log1p(-p)) / k
        lw[tail_ids] <- log(qq + exp(cutoff))
      }
    }
  }
  lw[lw > 0] <- 0
  list(log_weights = lw + mx, pareto_k = k, tail_len = M)
}


## ---- psis_lfo -----------------------------------------------------------------

#' Approximate leave-future-out cross-validation (PSIS-LFO)
#'
#' Estimates the one-step-ahead leave-future-out expected log predictive
#' density
#' \deqn{\mathrm{elpd}_{\mathrm{LFO}} = \sum_{i=L}^{T-1}
#' \log p(y_{i+1} \mid y_{1:i}),}
#' where each term uses the posterior given ONLY \eqn{y_{1:i}}, without
#' refitting the model \eqn{T-L} times: the posterior draws are importance
#' weighted by \eqn{p(\theta\mid y_{1:i}) / p(\theta\mid y_{1:i^*})
#' \propto \exp\{\sum_{t \le i}\ell_t(\theta) - \sum_{t \le i^*}\ell_t(\theta)\}}
#' (\eqn{\ell_t} the per-period log contributions, \eqn{i^*} the sample of the
#' last fit) and the weights are Pareto smoothed (Buerkner, Gabry & Vehtari
#' 2020).
#'
#' @details
#' \strong{Direction.} The pass starts from the full-sample posterior
#' (\eqn{i^* = T}) that \code{result} already holds and moves BACKWARD,
#' \eqn{i = T-1, \dots, L}, so no fit to \eqn{y_{1:L}} is needed up front (the
#' paper's forward pass starts from one). The further \eqn{i} is from
#' \eqn{i^*}, the more the weights degenerate; the Pareto \eqn{\hat k} of each
#' step measures that.
#'
#' \strong{Refits.} When \eqn{\hat k >} \code{k_threshold} and
#' \code{refit_fn} is supplied, the model is refitted on \eqn{y_{1:i}}:
#' \code{refit_fn(i)} must return posterior draws given the first \eqn{i}
#' periods only (a draws matrix with the same named columns, or any object
#' \code{posterior_log_scores()} accepts as \code{result}); those draws become
#' the new reference (\eqn{i^* = i}, exact for that step). Their contributions
#' are still evaluated on the FULL sample by the same closure -- \eqn{\ell_t}
#' for \eqn{t \le i} does not depend on later data. Without \code{refit_fn},
#' every step uses PSIS from the full-sample fit and a warning lists the
#' predicted periods whose \eqn{\hat k} exceeds the threshold: those terms are
#' unreliable, and an \code{elpd_lfo} built on many of them should not be
#' trusted.
#'
#' \strong{PSIS} is implemented locally (the \code{loo} package's algorithm
#' with relative efficiency 1: generalised Pareto tail fit by Zhang & Stephens
#' 2009 with the weakly informative prior, smoothing of the largest
#' \eqn{\lceil\min(0.2S, 3\sqrt S)\rceil} weights, truncation at the largest raw
#' weight).
#'
#' \strong{Errors.} \code{se} is the usual data-uncertainty standard error
#' \eqn{\sqrt{n\,\mathrm{Var}(\mathrm{elpd}_i)}} used when comparing models;
#' \code{mc_se} is the Monte Carlo (importance-sampling) error of
#' \code{elpd_lfo} by the delta method on the self-normalised weights,
#' accounting for the correlation between steps that share draws.
#'
#' @param result Posterior (full-sample) result; see
#'   \code{\link{posterior_log_scores}}.
#' @param L Number of initial periods always conditioned on: predictions are
#'   made for periods \eqn{L+1, \dots, T}. An integer in \eqn{[1, T-1]}.
#' @param n_draws Number of posterior draws used (evenly spaced subset).
#'   Default 1000. Also applied to the draws a \code{refit_fn} returns.
#' @param k_threshold Pareto \eqn{\hat k} above which a step is refitted
#'   (with \code{refit_fn}) or flagged. Default 0.7.
#' @param refit_fn Optional \code{function(i)} returning posterior draws
#'   given \eqn{y_{1:i}}; see Details.
#' @param loglik_contrib_fn,... As in \code{\link{posterior_log_scores}}. The
#'   closure must return ONE contribution per period (per-observable matrices
#'   are summed over observables).
#' @return A list with \code{elpd_lfo}, \code{se}, \code{mc_se},
#'   \code{pointwise} (data frame: \code{i}, \code{period} \eqn{= i+1},
#'   \code{elpd}, \code{pareto_k}, \code{refit}, \code{mc_se}; its \code{elpd}
#'   column sums to \code{elpd_lfo}), \code{n_refits}, \code{high_k_periods}
#'   (predicted periods with \eqn{\hat k >} \code{k_threshold} that were NOT
#'   refitted), \code{L}, \code{k_threshold} and \code{n_draws}.
#' @references Buerkner, P.-C., Gabry, J., & Vehtari, A. (2020). Approximate
#'   leave-future-out cross-validation for Bayesian time series models.
#'   \emph{Journal of Statistical Computation and Simulation}, 90(14),
#'   2499-2523.
#'
#'   Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024). Pareto
#'   smoothed importance sampling. \emph{Journal of Machine Learning Research},
#'   25(72), 1-58.
#'
#'   Zhang, J., & Stephens, M. A. (2009). A new and efficient estimation
#'   method for the generalized Pareto distribution. \emph{Technometrics},
#'   51(3), 316-325.
#' @seealso \code{\link{posterior_log_scores}}, \code{\link{marginal_likelihoods}}
#' @examples
#' set.seed(1)
#' y <- as.numeric(stats::arima.sim(list(ar = 0.6), 60))
#' contrib <- function(theta)
#'   stats::dnorm(y[-1], theta[["rho"]] * y[-length(y)], 1, log = TRUE)
#' ## exact conjugate posterior of rho (flat prior, unit noise variance)
#' x <- y[-length(y)]
#' draws <- cbind(rho = stats::rnorm(2000, sum(x * y[-1]) / sum(x^2),
#'                                   1 / sqrt(sum(x^2))))
#' psis_lfo(draws, L = 40, loglik_contrib_fn = contrib)$elpd_lfo
#' @export
psis_lfo <- function(result, L, n_draws = 1000L, k_threshold = 0.7,
                     refit_fn = NULL, loglik_contrib_fn = NULL, ...) {
  where <- "psis_lfo"
  if (missing(L))
    .dynhr_abort(where, ": `L` (the number of initial periods always ",
                 "conditioned on) is required.",
                 class = "dynhr_error_predictive_input")
  if (length(k_threshold) != 1L || !is.finite(k_threshold))
    .dynhr_abort(where, ": `k_threshold` must be a single finite number.",
                 class = "dynhr_error_predictive_input")
  if (!is.null(refit_fn) && !is.function(refit_fn))
    .dynhr_abort(where, ": `refit_fn` must be NULL or a function(i).",
                 class = "dynhr_error_predictive_input")
  src <- .pc_source(result, where)
  if (is.null(src$draws) || nrow(src$draws) < 2L)
    .dynhr_abort(where, ": `result` holds fewer than 2 posterior draws.",
                 class = "dynhr_error_predictive_input")
  fn <- if (is.null(loglik_contrib_fn)) .pc_contrib_fn(src, list(...), where)
        else loglik_contrib_fn
  if (!is.function(fn))
    .dynhr_abort(where, ": `loglik_contrib_fn` must be a function(theta).",
                 class = "dynhr_error_predictive_input")

  ## Contributions at a (thinned) draw set; draws the model cannot solve at
  ## carry no posterior mass and would make every ratio NaN, so drop them.
  .ref_at <- function(draws, i_ref) {
    draws <- draws[.pc_thin_index(nrow(draws), n_draws, where), , drop = FALSE]
    ll <- .pc_loglik_eval(draws, fn, where)$ll
    ok <- apply(ll, 1L, function(r) all(is.finite(r)))
    if (!all(ok))
      .dynhr_warn(sprintf(paste0("%s: dropping %d of %d draws with a ",
                                 "non-finite likelihood contribution."),
                          where, sum(!ok), length(ok)),
                  class = "dynhr_warning_predictive_infeasible")
    if (sum(ok) < 2L)
      .dynhr_abort(where, ": fewer than 2 draws with a finite likelihood.",
                   class = "dynhr_error_predictive_input")
    ll  <- ll[ok, , drop = FALSE]
    cum <- matrix(t(apply(ll, 1L, cumsum)), nrow = nrow(ll))
    list(i = i_ref, ll = ll, cum = cum)
  }

  ref <- .ref_at(src$draws, NA_integer_)
  n_T <- ncol(ref$ll)
  ref$i <- n_T
  if (length(L) != 1L || is.na(L) || L != round(L) || L < 1 || L > n_T - 1)
    .dynhr_abort(sprintf("%s: `L` must be an integer in [1, %d] (T = %d).",
                         where, n_T - 1L, n_T),
                 class = "dynhr_error_predictive_input")
  L <- as.integer(L)

  i_seq <- seq(n_T - 1L, L)
  n_i   <- length(i_seq)
  elpd  <- numeric(n_i)
  khat  <- rep(NA_real_, n_i)
  refit <- logical(n_i)
  mcse  <- numeric(n_i)
  mc_var <- 0
  phi   <- numeric(nrow(ref$ll))

  for (j in seq_len(n_i)) {
    i <- i_seq[j]
    lr <- ref$cum[, i] - ref$cum[, ref$i]
    ps <- .psis_smooth(lr)
    khat[j] <- ps$pareto_k
    if (!is.null(refit_fn) && ps$pareto_k > k_threshold) {
      new_src <- .pc_source(refit_fn(i), where)
      if (is.null(new_src$draws) || nrow(new_src$draws) < 2L)
        .dynhr_abort(where, ": `refit_fn(", i, ")` returned fewer than 2 ",
                     "draws.", class = "dynhr_error_predictive_input")
      mc_var <- mc_var + sum(phi^2)
      ref <- .ref_at(new_src$draws, i)
      phi <- numeric(nrow(ref$ll))
      refit[j] <- TRUE
      ps <- list(log_weights = rep(0, nrow(ref$ll)))
    }
    w <- exp(ps$log_weights - max(ps$log_weights))
    w <- w / sum(w)
    lp_next <- ref$ll[, i + 1L]
    mx <- max(lp_next)
    p  <- exp(lp_next - mx)
    pbar <- sum(w * p)
    elpd[j] <- mx + log(pbar)
    infl <- w * (p / pbar - 1)
    mcse[j] <- sqrt(sum(infl^2))
    phi <- phi + infl
  }
  mc_var <- mc_var + sum(phi^2)

  ord <- rev(seq_len(n_i))
  pointwise <- data.frame(i = i_seq[ord], period = i_seq[ord] + 1L,
                          elpd = elpd[ord], pareto_k = khat[ord],
                          refit = refit[ord], mc_se = mcse[ord])
  high <- pointwise$period[!pointwise$refit &
                           pointwise$pareto_k > k_threshold]
  if (length(high) > 0L)
    .dynhr_warn(sprintf(paste0(
      "%s: Pareto k > %.2f at %d of %d predicted period(s) (%s%s); the PSIS ",
      "estimate of those terms is unreliable. Supply `refit_fn` to refit ",
      "the model there, or raise L."),
      where, k_threshold, length(high), n_i,
      paste(utils::head(high, 20L), collapse = ", "),
      if (length(high) > 20L) ", ..." else ""),
      class = "dynhr_warning_psis_high_k")

  list(elpd_lfo       = sum(elpd),
       se             = sqrt(n_i * stats::var(elpd)),
       mc_se          = sqrt(mc_var),
       pointwise      = pointwise,
       n_refits       = sum(refit),
       high_k_periods = high,
       L              = L,
       k_threshold    = k_threshold,
       n_draws        = nrow(ref$ll))
}


## ---- marginal_likelihoods -------------------------------------------------------

#' Log marginal likelihoods side by side (Laplace, MHM, THAMES, SMC)
#'
#' Collects every log marginal-likelihood (MDD) estimate a posterior result
#' supports into one table, the way Dynare reports the Laplace and modified
#' harmonic mean estimates together:
#' \describe{
#'   \item{\code{laplace}}{\code{\link{laplace_log_marglik}} at the mode
#'     (from \code{mode_result}, or the one the result carries); needs an
#'     exact Hessian. Deterministic, no s.e.}
#'   \item{\code{mhm}}{\code{\link{mdd_modified_harmonic_mean}} (Geweke 1999)
#'     on the draws and their log-posterior values.}
#'   \item{\code{thames}}{\code{\link{thames_mdd}} on the same pairs, with
#'     the support correction when the result records prior bounds.}
#'   \item{\code{smc}}{the SMC normalising-constant estimate, when an SMC
#'     sampler produced the draws.}
#' }
#' An estimator the result cannot support is reported as \code{NA} with the
#' reason in \code{note}; the function does not error for that.
#'
#' Disagreement between the draw-based estimators well beyond their
#' \code{se} is itself a diagnostic (non-Gaussian or multimodal posterior,
#' draws piled against a bound); do not pick the most favourable row.
#'
#' @param result A \code{dynhr_posterior_result}, a
#'   \code{dynhr_estimation_result} or a \code{dynhr_chains} object (the draws
#'   must be paired with their \code{post_logpost} values).
#' @param tau Truncation probability (or vector) for the modified harmonic
#'   mean; \code{seq(0.1, 0.9, 0.1)} reproduces Dynare's averaged estimate.
#'   Default 0.5.
#' @param mode_result Optional \code{dynhr_mode_result}
#'   (\code{\link{run_mode_finding}}) for the Laplace row when \code{result}
#'   does not carry one.
#' @param se_method \code{"batch"} (default; honest for autocorrelated MCMC
#'   draws) or \code{"iid"}, for the MHM and THAMES standard errors.
#' @return A data frame with columns \code{estimator}, \code{log_mdd},
#'   \code{se} and \code{note}, one row per estimator.
#' @seealso \code{\link{mdd_modified_harmonic_mean}}, \code{\link{thames_mdd}},
#'   \code{\link{laplace_log_marglik}}, \code{\link{posterior_log_scores}}
#' @examples
#' set.seed(1)
#' S <- matrix(c(1, 0.5, 0.5, 2), 2)
#' draws <- MASS::mvrnorm(4000, mu = c(1, -1), Sigma = S)
#' colnames(draws) <- c("a", "b")
#' ch <- list(chain = draws,
#'            post_logpost = -0.5 * mahalanobis(draws, c(1, -1), S))
#' marginal_likelihoods(ch, se_method = "iid")
#' @export
marginal_likelihoods <- function(result, tau = 0.5, mode_result = NULL,
                                 se_method = c("batch", "iid")) {
  se_method <- match.arg(se_method)
  src <- .pc_source(result, "marginal_likelihoods")

  rows <- list()
  .row <- function(est, val, se, note)
    data.frame(estimator = est, log_mdd = as.numeric(val),
               se = as.numeric(se), note = note, stringsAsFactors = FALSE)

  ## Laplace
  lap <- src$laplace
  if (!is.null(mode_result)) {
    lap <- if (is.numeric(mode_result$log_marglik_laplace) &&
               is.finite(mode_result$log_marglik_laplace))
             mode_result$log_marglik_laplace
           else if (!is.null(mode_result$hessian_exact))
             laplace_log_marglik(mode_result)
           else NA_real_
  }
  rows$laplace <- .row("laplace", lap, NA_real_,
                       if (is.finite(lap)) "mode + exact Hessian (deterministic)"
                       else "no exact Hessian at the mode")

  ## Draw-based estimators need (draw, log-posterior) pairs.
  draws <- src$draws
  lp    <- src$logpost
  d     <- if (is.null(draws)) 0L else ncol(draws)
  n_fin <- if (is.null(lp)) 0L else sum(is.finite(lp))
  why <- if (is.null(draws)) "no posterior draws"
         else if (is.null(lp)) "draws not paired with post_logpost values"
         else if (n_fin <= 2L * d) sprintf("too few finite draws (%d, d = %d)", n_fin, d)
         else NULL
  if (is.null(why)) {
    keep <- is.finite(lp)
    ev <- eigen(stats::cov(draws[keep, , drop = FALSE]), symmetric = TRUE,
                only.values = TRUE)$values
    if (!all(is.finite(ev)) || min(ev) <= max(ev) * 1e-14)
      why <- "singular covariance of the draws"
  }
  if (is.null(why)) {
    mhm <- mdd_modified_harmonic_mean(draws, lp, tau = tau, se_method = se_method)
    rows$mhm <- .row("mhm", mhm$log_mdd, mhm$se,
                     sprintf("Geweke (1999), tau = %s, %s s.e.",
                             paste(format(tau), collapse = "/"), se_method))
    th <- thames_mdd_from_chains(list(chain = draws, post_logpost = lp,
                                      support_lower = src$support_lower,
                                      support_upper = src$support_upper),
                                 se_method = se_method)
    rows$thames <- .row("thames", th$log_mdd, th$se,
                        if (is.finite(th$log_mdd))
                          sprintf("Metodiev et al. (2023), %s s.e.", se_method)
                        else "THAMES failed (see message)")
  } else {
    rows$mhm    <- .row("mhm", NA_real_, NA_real_, why)
    rows$thames <- .row("thames", NA_real_, NA_real_, why)
  }

  rows$smc <- .row("smc", src$smc, NA_real_,
                   if (is.finite(src$smc)) "SMC normalising constant"
                   else "not an SMC run")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
