## R/forecast-eval.R
## --------------------------------------------------------------------------
## Recursive (expanding-window) out-of-sample forecast evaluation.
##
## `score_forecast()` (R/diag-scoring-rules.R) scores ONE predictive against
## ONE realisation. This file is the driver that produces those pairs: it
## walks a set of forecast origins, re-estimates or re-filters at each one,
## forms the h-step predictive distribution, and scores it against what
## actually happened -- the standard pseudo-real-time backtest every policy
## shop runs before it trusts a model's forecasts.
##
## PREDICTIVE DISTRIBUTION.  For the linear Gaussian state space
##
##   s_t = T s_{t-1} + R eps_t,   y_t = Z s_{t-1} + D eps_t,  eps ~ N(0, Q)
##
## (dynhr Convention A: y_t loads on the LAGGED state) the h-step-ahead
## predictive given data through the origin t is exactly Gaussian:
##
##   m_0 = s_{t|t},                V_0 = P_{t|t}
##   m_k = T m_{k-1},              V_k = T V_{k-1} T' + R Q R'      (k >= 1)
##   E[y_{t+h} | F_t]   = Z m_{h-1}
##   Var[y_{t+h} | F_t] = Z V_{h-1} Z' + D Q D' + me_variance I
##
## The `+ D Q D'` term is the contemporaneous shock loading and the
## `V_{h-1}` term carries BOTH the accumulated future shocks and the
## filtering uncertainty about the origin state. Dropping P_{t|t} (a common
## shortcut: forecast from the point estimate of the terminal state) makes
## the predictive too sharp and the PITs systematically over-dispersed --
## which is exactly what the PIT calibration gate in
## tests/testthat/test-forecast-eval.R would catch.
##
## Because these moments are exact, the default path scores them in closed
## form (`n_draws = 0`): no Monte-Carlo noise in CRPS, the log score, the
## PITs or the interval coverage. `n_draws > 0` switches to an ensemble
## drawn from the same moments, which is what the multivariate energy and
## variogram scores need.
##
## BENCHMARK.  `model = "rw"` (or "ar1"/"mean") runs the same driver on a
## naive reduced-form benchmark instead of the structural model, so that
## `dm_test()` has something to compare against. A random walk is the
## benchmark a forecast evaluation is expected to beat.
##
## Public API: forecast_backtest(), plus summary/print/plot methods.
## --------------------------------------------------------------------------


## ---- Input plumbing ------------------------------------------------------

#' Coerce and validate the backtest data matrix
#'
#' @param data     T x n_obs numeric matrix or data.frame.
#' @param obs_vars Character vector of observable names, or NULL.
#' @return List with `Y` (T x n_obs matrix with column names) and `obs_vars`.
#' @noRd
.fbt_data <- function(data, obs_vars) {
  Y <- as.matrix(data)
  if (!is.numeric(Y))
    stop("forecast_backtest: 'data' must be numeric.", call. = FALSE)
  if (is.null(colnames(Y)) && !is.null(obs_vars) && ncol(Y) == length(obs_vars))
    colnames(Y) <- obs_vars
  if (is.null(obs_vars)) {
    if (is.null(colnames(Y)))
      stop("forecast_backtest: supply 'obs_vars', or give 'data' column names.",
           call. = FALSE)
    obs_vars <- colnames(Y)
  }
  ## SELECT BY NAME whenever the data is labelled -- never by position. A
  ## same-length obs_vars in a different order is a REORDER request, and
  ## silently relabelling the columns instead would feed each observable's
  ## series to the wrong equation while every dimension still checked out.
  if (!is.null(colnames(Y))) {
    miss <- setdiff(obs_vars, colnames(Y))
    if (length(miss))
      stop(sprintf(paste0("forecast_backtest: obs_vars not found in the data ",
                          "columns: %s (data has: %s)"),
                   paste(miss, collapse = ", "),
                   paste(colnames(Y), collapse = ", ")), call. = FALSE)
    Y <- Y[, obs_vars, drop = FALSE]
  } else {
    if (length(obs_vars) != ncol(Y))
      stop(sprintf(paste0("forecast_backtest: 'data' has no column names and ",
                          "%d column(s), but obs_vars has %d entries."),
                   ncol(Y), length(obs_vars)), call. = FALSE)
    colnames(Y) <- obs_vars
  }
  if (anyNA(Y))
    stop("forecast_backtest: 'data' contains NA. The recursive driver scores ",
         "realised values, so missing observations must be handled (or the ",
         "affected origins dropped) before calling it.", call. = FALSE)
  list(Y = Y, obs_vars = obs_vars)
}

#' Pull (model, compiled, dr, params) out of whatever the caller passed
#' @noRd
.fbt_parts <- function(model) {
  if (inherits(model, "dynhr_solved")) {
    return(list(model = model$model, compiled = model$compiled,
                dr = model$dr, params = model$params %||%
                  model$model$param_values))
  }
  if (is.list(model) && !is.null(model$model) && !is.null(model$dr)) {
    m <- model$model
    return(list(model = m, compiled = model$compiled,
                dr = model$dr, params = model$params %||% m$param_values))
  }
  stop("forecast_backtest: 'model' must be a dynhr_solved object (from ",
       "solve_model()), a list with $model/$compiled/$dr, or one of the ",
       "benchmark names \"rw\", \"ar1\", \"mean\".", call. = FALSE)
}

#' Re-solve the decision rule at a new parameter vector
#' @noRd
.fbt_resolve <- function(parts, params) {
  if (is.null(parts$compiled))
    stop("forecast_backtest: re-estimation needs the compiled model; pass a ",
         "dynhr_solved object from solve_model().", call. = FALSE)
  ssv <- solve_steady(parts$compiled, params,
                      endo_names = parts$model$var_names,
                      exo_names  = parts$model$varexo_names,
                      verbose = FALSE)
  solve_perturbation(parts$model, parts$compiled, ssv$values, params,
                     verbose = FALSE)
}

#' Per-observable measurement-error variance as an n_obs x T matrix
#' @noRd
.fbt_me_extra <- function(me_variance, n_obs, TT) {
  v <- as.numeric(me_variance)
  if (length(v) == 0L || all(v == 0)) return(NULL)
  if (any(!is.finite(v)) || any(v < 0))
    stop("forecast_backtest: me_variance must be finite and non-negative.",
         call. = FALSE)
  matrix(rep_len(v, n_obs), n_obs, TT)
}

#' Filtered state mean and covariance for every period of `Y`
#'
#' The Kalman forward pass at period t uses only y_1..y_t, so a SINGLE pass
#' over the longest window supplies s_{t|t} and P_{t|t} for every origin t at
#' once. That is what `refit = "once"` exploits; `refit = "each"` re-runs the
#' pass per origin, which is arithmetically the same recursion truncated
#' earlier and therefore returns identical numbers.
#'
#' `Y` is in LEVELS: the state space carries the observation intercept and the
#' smoother subtracts it, and `.fbt_predictive_ss()` adds it back so the
#' predictive mean is scored against the realised level.
#'
#' @return List with `s` (T x n_state) and `P` (n_state x n_state x T).
#' @noRd
.fbt_filter <- function(ss, Y, me_variance) {
  me <- .fbt_me_extra(me_variance, ncol(Y), nrow(Y))
  sm <- .kalman_smoother_ss(Y, ss, me_extra = me)
  list(s = sm$filtered_states, P = sm$filtered_cov, loglik = sm$loglik)
}

#' Exact h-step Gaussian predictive moments from a filtered origin state
#'
#' See the derivation in this file's header.
#'
#' @param ss          dsge_ss object (lagged-state timing).
#' @param s0,P0       Filtered state mean / covariance at the origin.
#' @param horizons    Integer vector of horizons (>= 1).
#' @param me_variance Scalar or per-observable measurement-error variance.
#' @return Named list (by horizon) of `list(mean =, cov =)`.
#' @noRd
.fbt_predictive_ss <- function(ss, s0, P0, horizons, me_variance = 0) {
  T_mat <- ss$T_mat; R_mat <- ss$R_mat
  Z_mat <- ss$Z_mat; D_mat <- ss$D_mat
  Q <- ss$Sigma_e
  if (is.null(Q)) {
    warning("forecast_backtest: state space has no Sigma_e; using Q = I, ",
            "which assumes every shock has stderr 1.", call. = FALSE)
    Q <- diag(ss$n_shock)
  }
  RQR <- R_mat %*% Q %*% t(R_mat)
  DQD <- D_mat %*% Q %*% t(D_mat)
  me  <- rep_len(as.numeric(me_variance), ss$n_obs)

  h_max <- max(horizons)
  m <- as.numeric(s0)
  V <- as.matrix(P0)
  out <- vector("list", length(horizons))
  names(out) <- as.character(horizons)
  ## k indexes the state one period BEHIND the observation, matching
  ## y_{t+h} = Z s_{t+h-1} + D eps_{t+h}: horizon h reads (m, V) after
  ## h-1 state steps.
  for (k in 0:(h_max - 1L)) {
    if (k > 0L) {
      m <- as.numeric(T_mat %*% m)
      V <- T_mat %*% V %*% t(T_mat) + RQR
      V <- 0.5 * (V + t(V))
    }
    h <- k + 1L
    if (h %in% horizons) {
      ## + d puts the predictive back in LEVELS, which is what the realised
      ## y_{t+h} the scoring rules see is in. The state pass is in deviations
      ## because the filter subtracted the same intercept on the way in.
      mu <- as.numeric(Z_mat %*% m)
      if (!is.null(ss$d) && !anyNA(ss$d)) mu <- mu + as.numeric(ss$d)
      Sg <- Z_mat %*% V %*% t(Z_mat) + DQD
      Sg <- 0.5 * (Sg + t(Sg))
      diag(Sg) <- diag(Sg) + me
      names(mu) <- ss$obs_names
      dimnames(Sg) <- list(ss$obs_names, ss$obs_names)
      out[[as.character(h)]] <- list(mean = mu, cov = Sg)
    }
  }
  out
}

#' Naive reduced-form benchmark predictives (random walk / AR(1) / mean)
#'
#' All three are the SAME recursion on a VAR(1)-in-levels representation
#'   y_{t+1} = c + B y_t + u_{t+1},  u ~ N(0, Su)
#' so the h-step moments come from one loop:
#'   m_0 = y_t (or the sample mean), V_0 = 0 (or the sample covariance);
#'   m_k = c + B m_{k-1},  V_k = B V_{k-1} B' + Su.
#'
#' \describe{
#'   \item{rw}{c = 0, B = I, Su = cov(diff(train)) -- the random walk.}
#'   \item{ar1}{Per-variable OLS AR(1) with intercept; B = diag(b), Su the
#'     residual covariance (cross-equation correlation retained).}
#'   \item{mean}{The unconditional sample mean and covariance at every
#'     horizon.}
#' }
#'
#' @param kind     "rw", "ar1" or "mean".
#' @param train    Training window (t x n_obs), origin = last row.
#' @param horizons Integer vector of horizons.
#' @param me_variance Added to the predictive variance diagonal, for parity
#'   with the structural path.
#' @return Named list (by horizon) of `list(mean =, cov =)`.
#' @noRd
.fbt_predictive_naive <- function(kind, train, horizons, me_variance = 0) {
  Y  <- as.matrix(train)
  n  <- nrow(Y); k <- ncol(Y)
  nm <- colnames(Y)
  if (n < 4L)
    stop(sprintf(paste0("forecast_backtest: the \"%s\" benchmark needs at ",
                        "least 4 training observations, got %d."), kind, n),
         call. = FALSE)

  if (kind == "rw") {
    cvec <- rep(0, k)
    B    <- diag(k)
    Su   <- stats::cov(diff(Y))
    m    <- Y[n, ]
    V    <- matrix(0, k, k)
  } else if (kind == "ar1") {
    x <- Y[-n, , drop = FALSE]
    z <- Y[-1L, , drop = FALSE]
    b <- numeric(k); a <- numeric(k)
    resid <- matrix(0, n - 1L, k)
    for (j in seq_len(k)) {
      fit <- stats::lm.fit(cbind(1, x[, j]), z[, j])
      a[j] <- fit$coefficients[1L]
      b[j] <- fit$coefficients[2L]
      if (!is.finite(b[j])) { a[j] <- mean(z[, j]); b[j] <- 0 }
      resid[, j] <- z[, j] - (a[j] + b[j] * x[, j])
    }
    cvec <- a
    B    <- diag(b, nrow = k)
    Su   <- crossprod(resid) / max(1L, nrow(resid) - 2L)
    m    <- Y[n, ]
    V    <- matrix(0, k, k)
  } else {                                   # "mean"
    ## B = 0 collapses the recursion to m_h = cvec and V_h = Su at every
    ## horizon, i.e. the unconditional sample moments -- so Su, not V0, is
    ## where the sample covariance has to go.
    cvec <- colMeans(Y)
    B    <- matrix(0, k, k)
    Su   <- stats::cov(Y)
    m    <- colMeans(Y)
    V    <- matrix(0, k, k)
  }

  me    <- rep_len(as.numeric(me_variance), k)
  h_max <- max(horizons)
  out <- vector("list", length(horizons)); names(out) <- as.character(horizons)
  for (h in seq_len(h_max)) {
    m <- as.numeric(cvec + B %*% m)
    V <- B %*% V %*% t(B) + Su
    V <- 0.5 * (V + t(V))
    if (h %in% horizons) {
      mu <- m; Sg <- V
      diag(Sg) <- diag(Sg) + me
      names(mu) <- nm; dimnames(Sg) <- list(nm, nm)
      out[[as.character(h)]] <- list(mean = mu, cov = Sg)
    }
  }
  out
}

#' Draw an ensemble from Gaussian predictive moments
#' @noRd
.fbt_draw <- function(pm, n_draws) {
  k  <- length(pm$mean)
  Sg <- 0.5 * (pm$cov + t(pm$cov))
  ev <- eigen(Sg, symmetric = TRUE)
  L  <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = k)
  Z  <- matrix(stats::rnorm(n_draws * k), n_draws, k)
  X  <- Z %*% t(L)
  X  <- sweep(X, 2L, pm$mean, `+`)
  colnames(X) <- names(pm$mean)
  X
}

#' Pool a list of per-parameter-draw predictives into one ensemble
#' @noRd
.fbt_draw_mixture <- function(pm_list, n_draws) {
  M   <- length(pm_list)
  per <- max(2L, ceiling(n_draws / M))
  do.call(rbind, lapply(pm_list, .fbt_draw, n_draws = per))
}


## ---- Estimation at an origin --------------------------------------------

#' Estimate the model on `Y[1:origin, ]` and return the parameter vector(s)
#'
#' @return List with `params` (a full parameter vector) and, for
#'   `estimator = "posterior"`, `theta_draws` (a thinned matrix of posterior
#'   draws) -- otherwise `NULL`.
#' @noRd
.fbt_estimate <- function(parts, Y_train, obs_vars, estimator, me_variance,
                          n_iter, n_post_draws, mcmc_draws, mcmc_warmup,
                          verbose) {
  if (estimator == "fixed") return(list(params = parts$params, theta_draws = NULL))

  solved <- structure(list(model = parts$model, compiled = parts$compiled,
                           dr = parts$dr, params = parts$params),
                      class = c("dynhr_solved", "list"))
  mode_res <- run_mode_finding(solved, Y_train, obs_vars = obs_vars,
                               n_iter = n_iter, me_variance = me_variance,
                               verbose = verbose)
  params <- apply_theta_to_params(parts$model, mode_res$theta_mode,
                                  params = parts$params)
  if (estimator == "mode")
    return(list(params = params, theta_draws = NULL))

  ## estimator == "posterior"
  ch <- mcmc(mode_res$log_post_fn, mode_res$theta_mode, mode_res$Sigma_prop,
             n_draws = mcmc_draws, n_warmup = mcmc_warmup)
  dr_mat <- ch$chain
  idx <- unique(round(seq(1, nrow(dr_mat), length.out = n_post_draws)))
  list(params = params, theta_draws = dr_mat[idx, , drop = FALSE])
}


## ---- The driver ----------------------------------------------------------

#' Recursive out-of-sample forecast backtest
#'
#' Walks an expanding window of forecast origins.  At each origin it
#' (re-)estimates the model on the data available up to that point, forms the
#' \eqn{h}-step-ahead predictive distribution for every requested horizon,
#' and scores it against what actually happened.  The result is a tidy table
#' of (origin, horizon, variable, score, value) rows plus \code{summary()},
#' \code{print()} and \code{plot()} methods.
#'
#' @section Predictive distribution:
#' For a linear Gaussian state space the \eqn{h}-step predictive is exactly
#' Gaussian and is computed in closed form from the filtered origin state
#' \eqn{s_{t|t}} and its covariance \eqn{P_{t|t}} (see the recursion in the
#' source header of \code{R/forecast-eval.R}).  Keeping \eqn{P_{t|t}} matters:
#' forecasting from the point estimate of the terminal state alone makes the
#' predictive too sharp and the PITs mis-calibrated.  With \code{n_draws = 0}
#' (the default) the scores are evaluated on those exact moments, so CRPS, the
#' log score, the PITs and the coverage indicators carry NO Monte-Carlo noise.
#' Set \code{n_draws > 0} to score an ensemble drawn from the same
#' distribution instead -- required for the multivariate \code{"energy"} and
#' \code{"variogram"} rules, and used automatically for
#' \code{estimator = "posterior"}, whose predictive is a mixture over
#' parameter draws.
#'
#' @section Refit policy:
#' \code{refit = "each"} re-estimates at every origin (pseudo-real-time, and
#' the honest but expensive choice).  \code{refit = "once"} estimates only at
#' the FIRST origin and thereafter merely re-filters the growing sample, which
#' isolates the contribution of new data from the contribution of new
#' parameter estimates.  Because the Kalman forward pass at period \eqn{t}
#' uses only \eqn{y_1 \ldots y_t}, \code{"once"} takes a single filter pass
#' over the longest window and slices it -- arithmetically identical to
#' re-filtering per origin, so with \code{estimator = "fixed"} the two
#' policies return the same numbers.
#'
#' @section Benchmarks:
#' Passing \code{model = "rw"} (random walk), \code{"ar1"} (per-variable OLS
#' AR(1) with intercept) or \code{"mean"} (unconditional sample moments) runs
#' the same driver on a naive reduced-form benchmark.  Score the structural
#' model and the benchmark over the same origins and hand the two score series
#' to \code{\link{dm_test}}.
#'
#' @param model A \code{dynhr_solved} object (from \code{\link{solve_model}}),
#'   a list carrying \code{$model}/\code{$compiled}/\code{$dr}, or one of the
#'   benchmark names \code{"rw"}, \code{"ar1"}, \code{"mean"}.
#' @param data T \eqn{\times} n_obs matrix (or data.frame) of observations in
#'   \strong{levels} -- the model's steady state is subtracted when filtering
#'   and added back to the predictive mean, so the scores compare like with
#'   like. Must not contain \code{NA}.
#' @param obs_vars Character vector of observable names.  Defaults to
#'   \code{colnames(data)}.
#' @param origins Integer vector of forecast origins (indices into the rows of
#'   \code{data}; the origin is the LAST observed period).  Default
#'   \code{NULL} uses every origin from \code{min_train} to
#'   \code{nrow(data) - max(horizons)}.
#' @param horizons Integer vector of forecast horizons, default \code{1:8}.
#' @param estimator \code{"mode"} (posterior mode via
#'   \code{\link{run_mode_finding}}), \code{"posterior"} (mode, then
#'   \code{\link{mcmc}}, integrating the predictive over parameter draws), or
#'   \code{"fixed"} (no estimation -- use the parameters already in
#'   \code{model}, or \code{theta}).
#' @param refit \code{"each"} or \code{"once"}; see the section above.
#' @param theta Optional named parameter vector applied with
#'   \code{\link{apply_theta_to_params}} before anything else.  Supplying it
#'   with \code{estimator = "fixed"} pins every window to the same parameters.
#' @param me_variance Measurement-error variance (scalar or one per
#'   observable), used in both the filter and the predictive variance.
#' @param rules Scoring rules passed to \code{\link{score_forecast}}.  Default
#'   \code{c("crps", "logs", "pit", "coverage")}.
#' @param coverage_levels Nominal central-interval levels, default
#'   \code{c(0.5, 0.9)}.
#' @param n_draws Ensemble size per (origin, horizon).  \code{0} (default)
#'   scores the exact Gaussian moments instead of drawing.
#' @param n_post_draws Number of thinned posterior parameter draws forming the
#'   predictive mixture when \code{estimator = "posterior"} (default 25).
#' @param min_train Shortest training window when \code{origins} is
#'   \code{NULL}.  Default \code{NULL} uses half the sample (at least 10).
#' @param n_iter Mode-finder iteration cap forwarded to
#'   \code{\link{run_mode_finding}}.
#' @param mcmc_draws,mcmc_warmup Chain length and warm-up for
#'   \code{estimator = "posterior"}.
#' @param verbose Print per-origin progress (default \code{FALSE}).
#' @param ... Reserved; must be empty.
#' @return An object of class \code{"dynhr_backtest"}: a list with
#'   \describe{
#'     \item{\code{scores}}{Tidy data.frame with columns \code{origin},
#'       \code{horizon}, \code{variable}, \code{score}, \code{value}.}
#'     \item{\code{predictive}}{List indexed by \code{"origin.horizon"} of
#'       \code{list(mean =, cov =)}.}
#'     \item{\code{realised}}{Matrix of the realised values scored.}
#'     \item{\code{params}}{Per-origin parameter vectors actually used.}
#'     \item{\code{meta}}{The settings of the run.}
#'   }
#' @export
#'
#' @seealso \code{\link{score_forecast}}, \code{\link{dm_test}},
#'   \code{\link{conditional_forecast}}
#'
#' @examples
#' ## Random-walk benchmark on simulated data -- no model solve needed.
#' set.seed(1)
#' Y <- matrix(cumsum(rnorm(120)), 120, 1, dimnames = list(NULL, "y"))
#' bt <- forecast_backtest("rw", Y, obs_vars = "y",
#'                         horizons = c(1L, 4L), estimator = "fixed")
#' summary(bt)
forecast_backtest <- function(model,
                              data,
                              obs_vars = NULL,
                              origins  = NULL,
                              horizons = 1:8,
                              estimator = c("mode", "posterior", "fixed"),
                              refit     = c("each", "once"),
                              theta = NULL,
                              me_variance = 0,
                              rules = c("crps", "logs", "pit", "coverage"),
                              coverage_levels = c(0.5, 0.9),
                              n_draws = 0L,
                              n_post_draws = 25L,
                              min_train = NULL,
                              n_iter = 200L,
                              mcmc_draws = 500L,
                              mcmc_warmup = 250L,
                              verbose = FALSE,
                              ...) {

  if (length(list(...)))
    stop("forecast_backtest: unused argument(s): ",
         paste(names(list(...)), collapse = ", "), call. = FALSE)

  estimator <- match.arg(estimator)
  refit     <- match.arg(refit)

  dd       <- .fbt_data(data, obs_vars)
  Y        <- dd$Y
  obs_vars <- dd$obs_vars
  TT       <- nrow(Y)

  horizons <- sort(unique(as.integer(horizons)))
  if (any(is.na(horizons)) || any(horizons < 1L))
    stop("forecast_backtest: horizons must be positive integers.", call. = FALSE)
  h_max <- max(horizons)

  is_naive <- is.character(model) && length(model) == 1L
  if (is_naive) {
    model <- match.arg(model, c("rw", "ar1", "mean"))
    if (estimator != "fixed") estimator <- "fixed"
  }

  ## ---- Origins ----------------------------------------------------------
  if (is.null(origins)) {
    if (is.null(min_train)) min_train <- max(10L, ceiling(TT / 2))
    min_train <- as.integer(min_train)
    origins <- seq.int(min_train, TT - h_max)
  } else {
    origins <- sort(unique(as.integer(origins)))
  }
  if (length(origins) == 0L || any(origins < 2L) || any(origins + h_max > TT))
    stop(sprintf(paste0("forecast_backtest: no usable forecast origin. With ",
                        "T = %d and max(horizons) = %d, origins must lie in ",
                        "[2, %d]."), TT, h_max, TT - h_max), call. = FALSE)

  n_draws <- as.integer(n_draws)
  needs_draws <- any(c("energy", "variogram") %in% rules) ||
    estimator == "posterior"
  if (needs_draws && n_draws < 2L) n_draws <- 1000L

  ## ---- Structural set-up ------------------------------------------------
  parts <- NULL; ss_cached <- NULL; filt_cached <- NULL
  if (!is_naive) {
    parts <- .fbt_parts(model)
    if (!is.null(theta))
      parts$params <- apply_theta_to_params(parts$model, theta,
                                            params = parts$params)
    if (!is.null(theta) && estimator == "fixed")
      parts$dr <- .fbt_resolve(parts, parts$params)
  } else if (!is.null(theta)) {
    stop("forecast_backtest: 'theta' is meaningless for a naive benchmark.",
         call. = FALSE)
  }

  ## `refit = "once"`: estimate at the FIRST origin only, then reuse the state
  ## space and take ONE filter pass over the longest window.
  params_once <- NULL; theta_draws_once <- NULL
  if (!is_naive && refit == "once") {
    est <- .fbt_estimate(parts, Y[seq_len(origins[1L]), , drop = FALSE],
                         obs_vars, estimator, me_variance, n_iter,
                         n_post_draws, mcmc_draws, mcmc_warmup, verbose)
    params_once      <- est$params
    theta_draws_once <- est$theta_draws
    dr_once   <- if (estimator == "fixed") parts$dr
                 else .fbt_resolve(parts, params_once)
    ss_cached <- build_dsge_state_space(parts$model, dr_once, obs_vars,
                                        verbose = FALSE, params = params_once)
    filt_cached <- .fbt_filter(ss_cached,
                               Y[seq_len(max(origins)), , drop = FALSE],
                               me_variance)
  }

  ## ---- Walk the origins -------------------------------------------------
  rows <- vector("list", length(origins))
  pred_store <- list()
  params_store <- list()

  for (i in seq_along(origins)) {
    t0    <- origins[i]
    train <- Y[seq_len(t0), , drop = FALSE]
    if (verbose)
      message(sprintf("forecast_backtest: origin %d/%d (t = %d)",
                      i, length(origins), t0))

    if (is_naive) {
      pm_h <- .fbt_predictive_naive(model, train, horizons, me_variance)
      pm_mix <- NULL
      params_store[[as.character(t0)]] <- NA
    } else if (refit == "once") {
      pm_h <- .fbt_predictive_ss(ss_cached, filt_cached$s[t0, ],
                                 filt_cached$P[, , t0], horizons, me_variance)
      pm_mix <- if (is.null(theta_draws_once)) NULL
                else .fbt_mixture(parts, theta_draws_once, train, obs_vars,
                                  horizons, me_variance)
      params_store[[as.character(t0)]] <- params_once
    } else {
      est <- .fbt_estimate(parts, train, obs_vars, estimator, me_variance,
                           n_iter, n_post_draws, mcmc_draws, mcmc_warmup,
                           verbose)
      dr_t <- if (estimator == "fixed") parts$dr
              else .fbt_resolve(parts, est$params)
      ss_t <- build_dsge_state_space(parts$model, dr_t, obs_vars,
                                     verbose = FALSE, params = est$params)
      ft   <- .fbt_filter(ss_t, train, me_variance)
      pm_h <- .fbt_predictive_ss(ss_t, ft$s[t0, ], ft$P[, , t0],
                                 horizons, me_variance)
      pm_mix <- if (is.null(est$theta_draws)) NULL
                else .fbt_mixture(parts, est$theta_draws, train, obs_vars,
                                  horizons, me_variance)
      params_store[[as.character(t0)]] <- est$params
    }

    rows[[i]] <- .fbt_score_origin(t0, horizons, pm_h, pm_mix, Y, obs_vars,
                                   rules, coverage_levels, n_draws)
    for (h in horizons)
      pred_store[[sprintf("%d.%d", t0, h)]] <- pm_h[[as.character(h)]]
  }

  scores <- do.call(rbind, rows)
  rownames(scores) <- NULL

  structure(list(
    scores     = scores,
    predictive = pred_store,
    realised   = Y,
    params     = params_store,
    meta = list(origins = origins, horizons = horizons, obs_vars = obs_vars,
                estimator = estimator, refit = refit, rules = rules,
                coverage_levels = coverage_levels, n_draws = n_draws,
                me_variance = me_variance,
                model = if (is_naive) model else "dsge",
                n_origins = length(origins))),
    class = c("dynhr_backtest", "list"))
}

#' Predictive mixture over posterior parameter draws
#' @noRd
.fbt_mixture <- function(parts, theta_draws, train, obs_vars, horizons,
                         me_variance) {
  t0 <- nrow(train)
  out <- vector("list", nrow(theta_draws))
  keep <- logical(nrow(theta_draws))
  for (m in seq_len(nrow(theta_draws))) {
    th <- theta_draws[m, ]
    names(th) <- colnames(theta_draws)
    pr <- tryCatch({
      pp   <- apply_theta_to_params(parts$model, th, params = parts$params)
      dr_m <- .fbt_resolve(parts, pp)
      ss_m <- build_dsge_state_space(parts$model, dr_m, obs_vars,
                                     verbose = FALSE, params = pp)
      fm   <- .fbt_filter(ss_m, train, me_variance)
      .fbt_predictive_ss(ss_m, fm$s[t0, ], fm$P[, , t0], horizons, me_variance)
    }, error = function(e) NULL)
    if (!is.null(pr)) { out[[m]] <- pr; keep[m] <- TRUE }
  }
  out <- out[keep]
  if (length(out) == 0L)
    stop("forecast_backtest: every posterior draw failed to solve; the ",
         "predictive mixture is empty.", call. = FALSE)
  ## Re-index: horizon -> list of per-draw moments.
  res <- vector("list", length(horizons))
  names(res) <- as.character(horizons)
  for (h in horizons)
    res[[as.character(h)]] <- lapply(out, function(p) p[[as.character(h)]])
  res
}

#' Score one origin across all horizons; returns a tidy data.frame
#' @noRd
.fbt_score_origin <- function(t0, horizons, pm_h, pm_mix, Y, obs_vars,
                              rules, coverage_levels, n_draws) {
  out <- vector("list", length(horizons))
  for (idx in seq_along(horizons)) {
    h  <- horizons[idx]
    y  <- Y[t0 + h, ]
    hk <- as.character(h)
    if (!is.null(pm_mix)) {
      X  <- .fbt_draw_mixture(pm_mix[[hk]], n_draws)
      sc <- score_forecast(X, y, rules = rules, obs_vars = obs_vars,
                           coverage_levels = coverage_levels)
    } else if (n_draws >= 2L) {
      X  <- .fbt_draw(pm_h[[hk]], n_draws)
      sc <- score_forecast(X, y, rules = rules, obs_vars = obs_vars,
                           coverage_levels = coverage_levels)
    } else {
      sc <- score_forecast(NULL, y, rules = rules, obs_vars = obs_vars,
                           predictive_moments = pm_h[[hk]],
                           coverage_levels = coverage_levels)
    }
    out[[idx]] <- .fbt_tidy_scores(t0, h, sc, obs_vars)
  }
  do.call(rbind, out)
}

#' Flatten one `dynhr_forecast_scores` into tidy rows
#' @noRd
.fbt_tidy_scores <- function(t0, h, sc, obs_vars) {
  mk <- function(variable, score, value)
    data.frame(origin = as.integer(t0), horizon = as.integer(h),
               variable = variable, score = score, value = as.numeric(value),
               stringsAsFactors = FALSE)
  parts <- list()
  if (!is.null(sc$crps))
    parts[[length(parts) + 1L]] <- mk(obs_vars, "crps", sc$crps)
  if (!is.null(sc$logs)) {
    parts[[length(parts) + 1L]] <- mk(obs_vars, "logs", sc$logs)
    if (!is.null(sc$logs_joint) && is.finite(sc$logs_joint))
      parts[[length(parts) + 1L]] <- mk("(joint)", "logs_joint", sc$logs_joint)
  }
  if (!is.null(sc$pit))
    parts[[length(parts) + 1L]] <- mk(obs_vars, "pit", sc$pit)
  if (!is.null(sc$energy))
    parts[[length(parts) + 1L]] <- mk("(joint)", "energy", sc$energy)
  if (!is.null(sc$variogram) && is.finite(sc$variogram))
    parts[[length(parts) + 1L]] <- mk("(joint)", "variogram", sc$variogram)
  if (!is.null(sc$coverage)) {
    cv <- sc$coverage
    for (l in seq_along(cv$levels))
      parts[[length(parts) + 1L]] <-
        mk(obs_vars, sprintf("coverage_%g", 100 * cv$levels[l]),
           cv$covered[l, ])
  }
  do.call(rbind, parts)
}


## ---- Methods -------------------------------------------------------------

#' Summarise a recursive forecast backtest
#'
#' Averages the per-origin scores by horizon and variable, turns the coverage
#' indicators into empirical coverage rates against their nominal levels, and
#' runs a Kolmogorov-Smirnov uniformity test on the PITs (the standard density
#' calibration check: a small p-value says the predictive is mis-calibrated,
#' most often over-confident).
#'
#' @param object A \code{dynhr_backtest} from \code{\link{forecast_backtest}}.
#' @param ... Ignored.
#' @return An object of class \code{"summary.dynhr_backtest"}: a list with
#'   \code{scores} (mean score by horizon/variable/rule), \code{coverage}
#'   (empirical vs nominal rates) and \code{pit} (KS statistic and p-value
#'   per horizon and variable).
#' @export
summary.dynhr_backtest <- function(object, ...) {
  s <- object$scores
  is_cov <- grepl("^coverage_", s$score)

  mean_tab <- NULL
  sm <- s[!is_cov & s$score != "pit", , drop = FALSE]
  if (nrow(sm)) {
    ag <- stats::aggregate(list(mean = sm$value),
                           by = list(horizon = sm$horizon,
                                     variable = sm$variable,
                                     score = sm$score),
                           FUN = mean)
    n_ag <- stats::aggregate(list(n = sm$value),
                             by = list(horizon = sm$horizon,
                                       variable = sm$variable,
                                       score = sm$score),
                             FUN = length)
    ag$n <- n_ag$n
    mean_tab <- ag[order(ag$score, ag$variable, ag$horizon), ]
    rownames(mean_tab) <- NULL
  }

  cov_tab <- NULL
  sc <- s[is_cov, , drop = FALSE]
  if (nrow(sc)) {
    ag <- stats::aggregate(list(empirical = sc$value),
                           by = list(horizon = sc$horizon,
                                     variable = sc$variable,
                                     score = sc$score),
                           FUN = mean)
    ag$nominal <- as.numeric(sub("^coverage_", "", ag$score)) / 100
    ag$n <- stats::aggregate(list(n = sc$value),
                             by = list(horizon = sc$horizon,
                                       variable = sc$variable,
                                       score = sc$score),
                             FUN = length)$n
    ag <- ag[, c("horizon", "variable", "nominal", "empirical", "n")]
    cov_tab <- ag[order(ag$nominal, ag$variable, ag$horizon), ]
    rownames(cov_tab) <- NULL
  }

  pit_tab <- NULL
  sp <- s[s$score == "pit", , drop = FALSE]
  if (nrow(sp)) {
    key <- unique(sp[, c("horizon", "variable")])
    res <- lapply(seq_len(nrow(key)), function(i) {
      u <- sp$value[sp$horizon == key$horizon[i] &
                      sp$variable == key$variable[i]]
      kt <- suppressWarnings(stats::ks.test(u, "punif"))
      data.frame(horizon = key$horizon[i], variable = key$variable[i],
                 n = length(u), mean_pit = mean(u),
                 ks_stat = unname(kt$statistic), ks_p = unname(kt$p.value),
                 stringsAsFactors = FALSE)
    })
    pit_tab <- do.call(rbind, res)
    pit_tab <- pit_tab[order(pit_tab$variable, pit_tab$horizon), ]
    rownames(pit_tab) <- NULL
  }

  structure(list(scores = mean_tab, coverage = cov_tab, pit = pit_tab,
                 meta = object$meta),
            class = c("summary.dynhr_backtest", "list"))
}

#' @param x A \code{summary.dynhr_backtest} object.
#' @param ... Ignored.
#' @rdname summary.dynhr_backtest
#' @export
print.summary.dynhr_backtest <- function(x, ...) {
  m <- x$meta
  cat(sprintf("dynhr forecast backtest summary\n"))
  cat(sprintf("  model: %s | estimator: %s | refit: %s\n",
              m$model, m$estimator, m$refit))
  cat(sprintf("  %d origins (%d..%d), horizons %s, %d observable(s)\n",
              m$n_origins, min(m$origins), max(m$origins),
              paste(m$horizons, collapse = ","), length(m$obs_vars)))
  cat(sprintf("  predictive: %s\n",
              if (m$n_draws >= 2L) sprintf("%d draws", m$n_draws)
              else "exact Gaussian moments"))
  if (!is.null(x$scores)) {
    cat("\n-- Mean score by horizon (lower is better) --\n")
    print(x$scores, row.names = FALSE)
  }
  if (!is.null(x$coverage)) {
    cat("\n-- Interval coverage (empirical vs nominal) --\n")
    print(x$coverage, row.names = FALSE)
  }
  if (!is.null(x$pit)) {
    cat("\n-- PIT uniformity (KS test; small p = mis-calibrated) --\n")
    print(x$pit, row.names = FALSE)
  }
  invisible(x)
}

#' Print a recursive forecast backtest
#'
#' @param x A \code{dynhr_backtest} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.dynhr_backtest <- function(x, ...) {
  m <- x$meta
  cat(sprintf("<dynhr_backtest: %s, %d origins x %d horizons x %d obs>\n",
              m$model, m$n_origins, length(m$horizons), length(m$obs_vars)))
  cat(sprintf("  estimator = %s, refit = %s, rules = %s\n",
              m$estimator, m$refit, paste(m$rules, collapse = ",")))
  cat(sprintf("  %d score rows; use summary() for the tables\n",
              nrow(x$scores)))
  invisible(x)
}

#' Plot a recursive forecast backtest
#'
#' Two base-graphics panels: a PIT histogram with the uniform reference line
#' (bars systematically above the line in the middle mean the predictive is
#' too WIDE; a U shape means it is too NARROW), and mean score by horizon,
#' one line per observable.
#'
#' @param x A \code{dynhr_backtest} object.
#' @param score Which score to plot by horizon; default the first of
#'   \code{"crps"}, \code{"logs"}, \code{"energy"} present.
#' @param pit_horizon Horizon whose PITs feed the histogram; default the
#'   shortest.
#' @param n_bins Number of PIT histogram bins (default 10).
#' @param ... Passed to the underlying \code{plot()} calls.
#' @return \code{x}, invisibly.
#' @export
plot.dynhr_backtest <- function(x, score = NULL, pit_horizon = NULL,
                                n_bins = 10L, ...) {
  s <- x$scores
  if (is.null(score)) {
    avail <- intersect(c("crps", "logs", "energy"), unique(s$score))
    if (length(avail) == 0L)
      stop("plot.dynhr_backtest: no plottable score in this backtest.",
           call. = FALSE)
    score <- avail[1L]
  }
  has_pit <- any(s$score == "pit")
  op <- graphics::par(mfrow = if (has_pit) c(1L, 2L) else c(1L, 1L))
  on.exit(graphics::par(op), add = TRUE)

  if (has_pit) {
    if (is.null(pit_horizon)) pit_horizon <- min(x$meta$horizons)
    u <- s$value[s$score == "pit" & s$horizon == pit_horizon]
    graphics::hist(u, breaks = seq(0, 1, length.out = n_bins + 1L),
                   freq = FALSE, col = "grey85", border = "white",
                   xlab = "PIT", main = sprintf("PIT, h = %d", pit_horizon),
                   ...)
    graphics::abline(h = 1, lty = 2)
  }

  ss <- s[s$score == score, , drop = FALSE]
  ag <- stats::aggregate(list(mean = ss$value),
                         by = list(horizon = ss$horizon,
                                   variable = ss$variable), FUN = mean)
  vars <- unique(ag$variable)
  ylim <- range(ag$mean, finite = TRUE)
  plot(NA, xlim = range(ag$horizon), ylim = ylim, xlab = "horizon",
       ylab = sprintf("mean %s", score),
       main = sprintf("%s by horizon", score), ...)
  for (j in seq_along(vars)) {
    a <- ag[ag$variable == vars[j], , drop = FALSE]
    a <- a[order(a$horizon), ]
    graphics::lines(a$horizon, a$mean, type = "b", pch = 19, col = j, lty = j)
  }
  if (length(vars) > 1L)
    graphics::legend("topleft", legend = vars, col = seq_along(vars),
                     lty = seq_along(vars), pch = 19, bty = "n", cex = 0.8)
  invisible(x)
}
