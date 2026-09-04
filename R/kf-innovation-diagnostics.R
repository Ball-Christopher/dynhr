## R/kf-innovation-diagnostics.R
## --------------------------------------------------------------------------
## Innovation-whiteness diagnostic for the Kalman filter.
##
## Theorem: at the TRUE model (and true parameters), the standardized
## one-step-ahead prediction errors (innovations) of a correctly specified
## Kalman filter are white noise with unit variance:
##   z_{i,t} = v_{i,t} / sqrt(F_{ii,t})  ~ iid N(0, 1)   (per observable i)
## This is a cheap, decisive, model-free oracle: any likelihood-evaluation
## bug (wrong Sigma_e, wrong ZZ/DD, a mis-routed me_extra/shock_scale column,
## a stale steady state, ...) tends to show up as either a standardized
## variance far from 1 or a nonzero lag-1 autocorrelation, well before it is
## visible in the loglik level itself.
##
## kalman_filter() (R/kalman-filter.R) does not return per-step innovations
## or forecast variances to the caller on any path -- only `filtered_states`
## and (optionally) `loglik_contrib`. Rather than edit that file, this
## diagnostic runs its own thin per-step Kalman recursion (mirroring the
## "dare"/.kf_step textbook formulas) so it can capture v_t and F_t directly.
## It reuses the package's own model-setup helpers (.get_shock_cov,
## solve_lyapunov) so it always agrees with kalman_filter()'s conventions.
## --------------------------------------------------------------------------

#' Kalman-filter innovation whiteness diagnostics
#'
#' At the true model and true parameters, the one-step-ahead Kalman-filter
#' innovations, standardized by their forecast standard deviation, are white
#' noise with unit variance. This function runs a textbook (per-step
#' Riccati) Kalman filter, forms the standardized innovations
#' \eqn{z_{i,t} = v_{i,t} / \sqrt{F_{ii,t}}} for each observable \eqn{i}, and
#' reports two per-observable diagnostics against their asymptotic null:
#' \itemize{
#'   \item the sample variance of \eqn{z_{i,\cdot}} vs. 1, using the
#'     asymptotic standard deviation \eqn{\sqrt{2/T_i}} of a chi-square-based
#'     variance estimator;
#'   \item the lag-1 sample autocorrelation of \eqn{z_{i,\cdot}} vs. 0, using
#'     the standard Bartlett-type standard deviation \eqn{1/\sqrt{T_i}}.
#' }
#' A joint (portmanteau-style) summary combining all observables is also
#' returned. This is a cheap, model-free correctness oracle: likelihood bugs
#' (wrong shock covariance, mis-routed per-period tunes, a stale steady
#' state, ...) typically produce a standardized variance far from 1 or a
#' nonzero lag-1 autocorrelation well before they are visible in the raw
#' log-likelihood.
#'
#' @param data observation matrix (\code{n_obs x T}), or a data frame /
#'   numeric vector coercible to one; \code{NA} marks a missing observation
#'   at that observable/period (see Details).
#' @param dr decision rule (output of \code{\link{solve_perturbation}}).
#' @param model compiled model object (output of \code{\link{compile_model}}).
#' @param params named numeric vector of parameter values.
#' @param obs_vars character vector of observed variable names.
#' @param lik_init character; \code{"stationary"} (default) initializes
#'   \code{P0} via the discrete Lyapunov equation (matching
#'   \code{\link{kalman_filter}}'s \code{lik_init = "stationary"} /
#'   \code{"auto"} on a stable model). \code{"diffuse"} is not supported by
#'   this thin diagnostic; use a burn-in and stick with \code{"stationary"}
#'   on near-unit-root models, or pre-filter with \code{\link{kalman_filter}}
#'   and inspect its residuals directly.
#' @param me_variance scalar measurement-error jitter added to the
#'   innovation covariance diagonal, matching \code{\link{kalman_filter}}'s
#'   \code{me_variance} convention (default \code{0}).
#'
#' @details
#' \strong{Missing data.} When \code{Y[i, t]} is \code{NA}, observable
#' \code{i} is simply skipped at period \code{t}: no innovation is formed
#' for it and the Kalman update at \code{t} uses only the non-missing rows
#' of \code{Y[, t]} (the standard multivariate missing-data treatment: drop
#' rows of \code{ZZ}/\code{DD} and skip the corresponding entries of
#' \code{v}/\code{F} for that observable at that period). Each observable's
#' whiteness statistics are then computed only over its own non-missing
#' periods (its own effective \eqn{T_i}), and the lag-1 autocorrelation
#' skips any pair that straddles a gap (i.e. it is computed only over
#' consecutive non-missing periods, not periods \code{t} apart in wall-clock
#' time). This is a diagnostic simplification, not the statistically exact
#' univariate/sequential treatment that \code{\link{kalman_filter}} uses
#' internally for irregular missingness; for data with substantial gaps
#' prefer running \code{\link{kalman_filter}} directly and treat this
#' diagnostic's output as indicative.
#'
#' @return A list with class \code{"kf_innovation_diagnostics"}:
#' \describe{
#'   \item{by_obs}{data.frame, one row per observable, with columns
#'     \code{obs_var}, \code{n_used} (non-missing periods), \code{var_z}
#'     (sample variance of standardized innovations), \code{z_var}
#'     (z-stat of \code{var_z} against the null of 1, using sd
#'     \eqn{\sqrt{2/n\_used}}), \code{acf1} (lag-1 autocorrelation of
#'     standardized innovations), \code{z_acf1} (z-stat of \code{acf1}
#'     against 0, using sd \eqn{1/\sqrt{n\_used}}).}
#'   \item{joint}{list with \code{max_abs_z} (largest \code{|z_var|} or
#'     \code{|z_acf1|} across all observables) and \code{n_flagged} (count
#'     of individual z-stats with \code{|z| > 4}, a conservative multiple-
#'     comparisons-aware threshold).}
#'   \item{z}{the \code{n_obs x T} matrix of standardized innovations
#'     (\code{NA} where \code{data} was missing).}
#' }
#'
#' @examples
#' \donttest{
#' txt <- "
#'   var c k;
#'   varexo eps eta;
#'   parameters rho;
#'   rho = 0.9;
#'   model;
#'     c = rho*c(-1) + eps;
#'     k = 0.5*k(-1) + eta;
#'   end;
#'   initval; c = 0; k = 0; end;
#'   shocks; var eps; stderr 0.1; var eta; stderr 0.1; end;
#' "
#' m  <- parse_mod(txt, verbose = FALSE)
#' cm <- compile_model(m, verbose = FALSE)
#' ss <- solve_steady(cm, m$param_values, endo_names = m$var_names,
#'                    exo_names = m$varexo_names, verbose = FALSE)
#' dr <- solve_perturbation(m, cm, ss$values, m$param_values, verbose = FALSE)
#' set.seed(1)
#' sim <- simulate_model(dr, n_periods = 500, model = m, burn_in = 50)
#' Y <- t(as.matrix(sim[, c("c", "k")]))
#' diag <- kf_innovation_diagnostics(Y, dr, m, m$param_values,
#'                                   obs_vars = c("c", "k"))
#' diag$by_obs
#' }
#' @export
kf_innovation_diagnostics <- function(data, dr, model, params, obs_vars,
                                      lik_init = c("stationary", "diffuse"),
                                      me_variance = 0) {
  lik_init <- match.arg(lik_init)
  if (identical(lik_init, "diffuse"))
    stop("kf_innovation_diagnostics: lik_init = \"diffuse\" is not ",
         "supported by this thin diagnostic (no exact-diffuse phase is ",
         "implemented here). Use lik_init = \"stationary\", or run ",
         "kalman_filter() directly and inspect its residuals.",
         call. = FALSE)

  state_idx <- dr$state_idx
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  n_state   <- length(state_idx)
  n_obs     <- length(obs_vars)

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("kf_innovation_diagnostics: observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "), call. = FALSE)

  ghx <- dr$ghx; ghu <- dr$ghu
  TT  <- ghx[state_idx, , drop = FALSE]
  RR  <- ghu[state_idx, , drop = FALSE]
  ZZ  <- ghx[obs_idx,   , drop = FALSE]
  DD  <- ghu[obs_idx,   , drop = FALSE]
  d   <- dr$ys[obs_vars]

  Sigma_e <- .get_shock_cov(model, exo, params)
  QQ      <- tcrossprod(RR %*% Sigma_e, RR)
  HH      <- tcrossprod(DD %*% Sigma_e, DD)
  SS      <- RR %*% Sigma_e %*% t(DD)
  me_diag <- me_variance * diag(n_obs)

  if (is.null(dim(data))) data <- matrix(data, nrow = n_obs)
  data <- as.matrix(data)
  if (nrow(data) != n_obs) data <- t(data)
  if (nrow(data) != n_obs)
    stop(sprintf("kf_innovation_diagnostics: Y must have %d rows/cols ",
                 "matching obs_vars; got %d x %d.",
                 n_obs, nrow(data), ncol(data)), call. = FALSE)
  n_T <- ncol(data)

  Y_minus_d <- data - d

  P0 <- solve_lyapunov(TT, QQ)
  if (anyNA(P0))
    stop("kf_innovation_diagnostics: solve_lyapunov() returned NaN -- TT ",
         "has unit-root eigenvalues; lik_init = \"stationary\" is not valid ",
         "for this model.", call. = FALSE)

  s <- numeric(n_state)
  P <- P0

  ## Standardized innovations, NA where the observable is missing at t.
  z_mat <- matrix(NA_real_, n_obs, n_T)

  for (t in seq_len(n_T)) {
    obs_ok <- which(!is.na(Y_minus_d[, t]))

    ## Prediction step (state only; shared regardless of which rows are
    ## observed this period).
    v_full <- Y_minus_d[, t] - as.numeric(ZZ %*% s)

    if (length(obs_ok) == 0L) {
      ## Nothing observed this period: pure prediction, no update.
      s <- as.numeric(TT %*% s)
      P <- tcrossprod(TT %*% P, TT) + QQ
      P <- (P + t(P)) * 0.5
      next
    }

    Zo  <- ZZ[obs_ok, , drop = FALSE]
    Do  <- DD[obs_ok, , drop = FALSE]
    v_o <- v_full[obs_ok]

    PZo <- P %*% t(Zo)
    HHo <- tcrossprod(Do %*% Sigma_e, Do) + me_diag[obs_ok, obs_ok, drop = FALSE]
    Fo  <- Zo %*% PZo + HHo
    Fo  <- (Fo + t(Fo)) * 0.5

    Fc <- tryCatch(chol(Fo), error = function(e) NULL)
    if (is.null(Fc))
      stop("kf_innovation_diagnostics: non-positive-definite innovation ",
           "covariance at period ", t, "; cannot standardize innovations.",
           call. = FALSE)
    Fi <- chol2inv(Fc)

    ## Standardized innovations: z_i = v_i / sqrt(F_ii). Uses the marginal
    ## variance F_ii (diagonal), NOT the full whitening transform F^{-1/2}v,
    ## so that per-observable statistics stay interpretable one series at a
    ## time (matching the CLAUDE.md per-observable convention); cross-
    ## observable correlation in F is not tested by this function.
    z_mat[obs_ok, t] <- v_o / sqrt(diag(Fo))

    SSo  <- RR %*% Sigma_e %*% t(Do)
    K_o  <- (TT %*% PZo + SSo) %*% Fi
    s    <- as.numeric(TT %*% s) + drop(K_o %*% v_o)
    TmKZ <- TT - K_o %*% Zo
    RmKD <- RR - K_o %*% Do
    P    <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Sigma_e, RmKD)
    P    <- (P + t(P)) * 0.5
  }

  by_obs <- data.frame(
    obs_var = obs_vars,
    n_used  = integer(n_obs),
    var_z   = numeric(n_obs),
    z_var   = numeric(n_obs),
    acf1    = numeric(n_obs),
    z_acf1  = numeric(n_obs),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_obs)) {
    zi <- z_mat[i, ]
    ok <- which(!is.na(zi))
    n_i <- length(ok)
    by_obs$n_used[i] <- n_i

    if (n_i < 2L) {
      by_obs$var_z[i]  <- NA_real_
      by_obs$z_var[i]  <- NA_real_
      by_obs$acf1[i]   <- NA_real_
      by_obs$z_acf1[i] <- NA_real_
      next
    }

    zi_ok <- zi[ok]
    var_i <- stats::var(zi_ok)
    by_obs$var_z[i] <- var_i
    by_obs$z_var[i] <- (var_i - 1) / sqrt(2 / n_i)

    ## Lag-1 autocorrelation over CONSECUTIVE non-missing periods only
    ## (skip any pair straddling a gap; see Details).
    consec <- ok[which(diff(ok) == 1L)]
    if (length(consec) >= 2L) {
      x0 <- zi[consec]
      x1 <- zi[consec + 1L]
      n1 <- length(consec)
      ## Standard (biased, mean-0-known) lag-1 autocorrelation estimator.
      acf1 <- sum(x0 * x1) / sum(zi_ok^2)
      by_obs$acf1[i]   <- acf1
      by_obs$z_acf1[i] <- acf1 / sqrt(1 / n1)
    } else {
      by_obs$acf1[i]   <- NA_real_
      by_obs$z_acf1[i] <- NA_real_
    }
  }

  all_z <- c(by_obs$z_var, by_obs$z_acf1)
  all_z <- all_z[is.finite(all_z)]
  joint <- list(
    max_abs_z = if (length(all_z)) max(abs(all_z)) else NA_real_,
    n_flagged = if (length(all_z)) sum(abs(all_z) > 4) else NA_integer_
  )

  structure(
    list(by_obs = by_obs, joint = joint, z = z_mat),
    class = "kf_innovation_diagnostics"
  )
}

#' @export
print.kf_innovation_diagnostics <- function(x, ...) {
  cat("KF innovation whiteness diagnostics\n")
  print(x$by_obs, row.names = FALSE)
  cat(sprintf("\nJoint: max |z| = %.3f, %d / %d z-stats flagged (|z| > 4)\n",
              x$joint$max_abs_z, x$joint$n_flagged, 2L * nrow(x$by_obs)))
  invisible(x)
}
