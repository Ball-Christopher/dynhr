## R/diag-deep-d35-softness.R
## --------------------------------------------------------------------------
## D35. Misspecification softness (sandwich / information-matrix equality).
##
## "Robust?" is the fourth question a parameter must answer to be called deep:
## would the estimate survive a small, plausible misspecification of the model?
## A parameter whose apparent precision evaporates once we stop assuming the
## model is exactly right is "soft" -- an artefact of the specification, not a
## structural fact.
##
## We operationalise this through the information-matrix equality. Under correct
## specification the observed information H = -d2 loglik / dtheta2 (the "bread")
## equals the variance of the summed score Omega (the "meat"); White (1982).
## For time series the meat is the Bartlett HAC long-run covariance of the
## per-period scores (Andrews 1991 bandwidth), so that misspecified DYNAMICS
## (serially correlated scores) also show up; it reduces to the outer product
## sum_t s_t s_t' when the scores are a martingale difference sequence. Their
## divergence is misspecification, and the misspecification-robust
## ("sandwich") covariance is
##     V_sand = H^{-1} Omega H^{-1}        (vs the nominal V = H^{-1}).
## The per-parameter SOFTNESS is the ratio of the robust to the nominal
## standard error,
##     rho_j = se_robust_j / se_nominal_j ,
## which is the Andrews-Mikusheva (2014) "discrepancy between two estimates of
## Fisher information" read parameter-by-parameter. rho_j ~ 1 = firm;
## rho_j >> 1 = soft. (rho_j < 1 is also a failure of the IM equality, e.g.
## thin tails, but it makes the nominal band conservative, so it is not
## flagged.)
##
## Input: a per-period log-likelihood function `loglik_contrib_fn(theta)` that
## returns the length-T vector of contributions (the Kalman prediction-error
## decomposition). D35 finite-differences it for the scores and the Hessian,
## so it is decoupled from how the likelihood is computed.
##
## Exposes `$result$passport_axis` (named logical, TRUE = robust/firm) for the
## Deep-Parameter Passport's "robust" column.
##
## References:
##   White, H. (1982). Maximum likelihood estimation of misspecified models.
##     Econometrica, 50(1), 1-25.
##   Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation
##     consistent covariance matrix estimation. Econometrica, 59(3), 817-858.
##   Andrews, I., & Mikusheva, A. (2014). Weak identification in maximum
##     likelihood: a question of information. AER P&P, 104(5), 195-199.
##   Bonhomme, S., & Weidner, M. (2022). Minimizing sensitivity to model
##     misspecification. Quantitative Economics, 13(3), 907-954.
## --------------------------------------------------------------------------


# Per-period score matrix S (T x k) by central finite differences of the
# log-likelihood contributions: S[t, j] = d l_t / d theta_j (the sign is that
# of the log-likelihood, so colSums(S) = 0 at an interior MLE). The step is
# relative, h_j = max(|theta_j|, 1e-3) * eps. A column whose perturbed
# evaluations are non-finite (or of the wrong length) is retried once with half
# the step and otherwise left NA. Errors raised by loglik_contrib_fn propagate:
# a callback that cannot evaluate a point should return non-finite values.
.d35_score_matrix <- function(loglik_contrib_fn, theta, eps = 1e-4) {
  k  <- length(theta)
  h0 <- pmax(abs(theta), 1e-3) * eps
  ll0 <- loglik_contrib_fn(theta)
  Tn <- length(ll0)
  S  <- matrix(NA_real_, Tn, k)
  for (j in seq_len(k)) {
    hj <- h0[j]
    for (attempt in 1:2) {
      tp <- theta; tp[j] <- tp[j] + hj
      tm <- theta; tm[j] <- tm[j] - hj
      lp <- loglik_contrib_fn(tp)
      lm <- loglik_contrib_fn(tm)
      if (length(lp) == Tn && length(lm) == Tn &&
          all(is.finite(lp)) && all(is.finite(lm))) {
        S[, j] <- (lp - lm) / (2 * hj)
        break
      }
      hj <- hj / 2
    }
  }
  S
}


# Long-run ("meat") covariance of the per-period scores, on the SUM scale so it
# is directly comparable with the observed information H = -d2 sum(l_t):
#     Omega = T * LRV(s_t),
# with LRV the Bartlett (Newey-West) estimator .newey_west() and the lag chosen
# by Andrews' (1991) AR(1) plug-in bandwidth S_T = 1.1447 (alpha T)^(1/3) --
# the rule D20 uses in .d20_moment_sampling_cov(), but applied to the
# STANDARDISED score columns (Andrews' unit weights would otherwise let the
# largest-scale score choose the lag for all of them).
#
# Why HAC and not the plain OPG crossprod(S): under a correctly specified
# prediction-error likelihood the scores are a martingale difference sequence
# and the two agree, but a misspecified DYNAMIC model leaves serially
# correlated scores, whose variance the plain OPG misses entirely (an
# iid-Gaussian fit to AR(1) data with rho = 0.5 has se_robust/se_nominal =
# sqrt(3) for the mean, which the OPG reports as 1). hac = FALSE gives the
# White (1982) OPG. Returns the k x k matrix with attribute "bandwidth"
# (number of Bartlett lags; 0 = none).
.d35_score_lrv <- function(S, hac = TRUE) {
  S  <- scale(as.matrix(S), scale = FALSE)
  Tn <- nrow(S)
  if (!isTRUE(hac) || Tn < 3L) {
    Om <- crossprod(S)
    attr(Om, "bandwidth") <- 0L
    return(Om)
  }
  num <- 0; den <- 0
  for (a in seq_len(ncol(S))) {
    v <- mean(S[, a]^2)
    if (!(v > 0)) next
    x1 <- S[-1, a]; x0 <- S[-Tn, a]
    r  <- min(max(sum(x1 * x0) / sum(x0^2), -0.99), 0.99)
    s2 <- mean((x1 - r * x0)^2) / v
    num <- num + 4 * r^2 * s2^2 / ((1 - r)^6 * (1 + r)^2)
    den <- den + s2^2 / (1 - r)^4
  }
  S_T <- if (den > 0) 1.1447 * (num / den * Tn)^(1 / 3) else 1
  L <- as.integer(min(max(ceiling(S_T) - 1, 0), Tn - 1))
  Om <- Tn * .newey_west(S, max_lag = L)
  Om <- (Om + t(Om)) / 2
  attr(Om, "bandwidth") <- L
  Om
}


# ---------------------------------------------------------------------------
#' D35. Misspecification softness (sandwich / IM-equality)
#'
#' Per-parameter ratio of the misspecification-robust (sandwich) standard error
#' to the nominal (inverse observed information) one,
#' \eqn{\rho_j = \sqrt{[H^{-1}\Omega H^{-1}]_{jj} / [H^{-1}]_{jj}}}, where
#' \eqn{H} is the finite-difference observed information of the total
#' log-likelihood and \eqn{\Omega} is the long-run covariance (sum scale) of
#' the finite-difference per-period scores. Under correct specification the
#' information-matrix equality gives \eqn{\rho_j \to 1}.
#'
#' @param theta       Named numeric vector -- the estimate (mode) to assess.
#' @param loglik_contrib_fn Function: \code{theta -> } numeric vector of the
#'   \eqn{T} per-period log-likelihood contributions (e.g. the Kalman
#'   prediction-error decomposition). Likelihood only (no prior term).
#' @param param_names Optional names (defaults to \code{names(theta)}).
#' @param prior_info  Optional \eqn{k \times k} prior precision matrix added to
#'   the bread \code{H} (use when assessing a posterior mode rather than an MLE).
#' @param softness_threshold \eqn{\rho} above which a parameter is flagged soft
#'   (default 2, i.e. the robust variance exceeds 4x the nominal one).
#' @param hac Logical. \code{TRUE} (default): \eqn{\Omega} is the Bartlett HAC
#'   long-run covariance of the scores (Andrews plug-in bandwidth), which also
#'   detects misspecified dynamics (serially correlated scores).
#'   \code{FALSE}: the White (1982) outer product of the scores.
#' @param eps         Finite-difference step scale (default 1e-4).
#' @param meta        Optional \code{\link{diag_meta}} provenance descriptor.
#' @return A \code{dynhr_diagnostic}; \code{$result$passport_axis} is a named
#'   logical (TRUE = robust/firm) for the Passport.
#' @noRd
d35_misspecification_softness <- function(theta,
                                          loglik_contrib_fn,
                                          param_names = NULL,
                                          prior_info  = NULL,
                                          softness_threshold = 2,
                                          hac  = TRUE,
                                          eps  = 1e-4,
                                          meta = NULL) {
  if (!is.function(loglik_contrib_fn))
    .dynhr_abort("D35: `loglik_contrib_fn` must be a function.")
  if (!is.numeric(theta) || length(theta) < 1L || !all(is.finite(theta)))
    .dynhr_abort("D35: `theta` must be a non-empty finite numeric vector.")
  if (!is.numeric(softness_threshold) || length(softness_threshold) != 1L ||
      !is.finite(softness_threshold) || softness_threshold < 1)
    .dynhr_abort("D35: `softness_threshold` must be a single number >= 1.")
  if (is.null(param_names))
    param_names <- names(theta) %||% paste0("theta_", seq_along(theta))
  if (length(param_names) != length(theta))
    .dynhr_abort("D35: `param_names` must have one entry per element of `theta`.")
  # Keep theta NAMED: a real loglik_contrib_fn keys off names(theta) to know
  # which parameters to set, so stripping names (as.numeric) silently freezes
  # it at the baseline and yields zero scores.
  theta <- stats::setNames(as.numeric(theta), param_names)
  k <- length(theta)
  if (!is.null(prior_info)) {
    prior_info <- as.matrix(prior_info)
    if (!identical(dim(prior_info), c(k, k)) || !all(is.finite(prior_info)))
      .dynhr_abort("D35: `prior_info` must be a finite ", k, " x ", k, " matrix.")
  }

  ll0 <- loglik_contrib_fn(theta)
  if (!is.numeric(ll0) || length(ll0) < 3L || !all(is.finite(ll0)))
    return(.make_result(pass = NA,
      summary = "D35 softness: loglik_contrib_fn did not return a finite per-period vector (length >= 3)."))

  # per-period scores
  S <- .d35_score_matrix(loglik_contrib_fn, theta, eps = eps)
  ok_cols <- which(apply(S, 2, function(c) all(is.finite(c))))
  dropped <- param_names[setdiff(seq_len(k), ok_cols)]
  if (length(ok_cols) < 1)
    return(.make_result(pass = NA,
      summary = "D35 softness: per-period scores could not be computed."))
  S  <- S[, ok_cols, drop = FALSE]
  pn <- param_names[ok_cols]
  th <- theta[ok_cols]

  # meat: long-run covariance of the scores (sum scale)
  Omega <- .d35_score_lrv(S, hac = hac)
  bandwidth <- attr(Omega, "bandwidth")
  attr(Omega, "bandwidth") <- NULL
  dimnames(Omega) <- list(pn, pn)

  # bread: observed information = -Hessian of the total loglik
  total_ll <- function(p) {
    full <- theta; full[ok_cols] <- p
    sum(loglik_contrib_fn(full))
  }
  H <- .deep_observed_information(total_ll, th, eps = max(eps, 1e-4))
  if (is.null(H))
    return(.make_result(pass = NA,
      summary = "D35 softness: observed-information Hessian could not be computed."))
  if (!is.null(prior_info))
    H <- H + prior_info[ok_cols, ok_cols, drop = FALSE]
  H <- (H + t(H)) / 2
  dimnames(H) <- list(pn, pn)

  # (pseudo-)inverse of the bread via its eigen-decomposition; the bread must
  # be positive definite for the sandwich to mean anything.
  eg  <- eigen(H, symmetric = TRUE)
  tol <- length(pn) * max(abs(eg$values)) * .Machine$double.eps
  nonpd <- !all(eg$values > tol)
  d_inv <- ifelse(abs(eg$values) > tol, 1 / eg$values, 0)
  Hinv  <- eg$vectors %*% (d_inv * t(eg$vectors))

  V_sand <- Hinv %*% Omega %*% Hinv
  dn <- diag(Hinv); ds <- diag(V_sand)
  se_nom <- ifelse(dn > 0, sqrt(pmax(dn, 0)), NA_real_)
  se_rob <- ifelse(ds >= 0, sqrt(pmax(ds, 0)), NA_real_)
  softness <- se_rob / se_nom
  if (nonpd) softness[] <- NA_real_         # sandwich undefined off a clean optimum
  names(softness) <- pn

  # global information-matrix-equality discrepancy
  im_trace_ratio <- sum(diag(Omega %*% Hinv)) / length(pn)   # ~1 if correct
  fro <- function(M) sqrt(sum(M^2))
  im_frob_gap <- fro(Omega - H) / max(fro(H), .Machine$double.eps)

  is_soft <- is.finite(softness) & softness > softness_threshold
  soft_params <- pn[is_soft]
  passport_axis <- stats::setNames(
    is.finite(softness) & softness <= softness_threshold, pn)

  pass <- if (nonpd || length(soft_params) > 0) FALSE
          else if (!all(is.finite(softness))) NA
          else TRUE

  tab <- data.frame(
    param      = pn,
    se_nominal = unname(se_nom),
    se_robust  = unname(se_rob),
    softness   = unname(softness),
    verdict    = ifelse(!is.finite(softness), "n/a",
                        ifelse(is_soft, "soft", "firm")),
    stringsAsFactors = FALSE)
  tab <- tab[order(-tab$softness, na.last = TRUE), , drop = FALSE]
  rownames(tab) <- NULL

  meat_lbl <- if (!isTRUE(hac)) "OPG meat"
              else sprintf("HAC meat, %d Bartlett lag%s", bandwidth,
                           if (bandwidth == 1L) "" else "s")

  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE))
    plots$softness <- .plot_d35_softness(tab, softness_threshold, meat_lbl,
                                         nonpd, meta)

  drop_txt <- if (length(dropped))
    sprintf(" Not assessed (non-finite scores): %s.", paste(dropped, collapse = ", "))
  else ""
  verdict_txt <-
    if (nonpd) "Information matrix not positive definite (not at a clean optimum); softness not assessed."
    else if (length(soft_params) == 0) sprintf("All estimates firm (rho<=%.2f).", softness_threshold)
    else sprintf("SOFT (rho>%.2f): %s (precision is a specification artefact).",
                 softness_threshold, paste(soft_params, collapse = ", "))
  summary_txt <- sprintf(
    "D35 Misspecification softness: %d params (%s), IM-ratio tr(Omega H^-1)/k=%.2f (1=correct). %s%s",
    length(pn), meat_lbl, im_trace_ratio, verdict_txt, drop_txt)

  badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
  worst <- utils::head(tab, 4)
  action <-
    if (isTRUE(pass))
      "Robust and nominal SEs agree; estimates survive small misspecification (IM equality ~holds)."
    else if (nonpd) "Re-locate the mode before reading sandwich SEs."
    else if (is.na(pass)) "Could not assess every parameter (see summary)."
    else sprintf("%s: robust SE >> nominal SE -- the data's score variability says these are far less pinned down than the model claims. Treat with a sandwich/robust band; suspect a misspecified equation or dynamics.",
                 paste(utils::head(soft_params, 3), collapse = ", "))
  llm <- paste(c(
    sprintf("D35 | Misspecification Softness (sandwich) | %s", badge),
    sprintf("  params=%d soft=%d threshold=rho>%.2f meat=%s im_ratio=%.2f im_frob_gap=%.2f",
            length(pn), length(soft_params), softness_threshold, meat_lbl,
            im_trace_ratio, im_frob_gap),
    sprintf("  softness(rho=se_robust/se_nominal): %s",
            paste(sprintf("%s=%.2f", worst$param, worst$softness), collapse = ", ")),
    if (length(soft_params))
      sprintf("  soft: %s", paste(soft_params, collapse = ", ")),
    if (length(dropped))
      sprintf("  not assessed: %s", paste(dropped, collapse = ", ")),
    if (nonpd) "  warning: H not positive definite (re-check the mode)." else NULL,
    sprintf("  action: %s", action)
  ), collapse = "\n")

  .make_result(
    result = list(table = tab, softness = softness,
                  Omega = Omega, H = H, bandwidth = bandwidth,
                  im_trace_ratio = im_trace_ratio, im_frob_gap = im_frob_gap,
                  soft_params = soft_params, dropped_params = dropped,
                  nonpd = nonpd, passport_axis = passport_axis),
    pass    = pass,
    plots   = plots,
    summary = summary_txt,
    llm_summary = llm)
}


# Bars of softness rho with the IM-equality line and the soft threshold.
.plot_d35_softness <- function(tab, threshold, meat_lbl, nonpd, meta) {
  df <- tab
  df$param <- factor(df$param, levels = rev(df$param))
  df$Verdict <- factor(df$verdict, levels = c("firm", "soft", "n/a"))
  fin <- is.finite(df$softness)
  df$softness_plot <- ifelse(fin, df$softness, 0)
  df$label <- ifelse(fin, sprintf("%.2f", df$softness), "not assessed")
  xmax <- max(c(df$softness_plot, threshold), na.rm = TRUE) * 1.15

  p <- ggplot2::ggplot(df, ggplot2::aes(x = softness_plot, y = param)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dotted",
                        colour = "grey35", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = dynhr_colours$red, linewidth = 0.5) +
    ggplot2::geom_col(ggplot2::aes(fill = Verdict), width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.15, size = 3.2,
                       colour = "grey20") +
    ggplot2::scale_fill_manual(
      values = c(firm = dynhr_colours$teal, soft = dynhr_colours$red,
                 `n/a` = dynhr_na_fill),
      drop = TRUE, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, xmax),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title    = "D35: Misspecification softness (sandwich vs Hessian SE)",
      subtitle = sprintf(paste0("rho = robust / nominal SE, %s\n",
                                "dotted: rho = 1 (IM equality holds);  dashed: rho = %.2f (soft threshold)"),
                         meat_lbl, threshold),
      caption = if (nonpd)
        "Observed information is not positive definite (not at a clean optimum): softness not assessed."
      else NULL,
      x = "rho = robust (sandwich) SE / nominal (inverse-Hessian) SE", y = NULL)
  p <- p + theme_dynhr()
  if (!any(fin)) p <- p + ggplot2::theme(legend.position = "none")
  .apply_meta(p, meta)
}
