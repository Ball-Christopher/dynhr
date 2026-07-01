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
## equals the outer product of the per-period scores Omega = sum_t s_t s_t'
## (the "meat"); White (1982). Their divergence is misspecification, and the
## misspecification-robust ("sandwich") covariance is
##     V_sand = H^{-1} Omega H^{-1}        (vs the nominal V = H^{-1}).
## The per-parameter SOFTNESS is the ratio of the robust to the nominal
## standard error,
##     rho_j = se_robust_j / se_nominal_j ,
## which is exactly the Andrews-Mikusheva (2014) "discrepancy between two
## estimates of Fisher information" read parameter-by-parameter, and the Mueller
## (2012) robustness spirit. rho_j ~ 1 = firm; rho_j >> 1 = soft.
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
##   Andrews, I., & Mikusheva, A. (2014). Weak identification in maximum
##     likelihood: a question of information. AER P&P, 104(5), 195-199.
##   Mueller, U. K. (2012). Measuring prior sensitivity and prior
##     informativeness in large Bayesian models. JME, 59(6), 581-597.
##   Bonhomme, S., & Weidner, M. (2022). Minimizing sensitivity to model
##     misspecification. Quantitative Economics, 13(3), 907-954.
## --------------------------------------------------------------------------


# Per-period score matrix S (T x k) by central finite differences of the
# log-likelihood contributions. Robust to occasional non-finite evaluations
# (halves the step once before giving up on a column).
.d35_score_matrix <- function(loglik_contrib_fn, theta, eps = 1e-4) {
  k  <- length(theta)
  h0 <- pmax(abs(theta), 1e-3) * eps
  ll0 <- loglik_contrib_fn(theta)
  Tn <- length(ll0)
  S  <- matrix(NA_real_, Tn, k)
  for (j in seq_len(k)) {
    hj <- h0[j]
    col <- NULL
    for (attempt in 1:2) {
      tp <- theta; tp[j] <- tp[j] + hj
      tm <- theta; tm[j] <- tm[j] - hj
      lp <- tryCatch(loglik_contrib_fn(tp), error = function(e) NULL)
      lm <- tryCatch(loglik_contrib_fn(tm), error = function(e) NULL)
      if (!is.null(lp) && !is.null(lm) &&
          length(lp) == Tn && length(lm) == Tn &&
          all(is.finite(lp)) && all(is.finite(lm))) {
        col <- (lp - lm) / (2 * hj); break
      }
      hj <- hj / 2
    }
    if (!is.null(col)) S[, j] <- col
  }
  S
}


# ---------------------------------------------------------------------------
#' D35. Misspecification softness (sandwich / IM-equality)
#'
#' @param theta       Named numeric vector -- the estimate (mode) to assess.
#' @param loglik_contrib_fn Function: \code{theta -> } numeric vector of the
#'   \eqn{T} per-period log-likelihood contributions (e.g. the Kalman
#'   prediction-error decomposition).
#' @param param_names Optional names (defaults to \code{names(theta)}).
#' @param prior_info  Optional \eqn{k \times k} prior precision matrix added to
#'   the bread \code{H} (use when assessing a posterior mode rather than an MLE).
#' @param softness_threshold \eqn{\rho} above which a parameter is flagged soft
#'   (default 2).
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
                                          eps  = 1e-4,
                                          meta = NULL) {
  tryCatch({
    if (is.null(param_names))
      param_names <- names(theta) %||% paste0("theta_", seq_along(theta))
    # Keep theta NAMED: a real loglik_contrib_fn keys off names(theta) to know
    # which parameters to set, so stripping names (as.numeric) silently freezes
    # it at the baseline and yields zero scores.
    theta <- stats::setNames(as.numeric(theta), param_names)
    k <- length(theta)

    ll0 <- tryCatch(loglik_contrib_fn(theta), error = function(e) NULL)
    if (is.null(ll0) || !all(is.finite(ll0)) || length(ll0) < 2)
      return(.make_result(pass = NA,
        summary = "D35 softness: loglik_contrib_fn did not return a finite per-period vector."))

    # meat: outer product of per-period scores
    S <- .d35_score_matrix(loglik_contrib_fn, theta, eps = eps)
    ok_cols <- which(apply(S, 2, function(c) all(is.finite(c))))
    if (length(ok_cols) < 1)
      return(.make_result(pass = NA,
        summary = "D35 softness: per-period scores could not be computed."))
    S   <- S[, ok_cols, drop = FALSE]
    pn  <- param_names[ok_cols]
    th  <- theta[ok_cols]
    Omega <- crossprod(S)                      # k x k

    # bread: observed information = -Hessian of the total loglik
    total_ll <- function(p) {
      full <- theta; full[ok_cols] <- p
      sum(loglik_contrib_fn(full))
    }
    H <- .deep_observed_information(total_ll, th, eps = max(eps, 1e-4))
    if (is.null(H))
      return(.make_result(pass = NA,
        summary = "D35 softness: observed-information Hessian could not be computed."))
    if (!is.null(prior_info)) {
      pi_sub <- prior_info[ok_cols, ok_cols, drop = FALSE]
      H <- H + pi_sub
    }
    H <- (H + t(H)) / 2

    # robust (pseudo-)inverse of the bread
    Hinv <- tryCatch(solve(H), error = function(e) NULL)
    if (is.null(Hinv)) {
      sv <- svd(H); tol <- max(dim(H)) * max(sv$d) * .Machine$double.eps
      d_inv <- ifelse(sv$d > tol, 1 / sv$d, 0)
      Hinv <- sv$v %*% (d_inv * t(sv$u))
    }

    V_nom  <- Hinv
    V_sand <- Hinv %*% Omega %*% Hinv

    dn <- pmax(diag(V_nom),  0)
    ds <- pmax(diag(V_sand), 0)
    se_nom  <- sqrt(dn)
    se_rob  <- sqrt(ds)
    softness <- ifelse(se_nom > 0, se_rob / se_nom, NA_real_)
    names(softness) <- pn

    nonpd <- any(diag(V_nom) <= 0)             # bread not positive definite

    # global information-matrix-equality discrepancy
    im_trace_ratio <- sum(diag(Omega %*% Hinv)) / length(ok_cols)   # ~1 if correct
    fro <- function(M) sqrt(sum(M^2))
    im_frob_gap <- fro(Omega - H) / max(fro(H), .Machine$double.eps)

    soft_params <- pn[is.finite(softness) & softness > softness_threshold]
    passport_axis <- stats::setNames(
      is.finite(softness) & softness <= softness_threshold, pn)

    pass <- length(soft_params) == 0 && !nonpd

    tab <- data.frame(
      param    = pn,
      se_nominal = se_nom,
      se_robust  = se_rob,
      softness   = as.numeric(softness),
      verdict    = ifelse(is.finite(softness) & softness > softness_threshold,
                          "soft", "firm"),
      stringsAsFactors = FALSE)
    tab <- tab[order(-tab$softness), , drop = FALSE]

    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE))
      plots$softness <- .plot_d35_softness(tab, softness_threshold, meta)

    summary_txt <- sprintf(
      "D35 Misspecification softness: %d params, IM-ratio tr(OmegaH^-1)/k=%.2f (1=correct). %s%s",
      length(pn), im_trace_ratio,
      if (length(soft_params) == 0) "All estimates robust (rho<=threshold)."
      else sprintf("SOFT: %s (precision is a specification artefact).",
                   paste(soft_params, collapse = ", ")),
      if (nonpd) " WARNING: information matrix not positive definite (not at a clean optimum)." else "")

    badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
    llm <- paste(c(
      sprintf("D35 | Misspecification Softness (sandwich) | %s", badge),
      sprintf("  params=%d soft=%d threshold=rho>%.1f im_ratio=%.2f im_frob_gap=%.2f",
              length(pn), length(soft_params), softness_threshold,
              im_trace_ratio, im_frob_gap),
      {
        worst <- utils::head(tab, 4)
        sprintf("  softness(rho=se_robust/se_nominal): %s",
                paste(sprintf("%s=%.2f", worst$param, worst$softness), collapse = ", "))
      },
      if (length(soft_params))
        sprintf("  soft: %s", paste(soft_params, collapse = ", ")),
      if (nonpd) "  warning: H not positive definite (re-check the mode)." else NULL,
      sprintf("  action: %s",
              if (isTRUE(pass))
                "Robust and nominal SEs agree; estimates survive small misspecification (IM equality ~holds)."
              else if (is.na(pass)) "Could not assess (see summary)."
              else sprintf("%s: robust SE >> nominal SE -- the data's score variability says these are far less pinned down than the model claims. Treat with a sandwich/robust band; suspect a misspecified equation.",
                           paste(utils::head(soft_params, 3), collapse = ", ")))
    ), collapse = "\n")

    .make_result(
      result = list(table = tab, softness = softness,
                    Omega = Omega, H = H,
                    im_trace_ratio = im_trace_ratio, im_frob_gap = im_frob_gap,
                    soft_params = soft_params, passport_axis = passport_axis),
      pass    = pass,
      plots   = plots,
      summary = summary_txt,
      llm_summary = llm)
  }, error = function(e) {
    .make_result(pass = NA, errored = TRUE,
                 summary = paste("D35 misspecification softness: ERROR --",
                                 conditionMessage(e)))
  })
}


# Bars of softness rho with nominal/robust SE context and the threshold.
.plot_d35_softness <- function(tab, threshold, meta) {
  df <- tab
  df$param <- factor(df$param, levels = rev(df$param))
  df$Verdict <- factor(df$verdict, levels = c("firm", "soft"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = softness, y = param, fill = Verdict)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dotted",
                        colour = dynhr_colours$grey, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = dynhr_colours$red, linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c(firm = dynhr_colours$teal, soft = dynhr_colours$red),
      drop = FALSE, name = NULL) +
    ggplot2::labs(
      title    = "D35: Misspecification softness (sandwich / IM-equality)",
      subtitle = sprintf("rho = robust SE / nominal SE; rho>%.1f (dashed) = soft; rho=1 (dotted) = IM equality",
                         threshold),
      x = expression(rho == "se"[robust] / "se"[nominal]), y = NULL)
  p <- p + theme_dynhr()
  .apply_meta(p, meta)
}
