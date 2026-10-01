## R/diag-mcmc-d6-prior-posterior.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R; ridge-density rewrite Phase-3+.
##
## D6 posterior vs prior -- small-multiple density overlay
##
## 0.9.4 refresh (adversarial review):
##  * The prior curve used to be RENORMALISED over a grid spanning only the
##    posterior's range (+/-15%). That turns a diffuse prior into a uniform
##    on the posterior's support and inflates the overlap coefficient:
##    N(0,100) prior vs N(0,1) posterior scored 0.40 (true OVL 0.027); a
##    posterior sitting in the prior's tail (N(3,0.3) vs N(0,1)) scored 0.094
##    (true 0.018). The prior is now used as the proper density it is.
##  * The posterior auto-dispatch built its density from log_prior_density(),
##    not from log_prior() (the density the samplers actually use), fed NA
##    bounds into a scalar `if` (error -> NA prior -> overlap 0 -> silent
##    PASS), silently substituted N(0,1) for an unmatched parameter name, and
##    never normalised a truncated prior. See .d6_prior_density_fn().
##  * Draw columns are matched to param_names BY NAME (they used to be
##    labelled positionally).
##  * A parameter with no usable prior density gets overlap NA and is
##    reported, instead of scoring 0 and counting towards PASS.
## --------------------------------------------------------------------------

## TRUE when (dist, p1..p4) name a proper law that .prior_law() can build.
## .prior_law() ABORTS on an improper / impossible shape (it is the sampler's
## guard); D6 wants NA for such a row instead, so the same conditions are
## tested here up front rather than catching the abort.
.diag_prior_law_proper <- function(d, p1, p2, p3 = NA_real_, p4 = NA_real_) {
  ## Uniform (Dynare convention): support [p3, p4] if both given -- p1/p2 may
  ## then be NA -- else p1 -+ sqrt(3) p2 (mean/sd).
  if (identical(d, "uniform")) {
    ab <- .uniform_ab(p1, p2, p3, p4)
    return(is.numeric(ab) && all(is.finite(ab)) && ab[1] < ab[2])
  }
  if (!is.numeric(p1) || !is.numeric(p2) || is.na(p1) || is.na(p2))
    return(FALSE)
  s <- .prior_gen_ab(p3, NA_real_)[1]
  switch(d,
    "normal"  = is.finite(p1) && is.finite(p2) && p2 > 0,
    "beta"    = {
      ab <- .prior_gen_ab(p3, p4)
      sh <- .beta_shapes(p1, p2, ab[1], ab[2])
      ab[2] > ab[1] && all(is.finite(sh)) && all(sh > 0)
    },
    "gamma"   = is.finite(p1 - s) && p1 - s > 0 && is.finite(p2) && p2 > 0,
    ## sd = Inf is a PROPER inverse gamma in both conventions (IG1
    ## alpha -> 1, IG2 shape = 2), so only p2 > 0 is required.
    "inv_gamma" =, "inv_gamma1" =, "inv_gamma2" =
      is.finite(p1 - s) && p1 - s > 0 && p2 > 0,
    FALSE)
}

## CDF of one prior, as the law log_prior() scores (R/prior-density.R
## .prior_law(): generalised beta on [p3, p4], shift p3 for gamma / inverse
## gammas), used only to get the truncation normaliser F(upper) - F(lower).
## Returns NA when the prior is improper or the distribution is unknown.
.d6_prior_cdf <- function(q, dist, p1, p2, p3 = NA_real_, p4 = NA_real_) {
  d <- .normalize_dist(dist)
  if (!.diag_prior_law_proper(d, p1, p2, p3, p4)) return(NA_real_)
  law <- .prior_law(d, p1, p2, p3, p4)
  ## `law$p` is only defined strictly inside the natural support.
  if (q <= law$lo) 0 else if (q >= law$hi) 1 else law$p(q)
}

## Build the D6 prior density from a prior_spec data.frame (name,
## distribution, p1, p2, lower, upper). The density is exp(log_prior()) on a
## one-row spec -- the exact function the posterior closures evaluate, bounds
## included -- divided by the prior mass inside [lower, upper] so a truncated
## prior is a proper density. Returns NA (not a stand-in density) for a name
## not in the spec or an improper prior.
.d6_prior_density_fn <- function(priors) {
  priors <- as.data.frame(priors, stringsAsFactors = FALSE)
  norm_const <- vapply(seq_len(nrow(priors)), function(i) {
    lo <- suppressWarnings(as.numeric(priors$lower[i]))
    hi <- suppressWarnings(as.numeric(priors$upper[i]))
    if (is.na(lo)) lo <- -Inf
    if (is.na(hi)) hi <-  Inf
    p1 <- as.numeric(priors$p1[i]); p2 <- as.numeric(priors$p2[i])
    ## Dynare p3/p4: generalised beta support / shift.
    p3 <- if (is.null(priors$p3)) NA_real_ else as.numeric(priors$p3[i])
    p4 <- if (is.null(priors$p4)) NA_real_ else as.numeric(priors$p4[i])
    Fh <- .d6_prior_cdf(hi, priors$distribution[i], p1, p2, p3, p4)
    Fl <- .d6_prior_cdf(lo, priors$distribution[i], p1, p2, p3, p4)
    z <- Fh - Fl
    if (is.finite(z) && z > 0) z else NA_real_
  }, numeric(1))

  function(x, param_name) {
    idx <- match(param_name, priors$name)
    if (is.na(idx) || is.na(norm_const[idx])) return(rep(NA_real_, length(x)))
    spec <- priors[idx, , drop = FALSE]
    vapply(x, function(xx) {
      lp <- log_prior(stats::setNames(xx, param_name), spec)
      exp(lp) / norm_const[idx]
    }, numeric(1))
  }
}

#' D6. Posterior vs prior updating
#'
#' Overlays each parameter's prior density and posterior kernel density in a
#' small-multiple grid (one panel per parameter, free axes) and scores the
#' overlap coefficient \eqn{\mathrm{OVL} = \int \min(p(\theta \mid y),
#' p(\theta))\,d\theta}.  OVL is 1 when the data leave the prior unchanged
#' and near 0 when they move or sharpen it substantially.  Parameters with
#' OVL > 0.80 are flagged (weakly updated) and sorted first.
#'
#' @param draws           Matrix (n_draws x n_params) -- posterior draws.
#'   When it has column names they are matched to \code{param_names} by name.
#' @param prior_density_fn Function: (x, param_name) -> prior density value(s).
#'   Must be the PROPER (normalised, truncation-adjusted) prior density used in
#'   estimation; it is not renormalised here. Return \code{NA} for a parameter
#'   with no known prior (its OVL is then \code{NA} and it is excluded from the
#'   badge). Called elementwise, so it need not be vectorised.
#' @param param_names     Character vector (optional; taken from \code{colnames(draws)})
#' @param meta            A \code{\link{diag_meta}} object for plot provenance
#'   (model name, data hash, date caption).  Pass \code{NULL} to suppress.
#' @param threshold       OVL above which a parameter is flagged (default 0.80).
#' @return dynhr_diagnostic list; \code{result$overlap_scores} (named, may
#'   contain NA), \code{result$uninformative}, \code{result$no_prior}.
#' @references Inman, H. F. and Bradley, E. L. (1989). The overlapping
#'   coefficient as a measure of agreement between probability distributions.
#'   \emph{Communications in Statistics - Theory and Methods} 18(10), 3851-3874.
#'
#' @noRd
d6_posterior_vs_prior <- function(draws,
                                  prior_density_fn,
                                  param_names = NULL,
                                  meta        = NULL,
                                  threshold   = 0.80) {

    draws <- as.matrix(draws)
    if (is.null(param_names)) {
      param_names <- if (!is.null(colnames(draws))) colnames(draws)
      else paste0("theta_", seq_len(ncol(draws)))
    }
    if (!is.null(colnames(draws))) {
      missing_cols <- setdiff(param_names, colnames(draws))
      if (length(missing_cols) > 0)
        .dynhr_abort("D6: draws have no column(s) for: ",
                     paste(missing_cols, collapse = ", "))
      draws <- draws[, param_names, drop = FALSE]
    } else if (ncol(draws) != length(param_names)) {
      .dynhr_abort(sprintf("D6: draws have %d columns but %d param_names.",
                           ncol(draws), length(param_names)))
    }
    n_par <- length(param_names)

    # ---- Per-parameter KDE and overlap score --------------------------------
    overlap_scores <- setNames(rep(NA_real_, n_par), param_names)
    density_rows   <- vector("list", n_par)

    for (j in seq_len(n_par)) {
      xj <- draws[is.finite(draws[, j]), j]
      if (length(xj) < 2L || stats::sd(xj) == 0) next
      post_kde   <- stats::density(xj, n = 512)
      dx         <- diff(post_kde$x[1:2])

      # Modestly extend the display grid (+/-15% of KDE range) so the prior is
      # not hard-clipped to the posterior support.
      kde_range  <- diff(range(post_kde$x))
      prior_x    <- seq(min(post_kde$x) - 0.15 * kde_range,
                        max(post_kde$x) + 0.15 * kde_range, length.out = 512L)
      .pf <- function(xs) vapply(xs, function(xx) {
        v <- prior_density_fn(xx, param_names[j])
        if (length(v) != 1L || is.na(v) || is.nan(v) || v < 0) NA_real_
        else if (!is.finite(v)) NA_real_ else v
      }, numeric(1))
      prior_wide <- .pf(prior_x)
      prior_kde  <- .pf(post_kde$x)

      has_prior <- any(is.finite(prior_kde)) && any(prior_kde > 0, na.rm = TRUE)
      if (has_prior) {
        pk <- prior_kde; pk[!is.finite(pk)] <- 0
        # Overlap coefficient: integral of min(posterior, prior).  min() is
        # bounded by the posterior, so the posterior KDE grid (which carries
        # essentially all posterior mass) is the integration domain.
        overlap_scores[j] <- sum(pmin(post_kde$y, pk)) * dx
      }

      density_rows[[j]] <- data.frame(
        x            = c(post_kde$x, if (has_prior) prior_x),
        density      = c(post_kde$y, if (has_prior) prior_wide),
        Distribution = rep(c("Posterior", "Prior"),
                           c(length(post_kde$x), if (has_prior) length(prior_x) else 0L)),
        Parameter    = param_names[j],
        stringsAsFactors = FALSE
      )
    }

    uninformative <- param_names[!is.na(overlap_scores) & overlap_scores > threshold]
    no_prior      <- param_names[is.na(overlap_scores)]
    scored        <- overlap_scores[!is.na(overlap_scores)]

    # ---- Small-multiple layout ----------------------------------------------
    plots <- list()
    all_df <- do.call(rbind, density_rows)

    if (!is.null(all_df) && requireNamespace("ggplot2", quietly = TRUE)) {
      # Flagged (weakly updated) first; within a group, highest OVL first.
      ord <- order(!(param_names %in% uninformative),
                   -ifelse(is.na(overlap_scores), -1, overlap_scores))
      ordered_params <- param_names[ord]
      strip_lab <- setNames(
        sprintf("%s%s  (OVL %s)",
                ifelse(param_names %in% uninformative, "* ", ""),
                param_names,
                ifelse(is.na(overlap_scores), "n/a: no prior",
                       sprintf("%.2f", overlap_scores))),
        param_names)
      all_df <- all_df[all_df$Parameter %in% ordered_params, ]
      all_df$Panel <- factor(strip_lab[all_df$Parameter],
                             levels = unname(strip_lab[ordered_params]))
      all_df$Distribution <- factor(all_df$Distribution,
                                    levels = c("Prior", "Posterior"))

      flag_label <- if (length(uninformative) > 0)
        sprintf("* weakly updated (OVL > %.2f): %s", threshold,
                paste(uninformative, collapse = ", "))
      else if (length(scored) == 0) "No parameter could be scored"
      else sprintf("No parameter has OVL > %.2f: all scored posteriors moved away from the prior",
                   threshold)
      if (length(no_prior) > 0)
        flag_label <- paste0(flag_label, "\nNo usable prior density: ",
                             paste(no_prior, collapse = ", "))

      cols <- c("Prior" = dynhr_colours$orange, "Posterior" = dynhr_colours$dark_blue)
      p <- ggplot2::ggplot(all_df, ggplot2::aes(x = x, group = Distribution)) +
        ggplot2::geom_area(ggplot2::aes(y = density, fill = Distribution),
                           position = "identity", alpha = 0.35, colour = NA) +
        ggplot2::geom_line(ggplot2::aes(y = density, colour = Distribution),
                           linewidth = 0.55) +
        ## 4 columns up to 24 parameters, 6 beyond: keeps a 68-parameter grid
        ## within one page height when the report caps the figure.
        ggplot2::facet_wrap(~ Panel, ncol = min(if (n_par > 24L) 6L else 4L, n_par),
                            scales = "free") +
        ggplot2::scale_x_continuous(
          n.breaks = 3L,
          labels   = function(x) formatC(x, format = "g", digits = 3)
        ) +
        ggplot2::scale_fill_manual(values = cols, name = NULL) +
        ggplot2::scale_colour_manual(values = cols, name = NULL) +
        theme_dynhr_compact() +
        ggplot2::theme(panel.spacing.x = ggplot2::unit(1.4, "lines")) +
        ggplot2::labs(
          title    = "D6: Prior vs posterior marginal densities",
          subtitle = flag_label,
          x        = "Parameter value (estimation scale)",
          y        = "Density",
          caption  = "OVL = overlap coefficient, integral of min(prior, posterior); 1 = data did not update the prior"
        )

      p <- .apply_meta(p, meta)
      ncol_pp <- min(if (n_par > 24L) 6L else 4L, n_par)
      attr(p, "dynhr_fig_height") <- min(30, max(5, ceiling(n_par / ncol_pp) * 1.9 + 1.2))
      plots$prior_posterior <- p
    }

    pass <- if (length(scored) == 0) NA else length(uninformative) == 0
    n_uninformative <- length(uninformative)
    rng <- if (length(scored)) sprintf("[%.2f, %.2f]", min(scored), max(scored)) else "n/a"
    badge <- if (is.na(pass)) "N/A" else if (pass) "PASS" else "FAIL"

    .make_result(
      result  = list(overlap_scores = overlap_scores,
                     uninformative  = uninformative,
                     no_prior       = no_prior,
                     threshold      = threshold),
      pass    = pass,
      plots   = plots,
      summary = paste0(
        sprintf("D6 Posterior vs prior: %d params. Overlap range %s. ", n_par, rng),
        if (is.na(pass)) "N/A -- no parameter has a usable prior density."
        else if (pass) sprintf("PASS -- no parameter has OVL > %.2f.", threshold)
        else sprintf("FAIL -- %d weakly updated (OVL > %.2f): %s",
                     n_uninformative, threshold, paste(uninformative, collapse = ", ")),
        if (length(no_prior) > 0)
          sprintf(" No usable prior for: %s.", paste(no_prior, collapse = ", "))
      ),
      llm_summary = {
        worst5_uninf <- head(sort(scored, decreasing = TRUE), min(5, n_uninformative))
        most_inf     <- head(sort(scored), 5)
        paste(c(
          sprintf("D6 | Prior vs Posterior | %s", badge),
          sprintf("  params=%d scored=%d uninformative_gt%.2f=%d/%d",
                  n_par, length(scored), threshold, n_uninformative, length(scored)),
          if (length(scored) > 0)
            sprintf("  overlap_min=%.2f (%s) overlap_max=%.2f (%s) overlap_median=%.2f",
                    min(scored), names(which.min(scored)),
                    max(scored), names(which.max(scored)),
                    stats::median(scored)),
          if (length(most_inf) > 0)
            sprintf("  most_informative: %s",
                    paste(sprintf("%s=%.2f", names(most_inf), most_inf),
                          collapse = ", ")),
          if (n_uninformative > 0)
            sprintf("  uninformative: %s",
                    paste(sprintf("%s=%.2f", names(worst5_uninf), worst5_uninf),
                          collapse = ", ")),
          if (length(no_prior) > 0)
            sprintf("  no_prior_density: %s", paste(no_prior, collapse = ", ")),
          sprintf("  action: %s",
                  if (is.na(pass))
                    "No prior density could be evaluated; check prior_spec names/distributions."
                  else if (n_uninformative == 0)
                    "All scored parameters updated by data."
                  else
                    sprintf("%s: data barely move the prior (overlap>%.2f). Check identification or the prior.",
                            paste(names(worst5_uninf)[seq_len(min(3, n_uninformative))],
                                  collapse = ", "), threshold))
        ), collapse = "\n")
      }
    )
}
