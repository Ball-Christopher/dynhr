## R/diag-pre-prior-sensitivity.R
## --------------------------------------------------------------------------
## Prior Sensitivity Diagnostic — compares mode estimates under informative
## priors vs flat (max-entropy / uniform) priors over the same support.
##
## Parameters whose posterior mode shifts substantially under flat priors
## are "prior-driven" — the data alone does not pin them down.
##
## Usage:
##   Added to run_all_diagnostics() when the model has an estimated_params
##   block and the data is available.
## --------------------------------------------------------------------------
##
## Flat-prior construction (0.9.4 bug fix PS1).
##
## The flat-prior model is built IN MEMORY from the ALREADY-PARSED
## `model$estimated_params` data frame, never by re-reading and re-splitting
## the .mod source text. The old implementation re-read `source_file`, located
## `estimated_params;` by regex, split each line on "," and took the LAST TWO
## comma fields as (lower, upper). That is wrong for every line carrying a
## trailing `%` / `//` comment that contains a comma -- e.g. NZSIM's
##
##   sigma,  GAMMA_PDF,  1,  0.04;   %In original NZSIMV6 coding,[2,1,0.04];
##
## whose "bounds" were read out of the comment, yielding
## `uniform_pdf, , , 1, 1` and an `extract_prior_spec()` abort ("uniform prior
## 'sigma' needs lower < upper"). It was equally wrong for the long Dynare form
## NAME, INIT, LB, UB, SHAPE, P1, P2, for /* */ blocks and for multi-statement
## lines. `parse_estimated_params_block()` (R/parse-blocks.R) already handles
## all of those, and priors enter neither the compiled model nor the steady
## state, so nothing needs re-solving.
## --------------------------------------------------------------------------

#' Quantile of ONE prior, in the SAME parameterisation as `log_prior_density()`
#'
#' `p1`/`p2` are the (mean, sd) hyperparameters a Dynare-style spec carries,
#' and `p3`/`p4` Dynare's generalisation parameters (beta support [p3, p4];
#' shift p3 for gamma and both inverse gammas). The quantile is `$q` of
#' `.prior_law()` (R/prior-density.R) -- the ONE (dist, p1..p4) -> law mapping
#' the density `.lp_dist1()` and the prior sampler use -- so the quantile and
#' the density cannot drift apart (they did: the shift and the generalised
#' beta support were ignored here).
#'
#' IG1 (Dynare.jl `inv_gamma`): X - s = sqrt(Y), Y ~ InvGamma(alpha, theta).
#' IG2 (standard): X - s ~ InvGamma(shape, scale).
#'
#' @param p     Probability in (0, 1).
#' @param dist  Distribution name (canonical or Dynare `*_pdf`).
#' @param p1,p2 Hyperparameters, as in the prior spec.
#' @param p3,p4 Dynare generalisation parameters (NA = default).
#' @return Scalar quantile (`NA_real_` for a degenerate parameterisation).
#' @noRd
.ps_prior_quantile <- function(p, dist, p1, p2, p3 = NA_real_, p4 = NA_real_) {
  d <- .normalize_dist(dist)
  ## Same properness screen as D6's CDF (.prior_law() aborts on an improper
  ## shape; here a degenerate row just yields NA).
  if (!.diag_prior_law_proper(d, p1, p2, p3, p4)) return(NA_real_)
  .prior_law(d, p1, p2, p3, p4)$q(p)
}

#' Finite flat-prior support for every row of an informative prior spec
#'
#' Rule, per parameter:
#'   * both `lower` and `upper` finite -> uniform on exactly [lower, upper]
#'     (rule `"bounds"`);
#'   * otherwise -> the informative prior's `q_lo` / `q_hi` quantiles,
#'     INTERSECTED with whichever side of [lower, upper] is finite
#'     (rule `"quantile"`).
#'
#' @param spec  Informative prior spec from `extract_prior_spec()`.
#' @param q_lo,q_hi Tail probabilities for the quantile rule.
#' @return data.frame(name, lower, upper, rule), one row per spec row.
#' @noRd
.ps_flat_support <- function(spec, q_lo = 0.005, q_hi = 0.995) {
  n  <- nrow(spec)
  lo <- numeric(n); hi <- numeric(n); rule <- character(n)
  for (i in seq_len(n)) {
    b_lo <- suppressWarnings(as.numeric(spec$lower[i]))
    b_hi <- suppressWarnings(as.numeric(spec$upper[i]))
    if (is.na(b_lo)) b_lo <- -Inf
    if (is.na(b_hi)) b_hi <- Inf
    if (is.finite(b_lo) && is.finite(b_hi)) {
      lo[i] <- b_lo; hi[i] <- b_hi; rule[i] <- "bounds"
    } else {
      p3 <- if (is.null(spec$p3)) NA_real_ else as.numeric(spec$p3[i])
      p4 <- if (is.null(spec$p4)) NA_real_ else as.numeric(spec$p4[i])
      ql <- .ps_prior_quantile(q_lo, spec$distribution[i],
                               spec$p1[i], spec$p2[i], p3, p4)
      qh <- .ps_prior_quantile(q_hi, spec$distribution[i],
                               spec$p1[i], spec$p2[i], p3, p4)
      lo[i] <- max(ql, b_lo); hi[i] <- min(qh, b_hi)
      rule[i] <- "quantile"
    }
  }
  data.frame(name = as.character(spec$name), lower = lo, upper = hi,
             rule = rule, stringsAsFactors = FALSE)
}

#' Measurement-error variance for the flat-prior re-run
#'
#' The flat run must be fitted with the SAME measurement-error variance as the
#' informative run; anything else conflates the prior change with an ME change
#' and every parameter looks prior-driven. `fallback` (the orchestrator's
#' `0.07 * var(data)` heuristic) is used only when the informative mode result
#' carries no `me_variance` at all.
#'
#' @param mode_result The informative `dynhr_mode_result`.
#' @param fallback    Per-observable ME variance to use if the mode result has
#'   none.
#' @return Scalar or per-observable numeric vector.
#' @noRd
.ps_me_variance <- function(mode_result, fallback) {
  mv <- if (is.list(mode_result)) mode_result$me_variance else NULL
  if (is.null(mv) || length(mv) == 0L || !any(is.finite(mv))) fallback else mv
}

#' Named mode vector out of whatever shape the caller's mode result carries
#'
#' `run_mode_finding()` / `run_full_estimation()` / `run_posterior_estimation()`
#' all put `theta_mode` at the TOP level of the result AND under `$mode`, and
#' different orchestrator paths hand this diagnostic different ones. Prefer
#' `$mode$theta_mode` (the optimiser's own record) and fall back to the
#' top-level copy.
#'
#' @param mode_result A mode/estimation result list.
#' @return Named numeric vector, or `NULL` when neither shape is usable.
#' @noRd
.ps_theta_mode <- function(mode_result) {
  if (!is.list(mode_result)) return(NULL)
  for (cand in list(mode_result$mode$theta_mode, mode_result$theta_mode)) {
    if (is.numeric(cand) && length(cand) > 0L && !is.null(names(cand)) &&
        all(nzchar(names(cand))) && all(is.finite(cand)))
      return(cand)
  }
  NULL
}

#' Rewrite an `estimated_params` data frame so every prior becomes uniform
#'
#' Row order, `type`, `name`, `name2`, `init`, `lb` and `ub` are preserved, so
#' `extract_prior_spec()` rebuilds the SAME canonical keys from the flat model
#' (`corr a,b`, `skew x`, the `eps_* -> sig_*` rename).
#'
#' @param ep      `model$estimated_params`.
#' @param support data.frame from `.ps_flat_support()`, row-aligned with `ep`.
#' @return The rewritten data frame.
#' @noRd
.ps_flatten_estimated_params <- function(ep, support) {
  ep$prior <- "uniform_pdf"
  ## The Dynare BOUND form `uniform_pdf, , , lower, upper`: P1/P2 (mean/sd)
  ## empty, the support in P3/P4. Convention-proof -- it does not depend on
  ## how P1/P2 of a uniform are read (Dynare: mean/sd), and it never trips
  ## extract_prior_spec()'s bounds-style-row warning.
  ep$p1 <- NA_real_
  ep$p2 <- NA_real_
  ep$p3 <- support$lower
  ep$p4 <- support$upper
  ep
}


#' Prior sensitivity: informative vs flat (uniform) priors
#'
#' Re-runs mode-finding with uniform priors over the same support as the
#' original informative priors, then compares the two mode vectors.
#' Parameters whose relative shift exceeds \code{threshold} are flagged as
#' "prior-driven."
#'
#' The flat-prior model is built in memory from the parsed
#' \code{model$estimated_params} rows -- no .mod file is rewritten and the model is
#' not re-solved (priors affect neither the compiled model nor the steady
#' state). Each parameter's uniform support is its informative prior's
#' \code{[lower, upper]} when both are finite, and otherwise that prior's 0.5\% /
#' 99.5\% quantiles intersected with whichever bound is finite; the rule used is
#' reported in \code{result$support}.
#'
#' The flat-prior search STARTS AT THE INFORMATIVE MODE (passed through
#' \code{run_mode_finding(theta_init = )}), not at the prior means -- which for
#' a uniform prior are the midpoint of its support. The question asked here is
#' local ("does the mode move when the prior is flattened?"), so a far-away
#' start turns it into a global-optimisation problem and reports the
#' optimiser's failure to converge as evidence about the priors. The flat
#' log-posterior at that start is recorded as \code{result$logpost_flat_at_inf};
#' when the flat search does not beat it (\code{result$improved} is
#' \code{FALSE}) the diagnostic returns INCONCLUSIVE (\code{pass = NA},
#' \code{errored = FALSE}) rather than a mode comparison.
#'
#' @param solved       A \code{dynhr_solved} object.
#' @param data         Observation matrix (T x n_obs).
#' @param obs_vars     Observable variable names.
#' @param mode_inf     The \code{dynhr_mode_result} from informative-prior run.
#' @param me_variance  Measurement error variance (same as used for informative).
#' @param threshold    Relative shift threshold for flagging (default 0.50,
#'   i.e. a 50% change in parameter value relative to its informative-mode
#'   magnitude is flagged).
#' @param n_iter       Optimizer iterations for flat-prior mode finding
#'   (default 500).  Fewer iterations are acceptable since we only need
#'   a rough comparison.
#' @param verbose      Print progress messages.
#' @param meta         Optional \code{\link{diag_meta}()} list applied to the plots
#'   (model name, sample size); \code{NULL} derives a default from the model.
#' @return A \code{dynhr_diagnostic} list.
#' @export
diag_prior_sensitivity <- function(solved,
                                    data,
                                    obs_vars,
                                    mode_inf,
                                    me_variance    = 0,
                                    threshold      = 0.50,
                                    n_iter         = 500L,
                                    verbose        = TRUE,
                                    meta           = NULL) {

  .vcat <- function(...) if (verbose) .dynhr_cat(...)

    .vcat("[prior_sensitivity] Building flat-prior model...\n")

    ep <- solved$model$estimated_params
    if (is.null(ep) || nrow(ep) == 0L) {
      return(.make_result(
        pass = NA,
        summary = "prior_sensitivity: no estimated_params block found.",
        errored = TRUE
      ))
    }

    # ---- 1. Flat-prior support, per parameter ----
    ## The informative spec comes from the package's own parser -- the same one
    ## run_mode_finding() uses -- so `corr` / `skew` / `stderr` keys and the
    ## eps_* -> sig_* rename are already canonical here.
    spec_inf <- extract_prior_spec(solved$model, verbose = FALSE)
    support  <- .ps_flat_support(spec_inf)

    bad <- !is.finite(support$lower) | !is.finite(support$upper) |
      support$lower >= support$upper
    if (any(bad)) {
      return(.make_result(
        pass = NA,
        summary = sprintf(
          paste0("prior_sensitivity: cannot build a finite flat support for ",
                 "parameter(s): %s."),
          paste(support$name[bad], collapse = ", ")),
        errored = TRUE
      ))
    }

    # ---- 2. Flat-prior model, IN MEMORY (no file rewrite, no re-solve) ----
    ## Priors touch neither the compiled model nor the steady state, so the
    ## solved object is reused verbatim apart from the estimated_params rows.
    ## run_mode_finding() re-derives its prior spec from `solved$model` on every
    ## call (nothing on `solved` caches one), so this edit is sufficient.
    solved_flat <- solved
    solved_flat$model$estimated_params <-
      .ps_flatten_estimated_params(ep, support)

    spec_flat <- extract_prior_spec(solved_flat$model, verbose = FALSE)
    .vcat(sprintf(
      "[prior_sensitivity] Flat priors: %d uniform (%d from bounds, %d from prior quantiles)\n",
      nrow(spec_flat), sum(support$rule == "bounds"),
      sum(support$rule == "quantile")))

    ## ---- 2b. Start the flat search AT the informative mode ----
    ## The question this diagnostic asks is LOCAL: "does the mode MOVE when the
    ## prior is flattened?" Starting the flat search from `priors$mean` -- for a
    ## uniform, the MIDPOINT of each support -- asks a different and much harder
    ## question: can Nelder-Mead find the global mode of a 68-dimensional
    ## posterior from a point far away? On NZSIM it could not (flat logpost 6451
    ## vs 9060 informative, 30/68 parameters "prior-driven"), and the answer was
    ## reported as evidence about the priors. Starting at the informative mode
    ## makes the comparison the intended one and makes the flat search cheap.
    theta_start <- .ps_theta_mode(mode_inf)
    if (is.null(theta_start)) {
      return(.make_result(
        pass = NA,
        summary = paste0("prior_sensitivity: the informative mode result ",
                         "carries no named `theta_mode` to start the ",
                         "flat-prior search from."),
        errored = TRUE
      ))
    }

    .vcat("[prior_sensitivity] Running mode-finding with flat priors ",
          "(starting at the informative mode)...\n")
    mode_flat <- run_mode_finding(solved_flat, data, obs_vars = obs_vars,
                                   n_iter = n_iter, method = "nmkb",
                                   me_variance = me_variance,
                                   theta_init = theta_start,
                                   verbose = FALSE)
    if (is.null(mode_flat) || is.null(mode_flat$mode) ||
        !is.finite(mode_flat$mode$logpost %||% NA_real_)) {
      return(.make_result(
        pass = NA,
        summary = "prior_sensitivity: flat-prior mode finding returned invalid result.",
        errored = TRUE
      ))
    }
    .vcat(sprintf("[prior_sensitivity] Flat mode logpost = %.2f\n",
                  mode_flat$mode$logpost))

    ## The flat log-posterior AT the informative mode: the flat search's own
    ## starting value, and hence the floor any converged flat run must clear.
    ## run_mode_finding() clips the start inside the flat bounds, so evaluate at
    ## the vector it ACTUALLY started from (`mode_flat$theta_init`) rather than
    ## at the raw informative mode -- otherwise a parameter sitting on a flat
    ## bound would make `improved` compare two different points.
    lp_start <- mode_flat$log_post_fn(mode_flat$theta_init)
    logpost_flat_at_inf <- lp_start$logpost %||% NA_real_
    improved <- isTRUE(is.finite(logpost_flat_at_inf) &&
                         mode_flat$mode$logpost > logpost_flat_at_inf)

    # ---- 3. Compare modes ----
    theta_inf <- mode_inf$mode$theta_mode
    theta_flat <- mode_flat$mode$theta_mode
    all_params <- union(names(theta_inf), names(theta_flat))

    comparisons <- data.frame(
      param = character(),
      informative = numeric(),
      flat = numeric(),
      abs_shift = numeric(),
      rel_shift = numeric(),
      flagged = logical(),
      stringsAsFactors = FALSE
    )

    for (p in all_params) {
      v_inf <- if (p %in% names(theta_inf)) theta_inf[p] else NA_real_
      v_flat <- if (p %in% names(theta_flat)) theta_flat[p] else NA_real_
      abs_shift <- abs(v_inf - v_flat)
      ref_val <- max(abs(v_inf), 1e-10)
      rel_shift <- abs_shift / ref_val
      flagged <- rel_shift > threshold
      comparisons <- rbind(comparisons, data.frame(
        param = p,
        informative = round(v_inf, 6),
        flat = round(v_flat, 6),
        abs_shift = round(abs_shift, 6),
        rel_shift = round(rel_shift, 3),
        flagged = flagged,
        stringsAsFactors = FALSE
      ))
    }

    n_flagged <- sum(comparisons$flagged)
    all_pass <- n_flagged == 0

    # ---- 4. Build result ----
    result <- list(
      comparisons = comparisons,
      ## Which support each flat prior was given, and by which rule -- the
      ## report shows this so a reader can tell a bound-driven flat prior from
      ## a quantile-driven one.
      support = support,
      prior_spec_inf = spec_inf,
      prior_spec_flat = spec_flat,
      ## Recorded so a caller can verify the flat run was fitted with the same
      ## measurement-error variance as the informative one.
      me_variance = me_variance,
      theta_inf = theta_inf,
      theta_flat = theta_flat,
      n_flagged = n_flagged,
      threshold = threshold,
      logpost_inf = mode_inf$mode$logpost,
      logpost_flat = mode_flat$mode$logpost,
      ## The flat search's own starting point (the informative mode, clipped
      ## inside the flat supports) and its flat log-posterior there, so a
      ## reader can see whether the flat optimiser actually got anywhere.
      theta_flat_init = mode_flat$theta_init,
      logpost_flat_at_inf = logpost_flat_at_inf,
      improved = improved
    )

    ## ---- 3b. Inconclusive: the flat search never beat its own start ----
    ## Nothing can be said about prior-drivenness from a search that did not
    ## move: the "flat mode" is just the informative mode, so every shift is
    ## zero for a reason that has nothing to do with the priors. INFO (pass =
    ## NA, errored = FALSE) is this package's inconclusive level -- see
    ## .badge_str() in R/diag-result.R.
    if (!improved) {
      return(.make_result(
        result  = result,
        pass    = NA,
        errored = FALSE,
        summary = sprintf(paste0(
          "Prior sensitivity: INCONCLUSIVE -- the flat-prior search did not ",
          "converge. Started at the informative mode (flat logpost %.4f) and ",
          "ended at %.4f, no improvement, so the mode comparison is not ",
          "informative about the priors. Raise n_iter (currently %d) or use a ",
          "stronger optimiser."),
          logpost_flat_at_inf, mode_flat$mode$logpost, as.integer(n_iter)),
        llm_summary = sprintf(
          "[INFO] PriorSensitivity inconclusive: flat search did not converge (start=%.4f end=%.4f n_iter=%d)",
          logpost_flat_at_inf, mode_flat$mode$logpost, as.integer(n_iter))
      ))
    }

    # Summary
    if (all_pass) {
      summary_str <- sprintf(
        "Prior sensitivity: %d/%d parameters within threshold (%.0f%%). All parameters are data-driven.",
        nrow(comparisons) - n_flagged, nrow(comparisons), threshold * 100)
    } else {
      flagged_params <- comparisons$param[comparisons$flagged]
      summary_str <- sprintf(
        "Prior sensitivity: %d/%d parameters flagged (threshold=%.0f%%). Prior-driven: %s.",
        n_flagged, nrow(comparisons), threshold * 100,
        paste(flagged_params, collapse = ", "))
    }

    # ---- 5. Plot ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      comp_plot <- comparisons
      comp_plot$param <- factor(comp_plot$param, levels = rev(comp_plot$param))
      comp_plot$status <- ifelse(comp_plot$flagged, "Prior-driven", "Data-driven")

      p <- ggplot2::ggplot(comp_plot) +
        ggplot2::geom_segment(
          ggplot2::aes(x = informative, y = param, 
                       xend = flat, yend = param,
                       colour = status),
          linewidth = 1.2, alpha = 0.8) +
        ggplot2::geom_point(
          ggplot2::aes(x = informative, y = param, colour = "Informative"),
          size = 2.5) +
        ggplot2::geom_point(
          ggplot2::aes(x = flat, y = param, colour = "Flat (uniform)"),
          size = 2.5, shape = 17) +
        ggplot2::scale_colour_manual(
          values = c("Informative" = "#0077BB",
                     "Flat (uniform)" = "#EE7733",
                     "Data-driven" = "#009988",
                     "Prior-driven" = "#CC3311"),
          name = NULL) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "Prior Sensitivity: Informative vs Flat Priors",
          subtitle = sprintf("Flags parameters with |shift| > %.0f%% of informative-mode magnitude",
                             threshold * 100),
          x = "Parameter value (mode)", y = NULL)

      plots$mode_comparison <- .apply_meta(p,
        meta %||% diag_meta(model_name = solved$model$source_file %||% "model"))

      # Relative shift bar chart
      comp_plot2 <- comp_plot
      comp_plot2$param <- factor(comp_plot2$param, levels = comp_plot2$param[order(comp_plot2$rel_shift)])
      p2 <- ggplot2::ggplot(comp_plot2) +
        ggplot2::geom_col(ggplot2::aes(x = rel_shift, y = param, fill = status),
                          alpha = 0.8) +
        ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                           colour = "#CC3311", linewidth = 0.5) +
        ggplot2::scale_fill_manual(
          values = c("Data-driven" = "#009988", "Prior-driven" = "#CC3311"),
          name = NULL) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "Prior Sensitivity: Relative Shift by Parameter",
          subtitle = sprintf("Threshold = %.0f%% (dashed line). Bars to the right are prior-driven.",
                             threshold * 100),
          x = sprintf("|shift| / |informative mode|"), y = NULL)

      plots$relative_shift <- .apply_meta(p2,
        meta %||% diag_meta(model_name = solved$model$source_file %||% "model"))
    }

    .make_result(
      result  = result,
      pass    = all_pass,
      plots   = plots,
      summary = summary_str,
      llm_summary = sprintf(
        "[%s] PriorSensitivity flagged=%d/%d threshold=%.2f",
        if (all_pass) "PASS" else "FAIL",
        n_flagged, nrow(comparisons), threshold)
    )
}
