## R/estimation-passport.R
## --------------------------------------------------------------------------
## Estimation Passport -- a thin, decision-oriented aggregation layer.
##
## This file adds NOTHING to the statistical machinery. dynhr already
## computes:
##   * per-parameter A-F deepness grades via deep_parameter_passport()
##     (R/diag-deep-passport.R), synthesising D1/D6/D20/D23/D33-36;
##   * a model-level FAIL/WARN/PASS executive verdict via
##     format_executive_summary() (R/diag-summary.R).
##
## estimation_passport() is the missing *decision* layer on top of those two
## deliverables. It answers the practitioner's question -- "can I trust this
## posterior for this policy question, and if not, what should I do next?" --
## by:
##   1. SURFACING (never recomputing) the existing per-parameter grades and the
##      existing model verdict;
##   2. ADDING a "suggested next experiment" recommendation engine that maps
##      *failed* axes / diagnostics to concrete actions (add an observable,
##      tighten a prior, reparameterise, check the likelihood, switch sampler);
##   3. Marking any axis / diagnostic that was not run as "not assessed" --
##      never fabricating a grade or a recommendation for an absent check.
##
## It deliberately does not touch any diag-* file or hot path: it reads their
## published outputs ($result$passport, $result$not_assessed, the diagnostic
## badges) and decides.
## --------------------------------------------------------------------------


# ---- recommendation rule book ---------------------------------------------
#
# Each rule is keyed by the *symptom* it detects. A rule fires ONLY when the
# evidence it relies on is actually present in the inputs; an absent diagnostic
# never produces a recommendation (it produces a "not assessed" note instead).
# The mapping is intentionally small and auditable -- one symptom -> one action.

# Failing deepness axis -> a concrete experiment. The axes are exactly the six
# columns deep_parameter_passport() exposes; see R/diag-deep-passport.R.
.passport_axis_actions <- list(
  identified = paste0(
    "weakly identified: add an observable informative for this parameter, ",
    "tighten its prior, or reparameterise to a better-identified combination"),
  informed = paste0(
    "data barely moves the prior: the sample is uninformative here \u2014 ",
    "tighten the prior or add data that loads on this parameter"),
  structural = paste0(
    "value appears borrowed / not primitive: re-derive it structurally or ",
    "declare its reduced-form mapping"),
  invariant = paste0(
    "not policy-invariant across regimes: do not use this parameter for the ",
    "counterfactual \u2014 re-estimate within regime or model the shift"),
  robust = paste0(
    "soft / sandwich-vs-Hessian mismatch: possible misspecification \u2014 ",
    "check the likelihood and measurement-error structure"),
  calibrated = paste0(
    "calibration tension: the fixed value conflicts with the data \u2014 ",
    "free it for estimation or revisit the calibration target"))


# ---- convergence (D5) reader ----------------------------------------------
# Reads the published d5$result$convergence data.frame (cols: param, rhat,
# ess_bulk, ess_tail) and the nuts divergence count. Returns NULL when D5 is
# absent so the caller can mark convergence "not assessed".
.read_convergence <- function(suite,
                              ess_target = 1000,
                              ess_tail_target = 400,
                              rhat_target = 1.01) {
  d5 <- suite[["d5"]]
  conv <- d5$result$convergence
  if (is.null(conv) || !is.data.frame(conv) || nrow(conv) == 0L) return(NULL)
  rhat_bad <- if ("rhat" %in% names(conv))
    conv$param[!is.na(conv$rhat) & conv$rhat > rhat_target] else character(0)
  ess_bad <- if ("ess_bulk" %in% names(conv))
    conv$param[!is.na(conv$ess_bulk) & conv$ess_bulk < ess_target] else character(0)
  ess_tail_bad <- if ("ess_tail" %in% names(conv))
    conv$param[!is.na(conv$ess_tail) & conv$ess_tail < ess_tail_target] else character(0)
  n_div <- d5$result$nuts_meta$n_divergent %||% 0L
  list(
    rhat_bad     = unique(stats::na.omit(rhat_bad)),
    ess_bad      = unique(stats::na.omit(c(ess_bad, ess_tail_bad))),
    n_divergent  = as.integer(n_div),
    rhat_target  = rhat_target,
    ess_target   = ess_target)
}


# ---- the model verdict, surfaced (not recomputed) -------------------------
# Reuses format_executive_summary() read-only and lifts the verdict word it
# already computed. We do NOT re-derive the FAIL/WARN/PASS rule here.
.surface_model_verdict <- function(suite, model_name = NULL) {
  out <- list(verdict = NA_character_, n_pass = NA_integer_,
              n_fail = NA_integer_, n_error = NA_integer_,
              n_info = NA_integer_, available = FALSE)
  has_diag <- any(vapply(suite, function(r) inherits(r, "dynhr_diagnostic"),
                         logical(1)))
  if (!has_diag) return(out)
  lines <- tryCatch(
    format_executive_summary(suite, model_name = model_name),
    error = function(e) NULL)
  if (is.null(lines)) return(out)
  vline <- grep("^OVERALL VERDICT:", lines, value = TRUE)
  if (length(vline)) {
    out$verdict <- if (grepl("FAIL", vline)) "FAIL"
                   else if (grepl("WARN", vline)) "WARN"
                   else if (grepl("PASS", vline)) "PASS"
                   else NA_character_
    out$available <- TRUE
  }
  cline <- grep("PASS .* FAIL .* ERROR .* INFO", lines, value = TRUE)
  if (length(cline)) {
    nums <- as.integer(regmatches(cline[1],
              gregexpr("[0-9]+", cline[1]))[[1]])
    if (length(nums) >= 4L) {
      out$n_pass  <- nums[1]; out$n_fail <- nums[2]
      out$n_error <- nums[3]; out$n_info <- nums[4]
    }
  }
  out
}


# ---------------------------------------------------------------------------
#' Estimation passport: a decision-oriented diagnostic aggregation
#'
#' Synthesises dynhr's existing diagnostic deliverables into a single,
#' decision-oriented answer to "can I trust this posterior, and what should I
#' do next?". It is a thin aggregation layer: it surfaces (never recomputes)
#' the per-parameter A-F grades from \code{\link{deep_parameter_passport}} and
#' the model-level FAIL/WARN/PASS verdict from the executive summary, then adds
#' a recommendation engine that maps the checks that *failed* to concrete next
#' experiments.
#'
#' @param suite A \code{dynhr_diagnostic_suite} (the output of
#'   \code{run_all_diagnostics()}). Optional if \code{passport} is supplied.
#' @param passport A \code{\link{deep_parameter_passport}} result (a
#'   \code{dynhr_diagnostic} whose \code{$result$passport} is the per-parameter
#'   scorecard). If \code{NULL} and a passport is present inside \code{suite}
#'   (element \code{deep_passport} / \code{passport}), it is taken from there.
#' @param model_name Optional character label for the model.
#' @param ess_target,ess_tail_target,rhat_target Convergence thresholds used
#'   only to *read* the published D5 result (no MCMC is rerun). Defaults: 1000
#'   bulk-ESS, 400 tail-ESS, R-hat 1.01.
#' @return An object of class \code{dynhr_estimation_passport}: a list with
#'   \describe{
#'     \item{model_verdict}{the surfaced FAIL/WARN/PASS verdict and counts.}
#'     \item{param_grades}{the per-parameter grade table (subset of the deep
#'       passport scorecard), or \code{NULL} if no passport was available.}
#'     \item{warnings}{model-level warnings (rank deficiency, convergence,
#'       prior-posterior conflict, softness) derived only from present checks.}
#'     \item{recommendations}{the suggested next experiments, each tagged with
#'       the symptom that triggered it.}
#'     \item{not_assessed}{axes / diagnostics that were not run.}
#'   }
#' @seealso \code{\link{deep_parameter_passport}},
#'   \code{summary.dynhr_diagnostic_suite}. For the complementary RUNNER that
#'   computes the paper's seven pre-flight checks directly from
#'   \code{(model, data)} instead of surfacing a pre-existing suite, see
#'   \code{\link{run_estimation_passport}}.
#' @export
estimation_passport <- function(suite = NULL,
                                passport = NULL,
                                model_name = NULL,
                                ess_target = 1000,
                                ess_tail_target = 400,
                                rhat_target = 1.01) {

  # Accept a bare deep_parameter_passport as the sole argument.
  if (is.null(passport) && inherits(suite, "dynhr_diagnostic") &&
      !is.null(suite$result$passport)) {
    passport <- suite
    suite <- NULL
  }
  if (is.null(suite)) suite <- list()
  # A passport may travel inside the suite under a conventional name.
  if (is.null(passport)) {
    for (nm in c("deep_passport", "passport", "deep_parameter_passport")) {
      cand <- suite[[nm]]
      if (inherits(cand, "dynhr_diagnostic") && !is.null(cand$result$passport)) {
        passport <- cand; break
      }
    }
  }

  recommendations <- list()
  warnings        <- character(0)
  not_assessed    <- character(0)

  add_rec <- function(symptom, action, params = NULL) {
    recommendations[[length(recommendations) + 1L]] <<- list(
      symptom = symptom, action = action,
      params  = if (is.null(params)) character(0) else as.character(params))
  }

  # ---- 1. surface the per-parameter grades (no recompute) ----------------
  param_grades <- NULL
  if (!is.null(passport) && !is.null(passport$result$passport)) {
    pg <- passport$result$passport
    keep <- intersect(c("param", "class", "is_deep", "identified", "informed",
                        "structural", "invariant", "robust", "calibrated",
                        "grade"), names(pg))
    param_grades <- pg[, keep, drop = FALSE]
    rownames(param_grades) <- NULL

    # axes the passport itself reported as never assessed
    pna <- passport$result$not_assessed
    if (!is.null(pna)) not_assessed <- union(not_assessed, pna)

    # ---- per-axis recommendations, only for deep params that FAILED -------
    deep_rows <- param_grades[isTRUE_vec(param_grades$is_deep), , drop = FALSE]
    for (axis in names(.passport_axis_actions)) {
      if (!axis %in% names(deep_rows)) next
      failed <- deep_rows$param[identical_false(deep_rows[[axis]])]
      failed <- failed[!is.na(failed)]
      if (length(failed)) {
        add_rec(sprintf("axis:%s", axis), .passport_axis_actions[[axis]], failed)
      }
    }
  } else {
    not_assessed <- union(not_assessed, "per-parameter grades (no passport)")
  }

  # ---- 2. model-level verdict (surfaced) ---------------------------------
  model_verdict <- .surface_model_verdict(suite, model_name = model_name)
  if (!model_verdict$available)
    not_assessed <- union(not_assessed, "model verdict (no diagnostic suite)")

  # ---- 3. model-level warnings + recommendations from present checks -----

  # Rank / local identification (D0 equation rank, D1 local ID weak_params)
  d0 <- suite[["d0"]]
  if (inherits(d0, "dynhr_diagnostic") && isFALSE(d0$pass)) {
    warnings <- c(warnings, "rank deficiency flagged by D0 (equation rank)")
    add_rec("model:rank_deficient",
            "rank-deficient system: re-check the equation/variable count and identifiability before trusting any estimate")
  }
  d1 <- suite[["d1"]]
  weak_id <- unique(c(suite[["d1"]]$result$weak_params,
                      suite[["d20"]]$result$weak_params))
  weak_id <- weak_id[!is.na(weak_id)]
  if (length(weak_id) && is.null(passport)) {
    # If we have no passport to carry the per-param axis recommendation, raise
    # the identification action at model level so it is never silently dropped.
    add_rec("model:weak_identification",
            .passport_axis_actions$identified, weak_id)
  }

  # Convergence (D5)
  conv <- .read_convergence(suite, ess_target, ess_tail_target, rhat_target)
  if (is.null(conv)) {
    not_assessed <- union(not_assessed, "convergence (D5)")
  } else {
    if (length(conv$rhat_bad)) {
      warnings <- c(warnings, sprintf(
        "R-hat > %.2f for: %s", conv$rhat_target,
        paste(conv$rhat_bad, collapse = ", ")))
    }
    if (length(conv$ess_bad)) {
      warnings <- c(warnings, sprintf(
        "low ESS (< %d bulk) for: %s", conv$ess_target,
        paste(conv$ess_bad, collapse = ", ")))
    }
    if (conv$n_divergent > 0L) {
      warnings <- c(warnings, sprintf(
        "%d divergent transition(s) in the sampler", conv$n_divergent))
    }
    if (length(conv$rhat_bad) || length(conv$ess_bad) || conv$n_divergent > 0L) {
      add_rec("model:convergence",
              "MCMC not converged: increase draws, switch/retune the sampler, or reparameterise the stiff directions",
              unique(c(conv$rhat_bad, conv$ess_bad)))
    }
  }

  # Prior-posterior conflict (D6 overlap) at the model level
  uninf <- suite[["d6"]]$result$uninformative
  if (is.null(suite[["d6"]]$result$overlap_scores)) {
    not_assessed <- union(not_assessed, "prior-posterior overlap (D6)")
  } else if (length(uninf) && is.null(passport)) {
    add_rec("model:uninformative", .passport_axis_actions$informed, uninf)
  }

  # Softness / robustness (D35 sandwich-vs-Hessian) at the model level
  soft <- suite[["d35"]]$result$soft_params
  if (!is.null(soft) && length(soft) && is.null(passport)) {
    warnings <- c(warnings, sprintf(
      "robustness (D35) flagged soft parameters: %s",
      paste(soft, collapse = ", ")))
    add_rec("model:softness", .passport_axis_actions$robust, soft)
  }

  # ---- 4. overall verdict line + headline --------------------------------
  n_fail_grades <- if (!is.null(param_grades))
    sum(param_grades$grade %in% c("D", "E", "F")) else 0L

  # "positive evidence" = something was actually assessed and came back clean.
  # A PASS suite verdict counts; so does a passport that graded >= 1 deep
  # parameter (i.e. at least one axis ran) with no D/E/F. Without either we
  # have nothing to vouch for and must say so.
  graded_something <- !is.null(param_grades) &&
    any(param_grades$grade %in% c("A", "B", "C", "D", "E", "F"))
  positive_evidence <- isTRUE(model_verdict$verdict == "PASS") || graded_something

  overall <- if (isTRUE(model_verdict$verdict == "FAIL") ||
                 any(grepl("rank deficiency", warnings))) {
    "DO NOT TRUST"
  } else if (isTRUE(model_verdict$verdict == "WARN") ||
             length(recommendations) > 0L || n_fail_grades > 0L ||
             length(warnings) > 0L) {
    "TRUST WITH CAVEATS"
  } else if (positive_evidence) {
    "TRUSTWORTHY"
  } else {
    "INSUFFICIENT EVIDENCE"
  }

  headline <- switch(overall,
    "DO NOT TRUST" =
      "A pipeline-blocking check failed; the posterior is not usable for policy as-is.",
    "TRUST WITH CAVEATS" =
      "Usable, but address the flagged parameters / recommendations before policy use.",
    "TRUSTWORTHY" =
      "No failures detected across the assessed checks; the posterior looks trustworthy.",
    "INSUFFICIENT EVIDENCE" =
      "Too few diagnostics were available to render a verdict; run more checks.")

  if (overall == "TRUSTWORTHY" && length(recommendations) == 0L) {
    # explicit, honest "nothing to do" rather than an empty list
    add_rec("none",
            "no failing checks: nothing to do \u2014 the estimate looks trustworthy on the assessed axes")
  }

  out <- list(
    overall         = overall,
    headline        = headline,
    model_verdict   = model_verdict,
    param_grades    = param_grades,
    n_fail_grades   = n_fail_grades,
    warnings        = warnings,
    recommendations = recommendations,
    not_assessed    = sort(unique(not_assessed)),
    model_name      = model_name)
  class(out) <- "dynhr_estimation_passport"
  out
}


# small vectorised "is identically TRUE" tolerant of NA / non-logical
isTRUE_vec <- function(v) {
  if (is.null(v)) return(logical(0))
  !is.na(v) & (v == TRUE)
}


# ---------------------------------------------------------------------------
#' Print an estimation passport
#'
#' @param x A \code{dynhr_estimation_passport}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @method print dynhr_estimation_passport
#' @export
print.dynhr_estimation_passport <- function(x, ...) {
  bar <- strrep("=", 72L)
  cat(bar, "\n")
  cat("  ESTIMATION PASSPORT",
      if (!is.null(x$model_name)) sprintf("  (%s)", x$model_name) else "", "\n",
      sep = "")
  cat(bar, "\n")

  mv <- x$model_verdict$verdict
  mv_str <- if (is.na(mv)) "n/a" else mv
  cat(sprintf("  Overall:        %s\n", x$overall))
  cat(sprintf("  Model verdict:  %s", mv_str))
  if (isTRUE(x$model_verdict$available))
    cat(sprintf("   (%d PASS / %d FAIL / %d ERROR / %d INFO)",
                x$model_verdict$n_pass %||% NA, x$model_verdict$n_fail %||% NA,
                x$model_verdict$n_error %||% NA, x$model_verdict$n_info %||% NA))
  cat("\n")
  cat("  ", x$headline, "\n", sep = "")
  cat("\n")

  # ---- per-parameter grades ----
  pg <- x$param_grades
  if (is.null(pg)) {
    cat("  Per-parameter grades: not assessed (no passport supplied)\n\n")
  } else {
    deep <- pg[isTRUE_vec(pg$is_deep), , drop = FALSE]
    cat("  PER-PARAMETER GRADES (deep parameters)\n")
    if (nrow(deep) == 0L) {
      cat("    (no deep parameters)\n")
    } else {
      w <- max(nchar(deep$param), 8L)
      for (i in seq_len(nrow(deep))) {
        flags <- character(0)
        for (a in c("identified", "informed", "structural", "invariant",
                    "robust", "calibrated")) {
          if (a %in% names(deep) && identical_false(deep[[a]][i]))
            flags <- c(flags, a)
        }
        cat(sprintf("    %-*s  [%s]%s\n", w, deep$param[i], deep$grade[i],
                    if (length(flags))
                      sprintf("  failed: %s", paste(flags, collapse = ", "))
                    else ""))
      }
    }
    cat("\n")
  }

  # ---- warnings ----
  if (length(x$warnings)) {
    cat("  MODEL WARNINGS\n")
    for (w in x$warnings) cat("    - ", w, "\n", sep = "")
    cat("\n")
  }

  # ---- recommendations ----
  cat("  SUGGESTED NEXT EXPERIMENT(S)\n")
  if (length(x$recommendations) == 0L) {
    cat("    (none)\n")
  } else {
    for (r in x$recommendations) {
      tag <- if (length(r$params))
        sprintf(" [%s]", paste(utils::head(r$params, 8), collapse = ", "))
      else ""
      cat("    - ", r$action, tag, "\n", sep = "")
    }
  }
  cat("\n")

  # ---- not assessed ----
  if (length(x$not_assessed)) {
    cat("  NOT ASSESSED: ", paste(x$not_assessed, collapse = "; "), "\n",
        sep = "")
  }
  cat(bar, "\n")
  invisible(x)
}
