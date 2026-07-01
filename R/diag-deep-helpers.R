## R/diag-deep-helpers.R
## --------------------------------------------------------------------------
## Deep-Parameter diagnostics — shared helpers.
##
## Provides the parameter taxonomy (`build_deep_spec()`), the structural-class
## auto-classifier, and the A-F "deepness" grade rubric consumed by the
## Deep-Parameter Passport (`diag-deep-passport.R`) and by the genuinely-new
## deep-parameter diagnostics D33 (structural-vs-reduced-form), D34 (policy
## invariance) and D35 (misspecification softness).
##
## See DEEP_PARAMETER_DIAGNOSTICS_PLAN.md for the design and the reconciliation
## against the existing identification block (D1/D20-D30).
## --------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Resolve the @dynhr:deep block from a model object.
# Mirrors diag_expectations()'s triple fallback so the taxonomy is found
# whether or not a pipeline has attached metadata to the model.
# ---------------------------------------------------------------------------
.resolve_deep_block <- function(model) {
  if (is.null(model)) return(list())
  if (!is.null(model$metadata$deep) && length(model$metadata$deep) > 0)
    return(model$metadata$deep)
  if (!is.null(model$deep) && length(model$deep) > 0)
    return(model$deep)
  sf <- model$source_file
  if (!is.null(sf) && length(sf) == 1 && !is.na(sf) && file.exists(sf)) {
    md <- tryCatch(extract_mod_metadata(sf), error = function(e) NULL)
    if (!is.null(md$deep) && length(md$deep) > 0) return(md$deep)
  }
  list()
}


# ---------------------------------------------------------------------------
# Heuristic structural-class auto-classifier (used only when a parameter has
# no declared @dynhr:deep entry).
#
# Honest about its limits: shock standard deviations and AR(1) persistences are
# detectable by naming convention, but a policy-smoothing coefficient such as
# `rho_r` is indistinguishable by name from a shock persistence and *will* be
# mis-tagged `shock` unless declared. The classifier therefore marks every
# auto-decision with `source = "auto"` so the passport can warn.
# ---------------------------------------------------------------------------
.auto_classify_param <- function(pname) {
  p <- tolower(pname)
  # Shock standard deviations / measurement-error scales. Note the trailing
  # underscores: `sig_a`/`sigma_e` are shock sds, but bare `sigma` (risk
  # aversion) and `sd` would be deep primitives, so we do NOT match `^sig`/`^sd`.
  if (grepl("^(sig_|sigma_|sd_|std_|se_|me_)", p) ||
      grepl("(_sd|_std|_se|_sig|_vol)$", p) ||
      grepl("^stderr", p)) {
    return(list(class = "shock", role = "auxiliary", source = "auto"))
  }
  # Shock AR(1) persistences.
  if (grepl("^rho_", p) || grepl("(_rho|_ar1?|_persistence)$", p)) {
    return(list(class = "shock", role = "auxiliary", source = "auto"))
  }
  list(class = "unknown", role = "primitive", source = "auto")
}


# ---------------------------------------------------------------------------
#' Build a deep-parameter taxonomy
#'
#' Normalises the \code{@dynhr:deep} metadata block (or an auto-classified
#' fallback) into a tidy table that the deep-parameter diagnostics and the
#' Deep-Parameter Passport consume.  Each row classifies one parameter by
#' structural \code{class}, \code{role} (primitive / reduced-form / auxiliary),
#' its reduced-form image (if any), whether it is \code{is_deep}, and its
#' \code{partition} (policy / private / auxiliary) for the Lucas-critique
#' invariance test.
#'
#' @param model        A parsed model (from \code{parse_mod}); used to find the
#'   \code{@dynhr:deep} block and the parameter names.  Optional if
#'   \code{param_names} and \code{block} are supplied directly.
#' @param param_names  Character vector of parameter names to classify.
#'   Defaults to \code{model$param_names}, else the names declared in the block.
#' @param block        A pre-parsed deep block (named list keyed by parameter).
#'   Defaults to the block resolved from \code{model}.
#' @return A \code{data.frame} of class \code{"dynhr_deep_spec"} with columns
#'   \code{param, class, role, reduced_form, from, source, is_deep, partition}.
#' @examples
#' \dontrun{
#'   ds <- build_deep_spec(model)
#'   print(ds)
#' }
#' @export
build_deep_spec <- function(model = NULL, param_names = NULL, block = NULL) {
  if (is.null(block)) block <- .resolve_deep_block(model)
  if (is.null(block)) block <- list()

  if (is.null(param_names)) {
    param_names <- if (!is.null(model$param_names) &&
                       length(model$param_names) > 0) model$param_names
                   else names(block)
  }
  param_names <- unique(param_names[nzchar(param_names)])
  if (length(param_names) == 0) {
    out <- data.frame(param = character(0), class = character(0),
                      role = character(0), reduced_form = character(0),
                      from = character(0), source = character(0),
                      is_deep = logical(0), partition = character(0),
                      stringsAsFactors = FALSE)
    class(out) <- c("dynhr_deep_spec", "data.frame")
    return(out)
  }

  rows <- lapply(param_names, function(pn) {
    sp <- block[[pn]]
    if (is.null(sp)) {
      a <- .auto_classify_param(pn)
      sp <- list(class = a$class, role = a$role,
                 reduced_form = NA_character_, from = NA_character_,
                 source = "auto")
    } else {
      sp$source       <- "declared"
      sp$class        <- sp$class %||% "unknown"
      sp$role         <- sp$role %||% "primitive"
      sp$reduced_form <- if (is.null(sp$reduced_form)) NA_character_ else sp$reduced_form
      sp$from         <- if (is.null(sp$from)) NA_character_ else sp$from
    }
    data.frame(
      param        = pn,
      class        = sp$class,
      role         = sp$role,
      reduced_form = sp$reduced_form,
      from         = sp$from,
      source       = sp$source,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, rows)
  df$is_deep   <- !(df$class %in% "shock") & !(df$role %in% "auxiliary")
  df$partition <- ifelse(df$class %in% "policy", "policy",
                         ifelse(df$is_deep, "private", "auxiliary"))
  rownames(df) <- NULL
  class(df) <- c("dynhr_deep_spec", "data.frame")
  df
}


#' Print a deep-parameter taxonomy
#'
#' @param x A \code{\link{build_deep_spec}} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.dynhr_deep_spec <- function(x, ...) {
  n_deep <- sum(x$is_deep)
  n_aux  <- sum(!x$is_deep)
  n_auto <- sum(x$source == "auto")
  cat(sprintf("dynhr deep-parameter taxonomy: %d params (%d deep, %d auxiliary)\n",
              nrow(x), n_deep, n_aux))
  if (n_auto > 0)
    cat(sprintf("  note: %d auto-classified (declare a @dynhr:deep block to override)\n",
                n_auto))
  by_part <- split(x$param, x$partition)
  for (p in names(by_part))
    cat(sprintf("  [%-9s] %s\n", p, paste(by_part[[p]], collapse = ", ")))
  invisible(x)
}


# ---------------------------------------------------------------------------
# A-F "deepness" grade rubric.
#
# A parameter starts at A and is demoted one letter for each *assessed* axis it
# fails. Axes (each TRUE = good, FALSE = fail, NA = not assessed yet):
#   identified  (Q1)  data pins it down, not just the prior   [D1/D20/D23/D24]
#   informed    (Q1)  posterior contracts with data           [D6/D21]
#   structural  (Q2)  a primitive, not a borrowed reduced form [D33]
#   invariant   (Q3)  stable when policy/regime shifts         [D34]
#   robust      (Q4)  survives small misspecification          [D35]
#   calibrated  (Q1)  fixed value is data-consistent/constant  [D36]
#
# (`calibrated` applies to fixed deep parameters; `identified`/`informed`/
# `invariant`/`robust` apply to estimated ones -- the relevant axes are simply
# NA for the other kind.)
#
# Auxiliary (non-deep) parameters are graded "-" (n/a) so the passport keeps
# deep and non-deep parameters visibly distinct.
# ---------------------------------------------------------------------------
.deepness_grade <- function(identified = NA, informed = NA, structural = NA,
                            invariant = NA, robust = NA, calibrated = NA,
                            is_deep = TRUE) {
  if (!isTRUE(is_deep)) return("-")
  vals <- list(identified, informed, structural, invariant, robust, calibrated)
  # A deep parameter with *no* assessed axis (e.g. calibrated and never touched
  # by the data) is graded "?" -- "unassessed" -- rather than a misleading "A".
  assessed <- sum(vapply(vals, function(v) !is.na(unname(v)[1]), logical(1)))
  if (assessed == 0L) return("?")
  # isFALSE() ignores names/attributes, so a named axis element such as
  # c(sigma = FALSE) is still recognised as a failure.
  fails <- sum(vapply(vals, function(v) isFALSE(unname(v)), logical(1)))
  grades <- c("A", "B", "C", "D", "E", "F")
  grades[min(length(grades), 1L + fails)]
}


# ---------------------------------------------------------------------------
# Numeric finite-difference Fisher information at theta from a scalar
# log-likelihood. Returns -Hessian (observed information). Used by the
# information-source axis and (later) by D33/D35. Falls back gracefully.
# ---------------------------------------------------------------------------
.deep_observed_information <- function(loglik_fn, theta, eps = 1e-4) {
  if (is.null(loglik_fn)) return(NULL)
  k <- length(theta)
  H <- matrix(NA_real_, k, k)
  f0 <- tryCatch(loglik_fn(theta), error = function(e) NA_real_)
  if (!is.finite(f0)) return(NULL)
  h <- pmax(abs(theta), 1) * eps
  ll <- function(p) tryCatch(loglik_fn(p), error = function(e) NA_real_)
  for (i in seq_len(k)) {
    for (j in i:k) {
      tpp <- theta; tpp[i] <- tpp[i] + h[i]; tpp[j] <- tpp[j] + h[j]
      tpm <- theta; tpm[i] <- tpm[i] + h[i]; tpm[j] <- tpm[j] - h[j]
      tmp <- theta; tmp[i] <- tmp[i] - h[i]; tmp[j] <- tmp[j] + h[j]
      tmm <- theta; tmm[i] <- tmm[i] - h[i]; tmm[j] <- tmm[j] - h[j]
      val <- (ll(tpp) - ll(tpm) - ll(tmp) + ll(tmm)) / (4 * h[i] * h[j])
      H[i, j] <- H[j, i] <- val
    }
  }
  if (anyNA(H)) return(NULL)
  -H
}


# ---------------------------------------------------------------------------
# Shared credible-interval summariser + overlap test.
#
# D16 (subsample stability) and D34 (policy invariance) both summarise per-
# regime posterior draws into (median, lo, hi) and test CI overlap. These
# helpers centralise that logic so the two diagnostics agree by construction
# and can consume the *same* draws (the pipeline re-estimates each regime once).
# ---------------------------------------------------------------------------
.summarise_draws_ci <- function(draws, label, param_names = NULL, ci_level = 0.90) {
  draws <- as.matrix(draws)
  if (is.null(param_names)) param_names <- colnames(draws)
  if (is.null(param_names)) param_names <- paste0("theta_", seq_len(ncol(draws)))
  a <- (1 - ci_level) / 2
  data.frame(
    param  = param_names,
    median = apply(draws, 2, stats::median),
    lo     = apply(draws, 2, stats::quantile, probs = a,     names = FALSE),
    hi     = apply(draws, 2, stats::quantile, probs = 1 - a, names = FALSE),
    sample = label,
    stringsAsFactors = FALSE)
}

# TRUE if the two intervals do not overlap.
.ci_disjoint <- function(lo1, hi1, lo2, hi2) (hi1 < lo2) || (hi2 < lo1)


# ---------------------------------------------------------------------------
#' Cheap (Laplace) posterior draws for a parameter block
#'
#' Mode-finds the supplied log-likelihood, forms the observed-information
#' Hessian, and draws from the implied Gaussian (Laplace) approximation. This
#' is the inexpensive way to produce the per-regime draws that D34 (policy
#' invariance) and D16 (subsample stability) consume -- re-estimating each
#' regime by a quick optimisation rather than a full MCMC run.
#'
#' @param loglik_fn Function: a named parameter vector -> scalar log-likelihood
#'   (or log-posterior). Should accept and return on the block being estimated.
#' @param start     Named numeric starting values (defines the block + names).
#' @param lower,upper Optional named numeric bounds (enables L-BFGS-B and clamps
#'   the draws).
#' @param n_draws   Number of draws (default 2000).
#' @param ridge     Tikhonov ridge added to the information matrix for a stable
#'   inverse (default 1e-6).
#' @param boundary  How to respect bounds. \code{"transform"} (default) runs the
#'   Laplace approximation in an unconstrained reparameterisation (logit for
#'   two-sided bounds, log for one-sided), so draws are always strictly interior
#'   and acquire an appropriate skew near a boundary instead of piling up on it
#'   -- the right behaviour for weakly-identified parameters that hit their
#'   bound (e.g. an IES at its upper limit). \code{"clamp"} keeps the legacy
#'   Gaussian-in-x approximation with the draws truncated to the bounds.
#' @return A list with \code{draws} (n_draws x k matrix, named columns),
#'   \code{mode}, \code{vcov} and \code{loglik_mode}.
#' @export
deep_laplace_draws <- function(loglik_fn, start, lower = NULL, upper = NULL,
                               n_draws = 2000, ridge = 1e-6,
                               boundary = c("transform", "clamp")) {
  boundary <- match.arg(boundary)
  k  <- length(start)
  nm <- names(start) %||% paste0("theta_", seq_len(k))
  start <- stats::setNames(as.numeric(start), nm)
  lo <- if (is.null(lower)) stats::setNames(rep(-Inf, k), nm) else lower[nm]
  hi <- if (is.null(upper)) stats::setNames(rep( Inf, k), nm) else upper[nm]
  has_bounds <- any(is.finite(lo)) || any(is.finite(hi))

  .laplace_gauss <- function(neg_fn, ll_fn, m0) {
    opt <- tryCatch(suppressWarnings(
      stats::optim(m0, neg_fn, method = "Nelder-Mead",
                   control = list(maxit = 1000))),
      error = function(e) NULL)
    mode <- if (is.null(opt)) m0 else opt$par
    H <- .deep_observed_information(ll_fn, mode)
    if (is.null(H)) H <- diag(length(mode))
    H <- (H + t(H)) / 2 + diag(ridge, length(mode))
    Sig <- tryCatch(solve(H), error = function(e) {
      sv <- svd(H); sv$v %*% ((1 / pmax(sv$d, ridge)) * t(sv$u)) })
    Sig <- (Sig + t(Sig)) / 2
    ev  <- eigen(Sig, symmetric = TRUE)
    Sig <- ev$vectors %*% (pmax(ev$values, ridge) * t(ev$vectors))
    R <- chol(Sig)
    Z <- matrix(stats::rnorm(n_draws * length(mode)), n_draws, length(mode))
    list(draws = sweep(Z %*% R, 2, mode, "+"), mode = mode, vcov = Sig,
         loglik_mode = if (is.null(opt)) NA_real_ else -opt$value)
  }

  if (boundary == "clamp" || !has_bounds) {
    # ---- legacy: Gaussian in x-space, draws truncated to the bounds ----
    g <- .laplace_gauss(
      neg_fn = function(th) { v <- tryCatch(loglik_fn(stats::setNames(th, nm)),
                                            error = function(e) NA_real_)
                              if (is.finite(v)) -v else 1e10 },
      ll_fn  = function(th) loglik_fn(stats::setNames(th, nm)),
      m0 = start)
    draws <- g$draws; colnames(draws) <- nm
    if (any(is.finite(lo))) draws <- sweep(draws, 2, ifelse(is.finite(lo), lo, -Inf), pmax)
    if (any(is.finite(hi))) draws <- sweep(draws, 2, ifelse(is.finite(hi), hi,  Inf), pmin)
    return(list(draws = draws, mode = stats::setNames(g$mode, nm),
                vcov = g$vcov, loglik_mode = g$loglik_mode))
  }

  # ---- boundary-aware: Laplace in an unconstrained reparameterisation ----
  tr <- .make_bound_transform(lo, hi)
  u0 <- tr$to_u(start)
  g  <- .laplace_gauss(
    neg_fn = function(u) { v <- tryCatch(loglik_fn(stats::setNames(tr$to_x(u), nm)),
                                         error = function(e) NA_real_)
                           if (is.finite(v)) -v else 1e10 },
    ll_fn  = function(u) loglik_fn(stats::setNames(tr$to_x(u), nm)),
    m0 = u0)
  draws_x <- t(apply(g$draws, 1, tr$to_x))
  if (ncol(draws_x) != k) draws_x <- matrix(draws_x, ncol = k, byrow = TRUE)
  colnames(draws_x) <- nm
  mode_x <- stats::setNames(tr$to_x(g$mode), nm)
  # delta-method x-space covariance at the mode (for reporting)
  J <- tr$dxdu(g$mode)
  list(draws = draws_x, mode = mode_x,
       vcov = (J * g$vcov) %*% diag(J, k), loglik_mode = g$loglik_mode)
}


# Per-parameter monotone map between a bounded box and unconstrained R^k.
# Two finite bounds -> scaled logit; one finite bound -> log; none -> identity.
.make_bound_transform <- function(lo, hi) {
  k <- length(lo)
  both  <- is.finite(lo) & is.finite(hi)
  lowr  <- is.finite(lo) & !is.finite(hi)
  uppr  <- !is.finite(lo) & is.finite(hi)
  eps   <- 1e-8
  to_u <- function(x) {
    u <- as.numeric(x)
    if (any(both)) {
      z <- (x[both] - lo[both]) / (hi[both] - lo[both])
      z <- pmin(pmax(z, eps), 1 - eps)
      u[both] <- stats::qlogis(z)
    }
    if (any(lowr)) u[lowr] <- log(pmax(x[lowr] - lo[lowr], eps))
    if (any(uppr)) u[uppr] <- log(pmax(hi[uppr] - x[uppr], eps))
    u
  }
  to_x <- function(u) {
    x <- as.numeric(u)
    if (any(both)) x[both] <- lo[both] + (hi[both] - lo[both]) * stats::plogis(u[both])
    if (any(lowr)) x[lowr] <- lo[lowr] + exp(u[lowr])
    if (any(uppr)) x[uppr] <- hi[uppr] - exp(u[uppr])
    x
  }
  dxdu <- function(u) {
    d <- rep(1, k)
    if (any(both)) { p <- stats::plogis(u[both]); d[both] <- (hi[both] - lo[both]) * p * (1 - p) }
    if (any(lowr)) d[lowr] <- exp(u[lowr])
    if (any(uppr)) d[uppr] <- exp(u[uppr])
    d
  }
  list(to_u = to_u, to_x = to_x, dxdu = dxdu)
}
