## R/ms-spec.R
## --------------------------------------------------------------------------
## Markov-switching DSGE specification (shock-variance switching only).
##
## ms_dsge_spec() -- constructor and validator for the MS setup:
##   - transition matrix P (h x h, rows sum to 1)
##   - per-regime shock-scale vectors (length n_exo each)
##   - initial regime distribution (default: ergodic dist of P)
##
## The coupled MS perturbation solver (structural-parameter switching) is
## deferred; only shock variances switch across regimes here.
## --------------------------------------------------------------------------


#' Create a Markov-switching DSGE specification
#'
#' Specifies a heteroskedastic Markov-switching model where only shock
#' variances differ across regimes.  The structural parameters (and therefore
#' \code{ghx}, \code{ghu}, \code{TT}, \code{ZZ}) are identical across regimes
#' and are computed once by the standard \code{solve_perturbation()} pipeline.
#' Only the per-regime shock covariance \eqn{\Sigma_e^{(s)}} differs.
#'
#' @param n_regimes  Integer >= 2.  Number of Markov regimes.
#' @param transition Numeric h x h matrix with non-negative entries whose rows
#'   sum to 1 (row-stochastic; entry \code{P[i,j]} is the probability of
#'   transitioning from regime \code{i} to regime \code{j}).
#' @param shock_scales Named list of length \code{n_regimes}.  Each element
#'   is a named numeric vector of length \code{n_exo} giving the shock
#'   standard-deviation scale factors for that regime.  Regime 1 is
#'   conventionally fixed to all 1s for identification, but this is not
#'   enforced (any strictly positive values are valid).  Names must match the
#'   model's exogenous variable names (checked when the filter is called).
#' @param pi0  Optional numeric vector of length \code{n_regimes}: initial
#'   regime probability distribution.  Defaults to the ergodic (stationary)
#'   distribution of \code{transition}.  Must sum to 1 and be non-negative.
#'
#' @return An object of class \code{"ms_dsge_spec"}.
#' @export
ms_dsge_spec <- function(n_regimes,
                          transition,
                          shock_scales,
                          pi0 = NULL) {

  ## ---- validate n_regimes ------------------------------------------------
  n_regimes <- as.integer(n_regimes)
  if (length(n_regimes) != 1L || is.na(n_regimes) || n_regimes < 2L)
    stop("ms_dsge_spec: n_regimes must be an integer >= 2.", call. = FALSE)

  ## ---- validate transition matrix ----------------------------------------
  if (!is.matrix(transition) ||
      nrow(transition) != n_regimes || ncol(transition) != n_regimes)
    stop(sprintf(
      "ms_dsge_spec: transition must be a %d x %d matrix.", n_regimes, n_regimes),
      call. = FALSE)
  if (any(transition < 0))
    stop("ms_dsge_spec: transition matrix must have non-negative entries.",
         call. = FALSE)
  row_sums <- rowSums(transition)
  if (any(abs(row_sums - 1) > 1e-10))
    stop(sprintf(
      "ms_dsge_spec: transition matrix rows must sum to 1; max deviation: %.2e",
      max(abs(row_sums - 1))), call. = FALSE)

  ## ---- validate shock_scales ---------------------------------------------
  if (!is.list(shock_scales) || length(shock_scales) != n_regimes)
    stop(sprintf(
      "ms_dsge_spec: shock_scales must be a list of length %d (one per regime).",
      n_regimes), call. = FALSE)
  n_exo <- length(shock_scales[[1L]])
  for (s in seq_len(n_regimes)) {
    ss <- shock_scales[[s]]
    if (!is.numeric(ss) || length(ss) != n_exo)
      stop(sprintf(
        "ms_dsge_spec: shock_scales[[%d]] must be a numeric vector of length %d.",
        s, n_exo), call. = FALSE)
    if (any(!is.finite(ss)) || any(ss <= 0))
      stop(sprintf(
        "ms_dsge_spec: shock_scales[[%d]] must have strictly positive finite values.",
        s), call. = FALSE)
  }

  ## ---- initial distribution pi0 (default: ergodic) -----------------------
  if (is.null(pi0)) {
    pi0 <- .ms_ergodic_dist(transition)
  } else {
    if (!is.numeric(pi0) || length(pi0) != n_regimes)
      stop(sprintf("ms_dsge_spec: pi0 must be a numeric vector of length %d.",
                   n_regimes), call. = FALSE)
    if (any(pi0 < 0) || !is.finite(sum(pi0)) || abs(sum(pi0) - 1) > 1e-10)
      stop("ms_dsge_spec: pi0 must be non-negative and sum to 1.", call. = FALSE)
    pi0 <- pi0 / sum(pi0)   # renormalise for floating-point tolerance
  }

  ## ---- build names for regimes -------------------------------------------
  regime_names <- names(shock_scales)
  if (is.null(regime_names))
    regime_names <- paste0("regime", seq_len(n_regimes))

  structure(
    list(
      n_regimes   = n_regimes,
      transition  = transition,
      shock_scales = shock_scales,
      pi0         = pi0,
      regime_names = regime_names,
      n_exo       = n_exo
    ),
    class = c("ms_dsge_spec", "list")
  )
}


#' @export
#' @noRd
print.ms_dsge_spec <- function(x, ...) {
  cat(sprintf("<ms_dsge_spec>  %d regimes, %d shocks\n",
              x$n_regimes, x$n_exo))
  cat("  Transition matrix (rows = from):\n")
  print(round(x$transition, 4))
  cat("  Ergodic dist:", paste(round(x$pi0, 4), collapse = ", "), "\n")
  cat("  Shock scales by regime:\n")
  for (s in seq_len(x$n_regimes)) {
    cat(sprintf("    %s: %s\n", x$regime_names[s],
                paste(round(x$shock_scales[[s]], 3), collapse = ", ")))
  }
  invisible(x)
}


## --------------------------------------------------------------------------
## ms_struct_spec() -- structural Markov-switching spec
##
## Carries the per-regime parameter lists, optional per-regime steady states,
## and transition matrix for a structural MS-DSGE model.  This is the
## companion spec to ms_dsge_spec (which handles shock-variance-only switching).
## Used by estimation_context() to route make_log_posterior() to the
## solve_ms_perturbation() + ms_kim_filter_struct() path.
## --------------------------------------------------------------------------

#' Create a structural Markov-switching DSGE specification
#'
#' Specifies a Markov-switching DSGE model where structural parameters (and
#' therefore the decision rules \code{ghx_s}, \code{ghu_s}) differ across
#' regimes.  Used with \code{\link{estimation_context}} to route the
#' estimation pipeline through \code{\link{solve_ms_perturbation}} and
#' \code{\link{ms_kim_filter_struct}} at each posterior draw.
#'
#' @param n_regimes      Integer >= 2.  Number of Markov regimes.
#' @param transition     Numeric h x h row-stochastic transition matrix.
#'   Entry \code{P[i,j]} is the probability of transitioning from regime
#'   \code{i} to regime \code{j}.
#' @param params_by_regime  List of length \code{n_regimes}.  Each element
#'   is a named numeric parameter vector for that regime.  Regime-specific
#'   parameters are used as the \code{params_by_regime} argument to
#'   \code{\link{solve_ms_perturbation}}.  Estimated parameters (from
#'   \code{estimated_params}) override the corresponding entry in each
#'   regime's vector at each posterior draw via \code{.apply_theta_to_params}.
#' @param ss_by_regime   Optional list of length \code{n_regimes}.  Each
#'   element is a named numeric steady-state vector (as returned by
#'   \code{solve_steady}).  When \code{NULL} (default), the steady state is
#'   solved at each posterior draw using the regime-specific parameters.
#' @param pi0            Optional numeric vector of length \code{n_regimes}:
#'   initial regime probability distribution.  Defaults to the ergodic
#'   (stationary) distribution of \code{transition}.  Must sum to 1 and be
#'   non-negative.
#'
#' @return An object of class \code{"ms_struct_spec"}.
#' @seealso \code{\link{estimation_context}}, \code{\link{solve_ms_perturbation}},
#'   \code{\link{ms_kim_filter_struct}}
#' @export
ms_struct_spec <- function(n_regimes,
                            transition,
                            params_by_regime,
                            ss_by_regime = NULL,
                            pi0 = NULL) {

  ## ---- validate n_regimes ------------------------------------------------
  n_regimes <- as.integer(n_regimes)
  if (length(n_regimes) != 1L || is.na(n_regimes) || n_regimes < 2L)
    stop("ms_struct_spec: n_regimes must be an integer >= 2.", call. = FALSE)

  ## ---- validate transition matrix ----------------------------------------
  if (!is.matrix(transition) ||
      nrow(transition) != n_regimes || ncol(transition) != n_regimes)
    stop(sprintf(
      "ms_struct_spec: transition must be a %d x %d matrix.", n_regimes, n_regimes),
      call. = FALSE)
  if (any(transition < 0))
    stop("ms_struct_spec: transition matrix must have non-negative entries.",
         call. = FALSE)
  row_sums <- rowSums(transition)
  if (any(abs(row_sums - 1) > 1e-10))
    stop(sprintf(
      "ms_struct_spec: transition matrix rows must sum to 1; max deviation: %.2e",
      max(abs(row_sums - 1))), call. = FALSE)

  ## ---- validate params_by_regime -----------------------------------------
  if (!is.list(params_by_regime) || length(params_by_regime) != n_regimes)
    stop(sprintf(
      "ms_struct_spec: params_by_regime must be a list of length %d.",
      n_regimes), call. = FALSE)
  for (s in seq_len(n_regimes)) {
    ps <- params_by_regime[[s]]
    if (!is.numeric(ps) || is.null(names(ps)))
      stop(sprintf(
        "ms_struct_spec: params_by_regime[[%d]] must be a named numeric vector.",
        s), call. = FALSE)
  }

  ## ---- validate ss_by_regime (optional) ----------------------------------
  if (!is.null(ss_by_regime)) {
    if (!is.list(ss_by_regime) || length(ss_by_regime) != n_regimes)
      stop(sprintf(
        "ms_struct_spec: ss_by_regime must be a list of length %d or NULL.",
        n_regimes), call. = FALSE)
    for (s in seq_len(n_regimes)) {
      ss_s <- ss_by_regime[[s]]
      if (!is.numeric(ss_s) || is.null(names(ss_s)))
        stop(sprintf(
          "ms_struct_spec: ss_by_regime[[%d]] must be a named numeric vector.",
          s), call. = FALSE)
    }
  }

  ## ---- initial distribution pi0 (default: ergodic) -----------------------
  if (is.null(pi0)) {
    pi0 <- .ms_ergodic_dist(transition)
  } else {
    if (!is.numeric(pi0) || length(pi0) != n_regimes)
      stop(sprintf("ms_struct_spec: pi0 must be a numeric vector of length %d.",
                   n_regimes), call. = FALSE)
    if (any(pi0 < 0) || !is.finite(sum(pi0)) || abs(sum(pi0) - 1) > 1e-10)
      stop("ms_struct_spec: pi0 must be non-negative and sum to 1.", call. = FALSE)
    pi0 <- pi0 / sum(pi0)
  }

  ## ---- build regime names ------------------------------------------------
  regime_names <- names(params_by_regime)
  if (is.null(regime_names))
    regime_names <- paste0("regime", seq_len(n_regimes))

  structure(
    list(
      n_regimes        = n_regimes,
      transition       = transition,
      params_by_regime = params_by_regime,
      ss_by_regime     = ss_by_regime,
      pi0              = pi0,
      regime_names     = regime_names
    ),
    class = c("ms_struct_spec", "list")
  )
}


#' @export
#' @noRd
print.ms_struct_spec <- function(x, ...) {
  cat(sprintf("<ms_struct_spec>  %d regimes (structural switching)\n", x$n_regimes))
  cat("  Transition matrix (rows = from):\n")
  print(round(x$transition, 4))
  cat("  Ergodic dist:", paste(round(x$pi0, 4), collapse = ", "), "\n")
  cat("  ss_by_regime:", if (is.null(x$ss_by_regime)) "NULL (solved per draw)"
                         else sprintf("pre-supplied (%d entries)", x$n_regimes), "\n")
  cat("  Regime parameter vectors:\n")
  for (s in seq_len(x$n_regimes)) {
    nms <- names(x$params_by_regime[[s]])
    cat(sprintf("    %s: %d params (%s)\n", x$regime_names[s], length(nms),
                paste(head(nms, 5), collapse = ", ")))
  }
  invisible(x)
}


## Internal: power-iteration ergodic distribution.
## (Replicates .compute_ergodic_dist_ms from R/diag-pre-d28-regime-ident.R
## but is self-contained here so ms-spec.R has no ordering dependency.)
## @noRd
.ms_ergodic_dist <- function(P, tol = 1e-14, max_iter = 10000L) {
  n  <- nrow(P)
  pi <- rep(1 / n, n)
  for (i in seq_len(max_iter)) {
    pi_new <- drop(pi %*% P)
    if (max(abs(pi_new - pi)) < tol) break
    pi <- pi_new
  }
  pi_new / sum(pi_new)
}


## Internal: build per-regime (QQ, HH, SS) covariance matrices.
## Given common (RR, DD, Sigma_e) and a list of shock-scale vectors,
## return a length-h list of lists(QQ, HH, SS, Sigma_e_s).
## @noRd
.ms_build_regime_covs <- function(RR, DD, Sigma_e, shock_scales) {
  h <- length(shock_scales)
  covs <- vector("list", h)
  for (s in seq_len(h)) {
    sc    <- shock_scales[[s]]
    Se_s  <- Sigma_e * outer(sc, sc)   # diag(sc) %*% Sigma_e %*% diag(sc)
    covs[[s]] <- list(
      QQ      = tcrossprod(RR %*% Se_s, RR),
      HH      = tcrossprod(DD %*% Se_s, DD),
      SS      = RR %*% Se_s %*% t(DD),
      Sigma_e = Se_s
    )
  }
  covs
}
