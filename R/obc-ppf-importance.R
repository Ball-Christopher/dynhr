## R/obc-ppf-importance.R
## --------------------------------------------------------------------------
## Phase C: PPF posterior importance re-weighting.
##
## Implements the Dynare 7 posterior importance re-weighting approach:
##   w_i = exp(loglik_PPF(theta_i) - loglik_PKF(theta_i))
## on a thinned subset of PKF posterior draws, returning normalized weights,
## ESS, and a verdict on PKF adequacy.
##
## Provides:
##   ppf_reweight_posterior()  -- importance re-weighting (exported)
## --------------------------------------------------------------------------


#' PPF posterior importance re-weighting
#'
#' Given existing PKF posterior draws (chains or a mode result), re-weights
#' them using PPF log-likelihoods. The importance weights are:
#'
#'   w_i = exp(loglik_PPF(theta_i) - loglik_PKF(theta_i))
#'
#' Normalised in log-space to avoid overflow. ESS / n near 1.0 indicates the
#' PKF approximation is adequate; low ESS/n signals PKF-induced posterior bias.
#'
#' PKF log-likelihoods are re-evaluated (chains do not carry per-draw logliks).
#'
#' @param chains         Matrix of posterior draws (n_draws x n_params) or a
#'                       named list with \code{$chains} or \code{$mode} component. If a
#'                       single-row matrix, treated as the mode.
#' @param model          dynhr_mod
#' @param compiled       dynhr_compiled
#' @param data           Observation matrix (T x n_obs or n_obs x T)
#' @param prior_spec     Prior specification data.frame
#' @param obs_vars       Character vector of observed variable names
#' @param specs          OBC spec list (from obc_parse_tags), or NULL to parse
#' @param n_thin         Thinning factor -- use every n_thin-th draw (default 10)
#' @param n_particles    Number of PPF particles per draw (default 2000)
#' @param me_variance    Measurement error variance (must be > 0; default 1e-4)
#' @param seed           Integer RNG seed (default 42L)
#' @return Named list with:
#'   \item{log_weights}{length-n vector of log importance weights}
#'   \item{weights}{length-n vector of normalized importance weights}
#'   \item{ess}{effective sample size}
#'   \item{ess_fraction}{ESS / n (near 1 = PKF adequate; < 0.3 = PKF biased)}
#'   \item{loglik_pkf}{length-n vector of re-evaluated PKF log-likelihoods}
#'   \item{loglik_ppf}{length-n vector of PPF log-likelihoods}
#'   \item{n_draws}{number of thinned draws used}
#'   \item{verdict}{character string diagnostic message}
#' @export
ppf_reweight_posterior <- function(chains, model, compiled, data,
                                    prior_spec, obs_vars,
                                    specs       = NULL,
                                    n_thin      = 10L,
                                    n_particles = 2000L,
                                    me_variance = 1e-4,
                                    seed        = 42L) {

  ## --- Guards ---------------------------------------------------------------
  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance <= 0)
    stop("ppf_reweight_posterior: me_variance must be a finite positive scalar.")

  ## --- Extract draw matrix --------------------------------------------------
  draw_mat <- .ppf_extract_draws(chains)
  if (is.null(draw_mat) || nrow(draw_mat) == 0L)
    stop("ppf_reweight_posterior: could not extract draws from 'chains'.")

  ## Thin
  n_all  <- nrow(draw_mat)
  idx    <- seq(1L, n_all, by = max(1L, as.integer(n_thin)))
  draws  <- draw_mat[idx, , drop = FALSE]
  n      <- nrow(draws)

  if (n == 0L)
    stop("ppf_reweight_posterior: no draws remaining after thinning.")

  ## --- Set up posterior functions -------------------------------------------
  if (is.null(specs)) specs <- obc_parse_tags(model)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  if (is.data.frame(data)) data <- as.matrix(data)
  Y <- if (nrow(data) == length(obs_vars)) data else t(data)

  lp_pkf <- make_log_posterior_obc_pkf(
    model, t(Y), prior_spec, obs_vars, compiled,
    specs = specs, me_variance = me_variance
  )

  ## --- Re-evaluate PKF and PPF at each thinned draw ------------------------
  loglik_pkf <- numeric(n)
  loglik_ppf <- numeric(n)

  for (i in seq_len(n)) {
    theta_i <- draws[i, ]
    if (!is.null(colnames(draws)))
      names(theta_i) <- colnames(draws)

    ## PKF
    res_pkf       <- lp_pkf(theta_i)
    loglik_pkf[i] <- if (!is.null(res_pkf$loglik)) res_pkf$loglik else -Inf

    ## PPF (fresh cache per draw, fixed seed offset by i for independent noise)
    if (is.finite(loglik_pkf[i])) {
      params <- .apply_theta_to_params(model, theta_i)

      ss_result <- solve_steady_state(model, compiled, params, verbose = FALSE)
      if (!is.null(ss_result) && ss_result$converged) {
        sys_cache_i <- cache_system_structure(compiled)
        ## Re-derive SSM-computed params for a consistent linearization point
        ## (no-op for non-SSM-parameter models; Tier 13 #1).
        params      <- ss_result$params %||% params
        sys_i       <- extract_system_matrices_fast(sys_cache_i, ss_result$ss, params)
        dr_i        <- .solve_from_system(sys_i, model, compiled, ss_result$ss, params, FALSE)
        if (!is.null(dr_i) && dr_i$bk_satisfied) {
          endo_i    <- model$var_names
          obs_idx_i <- match(obs_vars, endo_i)
          rc_i      <- new.env(parent = emptyenv(), hash = TRUE)
          obc_ensure_policy(0L, rc_i, sys_i, dr_i, specs, obs_idx_i)

          pf_i        <- ppf_likelihood(
            Y, dr_i, rc_i, sys_i, model, params, obs_vars, specs,
            obs_idx     = obs_idx_i,
            N           = n_particles,
            me_variance = me_variance,
            seed        = if (!is.null(seed)) seed + i else NULL
          )
          loglik_ppf[i] <- if (!is.null(pf_i$loglik)) pf_i$loglik else -Inf
        } else {
          loglik_ppf[i] <- -Inf
        }
      } else {
        loglik_ppf[i] <- -Inf
      }
    } else {
      loglik_ppf[i] <- -Inf
    }
  }

  ## --- Importance weights (log-space normalisation) -------------------------
  ## Restrict to draws where both are finite
  ok          <- is.finite(loglik_pkf) & is.finite(loglik_ppf)
  log_w_raw   <- rep(-Inf, n)
  log_w_raw[ok] <- loglik_ppf[ok] - loglik_pkf[ok]

  if (sum(ok) == 0L) {
    warning("ppf_reweight_posterior: no draws with finite PKF and PPF logliks.")
    return(list(
      log_weights  = log_w_raw,
      weights      = rep(1 / n, n),
      ess          = 1,
      ess_fraction = 1 / n,
      loglik_pkf   = loglik_pkf,
      loglik_ppf   = loglik_ppf,
      n_draws      = n,
      verdict      = "No valid draws -- cannot assess PKF adequacy."
    ))
  }

  ## Normalize in log-space
  log_w_max  <- max(log_w_raw[ok])
  log_w_norm <- log_w_raw - log_w_max
  w_unnorm   <- ifelse(is.finite(log_w_norm), exp(log_w_norm), 0)
  w_norm     <- w_unnorm / sum(w_unnorm)

  ess          <- .smc_ess(log_w_raw[ok])
  ess_fraction <- ess / sum(ok)

  verdict <- if (ess_fraction >= 0.9) {
    sprintf("PKF adequate (ESS/n = %.2f >= 0.90): PKF and PPF posteriors agree.", ess_fraction)
  } else if (ess_fraction >= 0.5) {
    sprintf("PKF moderate (ESS/n = %.2f): mild posterior discrepancy.", ess_fraction)
  } else {
    sprintf("PKF inadequate (ESS/n = %.2f < 0.50): PPF posterior differs materially.", ess_fraction)
  }

  list(
    log_weights  = log_w_raw,
    weights      = w_norm,
    ess          = ess,
    ess_fraction = ess_fraction,
    loglik_pkf   = loglik_pkf,
    loglik_ppf   = loglik_ppf,
    n_draws      = n,
    verdict      = verdict
  )
}


## Helper: extract draw matrix from various input formats
#' @noRd
.ppf_extract_draws <- function(chains) {
  if (is.matrix(chains) || is.data.frame(chains))
    return(as.matrix(chains))

  ## List input. Probe the CANONICAL draw locations first, with EXACT matching
  ## ([[ exact = TRUE ]] -- `$chain` would otherwise partial-match `$chains`):
  ##   dynhr_chains$chain ; dynhr_estimation_result$chains$chain ;
  ##   dynhr_posterior_result$pooled_draws. Legacy bare $chains/$draws matrices
  ## and a $mode vector are still honoured.
  if (is.list(chains)) {
    g <- function(x, nm) if (is.list(x)) x[[nm, exact = TRUE]] else NULL
    for (m in list(g(chains, "chain"),
                   g(g(chains, "chains"), "chain"),
                   g(chains, "pooled_draws"),
                   g(chains, "chains"),
                   g(chains, "draws")))
      if (is.matrix(m)) return(m)
    mo <- g(chains, "mode")
    if (is.numeric(mo))
      return(matrix(mo, nrow = 1L, dimnames = list(NULL, names(mo))))
    ## Last resort: a SINGLE unambiguous matrix field. Refuse to silently
    ## rbind several (e.g. the draws AND the data matrix) into a fake cloud.
    mats <- Filter(is.matrix, chains)
    if (length(mats) == 1L) return(mats[[1L]])
    if (length(mats) > 1L)
      stop(".ppf_extract_draws: ambiguous draw input -- multiple matrix fields (",
           paste(names(mats), collapse = ", "),
           "). Pass the draw matrix (e.g. chains$chain) explicitly.", call. = FALSE)
  }
  NULL
}
