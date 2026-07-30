## R/sv-rbpf-sbc.R
## --------------------------------------------------------------------------
## Simulation-based calibration (SBC; Talts et al. 2018) for the SV-on-shocks
## RB-PF posterior. For each replication: draw the SV hyperparameters from
## their priors, simulate volatility paths + data from the SV model at that
## draw, sample the posterior with a pseudo-marginal RWMH over the
## `likelihood = "sv_rbpf"` closure, and record the rank of the true value
## among thinned posterior draws. If the whole stack (parser -> factory ->
## RB-PF marginal-likelihood estimate -> pseudo-marginal MCMC) is correctly
## calibrated, the ranks are uniform (tested with sbc_uniformity_test()).
##
## Pseudo-marginal correctness notes (both are load-bearing):
##   - the likelihood closure is built with seed = NULL, so every evaluation
##     draws a fresh particle cloud (a fixed seed would make the estimate
##     deterministic and break the PMMH invariance -- the factory warns);
##   - the CURRENT draw's likelihood ESTIMATE is retained across iterations
##     and never re-evaluated (re-evaluating at the current point targets the
##     wrong stationary distribution).
## --------------------------------------------------------------------------


#' Simulation-based calibration of the SV-on-shocks RB-PF posterior
#'
#' Runs \code{n_repl} SBC replications on a canonical 2-shock linear panel
#' (AR(1) observables with calibrated persistence/scale) whose six SV
#' hyperparameters \code{mu/rho/sigma_eta} per shock are estimated with
#' natural priors (normal / beta / inverse-gamma). Returns per-parameter rank
#' statistics and the package's composite uniformity verdict.
#'
#' Runtime is dominated by \code{n_repl * (n_draws + n_warmup)} RB-PF
#' evaluations (compiled kernel: ~10ms at the defaults), i.e. roughly
#' 40s per replication at the defaults. Use \code{cores > 1} to run
#' replications in parallel (forked; not on Windows).
#'
#' @param n_repl Number of SBC replications (default 100).
#' @param T_obs Sample length of each simulated panel (default 60).
#' @param n_particles RB-PF particles per likelihood evaluation (default 500).
#' @param n_draws,n_warmup Post-warmup / warmup pseudo-marginal RWMH draws
#'   (defaults 8000 / 1500; the warmup adapts a full proposal COVARIANCE
#'   (Haario-style) plus a global scale toward 25\% acceptance, both frozen
#'   for the sampling phase, and is discarded).
#' @param thin_L Number of (approximately independent) thinned draws used for
#'   the rank statistic (default 75; ranks take values 0..thin_L).
#' @param seed Base RNG seed; replication r uses \code{seed + r} (so results
#'   are reproducible and independent of \code{cores}).
#' @param cores Parallel replications via \code{parallel::mclapply}
#'   (default 1 = serial).
#' @return A list: \code{ranks} (\code{n_repl x 6} matrix, columns named by
#'   parameter), \code{uniformity} (from \code{sbc_uniformity_test()}),
#'   \code{accept_rates}, \code{loglik_sd_at_truth} (per-replication PF noise
#'   diagnostic), and \code{settings}.
#' @seealso \code{\link{make_log_posterior_sv_rbpf}}, \code{sbc_uniformity_test()}
#' @export
sv_rbpf_sbc <- function(n_repl = 100L, T_obs = 60L, n_particles = 500L,
                        n_draws = 8000L, n_warmup = 1500L, thin_L = 75L,
                        seed = 20260712L, cores = 1L) {
  ## ---- canonical model: 2 AR(1) observables, calibrated dynamics ---------
  phi <- c(0.8, 0.5); sig <- c(1.0, 0.6)
  mod_txt <- paste(
    "var y1 y2;", "varexo e1 e2;",
    "parameters phi1 phi2 mu_1 rho_1 seta_1 mu_2 rho_2 seta_2;",
    sprintf("phi1=%g; phi2=%g;", phi[1], phi[2]),
    "mu_1=-0.2; rho_1=0.85; seta_1=0.35; mu_2=-0.2; rho_2=0.85; seta_2=0.35;",
    "model(linear); y1 = phi1*y1(-1) + e1; y2 = phi2*y2(-1) + e2; end;",
    sprintf("shocks; var e1; stderr %g; var e2; stderr %g; end;", sig[1], sig[2]),
    "stochastic_volatility;",
    "  var e1; mu = mu_1; rho = rho_1; sigma_eta = seta_1;",
    "  var e2; mu = mu_2; rho = rho_2; sigma_eta = seta_2;",
    "end;",
    "varobs y1 y2;",
    "estimated_params;",
    "  mu_1, normal_pdf, -0.2, 0.35;   rho_1, beta_pdf, 0.85, 0.07;",
    "  seta_1, inv_gamma_pdf, 0.35, 4; mu_2, normal_pdf, -0.2, 0.35;",
    "  rho_2, beta_pdf, 0.85, 0.07;    seta_2, inv_gamma_pdf, 0.35, 4;",
    "end;", sep = "\n")
  mf <- tempfile(fileext = ".mod")
  writeLines(mod_txt, mf)
  model    <- parse_mod(mf)
  compiled <- compile_model(model)
  priors   <- extract_prior_spec(model)
  par_nm   <- priors$name
  n_par    <- length(par_nm)
  prior_sd <- priors$std
  prior_sampler <- .smc_make_prior_sampler(priors)
  obs_vars <- c("y1", "y2")

  simulate_panel <- function(theta, burn = 300L) {
    ## theta order: (mu_1, rho_1, seta_1, mu_2, rho_2, seta_2) via par_nm.
    ## BURN-IN so the retained sample starts from the model's stationary
    ## joint (h, y) distribution -- the filter's stationarity assumption.
    ## (Starting at y_0 = 0 makes the opening artificially quiet and shows
    ## up in SBC as a systematic mu-down / seta-up bias.)
    th <- setNames(theta, par_nm)
    mu_v  <- c(th[["mu_1"]],  th[["mu_2"]])
    rho_v <- c(th[["rho_1"]], th[["rho_2"]])
    se_v  <- c(th[["seta_1"]], th[["seta_2"]])
    n_tot <- burn + T_obs
    Y <- matrix(0, 2, n_tot); yp <- c(0, 0)
    h <- mu_v + rnorm(2, 0, se_v / sqrt(1 - rho_v^2))
    for (t in seq_len(n_tot)) {
      if (t > 1L) h <- mu_v + rho_v * (h - mu_v) + rnorm(2, 0, se_v)
      e  <- rnorm(2, 0, sig * exp(h / 2))
      yp <- phi * yp + e
      Y[, t] <- yp
    }
    out <- t(Y[, (burn + 1L):n_tot, drop = FALSE])
    colnames(out) <- obs_vars
    out
  }

  run_one <- function(r) {
    set.seed(seed + r)
    ## 1. theta* ~ prior (bounded-support draws are the sampler's job).
    theta_star <- prior_sampler()[par_nm]
    ## 2. data | theta*.
    Y <- simulate_panel(theta_star)
    ## 3. pseudo-marginal RWMH on the sv_rbpf posterior (seed = NULL closure).
    lp_fn <- make_log_posterior(model, Y, priors, obs_vars = obs_vars,
                                compiled = compiled, likelihood = "sv_rbpf",
                                n_particles = n_particles)
    ## PF-noise diagnostic at the truth.
    ll_reps <- vapply(1:4, function(k) lp_fn(theta_star)$loglik, numeric(1))
    ll_sd <- stats::sd(ll_reps)

    cur_th <- theta_star
    cur_lp <- lp_fn(cur_th)$logpost
    ## fall back to a prior-mean start if the truth-draw is degenerate
    if (!is.finite(cur_lp)) {
      cur_th <- setNames(priors$mean, par_nm)
      cur_lp <- lp_fn(cur_th)$logpost
    }
    ## Adaptive-covariance (Haario-style) warmup: an ISOTROPIC prior-scaled
    ## proposal mixes catastrophically here (measured IACT ~350 for mu/rho vs
    ## ~10 for seta -- the global scale is pinned by the tight directions), so
    ## the warmup learns the chain covariance and proposes with its Cholesky;
    ## the global scale keeps adapting toward 25% acceptance. Both cov and
    ## scale are FROZEN after warmup (diminishing adaptation is confined to
    ## warmup, so the sampling phase is a valid fixed-kernel chain).
    log_scale <- log(0.5)
    L_prop <- diag(prior_sd, n_par)          # initial proposal Cholesky
    warm_store <- matrix(NA_real_, n_warmup, n_par)
    n_tot <- n_warmup + n_draws
    keep_every <- max(1L, floor(n_draws / thin_L))
    kept <- matrix(NA_real_, nrow = ceiling(n_draws / keep_every), ncol = n_par)
    k_i <- 0L; n_acc <- 0L
    for (it in seq_len(n_tot)) {
      prop <- cur_th + exp(log_scale) * as.numeric(L_prop %*% rnorm(n_par))
      prop_lp <- lp_fn(prop)$logpost
      ## Pseudo-marginal accept: compare against the RETAINED estimate.
      if (is.finite(prop_lp) &&
          log(runif(1)) < (prop_lp - cur_lp)) {
        cur_th <- prop; cur_lp <- prop_lp
        if (it > n_warmup) n_acc <- n_acc + 1L
        acc <- 1
      } else acc <- 0
      if (it <= n_warmup) {
        warm_store[it, ] <- cur_th
        log_scale <- log_scale + (acc - 0.25) / sqrt(it)   # Robbins-Monro
        ## Re-estimate the proposal covariance twice during warmup.
        if (it == floor(n_warmup / 2) || it == n_warmup) {
          cv <- stats::cov(warm_store[seq_len(it), , drop = FALSE])
          cv <- cv + diag(1e-8 + 1e-4 * prior_sd^2, n_par)   # regularise
          Lc <- tryCatch(t(chol(cv)), error = function(e) NULL)
          if (!is.null(Lc)) {
            L_prop <- Lc * (2.38 / sqrt(n_par)) / exp(log_scale)
          }
        }
      } else {
        j <- it - n_warmup
        if (j %% keep_every == 0L) {
          k_i <- k_i + 1L
          kept[k_i, ] <- cur_th
        }
      }
    }
    kept <- kept[seq_len(k_i), , drop = FALSE]
    ranks <- vapply(seq_len(n_par),
                    function(j) sum(kept[, j] < theta_star[j]), numeric(1))
    list(ranks = ranks, accept = n_acc / n_draws, ll_sd = ll_sd, L = k_i)
  }

  reps <- if (cores > 1L) {
    parallel::mclapply(seq_len(n_repl), run_one, mc.cores = cores,
                       mc.preschedule = FALSE)
  } else {
    lapply(seq_len(n_repl), run_one)
  }
  bad <- vapply(reps, function(x) inherits(x, "try-error") || is.null(x$ranks),
                logical(1))
  if (any(bad))
    warning(sprintf("sv_rbpf_sbc: %d/%d replications failed and were dropped.",
                    sum(bad), n_repl), call. = FALSE)
  reps <- reps[!bad]

  ranks_mat <- do.call(rbind, lapply(reps, function(x) x$ranks))
  colnames(ranks_mat) <- par_nm
  list(ranks = ranks_mat,
       uniformity = sbc_uniformity_test(ranks_mat),
       accept_rates = vapply(reps, function(x) x$accept, numeric(1)),
       loglik_sd_at_truth = vapply(reps, function(x) x$ll_sd, numeric(1)),
       settings = list(n_repl = n_repl, T_obs = T_obs,
                       n_particles = n_particles, n_draws = n_draws,
                       n_warmup = n_warmup, thin_L = thin_L, seed = seed))
}
