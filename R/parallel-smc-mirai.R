## R/parallel-smc-mirai.R
## --------------------------------------------------------------------------
## run_smc_mirai() -- correct parallel SMC for compiled DSGE models.
##
## The OLD parallel path (.smc_pool_setup + .smc_pmap, in parallel-mirai.R)
## shipped the log-posterior CLOSURE, which captures the compiled model's
## non-serialisable C++ external pointers -- fine for a pure-R toy posterior,
## but it fails on a real model ("replacement has length zero"). Here the pool
## is provisioned with .mirai_pool_init() (each daemon loads dynhr, maps the
## shared Y, compiles the model once, and builds .worker_lp), and the
## embarrassingly-parallel work -- stage-0 particle evaluation and each stage's
## RWMH mutation -- is farmed out via mirai_map task bodies that fetch
## .worker_lp from the daemon globalenv. The SMC bookkeeping (tempering,
## weights, ESS/resampling, marginal likelihood, adaptive mutation scaling)
## runs on the host, reusing the helpers in sampler-smc.R.
##
## Package-option shipping: this file starts NO daemon pool of its own -- both
## provisioning paths below (.mirai_pool_closure / .mirai_pool_init, in
## parallel-mirai.R) replay the host's `dynhr_set_options()` state on every
## daemon via `.dynhr_daemon_state()` / `.dynhr_daemon_apply()`, so a
## power_posterior / me_variance / debug_kf_errors set here reaches the workers.
## Do not add a bare `mirai::daemons()` + `everywhere()` here without also
## shipping that snapshot -- `.dynhr_opts` is a namespace-private environment,
## not base `options()`, and does not cross a process boundary by itself.
## --------------------------------------------------------------------------


#' Run SMC in parallel on a mirai daemon pool (compile-per-daemon).
#'
#' @param parsed_model parsed dynare model.
#' @param Y observation matrix.
#' @param prior_spec prior spec data.frame.
#' @param obs_names observed variable names.
#' @param n_particles number of particles.
#' @param ess_target target ESS ratio for adaptive tempering / resampling.
#' @param n_mh_steps RWMH mutation steps per stage.
#' @param mh_scale_factor initial mutation scale (/sqrt(d)); adapts thereafter.
#' @param mut_target target mutation acceptance rate (Herbst-Schorfheide).
#' @param lambda_schedule optional fixed tempering schedule (NULL = adaptive).
#' @param mixture_weights NULL (default) or a numeric length-3 vector of
#'   non-negative mixture weights, normalised internally, for the
#'   Herbst & Schorfheide (2014 JAE; 2015 book ch. 5) / FRBNY SMC.jl
#'   3-component mutation proposal -- see dynhr_smc() for the component
#'   definitions and the asymmetric-proposal MH correction. NULL keeps the
#'   single full-covariance random-walk mutation exactly as before.
#' @param seed_base base RNG seed.
#' @param n_cores worker count (NULL = auto).
#' @param me_variance measurement-error variance.
#' @param me_extra n_obs x T matrix of per-period extra ME variances (filter_tunes).
#' @param log_post_fn optional pre-built log-posterior closure. When supplied,
#'   the pool ships this closure once via \code{\link{.mirai_pool_closure}}
#'   instead of recompiling a standard Gaussian posterior per daemon via
#'   \code{\link{.mirai_pool_init}} -- the path for OBC/PKF and cumulant
#'   models, where \code{parsed_model}/\code{Y}/\code{obs_names} may be NULL.
#' @param verbose print per-stage progress.
#' @return a list matching dynhr_smc() output, ready for new_dynhr_chains().
#' @noRd
run_smc_mirai <- function(
    parsed_model = NULL, Y = NULL, prior_spec, obs_names = NULL,
    n_particles     = 2000L,
    ess_target      = 0.5,
    n_mh_steps      = 1L,
    mh_scale_factor = 0.5,
    mut_target      = 0.25,
    lambda_schedule = NULL,
    mixture_weights = NULL,
    seed_base       = 1L,
    n_cores         = NULL,
    me_variance     = 0,
    me_extra        = NULL,
    shock_scale     = NULL,
    system_priors   = NULL,
    lik_init        = "auto",
    tpf_options     = list(),
    gradient_policy = "auto",
    log_post_fn     = NULL,
    ctx             = NULL,
    verbose         = TRUE
) {
  ## Unpack ctx fields when provided (ctx wins over individual args).
  if (!is.null(ctx) && inherits(ctx, "dynhr_estimation_context")) {
    me_variance     <- ctx$me_variance
    me_extra        <- ctx$me_extra
    shock_scale     <- ctx$shock_scale
    system_priors   <- ctx$system_priors
    lik_init        <- ctx$lik_init        %||% "auto"
    tpf_options     <- ctx$tpf_options     %||% list()
    gradient_policy <- ctx$gradient_policy %||% "auto"
  }
  ## Validate and normalise the Herbst-Schorfheide mixture weights (see
  ## dynhr_smc() for the component definitions). NULL keeps the original
  ## single-component full-covariance random-walk mutation exactly as before.
  if (!is.null(mixture_weights)) {
    if (!is.numeric(mixture_weights) || length(mixture_weights) != 3L)
      stop("mixture_weights must be NULL or a numeric vector of length 3")
    if (any(mixture_weights < 0) || !any(mixture_weights > 0))
      stop("mixture_weights must be non-negative with at least one positive entry")
    mixture_weights <- mixture_weights / sum(mixture_weights)
  }
  use_mixture <- !is.null(mixture_weights)

  n_cores <- .mirai_n_cores(n_cores, n_particles)
  if (verbose)
    cat(sprintf("  Parallel SMC (mirai): %d particles on %d daemons\n",
                n_particles, n_cores))

  t_init <- proc.time()
  sh <- NULL
  if (!is.null(log_post_fn)) {
    .mirai_pool_closure(n_cores, log_post_fn)
  } else {
    sh <- .mirai_pool_init(n_cores, parsed_model, Y, prior_spec, obs_names,
                           me_variance, me_extra = me_extra,
                           shock_scale = shock_scale,
                           system_priors = system_priors,
                           lik_init = lik_init,
                           tpf_options = tpf_options)
  }
  on.exit({ mirai::daemons(NULL); if (!is.null(sh)) rm(sh) }, add = TRUE)
  if (verbose)
    cat(sprintf("  Daemon init: %.1f sec (load + compile + lp_fn)\n",
                (proc.time() - t_init)[["elapsed"]]))

  prior_sampler <- .smc_make_prior_sampler(prior_spec)
  par_names <- prior_spec$name
  t_start <- Sys.time()
  n_eval  <- 0L
  ll_floor <- -1e300

  ## ---- Stage 0: draw particles on the host, evaluate on the daemons --------
  set.seed(seed_base)
  theta_mat <- t(vapply(seq_len(n_particles), function(i) prior_sampler(),
                        numeric(length(par_names))))
  colnames(theta_mat) <- par_names
  d <- ncol(theta_mat)

  ## theta_mat / par_names are FREE variables in eval_task -> pass via `...`.
  eval_task <- function(i) {
    lpf <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
    th  <- theta_mat[i, ]; names(th) <- par_names
    r   <- lpf(th)
    list(loglik = r$loglik, logprior = r$logprior)
  }
  ## Sever the task env so mirai_map does NOT serialise run_smc_mirai's frame
  ## (parsed_model / Y, etc.) with every one of the n_particles tasks. Bind the
  ## data the task needs into a fresh env parented on the dynhr namespace (data
  ## bindings shadow any same-named dynhr/base functions); .worker_lp is fetched
  ## from the daemon globalenv inside the task.
  .ev_env <- new.env(parent = asNamespace("dynhr"))
  list2env(list(theta_mat = theta_mat, par_names = par_names), envir = .ev_env)
  environment(eval_task) <- .ev_env
  ev <- mirai::mirai_map(seq_len(n_particles), eval_task)[]
  log_liks <- vapply(ev, function(z) z$loglik %||% ll_floor, numeric(1))
  log_pris <- vapply(ev, function(z) z$logprior %||% -Inf, numeric(1))
  log_liks[!is.finite(log_liks)] <- ll_floor
  n_eval <- n_eval + n_particles
  n_valid <- sum(log_liks > ll_floor)
  if (verbose)
    cat(sprintf("  SMC: %d/%d particles have finite likelihood\n",
                n_valid, n_particles))
  ## Fail loud rather than tempering from an all-infeasible population (lambda
  ## would jump to 1 -> meaningless posterior + log_mlik ~ -1e308). n_valid==0
  ## also catches the interface footgun: SMC tempers the LIKELIHOOD and reads
  ## `$loglik`, so a $logpost-only closure floors every particle to ll_floor.
  if (n_valid == 0L) {
    stop(sprintf("SMC: 0/%d initial particles have a finite likelihood -- ",
                 n_particles),
         "cannot temper from an all-infeasible population. Either the ",
         "prior draws are outside the feasible region, or `log_post_fn` does ",
         "not return a finite `$loglik` (SMC tempers prior * lik^lambda and ",
         "needs $loglik/$logprior, NOT a $logpost-only closure).", call. = FALSE)
  }

  ## ---- Tempering loop ------------------------------------------------------
  lambda_curr <- 0; log_marginal <- 0; stage <- 0L
  lambda_trace <- ess_trace <- accept_trace <- scale_trace <- numeric(0)
  log_w <- rep(0, n_particles)
  mh_c  <- mh_scale_factor / sqrt(d)

  while (lambda_curr < 1) {
    stage <- stage + 1L
    lambda_next <- if (!is.null(lambda_schedule) && stage <= length(lambda_schedule))
      lambda_schedule[stage]
    else .smc_next_lambda(log_liks, lambda_curr, ess_target, n_particles)
    lambda_next <- min(lambda_next, 1)
    dlambda <- lambda_next - lambda_curr

    inc_log_w <- dlambda * log_liks

    # p(Y|M) contribution = sum_i W_prev_i * exp(inc_log_w_i), where W_prev are
    # the NORMALISED weights coming INTO this stage (reduces to the old
    # logsumexp(inc_log_w) - log(N) when log_w is uniform, i.e. after a
    # resample; see sampler-smc.R for the non-tempered derivation).
    log_w_prev_norm <- log_w - .smc_log_sum_exp(log_w)
    log_marginal <- log_marginal + .smc_log_sum_exp(log_w_prev_norm + inc_log_w)

    log_w <- log_w + inc_log_w

    w_norm <- exp(log_w - max(log_w)); w_norm <- w_norm / sum(w_norm)
    ess <- 1 / sum(w_norm^2)
    lambda_trace <- c(lambda_trace, lambda_next)
    ess_trace    <- c(ess_trace, ess)

    if (ess < ess_target * n_particles) {
      idx <- .smc_systematic_resample(w_norm, n_particles)
      theta_mat <- theta_mat[idx, , drop = FALSE]
      log_liks  <- log_liks[idx]; log_pris <- log_pris[idx]
      log_w     <- rep(0, n_particles)
      w_norm    <- rep(1 / n_particles, n_particles)
    }

    # Weighted proposal covariance (uniform weights if a resample just
    # occurred, non-uniform otherwise).
    Sigma_hat <- stats::cov.wt(theta_mat, wt = w_norm, method = "ML")$cov + diag(1e-6, d)
    L_prop <- t(chol(Sigma_hat))
    scale  <- mh_c
    stage_seed <- seed_base + stage * n_particles

    ## Herbst-Schorfheide (2014 JAE; 2015 book ch. 5) / FRBNY SMC.jl 3-component
    ## mixture proposal machinery. Only built when mixture_weights is non-NULL,
    ## so the default path neither computes these nor draws any extra RNG (the
    ## per-particle RNG stream below is therefore unchanged when NULL).
    theta_bar <- chol_S <- chol_Dg <- NULL
    if (use_mixture) {
      theta_bar <- as.numeric(w_norm %*% theta_mat)
      chol_S    <- scale * L_prop  # == t(chol(scale^2 * Sigma_hat))
      Dg        <- diag(diag(Sigma_hat), d)
      chol_Dg   <- scale * t(chol(Dg))
    }

    ## All names referenced below (theta_mat, ll_vec, ...) are FREE variables in
    ## mutate_task and therefore passed through `...` of mirai_map.
    mutate_task <- function(i) {
      lpf <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      set.seed(stage_seed + i)
      theta_i <- theta_mat[i, ]; names(theta_i) <- par_names
      ll_i <- ll_vec[i]; lp_i <- lp_vec[i]
      tlp_i <- lp_i + lambda_next * ll_i
      acc <- 0L
      for (s in seq_len(n_mh_steps)) {
        if (use_mixture) {
          # Choose mixture component, then propose accordingly.
          u <- runif(1)
          z <- rnorm(d)
          if (u < mixture_weights[1]) {
            theta_prop <- theta_i + as.numeric(chol_S %*% z)
          } else if (u < mixture_weights[1] + mixture_weights[2]) {
            theta_prop <- theta_i + as.numeric(chol_Dg %*% z)
          } else {
            theta_prop <- theta_bar + as.numeric(chol_S %*% z)
          }
        } else {
          z <- rnorm(d)
          theta_prop <- theta_i + scale * as.numeric(L_prop %*% z)
        }
        names(theta_prop) <- par_names
        rp <- lpf(theta_prop)
        ll_p <- rp$loglik; lp_p <- rp$logprior
        if (!is.finite(ll_p)) ll_p <- -1e300
        if (!is.finite(lp_p)) lp_p <- -Inf
        tlp_p <- lp_p + lambda_next * ll_p
        log_alpha <- tlp_p - tlp_i
        if (use_mixture) {
          # Asymmetric proposal (component 3 is independence-style): correct
          # the MH ratio with log q(prop -> curr) - log q(curr -> prop).
          log_q_fwd <- dynhr:::.smc_mixture_logq(theta_i, theta_prop, theta_bar,
                                                  chol_S, chol_Dg, mixture_weights)
          log_q_rev <- dynhr:::.smc_mixture_logq(theta_prop, theta_i, theta_bar,
                                                  chol_S, chol_Dg, mixture_weights)
          log_alpha <- log_alpha + (log_q_rev - log_q_fwd)
        }
        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
          theta_i <- theta_prop; ll_i <- ll_p; lp_i <- lp_p
          tlp_i <- tlp_p; acc <- acc + 1L
        }
      }
      list(theta = theta_i, loglik = ll_i, logprior = lp_i, accepted = acc)
    }
    ## Sever env via a data-bound child of the dynhr namespace (see eval_task);
    ## note `scale` collides with base::scale, so the data binding is essential.
    .mut_env <- new.env(parent = asNamespace("dynhr"))
    list2env(list(theta_mat = theta_mat, ll_vec = log_liks, lp_vec = log_pris,
                  lambda_next = lambda_next, scale = scale, L_prop = L_prop,
                  n_mh_steps = n_mh_steps, d = d, par_names = par_names,
                  stage_seed = stage_seed, use_mixture = use_mixture,
                  mixture_weights = mixture_weights, theta_bar = theta_bar,
                  chol_S = chol_S, chol_Dg = chol_Dg), envir = .mut_env)
    environment(mutate_task) <- .mut_env
    mut <- mirai::mirai_map(seq_len(n_particles), mutate_task)[]

    n_acc <- 0L
    for (i in seq_len(n_particles)) {
      theta_mat[i, ] <- mut[[i]]$theta
      log_liks[i]    <- mut[[i]]$loglik
      log_pris[i]    <- mut[[i]]$logprior
      n_acc <- n_acc + mut[[i]]$accepted
    }
    n_eval <- n_eval + n_particles * n_mh_steps
    accept_rate <- n_acc / (n_particles * n_mh_steps)
    accept_trace <- c(accept_trace, accept_rate)
    scale_trace  <- c(scale_trace, mh_c)
    mh_c <- mh_c * .smc_scale_adjust(accept_rate, target = mut_target)

    if (verbose)
      cat(sprintf("  SMC stage %d: lambda=%.4f ESS=%.0f accept=%.0f%% log_mlik=%.2f\n",
                  stage, lambda_next, ess, 100 * accept_rate, log_marginal))
    lambda_curr <- lambda_next
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  colnames(theta_mat) <- par_names
  logpost_final <- log_liks + log_pris
  ## Final normalised weights: w_norm is still in scope from the tempering
  ## loop (uniform if the last stage resampled; non-uniform otherwise).
  ## Stored so as_posterior_draws() can detect and correct for the
  ## non-uniform case.
  smc_weights_final <- w_norm
  list(
    chain            = theta_mat,
    particles        = theta_mat,
    smc_weights      = smc_weights_final,
    log_liks         = log_liks,
    log_priors       = log_pris,
    logpost_trace    = logpost_final,
    post_logpost     = logpost_final,
    log_marginal_lik = log_marginal,
    lambda_schedule  = lambda_trace,
    ess_schedule     = ess_trace,
    accept_schedule  = accept_trace,
    scale_schedule   = scale_trace,
    n_stages         = stage,
    n_particles      = n_particles,
    n_draws          = n_particles,
    n_burn           = 0L,
    n_eval           = n_eval,
    elapsed_secs     = elapsed,
    sampler          = "smc"
  )
}
