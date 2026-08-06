## R/tpf-sbc.R
## --------------------------------------------------------------------------
## Simulation-based calibration (SBC; Talts et al. 2018) for the ORDER-3
## Tempered Particle Filter (TPF) likelihood shipped in ec7c155
## (.tpf_propagate3_R / .tpf_log_weights3_R / tpf_run_period3;
## make_log_posterior_tpf(order = 3L)).
##
## Exactly analogous in structure to sv_rbpf_sbc() (R/sv-rbpf-sbc.R): for
## each replication, draw theta* ~ prior, simulate data at theta* on the
## SAME pruned order-3 state space the filter itself uses
## (simulate_model_order3(), re-solved at theta*), sample the posterior with
## a pseudo-marginal RWMH over the likelihood = "tpf" / order = 3L closure,
## and record the rank of theta* among thinned posterior draws. If the
## stack (parser -> order-3 solve -> TPF marginal-likelihood estimate ->
## pseudo-marginal MCMC) is correctly calibrated, the ranks are uniform
## (tested with sbc_uniformity_test()).
##
## Pseudo-marginal correctness notes (both load-bearing, copied verbatim
## from sv_rbpf_sbc's discipline):
##   - the likelihood closure is built with seed = NULL, so every evaluation
##     draws a fresh particle cloud (checked directly against
##     make_log_posterior_tpf's source: it only calls set.seed() when
##     `seed` is non-NULL -- confirmed NOT to cache/fix randomness);
##   - the CURRENT draw's likelihood ESTIMATE is retained across iterations
##     and never re-evaluated at the current point.
##
## Model: a 2-shock RBC model (rbc2shock's equations, TFP shock eps_a with
## fixed stderr sig_a = 0.02, and an AR(1) discount-factor shock eps_b with
## DISTINCT stderr sig_b -- CLAUDE.md multi-value convention) compiled to
## max_order = 3L, genuinely order-3-nonlinear (see test-tpf-order3.R
## .rbc2_o3_fixture, whose hxx/hxxx/hxxu/hxuu/huuu/hxss/huss tensors are all
## confirmed nonzero on this equation set).
##
## Estimated parameters (2, chain kept cheap): rho_a (TFP persistence,
## beta_pdf 0.7/0.1) and sig_b (discount-factor shock stderr, `stderr eps_b`
## -> auto-renamed to the model parameter sig_b, inv_gamma_pdf 0.01/0.003).
## Priors were verified (see roxygen below) to keep the order-3 solve/BK
## condition satisfied on 500/500 prior draws; a residual failure at a
## drawn theta* (steady state or BK/order-3 solve) aborts that replication
## with a recorded reason rather than a silent redraw.
## --------------------------------------------------------------------------


## Internal: build the canonical order-3-nonlinear 2-shock RBC model +
## priors used by tpf_order3_sbc(). Kept as a helper so the model text is
## defined once and is trivially inspectable.
.tpf_o3_sbc_model <- function() {
  mod_txt <- paste(
    "var c y k l w r a b;",
    "varexo eps_a eps_b;",
    "parameters beta alpha delta sigma eta rho_a sig_a rho_b sig_b;",
    "beta = 0.99; alpha = 0.33; delta = 0.025; sigma = 1; eta = 1;",
    "rho_a = 0.7; sig_a = 0.02; rho_b = 0.90; sig_b = 0.01;",
    "model;",
    "  c^(-sigma) = beta * exp(b) * c(+1)^(-sigma) * (1 + r(+1) - delta);",
    "  y = a * k(-1)^alpha * l^(1 - alpha);",
    "  k = (1 - delta) * k(-1) + y - c;",
    "  w = eta * c^sigma * l;",
    "  r = alpha * a * k(-1)^(alpha - 1) * l^(1 - alpha);",
    "  w = (1 - alpha) * a * k(-1)^alpha * l^(-alpha);",
    "  log(a) = rho_a * log(a(-1)) + eps_a;",
    "  b = rho_b * b(-1) + eps_b;",
    "end;",
    "initval;",
    "  a = 1.0; b = 0.0; k = 24.7; l = 0.928; y = 2.78; c = 2.16;",
    "  w = 2.00; r = 0.0351;",
    "end;",
    "shocks;",
    "  var eps_a; stderr sig_a;",
    "  var eps_b; stderr sig_b;",
    "end;",
    "varobs c y;",
    "estimated_params;",
    "  rho_a, beta_pdf, 0.7, 0.1;",
    "  stderr eps_b, inv_gamma_pdf, 0.01, 0.003;",
    "end;",
    sep = "\n")
  mf <- tempfile(fileext = ".mod")
  writeLines(mod_txt, mf)
  model    <- parse_mod(mf, verbose = FALSE)
  compiled <- compile_model(model, verbose = FALSE, max_order = 3L)
  priors   <- extract_prior_spec(model, verbose = FALSE)
  obs_vars <- c("c", "y")

  ## ---- Reference simulated-observable variance (design-time constant) ---
  ## me_variance MUST be a FIXED, KNOWN constant, identical for every
  ## replication and for both the simulator and the filter -- if it were
  ## instead re-derived per replication from that replication's OWN
  ## simulated-data variance (as in test-tpf-order3.R's oracle (c), which is
  ## fine there because it is not an SBC context), the realized variance is
  ## informative about sig_b and leaks into the nominally-exogenous
  ## me_variance, confounding the SBC (the generative model would then
  ## differ from the model the filter assumes -- me_variance held fixed).
  ## Fixed here via a ONE-TIME reference panel at the PRIOR MEAN theta, under
  ## a harness-internal RNG seed local to this computation (restored
  ## afterwards) so the constant is reproducible and identical across
  ## replications, seeds, and cores.
  theta_ref  <- setNames(priors$mean, priors$name)
  params_ref <- apply_theta_to_params(model, theta_ref)
  ss_ref     <- solve_steady(compiled, params_ref, verbose = FALSE)
  dr3_ref    <- solve_perturbation(model, compiled, ss_ref$values, params_ref,
                                   order = 3L, verbose = FALSE)
  if (is.null(dr3_ref) || !isTRUE(dr3_ref$bk_satisfied))
    stop("tpf_order3_sbc: the order-3 solve at the PRIOR MEAN (used only to ",
        "fix the reference measurement-error variance) failed to satisfy ",
        "the BK condition -- widen/re-center the priors.", call. = FALSE)
  m_ref <- model
  m_ref$param_values <- params_ref
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(999080204L)   ## harness-internal only; NOT the replication seed
  sim_ref <- simulate_model_order3(dr3_ref, n_periods = 200L, burn_in = 300L,
                                   model = m_ref, pruning = TRUE)
  if (!is.null(old_seed))
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  else if (exists(".Random.seed", envir = .GlobalEnv))
    rm(".Random.seed", envir = .GlobalEnv)
  ref_var <- mean(apply(sim_ref[, obs_vars, drop = FALSE], 2, var))

  list(model = model, compiled = compiled, priors = priors,
       obs_vars = obs_vars, ref_var = ref_var)
}


## Classify raw per-replication results from run_one() (called via
## parallel::mclapply()/lapply() in tpf_order3_sbc()).  A killed forked
## child (e.g. an OOM kill, or hitting a resource limit under mc.cores)
## surfaces as a bare NULL list entry -- neither a "try-error" condition
## object nor a list(failed = TRUE, ...): `NULL$failed` is NULL and
## `isTRUE(NULL)` is FALSE, so the old inline classification silently
## treated a killed child as a SUCCESSFUL replication.  That NULL then
## flowed into `ok_reps` and crashed
## `vapply(ok_reps, function(x) x$accept, numeric(1))` -- typically well
## after a multi-hour parallel run had otherwise finished (#8b, adversarial
## review). Refactored out of tpf_order3_sbc() so it is unit-testable by
## injecting a NULL into a results list directly.
#' @noRd
.tpf_sbc_classify_reps <- function(reps, n_repl) {
  is_null <- vapply(reps, is.null, logical(1))
  if (any(is_null))
    reps[is_null] <- lapply(reps[is_null], function(x)
      list(failed = TRUE, reason = "null_result_killed_worker"))
  is_err <- vapply(reps, function(x) inherits(x, "try-error"), logical(1))
  reps[is_err] <- lapply(reps[is_err], function(x)
    list(failed = TRUE, reason = "mclapply_error"))
  failed <- vapply(reps, function(x) isTRUE(x$failed), logical(1))
  if (any(failed))
    warning(sprintf("tpf_order3_sbc: %d/%d replications aborted at theta* ",
                    sum(failed), n_repl),
           "and were excluded from the ranks (see $failure_reasons).",
           call. = FALSE)
  failure_reasons <- vapply(reps[failed], function(x) x$reason %||% "unknown",
                            character(1))
  list(reps = reps, failed = failed, failure_reasons = failure_reasons)
}


#' Simulation-based calibration of the order-3 TPF posterior
#'
#' Runs \code{n_repl} SBC replications on a small, genuinely
#' order-3-nonlinear 2-shock RBC model (rbc2shock's equations; see
#' \code{.tpf_o3_sbc_model()}), estimating 2 parameters -- the TFP
#' persistence \code{rho_a} (beta_pdf 0.7/0.1) and the discount-factor
#' shock stderr \code{sig_b} (inv_gamma_pdf 0.01/0.003) -- with the
#' \code{likelihood = "tpf"}, \code{order = 3L} pseudo-marginal posterior.
#'
#' Data are simulated with \code{\link{simulate_model_order3}} (the SAME
#' pruned order-3 recursion the TPF filters), re-solving the model at each
#' drawn theta*, with i.i.d. Gaussian measurement error of variance
#' \code{me_var_const} (a DESIGN-TIME CONSTANT, fixed and known, identical
#' across every replication and passed unchanged to both the simulator and
#' the likelihood closure -- NOT re-derived per replication from that
#' replication's own simulated-data variance, which would make the
#' nominally-exogenous me_variance secretly informative about \code{sig_b}
#' and confound the SBC) added on top. \code{me_var_const} is
#' \code{me_var_frac} times the mean observable variance of a ONE-TIME
#' reference panel simulated at the PRIOR MEAN theta (see
#' \code{.tpf_o3_sbc_model()}). A burn-in (default 200 periods, absorbed
#' into \code{simulate_model_order3}'s
#' own \code{burn_in} argument) starts the retained sample from the model's
#' stationary distribution -- a cold start biases SBC exactly as in
#' \code{sv_rbpf_sbc}'s \code{simulate_panel}.
#'
#' Priors were empirically verified (500/500 draws) to keep the per-theta
#' steady-state + order-3/BK solve succeeding; see the file header of
#' \code{R/tpf-sbc.R}. Any RESIDUAL failure at a replication's theta* (rare)
#' aborts that replication with a recorded reason (counted in
#' \code{n_failed_repl}); it is never silently redrawn.
#'
#' Runtime, measured at implementation time (this machine, pure-R order-3
#' TPF kernel, no compiled backend for order = 3L): ONE likelihood
#' evaluation at \code{T_obs = 30}, \code{n_particles = 200} (the brief's
#' reference sizing), including the per-theta order-3 re-solve, took
#' ~2.4-3.0s. At that cost \code{n_draws = 2000}/\code{n_warmup = 600}
#' (the sv_rbpf_sbc-style starting point) would take ~1.5-2 HOURS per
#' replication, far outside budget. A first shrink to \code{T_obs = 20},
#' \code{n_particles = 100}, \code{n_mh = 0} cut the per-eval cost to
#' ~0.3-0.75s, but a 12-replication pilot at those settings showed the TPF's
#' adaptive tempering estimator was FAR too noisy for the pseudo-marginal
#' chain to mix (\code{loglik_sd_at_truth} up to the hundreds of nats,
#' acceptance collapsing to ~0 in most replications) -- \code{n_mh = 0}
#' (no per-stage particle mutation) was the culprit, not particle count:
#' switching to \code{n_mh = 1} cut the noise sd by ~5x at essentially the
#' SAME per-eval cost (mutation is cheap; resampling degeneracy without it
#' is not). A tight \code{me_var_frac} (the tempering instrument) further
#' inflated the noise; \code{me_var_frac = 0.1} tamed it to single digits
#' in most replications. The shipped defaults (\code{T_obs = 20},
#' \code{n_particles = 150}, \code{n_mh = 1}, \code{me_var_frac = 0.1},
#' ~0.4-0.55s per evaluation) with a short chain (\code{n_draws = 100},
#' \code{n_warmup = 50}, \code{thin_L = 25}) measured ~55-90s per
#' replication serially (~154 likelihood evaluations including the 4
#' PF-noise diagnostic draws) -- see the brief
#' \code{brief-tpf-order3-sbc-2026-08-04.md} for the full sizing derivation
#' and the deviation from its 2000/600 starting point.
#'
#' KNOWN SCIENTIFIC RISK (documented, not fixed here): the TPF's ADAPTIVE
#' tempering schedule makes the likelihood estimator only approximately
#' unbiased (the number and placement of phi-stages depends on the data and
#' particle cloud at each evaluation). If the FULL SBC run comes back
#' miscalibrated, the first fallback is a FIXED tempering schedule and/or
#' more particles -- that judgement call is deliberately left to the
#' orchestrator; this harness only provides the machinery.
#'
#' @param n_repl Number of SBC replications (default 100).
#' @param T_obs Sample length of each simulated panel (default 20).
#' @param n_particles TPF particles per likelihood evaluation (default 150).
#' @param n_draws,n_warmup Post-warmup / warmup pseudo-marginal RWMH draws
#'   (defaults 100 / 50; Haario-style adaptive-covariance warmup + a global
#'   Robbins-Monro scale toward 25\% acceptance, both frozen for sampling).
#' @param thin_L Number of (approximately independent) thinned draws used
#'   for the rank statistic (default 25; ranks take values 0..thin_L).
#' @param n_mh TPF per-stage RWMH mutation steps (default 1L -- NOT 0:
#'   measured to cut the PF-noise sd by ~5x at essentially the same
#'   per-eval cost, see the file header; passed straight to
#'   \code{make_log_posterior_tpf}).
#' @param me_var_frac Fraction of the (design-time, prior-mean REFERENCE
#'   panel's) observable variance used as the fixed, known
#'   \code{me_var_const} applied identically to every replication (default
#'   0.1 -- a smaller fraction, e.g. 0.02, was measured to leave the TPF's
#'   adaptive tempering badly noisy at this cheap sizing, with
#'   per-replication loglik_sd_at_truth up to the hundreds of nats and
#'   pseudo-marginal acceptance collapsing to ~0; see the file header). NOT
#'   re-derived per replication -- see the note above.
#' @param burn_in Burn-in periods absorbed into \code{simulate_model_order3}
#'   before the retained \code{T_obs}-period sample (default 200).
#' @param burn_in_init Filter-side particle burn-in passed to
#'   \code{make_log_posterior_tpf} (default 50). The generative side starts
#'   from the full pruned stationary joint (via \code{burn_in}), so the
#'   filter must too: with the default Lyapunov-only init (\code{x2_0 = 0})
#'   the mis-stated initial level is absorbed as extra persistence — the
#'   2026-08-04 R = 100 certification run showed a rho_a-specific rank
#'   shift (mean_rank_z -3.0, sig_b clean) from exactly this, after the
#'   obs-tensor fix removed the dominant bias.
#' @param seed Base RNG seed; replication r uses \code{seed + r}.
#' @param cores Parallel replications via \code{parallel::mclapply}
#'   (default 1 = serial).
#' @return A list: \code{ranks} (matrix, rows = successful replications,
#'   columns named \code{rho_a}/\code{sig_b}), \code{uniformity} (from
#'   \code{sbc_uniformity_test()}), \code{accept_rates},
#'   \code{loglik_sd_at_truth} (per-replication PF noise diagnostic, K = 4
#'   evaluations at theta*), \code{n_failed_repl} (replications aborted at
#'   theta* -- steady-state or order-3/BK solve failure -- excluded from
#'   \code{ranks}), \code{failure_reasons}, and \code{settings}.
#' @seealso \code{\link{make_log_posterior_tpf}}, \code{\link{sv_rbpf_sbc}},
#'   \code{sbc_uniformity_test()}
#' @export
tpf_order3_sbc <- function(n_repl = 100L, T_obs = 20L, n_particles = 150L,
                           n_draws = 100L, n_warmup = 50L, thin_L = 25L,
                           n_mh = 1L, me_var_frac = 0.1, burn_in = 200L,
                           burn_in_init = 50L,
                           seed = 20260804L, cores = 1L) {
  fx       <- .tpf_o3_sbc_model()
  model    <- fx$model
  compiled <- fx$compiled
  priors   <- fx$priors
  par_nm   <- priors$name
  n_par    <- length(par_nm)
  prior_sd <- priors$std
  prior_sampler <- .smc_make_prior_sampler(priors)
  obs_vars <- fx$obs_vars

  ## DESIGN-TIME CONSTANT: the SAME me_variance is used to generate every
  ## replication's data AND to build every replication's likelihood closure
  ## (see .tpf_o3_sbc_model()'s reference-panel derivation). NOT re-derived
  ## from each replication's own simulated variance -- that would make the
  ## nominally-exogenous, nominally-known me_variance secretly informative
  ## about sig_b, an SBC confound.
  me_var_const <- me_var_frac * fx$ref_var
  if (!is.finite(me_var_const) || me_var_const <= 0)
    stop("tpf_order3_sbc: computed me_var_const is degenerate (", me_var_const,
        "); check me_var_frac and the reference panel.", call. = FALSE)

  run_one <- function(r) {
    set.seed(seed + r)
    ## 1. theta* ~ prior.
    theta_star <- prior_sampler()[par_nm]

    ## 2. Re-solve the model at theta* (steady state + order-3 perturbation)
    ##    -- the SAME solve the likelihood closure performs internally.  A
    ##    residual failure here (rare; priors were verified 500/500 at
    ##    design time) aborts the replication with a recorded reason rather
    ##    than a silent redraw.
    params <- apply_theta_to_params(model, theta_star)
    ss <- tryCatch(solve_steady(compiled, params, verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(ss) || !isTRUE(ss$converged))
      return(list(failed = TRUE, reason = "steady_state_at_theta_star"))
    dr3 <- tryCatch(
      solve_perturbation(model, compiled, ss$values, params, order = 3L,
                        verbose = FALSE),
      error = function(e) NULL)
    if (is.null(dr3) || !isTRUE(dr3$bk_satisfied))
      return(list(failed = TRUE, reason = "order3_bk_solve_at_theta_star"))

    ## 3. Data | theta*, on the SAME pruned order-3 state space the TPF
    ##    filters (simulate_model_order3 absorbs the burn-in internally --
    ##    a cold start would bias SBC, mirroring sv_rbpf_sbc's rationale).
    m_th <- model
    m_th$param_values <- params
    sim <- simulate_model_order3(dr3, n_periods = T_obs, burn_in = burn_in,
                                 model = m_th, pruning = TRUE)
    Yl <- sweep(sim[, obs_vars, drop = FALSE], 2, dr3$ys[obs_vars], "+")
    ## me_var_const is the FIXED, KNOWN measurement-error variance (identical
    ## across all replications -- see the note at its computation above);
    ## NOT re-derived from this replication's own simulated variance.
    Y_noisy <- Yl + matrix(rnorm(length(Yl), 0, sqrt(me_var_const)),
                           nrow(Yl), ncol(Yl))

    ## 4. Pseudo-marginal RWMH on the likelihood = "tpf", order = 3L
    ##    closure (seed = NULL -> fresh particle cloud every evaluation).
    lp_fn <- make_log_posterior_tpf(model, t(Y_noisy), priors, obs_vars,
                                    compiled, me_variance = me_var_const,
                                    n_particles = n_particles, n_mh = n_mh,
                                    seed = NULL, order = 3L,
                                    burn_in_init = burn_in_init)

    ## PF-noise diagnostic at the truth.
    ll_reps <- vapply(1:4, function(k) lp_fn(theta_star)$loglik, numeric(1))
    ll_sd <- stats::sd(ll_reps)

    cur_th <- theta_star
    cur_lp <- lp_fn(cur_th)$logpost
    if (!is.finite(cur_lp)) {
      cur_th <- setNames(priors$mean, par_nm)
      cur_lp <- lp_fn(cur_th)$logpost
    }
    if (!is.finite(cur_lp))
      return(list(failed = TRUE, reason = "degenerate_chain_start"))

    ## Adaptive-covariance (Haario-style) warmup, copied verbatim from
    ## sv_rbpf_sbc's run_one: an isotropic prior-scaled proposal mixes
    ## poorly when the parameters have very different curvature, so the
    ## warmup learns the chain covariance and proposes with its Cholesky;
    ## the global scale keeps adapting toward 25% acceptance. Both cov and
    ## scale are FROZEN after warmup.
    log_scale <- log(0.5)
    L_prop <- diag(prior_sd, n_par)
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
        log_scale <- log_scale + (acc - 0.25) / sqrt(it)
        if (it == floor(n_warmup / 2) || it == n_warmup) {
          cv <- stats::cov(warm_store[seq_len(it), , drop = FALSE])
          cv <- cv + diag(1e-8 + 1e-4 * prior_sd^2, n_par)
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
    list(failed = FALSE, ranks = ranks, accept = n_acc / n_draws,
        ll_sd = ll_sd, L = k_i)
  }

  reps <- if (cores > 1L) {
    parallel::mclapply(seq_len(n_repl), run_one, mc.cores = cores,
                       mc.preschedule = FALSE)
  } else {
    lapply(seq_len(n_repl), run_one)
  }

  cls <- .tpf_sbc_classify_reps(reps, n_repl)
  reps <- cls$reps
  failed <- cls$failed
  failure_reasons <- cls$failure_reasons
  ok_reps <- reps[!failed]

  ranks_mat <- do.call(rbind, lapply(ok_reps, function(x) x$ranks))
  colnames(ranks_mat) <- par_nm
  ## True rank support: every replication's keep_every/kept-draws
  ## computation is deterministic given the (n_draws, thin_L) config, so
  ## every successful replication kept the SAME number of draws (x$L) --
  ## pass that explicitly instead of letting sbc_uniformity_test infer the
  ## support from the observed ranks, which is silently wrong whenever the
  ## top rank never appears by chance (#7 one-liner, adversarial review).
  L_support <- if (length(ok_reps)) ok_reps[[1L]]$L else NULL
  list(ranks = ranks_mat,
       uniformity = sbc_uniformity_test(ranks_mat, L = L_support),
       accept_rates = vapply(ok_reps, function(x) x$accept, numeric(1)),
       loglik_sd_at_truth = vapply(ok_reps, function(x) x$ll_sd, numeric(1)),
       n_failed_repl = sum(failed),
       failure_reasons = failure_reasons,
       settings = list(n_repl = n_repl, T_obs = T_obs,
                       n_particles = n_particles, n_draws = n_draws,
                       n_warmup = n_warmup, thin_L = thin_L, n_mh = n_mh,
                       me_var_frac = me_var_frac, me_var_const = me_var_const,
                       burn_in = burn_in, burn_in_init = burn_in_init,
                       seed = seed))
}
