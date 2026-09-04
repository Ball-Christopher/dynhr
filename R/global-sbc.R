## R/global-sbc.R
## --------------------------------------------------------------------------
## F1-C: simulation-based calibration (SBC; Talts et al. 2018) for the
## GLOBAL-solution bootstrap particle filter shipped in E1-B
## (make_log_posterior_global_pf(); `likelihood = "global_pf"`).
##
## WHY A DEDICATED HARNESS.  dynhr_sbc()'s data-generating process is the
## order-1/2 PERTURBATION state space of a user-supplied model.  Running the
## global PF against that DGP would certify nothing about the filter: the
## generative law and the filtered law would be two DIFFERENT models, so any
## rank non-uniformity would be model mismatch, not miscalibration.  The
## global PF's law is the PROJECTION policy, so its SBC needs a projection
## DGP -- exactly the shape tpf_order3_sbc() (R/tpf-sbc.R) has for the
## order-3 TPF, and the reason no `"global_pf"` branch was added to
## dynhr_sbc(): that entry point takes the caller's model, while this harness
## owns its fixture (the projection solve is the expensive, model-class-
## restricted part).
##
## THE TIMING AUDIT (F1-C task 1; the finding is recorded here because the
## harness is what depends on it).  Two recursions existed in the package for
## advancing a GlobalSolution:
##
##   filter (.global_pf_loglik)   feed_t = s_{t-1} + psi*eps_t/rho
##                                y_t    = policy(feed_t)
##                                s_t    = y_t[state_names]
##
##   simulate.GlobalSolution      y_t       = policy(state_lag)
##   (BEFORE this wave)           state_lag = compute_next_lag(y_t, eps_t)
##
## They are the same LAW one period apart: simulate()'s row t + 1 equalled
## the filter recursion's row t, so simulate()'s FIRST row was deterministic
## (the shock-free policy value at the initial lag) and its LAST drawn shock
## never entered the output.  A user simulating T periods therefore got a
## sample whose first observation had zero shock variance -- and an SBC built
## on it would have been mis-specified relative to the likelihood in exactly
## the first period.  simulate.GlobalSolution() is FIXED in this wave to the
## filter's timing (shock applied in the CURRENT period); see the regression
## in test-global-sbc.R, which pins the two recursions against each other on
## a shared shock path.
##
## Both sides of this harness therefore go through `.gpf_feed_lag()` -- the
## FILTER's own feed helper -- so the DGP cannot drift away from the
## likelihood's timing.
##
## INITIALISATION.  The filter starts its particle cloud at
## N(ss[state_names], P0) with P0 the order-1 stationary state covariance
## (`kf_stationary_init`, reordered to the projection's state ordering).  The
## DGP draws its initial state lag from the SAME normal.  Starting the DGP at
## the steady state instead (the natural-looking choice, and what
## simulate.GlobalSolution() still does by default) would make the first
## observations mis-specified relative to the likelihood and bias the ranks.
##
## STATE DOMAIN (F1-C task 4).  predict.GlobalSolution() CLIPS its input to
## `state_domain` (pmax/pmin on the normalised coordinate), silently.  Two
## consequences are handled here:
##   * The default auto domain is +-3.5 STATIONARY sd of each AR(1) state,
##     but the policy is evaluated at the FEED point s_{t-1} + psi*eps_t/rho,
##     whose sd is larger by a factor sqrt(1 + (sqrt(1-rho^2)/rho)^2).  This
##     harness therefore builds its own domain, covering `domain_cover`
##     FEED standard deviations at a PRIOR-UPPER shock sd, and passes that
##     SAME explicit domain to the DGP's solve_global() and to the likelihood
##     closure.  Being theta-independent it also keeps the model class fixed
##     across replications.
##   * The DGP path is nevertheless CHECKED against the domain: a replication
##     whose simulated feed points leave the box is aborted LOUDLY with
##     reason "state_domain_clipping" and counted in `n_failed_repl`, never
##     silently clipped.  It is a tripwire, not a filter -- it is expected
##     never to fire, and if it starts firing the domain is too narrow.
## --------------------------------------------------------------------------


## The fixture: the TWO-shock RBC of test-global-likelihood.R section 5
## (technology z, rho = 0.9; discount-factor g, rho = 0.6) with the two shock
## standard deviations promoted to MODEL PARAMETERS and ESTIMATED.  Both
## stderrs are estimated with DISTINCT values (0.02 vs 0.01) -- the CLAUDE.md
## multi-value rule: a harness that estimated one shock sd, or two equal
## ones, would be blind to a recycled per-shock index in either the DGP or
## the filter.
##
## `sd_scale` multiplies both shock sds AND their priors, so `sd_scale = 0.1`
## is a NEAR-LINEAR calibration of the same model -- the regime in which the
## global PF and the univariate Kalman filter agree, and hence the regime the
## harness-checking design oracle (compare_kf = TRUE) runs in.
#' @noRd
.global_pf_sbc_model <- function(sd_scale = 1) {
  sz <- 0.02 * sd_scale
  sg <- 0.01 * sd_scale
  fmt <- function(x) format(x, scientific = FALSE, digits = 12)
  mod_txt <- paste0(
    "var c k z g;\nvarexo eps_z eps_g;\n",
    "parameters beta alpha rho_z rho_g sig_z sig_g;\n",
    "beta = 0.99; alpha = 0.33; rho_z = 0.9; rho_g = 0.6;\n",
    "sig_z = ", fmt(sz), "; sig_g = ", fmt(sg), ";\n",
    "model;\n",
    "  1/c = beta * exp(g) * (1/c(+1)) * alpha * exp(z(+1)) * k^(alpha-1);\n",
    "  c + k = exp(z) * k(-1)^alpha;\n",
    "  z = rho_z * z(-1) + eps_z;\n",
    "  g = rho_g * g(-1) + eps_g;\n",
    "end;\n",
    "initval;\n  z = 0; g = 0; k = (0.33*0.99)^(1/(1-0.33));",
    " c = (1-0.33*0.99)*k^0.33;\nend;\n",
    "shocks; var eps_z; stderr sig_z; var eps_g; stderr sig_g; end;\n",
    "estimated_params;\n",
    "  stderr eps_z, inv_gamma_pdf, ", fmt(sz), ", ", fmt(0.3 * sz), ";\n",
    "  stderr eps_g, inv_gamma_pdf, ", fmt(sg), ", ", fmt(0.3 * sg), ";\n",
    "end;\nvarobs c k;\n")

  model    <- suppressMessages(parse_mod(mod_txt, verbose = FALSE))
  compiled <- compile_model(model, verbose = FALSE, max_order = 1L)
  priors   <- suppressMessages(extract_prior_spec(model, verbose = FALSE))
  list(model = model, compiled = compiled, priors = priors,
       obs_vars = c("c", "k"), sd_scale = sd_scale)
}


## Prior-implied, theta-INDEPENDENT collocation domain.
##
## For each AR(1) state the policy is evaluated at the FEED point
## s_{t-1} + psi*eps/rho, whose stationary sd is
##   sqrt( (psi*sig/sqrt(1-rho^2))^2 + (psi*sig/rho)^2 ),
## i.e. STRICTLY wider than the stationary sd of the state itself, which is
## what `.auto_state_domain()` covered before F2-C.  `sig` is taken at a
## prior-UPPER value (mean + `n_prior_sd` prior sds) so one domain covers
## every theta* the prior can draw.  Since F2-C the arithmetic itself is
## `.auto_state_domain()`'s -- this helper only supplies the wider cover and
## the prior-upper sds.
#' @noRd
.global_sbc_domain <- function(fx, domain_cover = 4, n_prior_sd = 5,
                              n_quad = 5L) {
  model    <- fx$model
  compiled <- fx$compiled
  params   <- model$param_values
  ss       <- solve_steady(compiled, params, verbose = FALSE)
  if (!isTRUE(ss$converged))
    stop(".global_sbc_domain: the fixture's baseline steady state did not ",
         "converge.", call. = FALSE)
  ss_vals <- ss$values
  endo    <- compiled$dynamic$endo_names
  exo     <- compiled$dynamic$exo_names
  sn      <- .global_state_names(compiled, endo)
  pairing <- .global_shock_pairing(compiled, params, ss_vals, sn,
                                   context = ".global_sbc_domain")

  ## Prior-upper shock sd, per estimated stderr; a shock whose sd is NOT
  ## estimated keeps its calibrated value.
  sds <- .get_shock_sds(model, params)
  for (i in seq_len(nrow(fx$priors))) {
    nm <- fx$priors$name[i]
    hi <- fx$priors$mean[i] + n_prior_sd * fx$priors$std[i]
    ## `stderr eps_x` priors are renamed to the model parameter (sig_x);
    ## map back to the shock through the shocks block's stderr expression.
    for (s in exo)
      if (identical(.shock_stderr_param(model, s), nm)) sds[[s]] <- hi
  }

  ## F2-C: the feed rule AND the shock-aware capital-like rule now live in
  ## `.auto_state_domain()` (R/global-solve.R), so this harness cannot drift
  ## away from the default the solver actually ships.  The only differences
  ## that remain are the two arguments below: a wider `domain_cover` and the
  ## prior-UPPER shock sds, which is what makes ONE box valid for every theta*
  ## the prior can draw.
  dr1 <- tryCatch(
    solve_perturbation(model, compiled, ss_vals, params, verbose = FALSE),
    error = function(e) NULL)
  dom <- .auto_state_domain(
    state_names  = sn,
    ss_vals      = ss_vals,
    shock_sds    = sds,
    model        = model,
    pairing      = pairing,
    exo          = exo,
    domain_cover = domain_cover,
    ## `domain_cover/3.5` is the inflation this harness always applied to the
    ## capital-like branch; keeping it on the CAP means that where the cap
    ## binds (the large-shock calibrations) the box is exactly the one that
    ## shipped before F2-C, and the SBC's own DGP does not move.
    ## F3-B moved the capital-like default from 12 stationary sds to 6; the
    ## harness tracks the default rather than pinning the old constant.  On
    ## this fixture the change is INERT: at the prior-UPPER shock sds the
    ## capped half-width binds at either cover, so the SBC's DGP box is
    ## unchanged (which is the point of expressing it as an inflation of
    ## whatever the shipped default is).
    endo_cover     = 6 * domain_cover / 3.5,
    level_cap_frac = 0.4 * domain_cover / 3.5,
    level_cap_min  = 0.1 * domain_cover / 3.5,
    state_sd       = .dr_state_sd(dr1, model, params, shock_sds = sds),
    ## F3-B: the AR(1) branch floors the box at the Gauss-Hermite reach, so
    ## the harness must declare the SAME `n_quad` its solves use.
    n_quad         = n_quad)
  list(state_domain = dom, ss = ss, state_names = sn)
}


## Which model parameter supplies shock `s`'s standard deviation, or NA.
## Read from the shocks block's stderr EXPRESSION (the parse-time text), so
## `stderr sig_z` resolves to "sig_z" and a numeric literal to NA.
#' @noRd
.shock_stderr_param <- function(model, shock) {
  v <- model$shocks$variances
  if (is.null(v) || !nrow(v)) return(NA_character_)
  i <- match(shock, v$name)
  if (is.na(i)) return(NA_character_)
  ex <- v$stderr_expr[i]
  if (is.na(ex) || !nzchar(ex)) return(NA_character_)
  if (ex %in% names(model$param_values)) ex else NA_character_
}


## Simulate one panel under the FILTER's law.
##
## Deliberately built on `.gpf_feed_lag()` and `predict.GlobalSolution()` --
## the very functions `.global_pf_loglik()` uses -- so the DGP is the
## likelihood's own transition, not a re-derivation of it.  Returns the full
## endogenous path plus the number of feed COORDINATES that fell outside
## `state_domain` (each of which predict() would silently clip).
##
## @param eps Optional n_T x n_exo matrix of realised shocks (already scaled
##   by the shock covariance).  NULL draws them from N(0, tcrossprod(Le)).
#' @noRd
.global_sbc_simulate <- function(g, s0, Le, n_T, eps = NULL) {
  sn      <- g$state_names
  n_state <- length(sn)
  n_exo   <- length(g$shock_names)
  if (is.null(eps))
    eps <- matrix(stats::rnorm(n_T * n_exo), n_T, n_exo) %*% t(Le)

  lag <- matrix(s0, 1L, n_state, dimnames = list(NULL, sn))
  out <- matrix(NA_real_, n_T, length(g$all_endo_names),
                dimnames = list(NULL, g$all_endo_names))
  n_clip <- 0L
  for (t in seq_len(n_T)) {
    feed <- .gpf_feed_lag(g, lag, eps[t, , drop = FALSE])
    for (j in seq_len(n_state)) {
      d <- g$state_domain[[sn[j]]]
      if (feed[1L, j] < d[1L] || feed[1L, j] > d[2L]) n_clip <- n_clip + 1L
    }
    y <- predict(g, feed)
    if (!all(is.finite(y))) return(NULL)
    out[t, ] <- y
    lag <- y[, sn, drop = FALSE]
  }
  list(path = out, n_clip = n_clip, n_feed = n_T * n_state)
}


## Classify raw per-replication results (mirrors the discipline of
## `.tpf_sbc_classify_reps()`, with this harness's own message; a killed
## mclapply child surfaces as a bare NULL and must NOT be counted a success).
#' @noRd
.global_sbc_classify_reps <- function(reps, n_repl) {
  is_null <- vapply(reps, is.null, logical(1))
  if (any(is_null))
    reps[is_null] <- lapply(reps[is_null], function(x)
      list(failed = TRUE, reason = "null_result_killed_worker"))
  is_err <- vapply(reps, function(x) inherits(x, "try-error"), logical(1))
  reps[is_err] <- lapply(reps[is_err], function(x)
    list(failed = TRUE, reason = "mclapply_error"))
  failed <- vapply(reps, function(x) isTRUE(x$failed), logical(1))
  if (any(failed))
    warning(sprintf("global_pf_sbc: %d/%d replications aborted at theta* ",
                    sum(failed), n_repl),
            "and were excluded from the ranks (see $failure_reasons).",
            call. = FALSE)
  list(reps = reps, failed = failed,
       failure_reasons = vapply(reps[failed],
                                function(x) x$reason %||% "unknown",
                                character(1)))
}


## One adaptive-covariance (Haario-style) RWMH chain over `lp_fn`, returning
## the rank of `theta_star` among the thinned draws.  Shared by the
## global-PF chain and the Kalman-filter comparison chain so the design
## oracle compares LIKELIHOODS, not two different samplers.
#' @noRd
.global_sbc_chain <- function(lp_fn, theta_star, prior_mean, prior_sd,
                              n_warmup, n_draws, thin_L) {
  n_par  <- length(theta_star)
  cur_th <- theta_star
  cur_lp <- lp_fn(cur_th)
  if (!is.finite(cur_lp)) {
    cur_th <- prior_mean
    cur_lp <- lp_fn(cur_th)
  }
  if (!is.finite(cur_lp)) return(NULL)

  log_scale  <- log(0.5)
  L_prop     <- diag(prior_sd, n_par)
  warm_store <- matrix(NA_real_, n_warmup, n_par)
  n_tot      <- n_warmup + n_draws
  keep_every <- max(1L, floor(n_draws / thin_L))
  kept <- matrix(NA_real_, nrow = ceiling(n_draws / keep_every), ncol = n_par)
  k_i <- 0L; n_acc <- 0L
  for (it in seq_len(n_tot)) {
    prop <- cur_th + exp(log_scale) *
      as.numeric(L_prop %*% stats::rnorm(n_par))
    prop_lp <- lp_fn(prop)
    ## Pseudo-marginal accept: the CURRENT point's retained estimate is
    ## never re-evaluated (re-evaluating it destroys the exactness).
    if (is.finite(prop_lp) && log(stats::runif(1)) < (prop_lp - cur_lp)) {
      cur_th <- prop; cur_lp <- prop_lp
      if (it > n_warmup) n_acc <- n_acc + 1L
      acc <- 1
    } else acc <- 0
    if (it <= n_warmup) {
      warm_store[it, ] <- cur_th
      log_scale <- log_scale + (acc - 0.25) / sqrt(it)
      if (it == floor(n_warmup / 2) || it == n_warmup) {
        cv <- stats::cov(warm_store[seq_len(it), , drop = FALSE])
        cv <- cv + diag(1e-12 + 1e-4 * prior_sd^2, n_par)
        Lc <- tryCatch(t(chol(cv)), error = function(e) NULL)
        if (!is.null(Lc)) L_prop <- Lc * (2.38 / sqrt(n_par)) / exp(log_scale)
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
  list(ranks = vapply(seq_len(n_par),
                      function(j) sum(kept[, j] < theta_star[j]), numeric(1)),
       accept = n_acc / n_draws, L = k_i)
}


#' Simulation-based calibration of the global-PF posterior
#'
#' Rank-uniformity SBC (Talts et al. 2018) for
#' \code{\link{make_log_posterior_global_pf}} -- the bootstrap particle
#' filter run on a \code{\link{solve_global}} projection policy
#' (\code{likelihood = "global_pf"}).  Per replication:
#' \code{theta*} is drawn from the prior, the model is re-solved GLOBALLY at
#' \code{theta*}, observables are simulated from that projection policy
#' \emph{under the filter's own timing and initialisation}, i.i.d. Gaussian
#' measurement error of a fixed, known variance is added, a pseudo-marginal
#' RWMH chain is run over the \code{global_pf} closure (\code{seed = NULL},
#' so every evaluation draws a fresh particle cloud -- a seeded closure
#' would break pseudo-marginal validity), and the rank of \code{theta*}
#' among the thinned draws is recorded.  Uniform ranks certify the whole
#' stack; \code{sbc_uniformity_test()} supplies the verdict.
#'
#' The fixture is the two-shock RBC of \code{test-global-likelihood.R}
#' (technology \code{z} with \code{rho = 0.9} and a discount-factor shock
#' \code{g} with \code{rho = 0.6}) with BOTH shock standard deviations
#' estimated, at DISTINCT values (0.02 and 0.01) -- the multi-value rule of
#' \code{CLAUDE.md}, without which a recycled per-shock index in either the
#' DGP or the filter would be invisible.  \code{c} and \code{k} are
#' observed.
#'
#' \strong{Timing.}  The DGP advances with \code{.gpf_feed_lag()}, the
#' FILTER's own AR(1) feed helper: \code{y_t = policy(s_\{t-1\} + psi
#' eps_t / rho)}, \code{s_t = y_t[state_names]}.  Its initial state lag is
#' drawn from \code{N(ss, P0)} with \code{P0} the order-1 stationary state
#' covariance the filter itself initialises the particle cloud with, NOT
#' from the steady state (which would mis-specify the early observations).
#'
#' \strong{State domain.}  \code{predict.GlobalSolution()} clips silently to
#' \code{state_domain}, and the default auto domain covers the stationary sd
#' of each AR(1) state rather than the wider sd of the FEED point the policy
#' is actually evaluated at.  This harness builds one theta-independent
#' domain covering \code{domain_cover} feed standard deviations at a
#' prior-upper shock sd and passes it to both the DGP solve and the
#' likelihood closure; a replication whose simulated path still leaves the
#' box is aborted with reason \code{"state_domain_clipping"} and counted in
#' \code{n_failed_repl} rather than silently clipped.
#'
#' \strong{Cost.}  Every likelihood evaluation re-runs \code{solve_global()},
#' which dominates: measured ~0.77 s per evaluation at the shipped defaults
#' (\code{poly_degree = 2}, \code{n_quad = 2}, \code{n_nodes = 3},
#' \code{solve_tol = 1e-5}; three states, two shocks), against ~1 ms for the
#' particle filter itself at \code{T_obs = 60} / \code{n_particles = 600}.
#' The projection settings are therefore the cheapest ones that still
#' converge, and \code{T_obs} / \code{n_particles} are generous because they
#' are nearly free.  A certification chain
#' (\code{n_warmup = 200}, \code{n_draws = 400}) is ~600 evaluations,
#' i.e. ~8 min per replication serially.
#'
#' @param n_repl Number of SBC replications (default 100).
#' @param T_obs Sample length of each simulated panel (default 60).
#' @param n_particles Bootstrap-PF particles per evaluation (default 600;
#'   measured PF-noise sd ~1 nat at the shipped sizing).
#' @param n_draws,n_warmup Post-warmup / warmup pseudo-marginal RWMH draws
#'   (defaults 100 / 50 -- deliberately small, so the always-on machinery
#'   smoke is cheap; the certification passes 400 / 200).
#' @param thin_L Number of thinned draws used for the rank statistic
#'   (default 25; ranks take values \code{0..thin_L}).
#' @param me_var_frac Fraction of the reference observable variance used as
#'   the FIXED, KNOWN measurement-error variance applied identically to
#'   every replication and to every likelihood closure (default 0.05).  It
#'   is derived ONCE from a reference panel at the PRIOR MEAN theta: deriving
#'   it per replication from that replication's own simulated variance would
#'   make the nominally-exogenous \code{me_variance} secretly informative
#'   about the shock sds being estimated, an SBC confound.
#' @param poly_degree,n_quad,n_nodes,solve_tol,solve_max_iter Projection
#'   settings, passed unchanged to \code{\link{solve_global}} in the DGP and
#'   to the likelihood closure.
#' @param domain_cover Feed standard deviations covered by the fixed
#'   collocation domain (default 4).
#' @param sd_scale Multiplies the fixture's shock sds AND their priors
#'   (default 1).  \code{sd_scale = 0.1} is the near-linear regime used by
#'   the \code{compare_kf} design oracle.
#' @param compare_kf Also run, on the SAME replications (same \code{theta*},
#'   same data, same chain settings), a chain over an exact univariate
#'   Kalman-filter posterior, and return its ranks as \code{ranks_kf}.  In
#'   the near-linear regime the two rank vectors must agree in distribution
#'   -- a check of the HARNESS, not of the likelihood.
#' @param seed Base RNG seed; replication \code{r} uses \code{seed + r}.
#' @param cores Parallel replications via \code{parallel::mclapply}
#'   (default 1 = serial), matching \code{\link{tpf_order3_sbc}}.
#' @return A list with \code{ranks} (matrix, one row per successful
#'   replication, columns \code{sig_z}/\code{sig_g}), \code{uniformity}
#'   (from \code{sbc_uniformity_test()}), \code{accept_rates},
#'   \code{loglik_sd_at_truth} (per-replication PF-noise diagnostic, 4
#'   evaluations at \code{theta*}), \code{n_failed_repl},
#'   \code{failure_reasons}, optionally \code{ranks_kf}, and
#'   \code{settings}.
#' @seealso \code{\link{make_log_posterior_global_pf}},
#'   \code{\link{solve_global}}, \code{\link{tpf_order3_sbc}}
#' @export
global_pf_sbc <- function(n_repl = 100L, T_obs = 60L, n_particles = 600L,
                          n_draws = 100L, n_warmup = 50L, thin_L = 25L,
                          me_var_frac = 0.05,
                          poly_degree = 2L, n_quad = 2L, n_nodes = 3L,
                          solve_tol = 1e-5, solve_max_iter = 200L,
                          domain_cover = 4, sd_scale = 1,
                          compare_kf = FALSE,
                          seed = 20260903L, cores = 1L) {
  fx       <- .global_pf_sbc_model(sd_scale)
  model    <- fx$model
  compiled <- fx$compiled
  priors   <- fx$priors
  par_nm   <- priors$name
  n_par    <- length(par_nm)
  prior_sd <- priors$std
  prior_mean <- stats::setNames(priors$mean, par_nm)
  prior_sampler <- .smc_make_prior_sampler(priors)
  obs_vars <- fx$obs_vars

  dm  <- .global_sbc_domain(fx, domain_cover = domain_cover,
                            n_quad = n_quad)
  dom <- dm$state_domain

  ## Solve the whole per-theta stack once; used by the reference panel, by
  ## every replication's DGP, and (structurally) by the likelihood closure.
  solve_at <- function(params) {
    ss <- tryCatch(solve_steady(compiled, params, verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(ss) || !isTRUE(ss$converged))
      return(list(reason = "steady_state_at_theta_star"))
    dr <- tryCatch(solve_perturbation(model, compiled, ss$values, params,
                                      verbose = FALSE),
                   error = function(e) NULL)
    if (is.null(dr) || !isTRUE(dr$bk_satisfied))
      return(list(reason = "bk_solve_at_theta_star"))
    g <- tryCatch(
      solve_global(compiled, ss, params, poly_degree = poly_degree,
                   n_quad = n_quad, n_nodes = n_nodes, state_domain = dom,
                   tol = solve_tol, max_iter = solve_max_iter,
                   verbose = FALSE),
      error = function(e) NULL)
    if (is.null(g) || !isTRUE(g$converged))
      return(list(reason = "projection_solve_at_theta_star"))

    ## Stationary state-lag initialisation, identical to the filter's.
    Sigma_e <- .get_shock_cov(model, g$shock_names, params)
    Le <- tryCatch(t(chol(Sigma_e)), error = function(e) .tpf_psd_sqrt(Sigma_e))
    TT <- dr$ghx[dr$state_idx, , drop = FALSE]
    RR <- dr$ghu[dr$state_idx, , drop = FALSE]
    P0 <- tryCatch(kf_stationary_init(TT, RR, Sigma_e), error = function(e) NULL)
    if (is.null(P0) || !all(is.finite(P0)))
      return(list(reason = "stationary_init_at_theta_star"))
    pos <- match(g$state_names, dr$state_vars)
    P0  <- P0[pos, pos, drop = FALSE]
    L0  <- tryCatch(t(chol(P0)), error = function(e) .tpf_psd_sqrt(P0))
    if (is.null(L0) || !all(is.finite(L0)))
      return(list(reason = "stationary_init_at_theta_star"))
    list(ss = ss, dr = dr, g = g, Le = Le, L0 = L0,
         s0 = as.numeric(g$ss_vals[g$state_names]))
  }

  ## ---- DESIGN-TIME CONSTANT me_variance ---------------------------------
  ## One reference panel at the PRIOR MEAN theta, under a harness-internal
  ## seed restored afterwards, so the constant is identical across
  ## replications, seeds and cores.
  ref <- solve_at(apply_theta_to_params(model, prior_mean))
  if (is.null(ref$g))
    stop("global_pf_sbc: the reference solve at the PRIOR MEAN theta failed (",
         ref$reason, ") -- re-centre the priors or widen the projection ",
         "settings.", call. = FALSE)
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(920260903L)
  sim_ref <- .global_sbc_simulate(ref$g, ref$s0, ref$Le, 400L)
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  else if (exists(".Random.seed", envir = .GlobalEnv))
    rm(".Random.seed", envir = .GlobalEnv)
  if (is.null(sim_ref))
    stop("global_pf_sbc: the reference panel at the prior mean is not finite.",
         call. = FALSE)
  ref_var <- mean(apply(sim_ref$path[, obs_vars, drop = FALSE], 2, stats::var))
  me_var_const <- me_var_frac * ref_var
  if (!is.finite(me_var_const) || me_var_const <= 0)
    stop("global_pf_sbc: computed me_var_const is degenerate (", me_var_const,
         "); check me_var_frac and the reference panel.", call. = FALSE)

  run_one <- function(r) {
    set.seed(seed + r)
    theta_star <- prior_sampler()[par_nm]
    params <- apply_theta_to_params(model, theta_star)
    sol <- solve_at(params)
    if (is.null(sol$g)) return(list(failed = TRUE, reason = sol$reason))

    ## Initial state lag from the SAME N(ss, P0) the filter uses.
    s0 <- sol$s0 + as.numeric(sol$L0 %*% stats::rnorm(length(sol$s0)))
    sim <- .global_sbc_simulate(sol$g, s0, sol$Le, T_obs)
    if (is.null(sim))
      return(list(failed = TRUE, reason = "nonfinite_simulated_path"))
    ## LOUD, counted rejection -- never a silent clip (F1-C task 4).
    if (sim$n_clip > 0L)
      return(list(failed = TRUE, reason = "state_domain_clipping"))

    Yl <- t(sim$path[, obs_vars, drop = FALSE])          # n_obs x T, LEVELS
    Y  <- Yl + matrix(stats::rnorm(length(Yl), 0, sqrt(me_var_const)),
                      nrow(Yl), ncol(Yl))

    lp_pf <- make_log_posterior_global_pf(
      model, Y, priors, obs_vars, compiled,
      me_variance = me_var_const, n_particles = n_particles,
      poly_degree = poly_degree, n_quad = n_quad, n_nodes = n_nodes,
      state_domain = dom, solve_tol = solve_tol,
      solve_max_iter = solve_max_iter, seed = NULL)
    f_pf <- function(th) lp_pf(th)$logpost

    ll_reps <- vapply(1:4, function(k) lp_pf(theta_star)$loglik, numeric(1))

    ch <- .global_sbc_chain(f_pf, theta_star, prior_mean, prior_sd,
                            n_warmup, n_draws, thin_L)
    if (is.null(ch)) return(list(failed = TRUE, reason = "degenerate_chain_start"))

    out <- list(failed = FALSE, ranks = ch$ranks, accept = ch$accept,
                ll_sd = stats::sd(ll_reps), L = ch$L, n_clip = sim$n_clip)

    if (compare_kf) {
      ## HARNESS oracle: an EXACT likelihood on the same data. The
      ## univariate filter is kept as the comparator for continuity; since
      ## F3-D (2026-09-03) the multivariate filter is the same true-ME law
      ## and either would do (see the header of test-global-likelihood.R).
      f_kf <- function(th) {
        lpri <- log_prior(th, priors)
        if (!is.finite(lpri)) return(-Inf)
        pp <- apply_theta_to_params(model, th)
        s2 <- tryCatch(solve_steady(compiled, pp, verbose = FALSE),
                       error = function(e) NULL)
        if (is.null(s2) || !isTRUE(s2$converged)) return(-Inf)
        d2 <- tryCatch(solve_perturbation(model, compiled, s2$values, pp,
                                          verbose = FALSE),
                       error = function(e) NULL)
        if (is.null(d2) || !isTRUE(d2$bk_satisfied)) return(-Inf)
        m2 <- model; m2$param_values <- pp
        kf <- tryCatch(kalman_filter(Y, d2, m2, pp, obs_vars,
                                     me_variance = me_var_const,
                                     lik_init = "stationary",
                                     method = "univariate",
                                     me_floor_check = FALSE),
                       error = function(e) NULL)
        if (is.null(kf) || !is.finite(kf$loglik)) return(-Inf)
        lpri + kf$loglik
      }
      ck <- .global_sbc_chain(f_kf, theta_star, prior_mean, prior_sd,
                              n_warmup, n_draws, thin_L)
      out$ranks_kf <- if (is.null(ck)) rep(NA_real_, n_par) else ck$ranks
    }
    out
  }

  reps <- if (cores > 1L) {
    parallel::mclapply(seq_len(n_repl), run_one, mc.cores = cores,
                       mc.preschedule = FALSE)
  } else lapply(seq_len(n_repl), run_one)

  cls <- .global_sbc_classify_reps(reps, n_repl)
  ok_reps <- cls$reps[!cls$failed]

  ranks_mat <- do.call(rbind, lapply(ok_reps, function(x) x$ranks))
  if (!is.null(ranks_mat)) colnames(ranks_mat) <- par_nm
  L_support <- if (length(ok_reps)) ok_reps[[1L]]$L else NULL

  out <- list(
    ranks = ranks_mat,
    uniformity = sbc_uniformity_test(ranks_mat, L = L_support),
    accept_rates = vapply(ok_reps, function(x) x$accept, numeric(1)),
    loglik_sd_at_truth = vapply(ok_reps, function(x) x$ll_sd, numeric(1)),
    n_failed_repl = sum(cls$failed),
    failure_reasons = cls$failure_reasons,
    settings = list(n_repl = n_repl, T_obs = T_obs,
                    n_particles = n_particles, n_draws = n_draws,
                    n_warmup = n_warmup, thin_L = thin_L,
                    me_var_frac = me_var_frac, me_var_const = me_var_const,
                    poly_degree = poly_degree, n_quad = n_quad,
                    n_nodes = n_nodes, solve_tol = solve_tol,
                    domain_cover = domain_cover, sd_scale = sd_scale,
                    state_domain = dom, seed = seed))
  if (compare_kf) {
    rk <- do.call(rbind, lapply(ok_reps, function(x) x$ranks_kf))
    if (!is.null(rk)) colnames(rk) <- par_nm
    out$ranks_kf <- rk
  }
  out
}
