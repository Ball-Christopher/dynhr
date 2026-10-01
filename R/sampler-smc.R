## R/sampler-smc.R
## --------------------------------------------------------------------------
## Phase-2 split from smc-monolith.R.
##
## dynhr_smc()  -- Sequential Monte Carlo with likelihood tempering.
## Helpers: .smc_log_sum_exp, .smc_ess, .smc_systematic_resample,
##           .smc_next_lambda, .smc_make_prior_sampler
## Diagnostics: smc_summary, smc_plot_diagnostics
## --------------------------------------------------------------------------

## Back-compat alias: the stable log-sum-exp lives once in R/solve-helpers.R
## as .logsumexp().  Kept because tests and inst/scripts call it by name.
.smc_log_sum_exp <- function(x) .logsumexp(x)


#' Effective sample size from log weights
#' @noRd
.smc_ess <- function(log_w) {
  log_w <- log_w - max(log_w)
  w <- exp(log_w)
  w <- w / sum(w)
  1 / sum(w^2)
}


#' Multivariate normal log-density via a Cholesky factor (no extra deps)
#'
#' @param x Numeric vector, the evaluation point.
#' @param mean Numeric vector, the mean.
#' @param chol_S Lower-triangular Cholesky factor L such that L %*% t(L) = Sigma.
#' @return Scalar log-density of N(mean, Sigma) at x.
#' @noRd
.dmvnorm_chol <- function(x, mean, chol_S) {
  d <- length(x)
  z <- forwardsolve(chol_S, x - mean)
  -0.5 * d * log(2 * pi) - sum(log(diag(chol_S))) - 0.5 * sum(z^2)
}


#' Log-density of the Herbst-Schorfheide 3-component mixture proposal
#'
#' q(from -> to) = alpha1 * N(to; from, S) + alpha2 * N(to; from, Dg)
#'                 + alpha3 * N(to; theta_bar, S)
#'
#' where S = c^2 * Sigma_hat (full covariance) and Dg = c^2 * diag(diag(Sigma_hat))
#' (diagonalised covariance). The first two components are symmetric random
#' walks (centred on `from`); the third is an independence-style draw centred
#' on the weighted particle mean `theta_bar`, which is what makes the mixture
#' proposal asymmetric and requires the MH correction below
#' (Herbst & Schorfheide 2014 JAE; FRBNY SMC.jl `mvnormal_mixture_density`).
#'
#' @param from Numeric vector, the proposal's conditioning point.
#' @param to Numeric vector, the point at which to evaluate the density.
#' @param theta_bar Numeric vector, weighted particle mean (component 3 centre).
#' @param chol_S Cholesky factor of the full-covariance proposal S.
#' @param chol_Dg Cholesky factor of the diagonalised proposal Dg.
#' @param alphas Length-3 normalised mixture weights (alpha1, alpha2, alpha3).
#' @return Scalar log-density log q(to | from).
#' @noRd
.smc_mixture_logq <- function(from, to, theta_bar, chol_S, chol_Dg, alphas) {
  log_terms <- c(
    log(alphas[1]) + .dmvnorm_chol(to, from,      chol_S),
    log(alphas[2]) + .dmvnorm_chol(to, from,      chol_Dg),
    log(alphas[3]) + .dmvnorm_chol(to, theta_bar, chol_S)
  )
  .smc_log_sum_exp(log_terms)
}


#' Adaptive scaling adjustment for the SMC mutation step (Herbst & Schorfheide
#' 2014, eq. for c_n; FRBNY SMC.jl `update_c!`).
#'
#' Multiplies the proposal scale by a logistic function of the previous stage's
#' mutation acceptance rate, centred on the target rate (default 0.25): the
#' factor lies in [0.95, 1.05], so the scale grows by up to 5% when acceptance
#' is above target (proposal too timid) and shrinks by up to 5% when below
#' (proposal too bold). This is what rescues the mutation under a diffuse prior:
#' a wildly over-wide proposal (acceptance -> 0) is shrunk stage by stage until
#' particles can actually move and concentrate on the posterior.
#'
#' @param accept_rate previous stage's average mutation acceptance rate.
#' @param target target acceptance rate (default 0.25).
#' @return multiplicative adjustment factor in [0.95, 1.05].
#' @noRd
.smc_scale_adjust <- function(accept_rate, target = 0.25) {
  x <- 16 * (accept_rate - target)
  0.95 + 0.10 * (exp(x) / (1 + exp(x)))
}


#' Systematic resampling (O(N), lower variance than multinomial)
#'
#' @param weights Normalized probability weights (sum to 1)
#' @param N Number of samples to draw
#' @return Integer vector of resampled indices
#' @noRd
.smc_systematic_resample <- function(weights, N = length(weights)) {
  cw <- cumsum(weights)
  cw[length(cw)] <- 1  # ensure exact sum

  u <- (seq_len(N) - 1 + runif(1)) / N
  idx <- integer(N)
  j <- 1L
  for (i in seq_len(N)) {
    while (cw[j] < u[i]) j <- j + 1L
    idx[i] <- j
  }
  idx
}


#' Adaptive tempering: find next lambda such that ESS = target
#'
#' Uses bisection to find lambda_next in [lambda_curr, 1] such that the
#' COMBINED weights -- the incoming weights \code{log_w_prev} times the
#' incremental weights -- give ESS = ess_target * N.
#'
#' This used to target the ESS of the INCREMENTAL
#' weights alone. After a stage that did not resample, the incoming weights
#' are already non-uniform, so the combined ESS fell to roughly
#' ess_target^2 * N (0.25 N at the default 0.5) before the resample fired.
#' \code{log_w_prev = NULL} (uniform incoming weights) is the old rule exactly
#' and is what the TPF's internal tempering, which resamples every step,
#' still passes.
#'
#' @param log_liks Vector of per-particle tempering scores (log-likelihoods)
#' @param lambda_curr Current tempering parameter
#' @param ess_target Target ESS ratio (0-1)
#' @param N Number of particles
#' @param log_w_prev Optional unnormalised log weights coming INTO the stage
#'   (\code{NULL} = uniform). If their own ESS is already below the target the
#'   combined target is unreachable and the incremental rule is used instead.
#' @return Next lambda value
#' @noRd
.smc_next_lambda <- function(log_liks, lambda_curr, ess_target, N,
                             log_w_prev = NULL) {
  target_ess <- ess_target * N
  base <- if (is.null(log_w_prev) || .smc_ess(log_w_prev) < target_ess) 0
          else log_w_prev

  # Check if we can go straight to lambda=1
  dlam <- 1 - lambda_curr
  inc_w <- base + dlam * log_liks
  if (.smc_ess(inc_w) >= target_ess) return(1)

  # Bisection

  lo <- lambda_curr
  hi <- 1
  for (iter in 1:50) {
    mid <- (lo + hi) / 2
    dlam <- mid - lambda_curr
    inc_w <- base + dlam * log_liks
    ess <- .smc_ess(inc_w)
    if (ess > target_ess) {
      lo <- mid
    } else {
      hi <- mid
    }
    if (hi - lo < 1e-8) break
  }
  lo
}


#' What SMC tempers, for one log-posterior evaluation
#'
#' A4. SMC draws stage 0 from the PARAMETER prior
#' (\code{.smc_make_prior_sampler()}) with uniform weights, so the only
#' consistent tempered path is
#' \deqn{\pi_\lambda(\theta) \propto p(\theta)\,\exp\{\lambda\,\phi(\theta)\},
#'   \quad \phi = \log s(\theta) + \zeta \log L(\theta)}
#' where \eqn{p} is the parameter prior, \eqn{s} the system prior and
#' \eqn{\zeta} the \code{power_posterior} exponent: stage 0 is then exact and
#' \eqn{\lambda = 1} is the full target every other sampler uses. It used to
#' temper the RAW \code{$loglik} with \code{$logprior} as the base, which (a)
#' ignored \eqn{\zeta} and (b) dropped the system prior ("extra"-mode
#' factories) or sampled stage 0 from the wrong base ("lp" mode, where
#' \code{$logprior} includes it). The components come from
#' \code{.posterior_parts()}; a closure without them keeps the historical
#' contract (\eqn{\phi} = \code{$loglik}).
#'
#' @param res List returned by the log-posterior closure.
#' @param ll_floor Floor replacing a non-finite score (the particle's weight
#'   then vanishes at the first positive tempering increment).
#' @param missing_prior Value for a closure that returns no \code{$logprior}
#'   (stage 0 used 0, a mutation proposal -Inf).
#' @return Named numeric \code{c(logprior, phi, loglik, log_sysprior)}:
#'   \code{logprior} the parameter prior (-Inf off support), \code{phi} the
#'   tempered score and \code{loglik} the raw likelihood (both floored),
#'   \code{log_sysprior} the system prior (0 when absent).
#' @noRd
.smc_particle_parts <- function(res, ll_floor = -1e300, missing_prior = 0) {
  parts <- attr(res, "posterior_parts", exact = TRUE)
  if (is.null(parts)) {
    lpp <- res$logprior
    if (is.null(lpp)) lpp <- missing_prior
    ll  <- res$loglik
    if (is.null(ll)) ll <- -Inf
    sp  <- 0
    phi <- ll
  } else {
    lpp <- parts[["logprior_param"]]
    ll  <- parts[["loglik"]]
    sp  <- parts[["log_sysprior"]]
    ## sp = 0, power = 1 gives phi == ll bit-for-bit (IEEE: 0 + 1 * x == x).
    phi <- sp + parts[["power"]] * ll
  }
  ## unname(): a closure may return named scalars (c(theta = ...)), which c()
  ## would turn into "logprior.theta".
  lpp <- unname(lpp); ll <- unname(ll); phi <- unname(phi); sp <- unname(sp)
  if (!is.finite(ll))  ll  <- ll_floor
  if (!is.finite(phi)) phi <- ll_floor
  if (!is.finite(lpp)) lpp <- -Inf
  c(logprior = lpp, phi = phi, loglik = ll, log_sysprior = sp)
}


#' Model-tempering M0 score: log system prior + log L_M0, floored
#'
#' Vectorised. A floored M0 likelihood, or a non-finite sum (system prior
#' -Inf), gives \code{ll_floor} -- never -Inf, so \code{(1 - 1) * phi0}
#' at the last bridge stage is 0 rather than NaN.
#' @noRd
.smc_phi0 <- function(sp, ll0, ll_floor = -1e300) {
  v <- sp + ll0
  v[!is.finite(v) | ll0 <= ll_floor] <- ll_floor
  v
}


#' Equally-weighted SMC result with every per-particle field re-indexed
#'
#' B5. The estimation drivers resampled \code{$chain}
#' to equal weights but left \code{$post_logpost} (and the other per-particle
#' vectors) in the original particle order, so THAMES paired draw i with the
#' log-posterior of some other particle. The resampling indices are drawn
#' ONCE here (through \code{as_posterior_draws()}, so the ESS warning and the
#' systematic resampler are the shared ones) and applied to every
#' per-particle field. \code{$particles} and \code{$smc_weights} keep the
#' original weighted cloud; \code{$resample_idx} records the map.
#'
#' @param res A \code{dynhr_smc()} / \code{run_smc_mirai()} result.
#' @param seed Optional seed for the resampling RNG.
#' @return \code{res} with \code{chain}, \code{post_logpost},
#'   \code{logpost_trace}, \code{log_liks}, \code{log_priors},
#'   \code{log_phi} and \code{log_sysprior} re-indexed, \code{n_draws} reset
#'   and \code{resample_idx} added.
#' @noRd
.smc_equal_weight_result <- function(res, seed = NULL) {
  src <- res$particles %||% res$chain
  N   <- nrow(src)
  idx <- as.integer(as_posterior_draws(
    list(chain = matrix(seq_len(N), ncol = 1L), smc_weights = res$smc_weights),
    seed = seed)[, 1L])
  res$chain <- src[idx, , drop = FALSE]
  for (f in c("post_logpost", "logpost_trace", "log_liks", "log_priors",
              "log_phi", "log_sysprior"))
    if (length(res[[f]]) == N) res[[f]] <- res[[f]][idx]
  res$n_draws      <- nrow(res$chain)
  res$resample_idx <- idx
  res
}


#' Construct a prior sampler from dynhr prior_spec
#'
#' Every parameter is drawn by `.rprior_dist1()` (R/prior-density.R), the ONE
#' prior sampler, which draws EXACTLY the law `log_prior()` scores: the
#' `.lp_dist1()` density (incl. Dynare's generalised beta on [p3, p4] and the
#' gamma / inverse-gamma shift p3) truncated to `[lower, upper]` by inverse
#' CDF. This used to carry its own switch, which drew a
#' bounded beta as a GENERALISED beta rescaled to [lower, upper] (with a silent
#' `max(v, 2)` shape clamp) while `log_prior()` scored a TRUNCATED standard
#' beta, and dispatched on `tolower()` instead of `.normalize_dist()`, so
#' "beta_pdf" and every other unrecognised name fell through to `rnorm(p1, p2)`.
#'
#' prior_spec can be:
#'   - A data.frame (from extract_prior_spec()) with columns name,
#'     distribution, p1, p2, lower, upper [, p3, p4];
#'   - a named list of lists (or an unnamed one whose elements carry `$name`),
#'     each with `$dist`/`$distribution`, `$p1`/`$mean`, `$p2`/`$sd`,
#'     optional `$p3`/`$p4` and `$lb`/`$lower`, `$ub`/`$upper`.
#'
#' @param prior_spec Prior specification
#' @return function() -> named numeric vector drawn from prior
#' @noRd
.smc_make_prior_sampler <- function(prior_spec) {

  # Normalise to a list of lists
  if (is.data.frame(prior_spec)) {
    spec_list <- lapply(seq_len(nrow(prior_spec)), function(i) as.list(prior_spec[i, ]))
    names(spec_list) <- prior_spec$name
  } else if (is.list(prior_spec) && !is.null(names(prior_spec))) {
    spec_list <- prior_spec
  } else if (is.list(prior_spec) && is.null(names(prior_spec))) {
    # Unnamed list -- try to extract names from $name field
    nms <- vapply(prior_spec, function(x) x$name, character(1))
    spec_list <- prior_spec
    names(spec_list) <- nms
  } else {
    stop("prior_spec format not recognised. Provide a named list or data.frame.")
  }

  num1 <- function(v, default) {
    if (is.null(v) || length(v) == 0L) return(default)
    v <- suppressWarnings(as.numeric(v[[1L]]))
    if (is.na(v)) default else v
  }

  # Build individual samplers
  samplers <- list()
  for (nm in names(spec_list)) {
    sp <- spec_list[[nm]]
    ## Read p1/p2 DIRECTLY (before falling back to mean/sd aliases).  A
    ## prior_spec data.frame from extract_prior_spec() carries both $p1/$p2
    ## (raw parameters) AND $mean/$std (derived summaries).  For UNIFORM
    ## p1 = lower and p2 = upper whereas mean = (p1+p2)/2 and
    ## std = (p2-p1)/sqrt(12) -- reading mean/sd instead of p1/p2 would use the
    ## MIDPOINT as the lower bound (shadow bug) and the std-dev as the upper.
    fn <- local({
      dist <- as.character(sp$dist %||% sp$distribution %||% "normal")
      p1 <- num1(sp$p1 %||% sp$mean %||% 0, NA_real_)
      p2 <- num1(sp$p2 %||% sp$sd %||% 1, NA_real_)
      if (is.na(p1) || is.na(p2))
        .dynhr_abort("prior for \"", nm, "\": p1/p2 must be numbers (got ",
                     p1, ", ", p2, ").", class = "dynhr_error_prior_spec")
      p3 <- num1(sp$p3, NA_real_)
      p4 <- num1(sp$p4, NA_real_)
      lo <- num1(sp$lb %||% sp$lower, -Inf)
      hi <- num1(sp$ub %||% sp$upper,  Inf)
      nm_ <- nm
      ## Validate the name and the shape NOW (classed error at construction,
      ## not mid-run): .prior_law() aborts on an unknown distribution or a
      ## degenerate parameterisation.
      invisible(.prior_law(.normalize_dist(dist), p1, p2, p3, p4, label = nm_))
      function() .rprior_dist1(1L, dist, p1, p2, p3, p4, lo, hi, label = nm_)
    })
    samplers[[nm]] <- fn
  }

  # Return a function that draws all parameters
  function() {
    theta <- vapply(samplers, function(f) f(), numeric(1))
    names(theta) <- names(samplers)
    theta
  }
}


## ---- Argument checks for the SMC family -----------------------------------
## One rule set per sampler, called by the sampler at entry and by
## validate_spec(), so a bad spec fails when it is built and the two cannot
## drift. Each helper takes a named list keyed by the sampler function's own
## argument names (absent = the sampler's default, which is valid) and returns
## the problems as a character vector (character(0) when fine). A message may
## carry a condition class in names(); an empty name means the sampler's
## default class.

## A scalar whole number >= min, or NULL if it is one; else the message.
#' @noRd
.smc_whole_problem <- function(x, nm, min = 1L) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x) &&
      x >= min) return(NULL)
  paste0("`", nm, "` must be a whole number >= ", min, " (got ",
         paste(format(x), collapse = ", "), ").")
}

## A scalar number strictly inside (lo, hi), with the upper end optionally
## closed. NULL if it is one; else the message.
#' @noRd
.smc_unit_problem <- function(x, nm, lo, hi, hi_closed = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x > lo &&
    (if (hi_closed) x <= hi else x < hi)
  if (ok) return(NULL)
  paste0("`", nm, "` must be a single number in (", lo, ", ", hi,
         if (hi_closed) "]" else ")", " (got ",
         paste(format(x), collapse = ", "), ").")
}

## The tempering-ladder rule shared by dynhr_smc(), run_smc_mirai() and
## dynhr_dsmh(): a fixed schedule is finite and strictly increasing in (0, 1].
## `end_at_one` additionally requires the last rung to be 1 (dynhr_dsmh() runs
## exactly the supplied ladder); dynhr_smc() continues adaptively after a
## shorter ladder, so it does not. Returns the message or NULL.
#' @noRd
.smc_schedule_problem <- function(x, nm, end_at_one = FALSE) {
  lam <- if (is.numeric(x)) as.numeric(x) else NA_real_
  if (length(lam) >= 1L && all(is.finite(lam)) && all(lam > 0) &&
      all(lam <= 1) && all(diff(lam) > 0) &&
      (!end_at_one || lam[length(lam)] == 1)) return(NULL)
  paste0("`", nm, "` must be strictly increasing in (0, 1]",
         if (end_at_one) " and end at 1" else "", " (got ",
         paste(format(x), collapse = ", "), ").")
}

#' @noRd
.smc_args_problem <- function(args, n_par = NULL) {
  acc <- new.env(parent = emptyenv())
  acc$p <- character(0)
  add <- function(msg, cls = "") {
    if (!is.null(msg)) acc$p <- c(acc$p, stats::setNames(msg, cls))
    invisible(NULL)
  }
  has <- function(nm) !is.null(args[[nm]])
  if (has("n_particles"))
    add(.smc_whole_problem(args$n_particles, "n_particles", 2L))
  if (has("ess_target"))
    add(.smc_unit_problem(args$ess_target, "ess_target", 0, 1, hi_closed = TRUE))
  if (has("n_mh_steps"))
    add(.smc_whole_problem(args$n_mh_steps, "n_mh_steps", 0L))
  if (has("mut_target"))
    add(.smc_unit_problem(args$mut_target, "mut_target", 0, 1))
  if (has("mh_scale_factor")) {
    s <- args$mh_scale_factor
    if (!(is.numeric(s) && length(s) == 1L && is.finite(s) && s > 0))
      add(paste0("`mh_scale_factor` must be a single positive number (got ",
                 paste(format(s), collapse = ", "), ")."))
  }
  if (has("mixture_weights")) {
    w <- args$mixture_weights
    if (!is.numeric(w) || length(w) != 3L || any(!is.finite(w)))
      add(paste0("`mixture_weights` must be NULL or a finite numeric vector ",
                 "of length 3 (got ", paste(format(w), collapse = ", "), ")."))
    else if (any(w < 0) || !any(w > 0))
      add(paste0("`mixture_weights` must be non-negative with at least one ",
                 "positive entry (got ", paste(format(w), collapse = ", "),
                 ")."))
  }

  ## Fixed tempering ladders. Under model tempering (approx_loglik_fn) the
  ## bridge is governed by phi_schedule and a lambda_schedule has no stage to
  ## govern; under likelihood tempering phi_schedule is an alias for
  ## lambda_schedule, so two different ladders are ambiguous.
  model_temp <- has("approx_loglik_fn")
  if (model_temp && has("lambda_schedule"))
    add(paste0("lambda_schedule is not supported with approx_loglik_fn (model ",
               "tempering): the bridge from the approximating model is ",
               "governed by phi_schedule. Pass phi_schedule instead, or drop ",
               "approx_loglik_fn to use lambda_schedule with likelihood ",
               "tempering."), "dynhr_error_inapplicable_argument")
  else if (!model_temp && has("lambda_schedule") && has("phi_schedule") &&
           !identical(as.numeric(args$lambda_schedule),
                      as.numeric(args$phi_schedule)))
    add(paste0("both lambda_schedule and phi_schedule were supplied and they ",
               "differ. Without approx_loglik_fn (likelihood tempering) ",
               "phi_schedule is an alias for lambda_schedule; pass one of ",
               "them."), "dynhr_error_schedule_conflict")
  else {
    if (!model_temp && has("lambda_schedule"))
      add(.smc_schedule_problem(args$lambda_schedule, "lambda_schedule"))
    if (has("phi_schedule"))
      add(.smc_schedule_problem(args$phi_schedule, "phi_schedule"))
  }
  acc$p
}

## n_mh_steps = 0 is legal (a no-mutation run reweights and resamples the
## prior cloud only, which some diagnostics use on purpose) but as a
## posterior sampler the cloud degenerates to copies of a few prior draws,
## so say so instead of returning it silently.
#' @noRd
.smc_warn_no_mutation <- function(fn, n_mh_steps) {
  if (identical(as.numeric(n_mh_steps), 0))
    .dynhr_warn(fn, ": n_mh_steps = 0 runs no mutation step: the particles ",
                "are only reweighted and resampled, so the cloud degenerates ",
                "to copies of a few prior draws. Use n_mh_steps >= 1 to ",
                "sample the posterior.", class = "dynhr_warning_smc_no_mutation")
  invisible(NULL)
}

## Abort with the helper's problems; a classed message sets the class.
#' @noRd
.smc_abort_problems <- function(fn, p, default_class) {
  cls <- names(p)[nzchar(names(p))]
  .dynhr_abort(fn, ": ", paste(p, collapse = " "),
               class = if (length(cls)) cls[1L] else default_class)
}


# ============================================================================
# dynhr_smc() -- Sequential Monte Carlo with likelihood tempering
# ============================================================================

#' @param log_post_fn function(theta) -> list(logpost, loglik, logprior)
#' @param prior_spec Prior specification (for constructing prior sampler)
#' @param prior_sampler Optional: function() -> named numeric vector from prior.
#'   If NULL, constructed from prior_spec via .smc_make_prior_sampler().
#' @param n_particles Number of particles (1000-10000 recommended)
#' @param ess_target Target ESS ratio for adaptive tempering (0.5-0.99)
#' @param n_mh_steps RWMH mutation steps per tempering stage
#' @param mh_scale_factor Proposal scale = mh_scale_factor / sqrt(d)
#' @param lambda_schedule Fixed likelihood-tempering schedule (NULL = adaptive).
#'   Not supported with approx_loglik_fn (error class
#'   dynhr_error_inapplicable_argument): use phi_schedule there.
#' @param approx_loglik_fn Optional: function(theta) -> scalar; log L_{M0}(theta).
#'   When non-NULL, enables Mlikota & Schorfheide (2024) model tempering:
#'   the bridge pi_phi ∝ p(theta) * L_{M0}^{1-phi} * L_{M1}^{phi} is used
#'   instead of the standard likelihood-tempering schedule. The incremental
#'   weight is dphi * (phi_M1 - phi_M0), where phi_M1 = log s + zeta * log L_M1
#'   (see the Target section) and phi_M0 = log s + log L_M0: M0 shares M1's
#'   parameter AND system prior and is not power-tempered. A non-finite
#'   L_M0 is floored like L_M1, and a particle at the floor under either
#'   model is given zero weight. When NULL (default), likelihood tempering.
#' @param phi_schedule Fixed phi schedule (NULL = adaptive, same bisection as
#'   lambda_schedule). With approx_loglik_fn it is the model-tempering bridge
#'   schedule; without it (likelihood tempering) it is an alias for
#'   lambda_schedule (a one-time classed message says so), and supplying both
#'   with different values is a classed error
#'   (dynhr_error_schedule_conflict).
#' @param log_Z_approx Scalar: log Z_{M0}, the log marginal likelihood of the
#'   approximating model M0 (e.g. from a prior SMC run on M0). Default 0.
#'   When model tempering is active, the returned log_marginal_lik equals
#'   log(Z_M1/Z_M0) + log_Z_approx. With log_Z_approx = 0 (default), the
#'   output is log(Z_M1/Z_M0).
#' @param init_particles Optional n_particles x n_par matrix of stage-0
#'   particle positions, used only with model tempering. For an UNBIASED
#'   log(Z_M1/Z_M0) these MUST be drawn from the M0 *posterior*
#'   (pi_0 ∝ p(theta) * L_{M0}(theta)) -- e.g. the particle cloud from a prior
#'   SMC run on M0 (the Mlikota-Schorfheide workflow). The bridge then
#'   telescopes Z(pi_phi)/Z(pi_{phi-1}) to Z_M1/Z_M0 exactly. When NULL
#'   (default) the stage-0 cloud is drawn from the prior: the particle
#'   positions still converge to the M1 posterior (the mutation kernel targets
#'   pi_phi), but log_marginal_lik is a valid MDD ratio ONLY if the M0
#'   posterior coincides with the prior (e.g. approx_loglik_fn constant). The
#'   returned `marginal_valid` flag records which case applies.
#' @param mixture_weights NULL (default) or a numeric length-3 vector of
#'   non-negative mixture weights (alpha1, alpha2, alpha3), normalised
#'   internally, for the Herbst & Schorfheide (2014 JAE; 2015 book ch. 5)
#'   3-component mutation proposal (see also FRBNY's SMC.jl
#'   `mvnormal_mixture_density`):
#'     alpha1: N(theta, c^2 * Sigma_hat)              -- full covariance (as when NULL)
#'     alpha2: N(theta, c^2 * diag(diag(Sigma_hat)))  -- diagonalised covariance
#'     alpha3: N(theta_bar, c^2 * Sigma_hat)          -- independence draw centred
#'             on the weighted particle mean theta_bar
#'   When NULL, mutation is the original single full-covariance random walk
#'   and the mixture machinery is not touched at all (no extra RNG draws,
#'   bit-identical to pre-mixture behaviour under fixed seeds). When non-NULL,
#'   the third component is asymmetric, so the MH acceptance ratio includes
#'   the proposal-density correction log q(y->x) - log q(x->y).
#' @param parallel Use future.apply for parallel evaluation
#' @param verbose Print progress
#' @param progressor progressr callback or NULL
#'
#' @section Target:
#'   Stage 0 draws from the PARAMETER prior p(theta) and the tempered path is
#'   pi_lambda ∝ p(theta) * exp(lambda * phi(theta)) with
#'   phi = log s(theta) + zeta * log L(theta): s is the system prior and zeta
#'   the power_posterior exponent, both read from the closure's
#'   "posterior_parts" (see .posterior_parts()). lambda = 1 is therefore the
#'   same target every other sampler uses, and
#'   log_marginal_lik = log of the integral of p * s * L^zeta over theta.
#'   Without a system prior and with zeta = 1 that is log p(Y | M). With a
#'   system prior it is log p(Y | M) under the joint prior p * s / C_s PLUS
#'   log C_s, where C_s = the integral of p * s is the system prior's
#'   normalising constant (typically unknown -- the system prior is an
#'   unnormalised density; estimate C_s with a run whose likelihood is
#'   constant if a comparable MDD is needed). A closure without the attribute
#'   is tempered on its raw $loglik with $logprior as the base (historical
#'   contract). Returned $log_liks is the RAW likelihood, $log_priors the
#'   parameter prior, $log_phi the tempered score, $log_sysprior the system
#'   prior, and $post_logpost the full target log p + phi (-Inf when
#'   infeasible).
#' @return List with particles, weights, marginal likelihood, diagnostics.
#'   When approx_loglik_fn is non-NULL (model tempering), log_marginal_lik
#'   equals log(Z_M1/Z_M0) + log_Z_approx, VALID only when the stage-0 cloud is
#'   the M0 posterior (supply init_particles; see that argument). tempering_mode
#'   is "model" or "likelihood"; marginal_valid is TRUE when log_marginal_lik is
#'   a valid marginal/ratio (standard tempering, or model tempering with
#'   init_particles), FALSE for model tempering started from the prior.
#' @noRd
dynhr_smc <- function(
    log_post_fn,
    prior_spec        = NULL,
    prior_sampler     = NULL,
    n_particles       = 2000L,
    ess_target        = 0.5,
    n_mh_steps        = 1L,
    mh_scale_factor   = 0.5,
    mut_target        = 0.25,
    lambda_schedule   = NULL,
    approx_loglik_fn  = NULL,
    phi_schedule      = NULL,
    log_Z_approx      = 0,
    init_particles    = NULL,
    mixture_weights   = NULL,
    parallel          = FALSE,
    backend           = "mirai",
    seed_base         = 1L,
    verbose           = TRUE,
    progressor        = NULL
) {
  stopifnot(is.function(log_post_fn))
  if (is.null(prior_sampler) && is.null(prior_spec))
    stop("Must provide either prior_sampler or prior_spec")

  ## Argument rules (one helper, shared with validate_spec()). Normalise the
  ## Herbst-Schorfheide mixture weights afterwards; NULL keeps the original
  ## single-component random-walk mutation exactly as before.
  prob <- .smc_args_problem(list(
    n_particles = n_particles, ess_target = ess_target,
    n_mh_steps = n_mh_steps, mh_scale_factor = mh_scale_factor,
    mut_target = mut_target, lambda_schedule = lambda_schedule,
    phi_schedule = phi_schedule, mixture_weights = mixture_weights,
    approx_loglik_fn = approx_loglik_fn))
  if (length(prob)) .smc_abort_problems("smc", prob, "dynhr_error_smc_args")
  .smc_warn_no_mutation("smc", n_mh_steps)
  if (!is.null(mixture_weights))
    mixture_weights <- mixture_weights / sum(mixture_weights)
  use_mixture <- !is.null(mixture_weights)

  ## Model tempering (Mlikota & Schorfheide 2024): bridge M0 -> M1 via
  ## pi_phi ∝ p(theta) * L_{M0}^{1-phi} * L_{M1}^{phi}.
  ## All new code is guarded by use_model_tempering; the default path
  ## (approx_loglik_fn = NULL) is bit-identical to the previous behaviour.
  use_model_tempering <- !is.null(approx_loglik_fn)

  ## (.smc_args_problem() above refuses a lambda_schedule here: it is the
  ## LIKELIHOOD-tempering ladder and a model-tempering run has no stage for it.
  ## dynhr_smc_model_tempered() honours it on its likelihood-tempered Stage 1
  ## and withholds it from this bridge run.)

  ## phi_schedule under LIKELIHOOD tempering: an alias for
  ## lambda_schedule. Every release since 0.8.1 read phi_schedule only on the
  ## model-tempering branch, so smc(..., phi_schedule = phi) without
  ## approx_loglik_fn ran ADAPTIVE tempering and said nothing. It is what the
  ## caller meant, so honour it; two DIFFERENT fixed schedules are ambiguous.
  if (!use_model_tempering && !is.null(phi_schedule)) {
    .dynhr_inform(
      "smc: phi_schedule is used as lambda_schedule (likelihood tempering; ",
      "phi_schedule names the model-tempering schedule when approx_loglik_fn ",
      "is supplied).",
      once = TRUE, key = "smc_phi_schedule_alias",
      class = "dynhr_message_schedule_alias")
    lambda_schedule <- phi_schedule
  }

  ## Parallel backend: a single persistent mirai daemon pool serves every
  ## tempering stage (particle eval + mutation). The old future backend spun a
  ## plan up per call; mirai daemons persist across all stages. The future path
  ## is retained as a fallback (backend = "future").
  use_mirai <- isTRUE(parallel) && identical(backend, "mirai") &&
    requireNamespace("mirai", quietly = TRUE)
  ## Counter so each mutation stage draws a fresh, reproducible RNG block
  ## (stage k seeds particle i with seed_base + k * n_particles + i).
  .smc_mut_stage <- 0L
  if (use_mirai) {
    .smc_pool_setup(.mirai_n_cores(NULL, n_particles))
    on.exit(mirai::daemons(NULL), add = TRUE)
  }

  # --- Build prior sampler if needed ---
  if (is.null(prior_sampler)) {
    prior_sampler <- .smc_make_prior_sampler(prior_spec)
  }

  t_start <- Sys.time()
  n_eval  <- 0L
  ## Floor for a non-finite tempering score (M1 and, B7, M0 alike).
  ll_floor <- -1e300

  # =========================================================================
  # Stage 0: Initialise particles from prior
  # =========================================================================
  if (verbose) .dynhr_inform("SMC: Drawing initial particles from prior...")

  particles <- vector("list", n_particles)
  log_liks  <- numeric(n_particles)
  log_pris  <- numeric(n_particles)

  ## Interface check (fail loud, before any tempering). SMC tempers the
  ## LIKELIHOOD -- pi_lambda ∝ prior · lik^lambda -- so log_post_fn MUST return
  ## a list with a numeric `$loglik` (and `$logprior`). A `$logpost`-only
  ## closure, as accepted by RWMH/NUTS, would leave every particle's $loglik
  ## NULL -> floored identically to ll_floor -> no likelihood variation ->
  ## lambda jumps straight to 1 -> a SILENT GARBAGE posterior. Catch it now.
  .if_chk <- log_post_fn(prior_sampler())
  if (!is.list(.if_chk) || is.null(.if_chk$loglik) ||
        !is.numeric(.if_chk$loglik) || length(.if_chk$loglik) != 1L) {
    stop("SMC requires `log_post_fn` to return a list with a scalar numeric ",
         "`$loglik` (it tempers the likelihood: prior * lik^lambda), but the ",
         "supplied function returned ",
         if (is.list(.if_chk)) "a list without a usable `$loglik`."
         else paste0("a ", class(.if_chk)[1], "."),
         " Pass the raw make_posterior()/make_log_posterior() closure (which ",
         "returns $loglik, $logprior and $logpost) -- not a $logpost-only ",
         "wrapper.", call. = FALSE)
  }

  ## Validate approx_loglik_fn (model tempering): must return a scalar numeric.
  ## Draw one prior sample and check immediately, before any expensive work.
  if (use_model_tempering) {
    if (!is.function(approx_loglik_fn))
      stop("approx_loglik_fn must be a function.", call. = FALSE)
    .chk_approx <- approx_loglik_fn(prior_sampler())
    if (!is.numeric(.chk_approx) || length(.chk_approx) != 1L)
      stop("approx_loglik_fn must return a scalar numeric (log L_{M0}(theta)), ",
           "not a list. It is separate from log_post_fn.", call. = FALSE)
  }

  ## init_particles (model tempering): the stage-0 cloud. For an UNBIASED
  ## log(Z_M1/Z_M0) the cloud at phi = 0 must be the M0 *posterior*
  ## (pi_0 ∝ p(theta) * L_{M0}(theta)) -- the bridge then telescopes
  ## Z(pi_phi)/Z(pi_{phi-1}) to Z_M1/Z_M0 exactly. If instead the cloud is the
  ## prior (init_particles = NULL), the particle positions still converge to
  ## the M1 posterior (the mutation kernel targets pi_phi), but the marginal
  ## estimate is biased because the first increment is taken under the prior,
  ## not under pi_0. So model tempering with a valid MDD requires the caller to
  ## supply init_particles drawn from the M0 posterior (e.g. the particle cloud
  ## from a prior SMC run on M0 -- the Mlikota-Schorfheide workflow).
  init_par_names <- NULL
  if (!is.null(init_particles)) {
    if (!use_model_tempering)
      stop("init_particles is only meaningful with model tempering; ",
           "supply approx_loglik_fn too.", call. = FALSE)
    init_particles <- as.matrix(init_particles)
    if (nrow(init_particles) != n_particles)
      stop(sprintf(paste0("init_particles must have n_particles = %d rows ",
                          "(one M0-posterior draw per particle); got %d."),
                   n_particles, nrow(init_particles)), call. = FALSE)
    init_par_names <- colnames(init_particles)
    if (is.null(init_par_names)) init_par_names <- names(prior_sampler())
  }
  marginal_valid <- (!use_model_tempering) || !is.null(init_particles)
  if (use_model_tempering && is.null(init_particles) && isTRUE(verbose))
    .dynhr_inform("dynhr_smc: model tempering without init_particles -- the particle ",
            "cloud converges to the M1 posterior, but log_marginal_lik is a ",
            "VALID log(Z_M1/Z_M0) only if the M0 posterior ~ the prior. Supply ",
            "init_particles drawn from the M0 posterior for an unbiased MDD.")

  # Draw and evaluate in parallel or serial
  .eval_particle <- function(i) {
    theta <- if (!is.null(init_particles)) {
      th <- init_particles[i, ]
      names(th) <- init_par_names
      th
    } else prior_sampler()
    res   <- log_post_fn(theta)
    ## Guard: log_post_fn must return list(loglik=, logprior=).  A missing
    ## field is handled inside .smc_particle_parts() (loglik -> floor,
    ## logprior -> 0) rather than erroring with "replacement has length zero".
    pp <- .smc_particle_parts(res, ll_floor, missing_prior = 0)
    ## Model tempering: evaluate M0 alongside M1. A non-finite M0 is
    ## FLOORED exactly like M1 (it used to be -Inf, which gave an Inf
    ## increment and, at phi = 1, 0 * -Inf = NaN).
    ll0_i <- if (use_model_tempering) {
      v <- approx_loglik_fn(theta)
      if (!is.finite(v)) ll_floor else v
    } else NA_real_
    list(theta = theta, loglik = pp[["loglik"]], logprior = pp[["logprior"]],
         phi = pp[["phi"]], sp = pp[["log_sysprior"]], loglik0 = ll0_i)
  }

  if (use_mirai) {
    results <- .smc_pmap(seq_len(n_particles), .eval_particle,
                         seed_base = seed_base)
  } else if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    results <- future.apply::future_lapply(seq_len(n_particles), .eval_particle,
                                            future.seed = TRUE)
  } else {
    results <- lapply(seq_len(n_particles), .eval_particle)
  }

  ## log_liks0: M0 log-likelihoods; only allocated when model tempering active.
  ## Maintained in parallel with log_liks through resample and mutation.
  log_liks0 <- if (use_model_tempering) numeric(n_particles) else NULL
  ## `log_phi` is what is TEMPERED (log system prior + power * loglik, see
  ## .smc_particle_parts()); log_liks keeps the RAW likelihood for consumers,
  ## log_pris the PARAMETER prior, log_sys the system prior. All four (and
  ## log_liks0) move together through resampling and mutation.
  log_phi <- numeric(n_particles)
  log_sys <- numeric(n_particles)

  for (i in seq_len(n_particles)) {
    particles[[i]] <- results[[i]]$theta
    log_liks[i]    <- results[[i]]$loglik
    log_pris[i]    <- results[[i]]$logprior
    log_phi[i]     <- results[[i]]$phi
    log_sys[i]     <- results[[i]]$sp
    if (use_model_tempering) log_liks0[i] <- results[[i]]$loglik0
  }
  n_eval <- n_eval + n_particles

  d <- length(particles[[1]])
  par_names <- names(particles[[1]])

  ## Non-finite scores were floored at ll_floor in .smc_particle_parts().
  n_valid <- sum(log_phi > ll_floor)
  if (verbose) {
    .dynhr_inform(sprintf("SMC: %d/%d particles have finite likelihood", n_valid, n_particles))
  }
  ## Fail loud rather than temper from an all-infeasible population: with every
  ## loglik at the floor there is no variation, lambda jumps to 1, and the
  ## "posterior" + log_marginal are meaningless (log_mlik ~ -1e308).
  if (n_valid == 0L) {
    stop(sprintf("SMC: 0/%d initial particles have a finite likelihood -- ",
                 n_particles),
         "cannot temper from an all-infeasible population. Likely causes: the ",
         "prior_sampler draws outside the model's feasible region, or ",
         "log_post_fn never returns a finite `$loglik`. Check the prior bounds ",
         "and that the model solves at prior draws.", call. = FALSE)
  }

  # =========================================================================
  # Tempering loop
  # =========================================================================
  lambda_curr <- 0
  log_marginal <- 0
  stage <- 0L

  # Storage for diagnostics
  lambda_trace  <- numeric(0)
  ess_trace     <- numeric(0)
  accept_trace  <- numeric(0)
  scale_trace   <- numeric(0)
  log_w         <- rep(0, n_particles)  # current log weights (start uniform)

  ## Mutation proposal scale, adapted between stages toward `mut_target`
  ## acceptance (Herbst-Schorfheide). Initialised to the classic 2.38-style
  ## factor / sqrt(d); the first stage thus behaves as before, and the logistic
  ## feedback then rescues a too-wide proposal under a diffuse prior.
  mh_c <- mh_scale_factor / sqrt(d)

  while (lambda_curr < 1) {
    stage <- stage + 1L

    # --- Determine next lambda / phi ---
    ## When model tempering, .smc_next_lambda operates on the log-lik DIFFERENCE
    ## (log_liks - log_liks0) so that the ESS target governs the bridge step.
    ## The phi_schedule argument plays the same role as lambda_schedule.
    ## Score whose increment reweights the cloud. Likelihood tempering: phi.
    ## Model tempering: phi_M1 - phi_M0 with phi_M0 = log_sys + logL_M0 (M0
    ## shares the prior side, incl. the system prior, with M1). A particle
    ## at the floor under EITHER model has zero target density for phi < 1,
    ## so its score is the floor (it dies at the first positive increment)
    ## rather than floor - floor = 0 or finite + 1e300.
    if (use_model_tempering) {
      log_phi0 <- .smc_phi0(log_sys, log_liks0, ll_floor)
      score <- ifelse(log_phi <= ll_floor | log_phi0 <= ll_floor, ll_floor,
                      log_phi - log_phi0)
      schedule <- phi_schedule
    } else {
      score    <- log_phi
      schedule <- lambda_schedule
    }
    adaptive_step <- is.null(schedule) || stage > length(schedule)
    lambda_next <- if (adaptive_step)
      .smc_next_lambda(score, lambda_curr, ess_target, n_particles,
                       log_w_prev = log_w)
    else schedule[stage]
    lambda_next <- min(lambda_next, 1)
    dlambda <- lambda_next - lambda_curr

    # --- Incremental weights ---
    inc_log_w <- dlambda * score

    # --- Marginal likelihood contribution ---
    # p(Y|M) contribution = sum_i W_prev_i * exp(inc_log_w_i), where W_prev are
    # the NORMALISED weights coming INTO this stage. When the previous stage
    # resampled, log_w is uniform (all zero) and log_w_prev_norm = -log(N), so
    # this reduces to the old logsumexp(inc_log_w) - log(N) formula -- but when
    # the previous stage skipped resampling, the incoming weights are
    # non-uniform and must be folded in here (Herbst & Schorfheide 2014).
    log_w_prev_norm <- log_w - .smc_log_sum_exp(log_w)
    log_marginal <- log_marginal +
      .smc_log_sum_exp(log_w_prev_norm + inc_log_w)

    log_w <- log_w + inc_log_w

    # --- Normalise weights ---
    max_lw <- max(log_w)
    w_norm <- exp(log_w - max_lw)
    w_sum  <- sum(w_norm)
    ## Guard: if max_lw = -Inf (all particles have weight 0) or the
    ## normalising sum is non-finite/zero (numerical underflow), fall back to
    ## uniform weights.  This is the correct SMC behaviour when all particles
    ## are at the likelihood floor (every particle equally unlikely): uniform
    ## weights, full ESS, no resampling, mutation drives exploration.
    if (!is.finite(w_sum) || w_sum <= 0) {
      w_norm <- rep(1 / n_particles, n_particles)
    } else {
      w_norm <- w_norm / w_sum
    }

    # --- ESS ---
    ess <- 1 / sum(w_norm^2)
    lambda_trace <- c(lambda_trace, lambda_next)
    ess_trace    <- c(ess_trace, ess)

    if (verbose) {
      .dynhr_inform(sprintf("SMC stage %d: lambda=%.4f  ESS=%.0f/%d  log_mlik=%.2f",
                      stage, lambda_next, ess, n_particles, log_marginal))
    }

    # --- Resample if ESS below threshold ---
    ## An adaptive step that stopped short of lambda = 1 was chosen to put
    ## the COMBINED ESS exactly at the target, so the cloud has reached the
    ## resampling threshold: resample. (Not resampling would leave the next
    ## stage's bisection starting AT the target with no room to move.)
    if (ess < ess_target * n_particles || (adaptive_step && lambda_next < 1)) {
      idx <- .smc_systematic_resample(w_norm, n_particles)
      particles <- particles[idx]
      log_liks  <- log_liks[idx]
      log_pris  <- log_pris[idx]
      log_phi   <- log_phi[idx]
      log_sys   <- log_sys[idx]
      ## log_liks0 maintenance (staleness guard): resample at the SAME indices
      ## as log_liks so the M0/M1 pair stays in sync. (High-severity risk from
      ## the scope: a stale log_liks0 silently corrupts the bridge weights.)
      if (use_model_tempering) log_liks0 <- log_liks0[idx]
      log_w     <- rep(0, n_particles)  # reset weights after resampling
      w_norm    <- rep(1 / n_particles, n_particles)  # weights now uniform
    }

    # --- Mutation: RWMH with tempered target ---
    # Proposal covariance from the (resampled) particle cloud, WEIGHTED by the
    # current normalised weights w_norm (uniform if a resample just occurred,
    # non-uniform otherwise). The scale that multiplies it, mh_c, is what
    # adapts between stages (see below) -- that adaptive feedback is the
    # Herbst-Schorfheide fix for a proposal that is far too wide or too narrow
    # (e.g. under a diffuse prior).
    theta_mat <- do.call(rbind, particles)
    Sigma_hat <- stats::cov.wt(theta_mat, wt = w_norm, method = "ML")$cov
    # Regularise
    Sigma_hat <- Sigma_hat + diag(1e-6, d)
    L_prop <- t(chol(Sigma_hat))
    scale <- mh_c

    ## Herbst-Schorfheide (2014 JAE; 2015 book ch. 5) / FRBNY SMC.jl 3-component
    ## mixture proposal machinery. Only built when mixture_weights is non-NULL,
    ## so the default path neither computes these nor draws any extra RNG.
    if (use_mixture) {
      theta_bar <- as.numeric(w_norm %*% theta_mat)
      chol_S  <- L_prop  # t(chol(scale^2 * Sigma_hat)) == scale * L_prop
      chol_S  <- scale * chol_S
      Dg      <- diag(diag(Sigma_hat), d)
      chol_Dg <- scale * t(chol(Dg))
    }

    ## Tempered log-target (`phi` is the tempered score, `lp` the
    ## PARAMETER prior -- see .smc_particle_parts()).
    ## Likelihood tempering: log pi_lambda = lp + lambda * phi
    ## Model tempering:      log pi_phi    = lp + (1-phi) * (sp + ll0) + phi * phi1
    ## The ll0/sp arguments are only used when use_model_tempering = TRUE; ll0
    ## is FLOORED, so (1 - 1) * ll0 is 0 at phi = 1, never 0 * -Inf = NaN.
    .tempered_lp <- function(theta, phi, lp, ll0 = ll_floor, sp = 0) {
      if (use_model_tempering)
        lp + (1 - lambda_next) * .smc_phi0(sp, ll0, ll_floor) +
          lambda_next * phi
      else
        lp + lambda_next * phi
    }

    n_accepted_total <- 0L

    .mutate_one <- function(i) {
      theta_i <- particles[[i]]
      ll_i    <- log_liks[i]
      lp_i    <- log_pris[i]
      phi_i   <- log_phi[i]
      sp_i    <- log_sys[i]
      ## ll0_i: M0 log-likelihood for the current particle. Only used when
      ## model tempering is active.
      ll0_i   <- if (use_model_tempering) log_liks0[i] else ll_floor
      tlp_i   <- .tempered_lp(theta_i, phi_i, lp_i, ll0_i, sp_i)
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

        res_prop <- log_post_fn(theta_prop)
        ## NULL fields (missing $loglik / $logprior from a non-conforming
        ## log_post_fn) are treated as infeasible rather than crashing with
        ## "argument is of length zero" in is.finite().
        pp_prop  <- .smc_particle_parts(res_prop, ll_floor,
                                        missing_prior = -Inf)
        ll_prop  <- pp_prop[["loglik"]]
        lp_prop  <- pp_prop[["logprior"]]
        phi_prop <- pp_prop[["phi"]]
        sp_prop  <- pp_prop[["log_sysprior"]]

        ## Model tempering: evaluate M0 at the proposal. This is the one extra
        ## M0 evaluation per RWMH step (only paid when use_model_tempering).
        ll0_prop <- if (use_model_tempering) {
          v <- approx_loglik_fn(theta_prop)
          if (!is.finite(v)) ll_floor else v
        } else ll_floor

        tlp_prop <- .tempered_lp(theta_prop, phi_prop, lp_prop, ll0_prop,
                                 sp_prop)
        log_alpha <- tlp_prop - tlp_i

        if (use_mixture) {
          # Asymmetric proposal (component 3 is independence-style): correct
          # the MH ratio with log q(prop -> curr) - log q(curr -> prop).
          log_q_fwd <- .smc_mixture_logq(theta_i, theta_prop, theta_bar,
                                          chol_S, chol_Dg, mixture_weights)
          log_q_rev <- .smc_mixture_logq(theta_prop, theta_i, theta_bar,
                                          chol_S, chol_Dg, mixture_weights)
          log_alpha <- log_alpha + (log_q_rev - log_q_fwd)
        }

        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
          theta_i <- theta_prop
          ll_i    <- ll_prop
          lp_i    <- lp_prop
          phi_i   <- phi_prop
          sp_i    <- sp_prop
          ll0_i   <- ll0_prop   # log_liks0 maintenance: update on accept
          tlp_i   <- tlp_prop
          acc     <- acc + 1L
        }
      }

      list(theta = theta_i, loglik = ll_i, logprior = lp_i, phi = phi_i,
           sp = sp_i, loglik0 = ll0_i, accepted = acc)
    }

    if (use_mirai) {
      .smc_mut_stage <- .smc_mut_stage + 1L
      mut_results <- .smc_pmap(seq_len(n_particles), .mutate_one,
                               seed_base = seed_base +
                                 .smc_mut_stage * n_particles)
    } else if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
      mut_results <- future.apply::future_lapply(seq_len(n_particles), .mutate_one,
                                                  future.seed = TRUE)
    } else {
      mut_results <- lapply(seq_len(n_particles), .mutate_one)
    }

    for (i in seq_len(n_particles)) {
      particles[[i]] <- mut_results[[i]]$theta
      log_liks[i]    <- mut_results[[i]]$loglik
      log_pris[i]    <- mut_results[[i]]$logprior
      log_phi[i]     <- mut_results[[i]]$phi
      log_sys[i]     <- mut_results[[i]]$sp
      ## log_liks0 maintenance: update from the mutation result which carries
      ## the M0 value of the particle's final position (accepted or original).
      if (use_model_tempering) log_liks0[i] <- mut_results[[i]]$loglik0
      n_accepted_total <- n_accepted_total + mut_results[[i]]$accepted
    }
    n_eval <- n_eval + n_particles * n_mh_steps

    accept_rate <- n_accepted_total / (n_particles * n_mh_steps)
    accept_trace <- c(accept_trace, accept_rate)
    scale_trace  <- c(scale_trace, mh_c)

    # Adapt the proposal scale for the NEXT stage from this stage's acceptance
    # (Herbst-Schorfheide logistic feedback toward mut_target).
    mh_c <- mh_c * .smc_scale_adjust(accept_rate, target = mut_target)

    if (verbose) {
      .dynhr_inform(sprintf("  Mutation: accept=%.0f%%  scale=%.4g", accept_rate * 100, scale))
    }

    lambda_curr <- lambda_next

    # Progress
    if (!is.null(progressor)) {
      progressor(
        message = sprintf("SMC lambda=%.3f ESS=%.0f", lambda_next, ess),
        amount = 1
      )
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  # =========================================================================
  # Assemble output
  # =========================================================================
  theta_mat <- do.call(rbind, particles)
  colnames(theta_mat) <- par_names

  # Final log-posteriors: the FULL lambda = 1 target, parameter prior + phi
  # (= the closure's own $logpost on every feasible draw, including the
  # system prior and the power_posterior exponent -- what THAMES must pair
  # with each draw). A floored (infeasible) particle is -Inf, not
  # -1e300 + prior, so downstream estimators drop it rather than let an
  # exp(1e300) reciprocal swamp them.
  logpost_final <- log_pris + log_phi
  logpost_final[log_phi <= ll_floor] <- -Inf

  # Final normalised particle weights.  When the last stage resampled, these
  # are uniform (1/N); when the last stage skipped resampling, they are the
  # non-uniform w_norm from that stage.  Stored so as_posterior_draws() can
  # detect and correct for the non-uniform case.
  smc_weights_final <- w_norm  # w_norm is still in scope from the loop

  ## Model tempering MDD bookkeeping:
  ## The SMC identity accumulates log(Z_M1 / Z_M0) in log_marginal.
  ## Adding log_Z_approx (= log Z_M0 from a prior run, default 0) yields
  ## log Z_M1. When log_Z_approx = 0 the output is log(Z_M1/Z_M0).
  if (use_model_tempering) {
    log_marginal <- log_marginal + log_Z_approx
  }

  # For compatibility with rwmh() diagnostics, alias particles as "chain"
  out <- list(
    chain            = theta_mat,
    particles        = theta_mat,
    smc_weights      = smc_weights_final,
    log_liks         = log_liks,
    log_priors       = log_pris,
    log_phi          = log_phi,
    log_sysprior     = log_sys,
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
    sampler          = "smc",
    tempering_mode   = if (use_model_tempering) "model" else "likelihood",
    log_Z_approx     = log_Z_approx,
    marginal_valid   = marginal_valid
  )
  ## Prior-support box, read by thames_mdd_from_chains() for the B6
  ## ellipsoid-outside-support correction (absent without a prior_spec).
  bnds <- .smc_support_bounds(prior_spec, par_names)
  if (!is.null(bnds$lower)) {
    out$support_lower <- bnds$lower
    out$support_upper <- bnds$upper
  }
  out
}


#' Prior-support box of an SMC run, for THAMES's support correction
#'
#' @param prior_spec The prior specification (data.frame or list) or NULL.
#' @param par_names Parameter names, in column order.
#' @return \code{list(lower, upper)} named numeric vectors, or NULLs when the
#'   spec is unavailable or carries no bounds.
#' @noRd
.smc_support_bounds <- function(prior_spec, par_names) {
  none <- list(lower = NULL, upper = NULL)
  if (is.null(prior_spec) || is.null(par_names)) return(none)
  get1 <- function(nm, keys, default) {
    sp <- if (is.data.frame(prior_spec)) {
      i <- match(nm, prior_spec$name)
      if (is.na(i)) return(default)
      as.list(prior_spec[i, , drop = FALSE])
    } else if (!is.null(prior_spec[[nm]])) prior_spec[[nm]]
    else {
      hit <- Filter(function(x) identical(x$name, nm), prior_spec)
      if (!length(hit)) return(default)
      hit[[1L]]
    }
    for (k in keys) {
      v <- sp[[k]]
      if (!is.null(v) && length(v) >= 1L) {
        v <- suppressWarnings(as.numeric(v[[1L]]))
        if (!is.na(v)) return(v)
      }
    }
    default
  }
  lo <- vapply(par_names, get1, numeric(1), keys = c("lb", "lower"),
               default = -Inf)
  hi <- vapply(par_names, get1, numeric(1), keys = c("ub", "upper"),
               default = Inf)
  if (all(!is.finite(lo)) && all(!is.finite(hi))) return(none)
  list(lower = lo, upper = hi)
}


# ============================================================================
# dynhr_smc_model_tempered() -- Mlikota & Schorfheide (2024) two-stage wrapper
# ============================================================================

#' Two-stage SMC model-tempering wrapper (Mlikota & Schorfheide 2024)
#'
#' Automates the unbiased Bayes-factor workflow:
#'   Stage 1 -- SMC on M0.  Obtains log Z_M0 and the M0 posterior cloud.
#'   Stage 2 -- M0->M1 bridge.  Feeds the M0 cloud as init_particles and
#'              log Z_M0 as log_Z_approx to a model-tempered SMC run.
#'              The returned log_marginal_lik = log Z_M1 (not just the ratio).
#'
#' ## MDD bookkeeping
#'
#'   log Z_M0  <- stage1$log_marginal_lik  (standard likelihood-tempering MDD)
#'   log(Z_M1/Z_M0) <- stage2 accumulates this bridge ratio.
#'   log Z_M1  = log(Z_M1/Z_M0) + log Z_M0  (what is returned in
#'               \code{$log_marginal_lik}).
#'
#' For the MDD to be UNBIASED the stage-1 posterior cloud must be exact.
#' Use enough particles and MH steps for Stage 1 to converge (check
#' \code{$stage1$ess_schedule}).
#'
#' ## Supplying M0
#'
#' Two modes are supported:
#'
#' \describe{
#'   \item{Scalar \code{approx_loglik_fn} only (no \code{log_post_fn_M0}):}{
#'     The function \code{approx_loglik_fn(theta)} is a scalar log-likelihood
#'     for M0.  The wrapper constructs the full M0 log-posterior by combining
#'     it with the prior density supplied via \code{prior_sampler} /
#'     \code{prior_spec} (the SAME prior as M1).}
#'   \item{Full \code{log_post_fn_M0} provided:}{
#'     Must return \code{list(logpost, loglik, logprior)} exactly like the M1
#'     \code{log_post_fn}.  Use this when M0 has a different prior from M1, or
#'     when you want to provide an already-evaluated M0 posterior.}
#' }
#'
#' \strong{Cost.} Without \code{log_post_fn_M0} the Stage-1 M0 log-posterior
#' is built from \code{log_post_fn_M1} itself: every Stage-1 evaluation calls
#' \code{log_post_fn_M1(theta)} (which runs the full Kalman filter, or the
#' full likelihood, of M1) to obtain the prior and system-prior terms, and
#' then \code{approx_loglik_fn(theta)}. Stage 1 therefore costs at least as
#' much per evaluation as Stage 2, and the cheap approximating likelihood
#' saves nothing there. Supply \code{log_post_fn_M0} (prior plus cheap
#' likelihood only) to make Stage 1 cheap.
#'
#' @param log_post_fn_M1 function(theta) -> list(logpost, loglik, logprior).
#'   The M1 log-posterior.  Same interface as \code{\link{dynhr_smc}}'s
#'   \code{log_post_fn}.
#' @param approx_loglik_fn function(theta) -> scalar.  M0 log-likelihood
#'   (the approximating model).  Used both to build the M0 run (Stage 1) and
#'   as the bridge denominator in the M1 run (Stage 2).
#' @param prior_sampler function() -> named numeric vector.  Prior sampler
#'   used for both stages.  Required when \code{log_post_fn_M0 = NULL}.
#' @param prior_spec Prior specification (data.frame or named list).
#'   Alternative to \code{prior_sampler}; forwarded to both stages.
#' @param log_post_fn_M0 Optional.  If supplied, used directly as the Stage-1
#'   log-posterior (must return list(logpost, loglik, logprior)).  When NULL
#'   (default), the M0 log-posterior is constructed from
#'   \code{approx_loglik_fn} + the prior.
#' @param n_particles Number of particles (same for both stages).
#' @param ess_target ESS ratio target for adaptive tempering.
#' @param n_mh_steps RWMH mutation steps per stage.
#' @param seed_M0 Integer seed for Stage 1 (M0 run).  Default 1L.
#' @param seed_M1 Integer seed for Stage 2 (M1 bridge run).  Default 2L.
#'   Both stages are seeded locally (results are reproducible for fixed
#'   seeds) and the caller's global RNG stream is restored on exit.
#' @param verbose Print progress messages.
#' @param ... Additional arguments forwarded to BOTH \code{dynhr_smc} calls
#'   (e.g. \code{mh_scale_factor}, \code{mixture_weights}, \code{parallel},
#'   \code{backend}, \code{n_mh_steps}, \code{ess_target} -- note that the
#'   last two are also explicit arguments above and take precedence when
#'   supplied directly), except \code{phi_schedule}, which is the Stage-2
#'   bridge schedule and is withheld from Stage 1, and \code{lambda_schedule},
#'   which is the Stage-1 likelihood-tempering ladder and is withheld from
#'   Stage 2.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{log_marginal_lik}}{log Z_M1 (the unbiased marginal likelihood
#'     of M1, NOT just the ratio).  Equals Stage-2 log(Z_M1/Z_M0) +
#'     log Z_M0 from Stage 1.}
#'   \item{\code{chain}, \code{particles}, \code{smc_weights}, ...}{
#'     All fields from the Stage-2 SMC result, for compatibility with
#'     \code{new_dynhr_chains()}.}
#'   \item{\code{stage1}}{Full Stage-1 SMC result for inspection.}
#'   \item{\code{log_Z_M0}}{log Z_M0 from Stage 1.}
#'   \item{\code{log_ratio}}{log(Z_M1/Z_M0) from Stage 2 (before adding
#'     log Z_M0).}
#'   \item{\code{marginal_valid}}{Always \code{TRUE}: the M0 cloud is used as
#'     \code{init_particles}, so the bridge is unbiased.}
#'   \item{\code{tempering_mode}}{"model" (set by Stage 2).}
#' }
#' @noRd
dynhr_smc_model_tempered <- function(
    log_post_fn_M1,
    approx_loglik_fn,
    prior_sampler     = NULL,
    prior_spec        = NULL,
    log_post_fn_M0    = NULL,
    n_particles       = 2000L,
    ess_target        = 0.5,
    n_mh_steps        = 1L,
    seed_M0           = 1L,
    seed_M1           = 2L,
    verbose           = TRUE,
    ...
) {
  stopifnot(is.function(log_post_fn_M1))
  stopifnot(is.function(approx_loglik_fn))
  if (is.null(prior_sampler) && is.null(prior_spec) && is.null(log_post_fn_M0))
    stop("Supply prior_sampler or prior_spec (and optionally log_post_fn_M0).",
         call. = FALSE)

  ## ── Stage 1: SMC on M0 ──────────────────────────────────────────────────
  ## Build M0 log-posterior if the caller didn't supply a full one.
  ## We need list(logpost, loglik, logprior) from a scalar approx_loglik_fn.
  ## The prior sampler already knows how to evaluate the prior density (we
  ## reconstruct it via the same prior_spec / prior_sampler) but that only
  ## gives us a SAMPLER, not a density evaluator.  The cleanest way to build
  ## logprior(theta) without duplicating the prior_spec machinery is to wrap
  ## approx_loglik_fn as the loglik and let the M1 log_post_fn supply the
  ## prior density at the same theta.
  ##
  ## Trick: call log_post_fn_M1(theta) to get $logprior (it's the same prior),
  ## then swap in the M0 loglik.  This avoids reimplementing a prior-density
  ## evaluator and is exact as long as M0 and M1 share the same prior.
  if (is.null(log_post_fn_M0)) {
    ## The M0 target shares M1's PRIOR SIDE -- parameter prior and system
    ## prior -- and swaps only the likelihood (untempered). The components are
    ## passed on as "posterior_parts" so stage 1 samples p * s * L_M0 from the
    ## parameter prior, the same M0 the stage-2 bridge assumes
    ## (.smc_phi0()); reading M1's $logprior alone folded the system prior
    ## into the stage-0 base in "lp" mode and dropped it in "extra" mode.
    log_post_fn_M0 <- function(theta) {
      lp_M1 <- log_post_fn_M1(theta)
      pp1   <- .posterior_parts(lp_M1)
      lp0   <- lp_M1$logprior
      ll0   <- approx_loglik_fn(theta)
      if (!is.finite(ll0)) ll0 <- -1e300
      out <- list(logpost = lp0 + ll0, loglik = ll0, logprior = lp0)
      attr(out, "posterior_parts") <- c(
        logprior_param = pp1[["logprior_param"]],
        log_sysprior   = pp1[["log_sysprior"]],
        loglik         = ll0,
        power          = 1)
      out
    }
  }

  if (verbose) .dynhr_inform("SMC model tempering: Stage 1 -- M0 run")
  .local_seed(seed_M0)  # the caller's RNG stream is restored on exit
  ## phi_schedule is the Stage-2 BRIDGE schedule. Stage 1 is likelihood
  ## tempering, where dynhr_smc() now reads phi_schedule as an alias for
  ## lambda_schedule -- so it must not reach Stage 1.
  dots1 <- list(...)
  dots1$phi_schedule <- NULL
  stage1 <- do.call(dynhr_smc, c(list(
    log_post_fn   = log_post_fn_M0,
    prior_sampler = prior_sampler,
    prior_spec    = prior_spec,
    n_particles   = n_particles,
    ess_target    = ess_target,
    n_mh_steps    = n_mh_steps,
    seed_base     = seed_M0,
    verbose       = verbose), dots1))
  log_Z_M0  <- stage1$log_marginal_lik
  ## Extract the M0 posterior particle cloud (rows = particles, cols = params)
  init_cloud <- stage1$chain   # already a matrix (n_particles x d)

  ## ── Stage 2: M0->M1 bridge ──────────────────────────────────────────────
  if (verbose) {
    .dynhr_inform(sprintf(
      "SMC model tempering: Stage 2 -- M0->M1 bridge  (log Z_M0 = %.3f)",
      log_Z_M0
    ))
  }
  .local_seed(seed_M1)
  ## lambda_schedule governs the likelihood-tempered Stage 1 only; the bridge
  ## is governed by phi_schedule, so lambda_schedule is withheld here.
  dots2 <- list(...)
  dots2$lambda_schedule <- NULL
  stage2 <- do.call(dynhr_smc, c(list(
    log_post_fn      = log_post_fn_M1,
    prior_sampler    = prior_sampler,
    prior_spec       = prior_spec,
    n_particles      = n_particles,
    ess_target       = ess_target,
    n_mh_steps       = n_mh_steps,
    approx_loglik_fn = approx_loglik_fn,
    init_particles   = init_cloud,
    log_Z_approx     = log_Z_M0,
    seed_base        = seed_M1,
    verbose          = verbose), dots2))

  ## ── Assemble output ─────────────────────────────────────────────────────
  ## Expose log_ratio (the bridge accumulation before adding log_Z_M0)
  log_ratio <- stage2$log_marginal_lik - log_Z_M0

  ## Start from the stage2 result (contains all the M1 cloud fields)
  out <- stage2
  ## Override / add wrapper-specific fields
  out$log_marginal_lik <- stage2$log_marginal_lik   # = log Z_M1
  out$log_Z_M0         <- log_Z_M0
  out$log_ratio        <- log_ratio
  out$stage1           <- stage1
  out$marginal_valid   <- TRUE   # always: init_particles supplied from M0 posterior
  out$tempering_mode   <- "model"
  out
}


# ============================================================================
# Diagnostics
# ============================================================================

#' Print summary of SMC results
#'
#' @param result Output from dynhr_smc()
#' @param probs Quantiles to report
#' @noRd
smc_summary <- function(result, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  cat("\n=== SMC Summary ===\n")
  cat(sprintf("  Particles: %d\n", result$n_particles))
  cat(sprintf("  Tempering stages: %d\n", result$n_stages))
  cat(sprintf("  Lambda schedule: %s\n",
              paste(round(result$lambda_schedule, 3), collapse = " -> ")))
  cat(sprintf("  ESS range: %.0f - %.0f\n",
              min(result$ess_schedule), max(result$ess_schedule)))
  cat(sprintf("  Mutation acceptance: %.0f%% - %.0f%%\n",
              min(result$accept_schedule) * 100,
              max(result$accept_schedule) * 100))
  cat(sprintf("  Log marginal likelihood: %.2f\n", result$log_marginal_lik))
  cat(sprintf("  Total evaluations: %d\n", result$n_eval))
  cat(sprintf("  Wall time: %.1f sec\n", result$elapsed_secs))

  # Parameter summary table
  ch <- result$chain
  qmat <- t(apply(ch, 2, quantile, probs = probs))
  means <- colMeans(ch)
  sds   <- apply(ch, 2, sd)

  out <- data.frame(
    mean = round(means, 4),
    sd   = round(sds, 4),
    round(qmat, 4),
    check.names = FALSE
  )

  cat("\nParameter estimates:\n")
  print(out)
  invisible(out)
}


#' Plot SMC tempering diagnostics
#'
#' @param result Output from dynhr_smc()
#' @noRd
smc_plot_diagnostics <- function(result) {
  if (!requireNamespace("graphics", quietly = TRUE)) {
    .dynhr_inform("Base graphics required for smc_plot_diagnostics()")
    return(invisible(NULL))
  }

  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(op))

  # 1. Lambda schedule
  plot(seq_along(result$lambda_schedule), result$lambda_schedule,
       type = "b", pch = 16, col = "steelblue",
       xlab = "Stage", ylab = expression(lambda),
       main = "Tempering Schedule")
  abline(h = 1, lty = 2, col = "grey60")

  # 2. ESS over stages
  plot(seq_along(result$ess_schedule), result$ess_schedule,
       type = "b", pch = 16, col = "darkgreen",
       xlab = "Stage", ylab = "ESS",
       main = "Effective Sample Size")
  abline(h = result$n_particles * 0.5, lty = 2, col = "red")

  # 3. Mutation acceptance rates
  plot(seq_along(result$accept_schedule), result$accept_schedule * 100,
       type = "b", pch = 16, col = "darkorange",
       xlab = "Stage", ylab = "Acceptance %",
       main = "Mutation Acceptance Rate")
  abline(h = 25, lty = 2, col = "grey60")

  # 4. Log marginal likelihood accumulation
  # Reconstruct cumulative from lambda schedule
  .dynhr_cat(sprintf("  Final log marginal likelihood: %.2f\n", result$log_marginal_lik))

  # Histogram of final log-posteriors
  hist(result$post_logpost[is.finite(result$post_logpost)],
       breaks = 40, col = "lightblue", border = "white",
       main = "Final Log-Posterior", xlab = "log p(theta|Y)")

  invisible(NULL)
}
