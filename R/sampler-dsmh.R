## R/sampler-dsmh.R
## --------------------------------------------------------------------------
## dynhr_dsmh() -- Dynamic Striated Metropolis-Hastings (DSMH) sampler.
##
## Waggoner, Wu & Zha (2016), "Striated Metropolis-Hastings sampler for
## high-dimensional models", J. Econometrics 192(2): 406-420, as implemented
## in Dynare 7 (matlab/estimation/smc/dsmh.m, posterior_sampling_method =
## 'dsmh'). Helpers: .dsmh_tempered, .dsmh_prop_factor, .dsmh_tune_c,
## .dsmh_args_problem; per-chain tasks .dsmh_tune_task / .dsmh_mutate_task,
## dispatched in-session or on a mirai pool by .dsmh_run_groups (n_cores).
## --------------------------------------------------------------------------


#' Tempered log kernel of one DSMH state
#'
#' \eqn{\log p(\theta) + \lambda \phi(\theta)}, \eqn{-\infty} for an
#' infeasible state (parameter prior off support or the tempered score at
#' the floor).
#' @param lp Parameter log prior (vector).
#' @param phi Tempered score (vector; see \code{.smc_particle_parts()}).
#' @param lambda Tempering exponent.
#' @param ll_floor The floor that marks an infeasible score.
#' @return Numeric vector.
#' @noRd
.dsmh_tempered <- function(lp, phi, lambda, ll_floor = -1e300) {
  out <- lp + lambda * phi
  out[!is.finite(lp) | phi <= ll_floor] <- -Inf
  out
}


#' Square-root factor of the weighted particle covariance
#'
#' Returns F with F F' = Omega, Omega the importance-weighted covariance of
#' the particle cloud (Dynare: \code{Omegachol}). Built from the symmetric
#' eigendecomposition with eigenvalues floored at \code{1e-12} times the
#' largest (and at \code{1e-300}), so a degenerate cloud (e.g. every particle
#' on one point in some direction) still gives a usable proposal instead of a
#' failed Cholesky.
#' @param theta n x d particle matrix.
#' @param w Normalised weights (length n).
#' @return d x d matrix F.
#' @noRd
.dsmh_prop_factor <- function(theta, w) {
  mu <- colSums(theta * w)
  z  <- sweep(theta, 2L, mu)
  Om <- crossprod(z * sqrt(w))
  Om <- (Om + t(Om)) / 2
  e  <- eigen(Om, symmetric = TRUE)
  ev <- pmax(e$values, max(e$values, 0) * 1e-12, 1e-300)
  e$vectors %*% diag(sqrt(ev), nrow = length(ev))
}


#' Proposal step \eqn{F z} without BLAS
#'
#' \code{rowSums()} rather than a matrix product, so a chain's proposal is
#' bit-identical whether it runs in this session or on a mirai daemon whose
#' BLAS threading differs (the daemons pin single-threaded BLAS).
#' @param Fm d x d factor.
#' @param z Length-d vector.
#' @return Length-d vector.
#' @noRd
.dsmh_mv <- function(Fm, z) rowSums(Fm * rep(z, each = nrow(Fm)))


#' Density stratum (1..M) of tempered log densities
#' @param t Tempered log density(ies) under \eqn{\pi_{i-1}}.
#' @param cuts The \code{M - 1} stratum boundaries.
#' @noRd
.dsmh_stratum <- function(t, cuts) findInterval(t, cuts, left.open = TRUE) + 1L


#' One DSMH evaluation: parameter prior, tempered score, raw loglik
#'
#' A point outside the prior support bounds is rejected without calling
#' \code{lpf} (\code{evaluated = 0}).
#' @param theta Parameter vector.
#' @param lpf Log-posterior closure.
#' @param par_names Parameter names.
#' @param bnds \code{.smc_support_bounds()} result.
#' @param ll_floor Infeasibility floor.
#' @return \code{list(lp, phi, ll, sp, evaluated)}.
#' @noRd
.dsmh_eval <- function(theta, lpf, par_names, bnds, ll_floor = -1e300) {
  names(theta) <- par_names
  if (!is.null(bnds$lower) &&
      (any(theta < bnds$lower) || any(theta > bnds$upper)))
    return(list(lp = -Inf, phi = ll_floor, ll = ll_floor, sp = 0,
                evaluated = 0L))
  pp <- .smc_particle_parts(lpf(theta), ll_floor, missing_prior = -Inf)
  list(lp = pp[["logprior"]], phi = pp[["phi"]], ll = pp[["loglik"]],
       sp = pp[["log_sysprior"]], evaluated = 1L)
}


#' Evaluate \code{expr} leaving the global RNG state as it found it
#'
#' Restores the global \code{.Random.seed} -- or its absence -- on exit.
#' Wraps everything the chain-group backend does between two host draws
#' (a chain run in this session; \code{mirai::daemons()} and
#' \code{mirai::mirai_map()}, which can advance the host stream), so the
#' host stream -- and hence every resampling draw and chain seed -- is the
#' same whichever backend runs the chains.
#' @noRd
.dsmh_keep_rng <- function(expr) {
  ge  <- globalenv()
  had <- exists(".Random.seed", envir = ge, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = ge, inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = ge)
    } else if (exists(".Random.seed", envir = ge, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = ge)
    }
  }, add = TRUE)
  expr
}


#' Evaluate \code{expr} on a chain's own Mersenne-Twister stream
#'
#' \code{set.seed(seed)} with the kinds pinned (a mirai daemon defaults to
#' L'Ecuyer-CMRG) inside \code{.dsmh_keep_rng()}.
#' @noRd
.dsmh_with_seed <- function(seed, expr) {
  .dsmh_keep_rng({
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    expr
  })
}


#' Draw one chain seed per group from the host stream
#' @noRd
.dsmh_group_seeds <- function(G) floor(stats::runif(G) * 2147483646) + 1


#' The log-posterior closure a chain task evaluates
#'
#' \code{job$lpf} in-session. On the daemon pool \code{job$lpf} is
#' \code{NULL} (the closure is not re-shipped with every task) and the
#' closure installed once per daemon by \code{.mirai_pool_closure()} as
#' \code{.worker_lp} is used.
#' @noRd
.dsmh_task_lpf <- function(job) {
  if (!is.null(job$lpf)) return(job$lpf)
  lpf <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
  if (!is.function(lpf))
    .dynhr_abort("dynhr_dsmh: no `.worker_lp` on this worker.",
                 class = "dynhr_error_dsmh_worker")
  lpf
}


#' One c_i tuning chain: \code{K} random-walk steps from its resampled start
#' @param j Chain (group) index.
#' @param job Shared inputs: \code{th} (G x d starts), \code{tl} (their
#'   tempered log densities), \code{cc}, \code{Fm}, \code{lambda}, \code{K},
#'   \code{par_names}, \code{bnds}, \code{lpf}.
#' @param seeds Per-group seeds.
#' @return \code{list(n_acc, n_eval)}.
#' @noRd
.dsmh_tune_task <- function(j, job, seeds) {
  lpf <- .dsmh_task_lpf(job)
  .dsmh_with_seed(seeds[[j]], {
    th <- job$th[j, ]
    tl <- job$tl[[j]]
    d  <- length(th)
    n_acc  <- 0L
    n_eval <- 0L
    for (k in seq_len(job$K)) {
      cand <- th + sqrt(job$cc) * .dsmh_mv(job$Fm, stats::rnorm(d))
      ev   <- .dsmh_eval(cand, lpf, job$par_names, job$bnds)
      n_eval <- n_eval + ev$evaluated
      tlx  <- .dsmh_tempered(ev$lp, ev$phi, job$lambda)
      if (is.finite(tlx) && stats::runif(1) < exp(tlx - tl)) {
        th <- cand
        tl <- tlx
        n_acc <- n_acc + 1L
      }
    }
    list(n_acc = n_acc, n_eval = n_eval)
  })
}


#' One striated Metropolis chain of a DSMH stage
#'
#' Runs \code{n_keep * tau} steps from the chain's resampled start and keeps
#' every \code{tau}-th state. Each step is either a striated jump
#' (probability \code{p_jump}) to a uniformly drawn member of the chain's own
#' density stratum of the previous cloud, or a random-walk step.
#' @param j Chain (group) index.
#' @param job Shared inputs: \code{pool} (previous cloud: \code{theta},
#'   \code{lp}, \code{phi}, \code{ll}, \code{sp}), \code{members},
#'   \code{cuts}, the resampled starts \code{th0}, \code{lp0}, \code{phi0},
#'   \code{ll0}, \code{sp0}, \code{lam}, \code{lam_prev}, \code{cc},
#'   \code{Fm}, \code{p_jump}, \code{n_keep}, \code{tau}, \code{par_names},
#'   \code{bnds}, \code{lpf}.
#' @param seeds Per-group seeds.
#' @return The kept states (\code{theta}, \code{lp}, \code{phi}, \code{ll},
#'   \code{sp}) and the move counters.
#' @noRd
.dsmh_mutate_task <- function(j, job, seeds) {
  lpf <- .dsmh_task_lpf(job)
  .dsmh_with_seed(seeds[[j]], {
    pool <- job$pool
    lam  <- job$lam
    dlam <- lam - job$lam_prev
    th   <- job$th0[j, ]
    lp0  <- job$lp0[[j]]; phi0 <- job$phi0[[j]]; ll0 <- job$ll0[[j]]
    sp0  <- job$sp0[[j]]
    tc0  <- .dsmh_tempered(lp0, phi0, lam)
    d    <- length(th)
    n_keep <- job$n_keep
    tau    <- job$tau
    k_theta <- matrix(NA_real_, n_keep, d)
    k_lp <- k_phi <- k_ll <- k_sp <- numeric(n_keep)
    n_rw <- n_rw_acc <- n_jp <- n_jp_acc <- n_eval <- 0L
    for (l in seq_len(n_keep * tau)) {
      if (stats::runif(1) < job$p_jump) {
        n_jp <- n_jp + 1L
        k   <- .dsmh_stratum(.dsmh_tempered(lp0, phi0, job$lam_prev),
                             job$cuts)
        mem <- job$members[[k]]
        if (length(mem) > 0L) {
          q <- mem[ceiling(stats::runif(1) * length(mem))]
          if (stats::runif(1) < exp(dlam * (pool$phi[q] - phi0))) {
            th  <- pool$theta[q, ]
            lp0 <- pool$lp[q]; phi0 <- pool$phi[q]; ll0 <- pool$ll[q]
            sp0 <- pool$sp[q]
            tc0 <- .dsmh_tempered(lp0, phi0, lam)
            n_jp_acc <- n_jp_acc + 1L
          }
        }
      } else {
        n_rw <- n_rw + 1L
        cand <- th + sqrt(job$cc) * .dsmh_mv(job$Fm, stats::rnorm(d))
        ev   <- .dsmh_eval(cand, lpf, job$par_names, job$bnds)
        n_eval <- n_eval + ev$evaluated
        tcx  <- .dsmh_tempered(ev$lp, ev$phi, lam)
        if (is.finite(tcx) && stats::runif(1) < exp(tcx - tc0)) {
          th  <- cand
          lp0 <- ev$lp; phi0 <- ev$phi; ll0 <- ev$ll; sp0 <- ev$sp
          tc0 <- tcx
          n_rw_acc <- n_rw_acc + 1L
        }
      }
      if (l %% tau == 0L) {
        r <- l %/% tau
        k_theta[r, ] <- th
        k_lp[r] <- lp0; k_phi[r] <- phi0; k_ll[r] <- ll0; k_sp[r] <- sp0
      }
    }
    list(theta = k_theta, lp = k_lp, phi = k_phi, ll = k_ll, sp = k_sp,
         n_rw = n_rw, n_rw_acc = n_rw_acc, n_jp = n_jp, n_jp_acc = n_jp_acc,
         n_eval = n_eval)
  })
}


#' Run the chain groups of one DSMH step, in-session or on the daemon pool
#'
#' Every chain draws only from its own seeded stream (\code{seeds}, drawn on
#' the host), so the result is identical whichever way -- and on however
#' many daemons -- the groups run. On the pool \code{job$lpf} must be
#' \code{NULL}: the daemons evaluate their \code{.worker_lp}.
#' @param task \code{.dsmh_tune_task} or \code{.dsmh_mutate_task}.
#' @param job Shared task inputs.
#' @param seeds Per-group seeds (one task per seed).
#' @param parallel Dispatch with \code{mirai::mirai_map()}.
#' @return List of task results, in group order.
#' @noRd
.dsmh_run_groups <- function(task, job, seeds, parallel) {
  ## Force the seeds NOW: a lazy host draw first forced inside
  ## .dsmh_keep_rng() below would be rolled back with mirai's RNG use.
  force(seeds)
  if (!parallel)
    return(lapply(seq_along(seeds), task, job = job, seeds = seeds))
  res <- .dsmh_keep_rng(
    mirai::mirai_map(seq_along(seeds), task,
                     .args = list(job = job, seeds = seeds))[])
  for (j in seq_along(res)) {
    if (inherits(res[[j]], "miraiError") || inherits(res[[j]], "errorValue"))
      .dynhr_abort("dynhr_dsmh: chain group ", j, " failed on a mirai ",
                   "daemon: ", paste(as.character(res[[j]]), collapse = " "),
                   class = "dynhr_error_dsmh_worker")
  }
  res
}


#' Start the DSMH daemon pool
#'
#' \code{.mirai_pool_closure()}: \code{n_cores} daemons that load the
#' INSTALLED dynhr, pass its worker-version check (a skewed build aborts
#' with class \code{dynhr_error_worker_version_skew}), replay the host's
#' \code{dynhr_set_options()} state and hold \code{log_post_fn} as
#' \code{.worker_lp}. The caller tears the pool down.
#' @noRd
.dsmh_pool_start <- function(n_cores, log_post_fn) {
  .mirai_pool_closure(n_cores, log_post_fn)
}


#' Tune the DSMH random-walk scale c_i (Dynare dsmh.m `tune_c`)
#'
#' Resamples \code{G} starting points from the importance weights, runs
#' \code{K} random-walk MH steps from each with proposal
#' \eqn{\theta + \sqrt{c} F z}, and accepts \code{c} once the acceptance rate
#' lies in \code{[alpha0, alpha1]}. Otherwise \code{c} is divided by 5
#' (acceptance below \eqn{m^5}, \eqn{m = (\alpha_0 + \alpha_1)/2}), multiplied
#' by 5 (above \eqn{m^{1/5}}), or rescaled by \eqn{\log m / \log a} in
#' between -- Dynare's rule verbatim. These tuning chains are discarded.
#' Capped at \code{max_iter} tries (Dynare loops without bound); the last
#' \code{c} is then kept with a warning.
#' The resampling draw and the per-chain seeds come from the host stream;
#' each of the \code{G} tuning chains then runs on its own seeded stream
#' (\code{.dsmh_tune_task}), in-session or on the daemon pool.
#' @param job_base Shared chain-task inputs (\code{par_names}, \code{bnds},
#'   \code{lpf}).
#' @param parallel Run the chains on the mirai pool.
#' @return \code{list(c, accept, converged, n_eval)}.
#' @noRd
.dsmh_tune_c <- function(c0, pool, w, lambda, Fm, G, K,
                         accept_target, max_iter, verbose, job_base,
                         parallel) {
  a0  <- accept_target[1L]; a1 <- accept_target[2L]
  mid <- 0.5 * (a0 + a1)
  lower_prob <- mid^5
  upper_prob <- mid^(1 / 5)
  cc  <- c0
  n_eval <- 0L
  acpt <- NA_real_
  for (it in seq_len(max_iter)) {
    idx <- .smc_systematic_resample(w, G)
    job <- c(job_base, list(
      th = pool$theta[idx, , drop = FALSE],
      tl = .dsmh_tempered(pool$lp[idx], pool$phi[idx], lambda),
      cc = cc, Fm = Fm, lambda = lambda, K = K))
    seeds <- .dsmh_group_seeds(G)
    res <- .dsmh_run_groups(.dsmh_tune_task, job, seeds, parallel)
    n_acc  <- sum(vapply(res, function(r) r$n_acc, integer(1)))
    n_eval <- n_eval + sum(vapply(res, function(r) r$n_eval, integer(1)))
    acpt <- n_acc / (G * K)
    if (verbose)
      .dynhr_inform(sprintf("  DSMH tune: c = %.4g, acceptance = %.3f", cc, acpt))
    if (a0 <= acpt && acpt <= a1)
      return(list(c = cc, accept = acpt, converged = TRUE, n_eval = n_eval))
    cc <- if (acpt < lower_prob) cc / 5
          else if (acpt <= upper_prob) cc * log(mid) / log(acpt)
          else cc * 5
  }
  list(c = cc, accept = acpt, converged = FALSE, n_eval = n_eval)
}


## Every argument rule of dynhr_dsmh() that needs no likelihood evaluation, in
## one place: called by dynhr_dsmh() at entry and by validate_spec(). `args` is
## a named list keyed by dynhr_dsmh()'s argument names (absent = default =
## valid). Returns the problems as a character vector, character(0) if fine.
## The "no n_obs / lambda1 / lambda_schedule" rule stays in dynhr_dsmh(): the
## runner fills n_obs from the data, so a spec cannot be judged on it.
#' @noRd
.dsmh_args_problem <- function(args, n_par = NULL) {
  d <- formals(dynhr_dsmh)
  get <- function(nm) if (nm %in% names(args)) args[[nm]] else eval(d[[nm]])
  has <- function(nm) !is.null(args[[nm]])
  p <- character(0)
  for (nm in c("n_tune_steps", "thin", "max_tune_iter"))
    p <- c(p, .smc_whole_problem(get(nm), nm, 1L))
  if (has("n_cores"))
    p <- c(p, .smc_whole_problem(args$n_cores, "n_cores", 1L))
  size <- .dsmh_size_problem(get("n_particles"), get("n_groups"),
                             get("n_strata"))
  if (!is.null(size)) p <- c(p, size)
  at <- get("accept_target")
  if (!is.numeric(at) || length(at) != 2L || !all(is.finite(at)) ||
      at[1L] <= 0 || at[2L] >= 1 || at[1L] >= at[2L])
    p <- c(p, paste0("`accept_target` must be c(alpha0, alpha1) with ",
                     "0 < alpha0 < alpha1 < 1 (got ",
                     paste(format(at), collapse = ", "), ")."))
  if (has("lambda_schedule")) {
    p <- c(p, .smc_schedule_problem(args$lambda_schedule, "lambda_schedule",
                                         end_at_one = TRUE))
  } else {
    p <- c(p, .smc_whole_problem(get("n_stages"), "n_stages", 1L))
    if (has("lambda1")) {
      l1 <- args$lambda1
      if (!is.numeric(l1) || length(l1) != 1L || !is.finite(l1) || l1 <= 0 ||
          l1 > 1)
        p <- c(p, paste0("`lambda1` must be in (0, 1] (got ",
                         paste(format(l1), collapse = ", "), ")."))
    } else if (has("n_obs")) {
      p <- c(p, .smc_whole_problem(args$n_obs, "n_obs", 1L))
    }
  }
  p
}


## The cloud-size rules: n_particles, n_groups and n_strata whole numbers
## (n_particles >= 2), n_particles a multiple of n_groups (each of the G chains
## yields n_particles / G particles) and of n_strata. Returns the problem as a
## message, or NULL. Shared by dynhr_dsmh() and validate_spec(), which checks
## a dsmh sampler_spec when the spec is built instead of after the mode stage.
#' @noRd
.dsmh_size_problem <- function(n_particles, n_groups, n_strata) {
  whole <- function(x, min) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x == round(x) && x >= min
  for (a in list(list(n_particles, "n_particles", 2L),
                 list(n_groups, "n_groups", 1L), list(n_strata, "n_strata", 1L)))
    if (!whole(a[[1L]], a[[3L]]))
      return(paste0("`", a[[2L]], "` must be a whole number >= ", a[[3L]],
                    " (got ", paste(format(a[[1L]]), collapse = ", "), ")."))
  if (n_particles %% n_groups != 0)
    return(paste0("`n_particles` (", n_particles, ") must be a multiple of ",
                  "`n_groups` (", n_groups, ")."))
  if (n_particles %% n_strata != 0)
    return(paste0("`n_strata` (", n_strata, ") must divide `n_particles` (",
                  n_particles, ") evenly (Dynare: M must divide N*G)."))
  NULL
}


#' Dynamic Striated Metropolis-Hastings (DSMH) posterior sampler
#'
#' Samples the posterior with the Dynamic Striated Metropolis-Hastings
#' algorithm of Waggoner, Wu & Zha (2016), the \code{'dsmh'}
#' \code{posterior_sampling_method} of Dynare 7. Like \code{\link{smc}} it
#' moves a particle cloud from the prior to the posterior along a likelihood
#' tempering ladder \eqn{\pi_i \propto p(\theta)\,L(\theta)^{\lambda_i}},
#' \eqn{0 < \lambda_1 < \dots < \lambda_H = 1}; unlike SMC, each stage
#' regenerates the WHOLE cloud by running \code{n_groups} Metropolis chains
#' whose moves mix
#' \itemize{
#'   \item a random walk \eqn{\theta + \sqrt{c_i} F z} with \eqn{F F'} the
#'     importance-weighted particle covariance and \eqn{c_i} tuned to the
#'     acceptance band \code{accept_target}, and
#'   \item (probability \eqn{0.1/}\code{thin} per step) a \emph{striated}
#'     jump: the previous stage's cloud is cut into \code{n_strata} equally
#'     populated strata by its density level under \eqn{\pi_{i-1}}, and the
#'     chain jumps to a uniformly drawn particle of its OWN stratum,
#'     accepted with probability
#'     \eqn{\min\{1, \exp[(\lambda_i - \lambda_{i-1})(\log L' - \log L)]\}}.
#' }
#' A stratum collects points of equal \eqn{\pi_{i-1}} density wherever they
#' are, so the striated jump carries a chain between well-separated modes
#' that a random walk cannot cross -- the property that lets DSMH weight the
#' modes of a multimodal posterior correctly.
#'
#' @section Algorithm and Dynare parity:
#'   Stage \eqn{i}: (1) importance weights
#'   \eqn{w_j \propto L_j^{\lambda_i - \lambda_{i-1}}} on the cloud from stage
#'   \eqn{i - 1} (the prior at \eqn{i = 1}), which also update the log
#'   marginal likelihood estimate \eqn{\log \hat Z \mathrel{+}= \log
#'   \overline{w}}; (2) tune \eqn{c_i}; (3) resample \code{n_groups} chain
#'   starts from \eqn{w} (systematic resampling = Dynare's \code{kitagawa}),
#'   run \code{n_particles / n_groups * thin} steps and keep every
#'   \code{thin}-th state of every chain as the new cloud. The tempering
#'   ladder is Dynare's \eqn{\lambda_i = \lambda_1^{(H - i)/(H - 1)}} with
#'   \eqn{\lambda_1 = 1/(10\, d\, T)} (\eqn{d} parameters, \eqn{T}
#'   observations). Defaults are Dynare's (\code{H = 25}, \code{M = 20},
#'   \code{N * G = 200 * 10}, \code{K = 10}, \code{tau = 10},
#'   \code{alpha0/alpha1 = 0.2/0.3}). Deliberate deviations from
#'   \code{dsmh.m}: strata are assigned by density VALUE for both the chain
#'   state and the pool (Dynare indexes the pool by sorted position, which
#'   disagrees with the value rule on ties -- and MH output has ties); the
#'   striated draw is uniform over the whole stratum and uses a uniform
#'   independent of the acceptance uniform (Dynare reuses one uniform for
#'   both, and its index formula almost never picks a stratum's last
#'   member); the \eqn{c_i} search is capped at \code{max_tune_iter} tries.
#'   Stage 0 redraws prior draws with an infeasible likelihood (as Dynare
#'   does) and adds \eqn{\log} of the feasible fraction to \eqn{\log \hat Z},
#'   so the estimate stays the marginal likelihood under the full prior.
#'
#' @section Parallel runs and seeding:
#'   The \code{n_groups} chains of a stage are independent given the
#'   previous cloud, so with \code{n_cores >= 2} each stage's chains -- and
#'   each \eqn{c_i} tuning try's chains -- run on a pool of \code{n_cores}
#'   mirai daemons (capped at \code{n_groups}); the weights, strata,
#'   resampling and \eqn{\log \hat Z} stay on the host. Every chain runs on
#'   its own Mersenne-Twister stream, seeded from a per-chain seed that the
#'   host draws, so a seeded run gives IDENTICAL output for every
#'   \code{n_cores}, serial included. The daemons load the INSTALLED dynhr
#'   and receive \code{log_post_fn} once; a pool whose installed build
#'   differs from the session's is refused (class
#'   \code{dynhr_error_worker_version_skew}; see \code{\link{dynhr_sitrep}}).
#'   The closure must be self-contained and serialisable: on a daemon it
#'   cannot see this session's global variables (build it inside a function
#'   that captures its data), and a closure over a compiled model that holds
#'   external pointers may not survive the trip. Stage 0 (the prior draws)
#'   runs on the host.
#'
#' @section Target:
#'   As in \code{\link{smc}}, the tempered score is
#'   \eqn{\phi = \log s(\theta) + \zeta \log L(\theta)} (system prior
#'   \eqn{s}, \code{power_posterior} exponent \eqn{\zeta}), read from the
#'   closure's \code{"posterior_parts"}; stage 0 draws the PARAMETER prior.
#'   \eqn{\lambda = 1} is therefore the target of every other sampler.
#'
#' @param log_post_fn Log-posterior closure from \code{\link{make_posterior}}
#'   (or any \code{function(theta)} returning a list with scalar
#'   \code{$loglik} and \code{$logprior}).
#' @param prior_spec Prior specification (from \code{\link{prior_spec}});
#'   used to draw stage 0 and to pre-screen proposals against the support
#'   bounds.
#' @param prior_sampler Optional \code{function()} returning one named prior
#'   draw; overrides \code{prior_spec} for stage 0.
#' @param n_particles Cloud size \eqn{N G} (default 2000 = Dynare's
#'   \code{N = 200} x \code{G = 10}). Must be a multiple of \code{n_groups}
#'   and of \code{n_strata}.
#' @param n_groups Number of Metropolis chains per stage, \code{G}
#'   (default 10).
#' @param n_stages Number of tempering stages, \code{H} (default 25).
#'   Ignored when \code{lambda_schedule} is given.
#' @param n_strata Number of density strata, \code{M} (default 20).
#' @param n_tune_steps MH steps per chain in each \eqn{c_i} tuning try,
#'   \code{K} (default 10).
#' @param thin Thinning, \code{tau} (default 10): each chain runs
#'   \code{thin} steps per kept state, and a striated jump is attempted with
#'   probability \code{0.1 / thin} per step.
#' @param accept_target Acceptance band \code{c(alpha0, alpha1)} for tuning
#'   \eqn{c_i} (default \code{c(0.2, 0.3)}).
#' @param lambda1 First rung \eqn{\lambda_1} of the tempering ladder.
#'   \code{NULL} (default) uses Dynare's \eqn{1/(10\, d\, T)} with
#'   \eqn{T} = \code{n_obs}.
#' @param n_obs Number of observations \eqn{T}, used only for the default
#'   \code{lambda1}. One of \code{lambda1}, \code{n_obs} or
#'   \code{lambda_schedule} is required.
#' @param lambda_schedule Optional explicit ladder: a strictly increasing
#'   numeric vector in \eqn{(0, 1]} ending at 1. Overrides \code{n_stages}
#'   and \code{lambda1}.
#' @param max_tune_iter Maximum \eqn{c_i} tuning tries per stage
#'   (default 50).
#' @param seed Optional integer seed. The run is then reproducible (for any
#'   \code{n_cores}) and the caller's RNG stream is restored on exit.
#' @param n_cores Number of mirai daemons for the chain groups.
#'   \code{NULL} (default) or 1 runs in this session; \code{>= 2} starts a
#'   pool of \code{min(n_cores, n_groups)} daemons for the run (see the
#'   parallel section).
#' @param verbose Print per-stage progress.
#'
#' @return A \code{\link{dynhr_chains}} object (\code{sampler = "dsmh"}) with
#'   the equally-weighted final cloud in \code{$chain} (rows ordered by kept
#'   step, \code{n_groups} rows per step; \code{$group} gives each row's
#'   chain), \code{$log_marginal_lik} (the Dynare \code{zhat} estimate of
#'   \eqn{\log \int p\, s\, L^{\zeta}}), \code{$post_logpost} (full target
#'   log density of each draw), \code{$log_liks} (raw log likelihood),
#'   \code{$log_priors}, \code{$log_phi}, \code{$log_sysprior},
#'   \code{$lambda_schedule}, \code{$ess_schedule} (importance-weight ESS
#'   per stage), \code{$c_schedule}, \code{$accept_schedule} (mutation
#'   random-walk acceptance), \code{$jump_schedule} (striated-jump
#'   acceptance), \code{$feasible_fraction}, \code{$n_eval} and
#'   \code{$elapsed_secs}.
#'
#' @references
#'   Waggoner, D. F., Wu, H., & Zha, T. (2016). Striated Metropolis-Hastings
#'     sampler for high-dimensional models. \emph{Journal of Econometrics},
#'     192(2), 406-420. \doi{10.1016/j.jeconom.2016.02.007}
#'
#'   Dynare 7 reference manual, \code{posterior_sampling_method = 'dsmh'};
#'     source \code{matlab/estimation/smc/dsmh.m}.
#' @seealso \code{\link{smc}}, \code{\link{dime}}, \code{\link{dynhr_mcmc}},
#'   \code{\link{run_posterior_estimation}} (\code{methods = "DSMH"})
#' @examples
#' ## Conjugate normal mean: y_t ~ N(mu, 1), mu ~ N(0, 2^2).
#' set.seed(1)
#' y <- rnorm(20, mean = 0.5)
#' priors <- list(mu = list(dist = "normal", p1 = 0, p2 = 2))
#' log_post <- function(theta) {
#'   lp <- dnorm(theta[["mu"]], 0, 2, log = TRUE)
#'   ll <- sum(dnorm(y, theta[["mu"]], 1, log = TRUE))
#'   list(logpost = lp + ll, loglik = ll, logprior = lp)
#' }
#' ch <- dynhr_dsmh(log_post, priors, n_particles = 400L, n_groups = 10L,
#'                  n_stages = 8L, n_strata = 10L, thin = 2L,
#'                  n_obs = length(y), seed = 1L, verbose = FALSE)
#' c(mean(ch$chain[, "mu"]), ch$log_marginal_lik)
#' @export
dynhr_dsmh <- function(log_post_fn,
                       prior_spec      = NULL,
                       prior_sampler   = NULL,
                       n_particles     = 2000L,
                       n_groups        = 10L,
                       n_stages        = 25L,
                       n_strata        = 20L,
                       n_tune_steps    = 10L,
                       thin            = 10L,
                       accept_target   = c(0.2, 0.3),
                       lambda1         = NULL,
                       n_obs           = NULL,
                       lambda_schedule = NULL,
                       max_tune_iter   = 50L,
                       seed            = NULL,
                       n_cores         = NULL,
                       verbose         = TRUE) {
  if (!is.function(log_post_fn))
    .dynhr_abort("dynhr_dsmh: `log_post_fn` must be a function.",
                 class = "dynhr_error_dsmh_args")
  if (is.null(prior_sampler) && is.null(prior_spec))
    .dynhr_abort("dynhr_dsmh: supply `prior_spec` or `prior_sampler`.",
                 class = "dynhr_error_dsmh_args")
  prob <- .dsmh_args_problem(list(
    n_particles = n_particles, n_groups = n_groups, n_strata = n_strata,
    n_tune_steps = n_tune_steps, thin = thin, max_tune_iter = max_tune_iter,
    n_cores = n_cores, accept_target = accept_target, n_stages = n_stages,
    lambda1 = lambda1, n_obs = n_obs, lambda_schedule = lambda_schedule))
  if (length(prob))
    .dynhr_abort("dynhr_dsmh: ", paste(prob, collapse = " "),
                 class = "dynhr_error_dsmh_args")
  n_particles  <- as.integer(n_particles)
  G            <- as.integer(n_groups)
  M            <- as.integer(n_strata)
  K            <- as.integer(n_tune_steps)
  tau          <- as.integer(thin)
  max_tune_iter <- as.integer(max_tune_iter)
  n_workers    <- if (is.null(n_cores)) 1L
                  else .mirai_n_cores(as.integer(n_cores), G)
  N_per <- n_particles %/% G

  .local_seed(seed)

  if (is.null(prior_sampler))
    prior_sampler <- .smc_make_prior_sampler(prior_spec)

  ll_floor <- -1e300
  t_start  <- Sys.time()

  ## ---- Interface check (as dynhr_smc): a scalar numeric $loglik ----------
  theta_chk <- prior_sampler()
  chk <- log_post_fn(theta_chk)
  if (!is.list(chk) || is.null(chk$loglik) || !is.numeric(chk$loglik) ||
      length(chk$loglik) != 1L)
    .dynhr_abort("dynhr_dsmh: `log_post_fn` must return a list with a ",
                 "scalar numeric `$loglik` (DSMH tempers the likelihood). ",
                 "Pass the make_posterior() closure, not a $logpost-only ",
                 "wrapper.", class = "dynhr_error_dsmh_args")
  par_names <- names(theta_chk)
  d <- length(theta_chk)

  ## ---- Tempering ladder ---------------------------------------------------
  if (!is.null(lambda_schedule)) {
    lam <- as.numeric(lambda_schedule)
  } else {
    H <- as.integer(n_stages)
    if (is.null(lambda1)) {
      if (is.null(n_obs))
        .dynhr_abort("dynhr_dsmh: supply `n_obs` (Dynare's first rung is ",
                     "lambda1 = 1/(10 * n_par * n_obs)), `lambda1` or ",
                     "`lambda_schedule`.", class = "dynhr_error_dsmh_args")
      lambda1 <- 1 / (10 * d * as.integer(n_obs))
    }
    ## Dynare: lambda_i = exp((H - i)/(H - 1) * log(lambda1)), i = 1..H.
    lam <- if (H == 1L) 1 else exp((H - seq_len(H)) / (H - 1) * log(lambda1))
    lam[H] <- 1
  }
  H <- length(lam)

  ## Support bounds pre-screen every chain proposal (.dsmh_eval).
  bnds <- .smc_support_bounds(prior_spec, par_names)

  ## ---- Stage 0: feasible prior draws (Dynare redraws infeasible ones) -----
  if (verbose) .dynhr_inform("DSMH: drawing ", n_particles,
                             " initial particles from the prior...")
  theta <- matrix(NA_real_, n_particles, d, dimnames = list(NULL, par_names))
  lp <- phi <- ll <- sp <- numeric(n_particles)
  n_tried  <- 0L
  n_eval   <- 1L
  max_try  <- 1000L * n_particles
  filled   <- 0L
  while (filled < n_particles) {
    if (n_tried >= max_try)
      .dynhr_abort(sprintf(paste0("dynhr_dsmh: only %d of %d prior draws ",
                                  "had a finite likelihood after %d tries."),
                           filled, n_particles, n_tried),
                   class = "dynhr_error_dsmh_infeasible")
    th <- prior_sampler()
    n_tried <- n_tried + 1L
    pp <- .smc_particle_parts(log_post_fn(th), ll_floor, missing_prior = 0)
    n_eval <- n_eval + 1L
    if (pp[["phi"]] > ll_floor && is.finite(pp[["logprior"]])) {
      filled <- filled + 1L
      theta[filled, ] <- th
      lp[filled]  <- pp[["logprior"]]
      phi[filled] <- pp[["phi"]]
      ll[filled]  <- pp[["loglik"]]
      sp[filled]  <- pp[["log_sysprior"]]
    }
  }
  feasible_fraction <- n_particles / n_tried
  log_Z <- log(feasible_fraction)

  ## ---- Chain-group backend ------------------------------------------------
  ## Chain tasks read the closure from job$lpf in-session; on the pool it is
  ## installed once per daemon (.worker_lp) and job$lpf stays NULL.
  parallel <- n_workers > 1L
  job_base <- list(par_names = par_names, bnds = bnds,
                   lpf = if (parallel) NULL else log_post_fn)
  if (parallel) {
    if (verbose)
      .dynhr_inform(sprintf("DSMH: %d chain groups on %d mirai daemons.",
                            G, n_workers))
    on.exit(mirai::daemons(NULL), add = TRUE)
    .dsmh_keep_rng(.dsmh_pool_start(n_workers, log_post_fn))
  }

  ess_tr <- c_tr <- acc_tr <- jump_tr <- numeric(H)
  tune_ok <- logical(H)
  cc <- 1
  p_jump <- 0.1 / tau
  lam_prev <- 0

  for (i in seq_len(H)) {
    dlam <- lam[i] - lam_prev

    ## (1) importance weights and the marginal-likelihood increment
    lw <- dlam * phi
    log_Z <- log_Z + .logsumexp(lw) - log(n_particles)
    w <- exp(lw - max(lw))
    w <- w / sum(w)
    ess_tr[i] <- 1 / sum(w^2)
    pool <- list(theta = theta, lp = lp, phi = phi, ll = ll, sp = sp)
    Fm <- .dsmh_prop_factor(theta, w)

    ## (2) tune c_i
    tn <- .dsmh_tune_c(cc, pool, w, lam[i], Fm, G, K,
                       accept_target, max_tune_iter, verbose, job_base,
                       parallel)
    cc <- tn$c
    n_eval <- n_eval + tn$n_eval
    tune_ok[i] <- tn$converged
    if (!tn$converged)
      .dynhr_warn(sprintf(paste0("dynhr_dsmh: stage %d: c_i tuning did not ",
                                 "reach the acceptance band [%.2f, %.2f] in ",
                                 "%d tries (last acceptance %.3f); using ",
                                 "c = %.4g."), i, accept_target[1L],
                          accept_target[2L], max_tune_iter, tn$accept, cc),
                  class = "dynhr_warning_dsmh_tune")

    ## (3) striated MH. Strata: the pool cut into M equally populated
    ## density levels of pi_{i-1}, assigned by VALUE (state and pool alike).
    t_prev_pool <- .dsmh_tempered(lp, phi, lam_prev)
    cuts <- sort(t_prev_pool)[seq_len(M - 1L) * (n_particles %/% M)]
    members <- split(seq_len(n_particles),
                     factor(.dsmh_stratum(t_prev_pool, cuts),
                            levels = seq_len(M)))

    ## The G chains are independent given the pool: one task per chain,
    ## each on its own seeded stream (identical output for any n_cores).
    idx <- .smc_systematic_resample(w, G)
    job <- c(job_base, list(
      pool = pool, members = members, cuts = cuts,
      th0 = theta[idx, , drop = FALSE], lp0 = lp[idx], phi0 = phi[idx],
      ll0 = ll[idx], sp0 = sp[idx], lam = lam[i], lam_prev = lam_prev,
      cc = cc, Fm = Fm, p_jump = p_jump, n_keep = N_per, tau = tau))
    seeds <- .dsmh_group_seeds(G)
    res <- .dsmh_run_groups(.dsmh_mutate_task, job, seeds, parallel)

    ## Chain j's r-th kept state is row (r - 1) * G + j: rows ordered by
    ## kept step, G rows per step.
    group <- rep(seq_len(G), times = N_per)
    n_rw <- n_rw_acc <- n_jp <- n_jp_acc <- 0L
    for (j in seq_len(G)) {
      rows <- (seq_len(N_per) - 1L) * G + j
      rj <- res[[j]]
      theta[rows, ] <- rj$theta
      lp[rows] <- rj$lp; phi[rows] <- rj$phi; ll[rows] <- rj$ll
      sp[rows] <- rj$sp
      n_rw     <- n_rw + rj$n_rw
      n_rw_acc <- n_rw_acc + rj$n_rw_acc
      n_jp     <- n_jp + rj$n_jp
      n_jp_acc <- n_jp_acc + rj$n_jp_acc
      n_eval   <- n_eval + rj$n_eval
    }
    c_tr[i]    <- cc
    acc_tr[i]  <- if (n_rw > 0L) n_rw_acc / n_rw else NA_real_
    jump_tr[i] <- if (n_jp > 0L) n_jp_acc / n_jp else NA_real_
    lam_prev   <- lam[i]
    if (verbose)
      .dynhr_inform(sprintf(paste0("DSMH stage %d/%d: lambda=%.4g  ESS=%.0f/%d",
                                   "  c=%.4g  accept=%.3f  jump=%.3f  ",
                                   "log_mlik=%.3f"),
                            i, H, lam[i], ess_tr[i], n_particles, cc,
                            acc_tr[i], jump_tr[i], log_Z))
  }

  logpost <- .dsmh_tempered(lp, phi, 1)
  out <- list(
    chain             = theta,
    particles         = theta,
    group             = group,
    log_liks          = ll,
    log_priors        = lp,
    log_phi           = phi,
    log_sysprior      = sp,
    post_logpost      = logpost,
    logpost_trace     = logpost,
    log_marginal_lik  = log_Z,
    lambda_schedule   = lam,
    ess_schedule      = ess_tr,
    c_schedule        = c_tr,
    accept_schedule   = acc_tr,
    jump_schedule     = jump_tr,
    tune_converged    = tune_ok,
    feasible_fraction = feasible_fraction,
    n_stages          = H,
    n_particles       = n_particles,
    n_groups          = G,
    n_burn            = 0L,
    n_eval            = n_eval,
    elapsed_secs      = as.numeric(difftime(Sys.time(), t_start,
                                            units = "secs")),
    sampler           = "dsmh"
  )
  if (!is.null(bnds$lower)) {
    out$support_lower <- bnds$lower
    out$support_upper <- bnds$upper
  }
  new_dynhr_chains(out, "dsmh")
}
