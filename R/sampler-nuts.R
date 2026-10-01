## R/sampler-nuts.R
## --------------------------------------------------------------------------
## Phase-2 split from nuts-monolith.R.
##
## Internal: .nuts_build_tree()    -- recursive tree builder (Betancourt 2017
##                                    multinomial + generalized U-turn)
##           .nuts_warmup_windows() -- Stan-style windowed warmup schedule
## Sampler:  dynhr_nuts()          -- NUTS with dual averaging +
##                                    windowed mass-matrix adaptation
##
## Depends on HMC internals in sampler-hmc.R (loaded together in the package).
## --------------------------------------------------------------------------

# NUTS internals
# ============================================================================

#' Recursive tree building for NUTS (Betancourt 2017 multinomial sampling)
#'
#' Implements multinomial trajectory sampling from Betancourt (2017) "A
#' Conceptual Introduction to Hamiltonian Monte Carlo" (arXiv:1701.02434).
#' Key differences from Hoffman-Gelman 2014 (HG14) slice sampler:
#'
#' \itemize{
#'   \item Base case returns \code{log_weight = joint1} (the joint log density
#'     of the proposed leaf) rather than a slice-validity count. The slice
#'     variable \code{log_u} is no longer passed to \code{.nuts_build_tree} --
#'     divergence is the only energy-based stopping criterion.
#'   \item Recursive case combines subtrees via log-sum-exp and applies
#'     \emph{biased progressive multinomial sampling}: the new subtree's
#'     proposal replaces the current one with probability
#'     \code{exp(tree2$log_weight - log_sum_exp(tree$log_weight, tree2$log_weight))}.
#'   \item Generalized U-turn criterion: stopping fires if EITHER subtree
#'     alone would stop (recursive) OR the combined tree boundary would
#'     (original HG14 velocity-metric check). Because we short-circuit when
#'     the left subtree stops, the right subtree is never built in that case;
#'     log_weight then covers only the left subtree (benign -- the trajectory
#'     would have been rejected anyway).
#' }
#'
#' @param theta Position
#' @param r Momentum
#' @param v Direction (-1 or +1)
#' @param j Tree depth
#' @param eps Step size
#' @param lp_fn Scalar log-posterior function
#' @param grad_fn Gradient function
#' @param M_inv_diag Inverse mass diagonal (numeric vector; diagonal path)
#' @param M_inv Dense inverse mass matrix (d×d; dense path; NULL = diagonal)
#' @param joint0 Initial joint log density (for acceptance stats and
#'   divergence detection: \code{joint0 - joint1 > delta_max} is divergent)
#' @param delta_max Maximum energy error before flagging divergence (default 1000)
#' @param g The gradient at \code{theta}, carried from the evaluation that
#'   produced it (fused path; the separate-call path too); NULL:
#'   computed with \code{grad_fn}.
#' @param vg_fn Fused path: \code{function(theta) -> list(lp, grad)}
#'   (\code{.hmc_fused_target()}); NULL = separate \code{lp_fn} /
#'   \code{grad_fn} calls.
#'
#' @return list with theta_minus, r_minus, theta_plus, r_plus,
#'         theta_prime, log_weight, stop, sum_alpha, n_leaves, divergent,
#'         g_minus, g_plus (gradients at the two ends) and lp_prime, g_prime
#'         (log density and gradient at theta_prime)
#' @noRd
.nuts_build_tree <- function(theta, r, v, j, eps,
                              lp_fn, grad_fn, M_inv_diag, joint0,
                              delta_max = 1000, M_inv = NULL,
                              g = NULL, vg_fn = NULL) {
  par_names <- names(theta)

  if (j == 0L) {
    # --- Base case: single leapfrog step ---
    step <- .hmc_leapfrog(theta, r, v * eps, grad_fn, M_inv_diag, M_inv = M_inv,
                          g0 = g, vg_fn = vg_fn)

    if (is.null(step) || any(!is.finite(step$theta)) || any(!is.finite(step$r))) {
      # Divergent step: log_weight = -Inf (zero probability mass)
      return(list(
        theta_minus = theta, r_minus = r,
        theta_plus  = theta, r_plus  = r,
        theta_prime = theta,
        log_weight = -Inf, stop = TRUE,
        sum_alpha = 0, n_leaves = 1L,
        divergent = TRUE,
        g_minus = g, g_plus = g, lp_prime = NULL, g_prime = g
      ))
    }

    theta1 <- step$theta
    r1     <- step$r
    names(theta1) <- par_names

    ## Fused: the leapfrog's end-point evaluation already carries the value.
    ## Either path: its gradient is carried to the next step from here.
    lp1    <- if (is.null(vg_fn)) lp_fn(theta1) else step$lp
    g1     <- step$g
    joint1 <- lp1 - .hmc_kinetic(r1, M_inv_diag, M_inv = M_inv)

    # Multinomial log-weight: unnormalised log probability of this leaf.
    # No slice-validity filter; divergence is handled below.
    log_weight <- if (is.finite(joint1)) joint1 else -Inf

    # Stop criterion: energy excursion too large (divergence)
    s_ok <- is.finite(joint1) && (joint0 - joint1 < delta_max)

    # Acceptance statistic for dual averaging (unchanged from HG14)
    alpha1 <- min(1, exp(joint1 - joint0))
    if (!is.finite(alpha1)) alpha1 <- 0

    list(
      theta_minus = theta1, r_minus = r1,
      theta_plus  = theta1, r_plus  = r1,
      theta_prime = theta1,
      log_weight  = log_weight,
      stop        = !s_ok,
      sum_alpha   = alpha1, n_leaves = 1L,
      divergent   = (joint0 - joint1 > delta_max),
      g_minus = g1, g_plus = g1, lp_prime = lp1, g_prime = g1
    )
  } else {
    # --- Recursion ---
    tree <- .nuts_build_tree(theta, r, v, j - 1L, eps,
                              lp_fn, grad_fn, M_inv_diag, joint0, delta_max,
                              M_inv = M_inv, g = g, vg_fn = vg_fn)
    # Early exit: if left subtree already stops, combined tree stops too.
    # The right subtree is NOT built (log_weight covers only the left
    # subtree, which is correct: the trajectory would be rejected anyway).
    if (tree$stop) return(tree)

    if (v == -1) {
      tree2 <- .nuts_build_tree(tree$theta_minus, tree$r_minus, v,
                                 j - 1L, eps, lp_fn, grad_fn, M_inv_diag,
                                 joint0, delta_max, M_inv = M_inv,
                                 g = tree$g_minus, vg_fn = vg_fn)
      theta_minus <- tree2$theta_minus
      r_minus     <- tree2$r_minus
      g_minus     <- tree2$g_minus
      theta_plus  <- tree$theta_plus
      r_plus      <- tree$r_plus
      g_plus      <- tree$g_plus
    } else {
      tree2 <- .nuts_build_tree(tree$theta_plus, tree$r_plus, v,
                                 j - 1L, eps, lp_fn, grad_fn, M_inv_diag,
                                 joint0, delta_max, M_inv = M_inv,
                                 g = tree$g_plus, vg_fn = vg_fn)
      theta_minus <- tree$theta_minus
      r_minus     <- tree$r_minus
      g_minus     <- tree$g_minus
      theta_plus  <- tree2$theta_plus
      r_plus      <- tree2$r_plus
      g_plus      <- tree2$g_plus
    }

    # Multinomial biased progressive sampling (Betancourt 2017, Alg 6):
    # combine log-weights via log-sum-exp and accept the new subtree's
    # proposal with probability exp(tree2$log_weight - log_w_total).
    log_w_total <- if (is.finite(tree$log_weight) && is.finite(tree2$log_weight)) {
      w_max <- max(tree$log_weight, tree2$log_weight)
      w_max + log(exp(tree$log_weight - w_max) + exp(tree2$log_weight - w_max))
    } else if (is.finite(tree2$log_weight)) {
      tree2$log_weight
    } else {
      tree$log_weight   # -Inf if both are -Inf
    }

    theta_prime <- tree$theta_prime
    lp_prime    <- tree$lp_prime
    g_prime     <- tree$g_prime
    if (is.finite(tree2$log_weight) && is.finite(log_w_total)) {
      accept_sub <- exp(tree2$log_weight - log_w_total)
      if (runif(1) < accept_sub) {
        theta_prime <- tree2$theta_prime
        lp_prime    <- tree2$lp_prime
        g_prime     <- tree2$g_prime
      }
    }

    # Generalized U-turn criterion: stop if EITHER subtree stops OR the
    # combined tree boundary shows a U-turn (velocity-metric, see below).
    # The velocity-metric dot product (M^-1 r, not raw r) is critical:
    # with a non-identity mass matrix the raw-momentum test almost never
    # fires on a strongly anisotropic posterior, causing trees to grow to
    # max depth.
    dtheta <- theta_plus - theta_minus
    if (is.null(M_inv)) {
      ## Diagonal path (bit-identical to pre-dense behaviour)
      u_turn_full <- (sum(dtheta * (M_inv_diag * r_minus)) < 0) ||
                     (sum(dtheta * (M_inv_diag * r_plus))  < 0)
    } else {
      ## Dense path: velocity = M⁻¹ r
      u_turn_full <- (sum(dtheta * .metric_apply_inv(M_inv, r_minus)) < 0) ||
                     (sum(dtheta * .metric_apply_inv(M_inv, r_plus))  < 0)
    }
    stop_flag <- tree2$stop || u_turn_full

    list(
      theta_minus = theta_minus, r_minus = r_minus,
      theta_plus  = theta_plus,  r_plus  = r_plus,
      theta_prime = theta_prime,
      log_weight  = log_w_total,
      stop        = stop_flag,
      sum_alpha   = tree$sum_alpha + tree2$sum_alpha,
      n_leaves    = tree$n_leaves + tree2$n_leaves,
      divergent   = tree$divergent || tree2$divergent,
      g_minus = g_minus, g_plus = g_plus, lp_prime = lp_prime, g_prime = g_prime
    )
  }
}


# ============================================================================
# .nuts_warmup_windows() -- Stan-style windowed warmup schedule
# ============================================================================

#' Compute windowed warmup schedule for NUTS mass-matrix adaptation
#'
#' Produces a three-phase warmup schedule modelled on Stan's windowed
#' adaptation (Stan reference manual, sec. "HMC algorithm"):
#' \itemize{
#'   \item \strong{Fast-init} \code{[1, init_buffer]}: step size only (no
#'     mass-matrix update). Default: \code{min(75, floor(0.15 * n_warmup))}.
#'   \item \strong{Slow windows} \code{[init_buffer+1, n_warmup-term_buffer]}:
#'     a sequence of expanding windows whose widths double from
#'     \code{base_window}. At the end of each window the diagonal mass matrix
#'     is updated from a memoryless (within-window) variance estimate, the
#'     step size is re-found, and the dual-averaging state is reset (this
#'     must happen at EVERY window boundary, not just the first, so the
#'     averager does not accumulate stats from the wrong mass geometry).
#'   \item \strong{Fast-final} \code{[n_warmup-term_buffer+1, n_warmup]}:
#'     step-size only. Default: \code{min(50, floor(0.10 * n_warmup))}.
#' }
#'
#' For small \code{n_warmup} the buffers are shrunk proportionally so that the
#' three phases always fit; if the slow window budget is < \code{base_window}
#' a single combined-window fallback is used.
#'
#' @param n_warmup   Number of warmup iterations.
#' @param init_buffer  Fast-init window width (NULL = auto).
#' @param base_window  Minimum slow-window width (default 25L).
#' @param term_buffer  Fast-final window width (NULL = auto).
#'
#' @return \code{data.frame(start, end, type)} where \code{type} is
#'   \code{"fast_init"}, \code{"slow"}, or \code{"fast_final"}. The
#'   \code{"slow"} rows correspond to individual expanding windows; the
#'   mass-matrix update fires at each \code{end} of a slow row. Returns a
#'   single row of type \code{"fast_init"} covering all of \code{n_warmup}
#'   when \code{n_warmup < 10}.
#' @noRd
.nuts_warmup_windows <- function(n_warmup, init_buffer = NULL,
                                  base_window = 25L, term_buffer = NULL) {
  n_warmup    <- as.integer(n_warmup)
  base_window <- as.integer(base_window)

  if (n_warmup < 10L) {
    # Degenerate case: skip all adaptation, single fast-init window
    return(data.frame(start = 1L, end = n_warmup, type = "fast_init",
                      stringsAsFactors = FALSE))
  }

  # Auto-size the fast buffers (proportional shrink for small n_warmup)
  if (is.null(init_buffer)) {
    init_buffer <- as.integer(min(75L, max(1L, floor(0.15 * n_warmup))))
  } else {
    init_buffer <- as.integer(init_buffer)
  }
  if (is.null(term_buffer)) {
    term_buffer <- as.integer(min(50L, max(1L, floor(0.10 * n_warmup))))
  } else {
    term_buffer <- as.integer(term_buffer)
  }

  # Ensure the three phases actually fit; shrink buffers if necessary
  slack <- n_warmup - init_buffer - term_buffer
  if (slack < base_window) {
    # Compress proportionally, leaving at least 1 iteration each
    total_buf <- init_buffer + term_buffer
    if (total_buf > 0) {
      init_buffer <- max(1L, as.integer(floor(init_buffer / total_buf * (n_warmup - 1L))))
      term_buffer <- max(1L, n_warmup - 1L - init_buffer)
    }
    slack <- n_warmup - init_buffer - term_buffer
  }

  rows <- list()

  # Phase 1: fast-init
  rows[[length(rows) + 1L]] <- data.frame(
    start = 1L, end = init_buffer, type = "fast_init",
    stringsAsFactors = FALSE
  )

  # Phase 2: slow windows with doubling widths
  slow_start <- init_buffer + 1L
  slow_end   <- n_warmup - term_buffer

  if (slow_end >= slow_start) {
    w <- base_window
    cur <- slow_start
    while (cur <= slow_end) {
      win_end <- min(cur + w - 1L, slow_end)
      rows[[length(rows) + 1L]] <- data.frame(
        start = cur, end = win_end, type = "slow",
        stringsAsFactors = FALSE
      )
      cur <- win_end + 1L
      w   <- w * 2L
    }
  }

  # Phase 3: fast-final
  if (term_buffer > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      start = n_warmup - term_buffer + 1L, end = n_warmup, type = "fast_final",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}


# ============================================================================
# .nuts_transition() -- one NUTS transition at a fixed step size and metric
# ============================================================================

#' One NUTS transition (momentum refresh + multinomial doubling tree)
#'
#' The loop body of \code{dynhr_nuts()}, factored out so the pooled
#' multi-chain warmup (\code{.nuts_pooled_warmup()}, sampler-nuts-pooled.R)
#' runs the identical kernel. RNG consumption order is unchanged: one momentum
#' draw, then per doubling one direction draw and the tree's multinomial draws.
#'
#' @param theta Current position (sampler space), named.
#' @param lp_curr Log density at \code{theta}.
#' @param eps Step size.
#' @param lp_fn,grad_fn Scalar log density and its gradient.
#' @param M_diag,M_inv_diag Diagonal mass / inverse mass (diagonal path).
#' @param M_inv,chol_M Dense or low-rank metric (NULL = diagonal path).
#' @param max_treedepth,delta_max As in \code{dynhr_nuts()}.
#' @param g_curr The gradient at \code{theta} when already known (NULL: one
#'   \code{grad_fn} call here); carried along the trajectory on both paths
#'.
#' @param vg_fn Fused path: the value-and-gradient function
#'   (\code{.hmc_fused_target()}); NULL = separate \code{lp_fn} /
#'   \code{grad_fn} calls.
#' @return list(theta, joint0, depth, divergent, alpha_sum, n_alpha, n_leaves,
#'   lp, g) -- \code{theta} is the selected proposal (unnamed-safe: callers
#'   re-apply names), \code{alpha_sum / n_alpha} the trajectory-averaged
#'   acceptance statistic, \code{n_leaves} the new-point evaluations spent;
#'   \code{g} is the gradient at \code{theta}, and on the fused path
#'   \code{lp} its log density (NULL otherwise).
#' @noRd
.nuts_transition <- function(theta, lp_curr, eps, lp_fn, grad_fn,
                             M_diag, M_inv_diag, M_inv, chol_M,
                             max_treedepth, delta_max,
                             g_curr = NULL, vg_fn = NULL) {
  d <- length(theta)
  # --- Sample momentum ---
  r0 <- .hmc_sample_momentum(d, M_diag, chol_M = chol_M)
  ## The gradient at the start point, once for both directions (the
  ## separate-call path used to recompute it at the first leaf of each).
  if (is.null(g_curr) && max_treedepth > 0L) g_curr <- grad_fn(theta)

  # --- Joint log density (multinomial: no slice variable needed) ---
  joint0 <- lp_curr - .hmc_kinetic(r0, M_inv_diag, M_inv = M_inv)

  # --- Initialize tree ---
  # Multinomial: initial point has log_weight = joint0 (the reference mass)
  theta_minus <- theta
  theta_plus  <- theta
  r_minus     <- r0
  r_plus      <- r0
  theta_m     <- theta          # current proposal
  lp_m        <- lp_curr        # ... its log density and gradient
  g_m         <- g_curr
  g_minus     <- g_curr         # gradients at the two trajectory ends
  g_plus      <- g_curr
  log_w_total <- joint0         # cumulative log-weight of all leaves so far
  j <- 0L
  stop <- FALSE
  any_divergent <- FALSE
  n_leaves <- 0L
  # Dual-averaging acceptance statistic, accumulated over the WHOLE trajectory
  # (Hoffman & Gelman 2014, Alg. 6 -- alpha / n_alpha). Using only the final
  # doubling makes the statistic swing between 0 and 1 iteration-to-iteration
  # and destabilises step-size adaptation.
  alpha_sum <- 0
  n_alpha   <- 0L

  while (!stop && j < max_treedepth) {
    # Choose direction
    v <- sample(c(-1L, 1L), 1)

    if (v == -1L) {
      tree <- .nuts_build_tree(theta_minus, r_minus, v, j, eps,
                                lp_fn, grad_fn, M_inv_diag, joint0, delta_max,
                                M_inv = M_inv, g = g_minus, vg_fn = vg_fn)
      theta_minus <- tree$theta_minus
      r_minus     <- tree$r_minus
      g_minus     <- tree$g_minus
    } else {
      tree <- .nuts_build_tree(theta_plus, r_plus, v, j, eps,
                                lp_fn, grad_fn, M_inv_diag, joint0, delta_max,
                                M_inv = M_inv, g = g_plus, vg_fn = vg_fn)
      theta_plus <- tree$theta_plus
      r_plus     <- tree$r_plus
      g_plus     <- tree$g_plus
    }

    # BIASED progressive multinomial (Betancourt 2017, Stan): accept the new
    # subtree's proposal with probability min(1, W_new / W_old) -- the new
    # subtree's weight relative to the trajectory weight BEFORE merging.
    # This deliberately over-weights the fresh subtree (anti-correlated
    # exploration); it remains an exact sampler because the within-subtree
    # selection (in .nuts_build_tree) is unbiased multinomial. (The initial
    # point's weight joint0 is already in log_w_total at the first doubling.)
    if (!tree$stop && is.finite(tree$log_weight)) {
      accept_prob <- min(1, exp(tree$log_weight - log_w_total))
      if (runif(1) < accept_prob) {
        theta_m <- tree$theta_prime
        lp_m    <- tree$lp_prime
        g_m     <- tree$g_prime
      }
      log_w_total <- if (is.finite(log_w_total)) {
        w_max <- max(log_w_total, tree$log_weight)
        w_max + log(exp(log_w_total - w_max) + exp(tree$log_weight - w_max))
      } else {
        tree$log_weight
      }
    }

    alpha_sum <- alpha_sum + tree$sum_alpha
    n_alpha   <- n_alpha + tree$n_leaves
    if (tree$divergent) any_divergent <- TRUE

    # Outer U-turn check on the full tree boundary (velocity metric).
    # The generalized check is already embedded in .nuts_build_tree
    # (recursive subtrees); this outer check guards the full doublings.
    dtheta  <- theta_plus - theta_minus
    if (is.null(M_inv)) {
      ## Diagonal path (bit-identical to pre-dense behaviour)
      v_minus <- M_inv_diag * r_minus
      v_plus  <- M_inv_diag * r_plus
    } else {
      ## Dense path: velocity = M⁻¹ r
      v_minus <- .metric_apply_inv(M_inv, r_minus)
      v_plus  <- .metric_apply_inv(M_inv, r_plus)
    }
    stop <- tree$stop ||
      (sum(dtheta * v_minus) < 0) ||
      (sum(dtheta * v_plus) < 0)

    n_leaves <- n_leaves + tree$n_leaves

    j <- j + 1L
  }

  list(theta = theta_m, joint0 = joint0, depth = j, divergent = any_divergent,
       alpha_sum = alpha_sum, n_alpha = n_alpha, n_leaves = n_leaves,
       lp = if (is.null(vg_fn)) NULL else lp_m,
       g  = g_m)
}


# ============================================================================
# dynhr_nuts() -- NUTS with dual averaging + windowed mass adaptation
# ============================================================================

## --------------------------------------------------------------------------
## NUTS argument checker
## --------------------------------------------------------------------------

#' Problems with the arguments of dynhr_nuts()
#'
#' Rules and why: the tree depth, the warmup window sizes and the draw counts
#' are counts (a depth of 0 never moves the chain; a zero-width slow window
#' never advances the warmup schedule); the step size and the divergence
#' threshold are positive; the target acceptance is a probability strictly
#' inside (0, 1) (dual averaging drives the step size to 0 or infinity at the
#' ends); the mass diagonal / dense metric are positive / positive definite
#' and match the parameter count.
#'
#' @param args Named list of the supplied dynhr_nuts() arguments.
#' @param n_par Number of estimated parameters, or NULL when not known.
#' @return Character vector of problems; \code{character(0)} when fine.
#' @noRd
.nuts_args_problem <- function(args, n_par = NULL) {
  lr_names <- c("cutoff", "max_rank", "gamma")
  rules <- c(list(
    n_draws       = .mcmc_r_whole(1L),
    n_warmup      = .mcmc_r_whole(0L),
    step_size     = .mcmc_r_pos(),
    max_treedepth = .mcmc_r_whole(1L),
    target_accept = .mcmc_r_range(0, 1, TRUE, TRUE),
    adapt_mass    = .mcmc_r_flag(),
    init_buffer   = .mcmc_r_whole(0L),
    base_window   = .mcmc_r_whole(1L),
    term_buffer   = .mcmc_r_whole(0L),
    delta_max     = .mcmc_r_pos(allow_inf = TRUE),
    metric        = .mcmc_r_choice(c("diagonal", "warmup_dense", "lowrank", "fisher_diag")),
    lowrank_control = function(v, nm, n_par) {
      if (is.list(v) && (length(v) == 0L || !is.null(names(v))) &&
          all(names(v) %in% lr_names)) return(NULL)
      paste0("`lowrank_control` must be a named list with elements in: ",
             paste(lr_names, collapse = ", "), ".")
    }),
    .mcmc_metric_rules(dynhr_nuts))
  .mcmc_check_args(args, dynhr_nuts, "dynhr_nuts", rules, n_par)
}


#' @param log_post_fn function(theta) -> list(logpost, loglik, logprior)
#' @param theta_init Named numeric vector of starting parameter values
#' @param n_draws Post-warmup draws to retain
#' @param n_warmup Warmup iterations (discarded; used for adaptation)
#' @param step_size Initial step size (NULL = auto-find)
#' @param max_treedepth Maximum tree depth (2^max_treedepth leapfrog steps).
#'   Default 8. With the numerical-gradient log-posterior each leapfrog step
#'   costs ~2*d evaluations, so depth 10 (1023 steps) is a wall-clock cliff;
#'   8 (255 steps) caps the worst-case iteration cost while leaving normal
#'   depth-3-5 trajectories untouched. Raise it for cheap/analytic gradients.
#' @param target_accept Target acceptance rate for dual averaging (0.6-0.95)
#' @param adapt_mass Adapt diagonal mass matrix during warmup using Stan-style
#'   windowed adaptation. When TRUE (default), the warmup is split into
#'   three phases: fast-init (step size only), expanding slow windows (mass
#'   update from within-window variance + step-size re-find + dual-averaging
#'   reset at each window boundary), and fast-final (step size only).
#'   Convention (Stan's): the momentum is \code{p ~ N(0, M)} and the adapted
#'   quantity is the INVERSE mass \code{M^{-1}}, set to the within-window
#'   posterior variance (\code{metric = "diagonal"}), so each coordinate is
#'   rescaled to unit posterior scale. Before 0.9.3.50 the mass itself was set
#'   to the variance (inverted), which squared the target's conditioning
#'   instead of removing it; results remained valid at convergence but finite
#'   runs could be badly under-mixed. \code{mass_diag} is the MASS diagonal
#'   \code{M} (i.e. \code{1 / variance}), and the returned \code{mass_matrix}
#'   is the diagonal mass \code{M} on the diagonal paths and the dense
#'   INVERSE mass \code{M^{-1}} on the dense / low-rank paths.
#' @param init_buffer Fast-init window width (NULL = auto: min(75, 15% of
#'   n_warmup)).
#' @param base_window Minimum slow-window width for windowed mass adaptation
#'   (default 25L); windows double in width.
#' @param term_buffer Fast-final window width (NULL = auto: min(50, 10% of
#'   n_warmup)).
#' @param grad_fn Optional analytical gradient function(theta) -> numeric.
#'   When it carries a fused log-posterior (\code{attr(grad_fn,
#'   "logpost_grad")}, as \code{make_posterior_grad()}'s Gaussian closures do)
#'   that equals \code{log_post_fn} at \code{theta_init} to 1e-10 relative,
#'   every new leaf costs ONE value-and-gradient evaluation (gradients carried
#'   along the trajectory; \code{.hmc_fused_target()}) and \code{log_post_fn}
#'   is not called after that check.
#' @param grad_method "simple" or "Richardson" for numerical gradients
#' @param delta_max Maximum energy error before marking divergence
#' @param metric Mass-matrix adaptation: \code{"diagonal"} (default; inverse
#'   mass = within-window variance, Stan's rule),
#'   \code{"warmup_dense"} (Ledoit-Wolf dense at the final slow window),
#'   \code{"fisher_diag"} (opt-in; Seyboldt, Carlson & Carpenter 2026,
#'   arXiv:2603.18845, Theorem 2.2): same Stan window schedule as
#'   \code{"diagonal"}, but at each slow-window end the inverse mass is
#'   \code{sigma^2 = sqrt(var(x) / var(g))} from the window's draws x AND
#'   log-density gradients g (\code{.fisher_diag_estimate()}), the diagonal
#'   minimiser of the sample Fisher divergence to N(0, I). Costs one extra
#'   gradient per warmup iteration inside a slow window (counted in
#'   \code{n_grad_evals}); sampling is unchanged. Or
#'   \code{"lowrank"} (opt-in; Lao 2026, arXiv:2607.23788): starts diagonal
#'   and at every window endpoint re-estimates a low-rank-plus-diagonal
#'   inverse mass \code{G = diag(sd) (I + U (Lambda - I) U') diag(sd)} from
#'   the window's draws AND log-density gradients (the Fisher-divergence
#'   estimator of Seyboldt, Carlson & Carpenter 2026, arXiv:2603.18845,
#'   Algorithm 1; see \code{.lowrank_estimate()}), promoting from diagonal
#'   (rank 0) to rank
#'   \code{k >= 1} when eigenvalues clear the threshold. Uses its own
#'   schedule (\code{.lowrank_warmup_windows()}: 1-step init, first window
#'   \code{8 (k_cap + 1)}, 1.5x growth, final 15\% step size only), so
#'   \code{init_buffer}/\code{base_window}/\code{term_buffer} are ignored.
#'   Costs one extra gradient per warmup iteration inside a slow window
#'   (counted in \code{n_grad_evals}). The result gains a \code{lowrank}
#'   element (final rank, per-window history) and \code{mass_matrix} is the
#'   dense inverse mass \code{G}. Ignored when a dense \code{M_inv} is supplied.
#' @param lowrank_control Optional list for \code{metric = "lowrank"}:
#'   \code{cutoff} (eigenvalue-ratio promotion threshold, default 2),
#'   \code{max_rank} (rank cap below \code{k_cap = min(50, floor(d/2))};
#'   0 = Fisher diagonal only), \code{gamma} (ridge, default 1e-5).
#' @param verbose Print progress
#' @param progressor progressr callback or NULL
#' @param chain_id Label for progress messages
#' @param transform Optional "dynhr_param_transform" object (from
#'   \code{\link{build_param_transform}}). When non-NULL (opt-in), NUTS runs
#'   in UNCONSTRAINED eta-space: `theta_init` is mapped to
#'   `eta0 = transform$to_unconstrained(theta_init)`, the target is
#'   `make_transformed_logpost(log_post_fn, transform, include_jacobian =
#'   TRUE)`, and `grad_fn` (if supplied) is wrapped via
#'   `make_transformed_grad()`; if `grad_fn` is NULL the numerical gradient
#'   differentiates the wrapped (Jacobian-included) eta-space target
#'   directly. `mass_diag` is interpreted in eta-space (the orchestrator
#'   converts a theta-space inverse-Hessian diagonal via the delta method).
#'   The returned `chain`/`full_chain` are mapped back to theta-space, and
#'   `logpost_trace`/`post_logpost`/`energy_trace` store THETA-SPACE
#'   log-posterior (Jacobian subtracted back out). Divergence/acceptance
#'   bookkeeping is unchanged (computed in eta-space, where the sampler
#'   actually moves). When NULL (default), behaviour is bit-identical to
#'   before.
#' @param checkpoint Optional list for memory-streamed restartable
#'   checkpointing (opt-in; NULL = disabled). When supplied, draws are
#'   streamed to per-chain binary files in \code{flush_every}-row chunks so
#'   in-RAM draw storage is bounded by \code{flush_every * d} instead of
#'   \code{n_total * d}. A restart state (position, log-posterior, frozen
#'   step size, frozen mass matrix, RNG, draw count) is saved after every
#'   flush.  Resume (\code{resume = TRUE}) continues from the saved state
#'   POST-WARMUP: warmup is not replayed; the frozen adaptation (step size +
#'   mass matrix) is restored and new sampling draws are appended.  Fields:
#'   \describe{
#'     \item{dir}{Directory for checkpoint files (must exist or be creatable).}
#'     \item{flush_every}{Rows per flush (default 1000L).}
#'     \item{resume}{Logical; TRUE to continue a prior run.}
#'     \item{fingerprint}{Config fingerprint list for mismatch detection.}
#'     \item{write_meta}{Logical; FALSE suppresses meta.rds write (parallel
#'       orchestrator writes it once to avoid races; default TRUE).}
#'     \item{return_chain}{Logical; FALSE skips the end read-back (draws
#'       remain on disk only; default TRUE).}
#'   }
#'
#' @return List compatible with rwmh() output + NUTS-specific diagnostics.
#'   \code{acceptance_rate} and \code{accept_stat} are the mean post-warmup
#'   NUTS acceptance statistic (the trajectory-averaged Metropolis alpha that
#'   \code{target_accept} targets; Stan's mean \code{accept_stat__}),
#'   \code{accept_stats} its per-iteration values, and \code{move_rate} the
#'   fraction of post-warmup iterations whose draw differs from the previous
#'   one.
#' @noRd
dynhr_nuts <- function(
    log_post_fn,
    theta_init,
    n_draws       = 2000L,
    n_warmup      = 1000L,
    step_size     = NULL,
    max_treedepth = 8L,
    target_accept = 0.80,
    adapt_mass    = TRUE,
    init_buffer   = NULL,
    base_window   = 25L,
    term_buffer   = NULL,
    grad_fn       = NULL,
    grad_method   = "forward",
    delta_max     = 1000,
    mass_diag     = NULL,
    metric        = c("diagonal", "warmup_dense", "lowrank", "fisher_diag"),
    verbose       = TRUE,
    progressor    = NULL,
    chain_id      = NULL,
    transform     = NULL,
    M_inv         = NULL,
    chol_M        = NULL,
    checkpoint    = NULL,
    lowrank_control = NULL
) {
  stopifnot(is.function(log_post_fn), is.numeric(theta_init))
  .mcmc_abort_if_problems("dynhr_nuts", .nuts_args_problem(
    list(n_draws = n_draws, n_warmup = n_warmup, step_size = step_size,
         max_treedepth = max_treedepth, target_accept = target_accept,
         adapt_mass = adapt_mass, init_buffer = init_buffer,
         base_window = base_window, term_buffer = term_buffer,
         grad_fn = grad_fn, delta_max = delta_max, mass_diag = mass_diag,
         metric = metric, M_inv = M_inv, chol_M = chol_M,
         lowrank_control = lowrank_control),
    n_par = length(theta_init)))
  metric <- match.arg(metric)
  lr_ctrl <- list(cutoff = 2, max_rank = NULL, gamma = 1e-5)
  if (!is.null(lowrank_control))
    lr_ctrl[names(lowrank_control)] <- lowrank_control
  d <- length(theta_init)
  par_names <- names(theta_init)
  n_total <- n_draws + n_warmup

  # ---- Checkpoint / streaming (opt-in). When `checkpoint` is a list carrying
  # a `dir`, draws (theta-space chain) + logpost are streamed to per-chain
  # files in flush_every-row chunks (RAM bounded by flush_every * d, not
  # n_total * d). A restart state is saved after every flush. Resume
  # (`checkpoint$resume = TRUE`) loads the frozen post-warmup adaptation
  # (step size, mass matrix) and RNG state and continues adding sampling
  # draws -- warmup is NOT replayed on resume.
  ckpt        <- !is.null(checkpoint)
  ckpt_resume <- ckpt && isTRUE(checkpoint$resume)
  flush_every <- if (ckpt) as.integer(checkpoint$flush_every %||% 1000L) else NA_integer_
  ckpt_paths  <- if (ckpt) .ckpt_paths(checkpoint$dir, chain_id) else NULL

  # ---- Transformed (opt-in) target: operate on eta = to_unconstrained(theta)
  if (!is.null(transform)) {
    target_fn  <- make_transformed_logpost(log_post_fn, transform, include_jacobian = TRUE)
    state_init <- transform$to_unconstrained(theta_init)
    names(state_init) <- par_names
  } else {
    target_fn  <- log_post_fn
    state_init <- theta_init
  }

  # --- Scalar log-posterior wrapper (operates on eta when transformed) ---
  .lp_scalar <- function(theta) {
    names(theta) <- par_names
    res <- target_fn(theta)
    val <- if (is.list(res)) res$logpost else res
    if (!is.finite(val)) -1e300 else val
  }

  # --- Gradient function ---
  # kernel_stats: make_posterior_grad()'s "implicit"/"adjoint"/
  # "adjoint_solution" closures attach fallback-usage counters as
  # attr(grad_fn, "kernel_stats") -- captured from the caller-supplied
  # grad_fn BEFORE any transform-wrapping (which would drop the attribute)
  # and surfaced in the returned diagnostics list below.
  grad_fn_kernel_stats <- attr(grad_fn, "kernel_stats")
  if (is.null(grad_fn)) {
    .grad <- function(theta) .hmc_gradient(.lp_scalar, theta, method = grad_method)
  } else if (!is.null(transform)) {
    .grad <- make_transformed_grad(grad_fn, transform)
  } else {
    .grad <- grad_fn
  }
  # --- Fused value + gradient: one evaluation per new leaf, gradients
  # carried along the trajectory; NULL = separate calls exactly as before.
  # (A checkpoint resume keeps its saved lp_curr; the gradient at the saved
  # position is re-evaluated once below.)
  fz <- .hmc_fused_target(log_post_fn, grad_fn, state_init, par_names,
                          transform = transform, verbose = verbose,
                          sampler = "NUTS")
  vg_fn <- if (is.null(fz)) NULL else fz$vg
  ## gradient at state_init when known (fused; also the separate-call
  ## path once the step-size search has taken it)
  g_init <- if (is.null(fz)) NULL else fz$g0

  # --- Mass matrix ---
  # Dense path: M_inv (d×d) and chol_M (upper Cholesky of M = solve(M_inv))
  # are pre-supplied; diagonal mass adaptation is disabled so the caller's
  # fixed dense metric is preserved throughout warmup.
  # Diagonal path: mass_diag = 1/diag(Sigma_prop) pre-conditions NUTS to the
  # posterior scale; the windowed mass adaptation (adapt_mass=TRUE) will
  # overwrite this starting point as warmup progresses.
  use_dense <- !is.null(M_inv) && is.matrix(M_inv)
  if (use_dense) {
    ## Dense path: sentinel diagonals (never used in the hot path)
    M_diag     <- rep(1, d)
    M_inv_diag <- NULL  ## signals: ignore diagonal path
  } else if (!is.null(mass_diag) && length(mass_diag) == d) {
    M_diag     <- pmax(as.numeric(mass_diag), 1e-12)
    M_inv_diag <- 1 / M_diag
  } else {
    M_diag     <- rep(1, d)
    M_inv_diag <- rep(1, d)
  }

  # --- Find initial step size ---
  # On resume the frozen step size comes from the saved state; skip the
  # expensive step-size search entirely (it would be thrown away anyway).
  if (ckpt_resume) {
    eps0 <- 1  # placeholder; overwritten from saved state before the loop
  } else if (is.null(step_size)) {
    ## The separate-call path takes the start gradient once, for the
    ## search and the first transition alike
    if (is.null(fz)) g_init <- .grad(state_init)
    eps0 <- .hmc_find_stepsize(state_init, .lp_scalar, .grad,
                                M_inv_diag, M_diag,
                                M_inv = M_inv, chol_M = chol_M,
                                vg_fn = vg_fn, lp0 = fz$lp0, g0 = g_init)
  } else {
    eps0 <- step_size
  }
  if (verbose && !ckpt_resume) .dynhr_inform(sprintf("NUTS: initial step_size = %.4e", eps0))

  # --- Dual averaging parameters ---
  mu      <- log(10 * eps0)   # target log step size
  eps_bar <- 1                # averaged step size (log scale tracking)
  H_bar   <- 0                # averaged acceptance stat
  gamma_da <- 0.05
  t0_da    <- 10
  kappa_da <- 0.75
  eps_m    <- eps0
  # da_m tracks the iteration count for dual averaging; reset to 0 at each
  # mass-matrix update so the averager sees a fresh trajectory.
  da_m <- 0L

  # --- Windowed warmup schedule (Stan-style) ---
  # metric = "lowrank" (Lao 2026) is active only when adapting and no fixed
  # dense metric was supplied; it uses its own dimension-derived schedule.
  use_lowrank <- identical(metric, "lowrank") && !use_dense && adapt_mass
  # metric = "fisher_diag" (Seyboldt et al. 2026) keeps the Stan schedule and
  # the diagonal hot path; only the slow-window estimator differs (draws AND
  # scores), so it needs the per-iteration gradients too.
  use_fisher_diag <- identical(metric, "fisher_diag") && !use_dense && adapt_mass
  use_scores      <- use_lowrank || use_fisher_diag
  warmup_windows <- if (adapt_mass && n_warmup >= 10L) {
    if (use_lowrank) {
      .lowrank_warmup_windows(n_warmup, d)
    } else {
      .nuts_warmup_windows(n_warmup, init_buffer = init_buffer,
                            base_window = base_window, term_buffer = term_buffer)
    }
  } else {
    NULL
  }
  # Low-rank / Fisher-diagonal bookkeeping: gradients at the warmup states (the
  # Fisher estimators need scores as well as draws) and the per-window
  # decision history (low-rank only).
  grad_chain <- if (use_scores) matrix(NA_real_, n_warmup, d) else NULL
  lr_history <- list()
  lr_last    <- NULL

  # Build a lookup: for each warmup iteration m (1-indexed in the loop),
  # is it the END of a slow window?  We precompute slow-window end indices.
  slow_window_ends <- if (!is.null(warmup_windows)) {
    warmup_windows$end[warmup_windows$type == "slow"]
  } else {
    integer(0)
  }
  # Also track which slow-window a given iteration falls in, for
  # within-window accumulation.
  # (We use a running window-start/end updated when we hit a slow boundary.)
  cur_slow_start <- if (length(slow_window_ends) > 0L) {
    warmup_windows$start[warmup_windows$type == "slow"][1L]
  } else {
    NA_integer_
  }
  cur_slow_end <- if (length(slow_window_ends) > 0L) {
    slow_window_ends[1L]
  } else {
    NA_integer_
  }
  slow_win_idx <- 1L  # which slow window we are currently in

  # --- Storage ---
  # `chain`/`full_chain` are always THETA-SPACE (mapped back via
  # to_constrained() when transformed). `state_chain` tracks the sampler
  # state (eta-space when transformed, theta-space otherwise) and is used
  # only for mass-matrix adaptation (only the current slow-window rows are
  # ever needed; in checkpoint mode we bound it to the adaptation window
  # rather than n_total rows).
  #
  # In checkpoint mode: chain/logpost_trace live on disk; we hold only a
  # flush-sized buffer in RAM. Small bookkeeping vectors (treedepths,
  # divergences, energy_trace) are retained in full because they are small
  # (integer/logical/double, n_total elements) and are needed for diagnostics.
  # On resume these are initialised to their neutral values for the resumed
  # draws only (warmup rows are zeros/FALSE -- diagnostics cover post-warmup).
  if (ckpt) {
    buf       <- matrix(NA_real_, nrow = flush_every, ncol = d)
    buf_lp    <- numeric(flush_every)
    buf_i     <- 0L
    chain     <- NULL
    logpost_trace <- NULL
    # state_chain: bound to the largest possible slow window. The maximum
    # window size is at most n_warmup rows; preallocate that but no more.
    # (Non-checkpoint path keeps n_total rows for simplicity / bit-identity.)
    state_chain <- matrix(NA_real_, n_warmup, d, dimnames = list(NULL, par_names))
  } else {
    chain         <- matrix(NA_real_, n_total, d, dimnames = list(NULL, par_names))
    logpost_trace <- numeric(n_total)
    state_chain   <- matrix(NA_real_, n_total, d, dimnames = list(NULL, par_names))
  }
  treedepths    <- integer(n_total)
  divergences   <- logical(n_total)
  energy_trace  <- numeric(n_total)  # Hamiltonian H = logpost - kinetic (for BFMI)
  # Per-iteration NUTS acceptance statistic (the trajectory-averaged alpha the
  # dual averaging targets) and whether the draw moved. NA = no transition
  # recorded in this run (draw 1; the pre-resume rows of a resumed run), so
  # the reported means cover exactly the transitions this run made.
  accept_stats  <- rep(NA_real_, n_total)
  moved         <- rep(NA, n_total)
  n_grad_evals  <- 0L

  theta   <- state_init
  ## Fused: the start point's value comes from the same function as
  ## every leaf's; g_curr is the gradient there.
  lp_curr <- if (is.null(fz)) .lp_scalar(theta) else fz$lp0
  g_curr  <- g_init   # Carried on the separate-call path too
  trace_lp_curr <- if (!is.null(transform)) {
    lp_curr - transform$log_jacobian(theta)
  } else {
    lp_curr
  }
  n_divergent_total <- 0L   # running count (supplement divergences[] for resume)
  i_start <- 1L             # loop starts at i_start + 1 (i.e., i=2 on a fresh run)

  if (ckpt_resume) {
    # ---- Continue a saved run. Refuse mismatched configuration, restore the
    # frozen post-warmup adaptation state (step size, mass matrix, dense path
    # matrices), position, log-posterior, accept count, draw count, and --
    # critically -- the RNG state, so the continuation is bit-identical to a
    # single longer run.
    .ckpt_meta_verify(ckpt_paths$meta, "nuts", checkpoint$fingerprint)
    st <- .ckpt_load_state(ckpt_paths$state)
    if (st$n_done < st$n_warmup)
      stop("checkpoint resume: saved state is mid-warmup (n_done=", st$n_done,
           " < n_warmup=", st$n_warmup, "). NUTS checkpoints are only resumable ",
           "after warmup completes. Start a fresh run.", call. = FALSE)
    theta         <- st$theta
    lp_curr       <- st$lp_curr
    trace_lp_curr <- st$trace_lp_curr
    ## fused: the gradient at the saved position (a deterministic function of
    ## it, so the continuation matches the uninterrupted run)
    g_curr <- if (!is.null(vg_fn)) vg_fn(theta)$grad else NULL
    eps_m         <- st$eps_m
    M_diag        <- st$M_diag
    M_inv_diag    <- st$M_inv_diag
    if (!is.null(st$M_inv))    M_inv    <- st$M_inv
    if (!is.null(st$chol_M))   chol_M   <- st$chol_M
    n_warmup      <- st$n_warmup        # original warmup fixes the retained set
    n_divergent_total <- st$n_divergent_total
    i_start       <- st$n_done
    .ckpt_truncate(ckpt_paths, st$n_done, d)  # drop any post-state partial flush
    assign(".Random.seed", st$rng, envir = .GlobalEnv)
    # Recompute n_total from the ORIGINAL n_draws for this run.
    # On resume the caller passes the EXTENDED n_draws (original + extra).
    n_total <- n_draws + n_warmup
    # Extend the small bookkeeping vectors if the new n_total is larger.
    if (length(treedepths) < n_total) {
      treedepths    <- c(treedepths,   integer(n_total - length(treedepths)))
      divergences   <- c(divergences,  logical(n_total - length(divergences)))
      energy_trace  <- c(energy_trace, numeric(n_total - length(energy_trace)))
      accept_stats  <- c(accept_stats, rep(NA_real_, n_total - length(accept_stats)))
      moved         <- c(moved, rep(NA, n_total - length(moved)))
    }
  } else {
    # ---- Fresh run: record draw 1 (to the streaming buffer or in-RAM chain).
    stored1 <- if (!is.null(transform)) transform$to_constrained(theta) else theta
    if (ckpt) {
      unlink(c(ckpt_paths$draws, ckpt_paths$lp))   # clear any stale fresh-run files
      if (!isFALSE(checkpoint$write_meta))
        .ckpt_meta_write(ckpt_paths$meta, "nuts", checkpoint$fingerprint)
      buf_i <- 1L; buf[1, ] <- stored1; buf_lp[1] <- trace_lp_curr
    } else {
      chain[1, ]       <- stored1
      logpost_trace[1] <- trace_lp_curr
    }
    state_chain[1, ]  <- theta
    treedepths[1]     <- 0L
    divergences[1]    <- FALSE
    # energy_trace[1] uses r=0 (kinetic=0)
    energy_trace[1]   <- if (!is.null(transform)) {
      lp_curr - transform$log_jacobian(theta)
    } else {
      lp_curr
    }
  }
  t_start <- Sys.time()

  for (m in i_start + seq_len(max(0L, n_total - i_start))) {
    # --- One NUTS transition (momentum draw + multinomial tree) ---
    tr <- .nuts_transition(theta, lp_curr, eps_m, .lp_scalar, .grad,
                           M_diag, M_inv_diag, M_inv, chol_M,
                           max_treedepth, delta_max,
                           g_curr = g_curr, vg_fn = vg_fn)
    theta_m       <- tr$theta
    j             <- tr$depth
    any_divergent <- tr$divergent
    alpha_sum     <- tr$alpha_sum
    n_alpha       <- tr$n_alpha
    n_grad_evals  <- n_grad_evals + tr$n_leaves
    # Acceptance stat averaged over the whole trajectory (see the accumulator
    # in .nuts_transition()): recorded for every iteration, fed to the dual
    # averaging during warmup.
    alpha_m         <- if (n_alpha > 0) alpha_sum / n_alpha else 0
    accept_stats[m] <- alpha_m
    moved[m]        <- any(theta_m != theta)
    # energy_trace stores the Hamiltonian using the THETA-SPACE logpost (for
    # BFMI), i.e. with the eta-space Jacobian subtracted back out -- the
    # NUTS dynamics themselves (tree building, U-turn checks) use the
    # eta-space `joint0`, unchanged.
    energy_trace[m] <- if (!is.null(transform)) {
      tr$joint0 - transform$log_jacobian(theta)
    } else {
      tr$joint0
    }

    # --- Update state ---
    names(theta_m) <- par_names
    theta   <- theta_m
    ## Fused: the selected point's value and gradient came with it (the leaf
    ## evaluation, or the unchanged start point's); no re-evaluation.
    ## The separate-call path carries the gradient too (same values).
    lp_curr <- if (is.null(vg_fn)) .lp_scalar(theta) else tr$lp
    g_curr  <- tr$g
    trace_lp_curr <- if (!is.null(transform)) {
      lp_curr - transform$log_jacobian(theta)
    } else {
      lp_curr
    }

    stored_m <- if (!is.null(transform)) transform$to_constrained(theta) else theta
    if (ckpt) {
      # Checkpoint path: buffer the draw; flush to disk every flush_every rows.
      buf_i <- buf_i + 1L
      buf[buf_i, ]  <- stored_m
      buf_lp[buf_i] <- trace_lp_curr
      if (buf_i >= flush_every || m == n_total) {
        .ckpt_append_draws(ckpt_paths$draws, buf[seq_len(buf_i), , drop = FALSE])
        .ckpt_append_lp(ckpt_paths$lp, buf_lp[seq_len(buf_i)])
        .ckpt_save_state(ckpt_paths$state, list(
          theta         = theta,
          lp_curr       = lp_curr,
          trace_lp_curr = trace_lp_curr,
          eps_m         = eps_m,
          M_diag        = M_diag,
          M_inv_diag    = M_inv_diag,
          M_inv         = M_inv,
          chol_M        = chol_M,
          n_done        = m,
          n_warmup      = n_warmup,
          n_draws_target = n_draws,
          n_divergent_total = n_divergent_total + as.integer(any_divergent),
          rng           = get(".Random.seed", envir = .GlobalEnv)
        ))
        buf_i <- 0L
      }
    } else {
      chain[m, ]       <- stored_m
      logpost_trace[m] <- trace_lp_curr
    }
    # state_chain: only needed during warmup for mass adaptation.
    # In checkpoint mode state_chain is bounded to n_warmup rows (1-indexed
    # directly by m during warmup).  In non-checkpoint mode it is n_total rows
    # indexed by m.  Either way we only write during warmup.
    if (m <= n_warmup) state_chain[m, ] <- theta
    if (use_scores && m <= n_warmup && !is.na(cur_slow_start) && m >= cur_slow_start) {
      if (is.null(g_curr)) {
        g_curr          <- .grad(theta)
        n_grad_evals    <- n_grad_evals + 1L
      }
      grad_chain[m, ] <- g_curr   # carried with the state: no extra call
    }
    treedepths[m]  <- j
    divergences[m] <- any_divergent
    if (any_divergent) n_divergent_total <- n_divergent_total + 1L

    # --- Dual averaging (step size adaptation during warmup) ---
    if (m <= n_warmup) {
      da_m <- da_m + 1L
      w <- 1 / (da_m + t0_da)
      H_bar <- (1 - w) * H_bar + w * (target_accept - alpha_m)
      log_eps_m <- mu - (sqrt(da_m) / gamma_da) * H_bar
      eps_m <- exp(log_eps_m)
      m_kappa <- da_m^(-kappa_da)
      eps_bar <- exp(m_kappa * log_eps_m + (1 - m_kappa) * log(eps_bar))
    }

    # --- Windowed mass-matrix adaptation ---
    # At the END of each slow window: compute within-window variance (diagonal)
    # or within-window covariance (warmup_dense), update the mass matrix,
    # re-find step size, and RESET dual averaging.
    # The reset at EVERY slow-window boundary ensures the dual averager does not
    # accumulate stats from the wrong mass geometry.
    # Disabled when use_dense=TRUE: the caller's fixed dense metric is preserved
    # (and once warmup_dense succeeds, use_dense is set TRUE to freeze it).
    if (!use_dense && adapt_mass && m <= n_warmup && !is.na(cur_slow_end) && m == cur_slow_end) {
      # Within-window samples (in sampler-state space, i.e. eta when transformed)
      # In checkpoint mode state_chain is bounded to n_warmup rows indexed by m
      # directly (same as non-checkpoint path during warmup).
      win_rows <- seq.int(cur_slow_start, cur_slow_end)
      win_rows <- win_rows[win_rows >= 1L & win_rows <= m]
      if (use_lowrank && length(win_rows) < 3L) {
        ## Window too short for the low-rank estimator: keep the current metric.
        lr_history[[length(lr_history) + 1L]] <- data.frame(
          start = cur_slow_start, end = cur_slow_end, n = length(win_rows),
          rank = NA_integer_, threshold = NA_real_, step_size = eps_m)
      } else if (use_lowrank) {
        ## Low-rank-plus-diagonal (Lao 2026 controller, single chain): the
        ## window's draws + scores give the Fisher estimate; rank 0 = the
        ## diagonal route, rank >= 1 = promotion. Every window is memoryless
        ## and re-estimated; step size re-found and dual averaging reset.
        lr_est <- .lowrank_estimate(state_chain[win_rows, , drop = FALSE],
                                    grad_chain[win_rows, , drop = FALSE],
                                    cutoff   = lr_ctrl$cutoff,
                                    max_rank = lr_ctrl$max_rank,
                                    gamma    = lr_ctrl$gamma)
        lr_last    <- lr_est
        M_inv      <- lr_est$metric
        chol_M     <- lr_est$metric
        M_diag     <- rep(1, d)   ## sentinel (low-rank path active)
        M_inv_diag <- NULL        ## signals: ignore diagonal path
        eps0 <- .hmc_find_stepsize(theta, .lp_scalar, .grad,
                                   M_inv_diag, M_diag,
                                   M_inv = M_inv, chol_M = chol_M,
                                   vg_fn = vg_fn, lp0 = lp_curr, g0 = g_curr)
        mu      <- log(10 * eps0)
        eps_bar <- 1
        H_bar   <- 0
        eps_m   <- eps0
        da_m    <- 0L
        lr_history[[length(lr_history) + 1L]] <- data.frame(
          start = cur_slow_start, end = cur_slow_end, n = length(win_rows),
          rank = lr_est$rank, threshold = lr_est$threshold, step_size = eps0)
        if (verbose) .dynhr_inform(sprintf(
          "NUTS: slow window [%d,%d] -> lowrank metric rank %d (k_max %d), reset step_size = %.4e",
          cur_slow_start, cur_slow_end, lr_est$rank, lr_est$k_max, eps0))
      } else if (length(win_rows) >= 2L) {
        if (identical(metric, "warmup_dense")) {
          ## Dense path: apply Ledoit-Wolf shrunken covariance -> M_inv.
          ## For all-but-last slow windows: update diagonal (same as default)
          ## so step-size adaptation benefits from improving variance estimates.
          ## At the FINAL slow window only: attempt dense; freeze if successful;
          ## fall back to the already-updated diagonal on failure.
          ## This mirrors Stan's logic: diagonal updates each window to keep DA
          ## well-conditioned, then a single dense update at the final window.
          is_last_slow_window <- (cur_slow_end == slow_window_ends[length(slow_window_ends)])

          if (is_last_slow_window) {
            win_draws  <- state_chain[win_rows, , drop = FALSE]
            dense_res  <- .warmup_dense_metric(win_draws, verbose = verbose)
            if (!is.null(dense_res)) {
              M_inv      <- dense_res$M_inv
              chol_M     <- dense_res$chol_M
              M_diag     <- rep(1, d)   ## sentinel (dense path active)
              M_inv_diag <- NULL        ## signals: ignore diagonal path
              use_dense  <- TRUE        ## freeze: no further mass adaptation
              eps0 <- .hmc_find_stepsize(theta, .lp_scalar, .grad,
                                         M_inv_diag, M_diag,
                                         M_inv = M_inv, chol_M = chol_M,
                                         vg_fn = vg_fn, lp0 = lp_curr,
                                         g0 = g_curr)
              mu      <- log(10 * eps0)
              eps_bar <- 1
              H_bar   <- 0
              eps_m   <- eps0
              da_m    <- 0L
              if (verbose) .dynhr_inform(sprintf(
                "NUTS: slow window [%d,%d] (final) -> warmup_dense mass set (lambda=%.3f), reset step_size = %.4e",
                cur_slow_start, cur_slow_end, dense_res$lambda, eps0))
            } else {
              ## Fallback to diagonal for the final window
              vars <- apply(state_chain[win_rows, , drop = FALSE], 2, var)
              vars[!is.finite(vars) | vars < 1e-12] <- 1
              ## Stan's rule: the INVERSE mass is the posterior variance. This used
              ## to set the MASS to the variance (inverted), squaring the problem's
              ## conditioning instead of removing it.
              M_inv_diag <- vars
              M_diag     <- 1 / vars
              eps0 <- .hmc_find_stepsize(theta, .lp_scalar, .grad, M_inv_diag, M_diag,
                                         vg_fn = vg_fn, lp0 = lp_curr, g0 = g_curr)
              mu      <- log(10 * eps0)
              eps_bar <- 1
              H_bar   <- 0
              eps_m   <- eps0
              da_m    <- 0L
              if (verbose) .dynhr_inform(sprintf(
                "NUTS: slow window [%d,%d] (final) -> warmup_dense fallback to diagonal, reset step_size = %.4e",
                cur_slow_start, cur_slow_end, eps0))
            }
          } else {
            ## Intermediate slow window: update diagonal to improve step-size DA
            vars <- apply(state_chain[win_rows, , drop = FALSE], 2, var)
            vars[!is.finite(vars) | vars < 1e-12] <- 1
            ## Stan's rule: the INVERSE mass is the posterior variance. This used
            ## to set the MASS to the variance (inverted), squaring the problem's
            ## conditioning instead of removing it.
            M_inv_diag <- vars
            M_diag     <- 1 / vars
            eps0 <- .hmc_find_stepsize(theta, .lp_scalar, .grad, M_inv_diag, M_diag,
                                       vg_fn = vg_fn, lp0 = lp_curr, g0 = g_curr)
            mu      <- log(10 * eps0)
            eps_bar <- 1
            H_bar   <- 0
            eps_m   <- eps0
            da_m    <- 0L
            if (verbose) .dynhr_inform(sprintf(
              "NUTS: slow window [%d,%d] -> diag update (warmup_dense pending final), reset step_size = %.4e",
              cur_slow_start, cur_slow_end, eps0))
          }
        } else {
          ## Default diagonal adaptation
          if (use_fisher_diag) {
            ## Fisher-divergence diagonal (Seyboldt et al. 2026, Thm 2.2):
            ## inverse mass sqrt(var(x) / var(g)) from draws AND scores. Rows
            ## without a stored score (a window starting at the unrecorded
            ## first state) are dropped.
            fd_rows <- win_rows[stats::complete.cases(grad_chain[win_rows, , drop = FALSE])]
            vars <- if (length(fd_rows) >= 2L) {
              .fisher_diag_estimate(state_chain[fd_rows, , drop = FALSE],
                                    grad_chain[fd_rows, , drop = FALSE])
            } else {
              apply(state_chain[win_rows, , drop = FALSE], 2, var)
            }
          } else {
            vars <- apply(state_chain[win_rows, , drop = FALSE], 2, var)
          }
          vars[!is.finite(vars) | vars < 1e-12] <- 1
          ## Stan's rule: the INVERSE mass is the posterior variance. This used
          ## to set the MASS to the variance (inverted), squaring the problem's
          ## conditioning instead of removing it.
          M_inv_diag <- vars
          M_diag     <- 1 / vars
          # Re-find step size and reset dual averaging
          eps0 <- .hmc_find_stepsize(theta, .lp_scalar, .grad, M_inv_diag, M_diag,
                                     vg_fn = vg_fn, lp0 = lp_curr, g0 = g_curr)
          mu      <- log(10 * eps0)
          eps_bar <- 1
          H_bar   <- 0
          eps_m   <- eps0
          da_m    <- 0L
          if (verbose) {
            .dynhr_inform(sprintf("NUTS: slow window [%d,%d] -> %s mass update, reset step_size = %.4e",
                            cur_slow_start, cur_slow_end,
                            if (use_fisher_diag) "fisher_diag" else "diagonal", eps0))
          }
        }
      }

      # Advance to the next slow window
      slow_win_idx <- slow_win_idx + 1L
      if (slow_win_idx <= length(slow_window_ends)) {
        cur_slow_start <- warmup_windows$start[warmup_windows$type == "slow"][slow_win_idx]
        cur_slow_end   <- slow_window_ends[slow_win_idx]
      } else {
        cur_slow_start <- NA_integer_
        cur_slow_end   <- NA_integer_
      }
    }

    # --- Fix step size at end of warmup (use dual-averaged value) ---
    if (m == n_warmup) {
      eps_m <- eps_bar
      if (verbose) .dynhr_inform(sprintf("NUTS: warmup complete, final step_size = %.4e", eps_m))
    }

    # --- Progress ---
    if (m %% 200 == 0 || m == n_total) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / m * (n_total - m)
      ch_lab  <- if (is.null(chain_id)) "?" else as.character(chain_id)
      n_div   <- sum(divergences[seq_len(m)])
      phase   <- if (m <= n_warmup) "warmup" else "sample"
      msg <- sprintf("NUTS Ch%s [%s] %d/%d depth=%.1f lp=%.1f eps=%.3e div=%d ETA=%.0fs",
                      ch_lab, phase, m, n_total,
                      mean(treedepths[max(1, m - 199):m]),
                      trace_lp_curr, eps_m, n_div, eta)
      if (!is.null(progressor)) {
        progressor(message = msg, amount = 1)
      } else if (verbose) {
        .dynhr_inform(msg)
      }
    }
  }

  if (ckpt) {
    # Materialize the full chain from the streamed files for the return value.
    # checkpoint$return_chain = FALSE skips this for very long runs (the draws
    # remain on disk; read them with .ckpt_read_draws() or resume to extend).
    logpost_trace <- .ckpt_read_lp(ckpt_paths$lp)
    chain <- if (isFALSE(checkpoint$return_chain)) NULL else
      .ckpt_read_draws(ckpt_paths$draws, d, par_names)
  }

  elapsed        <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  post_chain     <- if (is.null(chain)) NULL else
    chain[(n_warmup + 1):n_total, , drop = FALSE]
  post_logpost   <- logpost_trace[(n_warmup + 1):n_total]
  post_depths    <- treedepths[(n_warmup + 1):n_total]
  post_divs      <- divergences[(n_warmup + 1):n_total]
  post_energy    <- energy_trace[(n_warmup + 1):n_total]
  post_accept    <- accept_stats[(n_warmup + 1):n_total]
  post_moved     <- moved[(n_warmup + 1):n_total]
  accept_stat    <- if (any(!is.na(post_accept))) mean(post_accept, na.rm = TRUE)
                    else NA_real_
  move_rate      <- if (any(!is.na(post_moved))) mean(post_moved, na.rm = TRUE)
                    else NA_real_

  out <- list(
    chain           = post_chain,
    full_chain      = chain,
    logpost_trace   = logpost_trace,
    post_logpost    = post_logpost,
    ## acceptance_rate is the mean post-warmup NUTS acceptance statistic
    ## (Stan's mean accept_stat__; what target_accept targets). Before
    ## 0.9.3.127 it was 1 - mean(treedepth == 0), identically 1.
    acceptance_rate = accept_stat,
    accept_stat     = accept_stat,
    accept_stats    = post_accept,
    move_rate       = move_rate,
    step_size       = eps_m,
    mass_matrix     = if (inherits(M_inv, "dynhr_lowrank_metric")) .lowrank_dense_inv(M_inv)
                      else if (use_dense && !is.null(M_inv)) M_inv else M_diag,
    treedepths      = post_depths,
    divergences     = post_divs,
    n_divergent     = sum(post_divs),
    mean_treedepth  = mean(post_depths),
    energy_trace    = post_energy,
    n_draws         = as.integer(n_draws),
    n_burn          = as.integer(n_warmup),
    n_grad_evals    = n_grad_evals,
    elapsed_secs    = elapsed,
    sampler         = "nuts",
    checkpoint_dir  = if (ckpt) checkpoint$dir else NULL,
    kernel_stats    = if (is.null(grad_fn_kernel_stats)) NULL
                       else as.list(grad_fn_kernel_stats)
  )
  if (use_lowrank) {
    out$lowrank <- list(
      rank        = if (is.null(lr_last)) 0L else lr_last$rank,
      promoted    = any(vapply(lr_history, function(h) isTRUE(h$rank >= 1L), logical(1))),
      eigenvalues = if (is.null(lr_last)) NULL else lr_last$eigenvalues,
      threshold   = if (is.null(lr_last)) NA_real_ else lr_last$threshold,
      metric      = if (inherits(M_inv, "dynhr_lowrank_metric")) M_inv else NULL,
      history     = if (length(lr_history)) do.call(rbind, lr_history) else NULL
    )
  }
  out
}


# ============================================================================
# Diagnostics
# ============================================================================

#' Print summary of HMC/NUTS results
#'
#' @param result Output from dynhr_hmc() or dynhr_nuts()
#' @param probs Quantiles to report
#' @noRd
