## R/sampler-rwmh.R
## --------------------------------------------------------------------------
## Phase-2 split from estimation-monolith.R.
##
## Random Walk Metropolis-Hastings (RWMH) sampler.
## Called by estimate_model() in estimate-monolith.R and by run_mcmc_parallel().
## --------------------------------------------------------------------------

## --------------------------------------------------------------------------
## Argument checkers shared by the MCMC samplers
##
## Every sampler keeps ONE `.<sampler>_args_problem(args, n_par = NULL)`: a
## pure function of the arguments (no RNG, no likelihood evaluation) that
## returns the problems as a character vector (`character(0)` = fine). The
## sampler calls it at entry and aborts with the joined messages; a caller
## that already knows the arguments (a spec being validated) calls the same
## function, so a bad value fails when the spec is built rather than after
## the mode stage, and the two cannot drift apart.
##
## `args` is a named list keyed by the SAMPLER FUNCTION's formal argument
## names and holds only what the caller supplies: an absent (or NULL) entry
## is the sampler's default, which is valid by construction. `n_par` is the
## number of estimated parameters when known; rules that need it are skipped
## when it is NULL.
## --------------------------------------------------------------------------

## Short rendering of an offending value for a message.
.mcmc_fmt <- function(v) {
  if (is.numeric(v) && length(v) == 1L) format(v, digits = 6)
  else if (is.character(v) && length(v) == 1L) paste0("\"", v, "\"")
  else paste0("<", class(v)[1L], " of length ", length(v), ">")
}

.mcmc_scalar <- function(v) is.numeric(v) && length(v) == 1L && !is.na(v)

## A single finite whole number.
.mcmc_whole_ok <- function(v) .mcmc_scalar(v) && is.finite(v) && v == round(v)

## Rule constructors. A rule is function(v, nm, n_par) -> NULL | message.
.mcmc_r_whole <- function(min) function(v, nm, n_par) {
  if (.mcmc_whole_ok(v) && v >= min) return(NULL)
  sprintf("`%s` must be a whole number >= %s (got %s).", nm, min, .mcmc_fmt(v))
}

.mcmc_r_pos <- function(allow_inf = FALSE) function(v, nm, n_par) {
  if (.mcmc_scalar(v) && v > 0 && (allow_inf || is.finite(v))) return(NULL)
  sprintf("`%s` must be a single positive %snumber (got %s).", nm,
          if (allow_inf) "" else "finite ", .mcmc_fmt(v))
}

## Number in an interval; `lo_open` / `hi_open` say whether the end is excluded.
.mcmc_r_range <- function(lo, hi, lo_open, hi_open) function(v, nm, n_par) {
  ok <- .mcmc_scalar(v) && is.finite(v) &&
    (if (lo_open) v > lo else v >= lo) && (if (hi_open) v < hi else v <= hi)
  if (ok) return(NULL)
  sprintf("`%s` must be a single number in %s%s, %s%s (got %s).", nm,
          if (lo_open) "(" else "[", lo, hi, if (hi_open) ")" else "]",
          .mcmc_fmt(v))
}

.mcmc_r_flag <- function() function(v, nm, n_par) {
  if (isTRUE(v) || isFALSE(v)) return(NULL)
  sprintf("`%s` must be TRUE or FALSE (got %s).", nm, .mcmc_fmt(v))
}

.mcmc_r_fun <- function() function(v, nm, n_par) {
  if (is.function(v)) return(NULL)
  sprintf("`%s` must be NULL or a function (got %s).", nm, .mcmc_fmt(v))
}

## `choices` are the sampler's own match.arg() choices; the unevaluated
## default vector (all choices) is accepted as "use the default".
.mcmc_r_choice <- function(choices) function(v, nm, n_par) {
  if (is.character(v) && (identical(v, choices) ||
                          (length(v) == 1L && v %in% choices))) return(NULL)
  sprintf("`%s` must be one of %s (got %s).", nm,
          paste0("\"", choices, "\"", collapse = ", "), .mcmc_fmt(v))
}

## Diagonal of the mass matrix: one finite positive entry per parameter.
.mcmc_r_mass <- function() function(v, nm, n_par) {
  if (!is.numeric(v) || !length(v) || any(!is.finite(v)) || any(v <= 0))
    return(sprintf("`%s` must be a numeric vector of finite positive entries.", nm))
  if (!is.null(n_par) && length(v) != n_par)
    return(sprintf("`%s` must have one entry per estimated parameter (%d; got %d).",
                   nm, n_par, length(v)))
  NULL
}

## A square numeric matrix of the right order (a Cholesky factor).
.mcmc_r_square <- function() function(v, nm, n_par) {
  if (!is.matrix(v) || !is.numeric(v) || nrow(v) != ncol(v) || any(!is.finite(v)))
    return(sprintf("`%s` must be a finite square numeric matrix.", nm))
  if (!is.null(n_par) && nrow(v) != n_par)
    return(sprintf("`%s` must be %d x %d, one row/column per estimated parameter (got %d x %d).",
                   nm, n_par, n_par, nrow(v), ncol(v)))
  NULL
}

## Symmetric positive definite matrix of order n_par. `strict` is for a
## matrix that is factored (chol) or inverted: its smallest eigenvalue must
## clear the numerical-rank tolerance. Without it (a proposal covariance,
## which the sampler regularises itself) a numerically singular positive
## semi-definite matrix is accepted, but an indefinite / negative-definite
## one, or one with a zero or negative variance (a frozen coordinate), is not.
.mcmc_r_spd <- function(strict = FALSE) function(v, nm, n_par) {
  if (!is.matrix(v) || !is.numeric(v) || nrow(v) != ncol(v))
    return(sprintf("`%s` must be a square numeric matrix (got %s).", nm, .mcmc_fmt(v)))
  if (!is.null(n_par) && nrow(v) != n_par)
    return(sprintf("`%s` must be %d x %d, one row/column per estimated parameter (got %d x %d).",
                   nm, n_par, n_par, nrow(v), ncol(v)))
  if (any(!is.finite(v)))
    return(sprintf("`%s` must be finite.", nm))
  sc <- max(abs(v))
  if (sc == 0)
    return(sprintf("`%s` must be positive definite (it is all zeros).", nm))
  asym <- max(abs(v - t(v)))
  if (asym > 1e-6 * sc)
    return(sprintf("`%s` must be symmetric (largest asymmetry %s).", nm,
                   format(asym, digits = 3)))
  ev <- eigen((v + t(v)) / 2, symmetric = TRUE, only.values = TRUE)$values
  lo <- min(ev); hi <- max(ev)
  bad <- if (strict) lo <= .Machine$double.eps * nrow(v) * hi
         else (hi <= 0 || lo < -1e-8 * hi || any(diag(v) <= 0))
  if (bad)
    return(sprintf("`%s` must be positive definite (smallest eigenvalue %s%s).",
                   nm, format(lo, digits = 3),
                   if (!strict && any(diag(v) <= 0))
                     "; every diagonal entry must be > 0" else ""))
  NULL
}

## Run `rules` (named list of rule functions) over the supplied `args`; also
## report names that are not formals of `fn`.
.mcmc_check_args <- function(args, fn, fn_name, rules, n_par = NULL) {
  if (is.null(args)) args <- list()
  nms <- names(args)
  if (length(args) && (is.null(nms) || any(!nzchar(nms))))
    return("arguments must be a named list.")
  problems <- character(0)
  ## A function that takes `...` (a wrapper) accepts any name.
  unknown <- if ("..." %in% names(formals(fn))) character(0)
             else setdiff(nms, names(formals(fn)))
  if (length(unknown))
    problems <- sprintf("unknown argument%s %s (%s() accepts: %s).",
                        if (length(unknown) > 1L) "s" else "",
                        paste0("`", unknown, "`", collapse = ", "), fn_name,
                        paste0("`", setdiff(names(formals(fn)),
                                            c("log_post_fn", "theta_init", "theta0")),
                               "`", collapse = ", "))
  for (nm in intersect(names(rules), nms)) {
    v <- args[[nm]]
    if (is.null(v)) next
    msg <- rules[[nm]](v, nm, n_par)
    if (!is.null(msg)) problems <- c(problems, msg)
  }
  problems
}

## The value of `nm`: the supplied one, else the sampler's own default.
.mcmc_arg <- function(args, nm, fn) {
  v <- args[[nm]]
  if (!is.null(v)) v else eval(formals(fn)[[nm]])
}

## Sampler entry: abort with the joined problems.
.mcmc_abort_if_problems <- function(fn_name, problems) {
  if (length(problems))
    .dynhr_abort(fn_name, ": ", paste(problems, collapse = " "),
                 class = "dynhr_error_invalid_argument")
  invisible(NULL)
}

## A metric argument may also be the low-rank metric object (class
## dynhr_lowrank_metric, R/sampler-metric.R) that pooled NUTS passes as
## M_inv / chol_M: its constructor validated it; check only its order.
.mcmc_r_or_lowrank <- function(rule) function(v, nm, n_par) {
  if (!inherits(v, "dynhr_lowrank_metric")) return(rule(v, nm, n_par))
  if (!is.null(n_par) && length(v$sd) != n_par)
    return(sprintf("`%s` (low-rank metric) must have order %d, one per estimated parameter (got %d).",
                   nm, n_par, length(v$sd)))
  NULL
}

## Metric rules shared by the gradient samplers that take them.
.mcmc_metric_rules <- function(fn) {
  r <- list(mass_diag = .mcmc_r_mass(),
            M_inv     = .mcmc_r_or_lowrank(.mcmc_r_spd(strict = FALSE)),
            chol_M    = .mcmc_r_or_lowrank(.mcmc_r_square()),
            grad_fn   = .mcmc_r_fun())
  r[intersect(names(r), names(formals(fn)))]
}


## --------------------------------------------------------------------------
## RWMH argument checker
## --------------------------------------------------------------------------

#' Problems with the arguments of rwmh()
#'
#' @param args Named list of the supplied rwmh() arguments (absent = default).
#' @param n_par Number of estimated parameters, or NULL when not known.
#'   \code{n_draws} counts burn-in draws too, so it must exceed \code{n_burn}.
#' @return Character vector of problems; \code{character(0)} when fine.
#' @noRd
.rwmh_args_problem <- function(args, n_par = NULL) {
  rules <- list(
    Sigma_prop  = .mcmc_r_spd(strict = FALSE),
    n_draws     = .mcmc_r_whole(1L),
    n_burn      = .mcmc_r_whole(0L),
    scale       = .mcmc_r_pos(),
    target_rate = .mcmc_r_range(0, 1, TRUE, TRUE),
    adapt_every = .mcmc_r_whole(1L),
    adapt_cov   = .mcmc_r_flag(),
    n_blocks    = function(v, nm, n_par) {
      ok <- .mcmc_whole_ok(v) && v >= 1 && (is.null(n_par) || v <= n_par)
      if (ok) return(NULL)
      sprintf("`n_blocks` must be a single integer in [1, %s] (got %s).",
              if (is.null(n_par)) "n_par" else n_par, .mcmc_fmt(v))
    })
  p <- .mcmc_check_args(args, rwmh, "rwmh", rules, n_par)
  nd <- .mcmc_arg(args, "n_draws", rwmh)
  nb <- .mcmc_arg(args, "n_burn", rwmh)
  if (.mcmc_whole_ok(nd) && .mcmc_whole_ok(nb) && nb >= 0 && nd >= 1 && nb >= nd)
    p <- c(p, sprintf(paste0("`n_draws` counts the burn-in draws too, so it must exceed ",
                             "`n_burn` (got n_draws = %s, n_burn = %s): no draws would be kept."),
                      nd, nb))
  p
}


#' Random Walk Metropolis-Hastings sampler
#'
#' @param log_post_fn  Function(theta) -> list(logpost, loglik, logprior)
#' @param theta0       Named initial parameter vector
#' @param Sigma_prop   Proposal covariance matrix (n_par x n_par). When
#'   `transform` is non-NULL this is interpreted as an ETA-SPACE (not
#'   theta-space) covariance -- see `transform` below.
#' @param n_draws      Total draws (including burn-in)
#' @param n_burn       Burn-in draws to discard from the returned chain
#' @param scale        Initial scaling factor for the proposal
#' @param target_rate  Target acceptance rate for adaptive scaling
#' @param adapt_every  Adapt scale every N draws (during burn-in)
#' @param verbose      Print progress messages
#' @param progressor   Optional progressr function (or NULL)
#' @param chain_id     Chain label for progress messages (or NULL)
#' @param transform    Optional "dynhr_param_transform" object (from
#'   \code{\link{build_param_transform}}). When non-NULL (opt-in), the chain
#'   is run in UNCONSTRAINED eta-space:
#'   \itemize{
#'     \item `eta0 = transform$to_unconstrained(theta0)` is the starting
#'       state; proposals/adaptation operate on eta.
#'     \item the target is `lp_eta = make_transformed_logpost(log_post_fn,
#'       transform, include_jacobian = TRUE)`, i.e. the theta-space
#'       log-posterior PLUS the change-of-variables Jacobian, so eta-space
#'       draws are correctly distributed (Metropolis-Hastings ratios cancel
#'       the proposal density, so a symmetric eta-space random walk needs no
#'       extra correction beyond the Jacobian already folded into `lp_eta`).
#'     \item `Sigma_prop` (and its adaptive `scale`) are interpreted as an
#'       ETA-SPACE covariance -- the caller is responsible for converting a
#'       theta-space covariance via the delta method (see
#'       \code{run_posterior_estimation}'s `transform_params` path).
#'     \item each STORED chain row is `to_constrained(eta)` (theta-space, as
#'       always); `logpost_trace` stores the THETA-SPACE log-posterior (the
#'       Jacobian term is subtracted back out via one extra cheap call to
#'       `transform$log_jacobian()`), so traces remain directly comparable
#'       with the untransformed sampler's output.
#'   }
#'   When NULL (default), behaviour is bit-identical to before.
#' @param adapt_cov    Opt-in (default \code{FALSE}): Haario et al. (2001)
#'   adaptive proposal covariance. When \code{TRUE}, a history of
#'   SAMPLING-SPACE states (eta-space when `transform` is set, theta-space
#'   otherwise) is kept for the burn-in draws. Every `adapt_every` draws
#'   during burn-in, once at least `max(200, 10 * n_par)` states have been
#'   observed, the empirical covariance of the history so far is used to
#'   recompute the proposal Cholesky factor:
#'   \code{Sigma_emp = cov(state_hist) + 1e-8 * diag(n_par)} (Haario's
#'   epsilon regularisation, which keeps `Sigma_emp` non-degenerate even when
#'   a coordinate has not yet moved) and `L <- .robust_chol(Sigma_emp,
#'   n_par)`. The existing scalar `scale` adaptation continues to run on top
#'   of the adapted shape -- it tunes the overall magnitude, so Haario's
#'   classic fixed `2.38^2 / d` factor is subsumed by `scale` and is not
#'   applied separately. Adaptation stops at the end of burn-in (the
#'   covariance and `L` are frozen for the retained draws), so detailed
#'   balance holds for the post-burn-in chain and no diminishing-adaptation
#'   condition is required. When \code{FALSE} (default), no state history is
#'   allocated and no extra random numbers are drawn; behaviour is
#'   bit-identical to before.
#' @param n_blocks     Opt-in (default \code{1L}): randomized parameter
#'   blocking (Chib & Ramamurthy 2010; implementation follows Herbst &
#'   Schorfheide 2015, ch. 4). When `n_blocks > 1`, each iteration draws a
#'   random permutation of `1:n_par` (one `sample()` call) and splits it into
#'   `n_blocks` contiguous chunks (as equal in size as possible). For each
#'   block in turn, only those coordinates are proposed
#'   (`theta_prop[block] = theta_curr[block] + scale * L_b %*% z`, with `L_b
#'   = .robust_chol(Sigma_curr[block, block], length(block))`, recomputed
#'   per block per iteration -- cheap relative to a posterior evaluation),
#'   the FULL posterior is evaluated at the block-updated point, and the
#'   point is accepted/rejected; the state evolves WITHIN the iteration (a
#'   later block proposes from the post-earlier-block state). Each draw thus
#'   costs `n_blocks` posterior evaluations. `accepted[i]` records whether
#'   ANY block moved; the scalar `scale` adaptation instead uses the MEAN
#'   per-block acceptance rate (the any-block rate saturates near 1 as
#'   `n_blocks` grows, which would otherwise make `scale` explode). The mean
#'   block acceptance rate is returned as `block_accept_rate` (`NA` when
#'   `n_blocks == 1`). Blocks index whatever space the sampler is operating
#'   in (eta-space when `transform` is set) and, with `adapt_cov = TRUE`, use
#'   submatrices of the adapted `Sigma_curr`. When \code{1L} (default),
#'   behaviour is bit-identical to before.
#' @return List: chain, full_chain, logpost_trace, post_logpost,
#'         acceptance_rate, scale, n_draws, n_burn, block_accept_rate
#' @noRd
rwmh <- function(log_post_fn, theta0, Sigma_prop,
                 n_draws = 10000L, n_burn = 5000L,
                 scale = 0.5, target_rate = 0.25,
                 adapt_every = 100L, verbose = TRUE,
                 progressor = NULL, chain_id = NULL,
                 transform = NULL,
                 adapt_cov = FALSE, n_blocks = 1L,
                 checkpoint = NULL) {

  n_par     <- length(theta0)
  par_names <- names(theta0)

  .mcmc_abort_if_problems("rwmh", .rwmh_args_problem(
    list(Sigma_prop = Sigma_prop, n_draws = n_draws, n_burn = n_burn,
         scale = scale, target_rate = target_rate, adapt_every = adapt_every,
         adapt_cov = adapt_cov, n_blocks = n_blocks), n_par = n_par))
  n_blocks <- as.integer(n_blocks)

  Sigma_curr <- Sigma_prop
  L <- .robust_chol(Sigma_curr, n_par)

  # ---- Checkpoint / streaming (opt-in). When `checkpoint` is a list carrying a
  # `dir`, draws are streamed to per-chain files in flush_every-row chunks (RAM
  # bounded by flush_every * n_par, not n_draws * n_par) and a restart state is
  # saved after every flush. `checkpoint$resume = TRUE` continues a prior run
  # from its saved state -- exactly (RNG, position, scale, proposal covariance),
  # so resume(N1)+resume(N2) is bit-identical to a single run of N1+N2.
  ckpt        <- !is.null(checkpoint)
  ckpt_resume <- ckpt && isTRUE(checkpoint$resume)
  flush_every <- if (ckpt) as.integer(checkpoint$flush_every %||% 1000L) else NA_integer_
  ckpt_paths  <- if (ckpt) .ckpt_paths(checkpoint$dir, chain_id) else NULL

  # ---- Transformed (opt-in) target: operate on eta = to_unconstrained(theta)
  if (!is.null(transform)) {
    lp_target  <- make_transformed_logpost(log_post_fn, transform, include_jacobian = TRUE)
    state0     <- transform$to_unconstrained(theta0)
    names(state0) <- par_names
  } else {
    lp_target  <- log_post_fn
    state0     <- theta0
  }

  # In checkpoint mode only a flush-sized buffer lives in RAM; the full chain is
  # read back from disk at the end for the return value. Otherwise pre-allocate
  # the full chain exactly as before (non-checkpoint path is byte-identical).
  if (ckpt) {
    buf    <- matrix(NA_real_, nrow = flush_every, ncol = n_par)
    buf_lp <- numeric(flush_every)
    buf_i  <- 0L
    chain  <- NULL
    logpost_trace <- NULL
  } else {
    chain         <- matrix(NA_real_, nrow = n_draws, ncol = n_par)
    colnames(chain) <- par_names
    logpost_trace <- numeric(n_draws)
  }
  accepted      <- logical(n_draws)

  # ---- Haario adaptive covariance (opt-in): history of SAMPLING-SPACE
  # states (eta-space when `transform` is set), preallocated to n_burn rows
  # since adaptation only runs during burn-in. When FALSE, no history is
  # allocated and no extra RNG draws occur anywhere below.
  state_hist <- if (adapt_cov) matrix(NA_real_, nrow = n_burn, ncol = n_par) else NULL
  adapt_min_n <- max(200L, 10L * n_par)

  # ---- Randomized blocking (opt-in): per-block acceptance bookkeeping.
  # `*_window` counters drive the scale adaptation (reset every adapt_every
  # draws); `*_total` counters accumulate over the whole run for the
  # reported `block_accept_rate`.
  block_accept_window <- 0L
  block_total_window  <- 0L
  block_accept_total  <- 0L
  block_total_total   <- 0L

  state_curr  <- state0
  result_curr <- lp_target(state_curr)
  lp_curr     <- result_curr$logpost
  if (!is.finite(lp_curr))
    stop("Initial parameter vector has -Inf log-posterior. Check starting values.")

  # Theta-space logpost for the trace (subtract the Jacobian back out when
  # transformed; identical to lp_curr otherwise).
  trace_lp_curr <- if (!is.null(transform)) {
    lp_curr - transform$log_jacobian(state_curr)
  } else {
    lp_curr
  }

  n_accept <- 0L
  i_start  <- 1L
  t_start  <- Sys.time()

  if (ckpt_resume) {
    # ---- Continue a saved run. First refuse a mismatched configuration (same
    # parameters are forced), then restore position, scale, proposal covariance,
    # accept count, burn-in and draw count, and -- crucially -- the RNG state,
    # so the continuation draws the SAME random stream a single long run would.
    .ckpt_meta_verify(ckpt_paths$meta, "rwmh", checkpoint$fingerprint)
    st <- .ckpt_load_state(ckpt_paths$state)
    state_curr    <- st$state_curr
    lp_curr       <- st$lp_curr
    trace_lp_curr <- st$trace_lp_curr
    scale         <- st$scale
    Sigma_curr    <- st$Sigma_curr
    L             <- .robust_chol(Sigma_curr, n_par)
    n_accept      <- st$n_accept
    n_burn        <- st$n_burn          # original burn-in fixes the retained set
    i_start       <- st$n_done
    .ckpt_truncate(ckpt_paths, st$n_done, n_par)  # drop any post-state partial flush
    assign(".Random.seed", st$rng, envir = .GlobalEnv)
  } else {
    # ---- Fresh run: record draw 1 (to the streaming buffer or the in-RAM chain).
    stored1 <- if (!is.null(transform)) transform$to_constrained(state_curr) else state_curr
    if (ckpt) {
      unlink(c(ckpt_paths$draws, ckpt_paths$lp))    # clear any stale fresh-run files
      # The shared meta.rds is written once by the orchestrator on the parallel
      # path (write_meta = FALSE there) to avoid a multi-daemon write race.
      if (!isFALSE(checkpoint$write_meta))
        .ckpt_meta_write(ckpt_paths$meta, "rwmh", checkpoint$fingerprint)
      buf_i <- 1L; buf[1, ] <- stored1; buf_lp[1] <- trace_lp_curr
    } else {
      chain[1, ]       <- stored1
      logpost_trace[1] <- trace_lp_curr
    }
    accepted[1] <- TRUE
    if (adapt_cov) state_hist[1, ] <- state_curr
  }

  # Fresh: i = 2..n_draws (draw 1 already recorded). Resume: i = n_done+1..n_draws
  # (seq_len(0) when no additional draws are requested -> loop is skipped).
  for (i in i_start + seq_len(max(0L, n_draws - i_start))) {

    if (n_blocks == 1L) {
      ## ---- Unblocked path: identical to the pre-blocking code -----------
      z          <- rnorm(n_par)
      state_prop <- state_curr + scale * as.numeric(L %*% z)
      names(state_prop) <- par_names

      lp_prop   <- lp_target(state_prop)$logpost
      log_alpha <- lp_prop - lp_curr
      if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        state_curr <- state_prop; lp_curr <- lp_prop
        trace_lp_curr <- if (!is.null(transform)) {
          lp_curr - transform$log_jacobian(state_curr)
        } else {
          lp_curr
        }
        n_accept <- n_accept + 1L; accepted[i] <- TRUE
      }
    } else {
      ## ---- Randomized blocking (Chib & Ramamurthy 2010 / Herbst &
      ## Schorfheide 2015 ch. 4): one random permutation of 1:n_par per
      ## iteration, split into n_blocks contiguous (as-equal-as-possible)
      ## chunks. Each block proposes a move of only its coordinates against
      ## the FULL posterior; the state evolves within the iteration.
      perm   <- sample.int(n_par)
      ## Contiguous, as-equal-as-possible split of 1:n_par into n_blocks
      ## chunks (the first `n_par %% n_blocks` chunks get one extra element).
      base_size <- n_par %/% n_blocks
      n_extra   <- n_par %% n_blocks
      block_sizes <- rep(base_size, n_blocks) + c(rep(1L, n_extra), rep(0L, n_blocks - n_extra))
      block_grp   <- rep.int(seq_len(n_blocks), block_sizes)
      blocks      <- split(perm, block_grp)

      any_accept <- FALSE
      for (idx in blocks) {
        nb  <- length(idx)
        L_b <- .robust_chol(Sigma_curr[idx, idx, drop = FALSE], nb)
        z   <- rnorm(nb)

        state_prop <- state_curr
        state_prop[idx] <- state_curr[idx] + scale * as.numeric(L_b %*% z)
        names(state_prop) <- par_names

        lp_prop   <- lp_target(state_prop)$logpost
        log_alpha <- lp_prop - lp_curr

        block_total_window <- block_total_window + 1L
        block_total_total  <- block_total_total + 1L
        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
          state_curr <- state_prop; lp_curr <- lp_prop
          trace_lp_curr <- if (!is.null(transform)) {
            lp_curr - transform$log_jacobian(state_curr)
          } else {
            lp_curr
          }
          any_accept <- TRUE
          block_accept_window <- block_accept_window + 1L
          block_accept_total  <- block_accept_total + 1L
        }
      }

      if (any_accept) { n_accept <- n_accept + 1L; accepted[i] <- TRUE }
    }

    stored_i <- if (!is.null(transform)) transform$to_constrained(state_curr) else state_curr
    if (ckpt) {
      buf_i <- buf_i + 1L
      buf[buf_i, ]  <- stored_i
      buf_lp[buf_i] <- trace_lp_curr
      if (buf_i >= flush_every || i == n_draws) {
        # Flush the buffer to disk and persist the restart state (state written
        # atomically AFTER the draws, so n_done never exceeds what is on disk).
        .ckpt_append_draws(ckpt_paths$draws, buf[seq_len(buf_i), , drop = FALSE])
        .ckpt_append_lp(ckpt_paths$lp, buf_lp[seq_len(buf_i)])
        .ckpt_save_state(ckpt_paths$state, list(
          state_curr = state_curr, lp_curr = lp_curr, trace_lp_curr = trace_lp_curr,
          scale = scale, Sigma_curr = Sigma_curr, n_done = i, n_burn = n_burn,
          n_accept = n_accept, n_draws_target = n_draws,
          rng = get(".Random.seed", envir = .GlobalEnv)))
        buf_i <- 0L
      }
    } else {
      chain[i, ]       <- stored_i
      logpost_trace[i] <- trace_lp_curr
    }

    if (adapt_cov && i <= n_burn) state_hist[i, ] <- state_curr

    # ---- Haario adaptive covariance update (burn-in only; frozen
    # thereafter so the retained draws satisfy detailed balance under a
    # fixed proposal -- no diminishing-adaptation condition is needed).
    if (adapt_cov && i %% adapt_every == 0 && i <= n_burn && i >= adapt_min_n) {
      Sigma_emp <- stats::cov(state_hist[1:i, , drop = FALSE]) +
        1e-8 * diag(n_par)
      Sigma_curr <- Sigma_emp
      L <- .robust_chol(Sigma_curr, n_par)
    }

    # ---- Scalar scale adaptation. With blocking, the any-block acceptance
    # rate saturates near 1 as n_blocks grows (at least one of n_blocks
    # proposals is very likely to be accepted), which would drive `scale`
    # to explode; use the mean PER-BLOCK acceptance rate instead.
    if (i %% adapt_every == 0 && i <= n_burn) {
      if (n_blocks == 1L) {
        recent_rate <- mean(accepted[max(1, i - adapt_every + 1):i])
      } else {
        recent_rate <- if (block_total_window > 0L)
          block_accept_window / block_total_window else 0
        block_accept_window <- 0L
        block_total_window  <- 0L
      }
      if (recent_rate > target_rate + 0.05) scale <- scale * 1.1
      else if (recent_rate < target_rate - 0.05) scale <- scale * 0.9
    }

    if (i %% 1000 == 0 || i == n_draws) {
      rate    <- n_accept / (i - 1)
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / i * (n_draws - i)
      label   <- if (is.null(chain_id)) "?" else as.character(chain_id)
      msg <- sprintf("Ch%s %d/%d accept=%.0f%% lp=%.1f scale=%.3f ETA=%.0fs",
                     label, i, n_draws, rate * 100, trace_lp_curr, scale, eta)
      if (!is.null(progressor)) progressor(message = msg, amount = 1)
      else if (verbose) .dynhr_cat("  ", msg, "\n")
    }
  }

  if (ckpt) {
    # Materialize the full chain from the streamed files for the return value.
    # checkpoint$return_chain = FALSE skips this for very long runs (the draws
    # remain on disk; read them with .ckpt_read_draws() or resume to extend).
    logpost_trace <- .ckpt_read_lp(ckpt_paths$lp)
    chain <- if (isFALSE(checkpoint$return_chain)) NULL else
      .ckpt_read_draws(ckpt_paths$draws, n_par, par_names)
  }
  post_chain   <- if (is.null(chain)) NULL else
    chain[(n_burn + 1):n_draws, , drop = FALSE]
  post_logpost <- logpost_trace[(n_burn + 1):n_draws]

  # ---- Mean per-block acceptance rate over the WHOLE run (for reporting);
  # NA when n_blocks == 1.
  block_accept_rate <- if (n_blocks == 1L) {
    NA_real_
  } else {
    block_accept_total / block_total_total
  }

  elapsed_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  list(
    chain             = post_chain,
    full_chain        = chain,
    logpost_trace     = logpost_trace,
    post_logpost      = post_logpost,
    acceptance_rate   = n_accept / (n_draws - 1),
    scale             = scale,
    n_draws           = n_draws,
    n_burn            = n_burn,
    block_accept_rate = block_accept_rate,
    checkpoint_dir    = if (ckpt) checkpoint$dir else NULL,
    elapsed_secs      = elapsed_secs,
    sampler           = "rwmh"
  )
}


## ---------------------------------------------------------------------------
## Correlated Pseudo-Marginal (CPM) RWMH sampler
##
## Implements Algorithm 1 of Deligiannidis, Doucet & Pitt (2018,
## J. Econometrics 206:799-829). At each step, the latent auxiliary variable
## U (the n_e x N standard-normal draws fed to the particle filter per period)
## is updated via an AR(1) (Crank-Nicolson) step jointly with theta, so that
## consecutive loglik evaluations are highly correlated. This dramatically
## reduces PMCMC loglik variance without changing the number of particles.
##
## Correctness: U is part of the JOINT state (theta, U). On acceptance BOTH
## theta and U move; on rejection BOTH stay. This ensures the chain targets
## the correct joint distribution whose theta-marginal is the posterior.
## ---------------------------------------------------------------------------

#' Correlated Pseudo-Marginal RWMH sampler (serial path only)
#'
#' @param log_post_fn  Function(theta, U_list = NULL) -> list(logpost, U_list).
#'   Must be the inner function returned by \code{make_log_posterior_tpf}.
#' @param theta0       Named initial parameter vector.
#' @param Sigma_prop   Proposal covariance (n_par x n_par): an ETA-SPACE
#'   covariance when \code{transform} is non-NULL, theta-space otherwise.
#' @param n_draws      Total draws (including burn-in).
#' @param n_burn       Burn-in draws to discard.
#' @param scale        Initial proposal scale factor.
#' @param target_rate  Target acceptance rate for adaptive scaling.
#' @param adapt_every  Adapt scale every N draws (during burn-in).
#' @param transform    Optional "dynhr_param_transform" (from
#'   \code{build_param_transform}), as in \code{rwmh()}: the random walk runs
#'   on \code{eta = to_unconstrained(theta)} with the target
#'   \code{logpost(to_constrained(eta)) + log_jacobian(eta)}; stored chain
#'   rows are \code{to_constrained(eta)} and \code{logpost_trace} /
#'   \code{post_logpost} the theta-space log-posterior. Before this
#'   argument existed, \code{run_posterior_estimation()}'s transform_params
#'   path handed the eta-space proposal covariance to a THETA-space walk.
#' @param rho_u        AR(1) correlation for U: U' = rho_u * U + sqrt(1 - rho_u^2) * Z.
#'   Range (0, 1). rho_u = 0 is independent (standard PMCMC); rho_u -> 1
#'   gives maximum correlation. Default 0.99.
#' @param verbose      Print progress messages.
#' @param progressor   Optional progressr function.
#' @param chain_id     Chain label for progress messages.
#' @return List: chain, full_chain, logpost_trace, post_logpost,
#'         acceptance_rate, scale, n_draws, n_burn, block_accept_rate.
#' @noRd
rwmh_cpm <- function(log_post_fn, theta0, Sigma_prop,
                     n_draws = 10000L, n_burn = 5000L,
                     scale = 0.5, target_rate = 0.25,
                     adapt_every = 100L,
                     rho_u = 0.99,
                     verbose = TRUE,
                     progressor = NULL, chain_id = NULL,
                     transform = NULL) {

  ## Input validation
  if (!is.numeric(rho_u) || length(rho_u) != 1L || !is.finite(rho_u) ||
      rho_u <= 0 || rho_u >= 1)
    stop("rwmh_cpm: rho_u must be a single numeric value in (0, 1).",
         " Got: ", rho_u, call. = FALSE)

  n_par     <- length(theta0)
  par_names <- names(theta0)

  L <- .robust_chol(Sigma_prop, n_par)

  chain         <- matrix(NA_real_, nrow = n_draws, ncol = n_par)
  colnames(chain) <- par_names
  logpost_trace <- numeric(n_draws)
  accepted      <- logical(n_draws)

  ## --- Sampling space: eta (transform) or theta ------------------------------
  ## The walk runs on `state`; `.cpm_eval()` returns the log_post_fn result
  ## with $logpost the SAMPLING-space target (+ log|dtheta/deta| under the
  ## transform) and $lp_theta the theta-space log-posterior for the trace.
  .cpm_theta <- function(state) {
    th <- if (is.null(transform)) state else transform$to_constrained(state)
    names(th) <- par_names
    th
  }
  .cpm_eval <- function(state, U_list) {
    r <- log_post_fn(.cpm_theta(state), U_list = U_list)
    r$lp_theta <- r$logpost
    if (!is.null(transform) && !is.null(r$logpost))
      r$logpost <- r$logpost + transform$log_jacobian(state)
    r
  }
  state0 <- if (is.null(transform)) theta0 else transform$to_unconstrained(theta0)
  names(state0) <- par_names

  ## --- Prime U_curr and lp_curr from theta0 --------------------------------
  res0   <- .cpm_eval(state0, U_list = NULL)
  lp_curr <- res0$logpost
  lp_theta_curr <- res0$lp_theta
  if (!is.finite(lp_curr))
    stop("rwmh_cpm: Initial parameter vector has -Inf log-posterior. ",
         "Check starting values.", call. = FALSE)
  U_curr <- res0$U_list

  if (is.null(U_curr))
    stop("rwmh_cpm: log_post_fn did not return U_list. ",
         "Use make_log_posterior_tpf (TPF likelihood only).", call. = FALSE)

  ## Infer dimensions from U_curr.
  ## U_curr layout (3T+1 entries; 2T+1 for the legacy layout):
  ##   [[1]]:            init normals (n_s x N matrix)
  ##   [[2]].[[T+1]]:    per-period shock normals (n_e x N matrix each)
  ##   [[T+2]].[[2T+1]]: per-period phi=1 resampling z scalars (length-1 numeric each)
  ##   [[2T+2]].[[3T+1]]: per-period mid-stage z K-vectors (length-K numeric each)
  ## Note: legacy callers may supply only T+1 or 2T+1 entries; those are
  ## handled gracefully (missing slots treated as uncorrelated / fresh draws).
  n_U   <- length(U_curr)  # 3T+1 (or 2T+1 or T+1 for legacy)
  N     <- ncol(U_curr[[1L]])
  rho_c <- sqrt(1 - rho_u^2)  # complementary coefficient

  ## Classify each U slot by type for the AR(1) update:
  ##   "matrix" — init normals or shock normals (nrow x N)
  ##   "scalar" — phi=1 resampling z (length-1 numeric)
  ##   "kvec"   — mid-stage z K-vector (length-K numeric, K > 1)
  U_type <- vapply(U_curr, function(u) {
    if (is.matrix(u))                              "matrix"
    else if (is.numeric(u) && length(u) == 1L)    "scalar"
    else                                           "kvec"    # K-vector
  }, character(1L))

  ## The U_list=NULL priming call returns NA in the scalar/kvec z slots (the
  ## filter drew its uniforms internally and did not record a z). NA would
  ## propagate through the AR(1) update, silently disabling sorted-CPM resampling
  ## for the whole chain — draw fresh z ~ N(0,1) now for each NA-filled slot.
  U_curr <- lapply(U_curr, function(u) {
    if (!is.matrix(u) && is.numeric(u) && any(is.na(u)))
      rnorm(length(u))
    else u
  })

  ## Pre-record per-slot dimensions for the AR(1) update (only for matrix slots)
  U_dims <- lapply(U_curr, function(U_t) {
    if (is.matrix(U_t)) c(nrow(U_t), N) else c(length(U_t), 1L)
  })

  state_curr <- state0
  chain[1, ]       <- .cpm_theta(state_curr)
  logpost_trace[1] <- lp_theta_curr
  accepted[1]      <- TRUE
  n_accept         <- 0L
  t_start          <- Sys.time()

  for (i in 2:n_draws) {

    ## -- Propose theta ---------------------------------------------------------
    z          <- rnorm(n_par)
    state_prop <- state_curr + scale * as.numeric(L %*% z)
    names(state_prop) <- par_names

    ## -- Propose U via AR(1) (pCN / autoregressive) update --------------------
    ## Matrix slots (init normals + shock normals): standard AR(1) on each element.
    ## Scalar slots (phi=1 resampling z): AR(1) on the scalar z ~ N(0,1) stored there;
    ##   the uniform used for resampling is u = pnorm(z), so correlating z gives
    ##   correlated resampling uniforms (Deligiannidis et al. sorted CPM).
    ## K-vector slots (mid-stage resampling z's): AR(1) on each element independently;
    ##   same N(0,1) marginal is preserved by the AR(1) update.
    U_prop <- mapply(function(U_t, dims, utype) {
      if (identical(utype, "matrix")) {
        nr <- dims[1L]
        rho_u * U_t + rho_c * matrix(rnorm(nr * N), nr, N)
      } else {
        ## Scalar z OR K-vector z: AR(1) update preserving N(0,1) marginal
        ## (length(U_t) == 1 for scalar, length(U_t) == K for kvec)
        rho_u * U_t + rho_c * rnorm(length(U_t))
      }
    }, U_curr, U_dims, U_type, SIMPLIFY = FALSE)

    ## -- Evaluate loglik at (theta_prop, U_prop) --------------------------------
    res_prop <- .cpm_eval(state_prop, U_list = U_prop)
    lp_prop  <- res_prop$logpost

    ## -- Joint MH accept/reject (theta AND U move together) --------------------
    log_alpha <- lp_prop - lp_curr
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      state_curr <- state_prop
      lp_curr    <- lp_prop
      lp_theta_curr <- res_prop$lp_theta
      U_curr     <- res_prop$U_list  # CRITICAL: accept U_prop (not U_curr!)
      n_accept   <- n_accept + 1L
      accepted[i] <- TRUE
    }
    ## On rejection: state_curr, lp_curr, U_curr all stay unchanged.
    ## U_prop is discarded. This is the correct CPM stationary distribution.

    chain[i, ]       <- .cpm_theta(state_curr)
    logpost_trace[i] <- lp_theta_curr

    ## -- Adaptive scale (during burn-in) ---------------------------------------
    if (i %% adapt_every == 0L && i <= n_burn) {
      recent_rate <- mean(accepted[max(1L, i - adapt_every + 1L):i])
      if (recent_rate > target_rate + 0.05) scale <- scale * 1.1
      else if (recent_rate < target_rate - 0.05) scale <- scale * 0.9
    }

    if (i %% 1000 == 0 || i == n_draws) {
      rate    <- n_accept / (i - 1)
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      eta     <- elapsed / i * (n_draws - i)
      label   <- if (is.null(chain_id)) "?" else as.character(chain_id)
      msg <- sprintf("CPM Ch%s %d/%d accept=%.0f%% lp=%.1f scale=%.3f rho_u=%.3f ETA=%.0fs",
                     label, i, n_draws, rate * 100, lp_theta_curr, scale, rho_u, eta)
      if (!is.null(progressor)) progressor(message = msg, amount = 1)
      else if (verbose) .dynhr_cat("  ", msg, "\n")
    }
  }

  post_chain   <- chain[(n_burn + 1):n_draws, , drop = FALSE]
  post_logpost <- logpost_trace[(n_burn + 1):n_draws]

  elapsed_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  list(
    chain             = post_chain,
    full_chain        = chain,
    logpost_trace     = logpost_trace,
    post_logpost      = post_logpost,
    acceptance_rate   = n_accept / (n_draws - 1),
    scale             = scale,
    n_draws           = n_draws,
    n_burn            = n_burn,
    block_accept_rate = NA_real_,
    elapsed_secs      = elapsed_secs,
    sampler           = "rwmh"
  )
}
