## R/hank-het-block2.R
## --------------------------------------------------------------------------
## Two-asset (liquid/illiquid) "het2 block": packages a steady-state two-asset
## household solution (hank_egm2_solve policies + Young's-method joint
## distribution) with the operations its sequence-space Jacobian and nonlinear
## transitions need. The two-asset counterpart of R/hank-het-block.R.
##
## Aggregate block OUTPUTS, all aggregated with the BEGINNING-of-period-t
## distribution D_t:
##   B_t   = sum_x D_t(x) b'(x)    aggregate liquid holdings
##   A_t   = sum_x D_t(x) a'(x)    aggregate illiquid holdings
##   C_t   = sum_x D_t(x) c(x)     aggregate consumption
##   CHI_t = sum_x D_t(x) Psi(x)   aggregate adjustment cost (a real resource
##                                 cost, so the GE resource constraint needs it)
## Block INPUTS are the aggregate paths {rb_t, ra_t, w_t}, plus -- when built
## with a Pi_fn/Pi_inputs pair -- named transition-probability inputs, exactly
## as in the one-asset block.
##
## CLASS: "hank_het2_block", which DELIBERATELY does not inherit
## "hank_het_block". Four existing gates type-check the one-asset class
## (hank_mixture_block_spec, validate_hank_block, hank_hh_hmm_loglik,
## hank_sam_reiter_linearize) and every one of them is one-asset-specific --
## a non-inheriting class makes each fail loudly instead of silently
## misinterpreting a 3-D policy as a 2-D one.
## --------------------------------------------------------------------------


#' Reject a two-asset block at a one-asset entry point, loudly
#'
#' Every one-asset routine assumes \code{n_e x n_a} policies and a single
#' marginal value.  Handed a \code{\link{hank_het2_block}} they all fail --
#' but they fail CRYPTICALLY ("argument is not a matrix", "non-conformable
#' arrays", or, worse, a \code{matrix()} recycling warning), because they trip
#' over a shape rather than a contract.  This turns that into one clear
#' message naming the two-asset equivalent, or saying plainly that the path is
#' not implemented for two assets.
#'
#' Not decorative.  The non-inheriting class protects the entry points that
#' type-check, but several do not: \code{hank_sam_reiter_linearize()} guards
#' with \code{!inherits(block, "hank_het_block") && is.null(block$Pi_fn)} --
#' an AND, so a two-asset block that HAS a \code{Pi_fn} sails past the class
#' test and dies downstream.  Reiter is exactly the path
#' \code{briefs/19-twoasset-hank-scope.md} section 7 defers, so a future caller
#' reaching for it must be told it does not exist rather than shown a shape
#' error.
#'
#' @param x The object handed to the caller.
#' @param caller Character: the calling function, for the message.
#' @param use Character or \code{NULL}: the two-asset equivalent to point at.
#'   \code{NULL} means no equivalent exists (the path is not implemented for
#'   two assets).
#' @return Invisibly \code{NULL}; called for the error.
#' @keywords internal
.hank_reject_het2 <- function(x, caller, use = NULL) {
  if (inherits(x, "hank_het2d_block"))
    stop(caller, "(): this is a DISCRETE-ADJUSTMENT two-asset block ",
         "(hank_het2d_block); use the het2d machinery (hank_het2d_jacobian, ",
         "hank_td2d_nonlinear, hank_het2d_block_spec).", call. = FALSE)
  if (!inherits(x, "hank_het2_block")) return(invisible(NULL))
  stop(caller, "(): this is a TWO-asset block (hank_het2_block), whose ",
       "policies are n_e x n_b x n_a arrays over a joint (e, b, a) cell space ",
       "and which carries two marginal values (Vb, Va). ",
       if (!is.null(use))
         paste0("Use ", use, "() instead.")
       else
         paste0("There is no two-asset equivalent: this path is not ",
                "implemented for two assets (see ",
                "briefs/19-twoasset-hank-scope.md). Use the sequence-space ",
                "route -- hank_het2_jacobian() -> hank_model() -> ",
                "hank_state_space() -- which is."),
       call. = FALSE)
}


#' Construct a two-asset heterogeneous-agent household block at steady state
#'
#' Solves the two-asset household problem (\code{\link{hank_egm2_solve}}) and
#' its stationary joint distribution at fixed aggregate prices
#' \code{(rb, ra, w)}, and stores everything the sequence-space Jacobian and
#' nonlinear-transition routines need.  The two-asset counterpart of
#' \code{\link{hank_het_block}}.
#'
#' With \code{ra > rb} and a nontrivial adjustment cost this block produces
#' \emph{wealthy hand-to-mouth} households -- constrained in the liquid asset
#' while holding substantial illiquid wealth -- which a one-asset block cannot
#' represent (there, constrained households are necessarily poor).
#'
#' @param b_grid Numeric: increasing LIQUID grid; \code{b_grid[1]} is the
#'   liquid floor.
#' @param a_grid Numeric: increasing ILLIQUID grid.
#' @param Pi Numeric \code{n_e x n_e}: row-stochastic income transition matrix.
#' @param e Numeric length-\code{n_e}: income levels (see
#'   \code{\link{hank_income_rouwenhorst}}).
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param rb,ra Steady-state liquid and illiquid returns.
#' @param w Steady-state wage (income is \code{y = w * e}).
#' @param chi0,chi1,chi2 Adjustment-cost parameters (see \code{\link{.hank_psi}};
#'   reference calibration \code{0.25 / 6.5 / 2}).
#' @param n_k,k_max Size and top of the multiplier grid for the
#'   liquid-constrained branch.
#' @param tol,maxit Passed to \code{\link{hank_egm2_solve}}.
#' @param dist_tol,dist_maxit Tolerance and iteration cap for the stationary
#'   distribution.
#' @param backend Character: \code{"cpp"} (default) or \code{"R"}.  Defaults to
#'   \code{getOption("dynhr.hank_backend", "cpp")}.  \code{Lambda} is always
#'   built via \code{\link{hank_forward_operator2}} for the returned object's
#'   contract, regardless of backend.
#' @param Vb_init,Va_init Optional \code{n_e x n_b x n_a} initial marginal
#'   values, passed straight through to \code{\link{hank_egm2_solve}}.
#'   \code{NULL} (default) uses that solver's own guess.
#' @param Tr Uniform lump-sum transfer added to every state's income
#'   (\code{y = w*e + Tr}); default \code{0} reproduces the transfer-free
#'   household byte-for-byte.
#' @param Pi_fn,Pi_inputs Optional transition-probability machinery, exactly as
#'   in \code{\link{hank_het_block}}: supplying them makes the named inputs
#'   perturbable aggregate inputs alongside \code{(rb, ra, w)}.  Names must not
#'   collide with \code{"rb"}/\code{"ra"}/\code{"w"}.
#' @param theta_coll Scalar LTV in \eqn{[0, 1)} on end-of-period illiquid
#'   collateral (D1): households may borrow down to
#'   \eqn{b' \ge b_{grid}[1] - \theta a'}.  Solved in the gap coordinate
#'   \eqn{x = b + \theta a} (see \code{\link{hank_egm2_solve}}), so the
#'   block's grids, distribution and \code{b}/\code{Vb} policies live in
#'   \eqn{x}, while \code{b_liq} and the aggregate \code{B} report the TRUE
#'   liquid position \eqn{b = x - \theta a}.  \code{theta_coll} is also a
#'   perturbable aggregate input (the macroprudential/LTV dial).
#' @param threads Worker threads for the compiled household solve, forwarded to
#'   \code{\link{hank_egm2_solve}}; \code{NULL} (default) resolves from
#'   \code{getOption("dynhr.hank3_threads")} then a machine default, and
#'   \code{1} forces the serial path.  Affects only the steady-state EGM solve,
#'   and does so bit-identically (see \code{\link{hank_egm2_solve}}); the
#'   distribution iteration and the sequence-space Jacobians are single-
#'   threaded, and \code{\link{hank_td2_nonlinear}} runs on the R reference
#'   step, which the compiled kernel's threading does not reach at all.
#'
#' @return An object of class \code{hank_het2_block} with the steady-state
#'   policies (\code{b}, \code{a}, \code{c}), per-cell adjustment cost
#'   \code{chi}, marginal values \code{Vb}/\code{Va} (each
#'   \code{n_e x n_b x n_a}), the joint distribution \code{D} (vector) and
#'   forward operator \code{Lambda}, aggregates \code{B}, \code{A}, \code{C},
#'   \code{CHI}, the calibration, and (when supplied) \code{Pi_fn}/\code{Pi_inputs}.
#' @seealso \code{\link{hank_het_block}} (one-asset),
#'   \code{\link{hank_egm2_solve}}, \code{\link{hank_forward_operator2}}
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' blk <- hank_het2_block(hank_asset_grid(10, 12, 0), hank_asset_grid(40, 12, 0),
#'                        inc$Pi, inc$e, beta = 0.95, eis = 0.5,
#'                        rb = 0.005, ra = 0.02, w = 1,
#'                        chi0 = 0.25, chi1 = 6.5, chi2 = 2, n_k = 10)
#' c(B = blk$B, A = blk$A, C = blk$C)
#' @export
hank_het2_block <- function(b_grid, a_grid, Pi, e, beta, eis, rb, ra, w,
                            chi0 = 0.25, chi1 = 6.5, chi2 = 2,
                            n_k = 50L, k_max = 1,
                            tol = 1e-10, maxit = 5000L,
                            dist_tol = 1e-13, dist_maxit = 200000L,
                            backend = getOption("dynhr.hank_backend", "cpp"),
                            Pi_fn = NULL, Pi_inputs = NULL, Tr = 0,
                            theta_coll = 0,
                            Vb_init = NULL, Va_init = NULL,
                            threads = NULL) {
  backend <- match.arg(backend, c("R", "cpp"))
  if (!is.numeric(Tr) || length(Tr) != 1L || !is.finite(Tr))
    stop("hank_het2_block: 'Tr' must be a finite numeric scalar (the uniform ",
         "lump-sum transfer; 0 restores the transfer-free household).")
  if (!is.numeric(e) || !all(is.finite(e)))
    stop("hank_het2_block: 'e' must be a finite numeric vector.")
  if (!is.numeric(w) || length(w) != 1L || !is.finite(w))
    stop("hank_het2_block: 'w' must be a finite numeric scalar.")
  ## Pi_fn / Pi_inputs contract -- same rules as hank_het_block, with the
  ## two-asset price names.
  if (xor(is.null(Pi_fn), is.null(Pi_inputs)))
    stop("hank_het2_block: supply Pi_fn and Pi_inputs together (or neither).")
  if (!is.null(Pi_fn)) {
    if (!is.function(Pi_fn))
      stop("hank_het2_block: Pi_fn must be a function.")
    nm <- names(Pi_inputs)
    if (!is.list(Pi_inputs) || length(Pi_inputs) == 0L ||
        is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
      stop("hank_het2_block: Pi_inputs must be a non-empty, fully and ",
           "uniquely named list of steady-state transition-input values.")
    if (any(nm %in% c("rb", "ra", "w", "Tr", "theta_coll")))
      stop("hank_het2_block: Pi_inputs names must not collide with the ",
           "aggregate inputs ('rb', 'ra', 'w', 'Tr', 'theta_coll').")
    Pi_check <- do.call(Pi_fn, Pi_inputs)
    if (!isTRUE(all.equal(unname(as.matrix(Pi_check)), unname(as.matrix(Pi)),
                          tolerance = 1e-10)))
      stop("hank_het2_block: Pi_fn evaluated at Pi_inputs does not reproduce ",
           "the steady-state Pi (the transition inputs and Pi are ",
           "inconsistent).")
  }

  hh <- hank_egm2_solve(b_grid, a_grid, y = w * e + Tr, rb = rb, ra = ra,
                        beta = beta, eis = eis, chi0 = chi0, chi1 = chi1,
                        chi2 = chi2, Pi = Pi, n_k = n_k, k_max = k_max,
                        tol = tol, maxit = maxit, backend = backend,
                        theta_coll = theta_coll,
                        Vb_init = Vb_init, Va_init = Va_init,
                        threads = threads)
  if (!hh$converged)
    warning("hank_het2_block: household EGM did not converge at steady state")

  ## The distribution lives on the SOLVER's state coordinates -- under
  ## collateral that is the gap x = b + theta*a (hh$b is the x-policy), so the
  ## forward operator takes hh$b as-is.  The LIQUID aggregate, however, must
  ## be the TRUE position: B = E[b'] = E[x' - theta*a'] via hh$b_liq (which IS
  ## hh$b at theta_coll = 0).
  ## elapsed_dist times the DISTRIBUTION only -- the power iteration, not the
  ## Lambda build, matching how the three-asset block splits its two windows.
  ## Reporting one number for both would make the two families' manifest rows
  ## silently non-comparable.
  Lam <- hank_forward_operator2(hh$b, hh$a, b_grid, a_grid, Pi)
  .td0 <- proc.time()[["elapsed"]]
  sd  <- hank_stationary_dist(Lam, tol = dist_tol, maxit = dist_maxit,
                              backend = backend)
  .elapsed_dist <- proc.time()[["elapsed"]] - .td0
  D   <- sd$d
  structure(
    list(b_grid = b_grid, a_grid = a_grid, k_grid = hh$k_grid,
         Pi = Pi, e = e, beta = beta, eis = eis,
         rb = rb, ra = ra, w = w, Tr = Tr,
         theta_coll = theta_coll,
         chi0 = chi0, chi1 = chi1, chi2 = chi2,
         b = hh$b, a = hh$a, c = hh$c, chi = hh$chi, b_liq = hh$b_liq,
         Vb = hh$Vb, Va = hh$Va,
         Lambda = Lam, D = D,
         B = hank_aggregate2(D, hh$b_liq),
         A = hank_aggregate2(D, hh$a),
         C = hank_aggregate2(D, hh$c),
         CHI = hank_aggregate2(D, hh$chi),
         n_e = length(e), n_b = length(b_grid), n_a = length(a_grid),
         n_k = length(hh$k_grid),
         dist_converged = sd$converged,
         ## Run metadata hank_het2_manifest() reads, mirroring the one- and
         ## three-asset blocks so a driver can rbind manifest rows across
         ## families. All of it comes from the solve that actually ran -- none
         ## is re-derived, and last_value_gap stays NA because the two-asset
         ## solver's convergence test is on the POLICY gap, so reporting a
         ## value gap here would be inventing a quantity it never computed.
         backend = hh$backend, threads = hh$threads,
         iterations = hh$iterations, converged = hh$converged,
         last_policy_gap = hh$last_policy_gap,
         last_value_gap = NA_real_,
         elapsed_solve = hh$elapsed, elapsed_dist = .elapsed_dist,
         Pi_fn = Pi_fn, Pi_inputs = Pi_inputs),
    class = c("hank_het2_block", "hank_block"))
}


#' One nonlinear two-asset EGM backward step at given aggregate prices
#'
#' Thin wrapper mapping block inputs \code{(rb, ra, w)} to the two-asset EGM
#' step; the two-asset counterpart of \code{\link{.hank_block_step}}, and the
#' primitive the fake-news sweep and nonlinear transition are built from.
#'
#' Note the \code{Psi1_grid} argument depends on \code{ra} (not only on
#' \code{a_grid}), so it must be rebuilt whenever \code{ra} is perturbed --
#' passing a stale grid would silently mis-evaluate the illiquid FOC.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param Vb_p,Va_p Next-period marginal values (\code{n_e x n_b x n_a}).
#' @param rb,ra,w Aggregate returns and wage this period.
#' @param Pi Optional income transition matrix overriding \code{block$Pi} for
#'   this period (time-varying transition probabilities).
#' @param Tr Optional lump-sum transfer overriding the block's steady-state one.
#' @param theta_coll,dtheta_next Collateral LTV and its one-period increment;
#'   see \code{\link{.hank_egm2_step}}.
#' @param backend Character: \code{"cpp"} (default, from
#'   \code{getOption("dynhr.hank_backend")}) or \code{"R"}.  The pure-R
#'   \code{\link{.hank_egm2_step}} remains the REFERENCE implementation and the
#'   other half of the exact cross-oracle pair pinned by
#'   \code{test-hank-egm2-cpp-parity.R}; it is selected here (and by the option
#'   globally) rather than deleted.  The compiled kernel does not implement the
#'   collateral gap coordinate, so \code{theta_coll != 0} or
#'   \code{dtheta_next != 0} FORCES the R path regardless of this argument --
#'   silently dropping collateral would be a wrong answer, not a slow one.
#' @param threads Worker threads for the compiled backend, or \code{NULL} to
#'   resolve via \code{\link{hank_resolve_threads}}.  Callers that iterate the
#'   step (the nonlinear transition, the fake-news sweep) should resolve ONCE
#'   and pass the integer down.  The kernel is bit-identical across thread
#'   counts by construction, so this never changes the answer.
#' @return List with \code{Vb}, \code{Va}, \code{b}, \code{a}, \code{c},
#'   \code{chi}.
#' @keywords internal
.hank_block_step2 <- function(block, Vb_p, Va_p, rb, ra, w, Pi = NULL,
                              Tr = NULL, theta_coll = NULL, dtheta_next = 0,
                              backend = getOption("dynhr.hank_backend", "cpp"),
                              threads = NULL) {
  backend <- match.arg(backend, c("R", "cpp"))
  if (is.null(Pi)) Pi <- block$Pi
  if (is.null(Tr)) Tr <- .hank_block_tr(block)
  if (is.null(theta_coll))
    theta_coll <- if (is.null(block$theta_coll)) 0 else block$theta_coll
  n_a <- block$n_a
  ## Psi1 depends on ra, so it cannot be cached on the block across an
  ## ra perturbation.
  Psi1 <- .hank_psi(matrix(block$a_grid, n_a, n_a),
                    matrix(block$a_grid, n_a, n_a, byrow = TRUE),
                    ra, block$chi0, block$chi1, block$chi2)$Psi1
  y <- w * block$e + Tr
  ## Collateral is R-only (the compiled kernel has no theta_coll/dtheta_next
  ## arguments at all), so route it there rather than silently ignoring it.
  if (backend == "cpp" && theta_coll == 0 && dtheta_next == 0)
    return(hank_egm2_step_cpp(Vb_p, Va_p, block$b_grid, block$a_grid,
                              block$k_grid, y, rb, ra, block$beta, block$eis,
                              block$chi0, block$chi1, block$chi2,
                              as.matrix(Pi), Psi1,
                              hank_resolve_threads(threads)))
  .hank_egm2_step(Vb_p, Va_p, block$b_grid, block$a_grid, block$k_grid,
                  y = y, rb = rb, ra = ra, beta = block$beta,
                  eis = block$eis, chi0 = block$chi0, chi1 = block$chi1,
                  chi2 = block$chi2, Pi = Pi, Psi1_grid = Psi1,
                  theta_coll = theta_coll, dtheta_next = dtheta_next)
}


#' Nonlinear perfect-foresight transition of a two-asset het block
#'
#' Given aggregate input PATHS over horizon \code{T_h} (with terminal
#' conditions returning to steady state), computes the aggregate output paths
#' by a backward household solve carrying the \code{(Vb, Va)} pair followed by
#' a forward distribution simulation.  The block's nonlinear map from input
#' paths to output paths, and the brute-force reference used to validate the
#' sequence-space Jacobian.  The two-asset counterpart of
#' \code{\link{hank_td_nonlinear}}.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param rb_path,ra_path,w_path Numeric length-\code{T_h} input paths (levels).
#'   Missing entries default to the steady-state value.
#' @param T_h Integer horizon (default: the length of the longest supplied path).
#' @param keep_policies Logical (default \code{FALSE}). When \code{TRUE}, also
#'   return the per-period \code{b_pol}/\code{a_pol}/\code{c_pol} lists and the
#'   per-period sparse \code{Lambda} list (same no-transpose convention as
#'   \code{block$Lambda}).
#' @param Tr_path Optional length-\code{T_h} LEVEL path of the uniform
#'   lump-sum transfer; default holds it at the block's steady-state
#'   \code{Tr}.
#' @param pi_input_paths Optional named list of length-\code{T_h} LEVEL paths
#'   for (a subset of) the block's transition-probability inputs.  Requires a
#'   block built with \code{Pi_fn}/\code{Pi_inputs}.
#' @param D0 Optional initial (beginning-of-period-1) distribution, a
#'   length-\code{n_e*n_b*n_a} nonnegative vector summing to 1 in the package's
#'   two-asset order.  Default \code{NULL} starts from the steady-state
#'   \code{block$D}; supply it for state-dependence experiments.
#' @param theta_path Optional length-\code{T_h} LEVEL path of the collateral
#'   LTV \code{theta_coll}; default holds it at the block's steady-state
#'   value.  A moving path re-bases the gap coordinate each period: the
#'   backward steps receive \eqn{\theta_{t+1}-\theta_t} to pre-shift the
#'   continuation marginals, and the forward push shifts the liquid policy by
#'   the same amount (arrival \eqn{x_{t+1} = x' + (\theta_{t+1}-\theta_t)a'});
#'   the reported \code{B} is the true liquid aggregate
#'   \eqn{E[x' - \theta_t a']} throughout.
#' @param backend Character: \code{"cpp"} (default, from
#'   \code{getOption("dynhr.hank_backend")}) or \code{"R"} for the pure-R
#'   reference backward step.  A non-zero \code{theta_path} forces \code{"R"}
#'   (the compiled kernel has no collateral coordinate).
#' @param threads Worker threads for the compiled backward step, or \code{NULL}
#'   (default) to resolve from \code{getOption("dynhr.hank_threads")} and then a
#'   machine default (see \code{\link{hank_resolve_threads}}); \code{1} forces
#'   the serial kernel.  The kernel is bit-identical across thread counts, so
#'   this is a speed knob only.  Ignored by \code{backend = "R"}.
#'
#' @return A list with numeric length-\code{T_h} paths \code{B}, \code{A},
#'   \code{C}, \code{CHI}, and \code{Dpath} (\code{n_cell x T_h}): column
#'   \code{t} is the beginning-of-period-\code{t} distribution over which the
#'   date-\code{t} aggregates are taken.
#' @seealso \code{\link{hank_td_nonlinear}} (one-asset),
#'   \code{\link{hank_het2_block}}
#' @export
hank_td2_nonlinear <- function(block, rb_path = NULL, ra_path = NULL,
                               w_path = NULL, T_h = NULL,
                               keep_policies = FALSE, pi_input_paths = NULL,
                               D0 = NULL, Tr_path = NULL,
                               theta_path = NULL,
                               backend = getOption("dynhr.hank_backend", "cpp"),
                               threads = NULL) {
  backend <- match.arg(backend, c("R", "cpp"))
  ## Resolve ONCE for the whole backward pass: the resolver hits
  ## parallel::detectCores() and the once-per-regime reporter, neither of which
  ## belongs inside a T_h-long loop.
  threads <- hank_resolve_threads(threads)
  if (!inherits(block, "hank_het2_block"))
    stop("hank_td2_nonlinear: 'block' must be a hank_het2_block.")
  if (!is.null(D0)) {
    D0 <- as.numeric(D0)
    if (length(D0) != length(block$D) || any(D0 < -1e-12) ||
        abs(sum(D0) - 1) > 1e-8)
      stop("hank_td2_nonlinear: D0 must be a nonnegative length-",
           length(block$D), " distribution summing to 1 (package two-asset ",
           "order: illiquid index fastest).", call. = FALSE)
  }
  if (is.null(T_h)) {
    lens <- c(length(rb_path), length(ra_path), length(w_path),
              length(Tr_path), length(theta_path),
              if (!is.null(pi_input_paths))
                vapply(pi_input_paths, length, integer(1)))
    T_h <- max(lens, if (all(lens == 0L)) 1L else 0L)
  }
  if (is.null(rb_path)) rb_path <- rep(block$rb, T_h)
  if (is.null(ra_path)) ra_path <- rep(block$ra, T_h)
  if (is.null(w_path))  w_path  <- rep(block$w,  T_h)
  if (is.null(Tr_path)) Tr_path <- rep(.hank_block_tr(block), T_h)
  if (is.null(theta_path))
    theta_path <- rep(if (is.null(block$theta_coll)) 0
                      else block$theta_coll, T_h)
  stopifnot(length(rb_path) == T_h, length(ra_path) == T_h,
            length(w_path) == T_h, length(Tr_path) == T_h,
            length(theta_path) == T_h)
  Pi_path <- .hank2_pi_path(block, pi_input_paths, T_h)
  Pi_at   <- function(t) if (is.null(Pi_path)) block$Pi else Pi_path[[t]]

  ## Backward: terminal (Vb, Va) at steady state.  Under a theta PATH the gap
  ## coordinate x = b + theta_t*a is defined with the CONTEMPORANEOUS theta,
  ## so each step also receives dtheta_next = theta_{t+1} - theta_t (steady
  ## state beyond the horizon) to pre-shift the continuation marginals; see
  ## .hank_egm2_step.  All-zero at a constant path: bit-exact back-compat.
  th_ss <- if (is.null(block$theta_coll)) 0 else block$theta_coll
  dth   <- c(diff(theta_path), th_ss - theta_path[T_h])
  b_pol <- vector("list", T_h); a_pol <- vector("list", T_h)
  c_pol <- vector("list", T_h); chi_p <- vector("list", T_h)
  Vb <- block$Vb; Va <- block$Va
  for (t in T_h:1L) {
    step <- .hank_block_step2(block, Vb, Va, rb_path[t], ra_path[t],
                              w_path[t], Pi = Pi_at(t), Tr = Tr_path[t],
                              theta_coll = theta_path[t],
                              dtheta_next = dth[t],
                              backend = backend, threads = threads)
    b_pol[[t]] <- step$b; a_pol[[t]] <- step$a
    c_pol[[t]] <- step$c; chi_p[[t]] <- step$chi
    Vb <- step$Vb; Va <- step$Va
  }

  ## Forward: beginning distribution D_1 = D_ss (or the caller's D0).
  ## Collateral bookkeeping mirrors the backward pass: the liquid AGGREGATE is
  ## the true position b' = x' - theta_t*a', while the mass chosen at x' this
  ## period ARRIVES at x' + dtheta_next*a' tomorrow (the coordinate re-basing),
  ## so the operator gets the shifted liquid policy.  Both are exact no-ops on
  ## a constant path.
  coll <- any(theta_path != 0) || th_ss != 0
  D <- if (is.null(D0)) block$D else D0
  ## Date-1 coordinate re-basing: D (default block$D, or a caller D0) is the
  ## PHYSICAL beginning-of-period-1 distribution expressed in STEADY-STATE
  ## coordinates x = b + theta_ss*a.  If theta_1 differs (an unanticipated LTV
  ## move at date 1), the same physical states sit at x = b + theta_1*a --
  ## every collateral holder shifts by (theta_1 - theta_ss)*a.  Omitting this
  ## silently confiscates wealth from them (measured: it flipped the impact
  ## sign of J[C][theta]).  Same device as the Fisher-channel D0 shift.
  if (coll && theta_path[1L] != th_ss)
    D <- .hank2_shift_x(D, theta_path[1L] - th_ss, block)
  B <- numeric(T_h); A <- numeric(T_h); C <- numeric(T_h); CHI <- numeric(T_h)
  Dpath <- matrix(0, length(D), T_h)
  Lam_keep <- if (keep_policies) vector("list", T_h) else NULL
  for (t in seq_len(T_h)) {
    Dpath[, t] <- D
    B[t]   <- hank_aggregate2(D, if (coll) b_pol[[t]] - theta_path[t] * a_pol[[t]]
                                 else b_pol[[t]])
    A[t]   <- hank_aggregate2(D, a_pol[[t]])
    C[t]   <- hank_aggregate2(D, c_pol[[t]])
    CHI[t] <- hank_aggregate2(D, chi_p[[t]])
    Lam <- hank_forward_operator2(if (coll && dth[t] != 0)
                                    b_pol[[t]] + dth[t] * a_pol[[t]]
                                  else b_pol[[t]],
                                  a_pol[[t]], block$b_grid,
                                  block$a_grid, Pi_at(t))
    if (keep_policies) Lam_keep[[t]] <- Lam
    D <- as.numeric(Matrix::t(Lam) %*% D)
  }
  out <- list(B = B, A = A, C = C, CHI = CHI, Dpath = Dpath)
  if (keep_policies) {
    out$b_pol <- b_pol; out$a_pol <- a_pol; out$c_pol <- c_pol
    out$chi_pol <- chi_p; out$Lambda <- Lam_keep
  }
  out
}


#' Shift a two-asset distribution along the liquid axis by delta * a
#'
#' The collateral coordinate re-basing (D1, brief 19 section 9.12): mass at
#' \eqn{(e, x, a)} moves to \eqn{(e, x + \delta a, a)}, per a-slice, via the
#' VALIDATED joint-lottery forward operator with an identity income transition
#' (the a-component queries grid knots exactly, so only the liquid lottery is
#' live; edge mass clamps -- the forced-deleveraging convention).  Used to
#' re-base the initial distribution when \eqn{\theta_1 \ne \theta_{ss}} and to
#' build the distribution derivative in the theta column of the fake-news
#' Jacobian.
#' @keywords internal
.hank2_shift_x <- function(D, delta, block) {
  Xq <- .hank2_bcast_mid(block$b_grid, block$n_e, block$n_b, block$n_a) +
    delta * .hank2_bcast_a(block$a_grid, block$n_e, block$n_b, block$n_a)
  Aq <- .hank2_bcast_a(block$a_grid, block$n_e, block$n_b, block$n_a)
  S <- hank_forward_operator2(Xq, Aq, block$b_grid, block$a_grid,
                              diag(block$n_e))
  as.numeric(Matrix::t(S) %*% D)
}


#' Per-period transition-matrix path for a two-asset block
#'
#' Two-asset sibling of \code{\link{.hank_pi_path}} (kept separate rather than
#' generalising the validated one-asset helper).
#' @keywords internal
.hank2_pi_path <- function(block, pi_input_paths, T_h) {
  if (is.null(pi_input_paths) || length(pi_input_paths) == 0L) return(NULL)
  if (is.null(block$Pi_fn))
    stop(".hank2_pi_path: transition-input paths supplied, but the block has ",
         "no Pi_fn/Pi_inputs (rebuild it via hank_het2_block(..., Pi_fn, ",
         "Pi_inputs)).")
  nm  <- names(pi_input_paths)
  bad <- setdiff(nm, names(block$Pi_inputs))
  if (length(bad))
    stop(".hank2_pi_path: unknown transition input(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; the block's Pi_inputs are ",
         paste0("'", names(block$Pi_inputs), "'", collapse = ", "), ".")
  ok <- vapply(pi_input_paths, function(p) length(p) == T_h, logical(1))
  if (!all(ok))
    stop(".hank2_pi_path: every transition-input path must have length T_h (",
         T_h, ").")
  lapply(seq_len(T_h), function(t) {
    args <- block$Pi_inputs
    for (k in nm) args[[k]] <- pi_input_paths[[k]][t]
    do.call(block$Pi_fn, args)
  })
}


#' Decompose a two-asset consumption response into cash-flow and revaluation
#'
#' The anti-conflation instrument (briefs/19 sections 9.6-9.7): runs the SAME
#' household block through \code{\link{hank_td2_nonlinear}} on three
#' illiquid-return paths -- the FULL path, a caller-supplied CASH-FLOW-ONLY
#' counterfactual (price frozen; e.g. the dividend-yield component for equity,
#' or a constant for a geometric bond, whose coupons are fixed so its entire
#' return movement is revaluation), and the implied REVALUATION-ONLY path
#' \code{ra_ss + (ra_full - ra_cashflow)} -- and reports the channel
#' contributions under names that cannot silently merge.  The nonlinear
#' channels need not add up exactly; the \code{interaction} residual is
#' reported rather than hidden.
#'
#' Incidence: the impact-date consumption change per household is aggregated
#' within mass-weighted illiquid-wealth groups, giving the
#' \code{groups x channels} incidence table that answers "WHO is affected
#' through WHICH channel" -- the question conflating wealth with cash-flow
#' gets wrong.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param ra_full Numeric length-\code{T_h}: the full illiquid-return path
#'   (cash flow + revaluation), in levels.
#' @param ra_cashflow Numeric length-\code{T_h}: the cash-flow-only
#'   counterfactual, in levels.  For a \code{\link{hank_bond_block}} asset this
#'   is \code{rep(block$ra, T_h)}: fixed coupons mean NO cash-flow channel.
#' @param rb_path,w_path,Tr_path Optional other input paths (levels), shared
#'   across all three runs; defaults hold them at steady state.
#' @param n_groups Integer: number of mass-weighted illiquid-wealth groups for
#'   the incidence table.
#'
#' @return A list with \code{paths} (a data frame: \code{C_full},
#'   \code{C_cashflow}, \code{C_revaluation}, \code{C_interaction}, and the
#'   same for \code{A} and \code{B}, all as DEVIATIONS from steady state) and
#'   \code{incidence} (impact-date per-household consumption change by
#'   illiquid group and channel, with the group boundaries as attributes).
#' @seealso \code{\link{hank_bond_block}}, \code{\link{hank_td2_nonlinear}}
#' @export
hank_td2_reval_decompose <- function(block, ra_full, ra_cashflow,
                                     rb_path = NULL, w_path = NULL,
                                     Tr_path = NULL, n_groups = 3L) {
  if (!inherits(block, "hank_het2_block"))
    stop("hank_td2_reval_decompose: 'block' must be a hank_het2_block.")
  T_h <- length(ra_full)
  if (length(ra_cashflow) != T_h)
    stop("hank_td2_reval_decompose: 'ra_full' and 'ra_cashflow' must have ",
         "the same length.")
  ra_reval <- block$ra + (ra_full - ra_cashflow)

  run <- function(ra_p)
    hank_td2_nonlinear(block, rb_path = rb_path, ra_path = ra_p,
                       w_path = w_path, T_h = T_h, Tr_path = Tr_path,
                       keep_policies = TRUE)
  td_f <- run(ra_full); td_c <- run(ra_cashflow); td_r <- run(ra_reval)

  dev <- function(td, f) td[[f]] - block[[f]]
  paths <- data.frame(t = seq_len(T_h))
  for (f in c("C", "A", "B")) {
    paths[[paste0(f, "_full")]]        <- dev(td_f, f)
    paths[[paste0(f, "_cashflow")]]    <- dev(td_c, f)
    paths[[paste0(f, "_revaluation")]] <- dev(td_r, f)
    paths[[paste0(f, "_interaction")]] <-
      dev(td_f, f) - dev(td_c, f) - dev(td_r, f)
  }

  ## Impact-date incidence by mass-weighted illiquid-wealth group.
  n_e <- block$n_e; n_b <- block$n_b; n_a <- block$n_a
  A_cell <- .hank2_arr_to_vec(.hank2_bcast_a(block$a_grid, n_e, n_b, n_a))
  D <- block$D
  ## group boundaries: mass-weighted quantiles of the illiquid position
  o <- order(A_cell); cw <- cumsum(D[o])
  cuts <- vapply(seq_len(n_groups - 1L) / n_groups,
                 function(qq) A_cell[o][which(cw >= qq)[1L]], numeric(1))
  cuts <- unique(c(-Inf, cuts, Inf))
  gid <- cut(A_cell, cuts, labels = FALSE)
  inc <- matrix(0, length(unique(gid)), 3L,
                dimnames = list(paste0("illiq_g", sort(unique(gid))),
                                c("full", "cashflow", "revaluation")))
  for (ch in seq_len(3L)) {
    td <- list(td_f, td_c, td_r)[[ch]]
    dc <- .hank2_arr_to_vec(td$c_pol[[1L]] - block$c)
    for (g in sort(unique(gid))) {
      m <- gid == g
      inc[g, ch] <- sum(D[m] * dc[m]) / sum(D[m])
    }
  }
  attr(inc, "group_cuts") <- cuts
  list(paths = paths, incidence = inc)
}
