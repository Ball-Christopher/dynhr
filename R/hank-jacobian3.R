## The shared input/output vocabulary of BOTH three-asset Jacobians. Defined
## here, above the first roxygen block, because anything placed between a
## roxygen block and its function silently steals that documentation (and its
## @export) -- which is how hank_het3_jacobian_nd briefly lost its NAMESPACE
## entry while this was being written.
.hank3_jac_prices  <- c("rd","rf","ra","w","px","Tr")
.hank3_jac_outputs <- c("D","F","A","C","CHI","PHI")

## The DEFAULT input set of both three-asset Jacobians: the prices, plus any
## transition-probability inputs the block was built with (hank_het3_block's
## Pi_fn/Pi_inputs -- e.g. the job-finding rate f and separation rate s of
## hank_employment_income). Kept a function of the block rather than a
## constant so a Pi-carrying block gets its f/s columns by DEFAULT: a caller
## who built the block with those inputs asked for them, and a default that
## silently dropped them is how a household ends up estimated as if
## employment risk were fixed.
.hank3_jac_inputs <- function(block)
  c(.hank3_jac_prices, names(block$Pi_inputs))

## Is `input` a transition-probability input of this block (as opposed to a
## price)? These columns need the joint (policy, Pi) forward derivative.
.hank3_is_pi_input <- function(block, input)
  !is.null(block$Pi_inputs) && input %in% names(block$Pi_inputs)

#' Validate requested three-asset Jacobian inputs
#'
#' Shared by \code{\link{hank_het3_jacobian}} and
#' \code{\link{hank_het3_jacobian_nd}} so both routes accept exactly the same
#' vocabulary. Counterpart of \code{.hank_het_check_inputs} for the one-asset
#' block.
#' @param block A \code{\link{hank_het3_block}}.
#' @param inputs Character vector of requested inputs.
#' @param caller Name of the calling function, for error messages.
#' @return \code{inputs}, validated.
#' @keywords internal
.hank3_check_inputs <- function(block, inputs, caller) {
  allowed <- .hank3_jac_inputs(block)
  bad <- setdiff(inputs, allowed)
  if (length(bad))
    stop(caller, ": unsupported input(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; this block supports ",
         paste0("'", allowed, "'", collapse = ", "),
         if (is.null(block$Pi_inputs))
           " (build the block with Pi_fn/Pi_inputs to add transition-probability inputs)"
         else "", ".")
  inputs
}

# The singleton-f_grid reduction pins f' = f = 0 and delegates to
# hank_egm2_solve, which has never heard of px. Every output's derivative with
# respect to px is therefore EXACTLY zero -- px multiplies zero on both sides
# of the budget -- by the same argument that makes the rf column exactly zero
# there, which test-hank-egm3.R already asserts.
#
# The zero is returned rather than computed because hank_egm3_solve REFUSES
# px != 1 on this path: that guard exists so a user cannot believe they are
# using the valuation channel when nothing is being valued, and it is right to
# keep it, but an internal FD sweep must not trip over it. Erroring instead
# would also make `px` in the default `inputs` break every reduction block.
.hank3_px_is_inert <- function(block) length(block$f_grid) == 1L


#' Numerical-differentiation Jacobian for the three-asset reference block
#'
#' Brute-force sequence-space Jacobian: perturbs each aggregate input at each
#' date and re-runs the full nonlinear transition
#' (\code{\link{hank_td3_nonlinear}}), central-differencing the aggregate
#' outputs. Correct but \eqn{O(T)} solves per (input, date); this is the
#' mandatory acceptance oracle for the fake-news Jacobian
#' \code{\link{hank_het3_jacobian}}, the role \code{\link{hank_het_jacobian_nd}}
#' plays for the one-asset block.
#'
#' @param block A \code{\link{hank_het3_block}}.
#' @param T_h Integer horizon.
#' @param inputs Character subset of \code{c("rd", "rf", "ra", "w", "px", "Tr")}
#'   plus, for a block built with \code{Pi_fn}/\code{Pi_inputs}, that block's
#'   named transition-probability inputs (\code{names(block$Pi_inputs)});
#'   \code{NULL} (the default) takes all of them. Transition inputs are swept
#'   as \code{pi_input_paths} of \code{\link{hank_td3_nonlinear}}, prices as
#'   price paths -- both through the same nonlinear transition, so this stays
#'   a genuinely independent oracle for either kind of column.
#' @param outputs Character subset of \code{c("D", "F", "A", "C", "CHI",
#'   "PHI")} (aggregate liquid, foreign, capital holdings, consumption, and
#'   capital/foreign adjustment resources).
#' @param delta Numeric FD step for the input perturbation.
#' @param threads Worker threads for the backward policy step, forwarded to
#'   every \code{\link{hank_td3_nonlinear}} call this sweep makes.
#'   \code{NULL} (the default) resolves via
#'   \code{getOption("dynhr.hank3_threads")} and then a machine-derived
#'   default, exactly as \code{\link{hank_het3_jacobian}} does; \code{1}
#'   forces the serial path. Output is bit-identical at every thread count, so
#'   this only affects wall time. Note that \code{R CMD check} sets
#'   \code{_R_CHECK_LIMIT_CORES_}, which clamps the resolved count to 2 --
#'   set \code{options(dynhr.hank_report_threads = TRUE)} if you need to see
#'   what a long batch actually resolved to.
#'
#' @return Nested list \code{J[[output]][[input]]}, each a \code{T_h x T_h}
#'   matrix with \code{[t, s] = dO_t/dI_s} (central difference).
#' @seealso \code{\link{hank_het3_jacobian}} (the fake-news Jacobian this
#'   validates), \code{\link{hank_het_jacobian_nd}} (one-asset),
#'   \code{\link{hank_td3_nonlinear}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' J <- hank_het3_jacobian_nd(blk, T_h = 2, delta = 3e-6)
#' dim(J$D$rd)
#' @keywords internal
#' @export
hank_het3_jacobian_nd <- function(block,T_h,inputs=NULL,outputs=.hank3_jac_outputs,delta=1e-5,
                                  threads=NULL) {
  if(!inherits(block,"hank_het3_block"))stop("hank_het3_jacobian_nd: block must be hank_het3_block")
  if(!is.numeric(T_h)||length(T_h)!=1L||T_h<1||!is.finite(T_h))stop("hank_het3_jacobian_nd: T_h must be positive")
  if(is.null(inputs))inputs<-.hank3_jac_inputs(block)
  inputs<-.hank3_check_inputs(block,inputs,"hank_het3_jacobian_nd")
  bad<-setdiff(outputs,.hank3_jac_outputs);if(length(bad))stop("unsupported output(s): ",paste(bad,collapse=", "))
  base<-list(rd=rep(block$rd,T_h),rf=rep(block$rf,T_h),ra=rep(block$ra,T_h),w=rep(block$w,T_h),px=rep(if(is.null(block$px))1 else block$px,T_h),Tr=rep(.hank_block_tr(block),T_h));J<-setNames(lapply(outputs,function(o)setNames(lapply(inputs,function(i)matrix(0,T_h,T_h)),inputs)),outputs)
  # Name the path arguments in FULL. These used to reach hank_td3_nonlinear by
  # R's partial matching ("rd" -> "rd_path"), which held only while no other
  # formal shared the prefix -- px_path/keep_policies made that a coin flip
  # waiting to be tossed.
  as_paths<-function(p)setNames(p,paste0(names(p),"_path"))
  ## A transition-probability input is perturbed as a pi_input_path, not as a
  ## price path: it moves Pi_s rather than the budget. Everything else about
  ## the sweep -- the central difference, the horizon, the outputs -- is
  ## identical, which is the point of routing both through one nonlinear
  ## transition.
  ## `threads` is FORWARDED, not dropped. It used to be absent from this
  ## function's formals entirely, so every ND battery silently resolved its own
  ## thread count from getOption("dynhr.hank3_threads") / a machine default --
  ## while its sibling hank_het3_jacobian() took an explicit `threads` and
  ## hank_td3_nonlinear() (called right here) has always accepted one. A
  ## downstream paper's multi-hour ND gate was running at the machine default
  ## for exactly this reason, and under _R_CHECK_LIMIT_CORES_ that default
  ## clamps to 2. Output is bit-identical at every thread count, so this is a
  ## wall-clock and reproducibility-of-EFFORT fix, not a numerical one.
  run<-function(p,pi_paths)do.call(hank_td3_nonlinear,
    c(list(block=block,T_h=T_h,pi_input_paths=pi_paths,threads=threads),as_paths(p)))
  for(i in inputs){
    if(i=="px"&&.hank3_px_is_inert(block))next   # exact zero column; see above
    is_pi<-.hank3_is_pi_input(block,i)
    for(s in seq_len(T_h)){
      p<-m<-base;pip<-pim<-NULL
      if(is_pi){x0<-rep(block$Pi_inputs[[i]],T_h);xp<-x0;xm<-x0;xp[s]<-xp[s]+delta;xm[s]<-xm[s]-delta
        pip<-setNames(list(xp),i);pim<-setNames(list(xm),i)
      } else {p[[i]][s]<-p[[i]][s]+delta;m[[i]][s]<-m[[i]][s]-delta}
      op<-run(p,pip);om<-run(m,pim)
      for(o in outputs)J[[o]][[i]][,s]<-(op[[o]]-om[[o]])/(2*delta)}
  }
  J
}

.hank3_arr_to_vec <- function(x) as.vector(aperm(x, c(4L, 3L, 2L, 1L)))

## Column-block width for the fake-news assembly (see hank_het3_jacobian).
## Chosen from a memory budget rather than fixed: the buffer is n_cell x B
## doubles, so B scales down as the grid grows. 512 MB by default, overridable
## via getOption("dynhr.hank3_fn_block_mb"); capped at 256 because the measured
## return flattens well before that, and floored at 1 so a huge grid still
## works (degenerating to the old column-at-a-time behaviour).
.hank3_fn_block <- function(n_cell,
                            budget_mb = getOption("dynhr.hank3_fn_block_mb", 512)) {
  if (!is.numeric(budget_mb) || length(budget_mb) != 1L || !is.finite(budget_mb) ||
      budget_mb <= 0)
    stop("hank_het3_jacobian: dynhr.hank3_fn_block_mb must be a positive number of MB.")
  max(1L, min(256L, as.integer(budget_mb * 1024^2 / (8 * max(1, n_cell)))))
}

.hank_curly_sweep3 <- function(block, T_h, input, outputs,
                               delta_in, delta_v, delta_d, threads = NULL,
                               backend = getOption("dynhr.hank3_backend", "cpp")) {
  D_ss <- block$D
  n_cell <- length(D_ss)
  state_dim <- dim(block$d)
  n_asset <- prod(state_dim[-1L])
  e_policy <- array(rep(block$e, times = n_asset), state_dim)
  omega_policy <- array(rep(.hank_block_omega(block), times = n_asset),
                        state_dim)

  is_pi <- .hank3_is_pi_input(block, input)

  aggregate <- function(x) sum(D_ss * .hank3_arr_to_vec(x))
  distribution_response <- function(step_p, step_m, scale,
                                     Pi_p = block$Pi, Pi_m = block$Pi) {
    .hank_forward_legs3(block, step_p, step_m, scale, Pi_p, Pi_m)
  }
  one_step <- function(Vd, Vf, Va, prices, Pi = block$Pi) {
    if (identical(backend, "cpp")) {
      .hank_egm3_step(
        block$d_grid, block$f_grid, block$a_grid,
        prices$w * block$e + prices$Tr * .hank_block_omega(block), Pi,
        prices$rd, prices$rf, prices$ra,
        block$beta, block$eis,
        block$chi0, block$chi1, block$chi2,
        block$phi0, block$phi1, block$phi2,
        Vd, Vf, Va, prices$px, threads
      )
    } else {
      ## backend = "R": .hank_egm3_step() unconditionally routes to the
      ## COMPILED kernel (hank_egm3_step_cpp) whenever the block is not the
      ## singleton-f reduction, ignoring dynhr.hank3_backend entirely -- so
      ## the s >= 2 per-step loop below (the "documented pure-R diagnostic
      ## path") silently ran the compiled step even under
      ## options(dynhr.hank3_backend = "R") (#10, adversarial review). Call
      ## .hank_egm3_step_r() directly instead, which threads `backend`
      ## through to hank_egm3_solve()'s own R implementation (capped at 250
      ## states, same as every other R-backend three-asset path).
      .hank_egm3_step_r(
        block$d_grid, block$f_grid, block$a_grid,
        prices$w * block$e + prices$Tr * .hank_block_omega(block), Pi,
        prices$rd, prices$rf, prices$ra,
        block$beta, block$eis,
        block$chi0, block$chi1, block$chi2,
        block$phi0, block$phi1, block$phi2,
        Vd, Vf, Va, prices$px, backend = backend
      )
    }
  }
  record <- function(step_p, step_m, scale, s, curlyY, curlyD,
                     Pi_p = block$Pi, Pi_m = block$Pi) {
    dd <- (step_p$d - step_m$d) / (2 * scale)
    df <- (step_p$f - step_m$f) / (2 * scale)
    da <- (step_p$a - step_m$a) / (2 * scale)
    dc <- (step_p$c - step_m$c) / (2 * scale)
    if ("D" %in% outputs) curlyY$D[s] <- aggregate(dd)
    if ("F" %in% outputs) curlyY$F[s] <- aggregate(df)
    if ("A" %in% outputs) curlyY$A[s] <- aggregate(da)
    if ("C" %in% outputs) curlyY$C[s] <- aggregate(dc)
    # The adjustment costs are differentiated as the kernel REPORTS them, not
    # re-derived from da/df here: Psi/Phi depend on the current asset state as
    # well as the policy, so rebuilding them outside the kernel would silently
    # drop the Psi2 term and produce a Jacobian that only ND could catch.
    if ("CHI" %in% outputs)
      curlyY$CHI[s] <- aggregate((step_p$chi - step_m$chi) / (2 * scale))
    if ("PHI" %in% outputs)
      curlyY$PHI[s] <- aggregate((step_p$phi - step_m$phi) / (2 * scale))
    ## Internal labour/transfer income output used to differentiate the exact
    ## household budget. Current income moves directly only for a current wage
    ## or transfer shock; all later income effects come from curly-D acting on
    ## the steady income state vector in the fake-news assembly.
    if ("Y" %in% outputs)
      curlyY$Y[s] <- if (s == 1L && input == "w") aggregate(e_policy) else
        if (s == 1L && input == "Tr") aggregate(omega_policy) else 0
    curlyD[, s] <- distribution_response(step_p, step_m, scale, Pi_p, Pi_m)
    list(
      curlyY = curlyY,
      curlyD = curlyD,
      dVd = (step_p$Vd - step_m$Vd) / (2 * scale),
      dVf = (step_p$Vf - step_m$Vf) / (2 * scale),
      dVa = (step_p$Va - step_m$Va) / (2 * scale)
    )
  }

  curlyY <- setNames(lapply(outputs, function(x) numeric(T_h)), outputs)
  curlyD <- matrix(0, n_cell, T_h)
  prices <- list(rd = block$rd, rf = block$rf,
                 ra = block$ra, w = block$w,
                 px = if (is.null(block$px)) 1 else block$px,
                 Tr = .hank_block_tr(block))
  pp <- pm <- prices
  ## s = 1: the direct input shock at the current date.
  ##
  ## A PRICE moves the budget, and the shock date's forward push is the
  ## policy-only derivative. A TRANSITION-PROBABILITY input instead moves
  ## Pi_1, which enters in TWO places (the file header of R/hank-jacobian.R
  ## sets out the same decomposition for the one-asset block):
  ##   (a) the date-1 backward step, whose expectation uses Pi_1 -- so the
  ##       policies react at all t <= 1 exactly as they do for a price, and
  ##       the anticipation for s >= 2 rides the same (dVd, dVf, dVa)
  ##       recursion with NO further change;
  ##   (b) Lambda_1 DIRECTLY, so curly-D at the shock date is the joint
  ##       (policy, Pi) directional derivative, taken with one shared step in
  ##       the input so both legs describe the same perturbation.
  ## Omitting (b) is a silent, plausible-looking error: the column keeps the
  ## right sign and shape and is simply wrong in magnitude, which is why the
  ## ND oracle is the acceptance gate for these columns and not a formality.
  Pi_p <- Pi_m <- block$Pi
  if (is_pi) {
    Pi_in_p <- .hank_pi_perturb(block, input, +delta_in)
    Pi_in_m <- .hank_pi_perturb(block, input, -delta_in)
    ## Use the SAME actual Pi legs as the policy solve. Reconstructing a second
    ## delta_d pair would no longer describe the nonlinear legs passed to the
    ## accounting-preserving forward derivative.
    Pi_p <- Pi_in_p
    Pi_m <- Pi_in_m
  } else {
    pp[[input]] <- pp[[input]] + delta_in
    pm[[input]] <- pm[[input]] - delta_in
  }
  first <- record(
    one_step(block$Vd, block$Vf, block$Va, pp,
             if (is_pi) Pi_in_p else block$Pi),
    one_step(block$Vd, block$Vf, block$Va, pm,
             if (is_pi) Pi_in_m else block$Pi),
    delta_in, 1L, curlyY, curlyD, Pi_p, Pi_m
  )
  curlyY <- first$curlyY
  curlyD <- first$curlyD
  dVd <- first$dVd
  dVf <- first$dVf
  dVa <- first$dVa

  ## FUSED COMPILED SWEEP. hank_curly_sweep3_cpp() runs the whole s = 2..T_h
  ## recursion -- both EGM step legs per date, record()'s differencing and
  ## aggregation, and the curly-D leg accumulation -- in one compiled call
  ## under ONE worker pool, instead of ~13 interpreted n_cell passes plus two
  ## pool spawns per date. It reuses the step and forward-legs kernel bodies
  ## (egm3_expect/egm3_run_tasks/forward_legs3_core), so the arithmetic is the
  ## same operations in the same order and the contract is BIT-IDENTITY with
  ## the loop below, which test-hank-jacobian3-fused.R asserts with
  ## expect_identical().
  ##
  ## GATE. The singleton-f reduction routes .hank_egm3_step to the two-asset
  ## reduction solver (see R/hank-egm3.R), which the fused kernel does not
  ## implement -- it must FALL THROUGH to the loop below, like the two-asset
  ## sweep's collateral gate. The s = 1 term above stays in R for every path:
  ## that is where the input-type (price vs Pi) dispatch lives.
  ## options(dynhr.hank3_fused_sweep = FALSE) is the escape hatch that forces
  ## the per-step loop back on with everything else unchanged -- the
  ## equivalence tests use it, and it is the first thing to try if a
  ## three-asset Jacobian ever looks wrong.
  ## `identical(backend, "cpp")` mirrors the two-asset gate
  ## (R/hank-jacobian2.R's `use_fused <- identical(backend, "cpp") && ...`):
  ## the fused kernel IS the compiled backend, so it must not run under
  ## options(dynhr.hank3_backend = "R") (#10, adversarial review) -- that
  ## falls through to the per-step loop below, which now itself honours
  ## `backend` via `one_step()`.
  use_fused <- identical(backend, "cpp") && T_h >= 2L &&
    length(block$f_grid) > 1L &&
    isTRUE(getOption("dynhr.hank3_fused_sweep", TRUE))
  if (use_fused) {
    fuse_outputs <- intersect(outputs, c("D", "F", "A", "C", "CHI", "PHI"))
    ## The internal "Y" output is identically zero at s >= 2 (record() only
    ## sets it at the shock date for w/Tr), so curlyY$Y keeps its zero
    ## initialization; the compiled sweep carries the variable output set.
    fs <- hank_curly_sweep3_cpp(
      block$Vd, block$Vf, block$Va, dVd, dVf, dVa,
      block$d_grid, block$f_grid, block$a_grid,
      prices$w * block$e + prices$Tr * .hank_block_omega(block), block$Pi,
      prices$rd, prices$rf, prices$ra, block$beta, block$eis,
      block$chi0, block$chi1, block$chi2,
      block$phi0, block$phi1, block$phi2, prices$px,
      block$D, fuse_outputs, delta_v, T_h,
      hank_resolve_threads(threads))
    tail <- seq_len(T_h - 1L) + 1L
    for (o in fuse_outputs) curlyY[[o]][tail] <- fs$curlyY[[o]]
    curlyD[, tail] <- fs$curlyD
    return(list(curlyY = curlyY, curlyD = curlyD))
  }

  if (T_h >= 2L) for (s in 2L:T_h) {
    # All three marginal values describe one perturbation direction.  A single
    # shared scale is essential: scaling them separately changes that direction.
    h <- delta_v / max(1, max(abs(dVd)), max(abs(dVf)), max(abs(dVa)))
    next_step <- record(
      one_step(block$Vd + h * dVd, block$Vf + h * dVf,
               block$Va + h * dVa, prices),
      one_step(block$Vd - h * dVd, block$Vf - h * dVf,
               block$Va - h * dVa, prices),
      h, s, curlyY, curlyD
    )
    curlyY <- next_step$curlyY
    curlyD <- next_step$curlyD
    dVd <- next_step$dVd
    dVf <- next_step$dVf
    dVa <- next_step$dVa
  }
  list(curlyY = curlyY, curlyD = curlyD)
}

#' Fake-news Jacobian for the three-asset reference block
#'
#' Sequence-space Jacobian of the three-asset block via the fake-news
#' algorithm (Auclert, Bardoczy, Rognlie & Straub 2021), the three-asset
#' counterpart of \code{\link{hank_het_jacobian}} /
#' \code{\link{hank_het2_jacobian}}. Because the two costly assets are chosen
#' JOINTLY, a single perturbation direction must move all three marginal
#' values \code{(Vd, Vf, Va)} together -- \code{.hank_curly_sweep3}
#' rescales them by one shared step size for exactly this reason, since
#' scaling them separately would change the direction being differentiated --
#' and the forward-operator derivative (\code{.hank_forward_direction3})
#' perturbs all three policies at once. The local policy derivative follows
#' the steady-state active set (which asset, if any, sits at a corner);
#' \code{\link{hank_het3_jacobian_nd}} is the numerical acceptance oracle at
#' any kink where that could go wrong.
#'
#' @param block A \code{\link{hank_het3_block}}.
#' @param T_h Integer horizon.
#' @param inputs Character subset of \code{c("rd", "rf", "ra", "w", "px", "Tr")}
#'   plus, for a block built with \code{Pi_fn}/\code{Pi_inputs}, that block's
#'   named transition-probability inputs (\code{names(block$Pi_inputs)}, e.g.
#'   the job-finding rate \code{"f"} and separation rate \code{"s"} of
#'   \code{\link{hank_employment_income}}). \code{NULL} (the default) takes
#'   all of them. A transition-probability column is NOT a relabelled price
#'   column: the input moves the income transition matrix, so it enters the
#'   distribution's law of motion directly as well as through the policy --
#'   see the \code{s = 1} branch of \code{.hank_curly_sweep3}.
#'   \code{"px"} is the foreign-valuation column the Stage-4 contract's
#'   P-Jacobian gate requires. It is a genuine extra column rather than a
#'   rescaling of \code{"rf"}, because \code{px} enters the PURCHASE side of
#'   the budget as well as the payoff side (see \code{\link{hank_egm3_solve}});
#'   had it only revalued the predetermined stock it would be exactly
#'   \code{(1 + rf)} times the \code{rf} column, and matching ND on it would
#'   prove nothing.
#' @param outputs Character subset of \code{c("D", "F", "A", "C", "CHI",
#'   "PHI")} (aggregate liquid, foreign, capital holdings, consumption, and
#'   capital/foreign adjustment resources).
#' @param threads Worker threads for the backward policy step, which is
#'   \strong{87.8\%} of this function's cost (measured by \code{Rprof} at the
#'   32,256-state rung: 25.44 s of 28.98 s -- the step dominates the sweep, and
#'   the sweep dominates the Jacobian). \code{NULL} (the default) resolves via
#'   \code{getOption("dynhr.hank3_threads")} and then a machine-derived
#'   default, exactly as \code{\link{hank_egm3_solve}} does; \code{1} forces
#'   the serial path. The compiled step is BIT-IDENTICAL at every thread count,
#'   so this is a throughput knob and never a source of numerical difference.
#' @param delta_in Numeric FD step for the input (rd/rf/ra/w/px) perturbation
#'   in the \code{s = 1} (contemporaneous) term.
#' @param delta_v Numeric FD step for the backward joint-marginal-value
#'   propagation used in the \code{s >= 2} anticipation terms; rescaled
#'   internally by the current magnitude of \code{(dVd, dVf, dVa)} so the step
#'   stays a small RELATIVE perturbation of the propagated direction.
#' @param delta_d Retained for API/cache compatibility. The distributional
#'   (curly-D) response now differences the actual plus/minus policy legs at
#'   their generating \code{delta_in} or propagated-value step, so it does not
#'   construct a second \code{delta_d} perturbation.
#' @param backend Character: \code{"cpp"} (default, from
#'   \code{getOption("dynhr.hank3_backend")}) or \code{"R"} for the pure-R
#'   diagnostic path (capped at 250 states per \code{\link{hank_egm3_solve}}).
#'   Mirrors \code{\link{hank_het2_jacobian}}'s \code{backend} argument: under
#'   \code{"R"} the fused compiled sweep is skipped (it IS the compiled
#'   backend), and the per-step \code{s >= 2} loop also routes every backward
#'   step through the R implementation, not just the \code{s = 1} term.
#'
#' @return Nested list \code{J[[output]][[input]]}, each a \code{T_h x T_h}
#'   matrix with \code{[t, s] = dO_t/dI_s}. (The date-0 response vectors are
#'   the first Jacobian row, \code{J[[o]][[i]][1, ]}.)
#' @seealso \code{\link{hank_het3_jacobian_nd}} (the numerical oracle this is
#'   validated against), \code{\link{hank_het_jacobian}} (one-asset),
#'   \code{\link{hank_het2_jacobian}} (two-asset), \code{\link{hank_het3_block}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' J <- hank_het3_jacobian(blk, T_h = 2)
#' dim(J$D$rd)
#' @export
hank_het3_jacobian <- function(block, T_h,
                               inputs = NULL,
                               outputs = .hank3_jac_outputs,
                               delta_in = 1e-5, delta_v = 1e-6,
                               delta_d = 1e-6, threads = NULL,
                               backend = getOption("dynhr.hank3_backend", "cpp")) {
  backend <- match.arg(backend, c("cpp", "R"))
  if (!inherits(block, "hank_het3_block"))
    stop("hank_het3_jacobian: block must be hank_het3_block")
  if (!is.numeric(T_h) || length(T_h) != 1L || T_h < 1 || !is.finite(T_h))
    stop("hank_het3_jacobian: T_h must be positive")
  if (is.null(inputs)) inputs <- .hank3_jac_inputs(block)
  inputs <- .hank3_check_inputs(block, inputs, "hank_het3_jacobian")
  bad <- setdiff(outputs, .hank3_jac_outputs)
  if (length(bad)) stop("unsupported output(s): ", paste(bad, collapse = ", "))

  output_policy <- list(D = block$d, F = block$f,
                        A = block$a, C = block$c,
                        CHI = block$chi, PHI = block$phi)
  # Blocks saved by a dynhr that predates the adjustment-resource outputs carry
  # no chi/phi, and CHI/PHI are now in the DEFAULT outputs -- so say what is
  # wrong instead of failing on a NULL policy deep inside the sweep. Unlike px
  # (which defaults to 1, the pre-A4 kernel exactly) there is no safe default
  # here: silently recomputing the costs would hide the version mismatch.
  stale <- intersect(outputs, names(Filter(is.null, output_policy)))
  if (length(stale))
    stop("hank_het3_jacobian: this block carries no ",
         paste(tolower(stale), collapse = "/"),
         " policy, so output(s) ", paste(stale, collapse = ", "),
         " cannot be formed. It was built by a dynhr predating the ",
         "adjustment-resource outputs; rebuild it with hank_het3_block().")

  ## Consumption is pinned by the kernel's exact statewise budget. Computing
  ## its finite-step policy derivative separately from the asset/cost and
  ## distribution derivatives leaves small active-set/product-rule residuals
  ## (especially px*F). Carry the steady income state as an internal output and
  ## differentiate that exact identity after assembling the other outputs.
  requested_outputs <- outputs
  work_outputs <- outputs
  budget_close <- "C" %in% outputs &&
    !any(vapply(output_policy[c("D", "F", "A", "CHI", "PHI")],
                is.null, logical(1L)))
  ## C is DROPPED from the swept set when the budget closes it, not merely
  ## added-to. Sweeping it would cost a full extra output -- its curly-Y
  ## aggregate per input plus, in the assembly below, an E_s = Lambda^s y
  ## streaming pass over every cell -- and then be overwritten by the identity.
  ## That is ~1/7 of the output work on a seven-output build, and the streaming
  ## pass is the part that scales with the grid, so on a million-state block it
  ## is the whole saving. J$C is created by the identity block itself.
  if (budget_close)
    work_outputs <- setdiff(unique(c(outputs, "D", "F", "A", "CHI", "PHI", "Y")),
                            "C")
  state_dim <- dim(block$d)
  n_asset <- prod(state_dim[-1L])
  if (budget_close)
    output_policy$Y <- array(rep(block$y, times = n_asset), state_dim)

  J <- setNames(lapply(work_outputs, function(output)
    setNames(lapply(inputs, function(input) matrix(0, T_h, T_h)), inputs)),
    work_outputs)

  # Exact-zero columns on the reduction are never swept; see .hank3_px_is_inert.
  valid_inputs <- inputs[!(inputs == "px" & .hank3_px_is_inert(block))]

  ## E-HOISTING. The expectation stream E_s = Lambda^s y depends only on
  ## (block, output policy) -- NOT on the input being differentiated -- so
  ## the streaming loop below was regenerating the identical E_s sequence
  ## once per input, n_inputs times over. When every input's curlyD (n_cell
  ## x T_h) fits in memory SIMULTANEOUSLY, run every input's backward sweep
  ## first, then loop outputs on the OUTSIDE and generate each E block ONCE,
  ## crossprod-ing it against every input's curlyD in turn. Above the memory
  ## budget (`getOption("dynhr.hank3_jac_e_hoist_bytes", 2e9)`), fall back to
  ## the original per-input streaming path unchanged -- that is the only path
  ## that works at the paper's 1.29M-state scale, where holding every input's
  ## curlyD at once is not an option.
  ##
  ## Bit-identity is the contract: the hoisted path performs the SAME
  ## per-(input, output) crossprod(curlyD, E_block) calls, in the same block
  ## widths and the same sequential E_s generation order, as the streaming
  ## path -- only the loop nesting changes, never any accumulation's
  ## association. E_s is a pure function of (block, output), so which input's
  ## loop it happens to be generated under cannot change its bits.
  n_cell <- length(block$D)
  e_hoist_bytes <- getOption("dynhr.hank3_jac_e_hoist_bytes", 2e9)
  use_hoist <- length(valid_inputs) > 0L &&
    as.numeric(length(valid_inputs)) * n_cell * T_h * 8 <= e_hoist_bytes

  .assemble_jacobian <- function(fake_news) {
    jacobian <- matrix(0, T_h, T_h)
    jacobian[1L, ] <- fake_news[1L, ]
    if (T_h >= 2L) for (t in 2L:T_h) {
      jacobian[t, 1L] <- fake_news[t, 1L]
      jacobian[t, 2L:T_h] <-
        jacobian[t - 1L, 1L:(T_h - 1L)] + fake_news[t, 2L:T_h]
    }
    jacobian
  }

  if (use_hoist) {
    sweeps <- setNames(lapply(valid_inputs, function(input)
      .hank_curly_sweep3(block, T_h, input, work_outputs,
                         delta_in, delta_v, delta_d, threads, backend)),
      valid_inputs)
    for (output in work_outputs) {
      fake_news <- setNames(
        lapply(valid_inputs, function(i) matrix(0, T_h, T_h)), valid_inputs)
      for (input in valid_inputs)
        fake_news[[input]][1L, ] <- sweeps[[input]]$curlyY[[output]]
      expectation <- .hank3_arr_to_vec(output_policy[[output]])
      if (T_h >= 2L) {
        blk_w <- .hank3_fn_block(n_cell)
        t <- 2L
        while (t <= T_h) {
          nb <- min(blk_w, T_h - t + 1L)
          E <- matrix(0, n_cell, nb)
          for (j in seq_len(nb)) {
            E[, j] <- expectation
            expectation <- .hank_forward_apply3(block, expectation)
          }
          for (input in valid_inputs)
            fake_news[[input]][t:(t + nb - 1L), ] <-
              t(crossprod(sweeps[[input]]$curlyD, E))
          t <- t + nb
        }
      }
      for (input in valid_inputs)
        J[[output]][[input]] <- .assemble_jacobian(fake_news[[input]])
    }
  } else {
    for (input in valid_inputs) {
      sweep <- .hank_curly_sweep3(
        block, T_h, input, work_outputs, delta_in, delta_v, delta_d, threads,
        backend
      )
      for (output in work_outputs) {
        fake_news <- matrix(0, T_h, T_h)
        fake_news[1L, ] <- sweep$curlyY[[output]]
        # Stream E_s = Lambda^s y rather than retaining T full state vectors for
        # every output, but multiply in COLUMN BLOCKS.
        #
        # One row at a time, `crossprod(curlyD, expectation)` is a dgemv that
        # streams the whole n_cell x T curlyD for EVERY period, so the assembly
        # moves n_cell*T^2*8 bytes and is purely bandwidth-bound. Buffering B
        # expectation vectors and issuing one crossprod(curlyD, E_block) is the
        # SAME arithmetic, but reads curlyD once per block -- traffic falls by a
        # factor of B, and a dgemm is compute-bound and BLAS-threaded where a
        # dgemv is not. Measured 16x at B=32 and 27x at B=128 on the assembly
        # itself (n_cell = 322,560, T = 400).
        #
        # This is the O(T^2) term, and it is invisible until T is large: 0.14% of
        # a T=20 Jacobian, but 34% of a projected T=400 one -- a share that ROSE
        # when the policy step was threaded, because that cut the linear term ~4x.
        #
        # The expectations are still generated strictly sequentially, so nothing
        # about the recursion changes; only B of them are held at once. Results
        # agree with the column-at-a-time form to BLAS accumulation order
        # (~1e-13 relative), not bitwise -- an unavoidable consequence of letting
        # the BLAS choose its summation order, and far inside every gate here.
        expectation <- .hank3_arr_to_vec(output_policy[[output]])
        if (T_h >= 2L) {
          blk_w <- .hank3_fn_block(n_cell)
          t <- 2L
          while (t <= T_h) {
            nb <- min(blk_w, T_h - t + 1L)
            E <- matrix(0, n_cell, nb)
            for (j in seq_len(nb)) {
              E[, j] <- expectation
              expectation <- .hank_forward_apply3(block, expectation)
            }
            ## crossprod gives T_h x nb; column j is the fake-news row for
            ## period t + j - 1.
            fake_news[t:(t + nb - 1L), ] <- t(crossprod(sweep$curlyD, E))
            t <- t + nb
          }
        }
        J[[output]][[input]] <- .assemble_jacobian(fake_news)
      }
    }
  }
  if (budget_close) {
    I0 <- diag(T_h)
    L <- if (T_h == 1L) matrix(0, 1L, 1L) else
      rbind(0, cbind(diag(T_h - 1L), 0))
    px <- if (is.null(block$px)) 1 else block$px
    ## Income depends on the INCOME STATE ALONE -- hank_het3_block() builds
    ## y = w*e + Tr*omega, both functions of e -- so with Pi fixed an asset
    ## lottery cannot move the income-state marginal, and every fixed-Pi price
    ## input has an EXACTLY zero income-distribution response. Enforcing that
    ## analytically avoids retaining an O(N)/step accumulation residue from
    ## curly-D, which on million-state blocks is the dominant error term.
    ##
    ## The precondition is that y has that form, NOT that the transfer is zero.
    ## An earlier version tested `identical(y, w*e)`, which is the Tr == 0 case
    ## only: on a transfer block it fell through to the accumulated path and
    ## the budget residual degraded from ~2e-16 to ~2e-11 on a 54-cell block --
    ## silently, and growing with the grid. Compared with a tolerance rather
    ## than `identical()` because a block carrying a y built in a different
    ## association order is still structurally the same income function.
    y_struct <- block$w * block$e +
      .hank_block_tr(block) * .hank_block_omega(block)
    if (length(block$y) == length(y_struct) &&
        isTRUE(all.equal(as.numeric(block$y), as.numeric(y_struct),
                         tolerance = 1e-12))) {
      pi_inputs <- if (is.null(block$Pi_inputs)) character(0) else
        names(block$Pi_inputs)
      for (input in setdiff(inputs, c(pi_inputs, "w", "Tr")))
        J$Y[[input]][] <- 0
      if ("w" %in% inputs) {
        e_state <- .hank3_arr_to_vec(array(rep(block$e, times = n_asset),
                                           state_dim))
        J$Y$w <- sum(block$D * e_state) * I0
      }
      if ("Tr" %in% inputs) {
        omega_state <- .hank3_arr_to_vec(array(
          rep(.hank_block_omega(block), times = n_asset), state_dim))
        J$Y$Tr <- sum(block$D * omega_state) * I0
      }
    }
    for (input in inputs) {
      C <- J$Y[[input]] +
        (1 + block$rd) * L %*% J$D[[input]] +
        px * (1 + block$rf) * L %*% J$F[[input]] +
        (1 + block$ra) * L %*% J$A[[input]] -
        J$D[[input]] - px * J$F[[input]] - J$A[[input]] -
        J$CHI[[input]] - J$PHI[[input]]
      if (input == "rd") C <- C + block$D_agg * I0
      if (input == "rf") C <- C + px * block$F_agg * I0
      if (input == "ra") C <- C + block$A_agg * I0
      if (input == "px") C <- C + block$rf * block$F_agg * I0
      J$C[[input]] <- C
    }
  }
  J[requested_outputs]
}
