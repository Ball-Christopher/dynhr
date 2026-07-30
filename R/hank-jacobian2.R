## R/hank-jacobian2.R
## --------------------------------------------------------------------------
## Sequence-space Jacobian of the TWO-ASSET household block via the fake-news
## algorithm (Auclert, Bardoczy, Rognlie & Straub 2021). The two-asset
## counterpart of R/hank-jacobian.R.
##
## The algorithm is unchanged from the one-asset case -- backward sweep to
## curly-Y / curly-D, expectation vectors E_s = Lambda_ss^s y^o, fake-news
## matrix, diagonal cumulation -- because steps 2-4 are steady-state objects
## that only care about the cell space, not about how many assets built it.
## What changes is step 1, in three specific ways:
##
##   (a) The household carries TWO marginal values (Vb, Va), so the
##       anticipation recursion propagates the PAIR. They are scaled by ONE
##       shared step h: they are two derivatives of a single household
##       problem, and rescaling them independently would step along a
##       different direction than (dVb, dVa) and give the wrong directional
##       derivative.
##   (b) curly-D must perturb BOTH policies together, since the joint forward
##       operator is a function of (b', a') jointly -- perturbing one at a time
##       and adding would drop the interaction through the product lottery.
##   (c) There are three price inputs (rb, ra, w) rather than two, and `ra`
##       additionally enters the ADJUSTMENT COST -- so a date-s perturbation of
##       ra moves Psi and Psi1 as well as the budget. .hank_block_step2()
##       rebuilds Psi1 from the perturbed ra for exactly this reason; a cached
##       Psi1 would silently make dJ/dra wrong while leaving every other column
##       right, which the ND oracle is what catches.
##
## This file does NOT touch .hank_curly_sweep(): the one-asset sweep is
## validated against its own brute-force oracle and is load-bearing for the
## whole one-asset HANK stack. The ~40 lines of shared cumulation are
## duplicated deliberately (extracting a common helper is an optional later
## cleanup, not part of this wave).
##
## The mandatory oracle is hank_het2_jacobian_nd(): brute-force numerical
## differentiation of hank_td2_nonlinear(). See test-hank-jacobian2.R.
## --------------------------------------------------------------------------


#' Validate requested two-asset het-block Jacobian inputs
#'
#' The admissible aggregate inputs are the prices \code{c("rb", "ra", "w")}
#' plus, for a block built with \code{Pi_fn}/\code{Pi_inputs}, the block's named
#' transition-probability inputs.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param inputs Character vector of requested inputs.
#' @return \code{inputs}, validated.
#' @keywords internal
.hank_het2_check_inputs <- function(block, inputs) {
  if (!inherits(block, "hank_het2_block"))
    stop("hank_het2_jacobian: 'block' must be a hank_het2_block.")
  allowed <- c("rb", "ra", "w", "Tr", "theta_coll", names(block$Pi_inputs))
  bad <- setdiff(inputs, allowed)
  if (length(bad))
    stop("unsupported two-asset het-block input(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; this block supports ",
         paste0("'", allowed, "'", collapse = ", "),
         if (is.null(block$Pi_inputs))
           " (build the block with Pi_fn/Pi_inputs to add transition-probability inputs)"
         else "", ".")
  inputs
}


#' Validate requested two-asset het-block outputs
#' @keywords internal
.hank_het2_check_outputs <- function(outputs) {
  allowed <- c("B", "A", "C", "CHI")
  bad <- setdiff(outputs, allowed)
  if (length(bad))
    stop("unsupported two-asset het-block output(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; supported outputs are ",
         paste0("'", allowed, "'", collapse = ", "), ".")
  outputs
}


#' Brute-force numerical-differentiation Jacobian of a two-asset het block
#'
#' Reference sequence-space Jacobian: perturbs each aggregate input at each
#' date and re-runs the full nonlinear transition
#' (\code{\link{hank_td2_nonlinear}}), central-differencing the aggregate
#' outputs.  Correct but \eqn{O(T)} solves per (input, date) -- this is the
#' mandatory oracle for \code{\link{hank_het2_jacobian}}, not a production
#' path.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param T_h Integer horizon.
#' @param inputs Character subset of \code{c("rb", "ra", "w", "Tr")} plus the block's
#'   transition-probability inputs.
#' @param outputs Character subset of \code{c("B", "A", "C", "CHI")}.
#' @param delta Numeric FD step for the input perturbation.
#' @param backend Character: \code{"cpp"} (default, from
#'   \code{getOption("dynhr.hank_backend")}) or \code{"R"} for the pure-R
#'   reference backward step.  Forwarded to \code{\link{hank_td2_nonlinear}}.
#' @param threads Worker threads for the compiled backward step, or \code{NULL}
#'   (default) to resolve (see \code{\link{hank_resolve_threads}}).  Forwarded
#'   to every \code{\link{hank_td2_nonlinear}} call; the kernel is
#'   bit-identical across thread counts, so this is a speed knob only.
#'
#' @return Nested list \code{J[[output]][[input]]}, each a \code{T_h x T_h}
#'   matrix with \code{[t, s] = dO_t/dI_s}.
#' @seealso \code{\link{hank_het2_jacobian}},
#'   \code{\link{hank_het_jacobian_nd}} (one-asset)
#' @export
hank_het2_jacobian_nd <- function(block, T_h,
                                  inputs = c("rb", "ra", "w"),
                                  outputs = c("B", "A", "C"),
                                  delta = 1e-5,
                                  backend = getOption("dynhr.hank_backend",
                                                      "cpp"),
                                  threads = NULL) {
  backend <- match.arg(backend, c("R", "cpp"))
  threads <- hank_resolve_threads(threads)
  inputs  <- .hank_het2_check_inputs(block, inputs)
  outputs <- .hank_het2_check_outputs(outputs)
  base <- list(rb = rep(block$rb, T_h), ra = rep(block$ra, T_h),
               w = rep(block$w, T_h),
               Tr = rep(.hank_block_tr(block), T_h),
               theta_coll = rep(if (is.null(block$theta_coll)) 0
                                else block$theta_coll, T_h))

  J <- setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) matrix(0, T_h, T_h)), inputs)), outputs)

  for (i in inputs) {
    for (s in seq_len(T_h)) {
      pp <- base; pm <- base
      pip_p <- NULL; pip_m <- NULL
      if (i %in% c("rb", "ra", "w", "Tr", "theta_coll")) {
        pp[[i]][s] <- pp[[i]][s] + delta
        pm[[i]][s] <- pm[[i]][s] - delta
      } else {
        x0 <- rep(block$Pi_inputs[[i]], T_h)
        xp <- x0; xp[s] <- xp[s] + delta
        xm <- x0; xm[s] <- xm[s] - delta
        pip_p <- setNames(list(xp), i)
        pip_m <- setNames(list(xm), i)
      }
      out_p <- hank_td2_nonlinear(block, rb_path = pp$rb, ra_path = pp$ra,
                                  w_path = pp$w, T_h = T_h,
                                  pi_input_paths = pip_p, Tr_path = pp$Tr,
                                  theta_path = pp$theta_coll,
                                  backend = backend, threads = threads)
      out_m <- hank_td2_nonlinear(block, rb_path = pm$rb, ra_path = pm$ra,
                                  w_path = pm$w, T_h = T_h,
                                  pi_input_paths = pip_m, Tr_path = pm$Tr,
                                  theta_path = pm$theta_coll,
                                  backend = backend, threads = threads)
      for (o in outputs)
        J[[o]][[i]][, s] <- (out_p[[o]] - out_m[[o]]) / (2 * delta)
    }
  }
  J
}


#' Backward sweep for the two-asset fake-news algorithm
#'
#' Two-asset counterpart of \code{\link{.hank_curly_sweep}}: curly-Y (per-output
#' date-0 outcome response) and curly-D (induced one-period-ahead distribution
#' change) to an anticipated shock to input \code{i} at horizon
#' \code{s = 1 .. T_h}.
#'
#' See the file header for the three ways this differs from the one-asset
#' sweep: the \code{(dVb, dVa)} pair propagates under ONE shared step; curly-D
#' perturbs BOTH policies jointly; and \code{ra} moves the adjustment cost as
#' well as the budget.
#'
#' The \code{s >= 2} loop is the expensive part (two backward steps per date).
#' Under \code{backend = "cpp"}, and only for a block and an input that carry no
#' collateral, it is delegated whole to \code{hank_curly_sweep2_cpp()}, which
#' holds ONE worker pool across all \code{2*(T_h-1)} steps instead of spawning
#' one per step. Everything else -- the \code{s = 1} term, curly-D and the
#' curlyY aggregation -- stays in R.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param T_h Integer horizon.
#' @param i Character name of the aggregate input being shocked (one of
#'   \code{"rb"}, \code{"ra"}, \code{"w"}, \code{"Tr"}, \code{"theta_coll"}, or
#'   one of the block's transition-probability inputs).
#' @param outputs Character subset of \code{c("B", "A", "C", "CHI")}.
#' @param delta_in FD step for the input perturbation in the \code{s = 1} term.
#' @param delta_va Relative FD step for the backward value-function propagation
#'   (applied to the \code{(Vb, Va)} pair jointly).
#' @param delta_d FD step for the distributional (curly-D) response.
#' @param backend Character \code{"cpp"}/\code{"R"} backward-step backend.
#' @param threads Resolved integer worker count for the compiled backend.
#' @return List with \code{curlyY} (named list over \code{outputs}, each
#'   length-\code{T_h}) and \code{curlyD} (\code{n_cell x T_h}).
#' @keywords internal
.hank_curly_sweep2 <- function(block, T_h, i, outputs,
                               delta_in, delta_va, delta_d,
                               backend = getOption("dynhr.hank_backend", "cpp"),
                               threads = NULL) {
  b_grid <- block$b_grid; a_grid <- block$a_grid
  Pi <- block$Pi; D_ss <- block$D
  Vb_ss <- block$Vb; Va_ss <- block$Va
  b_ss <- block$b; a_ss <- block$a
  is_pi_input <- !(i %in% c("rb", "ra", "w", "Tr", "theta_coll"))
  n_cell <- block$n_e * block$n_b * block$n_a
  ## COLLATERAL (D1, brief 19 section 9.12): policies live in the gap
  ## coordinate x = b + theta*a, so the LIQUID-output response is
  ## d(b_liq) = dB - theta*dA for EVERY input; the theta input additionally
  ## enters TWO adjacent steps (theta_today at the shock date, dtheta_next one
  ## period before -- the coordinate re-basing), contributes -a_ss to the
  ## date-s liquid output (B_s = E[x' - theta_s a']), and shifts the forward
  ## operator's liquid argument by -+a_ss at distances 1/2.  All corrections
  ## vanish identically at theta = 0.
  th   <- if (is.null(block$theta_coll)) 0 else block$theta_coll
  aggB <- function(dB, dA) hank_aggregate2(D_ss, if (th > 0) dB - th * dA
                                                 else dB)

  ## Distributional response (curly-D) to a JOINT policy change (dB, dA).
  ## Both policies move together: the joint operator is a function of (b', a')
  ## jointly, so perturbing them one at a time would drop the interaction
  ## through the product lottery.
  ##
  ## MATRIX-FREE (0.9.0.0026): the two sparse hank_forward_operator2() builds
  ## this used to do per call were 94-97% of the (build + matvec) cost across
  ## the whole size range measured, so the sweep spent nearly all its forward
  ## time assembling n_cell x n_cell operators for a product that never needs
  ## one. .hank_forward_push2() contracts the same identity (see its comment
  ## in R/hank-distribution2.R); agreement with the sparse route is
  ## round-off-level, not approximate. Note the 1/(2*delta_d) ~ 1e6
  ## amplification below: that is why the parity gate on the push is 1e-14
  ## RELATIVE and not merely "close".
  curlyD_from_pol <- function(dB, dA, Pi_p = Pi, Pi_m = Pi) {
    (.hank_forward_push2(b_ss + delta_d * dB, a_ss + delta_d * dA,
                         b_grid, a_grid, Pi_p, D_ss) -
       .hank_forward_push2(b_ss - delta_d * dB, a_ss - delta_d * dA,
                           b_grid, a_grid, Pi_m, D_ss)) / (2 * delta_d)
  }
  agg <- function(x) hank_aggregate2(D_ss, x)

  curlyY <- setNames(lapply(outputs, function(o) numeric(T_h)), outputs)
  curlyD <- matrix(0, n_cell, T_h)

  ## --- s = 1: the direct input shock at the current date --------------------
  px <- list(rb = block$rb, ra = block$ra, w = block$w,
             Tr = .hank_block_tr(block),
             theta_coll = if (is.null(block$theta_coll)) 0
                          else block$theta_coll)
  if (is_pi_input) {
    sp <- .hank_block_step2(block, Vb_ss, Va_ss, px$rb, px$ra, px$w,
                            Pi = .hank_pi_perturb(block, i, +delta_in),
                            backend = backend, threads = threads)
    sm <- .hank_block_step2(block, Vb_ss, Va_ss, px$rb, px$ra, px$w,
                            Pi = .hank_pi_perturb(block, i, -delta_in),
                            backend = backend, threads = threads)
  } else {
    pp <- px; pp[[i]] <- pp[[i]] + delta_in
    pm <- px; pm[[i]] <- pm[[i]] - delta_in
    ## NOTE for i == "ra": .hank_block_step2 rebuilds Psi1 from the ra it is
    ## handed, so the adjustment-cost channel is differenced here too.
    ## NOTE for i == "theta_coll": a one-period deviation at the shock date
    ## also re-bases tomorrow's coordinate, so the SAME step carries
    ## dtheta_next = theta_ss - theta_shock = -(the perturbation).
    sp <- .hank_block_step2(block, Vb_ss, Va_ss, pp$rb, pp$ra, pp$w,
                            Tr = pp$Tr, theta_coll = pp$theta_coll,
                            dtheta_next = if (i == "theta_coll") -delta_in
                                          else 0,
                            backend = backend, threads = threads)
    sm <- .hank_block_step2(block, Vb_ss, Va_ss, pm$rb, pm$ra, pm$w,
                            Tr = pm$Tr, theta_coll = pm$theta_coll,
                            dtheta_next = if (i == "theta_coll") +delta_in
                                          else 0,
                            backend = backend, threads = threads)
  }
  dB   <- (sp$b   - sm$b)   / (2 * delta_in)
  dA   <- (sp$a   - sm$a)   / (2 * delta_in)
  dC   <- (sp$c   - sm$c)   / (2 * delta_in)
  dCHI <- (sp$chi - sm$chi) / (2 * delta_in)
  dVb  <- (sp$Vb  - sm$Vb)  / (2 * delta_in)
  dVa  <- (sp$Va  - sm$Va)  / (2 * delta_in)
  if ("B" %in% outputs)
    curlyY[["B"]][1L] <- aggB(dB, dA) -
      (if (i == "theta_coll") hank_aggregate2(D_ss, a_ss) else 0)
  if ("A"   %in% outputs) curlyY[["A"]][1L]   <- agg(dA)
  if ("C"   %in% outputs) curlyY[["C"]][1L]   <- agg(dC)
  if ("CHI" %in% outputs) curlyY[["CHI"]][1L] <- agg(dCHI)
  if (is_pi_input) {
    ## Pi enters Lambda DIRECTLY on the shock date, so difference the joint
    ## (policy, Pi) direction with the same step.
    curlyD[, 1L] <- curlyD_from_pol(dB, dA,
                                    Pi_p = .hank_pi_perturb(block, i, +delta_d),
                                    Pi_m = .hank_pi_perturb(block, i, -delta_d))
  } else if (i == "theta_coll") {
    ## Mass chosen at x' on the shock date arrives at x' + (dtheta_next)a'
    ## with dtheta_next = -dtheta: an EXTRA liquid-policy shift of -a_ss per
    ## unit theta in the operator.
    curlyD[, 1L] <- curlyD_from_pol(dB - a_ss, dA)
  } else {
    curlyD[, 1L] <- curlyD_from_pol(dB, dA)
  }

  ## theta only: DATE-1 COORDINATE RE-BASING of the predetermined states.  An
  ## unanticipated theta_1 move re-expresses the physical (b, a) states as
  ## x = b + theta_1*a, i.e. dD_1 = d/dtheta[shift of D_ss by dtheta*a]
  ## (mirroring hank_td2_nonlinear).  It adds (i) dD_1-weighted steady-state
  ## outputs to every curlyY at s = 1 and (ii) its propagation through the
  ## steady-state operator to curlyD[, 1].
  if (i == "theta_coll") {
    dD1 <- (.hank2_shift_x(D_ss, +delta_d, block) -
              .hank2_shift_x(D_ss, -delta_d, block)) / (2 * delta_d)
    y1 <- list(B   = if (is.null(block$b_liq)) b_ss else block$b_liq,
               A   = a_ss, C = block$c, CHI = block$chi)
    for (o in outputs)
      curlyY[[o]][1L] <- curlyY[[o]][1L] +
        sum(dD1 * .hank2_arr_to_vec(y1[[o]]))
    curlyD[, 1L] <- curlyD[, 1L] +
      as.numeric(Matrix::t(block$Lambda) %*% dD1)
  }

  ## --- s >= 2: propagate the anticipation via the (Vb, Va) PAIR -------------
  ## seq_len(T_h - 1L) + 1L, not 2L:T_h -- the latter counts DOWN to c(2, 1)
  ## at T_h = 1 and indexes out of bounds.
  dVb_prev <- dVb; dVa_prev <- dVa

  ## FUSED COMPILED SWEEP (0.9.0.0038). The per-step compiled path spawns and
  ## joins a worker pool inside EVERY hank_egm2_step_cpp call -- 2*(T_h-1) = 98
  ## spawns at T_h = 50 around a sub-millisecond kernel, which is why the
  ## compiled port bought 3.2-3.5x from compilation alone and NOTHING from
  ## threading. hank_curly_sweep2_cpp() runs this whole loop under ONE pool,
  ## exactly as hank_egm2_solve_cpp() already does for the solve.
  ##
  ## GATE. The compiled kernel has no collateral coordinate (there is no
  ## theta_coll anywhere in src/hank_egm2.cpp), so a collateral block is
  ## excluded -- and so is the "theta_coll" INPUT even at theta_coll = 0,
  ## because its s = 2 step carries dtheta_next = +-h. Both fall through to the
  ## per-step loop below, which dispatches to the R reference step exactly as
  ## before. backend = "R" also falls through, keeping the pure-R sweep
  ## reachable as the reference.
  ## options(dynhr.hank_fused_sweep = FALSE) forces the per-step loop back on
  ## with everything else unchanged. That is the ESCAPE HATCH the equivalence
  ## tests use to difference the two compiled routes against each other, and
  ## the first thing to try if a two-asset Jacobian ever looks wrong.
  use_fused <- identical(backend, "cpp") && th == 0 &&
    i != "theta_coll" && T_h >= 2L &&
    isTRUE(getOption("dynhr.hank_fused_sweep", TRUE))
  if (use_fused) {
    n_a <- block$n_a
    Psi1 <- .hank_psi(matrix(block$a_grid, n_a, n_a),
                      matrix(block$a_grid, n_a, n_a, byrow = TRUE),
                      block$ra, block$chi0, block$chi1, block$chi2)$Psi1
    fs <- hank_curly_sweep2_cpp(Vb_ss, Va_ss, dVb, dVa,
                                block$b_grid, block$a_grid, block$k_grid,
                                px$w * block$e + px$Tr,
                                px$rb, px$ra, block$beta, block$eis,
                                block$chi0, block$chi1, block$chi2,
                                as.matrix(Pi), Psi1,
                                delta_va, T_h,
                                hank_resolve_threads(threads))
    dm <- c(block$n_e, block$n_b, block$n_a)
    as_arr <- function(v) array(v, dm)
    for (s in seq_len(T_h - 1L) + 1L) {
      j <- s - 1L
      dB <- as_arr(fs$dB[, j]); dA <- as_arr(fs$dA[, j])
      if ("B"   %in% outputs) curlyY[["B"]][s]   <- aggB(dB, dA)
      if ("A"   %in% outputs) curlyY[["A"]][s]   <- agg(dA)
      if ("C"   %in% outputs) curlyY[["C"]][s]   <- agg(as_arr(fs$dC[, j]))
      if ("CHI" %in% outputs) curlyY[["CHI"]][s] <- agg(as_arr(fs$dCHI[, j]))
      curlyD[, s] <- curlyD_from_pol(dB, dA)
    }
    return(list(curlyY = curlyY, curlyD = curlyD))
  }

  for (s in seq_len(T_h - 1L) + 1L) {
    ## ONE shared step for both marginals: they are two derivatives of a single
    ## household problem, so the directional derivative must move along
    ## (dVb, dVa) jointly.  For the theta input at DISTANCE 2 only, the step
    ## one period before the shock ALSO sees dtheta_next = +dtheta (the
    ## arrival re-basing into the shock date's coordinate); it must ride the
    ## same FD direction as the (dVb, dVa) pair, scaled by the same h.
    h <- delta_va / max(1, max(abs(dVb_prev)), max(abs(dVa_prev)))
    dth2 <- if (i == "theta_coll" && s == 2L) h else 0
    sp <- .hank_block_step2(block, Vb_ss + h * dVb_prev, Va_ss + h * dVa_prev,
                            px$rb, px$ra, px$w, Tr = px$Tr,
                            theta_coll = px$theta_coll, dtheta_next = +dth2,
                            backend = backend, threads = threads)
    sm <- .hank_block_step2(block, Vb_ss - h * dVb_prev, Va_ss - h * dVa_prev,
                            px$rb, px$ra, px$w, Tr = px$Tr,
                            theta_coll = px$theta_coll, dtheta_next = -dth2,
                            backend = backend, threads = threads)
    dB   <- (sp$b   - sm$b)   / (2 * h)
    dA   <- (sp$a   - sm$a)   / (2 * h)
    dC   <- (sp$c   - sm$c)   / (2 * h)
    dCHI <- (sp$chi - sm$chi) / (2 * h)
    dVb  <- (sp$Vb  - sm$Vb)  / (2 * h)
    dVa  <- (sp$Va  - sm$Va)  / (2 * h)
    if ("B"   %in% outputs) curlyY[["B"]][s]   <- aggB(dB, dA)
    if ("A"   %in% outputs) curlyY[["A"]][s]   <- agg(dA)
    if ("C"   %in% outputs) curlyY[["C"]][s]   <- agg(dC)
    if ("CHI" %in% outputs) curlyY[["CHI"]][s] <- agg(dCHI)
    ## The same arrival re-basing shifts the forward operator's liquid
    ## argument by +a_ss at distance 2 (mass leaving the pre-shock date). But
    ## curlyY[1] (hence J[1,1]) already carries +sum(dD1 . y) from the date-1
    ## rebasing above, and by translation invariance date-2's own-shock
    ## response is IDENTICAL to date-1's -- so that shared dD1 piece must
    ## cancel out of F[2,2] = J[2,2] - J[1,1], or it double-counts (measured:
    ## it inflated every diagonal J[t,t], t >= 2, by the same constant
    ## sum(dD1 . a_ss)). Subtracting dD1 here is that cancellation.
    curlyD[, s] <- if (i == "theta_coll" && s == 2L)
      curlyD_from_pol(dB + a_ss, dA) - dD1 else curlyD_from_pol(dB, dA)
    dVb_prev <- dVb; dVa_prev <- dVa
  }

  list(curlyY = curlyY, curlyD = curlyD)
}


#' Sequence-space Jacobian of a two-asset het block via the fake-news algorithm
#'
#' The two-asset counterpart of \code{\link{hank_het_jacobian}}: the block
#' Jacobian \eqn{J^{o,i}[t,s] = dO_t/dI_s} of aggregate outputs
#' (\code{B}, \code{A}, \code{C}, \code{CHI}) with respect to anticipated
#' aggregate-input paths (\code{rb}, \code{ra}, \code{w}, plus any
#' transition-probability inputs), computed by one backward sweep plus cheap
#' bookkeeping instead of \eqn{O(T)} nonlinear solves.
#'
#' Validated against \code{\link{hank_het2_jacobian_nd}} (brute-force numerical
#' differentiation of the nonlinear transition) -- see the PRIMARY GATE in
#' \code{test-hank-jacobian2.R}.
#'
#' @inheritParams hank_het2_jacobian_nd
#' @param delta_in FD step for the input perturbation in the \code{s = 1} term.
#' @param delta_va Relative FD step for the backward value-function propagation
#'   (applied to the \code{(Vb, Va)} pair jointly).
#' @param delta_d FD step for the distributional (curly-D) response.
#' @param backend Character: \code{"cpp"} (default, from
#'   \code{getOption("dynhr.hank_backend")}) or \code{"R"}.  The backward sweep
#'   is ~76-80\% of a two-asset Jacobian build, so this is the knob that
#'   matters; \code{"R"} keeps the pure-R reference step reachable for
#'   debugging and cross-checking.  A collateral block (\code{theta_coll > 0})
#'   forces \code{"R"} internally.
#' @param threads Worker threads for the compiled backward step, or \code{NULL}
#'   (default) to resolve (see \code{\link{hank_resolve_threads}}); \code{1}
#'   forces the serial kernel.  Bit-identical across thread counts.  Under
#'   \code{backend = "cpp"} on a collateral-free block the whole
#'   \code{s >= 2} backward sweep runs in ONE compiled call holding ONE worker
#'   pool, so these threads are spawned once per input rather than once per
#'   step; \code{options(dynhr.hank_fused_sweep = FALSE)} reverts to the
#'   per-step calls (bit-identical, only slower).  Note the backward sweep is
#'   no longer where a two-asset Jacobian spends its time -- once compiled it
#'   is ~15\% of a build, and \code{.hank_forward_push2} (curly-D, still in R)
#'   is ~76\% -- so this knob moves the total by ~5-10\%, not by its own
#'   speedup.
#'
#' @return Nested list \code{J[[output]][[input]]}, each a \code{T_h x T_h}
#'   matrix with \code{[t, s] = dO_t/dI_s}.
#' @seealso \code{\link{hank_het_jacobian}} (one-asset),
#'   \code{\link{hank_het2_block}}, \code{\link{hank_het2_jacobian_nd}}
#' @export
hank_het2_jacobian <- function(block, T_h,
                               inputs = c("rb", "ra", "w"),
                               outputs = c("B", "A", "C"),
                               delta_in = 1e-5, delta_va = 1e-6,
                               delta_d = 1e-6,
                               backend = getOption("dynhr.hank_backend", "cpp"),
                               threads = NULL) {
  backend <- match.arg(backend, c("R", "cpp"))
  threads <- hank_resolve_threads(threads)
  inputs  <- .hank_het2_check_inputs(block, inputs)
  outputs <- .hank_het2_check_outputs(outputs)
  Lam <- block$Lambda

  ## --- Step 2: expectation vectors E_s = Lambda^s y^o, s = 0 .. T-1 ---------
  ## Under collateral the B output is the TRUE liquid policy b_liq = x' -
  ## theta*a' (identical to block$b at theta = 0; pre-D1 blocks lack the
  ## field).
  y_out <- list(B   = .hank2_arr_to_vec(if (is.null(block$b_liq)) block$b
                                        else block$b_liq),
                A   = .hank2_arr_to_vec(block$a),
                C   = .hank2_arr_to_vec(block$c),
                CHI = .hank2_arr_to_vec(block$chi))
  Elist <- setNames(vector("list", length(outputs)), outputs)
  for (o in outputs) {
    E <- vector("list", T_h)
    E[[1L]] <- y_out[[o]]                        # E_0
    if (T_h >= 2L)
      for (s in 2L:T_h) E[[s]] <- as.numeric(Lam %*% E[[s - 1L]])
    Elist[[o]] <- E
  }

  J <- setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) matrix(0, T_h, T_h)), inputs)), outputs)

  for (i in inputs) {
    ## --- Step 1: backward sweep -> curlyY[[o]][s], curlyD[, s] --------------
    sweep  <- .hank_curly_sweep2(block, T_h, i, outputs,
                                 delta_in, delta_va, delta_d,
                                 backend = backend, threads = threads)
    curlyY <- sweep$curlyY
    curlyD <- sweep$curlyD

    ## --- Steps 3-4: fake-news matrix F then Jacobian J ----------------------
    ## Identical to the one-asset recursion (R/hank-jacobian.R): these are
    ## steady-state objects and do not care how many assets built the cell
    ## space.
    for (o in outputs) {
      Fm <- matrix(0, T_h, T_h)
      Fm[1L, ] <- curlyY[[o]]                       # F[0, s] = curlyY[s]
      E <- Elist[[o]]
      if (T_h >= 2L)
        for (tt in 2L:T_h)
          Fm[tt, ] <- as.numeric(crossprod(curlyD, E[[tt - 1L]]))
      Jm <- matrix(0, T_h, T_h)
      Jm[1L, ] <- Fm[1L, ]
      if (T_h >= 2L)
        for (tt in 2L:T_h) {
          Jm[tt, 1L] <- Fm[tt, 1L]
          Jm[tt, 2L:T_h] <- Jm[tt - 1L, 1L:(T_h - 1L)] + Fm[tt, 2L:T_h]
        }
      J[[o]][[i]] <- Jm
    }
  }
  J
}
