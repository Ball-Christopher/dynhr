## R/hank-jacobian2d.R
## --------------------------------------------------------------------------
## Fake-news sequence-space Jacobian for the DISCRETE-ADJUSTMENT two-asset
## household (hank_het2d_block). Same algorithm as the smooth two-asset case
## (R/hank-jacobian2.R) with two structural extensions:
##
##   (a) the anticipation recursion propagates the (dV, dVb, dVa) TRIPLE --
##       the value LEVEL is a state of the recursion here because the discrete
##       choice compares branch levels -- under ONE shared FD step (they are
##       three derivatives of a single household problem);
##   (b) curly-D perturbs the WHOLE transition quintuple (b_A, a_A, b_N, a_N,
##       P) jointly through the P-mixed operator: the choice probability is
##       part of the law of motion, so a price shock moves the distribution
##       through WHO adjusts as well as through where adjusters go. P is
##       deliberately NOT clamped to [0, 1] inside the central difference --
##       clamping would bias the derivative, and the mixture with P slightly
##       outside the unit interval is still a well-defined linear operator for
##       FD purposes.
##
## The taste-shock scale sigma_taste is what makes this Jacobian well-posed:
## it smooths the adjust/no-adjust boundary in prices, which is exactly why
## the deterministic (sigma -> 0) limit is not offered.
##
## Inputs: the aggregate prices (rb, ra, w), the transfer Tr, and -- for a
## block built with Pi_fn/Pi_inputs -- its named transition-probability
## inputs, which now reach this tier too (they were skipped in the original
## wave for want of a use case, not for any structural reason: the discrete
## adjust/no-adjust choice lives inside the backward step, and Pi enters that
## step and the forward operator exactly as it does at every other tier).
## Validated against hank_het2d_jacobian_nd -- brute-force differentiation of
## hank_td2d_nonlinear -- at the package's standard 1e-5 gate.
## --------------------------------------------------------------------------


#' Validate requested discrete-adjustment het-block Jacobian inputs
#' @keywords internal
.hank_het2d_check_inputs <- function(block, inputs) {
  if (!inherits(block, "hank_het2d_block"))
    stop("hank_het2d_jacobian: 'block' must be a hank_het2d_block (the ",
         "smooth two-asset block routes through hank_het2_jacobian).")
  allowed <- c("rb", "ra", "w", "Tr", names(block$Pi_inputs))
  bad <- setdiff(inputs, allowed)
  if (length(bad))
    stop("unsupported het2d-block input(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; supported inputs are ",
         paste0("'", allowed, "'", collapse = ", "),
         if (is.null(block$Pi_inputs))
           " (build the block with Pi_fn/Pi_inputs to add transition-probability inputs)"
         else "", ".")
  inputs
}


#' Validate requested discrete-adjustment het-block outputs
#' @keywords internal
.hank_het2d_check_outputs <- function(outputs) {
  allowed <- c("B", "A", "C", "CHI", "ADJ")
  bad <- setdiff(outputs, allowed)
  if (length(bad))
    stop("unsupported het2d-block output(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; supported outputs are ",
         paste0("'", allowed, "'", collapse = ", "), ".")
  outputs
}


#' Per-cell outcome vectors for one backward step of the het2d household
#' @keywords internal
.hank_het2d_outcomes <- function(st, F_adj) {
  list(B   = st$P * st$b_A + (1 - st$P) * st$b_N,
       A   = st$P * st$a_A + (1 - st$P) * st$a_N,
       C   = st$P * st$c_A + (1 - st$P) * st$c_N,
       CHI = st$P * (st$chi_A + F_adj),
       ADJ = st$P)
}


#' Backward sweep for the discrete-adjustment fake-news algorithm
#' @keywords internal
.hank_curly_sweep2d <- function(block, T_h, i, outputs,
                                delta_in, delta_va, delta_d) {
  D_ss <- block$D
  n_cell <- block$n_e * block$n_b * block$n_a
  F_adj <- block$F_adj

  ## Directional derivative of the forward update through the FULL quintuple.
  is_pi_input <- !(i %in% c("rb", "ra", "w", "Tr"))

  curlyD_from_step <- function(dst, Pi_p = block$Pi, Pi_m = block$Pi) {
    per <- function(sgn) list(
      b_A = block$b_A + sgn * delta_d * dst$b_A,
      a_A = block$a_A + sgn * delta_d * dst$a_A,
      b_N = block$b_N + sgn * delta_d * dst$b_N,
      a_N = block$a_N + sgn * delta_d * dst$a_N,
      P   = block$P   + sgn * delta_d * dst$P)
    ## MATRIX-FREE (0.9.0.0026): each hank_forward_operator2d() build is TWO
    ## sparse hank_forward_operator2() builds, so this call site was the
    ## heaviest instance of the "materialize an n_cell x n_cell operator for a
    ## single matvec" pattern in the package -- four builds per evaluation.
    ## .hank_forward_push2d() reweights the source mass by P instead (see its
    ## comment in R/hank-distribution2.R).
    (.hank_forward_push2d(per(+1), block$b_grid, block$a_grid, Pi_p, D_ss) -
       .hank_forward_push2d(per(-1), block$b_grid, block$a_grid, Pi_m,
                            D_ss)) / (2 * delta_d)
  }
  agg <- function(x) hank_aggregate2(D_ss, x)
  out_ss <- .hank_het2d_outcomes(block, F_adj)

  diff_step <- function(sp, sm, h2) {
    op <- .hank_het2d_outcomes(sp, F_adj); om <- .hank_het2d_outcomes(sm, F_adj)
    list(dY = lapply(setNames(outputs, outputs),
                     function(o) (op[[o]] - om[[o]]) / h2),
         dq = list(b_A = (sp$b_A - sm$b_A) / h2, a_A = (sp$a_A - sm$a_A) / h2,
                   b_N = (sp$b_N - sm$b_N) / h2, a_N = (sp$a_N - sm$a_N) / h2,
                   P   = (sp$P   - sm$P)   / h2),
         dV  = (sp$V  - sm$V)  / h2,
         dVb = (sp$Vb - sm$Vb) / h2,
         dVa = (sp$Va - sm$Va) / h2)
  }

  curlyY <- setNames(lapply(outputs, function(o) numeric(T_h)), outputs)
  curlyD <- matrix(0, n_cell, T_h)

  ## s = 1: the direct input shock this period.
  px <- list(rb = block$rb, ra = block$ra, w = block$w,
             Tr = .hank_block_tr(block))
  pp <- px; pm <- px
  Pi_in_p <- Pi_in_m <- block$Pi
  if (is_pi_input) {
    ## Pi_1 enters the date-1 backward step AND Lambda_1 directly, exactly as
    ## in the smooth two-asset and three-asset sweeps. The discrete
    ## adjust/no-adjust choice changes nothing about that decomposition: it
    ## lives inside the step, and the quintuple derivative already carries it.
    Pi_in_p <- .hank_pi_perturb(block, i, +delta_in)
    Pi_in_m <- .hank_pi_perturb(block, i, -delta_in)
  } else {
    pp[[i]] <- pp[[i]] + delta_in
    pm[[i]] <- pm[[i]] - delta_in
  }
  sp <- .hank_block_step2d(block, block$V, block$Vb, block$Va,
                           pp$rb, pp$ra, pp$w, Tr = pp$Tr, Pi = Pi_in_p)
  sm <- .hank_block_step2d(block, block$V, block$Vb, block$Va,
                           pm$rb, pm$ra, pm$w, Tr = pm$Tr, Pi = Pi_in_m)
  d1 <- diff_step(sp, sm, 2 * delta_in)
  for (o in outputs) curlyY[[o]][1L] <- agg(d1$dY[[o]])
  curlyD[, 1L] <- if (is_pi_input)
    curlyD_from_step(d1$dq,
                     Pi_p = .hank_pi_perturb(block, i, +delta_d),
                     Pi_m = .hank_pi_perturb(block, i, -delta_d))
  else curlyD_from_step(d1$dq)

  ## s >= 2: anticipation, propagated via the (dV, dVb, dVa) TRIPLE under one
  ## shared step -- three derivatives of a single household problem.
  dV_prev <- d1$dV; dVb_prev <- d1$dVb; dVa_prev <- d1$dVa
  for (s in seq_len(T_h - 1L) + 1L) {
    h <- delta_va / max(1, max(abs(dV_prev)), max(abs(dVb_prev)),
                        max(abs(dVa_prev)))
    sp <- .hank_block_step2d(block, block$V + h * dV_prev,
                             block$Vb + h * dVb_prev,
                             block$Va + h * dVa_prev,
                             px$rb, px$ra, px$w, Tr = px$Tr)
    sm <- .hank_block_step2d(block, block$V - h * dV_prev,
                             block$Vb - h * dVb_prev,
                             block$Va - h * dVa_prev,
                             px$rb, px$ra, px$w, Tr = px$Tr)
    ds <- diff_step(sp, sm, 2 * h)
    for (o in outputs) curlyY[[o]][s] <- agg(ds$dY[[o]])
    curlyD[, s] <- curlyD_from_step(ds$dq)
    dV_prev <- ds$dV; dVb_prev <- ds$dVb; dVa_prev <- ds$dVa
  }
  list(curlyY = curlyY, curlyD = curlyD)
}


#' Brute-force ND Jacobian of a discrete-adjustment het block
#'
#' Reference Jacobian by central differences of
#' \code{\link{hank_td2d_nonlinear}} -- \eqn{O(T)} full transitions per
#' (input, date), the mandatory oracle for \code{\link{hank_het2d_jacobian}}.
#'
#' @param block A \code{\link{hank_het2d_block}}.
#' @param T_h Horizon.
#' @param inputs Subset of \code{c("rb", "ra", "w", "Tr")}.
#' @param outputs Subset of \code{c("B", "A", "C", "CHI", "ADJ")}.
#' @param delta FD step.
#' @return Nested list \code{J[[output]][[input]]} of \code{T_h x T_h}
#'   matrices.
#' @seealso \code{\link{hank_het2d_jacobian}}
#' @keywords internal
#' @export
hank_het2d_jacobian_nd <- function(block, T_h,
                                   inputs = c("rb", "ra", "w"),
                                   outputs = c("B", "A", "C"),
                                   delta = 1e-5) {
  inputs  <- .hank_het2d_check_inputs(block, inputs)
  outputs <- .hank_het2d_check_outputs(outputs)
  base <- list(rb = rep(block$rb, T_h), ra = rep(block$ra, T_h),
               w = rep(block$w, T_h),
               Tr = rep(.hank_block_tr(block), T_h))
  J <- setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) matrix(0, T_h, T_h)), inputs)), outputs)
  for (i in inputs) for (s in seq_len(T_h)) {
    pp <- base; pm <- base; pip <- NULL; pim <- NULL
    if (i %in% c("rb", "ra", "w", "Tr")) {
      pp[[i]][s] <- pp[[i]][s] + delta
      pm[[i]][s] <- pm[[i]][s] - delta
    } else {
      ## Transition-probability input: swept as a pi_input_path, since it
      ## moves Pi_s rather than the budget.
      x0 <- rep(block$Pi_inputs[[i]], T_h)
      xp <- x0; xp[s] <- xp[s] + delta
      xm <- x0; xm[s] <- xm[s] - delta
      pip <- setNames(list(xp), i); pim <- setNames(list(xm), i)
    }
    out_p <- hank_td2d_nonlinear(block, rb_path = pp$rb, ra_path = pp$ra,
                                 w_path = pp$w, Tr_path = pp$Tr, T_h = T_h,
                                 pi_input_paths = pip)
    out_m <- hank_td2d_nonlinear(block, rb_path = pm$rb, ra_path = pm$ra,
                                 w_path = pm$w, Tr_path = pm$Tr, T_h = T_h,
                                 pi_input_paths = pim)
    for (o in outputs)
      J[[o]][[i]][, s] <- (out_p[[o]] - out_m[[o]]) / (2 * delta)
  }
  J
}


#' Fake-news Jacobian of a discrete-adjustment two-asset het block
#'
#' See the file header for the two structural extensions over the smooth
#' two-asset algorithm (the \code{(dV, dVb, dVa)} triple; the joint quintuple
#' perturbation in curly-D).  Steps 2-4 (expectation vectors, fake-news
#' cumulation) are the standard steady-state recursion and do not care about
#' the discrete choice.
#'
#' @inheritParams hank_het2d_jacobian_nd
#' @param delta_in,delta_va,delta_d FD steps (input, value-triple,
#'   distribution).
#' @return Nested list \code{J[[output]][[input]]} of \code{T_h x T_h}
#'   matrices with \code{[t, s] = dO_t/dI_s}.
#' @seealso \code{\link{hank_het2d_jacobian_nd}} (the oracle),
#'   \code{\link{hank_het2_jacobian}} (smooth)
#' @export
hank_het2d_jacobian <- function(block, T_h,
                                inputs = c("rb", "ra", "w"),
                                outputs = c("B", "A", "C"),
                                delta_in = 1e-5, delta_va = 1e-6,
                                delta_d = 1e-6) {
  inputs  <- .hank_het2d_check_inputs(block, inputs)
  outputs <- .hank_het2d_check_outputs(outputs)
  Lam <- block$Lambda

  y_out <- lapply(.hank_het2d_outcomes(block, block$F_adj), .hank2_arr_to_vec)
  Elist <- setNames(vector("list", length(outputs)), outputs)
  for (o in outputs) {
    E <- vector("list", T_h)
    E[[1L]] <- y_out[[o]]
    for (s in seq_len(T_h - 1L) + 1L)
      E[[s]] <- as.numeric(Lam %*% E[[s - 1L]])
    Elist[[o]] <- E
  }

  J <- setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) matrix(0, T_h, T_h)), inputs)), outputs)
  for (i in inputs) {
    sweep  <- .hank_curly_sweep2d(block, T_h, i, outputs,
                                  delta_in, delta_va, delta_d)
    for (o in outputs) {
      Fm <- matrix(0, T_h, T_h)
      Fm[1L, ] <- sweep$curlyY[[o]]
      E <- Elist[[o]]
      for (tt in seq_len(T_h - 1L) + 1L)
        Fm[tt, ] <- as.numeric(crossprod(sweep$curlyD, E[[tt - 1L]]))
      Jm <- matrix(0, T_h, T_h)
      Jm[1L, ] <- Fm[1L, ]
      for (tt in seq_len(T_h - 1L) + 1L) {
        Jm[tt, 1L] <- Fm[tt, 1L]
        Jm[tt, 2L:T_h] <- Jm[tt - 1L, 1L:(T_h - 1L)] + Fm[tt, 2L:T_h]
      }
      J[[o]][[i]] <- Jm
    }
  }
  J
}
