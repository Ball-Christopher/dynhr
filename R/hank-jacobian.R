## R/hank-jacobian.R
## --------------------------------------------------------------------------
## Sequence-space Jacobian of a heterogeneous-agent household block via the
## fake-news algorithm (Auclert, Bardoczy, Rognlie & Straub 2021, Econometrica).
##
## The block Jacobian J^{o,i}[t, s] = dO_t / dI_s is the response of aggregate
## output o at date t to a (perfect-foresight, anticipated) shock to aggregate
## input i at date s.  Computed naively it needs O(T) full backward+forward
## solves; the fake-news algorithm reduces it to a single backward sweep plus
## cheap bookkeeping (ABRS Prop. 1 / eq. 25):
##
##   1. Backward sweep -> curly-Y (date-0 outcome response to a shock s periods
##      ahead) and curly-D (the induced one-period-ahead distribution change).
##   2. Expectation vectors E_s = Lambda_ss^s y^o  (Lemma 3).
##   3. Fake-news matrix  F[0,s] = curlyY[s];  F[t,s] = <E_{t-1}, curlyD[s]>.
##   4. Jacobian          J[0,s] = F[0,s];  J[t,s] = J[t-1,s-1] + F[t,s].
##
## Per-step derivatives here are taken by central finite differences (the ABRS
## reference uses analytic within-step derivatives; econpizza uses autodiff --
## all three validate against the same brute-force numerical-differentiation
## reference, which is the mandatory oracle in test-hank-jacobian.R).
##
## Lambda_ss is the row-stochastic forward operator (hank_forward_operator), so
## E_s = Lambda_ss %*% E_{s-1} with NO transpose (E_s(x) = expected future
## outcome from state x), while distributions push forward as t(Lambda) %*% D.
## --------------------------------------------------------------------------


#' Brute-force numerical-differentiation Jacobian of a het block
#'
#' Reference sequence-space Jacobian: perturbs each aggregate input at each date
#' and re-runs the full nonlinear transition (\code{\link{hank_td_nonlinear}}),
#' central-differencing the aggregate outputs.  Correct but \eqn{O(T)} solves
#' per (input, date); used to validate \code{\link{hank_het_jacobian}}.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param T_h Integer horizon.
#' @param inputs Character subset of \code{c("r", "w")}.
#' @param outputs Character subset of \code{c("A", "C")}.
#' @param delta Numeric FD step for the input perturbation.
#'
#' @return Nested list \code{J[[output]][[input]]}, each a \code{T x T} matrix
#'   with \code{[t, s] = dO_t/dI_s}.
#' @export
hank_het_jacobian_nd <- function(block, T_h,
                                 inputs = c("r", "w"),
                                 outputs = c("A", "C"),
                                 delta = 1e-5) {
  inputs  <- match.arg(inputs, several.ok = TRUE)
  outputs <- match.arg(outputs, several.ok = TRUE)
  r0 <- rep(block$r, T_h); w0 <- rep(block$w, T_h)

  J <- setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) matrix(0, T_h, T_h)), inputs)), outputs)

  for (i in inputs) {
    for (s in seq_len(T_h)) {
      rp <- r0; wp <- w0; rm <- r0; wm <- w0
      if (i == "r") { rp[s] <- rp[s] + delta; rm[s] <- rm[s] - delta }
      else          { wp[s] <- wp[s] + delta; wm[s] <- wm[s] - delta }
      out_p <- hank_td_nonlinear(block, r_path = rp, w_path = wp, T_h = T_h)
      out_m <- hank_td_nonlinear(block, r_path = rm, w_path = wm, T_h = T_h)
      for (o in outputs)
        J[[o]][[i]][, s] <- (out_p[[o]] - out_m[[o]]) / (2 * delta)
    }
  }
  J
}


#' Brute-force numerical-differentiation DISTRIBUTION Jacobian of a het block
#'
#' Reference sequence-space distribution Jacobian \eqn{J^D[t,s,] = dD_t/dI_s}
#' (the change in the full \code{(n_e*n_a)}-vector cross-sectional distribution
#' at date \code{t} induced by a shock to aggregate input \code{i} at date
#' \code{s}), computed exactly like \code{\link{hank_het_jacobian_nd}} but
#' keeping the full \code{Dpath} instead of aggregating it into \code{A}/\code{C}.
#' Used to validate \code{\link{hank_het_dist_jacobian}}.
#'
#' @inheritParams hank_het_jacobian_nd
#'
#' @return Named list \code{JD_nd[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_a)} with \code{JD_nd[[i]][t, s, ] =
#'   dD_t/dI_s} (central difference).
#' @export
hank_het_dist_jacobian_nd <- function(block, T_h,
                                      inputs = c("r", "w"),
                                      delta = 1e-5) {
  inputs <- match.arg(inputs, several.ok = TRUE)
  r0 <- rep(block$r, T_h); w0 <- rep(block$w, T_h)
  n_cell <- block$n_e * block$n_a

  JD_nd <- setNames(lapply(inputs, function(i) array(0, c(T_h, T_h, n_cell))),
                    inputs)

  for (i in inputs) {
    for (s in seq_len(T_h)) {
      rp <- r0; wp <- w0; rm <- r0; wm <- w0
      if (i == "r") { rp[s] <- rp[s] + delta; rm[s] <- rm[s] - delta }
      else          { wp[s] <- wp[s] + delta; wm[s] <- wm[s] - delta }
      out_p <- hank_td_nonlinear(block, r_path = rp, w_path = wp, T_h = T_h)
      out_m <- hank_td_nonlinear(block, r_path = rm, w_path = wm, T_h = T_h)
      dD <- (out_p$Dpath - out_m$Dpath) / (2 * delta)   # (n_cell x T_h), col t
      for (tt in seq_len(T_h)) JD_nd[[i]][tt, s, ] <- dD[, tt]
    }
  }
  JD_nd
}


#' Backward sweep shared by \code{hank_het_jacobian} and
#' \code{hank_het_dist_jacobian}: curly-Y (per-output date-0 outcome response)
#' and curly-D (induced one-period-ahead distribution change) to an anticipated
#' shock to input \code{i} at horizon \code{s = 1 .. T_h} (s=1 is the direct
#' current-period shock; s>=2 propagate via the value-function derivative).
#'
#' @return List with \code{curlyY} (named list over \code{outputs}, each a
#'   length-\code{T_h} vector) and \code{curlyD} (\code{(n_e*n_a) x T_h}
#'   matrix, column \code{s}).
#' @keywords internal
.hank_curly_sweep <- function(block, T_h, i, outputs,
                              delta_in, delta_va, delta_d) {
  a_grid <- block$a_grid; Pi <- block$Pi; D_ss <- block$D
  Va_ss <- block$Va; a_ss <- block$a

  ## Helper: distributional response (curly-D) to a savings-policy change dA.
  curlyD_from_dA <- function(dA) {
    Lp <- hank_forward_operator(a_ss + delta_d * dA, a_grid, Pi)
    Lm <- hank_forward_operator(a_ss - delta_d * dA, a_grid, Pi)
    (as.numeric(Matrix::t(Lp) %*% D_ss) -
       as.numeric(Matrix::t(Lm) %*% D_ss)) / (2 * delta_d)
  }

  curlyY <- setNames(lapply(outputs, function(o) numeric(T_h)), outputs)
  curlyD <- matrix(0, block$n_e * block$n_a, T_h)

  ## s = 1: direct input shock at the current date.
  if (i == "r") {
    sp <- .hank_block_step(block, Va_ss, block$r + delta_in, block$w)
    sm <- .hank_block_step(block, Va_ss, block$r - delta_in, block$w)
  } else {
    sp <- .hank_block_step(block, Va_ss, block$r, block$w + delta_in)
    sm <- .hank_block_step(block, Va_ss, block$r, block$w - delta_in)
  }
  dA  <- (sp$a  - sm$a)  / (2 * delta_in)
  dC  <- (sp$c  - sm$c)  / (2 * delta_in)
  dVa <- (sp$Va - sm$Va) / (2 * delta_in)
  if ("A" %in% outputs) curlyY[["A"]][1L] <- hank_aggregate(D_ss, dA)
  if ("C" %in% outputs) curlyY[["C"]][1L] <- hank_aggregate(D_ss, dC)
  curlyD[, 1L] <- curlyD_from_dA(dA)

  ## s >= 2: propagate the anticipation via the value-function derivative.
  dVa_prev <- dVa
  for (s in 2L:T_h) {
    h  <- delta_va / max(1, max(abs(dVa_prev)))
    sp <- .hank_block_step(block, Va_ss + h * dVa_prev, block$r, block$w)
    sm <- .hank_block_step(block, Va_ss - h * dVa_prev, block$r, block$w)
    dA  <- (sp$a  - sm$a)  / (2 * h)
    dC  <- (sp$c  - sm$c)  / (2 * h)
    dVa <- (sp$Va - sm$Va) / (2 * h)
    if ("A" %in% outputs) curlyY[["A"]][s] <- hank_aggregate(D_ss, dA)
    if ("C" %in% outputs) curlyY[["C"]][s] <- hank_aggregate(D_ss, dC)
    curlyD[, s] <- curlyD_from_dA(dA)
    dVa_prev <- dVa
  }

  list(curlyY = curlyY, curlyD = curlyD)
}


#' Sequence-space Jacobian of a het block via the fake-news algorithm
#'
#' @inheritParams hank_het_jacobian_nd
#' @param delta_in FD step for the input (r/w) perturbation in the s=0 term.
#' @param delta_va Relative FD step for the backward value-function propagation.
#' @param delta_d FD step for the distributional (curly-D) response.
#'
#' @return Nested list \code{J[[output]][[input]]}, each a \code{T x T} matrix
#'   with \code{[t, s] = dO_t/dI_s}.  (The date-0 response vectors are the first
#'   Jacobian row, \code{J[[o]][[i]][1, ]}.)
#' @export
hank_het_jacobian <- function(block, T_h,
                              inputs = c("r", "w"),
                              outputs = c("A", "C"),
                              delta_in = 1e-5, delta_va = 1e-6,
                              delta_d = 1e-6) {
  inputs  <- match.arg(inputs, several.ok = TRUE)
  outputs <- match.arg(outputs, several.ok = TRUE)
  Lam <- block$Lambda

  ## --- Step 2: expectation vectors E_s = Lambda^s y^o, s = 0 .. T-1 ---
  y_out <- list(A = .hank_mat_to_vec(block$a),   # per-agent savings
                C = .hank_mat_to_vec(block$c))   # per-agent consumption
  Elist <- setNames(vector("list", length(outputs)), outputs)
  for (o in outputs) {
    E <- vector("list", T_h)
    E[[1L]] <- y_out[[o]]                        # E_0
    for (s in 2L:T_h) E[[s]] <- as.numeric(Lam %*% E[[s - 1L]])
    Elist[[o]] <- E
  }

  J <- setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) matrix(0, T_h, T_h)), inputs)), outputs)

  for (i in inputs) {
    ## --- Step 1: backward sweep -> curlyY[[o]][s], curlyD[, s] ---
    sweep  <- .hank_curly_sweep(block, T_h, i, outputs,
                                delta_in, delta_va, delta_d)
    curlyY <- sweep$curlyY
    curlyD <- sweep$curlyD

    ## --- Steps 3-4: fake-news matrix F then Jacobian J ---
    for (o in outputs) {
      Fm <- matrix(0, T_h, T_h)
      Fm[1L, ] <- curlyY[[o]]                       # F[0, s] = curlyY[s]
      E <- Elist[[o]]
      for (tt in 2L:T_h) {
        Etm1 <- E[[tt - 1L]]                         # E_{t-1}
        Fm[tt, ] <- as.numeric(crossprod(curlyD, Etm1))  # <E_{t-1}, curlyD[,s]>
      }
      ## Diagonal cumulative sum.
      Jm <- matrix(0, T_h, T_h)
      Jm[1L, ] <- Fm[1L, ]
      for (tt in 2L:T_h) {
        Jm[tt, 1L] <- Fm[tt, 1L]
        Jm[tt, 2L:T_h] <- Jm[tt - 1L, 1L:(T_h - 1L)] + Fm[tt, 2L:T_h]
      }
      J[[o]][[i]] <- Jm
    }
  }
  J
}


#' Sequence-space DISTRIBUTION Jacobian of a het block via the fake-news
#' algorithm
#'
#' Extends the fake-news algorithm (\code{\link{hank_het_jacobian}}) to expose
#' the full distributional response \eqn{J^D[t, s, ] = dD_t/dI_s} (the change
#' in the \code{(n_e*n_a)}-vector cross-sectional distribution at date \code{t}
#' from an anticipated shock to input \code{i} at date \code{s}), instead of
#' aggregating it into scalar outputs \code{A}/\code{C}.
#'
#' Reuses the same backward sweep (\code{curlyD}) as \code{hank_het_jacobian}.
#' Distributions push forward under the TRANSPOSE of the steady-state forward
#' operator (\code{t(Lambda) \%*\% d}; see the file header), so the
#' distribution fake-news matrix cumulates by repeatedly applying
#' \code{t(Lambda)} to \code{curlyD[,s]} rather than by dotting against an
#' expectation vector (the aggregate-output equivalent of that projection).
#'
#' TIMING: \code{D_t} is the distribution ENTERING period \code{t} (a
#' predetermined state), so \code{D_1 = D_ss} always and row \code{t = 1} of
#' \code{J^D} is identically zero for every shock date \code{s} -- no
#' anticipated shock can move the very first, already-fixed, initial
#' distribution.  \code{curlyD[, s]} is the response of the policy used IN
#' period \code{s} (see \code{\link{.hank_curly_sweep}}: its \code{s = 1} term
#' is built from a policy step at the current period, exactly as
#' \code{curlyY[1]} is in \code{hank_het_jacobian}), which the forward operator
#' then turns into a distribution change one calendar period later, at
#' \code{t = s + 1}.  So the fake-news distribution response first appears at
#' \code{t = s + 1}, not \code{t = s} as for the aggregate outputs -- i.e. the
#' whole cumulation is the aggregate-Jacobian recursion shifted down by one row
#' (see the ND-vs-fake-news diagnostic in test-hank-dist-jacobian.R, which
#' pinned this down to a clean off-by-one before this comment was written).
#'
#' @inheritParams hank_het_jacobian
#'
#' @return Named list \code{JD[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_a)} with \code{JD[[i]][t, s, ] = dD_t/dI_s}.
#' @export
hank_het_dist_jacobian <- function(block, T_h,
                                   inputs = c("r", "w"),
                                   delta_in = 1e-5, delta_va = 1e-6,
                                   delta_d = 1e-6) {
  inputs <- match.arg(inputs, several.ok = TRUE)
  Lam <- block$Lambda
  n_cell <- block$n_e * block$n_a
  P <- function(x) as.numeric(Matrix::t(Lam) %*% x)   ## distribution push

  JD <- setNames(lapply(inputs, function(i) array(0, c(T_h, T_h, n_cell))),
                inputs)
  if (T_h < 2L) return(JD)   ## row t=1 (D_ss, fixed) is the only row; all-zero

  for (i in inputs) {
    ## --- Step 1: backward sweep -> curlyD[, s] (curlyY unused here) ---
    sweep  <- .hank_curly_sweep(block, T_h, i, outputs = "A",
                                delta_in, delta_va, delta_d)
    curlyD <- sweep$curlyD

    ## --- Step 3: distribution fake-news F^D, indexed by CALENDAR date t ---
    ## FD[[2]][, s] = curlyD[, s]   (first possible response: t = s + 1 = 2
    ##                                relative offset 1, hit when s = 1)
    ## FD[[t]][, s] = P(FD[[t-1]][, s])   for t >= 3
    ## (FD[[1]] would be t=1, always zero -- omitted; loop starts at t=2.)
    FD <- vector("list", T_h)
    FD[[2L]] <- curlyD
    for (tt in seq_len(T_h - 2L) + 2L) {              # tt = 3 .. T_h, empty if T_h < 3
      prev <- FD[[tt - 1L]]
      cur  <- matrix(0, n_cell, T_h)
      for (s in seq_len(T_h)) cur[, s] <- P(prev[, s])
      FD[[tt]] <- cur
    }

    ## --- Step 4: diagonal cumulation, vector-valued per (t, s) ---
    ## Mirrors hank_het_jacobian's J[t,s] = J[t-1,s-1] + F[t,s] exactly, just
    ## shifted down one row (t=1 row is the fixed, unresponsive D_ss and stays
    ## all-zero; the recursion proper starts at t=2).
    JD[[i]][2L, 1L, ] <- FD[[2L]][, 1L]               # JD[2, 1, ] = FD[2][, 1]
    for (s in seq_len(T_h - 1L) + 1L)                 # s = 2 .. T_h
      JD[[i]][2L, s, ] <- FD[[2L]][, s]               # JD[2,s,]=FD[2][,s] (t-1=1 row is 0)
    for (tt in seq_len(T_h - 2L) + 2L) {               # tt = 3 .. T_h, empty if T_h < 3
      JD[[i]][tt, 1L, ] <- FD[[tt]][, 1L]             # JD[t, 1, ] = FD[t][, 1]
      for (s in seq_len(T_h - 1L) + 1L) {              # s = 2 .. T_h
        JD[[i]][tt, s, ] <- JD[[i]][tt - 1L, s - 1L, ] + FD[[tt]][, s]
      }
    }
  }
  JD
}
