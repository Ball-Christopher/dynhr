#' Sparse three-asset Young forward operator
#'
#' Three-asset counterpart of \code{\link{hank_forward_operator2}} /
#' \code{\link{hank_forward_operator}}: builds the row-stochastic transition
#' matrix \code{Lambda} on the joint \code{(e, d, f, a)} grid by composing the
#' income transition \code{Pi} with three independent trilinear lotteries
#' (product lotteries, one per asset), each splitting mass between the two
#' bracketing gridpoints in proportion to distance, exactly conserving total
#' mass. It depends only on POLICIES, not on the marginal values that produced
#' them, so it is independent of \code{\link{hank_egm3_solve}} and can gate
#' distribution accounting against the discrete
#' \code{\link{hank_egm3_prototype_dist}} reference before a Jacobian is
#' attempted.
#'
#' @param d_pol,f_pol,a_pol Numeric \code{n_e x n_d x n_f x n_a} arrays: the
#'   next-period liquid, foreign and capital policies (conformable).
#' @param d_grid,f_grid,a_grid Numeric: the increasing liquid, foreign and
#'   capital grids the policies were solved on.
#' @param Pi Numeric \code{n_e x n_e}: row-stochastic income transition matrix.
#'
#' @return A sparse \code{dgCMatrix} (\code{Matrix} package), \code{N x N}
#'   with \code{N = n_e * n_d * n_f * n_a}, row-stochastic. Cell order matches
#'   the two-asset convention extended to three assets: income \code{e} is
#'   slowest, then domestic liquid \code{d}, then foreign \code{f}, with
#'   capital \code{a} fastest --
#'   \code{index(e, d, f, a) = a + n_a*((f-1) + n_f*((d-1) + n_d*(e-1)))}.
#'   Distributions update as \code{d_next = t(Lambda) \%*\% d}.
#' @seealso \code{\link{hank_forward_operator2}} (two-asset),
#'   \code{\link{hank_forward_operator}} (one-asset),
#'   \code{\link{hank_egm3_solve}}, \code{\link{hank_stationary_dist}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' hh <- hank_egm3_solve(dg, fg, ag, e, Pi, rd = .01, rf = .015, ra = .02,
#'                       beta = .97, eis = .5, chi1 = .2, phi1 = .1,
#'                       tol = 1e-5, maxit = 250)
#' Q <- hank_forward_operator3(hh$d, hh$f, hh$a, dg, fg, ag, Pi)
#' range(Matrix::rowSums(Q))
#' @export
hank_forward_operator3 <- function(d_pol,f_pol,a_pol,d_grid,f_grid,a_grid,Pi) {
  if(!is.array(d_pol)||length(dim(d_pol))!=4L||!identical(dim(d_pol),dim(f_pol))||!identical(dim(d_pol),dim(a_pol)))stop("hank_forward_operator3: policies must be conformable e x d x f x a arrays")
  ne<-dim(d_pol)[1];nd<-dim(d_pol)[2];nf<-dim(d_pol)[3];na<-dim(d_pol)[4];.hank_check_markov(Pi,ne,caller="hank_forward_operator3")
  if(any(!is.finite(c(d_pol,f_pol,a_pol)))||any(diff(d_grid)<=0)||any(diff(f_grid)<=0)||any(diff(a_grid)<=0))stop("hank_forward_operator3: finite policies and increasing grids required")
  # Package cell order extends the two-asset convention: income is slowest,
  # then domestic liquid and foreign assets, with capital fastest.
  N<-ne*nd*nf*na;ix<-function(e,d,f,a)a+na*((f-1)+nf*((d-1)+nd*(e-1)))
  lot<-function(g,z){if(length(g)==1L)return(c(1,1,1,0));j<-max(1L,min(findInterval(z,g),length(g)-1L));p<-min(1,max(0,(g[j+1]-z)/(g[j+1]-g[j])));c(j,j+1,p,1-p)}
  from<-integer(N*ne*8);to<-integer(N*ne*8);val<-numeric(N*ne*8);k<-0L
  for(e in seq_len(ne))for(d in seq_len(nd))for(f in seq_len(nf))for(a in seq_len(na)){ld<-lot(d_grid,d_pol[e,d,f,a]);lf<-lot(f_grid,f_pol[e,d,f,a]);la<-lot(a_grid,a_pol[e,d,f,a]);for(ep in seq_len(ne))for(jd in 1:2)for(jf in 1:2)for(ja in 1:2){k<-k+1L;from[k]<-ix(e,d,f,a);to[k]<-ix(ep,ld[jd],lf[jf],la[ja]);val[k]<-Pi[e,ep]*ld[jd+2]*lf[jf+2]*la[ja+2]}}
  Matrix::sparseMatrix(i=from[seq_len(k)],j=to[seq_len(k)],x=val[seq_len(k)],dims=c(N,N))
}

.hank_forward_apply3 <- function(block, x, transpose = FALSE) {
  hank_forward_apply3_cpp(block$d, block$f, block$a,
    block$d_grid, block$f_grid, block$a_grid, block$Pi, x, transpose)
}

## Pi_p/Pi_m default to the block's steady-state Pi (policy-only derivative);
## pass the two perturbed transition matrices to get the JOINT (policy, Pi)
## derivative a transition-probability Jacobian column needs at its shock date.
.hank_forward_direction3 <- function(block, dd, df, da, delta,
                                     Pi_p = NULL, Pi_m = NULL) {
  hank_forward_direction3_cpp(block$d, block$f, block$a, dd, df, da,
    block$d_grid, block$f_grid, block$a_grid, block$Pi, block$D, delta,
    Pi_p, Pi_m)
}

## Central derivative using the ACTUAL nonlinear policy and transition legs.
## This is the accounting-preserving path used by the three-asset fake-news
## sweep: at an active bound the two legs need not be symmetric around the
## steady policy, so reconstructing them from a central policy derivative can
## clip one leg and lose the exact asset first moments.
.hank_forward_legs3 <- function(block, step_p, step_m, step,
                                Pi_p = block$Pi, Pi_m = block$Pi) {
  hank_forward_legs3_cpp(step_p$d, step_p$f, step_p$a,
    step_m$d, step_m$f, step_m$a,
    block$d_grid, block$f_grid, block$a_grid,
    Pi_p, Pi_m, block$D, step)
}


#' Sequence-space DISTRIBUTION Jacobian of the three-asset block via the
#' fake-news algorithm
#'
#' Three-asset counterpart of \code{\link{hank_het_dist_jacobian}} (one-asset)
#' and \code{\link{hank_het2_dist_jacobian}} (two-asset): extends
#' \code{\link{hank_het3_jacobian}}'s fake-news algorithm to expose the full
#' distributional response \eqn{J^D[t, s, ] = dD_t/dI_s} (the change in the
#' \code{(n_e*n_d*n_f*n_a)}-vector cross-sectional distribution at date
#' \code{t} induced by an anticipated shock to aggregate input \code{i} at
#' date \code{s}), instead of aggregating it into scalar outputs
#' \code{A}/\code{C}/....
#'
#' Reuses the identical backward sweep (\code{.hank_curly_sweep3},
#' hence \code{curlyD}) as \code{\link{hank_het3_jacobian}}. Distributions
#' push forward under the TRANSPOSE of the steady-state joint Young operator
#' (\code{.hank_forward_apply3} with \code{transpose = TRUE}; matches
#' \code{t(Lambda) \%*\% d} in \code{\link{hank_forward_operator3}}'s
#' documented convention), so the distribution fake-news matrix cumulates by
#' repeatedly applying that transpose to \code{curlyD[, s]} rather than by
#' dotting against an expectation vector (the aggregate-output equivalent of
#' that projection, used by \code{\link{hank_het3_jacobian}}).
#'
#' TIMING: exactly as for the one- and two-asset blocks, \code{D_t} is the
#' distribution ENTERING period \code{t} (a predetermined state), so
#' \code{D_1 = D_ss} always and row \code{t = 1} of \code{J^D} is identically
#' zero for every shock date \code{s}. \code{curlyD[, s]} is the response of
#' the policy USED in period \code{s} (\code{s = 1} is the direct current-
#' period shock; \code{s >= 2} anticipation terms propagate via the joint
#' \code{(Vd, Vf, Va)} derivative -- see \code{.hank_curly_sweep3}),
#' which the forward operator turns into a distribution change one calendar
#' period later, at \code{t = s + 1}. So the whole cumulation is the
#' aggregate-Jacobian recursion (\code{\link{hank_het3_jacobian}}'s
#' \code{.assemble_jacobian}) shifted down by one row.
#'
#' The singleton-\code{f_grid} reduction's \code{px} column is an EXACT zero
#' here too, for the same reason as in \code{\link{hank_het3_jacobian}} (see
#' \code{.hank3_px_is_inert}): the reduction delegates to a kernel that has
#' never heard of \code{px}, so its distributional response is zero by the
#' same argument that makes its aggregate response zero, and is returned
#' rather than swept.
#'
#' @inheritParams hank_het3_jacobian
#'
#' @return Named list \code{JD[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_d*n_f*n_a)} with \code{JD[[i]][t, s, ] =
#'   dD_t/dI_s}.
#' @seealso \code{\link{hank_het3_dist_jacobian_nd}} (the numerical oracle
#'   this is validated against), \code{\link{hank_het_dist_jacobian}}
#'   (one-asset), \code{\link{hank_het3_jacobian}} (the aggregate three-asset
#'   Jacobian sharing this function's backward sweep)
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' JD <- hank_het3_dist_jacobian(blk, T_h = 3, inputs = c("rd", "w"))
#' dim(JD$rd)
#' @export
hank_het3_dist_jacobian <- function(block, T_h,
                                    inputs = NULL,
                                    delta_in = 1e-5, delta_v = 1e-6,
                                    delta_d = 1e-6, threads = NULL) {
  if (!inherits(block, "hank_het3_block"))
    stop("hank_het3_dist_jacobian: block must be hank_het3_block")
  if (!is.numeric(T_h) || length(T_h) != 1L || T_h < 1 || !is.finite(T_h))
    stop("hank_het3_dist_jacobian: T_h must be positive")
  if (is.null(inputs)) inputs <- .hank3_jac_inputs(block)
  inputs <- .hank3_check_inputs(block, inputs, "hank_het3_dist_jacobian")

  n_cell <- length(block$D)
  P <- function(x) .hank_forward_apply3(block, x, transpose = TRUE)

  JD <- setNames(lapply(inputs, function(i) array(0, c(T_h, T_h, n_cell))),
                inputs)
  if (T_h < 2L) return(JD)   ## row t=1 (D_ss, fixed) is the only row; all-zero

  # Exact-zero column on the reduction is never swept; see .hank3_px_is_inert.
  valid_inputs <- inputs[!(inputs == "px" & .hank3_px_is_inert(block))]

  for (i in valid_inputs) {
    ## --- Step 1: backward sweep -> curlyD[, s] (curlyY not needed here) ---
    sweep  <- .hank_curly_sweep3(block, T_h, i, character(0),
                                 delta_in, delta_v, delta_d, threads)
    curlyD <- sweep$curlyD

    ## --- Step 3: distribution fake-news F^D, indexed by CALENDAR date t ---
    ## FD[[2]][, s] = curlyD[, s]   (first possible response: t = s + 1 = 2,
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
    ## Mirrors hank_het3_jacobian's aggregate assembly exactly, just shifted
    ## down one row (t=1 row is the fixed, unresponsive D_ss and stays
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


#' Brute-force numerical-differentiation DISTRIBUTION Jacobian of the
#' three-asset block
#'
#' Reference sequence-space distribution Jacobian \eqn{J^D[t,s,] = dD_t/dI_s}
#' for the three-asset block, computed exactly like
#' \code{\link{hank_het3_jacobian_nd}} but keeping the full \code{Dpath}
#' (already returned unconditionally by \code{\link{hank_td3_nonlinear}})
#' instead of aggregating it into \code{D}/\code{F}/\code{A}/\code{C}/....
#' Used to validate \code{\link{hank_het3_dist_jacobian}}.
#'
#' @inheritParams hank_het3_jacobian_nd
#'
#' @return Named list \code{JD_nd[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_d*n_f*n_a)} with \code{JD_nd[[i]][t, s, ] =
#'   dD_t/dI_s} (central difference).
#' @seealso \code{\link{hank_het3_dist_jacobian}} (the fake-news distribution
#'   Jacobian this validates), \code{\link{hank_het_dist_jacobian_nd}}
#'   (one-asset)
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' JD_nd <- hank_het3_dist_jacobian_nd(blk, T_h = 3, inputs = c("rd", "w"),
#'                                    delta = 3e-6)
#' dim(JD_nd$rd)
#' @export
hank_het3_dist_jacobian_nd <- function(block, T_h, inputs = NULL,
                                       delta = 1e-5, threads = NULL) {
  if (!inherits(block, "hank_het3_block"))
    stop("hank_het3_dist_jacobian_nd: block must be hank_het3_block")
  if (!is.numeric(T_h) || length(T_h) != 1L || T_h < 1 || !is.finite(T_h))
    stop("hank_het3_dist_jacobian_nd: T_h must be positive")
  if (is.null(inputs)) inputs <- .hank3_jac_inputs(block)
  inputs <- .hank3_check_inputs(block, inputs, "hank_het3_dist_jacobian_nd")

  n_cell <- length(block$D)
  base <- list(rd = rep(block$rd, T_h), rf = rep(block$rf, T_h),
              ra = rep(block$ra, T_h), w = rep(block$w, T_h),
              px = rep(if (is.null(block$px)) 1 else block$px, T_h),
              Tr = rep(.hank_block_tr(block), T_h))
  # Name the path arguments in FULL -- see hank_het3_jacobian_nd's identical
  # comment: partial matching ("rd" -> "rd_path") is a coin flip once other
  # formals share the prefix (px_path/keep_policies).
  as_paths <- function(p) setNames(p, paste0(names(p), "_path"))
  run <- function(p, pi_paths) do.call(hank_td3_nonlinear,
    c(list(block = block, T_h = T_h, pi_input_paths = pi_paths, threads = threads),
      as_paths(p)))

  JD_nd <- setNames(lapply(inputs, function(i) array(0, c(T_h, T_h, n_cell))),
                    inputs)
  for (i in inputs) {
    if (i == "px" && .hank3_px_is_inert(block)) next   # exact zero column
    is_pi <- .hank3_is_pi_input(block, i)
    for (s in seq_len(T_h)) {
      p <- m <- base; pip <- pim <- NULL
      if (is_pi) {
        x0 <- rep(block$Pi_inputs[[i]], T_h)
        xp <- x0; xp[s] <- xp[s] + delta
        xm <- x0; xm[s] <- xm[s] - delta
        pip <- setNames(list(xp), i); pim <- setNames(list(xm), i)
      } else {
        p[[i]][s] <- p[[i]][s] + delta
        m[[i]][s] <- m[[i]][s] - delta
      }
      op <- run(p, pip); om <- run(m, pim)
      dD <- (op$Dpath - om$Dpath) / (2 * delta)   # (n_cell x T_h), col t
      for (tt in seq_len(T_h)) JD_nd[[i]][tt, s, ] <- dD[, tt]
    }
  }
  JD_nd
}
