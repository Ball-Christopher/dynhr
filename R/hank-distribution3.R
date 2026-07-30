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
