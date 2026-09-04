#' Small-grid three-asset household prototype
#'
#' A deliberately SEPARATE reference solver for domestic liquid claims
#' \code{d}, gross foreign assets \code{f}, and illiquid domestic capital
#' \code{a} -- independent of \code{\link{hank_egm3_solve}}'s EGM machinery,
#' so it can act as a discrete-choice oracle for it. At every state it does an
#' exhaustive grid search over the discrete next-period liquid/foreign/capital
#' action space
#' (value function iteration, no interpolation, no first-order condition),
#' which is slow but has no smoothness or monotonicity assumption to get
#' wrong. Capped at 512 states: this establishes and sanity-checks the
#' three-asset budget/friction structure before the production EGM/Jacobian
#' implementation is trusted to get it right on a larger grid.
#'
#' @param d_grid Numeric increasing grid, length >= 2: DOMESTIC LIQUID claims.
#'   MAY include negative values (\code{d_grid[1]} is the borrowing/liquid
#'   floor) -- the only one of the three grids that may.
#' @param f_grid Numeric increasing grid, length >= 2 (or a single
#'   non-negative value to hold foreign assets fixed, e.g. the reduction
#'   oracle at \code{0} -- see \code{\link{hank_egm3_solve}}): gross FOREIGN
#'   assets. Must be non-negative (no short foreign position).
#' @param a_grid Numeric increasing grid, length >= 2, non-negative (no short
#'   position): illiquid domestic CAPITAL.
#' @param y Numeric length-\code{n_e}: labour income by income state.
#' @param Pi Numeric \code{n_e x n_e}: row-stochastic income transition matrix.
#' @param rd,rf,ra Numeric: returns on liquid, foreign and capital holdings.
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param chi0,chi1,chi2 Capital adjustment-cost parameters (see
#'   \code{\link{.hank_psi}}): \code{chi0 > 0}, \code{chi1 >= 0},
#'   \code{chi2 > 1}.
#' @param phi0,phi1,phi2 Foreign-portfolio adjustment-cost parameters, same
#'   functional form and constraints as \code{chi0}/\code{chi1}/\code{chi2}.
#' @param tol Convergence tolerance on \eqn{\max|V_{new} - V|} across a value
#'   iteration.
#' @param maxit Maximum number of value iterations.
#'
#' @return An object of class \code{hank_egm3_prototype}: a list with the
#'   echoed calibration (\code{d_grid}, \code{f_grid}, \code{a_grid},
#'   \code{y}, \code{Pi}, \code{rd}, \code{rf}, \code{ra}, \code{beta},
#'   \code{eis}, \code{chi0}, \code{chi1}, \code{chi2}, \code{phi0},
#'   \code{phi1}, \code{phi2}) plus:
#'   \describe{
#'     \item{\code{d}, \code{f}, \code{a}, \code{c}}{Policies (each
#'       \code{n_e x n_d x n_f x n_a}): the argmax next-period liquid, foreign
#'       and capital holdings and the implied consumption.}
#'     \item{\code{V}}{Converged value function, same shape.}
#'     \item{\code{iterations}}{Number of value iterations run.}
#'     \item{\code{converged}}{Logical: whether \code{tol} was met.}
#'     \item{\code{state_count}}{\code{n_e * n_d * n_f * n_a} (an error is
#'       raised above 512).}
#'   }
#' @seealso \code{\link{hank_egm3_prototype_dist}} (the matching stationary
#'   distribution), \code{\link{hank_egm3_solve}} (the production EGM solver
#'   this establishes the budget for)
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' x <- hank_egm3_prototype(dg, fg, ag, e, Pi, rd = .01, rf = .015, ra = .02,
#'                          beta = .97, eis = .5, chi1 = .2, phi1 = .1,
#'                          tol = 1e-5, maxit = 250)
#' x$converged
#' @keywords internal
#' @export
hank_egm3_prototype <- function(d_grid, f_grid, a_grid, y, Pi, rd, rf, ra,
                                beta, eis, chi0=.25, chi1=6.5, chi2=2,
                                phi0=.25, phi1=.5, phi2=2,
                                tol=1e-9, maxit=2000L) {
  ck <- function(x,nm,nn=FALSE,single=FALSE) { if(!is.numeric(x)||length(x)<if(single)1 else 2||any(!is.finite(x))||length(x)>1&&any(diff(x)<=0)||nn&&min(x)<0) stop("hank_egm3_prototype: invalid ",nm) }
  ck(d_grid,"d_grid");ck(f_grid,"f_grid",TRUE,TRUE);ck(a_grid,"a_grid",TRUE)
  if(length(f_grid)==1L && f_grid!=0)stop("hank_egm3_prototype: singleton foreign grid must be zero (reduction oracle)")
  .hank_check_markov(Pi,length(y),caller="hank_egm3_prototype")
  if(any(!is.finite(c(rd,rf,ra,beta,eis,chi0,chi1,chi2,phi0,phi1,phi2,tol,maxit)))||beta<=0||beta>=1||eis<=0||chi0<=0||chi1<0||chi2<=1||phi0<=0||phi1<0||phi2<=1||tol<=0||maxit<1)stop("hank_egm3_prototype: invalid parameters")
  ne<-length(y);nd<-length(d_grid);nf<-length(f_grid);na<-length(a_grid);ns<-ne*nd*nf*na
  if(ns>512L)stop("hank_egm3_prototype: prototype cap is 512 states; no paper-scale use")
  u <- function(c) if(eis==1)log(c) else (c^(1-1/eis)-1)/(1-1/eis)
  V<-array(0,c(ne,nd,nf,na)); D<-F<-A<-C<-array(NA_real_,dim(V)); ok<-FALSE
  for(it in seq_len(as.integer(maxit))) { W<-V
    for(e in seq_len(ne))for(id in seq_len(nd))for(jf0 in seq_len(nf))for(ia in seq_len(na)) {
      best<- -Inf;ans<-rep(NA_real_,4); r<-y[e]+(1+rd)*d_grid[id]+(1+rf)*f_grid[jf0]+(1+ra)*a_grid[ia]
      for(jd in seq_len(nd))for(jf in seq_len(nf))for(ja in seq_len(na)){ psi_a<-.hank_psi(a_grid[ja],a_grid[ia],ra,chi0,chi1,chi2)$Psi;psi_f<-.hank_psi(f_grid[jf],f_grid[jf0],rf,phi0,phi1,phi2)$Psi; cc<-r-d_grid[jd]-f_grid[jf]-a_grid[ja]-psi_a-psi_f; if(cc>0){v<-u(cc)+beta*sum(Pi[e,]*V[,jd,jf,ja]);if(v>best){best<-v;ans<-c(d_grid[jd],f_grid[jf],a_grid[ja],cc)}} }
      if(!is.finite(best))stop("hank_egm3_prototype: infeasible state");W[e,id,jf0,ia]<-best;D[e,id,jf0,ia]<-ans[1];F[e,id,jf0,ia]<-ans[2];A[e,id,jf0,ia]<-ans[3];C[e,id,jf0,ia]<-ans[4]
    }
    if(max(abs(W-V))<tol){V<-W;ok<-TRUE;break};V<-W
  }
  structure(list(d_grid=d_grid,f_grid=f_grid,a_grid=a_grid,y=y,Pi=Pi,rd=rd,rf=rf,ra=ra,beta=beta,eis=eis,chi0=chi0,chi1=chi1,chi2=chi2,phi0=phi0,phi1=phi1,phi2=phi2,d=D,f=F,a=A,c=C,V=V,iterations=it,converged=ok,state_count=ns),class=c("hank_egm3_prototype", "hank_block"))
}

#' Stationary distribution for the three-asset prototype
#'
#' Builds the DENSE \code{N x N} transition matrix for a
#' \code{\link{hank_egm3_prototype}} solve by explicit trilinear lottery
#' interpolation (each of the three assets splits mass between its two
#' bracketing gridpoints, composed with the income transition \code{Pi}), then
#' iterates it to the invariant distribution. Deliberately dense rather than
#' sparse -- it is an independent accounting oracle for
#' \code{\link{hank_forward_operator3}} / \code{\link{hank_stationary_dist}}'s
#' sparse production path, and is only ever run at the prototype's small
#' (<= 512-state) cap.
#'
#' @param x A \code{\link{hank_egm3_prototype}} solve.
#' @param tol Convergence tolerance on \eqn{\max|D_{new} - D|} across a power
#'   iteration.
#' @param maxit Maximum number of power iterations.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{D}}{Stationary distribution, length-\code{N} vector
#'       (\code{N = x$state_count}), sums to 1.}
#'     \item{\code{Q}}{The dense \code{N x N} row-stochastic transition
#'       matrix.}
#'     \item{\code{converged}}{Logical: whether \code{tol} was met.}
#'     \item{\code{iterations}}{Number of power iterations run.}
#'     \item{\code{mass}}{\code{sum(D)}, a conservation check (should equal
#'       1).}
#'     \item{\code{D_agg}, \code{F_agg}, \code{A_agg}, \code{C_agg}}{Aggregate
#'       liquid, foreign and capital holdings and consumption,
#'       \code{sum(D * policy)}.}
#'   }
#' @seealso \code{\link{hank_egm3_prototype}},
#'   \code{\link{hank_forward_operator3}}, \code{\link{hank_stationary_dist}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' x <- hank_egm3_prototype(dg, fg, ag, e, Pi, rd = .01, rf = .015, ra = .02,
#'                          beta = .97, eis = .5, chi1 = .2, phi1 = .1,
#'                          tol = 1e-5, maxit = 250)
#' z <- hank_egm3_prototype_dist(x)
#' z$mass
#' @keywords internal
#' @export
hank_egm3_prototype_dist <- function(x, tol=1e-12, maxit=100000L) {
  if(!inherits(x,"hank_egm3_prototype"))stop("hank_egm3_prototype_dist: x must be a prototype solve")
  ne<-length(x$y);nd<-length(x$d_grid);nf<-length(x$f_grid);na<-length(x$a_grid);N<-x$state_count
  ix <- function(e,d,f,a) e+ne*((d-1)+nd*((f-1)+nf*(a-1)))
  wt <- function(g,z) { if(length(g)==1L)return(c(1,1,1,0));j<-findInterval(z,g);j<-max(1,min(j,length(g)-1));q<-(z-g[j])/(g[j+1]-g[j]);c(j,j+1,1-q,q) }
  Q<-matrix(0,N,N)
  for(e in seq_len(ne))for(d in seq_len(nd))for(f in seq_len(nf))for(a in seq_len(na)) {
    wd<-wt(x$d_grid,x$d[e,d,f,a]);wf<-wt(x$f_grid,x$f[e,d,f,a]);wa<-wt(x$a_grid,x$a[e,d,f,a]);r<-ix(e,d,f,a)
    for(ep in seq_len(ne))for(jd in 1:2)for(jf in 1:2)for(ja in 1:2) Q[r,ix(ep,wd[jd],wf[jf],wa[ja])]<-Q[r,ix(ep,wd[jd],wf[jf],wa[ja])]+x$Pi[e,ep]*wd[jd+2]*wf[jf+2]*wa[ja+2]
  }
  D<-rep(1/N,N);ok<-FALSE
  for(it in seq_len(as.integer(maxit))){Dn<-as.numeric(D%*%Q);if(max(abs(Dn-D))<tol){D<-Dn;ok<-TRUE;break};D<-Dn}
  agg<-function(z)sum(D*as.vector(z))
  list(D=D,Q=Q,converged=ok,iterations=it,mass=sum(D),D_agg=agg(x$d),F_agg=agg(x$f),A_agg=agg(x$a),C_agg=agg(x$c))
}
