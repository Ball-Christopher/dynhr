#' Construct a three-asset heterogeneous-agent household block
#'
#' Solves the three-asset household problem (\code{\link{hank_egm3_solve}})
#' and its stationary joint distribution
#' (\code{\link{hank_forward_operator3}} + \code{\link{hank_stationary_dist}})
#' at fixed aggregate prices \code{(rd, rf, ra, w)}, and packages everything a
#' later nonlinear transition/Jacobian implementation needs. Experimental
#' Stage-4 R-reference counterpart of \code{\link{hank_het2_block}}, adding a
#' foreign asset alongside the domestic liquid/capital pair.
#'
#' @param d_grid Numeric: increasing DOMESTIC LIQUID grid, length >= 3;
#'   \code{d_grid[1]} is the liquid floor.
#' @param f_grid Numeric: increasing FOREIGN asset grid, length >= 3, or the
#'   singleton \code{0} to route to the exact two-asset reduction (see
#'   \code{\link{hank_egm3_solve}}).
#' @param a_grid Numeric: increasing domestic CAPITAL grid, length >= 3.
#' @param Pi Numeric \code{n_e x n_e}: row-stochastic income transition matrix.
#' @param e Numeric length-\code{n_e}: income levels (see
#'   \code{\link{hank_income_rouwenhorst}}).
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param rd,rf,ra Steady-state liquid, foreign and capital returns.
#' @param w Steady-state wage; income is \code{y = w * e}.
#' @param chi0,chi1,chi2 Capital adjustment-cost parameters (see
#'   \code{\link{.hank_psi}}; reference calibration \code{0.25 / 6.5 / 2}).
#' @param phi0,phi1,phi2 Foreign-portfolio adjustment-cost parameters, same
#'   functional form and constraints as \code{chi0}/\code{chi1}/\code{chi2}.
#' @param px Steady-state world price of foreign claims (default \code{1});
#'   see \code{\link{hank_egm3_solve}} for where it enters and why. It is a
#'   first-class INPUT here, not a fixed parameter: it is the level around
#'   which \code{\link{hank_td3_nonlinear}}'s \code{px_path} and the
#'   \code{"px"} column of \code{\link{hank_het3_jacobian}} are taken.
#' @param Pi_fn Optional function rebuilding the income transition matrix from
#'   named transition-probability inputs (e.g. the \code{Pi_fn(f, s)} returned
#'   by \code{\link{hank_employment_income}}, whose \code{f} and \code{s} are
#'   the job-finding and separation rates). Supplying it -- together with
#'   \code{Pi_inputs} -- makes those inputs first-class perturbable aggregate
#'   inputs of the block, i.e. extra columns of
#'   \code{\link{hank_het3_jacobian}} and extra paths of
#'   \code{\link{hank_td3_nonlinear}}, alongside the prices. \code{e} must stay
#'   FIXED as the transition inputs move (only \code{Pi} responds), which is
#'   the convention \code{\link{hank_employment_income}} is built to.
#' @param Pi_inputs Named list of the steady-state values of \code{Pi_fn}'s
#'   arguments (e.g. \code{list(f = f, s = s)}); required with \code{Pi_fn}
#'   (and only then). Names must not collide with the price inputs
#'   (\code{"rd"}, \code{"rf"}, \code{"ra"}, \code{"w"}, \code{"px"}), and
#'   \code{do.call(Pi_fn, Pi_inputs)} must reproduce \code{Pi} (checked here --
#'   an inconsistent pair would silently differentiate around a point the
#'   block is not sitting at).
#' @param Tr Finite numeric scalar: a lump-sum transfer entering income
#'   additively, \code{y = w * e + Tr * Tr_incidence}. Default \code{0}
#'   restores the transfer-free household BYTE-FOR-BYTE (adding an exact zero
#'   to income). This closes the last input asymmetry against the one- and
#'   two-asset tiers, which have carried a transfer for some time.
#' @param Tr_incidence Optional finite numeric length-\code{n_e} incidence
#'   weight, normalised here so its mean under the \code{Pi}-invariant
#'   distribution is 1 (see \code{\link{.hank_normalize_incidence}}; the
#'   normalisation is against the invariant distribution rather than the
#'   solved \code{D} precisely so it is not circular). \code{NULL} (default)
#'   is uniform incidence. Indexed by the INCOME state only. Note the weight
#'   is block STATE, not a perturbable input: \code{\link{hank_het3_jacobian}}
#'   differentiates the scalar \code{Tr} holding the rule fixed, so the
#'   \code{"Tr"} column is the aggregate response to a transfer distributed by
#'   this rule -- change the rule and you have a different block.
#' @param tol,maxit Passed to \code{\link{hank_egm3_solve}}.
#' @param dist_tol,dist_maxit Tolerance and iteration cap for the stationary
#'   distribution (\code{\link{hank_stationary_dist}}, R backend).
#' @param backend,threads Passed straight through to
#'   \code{\link{hank_egm3_solve}} (by name, not position -- see that
#'   function's own defaults and documentation). Defaults reproduce this
#'   function's pre-existing behaviour exactly.
#' @param relax Damping weight in \eqn{(0, 1]} on the marginal-value update,
#'   forwarded BY NAME to \code{\link{hank_egm3_solve}}; default \code{.5}
#'   reproduces this function's pre-existing behaviour exactly.
#'
#'   \strong{Convergence and \code{relax}: which direction actually helps.}
#'   Measured on a 14-income-state, \code{(n_d,n_f,n_a)=(8,4,8)} fixture
#'   (3,584 states): \code{relax=1} converges in 927 iterations,
#'   \code{relax=.75} in 1,236, \code{relax=.5} (the default) in 1,854, and
#'   \code{relax=.25}, \code{.1} and \code{.05} all FAIL to converge inside a
#'   2,000-iteration cap. Iterations scale as \code{1/relax}: the damped
#'   update \code{V <- (1-relax)*V + relax*V_new} has contraction modulus
#'   \code{1 - relax*(1-rho)} with \code{rho ~ beta}, so damping can only ever
#'   SLOW a monotone contraction -- it helps against oscillation, and on this
#'   class of fixture the oscillation is transient (an early rise in the value
#'   gap around iteration 400 as the active set shuffles, well before the tail
#'   collapse). \strong{A user reaching for \code{relax} to fix a
#'   non-convergence should raise it toward \code{1}, not lower it}: lowering
#'   it below ~0.5 can prevent convergence entirely rather than help. Warm
#'   continuation (see \code{Vd_init}/\code{Vf_init}/\code{Va_init} below) is
#'   worth a further ~3.5x and is what actually traverses a genuine
#'   continuation frontier where the cold solve does not converge inside the
#'   iteration cap at all.
#' @param Vd_init,Vf_init,Va_init Optional \code{n_e x n_d x n_f x n_a} initial
#'   marginal values, forwarded BY NAME to \code{\link{hank_egm3_solve}}.
#'   \code{NULL} (default) reproduces this function's pre-existing
#'   from-scratch initialization exactly. This is the warm-start API the
#'   Stage-4 calibration driver needs: the block already returns its own
#'   converged \code{Vd}/\code{Vf}/\code{Va}, so a later, nearby-parameter
#'   evaluation can continue from them instead of solving cold --
#'   \preformatted{
#'   blk1 <- hank_het3_block(d_grid, f_grid, a_grid, Pi, e, beta = b0, ...)
#'   blk2 <- hank_het3_block(d_grid, f_grid, a_grid, Pi, e, beta = b1, ...,
#'                           Vd_init = blk1$Vd, Vf_init = blk1$Vf,
#'                           Va_init = blk1$Va)
#'   }
#'   which converges in materially fewer iterations than a cold
#'   \code{blk2}, and can converge at all where a cold solve at \code{b1}
#'   would exhaust \code{maxit} (see \code{relax} above). Invalid arrays
#'   (wrong dimension, non-finite) fail with \code{\link{hank_egm3_solve}}'s
#'   own validation message -- this function does not re-validate them.
#' @param strict Logical, default \code{TRUE} (the pre-existing behaviour):
#'   a non-converged EGM solve or stationary distribution \code{stop()}s, with
#'   a message reporting \code{iterations}, \code{last_value_gap} (or the
#'   distribution iteration count), \code{relax}, \code{maxit}/\code{dist_maxit}
#'   and \code{tol}/\code{dist_tol}, plus the steer above (raise \code{relax}
#'   toward \code{1} and/or warm-start; do not lower \code{relax}). When
#'   \code{strict = FALSE}, the same failure instead returns an object of
#'   class \code{\link{hank_het3_block_failed}} -- deliberately NOT class
#'   \code{hank_het3_block} -- carrying the underlying household solve, the
#'   gaps, \code{relax}, \code{maxit}/\code{dist_maxit}, and (when cheaply
#'   obtainable) \code{\link{hank_euler3_residual}}-style active-regime
#'   counts. A calibration root-finder can inspect this to distinguish slow
#'   convergence from a broken active set, but the DIFFERENT class means a
#'   non-converged household can never be silently fed to
#'   \code{\link{hank_het3_jacobian}} or \code{\link{hank_td3_nonlinear}},
#'   both of which require class \code{hank_het3_block} and will error on it.
#'
#' @return An object of class \code{hank_het3_block}: everything
#'   \code{\link{hank_egm3_solve}} returns (echoed calibration, policies
#'   \code{d}/\code{f}/\code{a}/\code{c}, marginal values
#'   \code{Vd}/\code{Vf}/\code{Va}, \code{iterations}, \code{converged},
#'   \code{last_value_gap}, \code{last_policy_gap}, \code{state_count},
#'   \code{backend}, \code{threads}, \code{elapsed}), plus:
#'   \describe{
#'     \item{\code{e}, \code{w}}{The income levels and wage supplied.}
#'     \item{\code{Lambda}}{Sparse forward operator
#'       (\code{\link{hank_forward_operator3}}).}
#'     \item{\code{D}}{Stationary joint distribution (length
#'       \code{state_count} vector, sums to 1).}
#'     \item{\code{dist_converged}}{Logical: whether the distribution
#'       iteration met \code{dist_tol}.}
#'     \item{\code{D_agg}, \code{F_agg}, \code{A_agg}, \code{C}}{Aggregate
#'       liquid, foreign and capital holdings and consumption,
#'       \code{sum(D * policy)}.}
#'     \item{\code{CHI}, \code{PHI}}{Aggregate CAPITAL and FOREIGN adjustment
#'       resources, \code{sum(D * cost)}. These are the goods absorbed by
#'       rebalancing rather than consumed or invested; the goods-market
#'       identity the block feeds needs them explicitly, which is why the
#'       kernel now reports the per-cell costs. Counterpart of
#'       \code{\link{hank_het2_block}}'s \code{CHI}.}
#'     \item{\code{n_e}, \code{n_d}, \code{n_f}, \code{n_a}}{Grid sizes.}
#'     \item{\code{elapsed_solve}, \code{elapsed_dist}}{Wall-clock seconds
#'       (\code{proc.time()[["elapsed"]]} differences) spent in the household
#'       solve (\code{\link{hank_egm3_solve}}'s own \code{elapsed}) and in the
#'       stationary-distribution iteration (\code{\link{hank_stationary_dist}})
#'       respectively, for the run manifest
#'       (\code{\link{hank_het3_manifest}}).}
#'   }
#'   With \code{strict = TRUE} (default) errors (rather than warns) if either
#'   the household EGM or the stationary distribution fails to converge; with
#'   \code{strict = FALSE} returns a \code{\link{hank_het3_block_failed}}
#'   instead (see \code{strict} above).
#' @seealso \code{\link{hank_het2_block}} (two-asset),
#'   \code{\link{hank_het_block}} (one-asset), \code{\link{hank_egm3_solve}},
#'   \code{\link{hank_td3_nonlinear}}, \code{\link{hank_het3_jacobian}},
#'   \code{\link{hank_het3_block_failed}}, \code{\link{hank_egm3_regrid_values}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' c(D = blk$D_agg, F = blk$F_agg, A = blk$A_agg, C = blk$C)
#'
#' ## Warm-start idiom: reuse the converged marginal values at a NEARBY
#' ## parameter point instead of solving cold (see Vd_init/Vf_init/Va_init).
#' blk2 <- hank_het3_block(dg, fg, ag, Pi, e, beta = .971, eis = .5,
#'                         rd = .01, rf = .015, ra = .02, w = 1,
#'                         chi1 = .2, phi1 = .1, maxit = 250,
#'                         Vd_init = blk$Vd, Vf_init = blk$Vf, Va_init = blk$Va)
#' blk2$iterations <= blk$iterations
#' @export
hank_het3_block <- function(d_grid,f_grid,a_grid,Pi,e,beta,eis,rd,rf,ra,w,
                            chi0=.25,chi1=6.5,chi2=2,phi0=.25,phi1=.5,phi2=2,
                            px=1,
                            tol=1e-5,maxit=300,dist_tol=1e-12,dist_maxit=100000,
                            backend=getOption("dynhr.hank3_backend","cpp"),
                            threads=NULL,
                            relax=.5,Vd_init=NULL,Vf_init=NULL,Va_init=NULL,
                            strict=TRUE,Pi_fn=NULL,Pi_inputs=NULL,
                            Tr=0,Tr_incidence=NULL) {
  .hank_check_pi_fn(Pi_fn,Pi_inputs,Pi,
                    reserved=c("rd","rf","ra","w","px","Tr"),
                    caller="hank_het3_block")
  if(!is.numeric(Tr)||length(Tr)!=1L||!is.finite(Tr))
    stop("hank_het3_block: 'Tr' must be a finite numeric scalar (the ",
         "lump-sum transfer; 0 restores the transfer-free household).")
  omega<-.hank_normalize_incidence(Tr_incidence,e,Pi,"hank_het3_block")
  ## y + 0 is bit-identical to y and Tr * 1 to Tr, so the default path here is
  ## byte-for-byte the pre-transfer block -- the same argument the one-asset
  ## tier relies on (see .hank_block_omega).
  y<-w*e+Tr*omega
  ## backend/threads/relax/Vd_init/Vf_init/Va_init passed BY NAME:
  ## hank_egm3_solve is otherwise called positionally here, and partial
  ## argument matching on a positional call is exactly the trap that bit A4
  ## ("rd" partial-matched "rd_path"). Naming these keeps them immune to any
  ## future reordering of the positional block above.
  hh<-hank_egm3_solve(d_grid,f_grid,a_grid,y,Pi,rd,rf,ra,beta,eis,chi0,chi1,chi2,phi0,phi1,phi2,tol,maxit,
    relax=relax,Vd_init=Vd_init,Vf_init=Vf_init,Va_init=Va_init,
    px=px,backend=backend,threads=threads)
  if(!hh$converged){
    msg<-.hank3_nonconverge_msg("EGM",hh$iterations,maxit,relax,tol,
      sprintf("last_value_gap=%s, last_policy_gap=%s",
        format(hh$last_value_gap,digits=4),format(hh$last_policy_gap,digits=4)))
    if(strict)stop(msg)
    return(.hank3_block_failed(hh,msg,"egm",d_grid,f_grid,a_grid,Pi,e,w,relax,maxit,tol,
      NA_integer_,dist_maxit,dist_tol))
  }
  Q<-hank_forward_operator3(hh$d,hh$f,hh$a,d_grid,f_grid,a_grid,Pi)
  t0_dist<-proc.time()[["elapsed"]]
  sd<-hank_stationary_dist(Q,tol=dist_tol,maxit=dist_maxit,backend="R")
  elapsed_dist<-proc.time()[["elapsed"]]-t0_dist
  if(!sd$converged){
    msg<-.hank3_nonconverge_msg("stationary distribution",sd$iterations,dist_maxit,relax,dist_tol,NULL)
    if(strict)stop(msg)
    return(.hank3_block_failed(hh,msg,"stationary_dist",d_grid,f_grid,a_grid,Pi,e,w,relax,maxit,tol,
      sd$iterations,dist_maxit,dist_tol,Lambda=Q))
  }
  flat<-function(z)as.vector(aperm(z,c(4,3,2,1)))
  structure(c(hh,list(e=e,w=w,Tr=Tr,Tr_incidence=omega,Lambda=Q,D=sd$d,dist_converged=sd$converged,
    D_agg=sum(sd$d*flat(hh$d)),F_agg=sum(sd$d*flat(hh$f)),A_agg=sum(sd$d*flat(hh$a)),C=sum(sd$d*flat(hh$c)),
    CHI=sum(sd$d*flat(hh$chi)),PHI=sum(sd$d*flat(hh$phi)),
    n_e=length(e),n_d=length(d_grid),n_f=length(f_grid),n_a=length(a_grid),
    Pi_fn=Pi_fn,Pi_inputs=Pi_inputs,
    elapsed_solve=hh$elapsed,elapsed_dist=elapsed_dist)),class="hank_het3_block")
}

## A non-convergence message a calibration driver can act on: iterations,
## gaps, relax, maxit and tol, plus a one-line steer that matches the
## measured direction (raise relax toward 1 and/or warm-start; do NOT lower
## relax -- see hank_het3_block's own `relax` documentation for the
## measurement). `what` names which solve failed ("EGM" or "stationary
## distribution"); `extra` is an optional pre-formatted trailing clause
## (the EGM path adds the value/policy gaps that the distribution path does
## not have).
.hank3_nonconverge_msg <- function(what,iterations,cap,relax,tol,extra) {
  paste0("hank_het3_block: ",what," did not converge (iterations=",iterations,
         "/",cap,", tol=",format(tol,digits=4),", relax=",format(relax,digits=3),
         if(!is.null(extra))paste0(", ",extra) else "",
         "). Raise 'relax' toward 1 and/or warm-start from a neighbouring ",
         "converged solve's Vd/Vf/Va (Vd_init/Vf_init/Va_init) -- lowering ",
         "'relax' below ~0.5 can PREVENT convergence entirely on this class ",
         "of fixture (contraction modulus 1 - relax*(1-rho), rho ~ beta; ",
         "damping only ever slows a monotone contraction, it does not fix ",
         "one). See ?hank_het3_block.")
}

## Constructor for the strict=FALSE failure return. Carries the underlying
## household solve `hh` (whatever hank_egm3_solve produced at the point
## iteration stopped -- policies, marginal values, gaps), the failure stage,
## the diagnostic message, and -- when cheaply obtainable -- regime counts
## from hank_euler3_residual(), which only needs the policy/grid fields `hh`
## already carries regardless of convergence. Wrapped in tryCatch because a
## severely broken active set (the very thing this path exists to report) is
## exactly the case where a diagnostic built ON TOP of the failed policies
## might itself choke; a failed-object print should never itself error.
.hank3_block_failed <- function(hh,msg,stage,d_grid,f_grid,a_grid,Pi,e,w,relax,
                                maxit,tol,dist_iterations,dist_maxit,dist_tol,
                                Lambda=NULL) {
  regime_counts <- tryCatch(hank_euler3_residual(hh)$regime$counts,
                            error=function(err)NULL)
  structure(list(
    hh=hh,stage=stage,message=msg,
    d_grid=d_grid,f_grid=f_grid,a_grid=a_grid,Pi=Pi,e=e,w=w,
    iterations=hh$iterations,
    last_value_gap=hh$last_value_gap,last_policy_gap=hh$last_policy_gap,
    dist_iterations=dist_iterations,
    dist_converged=if(stage=="stationary_dist")FALSE else NA,
    relax=relax,maxit=maxit,tol=tol,dist_maxit=dist_maxit,dist_tol=dist_tol,
    Lambda=Lambda,regime_counts=regime_counts
  ),class="hank_het3_block_failed")
}

#' A failed three-asset household-block solve
#'
#' The object \code{\link{hank_het3_block}} returns, instead of erroring, when
#' called with \code{strict = FALSE} and either the household EGM or the
#' stationary distribution fails to converge. Deliberately NOT class
#' \code{hank_het3_block}: \code{\link{hank_het3_jacobian}} and
#' \code{\link{hank_td3_nonlinear}} both require that class and will error on
#' this one, so a non-converged household can never be silently consumed as
#' if it were a valid fixed point.
#'
#' @param x A \code{hank_het3_block_failed} object.
#' @param ... Unused; present for the generic \code{print} signature.
#' @return \code{x}, invisibly.
#' @details Fields: \code{hh} (the underlying \code{\link{hank_egm3_solve}}
#'   result at the point iteration stopped -- policies, marginal values,
#'   \code{iterations}, \code{last_value_gap}, \code{last_policy_gap}),
#'   \code{stage} (\code{"egm"} or \code{"stationary_dist"}: which solve
#'   failed), \code{message} (the diagnostic \code{\link{hank_het3_block}}
#'   would have \code{stop()}'d with under \code{strict = TRUE}),
#'   \code{iterations}, \code{last_value_gap}, \code{last_policy_gap},
#'   \code{dist_iterations} (\code{NA} when the EGM stage itself failed),
#'   \code{relax}, \code{maxit}, \code{tol}, \code{dist_maxit}, \code{dist_tol},
#'   \code{Lambda} (the forward operator, only when the failure was at the
#'   stationary-distribution stage), and -- when cheaply obtainable --
#'   \code{regime_counts}, a \code{\link{hank_euler3_residual}}-style
#'   floor/interior/ceiling count per asset at the last iterate. This last
#'   field is useful for distinguishing slow-but-healthy convergence (regime
#'   counts already settled, just more iterations needed) from a genuinely
#'   broken active set.
#' @aliases hank_het3_block_failed
#' @seealso \code{\link{hank_het3_block}}
#' @export
print.hank_het3_block_failed <- function(x, ...) {
  cat("<hank_het3_block_failed> stage: ", x$stage, "\n", sep = "")
  cat(x$message, "\n\n", sep = "")
  if (!is.null(x$regime_counts)) {
    cat("Regime counts (floor/interior/ceiling) at the last iterate:\n")
    print(x$regime_counts)
  }
  invisible(x)
}

#' Nonlinear perfect-foresight transition of a three-asset block
#'
#' Given aggregate input PATHS over horizon \code{T_h} (with terminal
#' conditions returning to steady state), computes the aggregate output paths
#' by a backward household solve carrying the \code{(Vd, Vf, Va)} triple,
#' followed by a forward distribution simulation. Three-asset counterpart of
#' \code{\link{hank_td2_nonlinear}}; this is the brute-force reference
#' \code{\link{hank_het3_jacobian_nd}} central-differences to validate the
#' fake-news Jacobian \code{\link{hank_het3_jacobian}}.
#'
#' @param block A \code{\link{hank_het3_block}}.
#' @param rd_path,rf_path,ra_path,w_path Numeric length-\code{T_h} input paths
#'   (levels): liquid return, foreign return, capital return, wage. A
#'   \code{NULL} path (default) holds the block's steady-state value constant.
#' @param px_path Numeric length-\code{T_h} path of the world price of foreign
#'   claims; \code{NULL} (default) holds \code{block$px}. This is where the
#'   valuation channel has its content: a PERMANENT \code{px} nearly cancels
#'   between purchase and payoff, whereas a date-\code{t} move revalues a stock
#'   bought at a different price, so the interesting object is the path, not
#'   the level. Every entry must be positive.
#' @param T_h Integer horizon; default the length of the longest supplied
#'   path (\code{1} if none is supplied).
#' @param keep_policies Logical (default \code{FALSE}). When \code{TRUE}, also
#'   return the per-period \code{d_pol}/\code{f_pol}/\code{a_pol}/\code{c_pol}
#'   lists and the per-period sparse \code{Lambda} list (same no-transpose
#'   convention as \code{block$Lambda}).
#' @param threads Worker threads for the backward policy step, the dominant
#'   cost of this transition. \code{NULL} (the default) resolves via
#'   \code{getOption("dynhr.hank3_threads")} and then a machine-derived
#'   default, exactly as \code{\link{hank_egm3_solve}} does; \code{1} forces
#'   the serial path. Bit-identical at every thread count.
#' @param pi_input_paths Optional named list of length-\code{T_h} LEVEL paths
#'   for (a subset of) the block's transition-probability inputs
#'   (\code{names(block$Pi_inputs)}, e.g. \code{list(f = f_path, s = s_path)});
#'   missing inputs stay at their steady-state values. Requires a block built
#'   with \code{Pi_fn}/\code{Pi_inputs}. \code{Pi_t = Pi_fn(inputs_t)} is the
#'   transition BETWEEN \code{t} and \code{t+1}: it enters the date-\code{t}
#'   backward expectation AND the date-\code{t} forward push, the same timing
#'   convention as the one-asset \code{\link{hank_td_nonlinear}}. Default
#'   \code{NULL} (constant steady-state \code{Pi}; existing calls are
#'   byte-identical).
#' @param Tr_path Optional length-\code{T_h} LEVEL path of the lump-sum
#'   transfer; \code{NULL} (default) holds the block's steady-state \code{Tr}.
#'   The block's own \code{Tr_incidence} weight distributes it across income
#'   states every period -- the path scales the aggregate, not the rule.
#' @param D0 Optional initial (beginning-of-period-1) distribution, a
#'   length-\code{n_e*n_d*n_f*n_a} nonnegative vector summing to 1 in the
#'   package's three-asset cell order (see
#'   \code{\link{hank_forward_operator3}}). Default \code{NULL} starts from
#'   the steady-state \code{block$D}; supply it for state-dependence
#'   experiments.
#'
#' @return A list with numeric length-\code{T_h} paths \code{D}, \code{F},
#'   \code{A}, \code{C} (aggregate liquid, foreign, capital holdings and
#'   consumption), \code{CHI} and \code{PHI} (aggregate capital and foreign
#'   adjustment resources), and \code{Dpath} (\code{n_cell x T_h}): column
#'   \code{t} is the beginning-of-period-\code{t} distribution over which the
#'   date-\code{t} aggregates are taken. When \code{keep_policies = TRUE},
#'   also \code{d_pol}, \code{f_pol}, \code{a_pol}, \code{c_pol},
#'   \code{chi_pol}, \code{phi_pol} (each a length-\code{T_h} list of
#'   \code{n_e x n_d x n_f x n_a} arrays) and \code{Lambda} (length-\code{T_h}
#'   list of sparse forward operators).
#' @seealso \code{\link{hank_td2_nonlinear}} (two-asset),
#'   \code{\link{hank_td_nonlinear}} (one-asset), \code{\link{hank_het3_block}},
#'   \code{\link{hank_het3_jacobian}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' td <- hank_td3_nonlinear(blk, T_h = 2)
#' td$D
#' @export
hank_td3_nonlinear <- function(block,rd_path=NULL,rf_path=NULL,ra_path=NULL,w_path=NULL,px_path=NULL,T_h=NULL,keep_policies=FALSE,D0=NULL,threads=NULL,pi_input_paths=NULL,Tr_path=NULL) {
  if(!inherits(block,"hank_het3_block"))stop("hank_td3_nonlinear: block must be hank_het3_block")
  if(is.null(T_h))T_h<-max(1L,length(rd_path),length(rf_path),length(ra_path),length(w_path),length(px_path),length(Tr_path),
    if(is.null(pi_input_paths))0L else vapply(pi_input_paths,length,integer(1)))
  rd_path<-if(is.null(rd_path))rep(block$rd,T_h)else rd_path;rf_path<-if(is.null(rf_path))rep(block$rf,T_h)else rf_path;ra_path<-if(is.null(ra_path))rep(block$ra,T_h)else ra_path;w_path<-if(is.null(w_path))rep(block$w,T_h)else w_path
  ## Blocks predating the transfer carry no Tr field; .hank_block_tr reads
  ## that as 0, and Tr*omega then adds an exact zero to income.
  Tr_path<-if(is.null(Tr_path))rep(.hank_block_tr(block),T_h)else Tr_path
  omega<-.hank_block_omega(block)
  ## Blocks built before px existed carry no field; 1 is the pre-A4 kernel.
  px_path<-if(is.null(px_path))rep(if(is.null(block$px))1 else block$px,T_h)else px_path
  if(any(vapply(list(rd_path,rf_path,ra_path,w_path,px_path,Tr_path),length,integer(1))!=T_h))stop("hank_td3_nonlinear: paths must have T_h entries")
  ## px is a PRICE and divides nothing downstream, but a non-positive entry
  ## would silently invert the foreign FOC rather than fail; the kernel rejects
  ## it per solve, and this catches the whole path in one place.
  if(!is.numeric(px_path)||any(!is.finite(px_path))||any(px_path<=0))stop("hank_td3_nonlinear: px_path entries must be finite and positive")
  ## Per-period transition matrices; NULL = steady-state Pi everywhere, which
  ## is the zero-allocation path every pre-existing call takes.
  Pi_path<-.hank_pi_path(block,pi_input_paths,T_h);Pi_at<-function(t)if(is.null(Pi_path))block$Pi else Pi_path[[t]]
  dp<-fp<-ap<-cp<-chp<-php<-vector("list",T_h);Vd<-block$Vd;Vf<-block$Vf;Va<-block$Va
  for(t in T_h:1L){s<-.hank_egm3_step(block$d_grid,block$f_grid,block$a_grid,w_path[t]*block$e+Tr_path[t]*omega,Pi_at(t),rd_path[t],rf_path[t],ra_path[t],block$beta,block$eis,block$chi0,block$chi1,block$chi2,block$phi0,block$phi1,block$phi2,Vd,Vf,Va,px_path[t],threads);dp[[t]]<-s$d;fp[[t]]<-s$f;ap[[t]]<-s$a;cp[[t]]<-s$c;chp[[t]]<-s$chi;php[[t]]<-s$phi;Vd<-s$Vd;Vf<-s$Vf;Va<-s$Va}
  flat<-function(z)as.vector(aperm(z,c(4,3,2,1)));D<-if(is.null(D0))block$D else D0;Dpath<-matrix(0,length(D),T_h);DA<-FA<-AA<-CA<-CHI<-PHI<-numeric(T_h);Lams<-if(keep_policies)vector("list",T_h)else NULL
  for(t in seq_len(T_h)){Dpath[,t]<-D;DA[t]<-sum(D*flat(dp[[t]]));FA[t]<-sum(D*flat(fp[[t]]));AA[t]<-sum(D*flat(ap[[t]]));CA[t]<-sum(D*flat(cp[[t]]));CHI[t]<-sum(D*flat(chp[[t]]));PHI[t]<-sum(D*flat(php[[t]]));L<-hank_forward_operator3(dp[[t]],fp[[t]],ap[[t]],block$d_grid,block$f_grid,block$a_grid,Pi_at(t));if(keep_policies)Lams[[t]]<-L;D<-as.numeric(Matrix::t(L)%*%D)}
  out<-list(D=DA,F=FA,A=AA,C=CA,CHI=CHI,PHI=PHI,Dpath=Dpath);if(keep_policies){out$d_pol<-dp;out$f_pol<-fp;out$a_pol<-ap;out$c_pol<-cp;out$chi_pol<-chp;out$phi_pol<-php;out$Lambda<-Lams};out
}


#' Wrap a three-asset household as a sequence-space DAG block
#'
#' Registers a \code{\link{hank_het3_block}} as a node of a general
#' sequence-space model (\code{\link{hank_model}}), the three-asset sibling of
#' \code{\link{hank_het_block_spec}} (one-asset) and
#' \code{\link{hank_het2_block_spec}} (two-asset). Once wrapped, the block's
#' fake-news Jacobian is accumulated along the DAG by the chain rule and the
#' equilibrium is solved as \eqn{dU = -H_U^{-1} H_Z dZ}; the same model object
#' then feeds \code{\link{hank_state_space}} and
#' \code{\link{hank_kalman_loglik}}.
#'
#' \strong{What this gives you, and what remains yours.} The package owns the
#' household block, its Jacobian, the DAG accumulation, the GE solve, and the
#' likelihood machinery downstream. It does NOT own -- and cannot infer -- the
#' MODEL: which aggregate blocks exist, which paths are unknowns, which
#' equations are the zero targets, what drives the system, and the steady state
#' those are taken around. Those are the modelling decisions, and they are
#' arguments to \code{\link{hank_model}}. In particular there is no built-in
#' three-asset closure the way \code{\link{hank_twoasset_model}} supplies a
#' demonstration two-asset one: a three-asset economy has to say what the
#' foreign asset IS, and that is a modelling choice, not a default.
#'
#' \strong{Transition-probability inputs} (\code{f}, \code{s} from
#' \code{\link{hank_employment_income}}, say) are opt-in here even when the
#' block carries them, unlike \code{\link{hank_het3_jacobian}} where they
#' default on. The reason is the DAG: every declared input must be produced by
#' some upstream block or be a source, so declaring \code{"f"} commits you to
#' supplying it. Name them in \code{inputs} once an upstream matching block
#' produces them.
#'
#' @param name Character block name.
#' @param block A \code{\link{hank_het3_block}} solved at steady state.
#' @param inputs Character: aggregate inputs. Default \code{c("rd", "rf",
#'   "ra", "w")} -- the always-present prices. Add \code{"px"} to wire the
#'   foreign-valuation channel, and any of \code{names(block$Pi_inputs)} to
#'   wire endogenous transition probabilities. Validated here, so a bad wiring
#'   fails at spec time rather than inside the Jacobian dispatch.
#' @param outputs Character: a subset of \code{c("D", "F", "A", "C", "CHI",
#'   "PHI")} (aggregate liquid, foreign and capital holdings, consumption, and
#'   the capital/foreign adjustment resources). \code{CHI} and \code{PHI} are
#'   real goods absorbed by rebalancing, so a resource constraint that nets
#'   them out must request them.
#'
#' @return An object of class \code{hank_block} (kind \code{"het3"}).
#' @seealso \code{\link{hank_het_block_spec}} (one-asset),
#'   \code{\link{hank_het2_block_spec}} (two-asset), \code{\link{hank_model}},
#'   \code{\link{hank_het3_jacobian}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' spec <- hank_het3_block_spec("hh", blk, outputs = c("D", "C"))
#' spec$kind
#' @export
hank_het3_block_spec <- function(name, block,
                                 inputs = c("rd", "rf", "ra", "w"),
                                 outputs = c("D", "F", "A", "C")) {
  ## Reject the lower tiers explicitly. Validation below is on input NAMES,
  ## and a het/het2 block would fail it -- but with an error about an
  ## unsupported input rather than about the block, which sends the reader to
  ## the wrong place. Same reasoning as hank_het_block_spec's het2 guard.
  if (inherits(block, "hank_het_block"))
    stop("hank_het3_block_spec(): this is a ONE-asset block ",
         "(hank_het_block); use hank_het_block_spec(), whose inputs are ",
         "('r', 'w') and outputs ('A', 'C').")
  if (inherits(block, "hank_het2_block"))
    stop("hank_het3_block_spec(): this is a TWO-asset block ",
         "(hank_het2_block); use hank_het2_block_spec(), whose inputs are ",
         "('rb', 'ra', 'w') and outputs ('B', 'A', 'C').")
  if (!inherits(block, "hank_het3_block"))
    stop("hank_het3_block_spec(): 'block' must be a hank_het3_block.")
  .hank3_check_inputs(block, inputs, "hank_het3_block_spec")
  bad <- setdiff(outputs, .hank3_jac_outputs)
  if (length(bad))
    stop("hank_het3_block_spec(): unsupported output(s) ",
         paste0("'", bad, "'", collapse = ", "), "; this block supports ",
         paste0("'", .hank3_jac_outputs, "'", collapse = ", "), ".")
  structure(list(name = name, kind = "het3", inputs = inputs,
                 outputs = outputs, block = block),
            class = "hank_block")
}
