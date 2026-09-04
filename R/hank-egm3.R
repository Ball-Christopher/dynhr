## Experimental R-only three-asset EGM reference (Stage 4).
## The foreign and capital FOCs are solved jointly; this is intentionally
## small-grid/interior-only until it passes the discrete-oracle gates.

.hank3_interp2 <- function(Z, xg, yg, x, y) {
  ix <- max(1L,min(findInterval(x,xg),length(xg)-1L)); iy <- max(1L,min(findInterval(y,yg),length(yg)-1L))
  wx <- (x-xg[ix])/(xg[ix+1]-xg[ix]); wy <- (y-yg[iy])/(yg[iy+1]-yg[iy])
  (1-wx)*(1-wy)*Z[ix,iy]+wx*(1-wy)*Z[ix+1,iy]+(1-wx)*wy*Z[ix,iy+1]+wx*wy*Z[ix+1,iy+1]
}


.hank3_active_foc <- function(foc, start, lo, hi, tol=1e-7) {
  valid <- function(z) { g<-foc(z); int<-z>lo+1e-6&z<hi-1e-6; all((int&abs(g)<1e-5)|(!int&(z<=lo+1e-6&g<=1e-5|z>=hi-1e-6&g>=-1e-5))) }
  cand <- list(); add <- function(z) { if(all(is.finite(z))) cand[[length(cand)+1L]] <<- pmin(hi,pmax(lo,z)) }
  raw <- try(nleqslv::nleqslv(start,foc,method="Broyden",global="dbldog",control=list(ftol=tol,maxit=80))$x,silent=TRUE)
  if(!inherits(raw,"try-error")) {
    raw <- pmin(hi,pmax(lo,raw))
    # Strict convexity makes a valid joint-KKT point unique. Returning it now
    # avoids eight redundant edge/corner evaluations in the common interior
    # case; the full active-set enumeration remains the boundary fallback.
    if(valid(raw))return(raw)
    add(raw)
  }
  # One active costly asset: minimise the remaining scalar Euler residual.
  for(ff in c(lo[1],hi[1])) { o<-stats::optimize(function(aa)foc(c(ff,aa))[2]^2,interval=c(lo[2],hi[2]),tol=tol);add(c(ff,o$minimum)) }
  for(aa in c(lo[2],hi[2])) { o<-stats::optimize(function(ff)foc(c(ff,aa))[1]^2,interval=c(lo[1],hi[1]),tol=tol);add(c(o$minimum,aa)) }
  for(ff in c(lo[1],hi[1]))for(aa in c(lo[2],hi[2]))add(c(ff,aa))
  keep<-Filter(valid,cand);if(!length(keep))stop("hank_egm3_solve: no active-set candidate satisfies joint KKT conditions")
  keep[[which.min(vapply(keep,function(z)sum(abs(foc(z))),numeric(1)))]]
}

# One exact backward update for a finite-horizon transition.  Keeping this as
# a wrapper around the tested stationary kernel prevents the active-set and
# liquid-floor logic from diverging between steady state and transition paths.
## One backward step via hank_egm3_solve's own machinery (tol below any
## reachable gap, maxit = 1, relax = 1 => exactly one undamped update).
## `backend` matters here: on the SINGLETON-f reduction this is the only step
## available, and pinning it to "R" both slowed it down and capped it at 250
## states. It defaults to the package backend so the reduction gets the
## compiled two-asset kernel, and "R" stays reachable for the reference.
.hank_egm3_step_r <- function(d_grid,f_grid,a_grid,y,Pi,rd,rf,ra,beta,eis,
                              chi0,chi1,chi2,phi0,phi1,phi2,Vd,Vf,Va,px=1,
                              backend=getOption("dynhr.hank3_backend","cpp")) {
  hank_egm3_solve(d_grid,f_grid,a_grid,y,Pi,rd,rf,ra,beta,eis,
    chi0,chi1,chi2,phi0,phi1,phi2,tol=1e-300,maxit=1L,relax=1,
    Vd_init=Vd,Vf_init=Vf,Va_init=Va,px=px,backend=backend)
}

## Shape validation shared by .hank_egm3_step() and the compiled kernel
## (src/hank_egm3.cpp::egm3_check_shapes). Names the offending argument and
## the shape it must have; the two implementations are kept deliberately
## parallel so a caller gets the same diagnosis whichever branch runs.
.hank_egm3_check_step_shapes <- function(d_grid, f_grid, a_grid, y, Pi,
                                         Vd, Vf, Va) {
  dm <- dim(Vd)
  if (is.null(dm) || length(dm) != 4L)
    stop(".hank_egm3_step: `Vd` must be a four-dimensional array ",
         "(n_e x n_d x n_f x n_a); dim(Vd) is ",
         if (is.null(dm)) "NULL" else paste(dm, collapse = " x "), ".",
         call. = FALSE)
  ne <- length(y); nd <- length(d_grid)
  nf <- length(f_grid); na <- length(a_grid)
  want <- c(ne, nd, nf, na)
  if (!identical(as.integer(dm), as.integer(want)))
    stop(".hank_egm3_step: dim(`Vd`) must be ",
         "length(y) x length(d_grid) x length(f_grid) x length(a_grid) = ",
         paste(want, collapse = " x "), "; got ",
         paste(dm, collapse = " x "), ".", call. = FALSE)
  for (nm in c("Vf", "Va")) {
    x <- if (nm == "Vf") Vf else Va
    if (is.null(dim(x)) || !identical(as.integer(dim(x)), as.integer(want)))
      stop(".hank_egm3_step: `", nm, "` must have the same dim as `Vd` (",
           paste(want, collapse = " x "), "); got ",
           if (is.null(dim(x))) paste0("a length-", length(x),
                                       " object with no dim")
           else paste(dim(x), collapse = " x "), ".", call. = FALSE)
  }
  if (!is.matrix(Pi) || nrow(Pi) != ne || ncol(Pi) != ne)
    stop(".hank_egm3_step: `Pi` must be an n_e x n_e matrix with ",
         "n_e = length(y) = ", ne, "; got ",
         if (is.matrix(Pi)) paste(dim(Pi), collapse = " x ")
         else paste0("a non-matrix of length ", length(Pi)), ".",
         call. = FALSE)
  invisible(TRUE)
}

## `threads` MUST be forwarded here. hank_egm3_step_cpp defaults it to 1, so
## every caller that omitted it ran the policy root single-threaded -- and this
## step is 87.8% of a fake-news Jacobian (measured by Rprof at the medium rung:
## 25.44 s of 28.98 s). A5's threading therefore reached hank_egm3_solve's own
## fixed-point loop and nothing else, which is exactly the part of a production
## run that is NOT the bottleneck at paper scale. Bit-identical at every thread
## count, as test-hank-egm3-threads.R asserts, so this is purely a throughput
## knob on both callers (the fake-news sweep and hank_td3_nonlinear).
.hank_egm3_step <- function(d_grid,f_grid,a_grid,y,Pi,rd,rf,ra,beta,eis,
                            chi0,chi1,chi2,phi0,phi1,phi2,Vd,Vf,Va,px=1,
                            threads=1L) {
  ## Shape gate (B4). The compiled kernel reads dim(Vd) and then indexes every
  ## grid with raw pointer arithmetic, so a dropped dim attribute or a grid
  ## whose length disagrees with the array used to be an out-of-bounds READ,
  ## not an error. The kernel now checks too -- these R checks exist so the
  ## message is identical on BOTH branches below, including the singleton-f
  ## reduction, which never reaches the kernel at all.
  .hank_egm3_check_step_shapes(d_grid, f_grid, a_grid, y, Pi, Vd, Vf, Va)
  if (length(f_grid) == 1L)
    ## The singleton-f reduction has no three-asset compiled step; it routes to
    ## hank_egm3_solve's reduction path, which delegates to the two-asset
    ## COMPILED kernel. That is what lifts the old 250-state ceiling on
    ## Jacobians and nonlinear transitions over the reduction control.
    return(.hank_egm3_step_r(d_grid,f_grid,a_grid,y,Pi,rd,rf,ra,beta,eis,
      chi0,chi1,chi2,phi0,phi1,phi2,Vd,Vf,Va,px))
  hank_egm3_step_cpp(Vd,Vf,Va,d_grid,f_grid,a_grid,y,Pi,rd,rf,ra,beta,eis,
    chi0,chi1,chi2,phi0,phi1,phi2,px,hank_resolve_threads(threads))
}

.hank3_lift2 <- function(x) {
  dx <- dim(x)
  array(x, c(dx[1L], dx[2L], 1L, dx[3L]))
}

.hank3_reduction_solve <- function(d_grid, f_grid, a_grid, y, Pi, rd, rf, ra,
                                   beta, eis, chi0, chi1, chi2,
                                   phi0, phi1, phi2, tol, maxit,
                                   Vd_init, Vf_init, Va_init,
                                   backend = getOption("dynhr.hank3_backend", "cpp")) {
  if (length(f_grid) != 1L || f_grid[1L] != 0)
    stop(".hank3_reduction_solve: foreign grid must be the singleton zero")
  ne <- length(y); nd <- length(d_grid); na <- length(a_grid)
  drop_f <- function(x, nm) {
    if (is.null(x)) return(NULL)
    want <- c(ne, nd, 1L, na)
    if (!is.array(x) || !identical(dim(x), want) || any(!is.finite(x)))
      stop("hank_egm3_solve: ", nm,
           " must be a finite e x d x 1 x a array on the reduction path")
    array(x, c(ne, nd, na))
  }
  t0 <- proc.time()[["elapsed"]]
  hh <- hank_egm2_solve(
    d_grid, a_grid, y = y, rb = rd, ra = ra, beta = beta, eis = eis,
    chi0 = chi0, chi1 = chi1, chi2 = chi2, Pi = Pi,
    tol = tol, maxit = maxit,
    Vb_init = drop_f(Vd_init, "Vd_init"),
    ## The reduction IS the two-asset problem, so it should reach the two-asset
    ## COMPILED kernel like any other two-asset solve. This used to be pinned to
    ## "R", which left the ladder's (14,72,1,80) reduction control running an
    ## interpreted kernel and -- because .hank_egm3_step routes through here --
    ## made a Jacobian over it impossible above 250 states. Honouring the
    ## caller's backend keeps the reduction a genuine drop-in for the general
    ## path at BOTH backends, which is how test-hank-egm3.R now gates it.
    Va_init = drop_f(Va_init, "Va_init"), backend = backend
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  invisible(drop_f(Vf_init, "Vf_init"))
  Dp <- .hank3_lift2(hh$b)
  Ap <- .hank3_lift2(hh$a)
  Cp <- .hank3_lift2(hh$c)
  Vd <- .hank3_lift2(hh$Vb)
  Va <- .hank3_lift2(hh$Va)
  Vf <- (1 + rf) * Cp^(-1 / eis)
  structure(
    list(
      d_grid = d_grid, f_grid = f_grid, a_grid = a_grid,
      y = y, Pi = Pi, rd = rd, rf = rf, ra = ra,
      beta = beta, eis = eis,
      chi0 = chi0, chi1 = chi1, chi2 = chi2,
      phi0 = phi0, phi1 = phi1, phi2 = phi2,
      px = 1,
      d = Dp, f = array(0, dim(Dp)), a = Ap, c = Cp,
      ## Capital adjustment comes straight from the two-asset solver, which has
      ## always reported it. Foreign adjustment is IDENTICALLY zero here rather
      ## than merely small: f' = f = 0 throughout, and Phi(0, 0) = 0 for every
      ## admissible phi2 > 1. Returning the array (not NULL) keeps the
      ## reduction a drop-in for the general path in the resource accounting.
      chi = .hank3_lift2(hh$chi), phi = array(0, dim(Dp)),
      Vd = Vd, Vf = Vf, Va = Va,
      iterations = hh$iterations, converged = hh$converged,
      last_value_gap = NA_real_, last_policy_gap = NA_real_,
      state_count = ne * nd * na, reduction = "hank_egm2_solve",
      ## What actually ran: hank_egm2_solve at the backend the caller asked for.
      backend = backend, elapsed = elapsed
    ),
    class = "hank_egm3_solve"
  )
}

## Smallest damping the auto-retry will fall back to. Below this the iteration
## is so slow (iterations scale as 1/relax) that failing loudly is more useful
## than crawling: at relax = 1/16 a fixture needing ~900 undamped iterations
## would need ~15,000.
.relax_floor <- function() getOption("dynhr.hank3_relax_floor", 0.0625)

#' Resolve a worker-thread count for the compiled HANK kernels
#'
#' Shared by the TWO- and THREE-asset compiled kernels (the two-asset one
#' gained real worker threading in 0.9.0.0032), which is why the name carries
#' no asset count.
#'
#' Resolution order: an explicit \code{threads} argument, then
#' \code{getOption("dynhr.hank_threads")}, then the legacy
#' \code{getOption("dynhr.hank3_threads")} (still honoured -- downstream code
#' sets it), then a machine-derived default.
#'
#' The default leaves headroom rather than claiming the machine: one core free
#' on a small box, two on anything larger (measured 2026-07-27, the EGM policy
#' root is 99.0\% of a block build, so the useful ceiling is set by hardware,
#' not by the algorithm). Note that efficiency cores count toward
#' \code{detectCores()} but do far less of this dense floating-point work than
#' performance cores, so the realised speed-up saturates below the returned
#' count -- treat it as a budget, not a prediction.
#'
#' Under \code{R CMD check} (\code{_R_CHECK_LIMIT_CORES_} set) this caps at 2,
#' per CRAN policy, whatever the machine or the option says.
#'
#' @param threads Explicit worker count, or \code{NULL} to resolve.
#' @param ncores Core count to derive the default from; defaults to
#'   \code{parallel::detectCores()}. Exposed for testing the formula.
#' @return A positive integer worker count.
#' @keywords internal
hank_resolve_threads <- function(threads = NULL, ncores = NULL) {
  source <- "argument"
  ## TWO option names, deliberately. `dynhr.hank_threads` is the current one --
  ## this resolver is shared by the two- AND three-asset kernels since
  ## 0.9.0.0032, so a name carrying a "3" is a misnomer. But
  ## `dynhr.hank3_threads` is USER-FACING and already set by downstream code
  ## (the NZ HANK paper's Stage 4 steady-state driver sets it explicitly), so
  ## it is still honoured rather than silently ignored -- silently dropping a
  ## thread-count option does not fail, it just runs at the wrong width, which
  ## is exactly the invisible-slowdown class of bug this whole area exists to
  ## prevent. New name wins when both are set.
  if (is.null(threads)) {
    threads <- getOption("dynhr.hank_threads", NULL)
    source <- "option(dynhr.hank_threads)"
  }
  if (is.null(threads)) {
    threads <- getOption("dynhr.hank3_threads", NULL)
    source <- "option(dynhr.hank3_threads) [legacy name]"
  }
  if (is.null(threads)) {
    if (is.null(ncores)) ncores <- parallel::detectCores(logical = TRUE)
    if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
        ncores < 1) ncores <- 1L
    ncores  <- as.integer(ncores)
    threads <- if (ncores <= 4L) ncores - 1L else ncores - 2L
    source  <- paste0("machine default (detectCores=", ncores, ")")
  }
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads))
    stop("hank_egm3_solve: `threads` must be a single finite number.")
  threads <- max(1L, as.integer(threads))
  ## CRAN forbids more than two cores in examples/tests/vignettes.
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) && threads > 2L) {
    threads <- 2L
    source  <- paste0(source, ", CLAMPED to 2 by _R_CHECK_LIMIT_CORES_")
  }
  .hank_report_threads(threads, source)
  threads
}


## Opt-in visibility for the resolved thread count, off by default.
##
## WHY IT REPORTS ON CHANGE, NOT PER CALL. This resolver runs once per EGM
## solve, and an ND Jacobian battery performs thousands of them -- a message
## on every call is the "per-CALL warning inside a closure a sampler hits 20k
## times" anti-pattern this codebase has already been bitten by: it floods the
## log, gets switched off, and then the diagnostic does not exist when it is
## needed. Reporting only when the resolved (count, source) pair DIFFERS from
## the last one reported gives exactly one line per regime change, which is
## the event worth seeing: a long batch that silently drops from 16 threads to
## 2 because _R_CHECK_LIMIT_CORES_ is set announces itself once, loudly.
##
## Enable with options(dynhr.hank_report_threads = TRUE).
.hank_thread_report <- new.env(parent = emptyenv())
.hank_thread_report$last <- NULL

.hank_report_threads <- function(threads, source) {
  if (!isTRUE(getOption("dynhr.hank_report_threads", FALSE))) return(invisible(NULL))
  key <- paste0(threads, "|", source)
  if (identical(.hank_thread_report$last, key)) return(invisible(NULL))
  .hank_thread_report$last <- key
  message("dynhr: HANK compiled kernel using ", threads,
          " thread", if (threads == 1L) "" else "s", " [", source, "]")
  invisible(NULL)
}

#' Reset the thread-count reporter (testing)
#'
#' Clears the memo used by the \code{dynhr.hank_report_threads} reporter so the
#' next resolution reports again.
#'
#' @return Invisibly, the previously-memoized report key (or \code{NULL}).
#' @keywords internal
.hank_reset_thread_report <- function() {
  old <- .hank_thread_report$last
  .hank_thread_report$last <- NULL
  invisible(old)
}


#' Experimental three-asset endogenous-grid solver
#'
#' Three-asset counterpart of \code{\link{hank_egm2_solve}}: a household holds
#' domestic liquid claims \code{d} (return \code{rd}, floor \code{d_grid[1]}),
#' gross FOREIGN assets \code{f} (return \code{rf}, convex portfolio cost
#' \eqn{\Phi(f', f)}) and domestic CAPITAL \code{a} (return \code{ra}, convex
#' adjustment cost \eqn{\Psi(a', a)}; see \code{\link{.hank_psi}}, whose
#' functional form both frictions reuse with the \code{phi*}/\code{chi*}
#' parameter sets respectively). The liquid Euler equation is solved by
#' inversion, as in the one- and two-asset solvers, while the two costly-asset
#' first-order conditions are solved JOINTLY at every liquid gridpoint by a
#' safeguarded Broyden root (\code{.hank3_active_foc}), because moving
#' either costly asset changes the marginal value of the other.
#'
#' This is a Stage-4 EXPERIMENTAL reference: on the R backend it is capped to
#' a small state space, and it has no dedicated borrowing-corner branch beyond
#' the single liquid-floor case handled explicitly below \code{min(de)}. It
#' exists to gate the discrete \code{\link{hank_egm3_prototype}} oracle before
#' a production implementation. When \code{f_grid} is the SINGLETON \code{0}
#' it instead delegates exactly to \code{\link{hank_egm2_solve}} (the foreign
#' asset held at zero throughout), which doubles as a correctness reduction
#' test and an escape hatch from the small-grid cap.
#'
#' @param d_grid Numeric: increasing DOMESTIC LIQUID grid, length >= 3 (MAY
#'   include negative values -- \code{d_grid[1]} is the borrowing/liquid
#'   floor).
#' @param f_grid Numeric: increasing FOREIGN asset grid, length >= 3, or the
#'   SINGLETON \code{0} to route to the exact two-asset reduction (see
#'   Details).
#' @param a_grid Numeric: increasing domestic CAPITAL grid, length >= 3.
#' @param y Numeric length-\code{n_e}: labour income by income state.
#' @param Pi Numeric \code{n_e x n_e}: row-stochastic income transition matrix.
#' @param rd,rf,ra Numeric: returns on liquid, foreign and capital holdings.
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param chi0,chi1,chi2 Capital adjustment-cost parameters, passed to
#'   \code{\link{.hank_psi}} for \eqn{\Psi(a', a)}: \code{chi0 > 0}
#'   (denominator shift), \code{chi1 >= 0} (scale), \code{chi2 > 1}
#'   (curvature).
#' @param phi0,phi1,phi2 Foreign-portfolio adjustment-cost parameters, passed
#'   to \code{\link{.hank_psi}} for \eqn{\Phi(f', f)} with the same functional
#'   form and constraints as \code{chi0}/\code{chi1}/\code{chi2}.
#' @param tol Convergence tolerance on \eqn{\max(|dV_d|, |dV_f|, |dV_a|)}
#'   across a backward iteration.
#' @param maxit Maximum number of backward iterations.
#' @param relax Damping weight in \eqn{(0, 1]} on the marginal-value update
#'   (\code{V <- (1 - relax) * V + relax * V_new}); \code{1} is undamped.
#'   Guards the joint active-set iteration against overshoot near a kink.
#'
#'   \strong{Which direction actually helps a non-convergence.} Measured on a
#'   14-income-state, \code{(n_d,n_f,n_a)=(8,4,8)} fixture (3,584 states):
#'   \code{relax=1} converges in 927 iterations, \code{relax=.75} in 1,236,
#'   \code{relax=.5} (the package default elsewhere, e.g.
#'   \code{\link{hank_het3_block}}) in 1,854, and \code{relax=.25}, \code{.1}
#'   and \code{.05} all FAIL inside a 2,000-iteration cap. Iterations scale as
#'   \code{1/relax}: the damped update's contraction modulus is
#'   \code{1 - relax*(1-rho)} with \code{rho ~ beta}, so damping can only ever
#'   SLOW a monotone contraction -- it helps against oscillation, and on this
#'   class of fixture the oscillation is transient (an early rise in the value
#'   gap as the active set shuffles, well before the tail collapse). A caller
#'   reaching for \code{relax} to fix a non-convergence should therefore raise
#'   it TOWARD \code{1}, not lower it: lowering \code{relax} below ~0.5 can
#'   prevent convergence entirely on this class of fixture rather than help.
#'   Warm-starting from \code{Vd_init}/\code{Vf_init}/\code{Va_init} at a
#'   nearby already-converged solve is worth a further ~3.5x and is what
#'   actually traverses a genuine continuation frontier where even
#'   \code{relax=1} does not converge inside the iteration cap.
#' @param Vd_init,Vf_init,Va_init Optional \code{n_e x n_d x n_f x n_a} initial
#'   marginal values. \code{NULL} (default) starts every one from the
#'   positive, finite marginal utility of the stay-put budget.
#' @param px World price of foreign claims (default \code{1}, the pre-A4
#'   kernel exactly). Foreign claims are a QUANTITY traded at \code{px}, so it
#'   multiplies the foreign position on BOTH sides of the budget --
#'   \eqn{c + d' + p_x f' + a' + \Psi + \Phi = y + (1+r_d) d + p_x (1+r_f) f +
#'   (1+r_a) a} -- and correspondingly appears in the foreign FOC as
#'   \eqn{p_x + \Phi_1} and in the envelope as \eqn{p_x (1+r_f) - \Phi_2}.
#'   Making it revalue only the predetermined stock was rejected: it would be
#'   exactly collinear with \code{rf} (both multiplying \eqn{f} alone), and a
#'   price that never touches the purchase side cannot clear a market for new
#'   purchases, which is what the Stage-4 contract asks of it. The adjustment
#'   technology is a domestic-goods resource cost and is NOT scaled by
#'   \code{px}. A PERMANENT \code{px} is close to neutral in steady state
#'   (\eqn{p_x(1+r_f)/p_x = 1+r_f}); the channel's content is dynamic, which is
#'   why \code{\link{hank_td3_nonlinear}} takes it as a path. Rejected on the
#'   singleton-\code{f_grid} reduction, where the foreign position is pinned at
#'   zero and there is nothing to price.
#' @param backend Character: \code{"cpp"} (default, via
#'   \code{getOption("dynhr.hank3_backend", "cpp")}) for the compiled
#'   active-set/fused solve, or \code{"R"} for the capped (<= 250 states)
#'   correctness reference. Both backends share validation, so they reject
#'   identical inputs. The singleton-\code{f_grid} reduction runs neither
#'   backend directly -- it calls \code{\link{hank_egm2_solve}}, which honours
#'   its own default backend.
#' @param auto_relax Logical (default \code{FALSE}). The endogenous-grid
#'   safeguard fires when an update overshoots into a region where the implied
#'   liquid grid is non-monotone. That is a property of the STEP, not of the
#'   starting values, so the compiled loop now returns the last good marginal
#'   values and this argument decides what happens next: \code{TRUE} halves
#'   \code{relax} and continues from that iterate (emitting a
#'   \code{message()} each time, down to
#'   \code{getOption("dynhr.hank3_relax_floor", 0.0625)}); \code{FALSE}
#'   (the default) reports the failure immediately.
#'
#'   \strong{Default \code{FALSE} on measured evidence: damping DELAYS this
#'   failure rather than curing it.} On a cross-rung continuation from
#'   \code{(14,24,8,24)} to \code{(14,72,4,80)}, the safeguard fired at
#'   iteration 197 at \code{relax = 1}, at 689 at \code{0.5}, and at 1368 at
#'   \code{0.25} -- roughly in proportion to \code{1/relax}, i.e. the same
#'   trajectory walked more slowly into the same region. At \code{0.125} it
#'   survived 2,000 iterations without converging (gap 0.399). The whole
#'   ladder cost 17 minutes and still failed, which inside a calibration loop
#'   is worse than failing in 196 iterations with a message that says why. Turn
#'   it on deliberately when you want the ladder tried; do not rely on it to
#'   rescue a grid that cannot be solved.
#' @param threads Worker threads for the compiled backend's policy root, which
#'   is 99\% of a block build. \code{NULL} (the default) resolves via
#'   \code{getOption("dynhr.hank3_threads")} and then a machine-derived
#'   default (see \code{hank_resolve_threads}); \code{1} forces the serial
#'   path. Output is BIT-IDENTICAL at every thread count -- each state's
#'   policies come from the same operations in the same order however the work
#'   is divided -- so this is purely a throughput knob and never a source of
#'   numerical difference. Ignored by \code{backend = "R"} and by the
#'   singleton-\code{f_grid} reduction.
#'
#' @return An object of class \code{hank_egm3_solve}: a list with the echoed
#'   calibration (\code{d_grid}, \code{f_grid}, \code{a_grid}, \code{y},
#'   \code{Pi}, \code{rd}, \code{rf}, \code{ra}, \code{beta}, \code{eis},
#'   \code{chi0}, \code{chi1}, \code{chi2}, \code{phi0}, \code{phi1},
#'   \code{phi2}, \code{px}) plus:
#'   \describe{
#'     \item{\code{d}, \code{f}, \code{a}, \code{c}}{Policies (each
#'       \code{n_e x n_d x n_f x n_a}): next-period liquid, foreign and
#'       capital holdings and current consumption.}
#'     \item{\code{chi}, \code{phi}}{Per-cell CAPITAL (\eqn{\Psi(a', a)}) and
#'       FOREIGN (\eqn{\Phi(f', f)}) adjustment costs at the returned policy,
#'       same shape. These are real resource costs, so they are the wedge
#'       between household financing and goods supply; they are reported
#'       separately (rather than summed) because the two frictions carry
#'       different parameter sets and enter the external accounts differently.
#'       Naming follows \code{\link{hank_egm2_solve}}'s \code{chi}, which is
#'       likewise the \code{chi*}-parameterised capital cost.}
#'     \item{\code{Vd}, \code{Vf}, \code{Va}}{Converged marginal values of
#'       liquid, foreign and capital wealth, same shape.}
#'     \item{\code{iterations}}{Number of backward iterations run.}
#'     \item{\code{converged}}{Logical: whether \code{tol} was met.}
#'     \item{\code{last_value_gap}}{Final \eqn{\max|dV|}; \code{NA} on the
#'       reduction path, which does not run this solver's own loop.}
#'     \item{\code{last_policy_gap}}{Final \eqn{\max|d(\mathrm{policy})|}
#'       between the last two iterations; \code{NA} on the reduction path.}
#'     \item{\code{state_count}}{\code{n_e * n_d * n_f * n_a}.}
#'     \item{\code{threads}}{Worker threads actually used (compiled backend
#'       only), for the run manifest.}
#'     \item{\code{backend}}{Character: the backend that actually ran
#'       (\code{"cpp"} or \code{"R"}). On the singleton-\code{f_grid}
#'       reduction this is always \code{"R"} -- what
#'       \code{\link{hank_egm2_solve}} was called with -- regardless of the
#'       \code{backend} argument supplied here, which that path ignores.}
#'     \item{\code{elapsed}}{Wall-clock seconds
#'       (\code{proc.time()[["elapsed"]]} difference) spent in the solve
#'       itself (the backward iteration / compiled call / delegated
#'       \code{hank_egm2_solve} call), not the whole function.}
#'     \item{\code{reduction}}{Present (\code{"hank_egm2_solve"}) only on the
#'       singleton-\code{f_grid} path: names the solver that actually produced
#'       the policies.}
#'   }
#' @seealso \code{\link{hank_egm2_solve}} (the two-asset solver this reduces
#'   to at \code{f_grid = 0}), \code{\link{hank_egm_solve}} (one-asset),
#'   \code{\link{hank_egm3_prototype}} (the discrete-choice oracle),
#'   \code{\link{hank_het3_block}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' hh <- hank_egm3_solve(dg, fg, ag, e, Pi, rd = .01, rf = .015, ra = .02,
#'                       beta = .97, eis = .5, chi1 = .2, phi1 = .1,
#'                       tol = 1e-5, maxit = 250)
#' hh$converged
#' @export
hank_egm3_solve <- function(d_grid, f_grid, a_grid, y, Pi, rd, rf, ra,
                            beta, eis, chi0=.25, chi1=6.5, chi2=2,
                            phi0=.25, phi1=.5, phi2=2,
                            tol=1e-7, maxit=500L, relax=.5,
                            Vd_init=NULL, Vf_init=NULL, Va_init=NULL,
                            px=1,
                            backend=getOption("dynhr.hank3_backend","cpp"),
                            threads=NULL, auto_relax=FALSE) {
  backend <- match.arg(backend,c("cpp","R"))
  nthreads <- hank_resolve_threads(threads)
  if(!is.numeric(px)||length(px)!=1L||!is.finite(px)||px<=0)
    stop("hank_egm3_solve: px (world price of foreign claims) must be a single finite positive number.")
  reduction <- length(f_grid) == 1L && isTRUE(f_grid[1L] == 0)
  if(length(d_grid)<3||length(a_grid)<3||(!reduction&&length(f_grid)<3))
    stop("hank_egm3_solve: d/a grids and every non-degenerate foreign grid need >=3 points")
  if(length(f_grid)==1L&&!reduction)
    stop("hank_egm3_solve: singleton foreign grid must be zero")
  if(!is.numeric(relax)||length(relax)!=1L||!is.finite(relax)||relax<=0||relax>1)stop("hank_egm3_solve: relax must be in (0,1]")
  ## The 250-state cap belongs to the THREE-ASSET R reference, which the
  ## reduction never runs -- it delegates wholesale to hank_egm2_solve. The
  ## check therefore has to come AFTER the reduction dispatch below, or a
  ## perfectly scalable two-asset problem gets rejected by a cap on a code
  ## path it does not touch. That is what blocked a Jacobian over the
  ## ladder's own (14,72,1,80) = 80,640-state reduction control.
  .hank_check_markov(Pi,length(y),caller="hank_egm3_solve")
  if(reduction&&px!=1)
    stop("hank_egm3_solve: the singleton-f_grid reduction holds the foreign ",
         "position at zero, so a world price px != 1 has nothing to price. ",
         "Supply a non-degenerate f_grid to use the valuation channel.")
  if(reduction)return(.hank3_reduction_solve(
    d_grid,f_grid,a_grid,y,Pi,rd,rf,ra,beta,eis,
    chi0,chi1,chi2,phi0,phi1,phi2,tol,maxit,Vd_init,Vf_init,Va_init,
    backend=backend))
  if(backend=="R"&&length(y)*length(d_grid)*length(f_grid)*length(a_grid)>250L)
    stop("hank_egm3_solve: R reference cap is 250 states")
  ne<-length(y);nd<-length(d_grid);nf<-length(f_grid);na<-length(a_grid)
  uc <- function(c)c^(-1/eis); invuc <- function(v)v^(-eis)
  # Positive, finite marginal-value initialisation from the stay-put budget.
  DD<-array(rep(rep(d_grid,each=ne),nf*na),c(ne,nd,nf,na));FF<-array(rep(rep(f_grid,each=ne*nd),na),c(ne,nd,nf,na));AA<-array(rep(a_grid,each=ne*nd*nf),c(ne,nd,nf,na));YY<-array(rep(y,nd*nf*na),c(ne,nd,nf,na))
  V0<-uc(pmax(YY+rd*DD+rf*FF+ra*AA,1e-4)); chkV<-function(V,nm){if(!is.null(V)&&(!is.array(V)||!identical(dim(V),dim(V0))||any(!is.finite(V))))stop("hank_egm3_solve: ",nm," must be a finite e x d x f x a array")};chkV(Vd_init,"Vd_init");chkV(Vf_init,"Vf_init");chkV(Va_init,"Va_init")
  Vd<-if(is.null(Vd_init))V0 else Vd_init;Vf<-if(is.null(Vf_init))V0 else Vf_init;Va<-if(is.null(Va_init))V0 else Va_init;ok<-FALSE;pol_gap<-Inf;Dold<-Fold<-Aold<-Cold<-NULL
  if(backend=="cpp") {
    t0<-proc.time()[["elapsed"]]
    ## A step-guard failure means THIS update overshot into an infeasible
    ## region -- not that the starting values are unusable. The compiled loop
    ## now returns the LAST GOOD marginal values with step_status = 1 instead
    ## of discarding them, so recovery is: damp and continue from there.
    ##
    ## Undamped iteration (relax = 1) is the fastest choice when it is stable
    ## and it is exactly what removes the safeguard when it is not. Measured on
    ## the paper's cross-rung continuation ((14,24,8,24) -> (14,72,4,80)):
    ## relax = 1 ran 196 clean iterations, the gap falling 30.5 -> 0.09 by
    ## iteration 12 and drifting back to 1.35 by 196, then overshot at 197.
    ## Every damped run cleared the same point. Halving on failure keeps the
    ## speed of relax = 1 where it works without the cliff where it does not.
    ## RESTART FROM THE ORIGINAL VALUES, not from the last good iterate.
    ## Damping changes how much of a step is absorbed, NOT the step itself, so
    ## resuming from the iterate that produced the bad step just recomputes it
    ## and fails again on the first try -- verified. What the measurements
    ## actually show is that a damped run traces a DIFFERENT trajectory that
    ## never reaches the offending region at all: from the same start,
    ## relax = 1 died at iteration 197 while 0.75 / 0.5 / 0.25 all walked past
    ## it. So the retry has to re-run from the beginning at lower damping.
    ## The completed iterations are lost; correctness first, and the message
    ## says how many so the cost is visible rather than mysterious.
    relax_used<-relax; retries<-0L; last_status<-0L; last_error<-""
    repeat {
      ans<-hank_egm3_solve_cpp(Vd,Vf,Va,d_grid,f_grid,a_grid,y,Pi,rd,rf,ra,
        beta,eis,chi0,chi1,chi2,phi0,phi1,phi2,tol,as.integer(maxit),
        relax_used,px,nthreads)
      last_status<-if(is.null(ans$step_status)) 0L else as.integer(ans$step_status)
      last_error<-if(is.null(ans$step_error)) "" else ans$step_error
      if(last_status!=1L||!isTRUE(auto_relax)||relax_used<=.relax_floor()) break
      relax_used<-max(.relax_floor(),relax_used/2)
      retries<-retries+1L
      message("hank_egm3_solve: the endogenous-grid safeguard fired after ",
              ans$iterations," iterations at relax = ",
              format(relax_used*2,digits=3),". This is an OVERSHOOT of the ",
              "update, not a bad starting point -- damping traces a different ",
              "trajectory, so restarting from the supplied values with ",
              "relax = ",format(relax_used,digits=3),
              ". Set auto_relax = FALSE to make this an error instead.")
    }
    ans$relax_used<-relax_used; ans$relax_retries<-retries
    if(last_status==1L)
      stop("hank_egm3_solve: ",last_error,
           " -- the endogenous-grid safeguard still fires at relax = ",
           format(relax_used,digits=3),
           " after ",retries," damped retr", if(retries==1L)"y" else "ies",
           ". This is an overshoot guard, ",
           "not a bad starting point: lower `relax` further, or coarsen the ",
           "grid. See ?hank_egm3_solve.",call.=FALSE)
    ans$step_status<-NULL; ans$step_error<-NULL
    elapsed<-proc.time()[["elapsed"]]-t0
    return(structure(c(list(d_grid=d_grid,f_grid=f_grid,a_grid=a_grid,y=y,
      Pi=Pi,rd=rd,rf=rf,ra=ra,beta=beta,eis=eis,chi0=chi0,chi1=chi1,
      chi2=chi2,phi0=phi0,phi1=phi1,phi2=phi2,px=px,
      state_count=ne*nd*nf*na,threads=nthreads,backend=backend,
      elapsed=elapsed),ans),
      class="hank_egm3_solve"))
  }
  t0_r<-proc.time()[["elapsed"]]
  for(it in seq_len(as.integer(maxit))) {
    Ed<-Ef<-Ea<-array(0,c(ne,nd,nf,na));for(e in seq_len(ne)){Ed[e,,,]<-apply(Vd,c(2,3,4),function(q)sum(Pi[e,]*q));Ef[e,,,]<-apply(Vf,c(2,3,4),function(q)sum(Pi[e,]*q));Ea[e,,,]<-apply(Va,c(2,3,4),function(q)sum(Pi[e,]*q))}
    Dp<-Fp<-Ap<-Cp<-array(NA_real_,c(ne,nd,nf,na))
    for(e in seq_len(ne))for(jf0 in seq_len(nf))for(ja0 in seq_len(na)) {
      de<-numeric(nd);fe<-ae<-ce<-numeric(nd)
      z_start<-c(f_grid[jf0],a_grid[ja0])
      for(jd in seq_len(nd)) {
        foc<-function(z){ff<-z[1];aa<-z[2];wd<-.hank3_interp2(Ed[e,jd,,],f_grid,a_grid,ff,aa);wf<-.hank3_interp2(Ef[e,jd,,],f_grid,a_grid,ff,aa);wa<-.hank3_interp2(Ea[e,jd,,],f_grid,a_grid,ff,aa);c(wf/wd-px-.hank_psi(ff,f_grid[jf0],rf,phi0,phi1,phi2)$Psi1,wa/wd-1-.hank_psi(aa,a_grid[ja0],ra,chi0,chi1,chi2)$Psi1)}
        z<-.hank3_active_foc(foc,z_start,c(min(f_grid),min(a_grid)),c(max(f_grid),max(a_grid)));z_start<-z
        wd<-.hank3_interp2(Ed[e,jd,,],f_grid,a_grid,z[1],z[2]);cc<-invuc(beta*wd);pf<-.hank_psi(z[1],f_grid[jf0],rf,phi0,phi1,phi2)$Psi;pa<-.hank_psi(z[2],a_grid[ja0],ra,chi0,chi1,chi2)$Psi
        de[jd]<-(cc+d_grid[jd]+px*z[1]+z[2]+pf+pa-y[e]-px*(1+rf)*f_grid[jf0]-(1+ra)*a_grid[ja0])/(1+rd);fe[jd]<-z[1];ae[jd]<-z[2];ce[jd]<-cc
      }
      if(any(diff(de)<=0))stop("hank_egm3_solve: non-monotone endogenous liquid grid; corner/iteration safeguard required")
      zz_start<-c(f_grid[jf0],a_grid[ja0])
      for(id in seq_len(nd)) {
        if(d_grid[id] < min(de)) {
          # Liquid constraint: d' is fixed at its floor, so uc(c)=beta*Wd is
          # an INEQUALITY.  The two costly-asset FOCs instead use the budget
          # marginal utility at this current d state.
          ffoc <- function(z) {
            pf<-.hank_psi(z[1],f_grid[jf0],rf,phi0,phi1,phi2);pa<-.hank_psi(z[2],a_grid[ja0],ra,chi0,chi1,chi2)
            cc<-y[e]+(1+rd)*d_grid[id]+px*(1+rf)*f_grid[jf0]+(1+ra)*a_grid[ja0]-d_grid[1]-px*z[1]-z[2]-pf$Psi-pa$Psi
            if(!is.finite(cc)||cc<=0)return(c(1e12,1e12))
            wf<-.hank3_interp2(Ef[e,1,,],f_grid,a_grid,z[1],z[2]);wa<-.hank3_interp2(Ea[e,1,,],f_grid,a_grid,z[1],z[2])
            c(beta*wf/uc(cc)-px-pf$Psi1,beta*wa/uc(cc)-1-pa$Psi1)
          }
          zz<-.hank3_active_foc(ffoc,zz_start,c(min(f_grid),min(a_grid)),c(max(f_grid),max(a_grid)));zz_start<-zz
          pf<-.hank_psi(zz[1],f_grid[jf0],rf,phi0,phi1,phi2)$Psi;pa<-.hank_psi(zz[2],a_grid[ja0],ra,chi0,chi1,chi2)$Psi
          Dp[e,id,jf0,ja0]<-d_grid[1];Fp[e,id,jf0,ja0]<-zz[1];Ap[e,id,jf0,ja0]<-zz[2];Cp[e,id,jf0,ja0]<-y[e]+(1+rd)*d_grid[id]+px*(1+rf)*f_grid[jf0]+(1+ra)*a_grid[ja0]-d_grid[1]-px*zz[1]-zz[2]-pf-pa
        } else {
          # c comes from the BUDGET at the interpolated portfolio, not from a
          # fourth independent interpolation of the endogenous-grid c: the
          # budget is nonlinear in (f',a') through Phi/Psi, so interpolating c
          # alongside them violates it off a source knot. Mirrors
          # hank_egm2_solve, and the compiled kernel does the same.
          dpi<-approx(de,d_grid,xout=d_grid[id],rule=2)$y
          fpi<-approx(de,fe,xout=d_grid[id],rule=2)$y
          api<-approx(de,ae,xout=d_grid[id],rule=2)$y
          pfi<-.hank_psi(fpi,f_grid[jf0],rf,phi0,phi1,phi2)$Psi
          pai<-.hank_psi(api,a_grid[ja0],ra,chi0,chi1,chi2)$Psi
          Dp[e,id,jf0,ja0]<-dpi;Fp[e,id,jf0,ja0]<-fpi;Ap[e,id,jf0,ja0]<-api
          Cp[e,id,jf0,ja0]<-y[e]+(1+rd)*d_grid[id]+px*(1+rf)*f_grid[jf0]+(1+ra)*a_grid[ja0]-dpi-px*fpi-api-pfi-pai
        }
      }
    }
    # Envelope: dc/df=(1+rf)-Psi_F,2 and dc/da=(1+ra)-Psi_A,2.
    # Psi2 is normally negative, so adding it would incorrectly *reduce* the
    # value of current costly wealth; this follows the validated two-asset
    # `(rho_a - Psi2) * uc` update exactly.
    psF<-.hank_psi(Fp,FF,rf,phi0,phi1,phi2);psA<-.hank_psi(Ap,AA,ra,chi0,chi1,chi2)
    # Same object that the budget above charged the household, evaluated AT the
    # final policy -- the adjustment technologies are real resource costs and
    # the block must be able to report them. Free here: the envelope already
    # needs Psi2 from these very calls.
    Phic<-psF$Psi;Chic<-psA$Psi
    Vdn<-uc(Cp)*(1+rd);Vfn<-uc(Cp)*(px*(1+rf)-psF$Psi2);Van<-uc(Cp)*((1+ra)-psA$Psi2)
    gap<-max(abs(Vdn-Vd),abs(Vfn-Vf),abs(Van-Va));if(!is.null(Dold))pol_gap<-max(abs(Dp-Dold),abs(Fp-Fold),abs(Ap-Aold),abs(Cp-Cold));Dold<-Dp;Fold<-Fp;Aold<-Ap;Cold<-Cp;Vd<-(1-relax)*Vd+relax*Vdn;Vf<-(1-relax)*Vf+relax*Vfn;Va<-(1-relax)*Va+relax*Van;if(gap<tol){ok<-TRUE;break}
  }
  elapsed_r<-proc.time()[["elapsed"]]-t0_r
  # px MUST be echoed on this path too, not only the compiled one. Downstream
  # consumers (hank_het3_block, and hank_egm3_diagnostics, which defaults a
  # missing px to 1 for pre-A4 objects) would otherwise silently diagnose an
  # R-backend px != 1 solve at px = 1 -- correct-looking at px = 1, wrong
  # exactly where the valuation channel has content.
  structure(list(d_grid=d_grid,f_grid=f_grid,a_grid=a_grid,y=y,Pi=Pi,rd=rd,rf=rf,ra=ra,beta=beta,eis=eis,chi0=chi0,chi1=chi1,chi2=chi2,phi0=phi0,phi1=phi1,phi2=phi2,px=px,d=Dp,f=Fp,a=Ap,c=Cp,chi=Chic,phi=Phic,Vd=Vd,Vf=Vf,Va=Va,iterations=it,converged=ok,last_value_gap=gap,last_policy_gap=pol_gap,state_count=ne*nd*nf*na,backend=backend,elapsed=elapsed_r),class="hank_egm3_solve")
}
