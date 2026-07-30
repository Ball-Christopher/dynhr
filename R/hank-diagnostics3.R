## R/hank-diagnostics3.R
## --------------------------------------------------------------------------
## Stage-4 gate: P-budget/Euler diagnostics for the three-asset HANK household
## (hank_egm3_solve / hank_het3_block). Three-asset counterpart of
## hank_euler2_residual (R/hank-egm2.R), reporting REPORTED RESIDUALS -- not
## just pass/fail assertions -- for a run manifest: the budget identity, all
## three Euler/FOC residuals (with complementarity handled at bounds), the
## active-regime classification, and (where a distribution is available or
## can be built) the boundary/interior mass and aggregate accounting.
## --------------------------------------------------------------------------


## Broadcast helpers on the (e, d, f, a)-shaped array, one per axis. Mirror
## the inline constructions in hank_egm3_solve() (R/hank-egm3.R, the DD/FF/AA/
## YY block) exactly, so the budget identity below is evaluated against the
## SAME broadcast convention the solver itself uses.
.hank3_bcast_e <- function(v, ne, nd, nf, na)
  array(rep(v, nd * nf * na), c(ne, nd, nf, na))
.hank3_bcast_d <- function(v, ne, nd, nf, na)
  array(rep(rep(v, each = ne), nf * na), c(ne, nd, nf, na))
.hank3_bcast_f <- function(v, ne, nd, nf, na)
  array(rep(rep(v, each = ne * nd), na), c(ne, nd, nf, na))
.hank3_bcast_a <- function(v, ne, nd, nf, na)
  array(rep(v, each = ne * nd * nf), c(ne, nd, nf, na))


## Interpolation coordinates against one increasing (or degenerate, length-1)
## axis, for the trilinear reconstruction below. Weights EXTRAPOLATE outside
## the grid (matching .hank2_bilin's convention in hank-egm2.R): an oracle
## that clamped to the edge instead would charge a spurious residual to
## exactly the cells where the policy legitimately overshoots the grid via
## linear interpolation upstream.
.hank3_axis_interp <- function(g, q) {
  n <- length(g)
  if (n == 1L) return(list(i = rep(1L, length(q)), p = rep(1, length(q))))
  i <- pmin(pmax(findInterval(q, g), 1L), n - 1L)
  p <- (g[i + 1L] - q) / (g[i + 1L] - g[i])
  list(i = i, p = p)
}


## Trilinear reconstruction of a (n_d x n_f x n_a) grid function at vectorized
## query points (dq, fq, aq). Degenerate (length-1) axes contribute weight 1
## at their sole index (the "hi" corner is clamped to the same index, and its
## coefficient (1 - p) is exactly 0), so this also handles the singleton
## foreign-grid reduction without a special case.
.hank3_interp3 <- function(Z, dg, fg, ag, dq, fq, aq) {
  nd <- length(dg); nf <- length(fg); na <- length(ag)
  di <- .hank3_axis_interp(dg, dq)
  fi <- .hank3_axis_interp(fg, fq)
  ai <- .hank3_axis_interp(ag, aq)
  dhi <- pmin(di$i + 1L, nd); fhi <- pmin(fi$i + 1L, nf); ahi <- pmin(ai$i + 1L, na)
  pd <- di$p; pf <- fi$p; pa <- ai$p
  pd * pf * pa * Z[cbind(di$i, fi$i, ai$i)] +
    (1 - pd) * pf * pa * Z[cbind(dhi, fi$i, ai$i)] +
    pd * (1 - pf) * pa * Z[cbind(di$i, fhi, ai$i)] +
    pd * pf * (1 - pa) * Z[cbind(di$i, fi$i, ahi)] +
    (1 - pd) * (1 - pf) * pa * Z[cbind(dhi, fhi, ai$i)] +
    (1 - pd) * pf * (1 - pa) * Z[cbind(dhi, fi$i, ahi)] +
    pd * (1 - pf) * (1 - pa) * Z[cbind(di$i, fhi, ahi)] +
    (1 - pd) * (1 - pf) * (1 - pa) * Z[cbind(dhi, fhi, ahi)]
}


## Regime classification of one policy array against its grid: -1 = at the
## lower bound (floor), +1 = at the upper bound (grid ceiling), 0 = interior,
## within `tol`. A degenerate (length-1) grid has no distinct bounds; every
## state is reported "0" (interior) there since the axis never truly moves
## (the singleton-foreign reduction, see hank_egm3_solve's Details).
.hank3_regime <- function(pol, grid, tol) {
  if (length(grid) == 1L) return(array(0L, dim(pol)))
  lo <- grid[1L]; hi <- grid[length(grid)]
  out <- array(0L, dim(pol))
  out[pol <= lo + tol] <- -1L
  out[pol >= hi - tol] <- 1L
  out
}


## Six boundary masses plus the fully-interior mass of a stationary
## distribution over the three-asset state space. ONE implementation, so the
## cell-order convention below is written down (and can go stale) in exactly
## one place; hank_euler3_residual() and hank_het3_manifest() both call this
## rather than each re-deriving it.
##
## @param D Numeric length-`ne * nd * nf * na` stationary distribution,
##   package cell order (see hank_forward_operator3(): e slowest, then d,
##   then f, with a fastest).
## @param ne,nd,nf,na Grid sizes.
## @return A named numeric vector `c(d_lo, d_hi, f_lo, f_hi, a_lo, a_hi,
##   interior)`; `interior` is `NA_real_` when any of the three non-degenerate
##   axes has fewer than 3 points (no distinct interior to sum).
.hank3_grid_mass <- function(D, ne, nd, nf, na) {
  ## D's cell order (hank_forward_operator3): e slowest, then d, then f,
  ## with a fastest -- i.e. array(D, c(na, nf, nd, ne)) reproduces it
  ## directly (first axis fastest-varying, matching hank_het3_block's own
  ## `flat <- function(z) as.vector(aperm(z, c(4,3,2,1)))` convention).
  ## The summation itself is the family-agnostic `.hank_grid_mass_axes()`
  ## (R/hank-manifest-common.R), shared with the one- and two-asset manifests;
  ## the cell-order convention above is the only three-asset part.
  .hank_grid_mass_axes(D, c(na, nf, nd, ne),
                       c("a", "f", "d", NA_character_))
}


#' First-order-condition and budget residuals of a converged three-asset policy
#'
#' Self-contained correctness/reporting oracle for \code{\link{hank_egm3_solve}}
#' and \code{\link{hank_het3_block}}: the three-asset counterpart of
#' \code{\link{hank_euler2_residual}}. Where that function gates the two-asset
#' solver with pass/fail assertions, this one is built to REPORT NUMBERS for a
#' run manifest -- the Stage-4 paper contract requires a P-budget/Euler gate
#' with residuals covering the three grid ceilings, each Euler equation, the
#' budget identity, and aggregate accounting (see the file header).
#'
#' \strong{Budget.} The per-state residual of
#' \deqn{c + d' + f' + a' + \Psi(a', a) + \Phi(f', f) -
#'       \left[y + (1+r_d) d + (1+r_f) f + (1+r_a) a\right]}
#' evaluated directly against the solved policy arrays (no reconstruction).
#' On a converged solve this is EXACT to floating point at every state
#' (measured 8.9e-16 on the 54-state reference, in both backends), because
#' \code{hank_egm3_solve} forms consumption FROM this identity once the
#' portfolio is fixed, on all three branches.
#'
#' That was not always so, and the history is why the loud warning below
#' exists. This diagnostic's first run found the unconstrained branch
## The primes below MUST stay \eqn, not \code: Rd's parser reads an apostrophe
## inside \code{} as OPENING a quoted string, which only closes on the next
## apostrophe. An ODD number of primes in one \code run (three, here) leaves it
## open and swallows the rest of the section, so checkRd reports every later
## tag as "invalid in a \code block" -- 40 lines away, naming innocent lines.
## \eqn is also what these are: math symbols, as \eqn{u'(c)} below already is.
#' interpolating \eqn{c} independently of \eqn{d'}, \eqn{f'}, \eqn{a'}
#' off the endogenous grid; since the budget is NONLINEAR in the portfolio
#' (through \eqn{\Psi}/\eqn{\Phi}), that fourth interpolation did not
#' preserve it between knots -- 1.98e-4 on the reference fixture, shrinking
#' only under refinement, i.e. a genuine discretization error rather than
#' solver noise, and a cap on any aggregate accounting built on it. The kernel
#' was brought onto \code{\link{hank_egm2_solve}}'s long-standing convention
#' (consumption from the budget after the assets are set) and the residual
#' fell to machine precision, grid-independently.
#'
#' So if the overall maximum exceeds the internal threshold,
#' \code{hank_euler3_residual} warns rather than silently reporting a number
#' that looks like machine precision but is not: it now means either a
#' regression to interpolated consumption, or an object built by an older
#' dynhr.
#'
#' \strong{Euler/FOC residuals.} Continuation marginal values are
#' RECONSTRUCTED FROM THE POLICIES (not the stored \code{Vd}/\code{Vf}/\code{Va}),
#' exactly as \code{\link{hank_euler2_residual}} does, for the same two
#' reasons: independence from the solver's own marginal-value bookkeeping, and
#' accuracy (interpolating the near-linear policies commits far less error
#' than interpolating the violently convex \eqn{u'(c)} composite). All three
#' conditions are written with a single sign convention -- \eqn{g_x} is the
#' marginal benefit of increasing \eqn{x'} net of its marginal cost, so
#' \eqn{g_x = 0} where the choice is interior, \eqn{g_x \le 0} is consistent
#' with sitting at a LOWER bound (the household would cut \eqn{x'} further if
#' it could), and \eqn{g_x \ge 0} is consistent with sitting at an UPPER bound:
#' \deqn{g_d = W_d(e, d', f', a') - u'(c)}
#' \deqn{g_f = W_f(e, d', f', a') - u'(c)\,(1 + \Phi_1(f', f))}
#' \deqn{g_a = W_a(e, d', f', a') - u'(c)\,(1 + \Psi_1(a', a))}
#' where \eqn{W_d = (1+r_d)\,E[u'(c')]}, \eqn{W_f = E[(1+r_f-\Phi_2(f'',f'))\,u'(c')]},
#' \eqn{W_a = E[(1+r_a-\Psi_2(a'',a'))\,u'(c')]}, the expectation is over
#' \code{Pi} at the chosen \eqn{(d',f',a')}, and \eqn{f''}/\eqn{a''} are the
#' interpolated NEXT foreign/capital choices (the envelope's second state).
#' The COMPLEMENTARITY-RESPECTING VIOLATION -- what should be approximately
#' zero everywhere -- is \eqn{|g_x|} where \eqn{x'} is interior, \eqn{\max(0,
#' g_x)} at a lower bound (a positive \eqn{g_x} there is the wrong sign: the
#' household would want to raise \eqn{x'} but is reported constrained) and
#' \eqn{\max(0, -g_x)} at an upper bound.
#'
#' \strong{Regime, grid mass, aggregates.} \code{regime} classifies every
#' state's CHOSEN policy against its own grid (floor / interior / ceiling);
#' \code{grid_mass} (when a stationary distribution is available or can be
#' built) reports the STATE distribution's mass at each of the six grid
#' boundaries plus the fully-interior mass; \code{aggregate} reconstructs
#' \code{D}/\code{F}/\code{A}/\code{C} from that distribution and, for a
#' \code{\link{hank_het3_block}}, the gap against the block's own
#' \code{D_agg}/\code{F_agg}/\code{A_agg}/\code{C}.
#'
#' @param hh A solved three-asset household from \code{\link{hank_egm3_solve}}
#'   or a \code{\link{hank_het3_block}}.
#' @param constraint_tol Cells within \code{constraint_tol} of a grid's lower
#'   or upper endpoint are classified as at that bound (both for the
#'   complementarity-respecting violation and for \code{regime}).
#'
#' @return A list:
#'   \describe{
#'     \item{\code{budget}}{List with \code{residual} (\code{n_e x n_d x n_f
#'       x n_a}) and \code{max_abs}.}
#'     \item{\code{euler_d}, \code{euler_f}, \code{euler_a}}{Each a list with
#'       \code{raw} (the signed \eqn{g_x} array), \code{violation} (the
#'       complementarity-respecting array described above, \code{NA} on a
#'       degenerate singleton grid) and \code{max_abs_violation}.}
#'     \item{\code{regime}}{List with \code{d}, \code{f}, \code{a}
#'       (\code{n_e x n_d x n_f x n_a} integer arrays, \code{-1}/\code{0}/
#'       \code{1} for floor/interior/ceiling) and \code{counts} (a
#'       \code{data.frame} of floor/interior/ceiling counts per asset,
#'       summing to \code{length(hh$c)} in each row).}
#'     \item{\code{grid_mass}}{\code{NULL} if no stationary distribution is
#'       available or can be built (i.e. \code{hh} lacks the fields
#'       \code{\link{hank_forward_operator3}} needs); otherwise a named
#'       numeric vector \code{d_lo}, \code{d_hi}, \code{f_lo}, \code{f_hi},
#'       \code{a_lo}, \code{a_hi}, \code{interior}, plus \code{dist_source}
#'       (\code{"hh$D"} or \code{"built"}) and \code{dist_converged}.
#'
#'       \strong{These masses OVERLAP ACROSS AXES and do not sum to 1.} Each is
#'       the marginal mass on one boundary of one axis, so a household sitting
#'       at both the liquid floor and the foreign floor is counted in
#'       \code{d_lo} AND \code{f_lo}. What holds is per-axis:
#'       \code{d_lo + d_hi <= 1}, and likewise for \code{f} and \code{a},
#'       with equality when no mass is interior on that axis.
#'       \code{interior} is the mass strictly inside on ALL THREE axes at
#'       once, so it is a lower bound on each axis's interior mass rather than
#'       a complement of the six. On the 54-state reference the seven sum to
#'       2.42. (An earlier version of this page claimed they sum to
#'       \code{<= 1}; the installed-package smoke test asserted that and
#'       failed, which is how it was found. Within one axis the boundaries are
#'       indeed disjoint whenever that grid has \code{>= 3} points -- the
#'       error was extending that across axes.)}
#'     \item{\code{aggregate}}{\code{NULL} under the same condition as
#'       \code{grid_mass}; otherwise a list with \code{D}, \code{F}, \code{A},
#'       \code{C} (reconstructed from the distribution) and, when \code{hh}
#'       carries its own \code{D_agg}/\code{F_agg}/\code{A_agg}/\code{C}
#'       (a \code{hank_het3_block}), \code{gap} (reconstructed minus
#'       reported, each \code{~0}).}
#'     \item{\code{max_abs}}{Compact named numeric vector for a run manifest:
#'       \code{budget}, \code{euler_d}, \code{euler_f}, \code{euler_a} (the
#'       last three are the complementarity-respecting violations).}
#'   }
#' @seealso \code{\link{hank_euler2_residual}} (two-asset),
#'   \code{\link{hank_egm3_solve}}, \code{\link{hank_het3_block}},
#'   \code{\link{hank_forward_operator3}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' hh <- hank_egm3_solve(dg, fg, ag, e, Pi, rd = .01, rf = .015, ra = .02,
#'                       beta = .97, eis = .5, chi1 = .2, phi1 = .1,
#'                       tol = 1e-5, maxit = 250)
#' diag3 <- hank_euler3_residual(hh)
#' diag3$max_abs
#' @export
hank_euler3_residual <- function(hh, constraint_tol = 1e-8) {
  need <- c("d", "f", "a", "c", "d_grid", "f_grid", "a_grid", "y", "Pi",
            "beta", "eis", "rd", "rf", "ra",
            "chi0", "chi1", "chi2", "phi0", "phi1", "phi2")
  if (!is.list(hh) || !all(need %in% names(hh)))
    stop("hank_euler3_residual: 'hh' must be a solved three-asset household ",
         "(from hank_egm3_solve) or a hank_het3_block.")
  if (!(is.numeric(constraint_tol) && length(constraint_tol) == 1L &&
        is.finite(constraint_tol) && constraint_tol >= 0))
    stop("hank_euler3_residual: 'constraint_tol' must be a finite ",
         "non-negative scalar.")

  d_grid <- hh$d_grid; f_grid <- hh$f_grid; a_grid <- hh$a_grid
  ne <- dim(hh$c)[1L]; nd <- dim(hh$c)[2L]; nf <- dim(hh$c)[3L]; na <- dim(hh$c)[4L]
  eis <- hh$eis; beta <- hh$beta
  uc_pow <- function(c) pmax(c, tiny_floor())^(-1 / eis)

  ## ---- 1. Budget identity -------------------------------------------------
  D_cur <- .hank3_bcast_d(d_grid, ne, nd, nf, na)
  F_cur <- .hank3_bcast_f(f_grid, ne, nd, nf, na)
  A_cur <- .hank3_bcast_a(a_grid, ne, nd, nf, na)
  Y_cur <- .hank3_bcast_e(hh$y, ne, nd, nf, na)
  Psi_a_now <- .hank_psi(hh$a, A_cur, hh$ra, hh$chi0, hh$chi1, hh$chi2)$Psi
  Psi_f_now <- .hank_psi(hh$f, F_cur, hh$rf, hh$phi0, hh$phi1, hh$phi2)$Psi
  ## World price of foreign claims (A4). Objects built before px existed do
  ## not carry the field; default to 1, which is the pre-px kernel exactly.
  ## px must appear in BOTH the budget and the f-FOC, or this oracle would
  ## report a spurious violation against a correct px != 1 solve -- and would
  ## look fine at px = 1, which is how such a gap survives.
  px <- if (is.null(hh$px)) 1 else hh$px

  budget_resid <- hh$c + hh$d + px * hh$f + hh$a + Psi_a_now + Psi_f_now -
    (Y_cur + (1 + hh$rd) * D_cur + px * (1 + hh$rf) * F_cur +
       (1 + hh$ra) * A_cur)
  budget_max_abs <- max(abs(budget_resid))
  ## Loud, not silent: a machine-precision-looking number that is not one
  ## means the policy arrays and the budget have drifted apart (file-header
  ## note above spells out the known unconstrained-branch cause).
  budget_loud_tol <- 1e-6
  if (budget_max_abs > budget_loud_tol)
    warning("hank_euler3_residual: budget residual max|.| = ",
            format(budget_max_abs, digits = 4), " exceeds ",
            format(budget_loud_tol, digits = 2), " (not machine precision). ",
            "The policy arrays and the budget identity have drifted apart on ",
            "at least one state. hank_egm3_solve() forms consumption FROM the ",
            "budget once the portfolio is fixed, so a converged solve should ",
            "be exact here; this diagnostic is what caught the kernel when it ",
            "did not (it once interpolated c independently of (d', f', a') ",
            "off the endogenous grid, which a nonlinear budget does not ",
            "preserve between knots -- 1.98e-4, fixed 2026-07-27). A non-zero ",
            "residual now means either a regression to that pattern, or an ",
            "object built by an older dynhr. See ?hank_euler3_residual.",
            call. = FALSE)

  ## ---- 2. Continuation marginal values, reconstructed from POLICIES -------
  Wd_at <- Wf_at <- Wa_at <- array(0, c(ne, nd, nf, na))
  for (e in seq_len(ne)) {
    dq <- as.numeric(hh$d[e, , , ])
    fq <- as.numeric(hh$f[e, , , ])
    aq <- as.numeric(hh$a[e, , , ])
    EVd <- EVf <- EVa <- numeric(nd * nf * na)
    for (ep in seq_len(ne)) {
      cp  <- .hank3_interp3(array(hh$c[ep, , , ], c(nd, nf, na)),
                            d_grid, f_grid, a_grid, dq, fq, aq)
      fpp <- .hank3_interp3(array(hh$f[ep, , , ], c(nd, nf, na)),
                            d_grid, f_grid, a_grid, dq, fq, aq)
      app <- .hank3_interp3(array(hh$a[ep, , , ], c(nd, nf, na)),
                            d_grid, f_grid, a_grid, dq, fq, aq)
      ucp <- uc_pow(cp)
      Psi2f_p <- .hank_psi(fpp, fq, hh$rf, hh$phi0, hh$phi1, hh$phi2)$Psi2
      Psi2a_p <- .hank_psi(app, aq, hh$ra, hh$chi0, hh$chi1, hh$chi2)$Psi2
      EVd <- EVd + hh$Pi[e, ep] * (1 + hh$rd) * ucp
      EVf <- EVf + hh$Pi[e, ep] * (px * (1 + hh$rf) - Psi2f_p) * ucp
      EVa <- EVa + hh$Pi[e, ep] * ((1 + hh$ra) - Psi2a_p) * ucp
    }
    Wd_at[e, , , ] <- beta * EVd
    Wf_at[e, , , ] <- beta * EVf
    Wa_at[e, , , ] <- beta * EVa
  }

  ## ---- 3. Euler/FOC residuals, single sign convention ---------------------
  uc_now <- uc_pow(hh$c)
  Psi1_f_now <- .hank_psi(hh$f, F_cur, hh$rf, hh$phi0, hh$phi1, hh$phi2)$Psi1
  Psi1_a_now <- .hank_psi(hh$a, A_cur, hh$ra, hh$chi0, hh$chi1, hh$chi2)$Psi1
  g_d <- Wd_at - uc_now
  g_f <- Wf_at - uc_now * (px + Psi1_f_now)
  g_a <- Wa_at - uc_now * (1 + Psi1_a_now)

  regime_d <- .hank3_regime(hh$d, d_grid, constraint_tol)
  regime_f <- .hank3_regime(hh$f, f_grid, constraint_tol)
  regime_a <- .hank3_regime(hh$a, a_grid, constraint_tol)

  violation_of <- function(g, regime, grid) {
    if (length(grid) == 1L) return(array(NA_real_, dim(g)))
    v <- abs(g)
    v[regime == -1L] <- pmax(0, g[regime == -1L])
    v[regime ==  1L] <- pmax(0, -g[regime ==  1L])
    v
  }
  viol_d <- violation_of(g_d, regime_d, d_grid)
  viol_f <- violation_of(g_f, regime_f, f_grid)
  viol_a <- violation_of(g_a, regime_a, a_grid)

  ## ---- 4. Regime counts ----------------------------------------------------
  count_row <- function(regime, nm) {
    data.frame(asset = nm,
               floor = sum(regime == -1L),
               interior = sum(regime == 0L),
               ceiling = sum(regime == 1L))
  }
  regime_counts <- do.call(rbind, list(count_row(regime_d, "d"),
                                        count_row(regime_f, "f"),
                                        count_row(regime_a, "a")))

  ## ---- 5. Stationary distribution (if carried or buildable) ---------------
  grid_mass <- NULL; aggregate <- NULL
  can_build <- all(is.finite(c(hh$d, hh$f, hh$a))) &&
    all(diff(d_grid) > 0) && all(diff(f_grid) > 0) && all(diff(a_grid) > 0)
  D <- NULL; dist_source <- NULL; dist_converged <- NA
  if (!is.null(hh$D) && is.numeric(hh$D) &&
      length(hh$D) == ne * nd * nf * na) {
    D <- hh$D; dist_source <- "hh$D"
    dist_converged <- if (!is.null(hh$dist_converged)) hh$dist_converged else NA
  } else if (can_build) {
    Q  <- hank_forward_operator3(hh$d, hh$f, hh$a, d_grid, f_grid, a_grid, hh$Pi)
    sd <- hank_stationary_dist(Q, tol = 1e-12, maxit = 100000L, backend = "R")
    D <- sd$d; dist_source <- "built"; dist_converged <- sd$converged
  }
  if (!is.null(D)) {
    grid_mass <- .hank3_grid_mass(D, ne, nd, nf, na)
    attr(grid_mass, "dist_source") <- dist_source
    attr(grid_mass, "dist_converged") <- dist_converged

    flat <- function(z) as.vector(aperm(z, c(4, 3, 2, 1)))
    D_recon <- sum(D * flat(hh$d))
    F_recon <- sum(D * flat(hh$f))
    A_recon <- sum(D * flat(hh$a))
    C_recon <- sum(D * flat(hh$c))
    aggregate <- list(D = D_recon, F = F_recon, A = A_recon, C = C_recon)
    if (!is.null(hh$D_agg) && !is.null(hh$F_agg) && !is.null(hh$A_agg) &&
        !is.null(hh$C))
      aggregate$gap <- c(D = D_recon - hh$D_agg, F = F_recon - hh$F_agg,
                         A = A_recon - hh$A_agg, C = C_recon - hh$C)
  }

  max_abs <- c(budget = budget_max_abs,
               euler_d = max(viol_d, na.rm = TRUE),
               euler_f = if (all(is.na(viol_f))) NA_real_ else max(viol_f, na.rm = TRUE),
               euler_a = if (all(is.na(viol_a))) NA_real_ else max(viol_a, na.rm = TRUE))

  list(
    budget = list(residual = budget_resid, max_abs = budget_max_abs),
    euler_d = list(raw = g_d, violation = viol_d,
                   max_abs_violation = max_abs[["euler_d"]]),
    euler_f = list(raw = g_f, violation = viol_f,
                   max_abs_violation = max_abs[["euler_f"]]),
    euler_a = list(raw = g_a, violation = viol_a,
                   max_abs_violation = max_abs[["euler_a"]]),
    regime = list(d = regime_d, f = regime_f, a = regime_a,
                  counts = regime_counts),
    grid_mass = grid_mass,
    aggregate = aggregate,
    max_abs = max_abs
  )
}
