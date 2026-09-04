## R/hank-sam-reiter.R
## --------------------------------------------------------------------------
## Finite-state Reiter emission with ENDOGENOUS transition probabilities
## (HANK+SAM, W7 of the RANK/HANK paper; design in
## .claude/orchestration/sam/DESIGN_endogenous_pi.md, E2(a)).
##
## The KS emission (R/hank-truncation.R) linearizes the household block only
## in (Va', r, w) and hard-codes the KS firm. Here the household's income
## TRANSITION MATRIX Pi_t = Pi_fn(f_t, s_t, ...) is itself an equilibrium
## object (job-finding/separation rates from a search-and-matching block), so
## the structural system needs two NEW derivative groups per transition input
## x (design note E2):
##   backward:  gVa_x, ga_x, gc_x = d(EGM step)/dx   via .hank_pi_perturb
##   forward:   D_x = d(Lambda(a_ss, Pi(x))' D_ss)/dx  (the DIRECT Pi entry in
##              the law of motion; Da covers only the policy entry)
##
## Timing (pinned in R/hank-het-block.R and verified against the fake-news
## sweep + the sequence-space GE oracle in test-hank-sam-reiter.R): Pi_t is
## the transition BETWEEN t and t+1 and is dated by x_t — it enters the
## date-t backward step (with Va_{t+1}) AND the date-t forward push. So all
## x-channel terms are date-t (B-side) entries, with anticipation of future
## x carried through the G_V recursion exactly as for prices.
##
## v1 closure (the toy SAM DAG proven at the sequence-space layer,
## test-hank-employment-ge.R): an exogenous AR(1) matching efficiency
##   dzm_{t+1} = rho_zm dzm_t + eps_{t+1}
## drives the transition inputs with fixed steady-state elasticities
## (input_derivs, e.g. c(f = eta_m * f_ss) from f_t = f_ss * Zm_t^eta_m), the
## wage is FIXED, and a bond market with fixed supply clears on r_t (a STATIC
## jump: its A-column is zero and Klein/QZ absorbs it as an infinite
## generalized eigenvalue, like the singular G_V rows).
##
## System: A x_{t+1} = B x_t,  x = [dD (n); dzm; dVa (n); dr],
## predetermined (dD, dzm), jumps (dVa, dr):
##   dD_{t+1} = Lambda' dD_t + Da da_t + sum_x D_x dx_t
##   dzm_{t+1} = rho_zm dzm_t (+ eps)
##   dVa_t    = G_V dVa_{t+1} + gVa_r dr_t + sum_x gVa_x dx_t
##   0        = va.dD_t + D_ss.da_t                       (asset clearing)
## with da_t = ga_V dVa_{t+1} + ga_r dr_t + sum_x ga_x dx_t and
## dx_t = input_derivs[x] * dzm_t.
## --------------------------------------------------------------------------


#' FD linearization of a transition-input het block for the SAM Reiter emission
#'
#' The expensive, \code{rho_zm}-independent half of
#' \code{\link{hank_sam_reiter_statespace}}: the KS derivative groups
#' (\code{G_V}, \code{ga_V}, \code{gc_V}, r-channel, \code{Da}) plus, for each
#' transition-probability input of the block, the backward
#' (\code{gVa_x, ga_x, gc_x}) and direct-forward (\code{D_x}) channels.
#'
#' @param block A \code{\link{hank_het_block}} built WITH
#'   \code{Pi_fn}/\code{Pi_inputs}, solved at the bond-clearing steady state
#'   (for the anchored coarse block see \code{\link{hank_coarse_anchored}}).
#' @param inputs Which transition inputs to differentiate (default: all of
#'   \code{names(block$Pi_inputs)}).
#' @param delta_fd Relative finite-difference step.
#' @return An object of class \code{hank_sam_reiter_lin}.
#' @seealso \code{\link{hank_sam_reiter_statespace}}
#' @export
hank_sam_reiter_linearize <- function(block, inputs = NULL, delta_fd = 1e-6) {
  ## NOTE the guard below is an AND, so it does NOT fire for a block that has a
  ## Pi_fn but is not a hank_het_block -- which is exactly a two-asset block
  ## built for the HANK+SAM route. Reject that explicitly first: two-asset
  ## Reiter is deliberately not implemented (briefs/19, section 7), and without
  ## this the call dies downstream on "argument is not a matrix".
  .hank_reject_het2(block, "hank_sam_reiter_linearize", use = NULL)
  .hank_reject_wedge(block, "hank_sam_reiter_linearize")
  if (!inherits(block, "hank_het_block") && is.null(block$Pi_fn))
    stop("hank_sam_reiter_linearize(): `block` must be a hank_het_block ",
         "built with Pi_fn/Pi_inputs.", call. = FALSE)
  if (is.null(block$Pi_fn) || is.null(block$Pi_inputs))
    stop("hank_sam_reiter_linearize(): the block has no Pi_fn/Pi_inputs -- ",
         "endogenous transition probabilities require a block built with ",
         "them (see hank_employment_income()).", call. = FALSE)
  if (is.null(inputs)) inputs <- names(block$Pi_inputs)
  bad <- setdiff(inputs, names(block$Pi_inputs))
  if (length(bad) > 0L)
    stop("hank_sam_reiter_linearize(): unknown transition input(s): ",
         paste(bad, collapse = ", "), ".", call. = FALSE)

  blk <- block
  n_e <- blk$n_e; n_a <- blk$n_a; n <- n_e * n_a
  Pi <- blk$Pi; r <- blk$r; w <- blk$w

  vc   <- .hank_mat_to_vec
  unvc <- function(v) matrix(v, n_e, byrow = TRUE)
  Va_ss <- blk$Va; a_ss <- blk$a; c_ss <- blk$c; D_ss <- blk$D
  a_grid <- blk$a_grid; amin <- a_grid[1L]

  ## ---- G_V, ga_V, gc_V (identical shape to the KS linearization) ----------
  G_V <- matrix(0, n, n); ga_V <- matrix(0, n, n); gc_V <- matrix(0, n, n)
  vVa <- vc(Va_ss)
  for (j in seq_len(n)) {
    h <- delta_fd * (1 + abs(vVa[j]))
    Vp <- vVa; Vp[j] <- Vp[j] + h
    Vm <- vVa; Vm[j] <- Vm[j] - h
    sp <- .hank_block_step(blk, unvc(Vp), r, w)
    sm <- .hank_block_step(blk, unvc(Vm), r, w)
    G_V[, j]  <- (vc(sp$Va) - vc(sm$Va)) / (2 * h)
    ga_V[, j] <- (vc(sp$a)  - vc(sm$a))  / (2 * h)
    gc_V[, j] <- (vc(sp$c)  - vc(sm$c))  / (2 * h)
  }

  ## ---- r channel (the market-clearing unknown; w is FIXED in v1) ----------
  hp <- delta_fd
  sp <- .hank_block_step(blk, Va_ss, r + hp, w)
  sm <- .hank_block_step(blk, Va_ss, r - hp, w)
  gVa_r <- (vc(sp$Va) - vc(sm$Va)) / (2 * hp)
  ga_r  <- (vc(sp$a)  - vc(sm$a))  / (2 * hp)
  gc_r  <- (vc(sp$c)  - vc(sm$c))  / (2 * hp)

  ## ---- transition-input channels: backward AND direct-forward -------------
  ## MATRIX-FREE (0.9.0.0026), same reasoning as R/hank-truncation.R's Da loop:
  ## the Da loop below calls push() up to 2n times and the D_x loop twice per
  ## transition input, each of which built an n x n sparse Lambda for one
  ## matvec. .hank_forward_push() contracts it instead; agreement is
  ## round-off-level (test-hank-forward-push.R). Note the Pi argument matters
  ## here in a way it does not in the KS linearization: the D_x channel pushes
  ## the SAME policy through PERTURBED transition matrices, and the push mixes
  ## over income with crossprod(Pi, .) after the asset scatter.
  push <- function(A_mat, Pi_use)
    .hank_forward_push(A_mat, a_grid, Pi_use, D_ss)
  gVa_x <- list(); ga_x <- list(); gc_x <- list(); D_x <- list()
  for (x in inputs) {
    hx <- delta_fd * (1 + abs(blk$Pi_inputs[[x]]))
    Pi_p <- .hank_pi_perturb(blk, x, +hx)
    Pi_m <- .hank_pi_perturb(blk, x, -hx)
    ## backward: the date-t EGM step with the perturbed Pi_t
    sp <- .hank_block_step(blk, Va_ss, r, w, Pi = Pi_p)
    sm <- .hank_block_step(blk, Va_ss, r, w, Pi = Pi_m)
    gVa_x[[x]] <- (vc(sp$Va) - vc(sm$Va)) / (2 * hx)
    ga_x[[x]]  <- (vc(sp$a)  - vc(sm$a))  / (2 * hx)
    gc_x[[x]]  <- (vc(sp$c)  - vc(sm$c))  / (2 * hx)
    ## forward: the DIRECT Pi entry in Lambda(a_ss, Pi(x))' D_ss
    D_x[[x]] <- (push(a_ss, Pi_p) - push(a_ss, Pi_m)) / (2 * hx)
  }

  ## ---- Da (policy entry in the forward push; up-step at constrained nodes) --
  p0 <- push(a_ss, Pi)
  va <- vc(a_ss)
  Da <- matrix(0, n, n)
  for (k in seq_len(n)) {
    hk <- delta_fd * (1 + abs(va[k]))
    ap <- va; ap[k] <- ap[k] + hk
    if (va[k] <= amin + 1e-12) {
      Da[, k] <- (push(unvc(ap), Pi) - p0) / hk
    } else {
      am <- va; am[k] <- am[k] - hk
      Da[, k] <- (push(unvc(ap), Pi) - push(unvc(am), Pi)) / (2 * hk)
    }
  }

  structure(
    list(G_V = G_V, ga_V = ga_V, gc_V = gc_V,
         gVa_r = gVa_r, ga_r = ga_r, gc_r = gc_r,
         gVa_x = gVa_x, ga_x = ga_x, gc_x = gc_x, D_x = D_x,
         Da = Da, va = va, inputs = inputs,
         n = n, n_e = n_e, n_a = n_a, delta_fd = delta_fd, block = blk),
    class = c("hank_sam_reiter_lin", "hank_block"))
}


#' Reiter state space of a bond-economy HANK with endogenous transition risk
#'
#' Assembles and Klein-solves the finite-state Reiter system for a household
#' block whose income transition matrix responds to an exogenous AR(1)
#' matching-efficiency state \code{zm} (the toy SAM closure: fixed wage, bond
#' market clears on \code{r}). Returns the same
#' \code{(T_mat, R_mat, Z_mat, sigma_z)} contract as
#' \code{\link{hank_reiter_statespace}}, so
#' \code{\link{hank_reiter_irf}} and \code{\link{hank_reiter_kalman_loglik}}
#' apply unchanged — this is the Kalman-filterable HANK+SAM emission.
#'
#' @param x A transition-input \code{\link{hank_het_block}} (see
#'   \code{hank_sam_reiter_linearize}) or a precomputed
#'   \code{hank_sam_reiter_lin} (pass the latter per posterior draw: only the
#'   Klein solve depends on \code{rho_zm}).
#' @param input_derivs Named numeric vector: steady-state derivative of each
#'   responding transition input w.r.t. \code{zm} (e.g.
#'   \code{c(f = eta_m * f_ss)} for \code{f_t = f_ss * Zm_t^eta_m} with
#'   \code{zm = log Zm}). Inputs not named are held fixed.
#' @param rho_zm AR(1) persistence of the matching-efficiency state.
#' @param sigma_zm Innovation standard deviation (stored for the KF).
#' @param drop_dist_coord Drop the redundant mass-conservation distribution
#'   coordinate (default \code{TRUE}; keep it on for filtering).
#' @param u_cells Optional integer vector of distribution cells whose mass is
#'   reported as the extra observation row \code{"u"} (e.g. the unemployment
#'   cells \code{(n_prod*n_a+1):(2*n_prod*n_a)} of an
#'   \code{\link{hank_employment_income}} block).
#' @param delta_fd FD step (ignored when \code{x} is a linearization).
#'
#' @return An object of classes \code{hank_sam_reiter_ss} and
#'   \code{hank_reiter_ss}: \code{T_mat}, \code{R_mat}, \code{Z_mat} (rows
#'   \code{A}, \code{C}, \code{r}, \code{zm}, plus each responding input and
#'   optionally \code{u}), \code{state_names}, \code{spectral_radius},
#'   \code{Qj} (jump loading, rows = (dVa, dr)), \code{rho_z = rho_zm},
#'   \code{sigma_z = sigma_zm}.
#' @export
hank_sam_reiter_statespace <- function(x, input_derivs, rho_zm,
                                       sigma_zm = 0.01,
                                       drop_dist_coord = TRUE,
                                       u_cells = NULL,
                                       delta_fd = 1e-6) {
  if (!is.numeric(input_derivs) || is.null(names(input_derivs)) ||
      any(!nzchar(names(input_derivs))))
    stop("hank_sam_reiter_statespace(): `input_derivs` must be a fully named ",
         "numeric vector, e.g. c(f = eta_m * f_ss).", call. = FALSE)
  lin <- if (inherits(x, "hank_sam_reiter_lin")) x
         else hank_sam_reiter_linearize(x, inputs = names(input_derivs),
                                        delta_fd = delta_fd)
  bad <- setdiff(names(input_derivs), lin$inputs)
  if (length(bad) > 0L)
    stop("hank_sam_reiter_statespace(): input_derivs name(s) not in the ",
         "linearization: ", paste(bad, collapse = ", "), ".", call. = FALSE)
  .hank_sam_reiter_assemble(lin, input_derivs = input_derivs,
                            rho_zm = rho_zm, sigma_zm = sigma_zm,
                            drop_dist_coord = drop_dist_coord,
                            u_cells = u_cells)
}


#' @keywords internal
.hank_sam_reiter_assemble <- function(lin, input_derivs, rho_zm, sigma_zm,
                                      drop_dist_coord, u_cells = NULL) {
  blk <- lin$block
  n <- lin$n
  G_V <- lin$G_V; ga_V <- lin$ga_V; gc_V <- lin$gc_V
  gVa_r <- lin$gVa_r; ga_r <- lin$ga_r; gc_r <- lin$gc_r
  Da <- lin$Da; va <- lin$va
  D_ss <- blk$D; c_ss <- blk$c
  vc <- .hank_mat_to_vec

  ## Collapse the responding transition inputs onto the single zm state:
  ## d(channel)/dzm = sum_x d(channel)/dx * input_derivs[x].
  xs <- names(input_derivs)
  gVa_z <- Reduce(`+`, lapply(xs, function(x) lin$gVa_x[[x]] * input_derivs[[x]]))
  ga_z  <- Reduce(`+`, lapply(xs, function(x) lin$ga_x[[x]]  * input_derivs[[x]]))
  gc_z  <- Reduce(`+`, lapply(xs, function(x) lin$gc_x[[x]]  * input_derivs[[x]]))
  D_z   <- Reduce(`+`, lapply(xs, function(x) lin$D_x[[x]]   * input_derivs[[x]]))

  ## ---- structural system A x_{t+1} = B x_t, x = [dD; dzm; dVa; dr] --------
  m  <- 2L * n + 2L
  iD <- 1:n; iz <- n + 1L; iV <- n + 1L + 1:n; ir <- 2L * n + 2L
  npred <- n + 1L
  A <- matrix(0, m, m); B <- matrix(0, m, m)
  ## distribution law of motion (rows 1..n):
  ## dD_{t+1} = Lambda' dD_t + Da (ga_V dVa_{t+1} + ga_r dr_t + ga_z dzm_t)
  ##            + D_z dzm_t
  A[iD, iD] <- diag(n)
  A[iD, iV] <- -Da %*% ga_V
  B[iD, iD] <- as.matrix(Matrix::t(blk$Lambda))
  B[iD, ir] <- Da %*% ga_r
  B[iD, iz] <- Da %*% ga_z + D_z
  ## exogenous matching efficiency (row n+1):
  A[iz, iz] <- 1
  B[iz, iz] <- rho_zm
  ## backward/Euler block (rows n+2 .. 2n+1):
  ## dVa_t = G_V dVa_{t+1} + gVa_r dr_t + gVa_z dzm_t
  A[iV, iV] <- G_V
  B[iV, iV] <- diag(n)
  B[iV, ir] <- -gVa_r
  B[iV, iz] <- -gVa_z
  ## asset-market clearing (row 2n+2, STATIC in dr -- its A-column is zero,
  ## Klein absorbs the infinite eigenvalue):
  ## 0 = va.dD_t + D_ss.(ga_V dVa_{t+1} + ga_r dr_t + ga_z dzm_t)
  A[ir, iV] <- as.numeric(D_ss %*% ga_V)
  B[ir, iD] <- -va
  B[ir, ir] <- -sum(D_ss * ga_r)
  B[ir, iz] <- -sum(D_ss * ga_z)

  ## ---- Klein (2000) via QZ, identical bookkeeping to the KS assemble ------
  qzd <- QZ::qz.dgges(B, A)
  mu  <- complex(real = qzd$ALPHAR, imaginary = qzd$ALPHAI) / qzd$BETA
  sel <- is.finite(Mod(mu)) & (Mod(mu) < 1 + 1e-6)
  n_s <- sum(sel)
  if (n_s != npred)
    stop(sprintf(paste0("hank_sam_reiter_statespace(): Blanchard-Kahn failure ",
                        "-- %d stable generalized eigenvalues, need %d ",
                        "(predetermined states)."), n_s, npred))
  o <- QZ::qz.dtgsen(qzd$S, qzd$T, qzd$Q, qzd$Z, sel)
  Z11 <- o$Z[1:npred, 1:npred]
  Z21 <- o$Z[npred + 1:(n + 1L), 1:npred]
  S11 <- o$S[1:npred, 1:npred]
  T11 <- o$T[1:npred, 1:npred]
  Qj   <- Z21 %*% solve(Z11)                       # (dVa_t, dr_t) = Qj s_t
  Tmat <- Z11 %*% solve(T11, S11) %*% solve(Z11)   # s_{t+1} = Tmat s_t
  Rvec <- c(rep(0, n), 1)                          # shock enters dzm

  ## ---- observation rows (functions of s_t = (dD_t, dzm_t)) ----------------
  Qj_V <- Qj[seq_len(n), , drop = FALSE]           # dVa_t loading
  Qj_r <- Qj[n + 1L, ]                             # dr_t loading
  selz <- c(rep(0, n), 1)
  selD <- rbind(diag(n), 0)
  ## da_t = ga_V dVa_{t+1} + ga_r dr_t + ga_z dzm_t; dVa_{t+1} = Qj_V Tmat s_t
  da_map <- ga_V %*% (Qj_V %*% Tmat) + outer(ga_r, Qj_r) + outer(ga_z, selz)
  dc_map <- gc_V %*% (Qj_V %*% Tmat) + outer(gc_r, Qj_r) + outer(gc_z, selz)
  ZA <- as.numeric(va %*% t(selD)) + as.numeric(D_ss %*% da_map)
  ZC <- as.numeric(vc(c_ss) %*% t(selD)) + as.numeric(D_ss %*% dc_map)
  Zmat <- rbind(A = ZA, C = ZC, r = Qj_r, zm = selz)
  for (x in names(input_derivs))
    Zmat <- rbind(Zmat,
                  matrix(input_derivs[[x]] * selz, nrow = 1L,
                         dimnames = list(x, NULL)))
  if (!is.null(u_cells)) {
    zu <- numeric(npred); zu[u_cells] <- 1
    Zmat <- rbind(Zmat, u = zu)
  }
  state_names <- c(paste0("dD", seq_len(n)), "dzm")

  ## ---- drop the redundant mass-conservation distribution coordinate -------
  if (drop_dist_coord) {
    keep <- c(seq_len(n - 1L), n + 1L)             # drop dD_n
    E <- matrix(0, npred, npred - 1L)              # embed: dD_n = -sum(others)
    E[keep, ] <- diag(npred - 1L)
    E[n, 1:(n - 1L)] <- -1
    Tmat <- (Tmat %*% E)[keep, , drop = FALSE]
    Rvec <- Rvec[keep]
    Zmat <- Zmat %*% E
    Qj   <- Qj %*% E
    state_names <- state_names[keep]
  }
  sr <- max(Mod(eigen(Tmat, only.values = TRUE)$values))
  if (drop_dist_coord && sr >= 1 - 1e-10)
    warning(sprintf(paste0("hank_sam_reiter_statespace(): spectral radius ",
                           "%.8f >= 1 after dropping the distribution ",
                           "coordinate -- state space is not stationary."), sr))

  structure(
    list(T_mat = Tmat, R_mat = matrix(Rvec, ncol = 1L), Z_mat = Zmat,
         state_names = state_names, obs_names = rownames(Zmat),
         spectral_radius = sr, n_stable = n_s,
         n_infinite = sum(abs(qzd$BETA) < 1e-12),
         gen_eig_mod = Mod(mu), Qj = Qj,
         rho_z = rho_zm, sigma_z = sigma_zm,
         input_derivs = input_derivs,
         calibration = list(r = blk$r, w = blk$w, beta = blk$beta,
                            eis = blk$eis, A = blk$A, C = blk$C)),
    class = c("hank_sam_reiter_ss", "hank_reiter_ss", "hank_block"))
}
