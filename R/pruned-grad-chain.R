## R/pruned-grad-chain.R
## --------------------------------------------------------------------------
## C1 (phase b): chain the order-2 solution-derivative layer
## (solution_derivatives_order2(), R/gradient-solution-deriv-order2.R) through
## the AFVRR order-2 augmented-state-space assembly (.order2_aug_system(),
## R/stochsimul-monolith.R:1331, mirrored by pruned_ss_loglik(),
## R/pruned-state-space.R:278) to get d(Tlin, ZZ, d_y, c_drift, QQ, HH, SS,
## mu0, Sxi0)/d(theta_j) for every theta_j solution_derivatives_order2 covers.
##
## This is a pure Kronecker product-rule pass: every augmented-system matrix
## is an explicit (multi-)linear function of hx, hu, hxx, hxu, huu, hss,
## Sigma_e (all of which have known d/dtheta_j from the order-1/2 solution-
## derivative layers), so no new implicit-differentiation solve is needed
## here EXCEPT for the two stationary Lyapunov fixed points (Sigma_x, already
## differentiated by solution_derivatives_order2()'s d_Sigma_x, and the
## augmented Sxi0), whose derivative is itself a Lyapunov solve with a
## perturbed RHS -- exactly the pattern already used by `.o2sd_dSigma_x()`
## (R/gradient-solution-deriv-order2.R:191) and `solve_lyapunov()`
## (R/stochsimul-monolith.R:474): if X = A X A' + B then
## dX = A dX A' + (dA X A' + A X dA' + dB), i.e.
## dX = solve_lyapunov(A, dA X A' + A X dA' + dB).
##
## Convention: gradients are dL/dX_ij treating every matrix entry as
## independent -- matching .pruned_kf_correlated_adjoint()'s convention
## exactly (R/pruned-kf-adjoint.R), so the two layers compose by a plain
## Frobenius inner product (sum(bar_X * dX)) with no extra symmetrization
## Jacobian.
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## d(4th-moment Gaussian tensor)/dtheta, mirroring .fourth_moment_gaussian()
## (R/cumulant-likelihood.R:677) index-for-index:
##   M4[ij,kl] = S[i,j]S[k,l] + S[i,k]S[j,l] + S[i,l]S[j,k]
##   with i = (ij-1)%/%n + 1, j = (ij-1)%%n + 1 (and likewise k,l from kl).
## Plain product rule per term.
## ---------------------------------------------------------------------------
.pgc_d_fourth_moment_gaussian <- function(Sigma_e, dSigma_e) {
  n_u  <- nrow(Sigma_e)
  n_u2 <- n_u * n_u
  dM4  <- matrix(0, n_u2, n_u2)
  for (ij in seq_len(n_u2)) {
    i <- ((ij - 1L) %/% n_u) + 1L
    j <- ((ij - 1L) %% n_u) + 1L
    for (kl in seq_len(n_u2)) {
      k <- ((kl - 1L) %/% n_u) + 1L
      l <- ((kl - 1L) %% n_u) + 1L
      dM4[ij, kl] <-
        dSigma_e[i, j] * Sigma_e[k, l] + Sigma_e[i, j] * dSigma_e[k, l] +
        dSigma_e[i, k] * Sigma_e[j, l] + Sigma_e[i, k] * dSigma_e[j, l] +
        dSigma_e[i, l] * Sigma_e[j, k] + Sigma_e[i, l] * dSigma_e[j, k]
    }
  }
  dM4
}

#' d(Cov(r_t))/dtheta at stationarity (a = 0, P = Sigma_x), mirroring
#' .order2_cov_r(0, Sigma_x, Sigma_e) block-for-block
#' (R/stochsimul-monolith.R:1304-1327).  a = 0 identically at every theta
#' (the forward code always evaluates Cr0 with a fixed zero mean), so every
#' block that is linear in "a" alone (Cr[j1,j2], Cr[j1,j3]) is identically
#' zero and stays zero under differentiation; only the Sigma_e/M(=Sigma_x)
#' blocks move.
#' @noRd
.pgc_d_cov_r0 <- function(Sigma_e, dSigma_e, Sigma_x, dSigma_x) {
  n_u <- nrow(Sigma_e); n_s <- nrow(Sigma_x)
  M  <- Sigma_x
  dM <- dSigma_x

  vecSe  <- as.numeric(Sigma_e)
  dvecSe <- as.numeric(dSigma_e)
  dC_ee_ee <- .pgc_d_fourth_moment_gaussian(Sigma_e, dSigma_e) -
    (outer(dvecSe, vecSe) + outer(vecSe, dvecSe))

  ## C_ex1_x1e[(i-1)*n_s+j, (k-1)*n_u+l] = Sigma_e[i,l] * M[j,k]
  dC_ex1_x1e <- matrix(0, n_u * n_s, n_s * n_u)
  for (i in seq_len(n_u)) for (j in seq_len(n_s))
    for (k in seq_len(n_s)) for (l in seq_len(n_u))
      dC_ex1_x1e[(i - 1L) * n_s + j, (k - 1L) * n_u + l] <-
        dSigma_e[i, l] * M[j, k] + Sigma_e[i, l] * dM[j, k]

  d1 <- n_u; d2 <- n_u * n_s; d3 <- n_s * n_u; d4 <- n_u * n_u
  D  <- d1 + d2 + d3 + d4
  dCr <- matrix(0, D, D)
  i1 <- seq_len(d1); i2 <- (d1 + 1L):(d1 + d2)
  i3 <- (d1 + d2 + 1L):(d1 + d2 + d3); i4 <- (d1 + d2 + d3 + 1L):D

  dCr[i1, i1] <- dSigma_e
  ## a = 0 identically => Cr[i1,i2], Cr[i1,i3] and their derivatives are 0.
  dCr[i2, i2] <- kronecker(dSigma_e, M) + kronecker(Sigma_e, dM)
  dCr[i2, i3] <- dC_ex1_x1e
  dCr[i3, i3] <- kronecker(dM, Sigma_e) + kronecker(M, dSigma_e)
  dCr[i4, i4] <- dC_ee_ee
  dCr[lower.tri(dCr)] <- t(dCr)[lower.tri(dCr)]
  dCr
}


#' d(augmented order-2 system matrices)/dtheta_j for ONE parameter
#'
#' Chains the order-1/2 solution derivatives through the AFVRR augmented-
#' state assembly (\code{.order2_aug_system()}) and the stationary noise-
#' covariance / initial-condition construction in \code{pruned_ss_loglik()},
#' producing d(Tlin, ZZ, d_y, c_drift, QQ, HH, SS, mu0, Sxi0)/dtheta_j -- the
#' exact nine inputs \code{.pruned_kf_correlated_adjoint()} differentiates
#' against.
#'
#' @param sys       Augmented system list from \code{.order2_aug_system()}
#'   (\code{pss$sys}); carries the BASE hx, hu, Sigma_e, Tlin, G, Dxi, Gv,
#'   ix1/ix2/ik, n_s/n_u/n_endo/d.
#' @param dr2       The order-2 \code{DecisionRules2} (for base ghu/ghxu/ghuu
#'   full matrices -- \code{.order2_aug_system()} only keeps the state-row
#'   slices \code{hx}/\code{hu} in \code{sys}, not the full ghu/ghxu/ghuu, so
#'   we re-read them from \code{dr2} directly, matching what
#'   \code{.order2_aug_system()} itself reads).
#' @param obs_idx   Integer index of observed variables in \code{pss$endo_names}.
#' @param dG,dH     d(ghx)/dtheta_j, d(ghu)/dtheta_j (n_endo x n_s / n_u;
#'   \code{first$derivs[[pnm]]$dG}/\code{dH} from \code{solution_derivatives()}).
#' @param d_ghxx,d_ghxu,d_ghuu,d_ghss  Order-2 solution derivatives for this
#'   parameter (\code{derivs[[pnm]]} fields from
#'   \code{solution_derivatives_order2()}; all full n_endo-row).
#' @param dSigma_e   d(Sigma_e)/dtheta_j (n_u x n_u).
#' @param Sigma_x    BASE Sigma_x (= \code{solve_lyapunov(hx, hu Sigma_e hu')}).
#' @param d_Sigma_x  d(Sigma_x)/dtheta_j, from
#'   \code{solution_derivatives_order2()}'s \code{derivs[[pnm]]$d_Sigma_x}
#'   (R/gradient-solution-deriv-order2.R's \code{.o2sd_dSigma_x()}, which
#'   since C5 includes the \code{hu*dSigma_e*hu'} shock-covariance channel,
#'   so no local recomputation is needed here any more).
#' @return A list with \code{dTlin, dZZ, d_dy, d_c_drift, dQQ, dHH, dSS, dmu0,
#'   dSxi0} matching the shapes of the nine \code{.pruned_kf_correlated_adjoint}
#'   inputs.
#' @noRd
.pruned_aug_system_deriv <- function(sys, dr2, obs_idx, dG, dH,
                                     d_ghxx, d_ghxu, d_ghuu, d_ghss,
                                     dSigma_e, Sigma_x, d_Sigma_x) {
  n_s   <- sys$n_s
  n_u   <- sys$n_u
  d_dim <- sys$d
  ix1 <- sys$ix1; ix2 <- sys$ix2; ik <- sys$ik
  sidx <- dr2$state_idx

  Sigma_e <- sys$Sigma_e
  hx <- sys$hx; hu <- sys$hu           # n_s x n_s, n_s x n_u  (BASE)

  ## Base FULL matrices (.order2_aug_system reads these straight off dr2).
  ghu_f  <- dr2$ghu                     # n_endo x n_u
  ghxu_f <- dr2$ghxu                    # n_endo x (n_s*n_u)
  ghuu_f <- dr2$ghuu                    # n_endo x n_u^2
  huu    <- ghuu_f[sidx, , drop = FALSE]  # n_s x n_u^2  (BASE)

  dhx <- dG[sidx, , drop = FALSE]
  dhu <- dH[sidx, , drop = FALSE]
  dhxx <- d_ghxx[sidx, , drop = FALSE]
  dhuu <- d_ghuu[sidx, , drop = FALSE]

  ## ---- d(Sigma_x): taken directly from solution_derivatives_order2()'s
  ## d_Sigma_x (see @param above) -- since C5 fixed .o2sd_dSigma_x() at the
  ## source to include the hu*dSigma_e*hu' shock-covariance channel, the
  ## local Lyapunov recomputation this file used to carry (C1's workaround
  ## for the then-missing channel) is now redundant and would double-count.

  ## ---- d(Tlin) ------------------------------------------------------------
  dTlin <- matrix(0, d_dim, d_dim)
  dTlin[ix1, ix1] <- dhx
  dTlin[ix2, ix2] <- dhx
  dTlin[ix2, ik]  <- 0.5 * dhxx
  dTlin[ik, ik]   <- kronecker(dhx, hx) + kronecker(hx, dhx)

  ## ---- d(cc), d(c_u) --------------------------------------------------------
  dcc <- numeric(d_dim); dcc[ix2] <- 0.5 * d_ghss[sidx]

  vecSe  <- as.numeric(Sigma_e)
  dvecSe <- as.numeric(dSigma_e)
  dc_u <- numeric(d_dim)
  dc_u[ix2] <- 0.5 * as.numeric(dhuu %*% vecSe + huu %*% dvecSe)
  dc_u[ik]  <- as.numeric(
    (kronecker(dhu, hu) + kronecker(hu, dhu)) %*% vecSe +
    kronecker(hu, hu) %*% dvecSe
  )
  d_c_drift <- dcc + dc_u

  ## ---- G, Dxi, Gv (base + derivative) --------------------------------------
  D1 <- n_u; D2 <- n_u * n_s; D3 <- n_s * n_u; D4 <- n_u * n_u
  Dr <- D1 + D2 + D3 + D4
  j1 <- seq_len(D1); j2 <- (D1 + 1L):(D1 + D2)
  j3 <- (D1 + D2 + 1L):(D1 + D2 + D3); j4 <- (D1 + D2 + D3 + 1L):Dr

  G_base <- matrix(0, d_dim, Dr)
  G_base[ix1, j1] <- hu
  G_base[ix2, j2] <- ghxu_f[sidx, , drop = FALSE]
  G_base[ix2, j4] <- 0.5 * huu
  G_base[ik, j2]  <- kronecker(hu, hx)
  G_base[ik, j3]  <- kronecker(hx, hu)
  G_base[ik, j4]  <- kronecker(hu, hu)

  dGmat <- matrix(0, d_dim, Dr)
  dGmat[ix1, j1] <- dhu
  dGmat[ix2, j2] <- d_ghxu[sidx, , drop = FALSE]
  dGmat[ix2, j4] <- 0.5 * dhuu
  dGmat[ik, j2]  <- kronecker(dhu, hx) + kronecker(hu, dhx)
  dGmat[ik, j3]  <- kronecker(dhx, hu) + kronecker(hx, dhu)
  dGmat[ik, j4]  <- kronecker(dhu, hu) + kronecker(hu, dhu)

  Gv_base <- matrix(0, sys$n_endo, Dr)
  Gv_base[, j1] <- ghu_f
  Gv_base[, j2] <- ghxu_f
  Gv_base[, j4] <- 0.5 * ghuu_f

  dGv <- matrix(0, sys$n_endo, Dr)
  dGv[, j1] <- dH
  dGv[, j2] <- d_ghxu
  dGv[, j4] <- 0.5 * d_ghuu

  dDxi <- cbind(dG, dG, 0.5 * d_ghxx)

  r_mean  <- numeric(Dr); r_mean[j4]  <- vecSe
  dr_mean <- numeric(Dr); dr_mean[j4] <- dvecSe
  dc_v <- as.numeric(dGv %*% r_mean + Gv_base %*% dr_mean)

  ## ---- Observation-side (subset to obs_idx): ZZ, d_y -----------------------
  ## d_y = ys[obs_vars] + 0.5*ghss[obs_idx] + c_v[obs_idx]; d(ys)/dtheta_j
  ## (the steady-state derivative) is added by the caller via d1$dys, since
  ## it is NOT part of the order-2 augmented-system chain proper.
  dZZ  <- dDxi[obs_idx, , drop = FALSE]
  d_dy_no_ys <- 0.5 * d_ghss[obs_idx] + dc_v[obs_idx]

  ## ---- Observation-subset Gv (n_obs x Dr): HH/SS use ONLY the observed
  ## rows of Gv (mirrors pruned_ss_loglik's `Gv <- Gv_f[obs_idx, ]`), while
  ## c_v above correctly used the FULL Gv_base/dGv (c_v is built over all
  ## endo, then subset to obs_idx afterward -- matching
  ## .order2_aug_system()/pruned_ss_loglik() exactly).
  Gv_obs  <- Gv_base[obs_idx, , drop = FALSE]
  dGv_obs <- dGv[obs_idx, , drop = FALSE]

  ## ---- Cr0 and its derivative -----------------------------------------------
  Cr0  <- .order2_cov_r(numeric(n_s), Sigma_x, Sigma_e)
  dCr0 <- .pgc_d_cov_r0(Sigma_e, dSigma_e, Sigma_x, d_Sigma_x)

  ## ---- d(QQ) = d(G Cr0 G') ---------------------------------------------------
  dQQ <- dGmat %*% Cr0 %*% t(G_base) + G_base %*% dCr0 %*% t(G_base) +
         G_base %*% Cr0 %*% t(dGmat)
  dQQ <- (dQQ + t(dQQ)) * 0.5

  ## ---- d(HH) = d(Gv Cr0 Gv')  (me_variance jitter is theta-independent) -----
  dHH <- dGv_obs %*% Cr0 %*% t(Gv_obs) + Gv_obs %*% dCr0 %*% t(Gv_obs) +
         Gv_obs %*% Cr0 %*% t(dGv_obs)
  dHH <- (dHH + t(dHH)) * 0.5

  ## ---- d(SS) = d(G Cr0 Gv') ---------------------------------------------------
  dSS <- dGmat %*% Cr0 %*% t(Gv_obs) + G_base %*% dCr0 %*% t(Gv_obs) +
         G_base %*% Cr0 %*% t(dGv_obs)

  list(
    dTlin = dTlin, dZZ = dZZ, d_dy_no_ys = d_dy_no_ys, d_c_drift = d_c_drift,
    dQQ = dQQ, dHH = dHH, dSS = dSS,
    G_base = G_base, Gv_base = Gv_base, Cr0 = Cr0, dCr0 = dCr0
  )
}


#' d(mu0)/dtheta_j and d(Sxi0)/dtheta_j
#'
#' mu0 solves the linear system \code{(I - Tlin) mu0 = c_drift}, so
#' \code{d mu0 = solve(I - Tlin, d_c_drift + dTlin \%*\% mu0)}.
#'
#' Sxi0 solves the discrete Lyapunov fixed point \code{Sxi0 = Tlin Sxi0
#' Tlin' + QQ}; by the same pattern as \code{.o2sd_dSigma_x()}
#' (R/gradient-solution-deriv-order2.R:191), differentiating gives another
#' Lyapunov solve with a perturbed RHS:
#' \code{d Sxi0 = solve_lyapunov(Tlin, dTlin Sxi0 Tlin' + Tlin Sxi0 dTlin' + dQQ)}.
#'
#' @noRd
.pruned_stationary_init_deriv <- function(Tlin, QQ, mu0, Sxi0, dTlin, d_c_drift, dQQ) {
  d_dim <- nrow(Tlin)
  dmu0 <- as.numeric(solve(diag(d_dim) - Tlin, d_c_drift + as.numeric(dTlin %*% mu0)))

  dB <- dTlin %*% Sxi0 %*% t(Tlin) + Tlin %*% Sxi0 %*% t(dTlin) + dQQ
  dSxi0 <- solve_lyapunov(Tlin, dB)
  dSxi0 <- (dSxi0 + t(dSxi0)) * 0.5

  list(dmu0 = dmu0, dSxi0 = dSxi0)
}


## NOTE (C5): .pruned_d_ghss_sigma_channel() used to live here -- a local,
## closed-form re-solve (solve_perturbation_order2(..., Sigma_e = dSigma_e))
## that supplied the Sigma_e channel of d(ghss)/dtheta_j because
## solution_derivatives_order2()'s own d_ghss held Sigma_e fixed while
## differentiating. C5 fixed that gap AT THE SOURCE
## (R/gradient-solution-deriv-order2.R's d(ghss)/dtheta now includes the
## fp*ghuu*dvSe + H2(T_up,T_up)*dvSe terms), so sd2$derivs[[pnm]]$d_ghss
## already carries this channel and the local workaround has been removed
## (it would otherwise double-count).


#' End-to-end d(pruned loglik)/d(theta_j) for ONE parameter, chaining phase-a
#' (\code{.pruned_kf_correlated_adjoint}) through phase-b (this file).
#'
#' Builds the nine base augmented-system matrices exactly as
#' \code{pruned_ss_loglik()} does (mirroring its construction verbatim),
#' calls the phase-a adjoint ONCE to get \code{loglik} and
#' \code{grad$<matrix>} (dL/dX_ij), then contracts each matrix gradient
#' against this file's \code{d(X)/dtheta_j} via a Frobenius inner product
#' (\code{sum(grad$X * dX)}), summing over all nine inputs.
#'
#' @param pss        \code{pruned_ss} object.
#' @param Y          n_obs x T (or T x n_obs) data matrix.
#' @param obs_vars   Character vector of observed variable names.
#' @param sd2        Output of \code{solution_derivatives_order2(...,
#'   param_names)} at the SAME base point as \code{pss}.
#' @param dSigma_e_list  Named list of dSigma_e (n_u x n_u) per parameter in
#'   \code{sd2$param_names} (the caller supplies these; the shock-covariance
#'   construction lives outside \code{sd2}'s scope).
#' @param me_variance Scalar measurement-error jitter (default 0).
#' @param model,compiled,ss,dr1,params  Unused (kept for call-site
#'   compatibility). Previously supplied the base first-order solve inputs
#'   needed for a local, closed-form re-solve of the Sigma_e channel of
#'   \code{d(ghss)/dtheta_j} (\code{.pruned_d_ghss_sigma_channel()}); since
#'   C5 fixed that channel AT THE SOURCE in
#'   \code{solution_derivatives_order2()} (R/gradient-solution-deriv-order2.R),
#'   \code{sd2$derivs[[pnm]]$d_ghss} already carries it and no local
#'   correction/re-solve is needed here any more.
#' @return A list: \code{loglik} (scalar), \code{grad} (named numeric vector,
#'   one entry per \code{sd2$param_names}; \code{NA} for any parameter whose
#'   upstream solution derivative has \code{ok == FALSE}).
#' @noRd
.pruned_ss_loglik_grad_chain <- function(pss, Y, obs_vars, sd2, dSigma_e_list,
                                         me_variance = 0,
                                         model = NULL, compiled = NULL,
                                         ss = NULL, dr1 = NULL, params = NULL) {
  stopifnot(inherits(pss, "pruned_ss"))
  sys <- pss$sys
  dr2 <- pss$dr

  obs_idx <- match(obs_vars, pss$endo_names)
  if (any(is.na(obs_idx)))
    stop(".pruned_ss_loglik_grad_chain: obs_vars not found in pss: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))
  n_obs <- length(obs_vars)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)

  Tlin  <- sys$Tlin
  Dxi_f <- sys$Dxi
  Gv_f  <- sys$Gv
  d_dim <- sys$d

  ZZ  <- Dxi_f[obs_idx, , drop = FALSE]
  Gv  <- Gv_f[obs_idx,  , drop = FALSE]
  d_y <- pss$ys[obs_vars] + 0.5 * sys$ghss[obs_idx] + sys$c_v[obs_idx]
  c_drift <- sys$cc + sys$c_u

  Sigma_x <- solve_lyapunov(sys$hx, sys$hu %*% sys$Sigma_e %*% t(sys$hu))
  Cr0 <- .order2_cov_r(numeric(sys$n_s), Sigma_x, sys$Sigma_e)

  QQ <- sys$G %*% Cr0 %*% t(sys$G); QQ <- (QQ + t(QQ)) * 0.5
  HH <- Gv %*% Cr0 %*% t(Gv);       HH <- (HH + t(HH)) * 0.5
  if (me_variance > 0) HH <- HH + me_variance * diag(n_obs)
  SS <- sys$G %*% Cr0 %*% t(Gv)

  mu0  <- as.numeric(solve(diag(d_dim) - Tlin, c_drift))
  Sxi0 <- solve_lyapunov(Tlin, QQ); Sxi0 <- (Sxi0 + t(Sxi0)) * 0.5

  ad <- .pruned_kf_correlated_adjoint(Y, Tlin, ZZ, as.numeric(d_y), c_drift,
                                      QQ, HH, SS, mu0, Sxi0)

  param_names <- sd2$param_names
  grad <- setNames(rep(NA_real_, length(param_names)), param_names)

  for (pnm in param_names) {
    d1 <- sd2$first$derivs[[pnm]]
    d2 <- sd2$derivs[[pnm]]
    if (is.null(d1) || !isTRUE(d1$ok) || is.null(d2) || !isTRUE(d2$ok)) next
    dSigma_e <- dSigma_e_list[[pnm]]
    if (is.null(dSigma_e)) dSigma_e <- matrix(0, sys$n_u, sys$n_u)

    chn <- tryCatch(
      .pruned_aug_system_deriv(sys, dr2, obs_idx, d1$dG, d1$dH,
                               d2$d_ghxx, d2$d_ghxu, d2$d_ghuu, d2$d_ghss,
                               dSigma_e, Sigma_x, d2$d_Sigma_x),
      error = function(e) NULL)
    if (is.null(chn)) next

    d_dy_full <- chn$d_dy_no_ys + d1$dys[obs_vars]

    st <- .pruned_stationary_init_deriv(Tlin, QQ, mu0, Sxi0, chn$dTlin,
                                        chn$d_c_drift, chn$dQQ)

    gj <- sum(ad$grad$Tlin    * chn$dTlin) +
          sum(ad$grad$ZZ      * chn$dZZ) +
          sum(ad$grad$d_y     * d_dy_full) +
          sum(ad$grad$c_drift * chn$d_c_drift) +
          sum(ad$grad$QQ      * chn$dQQ) +
          sum(ad$grad$HH      * chn$dHH) +
          sum(ad$grad$SS      * chn$dSS) +
          sum(ad$grad$mu0     * st$dmu0) +
          sum(ad$grad$Sxi0    * st$dSxi0)

    grad[pnm] <- gj
  }

  list(loglik = ad$loglik, grad = grad)
}
