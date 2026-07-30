## R/metric-whittle-fim.R
## --------------------------------------------------------------------------
## Stage 2 (manifold-MCMC roadmap): Whittle (frequency-domain) Fisher
## information matrix as a position-dependent metric for smMALA.
##
## Theory:
##   The asymptotic Fisher information of a stationary Gaussian series is
##   I(theta)_{jk} = (T / 4*pi) * int_{-pi}^{pi}
##                     Re tr[ f(w)^{-1} d_j f(w) f(w)^{-1} d_k f(w) ] dw
##
##   f(w) is the spectral density; d_k f(w) is its parameter derivative.
##
## Quadrature convention:
##   Use Fourier grid omega_j = 2*pi*j/T over the half-period [0, pi]
##   (symmetry of the real Gaussian-process spectral density: the integral over
##   [-pi,pi] is 2*integral_0^pi). The TRAPEZOIDAL rule over [0,pi] gives
##     (T/4*pi)*2*(2*pi/T)*[ (1/2) g(0) + sum_{interior} g + (1/2) g(pi) ]
##       = (1/2) g(0) + sum_{interior} g + (1/2) g(pi).
##   So interior frequencies carry weight 1, but the endpoints omega=0 and the
##   Nyquist omega=pi each carry weight 1/2. A naive "weight 1 everywhere" sum
##   that drops omega=0 and full-weights pi leaves an O(1) endpoint error
##   (1/2)[g(0)-g(pi)] that does NOT vanish with T and blows up near a unit
##   root -- see whittle_fim() for the endpoint handling.
##
## Implementation:
##   Uses the STANDARD spectral density path (not the debiased EI path):
##   - f(w): .spectral_density_core() / .whittle_spectral_density()
##   - d_k f(w): computed via the same formula as .whittle_loglik_grad
##     standard path (lines 757-827 of whittle-likelihood.R):
##     dH = dZZ*B + ZZ*z*A^{-1}*(dTT*B + dRR) + dDD
##     dS = dH * Se * H^H + H * dSe * H^H + H * Se * dH^H
##   This is the EXACT spectral derivative (not Fejer-approximated), giving
##   the correct quadrature accuracy for a finite grid.
##
## Unit-root fallback:
##   When spectral radius >= 1, A=(I-TT*z) is near-singular.  In that case
##   whittle_fim() returns a valid metric via the fallback_metric argument
##   (default: prior Hessian or scaled identity via softabs_metric).
## --------------------------------------------------------------------------


#' Whittle Fisher information matrix (internal)
#'
#' Assembles the asymptotic frequency-domain Fisher information matrix from
#' the spectral density and its parameter derivatives.  No new derivatives
#' required beyond the existing solution-mover derivatives.
#'
#' @details
#' \strong{Formula:}
#' \deqn{
#'   I(\theta)_{jk} = \frac{T}{4\pi}
#'     \int_{-\pi}^{\pi} \mathrm{Re}\,\mathrm{tr}\bigl[
#'       f(\omega)^{-1}\,\partial_j f(\omega)\,
#'       f(\omega)^{-1}\,\partial_k f(\omega)
#'     \bigr]\,d\omega
#' }
#' Approximated by the trapezoidal rule on the Fourier grid
#' \eqn{\omega_j = 2\pi j / T} over the half-period \eqn{[0,\pi]}: interior
#' frequencies carry weight 1, while the endpoints \eqn{\omega=0} and the
#' Nyquist \eqn{\omega=\pi} each carry weight \eqn{1/2} (dropping \eqn{\omega=0}
#' or full-weighting \eqn{\pi} introduces an O(1) endpoint bias that grows near
#' a unit root). Recovers the analytic AR(1) Fisher information to machine
#' precision away from the unit-root boundary.
#'
#' @param TT       n_state x n_state state-transition matrix.
#' @param RR       n_state x n_shock shock-impact matrix.
#' @param ZZ       n_obs   x n_state observation matrix (lagged convention).
#' @param DD       n_obs   x n_shock direct-impact matrix.
#' @param Sigma_e  n_shock x n_shock shock covariance.
#' @param dss_list Named list; each element is a list with fields
#'   \code{dTT}, \code{dRR}, \code{dZZ}, \code{dDD}, \code{dSigma_e}
#'   (the per-parameter derivatives of the SSM matrices).
#' @param omega_grid Numeric vector of positive angular frequencies in (0, pi].
#'   Default: \code{2*pi*(1:floor(T_obs/2))/T_obs}.
#' @param T_obs    Sample length (integer).
#' @param prior_hess Optional n_par x n_par matrix: the \emph{negative}
#'   log-prior Hessian (\eqn{-\partial^2 \log p(\theta) / \partial\theta^2}).
#'   When supplied, added to the Whittle FIM to form the posterior metric.
#' @param fallback_metric Optional \code{list(G, G_inv, L, logdet)} used when
#'   the Whittle FIM is degenerate (unit-root boundary).
#' @param softabs_alpha Alpha for \code{softabs_metric} regularisation (1e6).
#' @param me_variance Scalar >= 0; ME variance added to spectral density diagonal.
#'
#' @return A named list:
#'   \item{G}{n_par x n_par positive-definite metric matrix}
#'   \item{G_inv}{Inverse of G}
#'   \item{L}{Upper Cholesky factor: \code{chol(G)}}
#'   \item{logdet}{log|G|}
#' @noRd
whittle_fim <- function(TT, RR, ZZ, DD, Sigma_e,
                         dss_list,
                         omega_grid,
                         T_obs,
                         prior_hess    = NULL,
                         fallback_metric = NULL,
                         softabs_alpha = 1e6,
                         me_variance   = 0) {

  n_par   <- length(dss_list)
  np      <- n_par
  nm      <- names(dss_list)
  ## Trapezoidal quadrature over the half-period [0, pi]. The positive-half
  ## Fourier sum approximates (T/4pi) int_{-pi}^{pi} ONLY with the correct
  ## endpoint weights: omega = 0 carries weight 1/2 and the Nyquist point
  ## omega = pi carries weight 1/2 (not 1); interior points carry weight 1.
  ## Dropping omega = 0 and giving pi full weight (the naive "weight 1 per
  ## positive Fourier frequency") leaves an O(1) endpoint error
  ## = (1/2)[g(0) - g(pi)] that does NOT vanish with T and BLOWS UP near a
  ## unit root (g(0) ~ 4/(1-rho)^2). We therefore augment the grid with
  ## omega = 0 if absent and weight both endpoints by 1/2.
  om_use <- if (any(abs(omega_grid) < 1e-12)) omega_grid else c(0, omega_grid)
  w_use  <- rep(1, length(om_use))
  w_use[abs(om_use)      < 1e-12] <- 0.5   # omega = 0 endpoint
  w_use[abs(om_use - pi) < 1e-9]  <- 0.5   # omega = pi (Nyquist) endpoint
  J       <- length(om_use)
  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)
  n_exo   <- ncol(RR)

  ## --- Zero-pad missing derivative components ---------------------------------
  zero_TT <- matrix(0, n_state, n_state)
  zero_RR <- matrix(0, n_state, n_exo)
  zero_ZZ <- matrix(0, n_obs, n_state)
  zero_DD <- matrix(0, n_obs, n_exo)
  zero_Se <- matrix(0, n_exo, n_exo)

  get_d <- function(d, field, zero) {
    if (is.null(d)) return(zero)
    v <- d[[field]]
    if (is.null(v)) zero else v
  }

  ## Extract derivative matrices for all params upfront
  dTT_arr <- lapply(dss_list, function(d) get_d(d, "dTT",      zero_TT))
  dRR_arr <- lapply(dss_list, function(d) get_d(d, "dRR",      zero_RR))
  dZZ_arr <- lapply(dss_list, function(d) get_d(d, "dZZ",      zero_ZZ))
  dDD_arr <- lapply(dss_list, function(d) get_d(d, "dDD",      zero_DD))
  dSe_arr <- lapply(dss_list, function(d) get_d(d, "dSigma_e", zero_Se))

  ## --- Accumulate FIM via spectral quadrature ---------------------------------
  ## Quadrature weight per frequency = 1 (see derivation above).
  ## Inner loop: per-frequency, compute Si_j = f^{-1}(w_j), then per-param
  ## dS_k, and accumulate FIM_{jk} += Re tr(Si dS_j Si dS_k).
  ##
  ## This mirrors .whittle_loglik_grad standard path (whittle-likelihood.R:757-827)
  ## but accumulates the outer product of gradients rather than summing them.

  FIM <- matrix(0, np, np)

  for (idx in seq_len(J)) {
    om_j <- om_use[idx]
    w_j  <- w_use[idx]
    z_j  <- exp(-1i * om_j)

    ## A_j = I - TT * z_j
    A_j     <- diag(n_state) - TT * z_j
    rcond_A <- rcond(A_j)

    if (!is.finite(rcond_A) || rcond_A <= .Machine$double.eps * 1e4) {
      ## Near-unit-root: skip this frequency
      next
    }

    ## B_j = z_j * A^{-1} * RR  (n_state x n_exo)
    B_j <- z_j * solve(A_j, RR)
    ## H_j = ZZ * B_j + DD  (n_obs x n_exo)
    H_j <- ZZ %*% B_j + DD

    ## S_j = H_j * Sigma_e * H_j^H  [+ me_variance * I]  (n_obs x n_obs complex)
    S_j <- H_j %*% Sigma_e %*% Conj(t(H_j))
    if (me_variance > 0) S_j <- S_j + me_variance * diag(n_obs)

    ## Invert S_j (eigendecompose for robustness, same as .whittle_loglik;
    ## .eigen_hermitian_safe avoids the zheev segfault under vecLib BLAS)
    ev_j <- tryCatch(.eigen_hermitian_safe(S_j), error = function(e) NULL)
    if (is.null(ev_j)) next
    ev_max_j <- max(ev_j$values)
    if (ev_max_j <= 0) next
    ev_clamped_j <- pmax(ev_j$values, .Machine$double.eps * ev_max_j)
    V_j  <- ev_j$vectors
    Si_j <- V_j %*% diag(1 / ev_clamped_j, nrow = length(ev_clamped_j)) %*% Conj(t(V_j))

    ## Per-param: compute SidS_k = Si_j %*% dS_k  (n_obs x n_obs complex)
    ## Store them to accumulate the outer product FIM_{jk} += Re tr(SidS_j SidS_k)
    SidS <- vector("list", np)
    for (k in seq_len(np)) {
      dTT_k <- dTT_arr[[k]]
      dRR_k <- dRR_arr[[k]]
      dZZ_k <- dZZ_arr[[k]]
      dDD_k <- dDD_arr[[k]]
      dSe_k <- dSe_arr[[k]]

      ## dH_k (same formula as .whittle_loglik_grad standard path, line 808-810):
      ##   inner = dTT_k * B_j + dRR_k
      ##   dG = dZZ_k * B_j + ZZ * z_j * A^{-1} * inner
      ##   dH = dG + dDD_k
      inner_k <- dTT_k %*% B_j + dRR_k
      dG_k    <- dZZ_k %*% B_j + ZZ %*% (z_j * solve(A_j, inner_k))
      dH_k    <- dG_k + dDD_k

      ## dS_k = dH_k Se H^H + H dSe_k H^H + H Se dH_k^H
      dS_k <- dH_k %*% Sigma_e  %*% Conj(t(H_j)) +
              H_j  %*% dSe_k    %*% Conj(t(H_j)) +
              H_j  %*% Sigma_e  %*% Conj(t(dH_k))

      SidS[[k]] <- Si_j %*% dS_k
    }

    ## Accumulate FIM: FIM[j_par, k_par] += Re tr(SidS_j SidS_k^T)
    ## = Re sum_{ab} (SidS_j)_{ab} * (SidS_k)_{ba}  [trace of product]
    ## For efficiency: Re tr(A B) = Re sum(A * t(B))   [elementwise Frobenius]
    for (j_par in seq_len(np)) {
      for (k_par in j_par:np) {
        val <- w_j * Re(sum(SidS[[j_par]] * t(SidS[[k_par]])))
        FIM[j_par, k_par] <- FIM[j_par, k_par] + val
        if (j_par != k_par)
          FIM[k_par, j_par] <- FIM[k_par, j_par] + val
      }
    }
  }

  ## Force symmetry
  FIM <- (FIM + t(FIM)) / 2

  ## Set names
  if (!is.null(nm)) {
    rownames(FIM) <- nm
    colnames(FIM) <- nm
  }

  ## --- Add prior Hessian (posterior metric) -----------------------------------
  if (!is.null(prior_hess)) {
    stopifnot(identical(dim(prior_hess), dim(FIM)))
    prior_sym <- (prior_hess + t(prior_hess)) / 2
    FIM <- FIM + prior_sym
  }

  ## --- Regularise via SoftAbs -------------------------------------------------
  m <- tryCatch(
    softabs_metric(FIM, alpha = softabs_alpha),
    error = function(e) NULL
  )

  if (is.null(m)) {
    return(.whittle_fim_fallback(np, nm, fallback_metric, prior_hess, softabs_alpha))
  }

  ## Restore parameter names (softabs_metric drops dimnames via eigen reassembly)
  if (!is.null(nm)) {
    rownames(m$G)     <- nm;  colnames(m$G)     <- nm
    rownames(m$G_inv) <- nm;  colnames(m$G_inv) <- nm
    rownames(m$L)     <- nm;  colnames(m$L)     <- nm
  }

  m
}


#' Assemble a dss_list (per-parameter SSM derivatives) at a point (internal)
#'
#' Builds the \code{dss_list} argument expected by \code{whittle_fim()}:
#' a named list, one entry per estimated parameter, each holding
#' \code{dTT}/\code{dRR}/\code{dZZ}/\code{dDD}/\code{dSigma_e} -- the
#' derivatives of the solved state-space matrices w.r.t. that parameter,
#' evaluated at \code{theta}.
#'
#' @details
#' Reuses the existing solution-derivative machinery
#' (\code{solution_derivatives()}, Tier 11/12 implicit-differentiation path)
#' for parameters that move the decision rule, exactly as the "implicit"
#' branch of \code{make_posterior_grad()} does (see
#' R/analytic-gradient.R around the \code{d_ss_list} construction). Any
#' parameter whose solution-derivative block is unusable
#' (\code{solution_derivatives()} failed or returned \code{ok = FALSE}) falls
#' back to a central finite difference on the solved SSM matrices directly
#' (central FD is acceptable here: this assembler runs once at the mode).
#' \code{dSigma_e} is always obtained by central FD of \code{.get_shock_cov()}
#' (cheap, and exact-zero for parameters absent from the shocks block),
#' mirroring \code{.dSigma_e_fd()} in \code{make_posterior_grad()}.
#'
#' @param model    dynhr_mod.
#' @param compiled dynhr_compiled.
#' @param theta    Named numeric parameter vector (typically the mode).
#' @param obs_vars Character vector of observed variable names.
#' @return list(dss_list = <named list>, TT=, RR=, ZZ=, DD=, Sigma_e=,
#'   ok = <logical, TRUE iff the base solve succeeded>).
#' @noRd
.assemble_dss_list_at_mode <- function(model, compiled, theta, obs_vars) {
  par_names <- names(theta)
  np        <- length(par_names)
  exo       <- model$varexo_names
  sys_cache <- cache_system_structure(compiled)

  .solve_dr <- function(th) {
    params <- .apply_theta_to_params(model, th)
    ss <- solve_steady_state(model, compiled, params, verbose = FALSE)
    if (is.null(ss) || !isTRUE(ss$converged)) return(NULL)
    params <- ss$params %||% params
    sys <- extract_system_matrices_fast(sys_cache, ss$ss, params)
    dr  <- .solve_from_system(sys, model, compiled, ss$ss, params, FALSE)
    if (is.null(dr) || !isTRUE(dr$bk_satisfied)) return(NULL)
    list(dr = dr, params = params)
  }

  base <- .solve_dr(theta)
  if (is.null(base))
    return(list(dss_list = NULL, ok = FALSE))

  dr      <- base$dr
  params  <- base$params
  si      <- dr$state_idx
  oi      <- match(obs_vars, dr$endo_names)
  TT      <- dr$ghx[si, , drop = FALSE]
  RR      <- dr$ghu[si, , drop = FALSE]
  ZZ      <- dr$ghx[oi, , drop = FALSE]
  DD      <- dr$ghu[oi, , drop = FALSE]
  Sigma_e <- .get_shock_cov(model, exo, params)

  ## dSigma_e/dtheta_nm by central FD of .get_shock_cov (mirrors
  ## .dSigma_e_fd() in make_posterior_grad(); exact-zero for parameters
  ## absent from the shocks block).
  .dSigma_e_fd1 <- function(nm) {
    h  <- 1e-6 * max(abs(theta[[nm]]), 1e-3)
    tp <- theta; tm <- theta
    tp[nm] <- tp[nm] + h; tm[nm] <- tm[nm] - h
    pp <- .apply_theta_to_params(model, tp, params)
    pm <- .apply_theta_to_params(model, tm, params)
    (.get_shock_cov(model, exo, pp) - .get_shock_cov(model, exo, pm)) / (2 * h)
  }

  ## Central FD on the solved SSM matrices directly (fallback path; also the
  ## sole path when solution_derivatives() is unavailable/fails for a param).
  .ssm_fd1 <- function(nm) {
    h  <- 1e-5 * max(abs(theta[[nm]]), 1e-3)
    tp <- theta; tm <- theta
    tp[nm] <- tp[nm] + h; tm[nm] <- tm[nm] - h
    dp <- .solve_dr(tp); dm <- .solve_dr(tm)
    if (is.null(dp) || is.null(dm)) return(NULL)
    list(
      dTT = (dp$dr$ghx[si, , drop = FALSE] - dm$dr$ghx[si, , drop = FALSE]) / (2 * h),
      dRR = (dp$dr$ghu[si, , drop = FALSE] - dm$dr$ghu[si, , drop = FALSE]) / (2 * h),
      dZZ = (dp$dr$ghx[oi, , drop = FALSE] - dm$dr$ghx[oi, , drop = FALSE]) / (2 * h),
      dDD = (dp$dr$ghu[oi, , drop = FALSE] - dm$dr$ghu[oi, , drop = FALSE]) / (2 * h)
    )
  }

  ## Classify parameters exactly as make_posterior_grad(): "sigma-like" params
  ## leave the decision rule unchanged (pure certainty-equivalence movers).
  is_sigma <- rep(FALSE, np); names(is_sigma) <- par_names
  g0 <- list(ghx = dr$ghx, ghu = dr$ghu)
  for (k in seq_len(np)) {
    th <- theta; h <- 1e-5 * max(abs(th[k]), 1e-3)
    th[k] <- th[k] + h
    d2 <- .solve_dr(th)
    if (!is.null(d2)) {
      dchg <- max(abs(d2$dr$ghx - g0$ghx), abs(d2$dr$ghu - g0$ghu))
      is_sigma[k] <- dchg < 1e-10
    }
  }
  sig_names <- par_names[is_sigma]
  num_names <- par_names[!is_sigma]

  dss_list <- vector("list", np)
  names(dss_list) <- par_names

  for (nm in sig_names)
    dss_list[[nm]] <- list(dSigma_e = .dSigma_e_fd1(nm))

  sd_res <- NULL
  if (length(num_names) > 0) {
    sd_res <- tryCatch(
      solution_derivatives(model, compiled, dr, params,
                            param_names = num_names, obs_vars = obs_vars),
      error = function(e) NULL
    )
    for (nm in num_names) {
      d <- if (!is.null(sd_res)) sd_res$derivs[[nm]] else NULL
      if (!is.null(d) && isTRUE(d$ok)) {
        dss_list[[nm]] <- list(dTT = d$dTT, dRR = d$dRR, dZZ = d$dZZ,
                                dDD = d$dDD,
                                dSigma_e = .dSigma_e_fd1(nm))
      } else {
        ## Fall back to direct central FD on the solved SSM matrices.
        fd <- .ssm_fd1(nm)
        dss_list[[nm]] <- if (!is.null(fd)) {
          c(fd, list(dSigma_e = .dSigma_e_fd1(nm)))
        } else {
          list(dSigma_e = .dSigma_e_fd1(nm))   ## last resort: zero solution-block
        }
      }
    }
  }

  list(dss_list = dss_list, TT = TT, RR = RR, ZZ = ZZ, DD = DD,
       Sigma_e = Sigma_e, ok = TRUE)
}


## Helper: return a fallback metric when the Whittle FIM is degenerate.
## Priority: supplied fallback_metric -> prior Hessian floor -> identity.
#' @noRd
.whittle_fim_fallback <- function(np, nm, fallback_metric, prior_hess,
                                   softabs_alpha) {
  if (!is.null(fallback_metric)) {
    return(fallback_metric)
  }
  ## Fallback: prior Hessian floor, or identity
  if (!is.null(prior_hess)) {
    m <- tryCatch(
      softabs_metric(prior_hess + diag(1e-4, np), alpha = softabs_alpha),
      error = function(e) NULL
    )
    if (!is.null(m)) {
      if (!is.null(nm)) {
        rownames(m$G)     <- nm;  colnames(m$G)     <- nm
        rownames(m$G_inv) <- nm;  colnames(m$G_inv) <- nm
        rownames(m$L)     <- nm;  colnames(m$L)     <- nm
      }
      return(m)
    }
  }
  ## Last resort: identity
  G <- diag(np)
  if (!is.null(nm)) { rownames(G) <- nm; colnames(G) <- nm }
  list(G = G, G_inv = G, L = diag(np), logdet = 0)
}
