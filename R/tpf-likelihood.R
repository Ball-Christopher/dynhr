## R/tpf-likelihood.R
## --------------------------------------------------------------------------
## Tempered Particle Filter (TPF) likelihood for nonlinear DSGE models.
##
## Reference: Herbst & Schorfheide (2019), "Tempered Particle Filtering",
##   Journal of Econometrics 210(1):26-44.
##
## The filter operates on the pruned order-2 state space (Andreasen,
## Fernandez-Villaverde & Rubio-Ramirez 2018, RES 85:1-49). The state vector
## is s_t = (x1_t, x2_t), dimension 2*n_state, where x1 is the first-order
## state and x2 is the second-order correction.
##
## Tempering instrument: measurement error variance is inflated to
## me_variance/phi at stage n (phi in [0,1]) within each period. The phi
## schedule is computed by adaptive bisection (reusing .smc_next_lambda from
## R/sampler-smc.R). Bootstrap proposal (state-transition density) is used.
## --------------------------------------------------------------------------


## ---- C++ dispatch checks --------------------------------------------------

#' Check whether the Rcpp TPF propagation kernels are available
#' @noRd
.HAS_RCPP_TPF <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("tpf_propagate_particles", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

#' Check whether the C++ per-period tempering loop is available
#' @noRd
.HAS_RCPP_TPF_PERIOD <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("tpf_run_period_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}


## ---- Bit-identical vectorized Kronecker helper ----------------------------

#' Row-index pair for a per-column pairwise Kronecker product
#'
#' For matrices A (na x N) and B (nb x N), \code{kronecker(A[,i], B[,i])}
#' for each column i is exactly \code{A[ia, i] * B[ib, i]} with
#' \code{ia = rep(seq_len(na), each = nb)}, \code{ib = rep(seq_len(nb), times
#' = na)} (A is the SLOW/outer factor, B is the FAST/inner factor — this
#' matches R's \code{kronecker()} convention). Verified \code{identical()}
#' TRUE against a \code{for (i in seq_len(N)) kronecker(A[,i], B[,i])} loop;
#' ~230x faster at N=150, n=2. Triple Kronecker products MUST group
#' inner-first (\code{a %x% (b %x% c)}, not left-to-right) to stay
#' bit-identical to \code{kronecker(a, kronecker(b, c))} (floating-point
#' reassociation is NOT identical).
#'
#' @param na,nb Integer row-counts of A, B.
#' @return list(ia, ib) index vectors of length na*nb.
#' @noRd
.tpf_kron_idx <- function(na, nb) {
  list(ia = rep(seq_len(na), each = nb), ib = rep(seq_len(nb), times = na))
}

#' Vectorized per-column pairwise Kronecker product: A[,i] %x% B[,i]
#' @noRd
.tpf_kron2 <- function(A, B, idx) {
  A[idx$ia, , drop = FALSE] * B[idx$ib, , drop = FALSE]
}


## ---- Pure-R propagation kernels ------------------------------------------

#' Propagate N particles one period forward (pure-R fallback)
#'
#' Implements the pruned second-order state transition for all N particles
#' simultaneously. Kronecker convention: columns of hxx are ordered
#' (state FAST x state SLOW), matching (x1 %x% x1); columns of hxu are
#' (state FAST x exo SLOW), matching (e %x% x1). This replicates exactly
#' simulate_model_order2() lines 919-923 in R/solve-perturbation-order2.R.
#'
#' @param particles  Matrix (2*n_s) x N; rows [1:n_s]=x1, [n_s+1:2*n_s]=x2.
#' @param shocks     Matrix n_e x N of drawn shocks.
#' @param hx,hu      State-row slices of ghx, ghu (n_s x n_s, n_s x n_e).
#' @param hxx,hxu,huu Second-order coefficient matrices (may be all-zero).
#' @param hss        n_s-vector (second-order sigma term).
#' @return Updated particles matrix (2*n_s) x N.
#' @noRd
.tpf_propagate_R <- function(particles, shocks, hx, hu,
                              hxx, hxu, huu, hss) {
  n_s <- nrow(hx)
  N   <- ncol(particles)
  n_e <- ncol(hu)

  x1 <- particles[seq_len(n_s), , drop = FALSE]          # n_s x N
  x2 <- particles[seq_len(n_s) + n_s, , drop = FALSE]    # n_s x N

  ## First-order state update
  x1_new <- hx %*% x1 + hu %*% shocks   # n_s x N

  ## Second-order state update (pruned):
  ##   x2_new = hx x2_prev
  ##           + 0.5 hxx (x1_prev %x% x1_prev)
  ##           + hxu (e %x% x1_prev)
  ##           + 0.5 huu (e %x% e)
  ##           + 0.5 hss
  x2_new <- hx %*% x2

  if (any(hxx != 0)) {
    ## (x1_prev %x% x1_prev): fast-first Kronecker product, n_s^2 x N
    idx_xx  <- .tpf_kron_idx(n_s, n_s)
    kron_xx <- .tpf_kron2(x1, x1, idx_xx)
    x2_new <- x2_new + 0.5 * hxx %*% kron_xx
  }
  if (any(hxu != 0)) {
    ## (e %x% x1_prev): (exo SLOW x state FAST) — matches hxu col ordering
    ## hxu cols: (state FAST, exo SLOW) -> vector is (e %x% x1)
    idx_ex  <- .tpf_kron_idx(n_e, n_s)
    kron_ex <- .tpf_kron2(shocks, x1, idx_ex)
    x2_new <- x2_new + hxu %*% kron_ex
  }
  if (any(huu != 0)) {
    idx_ee  <- .tpf_kron_idx(n_e, n_e)
    kron_ee <- .tpf_kron2(shocks, shocks, idx_ee)
    x2_new <- x2_new + 0.5 * huu %*% kron_ee
  }
  if (any(hss != 0)) {
    x2_new <- x2_new + 0.5 * hss  # n_s vector broadcast to n_s x N
  }

  rbind(x1_new, x2_new)
}


#' Compute particle log observation weights (pure-R)
#'
#' DSGE state-dating convention: y_t = ZZ*(x1_{t-1} + x2_{t-1}) + DD*e_t
#'   + d_obs + ghss_obs (matches Dynare/dynhr kalman_filter convention;
#'   see R/kalman-filter.R:669-671).
#'
#' \code{particles} here is the PRE-propagation period t-1 state; the
#' period-t shocks \code{shocks_t} are the draws used in this period's
#' bootstrap proposal.  For the RWMH mutation pass, \code{shocks_t = NULL}
#' and \code{particles} is the post-mutation state with shocks already folded
#' in via \code{DD \%*\% shocks_t} treated as a separate additive term stored
#' in \code{DD_shocks_mean} (pre-computed by the caller).
#'
#' @param particles    Matrix (2*n_s) x N (pre-propagation states).
#' @param y_t          Numeric vector n_obs.
#' @param ZZ           Matrix n_obs x n_s.
#' @param DD           Matrix n_obs x n_e (shock-to-obs loading, = ghu[obs_idx,]).
#' @param shocks_t     Matrix n_e x N drawn shocks (NULL for mutation-only call).
#' @param d_obs        Numeric vector n_obs (SS mean).
#' @param ghss_obs     Numeric vector n_obs (= 0.5 * ghss[obs_idx]).
#' @param me_variance  Scalar > 0 (nominal).
#' @param phi          Tempering level in (0, 1].
#' @param DD_shock_mean n_obs-length pre-computed DD*e_t column (used when
#'   shocks_t is NULL and DD contribution is constant across mutation steps).
#' @param hxx_obs n_obs x n_s^2 obs-row slice of ghxx (may be NULL/all-zero).
#' @param hxu_obs n_obs x (n_e*n_s) obs-row slice of ghxu (may be NULL/all-zero).
#' @param huu_obs n_obs x n_e^2 obs-row slice of ghuu (may be NULL/all-zero).
#' @return Numeric vector of length N.
#' @noRd
.tpf_log_weights_R <- function(particles, y_t, ZZ, DD, shocks_t,
                                d_obs, ghss_obs, me_variance, phi,
                                DD_shock_mean = NULL,
                                hxx_obs = NULL, hxu_obs = NULL, huu_obs = NULL,
                                idx = NULL) {
  n_s <- ncol(ZZ)
  N   <- ncol(particles)
  n_e <- ncol(DD)

  x1 <- particles[seq_len(n_s), , drop = FALSE]
  x2 <- particles[seq_len(n_s) + n_s, , drop = FALSE]

  ## Observation mean: ZZ*(x1+x2) + DD*e_t + d_obs + ghss_obs  (n_obs x N)
  mu <- ZZ %*% (x1 + x2)

  if (!is.null(shocks_t)) {
    mu <- mu + DD %*% shocks_t        # n_obs x N
  } else if (!is.null(DD_shock_mean)) {
    mu <- mu + DD_shock_mean          # broadcast single column
  }

  ## Nonlinear obs terms (2026-08-04 obs-tensor fix): the pruned-model
  ## observation mean also carries the quadratic terms in (x1_prev, e_t) that
  ## simulate_model_order2's observable reconstruction keeps (mirrors the
  ## order-3 ground truth at solve-perturbation-order3.R:1367-1373 with
  ## x2_prev/x3_prev == 0): 0.5*ghxx_obs(x1(x)x1) + ghxu_obs(e(x)x1) +
  ## 0.5*ghuu_obs(e(x)e). Vectorized per-column, same pattern as
  ## .tpf_propagate_R's state-row kron products; guarded so LINEAR models
  ## (all-zero obs tensors) stay bit-identical to the pre-fix nesting pin.
  ## `idx` (optional): caller-precomputed .tpf_kron_idx() index lists,
  ## keyed $xx/$ex/$ee (period-constant -- depend only on n_s/n_e, not on
  ## the particle/shock draws). Hot callers (tpf_run_period's mutation
  ## micro-calls, N=1 columns) pass these to avoid rebuilding the same
  ## rep()-index vectors on every call; direct/test callers that omit `idx`
  ## get the identical values computed inline (bit-identical either way).
  if (!is.null(hxx_obs) && any(hxx_obs != 0)) {
    idx_xx  <- if (!is.null(idx)) idx$xx else .tpf_kron_idx(n_s, n_s)
    kron_xx <- .tpf_kron2(x1, x1, idx_xx)
    mu <- mu + 0.5 * hxx_obs %*% kron_xx
  }
  if (!is.null(shocks_t) && !is.null(hxu_obs) && any(hxu_obs != 0)) {
    idx_ex  <- if (!is.null(idx)) idx$ex else .tpf_kron_idx(n_e, n_s)
    kron_ex <- .tpf_kron2(shocks_t, x1, idx_ex)
    mu <- mu + hxu_obs %*% kron_ex
  }
  if (!is.null(shocks_t) && !is.null(huu_obs) && any(huu_obs != 0)) {
    idx_ee  <- if (!is.null(idx)) idx$ee else .tpf_kron_idx(n_e, n_e)
    kron_ee <- .tpf_kron2(shocks_t, shocks_t, idx_ee)
    mu <- mu + 0.5 * huu_obs %*% kron_ee
  }

  mu <- mu + d_obs + ghss_obs         # broadcast n_obs vectors over N columns

  sd_phi <- sqrt(me_variance / phi)
  lw <- dnorm(as.numeric(y_t), mean = mu, sd = sd_phi, log = TRUE)
  ## dnorm drops the dim attribute when its first argument is as long as
  ## `mu` (e.g. the N = 1 single-particle calls in the mutation step);
  ## restore it so colSums always sees a matrix.
  dim(lw) <- dim(mu)
  colSums(lw)
}


## ---- Pure-R propagation kernels: pruned ORDER 3 --------------------------

#' Propagate N particles one period forward on the pruned ORDER-3 state
#'
#' Implements the pruned third-order (AFVRR 2018) state transition for all N
#' particles simultaneously.  GROUND TRUTH: simulate_model_order3() in
#' R/solve-perturbation-order3.R lines 1322-1360 — every term below is a
#' direct vectorized transcription of that loop body (which the order-3
#' pruned-KF file R/pruned-state-space-order3.R also designates as the
#' authoritative recursion; where the AFVRR paper's eq (14) transcription
#' differs, THE CODE WINS — e.g. hxx %*% (x1 %x% x2) carries coefficient 1).
#'
#' Kronecker column-ordering conventions (R/solve-perturbation-order3.R
#' lines 1310-1320, inherited from solve-perturbation-order2.R):
#'   hxx  cols (state FAST x state SLOW)          -> vector (x1 %x% x2) etc.
#'   hxu  cols (state FAST x exo SLOW)             -> vector (e %x% x)
#'   hxxu cols (state1 FASTEST, state2, exo SLOW)  -> vector (e %x% x1 %x% x1)
#'   hxuu cols (state FASTEST, exo1, exo2 SLOW)    -> vector (e %x% e %x% x1)
#'   hxxx cols (s1 FAST .. s3 SLOW)                -> vector (x1 %x% x1 %x% x1)
#'   huuu cols (e1 FAST .. e3 SLOW)                -> vector (e %x% e %x% e)
#'
#' Shock-covariance convention: shocks arrive here already scaled by
#' L_e = chol(Sigma_e) (drawn by the caller), exactly as in the order-2
#' kernel — hu/huu/... are ghu-convention matrices that EXCLUDE Sigma_e.
#'
#' @param particles Matrix (3*n_s) x N; rows [1:n_s]=x1, [n_s+1:2n_s]=x2,
#'   [2n_s+1:3n_s]=x3.
#' @param shocks    Matrix n_e x N of drawn shocks (model units).
#' @param hx,hu     State-row slices of ghx, ghu.
#' @param hxx,hxu,huu Second-order state-row coefficient matrices.
#' @param hss       n_s-vector (second-order sigma term).
#' @param hxxx,hxxu,hxuu,huuu Third-order state-row coefficient matrices.
#' @param hxss,huss n_s x n_s / n_s x n_e sigma-cross matrices (or NULL).
#' @param hs3       n_s-vector third-order sigma term (or NULL).
#' @return Updated particles matrix (3*n_s) x N.
#' @noRd
.tpf_propagate3_R <- function(particles, shocks, hx, hu,
                               hxx, hxu, huu, hss,
                               hxxx, hxxu, hxuu, huuu,
                               hxss = NULL, huss = NULL, hs3 = NULL) {
  n_s <- nrow(hx)
  N   <- ncol(particles)
  n_e <- ncol(hu)

  x1 <- particles[seq_len(n_s), , drop = FALSE]           # n_s x N
  x2 <- particles[seq_len(n_s) + n_s, , drop = FALSE]     # n_s x N
  x3 <- particles[seq_len(n_s) + 2L * n_s, , drop = FALSE]# n_s x N

  ## (A) First-order update — solve-perturbation-order3.R:1329
  x1_new <- hx %*% x1 + hu %*% shocks

  ## (B) Second-order update — solve-perturbation-order3.R:1333-1339
  ## (identical to the order-2 kernel .tpf_propagate_R; the conditional-add
  ## structure is kept IDENTICAL so a model with zero higher tensors runs
  ## the same floating-point operation sequence as the order-2 kernel).
  x2_new <- hx %*% x2
  if (any(hxx != 0)) {
    kron_xx <- .tpf_kron2(x1, x1, .tpf_kron_idx(n_s, n_s))
    x2_new <- x2_new + 0.5 * hxx %*% kron_xx
  }
  if (any(hxu != 0)) {
    kron_ex <- .tpf_kron2(shocks, x1, .tpf_kron_idx(n_e, n_s))
    x2_new <- x2_new + hxu %*% kron_ex
  }
  if (any(huu != 0)) {
    kron_ee <- .tpf_kron2(shocks, shocks, .tpf_kron_idx(n_e, n_e))
    x2_new <- x2_new + 0.5 * huu %*% kron_ee
  }
  if (any(hss != 0)) x2_new <- x2_new + 0.5 * hss

  ## (C) Third-order update — solve-perturbation-order3.R:1344-1356
  ##   x3_new = hx x3_prev
  ##          + hxx (x1_prev %x% x2_prev)          [coefficient 1, see :1342-1343]
  ##          + hxu (e %x% x2_prev)
  ##          + 0.5  hxxu (e %x% x1_prev %x% x1_prev)
  ##          + 0.5  hxuu (e %x% e %x% x1_prev)
  ##          + 1/6  hxxx (x1_prev %x% x1_prev %x% x1_prev)
  ##          + 1/6  huuu (e %x% e %x% e)
  ##          [+ 0.5 hxss x1_prev + 0.5 huss e + 1/6 hs3 when present]
  x3_new <- hx %*% x3
  if (any(hxx != 0)) {
    kron_x12 <- .tpf_kron2(x1, x2, .tpf_kron_idx(n_s, n_s))
    x3_new <- x3_new + hxx %*% kron_x12
  }
  if (any(hxu != 0)) {
    kron_ex2 <- .tpf_kron2(shocks, x2, .tpf_kron_idx(n_e, n_s))
    x3_new <- x3_new + hxu %*% kron_ex2
  }
  if (any(hxxu != 0)) {
    ## kronecker(e, kronecker(x1, x1)) -- grouped inner-first for bit-identity
    kron_xx_i <- .tpf_kron2(x1, x1, .tpf_kron_idx(n_s, n_s))
    kron_exx  <- .tpf_kron2(shocks, kron_xx_i, .tpf_kron_idx(n_e, n_s * n_s))
    x3_new <- x3_new + 0.5 * hxxu %*% kron_exx
  }
  if (any(hxuu != 0)) {
    ## kronecker(e, kronecker(e, x1)) -- grouped inner-first
    kron_ex_i <- .tpf_kron2(shocks, x1, .tpf_kron_idx(n_e, n_s))
    kron_eex  <- .tpf_kron2(shocks, kron_ex_i, .tpf_kron_idx(n_e, n_e * n_s))
    x3_new <- x3_new + 0.5 * hxuu %*% kron_eex
  }
  if (any(hxxx != 0)) {
    ## kronecker(x1, kronecker(x1, x1)) -- grouped inner-first
    kron_xx_i2 <- .tpf_kron2(x1, x1, .tpf_kron_idx(n_s, n_s))
    kron_xxx   <- .tpf_kron2(x1, kron_xx_i2, .tpf_kron_idx(n_s, n_s * n_s))
    x3_new <- x3_new + (1 / 6) * hxxx %*% kron_xxx
  }
  if (any(huuu != 0)) {
    ## kronecker(e, kronecker(e, e)) -- grouped inner-first
    kron_ee_i <- .tpf_kron2(shocks, shocks, .tpf_kron_idx(n_e, n_e))
    kron_eee  <- .tpf_kron2(shocks, kron_ee_i, .tpf_kron_idx(n_e, n_e * n_e))
    x3_new <- x3_new + (1 / 6) * huuu %*% kron_eee
  }
  if (!is.null(hxss) && any(hxss != 0)) x3_new <- x3_new + 0.5 * hxss %*% x1
  if (!is.null(huss) && any(huss != 0)) x3_new <- x3_new + 0.5 * huss %*% shocks
  if (!is.null(hs3)  && any(hs3  != 0)) x3_new <- x3_new + (1 / 6) * hs3

  rbind(x1_new, x2_new, x3_new)
}


#' Compute particle log observation weights on the pruned ORDER-3 state
#'
#' Order-3 analogue of .tpf_log_weights_R.  Measurement design mirrors the
#' order-2 TPF (LINEAR observation in the pruned layers, quadratic/cubic
#' obs terms dropped) extended with the order-3 LINEAR and CONSTANT pieces
#' the pruned-KF3 keeps (R/pruned-state-space-order3.R:700 folds 0.5*ghxss
#' into the x1 loading; :719/:1097 fold (1/6)*ghs3 into the obs constant):
#'
#'   y_t = ZZ*(x1+x2+x3)_{t-1} + ZZ_xss*x1_{t-1} + DD*e_t
#'         + d_obs + ghss_obs + ghs3_obs
#'
#' where ZZ_xss = 0.5*ghxss[obs_idx,] and ghs3_obs = (1/6)*ghs3[obs_idx].
#' The nonlinear obs terms (ghxx/ghxu/ghuu/ghxxu/... at obs rows) are
#' dropped exactly as the order-2 TPF drops the quadratic obs terms.
#'
#' @param particles (3*n_s) x N pre-propagation states.
#' @param ZZ_xss    n_obs x n_s extra x1 loading (may be all-zero).
#' @param ghs3_obs  n_obs vector (1/6 * ghs3 at obs rows; may be zero).
#' @param hxx_obs  n_obs x n_s^2 obs-row slice of ghxx (NULL/all-zero ok).
#' @param hxu_obs  n_obs x (n_e*n_s) obs-row slice of ghxu (NULL/all-zero ok).
#' @param huu_obs  n_obs x n_e^2 obs-row slice of ghuu (NULL/all-zero ok).
#' @param hxxu_obs n_obs x (n_e*n_s^2) obs-row slice of ghxxu (NULL/all-zero ok).
#' @param hxuu_obs n_obs x (n_e^2*n_s) obs-row slice of ghxuu (NULL/all-zero ok).
#' @param hxxx_obs n_obs x n_s^3 obs-row slice of ghxxx (NULL/all-zero ok).
#' @param huuu_obs n_obs x n_e^3 obs-row slice of ghuuu (NULL/all-zero ok).
#' @param huss_obs n_obs x n_e obs-row slice of ghuss (NULL/all-zero ok).
#' @inheritParams .tpf_log_weights_R
#' @return Numeric vector of length N.
#' @noRd
.tpf_log_weights3_R <- function(particles, y_t, ZZ, DD, shocks_t,
                                 d_obs, ghss_obs, ZZ_xss, ghs3_obs,
                                 me_variance, phi,
                                 DD_shock_mean = NULL,
                                 hxx_obs = NULL, hxu_obs = NULL, huu_obs = NULL,
                                 hxxu_obs = NULL, hxuu_obs = NULL,
                                 hxxx_obs = NULL, huuu_obs = NULL,
                                 huss_obs = NULL,
                                 idx = NULL) {
  ## `idx` (optional): caller-precomputed .tpf_kron_idx() lists, keyed
  ## $xx/$ex/$ee (pairwise) and $e_xx/$e_ex/$x_xx/$e_ee (triple, inner-first
  ## grouped) -- period-constant, so hot callers (tpf_run_period3's mutation
  ## micro-calls) pass these once instead of rebuilding them on every N=1
  ## call. Omitted -> computed inline (bit-identical either way).
  n_s <- ncol(ZZ)
  N   <- ncol(particles)
  n_e <- ncol(DD)

  x1 <- particles[seq_len(n_s), , drop = FALSE]
  x2 <- particles[seq_len(n_s) + n_s, , drop = FALSE]
  x3 <- particles[seq_len(n_s) + 2L * n_s, , drop = FALSE]

  ## NOTE (nesting exactness): (x1 + x2 + x3) evaluates as ((x1+x2)+x3), so
  ## when x3 == 0 identically this is bit-identical to the order-2 weights.
  mu <- ZZ %*% (x1 + x2 + x3)
  if (any(ZZ_xss != 0)) mu <- mu + ZZ_xss %*% x1

  if (!is.null(shocks_t)) {
    mu <- mu + DD %*% shocks_t
  } else if (!is.null(DD_shock_mean)) {
    mu <- mu + DD_shock_mean
  }

  ## Nonlinear obs terms (2026-08-04 obs-tensor fix): EXACT transcription of
  ## simulate_model_order3's y2/y3 observable reconstruction
  ## (solve-perturbation-order3.R:1362-1388), using the SAME (x1_prev,
  ## x2_prev, e) arguments the state recursion uses -- `particles` here ARE
  ## the pre-propagation x1_prev/x2_prev/x3_prev layers and `shocks_t` IS
  ## e_t.  Each term is guarded by the existing `if (any(tensor != 0))`
  ## conditional-add pattern (mirrors .tpf_propagate3_R) so a model with all
  ## higher tensors zero reproduces the pre-fix nesting pin bit-for-bit.
  have_e <- !is.null(shocks_t)

  ## -- y2 pieces: 0.5*hxx(x1(x)x1), hxu(e(x)x1), 0.5*huu(e(x)e) -----------
  if (!is.null(hxx_obs) && any(hxx_obs != 0)) {
    idx_xx  <- if (!is.null(idx)) idx$xx else .tpf_kron_idx(n_s, n_s)
    kron_xx <- .tpf_kron2(x1, x1, idx_xx)
    mu <- mu + 0.5 * hxx_obs %*% kron_xx
  }
  if (have_e && !is.null(hxu_obs) && any(hxu_obs != 0)) {
    idx_ex  <- if (!is.null(idx)) idx$ex else .tpf_kron_idx(n_e, n_s)
    kron_ex <- .tpf_kron2(shocks_t, x1, idx_ex)
    mu <- mu + hxu_obs %*% kron_ex
  }
  if (have_e && !is.null(huu_obs) && any(huu_obs != 0)) {
    idx_ee  <- if (!is.null(idx)) idx$ee else .tpf_kron_idx(n_e, n_e)
    kron_ee <- .tpf_kron2(shocks_t, shocks_t, idx_ee)
    mu <- mu + 0.5 * huu_obs %*% kron_ee
  }

  ## -- y3 pieces: hxx(x1(x)x2), hxu(e(x)x2), 0.5*hxxu(e(x)x1(x)x1),
  ##    0.5*hxuu(e(x)e(x)x1), (1/6)*hxxx(x1^3), (1/6)*huuu(e^3),
  ##    0.5*huss(e) [the term that was missing before this fix] ------------
  if (!is.null(hxx_obs) && any(hxx_obs != 0)) {
    idx_xx2  <- if (!is.null(idx)) idx$xx else .tpf_kron_idx(n_s, n_s)
    kron_x12 <- .tpf_kron2(x1, x2, idx_xx2)
    mu <- mu + hxx_obs %*% kron_x12
  }
  if (have_e && !is.null(hxu_obs) && any(hxu_obs != 0)) {
    idx_ex2  <- if (!is.null(idx)) idx$ex else .tpf_kron_idx(n_e, n_s)
    kron_ex2 <- .tpf_kron2(shocks_t, x2, idx_ex2)
    mu <- mu + hxu_obs %*% kron_ex2
  }
  if (have_e && !is.null(hxxu_obs) && any(hxxu_obs != 0)) {
    ## kronecker(e, kronecker(x1, x1)) -- grouped inner-first
    idx_xx_i  <- if (!is.null(idx)) idx$xx   else .tpf_kron_idx(n_s, n_s)
    idx_e_xx  <- if (!is.null(idx)) idx$e_xx else .tpf_kron_idx(n_e, n_s * n_s)
    kron_xx_i <- .tpf_kron2(x1, x1, idx_xx_i)
    kron_exx  <- .tpf_kron2(shocks_t, kron_xx_i, idx_e_xx)
    mu <- mu + 0.5 * hxxu_obs %*% kron_exx
  }
  if (have_e && !is.null(hxuu_obs) && any(hxuu_obs != 0)) {
    ## kronecker(e, kronecker(e, x1)) -- grouped inner-first
    idx_ex_i  <- if (!is.null(idx)) idx$ex   else .tpf_kron_idx(n_e, n_s)
    idx_e_ex  <- if (!is.null(idx)) idx$e_ex else .tpf_kron_idx(n_e, n_e * n_s)
    kron_ex_i <- .tpf_kron2(shocks_t, x1, idx_ex_i)
    kron_eex  <- .tpf_kron2(shocks_t, kron_ex_i, idx_e_ex)
    mu <- mu + 0.5 * hxuu_obs %*% kron_eex
  }
  if (!is.null(hxxx_obs) && any(hxxx_obs != 0)) {
    ## kronecker(x1, kronecker(x1, x1)) -- grouped inner-first
    idx_xx_i2  <- if (!is.null(idx)) idx$xx   else .tpf_kron_idx(n_s, n_s)
    idx_x_xx   <- if (!is.null(idx)) idx$x_xx else .tpf_kron_idx(n_s, n_s * n_s)
    kron_xx_i2 <- .tpf_kron2(x1, x1, idx_xx_i2)
    kron_xxx   <- .tpf_kron2(x1, kron_xx_i2, idx_x_xx)
    mu <- mu + (1 / 6) * hxxx_obs %*% kron_xxx
  }
  if (have_e && !is.null(huuu_obs) && any(huuu_obs != 0)) {
    ## kronecker(e, kronecker(e, e)) -- grouped inner-first
    idx_ee_i <- if (!is.null(idx)) idx$ee   else .tpf_kron_idx(n_e, n_e)
    idx_e_ee <- if (!is.null(idx)) idx$e_ee else .tpf_kron_idx(n_e, n_e * n_e)
    kron_ee_i <- .tpf_kron2(shocks_t, shocks_t, idx_ee_i)
    kron_eee  <- .tpf_kron2(shocks_t, kron_ee_i, idx_e_ee)
    mu <- mu + (1 / 6) * huuu_obs %*% kron_eee
  }
  if (have_e && !is.null(huss_obs) && any(huss_obs != 0)) {
    mu <- mu + 0.5 * huss_obs %*% shocks_t
  }

  mu <- mu + d_obs + ghss_obs
  if (any(ghs3_obs != 0)) mu <- mu + ghs3_obs

  sd_phi <- sqrt(me_variance / phi)
  lw <- dnorm(as.numeric(y_t), mean = mu, sd = sd_phi, log = TRUE)
  dim(lw) <- dim(mu)
  colSums(lw)
}


## ---- Per-period TPF filter -----------------------------------------------

#' Run one period of the Tempered Particle Filter
#'
#' Implements H&S (2019) Algorithm 1 within-period tempering:
#' given period t-1 particles, propagate via the bootstrap proposal
#' (state-transition density), then adaptively temper phi from 0 to 1,
#' resampling and mutating as needed.
#'
#' Returns the period-t particles and the log-marginal-likelihood
#' contribution log p(y_t | Y_{1:t-1}).
#'
#' @param particles   Matrix (2*n_s) x N (period t-1 particles = s_{t-1}).
#' @param y_t         Numeric vector n_obs (period t observation).
#' @param dr2         DecisionRules2 object.
#' @param Sigma_e     Shock covariance (n_e x n_e).
#' @param L_e         Lower-triangular Cholesky factor of Sigma_e.
#' @param ZZ          Observation loading matrix (n_obs x n_s).
#' @param DD          Shock-to-obs matrix (n_obs x n_e) = ghu[obs_idx,].
#' @param d_obs       Obs steady-state mean (n_obs vector).
#' @param ghss_obs    0.5 * ghss[obs_idx] (n_obs vector).
#' @param hxx_obs     n_obs x n_s^2 obs-row slice of ghxx (NULL/all-zero ok).
#' @param hxu_obs     n_obs x (n_e*n_s) obs-row slice of ghxu (NULL/all-zero ok).
#' @param huu_obs     n_obs x n_e^2 obs-row slice of ghuu (NULL/all-zero ok).
#' @param me_variance Scalar > 0.
#' @param ess_target  ESS target fraction (default 0.5).
#' @param n_mh        Herbst & Schorfheide (2019) mutation steps per phi
#'   stage (default 1). Each step proposes a fresh period-t shock
#'   \code{e' ~ N(0, Sigma_e)} (an independence proposal) with the
#'   ancestor state \code{s_{t-1}} held FIXED, and accepts with ratio
#'   \code{phi_curr * (loglik(e') - loglik(e))}. The state is never moved
#'   by mutation.
#' @param mh_scale    Unused since the mutation step was corrected to fix
#'   the ancestor state (no state random walk is proposed, so no proposal
#'   covariance scale is needed). Kept for API/back-compat only.
#' @param use_rcpp    Logical; use C++ propagation if TRUE.
#' @param backend     Character; \code{"cpp"} (default, uses the C++ period
#'   loop when compiled) or \code{"R"} (pure-R reference path).  When
#'   \code{"cpp"} is requested but the compiled symbol is absent, falls back
#'   to \code{"R"} silently.
#' @return list(particles = (2*n_s) x N matrix, log_lik_contrib = scalar).
#' @noRd
tpf_run_period <- function(particles, y_t, dr2, Sigma_e, L_e,
                            ZZ, DD, d_obs, ghss_obs,
                            me_variance, ess_target = 0.5,
                            n_mh = 1L, mh_scale = 1.0,
                            use_rcpp = .HAS_RCPP_TPF(),
                            backend  = if (.HAS_RCPP_TPF_PERIOD()) "cpp" else "R",
                            hxx_obs = NULL, hxu_obs = NULL, huu_obs = NULL,
                            U_normals  = NULL,   # CPM: n_e x N standard normals, or NULL
                            U_resample = NULL,   # CPM: scalar z ~ N(0,1) for phi=1 uniform
                                                 #      (u = pnorm(z)); NULL -> draw from RNG
                            U_mid      = NULL,   # CPM legacy: K-vector z_k ~ N(0,1) for mid-stage
                                                 #      resampling uniforms; never fires in practice
                            U_mutation = NULL) { # CPM Option A: pre-drawn mutation noise.
                                                 #   (n_2s+n_e+1) x (max_stages_u*n_mh*N) matrix.
                                                 #   NULL -> draw from RNG (bit-identical to pre-Option-A)

  N    <- ncol(particles)
  n_s  <- ncol(ZZ)
  n_e  <- nrow(Sigma_e)
  n_2s <- 2L * n_s
  n_obs <- nrow(ZZ)

  state_idx <- dr2$state_idx

  ## Extract state-row slices
  hx  <- dr2$ghx[state_idx, , drop = FALSE]
  hu  <- dr2$ghu[state_idx, , drop = FALSE]
  hxx <- dr2$ghxx[state_idx, , drop = FALSE]
  hxu <- dr2$ghxu[state_idx, , drop = FALSE]
  huu <- dr2$ghuu[state_idx, , drop = FALSE]
  hss <- dr2$ghss[state_idx]

  ## ---- C++ per-period loop dispatch ----------------------------------------
  if (identical(backend, "cpp") && .HAS_RCPP_TPF_PERIOD()) {
    res_cpp <- tpf_run_period_cpp(
      particles   = particles,
      y_t         = as.numeric(y_t),
      L_e         = L_e,
      hx          = hx,
      hu          = hu,
      hxx         = if (is.null(hxx)) matrix(0, n_s, n_s * n_s) else hxx,
      hxu         = if (is.null(hxu)) matrix(0, n_s, n_s * n_e) else hxu,
      huu         = if (is.null(huu)) matrix(0, n_s, n_e * n_e) else huu,
      hss         = if (is.null(hss)) numeric(n_s) else as.numeric(hss),
      ZZ          = ZZ,
      DD          = DD,
      d_obs       = if (is.null(d_obs)) numeric(nrow(ZZ)) else as.numeric(d_obs),
      ghss_obs    = if (is.null(ghss_obs)) numeric(nrow(ZZ)) else as.numeric(ghss_obs),
      hxx_obs     = if (is.null(hxx_obs)) matrix(0, n_obs, n_s * n_s) else hxx_obs,
      hxu_obs     = if (is.null(hxu_obs)) matrix(0, n_obs, n_s * n_e) else hxu_obs,
      huu_obs     = if (is.null(huu_obs)) matrix(0, n_obs, n_e * n_e) else huu_obs,
      me_variance = me_variance,
      ess_target  = ess_target,
      n_mh        = as.integer(n_mh),
      mh_scale    = mh_scale,
      max_stages  = 200L,
      U_normals   = U_normals,  # CPM: NULL -> draw internally (bit-identical to pre-CPM)
      U_resample  = U_resample, # CPM: NULL -> draw from RNG (bit-identical to pre-CPM)
      U_mid       = U_mid,      # CPM legacy: NULL -> draw from RNG for mid-stage (bit-identical)
      U_mutation  = U_mutation  # CPM Option A: NULL -> draw from RNG (bit-identical to pre-Option-A)
    )
    ## z_resample_used: scalar (NaN when non-CPM, actual z when CPM)
    z_res <- res_cpp$z_resample_used
    if (!is.null(z_res) && length(z_res) == 1L && is.nan(z_res))
      z_res <- NA_real_
    return(list(particles        = res_cpp$particles,
                log_lik_contrib  = res_cpp$log_lik_contrib,
                U_used           = res_cpp$U_used,
                z_resample_used  = z_res,
                u_mid_slots_used = res_cpp$u_mid_slots_used))
  }

  ## ---- Pure-R fallback (reference path) ------------------------------------

  ## Hoist the "any(tensor != 0)" obs-tensor scans OUT of the per-call weight
  ## kernel: these tensors are period-constant (they don't depend on the
  ## particle/shock draws), but .tpf_log_weights_R is called on every
  ## mutation micro-call (2x per particle per MH step per stage), and its
  ## `!is.null(x) && any(x != 0)` guard would otherwise re-scan the full
  ## obs-tensor matrix every time. Passing NULL for an all-zero tensor hits
  ## that SAME guard identically (any(NULL-branch) short-circuits the same
  ## way any(x != 0) == FALSE would), so this is bit-identical -- just
  ## computed once instead of on every call.
  hxx_obs_h <- if (!is.null(hxx_obs) && any(hxx_obs != 0)) hxx_obs else NULL
  hxu_obs_h <- if (!is.null(hxu_obs) && any(hxu_obs != 0)) hxu_obs else NULL
  huu_obs_h <- if (!is.null(huu_obs) && any(huu_obs != 0)) huu_obs else NULL

  ## Precompute the Kronecker row-index vectors ONCE per period (they depend
  ## only on n_s/n_e, not on the particle/shock draws) and pass them down to
  ## every .tpf_log_weights_R call, including the N=1 mutation micro-calls
  ## (the dominant remaining cost after the kron vectorization: rebuilding
  ## rep()-index vectors on every single-column call).
  idx_kw <- list(xx = .tpf_kron_idx(n_s, n_s),
                 ex = .tpf_kron_idx(n_e, n_s),
                 ee = .tpf_kron_idx(n_e, n_e))

  ## -- Step 1: draw shocks and propagate the state -------------------------
  ## State-dating convention (matches Dynare / dynhr KF):
  ##   y_t = ZZ * (x1_{t-1} + x2_{t-1}) + DD * e_t + d_obs + ghss_obs
  ##   s_t = transition(s_{t-1}, e_t)
  ## So we compute observation likelihoods from the PRE-propagation particles
  ## plus the drawn shocks, THEN propagate to get the period-t particles.
  ## CPM: use pre-drawn standard normals when supplied; otherwise draw fresh.
  z_mat_R <- if (!is.null(U_normals)) U_normals else matrix(rnorm(n_e * N), nrow = n_e)
  shocks <- L_e %*% z_mat_R

  ## Full measurement log-likelihoods (phi=1) for CURRENT draw of shocks.
  ## Uses PRE-propagation particles (= s_{t-1}) and current shocks (= e_t).
  log_liks <- .tpf_log_weights_R(particles, y_t, ZZ, DD, shocks,
                                  d_obs, ghss_obs, me_variance, phi = 1.0,
                                  hxx_obs = hxx_obs_h, hxu_obs = hxu_obs_h,
                                  huu_obs = huu_obs_h, idx = idx_kw)

  ## -- Step 2: adaptive phi tempering loop (H&S 2019 Algorithm 1) ---------
  phi_curr          <- 0
  log_w             <- rep(0, N)   # log importance weights (uniform initially)
  z_resample_used_R <- NA_real_    # CPM: z used for phi=1 uniform (NA = non-CPM)
  log_lik_contrib <- 0           # accumulator for log p(y_t | Y_{1:t-1})
  u_mid_idx_R     <- 1L          # CPM: next U_mid slot index to consume (1-based in R)

  max_stages <- 200L
  for (stage in seq_len(max_stages)) {

    ## Find next phi by bisection: ESS(delta_phi * log_liks) = ess_target * N
    phi_next <- .smc_next_lambda(log_liks, phi_curr, ess_target, N)

    ## Incremental weight update
    delta_phi <- phi_next - phi_curr
    inc_log_w <- delta_phi * log_liks
    log_w_new <- log_w + inc_log_w

    ## Accumulate log normalisation constant:
    ## log p(y_t | Y_{1:t-1}) += log(sum exp(log_w_new)) - log(sum exp(log_w))
    log_lik_contrib <- log_lik_contrib +
      (.smc_log_sum_exp(log_w_new) - .smc_log_sum_exp(log_w))

    log_w    <- log_w_new
    phi_curr <- phi_next

    ## Normalise
    log_w_c <- log_w - max(log_w)
    w_norm  <- exp(log_w_c)
    w_norm  <- w_norm / sum(w_norm)

    ess <- .smc_ess(log_w)

    if (phi_curr >= 1 - 1e-10) {
      ## Reached phi=1: ALWAYS resample before returning. The caller treats
      ## the returned particles as an equally-weighted draw from the period-t
      ## filtering distribution; returning them unweighted while log_w is
      ## non-uniform silently drops the final-stage weights and biases every
      ## subsequent period's likelihood contribution (a bias that does NOT
      ## shrink with N).
      ##
      ## CPM sorted path: when U_resample is supplied, sort particles by first
      ## state dimension (x1[1]) and use the pre-drawn uniform u = pnorm(z).
      ## Non-CPM path: standard unsorted systematic resample (bit-identical).
      if (!is.null(U_resample)) {
        sort_ord  <- order(particles[1L, ])
        w_sorted  <- w_norm[sort_ord]
        u0        <- max(1e-15, min(1 - 1e-15, pnorm(U_resample)))
        cw_sorted <- cumsum(w_sorted)
        cw_sorted[N] <- 1.0
        idx_sorted <- integer(N)
        j <- 1L
        for (ii in seq_len(N)) {
          u_i <- (ii - 1L + u0) / N
          while (j < N && cw_sorted[j] < u_i) j <- j + 1L
          idx_sorted[ii] <- j
        }
        idx <- sort_ord[idx_sorted]
        z_resample_used_R <- U_resample
      } else {
        idx <- .smc_systematic_resample(w_norm, N)
        z_resample_used_R <- NA_real_
      }
      particles <- particles[, idx, drop = FALSE]
      shocks    <- shocks[, idx, drop = FALSE]
      log_liks  <- log_liks[idx]
      break
    }

    ## -- Resample if ESS below threshold -----------------------------------
    resampled <- FALSE
    if (ess < ess_target * N) {
      ## CPM path: use pre-drawn U_mid z slot for sorted systematic resample
      if (!is.null(U_mid) && u_mid_idx_R <= length(U_mid)) {
        z_k       <- U_mid[[u_mid_idx_R]]
        u_mid_idx_R <- u_mid_idx_R + 1L
        u_k       <- max(1e-15, min(1 - 1e-15, pnorm(z_k)))
        sort_ord  <- order(particles[1L, ])
        w_sorted  <- w_norm[sort_ord]
        cw_sorted <- cumsum(w_sorted)
        cw_sorted[N] <- 1.0
        idx_sorted <- integer(N)
        j <- 1L
        for (ii in seq_len(N)) {
          u_i <- (ii - 1L + u_k) / N
          while (j < N && cw_sorted[j] < u_i) j <- j + 1L
          idx_sorted[ii] <- j
        }
        idx <- sort_ord[idx_sorted]
      } else {
        ## Non-CPM path or U_mid exhausted: fresh RNG draw (bit-identical)
        if (!is.null(U_mid) && u_mid_idx_R > length(U_mid)) {
          warning("tpf_run_period R fallback: U_mid slots exhausted (K=",
                  length(U_mid), "); falling back to fresh RNG draw. ",
                  "Increase max_stages_u to suppress.")
        }
        idx <- .smc_systematic_resample(w_norm, N)
      }
      particles <- particles[, idx, drop = FALSE]
      shocks    <- shocks[, idx, drop = FALSE]
      log_liks  <- log_liks[idx]   ## resample log_liks to keep in sync
      log_w     <- rep(0, N)
      w_norm    <- rep(1 / N, N)
      resampled <- TRUE
    }

    ## -- RWMH mutation over (2*n_s + n_e)-dimensional (s_{t-1}, e_t) space -
    ## Target: phi_curr-tempered measurement density p_phi(y_t | s_{t-1}, e_t)
    ## Proposal: random walk on s_{t-1} (cloud covariance) + new e_t draw.
    ## For simplicity we resample the shock e_t independently from N(0,Sigma_e)
    ## at each MH step (this is valid since shocks are independent of the state
    ## and the target distribution over e_t is the mixture of the measurement
    ## likelihood and the prior N(0, Sigma_e)).
    ## Mutation is an MCMC move targeting the phi_curr-tempered density and
    ## is only valid when the particle set is equally weighted -- i.e.
    ## immediately after a resample (H&S 2019 Algorithm 1: resample, THEN
    ## mutate). Mutating a weighted cloud and resetting log_w to uniform
    ## would silently discard the weights.
    ## Herbst & Schorfheide (2019): mutate ONLY the period-t shock e_t with
    ## the ancestor state s_{t-1} FIXED. An independence proposal
    ## e' ~ N(0, Sigma_e) (via L_e below) has acceptance ratio exactly
    ## phi_curr * (loglik(e') - loglik(e)) -- the formula already used
    ## here. Randomly walking the state (the old behaviour) drops the
    ## filtering density p(s_{t-1}|Y_{1:t-1}) that lives only in the
    ## resampled cloud, breaking kernel invariance for the tempered target
    ## and biasing the likelihood estimate upward (increasing with n_mh).
    ## z_s is still drawn/consumed below (RNG/CPM stream compatibility)
    ## but is deliberately unused.
    if (resampled && n_mh > 0L && N > 1L) {
      ## Option A (common random numbers): when a U_mutation buffer is supplied,
      ## consume pre-drawn N(0,1) noise from it with the EXACT C++ layout
      ## (src/tpf_propagate.cpp RWMH block) so the R and C++ paths are
      ## bit-identical under a shared buffer.  Column index (0-based):
      ##   col = stage*(n_mh*N) + step*N + particle
      ## Rows (0-based): [0, n_2s) = z_s; [n_2s, n_2s+n_e) = z_e;
      ##   row n_2s+n_e = z_u ~ N(0,1) transformed to log_u = log(pnorm(z_u)).
      ## R indices are 1-based, so add 1 to columns/rows.  C++ stage is 0-based
      ## (= R `stage` - 1); step/particle likewise.  On buffer exhaustion or
      ## when U_mutation is NULL, fall back to fresh rnorm/runif (the legacy,
      ## production behaviour — bit-identical to the pre-Option-A path).
      has_u_mut     <- !is.null(U_mutation) && is.matrix(U_mutation)
      u_mut_ncols   <- if (has_u_mut) ncol(U_mutation) else 0L
      stage0        <- stage - 1L   # C++ 0-based stage index
      u_mut_warned  <- FALSE        # once-per-period exhaustion warning

      for (i in seq_len(N)) {
        s_i   <- particles[, i]
        e_i   <- shocks[, i]
        e_mat <- matrix(e_i, nrow = n_e, ncol = 1L)
        s_mat <- matrix(s_i, nrow = n_2s, ncol = 1L)
        tlp_i <- phi_curr *
          .tpf_log_weights_R(s_mat, y_t, ZZ, DD, e_mat,
                              d_obs, ghss_obs, me_variance, phi = 1.0,
                              hxx_obs = hxx_obs_h, hxu_obs = hxu_obs_h,
                              huu_obs = huu_obs_h, idx = idx_kw)[1L]

        for (step in seq_len(n_mh)) {
          ## 0-based C++ column index; +1 for R's 1-based column access.
          col_idx0  <- stage0 * (n_mh * N) + (step - 1L) * N + (i - 1L)
          use_u_mut <- has_u_mut && (col_idx0 < u_mut_ncols)

          if (use_u_mut) {
            col_v <- U_mutation[, col_idx0 + 1L]
            z_s   <- col_v[seq_len(n_2s)]
            z_e   <- col_v[n_2s + seq_len(n_e)]
            z_u   <- col_v[n_2s + n_e + 1L]
            log_u <- log(pnorm(z_u))
          } else {
            if (has_u_mut && !u_mut_warned) {
              warning("tpf_run_period R fallback: U_mutation columns exhausted ",
                       "(have ", u_mut_ncols, ", need col ", col_idx0,
                       "); falling back to fresh RNG draws. ",
                       "Increase max_stages_u or n_mh to suppress.")
              u_mut_warned <- TRUE
            }
            z_s   <- rnorm(n_2s)
            z_e   <- rnorm(n_e)
            log_u <- log(runif(1L))
          }

          ## Ancestor state fixed (H&S 2019); z_s drawn above but unused --
          ## keeps the mutation-buffer layout / RNG stream bit-identical.
          s_prop  <- s_i
          e_prop  <- as.numeric(L_e %*% z_e)   # resample shock
          e_mat_p <- matrix(e_prop, nrow = n_e, ncol = 1L)
          s_mat_p <- matrix(s_prop, nrow = n_2s, ncol = 1L)
          tlp_p   <- phi_curr *
            .tpf_log_weights_R(s_mat_p, y_t, ZZ, DD, e_mat_p,
                                d_obs, ghss_obs, me_variance, phi = 1.0,
                                hxx_obs = hxx_obs_h, hxu_obs = hxu_obs_h,
                                huu_obs = huu_obs_h, idx = idx_kw)[1L]
          if (is.finite(tlp_p) && log_u < tlp_p - tlp_i) {
            s_i   <- s_prop
            e_i   <- e_prop
            tlp_i <- tlp_p
          }
        }
        particles[, i] <- s_i
        shocks[, i]    <- e_i
      }

      ## Recompute full log-likelihoods after mutation
      log_liks <- .tpf_log_weights_R(particles, y_t, ZZ, DD, shocks,
                                      d_obs, ghss_obs, me_variance, phi = 1.0,
                                      hxx_obs = hxx_obs_h, hxu_obs = hxu_obs_h,
                                      huu_obs = huu_obs_h, idx = idx_kw)
      ## Reset weights to uniform (temper from phi_curr at next iteration)
      log_w  <- rep(0, N)
    }
  }  # end phi loop

  ## -- Step 3: propagate pre-state particles to period-t state -------------
  ## Use the (possibly resampled/mutated) (s_{t-1}, e_t) pairs to produce s_t.
  if (use_rcpp) {
    ## Guard against NULL second-order matrices (first-order models)
    particles_new <- tpf_propagate_particles(
      particles,
      shocks,
      hx, hu,
      if (is.null(hxx)) matrix(0, n_s, n_s * n_s) else hxx,
      if (is.null(hxu)) matrix(0, n_s, n_s * n_e) else hxu,
      if (is.null(huu)) matrix(0, n_s, n_e * n_e) else huu,
      if (is.null(hss)) numeric(n_s)  else as.numeric(hss))
  } else {
    particles_new <- .tpf_propagate_R(particles, shocks,
                                       hx, hu, hxx, hxu, huu, hss)
  }

  list(particles        = particles_new,
       log_lik_contrib  = log_lik_contrib,
       U_used           = z_mat_R,              # CPM: shock normals used (pre-L_e)
       z_resample_used  = z_resample_used_R,    # CPM: z for phi=1 uniform (NA if non-CPM)
       u_mid_slots_used = u_mid_idx_R - 1L)     # CPM: mid-stage slots consumed
}


#' Run one period of the Tempered Particle Filter on the pruned ORDER-3 state
#'
#' Order-3 sibling of \code{tpf_run_period}.  Pure R only (no C++ kernel in
#' this increment); the tempering/resampling/mutation control flow is a
#' line-for-line mirror of tpf_run_period's pure-R reference path with the
#' state enlarged to (x1, x2, x3) (dimension 3*n_s) and the propagation /
#' weight kernels swapped for their order-3 versions.  The order-2 path is
#' entirely untouched by this function's existence.
#'
#' @param particles (3*n_s) x N period t-1 particles.
#' @param dr3 DecisionRules3 object (state-row slices extracted here).
#' @param ZZ_xss n_obs x n_s extra x1 obs loading (0.5*ghxss at obs rows).
#' @param ghs3_obs n_obs vector ((1/6)*ghs3 at obs rows).
#' @param hxx_obs  n_obs x n_s^2 obs-row slice of ghxx (NULL/all-zero ok).
#' @param hxu_obs  n_obs x (n_e*n_s) obs-row slice of ghxu (NULL/all-zero ok).
#' @param huu_obs  n_obs x n_e^2 obs-row slice of ghuu (NULL/all-zero ok).
#' @param hxxu_obs n_obs x (n_e*n_s^2) obs-row slice of ghxxu (NULL/all-zero ok).
#' @param hxuu_obs n_obs x (n_e^2*n_s) obs-row slice of ghxuu (NULL/all-zero ok).
#' @param hxxx_obs n_obs x n_s^3 obs-row slice of ghxxx (NULL/all-zero ok).
#' @param huuu_obs n_obs x n_e^3 obs-row slice of ghuuu (NULL/all-zero ok).
#' @param huss_obs n_obs x n_e obs-row slice of ghuss (NULL/all-zero ok).
#' @inheritParams tpf_run_period
#' @return list(particles = (3*n_s) x N matrix, log_lik_contrib = scalar,
#'   U_used, z_resample_used, u_mid_slots_used) — same contract as
#'   tpf_run_period.
#' @noRd
tpf_run_period3 <- function(particles, y_t, dr3, Sigma_e, L_e,
                             ZZ, DD, d_obs, ghss_obs, ZZ_xss, ghs3_obs,
                             me_variance, ess_target = 0.5,
                             n_mh = 1L, mh_scale = 1.0,
                             hxx_obs = NULL, hxu_obs = NULL, huu_obs = NULL,
                             hxxu_obs = NULL, hxuu_obs = NULL,
                             hxxx_obs = NULL, huuu_obs = NULL,
                             huss_obs = NULL,
                             U_normals  = NULL,
                             U_resample = NULL,
                             U_mid      = NULL,
                             U_mutation = NULL) {

  N    <- ncol(particles)
  n_s  <- ncol(ZZ)
  n_e  <- nrow(Sigma_e)
  n_3s <- 3L * n_s

  state_idx <- dr3$state_idx

  ## State-row slices (see .tpf_propagate3_R for the recursion reference)
  hx   <- dr3$ghx  [state_idx, , drop = FALSE]
  hu   <- dr3$ghu  [state_idx, , drop = FALSE]
  hxx  <- dr3$ghxx [state_idx, , drop = FALSE]
  hxu  <- dr3$ghxu [state_idx, , drop = FALSE]
  huu  <- dr3$ghuu [state_idx, , drop = FALSE]
  hss  <- dr3$ghss [state_idx]
  hxxx <- dr3$ghxxx[state_idx, , drop = FALSE]
  hxxu <- dr3$ghxxu[state_idx, , drop = FALSE]
  hxuu <- dr3$ghxuu[state_idx, , drop = FALSE]
  huuu <- dr3$ghuuu[state_idx, , drop = FALSE]
  hxss <- if (!is.null(dr3$ghxss)) dr3$ghxss[state_idx, , drop = FALSE] else NULL
  huss <- if (!is.null(dr3$ghuss)) dr3$ghuss[state_idx, , drop = FALSE] else NULL
  hs3  <- if (!is.null(dr3$ghs3))  dr3$ghs3[state_idx]                  else NULL

  ## Hoist the "any(tensor != 0)" obs-tensor scans OUT of the per-call weight
  ## kernel (see tpf_run_period's identical hoist for the full rationale):
  ## these tensors are period-constant, but .tpf_log_weights3_R is called on
  ## every mutation micro-call. Passing NULL for an all-zero tensor hits the
  ## SAME `!is.null(x) && any(x != 0)` guard identically -- bit-identical,
  ## just computed once instead of on every call.
  hxx_obs_h  <- if (!is.null(hxx_obs)  && any(hxx_obs  != 0)) hxx_obs  else NULL
  hxu_obs_h  <- if (!is.null(hxu_obs)  && any(hxu_obs  != 0)) hxu_obs  else NULL
  huu_obs_h  <- if (!is.null(huu_obs)  && any(huu_obs  != 0)) huu_obs  else NULL
  hxxu_obs_h <- if (!is.null(hxxu_obs) && any(hxxu_obs != 0)) hxxu_obs else NULL
  hxuu_obs_h <- if (!is.null(hxuu_obs) && any(hxuu_obs != 0)) hxuu_obs else NULL
  hxxx_obs_h <- if (!is.null(hxxx_obs) && any(hxxx_obs != 0)) hxxx_obs else NULL
  huuu_obs_h <- if (!is.null(huuu_obs) && any(huuu_obs != 0)) huuu_obs else NULL
  huss_obs_h <- if (!is.null(huss_obs) && any(huss_obs != 0)) huss_obs else NULL

  ## Precompute the Kronecker row-index vectors ONCE per period (period-
  ## constant; see tpf_run_period's identical hoist) and pass them down to
  ## every .tpf_log_weights3_R call, including the N=1 mutation micro-calls.
  idx_kw <- list(xx   = .tpf_kron_idx(n_s, n_s),
                 ex   = .tpf_kron_idx(n_e, n_s),
                 ee   = .tpf_kron_idx(n_e, n_e),
                 e_xx = .tpf_kron_idx(n_e, n_s * n_s),
                 e_ex = .tpf_kron_idx(n_e, n_e * n_s),
                 x_xx = .tpf_kron_idx(n_s, n_s * n_s),
                 e_ee = .tpf_kron_idx(n_e, n_e * n_e))

  ## -- Step 1: draw shocks; observation weights from PRE-propagation state --
  z_mat_R <- if (!is.null(U_normals)) U_normals else matrix(rnorm(n_e * N), nrow = n_e)
  shocks <- L_e %*% z_mat_R

  log_liks <- .tpf_log_weights3_R(particles, y_t, ZZ, DD, shocks,
                                   d_obs, ghss_obs, ZZ_xss, ghs3_obs,
                                   me_variance, phi = 1.0,
                                   hxx_obs = hxx_obs_h, hxu_obs = hxu_obs_h,
                                   huu_obs = huu_obs_h, hxxu_obs = hxxu_obs_h,
                                   hxuu_obs = hxuu_obs_h, hxxx_obs = hxxx_obs_h,
                                   huuu_obs = huuu_obs_h, huss_obs = huss_obs_h,
                                   idx = idx_kw)

  ## -- Step 2: adaptive phi tempering loop (mirror of tpf_run_period) ------
  phi_curr          <- 0
  log_w             <- rep(0, N)
  z_resample_used_R <- NA_real_
  log_lik_contrib   <- 0
  u_mid_idx_R       <- 1L

  max_stages <- 200L
  for (stage in seq_len(max_stages)) {

    phi_next  <- .smc_next_lambda(log_liks, phi_curr, ess_target, N)
    delta_phi <- phi_next - phi_curr
    inc_log_w <- delta_phi * log_liks
    log_w_new <- log_w + inc_log_w

    log_lik_contrib <- log_lik_contrib +
      (.smc_log_sum_exp(log_w_new) - .smc_log_sum_exp(log_w))

    log_w    <- log_w_new
    phi_curr <- phi_next

    log_w_c <- log_w - max(log_w)
    w_norm  <- exp(log_w_c)
    w_norm  <- w_norm / sum(w_norm)

    ess <- .smc_ess(log_w)

    if (phi_curr >= 1 - 1e-10) {
      ## phi = 1: ALWAYS resample before returning (see tpf_run_period).
      if (!is.null(U_resample)) {
        sort_ord  <- order(particles[1L, ])
        w_sorted  <- w_norm[sort_ord]
        u0        <- max(1e-15, min(1 - 1e-15, pnorm(U_resample)))
        cw_sorted <- cumsum(w_sorted)
        cw_sorted[N] <- 1.0
        idx_sorted <- integer(N)
        j <- 1L
        for (ii in seq_len(N)) {
          u_i <- (ii - 1L + u0) / N
          while (j < N && cw_sorted[j] < u_i) j <- j + 1L
          idx_sorted[ii] <- j
        }
        idx <- sort_ord[idx_sorted]
        z_resample_used_R <- U_resample
      } else {
        idx <- .smc_systematic_resample(w_norm, N)
        z_resample_used_R <- NA_real_
      }
      particles <- particles[, idx, drop = FALSE]
      shocks    <- shocks[, idx, drop = FALSE]
      log_liks  <- log_liks[idx]
      break
    }

    ## -- Resample if ESS below threshold -----------------------------------
    resampled <- FALSE
    if (ess < ess_target * N) {
      if (!is.null(U_mid) && u_mid_idx_R <= length(U_mid)) {
        z_k         <- U_mid[[u_mid_idx_R]]
        u_mid_idx_R <- u_mid_idx_R + 1L
        u_k         <- max(1e-15, min(1 - 1e-15, pnorm(z_k)))
        sort_ord    <- order(particles[1L, ])
        w_sorted    <- w_norm[sort_ord]
        cw_sorted   <- cumsum(w_sorted)
        cw_sorted[N] <- 1.0
        idx_sorted  <- integer(N)
        j <- 1L
        for (ii in seq_len(N)) {
          u_i <- (ii - 1L + u_k) / N
          while (j < N && cw_sorted[j] < u_i) j <- j + 1L
          idx_sorted[ii] <- j
        }
        idx <- sort_ord[idx_sorted]
      } else {
        if (!is.null(U_mid) && u_mid_idx_R > length(U_mid)) {
          warning("tpf_run_period3: U_mid slots exhausted (K=",
                  length(U_mid), "); falling back to fresh RNG draw.")
        }
        idx <- .smc_systematic_resample(w_norm, N)
      }
      particles <- particles[, idx, drop = FALSE]
      shocks    <- shocks[, idx, drop = FALSE]
      log_liks  <- log_liks[idx]
      log_w     <- rep(0, N)
      w_norm    <- rep(1 / N, N)
      resampled <- TRUE
    }

    ## -- RWMH mutation over (3*n_s + n_e)-dimensional (s_{t-1}, e_t) space -
    ## Same structure as tpf_run_period; only the state dimension changes
    ## (n_3s here vs n_2s there), so the U_mutation buffer layout is
    ## rows [0, n_3s) = z_s; [n_3s, n_3s+n_e) = z_e; row n_3s+n_e = z_u.
    ## Herbst & Schorfheide (2019): mutate ONLY the period-t shock e_t with
    ## the ancestor state s_{t-1} FIXED -- see the order-2 block above for
    ## the full rationale (kernel invariance / upward likelihood bias).
    ## z_s is still drawn/consumed below (RNG/CPM stream compatibility)
    ## but is deliberately unused.
    if (resampled && n_mh > 0L && N > 1L) {
      has_u_mut     <- !is.null(U_mutation) && is.matrix(U_mutation)
      u_mut_ncols   <- if (has_u_mut) ncol(U_mutation) else 0L
      stage0        <- stage - 1L
      u_mut_warned  <- FALSE        # once-per-period exhaustion warning

      for (i in seq_len(N)) {
        s_i   <- particles[, i]
        e_i   <- shocks[, i]
        e_mat <- matrix(e_i, nrow = n_e, ncol = 1L)
        s_mat <- matrix(s_i, nrow = n_3s, ncol = 1L)
        tlp_i <- phi_curr *
          .tpf_log_weights3_R(s_mat, y_t, ZZ, DD, e_mat,
                               d_obs, ghss_obs, ZZ_xss, ghs3_obs,
                               me_variance, phi = 1.0,
                               hxx_obs = hxx_obs_h, hxu_obs = hxu_obs_h,
                               huu_obs = huu_obs_h, hxxu_obs = hxxu_obs_h,
                               hxuu_obs = hxuu_obs_h, hxxx_obs = hxxx_obs_h,
                               huuu_obs = huuu_obs_h, huss_obs = huss_obs_h,
                               idx = idx_kw)[1L]

        for (step in seq_len(n_mh)) {
          col_idx0  <- stage0 * (n_mh * N) + (step - 1L) * N + (i - 1L)
          use_u_mut <- has_u_mut && (col_idx0 < u_mut_ncols)

          if (use_u_mut) {
            col_v <- U_mutation[, col_idx0 + 1L]
            z_s   <- col_v[seq_len(n_3s)]
            z_e   <- col_v[n_3s + seq_len(n_e)]
            z_u   <- col_v[n_3s + n_e + 1L]
            log_u <- log(pnorm(z_u))
          } else {
            if (has_u_mut && !u_mut_warned) {
              warning("tpf_run_period3 R fallback: U_mutation columns exhausted ",
                       "(have ", u_mut_ncols, ", need col ", col_idx0,
                       "); falling back to fresh RNG draws. ",
                       "Increase max_stages_u or n_mh to suppress.")
              u_mut_warned <- TRUE
            }
            z_s   <- rnorm(n_3s)
            z_e   <- rnorm(n_e)
            log_u <- log(runif(1L))
          }

          ## Ancestor state fixed (H&S 2019); z_s drawn above but unused --
          ## keeps the mutation-buffer layout / RNG stream bit-identical.
          s_prop  <- s_i
          e_prop  <- as.numeric(L_e %*% z_e)
          e_mat_p <- matrix(e_prop, nrow = n_e, ncol = 1L)
          s_mat_p <- matrix(s_prop, nrow = n_3s, ncol = 1L)
          tlp_p   <- phi_curr *
            .tpf_log_weights3_R(s_mat_p, y_t, ZZ, DD, e_mat_p,
                                 d_obs, ghss_obs, ZZ_xss, ghs3_obs,
                                 me_variance, phi = 1.0,
                                 hxx_obs = hxx_obs_h, hxu_obs = hxu_obs_h,
                                 huu_obs = huu_obs_h, hxxu_obs = hxxu_obs_h,
                                 hxuu_obs = hxuu_obs_h, hxxx_obs = hxxx_obs_h,
                                 huuu_obs = huuu_obs_h, huss_obs = huss_obs_h,
                                 idx = idx_kw)[1L]
          if (is.finite(tlp_p) && log_u < tlp_p - tlp_i) {
            s_i   <- s_prop
            e_i   <- e_prop
            tlp_i <- tlp_p
          }
        }
        particles[, i] <- s_i
        shocks[, i]    <- e_i
      }

      log_liks <- .tpf_log_weights3_R(particles, y_t, ZZ, DD, shocks,
                                       d_obs, ghss_obs, ZZ_xss, ghs3_obs,
                                       me_variance, phi = 1.0,
                                       hxx_obs = hxx_obs_h, hxu_obs = hxu_obs_h,
                                       huu_obs = huu_obs_h, hxxu_obs = hxxu_obs_h,
                                       hxuu_obs = hxuu_obs_h, hxxx_obs = hxxx_obs_h,
                                       huuu_obs = huuu_obs_h, huss_obs = huss_obs_h,
                                       idx = idx_kw)
      log_w  <- rep(0, N)
    }
  }  # end phi loop

  ## -- Step 3: propagate pre-state particles to period-t state -------------
  particles_new <- .tpf_propagate3_R(particles, shocks,
                                      hx, hu, hxx, hxu, huu, hss,
                                      hxxx, hxxu, hxuu, huuu,
                                      hxss, huss, hs3)

  list(particles        = particles_new,
       log_lik_contrib  = log_lik_contrib,
       U_used           = z_mat_R,
       z_resample_used  = z_resample_used_R,
       u_mid_slots_used = u_mid_idx_R - 1L)
}


## ---- PMCMC preflight --------------------------------------------------------

#' Evaluate TPF loglik variance at a fixed theta (internal helper)
#'
#' Calls \code{log_post_fn(theta)} K times with different RNG seeds to
#' measure the Monte Carlo variance of the particle-filter log-likelihood.
#' The closure MUST have been built with \code{seed = NULL} (or the caller
#' must set.seed() before each evaluation externally); otherwise the PF is
#' deterministic and SD = 0.
#'
#' CRITICAL (Landmine 1): if the closure was built with \code{seed = 42L}
#' (non-NULL), it calls \code{set.seed(42)} at each evaluation, making the
#' PF fully deterministic.  This function always builds a fresh
#' \code{seed = NULL} wrapper around the supplied closure; the external
#' \code{set.seed(k)} drives the randomness.
#'
#' @param log_post_fn  A \code{function(theta)} returning
#'   \code{list(logpost, loglik, logprior)}.  Built with \code{seed = NULL}.
#' @param theta_mode   Named numeric vector: the parameter point at which
#'   to evaluate variance (typically the posterior mode).
#' @param K            Integer; number of independent evaluations (default 30).
#' @param verbose      Print a progress message.
#' @return Named list:
#'   \itemize{
#'     \item \code{sd} -- SD of the K log-likelihood values
#'     \item \code{mean} -- mean of the K log-likelihood values
#'     \item \code{logliks} -- numeric(K) raw values
#'     \item \code{n_needed} -- estimated N for SD < 1
#'       (\eqn{N_{needed} \approx N_{current} \cdot \text{SD}^2})
#'     \item \code{accept_noise_factor} -- \code{exp(-sd^2/2)} (Sherlock et al.)
#'   }
#' @noRd
.tpf_pmcmc_preflight <- function(log_post_fn, theta_mode, K = 30L,
                                  verbose = FALSE) {
  K <- as.integer(K)
  if (K <= 0L) {
    return(list(sd = NA_real_, mean = NA_real_, logliks = numeric(0),
                n_needed = NA_integer_, accept_noise_factor = NA_real_,
                skipped = TRUE))
  }

  if (verbose) message(sprintf("  [TPF preflight] evaluating loglik SD (K = %d)...", K))

  logliks <- vapply(seq_len(K), function(k) {
    set.seed(k)   # vary RNG stream externally; closure must have seed = NULL
    res <- tryCatch(log_post_fn(theta_mode), error = function(e) NULL)
    if (is.null(res) || !is.finite(res$loglik)) NA_real_ else res$loglik
  }, numeric(1L))

  valid_ll <- logliks[is.finite(logliks)]
  if (length(valid_ll) < 2L) {
    warning(".tpf_pmcmc_preflight: fewer than 2 finite loglik evaluations; ",
            "SD cannot be estimated. Check that the model evaluates at theta_mode.",
            call. = FALSE)
    return(list(sd = NA_real_, mean = NA_real_, logliks = logliks,
                n_needed = NA_integer_, accept_noise_factor = NA_real_,
                skipped = FALSE))
  }

  sd_ll   <- sd(valid_ll)
  mean_ll <- mean(valid_ll)
  ## N_needed ~ N_current * sd_ll^2 (SD scales ~ 1/sqrt(N))
  ## The attribute n_particles is not available here; caller should multiply.
  ## We return the scaling factor; caller computes n_needed = n_particles * sd_ll^2.
  list(
    sd                  = sd_ll,
    mean                = mean_ll,
    logliks             = logliks,
    n_needed_factor     = sd_ll^2,   # multiply by current n_particles to get N_needed
    accept_noise_factor = exp(-sd_ll^2 / 2),  # Sherlock et al. 2015
    skipped             = FALSE
  )
}


#' Preflight check for PMCMC / particle-MCMC variance
#'
#' Evaluates the TPF log-likelihood K times with different RNG seeds at
#' \code{theta} (typically the posterior mode) to measure Monte Carlo
#' variance. A standard deviation > 1 (the Dynare/Herbst-Schorfheide
#' threshold) indicates that particle-MCMC acceptance will be dominated by
#' log-likelihood noise rather than posterior geometry, and more particles
#' are needed.
#'
#' @param log_post_fn  A log-posterior function built by
#'   \code{\link{make_log_posterior_tpf}} with \code{seed = NULL}.
#' @param theta        Named numeric vector at which to measure variance.
#' @param K            Number of independent evaluations (default 30).
#' @param n_particles  Current number of particles (for the N-recommendation
#'   message; does not affect computation).
#' @param verbose      Print status and recommendation (default TRUE).
#' @return A named list with fields \code{sd}, \code{mean}, \code{logliks},
#'   \code{n_needed}, \code{accept_noise_factor}.
#' @export
tpf_loglik_sd_preflight <- function(log_post_fn, theta, K = 30L,
                                     n_particles = NA_integer_,
                                     verbose = TRUE) {
  pf <- .tpf_pmcmc_preflight(log_post_fn, theta, K = K, verbose = verbose)
  if (isTRUE(pf$skipped)) return(pf)

  n_needed <- if (!is.na(pf$sd) && !is.na(n_particles) && n_particles > 0) {
    ceiling(n_particles * pf$n_needed_factor)
  } else NA_integer_

  pf$n_needed <- n_needed

  if (verbose && !is.na(pf$sd)) {
    message(sprintf("  TPF loglik SD at mode = %.3f  (K = %d evaluations)", pf$sd, K))
    if (pf$sd > 1) {
      message(sprintf(
        paste0(
          "  WARNING: SD > 1 (Dynare threshold). PMCMC acceptance will be dominated by\n",
          "  loglik noise. Current n_particles = %s; estimated n_particles needed for\n",
          "  SD < 1: ~%s (N_needed ~ N * SD^2)."),
        if (is.na(n_particles)) "unknown" else format(n_particles),
        if (is.na(n_needed))    "unknown" else format(n_needed)))
    }
    message(sprintf("  Acceptance noise factor exp(-SD^2/2) = %.3f",
                    pf$accept_noise_factor))
  }

  pf
}


## ---- Factory: make_log_posterior_tpf ------------------------------------

#' Construct the tempered particle filter log-posterior function
#'
#' Returns a \code{function(theta)} that evaluates the log-posterior using
#' the Herbst & Schorfheide (2019) Tempered Particle Filter on the pruned
#' order-2 state space (default), or on the pruned ORDER-3 (AFVRR 2018)
#' state space when \code{order = 3L}.
#'
#' @param model       dynhr_mod from \code{\link{parse_mod}}.
#' @param data        Observation matrix (n_obs x T); columns are time periods.
#' @param prior_spec  Prior specification from \code{extract_prior_spec}.
#' @param obs_vars    Character vector of observed variable names.
#' @param compiled    dynhr_compiled from \code{\link{compile_model}}.
#' @param me_variance Measurement error variance (scalar, must be > 0).
#'   This is the TPF tempering instrument; me_variance = 0 is not allowed.
#' @param n_particles Number of particles (default 1000).
#' @param ess_target  ESS resampling threshold as fraction of N (default 0.5).
#' @param n_mh        Herbst & Schorfheide (2019) mutation steps per phi
#'   stage (default 1); each step is an independence proposal on the
#'   period-t shock only (ancestor state fixed) -- see the internal
#'   \code{tpf_run_period()} mutation step.
#' @param mh_scale    Unused since the mutation step fix (ancestor state is
#'   never moved); kept for API/back-compat only.
#' @param seed        Integer RNG seed for reproducibility (NULL = no fixed seed).
#' @param system_priors Optional system priors list (see \code{\link{sp_irf}}).
#' @param max_stages_u  Maximum number of tempering stages per observation
#'   (default 16).  The filter raises phi from 0 to 1 in at most this many
#'   steps; a smaller value speeds up the filter at the cost of coarser
#'   tempering.
#' @param order  Integer, \code{2L} (default) or \code{3L}: perturbation order
#'   of the pruned state space the particles live on.  \code{2L} is the
#'   original (C++-accelerated) path and is byte-identical to previous
#'   releases; \code{3L} runs the AFVRR pruned third-order transition
#'   (pure-R kernel, see \code{tpf_run_period3}).  Selected via
#'   \code{pruned_order = 3L} from \code{make_log_posterior} /
#'   \code{estimation_context} when \code{likelihood = "tpf"}.
#' @param burn_in_init Integer >= 0 (default 50L, flipped 2026-08-05 from the
#'   prior default of 0). When positive, the initial particle cloud is
#'   propagated this many periods with fresh shocks before the first
#'   observation, so it starts from (approximately) the FULL pruned
#'   stationary joint instead of the first-order Lyapunov draw with
#'   \code{x2_0 = x3_0 = 0}. The Lyapunov-only init mis-states the initial
#'   LEVEL (the pruned model's stationary x2 mean is the nonzero
#'   risk-adjustment shift), which short samples absorb as extra
#'   persistence — caught by the 2026-08-04 order-3 rank-uniformity SBC as
#'   a rho-specific rank shift after the obs-tensor fix. ~50 is ample for
#'   business-cycle persistence and matches the SBC-certified order-3 TPF
#'   config, hence the new default. Pass \code{burn_in_init = 0} to
#'   reproduce pre-0.9.2 behavior (needed for fixed-seed pins keyed to the
#'   old Lyapunov-only init). Incompatible with an externally supplied
#'   \code{U_list} (CPM): callers that build a TPF closure for
#'   \code{rwmh_cpm}/\code{U_list} use MUST pass \code{burn_in_init = 0L}
#'   explicitly, since the CPM slot layout has no burn-in slots — the
#'   closure errors at the boundary if this is violated.
#' @section Missing data:
#' \code{data} may contain \code{NA}/non-finite entries, handled per the
#' package's usual per-element (univariate-KF-style) convention rather than
#' dropping whole periods:
#' \itemize{
#'   \item A period with SOME but not all observables missing is handled by
#'     subsetting the observation vector (and the corresponding rows of the
#'     loading matrices/obs-tensors) to the observed elements for that
#'     period only, then running the normal tempering step on the reduced
#'     observation.
#'   \item A period with ALL observables missing draws that period's shocks
#'     and propagates the particle cloud one step forward with no
#'     weighting/resampling/mutation (a pure predictive step), contributing
#'     0 to the log-likelihood -- it does NOT glue the periods on either
#'     side of the gap together with a single multi-period transition (the
#'     pre-0.9.2 behavior, which silently evaluated the likelihood of a
#'     different, gap-deleted model).
#'   \item Fully-observed periods are unaffected: identical code path and
#'     bit-identical results to previous releases.
#' }
#' @return A \code{function(theta)} returning
#'   \code{list(logpost, loglik, logprior)}.
#' @export
make_log_posterior_tpf <- function(model, data, prior_spec, obs_vars,
                                    compiled,
                                    me_variance,
                                    n_particles  = 1000L,
                                    ess_target   = 0.5,
                                    n_mh         = 1L,
                                    mh_scale     = 1.0,
                                    seed         = NULL,
                                    system_priors = NULL,
                                    max_stages_u  = 16L,
                                    order         = 2L,
                                    burn_in_init  = 50L) {

  ## Force promises for closure-capture safety (same pattern as cumulant branch)
  force(data); force(prior_spec); force(obs_vars); force(me_variance)
  force(n_particles); force(ess_target); force(n_mh); force(mh_scale)
  force(seed); force(system_priors); force(max_stages_u); force(order)
  force(burn_in_init)

  if (!is.numeric(burn_in_init) || length(burn_in_init) != 1L ||
      !is.finite(burn_in_init) || burn_in_init < 0 ||
      burn_in_init != as.integer(burn_in_init)) {
    stop("make_log_posterior_tpf: `burn_in_init` must be a non-negative ",
         "integer.", call. = FALSE)
  }
  burn_in_init <- as.integer(burn_in_init)

  if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
      !(as.integer(order) %in% c(2L, 3L)) || order != as.integer(order)) {
    stop("make_log_posterior_tpf: `order` must be 2 or 3.", call. = FALSE)
  }
  order  <- as.integer(order)
  order3 <- (order == 3L)
  ## Note: U_list is an argument to the INNER function (not captured here),
  ## which avoids rebuilding the outer closure on every CPM MCMC step.

  ## --- Hard stop: me_variance = 0 is degenerate (Landmine 1) --------------
  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance <= 0) {
    stop(
      "make_log_posterior_tpf: 'me_variance' must be a finite positive scalar.\n",
      "The TPF uses measurement error variance as its tempering instrument:\n",
      "effective variance at tempering stage phi is me_variance / phi.\n",
      "me_variance = 0 (or negative / non-finite) makes the filter degenerate.\n",
      "Choose me_variance > 0 (e.g. 1e-4 * var(data))."
    )
  }

  ## Warn on very small me_variance (weight-collapse risk, Landmine 3)
  if (me_variance < 1e-6) {
    warning(
      "make_log_posterior_tpf: me_variance = ", me_variance, " is very small. ",
      "The particle filter may suffer weight collapse (ESS -> 1). ",
      "Consider me_variance >= 1e-4 or at least 1% of data variance."
    )
  }

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  sys_cache <- cache_system_structure(compiled)

  ## data: n_obs x T (columns = time periods)
  T_obs <- ncol(data)

  ## Per-closure warm-start for steady-state solve
  ss_warm <- NULL

  function(theta, U_list = NULL) {
    if (!is.null(U_list) && burn_in_init > 0L) {
      stop("make_log_posterior_tpf: burn_in_init > 0 is incompatible with an ",
           "externally supplied U_list (the CPM slot layout has no burn-in ",
           "slots). Use burn_in_init = 0 with CPM/common-random-number paths.",
           call. = FALSE)
    }
    ## Fixed seed for reproducible filter evaluation.
    ## When seed is set: applies to ALL random draws (shock normals when
    ## U_list = NULL, plus resampling uniforms and RWMH mutation normals).
    ## When U_list is supplied: seed controls the non-U randomness (resampling
    ## uniforms); the shock normals come from U_list (seed has no effect on them).
    if (!is.null(seed)) set.seed(seed)

    lp <- log_prior(theta, prior_spec)
    if (!is.finite(lp))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    params <- .apply_theta_to_params(model, theta)

    ## ---- Steady state (warm-started) ------------------------------------
    ss_result <- solve_steady_state(model, compiled, params,
                                    y0 = ss_warm, verbose = FALSE)
    if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
      if (!is.null(ss_warm))
        ss_result <- solve_steady_state(model, compiled, params,
                                        verbose = FALSE)
      if (is.null(ss_result) || !isTRUE(ss_result$converged)) {
        ss_warm <<- NULL
        return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
      }
    }
    ss_warm <<- ss_result$ss

    ## ---- First-order perturbation ---------------------------------------
    ## Re-derive SSM-computed params for a consistent linearization point
    ## (no-op for non-SSM-parameter models; Tier 13 #1).
    params <- ss_result$params %||% params
    sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)
    dr1 <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
    if (is.null(dr1) || !isTRUE(dr1$bk_satisfied))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## Stationarity guard (mirrors R/posterior.R:150-158)
    ns <- length(dr1$state_idx)
    ev <- dr1$eigenvalues
    spectral_radius <- if (!is.null(ev) && length(ev) >= ns)
      max(Mod(ev[seq_len(ns)]))
    else
      max(Mod(eigen(dr1$ghx[dr1$state_idx, , drop = FALSE],
                    only.values = TRUE)$values))
    if (spectral_radius >= 1)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## ---- Higher-order perturbation (per-draw, expensive) ----------------
    Sigma_e <- .get_shock_cov(model, model$varexo_names, params)
    if (order3) {
      ## Pruned ORDER-3 state space: full third-order solve (includes the
      ## sigma-cross terms ghxss/ghuss via solve_perturbation's order-3 path).
      dr2 <- tryCatch(
        solve_perturbation(model, compiled, ss_result$ss, params,
                           order = 3L, Sigma_e = Sigma_e, verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(dr2) || !isTRUE(dr2$bk_satisfied))
        return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
    } else {
      dr2 <- tryCatch(
        solve_perturbation_order2(model, compiled, ss_result$ss, params,
                                   dr1, Sigma_e = Sigma_e, verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(dr2))
        return(list(logpost = -Inf, loglik = -Inf, logprior = lp))
    }

    ## ---- Observation mapping --------------------------------------------
    endo    <- dr2$endo_names
    obs_idx <- match(obs_vars, endo)
    if (any(is.na(obs_idx)))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ZZ       <- dr2$ghx[obs_idx, , drop = FALSE]   # n_obs x n_s
    DD       <- dr2$ghu[obs_idx, , drop = FALSE]   # n_obs x n_e
    d_obs    <- dr2$ys[obs_vars]                    # n_obs
    ## 0.5 * ghss at obs variables (Landmine 5)
    ghss_obs <- 0.5 * dr2$ghss[obs_idx]             # n_obs

    ## Nonlinear obs-row tensors (2026-08-04 obs-tensor fix): the pruned-model
    ## observation MEAN carries the same quadratic/cubic terms in (x1_prev,
    ## x2_prev, e_t) that the state recursion carries (simulate_model_order3's
    ## y2/y3 reconstruction, solve-perturbation-order3.R:1362-1388) -- these
    ## were previously DROPPED, making the TPF the likelihood of the wrong
    ## (linear-in-obs) model whenever the observed variables have nonzero
    ## quadratic/cubic obs tensors (SBC-caught).  Extract the obs_idx rows
    ## here (order-2 pieces always; order-3 pieces only when order3) and pass
    ## them down to tpf_run_period[3] so ALL THREE weight-evaluation call
    ## sites (initial weights, tempering re-weights, MH-mutation acceptance)
    ## use the identical observation mean.  NULL when the tensor itself is
    ## NULL (no Hessian/third-order arrays computed); the weight kernels
    ## treat NULL as an all-zero tensor and skip the term entirely.
    hxx_obs  <- if (!is.null(dr2$ghxx)) dr2$ghxx[obs_idx, , drop = FALSE] else NULL
    hxu_obs  <- if (!is.null(dr2$ghxu)) dr2$ghxu[obs_idx, , drop = FALSE] else NULL
    huu_obs  <- if (!is.null(dr2$ghuu)) dr2$ghuu[obs_idx, , drop = FALSE] else NULL
    hxxu_obs <- NULL; hxuu_obs <- NULL; hxxx_obs <- NULL
    huuu_obs <- NULL; huss_obs <- NULL
    if (order3) {
      hxxu_obs <- if (!is.null(dr2$ghxxu)) dr2$ghxxu[obs_idx, , drop = FALSE] else NULL
      hxuu_obs <- if (!is.null(dr2$ghxuu)) dr2$ghxuu[obs_idx, , drop = FALSE] else NULL
      hxxx_obs <- if (!is.null(dr2$ghxxx)) dr2$ghxxx[obs_idx, , drop = FALSE] else NULL
      huuu_obs <- if (!is.null(dr2$ghuuu)) dr2$ghuuu[obs_idx, , drop = FALSE] else NULL
      huss_obs <- if (!is.null(dr2$ghuss)) dr2$ghuss[obs_idx, , drop = FALSE] else NULL
    }
    ## Order-3 obs-side pieces, mirroring the pruned-KF3 treatment
    ## (R/pruned-state-space-order3.R:700 and :719/:1097): 0.5*ghxss folds
    ## into the x1 loading; (1/6)*ghs3 is the third-order obs constant.
    ## Both are exactly zero matrices/vectors when the model carries no
    ## sigma-cross terms, so they never perturb the nesting case.
    ZZ_xss <- NULL
    ghs3_obs <- NULL
    if (order3) {
      ZZ_xss <- if (!is.null(dr2$ghxss))
        0.5 * dr2$ghxss[obs_idx, , drop = FALSE]
      else matrix(0, nrow = length(obs_idx), ncol = ncol(ZZ))
      ghs3_obs <- if (!is.null(dr2$ghs3))
        (1 / 6) * dr2$ghs3[obs_idx]
      else numeric(length(obs_idx))
    }
    ## Tag observation matrices as a dsge_ss at the TPF boundary.
    ## TT_s/RR_s (state transition) are extracted below once n_s is known.
    ## ghss_obs has no home in dsge_ss (second-order term) — kept separate.
    ## Hot paths (tpf_run_period, C++ kernel) continue to receive bare matrices.
    TT_s <- dr2$ghx[dr2$state_idx, , drop = FALSE]  # n_s x n_s (needed below)
    RR_s <- dr2$ghu[dr2$state_idx, , drop = FALSE]  # n_s x n_e (needed below)
    ss_obs <- new_dsge_ss(
      T_mat   = TT_s,
      R_mat   = RR_s,
      Z_mat   = ZZ,
      D_mat   = DD,
      Sigma_e = Sigma_e,
      d       = d_obs,
      timing  = "lagged"
    )

    ## Cholesky of Sigma_e for shock sampling
    L_e <- tryCatch(t(chol(Sigma_e)), error = function(e) {
      eg   <- eigen(Sigma_e, symmetric = TRUE)
      vals <- pmax(eg$values, 0)
      eg$vectors %*% diag(sqrt(vals), nrow = length(vals))
    })

    n_s  <- length(dr2$state_idx)
    n_2s <- 2L * n_s
    ## Total particle state dimension: 2*n_s at order 2 (x1,x2), 3*n_s at
    ## order 3 (x1,x2,x3).  n_st == n_2s on the order-2 path, so every use
    ## below is value-identical to the pre-order-3 code there.
    n_st <- if (order3) 3L * n_s else n_2s
    n_e  <- nrow(Sigma_e)
    N    <- as.integer(n_particles)

    ## ---- Initialise particles from stationary distribution ---------------
    ## The Kalman filter initialises P_0 from the stationary Lyapunov equation
    ## P = TT * P * TT' + RR * Sigma_e * RR'.
    ## For comparability, draw x1_0^i ~ N(0, P0) (zero mean at SS deviations);
    ## x2_0 = 0 (no second-order prior correction at time 0).
    ## TT_s / RR_s already extracted above when building ss_obs.
    QQ_s <- tcrossprod(RR_s %*% Sigma_e, RR_s)       # n_s x n_s

    P0 <- tryCatch({
      ## Discrete Lyapunov equation: P = TT*P*TT' + QQ
      ## Use the dlyap-style iteration; Matrix::lyap is not always available.
      ## For small systems the direct solve is fine; 100 iterations converge
      ## quickly for stationary models (spectral radius < 1 guaranteed here).
      Pk <- QQ_s
      for (iter in seq_len(500L)) {
        Pk_new <- TT_s %*% Pk %*% t(TT_s) + QQ_s
        if (max(abs(Pk_new - Pk)) < 1e-12 * (1 + max(abs(Pk_new)))) break
        Pk <- Pk_new
      }
      Pk
    }, error = function(e) NULL)

    ## CPM U_list layout (3T+1 entries):
    ##   Slot 1:              init normals (n_s x N matrix)
    ##   Slots 2..(T+1):      per-period shock normals (n_e x N matrix each)
    ##   Slots (T+2)..(2T+1): per-period phi=1 resampling z scalars (scalar each,
    ##                         z ~ N(0,1); transformed to u = pnorm(z))
    ##   Slots (2T+2)..(3T+1): Option A mutation noise matrices per period.
    ##     When n_mh > 0: (n_2s+n_e+1) x (max_stages_u*n_mh*N) matrix of N(0,1) draws.
    ##     Column layout: stage*(n_mh*N) + step*N + particle (0-based).
    ##     Row layout: [0,n_2s) = z_s, [n_2s,n_2s+n_e) = z_e, n_2s+n_e = z_u (for log_u).
    ##     When n_mh == 0: NULL (mutation never fires, mutation buffer unused).
    ##   (Legacy: was K-vector z_k for mid-stage resamples; repurposed since
    ##    mid-stage resamples never fire in practice (u_mid_slots_used==0).)
    ## When U_list is NULL: all randomness drawn internally (standard path).
    ## Backward compatibility: if length(U_list) == 2T+1 (legacy Tier 9 layout),
    ## mutation buffer slots are treated as NULL (no CRN for mutation, fresh draws).
    U_init_normals <- if (!is.null(U_list)) U_list[[1L]] else NULL

    ## Mutation buffer dimensions (Option A). n_2s and n_e are known here.
    ## mut_buf_rows: number of rows per buffer column (z_s + z_e + 1 for z_u)
    ## mut_buf_cols: number of columns (max_stages_u * n_mh * N)
    ## When n_mh=0: no mutation buffer needed; slots will be NULL.
    ## (n_st = n_2s at order 2 — value-identical to the pre-order-3 code;
    ##  n_st = 3*n_s at order 3, whose mutation z_s spans (x1,x2,x3).)
    mut_buf_rows <- n_st + n_e + 1L
    mut_buf_cols <- as.integer(max_stages_u) * as.integer(n_mh) * N
    has_mutation <- (n_mh > 0L)

    if (is.null(P0) || !is.finite(max(abs(P0)))) {
      ## Fallback: start from zero (less accurate but still unbiased)
      particles <- matrix(0, nrow = n_st, ncol = N)
      ## Still need to consume / record the init normals slot
      U_init_used <- if (!is.null(U_init_normals)) U_init_normals else
                     matrix(0, nrow = n_s, ncol = N)
    } else {
      ## Draw x1_0^i ~ N(0, P0); x2_0^i = 0 (and x3_0^i = 0 at order 3)
      L_P0 <- tryCatch(t(chol(P0 + diag(1e-12, n_s))),
                       error = function(e) diag(sqrt(diag(P0) + 1e-12), n_s))
      z_init <- if (!is.null(U_init_normals)) U_init_normals else
                matrix(rnorm(n_s * N), nrow = n_s)
      U_init_used <- z_init
      x1_init     <- L_P0 %*% z_init
      particles   <- rbind(x1_init, matrix(0, nrow = n_st - n_s, ncol = N))
    }

    ## ---- Optional particle burn-in (burn_in_init > 0) --------------------
    ## The Lyapunov draw above puts x1_0 at the FIRST-order stationary
    ## distribution but x2_0 (and x3_0) at exactly 0 — while the pruned
    ## model's true stationary joint has a NONZERO x2 mean (the
    ## risk-adjustment level shift) and nonzero higher-layer dispersion.
    ## On short samples the resulting initial-level mismatch is absorbed as
    ## extra persistence: the 2026-08-04 order-3 SBC certification, after
    ## the obs-tensor fix removed the dominant bias, still showed a
    ## rho_a-specific rank shift (mean_rank_z -3.0) with sig_b clean.
    ## Propagating the cloud burn_in_init periods with fresh shocks before
    ## the first observation starts it from (approximately) the full pruned
    ## stationary joint instead — the same stationarity assumption the
    ## data-side simulate_model_* burn-in and the augmented-state pruned KF
    ## already make. Opt-in (default 0) so fixed-seed pins are unchanged;
    ## incompatible with an externally supplied U_list (the CPM slot layout
    ## has no burn-in slots), which is rejected at the factory boundary.
    if (burn_in_init > 0L) {
      st_b  <- dr2$state_idx
      hx_b  <- dr2$ghx [st_b, , drop = FALSE]
      hu_b  <- dr2$ghu [st_b, , drop = FALSE]
      hxx_b <- dr2$ghxx[st_b, , drop = FALSE]
      hxu_b <- dr2$ghxu[st_b, , drop = FALSE]
      huu_b <- dr2$ghuu[st_b, , drop = FALSE]
      hss_b <- dr2$ghss[st_b]
      for (b in seq_len(burn_in_init)) {
        shocks_b <- L_e %*% matrix(rnorm(n_e * N), nrow = n_e)
        particles <- if (order3) {
          .tpf_propagate3_R(
            particles, shocks_b, hx_b, hu_b, hxx_b, hxu_b, huu_b, hss_b,
            hxxx = dr2$ghxxx[st_b, , drop = FALSE],
            hxxu = dr2$ghxxu[st_b, , drop = FALSE],
            hxuu = dr2$ghxuu[st_b, , drop = FALSE],
            huuu = dr2$ghuuu[st_b, , drop = FALSE],
            hxss = if (!is.null(dr2$ghxss)) dr2$ghxss[st_b, , drop = FALSE],
            huss = if (!is.null(dr2$ghuss)) dr2$ghuss[st_b, , drop = FALSE],
            hs3  = if (!is.null(dr2$ghs3))  dr2$ghs3[st_b])
        } else {
          .tpf_propagate_R(particles, shocks_b, hx_b, hu_b,
                           hxx_b, hxu_b, huu_b, hss_b)
        }
      }
    }

    ## ---- State-row slices for FULLY-missing period propagation -----------
    ## Mirrors the identical extraction inside tpf_run_period[3] (see there
    ## for the recursion reference).  Needed directly here because a
    ## fully-missing period propagates the particle cloud one step WITHOUT
    ## going through tpf_run_period[3]'s tempering loop (no observation to
    ## temper against).
    st_idx_s <- dr2$state_idx
    hx_s  <- dr2$ghx [st_idx_s, , drop = FALSE]
    hu_s  <- dr2$ghu [st_idx_s, , drop = FALSE]
    hxx_s <- dr2$ghxx[st_idx_s, , drop = FALSE]
    hxu_s <- dr2$ghxu[st_idx_s, , drop = FALSE]
    huu_s <- dr2$ghuu[st_idx_s, , drop = FALSE]
    hss_s <- dr2$ghss[st_idx_s]
    hxxx_s <- NULL; hxxu_s <- NULL; hxuu_s <- NULL; huuu_s <- NULL
    hxss_s <- NULL; huss_s <- NULL; hs3_s  <- NULL
    if (order3) {
      hxxx_s <- dr2$ghxxx[st_idx_s, , drop = FALSE]
      hxxu_s <- dr2$ghxxu[st_idx_s, , drop = FALSE]
      hxuu_s <- dr2$ghxuu[st_idx_s, , drop = FALSE]
      huuu_s <- dr2$ghuuu[st_idx_s, , drop = FALSE]
      hxss_s <- if (!is.null(dr2$ghxss)) dr2$ghxss[st_idx_s, , drop = FALSE] else NULL
      huss_s <- if (!is.null(dr2$ghuss)) dr2$ghuss[st_idx_s, , drop = FALSE] else NULL
      hs3_s  <- if (!is.null(dr2$ghs3))  dr2$ghs3[st_idx_s]                  else NULL
    }

    ## ---- Main TPF loop over periods t = 1, ..., T -----------------------
    ## U_list layout: slot 1 = init normals; slots 2..(T+1) = shock normals;
    ## slots (T+2)..(2T+1) = phi=1 resampling z scalars (NULL -> draw from RNG);
    ## slots (2T+2)..(3T+1) = Option A mutation noise matrices (NULL when n_mh=0).
    loglik     <- 0
    ## U_realized: (3T+1)-element list; same layout as U_list.
    U_realized <- vector("list", 3L * T_obs + 1L)
    U_realized[[1L]] <- U_init_used

    for (t in seq_len(T_obs)) {
      y_t <- data[, t]
      obs_mask_t <- is.finite(y_t)
      n_obs_t    <- sum(obs_mask_t)

      ## Slot t+1 = period-t shock normals (n_e x N)
      U_normals_t <- if (!is.null(U_list)) U_list[[t + 1L]] else NULL

      if (n_obs_t == 0L) {
        ## -- FULLY-missing period: pure predictive step -------------------
        ## Draw (or, under CPM, consume the U_list slot for) this period's
        ## shocks and propagate every particle one step with the existing
        ## propagation kernels. No observation -> no weight/resample/mutate
        ## and the loglik contribution is exactly 0. This is the correct
        ## predictive behavior and keeps the CPM slot layout aligned by
        ## period (U_realized[[t+1]] still gets the shock normals used).
        ## The pre-fix behavior (`next`, skipping the period entirely) glued
        ## the periods on either side of the gap together with a single
        ## multi-period transition -- the likelihood of a different model.
        z_mat_na  <- if (!is.null(U_normals_t)) U_normals_t else
                     matrix(rnorm(n_e * N), nrow = n_e)
        shocks_na <- L_e %*% z_mat_na
        particles <- if (order3) {
          .tpf_propagate3_R(particles, shocks_na,
                             hx_s, hu_s, hxx_s, hxu_s, huu_s, hss_s,
                             hxxx_s, hxxu_s, hxuu_s, huuu_s,
                             hxss_s, huss_s, hs3_s)
        } else if (.HAS_RCPP_TPF()) {
          tpf_propagate_particles(
            particles, shocks_na, hx_s, hu_s,
            if (is.null(hxx_s)) matrix(0, n_s, n_s * n_s) else hxx_s,
            if (is.null(hxu_s)) matrix(0, n_s, n_s * n_e) else hxu_s,
            if (is.null(huu_s)) matrix(0, n_s, n_e * n_e) else huu_s,
            if (is.null(hss_s)) numeric(n_s) else as.numeric(hss_s))
        } else {
          .tpf_propagate_R(particles, shocks_na,
                            hx_s, hu_s, hxx_s, hxu_s, huu_s, hss_s)
        }
        U_realized[[t + 1L]]         <- z_mat_na
        U_realized[[T_obs + 1L + t]] <- NA_real_
        next
      }

      ## Slot T_obs+1+t = period-t phi=1 resampling z scalar
      U_resample_t <- if (!is.null(U_list) && length(U_list) > T_obs + 1L)
                        U_list[[T_obs + 1L + t]] else NULL
      ## Slot 2*T_obs+1+t = Option A mutation noise matrix.
      ## When U_list is NULL or absent: pre-draw fresh mutation buffer if n_mh > 0.
      ## When U_list has this slot: use it (CPM path, matrix with correlated draws).
      ## When n_mh == 0: NULL (mutation never fires, no draws needed).
      U_mutation_t <- if (!is.null(U_list) && length(U_list) > 2L * T_obs + 1L) {
        slot <- U_list[[2L * T_obs + 1L + t]]
        ## Legacy: if slot is a K-vector (not a matrix), treat as NULL (pre-Option-A)
        if (!is.null(slot) && is.matrix(slot)) slot else NULL
      } else if (has_mutation) {
        ## No U_list supplied: pre-draw fresh mutation buffer for this period.
        ## These become the "used" draws that get stored in U_realized and
        ## can be AR(1) updated by rwmh_cpm for the next CPM step.
        matrix(rnorm(mut_buf_rows * mut_buf_cols), nrow = mut_buf_rows)
      } else {
        NULL  ## n_mh=0: no mutation buffer
      }

      ## -- PARTIALLY-missing period: subset to observed rows only. --------
      ## `has_na_t == FALSE` (fully-observed) takes the ORIGINAL, un-copied
      ## matrices/vectors -- exactly the pre-fix call arguments -- so that
      ## branch is bit-identical to previous releases (requirement 3).
      has_na_t <- n_obs_t < length(y_t)
      if (has_na_t) {
        idx_obs      <- which(obs_mask_t)
        y_t_use      <- y_t[idx_obs]
        ZZ_use       <- ZZ[idx_obs, , drop = FALSE]
        DD_use       <- DD[idx_obs, , drop = FALSE]
        d_obs_use    <- d_obs[idx_obs]
        ghss_obs_use <- ghss_obs[idx_obs]
        hxx_obs_use  <- if (!is.null(hxx_obs)) hxx_obs[idx_obs, , drop = FALSE] else NULL
        hxu_obs_use  <- if (!is.null(hxu_obs)) hxu_obs[idx_obs, , drop = FALSE] else NULL
        huu_obs_use  <- if (!is.null(huu_obs)) huu_obs[idx_obs, , drop = FALSE] else NULL
        hxxu_obs_use <- NULL; hxuu_obs_use <- NULL; hxxx_obs_use <- NULL
        huuu_obs_use <- NULL; huss_obs_use <- NULL
        ZZ_xss_use   <- NULL; ghs3_obs_use <- NULL
        if (order3) {
          hxxu_obs_use <- if (!is.null(hxxu_obs)) hxxu_obs[idx_obs, , drop = FALSE] else NULL
          hxuu_obs_use <- if (!is.null(hxuu_obs)) hxuu_obs[idx_obs, , drop = FALSE] else NULL
          hxxx_obs_use <- if (!is.null(hxxx_obs)) hxxx_obs[idx_obs, , drop = FALSE] else NULL
          huuu_obs_use <- if (!is.null(huuu_obs)) huuu_obs[idx_obs, , drop = FALSE] else NULL
          huss_obs_use <- if (!is.null(huss_obs)) huss_obs[idx_obs, , drop = FALSE] else NULL
          ZZ_xss_use   <- if (!is.null(ZZ_xss))   ZZ_xss[idx_obs, , drop = FALSE]   else NULL
          ghs3_obs_use <- if (!is.null(ghs3_obs)) ghs3_obs[idx_obs]                 else NULL
        }
      } else {
        y_t_use      <- y_t
        ZZ_use       <- ZZ
        DD_use       <- DD
        d_obs_use    <- d_obs
        ghss_obs_use <- ghss_obs
        hxx_obs_use  <- hxx_obs
        hxu_obs_use  <- hxu_obs
        huu_obs_use  <- huu_obs
        hxxu_obs_use <- hxxu_obs
        hxuu_obs_use <- hxuu_obs
        hxxx_obs_use <- hxxx_obs
        huuu_obs_use <- huuu_obs
        huss_obs_use <- huss_obs
        ZZ_xss_use   <- ZZ_xss
        ghs3_obs_use <- ghs3_obs
      }

      res <- if (order3) {
        ## Pure-R order-3 period (no C++ kernel in this increment)
        tpf_run_period3(
          particles   = particles,
          y_t         = y_t_use,
          dr3         = dr2,           # DecisionRules3 (superset of dr2 fields)
          Sigma_e     = Sigma_e,
          L_e         = L_e,
          ZZ          = ZZ_use,
          DD          = DD_use,
          d_obs       = d_obs_use,
          ghss_obs    = ghss_obs_use,
          ZZ_xss      = ZZ_xss_use,
          ghs3_obs    = ghs3_obs_use,
          hxx_obs     = hxx_obs_use,
          hxu_obs     = hxu_obs_use,
          huu_obs     = huu_obs_use,
          hxxu_obs    = hxxu_obs_use,
          hxuu_obs    = hxuu_obs_use,
          hxxx_obs    = hxxx_obs_use,
          huuu_obs    = huuu_obs_use,
          huss_obs    = huss_obs_use,
          me_variance = me_variance,
          ess_target  = ess_target,
          n_mh        = n_mh,
          mh_scale    = mh_scale,
          U_normals   = U_normals_t,
          U_resample  = U_resample_t,
          U_mid       = NULL,
          U_mutation  = U_mutation_t
        )
      } else tpf_run_period(
        particles   = particles,
        y_t         = y_t_use,
        dr2         = dr2,
        Sigma_e     = Sigma_e,
        L_e         = L_e,
        ZZ          = ZZ_use,
        DD          = DD_use,
        d_obs       = d_obs_use,
        ghss_obs    = ghss_obs_use,
        hxx_obs     = hxx_obs_use,
        hxu_obs     = hxu_obs_use,
        huu_obs     = huu_obs_use,
        me_variance = me_variance,
        ess_target  = ess_target,
        n_mh        = n_mh,
        mh_scale    = mh_scale,
        use_rcpp    = .HAS_RCPP_TPF(),
        backend     = if (.HAS_RCPP_TPF_PERIOD()) "cpp" else "R",
        U_normals   = U_normals_t,
        U_resample  = U_resample_t,
        U_mid       = NULL,            # legacy: always NULL now (never fires)
        U_mutation  = U_mutation_t     # Option A: pre-drawn mutation buffer
      )
      particles                          <- res$particles
      loglik                             <- loglik + res$log_lik_contrib
      U_realized[[t + 1L]]               <- res$U_used          # n_e x N shock normals
      ## Phi=1 resampling z: assign via list() wrapping to avoid NULL-removal.
      ## In R, `x[[i]] <- NULL` REMOVES element i; wrapping in list() stores the value.
      ## For scalar z (non-CPM path returns NA_real_, not NULL), direct assignment is safe.
      U_realized[[T_obs + 1L + t]]       <- res$z_resample_used  # scalar z (NA if non-CPM)
      ## Mutation buffer slot: only assign when non-NULL to avoid list shrinkage.
      ## (In R, `x[[i]] <- NULL` removes element i, shrinking the pre-allocated list.)
      ## When n_mh=0, U_mutation_t is NULL and the slot stays as the pre-allocated NULL.
      if (!is.null(U_mutation_t))
        U_realized[[2L * T_obs + 1L + t]] <- U_mutation_t        # mutation matrix
    }

    if (!is.finite(loglik))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## ---- System priors --------------------------------------------------
    if (!is.null(system_priors)) {
      sp_lp <- .eval_system_priors(
        system_priors,
        list(theta   = theta,
             model   = model,
             dr      = dr1,
             Sigma_e = Sigma_e,
             params  = params))
      if (!is.finite(sp_lp))
        return(list(logpost = -Inf, loglik = loglik, logprior = lp))
      lp <- lp + sp_lp
    }

    list(logpost = .dynhr_opt("power_posterior", default = 1) * loglik + lp,
         loglik = loglik, logprior = lp,
         U_list = U_realized)  # CPM: standard normals used per period (n_e x N each)
  }
}
