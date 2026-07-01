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
    kron_xx <- matrix(0, nrow = n_s * n_s, ncol = N)
    for (i in seq_len(N)) {
      kron_xx[, i] <- kronecker(x1[, i], x1[, i])
    }
    x2_new <- x2_new + 0.5 * hxx %*% kron_xx
  }
  if (any(hxu != 0)) {
    ## (e %x% x1_prev): (exo SLOW x state FAST) — matches hxu col ordering
    ## hxu cols: (state FAST, exo SLOW) -> vector is (e %x% x1)
    kron_ex <- matrix(0, nrow = n_e * n_s, ncol = N)
    for (i in seq_len(N)) {
      kron_ex[, i] <- kronecker(shocks[, i], x1[, i])
    }
    x2_new <- x2_new + hxu %*% kron_ex
  }
  if (any(huu != 0)) {
    kron_ee <- matrix(0, nrow = n_e * n_e, ncol = N)
    for (i in seq_len(N)) {
      kron_ee[, i] <- kronecker(shocks[, i], shocks[, i])
    }
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
#' @return Numeric vector of length N.
#' @noRd
.tpf_log_weights_R <- function(particles, y_t, ZZ, DD, shocks_t,
                                d_obs, ghss_obs, me_variance, phi,
                                DD_shock_mean = NULL) {
  n_s <- ncol(ZZ)
  N   <- ncol(particles)

  x1 <- particles[seq_len(n_s), , drop = FALSE]
  x2 <- particles[seq_len(n_s) + n_s, , drop = FALSE]

  ## Observation mean: ZZ*(x1+x2) + DD*e_t + d_obs + ghss_obs  (n_obs x N)
  mu <- ZZ %*% (x1 + x2)

  if (!is.null(shocks_t)) {
    mu <- mu + DD %*% shocks_t        # n_obs x N
  } else if (!is.null(DD_shock_mean)) {
    mu <- mu + DD_shock_mean          # broadcast single column
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
#' @param me_variance Scalar > 0.
#' @param ess_target  ESS target fraction (default 0.5).
#' @param n_mh        RWMH mutation steps per phi stage (default 1).
#' @param mh_scale    RWMH proposal scale relative to cloud covariance.
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
                                  d_obs, ghss_obs, me_variance, phi = 1.0)

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
    if (resampled && n_mh > 0L && N > 1L) {
      ## Cloud covariance of the state component (weighted by w_norm)
      parts_T   <- t(particles)   # N x (2*n_s)
      w_bar     <- as.numeric(w_norm %*% parts_T)
      parts_c   <- sweep(parts_T, 2, w_bar)
      Sigma_hat <- crossprod(parts_c * w_norm, parts_c)
      Sigma_hat <- Sigma_hat + diag(1e-8, n_2s)

      L_prop <- tryCatch(
        t(chol(mh_scale^2 * Sigma_hat)),
        error = function(e) diag(mh_scale * 0.01, n_2s)
      )

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
      has_u_mut   <- !is.null(U_mutation) && is.matrix(U_mutation)
      u_mut_ncols <- if (has_u_mut) ncol(U_mutation) else 0L
      stage0      <- stage - 1L   # C++ 0-based stage index

      for (i in seq_len(N)) {
        s_i   <- particles[, i]
        e_i   <- shocks[, i]
        e_mat <- matrix(e_i, nrow = n_e, ncol = 1L)
        s_mat <- matrix(s_i, nrow = n_2s, ncol = 1L)
        tlp_i <- phi_curr *
          .tpf_log_weights_R(s_mat, y_t, ZZ, DD, e_mat,
                              d_obs, ghss_obs, me_variance, phi = 1.0)[1L]

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
            z_s   <- rnorm(n_2s)
            z_e   <- rnorm(n_e)
            log_u <- log(runif(1L))
          }

          s_prop  <- s_i + as.numeric(L_prop %*% z_s)
          e_prop  <- as.numeric(L_e %*% z_e)   # resample shock
          e_mat_p <- matrix(e_prop, nrow = n_e, ncol = 1L)
          s_mat_p <- matrix(s_prop, nrow = n_2s, ncol = 1L)
          tlp_p   <- phi_curr *
            .tpf_log_weights_R(s_mat_p, y_t, ZZ, DD, e_mat_p,
                                d_obs, ghss_obs, me_variance, phi = 1.0)[1L]
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
                                      d_obs, ghss_obs, me_variance, phi = 1.0)
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
        "  WARNING: SD > 1 (Dynare threshold). PMCMC acceptance will be dominated by\n",
        "  loglik noise. Current n_particles = %s; estimated n_particles needed for\n",
        "  SD < 1: ~%s (N_needed ~ N * SD^2).",
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
#' order-2 state space.
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
#' @param n_mh        RWMH mutation steps per phi stage (default 1).
#' @param mh_scale    RWMH proposal scale relative to particle cloud covariance.
#' @param seed        Integer RNG seed for reproducibility (NULL = no fixed seed).
#' @param system_priors Optional system priors list (see \code{\link{sp_irf}}).
#' @param max_stages_u  Maximum number of tempering stages per observation
#'   (default 16).  The filter raises phi from 0 to 1 in at most this many
#'   steps; a smaller value speeds up the filter at the cost of coarser
#'   tempering.
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
                                    max_stages_u  = 16L) {

  ## Force promises for closure-capture safety (same pattern as cumulant branch)
  force(data); force(prior_spec); force(obs_vars); force(me_variance)
  force(n_particles); force(ess_target); force(n_mh); force(mh_scale)
  force(seed); force(system_priors); force(max_stages_u)
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

    ## ---- Second-order perturbation (per-draw, expensive) ----------------
    Sigma_e <- .get_shock_cov(model, model$varexo_names, params)
    dr2 <- tryCatch(
      solve_perturbation_order2(model, compiled, ss_result$ss, params,
                                 dr1, Sigma_e = Sigma_e, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(dr2))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

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
    mut_buf_rows <- n_2s + n_e + 1L
    mut_buf_cols <- as.integer(max_stages_u) * as.integer(n_mh) * N
    has_mutation <- (n_mh > 0L)

    if (is.null(P0) || !is.finite(max(abs(P0)))) {
      ## Fallback: start from zero (less accurate but still unbiased)
      particles <- matrix(0, nrow = n_2s, ncol = N)
      ## Still need to consume / record the init normals slot
      U_init_used <- if (!is.null(U_init_normals)) U_init_normals else
                     matrix(0, nrow = n_s, ncol = N)
    } else {
      ## Draw x1_0^i ~ N(0, P0); x2_0^i = 0
      L_P0 <- tryCatch(t(chol(P0 + diag(1e-12, n_s))),
                       error = function(e) diag(sqrt(diag(P0) + 1e-12), n_s))
      z_init <- if (!is.null(U_init_normals)) U_init_normals else
                matrix(rnorm(n_s * N), nrow = n_s)
      U_init_used <- z_init
      x1_init     <- L_P0 %*% z_init
      particles   <- rbind(x1_init, matrix(0, nrow = n_s, ncol = N))
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
      if (any(!is.finite(y_t))) next

      ## Slot t+1 = period-t shock normals (n_e x N)
      U_normals_t <- if (!is.null(U_list)) U_list[[t + 1L]] else NULL
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

      res <- tpf_run_period(
        particles   = particles,
        y_t         = y_t,
        dr2         = dr2,
        Sigma_e     = Sigma_e,
        L_e         = L_e,
        ZZ          = ZZ,
        DD          = DD,
        d_obs       = d_obs,
        ghss_obs    = ghss_obs,
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
