## R/hank-kalman-ar-grad.R
## --------------------------------------------------------------------------
## Analytic / semi-analytic SCORE of the exact-AR sequence-space likelihood
## (R/hank-kalman-ar.R's hank_loglik_ar()), and a make_log_posterior_hank()-
## style gradient closure over its (rho, sigma) parameters.
##
## WHY.  dynhr's derivative stack (21 gradient/Hessian/adjoint files) serves
## ONLY the state-space/perturbation-decision-rule likelihood path
## (make_posterior_grad(), grad_method = "hybrid"/"implicit"/"adjoint"/
## "adjoint_solution"). hank_loglik_ar() -- the sequence-space exact-AR
## likelihood -- has none: every reference to hank_loglik_ar/autocov_slab/
## stacked/H_U across the whole stack is zero (verified by grep, not
## assumed). Raised from the NZ HANK paper
## (notes/dynhr_brief_seqspace_gradients.md in the paper repo); this file is
## the package-side port of that brief's validated prototype
## (paper repo: replication/R/hank_abrs_grad.R, hank_abrs_score() -- checked
## against numerical differentiation there to 3.3e-9 relative error on a
## synthetic fixture sized to the paper's own model).
##
## WHY make_posterior_grad() ISN'T extended instead (see briefs/20).
## make_posterior_grad() is built entirely around a PERTURBATION DECISION
## RULE (cache_system_structure(), .solve_dr(), an is_sigma classification by
## diffing ghx/ghu across a param move) -- none of that applies to a
## Theta_list. The natural fit mirrors how make_log_posterior_hank() ALREADY
## sits beside make_log_posterior() as a standalone parallel function rather
## than a branch inside it: make_posterior_grad_hank_ar() below does the same
## for the gradient.
##
## STRUCTURE (mirrors hank_loglik_ar() exactly; see that function's header):
##   G = sum_z sigma_z^2 A_z(Theta_z, rho_z) ; G[,,1] += diag(me_sd^2)
##   S = gather(G) via the fixed index built once per (Td, n_obs, NA-pattern)
##   loglik = -0.5 [ N log 2pi + log|S| + y' S^-1 y ]
## The gather is LINEAR, so dS/dphi = gather(dG/dphi), and with v = S^-1 y,
## Sinv = S^-1 (one Cholesky/inverse, shared by every parameter):
##   d loglik/d phi = -0.5 [ tr(Sinv dS) - v' dS v ]      (dS symmetric)
##
##   sigma_z : dG/dsigma_z = 2 sigma_z A_z -- FULLY ANALYTIC. A_z is exactly
##             the slab hank_loglik_ar() already builds and caches
##             (.hank_ar_autocov_slab()/cache$A), so a sigma-only move
##             recomputes NO new slabs.
##   me_sd_i : dG[,,1]/d(me_sd_i) is the rank-1 diagonal e_i e_i' * 2 me_sd_i
##             -- trivial closed form, no gather needed (see .hank_ar_dscore_me).
##   rho_z   : dA_z/drho_z has three channels (Theta_z's own geometric driving
##             path, the quasi-difference Psi = Th - rho*lag(Th), and the
##             weight kernels Wp/Wm). A fully analytic form is derivable but
##             error-prone (see the brief), so SEMI-ANALYTIC: central-
##             difference the SLAB ONLY via the existing (uncached, so FD taps
##             never evict a real cache entry) .hank_ar_autocov_kernel(), and
##             keep the expensive S-level algebra (Cholesky, inverse, traces)
##             analytic -- this avoids one Cholesky PER PARAMETER, which is
##             where a naive numerical gradient's cost actually lives.
##
## COST.  The contractions are linear in dS = gather(dG), so both collapse
## onto the compact G-space under the fixed gather index: one compiled pass
## (src/hank_ar_score.cpp) accumulates the weights w[k] = sum_{Sidx=k}
## (Sinv[r,c] - v_r v_c), after which each parameter is a n_obs^2*T_data dot
## product instead of an n_kept^2 gather + trace + matvec. Measured on the
## paper-sized fixture (T_h=400, 12 observables, 132 quarters, 11 shocks;
## installed -O2 build): sigma block 0.056 s vs 0.255 s for the per-parameter
## gather path and 0.424 s for numerical central differences (7.6x); whole
## shock block incl. the slab-bound rho FDs, 0.241 s vs 1.119 s (4.6x). The
## R-level rowsum() fallback lands in between (0.129 s) -- it is kept only as
## the parity reference (`use_cpp = FALSE`), matching the paper brief's
## finding that this win is contingent on the accumulation being compiled.
##
## THE KERNEL ADJOINT (rho_method = "adjoint", wrt = "theta").
## .hank_ar_slab_adjoint() differentiates the autocovariance kernel in reverse,
## giving dl/dTheta_z and the rho score with NO differencing anywhere. Measured
## on the same fixture (installed -O2):
##   whole shock block, rho via FD-the-slab : 0.194 s
##   whole shock block, rho via the adjoint : 0.215 s   (0.90x -- SLOWER)
##   dl/dTheta_z for all 11 shocks          : 0.112 s
## So the adjoint is NOT a speed win on the rho block: it buys exactness (no
## h_rho) and it needs dTheta_z/drho_z rather than two Theta re-derivations.
## The real payoff is dl/dTheta_z, which nothing else here can produce and
## which every STRUCTURAL parameter needs -- with it, a structural gradient
## costs one adjoint pass plus one linear-map application per parameter,
## instead of two cold likelihood evaluations (0.28 s) plus two model rebuilds
## PER PARAMETER. See briefs/21-structural-score-api-scope.md.
##
## Regenerating Theta_z at a perturbed rho needs a full sequence-space IRF
## (Theta_z's OWN definition is the response to the driving path rho_z^t), so
## unlike sigma/me the rho block cannot be computed from Y/Theta_list alone:
## hank_loglik_ar_grad() takes a `theta_fn(z, rho_z)` callback for this;
## make_posterior_grad_hank_ar() supplies one built from hank_model_irf()
## (cheap: model$G is the sequence-space Jacobian, precomputed once at
## hank_model() construction, so this is a small matvec, not a re-solve).
## --------------------------------------------------------------------------


## Symmetric bilinear contraction d(loglik)/d(phi) = -0.5*(tr(Sinv dS) -
## v' dS v), for dS already gathered to the stacked (kept) shape.
.hank_ar_dterm <- function(Sinv, v, dS) {
  dS <- (dS + t(dS)) / 2
  -0.5 * (sum(Sinv * dS) - as.numeric(crossprod(v, dS %*% v)))
}


## COMPACT ACCUMULATION (the compiled fast path; see src/hank_ar_score.cpp).
##
## Both contractions above are linear in dS = gather(dG) under the FIXED
## gather index, so they collapse onto the compact G-space (n_obs^2 * T_data
## entries, ~19k at paper scale vs ~2.5M stacked entries):
##
##   tr(Sinv gather(X)) - v' gather(X) v = sum_k w[k] X[k],
##   w[k] = sum_{(r,c): Sidx[r,c] = k} (Sinv[r,c] - v[r] v[c])
##
## Accumulating `w` costs ONE pass over the stacked entries, after which every
## parameter is a ~19k dot product (.hank_ar_dterm_compact) instead of a fresh
## N x N gather + N^2 trace + N^2 matvec. The R fallback below (rowsum) is the
## parity reference and is what the paper measured as slower than the naive
## gather -- the win is contingent on the compiled pass (briefs/20 5).
.hank_ar_score_weights <- function(Sinv, v, Sidx, n_g, use_cpp = TRUE) {
  if (isTRUE(use_cpp))
    return(hank_ar_score_weights_cpp(Sidx, Sinv, as.numeric(v),
                                     as.integer(n_g)))
  ww <- rowsum(as.numeric(Sinv - tcrossprod(as.numeric(v))),
               as.integer(Sidx), reorder = FALSE)
  w <- numeric(n_g)
  w[as.integer(rownames(ww))] <- ww[, 1L]
  w
}


## d(loglik)/d(phi) from the compact weights and dG/dphi in G-space (an
## n_obs x n_obs x T_data array, e.g. 2*sigma_z*A_z).
.hank_ar_dterm_compact <- function(w, dG) -0.5 * sum(w * as.numeric(dG))


## ADJOINT OF THE AUTOCOVARIANCE KERNEL (.hank_ar_autocov_kernel).
##
## Everything a structural parameter does to this likelihood, it does through
## the MA coefficients: theta -> Theta_z -> A_z -> S -> l. This routine
## supplies the middle link in reverse. Given the weight array `Wz` such that
## the scalar being differentiated is L = sum_k <Wz[,,k], A_z[,,k]> (for the
## loglik, Wz = -0.5 * sigma_z^2 * compact weights, reshaped), it returns
## dL/dTheta_z and the PARTIAL dL/drho_z at FIXED Theta.
##
## Derivation (notation as in .hank_ar_autocov_kernel):
##   A(k) = sum_d [ Wp[d,k] m_d + Wm[d,k] m_d' ] - Wp[0,k] m_0,
##   m_d  = sum_j Psi[j+d,]' Psi[j,],   Psi = Th - rho lag(Th)
## so L = sum_d <B_d, m_d> with
##   B_d = sum_k ( Wp[d,k] W_k + Wm[d,k] W_k' ) - [d = 0] sum_k Wp[0,k] W_k,
## which is two small matrix products against the SAME Wp/Wm the forward
## kernel builds (the reverse of its Gp/Gm folds). Differentiating m_d,
##   dL/dPsi[t,] = sum_{d<=t-1} B_d Psi[t-d,] + sum_{d<=q-t} B_d' Psi[t+d,]
## and the quasi-difference gives dL/dTh[t,] = dL/dPsi[t,] - rho dL/dPsi[t+1,].
## The rho partial has two channels: through Psi (dPsi/drho = -lag(Th)) and
## through the weight kernels themselves (dWp/drho, dWm/drho below).
##
## WHY THIS EXISTS.  The shipped rho block finite-differences the whole slab
## (two slab builds per shock) and therefore needs a `theta_fn` callback to
## re-derive Theta at a perturbed persistence; it also cannot serve structural
## parameters at all, since those need dL/dTheta itself. Both fall out of one
## adjoint pass here. The TOTAL rho derivative still needs dTheta_z/drho_z
## (Theta_z is DEFINED as the response to the rho_z^t driving path) -- but that
## is one linear-map application, not a re-solve, because Theta is linear in
## the driving path: see `dtheta_fn` in hank_loglik_ar_grad().
.hank_ar_slab_adjoint <- function(Wz, Th, rho, q, n_lags,
                                  want_theta = TRUE, want_rho = TRUE,
                                  use_cpp = TRUE) {
  n_obs <- ncol(Th)
  nn <- n_obs * n_obs
  Th <- Th[seq_len(q), , drop = FALSE]
  Psi <- Th - rho * rbind(rep(0, n_obs), Th[-q, , drop = FALSE])
  ds <- 0:(q - 1L); ks <- 0:n_lags
  omr <- 1 - rho^2
  Wp <- outer(ds, ks, function(d, k) rho^abs(k - d)) / omr
  Wm <- outer(ds, ks, function(d, k) rho^(k + d)) / omr

  ## col k = vec(W_k) and vec(W_k'): the two folds the forward kernel applies.
  Wmat  <- matrix(as.numeric(Wz), nn, n_lags + 1L)
  WmatT <- matrix(as.numeric(aperm(Wz, c(2L, 1L, 3L))), nn, n_lags + 1L)
  Bvec <- Wmat %*% t(Wp) + WmatT %*% t(Wm)                  # nn x q, col d+1
  Bvec[, 1L] <- Bvec[, 1L] - as.numeric(Wmat %*% Wp[1L, ])  # the -m_0 fold

  ## dL/dPsi: a block-Toeplitz contraction, both folds of every B_d. This is
  ## the adjoint's dominant cost (~q^2 n_obs^2), hence the compiled kernel; the
  ## R loop below is its parity reference (use_cpp = FALSE).
  GPsi <- if (isTRUE(use_cpp)) {
    hank_ar_slab_adjoint_psi_cpp(Bvec, Psi)
  } else {
    out <- matrix(0, q, n_obs)
    for (d in ds) {
      Bd <- matrix(Bvec[, d + 1L], n_obs, n_obs)
      if (all(Bd == 0)) next                   # rho^d underflow: nothing left
      idx <- seq_len(q - d)
      out[idx + d, ] <- out[idx + d, ] + Psi[idx, , drop = FALSE] %*% t(Bd)
      out[idx, ] <- out[idx, ] + Psi[idx + d, , drop = FALSE] %*% Bd
    }
    out
  }

  g_theta <- NULL
  if (isTRUE(want_theta)) {
    g_theta <- GPsi
    g_theta[-q, ] <- g_theta[-q, ] - rho * GPsi[-1L, , drop = FALSE]
  }

  g_rho <- NULL
  if (isTRUE(want_rho)) {
    ## channel 1: Psi = Th - rho lag(Th)
    dPsi <- -rbind(rep(0, n_obs), Th[-q, , drop = FALSE])
    g_rho <- sum(GPsi * dPsi)
    ## channel 2: the weight kernels. d/drho [rho^e / (1-rho^2)] =
    ## e rho^(e-1)/(1-rho^2) + (rho^e/(1-rho^2)) * 2 rho/(1-rho^2); the first
    ## term is identically 0 at e = 0 (never evaluate rho^-1).
    pow_d <- function(e) ifelse(e == 0, 0, e * rho^pmax(e - 1L, 0L))
    dWp <- outer(ds, ks, function(d, k) pow_d(abs(k - d))) / omr +
      Wp * (2 * rho / omr)
    dWm <- outer(ds, ks, function(d, k) pow_d(k + d)) / omr +
      Wm * (2 * rho / omr)
    dBvec <- Wmat %*% t(dWp) + WmatT %*% t(dWm)
    dBvec[, 1L] <- dBvec[, 1L] - as.numeric(Wmat %*% dWp[1L, ])
    M <- matrix(unlist(lapply(ds, function(d) {
      idx <- seq_len(q - d)
      crossprod(Psi[idx + d, , drop = FALSE], Psi[idx, , drop = FALSE])
    })), nrow = q, byrow = TRUE)                # row d+1 = vec(m_d)
    g_rho <- g_rho + sum(dBvec * t(M))
  }

  list(theta = g_theta, rho = g_rho)
}


## Diagonal me_sd score: dG[,,1]/d(me_sd_i) = 2*me_sd_i*e_i e_i', which is
## already diagonal in G, hence diagonal in the gather (a gathered entry
## S[r,c] is nonzero only at r == c with observable(r) == i). No new gather.
.hank_ar_dscore_me <- function(Sinv, v, me_sd, obs_of_row) {
  n_obs <- length(me_sd)
  vapply(seq_len(n_obs), function(i) {
    idx <- which(obs_of_row == i)
    if (!length(idx)) return(0)
    -me_sd[i] * (sum(diag(Sinv)[idx]) - sum(v[idx]^2))
  }, numeric(1))
}


#' Analytic / semi-analytic score of the exact-AR sequence-space likelihood
#'
#' Gradient counterpart of \code{\link{hank_loglik_ar}}: the sigma block is
#' exact-analytic (reuses the same cached autocovariance slabs the
#' likelihood itself builds -- a sigma-only move recomputes zero new slabs),
#' the me_sd block is a trivial closed form, and the rho block is
#' semi-analytic (central-differences the per-shock slab only, keeping the
#' expensive stacked-covariance algebra -- one Cholesky, shared across every
#' parameter -- fully analytic). See the file header for the derivation and
#' \code{briefs/20-seqspace-gradient-scope.md} for the validation history.
#'
#' @inheritParams hank_loglik_ar
#' @param wrt Character subset of \code{c("sigma", "rho", "me", "theta")}:
#'   which score blocks to compute. Requesting only \code{"sigma"} (the default
#'   use inside \code{\link{make_posterior_grad_hank_ar}} also asks for
#'   \code{"rho"}) skips the (relatively) more expensive rho block entirely.
#'   \code{"theta"} returns \code{d loglik / d Theta_z} (a \code{q x n_obs}
#'   matrix per shock) from the adjoint of the autocovariance kernel -- the
#'   quantity any upstream chain rule needs, e.g. to score STRUCTURAL
#'   parameters, which reach this likelihood only through \code{Theta}. It
#'   requires no callback and no model: hand it whatever produced
#'   \code{Theta_list}.
#' @param theta_fn Required when \code{"rho" \%in\% wrt}: a function
#'   \code{function(shock, rho_z)} returning the \code{T_h x n_obs}
#'   unit-innovation MA-coefficient matrix for \code{shock} at persistence
#'   \code{rho_z} (i.e. what \code{Theta_list[[shock]]} would be if the
#'   model were resolved at that persistence). \code{Theta_z} is DEFINED as
#'   the response to the driving path \code{rho_z^t}, so unlike the sigma/me
#'   blocks this cannot be derived from \code{Y}/\code{Theta_list} alone --
#'   see \code{\link{make_posterior_grad_hank_ar}} for the
#'   \code{hank_model_irf()}-based callback used when \code{ss} carries a
#'   \code{\link{hank_model}}.
#' @param h_rho Central-difference step for the rho-block slab FD
#'   (\code{rho_method = "fd_slab"} only).
#' @param rho_method How to obtain the rho score. \code{"fd_slab"} (default)
#'   central-differences the per-shock slab through \code{theta_fn}, keeping the
#'   stacked-covariance algebra analytic. \code{"adjoint"} is fully analytic:
#'   the kernel adjoint gives \code{d loglik/d Theta_z} and the partial
#'   \code{d loglik/d rho_z} at fixed \code{Theta}, and the total derivative
#'   adds \code{<d loglik/d Theta_z, d Theta_z/d rho_z>} using
#'   \code{dtheta_fn}. It needs no differencing step, but costs the
#'   block-Toeplitz adjoint pass; see the file header for measured costs.
#' @param dtheta_fn Required when \code{rho_method = "adjoint"} and
#'   \code{"rho" \%in\% wrt}: \code{function(shock, rho_z)} returning
#'   \code{d Theta_z / d rho_z} (\code{T_h x n_obs}). \code{Theta_z} is LINEAR
#'   in its driving path \code{rho_z^t}, so this is the same linear map applied
#'   to \code{d/d rho (rho^t) = t rho^(t-1)} -- one application, not a re-solve
#'   (see \code{\link{make_posterior_grad_hank_ar}}, which supplies it).
#' @param method \code{"compact"} (default) contracts every parameter against
#'   the compiled compact weights (one pass over the stacked entries, then a
#'   \code{n_obs^2 * T_data} dot product per parameter -- see
#'   \code{src/hank_ar_score.cpp}); \code{"gather"} is the original
#'   reference path (one \code{n_kept x n_kept} gather, trace and matvec per
#'   parameter). Mathematically identical; \code{"gather"} is retained as the
#'   parity oracle and costs ~4x more at paper scale.
#' @param use_cpp Logical: use the compiled weight accumulation
#'   (\code{method = "compact"} only). \code{FALSE} selects the
#'   \code{rowsum()}-based R fallback, which is the parity reference and much
#'   slower.
#'
#' @return A list: \code{loglik} (identical to what
#'   \code{\link{hank_loglik_ar}} would return, computed as part of the same
#'   call so the two are always consistent), and \code{sigma}/\code{rho}/
#'   \code{me}, each a named numeric vector (by shock, by shock, by
#'   observable respectively) or \code{NULL} if not requested via \code{wrt},
#'   plus \code{theta} (a named list of \code{q x n_obs} matrices) when
#'   \code{"theta" \%in\% wrt}.
#' @seealso \code{\link{hank_loglik_ar}}, \code{\link{make_posterior_grad_hank_ar}}
#' @export
hank_loglik_ar_grad <- function(Y, ss, rho = NULL, sigma = NULL, me_sd = 0,
                                q = NULL, check_boundary = TRUE,
                                boundary_tol = 1e-3, cache = NULL,
                                wrt = c("sigma", "rho", "me"),
                                theta_fn = NULL, h_rho = 1e-5,
                                method = c("compact", "gather"),
                                use_cpp = TRUE,
                                rho_method = c("fd_slab", "adjoint"),
                                dtheta_fn = NULL) {
  wrt <- match.arg(wrt, c("sigma", "rho", "me", "theta"), several.ok = TRUE)
  method <- match.arg(method)
  rho_method <- match.arg(rho_method)
  if ("rho" %in% wrt && rho_method == "fd_slab" && is.null(theta_fn))
    stop("hank_loglik_ar_grad: `theta_fn` is required when 'rho' %in% wrt ",
         "(Theta_z is DEFINED as the response to the rho_z^t driving path, ",
         "so the rho score needs a way to re-derive it at a perturbed ",
         "persistence -- see the `theta_fn` documentation).")
  if ("rho" %in% wrt && rho_method == "adjoint" && is.null(dtheta_fn))
    stop("hank_loglik_ar_grad: `dtheta_fn` is required when ",
         "rho_method = \"adjoint\" and 'rho' %in% wrt: the kernel adjoint ",
         "gives the rho score at FIXED Theta, and Theta_z itself depends on ",
         "rho_z (it is the response to the rho_z^t driving path), so the ",
         "total derivative needs dTheta_z/drho_z -- see the `dtheta_fn` ",
         "documentation.")
  if (is.null(cache)) cache <- new.env(parent = emptyenv())
  if (!is.environment(cache))
    stop("hank_loglik_ar_grad: `cache` must be an environment or NULL.")

  ## ONE preparation, shared verbatim with hank_loglik_ar(): it validates and
  ## normalizes (ss, rho, sigma, q) under that function's contract and error
  ## text, warms cache$A (every slab), and returns the summed autocovariance
  ## array. The loglik this IS the gradient of then comes from the SAME
  ## Cholesky the score blocks need, rather than from a second gather and
  ## factorization inside a separate hank_loglik_ar() call -- bit-identical by
  ## construction, since .hank_stacked_loglik() evaluates the same expression
  ## on the same (ch, yk).
  Y <- as.matrix(Y)
  prep <- .hank_ar_prepare(Y, ss, rho = rho, sigma = sigma, me_sd = me_sd,
                           q = q, check_boundary = check_boundary,
                           boundary_tol = boundary_tol, cache = cache)
  Theta_list <- prep$Theta_list
  rho <- prep$rho; sigma <- prep$sigma; me_sd <- prep$me_sd
  shocks <- prep$shocks
  Td <- nrow(Y); n_obs <- ncol(Y)
  q <- min(prep$q, nrow(as.matrix(Theta_list[[1L]])))
  G <- prep$G

  ## The gather index depends only on (T_data, n_obs, NA pattern), so a
  ## `cache` carried across calls keeps it through a STRUCTURAL move (which
  ## invalidates every slab but not this) -- see the cache documentation.
  ix <- .hank_stacked_index(Y, cache)
  keep <- ix$keep
  yk <- ix$yv[keep]; N <- length(yk)
  gath <- function(X) matrix(as.numeric(X)[ix$idx], N, N)
  S <- gath(G); S <- (S + t(S)) / 2
  ch <- chol(S)
  ll0 <- .hank_stacked_loglik_from_chol(ch, yk)
  z2 <- backsolve(ch, yk, transpose = TRUE)
  v <- backsolve(ch, z2)
  Sinv <- chol2inv(ch)

  ## One compact-weight pass replaces every per-parameter gather/trace/matvec
  ## (see .hank_ar_score_weights); dterm() below dispatches on `method` so the
  ## reference path stays available as a parity oracle.
  ## The adjoint blocks ("theta", and the rho score under rho_method =
  ## "adjoint") are contractions against the SAME weights, so build them
  ## whenever any consumer needs them -- not only when method == "compact".
  need_w <- method == "compact" || "theta" %in% wrt ||
    ("rho" %in% wrt && rho_method == "adjoint")
  w <- if (need_w)
    .hank_ar_score_weights(Sinv, v, ix$idx, length(G), use_cpp = use_cpp)
  else NULL
  dterm <- if (method == "compact")
    function(dG) .hank_ar_dterm_compact(w, dG)
  else function(dG) .hank_ar_dterm(Sinv, v, gath(dG))

  d_sigma <- d_rho <- d_me <- d_theta <- NULL

  ## ONE adjoint pass per shock, shared by the theta and rho blocks. The scalar
  ## it differentiates is sum_k <Wz[,,k], A_z[,,k]>, and loglik depends on A_z
  ## only through G += sigma_z^2 A_z, so Wz = -0.5 sigma_z^2 (compact weights).
  want_adj_theta <- "theta" %in% wrt
  want_adj_rho <- "rho" %in% wrt && rho_method == "adjoint"
  adj_list <- if (want_adj_theta || want_adj_rho)
    stats::setNames(lapply(shocks, function(z) {
      Theta <- as.matrix(Theta_list[[z]])
      .hank_ar_slab_adjoint(array(-0.5 * sigma[[z]]^2 * w, c(n_obs, n_obs, Td)),
                            Theta, rho[[z]], min(q, nrow(Theta)), Td - 1L,
                            want_theta = want_adj_theta || want_adj_rho,
                            want_rho = want_adj_rho, use_cpp = use_cpp)
    }), shocks)
  else NULL
  if (want_adj_theta)
    d_theta <- stats::setNames(lapply(shocks, function(z) adj_list[[z]]$theta),
                               shocks)

  if ("sigma" %in% wrt)
    d_sigma <- stats::setNames(vapply(shocks, function(z) {
      dterm(2 * sigma[[z]] * cache$A[[z]]$A)
    }, numeric(1)), shocks)

  if ("me" %in% wrt) {
    bI <- rep(seq_len(n_obs), Td)[keep]
    d_me <- stats::setNames(.hank_ar_dscore_me(Sinv, v, me_sd, bI),
                            colnames(Y) %||% paste0("obs", seq_len(n_obs)))
  }

  if ("rho" %in% wrt) {
    n_lags <- Td - 1L
    d_rho <- stats::setNames(rep(NA_real_, length(shocks)), shocks)
    for (z in shocks) {
      r0 <- rho[[z]]
      if (!is.finite(r0)) next
      if (rho_method == "adjoint") {
        ## total derivative = partial at fixed Theta (from the kernel adjoint)
        ## + the Theta channel, since Theta_z IS the response to rho_z^t.
        dTh <- as.matrix(dtheta_fn(z, r0))[seq_len(q), , drop = FALSE]
        d_rho[[z]] <- adj_list[[z]]$rho + sum(adj_list[[z]]$theta * dTh)
      } else {
        if (abs(r0) >= 1 - 2 * h_rho) next    # FD stencil cannot straddle 1
        Thp <- as.matrix(theta_fn(z, r0 + h_rho))[seq_len(q), , drop = FALSE]
        Thm <- as.matrix(theta_fn(z, r0 - h_rho))[seq_len(q), , drop = FALSE]
        Ap <- .hank_ar_autocov_kernel(Thp, r0 + h_rho, n_lags, q)
        Am <- .hank_ar_autocov_kernel(Thm, r0 - h_rho, n_lags, q)
        dA <- (Ap - Am) / (2 * h_rho)
        d_rho[[z]] <- dterm(sigma[[z]]^2 * dA)
      }
    }
  }

  list(loglik = ll0, sigma = d_sigma, rho = d_rho, me = d_me, theta = d_theta)
}


#' Score of the exact-AR likelihood w.r.t. STRUCTURAL model parameters
#'
#' Structural parameters (Taylor-rule coefficients, an NKPC slope, household
#' preferences) reach this likelihood only through the MA coefficients, so the
#' score factors as \eqn{dl/d\theta_k = \sum_z \langle dl/d\Theta_z,
#' d\Theta_z/d\theta_k\rangle}. The first factor is exact and shared by every
#' parameter (the kernel adjoint, \code{\link{hank_loglik_ar_grad}} with
#' \code{wrt = "theta"}); only the second is per-parameter, and it costs a
#' model rebuild rather than a likelihood evaluation.
#'
#' That split is the point: a plain central difference of the likelihood pays
#' two model rebuilds AND two \emph{cold} likelihood evaluations per parameter
#' (a structural move invalidates every cached autocovariance slab), while this
#' pays one adjoint pass for the whole gradient plus two rebuilds per
#' parameter -- and the rebuilds are themselves much cheaper than they look,
#' because the block-Jacobian cache skips the fake-news Jacobian for any
#' parameter that cannot touch the household. Measured scale and the full
#' rationale: \code{briefs/21-structural-score-api-scope.md}.
#'
#' @section Supplying derivatives, and why nothing is required: with
#'   \code{dtheta_fn = NULL} (default) \eqn{d\Theta_z/d\theta_k} is obtained by
#'   central differences of \code{model_fn} -- correct for every parameter,
#'   including ones entering a heterogeneous-agent block through a nonlinear
#'   steady state, at the cost of a differencing step. Supplying
#'   \code{dtheta_fn} makes it exact. A supplied \code{dtheta_fn} is VERIFIED
#'   against one finite-difference directional derivative per parameter unless
#'   \code{verify = FALSE}, because a misplaced derivative produces a
#'   plausible-but-wrong gradient rather than an error; pass
#'   \code{verify = FALSE} inside a sampler loop once it has passed, and the
#'   function says once, via \code{message()}, which parameters it is trusting.
#'
#' @section Reusing what the caller already has: the single largest cost in a
#'   structural gradient is not the adjoint or the chain rule but the ONE model
#'   rebuild this needs for the base \eqn{\Theta} (measured: 0.46 s at 11
#'   shocks x 12 observables, \code{T_h = 200}, against 0.008 s for the DAG
#'   propagation). A sampler or a block-coordinate optimizer already holds that
#'   model, so pass it: \code{model} (a \code{\link{hank_model}} at
#'   \code{theta}) or \code{Theta_list} (its MA coefficients, skipping the
#'   \code{.hank_theta_list} step as well) removes the rebuild entirely.
#'   \code{\link{hank_dtheta_fn}} takes the same \code{model}/\code{theta}
#'   seeding, so the exact route can be driven with ZERO rebuilds per gradient.
#'   Passing a \code{cache} environment across calls matters for the same
#'   reason on the likelihood side: the stacked gather index survives a
#'   structural move (only the autocovariance slabs are invalidated), and
#'   rebuilding it costs more than every slab put together.
#'
#'   Precedence and disagreement: \code{Theta_list} beats \code{model}, which
#'   beats \code{model_fn(theta)}. Supplying BOTH \code{model} and
#'   \code{Theta_list} is checked (one \code{.hank_theta_list} call, ~5% of a
#'   rebuild) and a mismatch is an error. A prebuilt base is NOT checked
#'   against \code{model_fn(theta)} by default, because that check costs
#'   exactly the rebuild being skipped; \code{verify = TRUE} does perform it,
#'   on the same "pay for it once, then turn it off" logic as the
#'   \code{dtheta_fn} verification.
#'
#' @param Y \code{T_data x n_obs} matrix of demeaned observations.
#' @param model_fn \code{function(theta)} returning a \code{\link{hank_model}}
#'   at the structural parameter vector \code{theta}. May be \code{NULL} only
#'   when a prebuilt base (\code{model} or \code{Theta_list}) and an exact
#'   \code{dtheta_fn} are supplied with \code{verify = FALSE} -- the two paths
#'   that would otherwise call it.
#' @param theta Named numeric vector of structural parameters.
#' @param model Optional prebuilt \code{\link{hank_model}} at \code{theta},
#'   used instead of calling \code{model_fn(theta)} for the base
#'   \eqn{\Theta_z}.
#' @param score Optional precomputed kernel adjoint: the list returned by
#'   \code{\link{hank_loglik_ar_grad}} with \code{"theta"} in \code{wrt}, at
#'   this same \code{(Theta, rho, sigma)}. The adjoint is shared by every
#'   structural parameter AND by the \code{rho}/\code{sigma} scores, so a
#'   caller estimating both blocks together has already paid for it; passing it
#'   here removes the second adjoint pass per gradient. Nothing can verify it
#'   corresponds to the current arguments -- that is the caller's obligation.
#' @param Theta_list Optional prebuilt named list of \code{T_h x n_obs} MA
#'   coefficient matrices at \code{theta} (what
#'   \code{model_fn(theta)} would produce for \code{observables} under
#'   \code{rho}), used instead of \code{model} / \code{model_fn}.
#' @param observables Character vector of observable names, in the column order
#'   of \code{Y}.
#' @param rho,sigma Named numeric vectors of per-shock persistence and
#'   innovation sd, held FIXED here (their own scores come from
#'   \code{\link{hank_loglik_ar_grad}}).
#' @param dtheta_fn Optional \code{function(theta, k)} returning
#'   \eqn{d\Theta_z/d\theta_k} for parameter name \code{k}: a named list of
#'   \code{T_h x n_obs} matrices, one per shock. \code{NULL} (default) obtains
#'   it by central differences of \code{model_fn}.
#' @param verify Logical: when \code{dtheta_fn} is supplied, check it against
#'   one central difference per parameter (relative tolerance \code{tol}) and
#'   error on a mismatch; and, when a prebuilt \code{model}/\code{Theta_list}
#'   is supplied alongside \code{model_fn}, check that base against
#'   \code{model_fn(theta)} (one extra rebuild). Ignored when \code{dtheta_fn}
#'   is \code{NULL} and no prebuilt base is supplied.
#' @param tol Relative tolerance for that verification.
#' @param h Central-difference step for \code{theta} (used by the FD route and
#'   by \code{verify}).
#' @inheritParams hank_loglik_ar_grad
#'
#' @return A list with \code{loglik}, \code{structural} (a named numeric vector
#'   matching \code{theta}), and \code{theta} (the per-shock
#'   \eqn{dl/d\Theta_z} matrices the score was contracted from).
#' @seealso \code{\link{hank_loglik_ar_grad}}, \code{\link{hank_loglik_ar}}
#' @export
hank_loglik_ar_structural_grad <- function(Y, model_fn, theta, observables,
                                           rho, sigma, me_sd = 0, q = NULL,
                                           dtheta_fn = NULL, verify = TRUE,
                                           tol = 1e-4, h = 1e-5, cache = NULL,
                                           check_boundary = FALSE,
                                           boundary_tol = 1e-3,
                                           model = NULL, Theta_list = NULL,
                                           score = NULL) {
  base_arg <- if (!is.null(Theta_list)) "Theta_list"
              else if (!is.null(model)) "model" else NULL
  have_base <- !is.null(base_arg)
  ## model_fn is needed by the FD-dTheta route and by every verification; it
  ## is dispensable only when neither can run.
  need_fn <- is.null(dtheta_fn) || isTRUE(verify) || !have_base
  if (is.null(model_fn) && !need_fn) model_fn <- function(th) NULL
  if (!is.function(model_fn))
    stop("hank_loglik_ar_structural_grad: `model_fn` must be a function of ",
         "theta returning a hank_model()", if (have_base)
           paste(" -- it stays required for the finite-difference dTheta",
                 "route and for verify = TRUE, so supply dtheta_fn with",
                 "verify = FALSE to drop it") else "", ".")
  nms <- names(theta)
  if (is.null(nms) || any(!nzchar(nms)))
    stop("hank_loglik_ar_structural_grad: `theta` must be a NAMED numeric ",
         "vector (the score is returned under the same names).")
  if (!is.null(model) && !inherits(model, "hank_model"))
    stop("hank_loglik_ar_structural_grad: `model` must be a hank_model() ",
         "object built at `theta`, or NULL.")
  Y <- as.matrix(Y)

  shock_specs <- stats::setNames(
    lapply(names(rho), function(z) list(rho = rho[[z]], sigma = sigma[[z]])),
    names(rho))
  theta_list_of <- function(mod) {
    if (!inherits(mod, "hank_model"))
      stop("hank_loglik_ar_structural_grad: `model_fn` must return a ",
           "hank_model() object.")
    .hank_theta_list(mod, shock_specs, observables)
  }
  theta_list_at <- function(th) theta_list_of(model_fn(th))

  ## Base Theta: Theta_list beats model beats model_fn(theta). Only the last
  ## costs a rebuild, which is the whole point of the first two.
  agree <- function(a, b, what) {
    ok <- identical(names(a), names(b)) && all(vapply(names(a), function(z) {
      A <- as.matrix(a[[z]]); B <- as.matrix(b[[z]])
      identical(dim(A), dim(B)) &&
        max(abs(A - B)) <= tol * max(max(abs(B)), 1e-8)
    }, logical(1)))
    if (!ok)
      stop("hank_loglik_ar_structural_grad: the supplied `", what[1L],
           "` disagrees with ", what[2L], " (relative tolerance ",
           format(tol), "). A base Theta from a different theta yields a ",
           "plausible-but-wrong gradient, so this is an error, not a warning.")
  }
  if (!is.null(Theta_list)) {
    if (!is.list(Theta_list) || !all(names(rho) %in% names(Theta_list)))
      stop("hank_loglik_ar_structural_grad: `Theta_list` must be a list with ",
           "an entry for every shock in {",
           paste(names(rho), collapse = ", "), "}.")
    Theta_list <- Theta_list[names(rho)]
    ## both supplied: cross-check, which costs one .hank_theta_list (~5% of a
    ## rebuild) -- cheap enough to be worth catching a mismatched pair.
    if (!is.null(model)) agree(Theta_list, theta_list_of(model),
                               c("Theta_list", "the model's own Theta"))
  } else if (!is.null(model)) {
    Theta_list <- theta_list_of(model)
  } else {
    Theta_list <- theta_list_at(theta)
  }
  ## A prebuilt base is trusted by default (checking it costs the rebuild it
  ## saves); verify = TRUE buys the check, like the dtheta_fn one.
  if (have_base && isTRUE(verify))
    agree(Theta_list, theta_list_at(theta), c(base_arg, "model_fn(theta)"))
  ## The kernel adjoint is shared by every structural parameter AND by the
  ## rho/sigma scores, so a caller estimating both blocks together (the runner
  ## in R/hank-run-ar.R) has already computed it: `score` lets it be reused
  ## instead of paying a second adjoint pass per gradient. It must be the
  ## output of hank_loglik_ar_grad() with "theta" in `wrt`, at THIS (Theta,
  ## rho, sigma) -- nothing here can check that, so it is documented as the
  ## caller's obligation and defaults to computing it.
  sc <- if (is.null(score)) {
    hank_loglik_ar_grad(Y, Theta_list, rho = rho, sigma = sigma,
                        me_sd = me_sd, q = q, cache = cache,
                        check_boundary = check_boundary,
                        boundary_tol = boundary_tol, wrt = "theta")
  } else {
    if (!is.list(score) || is.null(score$theta) || is.null(score$loglik))
      stop("hank_loglik_ar_structural_grad: `score` must be the list returned ",
           "by hank_loglik_ar_grad() with \"theta\" in `wrt` (it needs both ",
           "$theta and $loglik).")
    score
  }
  shocks <- names(sc$theta)
  q_eff <- nrow(sc$theta[[1L]])

  ## d Theta_z / d theta_k by central differences of the model builder. Both
  ## taps go through the block-Jacobian cache, so a parameter that leaves the
  ## het block untouched never re-runs the fake-news algorithm.
  fd_dtheta <- function(k) {
    tp <- theta; tp[[k]] <- tp[[k]] + h
    tm <- theta; tm[[k]] <- tm[[k]] - h
    Lp <- theta_list_at(tp); Lm <- theta_list_at(tm)
    stats::setNames(lapply(shocks, function(z)
      (as.matrix(Lp[[z]]) - as.matrix(Lm[[z]])) / (2 * h)), shocks)
  }
  contract <- function(dTh)
    sum(vapply(shocks, function(z)
      sum(sc$theta[[z]] * as.matrix(dTh[[z]])[seq_len(q_eff), , drop = FALSE]),
      numeric(1)))

  if (is.null(dtheta_fn)) {
    grad <- vapply(nms, function(k) contract(fd_dtheta(k)), numeric(1))
  } else {
    if (!isTRUE(verify))
      message("hank_loglik_ar_structural_grad: trusting the supplied ",
              "dtheta_fn unverified for ", paste(nms, collapse = ", "),
              " (verify = TRUE checks each against one central difference).")
    grad <- vapply(nms, function(k) {
      dTh <- dtheta_fn(theta, k)
      if (!is.list(dTh) || !all(shocks %in% names(dTh)))
        stop("hank_loglik_ar_structural_grad: dtheta_fn(theta, \"", k, "\") ",
             "must return a list with an entry for every shock in {",
             paste(shocks, collapse = ", "), "}.")
      if (isTRUE(verify)) {
        ref <- fd_dtheta(k)
        for (z in shocks) {
          a <- as.matrix(dTh[[z]])[seq_len(q_eff), , drop = FALSE]
          b <- ref[[z]][seq_len(q_eff), , drop = FALSE]
          rel <- max(abs(a - b)) / max(max(abs(b)), 1e-8)
          if (!is.finite(rel) || rel > tol)
            stop("hank_loglik_ar_structural_grad: the supplied dtheta_fn ",
                 "disagrees with a central difference of `model_fn` for ",
                 "parameter '", k, "', shock '", z, "' (relative ",
                 sprintf("%.2e", rel), " > tol ", format(tol), "). A ",
                 "misplaced derivative yields a plausible-but-wrong ",
                 "gradient, so this is an error, not a warning; pass ",
                 "verify = FALSE to bypass deliberately.")
        }
      }
      contract(dTh)
    }, numeric(1))
  }

  list(loglik = sc$loglik, structural = grad, theta = sc$theta)
}


#' Gradient of \code{\link{make_log_posterior_hank}}'s exact-AR posterior
#'
#' \code{\link{make_log_posterior_hank}}-shaped gradient closure: same
#' arguments, same \code{(rho_<shock>, sigma_<shock>)} parameterization, same
#' priors, but \code{likelihood = "exact_ar"} only (the Kalman path already
#' has a full gradient stack via \code{\link{make_posterior_grad}}) and
#' every call returns \code{list(logpost, loglik, logprior, grad)} where
#' \code{grad} is a named numeric vector matching \code{theta}.
#'
#' Reuses \code{\link{hank_loglik_ar_grad}} with a \code{theta_fn} built from
#' \code{\link{hank_model_irf}} for the rho block (cheap: \code{model$G} is
#' the sequence-space Jacobian, already assembled once by
#' \code{\link{hank_model}}, so this is a small matvec per FD tap, not a
#' re-solve -- the same cost \code{\link{hank_state_space}}'s own
#' \code{Theta_list} construction already pays per shock).
#'
#' @inheritParams make_log_posterior_hank
#' @param rho_method Passed to \code{\link{hank_loglik_ar_grad}}; with
#'   \code{"adjoint"} the required \code{dtheta_fn} is supplied here exactly
#'   (the IRF of the differentiated driving path).
#' @param boundary,boundary_tol Representability guard, as for
#'   \code{\link{make_log_posterior_hank}}. It matters MORE here: a
#'   gradient-based sampler follows the score, and the score of a
#'   terminal-boundary-contaminated \code{Theta} points somewhere, so nothing
#'   stops the chain walking further out. \code{"reject"} returns
#'   \code{logpost = -Inf} with \code{grad = NULL}, which every sampler in this
#'   package already handles as an infeasible draw.
#' @return A function \code{grad_fn(theta)} returning
#'   \code{list(logpost, loglik, logprior, grad)}; \code{grad} is
#'   \code{NULL} at an infeasible \code{theta} (mirrors
#'   \code{\link{make_log_posterior_hank}}'s \code{-Inf} early return).
#' @seealso \code{\link{make_log_posterior_hank}}, \code{\link{hank_loglik_ar_grad}}
#' @export
make_posterior_grad_hank_ar <- function(model, Y, observables, q = NULL,
                                        me_var = 0, me_sd = NULL,
                                        prior_rho_mean = 0.5, prior_rho_sd = 0.3,
                                        prior_sigma_sd = 0.05,
                                        rho_method = c("fd_slab", "adjoint"),
                                        boundary = c("warn", "reject",
                                                     "ignore"),
                                        boundary_tol = 1e-3) {
  rho_method <- match.arg(rho_method)
  boundary <- match.arg(boundary)
  if (!(is.numeric(boundary_tol) && length(boundary_tol) == 1L &&
        is.finite(boundary_tol) && boundary_tol > 0))
    stop("make_posterior_grad_hank_ar: `boundary_tol` must be a finite ",
         "positive scalar.")
  boundary_state <- new.env(parent = emptyenv())
  exo <- model$exogenous
  missing_obs <- setdiff(observables, names(model$G))
  if (length(missing_obs))
    stop("make_posterior_grad_hank_ar: observable(s) not produced by model: ",
         paste(missing_obs, collapse = ", "))
  if (is.null(me_sd)) me_sd <- sqrt(me_var)
  rep_named <- function(x, nm) {
    if (is.null(names(x))) stats::setNames(rep(x, length.out = length(nm)), nm)
    else x[nm]
  }
  rho_mean <- rep_named(prior_rho_mean, exo)
  rho_sd   <- rep_named(prior_rho_sd, exo)
  sig_sd   <- rep_named(prior_sigma_sd, exo)
  ar_cache <- new.env(parent = emptyenv())

  irf_of <- function(z, path) {
    irf <- hank_model_irf(model, stats::setNames(list(path), z))
    matrix(sapply(observables, function(o) irf[[o]]),
          model$T_h, length(observables), dimnames = list(NULL, observables))
  }
  theta_fn <- function(z, r) irf_of(z, r^(seq_len(model$T_h) - 1L))
  ## dTheta_z/drho_z for rho_method = "adjoint": Theta_z is LINEAR in its
  ## driving path, so differentiating the path is exact -- the same IRF applied
  ## to d/drho (rho^s) = s rho^(s-1), one application, not a re-solve.
  dtheta_fn <- function(z, r) {
    s <- seq_len(model$T_h) - 1L
    irf_of(z, ifelse(s == 0, 0, s * r^pmax(s - 1L, 0)))
  }

  function(theta) {
    rho_nm <- paste0("rho_", exo); sig_nm <- paste0("sigma_", exo)
    if (!all(c(rho_nm, sig_nm) %in% names(theta)))
      stop("make_posterior_grad_hank_ar: theta must have entries ",
           paste(c(rho_nm, sig_nm), collapse = ", "))
    rho   <- theta[rho_nm];   names(rho)   <- exo
    sigma <- theta[sig_nm];   names(sigma) <- exo

    if (any(rho <= -1 | rho >= 1) || any(sigma <= 0))
      return(list(logpost = -Inf, loglik = NA_real_, logprior = -Inf, grad = NULL))

    logprior <- sum(stats::dnorm(rho, rho_mean, rho_sd, log = TRUE)) +
      sum(stats::dnorm(sigma, 0, sig_sd, log = TRUE) + log(2))
    if (.hank_ar_boundary_gate(rho, model$T_h, boundary, boundary_tol,
                               "make_posterior_grad_hank_ar", boundary_state))
      return(list(logpost = -Inf, loglik = NA_real_, logprior = logprior,
                  grad = NULL))
    dlp_rho   <- -(rho - rho_mean) / rho_sd^2
    dlp_sigma <- -sigma / sig_sd^2

    shock_specs <- stats::setNames(
      lapply(exo, function(z) list(rho = rho[[z]], sigma = sigma[[z]])), exo)
    Theta_list <- .hank_theta_list(model, shock_specs, observables)

    sc <- tryCatch(
      hank_loglik_ar_grad(Y, Theta_list, rho = rho, sigma = sigma, me_sd = me_sd,
                          q = q, check_boundary = FALSE, cache = ar_cache,
                          wrt = c("sigma", "rho"), theta_fn = theta_fn,
                          rho_method = rho_method, dtheta_fn = dtheta_fn),
      error = function(e) NULL)
    if (is.null(sc) || !is.finite(sc$loglik))
      return(list(logpost = -Inf, loglik = if (is.null(sc)) NA_real_ else sc$loglik,
                  logprior = logprior, grad = NULL))

    grad <- stats::setNames(numeric(length(theta)), names(theta))
    grad[rho_nm] <- sc$rho[exo] + dlp_rho[exo]
    grad[sig_nm] <- sc$sigma[exo] + dlp_sigma[exo]

    list(logpost = sc$loglik + logprior, loglik = sc$loglik,
        logprior = logprior, grad = grad)
  }
}
