## R/svar-ica.R
## --------------------------------------------------------------------------
## Non-Gaussian SVAR identified by independent component analysis (ICA), and
## indirect inference (II) that matches its structural IRFs.
##
## svar_ica()             reduced-form VAR(p) by OLS (estimate_var()), whiten
##                        the residuals with the Cholesky factor of Sigma_u,
##                        then find the orthogonal rotation O whose rotated
##                        components are least dependent (distance-covariance
##                        ICA of Matteson & Tsay 2017) or maximally non-Gaussian
##                        (symmetric FastICA).  Impact matrix B = D O'.
## match_irfs_svar_ica()  the IRF-matching indirect-inference estimator of
##                        Martinoli, Di Francesco, Moneta & Seri (2026, LEM WP
##                        2026/19): the SAME SVAR-ICA auxiliary model is fitted
##                        to the data and to S model simulations (common random
##                        numbers, the model's own skewed CSN shocks), the
##                        simulated impact matrices are aligned to the data one
##                        by their MD1 minimum-distance index, and theta
##                        minimises || psi_hat - mean_s psi_s(theta) ||_W^2.
##
## Base R (stats) only; no ICA package.
## --------------------------------------------------------------------------


## --------------------------------------------------------------------------
## Rotation parameterisation: product of Givens rotations, one angle per pair
## (i < j).  K(K-1)/2 angles cover SO(K); reflections are irrelevant because
## the ICA contrast is invariant to column signs.
## --------------------------------------------------------------------------
.ica_givens <- function(omega, K) {
  O <- diag(K)
  k <- 0L
  for (i in seq_len(K - 1L)) {
    for (j in (i + 1L):K) {
      k <- k + 1L
      c0 <- cos(omega[k]); s0 <- sin(omega[k])
      G <- diag(K)
      G[i, i] <- c0; G[j, j] <- c0
      G[i, j] <- -s0; G[j, i] <- s0
      O <- G %*% O
    }
  }
  O
}


## --------------------------------------------------------------------------
## U-centred distance matrix (Szekely & Rizzo 2014; Huo & Szekely 2016):
##   A~_ij = a_ij - a_i./(n-2) - a_.j/(n-2) + a../((n-1)(n-2)),  A~_ii = 0.
## x is an n-vector or an n x m matrix (rows = observations).
## --------------------------------------------------------------------------
.dcov_ucentred <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (ncol(x) == 1L) {
    a <- abs(outer(x[, 1L], x[, 1L], "-"))
  } else {
    a <- as.matrix(stats::dist(x))
  }
  rs  <- rowSums(a)
  tot <- sum(rs)
  A <- a - outer(rs, rs, "+") / (n - 2) + tot / ((n - 1) * (n - 2))
  diag(A) <- 0
  A
}


## --------------------------------------------------------------------------
## Distance-covariance ICA contrast (Matteson & Tsay 2017, eq. (13) of
## Martinoli et al.): sum_{k=1}^{K-1} U_T(s_k, s_{(k+1):K}), with U_T the
## unbiased (U-statistic) squared distance covariance.  Zero in population iff
## the components are mutually independent.
## --------------------------------------------------------------------------
.ica_dcov_contrast <- function(S) {
  n <- nrow(S); K <- ncol(S)
  val <- 0
  for (k in seq_len(K - 1L)) {
    A <- .dcov_ucentred(S[, k])
    B <- .dcov_ucentred(S[, (k + 1L):K, drop = FALSE])
    val <- val + sum(A * B) / (n * (n - 3))
  }
  val
}


## --------------------------------------------------------------------------
## Symmetric FastICA (Hyvarinen 1999) on whitened data V (n x K), logcosh
## contrast g = tanh.  Returns the orthogonal unmixing rotation O (components
## S = V O').  Deterministic: starts from O_init.
## --------------------------------------------------------------------------
.ica_sym_decorrelate <- function(W) {
  e <- eigen(W %*% t(W), symmetric = TRUE)
  e$vectors %*% diag(1 / sqrt(e$values), nrow = length(e$values)) %*%
    t(e$vectors) %*% W
}

.ica_fastica <- function(V, O_init = diag(ncol(V)), maxit = 500L, tol = 1e-8) {
  n <- nrow(V)
  W <- .ica_sym_decorrelate(O_init)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    Y  <- V %*% t(W)
    g  <- tanh(Y)
    gp <- 1 - g^2
    W1 <- crossprod(g, V) / n - diag(colMeans(gp), nrow = ncol(V)) %*% W
    W1 <- .ica_sym_decorrelate(W1)
    delta <- max(abs(abs(diag(W1 %*% t(W))) - 1))
    W <- W1
    if (delta < tol) { converged <- TRUE; break }
  }
  list(O = W, converged = converged, iterations = it)
}


## --------------------------------------------------------------------------
## Distance-covariance rotation search.
##   K = 2: the contrast is pi/2-periodic in the single angle (a rotation by
##          pi/2 is a signed permutation), so a grid on [0, pi/2) followed by a
##          Brent refinement around the best node finds the global minimum.
##   K > 2: Nelder-Mead over the K(K-1)/2 Givens angles, relative to a set of
##          starting rotations: the FastICA solution plus `n_starts` random
##          Haar rotations (drawn under a local seed; the caller's RNG stream is
##          untouched).  The best local minimum is kept.
## --------------------------------------------------------------------------
.ica_dcov_rotation <- function(V, n_grid = 12L, n_starts = 3L, seed = 1L,
                               maxit = 400L) {
  K <- ncol(V)
  f_rot <- function(O) .ica_dcov_contrast(V %*% t(O))

  if (K == 2L) {
    grid <- seq(0, pi / 2, length.out = n_grid + 1L)[seq_len(n_grid)]
    vals <- vapply(grid, function(w) f_rot(.ica_givens(w, 2L)), numeric(1))
    step <- grid[2L] - grid[1L]
    w0   <- grid[which.min(vals)]
    opt  <- stats::optimize(function(w) f_rot(.ica_givens(w, 2L)),
                            lower = w0 - step, upper = w0 + step,
                            tol = 1e-6)
    best <- if (opt$objective < min(vals)) opt$minimum else w0
    return(list(O = .ica_givens(best, 2L), objective = min(opt$objective, min(vals))))
  }

  n_ang  <- K * (K - 1L) / 2L
  starts <- list(.ica_fastica(V)$O)
  if (n_starts > 0L) {
    rand <- .with_local_seed(seed, lapply(seq_len(n_starts), function(i) {
      q <- qr(matrix(stats::rnorm(K * K), K, K))
      Q <- qr.Q(q)
      Q %*% diag(sign(diag(qr.R(q))), nrow = K)
    }))
    starts <- c(starts, rand)
  }
  best <- list(O = starts[[1L]], objective = Inf)
  for (O0 in starts) {
    fn  <- function(w) f_rot(.ica_givens(w, K) %*% O0)
    opt <- stats::optim(rep(0, n_ang), fn, method = "Nelder-Mead",
                        control = list(maxit = maxit, reltol = 1e-10))
    if (opt$value < best$objective)
      best <- list(O = .ica_givens(opt$par, K) %*% O0, objective = opt$value)
  }
  best
}


## --------------------------------------------------------------------------
## Signed-permutation alignment.  For a K x K matrix G, the signed permutation
## C = P Lambda maximising tr(C G) -- equivalently minimising ||C G - I||_F,
## since ||C G - I||^2 = ||G||^2 + K - 2 tr(C G) for orthogonal C.  With
## C[j, perm[j]] = sign_j, tr(C G) = sum_j sign_j G[perm[j], j], so for a fixed
## permutation the optimal signs are sign(G[perm[j], j]) and only the K!
## permutations need enumerating (greedy assignment above K = 8).
## Returns list(C, perm, sign).
## --------------------------------------------------------------------------
.ica_best_signed_perm <- function(G) {
  K <- nrow(G)
  AG <- abs(G)
  if (K <= 8L) {
    perms <- .multiset_perms(seq_len(K))
    score <- vapply(perms, function(pp) sum(AG[cbind(pp, seq_len(K))]), numeric(1))
    perm  <- perms[[which.max(score)]]
  } else {
    perm <- integer(K)
    M <- AG
    for (r in seq_len(K)) {
      idx <- which(M == max(M), arr.ind = TRUE)[1L, ]
      perm[idx[2L]] <- idx[1L]
      M[idx[1L], ] <- -Inf; M[, idx[2L]] <- -Inf
    }
  }
  sg <- sign(G[cbind(perm, seq_len(K))])
  sg[sg == 0] <- 1
  C <- matrix(0, K, K)
  C[cbind(seq_len(K), perm)] <- sg
  list(C = C, perm = perm, sign = sg)
}


## --------------------------------------------------------------------------
## MD1 alignment (Martinoli et al., Definition 1): given an estimated impact
## matrix B (mixing, u = B e) and a reference impact matrix Psi_ref, choose the
## signed permutation C minimising ||C B^{-1} Psi_ref - I||_F and return the
## aligned impact matrix B C' (so that B C' ~ Psi_ref column by column) and
## the index D1 = min ||.||_F / sqrt(K - 1).
## --------------------------------------------------------------------------
.ica_align_md1 <- function(B, ref) {
  K <- nrow(B)
  G <- solve(B, ref)
  sp <- .ica_best_signed_perm(G)
  d1 <- sqrt(sum((sp$C %*% G - diag(K))^2)) / sqrt(max(K - 1L, 1L))
  list(B = B %*% t(sp$C), C = sp$C, D1 = d1)
}


## --------------------------------------------------------------------------
## Default normalisation (no reference): shock j is the one whose relative
## loading |B_jj| / ||B_.j|| is largest on variable j (the column permutation
## maximising sum_j |B_jj| / ||B_.j||), then each column's sign is set so that
## its diagonal entry is positive.
## --------------------------------------------------------------------------
.ica_normalise_default <- function(B) {
  Bn <- sweep(B, 2L, sqrt(colSums(B^2)), `/`)
  ## Column permutation: column perm[j] goes to position j, maximising
  ## sum_j |Bn[j, perm[j]]| = sum_j |t(Bn)[perm[j], j]| -- the signed-
  ## permutation problem above applied to t(Bn).
  sp <- .ica_best_signed_perm(t(Bn))
  Bp <- B[, sp$perm, drop = FALSE]
  s  <- sign(diag(Bp)); s[s == 0] <- 1
  sweep(Bp, 2L, s, `*`)
}


## --------------------------------------------------------------------------
## Per-component Jarque-Bera normality diagnostics on the structural shocks.
## --------------------------------------------------------------------------
.ica_gaussianity <- function(E, level) {
  n <- nrow(E)
  sk <- apply(E, 2L, function(v) { v <- v - mean(v); mean(v^3) / mean(v^2)^1.5 })
  ku <- apply(E, 2L, function(v) { v <- v - mean(v); mean(v^4) / mean(v^2)^2 - 3 })
  jb <- n / 6 * (sk^2 + ku^2 / 4)
  pv <- stats::pchisq(jb, df = 2, lower.tail = FALSE)
  data.frame(shock = colnames(E), skewness = sk, excess_kurtosis = ku,
             jarque_bera = jb, p_value = pv, gaussian = pv > level,
             row.names = NULL, stringsAsFactors = FALSE)
}


#' Non-Gaussian structural VAR identified by independent component analysis
#'
#' Fits a reduced-form VAR(\eqn{p}) by OLS (\code{\link{estimate_var}}) and
#' identifies the structural impact matrix statistically, by independent
#' component analysis (ICA) of the residuals, as in the SVAR-ICA literature
#' (Moneta et al. 2013; Gourieroux, Monfort & Renne 2017; Lanne, Meitz &
#' Saikkonen 2017) and the auxiliary model of Martinoli, Di Francesco, Moneta
#' & Seri (2026).  With \eqn{u_t = B \varepsilon_t}, \eqn{\varepsilon_t} having
#' mutually independent, unit-variance components of which at most one is
#' Gaussian, \eqn{B} is identified up to the order and signs of its columns.
#'
#' Algorithm: the residuals are whitened with the lower Cholesky factor
#' \eqn{D} of \eqn{\hat\Sigma_u} (\eqn{v_t = D^{-1}\hat u_t}; this imposes no
#' recursive structure), and the orthogonal rotation \eqn{O} is chosen so that
#' the components \eqn{O v_t} are least dependent.  Then \eqn{B = D O'},
#' \eqn{B B' = \hat\Sigma_u}, and the structural shocks are
#' \eqn{\hat\varepsilon_t = B^{-1}\hat u_t}.
#' \describe{
#'   \item{\code{method = "dcov"} (default)}{Minimises the distance-covariance
#'     contrast of Matteson & Tsay (2017),
#'     \eqn{\sum_{k=1}^{K-1} \mathcal U_T(\varepsilon_{k}, \varepsilon_{(k+1):K})}
#'     with \eqn{\mathcal U_T} the unbiased (U-statistic) squared distance
#'     covariance -- the estimator used by Martinoli et al.  Nonparametric: no
#'     assumption on the shape of the shock distributions.  The rotation is
#'     parameterised by \eqn{K(K-1)/2} Givens angles; for \eqn{K = 2} a grid
#'     plus Brent search finds the global minimum, for \eqn{K > 2} Nelder-Mead
#'     runs from the FastICA solution and \code{n_starts} random rotations.
#'     Cost is \eqn{O(T^2)} per contrast evaluation.}
#'   \item{\code{method = "fastica"}}{Symmetric FastICA (Hyvarinen 1999) with
#'     the log-cosh contrast: \eqn{O(T)} per iteration and much faster, but its
#'     contrast measures kurtosis, so it can lose power when the shocks are
#'     skewed yet close to mesokurtic.}
#' }
#'
#' @section Normalisation:
#' ICA identifies \eqn{B} only up to a signed permutation of its columns.
#' With \code{reference = NULL} the columns are ordered so that shock
#' \eqn{j} is the one loading most (relative to its column norm) on variable
#' \eqn{j} -- the permutation maximising
#' \eqn{\sum_j |B_{jj}| / \lVert B_{\cdot j}\rVert} -- and each column is
#' signed so that \eqn{B_{jj} > 0}.  With a \code{reference} impact matrix the
#' columns are instead aligned to it by the MD1 minimum-distance index of
#' Martinoli et al. (Definition 1): the signed permutation \eqn{C} minimising
#' \eqn{\lVert C B^{-1} \Psi_{ref} - I \rVert_F}; the index value is returned as
#' \code{md1}.
#'
#' @section Identification check:
#' Identification requires at most one Gaussian shock.  Each estimated
#' structural shock gets a Jarque-Bera test; when two or more are not
#' distinguishable from Gaussian at level \code{gauss_level}, the rotation
#' among them is not identified (every rotation of independent Gaussians is
#' again independent), \code{identified} is \code{FALSE}, and -- if
#' \code{warn_gaussian} -- a warning of class
#' \code{"dynhr_svar_ica_not_identified"} is raised.  This is a per-component
#' screen, not the rank test of Guay (2021) that Martinoli et al. recommend.
#'
#' @param data A \code{T x K} numeric matrix or data frame (rows = periods in
#'   order, columns = variables).
#' @param p VAR lag order.
#' @param horizon Number of IRF periods (index 1 = impact), as in
#'   \code{\link{var_irf}}.
#' @param type Deterministic terms of the VAR, as in \code{\link{estimate_var}}.
#' @param method \code{"dcov"} (default) or \code{"fastica"}; see Details.
#' @param reference Optional \code{K x K} reference impact matrix; when given,
#'   the columns of \eqn{B} are aligned to it by MD1 (see Normalisation).
#' @param gauss_level Level of the per-shock Jarque-Bera normality tests.
#'   Default \code{0.05}.
#' @param warn_gaussian Logical; warn when identification fails the Gaussianity
#'   screen.  Default \code{TRUE}.
#' @param n_starts Number of random starting rotations for the \code{"dcov"}
#'   search when \eqn{K > 2} (in addition to the FastICA start).  Default
#'   \code{3}.
#' @param seed Seed for those random starts (the caller's RNG stream is
#'   restored).  Default \code{1}.
#'
#' @return A list of class \code{svar_ica}:
#'   \describe{
#'     \item{B}{\code{K x K} structural impact matrix (\eqn{u_t = B\varepsilon_t};
#'       rows = variables, columns = shocks).}
#'     \item{rotation}{The orthogonal \eqn{O} with \eqn{B = D O'} (before the
#'       normalising signed permutation).}
#'     \item{chol}{The whitening Cholesky factor \eqn{D}.}
#'     \item{irfs}{Named list (\code{"shock1"}, ...) of \code{horizon x K}
#'       structural IRF matrices, the \code{\link{var_irf}} format.}
#'     \item{shocks}{\code{(T-p) x K} estimated structural shocks.}
#'     \item{objective}{ICA contrast at the optimum (distance covariance for
#'       \code{"dcov"}; the same distance-covariance contrast evaluated at the
#'       FastICA rotation for \code{"fastica"}).}
#'     \item{gaussianity}{Per-shock skewness, excess kurtosis, Jarque-Bera
#'       statistic and p-value.}
#'     \item{identified}{\code{TRUE} when at most one shock looks Gaussian.}
#'     \item{md1}{MD1 index to \code{reference} (\code{NA} without one).}
#'     \item{var_fit, method, p, horizon}{Bookkeeping.}
#'   }
#'
#' @references
#' Matteson, D. S. and Tsay, R. S. (2017). Independent component analysis via
#' distance covariance. \emph{Journal of the American Statistical Association},
#' 112(518), 623-637.
#'
#' Martinoli, M., Di Francesco, D., Moneta, A. and Seri, R. (2026). Estimation
#' of DSGE models by non-Gaussian vector autoregressions. LEM Working Paper
#' 2026/19, Sant'Anna School of Advanced Studies.
#'
#' Hyvarinen, A. (1999). Fast and robust fixed-point algorithms for independent
#' component analysis. \emph{IEEE Transactions on Neural Networks}, 10(3),
#' 626-634.
#'
#' @examples
#' set.seed(1)
#' Tn <- 400
#' e  <- cbind(rexp(Tn) - 1, sample(c(-1, 1), Tn, TRUE) * rexp(Tn) / sqrt(2))
#' B0 <- matrix(c(1, 0.5, -0.4, 1), 2, 2)
#' Y  <- matrix(0, Tn, 2)
#' for (t in 2:Tn) Y[t, ] <- 0.5 * Y[t - 1, ] + B0 %*% e[t, ]
#' fit <- svar_ica(Y, p = 1, horizon = 8)
#' fit$B          # ~ B0 up to column order / sign
#'
#' @seealso \code{\link{estimate_var}}, \code{\link{var_irf}},
#'   \code{\link{match_irfs_svar_ica}}
#' @export
svar_ica <- function(data, p = 1L, horizon = 20L,
                     type = c("const", "none", "trend"),
                     method = c("dcov", "fastica"),
                     reference = NULL, gauss_level = 0.05,
                     warn_gaussian = TRUE, n_starts = 3L, seed = 1L) {
  type   <- match.arg(type)
  method <- match.arg(method)
  horizon <- as.integer(horizon)
  if (length(horizon) != 1L || is.na(horizon) || horizon < 1L)
    .dynhr_abort("svar_ica: `horizon` must be a positive integer.",
                 class = "dynhr_input_error")
  Y <- as.matrix(data)
  if (!is.numeric(Y) || !all(is.finite(Y)))
    .dynhr_abort("svar_ica: `data` must be a finite numeric matrix.",
                 class = "dynhr_input_error")
  K <- ncol(Y)
  if (K < 2L)
    .dynhr_abort("svar_ica: ICA identification needs at least 2 variables.",
                 class = "dynhr_input_error")
  if (!is.null(reference)) {
    reference <- as.matrix(reference)
    if (!all(dim(reference) == c(K, K)) || !all(is.finite(reference)))
      .dynhr_abort(sprintf("svar_ica: `reference` must be a finite %d x %d matrix.",
                           K, K), class = "dynhr_input_error")
  }

  vf <- estimate_var(Y, p = p, type = type)
  U  <- vf$residuals
  if (nrow(U) < 5L)
    .dynhr_abort("svar_ica: too few VAR residuals for ICA.",
                 class = "dynhr_input_error")
  D  <- t(chol(vf$Sigma))
  V  <- t(forwardsolve(D, t(U)))           # whitened: crossprod(V)/n = I

  if (method == "dcov") {
    rot <- .ica_dcov_rotation(V, n_starts = as.integer(n_starts), seed = seed)
    O   <- rot$O
    obj <- rot$objective
  } else {
    O   <- .ica_fastica(V)$O
    obj <- .ica_dcov_contrast(V %*% t(O))
  }

  B <- D %*% t(O)
  md1 <- NA_real_
  if (is.null(reference)) {
    B <- .ica_normalise_default(B)
  } else {
    al  <- .ica_align_md1(B, reference)
    B   <- al$B
    md1 <- al$D1
  }
  shock_names <- paste0("shock", seq_len(K))
  dimnames(B) <- list(vf$var_names, shock_names)

  E <- t(solve(B, t(U)))
  colnames(E) <- shock_names
  gs <- .ica_gaussianity(E, gauss_level)
  identified <- sum(gs$gaussian) <= 1L
  if (!identified && isTRUE(warn_gaussian))
    .dynhr_warn(sprintf(paste0(
      "svar_ica: %d of %d estimated structural shocks are not distinguishable ",
      "from Gaussian (Jarque-Bera p > %g): the ICA rotation among them is not ",
      "identified and those columns of the impact matrix are arbitrary."),
      sum(gs$gaussian), K, gauss_level),
      class = "dynhr_svar_ica_not_identified")

  irfs <- .var_irf_from_coef(vf$A, B, horizon, vf$var_names)

  structure(
    list(B = B, rotation = O, chol = D, irfs = irfs, shocks = E,
         objective = obj, gaussianity = gs, identified = identified,
         md1 = md1, var_fit = vf, method = method, p = vf$p,
         horizon = horizon),
    class = "svar_ica")
}


#' @export
print.svar_ica <- function(x, ...) {
  cat(sprintf("<svar_ica>  K = %d, p = %d, method = %s, identified = %s\n",
              ncol(x$B), x$p, x$method, x$identified))
  cat("  impact matrix B (rows = variables, columns = shocks):\n")
  print(round(x$B, 4))
  invisible(x)
}


## --------------------------------------------------------------------------
## Stack a list of horizon x K IRF matrices (one per shock) into the vector
## psi = (vec(Psi_0)', ..., vec(Psi_H)')' of Martinoli et al., where
## Psi_l[i, j] = response of variable i to shock j at horizon l.
## --------------------------------------------------------------------------
.svar_ica_psi <- function(irfs) {
  H <- nrow(irfs[[1L]])
  K <- ncol(irfs[[1L]])
  unlist(lapply(seq_len(H), function(h)
    as.numeric(vapply(irfs, function(M) M[h, ], numeric(K)))),
    use.names = FALSE)
}


## --------------------------------------------------------------------------
## Model shocks from FIXED standard draws (common random numbers): the same
## CSN law as simulate_model()'s diagonal skew branch,
##   e_k = sigma_k (delta_k |z1_k| + sqrt(1 - delta_k^2) z2_k) - sigma_k delta_k sqrt(2/pi),
## delta_k = alpha_k / sqrt(1 + alpha_k^2), applied COLUMN-wise (per shock).
## `z1` holds |N(0,1)| draws, `z2` N(0,1) draws.  With alpha = 0 this is the
## Gaussian N(0, sigma_k^2).
## --------------------------------------------------------------------------
.svar_ica_csn_shocks <- function(z1, z2, sigma, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  std <- sweep(z1, 2L, delta, `*`) + sweep(z2, 2L, sqrt(1 - delta^2), `*`)
  sweep(sweep(std, 2L, sigma, `*`), 2L, sigma * delta * sqrt(2 / pi), `-`)
}


## --------------------------------------------------------------------------
## Weight matrix for the II criterion.
## --------------------------------------------------------------------------
.svar_ica_weight <- function(weight, q) {
  if (is.null(weight)) return(diag(q))
  if (is.matrix(weight)) {
    if (!all(dim(weight) == c(q, q)) || !all(is.finite(weight)))
      .dynhr_abort(sprintf("match_irfs_svar_ica: `weight` matrix must be finite %d x %d.",
                           q, q), class = "dynhr_input_error")
    return((weight + t(weight)) / 2)
  }
  if (is.numeric(weight) && (length(weight) == 1L || length(weight) == q)) {
    if (any(!is.finite(weight)) || any(weight < 0))
      .dynhr_abort("match_irfs_svar_ica: `weight` must be finite and non-negative.",
                   class = "dynhr_input_error")
    return(diag(rep_len(weight, q), nrow = q))
  }
  .dynhr_abort(sprintf(paste0(
    "match_irfs_svar_ica: `weight` must be NULL, \"inverse_variance\", a scalar, ",
    "a length-%d vector or a %d x %d matrix."), q, q, q),
    class = "dynhr_input_error")
}


#' Indirect inference by matching non-Gaussian SVAR-ICA impulse responses
#'
#' The impulse-response-matching indirect-inference estimator of Martinoli, Di
#' Francesco, Moneta & Seri (2026).  The auxiliary model is a structural VAR
#' identified by non-Gaussianity (\code{\link{svar_ica}}).  It is fitted once to
#' the data, giving the structural IRF vector \eqn{\hat\psi_T} and impact
#' matrix \eqn{\hat\Psi_T}, and at each candidate \eqn{\theta} to \eqn{S}
#' simulated samples \eqn{y^s(\theta)} of the same length, giving
#' \eqn{\tilde\psi^s_T(\theta)}; then
#' \deqn{\hat\theta = \arg\min_\theta \bigl(\hat\psi_T - \bar\psi_T(\theta)\bigr)'
#'   W \bigl(\hat\psi_T - \bar\psi_T(\theta)\bigr), \quad
#'   \bar\psi_T(\theta) = \frac1S \sum_{s=1}^S \tilde\psi^s_T(\theta),}
#' with \eqn{\psi = (\mathrm{vec}(\Psi_0)', \ldots, \mathrm{vec}(\Psi_H)')'}.
#' Data and simulations go through the SAME auxiliary estimator, so its
#' small-sample and lag-truncation biases cancel (the Sims--Cogley--Nason
#' argument).  Each simulated impact matrix is aligned to \eqn{\hat\Psi_T} by
#' the MD1 minimum-distance index (the signed permutation minimising
#' \eqn{\lVert C \tilde\Psi^{s,-1} \hat\Psi_T - I\rVert_F}), which resolves the
#' ICA order/sign indeterminacy without any economic labelling; by
#' Martinoli et al.'s Proposition 1 the estimator is invariant to the
#' normalisation of \eqn{\hat\Psi_T}.
#'
#' The simulated shocks are the model's own: each shock \eqn{k} is drawn from
#' its skew-normal (CSN) law with the \code{stderr} and \code{skew} of the
#' model's \code{shocks} block evaluated at \eqn{\theta} (so an estimated
#' shock stderr or skewness is honoured), built from standard draws fixed once
#' under \code{seed} -- common random numbers, which make the criterion a
#' smooth function of \eqn{\theta}.  ICA needs at least \eqn{K-1}
#' non-Gaussian shocks, so declare \code{skew} for them (or pass non-Gaussian
#' \code{shock_draws}).  Correlated shocks (\code{corr}) are rejected: the
#' auxiliary model assumes independent structural shocks.
#'
#' @param model A \code{dynhr_mod} or a path to a \code{.mod} file.
#' @param params_init Named numeric vector of starting values for
#'   \code{free_params}; other names are fixed calibration overrides.
#' @param data \code{T x K} matrix or data frame of observables, columns named
#'   after model variables.
#' @param free_params Character vector of parameters to estimate.
#' @param obs Model variables used as VAR observables (default: the column
#'   names of \code{data}).  Their number must equal the number of model
#'   shocks, so that the SVAR is square.
#' @param p VAR lag order.  Default \code{2}.
#' @param horizon Number of IRF periods matched, impact included (\eqn{H+1}).
#'   Default \code{5}.
#' @param n_sim Number of simulated samples \eqn{S}.  Default \code{20}.
#' @param n_periods Length of each simulated sample.  Default \code{nrow(data)}.
#' @param burn_in Periods discarded before each simulated sample.  Default
#'   \code{200}.
#' @param weight Weight matrix \eqn{W}: \code{NULL} (identity, as in
#'   Martinoli et al.'s Monte Carlo), a scalar, a length-\eqn{q} vector (a
#'   diagonal), a full \eqn{q \times q} matrix with \eqn{q = K^2 \cdot}
#'   \code{horizon}, or \code{"inverse_variance"}: the diagonal inverse of
#'   \eqn{(1 + 1/S)} times the across-simulation variance of
#'   \eqn{\tilde\psi^s} at \code{params_init} -- the diagonal of the
#'   II-optimal \eqn{\zeta_{II}^{-1}} when data and simulations are independent.
#' @param ica_method ICA estimator for data and simulations alike:
#'   \code{"dcov"} (default, as in Martinoli et al.) or \code{"fastica"}.
#' @param type VAR deterministic terms (see \code{\link{estimate_var}}).
#' @param order Perturbation order of the simulations (\code{1} or \code{2}).
#' @param compiled Optional pre-compiled model.
#' @param method Optimiser: any \code{\link[stats]{optim}} method.  For one
#'   free parameter \code{"Brent"} with finite \code{lower}/\code{upper} is
#'   the robust choice.
#' @param lower,upper Bounds for \code{"L-BFGS-B"} / \code{"Brent"}.
#' @param penalty Objective returned where the model cannot be solved or is
#'   not stationary.  Default \code{1e10}.
#' @param control Passed to \code{optim}.
#' @param seed Seed for the common random numbers (caller's RNG restored).
#' @param shock_draws Optional list of \code{n_sim} matrices, each
#'   \code{(n_periods + burn_in) x n_exo}, of STANDARDISED (mean 0, variance 1)
#'   shock draws replacing the CSN law, e.g. scaled Student-t draws as in
#'   Martinoli et al.'s Monte Carlo; column \eqn{k} is scaled by shock
#'   \eqn{k}'s stderr at \eqn{\theta}.
#' @param verbose Print the criterion at each evaluation.
#' @param ... Passed to \code{optim}.
#'
#' @return A list of class \code{c("irf_match_svar_ica", "irf_match_result")}:
#'   \describe{
#'     \item{par, params, objective, convergence, counts, message, optim}{As in
#'       \code{\link{match_irfs}}.}
#'     \item{target_psi, fitted_psi}{\eqn{\hat\psi_T} and
#'       \eqn{\bar\psi_T(\hat\theta)}.}
#'     \item{target_irf, fitted_irf}{The same as lists of \code{horizon x K}
#'       matrices keyed by shock (the data SVAR's normalisation).}
#'     \item{data_fit}{The \code{\link{svar_ica}} fit to the data.}
#'     \item{weight}{The \eqn{W} used.}
#'     \item{n_sim, n_periods, obs}{Bookkeeping.}
#'   }
#'
#' @references
#' Martinoli, M., Di Francesco, D., Moneta, A. and Seri, R. (2026). Estimation
#' of DSGE models by non-Gaussian vector autoregressions. LEM Working Paper
#' 2026/19, Sant'Anna School of Advanced Studies.
#'
#' @examples
#' \dontrun{
#' m <- parse_mod("
#'   var x pie; varexo e_d e_s; parameters rho kappa;
#'   rho = 0.7; kappa = 0.3;
#'   model(linear);
#'     x = rho*x(-1) + e_d - 0.5*e_s;
#'     pie = 0.99*pie(+1) + kappa*x + e_s;
#'   end;
#'   shocks; var e_d; stderr 1; var e_s; stderr 0.5;
#'     skew e_d = 8; skew e_s = -6; end;", verbose = FALSE)
#' ## ... data <- simulated or observed (x, pie) ...
#' fit <- match_irfs_svar_ica(m, c(kappa = 0.5), data, "kappa", p = 1,
#'                            horizon = 4, method = "Brent",
#'                            lower = 0.05, upper = 0.8)
#' fit$par
#' }
#'
#' @seealso \code{\link{svar_ica}}, \code{\link{match_irfs}},
#'   \code{\link{simulate_model}}
#' @export
match_irfs_svar_ica <- function(model, params_init, data, free_params,
                                obs = colnames(data), p = 2L, horizon = 5L,
                                n_sim = 20L, n_periods = NULL, burn_in = 200L,
                                weight = NULL,
                                ica_method = c("dcov", "fastica"),
                                type = c("const", "none", "trend"),
                                order = 1L, compiled = NULL,
                                method = "Nelder-Mead", lower = -Inf,
                                upper = Inf, penalty = 1e10,
                                control = list(), seed = 1L,
                                shock_draws = NULL, verbose = FALSE, ...) {
  ica_method <- match.arg(ica_method)
  type <- match.arg(type)

  ## ---- Model ----------------------------------------------------------
  if (is.character(model) && length(model) == 1L && file.exists(model))
    model <- parse_mod(model, verbose = FALSE)
  if (!inherits(model, "dynhr_mod"))
    .dynhr_abort("match_irfs_svar_ica: `model` must be a dynhr_mod or a path ",
                 "to a .mod file.", class = "dynhr_input_error")
  order <- as.integer(order)
  if (!order %in% c(1L, 2L))
    .dynhr_abort("match_irfs_svar_ica: `order` must be 1 or 2.",
                 class = "dynhr_input_error")

  ## ---- Free parameters / start ----------------------------------------
  if (!is.character(free_params) || length(free_params) == 0L ||
      anyDuplicated(free_params))
    .dynhr_abort("match_irfs_svar_ica: `free_params` must be a non-empty ",
                 "character vector without duplicates.",
                 class = "dynhr_input_error")
  unknown <- setdiff(free_params, model$param_names)
  if (length(unknown))
    .dynhr_abort(sprintf("match_irfs_svar_ica: free_params not in model$param_names: %s",
                         paste(unknown, collapse = ", ")),
                 class = "dynhr_input_error")
  if (is.null(names(params_init)) || !all(free_params %in% names(params_init)))
    .dynhr_abort("match_irfs_svar_ica: `params_init` must be a named vector ",
                 "with a value for every free parameter.",
                 class = "dynhr_input_error")
  theta0 <- setNames(as.numeric(params_init[free_params]), free_params)
  base_params <- model$param_values
  for (nm in setdiff(names(params_init), free_params))
    if (nm %in% names(base_params)) base_params[[nm]] <- params_init[[nm]]

  ## ---- Observables / data ---------------------------------------------
  exo <- model$varexo_names
  n_exo <- length(exo)
  if (is.null(obs))
    .dynhr_abort("match_irfs_svar_ica: supply `obs` or column names on `data`.",
                 class = "dynhr_input_error")
  if (!all(obs %in% model$var_names))
    .dynhr_abort(sprintf("match_irfs_svar_ica: obs not model variables: %s",
                         paste(setdiff(obs, model$var_names), collapse = ", ")),
                 class = "dynhr_input_error")
  K <- length(obs)
  if (K != n_exo)
    .dynhr_abort(sprintf(paste0(
      "match_irfs_svar_ica: %d observables but %d shocks; the ICA-identified ",
      "SVAR needs as many observables as structural shocks."), K, n_exo),
      class = "dynhr_input_error")
  Y <- as.matrix(data)
  if (!is.null(colnames(Y))) {
    if (!all(obs %in% colnames(Y)))
      .dynhr_abort("match_irfs_svar_ica: `data` lacks columns for all `obs`.",
                   class = "dynhr_input_error")
    Y <- Y[, obs, drop = FALSE]
  } else if (ncol(Y) == K) {
    colnames(Y) <- obs
  } else {
    .dynhr_abort("match_irfs_svar_ica: `data` needs one column per observable.",
                 class = "dynhr_input_error")
  }
  n_periods <- as.integer(n_periods %||% nrow(Y))
  burn_in <- as.integer(burn_in)
  n_sim <- as.integer(n_sim)
  if (n_sim < 1L)
    .dynhr_abort("match_irfs_svar_ica: `n_sim` must be >= 1.",
                 class = "dynhr_input_error")

  ## Independent shocks only.
  Sig0 <- .get_shock_cov(model, exo, base_params)
  if (any(abs(Sig0[upper.tri(Sig0)]) > 0))
    .dynhr_abort("match_irfs_svar_ica: the shocks block declares correlated ",
                 "shocks; the ICA auxiliary model assumes independent ",
                 "structural shocks.", class = "dynhr_input_error")

  ## ---- Auxiliary model on the data ------------------------------------
  data_fit <- svar_ica(Y, p = p, horizon = horizon, type = type,
                       method = ica_method, seed = seed)
  psi_hat <- .svar_ica_psi(data_fit$irfs)
  ref <- data_fit$B
  q <- length(psi_hat)

  ## ---- Common random numbers ------------------------------------------
  n_tot <- n_periods + burn_in
  if (is.null(shock_draws)) {
    crn <- .with_local_seed(seed, lapply(seq_len(n_sim), function(s)
      list(z1 = abs(matrix(stats::rnorm(n_tot * n_exo), n_tot, n_exo)),
           z2 = matrix(stats::rnorm(n_tot * n_exo), n_tot, n_exo))))
  } else {
    ok <- is.list(shock_draws) && length(shock_draws) == n_sim &&
      all(vapply(shock_draws, function(z) is.matrix(z) &&
                   all(dim(z) == c(n_tot, n_exo)) && all(is.finite(z)),
                 logical(1)))
    if (!ok)
      .dynhr_abort(sprintf(paste0(
        "match_irfs_svar_ica: `shock_draws` must be a list of %d finite ",
        "%d x %d matrices (n_periods + burn_in rows, one column per shock)."),
        n_sim, n_tot, n_exo), class = "dynhr_input_error")
  }

  if (is.null(compiled))
    compiled <- compile_model(model, verbose = FALSE, max_order = order)

  ## ---- Simulated auxiliary IRFs at theta ------------------------------
  ## Returns an n_sim x q matrix of aligned psi's, or NULL when theta is not
  ## admissible (no steady state, BK failure, non-stationary).
  .sim_psi <- function(theta) {
    params <- base_params
    params[free_params] <- theta
    sol <- .mom_solve_at(model, compiled, params, order)
    if (is.null(sol)) return(NULL)
    pr <- sol$params
    sig <- .get_shock_stderr(model, exo, pr)[exo]
    alp <- .get_shock_skewness(model, exo, pr)[exo]
    out <- matrix(NA_real_, n_sim, q)
    for (s in seq_len(n_sim)) {
      shocks <- if (is.null(shock_draws)) {
        .svar_ica_csn_shocks(crn[[s]]$z1, crn[[s]]$z2, sig, alp)
      } else {
        sweep(shock_draws[[s]], 2L, sig, `*`)
      }
      sim <- .simulate_dr_any_order(sol$dr, n_periods = n_periods,
                                    model = model, burn_in = burn_in,
                                    shocks = shocks)
      ys <- sim[, obs, drop = FALSE]
      if (!all(is.finite(ys))) return(NULL)
      fs <- svar_ica(ys, p = p, horizon = horizon, type = type,
                     method = ica_method, reference = ref,
                     warn_gaussian = FALSE, seed = seed)
      out[s, ] <- .svar_ica_psi(fs$irfs)
    }
    out
  }

  ## ---- Weight ---------------------------------------------------------
  if (identical(weight, "inverse_variance")) {
    ps0 <- .sim_psi(theta0)
    if (is.null(ps0) || n_sim < 2L)
      .dynhr_abort("match_irfs_svar_ica: weight = \"inverse_variance\" needs ",
                   "n_sim >= 2 and a solvable model at params_init.",
                   class = "dynhr_input_error")
    v <- apply(ps0, 2L, stats::var) * (1 + 1 / n_sim)
    v[v <= 0] <- min(v[v > 0])
    W <- diag(1 / v, nrow = q)
  } else {
    W <- .svar_ica_weight(weight, q)
  }

  ## ---- Criterion ------------------------------------------------------
  obj <- function(theta) {
    ps <- .sim_psi(theta)
    if (is.null(ps)) {
      if (verbose) .dynhr_cat(sprintf("  obj = %.6g  [penalty]\n", penalty))
      return(penalty)
    }
    d <- psi_hat - colMeans(ps)
    val <- as.numeric(crossprod(d, W %*% d))
    if (!is.finite(val)) val <- penalty
    if (verbose) .dynhr_cat(sprintf("  obj = %.6g\n", val))
    val
  }

  optim_args <- list(par = theta0, fn = obj, method = method,
                     control = control, ...)
  if (method %in% c("L-BFGS-B", "Brent")) {
    optim_args$lower <- lower
    optim_args$upper <- upper
  }
  opt <- do.call(stats::optim, optim_args)

  par_hat <- setNames(opt$par, free_params)
  params_hat <- base_params
  params_hat[free_params] <- par_hat
  ps_hat <- .sim_psi(par_hat)
  fitted_psi <- if (is.null(ps_hat)) NULL else colMeans(ps_hat)

  unstack <- function(psi) {
    if (is.null(psi)) return(NULL)
    arr <- array(psi, c(K, K, horizon))          # [variable, shock, horizon]
    out <- lapply(seq_len(K), function(j) {
      M <- t(arr[, j, , drop = TRUE])
      M <- matrix(M, horizon, K, dimnames = list(NULL, obs))
      M
    })
    names(out) <- paste0("shock", seq_len(K))
    out
  }

  structure(
    list(par = par_hat, params = params_hat, objective = opt$value,
         convergence = opt$convergence, counts = opt$counts,
         message = opt$message, optim = opt,
         target_psi = psi_hat, fitted_psi = fitted_psi,
         target_irf = unstack(psi_hat), fitted_irf = unstack(fitted_psi),
         data_fit = data_fit, weight = W, shock = "all (SVAR-ICA)",
         free_params = free_params, horizons = seq_len(horizon),
         n_sim = n_sim, n_periods = n_periods, obs = obs),
    class = c("irf_match_svar_ica", "irf_match_result"))
}
