## R/mom-se-bounds.R
## ---------------------------------------------------------------------------
## Standard-error BOUNDS for a moment-matching estimator whose moments come
## from different samples, so that only their marginal variances (and perhaps
## a few cross-moment correlations) are known.
##
## Estimand. A minimum-distance estimator with Jacobian G (p x k) and weight
## W has the first-order expansion
##     lambda' (theta_hat - theta) = x' (mu_hat - mu),
##     x = W G (G'WG)^{-1} lambda,
## so its variance is x' D Omega D x, D = diag(sd of mu_hat_j), Omega the
## (partly unknown) correlation matrix of the moment estimators.
##
## Diagonal-only information (Cha, Sasaki & Tan 2026, Theorem 3.1; the upper
## half is Cocci & Plagborg-Moller 2025, Lemma 1), with z_j = |x_j| s_j:
##     upper = sum_j z_j                     (all errors aligned),
##     lower = max(max_m z_m - sum_{j!=m} z_j, 0)   (maximal cancellation).
## The lower bound is zero exactly when the lengths z_j close into a polygon
## (the largest is no longer than the sum of the rest).
##
## Partial information (some Omega_jk known): the sharp bounds are the two
## semidefinite programs min / max  z' Omega z  s.t. Omega >= 0, diag = 1,
## Omega_jk = rho_jk on the known pairs (Cha, Sasaki & Tan 2026, Prop. 5.2;
## Cocci & Plagborg-Moller 2025, Sec. 3.3). They are solved here without an
## SDP dependency by a log-barrier Newton method on the DUAL
##     max b'y  s.t.  C - sum_k y_k A_k  >= 0,
## whose feasible points certify: the reported lower variance never exceeds
## the true minimum and the reported upper variance never falls below the
## true maximum (weak duality), to within the barrier gap `tol`.
##
## Known correlations of exactly +/-1 are merged before solving: two moment
## estimators with correlation +/-1 are (asymptotically) the same random
## variable up to sign and scale, so they collapse into one coordinate. This
## is exact, and it keeps the SDP's primal feasible set with a non-empty
## interior, which the barrier needs.
## ---------------------------------------------------------------------------


#' Standard-error bounds for moments estimated on different samples
#'
#' Bounds the standard error of a moment-matching (minimum-distance)
#' estimator when the empirical moments come from different samples, so that
#' only their marginal variances are known and the cross-moment correlations
#' are not. The worst case (upper bound) is Cocci & Plagborg-Moller's
#' (2025); the best case (lower bound) and the semidefinite program for
#' partially known correlations are Cha, Sasaki & Tan's (2026).
#'
#' For each target \eqn{\lambda' \theta} the estimator is asymptotically the
#' linear combination \eqn{x'\hat\mu} of the moments, with
#' \eqn{x = W G (G'WG)^{-1}\lambda}. With \eqn{s_j^2} the variance of
#' \eqn{\hat\mu_j} and \eqn{z_j = |x_j| s_j}:
#' \itemize{
#'   \item \strong{No correlation known.} The sharp bounds are closed form:
#'     \deqn{\underline{se} = \max\{\max_m z_m - \textstyle\sum_{j\ne m} z_j,
#'       \; 0\}, \qquad \overline{se} = \textstyle\sum_j z_j.}
#'     The upper bound is attained when all estimation errors are aligned,
#'     the lower bound when they cancel as far as possible. The lower bound
#'     is zero exactly when the largest \eqn{z_j} is no longer than the sum
#'     of the others: the lengths \eqn{z_j} then close into a polygon
#'     (column \code{polygon_closes}).
#'   \item \strong{Some correlations known} (\code{cor_known}). The bounds are
#'     the minimum and maximum of \eqn{z'\Omega z} over correlation matrices
#'     \eqn{\Omega \succeq 0} that match the known entries -- a pair of
#'     semidefinite programs. They are solved by a log-barrier Newton method
#'     on the dual problem, so no SDP package is needed. The dual values are
#'     reported, so the interval is conservative: the lower variance can be
#'     low, and the upper variance high, by at most \code{tol} times the
#'     squared diagonal-only upper bound.
#'   \item \strong{All correlations known.} Both bounds equal the usual
#'     sandwich standard error.
#' }
#' Known correlations of exactly \eqn{\pm 1} are merged before solving,
#' because two such moment estimators are the same random variable up to
#' sign and scale.
#'
#' The standard errors are on the scale of \code{var_moments}. Pass the
#' variances of the moment ESTIMATORS (for example squared reported standard
#' errors), not the variances of the underlying data.
#'
#' @param G Numeric \eqn{p \times k} Jacobian of the model moments with
#'   respect to the parameters (full column rank). Column names become the
#'   target names.
#' @param W Symmetric \eqn{p \times p} weight matrix, or \code{NULL} (the
#'   default) for the identity.
#' @param var_moments Length-\eqn{p} vector of the (non-negative) variances
#'   of the moment estimators.
#' @param cor_known Optional \eqn{p \times p} matrix of the known
#'   correlations between moment estimators. Use \code{NA} for unknown
#'   entries. The diagonal must be 1 or \code{NA}, and the known pattern must
#'   be symmetric. \code{NULL} (the default) means no correlation is known.
#' @param lambda Optional \eqn{k \times m} matrix (or length-\eqn{k} vector)
#'   whose columns are the linear combinations of the parameters to bound.
#'   The default, the identity, bounds each parameter.
#' @param tol Duality-gap tolerance of the semidefinite programs, relative to
#'   the squared diagonal-only upper bound (default \code{1e-10}).
#' @param max_iter Maximum total Newton steps per semidefinite program
#'   (default 500).
#'
#' @return An object of class \code{"dynhr_se_bounds"}: a list with
#'   \describe{
#'     \item{\code{bounds}}{Data frame with one row per target and the
#'       columns \code{target}, \code{se_lower}, \code{se_upper},
#'       \code{se_uncorrelated} (the standard error if every moment were
#'       uncorrelated with every other), \code{polygon_closes} (the
#'       diagonal-only zero-lower-bound condition), and \code{solver}
#'       (\code{"closed_form"}, \code{"sdp"} or \code{"known"}).}
#'     \item{\code{loadings}}{The \eqn{p \times m} matrix of loadings
#'       \eqn{x}.}
#'     \item{\code{cor_lower}, \code{cor_upper}}{Lists (one per target) of
#'       correlation matrices that attain the lower and upper bound. For the
#'       semidefinite-program case these are the solver's primal iterates,
#'       accurate to the tolerance.}
#'     \item{\code{var_moments}, \code{cor_known}}{The inputs.}
#'   }
#'
#' @references
#'   Cha, J., Sasaki, Y., & Tan, N. M. P. (2026). Bounds for standard errors
#'     in combined data. arXiv:2606.24867.
#'
#'   Cocci, M. D., & Plagborg-Moller, M. (2025). Standard errors for
#'     calibrated parameters. \emph{Review of Economic Studies}, 92(5),
#'     2952-2978. \doi{10.1093/restud/rdae099}
#'
#' @seealso \code{\link{method_of_moments}} (argument \code{se_bounds}).
#' @examples
#' ## Repeated measurements: two estimates of one parameter, equally weighted.
#' b <- mom_se_bounds(G = matrix(1, 2, 1, dimnames = list(NULL, "theta")),
#'                    var_moments = c(0.3, 0.1)^2)
#' b$bounds   # lower |0.3 - 0.1| / 2 = 0.1, upper (0.3 + 0.1) / 2 = 0.2
#'
#' ## Knowing the correlation pins the standard error down.
#' mom_se_bounds(G = matrix(1, 2, 1), var_moments = c(0.3, 0.1)^2,
#'               cor_known = matrix(c(1, 0.5, 0.5, 1), 2))$bounds
#' @export
mom_se_bounds <- function(G, W = NULL, var_moments, cor_known = NULL,
                          lambda = NULL, tol = 1e-10, max_iter = 500L) {
  if (is.data.frame(G)) G <- as.matrix(G)
  if (is.null(dim(G))) G <- matrix(G, ncol = 1L)
  if (!is.numeric(G) || !all(is.finite(G)))
    .dynhr_abort("mom_se_bounds(): `G` must be a finite numeric matrix.",
                 class = "dynhr_input_error")
  p <- nrow(G)
  k <- ncol(G)
  if (p < k)
    .dynhr_abort("mom_se_bounds(): `G` has ", p, " moments for ", k,
                 " parameters; need p >= k.", class = "dynhr_input_error")

  if (is.null(W)) W <- diag(p)
  W <- as.matrix(W)
  if (!is.numeric(W) || !all(dim(W) == c(p, p)) || !all(is.finite(W)))
    .dynhr_abort("mom_se_bounds(): `W` must be a finite ", p, " x ", p,
                 " matrix.", class = "dynhr_input_error")
  if (max(abs(W - t(W))) > 1e-10 * max(1, max(abs(W))))
    .dynhr_abort("mom_se_bounds(): `W` must be symmetric.",
                 class = "dynhr_input_error")

  var_moments <- as.numeric(var_moments)
  if (length(var_moments) != p || !all(is.finite(var_moments)) ||
      any(var_moments < 0))
    .dynhr_abort("mom_se_bounds(): `var_moments` must be ", p,
                 " finite non-negative variances (one per row of `G`).",
                 class = "dynhr_input_error")
  s <- sqrt(var_moments)

  pnames <- colnames(G)
  if (is.null(pnames)) pnames <- paste0("theta", seq_len(k))
  if (is.null(lambda)) {
    lambda <- diag(k)
    colnames(lambda) <- pnames
  }
  if (is.null(dim(lambda))) lambda <- matrix(lambda, ncol = 1L)
  lambda <- as.matrix(lambda)
  if (!is.numeric(lambda) || nrow(lambda) != k || !all(is.finite(lambda)))
    .dynhr_abort("mom_se_bounds(): `lambda` must be a finite matrix with ",
                 k, " rows (one per column of `G`).",
                 class = "dynhr_input_error")
  tnames <- colnames(lambda)
  if (is.null(tnames)) tnames <- paste0("comb", seq_len(ncol(lambda)))

  cor_known <- .sebnd_check_cor(cor_known, p)

  ## Loadings x = W G (G'WG)^{-1} lambda.
  A <- crossprod(G, W %*% G)
  qa <- qr(A)
  if (qa$rank < k)
    .dynhr_abort("mom_se_bounds(): G'WG is singular (rank ", qa$rank,
                 " < ", k, "); the parameters are not identified by the ",
                 "weighted moments.", class = "dynhr_input_error")
  X <- W %*% G %*% qr.solve(qa, lambda)
  dimnames(X) <- list(rownames(G), tnames)

  m <- ncol(X)
  rows <- vector("list", m)
  cor_lo <- vector("list", m)
  cor_up <- vector("list", m)
  for (i in seq_len(m)) {
    zs <- s * X[, i]                   # signed z: variance = zs' Omega zs
    z  <- abs(zs)
    res <- .sebnd_one(zs, cor_known, tol = tol, max_iter = as.integer(max_iter))
    zmax <- if (p) max(z) else 0
    rows[[i]] <- data.frame(
      target          = tnames[i],
      se_lower        = res$lower,
      se_upper        = res$upper,
      se_uncorrelated = sqrt(sum(z^2)),
      polygon_closes  = zmax <= sum(z) - zmax,
      solver          = res$solver,
      stringsAsFactors = FALSE)
    cor_lo[[i]] <- res$cor_lower
    cor_up[[i]] <- res$cor_upper
  }
  names(cor_lo) <- names(cor_up) <- tnames
  bounds <- do.call(rbind, rows)
  rownames(bounds) <- NULL

  structure(list(bounds = bounds, loadings = X,
                 cor_lower = cor_lo, cor_upper = cor_up,
                 var_moments = var_moments, cor_known = cor_known),
            class = "dynhr_se_bounds")
}


#' Print standard-error bounds
#'
#' @param x A \code{"dynhr_se_bounds"} object from
#'   \code{\link{mom_se_bounds}}.
#' @param digits Significant digits.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.dynhr_se_bounds <- function(x, digits = 4L, ...) {
  nk <- if (is.null(x$cor_known)) 0L else
    sum(!is.na(x$cor_known[upper.tri(x$cor_known)]))
  cat("=== dynhr standard-error bounds (combined data) ===\n")
  cat(sprintf("  Moments: %d | known cross-moment correlations: %d\n",
              length(x$var_moments), nk))
  tb <- x$bounds
  num <- c("se_lower", "se_upper", "se_uncorrelated")
  tb[num] <- lapply(tb[num], function(v) signif(v, digits))
  print(tb, row.names = FALSE)
  invisible(x)
}


## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

#' Validate `cor_known`: NULL, or a p x p matrix with NA for unknown entries,
#' diagonal 1/NA (set to 1), symmetric known pattern, |values| <= 1.
#' @noRd
.sebnd_check_cor <- function(cor_known, p) {
  if (is.null(cor_known)) return(NULL)
  cor_known <- as.matrix(cor_known)
  if (!all(dim(cor_known) == c(p, p)))
    .dynhr_abort("mom_se_bounds(): `cor_known` must be ", p, " x ", p, ".",
                 class = "dynhr_input_error")
  storage.mode(cor_known) <- "double"
  d <- diag(cor_known)
  if (any(!is.na(d) & abs(d - 1) > 1e-12))
    .dynhr_abort("mom_se_bounds(): the diagonal of `cor_known` must be 1 ",
                 "(or NA).", class = "dynhr_input_error")
  diag(cor_known) <- 1
  known <- !is.na(cor_known)
  if (!identical(known, t(known)))
    .dynhr_abort("mom_se_bounds(): the known pattern of `cor_known` must be ",
                 "symmetric.", class = "dynhr_input_error")
  if (any(!is.finite(cor_known[known])))
    .dynhr_abort("mom_se_bounds(): `cor_known` entries must be finite or NA.",
                 class = "dynhr_input_error")
  if (any(abs(cor_known[known] - t(cor_known)[known]) > 1e-12))
    .dynhr_abort("mom_se_bounds(): `cor_known` must be symmetric.",
                 class = "dynhr_input_error")
  if (any(abs(cor_known[known]) > 1 + 1e-12))
    .dynhr_abort("mom_se_bounds(): `cor_known` entries must lie in [-1, 1].",
                 class = "dynhr_input_error")
  cor_known[known] <- pmin(pmax(cor_known[known], -1), 1)
  cor_known
}


#' Bounds for one signed-z vector: returns list(lower, upper, solver,
#' cor_lower, cor_upper) with lower/upper on the SE (not variance) scale.
#' @noRd
.sebnd_one <- function(zs, cor_known, tol, max_iter) {
  p <- length(zs)
  ## Fully known: the ordinary standard error.
  if (!is.null(cor_known) && !anyNA(cor_known)) {
    ev <- eigen(cor_known, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8 * max(1, max(ev)))
      .dynhr_abort("mom_se_bounds(): `cor_known` is fully specified but not ",
                   "positive semidefinite (smallest eigenvalue ",
                   signif(min(ev), 3), ").",
                   class = "dynhr_se_bounds_infeasible")
    se <- sqrt(max(as.numeric(crossprod(zs, cor_known %*% zs)), 0))
    return(list(lower = se, upper = se, solver = "known",
                cor_lower = cor_known, cor_upper = cor_known))
  }
  if (is.null(cor_known)) {
    kn0 <- matrix(NA_real_, p, p)
    diag(kn0) <- 1
    red <- list(w = zs, P = diag(p), known = kn0)
  } else {
    red <- .sebnd_merge(zs, cor_known)
  }
  q <- length(red$w)
  kn_off <- !is.na(red$known) & upper.tri(red$known)
  if (!any(kn_off)) {
    cf <- .sebnd_closed_form(red$w)
    out <- list(lower = cf$lower, upper = cf$upper, solver = "closed_form",
                R_lo = cf$R_lower, R_up = cf$R_upper)
  } else if (!anyNA(red$known)) {
    ev <- eigen(red$known, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8 * max(1, max(ev)))
      .dynhr_abort("mom_se_bounds(): the known correlations admit no ",
                   "positive semidefinite completion.",
                   class = "dynhr_se_bounds_infeasible")
    se <- sqrt(max(as.numeric(crossprod(red$w, red$known %*% red$w)), 0))
    out <- list(lower = se, upper = se, solver = "known",
                R_lo = red$known, R_up = red$known)
  } else {
    scale <- sum(abs(red$w))
    if (scale == 0) {
      ## Nothing loads on a random moment. Any completion works; find one.
      sd0 <- .sebnd_sdp(matrix(0, q, q), red$known, tol, max_iter, sense = 1)
      out <- list(lower = 0, upper = 0, solver = "sdp",
                  R_lo = sd0$Omega, R_up = sd0$Omega)
    } else {
      wn <- red$w / scale
      C  <- tcrossprod(wn)
      lo <- .sebnd_sdp(C, red$known, tol, max_iter, sense = 1)
      up <- .sebnd_sdp(C, red$known, tol, max_iter, sense = -1)
      out <- list(lower = scale * sqrt(max(lo$value, 0)),
                  upper = scale * sqrt(max(up$value, 0)),
                  solver = "sdp", R_lo = lo$Omega, R_up = up$Omega)
    }
  }
  expand <- function(R) {
    Rf <- red$P %*% R %*% t(red$P)
    diag(Rf) <- 1
    (Rf + t(Rf)) * 0.5
  }
  list(lower = out$lower, upper = out$upper, solver = out$solver,
       cor_lower = expand(out$R_lo), cor_upper = expand(out$R_up))
}


#' Diagonal-only closed form (Cha-Sasaki-Tan Thm 3.1 / CPM Lemma 1) with
#' attaining correlation matrices for the signed vector w.
#' @noRd
.sebnd_closed_form <- function(w) {
  q <- length(w)
  sg <- ifelse(w < 0, -1, 1)
  z  <- abs(w)
  tot <- sum(z)
  jm <- which.max(z)
  zmax <- z[jm]
  lower <- max(zmax - (tot - zmax), 0)
  ## Upper: all sign-adjusted errors aligned, T = 1 1'.
  R_up <- tcrossprod(sg)
  ## Lower: unit vectors u_j with sum_j z_j u_j as short as possible.
  if (zmax >= tot - zmax) {
    ## No polygon: largest against all the rest (rank 1).
    t1 <- rep(-1, q)
    t1[jm] <- 1
    U <- matrix(t1, q, 1L)
  } else {
    ## Polygon closes: pack the lengths into three sides (largest-first into
    ## the currently shortest side; each side then is <= tot / 2, so the
    ## triangle inequalities hold) and lay the sides out as a triangle.
    ord  <- order(z, decreasing = TRUE)
    side <- integer(q)
    load <- numeric(3)
    for (j in ord) {
      b <- which.min(load)
      side[j] <- b
      load[b] <- load[b] + z[j]
    }
    ## Side 1 is the longest: it holds the largest element and at least
    ## as much as any other side after the first three assignments ...
    ## but to be safe, index sides by decreasing load.
    so <- order(load, decreasing = TRUE)
    a <- load[so[1]]; bb <- load[so[2]]; cc <- load[so[3]]
    cx <- (a^2 + cc^2 - bb^2) / (2 * a)
    cy <- sqrt(max(cc^2 - cx^2, 0))
    dir <- matrix(0, 3, 2)
    dir[so[1], ] <- c(1, 0)                              # A -> B
    vb <- c(cx - a, cy)                                  # B -> C
    vc <- c(-cx, -cy)                                    # C -> A
    dir[so[2], ] <- if (bb > 0) vb / sqrt(sum(vb^2)) else c(1, 0)
    dir[so[3], ] <- if (cc > 0) vc / sqrt(sum(vc^2)) else c(1, 0)
    U <- dir[side, , drop = FALSE]
  }
  Tm <- tcrossprod(U)
  diag(Tm) <- 1
  R_lo <- Tm * tcrossprod(sg)
  list(lower = lower, upper = tot, R_lower = R_lo, R_upper = R_up)
}


#' Merge moments linked by known correlations of exactly +/-1.
#'
#' Returns list(w, P, known): the reduced signed-z vector, the p x q
#' expansion matrix (P[j, g] = sign of j within its group g) and the q x q
#' reduced known-correlation matrix (NA = unknown). Aborts when the +/-1
#' links or the implied known entries contradict each other.
#' @noRd
.sebnd_merge <- function(zs, cor_known) {
  p <- length(zs)
  one <- !is.na(cor_known) & abs(abs(cor_known) - 1) <= 1e-12
  diag(one) <- FALSE
  grp <- rep(NA_integer_, p)
  sgn <- rep(0, p)
  ng <- 0L
  for (st in seq_len(p)) {
    if (!is.na(grp[st])) next
    ng <- ng + 1L
    grp[st] <- ng
    sgn[st] <- 1
    queue <- st
    while (length(queue)) {
      i <- queue[1L]
      queue <- queue[-1L]
      for (j in which(one[i, ])) {
        sj <- sgn[i] * sign(cor_known[i, j])
        if (is.na(grp[j])) {
          grp[j] <- ng
          sgn[j] <- sj
          queue <- c(queue, j)
        } else if (sgn[j] != sj) {
          .dynhr_abort("mom_se_bounds(): the +/-1 entries of `cor_known` ",
                       "contradict each other (moments ", i, " and ", j, ").",
                       class = "dynhr_se_bounds_infeasible")
        }
      }
    }
  }
  P <- matrix(0, p, ng)
  P[cbind(seq_len(p), grp)] <- sgn
  w <- as.numeric(crossprod(P, zs))
  known <- matrix(NA_real_, ng, ng)
  diag(known) <- 1
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    cij <- cor_known[i, j]
    if (is.na(cij)) next
    gi <- grp[i]; gj <- grp[j]
    v <- sgn[i] * sgn[j] * cij
    if (gi == gj) {
      if (abs(v - 1) > 1e-9)
        .dynhr_abort("mom_se_bounds(): known correlation of moments ", i,
                     " and ", j, " contradicts the +/-1 links in ",
                     "`cor_known`.", class = "dynhr_se_bounds_infeasible")
      next
    }
    old <- known[gi, gj]
    if (!is.na(old) && abs(old - v) > 1e-9)
      .dynhr_abort("mom_se_bounds(): known correlations in `cor_known` ",
                   "contradict the +/-1 links (moments ", i, " and ", j,
                   ").", class = "dynhr_se_bounds_infeasible")
    known[gi, gj] <- known[gj, gi] <- v
  }
  list(w = w, P = P, known = known)
}


#' Log-barrier dual solver for  min (sense = 1) / max (sense = -1)
#'   <C, Omega>  s.t.  Omega >= 0, Omega_ab = known_ab on known entries.
#'
#' Works on the dual  max b'y s.t. S(y) = sense*C - sum_k y_k A_k >= 0 with
#' A_k = e_a e_b' + e_b e_a' (off-diagonal) or e_a e_a' (diagonal), which is
#' always strictly feasible. Returns the certified value of the ORIGINAL
#' problem (dual bound, conservative in the direction of `sense`) and the
#' primal iterate Omega = mu S^{-1} rescaled to a unit diagonal.
#' @noRd
.sebnd_sdp <- function(C, known, tol, max_iter, sense) {
  q <- nrow(C)
  Cs <- sense * C
  idx <- which(!is.na(known) & upper.tri(known, diag = TRUE), arr.ind = TRUE)
  a <- idx[, 1L]
  b <- idx[, 2L]
  isd <- a == b
  wt <- ifelse(isd, 0.5, 1)
  bvec <- ifelse(isd, 1, 2 * known[cbind(a, b)])
  nk <- length(a)

  Amap <- function(y) {
    Y <- matrix(0, q, q)
    Y[cbind(a, b)] <- y
    Y[cbind(b, a)] <- y
    Y
  }
  eig_S <- function(y) eigen(Cs - Amap(y), symmetric = TRUE)

  ## Strictly feasible start: a large negative diagonal multiplier.
  y <- ifelse(isd, -(max(abs(Cs)) * q + 1), 0)
  ## Normalised objective scale: <C, Omega> lies in [-1, 1] for the
  ## normalised C used by the caller; anything past that is infeasibility.
  cap <- max(1, sum(abs(C))) * (1 + 1e-8) + 1e-8
  mu <- 1
  n_newton <- 0L
  exhausted <- FALSE
  repeat {
    ## Centering: minimise F(y) = -(1/mu) b'y - logdet S(y).
    stalled <- FALSE
    repeat {
      es <- eig_S(y)
      Q <- es$vectors %*% (t(es$vectors) / es$values)
      Q <- (Q + t(Q)) * 0.5
      tq <- 2 * wt * Q[cbind(a, b)]
      g <- bvec - mu * tq
      Hq <- 2 * outer(wt, wt) * (Q[a, b, drop = FALSE] * t(Q[a, b, drop = FALSE]) +
                                   Q[a, a, drop = FALSE] * Q[b, b, drop = FALSE])
      ## Jacobi-scaled eigen solve. As mu -> 0 the Newton matrix becomes
      ## ill-conditioned (S(y) nears singularity); that is the benign
      ## ill-conditioning of barrier methods, so drop the numerically null
      ## directions instead of failing.
      dh <- 1 / sqrt(pmax(diag(Hq), .Machine$double.xmin))
      eh <- eigen(Hq * outer(dh, dh), symmetric = TRUE)
      keep <- eh$values > 1e-13 * max(eh$values)
      if (!any(keep))
        .dynhr_abort("mom_se_bounds(): the semidefinite program's Newton ",
                     "system is degenerate.", class = "dynhr_se_bounds_solver")
      Vk <- eh$vectors[, keep, drop = FALSE]
      dlt <- dh * (Vk %*% (crossprod(Vk, dh * g) / eh$values[keep])) / mu
      dlt <- as.numeric(dlt)
      dec <- sum(g * dlt) / mu
      if (!is.finite(dec))
        .dynhr_abort("mom_se_bounds(): the semidefinite program diverged.",
                     class = "dynhr_se_bounds_solver")
      ## Centred enough: the dual point y is feasible whatever the centring,
      ## so b'y stays a valid bound; the centring only sets how tight it is.
      if (dec / 2 <= 1e-8) break
      if (stalled) break
      n_newton <- n_newton + 1L
      if (n_newton > max_iter) {
        exhausted <- TRUE
        break
      }
      ## Armijo on the CHANGE in F, formed without F itself: |F| ~ |b'y|/mu
      ## grows as mu -> 0 and its rounding would swamp the decrease. The
      ## change uses the REALISED step yn - y: when |y| is large a tiny step
      ## can round to no change at all, which must not count as progress.
      ld0 <- sum(log(es$values))
      st <- 1
      repeat {
        yn <- y + st * dlt
        ev <- eigen(Cs - Amap(yn), symmetric = TRUE, only.values = TRUE)$values
        ## A margin keeps S(yn) numerically PD for the full eigen solve.
        if (min(ev) > 1e-14 * max(abs(ev))) {
          dF <- -sum(bvec * (yn - y)) / mu - (sum(log(ev)) - ld0)
          if (dF <= -0.25 * st * dec && any(yn != y)) break
        }
        st <- st * 0.5
        if (st < 1e-10) break
      }
      ## No representable progress left at this mu: F is flat to rounding.
      if (st < 1e-10) break
      y <- yn
      ## Nearly centred but the full Newton step no longer passes the line
      ## search: the Newton direction is at its rounding floor (it degrades
      ## as S(y) nears singularity for small mu). Accept this centring after
      ## refreshing Q / dlt at the new y.
      stalled <- st < 1 && dec <= 1e-4
      ## Weak duality: b'y <= <sense*C, Omega> for every feasible Omega, and
      ## |<C, Omega>| <= cap. A dual value past the cap proves infeasibility.
      if (sum(bvec * y) > cap)
        .dynhr_abort("mom_se_bounds(): the known correlations in ",
                     "`cor_known` admit no positive semidefinite completion.",
                     class = "dynhr_se_bounds_infeasible")
    }
    if (exhausted || q * mu <= tol) break
    mu <- mu * 0.1
  }
  ## Primal from the last Newton system: Omega = mu (Q + Q A(dlt) Q) meets
  ## the known entries exactly (the constraints are linear in it) and is PSD
  ## up to the Newton decrement; clip, then restore the unit diagonal.
  Om <- mu * (Q + Q %*% Amap(dlt) %*% Q)
  Om <- (Om + t(Om)) * 0.5
  eo <- eigen(Om, symmetric = TRUE)
  Om <- eo$vectors %*% (t(eo$vectors) * pmax(eo$values, 0))
  dd <- sqrt(pmax(diag(Om), .Machine$double.eps))
  Om <- Om / outer(dd, dd)
  Om <- (Om + t(Om)) * 0.5
  ## The dual value is a valid bound whatever happened above; when the
  ## solver stopped early (a known block with no positive-definite
  ## completion makes the dual optimum unattained and the iteration stall)
  ## it is conservative but possibly not sharp. Say so.
  if (exhausted)
    .dynhr_warn("mom_se_bounds(): the semidefinite program stopped after ",
                max_iter, " Newton steps (duality gap ",
                signif(abs(sum(Cs * Om) - sum(bvec * y)), 3), " of the ",
                "squared worst case). The reported bound is still valid but ",
                "may not be sharp; this happens when the known correlations ",
                "admit only singular completions.",
                class = "dynhr_se_bounds_inexact", call. = FALSE)
  list(value = sense * sum(bvec * y), Omega = Om, mu = mu)
}
