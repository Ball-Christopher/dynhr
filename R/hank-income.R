## R/hank-income.R
## --------------------------------------------------------------------------
## Idiosyncratic income discretization for heterogeneous-agent (HANK) blocks.
##
## Discretizes a log-AR(1) idiosyncratic productivity process
##
##   log e_{t+1} = rho * log e_t + eps,   eps ~ N(0, sigma_eps^2)
##
## into a finite Markov chain on productivity *levels* e (normalized so the
## stationary-distribution mean of e is 1), returning the transition matrix
## Pi and the stationary distribution pi.
##
## The Rouwenhorst (1995) method is used by default: unlike Tauchen it matches
## the persistence and unconditional variance of the AR(1) exactly for any N,
## and it is the method used by the `sequence-jacobian` reference package.  The
## grid/scaling convention here is chosen to match that package so its published
## household-Jacobian goldens are reproducible (see
## tests/testthat/test-hank-household.R and briefs/08).
##
## CONVENTION (matches sequence-jacobian utilities.discretize.markov_rouwenhorst):
##   `sigma` is the UNCONDITIONAL standard deviation of log e, i.e.
##   sigma = sigma_eps / sqrt(1 - rho^2).  A symmetric grid on [-1, 1] is scaled
##   so its stationary variance equals sigma^2, then exponentiated and
##   mean-normalized: e = exp(s) / sum(pi * exp(s)).
## --------------------------------------------------------------------------


# =============================================================================
# Rouwenhorst transition matrix
# =============================================================================

#' Rouwenhorst transition matrix for a two-parameter symmetric chain
#'
#' Builds the N-state Rouwenhorst transition matrix from the single crossing
#' probability \code{p = (1 + rho) / 2} by the standard recursive doubling
#' construction.
#'
#' @param n Integer >= 2: number of states.
#' @param p Numeric in (0, 1): stay/crossing probability.
#'
#' @return An \code{n x n} row-stochastic transition matrix.
#' @keywords internal
.hank_rouwenhorst_Pi <- function(n, p) {
  if (n < 2L) stop("Rouwenhorst requires n >= 2")
  Pi <- matrix(c(p, 1 - p, 1 - p, p), 2L, 2L, byrow = TRUE)
  if (n == 2L) return(Pi)
  for (k in 3L:n) {
    m  <- k - 1L
    P0 <- matrix(0, k, k)
    P0[1:m, 1:m]           <- P0[1:m, 1:m]           + p       * Pi
    P0[1:m, 2:k]           <- P0[1:m, 2:k]           + (1 - p) * Pi
    P0[2:k, 1:m]           <- P0[2:k, 1:m]           + (1 - p) * Pi
    P0[2:k, 2:k]           <- P0[2:k, 2:k]           + p       * Pi
    ## interior rows were double-counted -> halve them
    P0[2:m, ] <- P0[2:m, ] / 2
    Pi <- P0
  }
  Pi
}


#' Stationary distribution of a finite Markov transition matrix
#'
#' Left eigenvector of \code{Pi} for eigenvalue 1, normalized to sum to one.
#' Falls back to iteration if the eigen solve is ill-conditioned.
#'
#' @param Pi A row-stochastic transition matrix.
#' @param tol Convergence tolerance for the iterative fallback.
#' @param maxit Maximum iterations for the fallback.
#'
#' @return A probability vector (the stationary distribution).
#' @keywords internal
.hank_stationary <- function(Pi, tol = 1e-14, maxit = 100000L) {
  n <- nrow(Pi)
  ## Try the eigen route: stationary pi solves pi Pi = pi  <=>  Pi' v = v.
  ev <- eigen(t(Pi))
  k  <- which.min(abs(ev$values - 1))
  v  <- Re(ev$vectors[, k])
  if (all(is.finite(v)) && abs(sum(v)) > 1e-12) {
    v <- v / sum(v)
    if (all(v >= -1e-10)) {
      v[v < 0] <- 0
      return(v / sum(v))
    }
  }
  ## Iterative fallback.
  d <- rep(1 / n, n)
  for (i in seq_len(maxit)) {
    d_new <- as.numeric(crossprod(Pi, d))  # Pi' d
    if (max(abs(d_new - d)) < tol) { d <- d_new; break }
    d <- d_new
  }
  d / sum(d)
}


# =============================================================================
# Public: discretize a log-AR(1) income process
# =============================================================================

#' Discretize an idiosyncratic log-AR(1) income process (Rouwenhorst)
#'
#' Discretizes \eqn{\log e_{t+1} = \rho \log e_t + \varepsilon} into an
#' \code{n}-state Markov chain on productivity levels, normalized so the
#' stationary mean of \eqn{e} is one.
#'
#' @param rho Numeric in (-1, 1): AR(1) persistence of log income.
#' @param sigma Numeric > 0: UNCONDITIONAL standard deviation of log income
#'   (i.e. \eqn{\sigma_\varepsilon / \sqrt{1 - \rho^2}}).  This matches the
#'   \code{sequence-jacobian} convention; pass the innovation std divided by
#'   \eqn{\sqrt{1-\rho^2}} if you have the conditional std instead.
#' @param n Integer >= 2: number of income states (default 7).
#' @param method Character: only \code{"rouwenhorst"} is supported (Tauchen may
#'   be added later); other values error.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{e}}{Numeric length-\code{n}: income levels, stationary-mean 1.}
#'     \item{\code{Pi}}{\code{n x n} row-stochastic transition matrix.}
#'     \item{\code{pi}}{Numeric length-\code{n}: stationary distribution.}
#'     \item{\code{rho}, \code{sigma}, \code{n}, \code{method}}{Echoed inputs.}
#'   }
#'
#' @examples
#' inc <- hank_income_rouwenhorst(rho = 0.966, sigma = 0.5, n = 7)
#' sum(inc$pi * inc$e)   # ~ 1 (mean-normalized)
#' @export
hank_income_rouwenhorst <- function(rho, sigma, n = 7L,
                                    method = c("rouwenhorst")) {
  method <- match.arg(method)
  n <- as.integer(n)
  if (!is.finite(rho) || abs(rho) >= 1) stop("rho must be in (-1, 1)")
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be > 0")

  p  <- (1 + rho) / 2
  Pi <- .hank_rouwenhorst_Pi(n, p)
  pi <- .hank_stationary(Pi)

  ## Symmetric base grid, scaled so its *stationary* variance equals sigma^2.
  s0     <- seq(-1, 1, length.out = n)
  mean0  <- sum(pi * s0)
  var0   <- sum(pi * (s0 - mean0)^2)
  scale  <- if (var0 > 0) sigma / sqrt(var0) else 0
  s      <- s0 * scale                 # log-income grid (mean 0 by symmetry)

  ## Levels, normalized so stationary mean of e is 1.
  e_raw <- exp(s)
  e     <- e_raw / sum(pi * e_raw)

  list(e = e, Pi = Pi, pi = pi,
       rho = rho, sigma = sigma, n = n, method = method)
}


# =============================================================================
# Public: employment-state extension (HANK+SAM household income process)
# =============================================================================

#' Employment-augmented idiosyncratic income process (HANK+SAM household)
#'
#' Extends a Rouwenhorst productivity chain with a two-point employment state
#' \eqn{m \in \{E, U\}} (Ravn-Sterk / den Haan-Rendahl-Riegler HANK+SAM
#' class): employed households earn \eqn{w e}, unemployed households earn
#' replacement income \eqn{b_{ui} w e} (replacement-rate parameterization on
#' the same productivity level). The employment transition is
#' \deqn{P(U \to E) = f, \qquad P(E \to U) = s,}
#' independent of the productivity transition, so the combined chain has the
#' Kronecker structure \eqn{\Pi = \Pi_m(f, s) \otimes \Pi_e} on the
#' \code{2 * n} combined states.
#'
#' STATE ORDERING: employment OUTER, productivity INNER -- combined state
#' \code{(m, j)} has index \code{(m - 1) * n + j} with \code{m = 1} employed,
#' \code{m = 2} unemployed. Rows \code{1..n} of a \code{(2n) x n_a} policy
#' matrix are therefore the employed states; the stationary unemployment rate
#' is the mass on rows \code{n+1..2n} and equals \eqn{s / (s + f)} exactly
#' (independence of the two chains).
#'
#' NORMALIZATION: the productivity levels \code{e_prod} keep the
#' \code{\link{hank_income_rouwenhorst}} convention (stationary mean 1), and
#' the EFFECTIVE levels \code{e = c(e_prod, b_ui * e_prod)} are NOT
#' re-normalized: their stationary mean is \eqn{(1 - u) + u\, b_{ui}} with
#' \eqn{u = s/(s+f)}. This is deliberate -- \code{e} must stay FIXED when
#' \eqn{(f_t, s_t)} move along a transition path (only \code{Pi} responds to
#' the aggregate inputs), so no quantity may bake the steady-state
#' \eqn{(f, s)} into the income levels.
#'
#' The returned \code{Pi_fn(f, s)} rebuilds the combined transition matrix at
#' arbitrary transition-probability inputs; pass it (with
#' \code{Pi_inputs = list(f = f, s = s)}) to \code{\link{hank_het_block}} to
#' make \code{f} and \code{s} perturbable aggregate inputs of the household
#' block (time-varying job-finding/separation risk driving the fake-news
#' Jacobian columns; see \code{\link{hank_het_jacobian}}).
#'
#' @param f Numeric in (0, 1]: steady-state job-finding rate, \eqn{P(U \to E)}.
#' @param s Numeric in (0, 1): steady-state separation rate, \eqn{P(E \to U)}.
#' @param rho Numeric in (-1, 1): AR(1) persistence of log productivity.
#' @param sigma Numeric > 0: UNCONDITIONAL standard deviation of log
#'   productivity (same convention as \code{\link{hank_income_rouwenhorst}}).
#' @param n Integer >= 2: number of productivity states (default 7).
#' @param b_ui Numeric > 0: unemployment-insurance replacement rate (share of
#'   the household's employed income received while unemployed; default 0.5).
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{e}}{Numeric length-\code{2n}: EFFECTIVE income levels,
#'       \code{c(e_prod, b_ui * e_prod)} (employed block first).}
#'     \item{\code{Pi}}{\code{2n x 2n} row-stochastic combined transition
#'       matrix, \code{kronecker(Pi_m(f, s), Pi_e)}.}
#'     \item{\code{pi}}{Numeric length-\code{2n}: stationary distribution,
#'       \code{kronecker(c(1 - u, u), pi_e)}.}
#'     \item{\code{Pi_fn}}{Function \code{(f, s) -> 2n x 2n} combined
#'       transition matrix (the object the transition-probability Jacobian
#'       columns perturb).}
#'     \item{\code{u}}{Stationary unemployment rate \code{s / (s + f)}.}
#'     \item{\code{idx_E}, \code{idx_U}}{Integer index vectors of the
#'       employed / unemployed combined states (rows \code{1..n} /
#'       \code{n+1..2n}).}
#'     \item{\code{e_prod}, \code{Pi_e}, \code{pi_e}}{The underlying
#'       productivity chain (\code{\link{hank_income_rouwenhorst}} output).}
#'     \item{\code{f}, \code{s}, \code{b_ui}, \code{rho}, \code{sigma},
#'       \code{n}, \code{n_m}, \code{method}}{Echoed inputs / metadata
#'       (\code{n_m = 2L}, \code{method = "rouwenhorst_employment"}).}
#'   }
#'
#' @examples
#' inc <- hank_employment_income(f = 0.7, s = 0.05, rho = 0.9, sigma = 0.6,
#'                               n = 3, b_ui = 0.5)
#' sum(inc$pi[inc$idx_U])         # unemployment rate = 0.05 / 0.75
#' inc$u
#' @export
hank_employment_income <- function(f, s, rho, sigma, n = 7L, b_ui = 0.5) {
  if (!is.finite(f) || f <= 0 || f > 1)
    stop("f (job-finding rate) must be in (0, 1]")
  if (!is.finite(s) || s <= 0 || s >= 1)
    stop("s (separation rate) must be in (0, 1)")
  if (!is.finite(b_ui) || b_ui <= 0)
    stop("b_ui (replacement rate) must be > 0 (zero income at the borrowing ",
         "constraint makes the household problem infeasible)")

  prod <- hank_income_rouwenhorst(rho = rho, sigma = sigma, n = n)
  n    <- prod$n
  Pi_e <- prod$Pi
  pi_e <- prod$pi

  ## Combined transition at arbitrary (f, s): employment outer, productivity
  ## inner. Closes over the FIXED productivity chain Pi_e -- only the
  ## employment margin responds to the aggregate inputs.
  Pi_fn <- function(f, s) {
    if (!is.finite(f) || f <= 0 || f > 1)
      stop("Pi_fn: f (job-finding rate) must be in (0, 1]")
    if (!is.finite(s) || s <= 0 || s >= 1)
      stop("Pi_fn: s (separation rate) must be in (0, 1)")
    Pi_m <- matrix(c(1 - s, s,
                     f, 1 - f), 2L, 2L, byrow = TRUE)
    kronecker(Pi_m, Pi_e)
  }

  u  <- s / (s + f)
  Pi <- Pi_fn(f, s)
  pi <- as.numeric(kronecker(c(1 - u, u), pi_e))

  list(e = c(prod$e, b_ui * prod$e), Pi = Pi, pi = pi,
       Pi_fn = Pi_fn, u = u,
       idx_E = seq_len(n), idx_U = n + seq_len(n),
       e_prod = prod$e, Pi_e = Pi_e, pi_e = pi_e,
       f = f, s = s, b_ui = b_ui,
       rho = rho, sigma = sigma, n = n, n_m = 2L,
       method = "rouwenhorst_employment")
}


# =============================================================================
# Non-Gaussian innovation presets (analytic mean/var/skew/ex-kurt)
# =============================================================================

#' Analytic moments and CDF for a mean-zero mixture-of-normals innovation
#'
#' Builds a two-component normal mixture with weights \code{(p, 1 - p)},
#' component means \code{(mu1, mu2)} satisfying \code{p*mu1 + (1-p)*mu2 = 0}
#' and \code{mu1 - mu2 = mix_mu_gap}, and component sds \code{(s1, s2)} with
#' \code{s1 / s2 = mix_s_ratio}, then rescales the whole mixture (means and
#' sds together) by a single factor so its variance is exactly
#' \code{sigma_eps^2}. Skewness/excess-kurtosis are scale-invariant so they
#' are computed pre-rescale.
#'
#' Sign convention: \code{mix_mu_gap > 0} together with \code{mix_p < 0.5}
#' puts the small-probability component at the HIGH mean
#' (\code{mu1 = (1 - p) * mix_mu_gap > 0} with weight \code{p < 0.5}),
#' which gives POSITIVELY skewed eps. Flipping the sign of \code{mix_mu_gap}
#' (or moving \code{p} across 0.5) flips the skew sign.
#'
#' @param sigma_eps Numeric > 0: target innovation standard deviation.
#' @param mix_p Numeric in (0, 1): weight on component 1.
#' @param mix_mu_gap Numeric: pre-scale mean gap \code{mu1 - mu2}.
#' @param mix_s_ratio Numeric > 0: pre-scale sd ratio \code{s1 / s2}.
#'
#' @return A list with \code{p, mu1, mu2, s1, s2} (post-rescale, ready to use
#'   in the mixture CDF) and \code{mean, var, skewness, ex_kurt} (analytic
#'   innovation moments; \code{mean = 0}, \code{var = sigma_eps^2} exactly).
#' @keywords internal
.hank_mixture_moments <- function(sigma_eps, mix_p, mix_mu_gap, mix_s_ratio) {
  if (!is.finite(mix_p) || mix_p <= 0 || mix_p >= 1)
    stop("mix_p must be in (0, 1)")
  if (!is.finite(mix_s_ratio) || mix_s_ratio <= 0)
    stop("mix_s_ratio must be > 0")

  p <- mix_p
  ## Pre-scale (arbitrary unit): s2_0 = 1, s1_0 = mix_s_ratio.
  mu1_0 <- (1 - p) * mix_mu_gap
  mu2_0 <- -p * mix_mu_gap
  s1_0  <- mix_s_ratio
  s2_0  <- 1

  ## Pre-scale central moments (mean is exactly 0 by construction).
  var0 <- p * (s1_0^2 + mu1_0^2) + (1 - p) * (s2_0^2 + mu2_0^2)
  mu3  <- p * (mu1_0^3 + 3 * mu1_0 * s1_0^2) +
    (1 - p) * (mu2_0^3 + 3 * mu2_0 * s2_0^2)
  mu4  <- p * (mu1_0^4 + 6 * mu1_0^2 * s1_0^2 + 3 * s1_0^4) +
    (1 - p) * (mu2_0^4 + 6 * mu2_0^2 * s2_0^2 + 3 * s2_0^4)

  skewness <- if (var0 > 0) mu3 / var0^1.5 else 0
  ex_kurt  <- if (var0 > 0) mu4 / var0^2 - 3 else 0

  ## Rescale means and sds together so the variance is exactly sigma_eps^2.
  scale <- if (var0 > 0) sigma_eps / sqrt(var0) else 0
  list(p = p, mu1 = mu1_0 * scale, mu2 = mu2_0 * scale,
       s1 = s1_0 * scale, s2 = s2_0 * scale,
       mean = 0, var = sigma_eps^2,
       skewness = skewness, ex_kurt = ex_kurt)
}


#' Analytic moments and scale for a Student-t innovation
#'
#' Scales a standard Student-t(\code{df}) variate by
#' \code{sigma_eps / sqrt(df / (df - 2))} so the result has variance exactly
#' \code{sigma_eps^2}; symmetric (skewness 0) with excess kurtosis
#' \code{6 / (df - 4)} for \code{df > 4} (infinite/undefined otherwise).
#'
#' @param sigma_eps Numeric > 0: target innovation standard deviation.
#' @param t_df Numeric > 2: Student-t degrees of freedom.
#'
#' @return A list with \code{df}, \code{scale} (multiplier on a standard
#'   Student-t variate), and \code{mean, var, skewness, ex_kurt}.
#' @keywords internal
.hank_student_t_moments <- function(sigma_eps, t_df) {
  if (!is.finite(t_df) || t_df <= 2)
    stop("t_df must be > 2 (Student-t variance is undefined/infinite otherwise)")
  if (t_df <= 4)
    warning("t_df <= 4: Student-t excess kurtosis is infinite/undefined")

  scale   <- sigma_eps / sqrt(t_df / (t_df - 2))
  ex_kurt <- if (t_df > 4) 6 / (t_df - 4) else Inf
  list(df = t_df, scale = scale,
       mean = 0, var = sigma_eps^2, skewness = 0, ex_kurt = ex_kurt)
}


# =============================================================================
# Public: discretize a log-AR(1) income process (Tauchen, non-Gaussian)
# =============================================================================

#' Discretize an idiosyncratic log-AR(1) income process with non-Gaussian
#' innovations (Tauchen)
#'
#' Discretizes \eqn{\log e_{t+1} = \rho \log e_t + \varepsilon} into an
#' \code{n}-state Markov chain on productivity levels via the Tauchen (1986)
#' method, generalized to a non-Gaussian innovation CDF \eqn{F}. This is an
#' alternative to \code{\link{hank_income_rouwenhorst}} for studying the
#' effect of skewed and/or fat-tailed idiosyncratic income risk.
#'
#' Unlike Rouwenhorst, Tauchen does NOT match the target persistence
#' \code{rho} and unconditional standard deviation \code{sigma} exactly at
#' small \code{n}: the discretization is exact only in the limit
#' \code{n -> Inf} (finer grid + wider truncation). The chain-implied moments
#' (persistence, sd, skewness, excess kurtosis of the stationary distribution
#' of log e) converge to their targets as \code{n} grows; \strong{n >= 21 is
#' recommended for production use}. \code{n = 7} is fine for cheap
#' plumbing/smoke checks where discretization fidelity does not matter.
#'
#' The innovation \code{eps} is scaled so
#' \code{Var(eps) = sigma_eps^2 = sigma^2 * (1 - rho^2)}, where \code{sigma}
#' follows the same convention as \code{\link{hank_income_rouwenhorst}}: the
#' UNCONDITIONAL standard deviation of log e. This rescaling holds exactly
#' for any innovation shape (it follows from linearity of the AR(1)); higher
#' stationary moments (skewness, kurtosis) of log e are shape-dependent and
#' are only recovered in the chain as \code{n} grows.
#'
#' @param rho Numeric in (-1, 1): AR(1) persistence of log income.
#' @param sigma Numeric > 0: UNCONDITIONAL standard deviation of log income
#'   (same convention as \code{\link{hank_income_rouwenhorst}}).
#' @param n Integer >= 2: number of income states (default 7).
#' @param m Numeric > 0: grid half-width in units of \code{sigma} (default 3);
#'   the grid is \code{seq(-m*sigma, m*sigma, length.out = n)}.
#' @param innov Character: innovation shape, one of \code{"gaussian"}
#'   (pure validation arm; matches \code{hank_income_rouwenhorst} in spirit),
#'   \code{"mixture"} (two-component normal mixture, can be skewed), or
#'   \code{"student_t"} (symmetric, fat-tailed).
#' @param mix_p Numeric in (0, 1): weight on mixture component 1 (used when
#'   \code{innov = "mixture"}).
#' @param mix_mu_gap Numeric: pre-scale mean gap between the two mixture
#'   components (used when \code{innov = "mixture"}). \code{mix_mu_gap > 0}
#'   with \code{mix_p < 0.5} gives POSITIVELY skewed eps (a small-probability
#'   high-mean component); see \code{\link{.hank_mixture_moments}}.
#' @param mix_s_ratio Numeric > 0: pre-scale sd ratio of the two mixture
#'   components (used when \code{innov = "mixture"}).
#' @param t_df Numeric > 2: Student-t degrees of freedom (used when
#'   \code{innov = "student_t"}); a warning is issued if \code{t_df <= 4}
#'   since kurtosis is then infinite/undefined.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{e}}{Numeric length-\code{n}: income levels, stationary-mean 1.}
#'     \item{\code{Pi}}{\code{n x n} row-stochastic transition matrix.}
#'     \item{\code{pi}}{Numeric length-\code{n}: stationary distribution.}
#'     \item{\code{rho}, \code{sigma}, \code{n}, \code{method}, \code{innov}}{Echoed inputs (\code{method = "tauchen_nongaussian"}).}
#'     \item{\code{innov_moments}}{List with the ANALYTIC innovation
#'       \code{mean}, \code{var}, \code{skewness}, \code{ex_kurt} of the
#'       chosen preset (not the chain-implied moments of log e).}
#'   }
#'
#' @examples
#' inc <- hank_income_nongaussian(rho = 0.9, sigma = 0.5, n = 21,
#'                                innov = "mixture", mix_p = 0.2,
#'                                mix_mu_gap = 2)
#' sum(inc$pi * inc$e)   # ~ 1 (mean-normalized)
#' inc$innov_moments$skewness   # > 0 (mix_p < 0.5, mix_mu_gap > 0)
#' @export
hank_income_nongaussian <- function(rho, sigma, n = 7L, m = 3,
                                    innov = c("gaussian", "mixture", "student_t"),
                                    mix_p = 0.2, mix_mu_gap = 2, mix_s_ratio = 1,
                                    t_df = 5) {
  innov <- match.arg(innov)
  n <- as.integer(n)
  if (!is.finite(rho) || abs(rho) >= 1) stop("rho must be in (-1, 1)")
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be > 0")
  if (n < 2L) stop("n must be >= 2")
  if (!is.finite(m) || m <= 0) stop("m must be > 0")

  sigma_eps <- sigma * sqrt(1 - rho^2)

  ## Innovation preset: CDF F(x) and analytic moments.
  if (innov == "gaussian") {
    F_eps <- function(x) stats::pnorm(x, mean = 0, sd = sigma_eps)
    innov_moments <- list(mean = 0, var = sigma_eps^2, skewness = 0, ex_kurt = 0)
  } else if (innov == "mixture") {
    mm <- .hank_mixture_moments(sigma_eps, mix_p, mix_mu_gap, mix_s_ratio)
    F_eps <- function(x) {
      mm$p * stats::pnorm(x, mean = mm$mu1, sd = mm$s1) +
        (1 - mm$p) * stats::pnorm(x, mean = mm$mu2, sd = mm$s2)
    }
    innov_moments <- list(mean = mm$mean, var = mm$var,
                          skewness = mm$skewness, ex_kurt = mm$ex_kurt)
  } else {  # student_t
    tm <- .hank_student_t_moments(sigma_eps, t_df)
    F_eps <- function(x) stats::pt(x / tm$scale, df = tm$df)
    innov_moments <- list(mean = tm$mean, var = tm$var,
                          skewness = tm$skewness, ex_kurt = tm$ex_kurt)
  }

  ## Tauchen grid and transition matrix.
  s <- seq(-m * sigma, m * sigma, length.out = n)
  d <- if (n > 1L) s[2L] - s[1L] else 0

  Pi <- matrix(0, n, n)
  for (i in seq_len(n)) {
    mu_i <- rho * s[i]
    for (j in seq_len(n)) {
      if (j == 1L) {
        Pi[i, j] <- F_eps(s[1L] - mu_i + d / 2)
      } else if (j == n) {
        Pi[i, j] <- 1 - F_eps(s[n] - mu_i - d / 2)
      } else {
        Pi[i, j] <- F_eps(s[j] - mu_i + d / 2) - F_eps(s[j] - mu_i - d / 2)
      }
    }
  }
  ## Numerical hygiene: clip tiny negatives from CDF roundoff, renormalize.
  Pi[Pi < 0] <- 0
  Pi <- Pi / rowSums(Pi)

  pi_stat <- .hank_stationary(Pi)

  ## Levels, normalized so stationary mean of e is 1.
  e_raw <- exp(s)
  e     <- e_raw / sum(pi_stat * e_raw)

  list(e = e, Pi = Pi, pi = pi_stat,
       rho = rho, sigma = sigma, n = n, method = "tauchen_nongaussian",
       innov = innov, innov_moments = innov_moments)
}
