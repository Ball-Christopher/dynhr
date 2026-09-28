## R/hank-limited-info.R
## --------------------------------------------------------------------------
## Limited-information (single-block) estimation of a heterogeneous-agent model
## block: Liu, Plagborg-Moller & Tan (2026), "Limited-Information Estimation of
## Heterogeneous Agent Models", arXiv:2608.13953.
##
## THE MOMENT CONDITION (their eq. 4 / Assumption 1).  Linearising the
## aggregated household outputs y_t in the macro "sufficient statistics" x_t
## gives y_t - y_ss = E_t sum_k J_k(theta) (x_{t+k} - x_ss) + xi_{t-1}(theta),
## where J_k(theta) = dy_0 / dx_k is the DATE-0 ROW of the block's
## sequence-space Jacobian (response of today's aggregate to news about the
## input k periods ahead, with the predetermined distribution held fixed; the
## distribution's effect is xi_{t-1}).  An instrument z_t orthogonal to
## xi_{t-1} and to households' forecast errors turns this into
##
##     Cov(y_t, z_t) = sum_k J_k(theta) Cov(x_{t+k}, z_t),
##
## which restricts NOTHING about how x_t is generated -- no GE closure.
##
## THE ESTIMATOR (Definition 1).  g(theta) = (2 pi / T) sum_j vec{ S_yz(w_j) -
## J(w_j; theta) S_xz(w_j) }, J(w; theta) = sum_k e^{i k w} J_k(theta), with
## the cross-periodogram of the demeaned data at the Fourier frequencies
## w_j = 2 pi j / T.  Summing over j is an exact reordering of
##
##     g(theta) = vec( A - sum_{m=0}^{T-1} Jt_m(theta) B_m ),
##     A   = (2 pi / T) sum_j S_yz(w_j)                   (= lag-0 cross-cov),
##     B_m = (2 pi / T) sum_j e^{i m w_j} S_xz(w_j)       (= CIRCULAR lag-m
##                                                          Cov(x_{t+m}, z_t)),
##     Jt_m = sum_{k = m mod T} J_k                        (aliasing, H > T),
##
## so A and B are computed ONCE from the periodogram and every theta costs one
## SSJ evaluation plus a (dy x dx T) %*% (dx T x dz) product.  The moments-
## mode input (population / pre-estimated covariances) is the time-domain eq.
## (5) with the same (A, B) structure, which is why one code path serves both.
##
## Inference: the HAC sandwich of Section 3.3 (Newey-West, L_T = ceil(2.24
## T^{1/3})), the over-identification statistic of Section 3.4, and the
## Gaussian multiplier bootstrap of the periodogram of Appendix A.
## --------------------------------------------------------------------------


#' Limited-information minimum-distance estimation of a heterogeneous-agent
#' block
#'
#' Estimates the structural parameters of ONE model block -- typically the
#' heterogeneous household block of a HANK model -- without specifying the
#' rest of the economy, following Liu, Plagborg-Moller & Tan (2026). The
#' block enters only through its sequence-space Jacobians (SSJs)
#' \eqn{J_k(\theta)} with respect to a vector of macro "sufficient statistics"
#' \eqn{x_t} (e.g. returns and after-tax earnings); identified shocks
#' \eqn{z_t} serve as instruments, and the estimator fits the moment
#' conditions
#' \deqn{\mathrm{Cov}(y_t, z_t) = \sum_{k \ge 0} J_k(\theta)\,
#'   \mathrm{Cov}(x_{t+k}, z_t)}{Cov(y_t, z_t) = sum_k J_k(theta) Cov(x_{t+k}, z_t)}
#' for the observed aggregated block outputs \eqn{y_t} (e.g. cross-sectional
#' moments of consumption and savings). No process for \eqn{x_t}, and no list
#' of the shocks driving the data, is assumed.
#'
#' @section Input modes:
#' \describe{
#'   \item{Data mode (\code{y}, \code{x}, \code{z})}{The paper's estimator
#'     (Definition 1): the moment function is built from the
#'     cross-periodogram of the demeaned series,
#'     \eqn{\hat g(\theta) = (2\pi/T)\sum_{j=0}^{T-1}
#'     \mathrm{vec}\{\hat S_{yz}(\omega_j) - J(\omega_j;\theta)\hat
#'     S_{xz}(\omega_j)\}}, \eqn{J(\omega;\theta) = \sum_k e^{\iota k\omega}
#'     J_k(\theta)}. Standard errors are HAC (Section 3.3) or bootstrap
#'     (Appendix A); the over-identification test uses the HAC long-run
#'     variance \eqn{\hat\Omega} (Section 3.4).}
#'   \item{Moments mode (\code{moments})}{The time-domain estimating equation
#'     (eq. 5) evaluated at user-supplied covariances --
#'     \code{moments$cov_yz} (\eqn{d_y \times d_z}), \code{moments$cov_xz}
#'     (\eqn{d_x \times d_z \times K}, slice \code{k + 1} =
#'     \eqn{\mathrm{Cov}(x_{t+k}, z_t)}), the asymptotic covariance
#'     \code{moments$Omega} of \eqn{\sqrt{T}\,\hat g(\theta_0)}
#'     (\eqn{d_g \times d_g}, \eqn{d_g = d_y d_z}, \code{vec} order: outputs
#'     fastest) and the sample size \code{moments$n_obs}. Useful when the
#'     impulse responses come from elsewhere (local projections, a published
#'     table). The bootstrap needs the data and is not available here.}
#' }
#'
#' @section The SSJ function:
#' \code{ssj_fn(theta)} returns a numeric array of dimension
#' \code{c(d_y, d_x, H)} whose slice \code{[, , k + 1]} is \eqn{J_k(\theta)},
#' the response of the date-0 block outputs to date-\eqn{k} news about the
#' inputs with the predetermined distribution held fixed -- i.e. the FIRST ROW
#' of the block's sequence-space Jacobian. \code{\link{hank_limited_info_ssj}}
#' builds one from a household-block constructor. When its dimnames and the
#' column names of \code{y}/\code{x} are both present they must agree in
#' ORDER (a transposed output or input silently matches the wrong moments).
#' Return \code{NULL} (or a non-finite array) for an infeasible \code{theta};
#' the optimiser then rejects that step. With \eqn{H > T} the tail is aliased
#' onto lag \eqn{k \bmod T}, exactly as \eqn{J(\omega_j;\theta)} does at the
#' Fourier frequencies.
#'
#' @section Derivatives:
#' The moment function is LINEAR in the SSJs, so
#' \eqn{\partial\hat g/\partial\theta_p = -\mathrm{vec}\sum_m
#' (\partial\tilde J_m/\partial\theta_p) B_m} is exact given
#' \eqn{\partial J_k/\partial\theta}. Supply \code{dssj_fn(theta)} (an array
#' \code{c(d_y, d_x, H, d_theta)}) when analytic parameter derivatives exist;
#' otherwise they are central finite differences of \code{ssj_fn} (the
#' one-asset/two-asset het blocks have no analytic parameter derivative of
#' their fake-news Jacobian; see \code{\link{hank_model_dtheta}}).
#'
#' @section Optimiser:
#' A box-projected Levenberg-Marquardt (Gauss-Newton) iteration on
#' \eqn{\hat g' W \hat g} using the Jacobian above. The paper uses L-BFGS-B
#' with multiple starts; the objective is the same, but LM exploits the
#' least-squares structure and converges to optimiser precision in a handful
#' of SSJ evaluations. It is a LOCAL method: for a multimodal criterion start
#' from several \code{theta0}.
#'
#' @param ssj_fn Function \code{theta -> array(d_y, d_x, H)} of SSJs (see the
#'   section).
#' @param theta0 Named numeric starting value.
#' @param y,x,z Data mode: \code{T x d_y} block outputs, \code{T x d_x} macro
#'   sufficient statistics, \code{T x d_z} instruments (matrices or data
#'   frames; demeaned internally).
#' @param moments Moments mode: a list with \code{cov_yz}, \code{cov_xz},
#'   \code{Omega}, \code{n_obs} (see the section). Give either
#'   \code{moments} or all of \code{y}, \code{x}, \code{z}.
#' @param lower,upper Optional box bounds (named or in \code{theta0} order).
#' @param weighting \code{"diagonal"} (default; the paper's scale-invariant
#'   \eqn{\hat V_z^{-1} \otimes \hat V_y^{-1}} in data mode, and
#'   \eqn{\mathrm{diag}(\Omega)^{-1}} in moments mode), \code{"optimal"}
#'   (\eqn{\hat\Omega^{-1}}; two-step from the diagonal estimate in data
#'   mode), \code{"identity"}, or a \eqn{d_g \times d_g} symmetric positive
#'   definite matrix. The paper cautions that \code{"optimal"} can behave
#'   poorly when \eqn{d_g} is large relative to \eqn{T}.
#' @param se \code{"hac"} (default) or \code{"bootstrap"} (data mode only).
#'   The analytic interval under-covered in the paper's simulation; the
#'   bootstrap interval and J-test p-value were accurate.
#' @param n_boot Bootstrap replications.
#' @param level Confidence level of the reported intervals.
#' @param bandwidth Newey-West truncation \eqn{L_T}; \code{NULL} uses the
#'   paper's \eqn{\lceil 2.24\,T^{1/3}\rceil}.
#' @param dssj_fn Optional analytic derivative (see the section).
#' @param fd_step Relative central-difference step for
#'   \eqn{\partial J_k/\partial\theta}
#'   (\eqn{h_p = } \code{fd_step} \eqn{\times \max(|\theta_p|, 10^{-3})}).
#' @param maxit Maximum LM iterations.
#' @param tol Convergence tolerance on the relative step and on the relative
#'   criterion decrease.
#' @param seed Bootstrap seed (the caller's RNG stream is restored);
#'   \code{NULL} to use the current stream.
#'
#' @return An object of class \code{"dynhr_limited_info"}: \code{estimate},
#'   \code{se}, \code{vcov} (\eqn{\hat\Sigma/T}), \code{ci} (a
#'   \eqn{d_\theta \times 2} matrix; bootstrap intervals are
#'   \eqn{[2\hat\theta - q_{1-\alpha/2}, 2\hat\theta - q_{\alpha/2}]}),
#'   \code{se_hac}, \code{j_stat}, \code{j_df}, \code{j_pvalue} (asymptotic
#'   \eqn{\chi^2}), \code{j_pvalue_boot} (bootstrap, or \code{NA}),
#'   \code{moment_fit} (a data frame: data vs fitted moments, their
#'   difference, HAC standard error and t-statistic), \code{G}, \code{W},
#'   \code{Omega}, \code{criterion} (\eqn{T\hat g'W\hat g}), \code{n_obs},
#'   \code{H}, \code{weighting}, \code{se_method}, \code{bandwidth},
#'   \code{convergence} (\code{"converged"}; otherwise \code{"boundary"},
#'   \code{"stalled"} or \code{"maxit"}, with a warning), \code{iterations},
#'   \code{n_ssj_eval}, \code{mode},
#'   and \code{boot} (the draws and bootstrap J statistics, or \code{NULL}).
#'
#' @references
#' Liu, L., Plagborg-Moller, M. & Tan, N. M. P. (2026). Limited-Information
#' Estimation of Heterogeneous Agent Models. arXiv:2608.13953 [econ.EM].
#'
#' Auclert, A., Bardoczy, B., Rognlie, M. & Straub, L. (2021). Using the
#' Sequence-Space Jacobian to Solve and Estimate Heterogeneous-Agent Models.
#' Econometrica 89(5), 2375-2408.
#'
#' Meyer, M. & Paparoditis, E. (2023). A frequency domain bootstrap for general
#' multivariate stationary processes. Bernoulli 29(3), 2367-2391.
#'
#' @seealso \code{\link{hank_limited_info_ssj}}, \code{\link{hank_het_jacobian}},
#'   \code{\link{method_of_moments}}, \code{\link{match_irfs}}
#' @examples
#' ## Example 1 of the paper: y_t = beta x_t + gamma E_t x_{t+1} + xi_{t-1},
#' ## so J_0 = beta, J_1 = gamma. Two instruments move x_t and x_{t+1}
#' ## differently, which identifies (beta, gamma).
#' set.seed(1)
#' T <- 400; e <- matrix(rnorm((T + 60) * 3), ncol = 3)
#' ma <- function(v, rho) stats::filter(v, rho, method = "recursive")
#' x_all <- ma(e[, 1], 0.9) + ma(e[, 2], 0.2) + ma(e[, 3], 0.5)
#' ex1 <- 0.9 * ma(e[, 1], 0.9) + 0.2 * ma(e[, 2], 0.2) + 0.5 * ma(e[, 3], 0.5)
#' y_all <- 0.6 * x_all + 0.3 * ex1 + 0.4 * c(0, e[-nrow(e), 3])
#' keep <- 61:(T + 60)
#' ssj <- function(th) array(c(th[["beta"]], th[["gamma"]]), c(1, 1, 2))
#' fit <- hank_limited_info_md(ssj, c(beta = 0.5, gamma = 0.5),
#'                             y = y_all[keep], x = x_all[keep],
#'                             z = e[keep, 1:2])
#' fit
#' @export
hank_limited_info_md <- function(ssj_fn, theta0, y = NULL, x = NULL, z = NULL,
                                 moments = NULL, lower = NULL, upper = NULL,
                                 weighting = c("diagonal", "optimal",
                                               "identity"),
                                 se = c("hac", "bootstrap"),
                                 n_boot = 199L, level = 0.95,
                                 bandwidth = NULL, dssj_fn = NULL,
                                 fd_step = 1e-4, maxit = 100L, tol = 1e-10,
                                 seed = 1L) {
  if (!is.function(ssj_fn))
    .dynhr_abort("hank_limited_info_md(): `ssj_fn` must be a function ",
                 "theta -> array(d_y, d_x, H).")
  if (!is.numeric(theta0) || !length(theta0) || any(!is.finite(theta0)) ||
      is.null(names(theta0)) || any(!nzchar(names(theta0))))
    .dynhr_abort("hank_limited_info_md(): `theta0` must be a finite, fully ",
                 "named numeric vector.")
  se <- match.arg(se)
  W_user <- NULL
  if (is.matrix(weighting)) {
    W_user <- weighting
    weighting <- "user"
  } else {
    weighting <- match.arg(weighting)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    .dynhr_abort("hank_limited_info_md(): `level` must be in (0, 1).")
  if (!is.null(dssj_fn) && !is.function(dssj_fn))
    .dynhr_abort("hank_limited_info_md(): `dssj_fn` must be NULL or a function.")

  p_names <- names(theta0)
  d_th <- length(theta0)
  lo <- .li_bound(lower, p_names, -Inf, "lower")
  up <- .li_bound(upper, p_names, Inf, "upper")
  if (any(lo >= up))
    .dynhr_abort("hank_limited_info_md(): every `lower` must be below `upper`.")
  theta0 <- pmin(pmax(theta0, lo), up)

  ## ---- Data -> (A, B) -----------------------------------------------------
  data_mode <- is.null(moments)
  if (data_mode) {
    if (is.null(y) || is.null(x) || is.null(z))
      .dynhr_abort("hank_limited_info_md(): give either `moments` or all of ",
                   "`y`, `x`, `z`.")
    y <- .li_as_matrix(y, "y")
    x <- .li_as_matrix(x, "x")
    z <- .li_as_matrix(z, "z")
    T_obs <- nrow(y)
    if (nrow(x) != T_obs || nrow(z) != T_obs)
      .dynhr_abort("hank_limited_info_md(): `y`, `x`, `z` must have the same ",
                   "number of rows (", T_obs, ", ", nrow(x), ", ", nrow(z), ").")
    if (T_obs < 8L)
      .dynhr_abort("hank_limited_info_md(): need at least 8 observations.")
    yc <- sweep(y, 2L, colMeans(y))
    xc <- sweep(x, 2L, colMeans(x))
    zc <- sweep(z, 2L, colMeans(z))
    AB <- .li_spectral_moments(yc, xc, zc)
    A <- AB$A
    B <- AB$B
    y_names <- colnames(y)
    x_names <- colnames(x)
    z_names <- colnames(z)
  } else {
    if (!is.null(y) || !is.null(x) || !is.null(z))
      .dynhr_abort("hank_limited_info_md(): give either `moments` or data ",
                   "(`y`, `x`, `z`), not both.")
    if (se == "bootstrap")
      .dynhr_abort("hank_limited_info_md(): se = \"bootstrap\" resamples the ",
                   "periodogram and needs the time series `y`, `x`, `z`; ",
                   "moments mode supports se = \"hac\" (with the supplied ",
                   "Omega) only.")
    mm <- .li_check_moments(moments)
    A <- mm$cov_yz
    B <- mm$cov_xz
    T_obs <- mm$n_obs
    y_names <- rownames(A)
    x_names <- dimnames(B)[[1L]]
    z_names <- colnames(A)
  }
  d_y <- nrow(A); d_z <- ncol(A); d_x <- dim(B)[1L]; K <- dim(B)[3L]
  if (dim(B)[2L] != d_z)
    .dynhr_abort("hank_limited_info_md(): cov_xz has ", dim(B)[2L],
                 " instrument columns but cov_yz has ", d_z, ".")
  d_g <- d_y * d_z
  if (d_g < d_th)
    .dynhr_abort("hank_limited_info_md(): ", d_g, " moments (d_y * d_z = ",
                 d_y, " * ", d_z, ") cannot identify ", d_th, " parameters ",
                 "(order condition; Section 2.3 of Liu, Plagborg-Moller & ",
                 "Tan 2026).")
  ## Names the USER gave are checked against the SSJ dimnames; the defaults
  ## below only label the output.
  y_chk <- y_names
  x_chk <- x_names
  if (is.null(y_names)) y_names <- paste0("y", seq_len(d_y))
  if (is.null(x_names)) x_names <- paste0("x", seq_len(d_x))
  if (is.null(z_names)) z_names <- paste0("z", seq_len(d_z))
  Bmat <- .li_bmat(B)

  ## ---- SSJ evaluation, cached by theta ------------------------------------
  cache <- new.env(parent = emptyenv())
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  ssj_at <- function(theta) {
    key <- paste(sprintf("%.17g", theta), collapse = "|")
    hit <- get0(key, envir = cache, inherits = FALSE)
    if (!is.null(hit)) return(hit$J)
    names(theta) <- p_names
    J <- ssj_fn(theta)
    counter$n <- counter$n + 1L
    J <- .li_check_ssj(J, d_y, d_x, y_chk, x_chk)
    assign(key, list(J = J), envir = cache)
    J
  }
  moment_at <- function(theta, AA = A, BB = Bmat) {
    J <- ssj_at(theta)
    if (is.null(J)) return(NULL)
    .li_g(AA, BB, .li_alias(J, K))
  }
  ## d vec(J_k)/d theta_p: analytic when supplied, else central FD in theta.
  dssj_at <- function(theta) {
    if (!is.null(dssj_fn)) {
      names(theta) <- p_names
      dJ <- dssj_fn(theta)
      H0 <- dim(ssj_at(theta))[3L]
      if (!is.numeric(dJ) || !identical(as.integer(dim(dJ)),
                                        c(d_y, d_x, H0, d_th)) ||
          any(!is.finite(dJ)))
        .dynhr_abort("hank_limited_info_md(): `dssj_fn` must return a finite ",
                     "array of dim c(", d_y, ", ", d_x, ", ", H0, ", ", d_th,
                     ").")
      return(dJ)
    }
    J0 <- ssj_at(theta)
    dJ <- array(0, c(dim(J0), d_th))
    for (p in seq_len(d_th)) {
      h <- fd_step * max(abs(theta[p]), 1e-3)
      tp <- theta; tm <- theta
      hp <- min(h, up[p] - theta[p]); hm <- min(h, theta[p] - lo[p])
      tp[p] <- theta[p] + hp; tm[p] <- theta[p] - hm
      Jp <- if (hp > 0) ssj_at(tp) else J0
      Jm <- if (hm > 0) ssj_at(tm) else J0
      if (is.null(Jp) || is.null(Jm) || hp + hm <= 0)
        .dynhr_abort("hank_limited_info_md(): `ssj_fn` is infeasible at a ",
                     "finite-difference neighbour of theta = (",
                     paste(format(theta), collapse = ", "), ") in '",
                     p_names[p], "'; shrink `fd_step` or tighten the bounds.")
      dJ[, , , p] <- (Jp - Jm) / (hp + hm)
    }
    dJ
  }
  G_at <- function(theta, BB = Bmat) {
    dJ <- dssj_at(theta)
    G <- matrix(0, d_g, d_th)
    for (p in seq_len(d_th))
      G[, p] <- -as.numeric(matrix(.li_alias(dJ[, , , p], K, d_y, d_x), d_y) %*%
                              BB)
    G
  }

  if (is.null(moment_at(theta0)))
    .dynhr_abort("hank_limited_info_md(): `ssj_fn` is infeasible at `theta0`.")
  H_ssj <- dim(ssj_at(theta0))[3L]

  ## ---- Scores and HAC long-run variance -----------------------------------
  if (data_mode) {
    L_T <- if (is.null(bandwidth)) as.integer(ceiling(2.24 * T_obs^(1 / 3)))
           else as.integer(bandwidth)
    if (length(L_T) != 1L || is.na(L_T) || L_T < 1L)
      .dynhr_abort("hank_limited_info_md(): `bandwidth` must be a positive ",
                   "integer.")
    omega_at <- function(theta, yy = yc, xx = xc, zz = zc)
      .li_hac(.li_scores(yy, xx, zz, ssj_at(theta)), L_T)
  } else {
    L_T <- NA_integer_
    ## Same formals as the data-mode closure above (R CMD check flags one
    ## local name bound to functions with different formals); the data
    ## arguments are unused (and never evaluated) in moment-matrix mode.
    omega_at <- function(theta, yy = yc, xx = xc, zz = zc) mm$Omega
  }

  ## ---- Weight matrix -------------------------------------------------------
  W <- switch(weighting,
    identity = diag(d_g),
    user = {
      if (!is.numeric(W_user) || !identical(dim(W_user), c(d_g, d_g)))
        .dynhr_abort("hank_limited_info_md(): a matrix `weighting` must be ",
                     d_g, " x ", d_g, ".")
      (W_user + t(W_user)) / 2
    },
    ## diagonal and the first step of optimal
    if (data_mode) kronecker(diag(1 / apply(zc, 2L, stats::var), d_z),
                             diag(1 / apply(yc, 2L, stats::var), d_y))
    else diag(1 / diag(mm$Omega), d_g))
  if (weighting == "optimal" && !data_mode)
    W <- .mom_omega_inverse(mm$Omega, ridge = 0)$W

  fit <- .li_lm(theta0, moment_at, G_at, W, lo, up, maxit, tol)
  steps <- 1L
  if (weighting == "optimal" && data_mode) {
    W <- .mom_omega_inverse(omega_at(fit$theta), ridge = 0)$W
    it1 <- fit$iterations
    fit <- .li_lm(fit$theta, moment_at, G_at, W, lo, up, maxit, tol)
    fit$iterations <- fit$iterations + it1
    steps <- 2L
  }
  theta_hat <- setNames(fit$theta, p_names)
  if (fit$convergence != "converged")
    .dynhr_warn("hank_limited_info_md(): the optimiser stopped with status '",
                fit$convergence, "' after ", fit$iterations, " iterations; ",
                "check the bounds and try other starting values.",
                class = "dynhr_limited_info_convergence")
  g_hat <- fit$g
  G <- G_at(theta_hat)
  Omega <- omega_at(theta_hat)
  mom_names <- paste0(rep(y_names, d_z), "|", rep(z_names, each = d_y))

  ## ---- Sandwich (Section 3.3) and J test (Section 3.4) --------------------
  inf <- .li_inference(G, W, Omega, g_hat, T_obs, d_th)
  if (!inf$identified)
    .dynhr_warn("hank_limited_info_md(): G'WG is singular at the estimate -- ",
                "the parameters are not locally identified by these ",
                "instruments (Assumption 6 of Liu, Plagborg-Moller & Tan ",
                "2026: the instruments must move the sufficient statistics ",
                "differently across horizons). Standard errors are NA.",
                class = "dynhr_limited_info_unidentified")
  vc <- inf$Sigma / T_obs
  dimnames(vc) <- list(p_names, p_names)
  se_hac <- setNames(sqrt(pmax(diag(vc), 0)), p_names)
  if (!inf$identified) se_hac[] <- NA_real_
  zq <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(theta_hat - zq * se_hac, theta_hat + zq * se_hac)

  ## ---- Bootstrap (Appendix A) ---------------------------------------------
  boot <- NULL
  j_p_boot <- NA_real_
  se_out <- se_hac
  if (se == "bootstrap") {
    .local_seed(seed)
    boot <- .li_bootstrap(yc, xc, zc, theta_hat, n_boot, W, ssj_at,
                          moment_at, G_at, omega_at, lo, up, maxit, tol,
                          T_obs, d_th, K)
    se_out <- setNames(apply(boot$theta, 2L, stats::sd), p_names)
    a <- 1 - level
    q_lo <- apply(boot$theta, 2L, stats::quantile, probs = a / 2, names = FALSE)
    q_hi <- apply(boot$theta, 2L, stats::quantile, probs = 1 - a / 2,
                  names = FALSE)
    ci <- cbind(2 * theta_hat - q_hi, 2 * theta_hat - q_lo)
    if (is.finite(inf$j_stat))
      j_p_boot <- mean(boot$j_stat >= inf$j_stat)
  }
  ci_names <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                            trim = TRUE), " %")
  dimnames(ci) <- list(p_names, ci_names)

  ## ---- Moment fit table ----------------------------------------------------
  data_m <- as.numeric(A)
  mom_se <- sqrt(pmax(diag(Omega), 0) / T_obs)
  moment_fit <- data.frame(
    moment = mom_names,
    output = rep(y_names, d_z),
    instrument = rep(z_names, each = d_y),
    data = data_m,
    fitted = data_m - g_hat,
    diff = g_hat,
    se = mom_se,
    t_stat = ifelse(mom_se > 0, g_hat / mom_se, NA_real_),
    stringsAsFactors = FALSE)
  G <- unname(G); W <- unname(W); Omega <- unname(Omega)
  dimnames(G) <- list(mom_names, p_names)
  dimnames(W) <- dimnames(Omega) <- list(mom_names, mom_names)

  out <- list(
    estimate = theta_hat, se = se_out, vcov = vc, ci = ci, se_hac = se_hac,
    j_stat = inf$j_stat, j_df = inf$j_df, j_pvalue = inf$j_p,
    j_pvalue_boot = j_p_boot, moment_fit = moment_fit,
    G = G, W = W, Omega = Omega, criterion = T_obs * fit$Q,
    n_obs = T_obs, H = H_ssj, n_moments = d_g, n_params = d_th,
    weighting = weighting, steps = steps, se_method = se,
    bandwidth = L_T, level = level, convergence = fit$convergence,
    iterations = fit$iterations, n_ssj_eval = counter$n,
    mode = if (data_mode) "data" else "moments", start = theta0,
    boot = boot)
  class(out) <- "dynhr_limited_info"
  out
}


#' @rdname hank_limited_info_md
#' @param digits Significant digits for \code{print}.
#' @param ... Ignored.
#' @export
print.dynhr_limited_info <- function(x, digits = 4L, ...) {
  cat("=== dynhr limited-information block estimate (Liu, Plagborg-Moller & Tan 2026) ===\n")
  cat(sprintf("  Mode: %s | T = %d | moments = %d | params = %d | H = %d\n",
              x$mode, as.integer(x$n_obs), x$n_moments, x$n_params,
              as.integer(x$H)))
  cat(sprintf("  Weighting: %s | SE: %s | convergence: %s (%d it, %d SSJ evals)\n",
              x$weighting, x$se_method, x$convergence,
              as.integer(x$iterations), as.integer(x$n_ssj_eval)))
  tab <- cbind(estimate = x$estimate, se = x$se, x$ci)
  print(signif(tab, digits))
  if (is.finite(x$j_stat)) {
    cat(sprintf("  Over-identification: J = %.4g, df = %d, p = %.4g",
                x$j_stat, as.integer(x$j_df), x$j_pvalue))
    if (is.finite(x$j_pvalue_boot))
      cat(sprintf(" (bootstrap p = %.4g)", x$j_pvalue_boot))
    cat("\n")
  } else {
    cat("  Over-identification: not available (exactly identified)\n")
  }
  invisible(x)
}


#' Limited-information SSJ function from a household-block constructor
#'
#' Wraps a function that builds a heterogeneous-agent block at a parameter
#' vector into the \code{ssj_fn} that \code{\link{hank_limited_info_md}}
#' needs: the block's sequence-space Jacobian is computed at horizon \code{H}
#' and its FIRST ROW is returned as the array \code{J[o, i, k + 1]} =
#' \eqn{\partial y^o_0 / \partial x^i_k} -- the response of today's aggregate
#' output to news about input \code{i} \eqn{k} periods ahead, holding the
#' predetermined distribution fixed, which is the \eqn{J_k(\theta)} of Liu,
#' Plagborg-Moller & Tan (2026, eq. 1). Later Jacobian rows (anticipation
#' filtered through the distribution) are the \eqn{\xi_{t-1}} term that the
#' instruments difference out, so they are not used.
#'
#' The inputs' DATING must match the data: \code{x[t, i]} is the input the
#' block sees at date \code{t} (e.g. the one-asset block's \code{r} is the
#' return paid on assets carried into \code{t}).
#'
#' @param block_fn Function \code{theta -> block} (e.g. a closure over
#'   \code{\link{hank_het_block}} holding the grids, income process and
#'   steady-state prices fixed). Return \code{NULL} when \code{theta} is
#'   infeasible.
#' @param H SSJ horizon (number of lags \eqn{k = 0, \ldots, H-1}); choose it
#'   so the SSJs have decayed (the paper suggests \eqn{H \gg T} is cheap).
#' @param inputs,outputs Input (sufficient-statistic) and output names, in
#'   the column order of \code{x} and \code{y}.
#' @param jacobian The block's SSJ routine, called as
#'   \code{jacobian(block, H, inputs = inputs, outputs = outputs)} and
#'   returning \code{J[[output]][[input]]} (\code{H x H}):
#'   \code{\link{hank_het_jacobian}} (default), \code{\link{hank_het2_jacobian}},
#'   \code{\link{hank_het3_jacobian}}, ...
#'
#' @return A function \code{theta -> array(length(outputs), length(inputs), H)}
#'   with dimnames \code{list(outputs, inputs, NULL)}.
#' @references Liu, L., Plagborg-Moller, M. & Tan, N. M. P. (2026).
#'   Limited-Information Estimation of Heterogeneous Agent Models.
#'   arXiv:2608.13953 [econ.EM].
#' @seealso \code{\link{hank_limited_info_md}}
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.6, 3)
#' ag  <- hank_asset_grid(60, 50, 0)
#' ssj <- hank_limited_info_ssj(function(th)
#'   hank_het_block(ag, inc$Pi, inc$e, beta = th[["beta"]],
#'                  eis = th[["eis"]], r = 0.01, w = 1),
#'   H = 40, inputs = c("r", "w"), outputs = c("A", "C"))
#' J <- ssj(c(beta = 0.95, eis = 0.5))
#' dim(J)        # 2 x 2 x 40
#' J["C", "r", 1:5]
#' @export
hank_limited_info_ssj <- function(block_fn, H, inputs, outputs,
                                  jacobian = hank_het_jacobian) {
  if (!is.function(block_fn) || !is.function(jacobian))
    .dynhr_abort("hank_limited_info_ssj(): `block_fn` and `jacobian` must be ",
                 "functions.")
  H <- as.integer(H)
  if (length(H) != 1L || is.na(H) || H < 1L)
    .dynhr_abort("hank_limited_info_ssj(): `H` must be a positive integer.")
  if (!is.character(inputs) || !length(inputs) || !is.character(outputs) ||
      !length(outputs))
    .dynhr_abort("hank_limited_info_ssj(): `inputs` and `outputs` must be ",
                 "non-empty character vectors.")
  function(theta) {
    blk <- block_fn(theta)
    if (is.null(blk)) return(NULL)
    Jl <- jacobian(blk, H, inputs = inputs, outputs = outputs)
    out <- array(0, c(length(outputs), length(inputs), H),
                 dimnames = list(outputs, inputs, NULL))
    for (o in outputs) for (i in inputs) out[o, i, ] <- Jl[[o]][[i]][1L, ]
    out
  }
}


## ---- internals ----------------------------------------------------------

#' Resolve a named/positional bound to a full named vector.
#' @noRd
.li_bound <- function(b, p_names, default, what) {
  out <- setNames(rep(default, length(p_names)), p_names)
  if (is.null(b)) return(out)
  if (!is.numeric(b))
    .dynhr_abort("hank_limited_info_md(): `", what, "` must be numeric.")
  if (is.null(names(b))) {
    if (length(b) != length(p_names))
      .dynhr_abort("hank_limited_info_md(): unnamed `", what, "` must have ",
                   "one entry per parameter.")
    out[] <- b
  } else {
    bad <- setdiff(names(b), p_names)
    if (length(bad))
      .dynhr_abort("hank_limited_info_md(): `", what, "` names unknown ",
                   "parameter(s) ", paste0("'", bad, "'", collapse = ", "), ".")
    out[names(b)] <- b
  }
  out
}

#' Coerce a data argument to a finite numeric matrix.
#' @noRd
.li_as_matrix <- function(v, what) {
  m <- if (is.data.frame(v)) as.matrix(v) else v
  if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
  if (!is.numeric(m) || length(dim(m)) != 2L || any(!is.finite(m)))
    .dynhr_abort("hank_limited_info_md(): `", what, "` must be a finite ",
                 "numeric matrix (T rows).")
  storage.mode(m) <- "double"
  m
}

#' Validate the moments-mode list.
#' @noRd
.li_check_moments <- function(m) {
  need <- c("cov_yz", "cov_xz", "Omega", "n_obs")
  if (!is.list(m) || length(setdiff(need, names(m))))
    .dynhr_abort("hank_limited_info_md(): `moments` must be a list with ",
                 paste0("`", need, "`", collapse = ", "), ".")
  A <- m$cov_yz
  if (is.null(dim(A))) A <- matrix(A, ncol = 1L)
  B <- m$cov_xz
  if (!is.numeric(A) || any(!is.finite(A)) || length(dim(A)) != 2L)
    .dynhr_abort("hank_limited_info_md(): moments$cov_yz must be a finite ",
                 "d_y x d_z matrix.")
  if (!is.numeric(B) || length(dim(B)) != 3L || any(!is.finite(B)))
    .dynhr_abort("hank_limited_info_md(): moments$cov_xz must be a finite ",
                 "d_x x d_z x K array (slice k+1 = Cov(x_{t+k}, z_t)).")
  d_g <- length(A)
  Om <- m$Omega
  if (!is.numeric(Om) || !identical(dim(Om), c(d_g, d_g)) ||
      any(!is.finite(Om)))
    .dynhr_abort("hank_limited_info_md(): moments$Omega must be a finite ",
                 d_g, " x ", d_g, " matrix (the asymptotic covariance of ",
                 "sqrt(T) * vec(moments), outputs fastest).")
  Om <- (Om + t(Om)) / 2
  if (min(eigen(Om, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    .dynhr_abort("hank_limited_info_md(): moments$Omega must be positive ",
                 "definite.")
  n <- m$n_obs
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0)
    .dynhr_abort("hank_limited_info_md(): moments$n_obs must be a positive ",
                 "number.")
  list(cov_yz = A, cov_xz = B, Omega = Om, n_obs = as.numeric(n))
}

#' Validate one SSJ array (NULL / non-finite -> NULL = infeasible).
#' @noRd
.li_check_ssj <- function(J, d_y, d_x, y_names, x_names) {
  if (is.null(J)) return(NULL)
  if (!is.numeric(J) || length(dim(J)) != 3L ||
      dim(J)[1L] != d_y || dim(J)[2L] != d_x)
    .dynhr_abort("hank_limited_info_md(): `ssj_fn` must return an array of ",
                 "dim c(", d_y, ", ", d_x, ", H) (outputs x inputs x lags).")
  if (any(!is.finite(J))) return(NULL)
  dn <- dimnames(J)
  chk <- function(have, want, what) {
    if (!is.null(have) && !is.null(want) && !identical(as.character(have),
                                                       as.character(want)))
      .dynhr_abort("hank_limited_info_md(): the SSJ ", what, " (",
                   paste(have, collapse = ", "), ") do not match the data ",
                   "columns (", paste(want, collapse = ", "), ") in order.",
                   class = "dynhr_limited_info_name_mismatch")
  }
  chk(dn[[1L]], y_names, "outputs")
  chk(dn[[2L]], x_names, "inputs")
  J
}

#' Moments A (d_y x d_z) and B (d_x x d_z x T) from the cross-periodogram.
#'
#' S_ab(w_j) = F_a(w_j) F_b(w_j)^* / (2 pi T) with F_a the DFT of the
#' demeaned series; A = (2 pi/T) sum_j S_yz(w_j) and
#' B_m = (2 pi/T) sum_j e^{i m w_j} S_xz(w_j), m = 0..T-1 (an inverse DFT).
#' R's fft indexes time from 0 rather than 1; the phase cancels in F_a F_b^*.
#' @noRd
.li_spectral_moments <- function(yc, xc, zc) {
  T_obs <- nrow(yc)
  Fy <- stats::mvfft(yc); Fx <- stats::mvfft(xc); Fz <- stats::mvfft(zc)
  d_x <- ncol(xc); d_z <- ncol(zc)
  A <- Re(crossprod(Fy, Conj(Fz))) / T_obs^2
  prod_xz <- Fx[, rep(seq_len(d_x), d_z), drop = FALSE] *
    Conj(Fz[, rep(seq_len(d_z), each = d_x), drop = FALSE])
  Bm <- Re(stats::mvfft(prod_xz, inverse = TRUE)) / T_obs^2   # T x (d_x d_z)
  B <- array(0, c(d_x, d_z, T_obs))
  for (b in seq_len(d_z)) for (i in seq_len(d_x))
    B[i, b, ] <- Bm[, (b - 1L) * d_x + i]
  dimnames(A) <- list(colnames(yc), colnames(zc))
  dimnames(B) <- list(colnames(xc), colnames(zc), NULL)
  list(A = A, B = B)
}

#' Flatten B (d_x x d_z x K) to the (d_x K) x d_z matrix whose row (i, m),
#' i fastest, matches matrix(Jt, d_y) columns.
#' @noRd
.li_bmat <- function(B) {
  d <- dim(B)
  matrix(aperm(B, c(1L, 3L, 2L)), d[1L] * d[3L], d[2L])
}

#' Alias an SSJ array (d_y x d_x x H) onto K lags: Jt_m = sum_{k = m mod K} J_k.
#' @noRd
.li_alias <- function(J, K, d_y = dim(J)[1L], d_x = dim(J)[2L]) {
  J <- array(J, c(d_y, d_x, length(J) / (d_y * d_x)))
  H <- dim(J)[3L]
  out <- array(0, c(d_y, d_x, K))
  for (k in seq_len(H)) {
    m <- ((k - 1L) %% K) + 1L
    out[, , m] <- out[, , m] + J[, , k]
  }
  out
}

#' Moment vector vec(A - sum_m Jt_m B_m).
#' @noRd
.li_g <- function(A, Bmat, Jt) {
  as.numeric(A - matrix(Jt, nrow(A)) %*% Bmat)
}

#' Per-observation scores psi_t = z_t (x) (y_t - sum_{k=0}^{T-t} J_k x_{t+k})
#' (Section 3.3; demeaned data).
#' @noRd
.li_scores <- function(yc, xc, zc, J) {
  T_obs <- nrow(yc)
  e <- yc
  for (k in seq_len(min(dim(J)[3L], T_obs)) - 1L) {
    rows <- seq_len(T_obs - k)
    e[rows, ] <- e[rows, , drop = FALSE] -
      xc[rows + k, , drop = FALSE] %*% t(matrix(J[, , k + 1L], dim(J)[1L]))
  }
  d_y <- ncol(yc); d_z <- ncol(zc)
  e[, rep(seq_len(d_y), d_z), drop = FALSE] *
    zc[, rep(seq_len(d_z), each = d_y), drop = FALSE]
}

#' Newey-West HAC: sum_{|h| <= L} (1 - |h|/L) Gamma(h), Gamma of demeaned
#' scores divided by T.
#' @noRd
.li_hac <- function(psi, L) {
  T_obs <- nrow(psi)
  pc <- sweep(psi, 2L, colMeans(psi))
  Om <- crossprod(pc) / T_obs
  for (h in seq_len(min(L, T_obs) - 1L)) {
    w <- 1 - h / L
    Gh <- crossprod(pc[(h + 1L):T_obs, , drop = FALSE],
                    pc[seq_len(T_obs - h), , drop = FALSE]) / T_obs
    Om <- Om + w * (Gh + t(Gh))
  }
  (Om + t(Om)) / 2
}

#' Sandwich variance and over-identification statistic.
#' @noRd
.li_inference <- function(G, W, Omega, g, T_obs, d_th) {
  Hm <- crossprod(G, W %*% G)
  Hm <- (Hm + t(Hm)) / 2
  ev <- eigen(Hm, symmetric = TRUE, only.values = TRUE)$values
  identified <- all(is.finite(ev)) && max(ev) > 0 && min(ev) > 1e-10 * max(ev)
  Sigma <- matrix(NA_real_, d_th, d_th)
  if (identified) {
    Hi <- solve(Hm)
    Sigma <- Hi %*% crossprod(G, W %*% Omega %*% W %*% G) %*% Hi
    Sigma <- (Sigma + t(Sigma)) / 2
  }
  oi <- .mom_omega_inverse(Omega, ridge = 0)
  j_df <- oi$rank - d_th
  j_stat <- NA_real_
  j_p <- NA_real_
  if (j_df > 0L && identified) {
    Oi <- oi$W
    M <- crossprod(G, Oi %*% G)
    P <- Oi - Oi %*% G %*% solve(M, crossprod(G, Oi))
    j_stat <- T_obs * as.numeric(crossprod(g, P %*% g))
    j_p <- stats::pchisq(j_stat, df = j_df, lower.tail = FALSE)
  }
  list(Sigma = Sigma, identified = identified, j_stat = j_stat,
       j_df = j_df, j_p = j_p)
}

#' Box-projected Levenberg-Marquardt on Q(theta) = g' W g.
#'
#' Stops when the Gauss-Newton predicted decrease grad' (G'WG)^{-1} grad falls
#' below `tol * Q` (the relative Newton decrement), when a step is below `tol`
#' relative, or when Q has fallen to 1e-20 of its start (noise-free moments).
#' The decrement test is what terminates a noisy fit: an SSJ is only as smooth
#' in theta as its steady-state solve (EGM tolerance), so the relative step
#' stalls near 1e-8 at the optimum while the decrement is ~1e-11 * Q. A fit
#' that can make no further progress for another reason reports "boundary"
#' (a parameter on a bound) or "stalled".
#' @noRd
.li_lm <- function(theta, moment_at, G_at, W, lo, up, maxit, tol) {
  Qf <- function(g) if (is.null(g)) Inf else as.numeric(crossprod(g, W %*% g))
  g <- moment_at(theta)
  Q <- Qf(g)
  Q_start <- Q
  lambda <- 1e-3
  conv <- "maxit"
  it <- 0L
  while (it < maxit) {
    if (Q <= 1e-20 * Q_start) { conv <- "converged"; break }
    it <- it + 1L
    G <- G_at(theta)
    Hm <- crossprod(G, W %*% G)
    grad <- as.numeric(crossprod(G, W %*% g))
    dH <- pmax(diag(Hm), 1e-12 * max(diag(Hm), 1e-300))
    decr <- as.numeric(crossprod(grad, solve(Hm + 1e-12 * diag(dH, length(dH)),
                                            grad)))
    if (is.finite(decr) && decr <= tol * Q) { conv <- "converged"; break }
    accepted <- FALSE
    while (lambda <= 1e12) {
      step <- -solve(Hm + lambda * diag(dH, length(dH)), grad)
      th_new <- pmin(pmax(theta + step, lo), up)
      g_new <- moment_at(th_new)
      Q_new <- Qf(g_new)
      if (Q_new < Q) { accepted <- TRUE; break }
      lambda <- lambda * 10
    }
    if (!accepted) {
      ## No damping reduces Q. If the undamped Gauss-Newton step is itself
      ## below sqrt(tol) relative, theta is stationary to the SSJ's own
      ## numerical noise (the noise-free-moments case, where the decrement
      ## stays ~Q at the floor): converged. Otherwise a genuine stall.
      gn <- -solve(Hm + 1e-12 * diag(dH, length(dH)), grad)
      conv <- if (max(abs(gn) / pmax(abs(theta), 1e-3)) < sqrt(tol)) "converged"
              else if (any(theta <= lo | theta >= up)) "boundary"
              else "stalled"
      break
    }
    rel_step <- max(abs(th_new - theta) / pmax(abs(theta), 1e-3))
    theta <- th_new; g <- g_new; Q <- Q_new
    lambda <- max(lambda / 10, 1e-12)
    if (rel_step < tol) { conv <- "converged"; break }
  }
  list(theta = theta, g = g, Q = Q, iterations = it, convergence = conv)
}

#' Epanechnikov-smoothed periodogram of zeta at the Fourier frequencies
#' (Appendix A; bandwidth B = T^-0.2, the zero frequency excluded).
#' @noRd
.li_smoothed_spectrum <- function(zeta) {
  T_obs <- nrow(zeta)
  d <- ncol(zeta)
  Fz <- stats::mvfft(zeta)
  S <- array(0 + 0i, c(d, d, T_obs))
  for (j in seq_len(T_obs))
    S[, , j] <- (Fz[j, ] %o% Conj(Fz[j, ])) / (2 * pi * T_obs)
  w <- 2 * pi * (seq_len(T_obs) - 1L) / T_obs
  Bw <- T_obs^(-0.2)
  f <- array(0 + 0i, c(d, d, T_obs))
  for (j in seq_len(T_obs)) {
    dist <- abs(w[j] - w[-1L])
    dist <- pmin(dist, 2 * pi - dist)
    kw <- pmax(1 - (dist / Bw)^2, 0)
    if (sum(kw) <= 0) kw[which.min(dist)] <- 1
    kw <- kw / sum(kw)
    acc <- matrix(0 + 0i, d, d)
    for (i in which(kw > 0)) acc <- acc + kw[i] * S[, , i + 1L]
    f[, , j] <- acc
  }
  f
}

#' Gaussian multiplier bootstrap of the periodogram (Appendix A): pseudo
#' Fourier transforms F_j ~ N_c(0, f(w_j)), recentred moments, re-optimised
#' estimate and over-identification statistic per draw.
#' @noRd
.li_bootstrap <- function(yc, xc, zc, theta_hat, n_boot, W, ssj_at,
                          moment_at, G_at, omega_at, lo, up, maxit, tol,
                          T_obs, d_th, K) {
  d_y <- ncol(yc); d_x <- ncol(xc); d_z <- ncol(zc)
  zeta <- cbind(yc, xc, zc)
  d <- ncol(zeta)
  iy <- seq_len(d_y); ix <- d_y + seq_len(d_x); iz <- d_y + d_x + seq_len(d_z)
  f <- .li_smoothed_spectrum(zeta)
  n_half <- (T_obs - 1L) %/% 2L
  ## Recentring term (2 pi / T) sum_j {f_yz - J(w_j; theta_hat) f_xz} over the
  ## frequencies the pseudo-periodogram populates (not 0 or pi), so the
  ## bootstrap moment has mean exactly zero at theta_hat.
  js <- c(seq_len(n_half), T_obs - seq_len(n_half)) + 1L
  f_yz <- matrix(0 + 0i, d_y, d_z)
  for (j in js) f_yz <- f_yz + f[iy, iz, j]
  A_f <- Re(f_yz) * 2 * pi / T_obs
  prod_xz <- matrix(0 + 0i, T_obs, d_x * d_z)
  for (j in js) prod_xz[j, ] <- as.vector(f[ix, iz, j])
  Bm_f <- Re(stats::mvfft(prod_xz, inverse = TRUE)) * 2 * pi / T_obs
  B_f <- array(0, c(d_x, d_z, T_obs))
  for (b in seq_len(d_z)) for (i in seq_len(d_x))
    B_f[i, b, ] <- Bm_f[, (b - 1L) * d_x + i]
  recentre <- A_f - matrix(.li_alias(ssj_at(theta_hat), K), d_y) %*%
    .li_bmat(B_f)
  ## F ~ N_c(0, f) (circularly symmetric, E[F F^*] = f, E[F F'] = 0) through
  ## its real 2d-vector form: (Re F, Im F) ~ N(0, [[R, -I], [I, R]] / 2) for
  ## f = R + iI Hermitian. Kept real on purpose: a complex Hermitian eigen()
  ## is avoidable and has crashed some LAPACK builds.
  roots <- lapply(seq_len(n_half) + 1L, function(j) {
    fj <- (f[, , j] + Conj(t(f[, , j]))) / 2
    R <- Re(fj); I <- Im(fj)
    M <- rbind(cbind(R, -I), cbind(I, R)) / 2
    eg <- eigen((M + t(M)) / 2, symmetric = TRUE)
    eg$vectors %*% diag(sqrt(pmax(eg$values, 0)), 2L * d)
  })
  theta_b <- matrix(NA_real_, n_boot, d_th, dimnames = list(NULL,
                                                             names(theta_hat)))
  j_b <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    Fd <- matrix(0 + 0i, T_obs, d)
    for (jj in seq_len(n_half)) {
      v <- as.vector(roots[[jj]] %*% stats::rnorm(2L * d))
      Fj <- complex(real = v[seq_len(d)], imaginary = v[d + seq_len(d)])
      ## DFT scale: periodogram = F_dft F_dft^* / (2 pi T) = Fj Fj^*.
      Fd[jj + 1L, ] <- Fj * sqrt(2 * pi * T_obs)
      Fd[T_obs - jj + 1L, ] <- Conj(Fd[jj + 1L, ])
    }
    zb <- Re(stats::mvfft(Fd, inverse = TRUE)) / T_obs
    ABb <- .li_spectral_moments(zb[, iy, drop = FALSE], zb[, ix, drop = FALSE],
                                zb[, iz, drop = FALSE])
    A_b <- ABb$A - recentre
    B_b <- .li_bmat(ABb$B)
    mom_b <- function(theta) moment_at(theta, AA = A_b, BB = B_b)
    G_b <- function(theta) G_at(theta, BB = B_b)
    fit_b <- .li_lm(theta_hat, mom_b, G_b, W, lo, up, maxit, tol)
    theta_b[b, ] <- fit_b$theta
    Om_b <- omega_at(fit_b$theta, zb[, iy, drop = FALSE],
                     zb[, ix, drop = FALSE], zb[, iz, drop = FALSE])
    j_b[b] <- .li_inference(G_b(fit_b$theta), W, Om_b, fit_b$g, T_obs,
                            d_th)$j_stat
  }
  list(theta = theta_b, j_stat = j_b)
}
