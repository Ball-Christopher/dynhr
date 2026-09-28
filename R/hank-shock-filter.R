## R/hank-shock-filter.R
## --------------------------------------------------------------------------
## Sequence-space least-squares shock filter (Rigato 2026, "A least-squares
## filter for sequence-space models", ECB Working Paper 3191).
##
## A linear model known only through its impulse responses
##   x^{ij}_k = d y^i_{t+k} / d eps^j_t,   k = 0, ..., T_IRF - 1,
## maps a stacked shock history e (shock j from t = -(T_IRF - 1) to T - 1,
## so every observation's full truncated MA window is covered) to the stacked
## observables y (observable i, t = 0, ..., T - 1) by the Toeplitz system
##   y = X e + u,   e ~ N(0, Sigma),   u ~ N(0, Omega),
## where X = [X^{ij}] and each X^{ij} (T x (T + T_IRF - 1)) carries the IRF
## along its rows, starting from the main diagonal (paper eq. 3). The most
## likely shock history is the least-squares solution
##   e_hat = Sigma X' (X Sigma X' + Omega)^{-1} y          (paper eq. 4)
## which is also E(e | y) under Gaussianity (eq. 5; the best linear predictor
## otherwise), with posterior covariance
##   Var(e | y) = Sigma - Sigma X' (X Sigma X' + Omega)^{-1} X Sigma.
## Heteroskedasticity is a time-varying diagonal Sigma, measurement error is
## Omega, and a missing observation deletes its row of X (paper section 2).
##
## Kalman equivalence (the test oracle): with the truncated-MA companion form
## of R/hank-kalman.R and q = T_IRF, the shift-register initial state is
## (eps_{-1}, ..., eps_{-q}) ~ N(0, sigma^2 I) -- exactly the paper's prior on
## the pre-sample shocks -- so the Kalman smoother's smoothed shocks equal
## e_hat and its log-likelihood equals the Gaussian log-density of the
## observed y returned here.
## --------------------------------------------------------------------------


#' Sequence-space least-squares shock filter
#'
#' Recovers the most likely history of structural shocks of a linear
#' (sequence-space) model from observed data, given only the model's impulse
#' responses, via the closed form of Rigato (2026). The stacked shock
#' history \eqn{e} -- each shock from \eqn{t = -(T_{IRF}-1)} to \eqn{T-1}, so
#' every observation's truncated moving-average window is covered -- and the
#' stacked observables \eqn{y} satisfy \eqn{y = X e + u}, where \eqn{X} is the
#' block-Toeplitz matrix of impulse responses, \eqn{e \sim N(0, \Sigma)}
#' and the measurement error \eqn{u \sim N(0, \Omega)}. The filter returns
#' \deqn{\hat e = \Sigma X' (X \Sigma X' + \Omega)^{-1} y = E(e \mid y),}
#' the posterior covariance
#' \eqn{\Sigma - \Sigma X' (X \Sigma X' + \Omega)^{-1} X \Sigma}, and the
#' Gaussian log-likelihood of the observed data. Missing observations
#' (\code{NA}) delete the corresponding rows of \eqn{X}; time-varying shock
#' variances make \eqn{\Sigma} time-varying.
#'
#' Under Gaussian shocks the output coincides with the Kalman smoother on the
#' truncated-MA companion form of \code{\link{hank_state_space}} (with
#' \code{q = T_irf}), whose initial shift-register state is the same prior on
#' the pre-sample shocks. Under non-Gaussian shocks with finite variance it
#' is the best linear predictor of the shocks. The cost is one Cholesky
#' factorisation of a matrix of size (number of observed data points).
#'
#' @param model_or_irfs One of:
#'   \itemize{
#'   \item a \code{\link{hank_model}}: the impulse response of each
#'     \code{model$exogenous} shock to a unit innovation of its AR(1) driving
#'     process (persistence \code{rho}) is computed with
#'     \code{\link{hank_model_irf}};
#'   \item a \code{\link{hank_state_space}} object: its per-shock MA
#'     coefficients (\code{Theta_list}, truncated at its \code{q}) are used;
#'   \item a named list of impulse-response matrices, one per shock, each
#'     \code{T_IRF x n_obs} with row \code{k + 1} holding the response of
#'     every observable \code{k} periods after a unit shock (a plain vector
#'     is accepted when there is one observable).
#'   }
#' @param data \code{T x n_obs} matrix or data frame of demeaned
#'   observations; \code{NA} marks a missing observation. When it has column
#'   names they are matched to the observable names.
#' @param shock_sd Innovation standard deviations: a scalar, a length-\code{J}
#'   vector (one per shock, matched by name when named), or a \code{J x T}
#'   matrix of per-period values for the sample periods (pre-sample shocks
#'   then use the first column), or a \code{J x (T + T_irf - 1)} matrix
#'   covering the full history from \eqn{t = -(T_{irf}-1)}.
#' @param me_variance Measurement-error VARIANCES (the package-wide
#'   convention; standard deviations are their square roots): a scalar, a
#'   length-\code{n_obs} vector (matched by name when named) or a
#'   \code{T x n_obs} matrix. Default \code{0} (no measurement error; then
#'   the observed data must not be stochastically singular, i.e. no more
#'   observables than shocks can explain).
#' @param pre_sample Logical. \code{TRUE} (default, the paper's filter)
#'   filters the \code{T_irf - 1} pre-sample shocks of every shock from their
#'   \eqn{N(0, \sigma^2)} prior; \code{FALSE} fixes them at zero, i.e. the
#'   economy starts at its steady state.
#' @param obs_vars Character vector of observable names. Required for a
#'   \code{hank_model} (defaults to \code{colnames(data)}); optional
#'   otherwise, where it selects/reorders columns of the impulse responses.
#' @param rho For a \code{hank_model} only: named numeric AR(1) persistence
#'   of each exogenous driving process (default \code{0} for every shock,
#'   i.e. i.i.d. shock paths). Missing names default to 0.
#' @param T_irf Integer truncation horizon of the impulse responses (default:
#'   all rows supplied, \code{model$T_h} for a \code{hank_model}).
#' @param return_cov Logical; also return the full posterior covariance of
#'   the stacked shock history (default \code{FALSE}).
#'
#' @return A list of class \code{hank_shock_filter} with
#'   \describe{
#'   \item{\code{shocks}}{\code{T x J} matrix of smoothed shocks
#'     \eqn{E(\varepsilon_t \mid y)}, \eqn{t = 0, \dots, T-1}.}
#'   \item{\code{shocks_sd}}{\code{T x J} posterior standard deviations (the
#'     square root of the posterior covariance diagonal).}
#'   \item{\code{pre_sample_shocks}, \code{pre_sample_sd}}{\code{(T_irf - 1) x
#'     J} matrices for \eqn{t = -(T_{irf}-1), \dots, -1} (zero-row matrices
#'     when \code{pre_sample = FALSE}).}
#'   \item{\code{fitted}}{\code{T x n_obs} model-implied observables
#'     \eqn{X \hat e} (also at missing observations).}
#'   \item{\code{residuals}}{\code{data - fitted} (the smoothed measurement
#'     error; \code{NA} where the data are missing).}
#'   \item{\code{loglik}}{Gaussian log-likelihood of the observed data.}
#'   \item{\code{cov}}{Full posterior covariance of the stacked history
#'     (shock-major, time-minor), or \code{NULL} unless \code{return_cov}.}
#'   \item{\code{shock_names}, \code{obs_names}, \code{T_irf},
#'     \code{pre_sample}, \code{n_observed}}{Metadata.}
#'   }
#' @references Rigato, R. D. (2026). A least-squares filter for
#'   sequence-space models. ECB Working Paper No. 3191.
#' @seealso \code{\link{hank_state_space}}, \code{\link{hank_kalman_loglik}},
#'   \code{\link{hank_model_irf}}
#' @examples
#' ## Two shocks, two observables, known impulse responses.
#' irfs <- list(a = cbind(y1 = 0.8^(0:19), y2 = 0.5 * 0.8^(0:19)),
#'              b = cbind(y1 = -0.3 * 0.6^(0:19), y2 = 0.6^(0:19)))
#' set.seed(1)
#' y <- matrix(rnorm(80), 40, 2, dimnames = list(NULL, c("y1", "y2")))
#' f <- hank_filter_shocks(irfs, y, shock_sd = c(a = 1, b = 0.5), me_variance = 0.01)
#' head(f$shocks)
#' @export
hank_filter_shocks <- function(model_or_irfs, data, shock_sd, me_variance = 0,
                               pre_sample = TRUE, obs_vars = NULL, rho = NULL,
                               T_irf = NULL, return_cov = FALSE) {
  what <- "hank_filter_shocks"
  data <- as.matrix(data)
  if (!is.numeric(data) && !all(is.na(data)))
    .dynhr_abort(what, ": `data` must be numeric.", class = "dynhr_error_input")
  storage.mode(data) <- "double"
  if (!is.logical(pre_sample) || length(pre_sample) != 1L || is.na(pre_sample))
    .dynhr_abort(what, ": `pre_sample` must be TRUE or FALSE.",
                 class = "dynhr_error_input")

  irf <- .hsf_irf_list(model_or_irfs, data, obs_vars, rho, what)
  shock_names <- names(irf)
  obs_names   <- colnames(irf[[1L]])
  J <- length(shock_names)
  I <- length(obs_names)

  ## Align the data columns with the observables.
  if (!is.null(colnames(data))) {
    miss <- setdiff(obs_names, colnames(data))
    if (length(miss))
      .dynhr_abort(what, ": observable(s) missing from `data` columns: ",
                   paste(miss, collapse = ", "), ".",
                   class = "dynhr_error_input")
    data <- data[, obs_names, drop = FALSE]
  } else if (ncol(data) != I) {
    .dynhr_abort(what, ": `data` has ", ncol(data), " column(s) but there ",
                 "are ", I, " observable(s).", class = "dynhr_error_input")
  }
  colnames(data) <- obs_names
  n_T <- nrow(data)
  if (n_T < 1L)
    .dynhr_abort(what, ": `data` has no rows.", class = "dynhr_error_input")

  L_max <- nrow(irf[[1L]])
  if (is.null(T_irf)) T_irf <- L_max
  if (length(T_irf) != 1L || !is.finite(T_irf) || T_irf < 1 ||
      T_irf != round(T_irf) || T_irf > L_max)
    .dynhr_abort(what, ": `T_irf` must be an integer in [1, ", L_max, "].",
                 class = "dynhr_error_input")
  T_irf <- as.integer(T_irf)

  n_pre <- if (pre_sample) T_irf - 1L else 0L
  n_per <- n_T + n_pre                      # periods per shock in e
  sd_mat <- .hsf_shock_sd(shock_sd, shock_names, n_T, T_irf, n_pre, what)
  me_mat <- sqrt(.hsf_me_var(me_variance, obs_names, n_T, what))

  ## Block-Toeplitz X: row (i, t) holds x^{ij}_k at the column of shock j at
  ## time t - k (paper eq. 3); a pre-sample time below -n_pre is dropped
  ## (pre_sample = FALSE fixes those shocks at zero).
  n_e <- J * n_per
  X <- matrix(0, I * n_T, n_e)
  tt <- rep(seq_len(n_T), times = T_irf)            # 1-based t + 1
  kk <- rep(seq_len(T_irf) - 1L, each = n_T)
  cc <- tt - kk + n_pre                            # column within shock block
  keep <- cc >= 1L
  tt <- tt[keep]; kk <- kk[keep]; cc <- cc[keep]
  for (i in seq_len(I)) {
    for (j in seq_len(J)) {
      X[cbind((i - 1L) * n_T + tt, (j - 1L) * n_per + cc)] <-
        irf[[j]][kk + 1L, i]
    }
  }

  sig2 <- as.vector(t(sd_mat))^2                   # shock-major stacking
  yv   <- as.vector(data)                          # observable-major stacking
  obs  <- !is.na(yv)
  n_o  <- sum(obs)

  if (n_o == 0L) {
    e_hat <- numeric(n_e)
    post_var <- sig2
    loglik <- 0
    cov <- if (isTRUE(return_cov)) diag(sig2, n_e) else NULL
  } else {
    Xo <- X[obs, , drop = FALSE]
    XS <- Xo * rep(sig2, each = n_o)               # X Sigma
    S  <- tcrossprod(XS, Xo)
    diag(S) <- diag(S) + as.vector(me_mat)[obs]^2
    ## Pivoted Cholesky never errors on a singular S; singularity is judged on
    ## the pivots (chol() succeeding is not a singularity test).
    R <- suppressWarnings(chol(S, pivot = TRUE))
    piv <- attr(R, "pivot")
    if (attr(R, "rank") < n_o || .kf_F_singular(R, S))
      .dynhr_abort(what, ": the covariance of the observed data is singular ",
                   "-- the observables are stochastically singular given the ",
                   "shocks (more observables than shocks, or an observable ",
                   "that is an exact combination of others). Set a positive ",
                   "`me_variance` or drop observables.",
                   class = "dynhr_error_singular")
    ## S[piv, piv] = R'R.
    z <- backsolve(R, yv[obs][piv], transpose = TRUE)
    alpha <- numeric(n_o)
    alpha[piv] <- backsolve(R, z)
    e_hat <- as.vector(crossprod(XS, alpha))       # Sigma X' S^{-1} y
    W <- backsolve(R, XS[piv, , drop = FALSE], transpose = TRUE)
    post_var <- pmax(sig2 - colSums(W^2), 0)
    loglik <- -0.5 * (n_o * log(2 * pi) + 2 * sum(log(diag(R))) + sum(z^2))
    cov <- if (isTRUE(return_cov)) diag(sig2, n_e) - crossprod(W) else NULL
  }

  e_mat  <- matrix(e_hat, n_per, J, dimnames = list(NULL, shock_names))
  sd_post <- matrix(sqrt(post_var), n_per, J, dimnames = list(NULL, shock_names))
  in_sample <- n_pre + seq_len(n_T)
  pre_idx   <- seq_len(n_pre)
  fitted <- matrix(as.vector(X %*% e_hat), n_T, I,
                   dimnames = list(rownames(data), obs_names))
  if (!is.null(cov)) {
    lab <- paste0(rep(shock_names, each = n_per), "[",
                  rep(seq_len(n_per) - 1L - n_pre, times = J), "]")
    dimnames(cov) <- list(lab, lab)
  }

  structure(list(
    shocks            = e_mat[in_sample, , drop = FALSE],
    shocks_sd         = sd_post[in_sample, , drop = FALSE],
    pre_sample_shocks = e_mat[pre_idx, , drop = FALSE],
    pre_sample_sd     = sd_post[pre_idx, , drop = FALSE],
    fitted            = fitted,
    residuals         = data - fitted,
    loglik            = loglik,
    cov               = cov,
    shock_names       = shock_names,
    obs_names         = obs_names,
    T_irf             = T_irf,
    pre_sample        = pre_sample,
    n_observed        = n_o),
    class = "hank_shock_filter")
}


## Impulse responses as a named list of T_IRF x n_obs matrices (one per
## shock, common column names = observables).
.hsf_irf_list <- function(x, data, obs_vars, rho, what) {
  if (inherits(x, "hank_model")) {
    if (is.null(obs_vars)) obs_vars <- colnames(data)
    if (is.null(obs_vars))
      .dynhr_abort(what, ": give `obs_vars` (or column names on `data`) for ",
                   "a hank_model.", class = "dynhr_error_input")
    miss <- setdiff(obs_vars, names(x$G))
    if (length(miss))
      .dynhr_abort(what, ": observable(s) not produced by the model: ",
                   paste(miss, collapse = ", "), ".",
                   class = "dynhr_error_input")
    exo <- x$exogenous
    if (!is.null(rho) && (is.null(names(rho)) ||
                          length(setdiff(names(rho), exo))))
      .dynhr_abort(what, ": `rho` must be named by model$exogenous = {",
                   paste(exo, collapse = ", "), "}.",
                   class = "dynhr_error_input")
    specs <- setNames(lapply(exo, function(z) {
      r <- if (!is.null(rho) && z %in% names(rho)) rho[[z]] else 0
      list(rho = r, sigma = 1)
    }), exo)
    return(.hank_theta_list(x, specs, obs_vars))
  }
  if (inherits(x, "dsge_ss") && !is.null(x$Theta_list)) {
    out <- lapply(x$Theta_list, function(m) m[seq_len(x$q), , drop = FALSE])
    return(.hsf_select_obs(out, obs_vars, what))
  }
  if (!is.list(x) || !length(x))
    .dynhr_abort(what, ": `model_or_irfs` must be a hank_model, a ",
                 "hank_state_space object, or a list of impulse-response ",
                 "matrices.", class = "dynhr_error_input")
  if (is.null(names(x)) || anyNA(names(x)) || any(names(x) == ""))
    names(x) <- paste0("shock", seq_along(x))
  x <- lapply(x, function(m) if (is.null(dim(m))) matrix(m, ncol = 1L) else
    as.matrix(m))
  d1 <- dim(x[[1L]])
  cn <- colnames(x[[1L]])
  for (m in x) {
    if (!is.numeric(m) || any(!is.finite(m)))
      .dynhr_abort(what, ": impulse responses must be finite numeric.",
                   class = "dynhr_error_input")
    if (!identical(dim(m), d1) || !identical(colnames(m), cn))
      .dynhr_abort(what, ": every shock's impulse-response matrix must have ",
                   "the same dimensions and column names.",
                   class = "dynhr_error_input")
  }
  if (is.null(cn)) {
    cn <- if (!is.null(colnames(data)) && length(colnames(data)) == d1[2L])
      colnames(data) else paste0("obs", seq_len(d1[2L]))
    x <- lapply(x, function(m) { colnames(m) <- cn; m })
  }
  .hsf_select_obs(x, obs_vars, what)
}

.hsf_select_obs <- function(irf, obs_vars, what) {
  if (is.null(obs_vars)) return(irf)
  miss <- setdiff(obs_vars, colnames(irf[[1L]]))
  if (length(miss))
    .dynhr_abort(what, ": observable(s) not in the impulse responses: ",
                 paste(miss, collapse = ", "), ".",
                 class = "dynhr_error_input")
  lapply(irf, function(m) m[, obs_vars, drop = FALSE])
}

## Shock standard deviations as a J x n_per matrix over the filtered periods
## (the last n_T columns are the sample periods).
.hsf_shock_sd <- function(shock_sd, shock_names, n_T, T_irf, n_pre, what) {
  J <- length(shock_names)
  n_full <- n_T + T_irf - 1L
  if (is.matrix(shock_sd)) {
    if (nrow(shock_sd) != J || !(ncol(shock_sd) %in% c(n_T, n_full)))
      .dynhr_abort(what, ": a matrix `shock_sd` must be ", J, " x ", n_T,
                   " (sample periods) or ", J, " x ", n_full,
                   " (full history from t = -(T_irf - 1)).",
                   class = "dynhr_error_input")
    if (!is.null(rownames(shock_sd))) {
      if (!setequal(rownames(shock_sd), shock_names))
        .dynhr_abort(what, ": rownames of `shock_sd` must be the shock names {",
                     paste(shock_names, collapse = ", "), "}.",
                     class = "dynhr_error_input")
      shock_sd <- shock_sd[shock_names, , drop = FALSE]
    }
    full <- if (ncol(shock_sd) == n_full) shock_sd else
      cbind(shock_sd[, rep(1L, T_irf - 1L), drop = FALSE], shock_sd)
  } else {
    v <- .hsf_named_vec(shock_sd, shock_names, "shock_sd", what)
    full <- matrix(v, J, n_full)
  }
  if (any(!is.finite(full)) || any(full < 0))
    .dynhr_abort(what, ": `shock_sd` must be finite and non-negative.",
                 class = "dynhr_error_input")
  full[, (n_full - n_T - n_pre + 1L):n_full, drop = FALSE]
}

## Measurement-error standard deviations as a T x n_obs matrix.
.hsf_me_var <- function(me_variance, obs_names, n_T, what) {
  I <- length(obs_names)
  if (is.matrix(me_variance)) {
    if (nrow(me_variance) != n_T || ncol(me_variance) != I)
      .dynhr_abort(what, ": a matrix `me_variance` must be ", n_T, " x ", I, ".",
                   class = "dynhr_error_input")
    if (!is.null(colnames(me_variance))) {
      if (!setequal(colnames(me_variance), obs_names))
        .dynhr_abort(what, ": colnames of `me_variance` must be the observable ",
                     "names.", class = "dynhr_error_input")
      me_variance <- me_variance[, obs_names, drop = FALSE]
    }
    out <- me_variance
  } else {
    v <- .hsf_named_vec(me_variance, obs_names, "me_variance", what)
    out <- matrix(v, n_T, I, byrow = TRUE)
  }
  if (any(!is.finite(out)) || any(out < 0))
    .dynhr_abort(what, ": `me_variance` must be finite and non-negative.",
                 class = "dynhr_error_input")
  out
}

## Scalar or per-name vector -> vector in `nms` order.
.hsf_named_vec <- function(v, nms, arg, what) {
  if (!is.numeric(v))
    .dynhr_abort(what, ": `", arg, "` must be numeric.",
                 class = "dynhr_error_input")
  if (length(v) == 1L) return(rep(as.numeric(v), length(nms)))
  if (length(v) != length(nms))
    .dynhr_abort(what, ": `", arg, "` must be a scalar or have length ",
                 length(nms), ".", class = "dynhr_error_input")
  if (!is.null(names(v))) {
    if (!setequal(names(v), nms))
      .dynhr_abort(what, ": names of `", arg, "` must be {",
                   paste(nms, collapse = ", "), "}.",
                   class = "dynhr_error_input")
    v <- v[nms]
  }
  as.numeric(v)
}
