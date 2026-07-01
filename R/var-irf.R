## R/var-irf.R
## --------------------------------------------------------------------------
## Reduced-form VAR(p) estimation + structural (Cholesky) IRFs + residual-
## bootstrap confidence bands.
##
## Purpose: generate the *empirical* IRF targets (and their sampling bands)
## that match_irfs() matches a structural model against, in the standard
## Ireland / RBC-IRF-matching workflow.  estimate_var() fits a reduced-form
## VAR by OLS equation-by-equation; var_irf() recovers structural IRFs via a
## Cholesky (recursive) identification of the residual covariance; and
## var_irf_bootstrap() produces percentile bands by residual resampling.
##
## Base R (stats) only.
## --------------------------------------------------------------------------


## --------------------------------------------------------------------------
## Internal: build the regressor matrix (lagged Y + deterministics) and the
## LHS for a VAR(p).  Y is T x K (rows = time, cols = variables).
## Returns list(X, Ylhs, K, p, type, det_cols).
## --------------------------------------------------------------------------
.var_design <- function(Y, p, type) {
  Y  <- as.matrix(Y)
  Tn <- nrow(Y); K <- ncol(Y)
  if (Tn <= p + 1L) stop("estimate_var: not enough observations for the lag order.")
  Ylhs <- Y[(p + 1L):Tn, , drop = FALSE]        # (T-p) x K
  neff <- nrow(Ylhs)

  ## Lag block: [y_{t-1}, ..., y_{t-p}], each y a K-vector.
  Xlags <- matrix(0, neff, K * p)
  for (lag in seq_len(p)) {
    Xlags[, ((lag - 1L) * K + 1L):(lag * K)] <-
      Y[(p + 1L - lag):(Tn - lag), , drop = FALSE]
  }

  ## Deterministic columns.
  det <- switch(type,
                none  = matrix(0, neff, 0),
                const = matrix(1, neff, 1, dimnames = list(NULL, "const")),
                trend = cbind(const = 1, trend = seq_len(neff)),
                stop("estimate_var: `type` must be 'const', 'none', or 'trend'."))

  X <- cbind(Xlags, det)
  list(X = X, Ylhs = Ylhs, K = K, p = p, type = type, ndet = ncol(det))
}


#' Estimate a reduced-form VAR(p) by OLS
#'
#' Fits a reduced-form vector autoregression
#' \deqn{y_t = c + A_1 y_{t-1} + \dots + A_p y_{t-p} + u_t,}
#' equation-by-equation by ordinary least squares (the OLS estimator is
#' efficient here because every equation shares the same regressors), and
#' returns the coefficient matrices together with the residual covariance.
#'
#' @param Y A \code{T x K} numeric matrix or data frame: rows are time periods
#'   (in order), columns are the \code{K} endogenous variables.  Column names,
#'   if present, are carried through to the IRFs.
#' @param p Lag order (positive integer).
#' @param type Deterministic terms: \code{"const"} (default, an intercept),
#'   \code{"none"}, or \code{"trend"} (intercept + linear trend).
#'
#' @return A list of class \code{var_fit} with elements:
#'   \describe{
#'     \item{A}{Length-\code{p} list of \code{K x K} autoregressive coefficient
#'       matrices \eqn{A_1,\dots,A_p} (row = equation, column = lagged variable).}
#'     \item{coef}{Full coefficient matrix \code{(K*p + ndet) x K} as returned
#'       by the OLS (rows = regressors, columns = equations).}
#'     \item{Sigma}{\code{K x K} residual covariance (MLE-style divisor
#'       \code{neff}).}
#'     \item{residuals}{\code{(T-p) x K} OLS residuals.}
#'     \item{p, type, K, var_names, nobs}{Bookkeeping.}
#'   }
#'
#' @examples
#' set.seed(1)
#' Tn <- 300
#' Y <- matrix(0, Tn, 2)
#' for (t in 2:Tn) Y[t, ] <- 0.6 * Y[t - 1, ] + rnorm(2)
#' fit <- estimate_var(Y, p = 1)
#' fit$A[[1]]   # ~ diag(0.6)
#'
#' @seealso \code{\link{var_irf}}, \code{\link{var_irf_bootstrap}}
#' @export
estimate_var <- function(Y, p = 1L, type = c("const", "none", "trend")) {
  type <- match.arg(type)
  p    <- as.integer(p)
  if (p < 1L) stop("estimate_var: `p` must be a positive integer.")

  var_names <- colnames(as.matrix(Y))
  d  <- .var_design(Y, p, type)
  X  <- d$X; Yl <- d$Ylhs; K <- d$K
  neff <- nrow(Yl)
  if (is.null(var_names)) var_names <- paste0("y", seq_len(K))

  ## OLS: B = (X'X)^{-1} X'Y, all equations at once.
  XtX  <- crossprod(X)
  coef <- solve(XtX, crossprod(X, Yl))            # (Kp+ndet) x K
  resid <- Yl - X %*% coef
  Sigma <- crossprod(resid) / neff                # K x K, divisor neff (MLE)

  colnames(coef) <- var_names
  dimnames(Sigma) <- list(var_names, var_names)

  ## Pull out the AR matrices A_1..A_p (rows of coef are stacked lags).
  A <- vector("list", p)
  for (lag in seq_len(p)) {
    block <- coef[((lag - 1L) * K + 1L):(lag * K), , drop = FALSE]  # K x K
    ## block[i, j] = effect of variable i (lag) on equation j; transpose so
    ## A_lag[eq, var] is the conventional orientation.
    A[[lag]] <- t(block)
    dimnames(A[[lag]]) <- list(var_names, var_names)
  }

  structure(
    list(A = A, coef = coef, Sigma = Sigma, residuals = resid,
         p = p, type = type, K = K, ndet = d$ndet,
         var_names = var_names, nobs = neff),
    class = "var_fit")
}


## --------------------------------------------------------------------------
## Internal: structural IRFs from AR matrices + a structural impact matrix B
## (such that u_t = B e_t, e_t unit-variance).  Returns a list keyed by shock
## name; each entry is a (horizon x K) matrix (rows = horizon 1..horizon,
## cols = variables) -- the format match_irfs() consumes.
##
## Recursion on the companion form: Phi_0 = I; Phi_h = sum_{j=1}^{min(h,p)}
## A_j Phi_{h-j}.  IRF of variable to structural shock s at horizon h is
## (Phi_h B)[, s].  Horizon index 1 == impact (Phi_0 B).
## --------------------------------------------------------------------------
.var_irf_from_coef <- function(A, B, horizon, var_names) {
  K <- nrow(B); p <- length(A)
  ## Reduced-form MA (Phi) matrices, h = 0 .. horizon-1.
  Phi <- vector("list", horizon)
  Phi[[1]] <- diag(K)                              # Phi_0
  for (h in 2:horizon) {
    acc <- matrix(0, K, K)
    for (j in seq_len(min(h - 1L, p)))
      acc <- acc + A[[j]] %*% Phi[[h - j]]
    Phi[[h]] <- acc
  }
  ## Structural responses: Theta_h = Phi_h %*% B.
  shock_names <- paste0("shock", seq_len(K))
  out <- vector("list", K)
  names(out) <- shock_names
  for (s in seq_len(K)) {
    M <- matrix(0, horizon, K, dimnames = list(NULL, var_names))
    for (h in seq_len(horizon))
      M[h, ] <- (Phi[[h]] %*% B)[, s]
    out[[s]] <- M
  }
  out
}


#' Structural impulse responses from a reduced-form VAR
#'
#' Recovers structural IRFs from a \code{\link{estimate_var}} fit using a
#' recursive (Cholesky) identification: the structural impact matrix is the
#' lower-triangular Cholesky factor \eqn{B} of the residual covariance, so
#' \eqn{u_t = B e_t} with mutually orthogonal unit-variance structural shocks
#' \eqn{e_t}.  Variable ordering therefore matters (first variable responds
#' only to its own shock on impact, etc.), exactly as in the standard
#' recursive SVAR.
#'
#' @param varfit A \code{var_fit} from \code{\link{estimate_var}}.
#' @param horizon Number of IRF periods to compute (horizon index 1 = impact).
#' @param identification Identification scheme; currently only
#'   \code{"cholesky"} (the default) is implemented.
#' @param shock Optional shock selector: \code{NULL} (default) returns all
#'   \code{K} structural shocks; an integer index or the variable name whose
#'   ordering position defines the shock returns just that shock's
#'   \code{horizon x K} matrix (the format \code{\link{match_irfs}} consumes
#'   directly).
#'
#' @return When \code{shock = NULL}, a named list (keyed \code{"shock1"},
#'   ..., \code{"shockK"}) of \code{horizon x K} IRF matrices (rows = horizon,
#'   columns = variables, named after the VAR variables).  When \code{shock}
#'   is given, the single \code{horizon x K} matrix for that shock.  The impact
#'   (horizon-1) responses equal the corresponding column of the Cholesky
#'   factor \eqn{B}.
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(400), 200, 2)
#' fit <- estimate_var(Y, p = 1)
#' irf <- var_irf(fit, horizon = 10)
#' irf$shock1[1, ]   # impact response to shock 1 = column 1 of chol(Sigma)
#'
#' @seealso \code{\link{estimate_var}}, \code{\link{var_irf_bootstrap}},
#'   \code{\link{match_irfs}}
#' @export
var_irf <- function(varfit, horizon = 20L,
                    identification = "cholesky", shock = NULL) {
  if (!inherits(varfit, "var_fit"))
    stop("var_irf: `varfit` must be a var_fit (from estimate_var()).")
  identification <- match.arg(identification, c("cholesky"))
  horizon <- as.integer(horizon)
  if (horizon < 1L) stop("var_irf: `horizon` must be >= 1.")

  ## Cholesky factor: lower triangular B with B B' = Sigma.
  B <- t(chol(varfit$Sigma))                       # chol() is upper; t() -> lower
  dimnames(B) <- list(varfit$var_names, varfit$var_names)

  irfs <- .var_irf_from_coef(varfit$A, B, horizon, varfit$var_names)

  if (is.null(shock)) return(irfs)

  ## Resolve a single-shock selector.
  idx <- if (is.character(shock)) {
    j <- match(shock, varfit$var_names)
    if (is.na(j)) stop(sprintf("var_irf: shock '%s' not a VAR variable.", shock))
    j
  } else {
    j <- as.integer(shock)
    if (j < 1L || j > varfit$K) stop("var_irf: `shock` index out of range.")
    j
  }
  irfs[[idx]]
}


#' Bootstrap confidence bands for VAR structural IRFs
#'
#' Residual-resampling (nonparametric) bootstrap of the structural IRFs from
#' \code{\link{var_irf}}.  For each replication the VAR residuals are resampled
#' with replacement, a bootstrap sample is regenerated recursively from the
#' fitted coefficients and the original initial conditions, the VAR is
#' re-estimated, and its Cholesky-identified IRFs are recomputed.  Percentile
#' bands are then read off the bootstrap distribution.  Fully deterministic
#' given \code{seed}.
#'
#' @param varfit A \code{var_fit} from \code{\link{estimate_var}}.
#' @param horizon Number of IRF periods.
#' @param Y The original \code{T x K} data used to fit \code{varfit} (needed to
#'   recover the \code{p} initial-condition rows for the recursive DGP).
#' @param n_boot Number of bootstrap replications.  Default \code{1000}.
#' @param ci Two-element vector of lower/upper percentiles for the band.
#'   Default \code{c(0.16, 0.84)} (the one-standard-deviation band).
#' @param shock Optional single-shock selector (see \code{\link{var_irf}}).
#'   \code{NULL} (default) bands every shock.
#' @param seed Integer RNG seed for reproducibility.  Default \code{1L}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{point}{The point-estimate IRF(s) (same shape as
#'       \code{\link{var_irf}}).}
#'     \item{lower, upper}{Lower/upper percentile bands, matching the shape of
#'       \code{point}.}
#'     \item{ci, n_boot}{The percentiles and replication count used.}
#'   }
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(400), 200, 2)
#' fit <- estimate_var(Y, p = 1)
#' bands <- var_irf_bootstrap(fit, horizon = 8, Y = Y, n_boot = 200, seed = 42)
#' bands$lower$shock1[1, ]   # lower band, impact, shock 1
#'
#' @seealso \code{\link{var_irf}}, \code{\link{match_irfs}}
#' @export
var_irf_bootstrap <- function(varfit, horizon = 20L, Y = NULL,
                              n_boot = 1000L, ci = c(0.16, 0.84),
                              shock = NULL, seed = 1L) {
  if (!inherits(varfit, "var_fit"))
    stop("var_irf_bootstrap: `varfit` must be a var_fit.")
  if (is.null(Y))
    stop("var_irf_bootstrap: supply the original data `Y` (for initial conditions).")
  if (length(ci) != 2L || any(ci < 0) || any(ci > 1) || ci[1] >= ci[2])
    stop("var_irf_bootstrap: `ci` must be c(lo, hi) with 0 <= lo < hi <= 1.")
  horizon <- as.integer(horizon)
  n_boot  <- as.integer(n_boot)

  Y  <- as.matrix(Y)
  K  <- varfit$K; p <- varfit$p; type <- varfit$type
  Tn <- nrow(Y)
  resid <- varfit$residuals
  neff  <- nrow(resid)

  ## Point estimate.
  point <- var_irf(varfit, horizon = horizon, shock = shock)

  ## Fitted-value reconstruction needs coef in the same layout estimate_var
  ## uses: regressors = [lag1(K), ..., lagp(K), det].  We rebuild each
  ## bootstrap path recursively from the original p initial rows.
  coef <- varfit$coef                              # (Kp+ndet) x K
  det_row <- function(t_idx) switch(type,
    none  = numeric(0),
    const = 1,
    trend = c(1, t_idx))

  one_rep <- function() {
    ## Resample residuals with replacement.
    idx <- sample.int(neff, neff, replace = TRUE)
    ub  <- resid[idx, , drop = FALSE]
    ## Regenerate Y* recursively from the original first p rows.
    Ys <- matrix(0, Tn, K)
    Ys[seq_len(p), ] <- Y[seq_len(p), , drop = FALSE]
    for (t in (p + 1L):Tn) {
      reg <- numeric(0)
      for (lag in seq_len(p)) reg <- c(reg, Ys[t - lag, ])
      reg <- c(reg, det_row(t - p))
      Ys[t, ] <- as.numeric(reg %*% coef) + ub[t - p, ]
    }
    fit_b <- tryCatch(estimate_var(Ys, p = p, type = type), error = function(e) NULL)
    if (is.null(fit_b)) return(NULL)
    fit_b$var_names <- varfit$var_names            # keep naming stable
    tryCatch(var_irf(fit_b, horizon = horizon, shock = shock),
             error = function(e) NULL)
  }

  ## Reproducible draws.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)

  reps <- vector("list", n_boot)
  kept <- 0L
  for (b in seq_len(n_boot)) {
    r <- one_rep()
    if (!is.null(r)) { kept <- kept + 1L; reps[[kept]] <- r }
  }
  reps <- reps[seq_len(kept)]
  if (kept == 0L) stop("var_irf_bootstrap: all bootstrap replications failed.")

  ## Quantiles, replicating `point`'s shape (single matrix or named list).
  q_of <- function(extract) {
    arr <- vapply(reps, extract, FUN.VALUE = extract(point))
    ## arr has dims (horizon, K, kept) for a matrix extractor.
    lo <- apply(arr, c(1, 2), stats::quantile, probs = ci[1], names = FALSE)
    hi <- apply(arr, c(1, 2), stats::quantile, probs = ci[2], names = FALSE)
    pt <- extract(point)
    dn <- dimnames(pt)
    dimnames(lo) <- dn; dimnames(hi) <- dn
    list(lower = lo, upper = hi)
  }

  if (is.null(shock)) {
    lower <- vector("list", length(point)); names(lower) <- names(point)
    upper <- lower
    for (nm in names(point)) {
      qq <- q_of(function(x) x[[nm]])
      lower[[nm]] <- qq$lower; upper[[nm]] <- qq$upper
    }
  } else {
    qq <- q_of(function(x) x)
    lower <- qq$lower; upper <- qq$upper
  }

  list(point = point, lower = lower, upper = upper,
       ci = ci, n_boot = kept)
}
