## R/lp-estimand.R
## --------------------------------------------------------------------------
## Population local-projection (LP) estimands of a solved dynhr model, and
## the gap between each LP specification and the model's true response.
##
## Reference: You, Z. (2026). How well are state-dependent local projections
##   capturing nonlinearities?  arXiv:2602.14455.
##
## The paper's QVAR(1,1) laboratory is "motivated by pruned second-order
## perturbation solutions"; here the DGP IS the AFVRR pruned order-2 solution
## of the model (order 1 is the special case ghxx = ghxu = ghuu = 0).
##
## Notation (all per UNIT of the observed shock e_i, deviations from det SS):
##   s     = x1_{t-1}, the first-order (pruned) state component at t-1
##   d     = Sigma_e[, i] / Sigma_e[i, i]  (E[e_t | e_it = 1]; = unit vector
##           e_i for a diagonal Sigma_e)
##   Psi_h = ghu d (h = 0),  ghx hx^{h-1} hu d (h >= 1)       -- linear IRF
##   M_h   = coefficient on e_it^2 of y_{t+h}: the deterministic pruned
##           impulse from s = 0 minus its linear part
##   Lam_h = n_endo x n_s loading of the shock-state interaction e_it * s:
##           Lam_0 = ghxu (d (x) I),
##           Lam_h = ghx B_{h-1} + 0.5 ghxx [kron(hx^h, r) + kron(r, hx^h)],
##           B_h   = hx  B_{h-1} + 0.5 hxx  [   same bracket              ],
##           B_0   = hxu (d (x) I),  r = hx^{h-1} hu d.
##
## Conditional average response (paper Prop. 1 / 7, generalised to a pruned
## DSGE with the contemporaneous ghuu term):
##   CAR_h(s, delta) = Psi_h delta + (Lam_h s) delta + M_h delta^2.
## (x2_{t-1} enters y_{t+h} linearly and never multiplies e_t, so it drops out
## of the CAR; only the first-order state component matters.)
##
## Population LP coefficients under Gaussian shocks (every E[odd power] = 0,
## s independent of e_t, E[s] = 0):
##   Linear  : beta_h = Psi_h                                     (Prop. 2/8)
##   AsymLP  : beta_h(+/-) = Psi_h +/- sigma_i m M_h,
##             m = sqrt(2/pi) / (1 - 2/pi)                          (Prop. 3)
##   LagLP   : theta_zu = Var(z)^{-1} Cov(z, s) Lam_h',
##             theta_u  = Psi_h - theta_zu' E[z]                    (Prop. 4)
##   Feas    : LagLP coefficients + M_h on e_it^2                   (Prop. 6)
##   Infeas  : (Psi_h, Lam_h, M_h) -- recovers the CAR exactly      (Prop. 9)
## Controls independent of e_it (lags of y, of the shocks, a constant) leave
## all of these unchanged, so they are the estimands of the usual lag-
## augmented regressions too.
##
## Conditional MSE given e_it = delta (Theorem 1, generalised with
## sigma2_s -> Lam Sigma_x Lam', sigma2_{s|y} -> Lam Sigma_{x|z} Lam'):
##   L_linear = delta^2 V_s + M^2 delta^4
##   L_asym   = delta^2 V_s + M^2 (delta^2 - sigma_i m |delta|)^2
##   L_lag    = delta^2 V_s|z + M^2 delta^4
##   L_feas   = delta^2 V_s|z
## and the unconditional distance (eq. 19) integrates these over
## e_it ~ N(0, sigma_i^2):  E[delta^2] = sigma_i^2, E[delta^4] = 3 sigma_i^4,
## E[(delta^2 - sigma_i m |delta|)^2] = nu_m sigma_i^4,
## nu_m = 3 - 4 m sqrt(2/pi) + m^2.
## --------------------------------------------------------------------------


#' Population local-projection estimands versus the model's true response
#'
#' For a solved dynhr model (order 1, or the pruned order-2 solution) and an
#' observed structural shock, computes in closed form what a population local
#' projection (LP) of each outcome on that shock recovers at horizons
#' \eqn{0, \ldots, H}, and compares it with the model's true conditional
#' average response (CAR).  This is the diagnostic of You (2026), who uses a
#' quadratic VAR -- the structure of a pruned second-order perturbation
#' solution -- to show which aspects of a nonlinear response each LP
#' specification can and cannot capture.
#'
#' The true response of outcome \eqn{y_j} to a shock of size \eqn{\delta} in
#' \eqn{e_i}, conditional on the first-order state \eqn{s = x^{(1)}_{t-1}}, is
#' \deqn{\mathrm{CAR}_h(s, \delta) = \Psi_h \delta + (\Lambda_h s)\,\delta +
#'   M_h \delta^2,}
#' a baseline (linear) effect, a state-dependent effect and a higher-order
#' effect.  At order 1, \eqn{\Lambda_h = 0} and \eqn{M_h = 0}.  The LP
#' specifications considered are those of the paper:
#' \itemize{
#'   \item \code{linear}: \eqn{y_{j,t+h} = \beta_h e_{it} + \pi' W_t + \epsilon}.
#'     Under Gaussian shocks its estimand is \eqn{\beta_h = \Psi_h}
#'     \emph{exactly}, at order 2 as well as order 1: the linear LP recovers
#'     the first-order IRF and misses both the state-dependent and the
#'     higher-order effect (the paper's headline result).
#'   \item \code{asym}: sign-interacted LP (separate intercept and slope for
#'     \eqn{e_{it} > 0} and \eqn{e_{it} \le 0});
#'     \eqn{\beta_h^{\pm} = \Psi_h \pm \sigma_i m M_h} with
#'     \eqn{m = \sqrt{2/\pi}/(1 - 2/\pi)}.
#'   \item \code{lag}: the shock interacted with lagged observable proxies
#'     \eqn{z_{t-1}} (default: the outcome's own lag);
#'     \eqn{\theta_{zu} = \mathrm{Var}(z)^{-1}\mathrm{Cov}(z, s)\Lambda_h'},
#'     \eqn{\theta_u = \Psi_h - \theta_{zu}' E[z]}.
#'   \item \code{feas}: the paper's proposed specification, \code{lag} plus a
#'     squared-shock regressor, whose coefficient is \eqn{M_h}.
#'   \item the infeasible specification that interacts the shock with the
#'     true state recovers the CAR exactly; its coefficients are
#'     \eqn{(\Psi_h, \Lambda_h, M_h)}, returned as \code{loading} and
#'     \code{q}.
#' }
#' Controls \eqn{W_t} independent of \eqn{e_{it}} (a constant, lags of the
#' outcomes and shocks) do not change any of these estimands.
#'
#' The approximation error of each specification is summarised by the
#' paper's conditional mean-squared error given \eqn{e_{it} = \delta}
#' (Theorem 1; \code{mse_*} columns) and by its unconditional distance
#' \eqn{D = (\sum_{h=0}^H E[(\mathrm{CAR}_h - \mathrm{IRF}_h)^2])^{1/2}}
#' over the joint distribution of state and shock (eq. 19;
#' \code{distance}).  The results require Gaussian shocks (a model declaring
#' skewed shocks gets a warning of class
#' \code{dynhr_warning_lp_gaussian}); for order 3 and
#' above the third-order terms enter the linear LP (the paper's Remark 2), so
#' those decision rules are rejected.
#'
#' \strong{Observed shock and correlated shocks.}  The LP regressor is the
#' observed shock \eqn{e_{it}} itself.  When \eqn{\Sigma_e} has off-diagonal
#' entries, the other shocks load on it through
#' \eqn{E[e_t | e_{it}] = d\, e_{it}}, \eqn{d = \Sigma_e[, i]/\Sigma_e[i, i]},
#' and every estimand above is the response to the impulse vector
#' \eqn{d\,\delta}.  For a diagonal \eqn{\Sigma_e}, \eqn{d} is the unit
#' vector; then \code{lp_linear} equals the \code{\link{compute_irfs}} IRF and
#' \code{car_mean} the \code{\link{compute_irfs_order2}} IRF, row for row.
#'
#' @param dr A \code{DecisionRules} (order 1) or \code{DecisionRules2}
#'   (order 2, pruned) object from \code{solve_perturbation()}.
#' @param model The parsed model (for the shock covariance).
#' @param shock Name of the observed shock (one of \code{dr$exo_names}).
#' @param outcome Character vector of endogenous variables; \code{NULL}
#'   (default) takes all.
#' @param horizon Maximum horizon \eqn{H}; horizons \eqn{0, \ldots, H} are
#'   returned (\eqn{h = 0} is the impact period, row 1 of
#'   \code{\link{compute_irfs}}).
#' @param shock_size Shock size \eqn{\delta} in standard deviations of the
#'   shock (default 1), the convention of \code{\link{compute_irfs}}.  Negative
#'   values give a negative shock, which matters for \code{asym} and for the
#'   higher-order effect.
#' @param proxy Optional character vector of endogenous variables whose
#'   \emph{lagged} values are the observable state proxies \eqn{z_{t-1}} of
#'   the \code{lag} and \code{feas} specifications, common to all outcomes.
#'   \code{NULL} (default) uses each outcome's own lag, as in the paper.
#' @param state Optional named numeric vector of the state variables
#'   (\code{dr$endo_names[dr$state_idx]}) at \eqn{t-1}, in LEVELS, read as the
#'   first-order state component.  When given, the true state-conditional
#'   response \code{car_state} is added.
#' @param proxy_value Optional named numeric vector of proxy values
#'   \eqn{z_{t-1}} in LEVELS (names covering the proxies used).  When given,
#'   the implied \code{lag} and \code{feas} IRFs at that value are added as
#'   \code{lp_lag_at} and \code{lp_feas_at}.
#' @param params Optional parameter vector; the shock covariance follows the
#'   rule of \code{\link{compute_irfs}} (\code{params}, else \code{dr$Sigma_e},
#'   else the model's calibration).
#' @return An object of class \code{"lp_estimand"}: a list with
#'   \describe{
#'     \item{table}{Data frame, one row per (outcome, horizon): \code{h},
#'       \code{outcome}; \code{lp_linear} (linear-LP IRF,
#'       \eqn{\Psi_h\delta}); \code{car_mean} (true response averaged over the
#'       state, \eqn{\Psi_h\delta + M_h\delta^2}; also the estimand of the LP
#'       with a squared-shock term); \code{higher_order}
#'       (\eqn{M_h\delta^2} = \code{car_mean - lp_linear}, the part the linear
#'       LP misses on average); \code{state_sd} (standard deviation across
#'       states of the state-dependent effect,
#'       \eqn{|\delta|(\Lambda_h\Sigma_x\Lambda_h')^{1/2}}); \code{lp_asym}
#'       (sign-interacted LP IRF for the sign of \eqn{\delta});
#'       \code{mse_linear}, \code{mse_asym}, \code{mse_lag}, \code{mse_feas}
#'       (conditional MSE given \eqn{e_{it} = \delta}); and, when requested,
#'       \code{car_state}, \code{lp_lag_at}, \code{lp_feas_at}.}
#'     \item{coefficients}{List of \code{(H+1) x n_outcome} matrices of
#'       population coefficients per unit of \eqn{e_{it}}: \code{beta_linear},
#'       \code{beta_pos}, \code{beta_neg}, \code{q} (squared-shock
#'       coefficient \eqn{M_h}), \code{theta_u}; plus per-outcome lists
#'       \code{theta_zu} (\code{(H+1) x n_proxy}) and \code{loading}
#'       (\code{(H+1) x n_s}, \eqn{\Lambda_h} on state deviations).}
#'     \item{distance}{Data frame of the unconditional distance \eqn{D} per
#'       outcome for \code{linear}, \code{asym}, \code{lag}, \code{feas}.}
#'     \item{shock, delta, shock_sd, direction, order, horizon, proxy}{The
#'       shock name, \eqn{\delta} in shock units, \eqn{\sigma_i}, the impulse
#'       direction \eqn{d}, the solution order, \eqn{H}, and the proxies used
#'       (a list per outcome).}
#'   }
#' @references You, Z. (2026). How well are state-dependent local projections
#'   capturing nonlinearities? arXiv:2602.14455.
#'
#'   Kolesar, M. and Plagborg-Moller, M. (2025). Dynamic causal effects in a
#'   nonlinear world: the good, the bad, and the ugly.
#' @seealso \code{\link{compute_irfs}}, \code{\link{compute_irfs_order2}},
#'   \code{\link{compute_girf}}, \code{\link{var_irf}}
#' @examples
#' \dontrun{
#' dr2 <- solve_perturbation(model, compiled, ss$values, model$param_values,
#'                           order = 2L)
#' lp <- lp_estimand(dr2, model, shock = "e_a", outcome = "y", horizon = 10)
#' lp$table[, c("h", "lp_linear", "car_mean", "higher_order")]
#' lp$distance
#' }
#' @export
lp_estimand <- function(dr, model, shock, outcome = NULL, horizon = 20L,
                        shock_size = 1, proxy = NULL, state = NULL,
                        proxy_value = NULL, params = NULL) {
  if (!inherits(dr, "DecisionRules"))
    .dynhr_abort("lp_estimand(): `dr` must be a DecisionRules object ",
                 "(solve_perturbation(), order 1 or 2).",
                 class = "dynhr_error_input")
  if (inherits(dr, "DecisionRules3"))
    .dynhr_abort("lp_estimand(): order >= 3 decision rules are not ",
                 "supported -- the closed forms need the order-2 (quadratic) ",
                 "structure; third-order terms enter the linear LP ",
                 "(You 2026, Remark 2).", class = "dynhr_error_input")

  endo <- dr$endo_names
  exo  <- dr$exo_names
  if (!is.character(shock) || length(shock) != 1L || !(shock %in% exo))
    .dynhr_abort("lp_estimand(): `shock` must be one of: ",
                 paste(exo, collapse = ", "), ".",
                 class = "dynhr_error_input")
  if (is.null(outcome)) outcome <- endo
  if (!is.character(outcome) || length(outcome) < 1L ||
      !all(outcome %in% endo))
    .dynhr_abort("lp_estimand(): `outcome` must name endogenous variables; ",
                 "unknown: ", paste(setdiff(outcome, endo), collapse = ", "),
                 class = "dynhr_error_input")
  if (!is.numeric(horizon) || length(horizon) != 1L || !is.finite(horizon) ||
      horizon < 0)
    .dynhr_abort("lp_estimand(): `horizon` must be a single integer >= 0.",
                 class = "dynhr_error_input")
  H <- as.integer(horizon)
  if (!is.numeric(shock_size) || length(shock_size) != 1L ||
      !is.finite(shock_size))
    .dynhr_abort("lp_estimand(): `shock_size` must be a finite number.",
                 class = "dynhr_error_input")
  if (!is.null(proxy) &&
      (!is.character(proxy) || length(proxy) < 1L || !all(proxy %in% endo)))
    .dynhr_abort("lp_estimand(): `proxy` must name endogenous variables.",
                 class = "dynhr_error_input")

  order2 <- inherits(dr, "DecisionRules2") && !is.null(dr$ghxx) &&
    nrow(dr$ghxx) == length(endo)

  ## ---- shock covariance, impulse direction, delta --------------------------
  Sigma_e <- .irf_shock_scale(dr, model, params)
  Sigma_e <- as.matrix(Sigma_e)
  i  <- match(shock, exo)
  s2 <- Sigma_e[i, i]
  if (!is.finite(s2) || s2 <= 0)
    .dynhr_abort("lp_estimand(): shock '", shock, "' has zero variance; an ",
                 "LP on it has no population estimand.",
                 class = "dynhr_error_input")
  ## The closed forms rest on E[odd power of e_t] = 0 (Gaussian shocks).  A
  ## skewed shock breaks them: the linear LP then also loads on the
  ## second-order effect (Kolesar & Plagborg-Moller 2025; You 2026, Sec. 2.2).
  alpha <- .get_shock_skewness(model, exo, params)
  if (any(alpha != 0))
    .dynhr_warn("lp_estimand(): the model declares skewed shocks (",
                paste(exo[alpha != 0], collapse = ", "), "); the population ",
                "LP formulas assume Gaussian shocks and are not exact here.",
                class = "dynhr_warning_lp_gaussian")
  sig   <- sqrt(s2)
  d     <- as.numeric(Sigma_e[, i]) / s2
  names(d) <- exo
  delta <- shock_size * sig
  m_c   <- sqrt(2 / pi) / (1 - 2 / pi)
  nu_m  <- 3 - 4 * m_c * sqrt(2 / pi) + m_c^2

  ## ---- decision-rule blocks -------------------------------------------------
  sidx <- dr$state_idx
  n_s  <- length(sidx)
  n_u  <- length(exo)
  n_e  <- length(endo)
  ghx  <- matrix(dr$ghx, n_e, n_s)
  ghu  <- matrix(dr$ghu, n_e, n_u)
  hx   <- ghx[sidx, , drop = FALSE]
  hu   <- ghu[sidx, , drop = FALSE]
  if (order2) {
    ghxx <- matrix(dr$ghxx, n_e, n_s * n_s)
    ghxu <- matrix(dr$ghxu, n_e, n_s * n_u)
    ghuu <- matrix(dr$ghuu, n_e, n_u * n_u)
  } else {
    ghxx <- matrix(0, n_e, n_s * n_s)
    ghxu <- matrix(0, n_e, n_s * n_u)
    ghuu <- matrix(0, n_e, n_u * n_u)
  }
  hxx <- ghxx[sidx, , drop = FALSE]
  hxu <- ghxu[sidx, , drop = FALSE]
  huu <- ghuu[sidx, , drop = FALSE]

  ## ---- horizon recursions: Psi_h d, M_h (d(x)d), Lam_h -----------------------
  ## ghxu / hxu columns are (state FAST, exo SLOW): the Kronecker vector is
  ## (e (x) x), so the bilinear map s -> ghxu (d (x) s) is ghxu (d (x) I_s).
  dI   <- kronecker(matrix(d, ncol = 1L), diag(n_s))            # n_u n_s x n_s
  Psi  <- matrix(0, H + 1L, n_e)
  Mq   <- matrix(0, H + 1L, n_e)
  Lam  <- array(0, dim = c(H + 1L, n_e, n_s))
  x1   <- as.numeric(hu %*% d)                                   # hx^0 hu d
  x2   <- 0.5 * as.numeric(huu %*% (d %x% d))
  B    <- hxu %*% dI                                             # n_s x n_s
  Psi[1L, ] <- as.numeric(ghu %*% d)
  Mq[1L, ]  <- 0.5 * as.numeric(ghuu %*% (d %x% d))
  Lam[1L, , ] <- ghxu %*% dI
  P <- diag(n_s)
  for (h in seq_len(H)) {
    P     <- hx %*% P                                            # hx^h
    r     <- matrix(x1, ncol = 1L)                               # hx^{h-1} hu d
    cross <- kronecker(P, r) + kronecker(r, P)                   # n_s^2 x n_s
    Psi[h + 1L, ] <- as.numeric(ghx %*% x1)
    Mq[h + 1L, ]  <- as.numeric(ghx %*% x2) +
      0.5 * as.numeric(ghxx %*% (x1 %x% x1))
    Lam[h + 1L, , ] <- ghx %*% B + 0.5 * ghxx %*% cross
    B  <- hx %*% B + 0.5 * hxx %*% cross
    x2 <- as.numeric(hx %*% x2) + 0.5 * as.numeric(hxx %*% (x1 %x% x1))
    x1 <- as.numeric(hx %*% x1)
  }

  ## ---- stationary moments: Var(x1), Cov(x1_{t-1}, y_{t-1}), Var(y), E[y] ---
  Sx <- if (n_s > 0L) solve_lyapunov(hx, hu %*% Sigma_e %*% t(hu)) else
    matrix(0, 0, 0)
  C_xy <- hx %*% Sx %*% t(ghx) + hu %*% Sigma_e %*% t(ghu)       # n_s x n_e
  if (order2) {
    st   <- .order2_stationary_moments(.order2_aug_system(dr, Sigma_e))
    Vy   <- st$var_cov
    mu_y <- as.numeric(st$mean)
  } else {
    Vy   <- ghx %*% Sx %*% t(ghx) + ghu %*% Sigma_e %*% t(ghu)
    mu_y <- as.numeric(dr$ys)
  }

  ## ---- optional evaluation state --------------------------------------------
  s_dev <- NULL
  if (!is.null(state)) {
    snames <- endo[sidx]
    if (!is.numeric(state) || is.null(names(state)) ||
        !all(snames %in% names(state)))
      .dynhr_abort("lp_estimand(): `state` must be a named numeric vector ",
                   "covering the state variables: ",
                   paste(snames, collapse = ", "), ".",
                   class = "dynhr_error_input")
    s_dev <- as.numeric(state[snames]) - as.numeric(dr$ys[sidx])
  }

  ## ---- per-outcome assembly -------------------------------------------------
  hh  <- 0:H
  mk  <- function() matrix(0, H + 1L, length(outcome),
                           dimnames = list(paste0("h", hh), outcome))
  beta_lin <- mk(); beta_pos <- mk(); beta_neg <- mk()
  q_mat <- mk(); theta_u <- mk()
  theta_zu <- list(); loading <- list(); proxies <- list()
  rows <- vector("list", length(outcome))
  dist <- data.frame(outcome = outcome, linear = 0, asym = 0, lag = 0,
                     feas = 0, stringsAsFactors = FALSE)
  sgn  <- if (delta > 0) 1 else if (delta < 0) -1 else 0

  for (k in seq_along(outcome)) {
    j   <- match(outcome[k], endo)
    prx <- if (is.null(proxy)) outcome[k] else proxy
    pidx <- match(prx, endo)
    Vz  <- Vy[pidx, pidx, drop = FALSE]
    ev  <- eigen((Vz + t(Vz)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 1e-12 * max(1, max(abs(ev))))
      .dynhr_abort("lp_estimand(): the proxy variance Var(",
                   paste(prx, collapse = ", "), ") is singular; choose ",
                   "proxies that move with the state.",
                   class = "dynhr_error_input")
    C_zx <- t(C_xy[, pidx, drop = FALSE])                         # n_z x n_s
    Vz_inv_Czx <- solve(Vz, C_zx)                                 # n_z x n_s
    Sx_z <- Sx - t(C_zx) %*% Vz_inv_Czx                           # Var(s | z)
    mu_z <- mu_y[pidx]

    lam <- matrix(Lam[, j, ], H + 1L, n_s)                        # (H+1) x n_s
    psi <- Psi[, j]
    q   <- Mq[, j]
    v_s  <- rowSums((lam %*% Sx) * lam)
    v_sz <- pmax(rowSums((lam %*% Sx_z) * lam), 0)
    tz  <- lam %*% t(Vz_inv_Czx)                                  # (H+1) x n_z
    tu  <- psi - as.numeric(tz %*% mu_z)
    colnames(tz)  <- prx
    colnames(lam) <- endo[sidx]
    rownames(tz) <- rownames(lam) <- paste0("h", hh)

    beta_lin[, k] <- psi
    beta_pos[, k] <- psi + sig * m_c * q
    beta_neg[, k] <- psi - sig * m_c * q
    q_mat[, k]    <- q
    theta_u[, k]  <- tu
    theta_zu[[outcome[k]]] <- tz
    loading[[outcome[k]]]  <- lam
    proxies[[outcome[k]]]  <- prx

    df <- data.frame(
      h            = hh,
      outcome      = outcome[k],
      lp_linear    = psi * delta,
      car_mean     = psi * delta + q * delta^2,
      higher_order = q * delta^2,
      state_sd     = abs(delta) * sqrt(pmax(v_s, 0)),
      lp_asym      = (psi + sgn * sig * m_c * q) * delta,
      mse_linear   = delta^2 * v_s + q^2 * delta^4,
      mse_asym     = delta^2 * v_s + q^2 * (delta^2 - sig * m_c * abs(delta))^2,
      mse_lag      = delta^2 * v_sz + q^2 * delta^4,
      mse_feas     = delta^2 * v_sz,
      stringsAsFactors = FALSE
    )
    if (!is.null(s_dev))
      df$car_state <- (psi + as.numeric(lam %*% s_dev)) * delta + q * delta^2
    if (!is.null(proxy_value)) {
      if (!is.numeric(proxy_value) || is.null(names(proxy_value)) ||
          !all(prx %in% names(proxy_value)))
        .dynhr_abort("lp_estimand(): `proxy_value` must be a named numeric ",
                     "vector covering the proxies: ",
                     paste(prx, collapse = ", "), ".",
                     class = "dynhr_error_input")
      zv <- as.numeric(proxy_value[prx])
      df$lp_lag_at  <- (tu + as.numeric(tz %*% zv)) * delta
      df$lp_feas_at <- df$lp_lag_at + q * delta^2
    }
    rows[[k]] <- df

    dist$linear[k] <- sqrt(sum(s2 * v_s + 3 * s2^2 * q^2))
    dist$asym[k]   <- sqrt(sum(s2 * v_s + nu_m * s2^2 * q^2))
    dist$lag[k]    <- sqrt(sum(s2 * v_sz + 3 * s2^2 * q^2))
    dist$feas[k]   <- sqrt(sum(s2 * v_sz))
  }

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  structure(
    list(table = tab,
         coefficients = list(beta_linear = beta_lin, beta_pos = beta_pos,
                             beta_neg = beta_neg, q = q_mat,
                             theta_u = theta_u, theta_zu = theta_zu,
                             loading = loading),
         distance = dist, shock = shock, delta = delta, shock_sd = sig,
         direction = d, order = if (order2) 2L else 1L, horizon = H,
         proxy = proxies),
    class = "lp_estimand")
}


#' Print an lp_estimand object
#'
#' @param x An \code{"lp_estimand"} object from \code{\link{lp_estimand}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.lp_estimand <- function(x, ...) {
  cat(sprintf(paste0("Population LP estimands (You 2026): shock '%s', ",
                     "delta = %.4g (%.4g sd), order %d, h = 0..%d\n"),
              x$shock, x$delta, x$delta / x$shock_sd, x$order, x$horizon))
  cat("Unconditional distance to the true response (eq. 19):\n")
  print(x$distance, row.names = FALSE, digits = 4)
  cat("Linear LP vs state-averaged true response (first rows):\n")
  print(utils::head(x$table[, c("outcome", "h", "lp_linear", "car_mean",
                                "higher_order", "state_sd")], 8L),
        row.names = FALSE, digits = 4)
  invisible(x)
}
