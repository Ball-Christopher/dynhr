## R/ms-filter.R
## --------------------------------------------------------------------------
## ms_kim_filter() -- Kim-Nelson (GPB(2)) filter for Markov-switching DSGE.
##
## Implements the Kim (1994) / Kim-Nelson (1999) filter for a state-space
## model where ONLY shock variances switch across regimes (structural
## parameters and TT, ZZ are common across regimes).
##
## Reference: Kim, C.-J. (1994), "Dynamic linear models with Markov-switching",
##   Journal of Econometrics 60(1-2), 1-22.
##   Kim, C.-J. & Nelson, C. R. (1999), "State-Space Models with Regime
##   Switching", MIT Press.
##
## dynhr state-space convention (lagged-state form):
##   s_t  = TT * s_{t-1} + RR * eps_t,   eps_t ~ N(0, I)   [ghu excl. Sigma_e]
##   y_t  = ZZ * s_{t-1} + DD * eps_t + d
##
## KEY INSIGHT: Because y_t depends on s_{t-1} (the PREVIOUS state), the
## innovation covariance uses the PRIOR covariance P_i, NOT the one-step-ahead
## predicted P_pred = TT*P_i*TT'+QQ. This is the fundamental difference from
## the standard (contemporaneous) state-space convention.
##
## Correct formulas for dynhr's lagged-state convention:
##   v_ij   = y_t - ZZ * b_i                      (innovation, uses prior b_i)
##   F_ij   = ZZ * P_i * ZZ' + HH_j               (uses PRIOR covariance P_i)
##   K_ij   = (TT * P_i * ZZ' + SS_j) * F_ij^{-1} (uses PRIOR P_i)
##   b_hat  = TT * b_i + K_ij * v_ij              (combined predict+update)
##   P_hat  = Joseph-form using PRIOR P_i
##
## Per-regime shock covariance:
##   Sigma_e^(j) = diag(scale_j) %*% Sigma_e %*% diag(scale_j)
##   QQ^(j)      = RR %*% Sigma_e^(j) %*% t(RR)
##   HH^(j)      = DD %*% Sigma_e^(j) %*% t(DD)
##   SS^(j)      = RR %*% Sigma_e^(j) %*% t(DD)
## --------------------------------------------------------------------------


#' Kim-Nelson filter for Markov-switching DSGE (shock-variance switching)
#'
#' Evaluates the log-likelihood of an MS-DSGE model where only shock variances
#' switch across regimes, using the Kim (1994) / GPB(2) filter.  The
#' structural parameters and decision rules (\code{ghx}, \code{ghu}) are
#' identical across regimes; only the per-regime shock covariance differs.
#'
#' @param data  Observation matrix (\code{n_obs x T}).  May contain \code{NA}s;
#'   missing observations at period \code{t} are handled by skipping the
#'   measurement update (propagating the state through transition only).
#' @param dr  Decision rule (output of \code{\link{solve_perturbation}}).
#' @param model  Compiled model object (output of \code{\link{compile_model}}).
#' @param params  Named numeric vector of parameter values.
#' @param obs_vars  Character vector of observed variable names.
#' @param ms_spec  An \code{\link{ms_dsge_spec}} object.
#' @param me_variance  Scalar variance of TRUE i.i.d. measurement error
#'   (default \code{0}).  See the \emph{Measurement error} section.
#' @param return_regime_probs  Logical; if \code{TRUE} return a
#'   \code{n_regimes x T} matrix of filtered regime probabilities
#'   \code{Pr[s_t = j | y_{1:t}]}.  Default \code{FALSE}.
#' @param return_state_path  Logical; if \code{TRUE} additionally return the
#'   per-regime collapsed filtered state moments, the per-path Durbin-Koopman
#'   backward blocks, and the state-space pieces the Kim smoother needs
#'   (\code{beta_filt}, \code{P_filt}, \code{dk_path}, \code{joint_filt},
#'   \code{cell_filt}, \code{TT}, \code{RR},
#'   \code{QQ_list}, \code{Sigma_e_list}, \code{P0_list}).  \code{joint_filt}
#'   is the \code{h x h x T} array of Hamilton-filter posterior JOINT regime
#'   probabilities \eqn{\Pr[s_{t-1} = i, s_t = j \mid y_{1:t}]};
#'   \code{cell_filt} (\code{"gpb3"} only) is the length-\code{T} list of
#'   UNCOLLAPSED per-cell posteriors, i.e. the per-TRIPLE filtered joint
#'   \eqn{\Pr[s_{t-2}, s_{t-1}, s_t \mid y_{1:t}]} the GPB(3) backward pass
#'   needs, and is \code{NULL} under \code{"gpb2"}.  Implies
#'   \code{return_regime_probs = TRUE}.  This is the interface
#'   \code{\link{ms_kim_smoother}} consumes; it exists so the smoother's
#'   forward pass IS the filter, never a second copy of the recursion.
#'   Default \code{FALSE}.
#' @param return_collapse_diag  Logical; if \code{TRUE} additionally return
#'   \code{collapse_diag} (\code{n_regimes x T}) and \code{collapse_max}, the
#'   GPB(2) COLLAPSE-QUALITY statistic.  See the "Collapse quality" section.
#'   Default \code{FALSE}.
#' @param lik_init  Character; state covariance initialisation: \code{"auto"},
#'   \code{"stationary"}, or \code{"kappa"}.  The exact-diffuse path is not
#'   yet supported for the MS filter; models with true unit roots should use
#'   \code{"kappa"}.  Each regime is initialised at its OWN unconditional
#'   state covariance (Kim & Nelson 1999); see \code{.ms_init_P0}.
#'
#' @param collapse  Character; the GPB collapse depth.  \code{"gpb2"}
#'   (default, and bit-identical to every release before this argument
#'   existed) is Kim's filter: it keeps \eqn{h} Gaussians, collapsing
#'   \eqn{h^2 \to h} each period on the LAST regime.  \code{"gpb3"} keeps
#'   \eqn{h^2} Gaussians indexed by the last TWO regimes
#'   \eqn{(s_{t-1}, s_t)}, collapsing \eqn{h^3 \to h^2}; it costs about
#'   \eqn{h} times as much per period and is exact whenever GPB(2)'s error
#'   is entirely due to the one-period collapse.  Use it when
#'   \code{collapse_max} from a \code{"gpb2"} run is O(1) nats or more (see
#'   the \emph{Collapse quality} section).  \code{return_state_path} is
#'   supported under \code{"gpb3"}, but the state path then carries the
#'   \eqn{h^2} pair-indexed components, which the Kim SMOOTHERS route to
#'   their own PAIR-INDEXED backward pass (F4-D): an adjoint per pair
#'   collapsed with the \eqn{h^3} smoothed joint
#'   \eqn{\Pr[s_{t-2}, s_{t-1}, s_t \mid y_{1:T}]}, whose weights come from
#'   \code{cell_filt}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{loglik}}{Total log-likelihood (scalar).}
#'     \item{\code{regime_probs}}{\code{n_regimes x T} matrix if
#'       \code{return_regime_probs = TRUE}, otherwise \code{NULL}.}
#'     \item{\code{n_obs}}{Number of observed variables.}
#'     \item{\code{n_T}}{Number of time periods.}
#'     \item{\code{beta_filt}, \code{P_filt}, \code{dk_path},
#'       \code{joint_filt}, \code{cell_filt}, \code{collapse}, \code{TT},
#'       \code{RR}, \code{QQ_list}, \code{Sigma_e_list}, \code{P0_list},
#'       \code{state_names}, \code{shock_names}, \code{lik_init}}{Present only
#'       when \code{return_state_path = TRUE} (\code{cell_filt} is
#'       \code{NULL} unless \code{collapse = "gpb3"}).}
#'     \item{\code{collapse_diag}, \code{collapse_max}}{Present only when
#'       \code{return_collapse_diag = TRUE}.}
#'   }
#'
#' @section Measurement error (TRUE noise, not a regulariser):
#' \code{me_variance} is the variance of a genuine i.i.d. observation noise
#' \eqn{u_t \sim N(0, \code{me\_variance} \cdot I)} appended to the
#' measurement equation,
#' \deqn{y_t = d_i + Z_i s_{t-1} + D_j \varepsilon_t + u_t ,}
#' so it enters BOTH the innovation covariance
#' \eqn{F_{ij} = Z P_i Z' + HH_j + \code{me\_variance} I} AND the Joseph
#' covariance update, which carries the extra \eqn{K_{ij} (\code{me\_variance}
#' I) K_{ij}'} term.  The resulting log-likelihood is the exact
#' joint-Gaussian likelihood of that model (verified against an
#' all-regime-path enumeration, and against
#' \code{kalman_filter(method = "univariate")} whenever the regime paths are
#' self-contained, e.g. \eqn{P = I}).  Every MS entry point
#' (\code{\link{ms_kim_filter}}, \code{\link{ms_kim_filter_struct}},
#' \code{\link{ms_kim_smoother}}, \code{\link{ms_kim_smoother_struct}}) uses
#' this same semantics, and it is the law the MS SBC data-generating
#' processes simulate from.
#'
#' This DIFFERS from the multivariate \code{\link{kalman_filter}} methods
#' (\code{"standard"}, \code{"dare"}, \code{"chandrasekhar"}), where
#' \code{me_variance} is documented as an \eqn{F}-only REGULARISER: it is
#' added to \eqn{F} but not carried through the covariance update, so it is
#' not an i.i.d. noise.  \code{kalman_filter(method = "univariate")} is the
#' true-ME law and is the single-regime oracle for the MS filters here.
#' The two conventions agree at \code{me_variance = 0} and differ by
#' \eqn{O(\code{me\_variance})} otherwise.
#'
#' @section Collapse quality (GPB(2) is an approximation, and it can fail):
#' Kim's filter keeps \eqn{h} Gaussians where the exact posterior is a mixture
#' over all \eqn{h^{T}} regime paths, collapsing \eqn{h^2 \to h} by
#' moment-matching each period.  \code{collapse_diag[t]} reports, in
#' LOG-LIKELIHOOD NATS, what the collapse at \eqn{t-1} cost at period
#' \eqn{t}: the period-\eqn{t} contribution recomputed against the \eqn{h^2}
#' UNCOLLAPSED components, minus the contribution the filter actually used,
#' \deqn{\log f_t^{*} - \log f_t .}
#' This is the one-step GPB(3)-vs-GPB(2) gap.  It is identically \code{0}
#' whenever the collapse is lossless --- in particular at every period under
#' \eqn{P = I}, and at \eqn{t = 1}, where nothing has been collapsed yet ---
#' and unlike a KL or a variance-share statistic it does not saturate when
#' the per-path covariances go singular, which is exactly the regime in which
#' the collapse fails.  \code{collapse_max} is \eqn{\max_t |collapse\_diag[t]|}.
#'
#' The collapse error is normally negligible, but it is NOT bounded.  On the
#' \code{rbc} two-structural-regime fixture with \code{n_obs = n_exo = 1} and
#' \code{me_variance = 0} the model is EXACTLY identified, so the per-path
#' covariances collapse toward singular.  Over 96 switching draws the filter
#' matches the exact all-path-enumeration log-likelihood to a median of
#' 3e-4 nats, but on one draw (seed 11) it misses by 48 nats --- while a full
#' GPB(3) reproduces the exact value on that same draw to \code{1e-15}.  That
#' contrast is the proof that the \eqn{h^2 \to h} collapse, and not the
#' recursion, is the entire source of the error.  Raising \code{me_variance}
#' off zero relieves the degeneracy: on that draw the GPB(2)-vs-GPB(3) gap
#' falls from 48 nats to 5.6 at \code{me_variance = 1e-2} and to 0.01 at
#' \code{1e-1}.  (Those three GPB(2)-vs-GPB(3) figures were measured before
#' F3-A made \code{me_variance} true measurement error; the qualitative
#' point --- adding observation noise relieves the degeneracy --- is
#' unchanged.  See the \emph{Measurement error} section.)
#'
#' Practical rule: \code{collapse_max} below ~0.1 nats means the collapse is
#' harmless on that sample.  If it runs to O(1) nats or more, treat the
#' log-likelihood as unreliable --- re-run with \code{collapse = "gpb3"}, add
#' measurement error, or observe fewer series than the model has shocks ---
#' rather than trusting the number.
#'
#' \code{collapse = "gpb3"} is that fix, at \eqn{h}x the cost: it keeps the
#' \eqn{h^2} components indexed by \eqn{(s_{t-1}, s_t)} and collapses
#' \eqn{h^3 \to h^2}.  On the seed-11 breakdown above it reproduces the exact
#' enumeration likelihood to \code{4e-16} where GPB(2) is 48 nats off, and
#' over the whole 96-case seed scan its worst error is \code{3e-14} against
#' GPB(2)'s 48 (its filtered regime probabilities likewise land within
#' \code{7e-16} of the exact posterior, where GPB(2) is 0.91 away on that
#' draw).  Under \code{"gpb3"} \code{collapse_diag} keeps its meaning ---
#' defer THIS filter's collapse by one period and rescore --- so it becomes
#' the GPB(4)-vs-GPB(3) gap, and it is \code{2e-14} on that same draw, which
#' is what says GPB(3) has CONVERGED here rather than merely improved.
#' Measured cost at \eqn{T = 300}, \eqn{h = 2}: 1.9x (structural filter),
#' 1.2x (reduced-form).  See \code{tests/testthat/test-ms-gpb3.R}.
#'
#' @seealso \code{\link{ms_kim_smoother}}, \code{\link{ms_irf}}
#' @export
ms_kim_filter <- function(data, dr, model, params, obs_vars, ms_spec,
                           me_variance = 0,
                           return_regime_probs = FALSE,
                           return_state_path = FALSE,
                           return_collapse_diag = FALSE,
                           lik_init = c("auto", "stationary", "kappa"),
                           collapse = c("gpb2", "gpb3")) {

  lik_init <- match.arg(lik_init)
  collapse <- match.arg(collapse)
  ## The Kim smoother needs the full per-regime forward moment path; it calls
  ## THIS function rather than re-deriving a second forward pass, so the
  ## smoother can never drift out of sync with the filter it smooths.
  if (isTRUE(return_state_path)) return_regime_probs <- TRUE

  ## ---- input checks -------------------------------------------------------
  if (!inherits(ms_spec, "ms_dsge_spec"))
    stop("ms_kim_filter: ms_spec must be an ms_dsge_spec object.", call. = FALSE)

  h <- ms_spec$n_regimes
  P <- ms_spec$transition    # h x h, rows sum to 1

  ## ---- extract state-space matrices (common across regimes) ---------------
  state_idx <- dr$state_idx
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  n_state   <- length(state_idx)
  n_exo     <- length(exo)
  n_obs     <- length(obs_vars)

  if (n_obs > n_exo)
    warning(sprintf("Stochastic singularity: %d obs but only %d shocks.", n_obs, n_exo))

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("ms_kim_filter: observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "), call. = FALSE)

  ## Validate that ms_spec shock-scale names match exo names (if names present)
  sc1 <- ms_spec$shock_scales[[1L]]
  if (!is.null(names(sc1))) {
    if (length(sc1) != n_exo ||
        !identical(sort(names(sc1)), sort(exo)))
      stop(sprintf(
        "ms_kim_filter: shock_scales names (%s) do not match model exo names (%s).",
        paste(names(sc1), collapse = ","), paste(exo, collapse = ",")),
        call. = FALSE)
    ## Reorder each regime's scale vector to match exo order
    ms_spec$shock_scales <- lapply(ms_spec$shock_scales, function(v) v[exo])
  } else if (length(sc1) != n_exo) {
    stop(sprintf(
      "ms_kim_filter: shock_scales length (%d) != n_exo (%d).",
      length(sc1), n_exo), call. = FALSE)
  }

  ghx <- dr$ghx; ghu <- dr$ghu
  TT  <- ghx[state_idx, , drop = FALSE]
  RR  <- ghu[state_idx, , drop = FALSE]
  ZZ  <- ghx[obs_idx,   , drop = FALSE]
  DD  <- ghu[obs_idx,   , drop = FALSE]
  d   <- dr$ys[obs_vars]
  tZZ <- t(ZZ)

  Sigma_e <- .get_shock_cov(model, exo, params)

  ## ---- per-regime covariance matrices (QQ, HH, SS, Sigma_e per regime) ---
  regime_covs <- .ms_build_regime_covs(RR, DD, Sigma_e, ms_spec$shock_scales)

  ## ---- prepare observations -----------------------------------------------
  if (is.null(dim(data))) data <- matrix(data, nrow = n_obs)
  if (nrow(data) != n_obs) data <- t(data)
  n_T <- ncol(data)
  Y_minus_d <- data - d

  ## ---- initialise state distributions (one per regime) --------------------
  ## PER-REGIME P_{0|0} (Kim & Nelson 1999, sec. 5.3): regime j starts at its
  ## OWN unconditional state covariance Lyapunov(TT, QQ_j), not at regime 1's.
  ##
  ## FIXED 2026-09-02 (E3-C).  The old code built ONE P0 from regime 1's QQ and
  ## handed the same matrix to every regime.  That is a bug, and the P = I
  ## oracle exposes it: with a degenerate transition matrix the Kim filter
  ## started in regime r must reproduce the fixed-regime-r Kalman filter
  ## exactly, but on the rbc 2-regime setup (scales 1 and 2, T = 150) it
  ## returned -74.0598 against the kalman_filter() oracle's +193.4512 -- a
  ## 267.51-nat gap that is CONSTANT in T (267.75 at T = 10, 267.5110 from
  ## T = 50 on), i.e. a pure initial-condition error, not a recursion error.
  QQ_list  <- lapply(regime_covs, `[[`, "QQ")
  P0_list  <- .ms_init_P0(TT, QQ_list, lik_init)
  lik_init <- attr(P0_list, "lik_init")

  ## ---- component storage ---------------------------------------------------
  ## The recursion carries a set of M Gaussian COMPONENTS.  Component m is the
  ## conditional law of s_{t-1} given y_{1:t-1} and a regime SUFFIX:
  ##
  ##   collapse = "gpb2"  M = h    component m = {s_{t-1} = m}
  ##   collapse = "gpb3"  M = h^2  component m = {s_{t-2} = i, s_{t-1} = j},
  ##                               m = i + (j - 1) h   (M = h at t = 1, where
  ##                               there is no s_{-1} to condition on)
  ##
  ##   Beta[[m]]      = E[s_{t-1} | suffix_m, y_{1:t-1}]   (FILTERED mean)
  ##   Pvar[[m]]      = Var[...]                           (FILTERED covariance)
  ##   comp_from[m]   = the regime of s_{t-1} in suffix m -- the FROM regime,
  ##                    the only part of the suffix the period-t update reads
  ##   comp_mass[m]   = Pr[suffix_m | y_{1:t-1}]
  ##
  ## GPB(2) and GPB(3) then differ in ONE thing: which coordinate of the
  ## post-update cell (m, k) is collapsed away.  GPB(2) keeps only the last
  ## regime k, so the h cells sharing a k merge; GPB(3) keeps the PAIR
  ## (s_{t-1}, s_t) = (comp_from[m], k), so only the cells sharing both merge.
  ## Everything else in the loop -- the Joseph update, the -Inf / singular-F /
  ## all-NA branches, the Durbin-Koopman blocks, the collapse diagnostic -- is
  ## identical, which is why `collapse = "gpb2"` is bit-identical to the
  ## pre-F3-C recursion: the same cell arithmetic is emitted in the same order.
  ##
  ## Initialised with s_0 = 0, P_0^{(j)} = P0_list[[j]].
  gpb3 <- identical(collapse, "gpb3")
  s0   <- numeric(n_state)
  Beta <- replicate(h, s0, simplify = FALSE)
  Pvar <- P0_list
  attr(Pvar, "lik_init") <- NULL

  comp_from <- seq_len(h)
  comp_mass <- ms_spec$pi0           # Pr[s_0 = j | y_{1:0}]
  ## Post-collapse component count and their FROM regimes.  Under gpb3 the new
  ## component (j, k) sits at index j + (k - 1) h, so the regime of s_t -- the
  ## FROM regime for the NEXT period -- is k.
  n_new    <- if (gpb3) h * h else h
  from_new <- if (gpb3) rep(seq_len(h), each = h) else seq_len(h)

  ## ---- constant for Gaussian likelihood -----------------------------------
  ll_const <- -0.5 * n_obs * log(2 * pi)
  me_diag  <- me_variance * diag(n_obs)

  loglik  <- 0
  ll_floor <- -1e300

  ## Optional: store filtered regime probs (h x T)
  if (return_regime_probs)
    reg_prob_out <- matrix(0, h, n_T)

  ## Optional: per-period collapse-quality statistic (length T, nats).
  ## Period 1 is 0 by definition -- nothing has been collapsed yet.  Under
  ## gpb3 this is the GPB(4)-vs-GPB(3) one-step gap, by exactly the same
  ## construction (defer THIS collapse by one period and rescore y_t).
  if (return_collapse_diag) {
    collapse_out <- numeric(n_T)
    cd_ZZ   <- rep(list(ZZ), h)
    cd_d    <- rep(list(numeric(n_obs)), h)   # Y is already offset-free here
    cd_HH   <- lapply(regime_covs, `[[`, "HH")
    prev_bh <- NULL
  }

  ## Optional: store the per-component COLLAPSED filtered moments
  ##   beta_path[, g, t] = beta_{t|t}^{(g)} ,  P_path[[t]][[g]] = P_{t|t}^{(g)}
  ## (consumed by ms_kim_smoother()).  Under gpb2 g indexes s_t (h of them);
  ## under gpb3 it indexes the PAIR (s_{t-1}, s_t) at j + (k - 1) h (h^2).
  ##
  ## Also store, per post-update CELL (component m -> to-regime k), the five
  ## blocks the Durbin-Koopman backward recursion needs.  They are all
  ## by-products of the measurement update that has already been done here, so
  ## the smoother never re-derives an innovation: it reuses this filter's own
  ## v / F^{-1} / K.
  ##   a  = Z_t' F^{-1} v          (n_state)          -> r_{t-1} forcing term
  ##   M  = Z_t' F^{-1} Z_t        (n_state x n_state) -> N_{t-1} forcing term
  ##   L  = TT - K Z_t             (n_state x n_state) -> the DK L_t
  ##   G  = (RR - K D_t)'          (n_shk x n_state)   -> shock adjoint loading
  ##   du = D_t' F^{-1} v          (n_shk)             -> shock forcing term
  ## with (a, M, du) = 0 and (L, G) = (TT, RR') on an all-missing or
  ## singular-F period, which is exactly kalman_smoother()'s no-observation
  ## branch (r_{t-1} = T' r_t, N_{t-1} = T' N_t T, eps = Q R' r).
  ## Cells are stored column-major in (m, k), i.e. at index m + (k - 1) M,
  ## which under gpb2 is the familiar i + (j - 1) h.
  ## joint_filt[i, j, t] = Pr[s_{t-1} = i, s_t = j | y_{1:t}] -- the Hamilton
  ## filter's own posterior JOINT, which the recursion already forms (and, up
  ## to F2-A, threw away after collapsing).  The joint regime smoother
  ## (ms_kim_smoother(regime_pass = "joint")) needs it: under dynhr's lag-1
  ## timing y_t loads s_{t-1} DIRECTLY, so this joint carries information about
  ## s_{t-1} that Kim's Pr[s_{t-1}=i | s_t=j, y_{1:t-1}] approximation
  ## discards.  Under gpb3 the SAME object is the collapsed group mass, which
  ## is the more accurate estimate of the same probability.
  ## cell_filt[[t]] is the FULL per-CELL posterior Pr[component_m at t-1,
  ## s_t = k | y_{1:t}] (M_t x h), i.e. under gpb3 the per-TRIPLE filtered
  ## joint Pr[s_{t-2} = i, s_{t-1} = j, s_t = k | y_{1:t}] laid out at
  ## [i + (j-1)h, k] (and the per-PAIR joint at t = 1, where there is no
  ## s_{-1}).  `joint_filt` is its collapse onto the surviving pair, so the
  ## gpb2 smoother needs nothing more; the GPB(3) BACKWARD pass does, because
  ## its triple pass conditions the collapsed-away s_{t-2} on y_{1:t} exactly
  ## as the gpb2 joint pass conditions s_{t-1} on y_{1:t} (F2-A, one level
  ## up).  Stored only under gpb3, so the gpb2 state path is untouched.
  if (return_state_path) {
    beta_path <- array(0, c(n_state, n_new, n_T))
    P_path    <- array(0, c(n_state, n_state, n_new, n_T))
    dk_path   <- vector("list", n_T)
    joint_filt <- array(0, c(h, h, n_T))
    cell_filt  <- if (gpb3) vector("list", n_T) else NULL
    .dk_null  <- list(a  = numeric(n_state),
                      M  = matrix(0, n_state, n_state),
                      L  = TT,
                      G  = t(RR),
                      du = numeric(n_exo))
  }

  ## ---- main Kim-Nelson loop -----------------------------------------------
  ##
  ## At each period t, for each (component m, to-regime k) CELL:
  ##   1. Predict: beta_pred = TT * b_m,  P_pred = TT * P_m * TT' + QQ_k
  ##      (one-step forecast of state, using regime-k shock variance)
  ##   2. Innovation (LAGGED-STATE): v = y_t - ZZ * b_m
  ##      (ZZ acts on s_{t-1}, so b_m is the correct base, NOT beta_pred)
  ##   3. Innovation covariance (LAGGED-STATE): F = ZZ * P_m * ZZ' + HH_k
  ##      (uses PRIOR P_m, not the predicted P_pred)
  ##   4. Kalman gain: K = (TT * P_m * ZZ' + SS_k) * F^{-1}
  ##   5. Updated state: b_hat = TT * b_m + K * v
  ##   6. Updated covariance: P_hat = Joseph-form
  ##   7. Likelihood weight: N(v; 0, F) * P[comp_from[m], k] * comp_mass[m]
  ##
  ## Collapsing (M*h cells -> n_new components):
  ##   b_g = sum_{cells c in g} wts[c] * b_hat_c
  ##   P_g = sum_c wts[c] * (P_hat_c + (b_hat_c - b_g)(b_hat_c - b_g)')
  ##   The outer-product cross-term is CRITICAL; omitting it is the most common
  ##   implementation bug (causes underestimated variance and -Inf loglik at t>1).
  ## The groups are CONTIGUOUS blocks of `n_coll` cells in the column-major
  ## (m, k) order, both under gpb2 (a whole k-column) and under gpb3 (the h
  ## cells of one (j, k) pair) -- which is what makes the two collapses the
  ## same code with a different block size.

  for (t in seq_len(n_T)) {
    y_t   <- Y_minus_d[, t]
    obs_ok <- is.finite(y_t)
    all_na <- !any(obs_ok)

    M      <- length(comp_mass)
    n_cell <- M * h
    n_coll <- n_cell %/% n_new     # cells merged per surviving component

    ## Storage for this period's cells.
    ## LOG joint density of each (component m -> to-regime k) cell: lp + log
    ## P[from_m, k] + log comp_mass[m]. Accumulated in log space (log-sum-exp
    ## below) so a tight regime whose Gaussian density underflows to 0 does NOT
    ## collapse f_y to 0 and spuriously return -Inf -- the failure mode when
    ## regimes have very different shock scales (e.g. a 0.5x volatility regime).
    log_lik_joint <- matrix(-Inf, M, h)   # [component m, to-regime k]
    beta_hat  <- matrix(0, n_state, n_cell)
    p_hat     <- vector("list", n_cell)
    ## DK blocks for this period, indexed dk_t[[m + (k - 1) * M]] (column-major,
    ## so the index matches log_lik_joint[m, k]).
    if (return_state_path) dk_t <- rep(list(.dk_null), n_cell)

    for (m in seq_len(M)) {
      i   <- comp_from[m]           # regime of s_{t-1} under this component
      b_i <- Beta[[m]]
      P_i <- Pvar[[m]]   # P_{t-1|t-1} for component m -- PRIOR covariance
      ## k-invariant blocks (depend only on the component) hoisted out of the
      ## k-loop: ZPZt_i is the state term of F; TPZt_i is the gain numerator.
      ## Only HH_k / SS_k / QQ_k carry the TO-regime k dependence.
      ZPZt_i <- ZZ %*% P_i %*% tZZ
      TPZt_i <- TT %*% P_i %*% tZZ

      for (j in seq_len(h)) {
        cell  <- m + (j - 1L) * M
        cov_j <- regime_covs[[j]]
        QQ_j  <- cov_j$QQ
        HH_j  <- cov_j$HH
        SS_j  <- cov_j$SS

        ## -- Prediction step ------------------------------------------------
        ## beta_pred = TT * b_i (one-step-ahead state forecast)
        ## P_pred    = TT * P_i * TT' + QQ_j (one-step-ahead cov forecast)
        ## P_pred is used for the UPDATED covariance, not for F_ij.
        beta_pred <- drop(TT %*% b_i)

        if (all_na) {
          ## No observations: skip measurement update; propagate state only.
          ## P_pred is referenced only here and on the singular-F fallback, so it
          ## is formed lazily -- the common (observed, non-singular) path, which
          ## uses the Joseph form on P_i directly, never computes it.
          P_pred <- tcrossprod(TT %*% P_i, TT) + QQ_j
          P_pred <- (P_pred + t(P_pred)) * 0.5
          beta_hat[, cell]    <- beta_pred
          p_hat[[cell]]       <- P_pred
          log_lik_joint[m, j] <- log(P[i, j] * comp_mass[m])  # density 1
        } else {
          ## -- Innovation (LAGGED-STATE convention) --------------------------
          ## v = y_t - ZZ * b_i   (b_i = s_{t-1|t-1}, not the predicted state)
          v_full <- y_t - as.numeric(ZZ %*% b_i)

          ## Handle partial NAs: zero out missing obs entries
          v_obs   <- v_full
          if (any(!obs_ok)) v_obs[!obs_ok] <- 0

          ## -- Innovation covariance F_ij (LAGGED-STATE) --------------------
          ## F_ij = ZZ * P_i * ZZ' + HH_j   (PRIOR P_i, not P_pred!)
          ## This is correct because y_t depends on s_{t-1}, so F is the
          ## variance of y_t - ZZ*s_{t-1|t-1} which uses P_{t-1|t-1} = P_i.
          F_ij <- ZPZt_i + HH_j + me_diag
          F_ij <- (F_ij + t(F_ij)) * 0.5

          ## Handle partial missing: project onto observed block
          if (any(!obs_ok)) {
            F_ij_obs <- F_ij[obs_ok, obs_ok, drop = FALSE]
            v_use    <- v_obs[obs_ok]
            n_obs_t  <- sum(obs_ok)
          } else {
            F_ij_obs <- F_ij
            v_use    <- v_obs
            n_obs_t  <- n_obs
          }

          Fc_ij <- tryCatch(chol(F_ij_obs), error = function(e) NULL)
          if (is.null(Fc_ij)) {
            ## Singular F for this cell: zero likelihood contribution; keep prior
            ## (lazy P_pred -- see the all_na branch above).
            P_pred <- tcrossprod(TT %*% P_i, TT) + QQ_j
            P_pred <- (P_pred + t(P_pred)) * 0.5
            log_lik_joint[m, j] <- -Inf
            beta_hat[, cell]    <- beta_pred
            p_hat[[cell]]       <- P_pred
            next
          }

          Fi_ij     <- chol2inv(Fc_ij)
          log_det_F <- 2 * sum(log(diag(Fc_ij)))
          ll_const_t <- -0.5 * n_obs_t * log(2 * pi)
          quad       <- drop(crossprod(v_use, Fi_ij %*% v_use))
          lp_ij      <- ll_const_t - 0.5 * (log_det_F + quad)

          log_lik_joint[m, j] <- lp_ij + log(P[i, j] * comp_mass[m])

          ## -- Kalman gain (LAGGED-STATE) ------------------------------------
          ## K_ij = (TT * P_i * ZZ' + SS_j) * F_ij^{-1}   (uses PRIOR P_i)
          if (any(!obs_ok)) {
            ZZ_obs  <- ZZ[obs_ok, , drop = FALSE]
            DD_obs  <- DD[obs_ok, , drop = FALSE]
            SS_j_obs <- RR %*% cov_j$Sigma_e %*% t(DD_obs)
            ## TT %*% P_i %*% t(ZZ_obs) is the observed-column subset of TPZt_i.
            K_ij    <- (TPZt_i[, obs_ok, drop = FALSE] + SS_j_obs) %*% Fi_ij
            b_upd   <- beta_pred + drop(K_ij %*% v_use)
            IKZ     <- TT - K_ij %*% ZZ_obs
            RmKD    <- RR - K_ij %*% DD_obs
          } else {
            K_ij <- (TPZt_i + SS_j) %*% Fi_ij
            b_upd <- beta_pred + drop(K_ij %*% v_use)
            IKZ   <- TT - K_ij %*% ZZ
            RmKD  <- RR - K_ij %*% DD
          }

          ## -- Updated covariance (Joseph-form) ---------------------------
          ## P_hat = (TT - K*ZZ) P_i (TT - K*ZZ)' + (RR - K*DD) Se (RR - K*DD)'
          ##         + K (me I) K'
          ## The last term is the measurement error propagated through the
          ## gain.  The one-step state error is
          ##   s_t - b_hat = (TT - K ZZ)(s_{t-1} - b_i) + (RR - K DD) eps - K u
          ## with u ~ N(0, me I) the TRUE i.i.d. observation noise, so leaving
          ## K me K' out makes me_variance an F-only regulariser instead of
          ## measurement error (that is the multivariate kalman_filter()
          ## convention, deliberately NOT this one -- see the roxygen).
          P_upd <- tcrossprod(IKZ %*% P_i, IKZ) +
                   tcrossprod(RmKD %*% cov_j$Sigma_e, RmKD)
          ## Guarded so me_variance = 0 stays bit-identical to the old code.
          if (me_variance > 0) P_upd <- P_upd + me_variance * tcrossprod(K_ij)
          P_upd <- (P_upd + t(P_upd)) * 0.5

          beta_hat[, cell]  <- b_upd
          p_hat[[cell]]     <- P_upd

          if (return_state_path) {
            ## ZZ_use / DD_use are the OBSERVED rows only; on a fully-observed
            ## period they are ZZ / DD themselves (no copy taken).
            ZZ_use <- if (any(!obs_ok)) ZZ[obs_ok, , drop = FALSE] else ZZ
            DD_use <- if (any(!obs_ok)) DD[obs_ok, , drop = FALSE] else DD
            Fv     <- Fi_ij %*% v_use
            dk_t[[cell]] <- list(
              a  = as.numeric(crossprod(ZZ_use, Fv)),
              M  = crossprod(ZZ_use, Fi_ij %*% ZZ_use),
              L  = IKZ,
              G  = t(RmKD),
              du = as.numeric(crossprod(DD_use, Fv))
            )
          }
        }
      }
    }

    ## -- Loglik contribution for period t (log-sum-exp, underflow-safe) ------
    mx <- max(log_lik_joint)
    if (!is.finite(mx)) {            # every path has -Inf log-density
      loglik <- -Inf
      break
    }
    log_f_y <- .logsumexp(log_lik_joint)
    loglik  <- loglik + log_f_y
    if (loglik < ll_floor) { loglik <- -Inf; break }

    ## -- Collapse-quality statistic: what the t-1 collapse cost THIS period --
    if (return_collapse_diag && !is.null(prev_bh) && !all_na) {
      lfd <- .ms_defer_logf(y_t, obs_ok, me_diag, prev_bh, prev_ph, prev_pj,
                            prev_from, P, cd_ZZ, cd_d, cd_HH)
      collapse_out[t] <- if (is.na(lfd)) NA_real_ else lfd - log_f_y
    }

    ## -- Hamilton filter: posterior cell probabilities -----------------------
    prob_joint <- exp(log_lik_joint - log_f_y)   # M x h, sums to 1

    ## Group mass Pr[surviving suffix g | y_{1:t}]: sum the cells of each
    ## contiguous block.  Under gpb2 a block IS a column of prob_joint, so this
    ## is the same colSums() the pre-F3-C code ran.
    group_mass_raw <- colSums(array(prob_joint, c(n_coll, n_new)))

    ## Guard: floor at 1e-300 to avoid division by zero in collapsing.  The RAW
    ## (unfloored) vector is kept: a group with EXACTLY zero posterior mass
    ## must be detected on the raw value, not on the floor (see the collapsing
    ## loop below).
    group_mass <- pmax(group_mass_raw, 1e-300)
    group_mass <- group_mass / sum(group_mass)

    ## -- Collapse: n_cell -> n_new ------------------------------------------
    ## Collapsed mean (Kim 1994, eq. 4, per group):
    ##   b_g = sum_{c in g} Pr[c | g, y_{1:t}] * b_hat_c
    ##
    ## Collapsed covariance with CROSS-TERM (Kim 1994, eq. 5):
    ##   P_g = sum_c wts[c] * (P_hat_c + (b_hat_c - b_g)(b_hat_c - b_g)')
    ##
    ## The outer-product term accounts for the variance INCREASE from merging
    ## several paths into one weighted mean.  Omitting it systematically
    ## understates the covariance, causing F to shrink below the true
    ## innovation variance in subsequent periods, which makes the filter assign
    ## overly high likelihoods to later innovations and eventually diverge to
    ## -Inf.
    Beta_new <- vector("list", n_new)
    Pvar_new <- vector("list", n_new)

    ## Fallback prior component for a group with EXACTLY zero mass: under gpb2
    ## the same-regime prior component (what the P = I / degenerate-pi0 oracle
    ## pins), under gpb3 the first prior component feeding the group, which
    ## carries the same FROM regime.
    fb <- if (gpb3) ((seq_len(n_new) - 1L) * n_coll) %% M + 1L else seq_len(h)

    for (g in seq_len(n_new)) {
      if (group_mass_raw[g] <= 0) {
        ## Group g is UNREACHABLE this period (exactly zero posterior mass --
        ## e.g. a degenerate pi0 under P = I).  Carry the prior moments forward.
        ## Testing the FLOORED probability instead would fall through here
        ## (1e-300 is not < 1e-300) and collapse with weights 0 / 1e-300 = 0,
        ## setting b_g = 0 and P_g = 0.  A zero covariance shrinks the dead
        ## group's F to HH alone, injecting a spurious high-density path back
        ## into the mixture at t+1 -- which is exactly what broke the
        ## P = I / degenerate-pi0 oracle.
        Beta_new[[g]] <- Beta[[fb[g]]]
        Pvar_new[[g]] <- Pvar[[fb[g]]]
        next
      }

      cells <- (g - 1L) * n_coll + seq_len(n_coll)
      wts   <- prob_joint[cells] / group_mass[g]   # sums to 1

      ## Collapsed mean
      b_j <- numeric(n_state)
      for (u in seq_len(n_coll)) b_j <- b_j + wts[u] * beta_hat[, cells[u]]
      Beta_new[[g]] <- b_j

      ## Collapsed covariance WITH CROSS-TERM
      P_j <- matrix(0, n_state, n_state)
      for (u in seq_len(n_coll)) {
        diff_ij <- beta_hat[, cells[u]] - b_j
        ## Critical: p_hat[[cell]] + outer-product cross-term
        P_j <- P_j + wts[u] * (p_hat[[cells[u]]] + tcrossprod(diff_ij))
      }
      Pvar_new[[g]] <- (P_j + t(P_j)) * 0.5

    }

    Beta      <- Beta_new
    Pvar      <- Pvar_new
    comp_mass <- group_mass
    comp_from <- from_new

    ## Filtered regime marginal Pr[s_t = k | y_{1:t}]: under gpb2 the group
    ## mass already IS that marginal; under gpb3 sum the (j, k) pairs over j.
    regime_prob <- if (gpb3) colSums(matrix(group_mass, h, h)) else group_mass

    ## Carry this period's PRE-collapse cells for the next period's
    ## collapse-quality statistic.  `ord` re-orders the cells so that gpb2
    ## reproduces the pre-F3-C (k outer, i inner) enumeration exactly, which
    ## keeps the diagnostic bit-identical as well as the loglik.
    if (return_collapse_diag) {
      ord <- as.integer(t(matrix(seq_len(n_cell), M, h)))
      prev_bh   <- beta_hat[, ord, drop = FALSE]
      prev_ph   <- p_hat[ord]
      prev_pj   <- as.numeric(prob_joint)[ord]
      prev_from <- rep(seq_len(h), each = M)[ord]
    }

    if (return_regime_probs)
      reg_prob_out[, t] <- regime_prob

    if (return_state_path) {
      for (g in seq_len(n_new)) {
        beta_path[, g, t]  <- Beta[[g]]
        P_path[, , g, t]   <- Pvar[[g]]
      }
      dk_path[[t]] <- dk_t
      joint_filt[, , t] <- if (gpb3) matrix(group_mass_raw, h, h) else prob_joint
      if (gpb3) cell_filt[[t]] <- prob_joint
    }
  }

  out <- list(
    loglik       = loglik,
    regime_probs = if (return_regime_probs) reg_prob_out else NULL,
    n_obs        = n_obs,
    n_T          = n_T
  )

  if (return_collapse_diag) {
    out$collapse_diag <- collapse_out
    out$collapse_max  <- suppressWarnings(max(abs(collapse_out), na.rm = TRUE))
    if (!is.finite(out$collapse_max)) out$collapse_max <- NA_real_
  }

  if (return_state_path) {
    ## Everything ms_kim_smoother() needs for the Kim (1994) backward pass:
    ## the per-regime filtered moments, the (common) transition matrix, the
    ## per-regime state-noise covariances, and the resolved initialisation.
    out$beta_filt   <- beta_path
    out$P_filt      <- P_path
    out$dk_path     <- dk_path
    out$joint_filt  <- joint_filt
    out$cell_filt   <- cell_filt
    out$collapse    <- collapse
    out$TT          <- TT
    out$RR          <- RR
    out$QQ_list     <- QQ_list
    out$Sigma_e_list <- lapply(regime_covs, `[[`, "Sigma_e")
    out$P0_list     <- P0_list
    out$state_names <- endo[state_idx]
    out$shock_names <- exo
    out$lik_init    <- lik_init
  }
  out
}


## ============================================================================
## Structural MS Kim-Nelson filter (regime-specific TT/ZZ/RR/DD)
## ============================================================================

#' Kim-Nelson filter for structural MS-DSGE (regime-specific decision rules)
#'
#' Evaluates the log-likelihood of a structural MS-DSGE model where BOTH
#' structural parameters AND shock variances can differ across regimes.
#' Each regime has its own state-space matrices \code{TT_s, RR_s, ZZ_s, DD_s}
#' extracted from an \code{MsDecisionRules} object (output of
#' \code{\link{solve_ms_perturbation}}).
#'
#' The shock covariance \eqn{\Sigma_e} is per-regime from the decision rules;
#' for shared-\eqn{\Sigma_e} with only structural params switching, supply the
#' same \code{Sigma_e} for all regimes (the common case).
#'
#' @param data          Observation matrix (n_obs x T). May contain \code{NA}s.
#' @param ms_dr      An \code{MsDecisionRules} object from
#'                   \code{\link{solve_ms_perturbation}}.
#' @param model      Compiled model object.
#' @param params     Named numeric parameter vector (used only to compute
#'                   \eqn{\Sigma_e} if \code{Sigma_e_by_regime} is \code{NULL}).
#' @param obs_vars   Character vector of observed variable names.
#' @param Sigma_e_by_regime  Optional list of length h giving per-regime
#'                   shock covariance matrices. If \code{NULL}, a common
#'                   \eqn{\Sigma_e} from \code{params} is used for all regimes.
#' @param me_variance  Scalar variance of TRUE i.i.d. measurement error
#'                   \eqn{u_t \sim N(0, \code{me\_variance} \cdot I)} added
#'                   to the measurement equation (default 0).  It enters both
#'                   \eqn{F_{ij}} and the Joseph covariance update (the
#'                   \eqn{K (\code{me\_variance} I) K'} term), so the
#'                   log-likelihood is the exact joint-Gaussian likelihood of
#'                   that model.  This is the same semantics as
#'                   \code{\link{ms_kim_filter}} --- see its
#'                   \emph{Measurement error} section --- and it DIFFERS from
#'                   the multivariate \code{\link{kalman_filter}} methods,
#'                   where \code{me_variance} is an \eqn{F}-only regulariser;
#'                   \code{kalman_filter(method = "univariate")} is the
#'                   single-regime true-ME oracle.
#' @param return_regime_probs  Logical; return \code{n_regimes x T} filtered
#'                   regime probability matrix (default FALSE).
#' @param return_state_path  Logical; if \code{TRUE} additionally return the
#'                   per-regime collapsed filtered moments, the per-path
#'                   Durbin-Koopman backward blocks and the state-space pieces
#'                   \code{\link{ms_kim_smoother_struct}} consumes
#'                   (\code{beta_filt}, \code{P_filt}, \code{dk_path},
#'                   \code{joint_filt}, \code{cell_filt},
#'                   \code{TT_list}, \code{RR_list}, \code{QQ_list},
#'                   \code{Sigma_e_list}, \code{P0_list}).  \code{joint_filt}
#'                   is the \code{h x h x T} array of Hamilton-filter posterior
#'                   JOINT regime probabilities
#'                   \eqn{\Pr[s_{t-1} = i, s_t = j \mid y_{1:t}]};
#'                   \code{cell_filt} (\code{"gpb3"} only) is the length-T
#'                   list of UNCOLLAPSED per-cell posteriors, the per-TRIPLE
#'                   filtered joint the GPB(3) backward pass needs.  Implies
#'                   \code{return_regime_probs = TRUE}.  The stored blocks are
#'                   by-products of the measurement update that has already
#'                   happened, so the flag is numerically inert: the
#'                   log-likelihood and the regime probabilities are
#'                   BIT-IDENTICAL with and without it (pinned in
#'                   \code{test-ms-smoother-struct.R}).  Default \code{FALSE}.
#' @param return_collapse_diag  Logical; if \code{TRUE} additionally return
#'                   \code{collapse_diag} (\code{n_regimes x T}) and
#'                   \code{collapse_max}, the GPB(2) COLLAPSE-QUALITY
#'                   statistic.  See the "Collapse quality" section of
#'                   \code{\link{ms_kim_filter}}; this filter is where it
#'                   matters most, because with regime-specific
#'                   \code{ZZ_s} / \code{d_s} the path means separate much
#'                   faster than in the shock-variance-only model.  Default
#'                   \code{FALSE}.
#' @param lik_init   State covariance initialisation: \code{"auto"},
#'                   \code{"stationary"}, or \code{"kappa"}.
#'
#' @param collapse  Character; the GPB collapse depth.  \code{"gpb2"}
#'   (default, and bit-identical to every release before this argument
#'   existed) is Kim's filter: it keeps \eqn{h} Gaussians, collapsing
#'   \eqn{h^2 \to h} each period on the LAST regime.  \code{"gpb3"} keeps
#'   \eqn{h^2} Gaussians indexed by the last TWO regimes
#'   \eqn{(s_{t-1}, s_t)}, collapsing \eqn{h^3 \to h^2}; it costs about
#'   \eqn{h} times as much per period and is exact whenever GPB(2)'s error
#'   is entirely due to the one-period collapse.  Use it when
#'   \code{collapse_max} from a \code{"gpb2"} run is O(1) nats or more (see
#'   the \emph{Collapse quality} section).  \code{return_state_path} is
#'   supported under \code{"gpb3"}, but the state path then carries the
#'   \eqn{h^2} pair-indexed components, which the Kim SMOOTHERS route to
#'   their own PAIR-INDEXED backward pass (F4-D): an adjoint per pair
#'   collapsed with the \eqn{h^3} smoothed joint
#'   \eqn{\Pr[s_{t-2}, s_{t-1}, s_t \mid y_{1:T}]}, whose weights come from
#'   \code{cell_filt}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{loglik}}{Total log-likelihood (scalar).}
#'     \item{\code{regime_probs}}{\code{n_regimes x T} matrix if
#'       \code{return_regime_probs = TRUE}, otherwise \code{NULL}.}
#'     \item{\code{n_obs}}{Number of observed variables.}
#'     \item{\code{n_T}}{Number of time periods.}
#'     \item{\code{beta_filt}, \code{P_filt}, \code{dk_path},
#'       \code{joint_filt}, \code{cell_filt}, \code{collapse},
#'       \code{TT_list},
#'       \code{RR_list}, \code{QQ_list}, \code{Sigma_e_list}, \code{P0_list},
#'       \code{d_list}, \code{state_names}, \code{shock_names},
#'       \code{lik_init}}{Present only when
#'       \code{return_state_path = TRUE} (\code{cell_filt} is \code{NULL}
#'       unless \code{collapse = "gpb3"}).}
#'     \item{\code{collapse_diag}, \code{collapse_max}}{Present only when
#'       \code{return_collapse_diag = TRUE}.}
#'   }
#' @seealso \code{\link{ms_kim_smoother_struct}}, \code{\link{ms_kim_filter}}
#' @export
ms_kim_filter_struct <- function(data, ms_dr, model, params, obs_vars,
                                  Sigma_e_by_regime = NULL,
                                  me_variance = 0,
                                  return_regime_probs = FALSE,
                                  return_state_path = FALSE,
                                  return_collapse_diag = FALSE,
                                  lik_init = c("auto", "stationary", "kappa"),
                                  collapse = c("gpb2", "gpb3")) {

  lik_init <- match.arg(lik_init)
  collapse <- match.arg(collapse)
  if (isTRUE(return_state_path)) return_regime_probs <- TRUE

  if (!inherits(ms_dr, "MsDecisionRules"))
    stop("ms_kim_filter_struct: ms_dr must be an MsDecisionRules object.", call. = FALSE)

  h <- length(ms_dr$dr)
  P <- ms_dr$P

  ## ---- extract per-regime state-space matrices ----------------------------
  ## Use regime 1 to get shared structural info (state_idx, endo, exo).
  dr1       <- ms_dr$dr[[1L]]
  state_idx <- dr1$state_idx
  endo      <- dr1$endo_names
  exo       <- dr1$exo_names
  n_state   <- length(state_idx)
  n_exo     <- length(exo)
  n_obs     <- length(obs_vars)

  if (n_obs > n_exo)
    warning(sprintf("Stochastic singularity: %d obs but only %d shocks.", n_obs, n_exo))

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("ms_kim_filter_struct: observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "), call. = FALSE)

  ## Per-regime TT_s, RR_s, ZZ_s, DD_s, d_s (observable SS mean)
  TT_list <- vector("list", h)
  RR_list <- vector("list", h)
  ZZ_list <- vector("list", h)
  DD_list <- vector("list", h)
  d_list  <- vector("list", h)

  for (s in seq_len(h)) {
    dr_s       <- ms_dr$dr[[s]]
    TT_list[[s]] <- dr_s$ghx[state_idx, , drop = FALSE]
    RR_list[[s]] <- dr_s$ghu[state_idx, , drop = FALSE]
    ZZ_list[[s]] <- dr_s$ghx[obs_idx,   , drop = FALSE]
    DD_list[[s]] <- dr_s$ghu[obs_idx,   , drop = FALSE]
    d_list[[s]]  <- dr_s$ys[obs_vars]
  }

  ## ---- per-regime shock covariance ----------------------------------------
  if (is.null(Sigma_e_by_regime)) {
    Sigma_e_common <- .get_shock_cov(model, exo, params)
    Sigma_e_by_regime <- replicate(h, Sigma_e_common, simplify = FALSE)
  }

  ## Build (QQ_s, HH_s, SS_s) per regime
  regime_covs <- vector("list", h)
  for (s in seq_len(h)) {
    RR_s  <- RR_list[[s]]
    DD_s  <- DD_list[[s]]
    Se_s  <- Sigma_e_by_regime[[s]]
    regime_covs[[s]] <- list(
      QQ      = tcrossprod(RR_s %*% Se_s, RR_s),
      HH      = tcrossprod(DD_s %*% Se_s, DD_s),
      SS      = RR_s %*% Se_s %*% t(DD_s),
      Sigma_e = Se_s
    )
  }

  ## ---- prepare observations -----------------------------------------------
  if (is.null(dim(data))) data <- matrix(data, nrow = n_obs)
  if (nrow(data) != n_obs) data <- t(data)
  n_T <- ncol(data)

  ## Observable SS mean for regime 1 (used as reference; each regime has its own d_s)
  ## For the innovation, we use the FROM-regime's d (= dr_s$ys[obs_vars]).

  ## ---- initialise state distributions (one per regime) --------------------
  ## PER-REGIME P_{0|0}: regime s starts at Lyapunov(TT_s, QQ_s), its OWN
  ## unconditional state covariance.  Here BOTH the transition matrix and the
  ## shock covariance are regime-specific, so sharing regime 1's P0 is even
  ## further off than in the shock-variance-only filter.  See the note in
  ## ms_kim_filter() for the P = I oracle that pins this.
  QQ_list  <- lapply(regime_covs, `[[`, "QQ")
  P0_list  <- .ms_init_P0(TT_list, QQ_list, lik_init)
  lik_init <- attr(P0_list, "lik_init")

  ## ---- component storage ---------------------------------------------------
  ## Same contract as ms_kim_filter(): the recursion carries M Gaussian
  ## COMPONENTS, component m being the law of s_{t-1} given y_{1:t-1} and a
  ## regime SUFFIX (the last one regime under gpb2, the last two under gpb3),
  ## with `comp_from[m]` the regime of s_{t-1} -- the only part of the suffix
  ## the period-t update reads -- and `comp_mass[m]` the suffix probability.
  ## GPB(2) and GPB(3) differ only in which coordinate of the post-update cell
  ## (m, k) is collapsed away, so `collapse = "gpb2"` is bit-identical to the
  ## pre-F3-C recursion.  See ms_kim_filter() for the full note.
  ##
  ## Initialise: s_0 = 0, P_0^{(s)} = P0_list[[s]]
  gpb3 <- identical(collapse, "gpb3")
  s0   <- numeric(n_state)
  Beta <- replicate(h, s0, simplify = FALSE)
  Pvar <- P0_list
  attr(Pvar, "lik_init") <- NULL

  comp_from <- seq_len(h)
  comp_mass <- ms_dr$pi0             # Pr[s_0 = i | y_{1:0}]
  n_new     <- if (gpb3) h * h else h
  from_new  <- if (gpb3) rep(seq_len(h), each = h) else seq_len(h)

  ## ---- Kim-Nelson loop (structural version) --------------------------------
  ##
  ## Key difference from shock-variance-only filter:
  ##   - Innovation uses FROM-regime i measurement matrices (ZZ_i, d_i)
  ##   - State prediction uses TO-regime j transition matrix (TT_j)
  ##   - Covariances use FROM-component P_m but TO-regime j QQ_j, HH_j, SS_j
  ##
  ## Following the comment block in ms-filter.R (dynhr lagged-state convention):
  ##   v_ij   = y_t - ZZ_i * b_i - d_i   (FROM-regime i measurement)
  ##   F_ij   = ZZ_i * P_i * ZZ_i' + HH_j (P_i from prior, HH_j to-regime)
  ##   K_ij   = (TT_j * P_i * ZZ_i' + SS_j) * F_ij^{-1}
  ##   b_hat  = TT_j * b_i + K_ij * v_ij
  ##   P_hat  = Joseph-form with TT_j - K_ij * ZZ_i and RR_j - K_ij * DD_j
  ##
  ## WHICH REGIME OWNS DD?  The law this filter assumes is
  ##   y_t = d_i + ZZ_i s_{t-1} + DD_j eps_t,   s_t = TT_j s_{t-1} + RR_j eps_t
  ## i.e. the state-dated blocks (ZZ_i, d_i) are FROM-regime and every block
  ## multiplying eps_t is TO-regime.  That is what F_ij already encodes
  ## (HH_j = DD_j Sigma_j DD_j'), so DD_j -- NOT DD_i -- must appear
  ## everywhere eps_t is projected out:
  ##   SS_j  = RR_j Sigma_j DD_j'         (cov(s_t, y_t) shock part)
  ##   P_hat = IKZ P_i IKZ' + (RR_j - K DD_j) Sigma_j (RR_j - K DD_j)'
  ## because s_t - b_hat = (TT_j - K ZZ_i)(s_{t-1} - b_i) + (RR_j - K DD_j) eps.
  ## Using DD_i in the Joseph term contradicts the F_ij used to build K and
  ## makes P_hat wrong on every OFF-DIAGONAL (i != j) path -- see F2-D.

  ll_const <- -0.5 * n_obs * log(2 * pi)
  me_diag  <- me_variance * diag(n_obs)
  loglik   <- 0
  ll_floor <- -1e300

  if (return_regime_probs)
    reg_prob_out <- matrix(0, h, n_T)

  ## Optional: per-period collapse-quality statistic (length T, nats).
  ## Period 1 is 0 by definition -- nothing has been collapsed yet.  Under
  ## gpb3 it is the GPB(4)-vs-GPB(3) one-step gap.
  if (return_collapse_diag) {
    collapse_out <- numeric(n_T)
    cd_HH   <- lapply(regime_covs, `[[`, "HH")
    prev_bh <- NULL
  }

  ## Optional state path for ms_kim_smoother_struct().  Same contract as
  ## ms_kim_filter(): per post-update CELL (component m -> to-regime k) store
  ## the five Durbin-Koopman blocks the backward pass needs, all of them
  ## by-products of the measurement update done below.  The regime dependence
  ## is BOTH-sided here:
  ##   a  = ZZ_i' F^{-1} v          (FROM-regime measurement loading)
  ##   M  = ZZ_i' F^{-1} ZZ_i
  ##   L  = TT_j - K ZZ_i           (TO-regime transition, FROM-regime Z)
  ##   G  = (RR_j - K DD_j)'      (TO-regime shock loading -- see above)
  ##   du = DD_j' F^{-1} v        (cov(eps_t, v) = Sigma_j DD_j')
  ## On an all-missing or singular-F period the blocks degenerate to
  ## (a, M, du) = 0 and (L, G) = (TT_j, RR_j'), i.e. the no-observation branch
  ## r_{t-1} = TT_j' r_t, N_{t-1} = TT_j' N_t TT_j, eps = Sigma_e RR_j' r.
  ## NOTE the null block is per TO-regime j, unlike the shock-variance-only
  ## filter where TT and RR are common.
  ## joint_filt[i, j, t] = Pr[s_{t-1} = i, s_t = j | y_{1:t}]; see the note in
  ## ms_kim_filter() -- consumed by the joint regime pass of
  ## ms_kim_smoother_struct().
  ## cell_filt[[t]]: the FULL per-CELL posterior (M_t x h) = the per-TRIPLE
  ## filtered joint Pr[s_{t-2}, s_{t-1}, s_t | y_{1:t}] under gpb3; see the
  ## note in ms_kim_filter().  Consumed by the GPB(3) backward pass; stored
  ## only under gpb3.
  if (return_state_path) {
    beta_path <- array(0, c(n_state, n_new, n_T))
    P_path    <- array(0, c(n_state, n_state, n_new, n_T))
    dk_path   <- vector("list", n_T)
    joint_filt <- array(0, c(h, h, n_T))
    cell_filt  <- if (gpb3) vector("list", n_T) else NULL
    .dk_null_j <- lapply(seq_len(h), function(j)
      list(a  = numeric(n_state),
           M  = matrix(0, n_state, n_state),
           L  = TT_list[[j]],
           G  = t(RR_list[[j]]),
           du = numeric(n_exo)))
  }

  for (t in seq_len(n_T)) {
    y_t    <- data[, t]
    obs_ok <- is.finite(y_t)
    all_na <- !any(obs_ok)

    M      <- length(comp_mass)
    n_cell <- M * h
    n_coll <- n_cell %/% n_new

    log_lik_joint <- matrix(-Inf, M, h)   # [component m, to-regime k]
    beta_hat  <- matrix(0, n_state, n_cell)
    p_hat     <- vector("list", n_cell)
    ## dk_t[[m + (k - 1) * M]] matches log_lik_joint[m, k] (column-major).
    if (return_state_path)
      dk_t <- .dk_null_j[rep(seq_len(h), each = M)]

    for (m in seq_len(M)) {
      i     <- comp_from[m]
      b_i   <- Beta[[m]]
      P_i   <- Pvar[[m]]
      TT_j_list <- TT_list   # reference to avoid repeated indexing

      ## FROM-regime i measurement matrices
      ZZ_i  <- ZZ_list[[i]]
      DD_i  <- DD_list[[i]]
      d_i   <- d_list[[i]]
      tZZ_i <- t(ZZ_i)

      for (j in seq_len(h)) {
        cell  <- m + (j - 1L) * M
        TT_j  <- TT_j_list[[j]]
        RR_j  <- RR_list[[j]]
        DD_j  <- DD_list[[j]]
        cov_j <- regime_covs[[j]]
        QQ_j  <- cov_j$QQ
        HH_j  <- cov_j$HH
        SS_j  <- cov_j$SS

        ## Prediction: TO-regime j transition
        beta_pred <- drop(TT_j %*% b_i)
        P_pred    <- tcrossprod(TT_j %*% P_i, TT_j) + QQ_j
        P_pred    <- (P_pred + t(P_pred)) * 0.5

        if (all_na) {
          beta_hat[, cell]    <- beta_pred
          p_hat[[cell]]       <- P_pred
          log_lik_joint[m, j] <- log(P[i, j] * comp_mass[m])
        } else {
          ## Innovation (FROM-regime i): v = y_t - ZZ_i * b_i - d_i
          v_full  <- y_t - as.numeric(ZZ_i %*% b_i) - d_i
          v_obs   <- v_full
          if (any(!obs_ok)) v_obs[!obs_ok] <- 0

          ## Innovation covariance (FROM-regime i P_i, TO-regime j HH_j)
          F_ij  <- ZZ_i %*% P_i %*% tZZ_i + HH_j + me_diag
          F_ij  <- (F_ij + t(F_ij)) * 0.5

          if (any(!obs_ok)) {
            F_ij_obs <- F_ij[obs_ok, obs_ok, drop = FALSE]
            v_use    <- v_obs[obs_ok]
            n_obs_t  <- sum(obs_ok)
          } else {
            F_ij_obs <- F_ij
            v_use    <- v_obs
            n_obs_t  <- n_obs
          }

          Fc_ij <- tryCatch(chol(F_ij_obs), error = function(e) NULL)
          if (is.null(Fc_ij)) {
            log_lik_joint[m, j] <- -Inf
            beta_hat[, cell]    <- beta_pred
            p_hat[[cell]]       <- P_pred
            next
          }

          Fi_ij     <- chol2inv(Fc_ij)
          log_det_F <- 2 * sum(log(diag(Fc_ij)))
          ll_t      <- -0.5 * n_obs_t * log(2 * pi)
          quad      <- drop(crossprod(v_use, Fi_ij %*% v_use))
          lp_ij     <- ll_t - 0.5 * (log_det_F + quad)

          log_lik_joint[m, j] <- lp_ij + log(P[i, j] * comp_mass[m])

          ## Kalman gain: K_ij = (TT_j * P_i * ZZ_i' + SS_j) * F_ij^{-1}
          if (any(!obs_ok)) {
            ZZ_i_obs  <- ZZ_i[obs_ok, , drop = FALSE]
            DD_j_obs  <- DD_j[obs_ok, , drop = FALSE]
            SS_j_obs  <- RR_j %*% cov_j$Sigma_e %*% t(DD_j_obs)
            K_ij      <- (TT_j %*% P_i %*% t(ZZ_i_obs) + SS_j_obs) %*% Fi_ij
            b_upd     <- beta_pred + drop(K_ij %*% v_use)
            IKZ       <- TT_j - K_ij %*% ZZ_i_obs
            RmKD      <- RR_j - K_ij %*% DD_j_obs
          } else {
            K_ij  <- (TT_j %*% P_i %*% tZZ_i + SS_j) %*% Fi_ij
            b_upd <- beta_pred + drop(K_ij %*% v_use)
            IKZ   <- TT_j - K_ij %*% ZZ_i
            RmKD  <- RR_j - K_ij %*% DD_j
          }

          ## Joseph-form updated covariance, INCLUDING the K (me I) K' term:
          ## me_variance is TRUE i.i.d. observation noise here, not an
          ## F-only regulariser (see ms_kim_filter()'s roxygen).
          P_upd <- tcrossprod(IKZ %*% P_i, IKZ) +
                   tcrossprod(RmKD %*% cov_j$Sigma_e, RmKD)
          ## Guarded so me_variance = 0 stays bit-identical to the old code.
          if (me_variance > 0) P_upd <- P_upd + me_variance * tcrossprod(K_ij)
          P_upd <- (P_upd + t(P_upd)) * 0.5

          beta_hat[, cell]  <- b_upd
          p_hat[[cell]]     <- P_upd

          if (return_state_path) {
            ZZ_use <- if (any(!obs_ok)) ZZ_i[obs_ok, , drop = FALSE] else ZZ_i
            DD_use <- if (any(!obs_ok)) DD_j[obs_ok, , drop = FALSE] else DD_j
            Fv     <- Fi_ij %*% v_use
            dk_t[[cell]] <- list(
              a  = as.numeric(crossprod(ZZ_use, Fv)),
              M  = crossprod(ZZ_use, Fi_ij %*% ZZ_use),
              L  = IKZ,
              G  = t(RmKD),
              du = as.numeric(crossprod(DD_use, Fv))
            )
          }
        }
      }
    }

    ## Log-sum-exp for period t
    mx <- max(log_lik_joint)
    if (!is.finite(mx)) {
      loglik <- -Inf
      break
    }
    log_f_y <- .logsumexp(log_lik_joint)
    loglik  <- loglik + log_f_y
    if (loglik < ll_floor) { loglik <- -Inf; break }

    ## Collapse-quality statistic: what the t-1 collapse cost THIS period.
    if (return_collapse_diag && !is.null(prev_bh) && !all_na) {
      lfd <- .ms_defer_logf(y_t, obs_ok, me_diag, prev_bh, prev_ph, prev_pj,
                            prev_from, P, ZZ_list, d_list, cd_HH)
      collapse_out[t] <- if (is.na(lfd)) NA_real_ else lfd - log_f_y
    }

    ## Hamilton filter: posterior cell probs, then group masses
    prob_joint     <- exp(log_lik_joint - log_f_y)
    group_mass_raw <- colSums(array(prob_joint, c(n_coll, n_new)))
    group_mass     <- pmax(group_mass_raw, 1e-300)
    group_mass     <- group_mass / sum(group_mass)

    ## Collapse: n_cell -> n_new (contiguous blocks of n_coll cells)
    Beta_new <- vector("list", n_new)
    Pvar_new <- vector("list", n_new)
    fb <- if (gpb3) ((seq_len(n_new) - 1L) * n_coll) %% M + 1L else seq_len(h)

    for (g in seq_len(n_new)) {
      ## Unreachable group: carry the prior moments forward (test the RAW
      ## mass, not the 1e-300 floor -- see ms_kim_filter() for why).
      if (group_mass_raw[g] <= 0) {
        Beta_new[[g]] <- Beta[[fb[g]]]
        Pvar_new[[g]] <- Pvar[[fb[g]]]
        next
      }
      cells <- (g - 1L) * n_coll + seq_len(n_coll)
      wts   <- prob_joint[cells] / group_mass[g]
      b_j <- numeric(n_state)
      for (u in seq_len(n_coll)) b_j <- b_j + wts[u] * beta_hat[, cells[u]]
      Beta_new[[g]] <- b_j

      P_j <- matrix(0, n_state, n_state)
      for (u in seq_len(n_coll)) {
        diff_ij <- beta_hat[, cells[u]] - b_j
        P_j <- P_j + wts[u] * (p_hat[[cells[u]]] + tcrossprod(diff_ij))
      }
      Pvar_new[[g]] <- (P_j + t(P_j)) * 0.5

    }

    Beta      <- Beta_new
    Pvar      <- Pvar_new
    comp_mass <- group_mass
    comp_from <- from_new

    regime_prob <- if (gpb3) colSums(matrix(group_mass, h, h)) else group_mass

    ## Carry this period's PRE-collapse cells for the next period's
    ## collapse-quality statistic (`ord` preserves the pre-F3-C enumeration
    ## order under gpb2, so the diagnostic is bit-identical too).
    if (return_collapse_diag) {
      ord <- as.integer(t(matrix(seq_len(n_cell), M, h)))
      prev_bh   <- beta_hat[, ord, drop = FALSE]
      prev_ph   <- p_hat[ord]
      prev_pj   <- as.numeric(prob_joint)[ord]
      prev_from <- rep(seq_len(h), each = M)[ord]
    }

    if (return_regime_probs)
      reg_prob_out[, t] <- regime_prob

    if (return_state_path) {
      for (g in seq_len(n_new)) {
        beta_path[, g, t] <- Beta[[g]]
        P_path[, , g, t]  <- Pvar[[g]]
      }
      dk_path[[t]] <- dk_t
      joint_filt[, , t] <- if (gpb3) matrix(group_mass_raw, h, h) else prob_joint
      if (gpb3) cell_filt[[t]] <- prob_joint
    }
  }

  out <- list(
    loglik       = loglik,
    regime_probs = if (return_regime_probs) reg_prob_out else NULL,
    n_obs        = n_obs,
    n_T          = n_T
  )

  if (return_collapse_diag) {
    out$collapse_diag <- collapse_out
    out$collapse_max  <- suppressWarnings(max(abs(collapse_out), na.rm = TRUE))
    if (!is.finite(out$collapse_max)) out$collapse_max <- NA_real_
  }

  if (return_state_path) {
    out$beta_filt    <- beta_path
    out$P_filt       <- P_path
    out$dk_path      <- dk_path
    out$joint_filt   <- joint_filt
    out$cell_filt    <- cell_filt
    out$collapse     <- collapse
    out$TT_list      <- TT_list
    out$RR_list      <- RR_list
    out$ZZ_list      <- ZZ_list
    out$DD_list      <- DD_list
    out$d_list       <- d_list
    out$QQ_list      <- QQ_list
    out$Sigma_e_list <- lapply(regime_covs, `[[`, "Sigma_e")
    out$P0_list      <- P0_list
    out$state_names  <- endo[state_idx]
    out$shock_names  <- exo
    out$lik_init     <- lik_init
  }
  out
}


## ============================================================================
## Shared initialisation helper
## ============================================================================

## Internal: PER-REGIME initial state covariance P_{0|0}^{(j)}.
##
## Kim & Nelson (1999) initialise every regime at its OWN unconditional state
## covariance, i.e. the solution of the regime-j Lyapunov equation
##   P_j = TT_j P_j TT_j' + QQ_j.
## Sharing regime 1's P0 across all regimes (the pre-2026-09-02 behaviour)
## breaks the P = I oracle: started in regime r with P = I, the Kim filter must
## reproduce the fixed-regime-r Kalman filter EXACTLY, and it cannot if it is
## handed the wrong regime's initial covariance.
##
## `lik_init` = "auto" is resolved ONCE, from regime 1, and the resulting
## family ("stationary" or "kappa") is then applied to EVERY regime: mixing a
## stationary P0 for one regime with a kappa-diffuse P0 for another would make
## the regime likelihoods incommensurable (the kappa offset is regime-specific)
## and the mixture weights meaningless.
##
## @param TT_list  A single n_state x n_state transition matrix (recycled to
##   all regimes) or a length-h list of per-regime transition matrices.
## @param QQ_list  Length-h list of per-regime state-noise covariances.
## @param lik_init "auto", "stationary" or "kappa".
## @return Length-h list of P_{0|0}^{(j)}, carrying the RESOLVED `lik_init`
##   as an attribute.
## @noRd
## Internal: GPB(2) COLLAPSE-QUALITY statistic, in LOG-LIKELIHOOD NATS.
##
## WHAT IS MEASURED.  At period t the Kim filter scores y_t against the h
## COLLAPSED priors (Beta[[i]], Pvar[[i]]), giving the period contribution
##   log f_t   = log sum_{i,j} Pr[s_{t-1}=i | y_{1:t-1}] P[i,j] N(y_t; ...).
## Had the h^2 -> h collapse at t-1 been DEFERRED by one period, the same
## contribution would instead be scored against the h^2 uncollapsed
## components (b_{ki}, P_{ki}) that the collapse merged:
##   log f_t^* = log sum_{k,i,j} Pr[s_{t-2}=k, s_{t-1}=i | y_{1:t-1}]
##                              P[i,j] N(y_t; ...).
## The statistic is the difference
##   collapse_diag[t] = log f_t^* - log f_t     (nats),
## i.e. exactly how many nats of period-t log-likelihood the previous
## period's collapse cost.  It is the GPB(3)-vs-GPB(2) one-step gap, so it is
## IDENTICALLY ZERO whenever the collapse is lossless -- in particular at
## every period under P = I, where each (k,i) mass is a point mass -- and it
## is measured in the only units that matter, those of the likelihood itself.
## Unlike a KL or a spread-share statistic it does not saturate when the
## per-path covariances are singular, which is precisely the regime in which
## the collapse fails.
##
## WHY IT EXISTS (F2-D).  Kim's collapse error is usually negligible, but it
## is NOT bounded.  On the rbc two-structural-regime fixture with
## n_obs = n_exo = 1 and me_variance = 0 the model is EXACTLY identified, the
## per-path covariances collapse toward singular, and on one draw (seed 11,
## P = .9/.9) the filter log-likelihood misses the exact all-path-enumeration
## value by 48 nats over 8 periods.  A full GPB(3) -- collapsing on the last
## TWO regimes throughout -- reproduces the exact value to 1e-15 on that same
## draw, which is what identifies the h^2 -> h collapse, and not the
## recursion, as the entire source of the error.  Raising me_variance off
## zero relieves the degeneracy (the GPB(2)-vs-GPB(3) gap on that draw falls
## from 48 nats to 5.6 at me_variance = 1e-2, and to 0.01 at 1e-1; those three
## figures were measured before F3-A made me_variance TRUE i.i.d. measurement
## error rather than an F-only regulariser).  There is
## no fix available inside GPB(2); this statistic exists so the breakdown is
## DETECTABLE instead of silent.
##
## Cost is h^3 Gaussian density evaluations per period (no gains, no Joseph
## form -- only the log density is needed), incurred only when the caller asks
## for it.  It is READ-ONLY with respect to the recursion, so the
## log-likelihood is bit-identical with and without the flag.
##
## @param y_t       Observation vector at t, ALREADY offset-free if the caller
##                  pre-subtracts d (the reduced-form filter does).
## @param obs_ok    Logical vector of non-missing rows of y_t.
## @param me_diag   n_obs x n_obs measurement-error matrix (added to F).
## It generalises verbatim to GPB(3): the caller hands it the PRE-collapse
## cells of period t-1 whatever their number, so under `collapse = "gpb3"` the
## same difference is the GPB(4)-vs-GPB(3) one-step gap -- the honest
## collapse-quality statistic for that filter, and expected to be ~0 wherever
## GPB(3) is already exact.
##
## @param y_t       Observation vector at t, ALREADY offset-free if the caller
##                  pre-subtracts d (the reduced-form filter does).
## @param obs_ok    Logical vector of non-missing rows of y_t.
## @param me_diag   n_obs x n_obs measurement-error matrix (added to F).
## @param beta_prev n_state x Mc matrix of PRE-collapse cell means from t-1.
## @param P_prev    Length-Mc list of PRE-collapse cell covariances from t-1.
## @param mass_prev Length-Mc vector of cell posterior masses at t-1.
## @param from_prev Length-Mc integer vector: the regime of s_{t-1} carried by
##                  each cell (the FROM regime of the period-t update).
## @param P         h x h transition matrix.
## @param ZZ_list   Length-h list of FROM-regime measurement loadings.
## @param d_list    Length-h list of FROM-regime observable intercepts.
## @param HH_list   Length-h list of TO-regime observation-noise covariances.
## @return Scalar log-density of y_t under the DEFERRED (uncollapsed) prior,
##   or `NA_real_` if no path is finite.
## @noRd
.ms_defer_logf <- function(y_t, obs_ok, me_diag, beta_prev, P_prev, mass_prev,
                           from_prev, P, ZZ_list, d_list, HH_list) {
  h  <- nrow(P)
  Mc <- length(mass_prev)
  lp <- rep(-Inf, Mc * h)
  n_use <- sum(obs_ok)
  idx <- 0L
  for (c in seq_len(Mc)) {
    m_ki <- mass_prev[c]
    i    <- from_prev[c]
    ZZ_i <- ZZ_list[[i]]; d_i <- d_list[[i]]
    b_ki <- beta_prev[, c]; P_ki <- P_prev[[c]]
    ZPZ  <- ZZ_i %*% P_ki %*% t(ZZ_i)
    v_f  <- y_t - as.numeric(ZZ_i %*% b_ki) - d_i
    for (j in seq_len(h)) {
      idx <- idx + 1L
      if (m_ki <= 0 || P[i, j] <= 0) next
      F_ij <- ZPZ + HH_list[[j]] + me_diag
      F_ij <- (F_ij + t(F_ij)) * 0.5
      Fc <- tryCatch(chol(F_ij[obs_ok, obs_ok, drop = FALSE]),
                     error = function(e) NULL)
      if (is.null(Fc)) next
      v <- v_f[obs_ok]
      q <- sum(backsolve(Fc, v, transpose = TRUE)^2)
      lp[idx] <- -0.5 * (n_use * log(2 * pi) + 2 * sum(log(diag(Fc))) + q) +
        log(m_ki * P[i, j])
    }
  }
  if (!any(is.finite(lp))) return(NA_real_)
  .logsumexp(lp)
}


.ms_init_P0 <- function(TT_list, QQ_list, lik_init) {
  h <- length(QQ_list)
  if (!is.list(TT_list)) TT_list <- replicate(h, TT_list, simplify = FALSE)

  if (identical(lik_init, "auto")) {
    P0_try <- tryCatch(solve_lyapunov(TT_list[[1L]], QQ_list[[1L]]),
                       error = function(e) NULL)
    ok_stat <- !is.null(P0_try) && all(is.finite(P0_try)) &&
      min(Re(eigen((P0_try + t(P0_try)) / 2, symmetric = TRUE,
                   only.values = TRUE)$values)) > -1e-8
    lik_init <- if (ok_stat) "stationary" else "kappa"
  }

  P0 <- vector("list", h)
  for (s in seq_len(h)) {
    TT_s <- TT_list[[s]]
    QQ_s <- QQ_list[[s]]
    if (identical(lik_init, "stationary")) {
      P_s <- tryCatch(solve_lyapunov(TT_s, QQ_s), error = function(e) NULL)
      if (is.null(P_s) || anyNA(P_s) || !all(is.finite(P_s)))
        P_s <- .build_P0(TT_s, QQ_s)
    } else {
      P_s <- .build_P0(TT_s, QQ_s)
    }
    P0[[s]] <- P_s
  }
  attr(P0, "lik_init") <- lik_init
  P0
}
