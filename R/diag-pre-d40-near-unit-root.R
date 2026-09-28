## R/diag-pre-d40-near-unit-root.R
## --------------------------------------------------------------------------
## D40 Near-unit-root diagnostic.
##
## Classifies every root of the first-order STATE transition
## TT = ghx[state_idx, ] (not the full QZ spectrum) as stationary / near-unit
## (half-life long relative to the sample) / unit / explosive, attributes the
## dominant root to state variables (participation factors) and to parameters
## (central finite differences of the spectral radius), and draws the roots in
## the complex plane with the unit circle and the half-life threshold circle.
## --------------------------------------------------------------------------

#' D40. Near-unit-root check on the state transition
#'
#' Solves the model at \code{theta} (or takes a supplied decision rule) and
#' inspects the eigenvalues of the first-order state transition
#' \code{TT = ghx[state_idx, ]} -- the matrix the Kalman filter propagates and
#' whose Lyapunov equation gives the stationary initial covariance \code{P0}.
#' Roots are judged by their MODULUS (so a complex pair counts as persistent
#' when \eqn{|\lambda|} is near 1, whatever its real part).
#'
#' @section Classification:
#' Each root \eqn{\lambda} is given a half-life
#' \eqn{h = \log(0.5) / \log|\lambda|} periods (0 for \eqn{\lambda = 0},
#' \code{Inf} for \eqn{|\lambda| \ge 1}) and a class:
#' \describe{
#'   \item{\code{explosive}}{\eqn{|\lambda| > 1 + } \code{unit_tol}: no stable
#'     solution along this direction.}
#'   \item{\code{unit}}{\eqn{|\lambda| \ge 1 - } \code{unit_tol}.  This is the
#'     same cut \code{kalman_filter(lik_init = "auto")} uses before trying the
#'     Lyapunov solve, and the cut \code{run_all_diagnostics()} uses to skip
#'     D41.  The stationary \code{P0} and unconditional moments are not
#'     defined (or not usable); the filter takes the exact-diffuse path and
#'     moment-based diagnostics are invalid.}
#'   \item{\code{near_unit}}{stationary but
#'     \eqn{h >} \code{halflife_frac * n_obs}: the root takes a sizeable share
#'     of the sample to decay, so the stationary \code{P0} is a very wide
#'     prior on the initial state and sample moments are poor estimates of
#'     the model's unconditional moments.  A root of 0.995 (h = 138) is
#'     near-unit at \code{n_obs = 200} but not at \code{n_obs = 20000}.  The
#'     equivalent modulus threshold is
#'     \eqn{0.5^{1 / (\mathrm{halflife\_frac} \cdot n_{obs})}}.}
#'   \item{\code{stationary}}{everything else.}
#' }
#' Badge: FAIL if any root is \code{explosive} or \code{unit} (the stationary
#' \code{P0} and the unconditional moments do not exist, so the moment-based
#' diagnostics downstream are invalid); \strong{WARN} if the worst root is only
#' \code{near_unit}; PASS otherwise.  Without \code{n_obs} the near-unit class
#' cannot be assigned: the badge is then FAIL for unit/explosive roots and
#' INFO otherwise.
#'
#' Near-unit is WARN, not FAIL, because long persistence is a legitimate model
#' feature, not a defect: Dynare (\code{qz_criterium},
#' \code{model_diagnostics}) and IRIS both route a persistent/unit-root system
#' to a diffuse-initialisation branch rather than rejecting it, and the
#' econometric response to a local-to-unity root is non-standard inference, not
#' a different model.  What changes is how you should filter and how much to
#' trust long-run responses -- exactly what a soft badge is for.
#'
#' @section Attribution:
#' \describe{
#'   \item{State participation}{For the dominant root with right eigenvector
#'     \eqn{w} and left eigenvector \eqn{v}, the participation factor of
#'     state \eqn{k} is \eqn{|v_k w_k| / \sum_j |v_j w_j|}
#'     (small-signal-stability convention); it says which state variables the
#'     persistent mode lives in.  For a complex pair the value is the same for
#'     both members.}
#'   \item{Parameter sensitivity}{\eqn{d\rho_{\max}/d\theta_i} by central
#'     finite differences (step \code{h}), re-solving the model at
#'     \eqn{\theta \pm h e_i}.  \code{NA} when a re-solve fails or violates
#'     Blanchard-Kahn, or when \code{compiled} is not supplied.}
#' }
#'
#' @param model        A \code{dynhr_mod} object (output of \code{parse_mod}).
#' @param compiled     A \code{dynhr_compiled} (output of
#'   \code{compile_model}).  Required for solving and for the parameter
#'   sensitivities; may be \code{NULL} when \code{dr} is supplied (then the
#'   sensitivities are skipped).
#' @param theta        Named numeric vector of parameter values.  If
#'   \code{NULL}, uses \code{model$param_values}.
#' @param params       Alias for \code{theta} (either may be supplied).
#' @param dr           Optional first-order decision rule
#'   (\code{solve_perturbation()} output) at \code{theta}; skips the baseline
#'   solve.
#' @param n_obs        Sample length \eqn{T} used for the half-life criterion.
#'   \code{NULL} (default) disables the \code{near_unit} class.
#' @param param_names  Character vector selecting which parameters to
#'   differentiate.  Defaults to all names in \code{theta};
#'   \code{character(0)} skips the sensitivities.
#' @param h            Finite-difference step size (default \code{1e-5}).
#' @param unit_tol     Tolerance for an exact unit root (default \code{1e-6},
#'   matching \code{kalman_filter()} and the D41 skip rule).
#' @param halflife_frac A root is \code{near_unit} when its half-life exceeds
#'   this fraction of \code{n_obs} (default \code{0.25}).  \strong{0.25 is a
#'   package choice with no literature source}: no paper or toolbox defines a
#'   half-life-versus-sample-length adequacy cutoff.  It is a heuristic for
#'   "this root takes a sizeable share of the sample to decay"; adjust it
#'   freely.
#' @param meta         Optional \code{diag_meta()} object for plot captions.
#' @param ...          Unused; reserved for future arguments.
#'
#' @return A \code{dynhr_diagnostic} list (see \code{.make_result}) with
#'   \code{$result} containing:
#'   \itemize{
#'     \item \code{roots} -- data.frame (one row per state root: \code{re},
#'       \code{im}, \code{modulus}, \code{half_life}, \code{class}), sorted by
#'       decreasing modulus
#'     \item \code{spectral_radius}, \code{distance_to_unit}
#'       (\code{1 - spectral_radius}), \code{critical_eigenvalue}
#'     \item \code{half_life} -- half-life of the dominant root (periods)
#'     \item \code{n_obs}, \code{halflife_frac}, \code{unit_tol},
#'       \code{modulus_threshold} (\code{NA} without \code{n_obs})
#'     \item \code{n_unit}, \code{n_explosive}, \code{n_near_unit}
#'     \item \code{stationary_p0_valid} -- \code{FALSE} iff some root has
#'       \eqn{|\lambda| > 1 - } \code{unit_tol} (the D41 skip rule)
#'     \item \code{state_participation} -- data.frame (\code{state},
#'       \code{participation}) for the dominant root
#'     \item \code{param_sensitivity} -- data.frame (\code{param},
#'       \code{d_rhomax_dparam}) sorted by decreasing absolute value
#'   }
#'   and \code{$plots$roots} (complex plane) plus, when any sensitivity is
#'   finite, \code{$plots$sensitivity}.
#'
#' @export
diag_near_unit_root <- function(model,
                                compiled      = NULL,
                                theta         = NULL,
                                params        = NULL,
                                dr            = NULL,
                                n_obs         = NULL,
                                param_names   = NULL,
                                h             = 1e-5,
                                unit_tol      = 1e-6,
                                halflife_frac = 0.25,
                                meta          = NULL,
                                ...) {

  ## ---- 0. Normalise inputs -----------------------------------------------
  if (!is.null(n_obs) &&
      (length(n_obs) != 1L || !is.finite(n_obs) || n_obs <= 0))
    .dynhr_abort("D40: `n_obs` must be a single positive number or NULL.")
  if (length(unit_tol) != 1L || !is.finite(unit_tol) || unit_tol < 0)
    .dynhr_abort("D40: `unit_tol` must be a single non-negative number.")
  if (length(halflife_frac) != 1L || !is.finite(halflife_frac) ||
      halflife_frac <= 0)
    .dynhr_abort("D40: `halflife_frac` must be a single positive number.")

  if (is.null(theta)) theta <- params
  if (is.null(theta)) theta <- model$param_values
  if (is.null(theta) && is.null(dr)) {
    return(.make_result(
      result = NULL, pass = NA,
      summary = "D40 Near-unit-root: supply theta, params or dr."
    ))
  }
  if (!is.null(theta)) theta <- unlist(theta)
  if (!is.null(theta) && is.null(names(theta)) && !is.null(model$param_names))
    names(theta) <- model$param_names
  if (is.null(param_names)) param_names <- names(theta)
  if (is.null(param_names)) param_names <- character(0)

  ## Re-solve at a parameter vector.  Returns NULL on failure or on a BK
  ## violation (the decision rule is then not a stable solution).  The
  ## steady state uses the model's own initval (no zero y0 override).
  .solve_at <- function(th) {
    params_full <- model$param_values
    for (nm in names(th)) params_full[[nm]] <- th[[nm]]
    ss <- tryCatch(solve_steady(compiled, params_full, verbose = FALSE),
                   error = function(e) .dynhr_reraise_bug(e, NULL))
    if (is.null(ss) || !isTRUE(ss$converged)) return(NULL)
    d <- tryCatch(suppressWarnings(
      solve_perturbation(model, compiled, ss = ss$values,
                         params = params_full, order = 1L, verbose = FALSE)),
      error = function(e) .dynhr_reraise_bug(e, NULL))
    if (is.null(d) || identical(d$bk_satisfied, FALSE)) return(NULL)
    d
  }

  ## ---- 1. Baseline decision rule ------------------------------------------
  if (is.null(dr)) {
    if (is.null(compiled)) {
      return(.make_result(
        result = NULL, pass = NA,
        summary = "D40 Near-unit-root: supply `compiled` (or a solved `dr`)."
      ))
    }
    dr <- .solve_at(theta)
    if (is.null(dr)) {
      return(.make_result(
        result = NULL, pass = FALSE,
        summary = paste0("D40 Near-unit-root: baseline solve failed (steady ",
                         "state did not converge or Blanchard-Kahn violated).")
      ))
    }
  }

  ## ---- 2. Roots of the STATE transition -----------------------------------
  T_mat <- .d40_state_transition(dr)
  n_s   <- nrow(T_mat)
  thr   <- if (is.null(n_obs)) NA_real_ else 0.5^(1 / (halflife_frac * n_obs))

  if (n_s == 0L) {
    return(.make_result(
      result = list(
        roots = .d40_classify(complex(0), n_obs, unit_tol, halflife_frac),
        spectral_radius = 0, distance_to_unit = 1,
        critical_eigenvalue = complex(0), half_life = 0,
        n_obs = n_obs, halflife_frac = halflife_frac, unit_tol = unit_tol,
        modulus_threshold = thr, n_unit = 0L, n_explosive = 0L,
        n_near_unit = 0L, stationary_p0_valid = TRUE,
        state_participation = data.frame(state = character(0),
                                         participation = numeric(0)),
        param_sensitivity = data.frame(param = param_names,
                                       d_rhomax_dparam = NA_real_,
                                       stringsAsFactors = FALSE)
      ),
      pass    = TRUE,
      summary = "D40 Near-unit-root: no state variables; nothing persists."
    ))
  }

  eg   <- eigen(T_mat, symmetric = FALSE)
  ev   <- eg$values
  roots <- .d40_classify(ev, n_obs, unit_tol, halflife_frac)
  crit  <- which.max(Mod(ev))[1L]
  crit_ev  <- ev[crit]
  rho_base <- Mod(crit_ev)
  hl_base  <- roots$half_life[1L]

  n_unit      <- sum(roots$class == "unit")
  n_explosive <- sum(roots$class == "explosive")
  n_near      <- sum(roots$class == "near_unit")

  ## ---- 3. State participation of the dominant root ------------------------
  state_names <- colnames(T_mat) %||% paste0("s", seq_len(n_s))
  part <- .d40_participation(T_mat, eg, crit)
  state_part <- data.frame(state = state_names, participation = part,
                           stringsAsFactors = FALSE)
  state_part <- state_part[order(state_part$participation,
                                 decreasing = TRUE), ]
  rownames(state_part) <- NULL

  ## ---- 4. Parameter sensitivity d(rho_max)/d(theta_i) ---------------------
  n_par <- length(param_names)
  dsens <- rep(NA_real_, n_par)
  if (!is.null(compiled)) {
    for (i in seq_len(n_par)) {
      pn <- param_names[i]
      if (!(pn %in% names(theta)) || !is.finite(theta[[pn]])) next
      th_p <- theta; th_p[[pn]] <- th_p[[pn]] + h
      th_m <- theta; th_m[[pn]] <- th_m[[pn]] - h
      dr_p <- .solve_at(th_p)
      dr_m <- .solve_at(th_m)
      if (is.null(dr_p) || is.null(dr_m)) next
      T_p <- .d40_state_transition(dr_p)
      T_m <- .d40_state_transition(dr_m)
      if (nrow(T_p) != n_s || nrow(T_m) != n_s) next
      dsens[i] <- (max(Mod(eigen(T_p, symmetric = FALSE,
                                 only.values = TRUE)$values)) -
                   max(Mod(eigen(T_m, symmetric = FALSE,
                                 only.values = TRUE)$values))) / (2 * h)
    }
  }
  param_sens_df <- data.frame(param = param_names, d_rhomax_dparam = dsens,
                              stringsAsFactors = FALSE)
  param_sens_df <- param_sens_df[order(abs(dsens), decreasing = TRUE,
                                       na.last = TRUE), , drop = FALSE]
  rownames(param_sens_df) <- NULL

  ## ---- 5. Badge + text ----------------------------------------------------
  ## Unit / explosive roots break the stationary P0 and the unconditional
  ## moments -> FAIL. A merely persistent (near-unit) root is a modelling
  ## fact, not an error -> WARN.
  n_hard <- n_unit + n_explosive
  pass <- if (n_hard > 0L) FALSE
          else if (n_near > 0L) TRUE
          else if (is.null(n_obs)) NA
          else TRUE
  warn <- (n_hard == 0L) && (n_near > 0L)
  badge <- .badge_str(list(pass = pass, errored = FALSE, warn = warn))

  fmt_hl <- function(x) if (is.finite(x)) sprintf("%.1f", x) else "Inf"
  ## 4 significant digits, not 6 decimal places: "spectral radius 0.850000"
  ## claims precision the eigensolver does not have, and a numerically-zero
  ## root printed as the signed zero "-0.000000".
  crit_str <- if (abs(Im(crit_ev)) > 1e-12)
    sprintf("%s%si", .fmt_sig(Re(crit_ev), zap = 1e-12),
            if (Im(crit_ev) >= 0) paste0("+", .fmt_sig(Im(crit_ev)))
            else .fmt_sig(Im(crit_ev)))
  else .fmt_sig(Re(crit_ev), zap = 1e-12)

  verdict <- if (n_explosive > 0L) {
    sprintf("%d explosive root(s): no stable solution", n_explosive)
  } else if (n_unit > 0L) {
    sprintf(paste0("%d unit root(s) (|lambda| >= 1 - %g): stationary P0 and ",
                   "unconditional moments undefined; the Kalman filter uses ",
                   "the exact-diffuse init and D41/moment diagnostics are ",
                   "invalid"), n_unit, unit_tol)
  } else if (n_near > 0L) {
    sprintf(paste0("%d near-unit root(s): half-life > %g x T = %.1f periods; ",
                   "stationary P0 is very diffuse and sample moments are ",
                   "unreliable"), n_near, halflife_frac, halflife_frac * n_obs)
  } else if (is.null(n_obs)) {
    "no unit root; supply n_obs for the half-life-vs-T check"
  } else {
    sprintf("all half-lives <= %g x T = %.1f periods", halflife_frac,
            halflife_frac * n_obs)
  }

  top3 <- utils::head(param_sens_df[is.finite(param_sens_df$d_rhomax_dparam), ,
                                    drop = FALSE], 3L)
  drivers <- if (nrow(top3) > 0L)
    paste(sprintf("%s=%.4g", top3$param, top3$d_rhomax_dparam),
          collapse = ", ")
  else "n/a"
  top_state <- paste(sprintf("%s=%.2f",
                             utils::head(state_part$state, 3L),
                             utils::head(state_part$participation, 3L)),
                     collapse = ", ")

  summary_text <- sprintf(
    paste0("D40 Near-unit-root: %d state root(s); spectral radius %s ",
           "(root %s, half-life %s periods)%s; %s. Dominant states: %s. ",
           "Top d(rho)/d(theta): %s."),
    n_s, .fmt_sig(rho_base), crit_str, fmt_hl(hl_base),
    if (is.null(n_obs)) "" else sprintf(", T = %g", n_obs),
    verdict, top_state, drivers)

  llm_text <- paste(c(
    sprintf("D40 | Near-unit-root check | %s", badge),
    sprintf(paste0("  spectral_radius=%.6f  half_life=%s  n_obs=%s  ",
                   "threshold_modulus=%s"),
            rho_base, fmt_hl(hl_base),
            if (is.null(n_obs)) "NA" else format(n_obs),
            if (is.na(thr)) "NA" else sprintf("%.6f", thr)),
    sprintf("  n_roots=%d  n_unit=%d  n_explosive=%d  n_near_unit=%d",
            n_s, n_unit, n_explosive, n_near),
    sprintf("  critical_root=%s  dominant_states: %s", crit_str, top_state),
    sprintf("  top_drivers: %s", drivers),
    sprintf("  action: %s",
            if (n_explosive > 0L)
              "Model has no stable solution at theta; check the parameterisation."
            else if (n_unit > 0L)
              "Unit root: use lik_init = 'auto'/'diffuse'; do not rely on unconditional moments."
            else if (n_near > 0L)
              "Persistence is long relative to T: treat moment-based diagnostics with caution; consider diffuse init; check the top drivers."
            else if (is.null(n_obs))
              "Pass n_obs to judge persistence against the sample length."
            else
              "No action needed.")
  ), collapse = "\n")

  ## ---- 6. Plots -----------------------------------------------------------
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots$roots <- .apply_meta(
      .d40_plot_roots(roots, thr, n_obs, halflife_frac, badge), meta)
    if (any(roots$modulus > 0))
      plots$half_life <- .apply_meta(
        .d40_plot_half_life(roots, n_obs, halflife_frac), meta)
    fin <- param_sens_df[is.finite(param_sens_df$d_rhomax_dparam), ,
                         drop = FALSE]
    if (nrow(fin) > 0L)
      plots$sensitivity <- .apply_meta(.d40_plot_sensitivity(fin, rho_base),
                                       meta)
  }

  .make_result(
    result = list(
      roots               = roots,
      spectral_radius     = rho_base,
      distance_to_unit    = 1 - rho_base,
      critical_eigenvalue = crit_ev,
      half_life           = hl_base,
      n_obs               = n_obs,
      halflife_frac       = halflife_frac,
      unit_tol            = unit_tol,
      modulus_threshold   = thr,
      n_unit              = n_unit,
      n_explosive         = n_explosive,
      n_near_unit         = n_near,
      stationary_p0_valid = !any(Mod(ev) > 1 - unit_tol),
      state_participation = state_part,
      param_sensitivity   = param_sens_df
    ),
    pass        = pass,
    warn        = warn,
    plots       = plots,
    summary     = summary_text,
    llm_summary = llm_text
  )
}

## Square state transition TT = ghx[state_idx, ] (n_state x n_state).
#' @noRd
.d40_state_transition <- function(dr) {
  ghx <- as.matrix(dr$ghx)
  si  <- dr$state_idx
  if (is.null(si)) {
    if (ncol(ghx) == 0L) return(matrix(0, 0, 0))
    .dynhr_abort("D40: decision rule has no `state_idx`.")
  }
  T_mat <- ghx[si, , drop = FALSE]
  if (nrow(T_mat) != ncol(T_mat))
    .dynhr_abort(sprintf("D40: state transition is %d x %d, not square.",
                         nrow(T_mat), ncol(T_mat)))
  T_mat
}

## One row per root, sorted by decreasing modulus.
#' @noRd
.d40_classify <- function(ev, n_obs, unit_tol, halflife_frac) {
  md <- Mod(ev)
  hl <- ifelse(md >= 1, Inf,
               ifelse(md == 0, 0, log(0.5) / log(pmax(md, 1e-300))))
  cls <- ifelse(md > 1 + unit_tol, "explosive",
                ifelse(md >= 1 - unit_tol, "unit", "stationary"))
  if (!is.null(n_obs))
    cls[cls == "stationary" & hl > halflife_frac * n_obs] <- "near_unit"
  ord <- order(md, decreasing = TRUE)
  out <- data.frame(re = Re(ev), im = Im(ev), modulus = md, half_life = hl,
                    class = cls, stringsAsFactors = FALSE)[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

## Participation factors |v_k w_k| / sum_j |v_j w_j| of root `k_root`.
## Left eigenvectors are the rows of W^{-1}; fall back to the eigenvectors of
## t(T) matched by eigenvalue when W is (numerically) singular.
#' @noRd
.d40_participation <- function(T_mat, eg, k_root) {
  W   <- eg$vectors
  w   <- W[, k_root]
  rc  <- if (nrow(W) > 0L) rcond(W) else 0
  if (is.finite(rc) && rc > 1e-12) {
    v <- solve(W)[k_root, ]
  } else {
    el <- eigen(t(T_mat), symmetric = FALSE)
    v  <- el$vectors[, which.min(Mod(el$values - eg$values[k_root]))[1L]]
  }
  p <- Mod(v * w)
  if (!is.finite(sum(p)) || sum(p) <= 0) return(rep(NA_real_, length(p)))
  p / sum(p)
}

## Roots in the complex plane: unit circle (solid), half-life threshold
## circle (dashed, when n_obs is known), points coloured by class.
#' @noRd
.d40_plot_roots <- function(roots, thr, n_obs, halflife_frac, badge) {
  ang <- seq(0, 2 * pi, length.out = 361L)
  circ <- data.frame(x = cos(ang), y = sin(ang), circle = "Unit circle")
  if (is.finite(thr))
    circ <- rbind(circ, data.frame(
      x = thr * cos(ang), y = thr * sin(ang),
      circle = sprintf("Near-unit threshold: half-life = %g x T (modulus %.4f)",
                       halflife_frac, thr)))
  circ$circle <- factor(circ$circle, levels = unique(circ$circle))
  n_circ <- nlevels(circ$circle)
  circ_cols <- c(dynhr_na_colour, unname(tol_vibrant["orange"]))[seq_len(n_circ)]

  cls_levels <- c("stationary", "near_unit", "unit", "explosive")
  cls_labels <- c(stationary = "Stationary", near_unit = "Near unit root",
                  unit = "Unit root", explosive = "Explosive")
  cls_cols <- c(stationary = unname(tol_vibrant["blue"]),
                near_unit  = unname(tol_vibrant["orange"]),
                unit       = unname(tol_vibrant["red"]),
                explosive  = unname(tol_vibrant["magenta"]))
  cls_shapes <- c(stationary = 16, near_unit = 17, unit = 15, explosive = 18)
  present <- cls_levels[cls_levels %in% roots$class]
  ## Draw the most severe class LAST so a unit root is not hidden under a
  ## near-unit root sitting next to it.
  pts <- roots[order(match(roots$class, cls_levels), roots$modulus), ,
               drop = FALSE]
  pts$class <- factor(pts$class, levels = cls_levels)

  lim <- max(1.1, max(roots$modulus, 0) * 1.05)
  top <- roots[1L, , drop = FALSE]
  subtitle <- if (nrow(roots) == 0L) "No state variables"
  else sprintf(paste0("%d state root(s); largest modulus %.4f, half-life %s ",
                      "periods\n%s. Badge: %s"),
               nrow(roots), top$modulus,
               if (is.finite(top$half_life)) sprintf("%.1f", top$half_life)
               else "Inf",
               if (is.null(n_obs)) "n_obs not given: no near-unit threshold"
               else sprintf("T = %g", n_obs),
               badge)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = dynhr_na_fill,
                        linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = dynhr_na_fill,
                        linewidth = 0.3) +
    ggplot2::geom_path(data = circ,
                       ggplot2::aes(x = x, y = y, linetype = circle,
                                    group = circle),
                       colour = rep(circ_cols, each = length(ang)),
                       linewidth = 0.6) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "dashed")[seq_len(n_circ)], name = NULL,
      guide = ggplot2::guide_legend(
        order = 2, override.aes = list(colour = circ_cols)))
  if (nrow(pts) > 0L) {
    p <- p +
      ggplot2::geom_point(data = pts,
                          ggplot2::aes(x = re, y = im, colour = class,
                                       shape = class),
                          size = 3) +
      ggplot2::scale_colour_manual(values = cls_cols[present],
                                   labels = cls_labels[present],
                                   breaks = present, name = NULL,
                                   guide = ggplot2::guide_legend(order = 1)) +
      ggplot2::scale_shape_manual(values = cls_shapes[present],
                                  labels = cls_labels[present],
                                  breaks = present, name = NULL,
                                  guide = ggplot2::guide_legend(order = 1))
  }
  p +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    theme_dynhr(base_size = 12) +
    ggplot2::theme(legend.box = "vertical") +
    ggplot2::labs(
      title    = "D40: Roots of the state transition ghx[state_idx, ]",
      subtitle = subtitle,
      x = "Real part", y = "Imaginary part"
    )
}

## Half-life of each non-zero root (log scale) against the half-life
## threshold halflife_frac * T; unit/explosive roots are pinned to the right
## edge (their half-life is infinite).
#' @noRd
.d40_plot_half_life <- function(roots, n_obs, halflife_frac) {
  d <- roots[roots$modulus > 0, , drop = FALSE]
  d <- utils::head(d, 25L)
  cut <- if (is.null(n_obs)) NA_real_ else halflife_frac * n_obs
  fin <- d$half_life[is.finite(d$half_life)]
  x_max <- max(c(fin, cut, 1), na.rm = TRUE) * 3
  d$x <- ifelse(is.finite(d$half_life), d$half_life, x_max)
  d$x <- pmax(d$x, 0.01)
  x_min <- min(0.5, min(d$x) / 2)
  ## Axis labels at 4 significant digits, and no "-0.000000": a root of
  ## -1.1e-17 is numerically zero and prints as "0".
  d$label <- ifelse(abs(d$im) > 1e-12,
                    paste0(.fmt_sig(d$re, zap = 1e-12),
                           ifelse(d$im >= 0, paste0("+", .fmt_sig(d$im)),
                                  .fmt_sig(d$im)), "i"),
                    .fmt_sig(d$re, zap = 1e-12))
  d$label <- make.unique(d$label, sep = " #")
  d$label <- factor(d$label, levels = rev(d$label))
  cls_cols <- c(stationary = unname(tol_vibrant["blue"]),
                near_unit  = unname(tol_vibrant["orange"]),
                unit       = unname(tol_vibrant["red"]),
                explosive  = unname(tol_vibrant["magenta"]))
  cls_labels <- c(stationary = "Stationary", near_unit = "Near unit root",
                  unit = "Unit root (half-life Inf)",
                  explosive = "Explosive (half-life Inf)")
  present <- names(cls_cols)[names(cls_cols) %in% d$class]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = label)) +
    ggplot2::geom_segment(ggplot2::aes(x = x_min, xend = x, yend = label),
                          colour = dynhr_na_fill, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = class), size = 3) +
    ggplot2::scale_colour_manual(values = cls_cols[present],
                                 labels = cls_labels[present],
                                 breaks = present, name = NULL) +
    ggplot2::scale_x_log10() +
    ggplot2::expand_limits(x = c(x_min, x_max))
  if (is.finite(cut))
    p <- p +
      ggplot2::geom_vline(xintercept = cut, linetype = "dashed",
                          colour = unname(tol_vibrant["orange"])) +
      ggplot2::annotate("text", x = cut, y = Inf, vjust = 1.2, hjust = 1.05,
                        size = 3.5, colour = unname(tol_vibrant["orange"]),
                        label = sprintf("%g x T = %.1f", halflife_frac, cut))
  p +
    theme_dynhr(base_size = 12) +
    ggplot2::labs(
      title = "D40: Half-life of each state root",
      subtitle = if (is.finite(cut))
        sprintf(paste0("Roots right of the dashed line decay too slowly for ",
                       "T = %g (near unit root)"), n_obs)
      else "Pass n_obs to draw the half-life threshold (fraction of T)",
      x = "Half-life, periods (log scale): log(0.5) / log(modulus)",
      y = "Root"
    )
}

## Bar chart of d(rho_max)/d(theta) (top 15 by magnitude).
#' @noRd
.d40_plot_sensitivity <- function(fin, rho_base) {
  fin <- utils::head(fin, 15L)
  fin$param <- factor(fin$param, levels = rev(fin$param))
  fin$direction <- ifelse(fin$d_rhomax_dparam > 0,
                          "Raises persistence", "Lowers persistence")
  ggplot2::ggplot(fin, ggplot2::aes(x = d_rhomax_dparam, y = param,
                                    fill = direction)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = dynhr_na_colour) +
    ggplot2::scale_fill_manual(
      values = c("Raises persistence" = unname(tol_vibrant["red"]),
                 "Lowers persistence" = unname(tol_vibrant["blue"])),
      name = NULL) +
    theme_dynhr(base_size = 12) +
    ggplot2::labs(
      title = "D40: Sensitivity of the spectral radius to parameters",
      subtitle = sprintf(
        "Central finite differences at theta; spectral radius = %.4f",
        rho_base),
      x = "d(spectral radius) / d(parameter)", y = NULL
    )
}
