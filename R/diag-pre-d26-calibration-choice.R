## R/diag-pre-d26-calibration-choice.R
## --------------------------------------------------------------------------
## D26 extension: WHICH parameters to calibrate and which to estimate
## (Alegre Canton 2026, arXiv:2606.25688).
##
## D26 answers "how much do my estimates move when a calibrated value is
## wrong?" for ONE calibration/estimation split. This file ranks EVERY split:
## for each admissible partition S (the estimated set) of the parameter vector
## eta = (alpha, beta), with alpha = eta_S estimated and beta = eta_{S^c} fixed,
##   D_ab   = -(G_a' W G_a)^-1 G_a' W G_b              (IFT, paper eq. 3)
##   D_gb   = Gamma_a D_ab + Gamma_b                    (chain rule, eq. 4)
##   K_S    = sqrt(|S^c|) * || D_gb Sigma_S ||_2        (Definition 1, eq. 5)
## with Sigma_S = diag(Delta_j, j in S^c) the normalisation of the calibration
## errors, and the selected partition is argmin_S K_S over the admissible set
## (Definition 2). K_S is exactly the worst-case first-order bias of the object
## of interest when every fixed parameter may be off by epsilon normalised
## units: WorstBias = epsilon * K_S + o(epsilon) (Lemma 2), and for a linear
## moment map the o(epsilon) term is zero.
##
## Everything is evaluated at ONE reference point from the Jacobians D26
## already has -- no re-estimation per partition (paper Section 5.2).
## --------------------------------------------------------------------------

#' Least-sensitive calibration/estimation partition (Alegre Canton 2026)
#'
#' Implements Algorithm 1 ("Least-Sensitive Partition") of Alegre Canton,
#' J. (2026), \emph{Choosing What to Calibrate and What to Estimate in
#' Structural Models}, arXiv:2606.25688 [econ.EM], June 2026. Given the
#' moment Jacobian \eqn{G_\eta = \partial g / \partial\eta'} of a
#' minimum-distance problem \eqn{\min_\alpha g(\alpha,\beta)' W g(\alpha,\beta)}
#' and the Jacobian \eqn{\Gamma_\eta} of an object of interest
#' \eqn{\gamma = \Gamma(\eta)}, both at a reference point \eqn{\eta^*} where the
#' model fits, every candidate estimated set \eqn{S} is scored by
#' \deqn{K_S = \sqrt{|S^c|}\;\big\| (\Gamma_\alpha D_{\alpha\beta} +
#'   \Gamma_\beta)\,\Sigma_S \big\|_2,\qquad
#'   D_{\alpha\beta} = -(G_\alpha' W G_\alpha)^{-1} G_\alpha' W G_\beta,}
#' with \eqn{\Sigma_S = \mathrm{diag}(\Delta_j)_{j \in S^c}} (the
#' \code{ranges}) and \eqn{\|\cdot\|_2} the spectral norm. \eqn{\epsilon K_S}
#' is the worst-case first-order bias of \eqn{\gamma} when each fixed
#' parameter is miscalibrated by \eqn{\epsilon} normalised units, i.e. over
#' \eqn{\|\Sigma_S^{-1}\delta\|_2 \le \epsilon\sqrt{|S^c|}} (paper Lemma 2;
#' exact for a linear moment map and linear \eqn{\Gamma}). Full estimation has
#' \eqn{K = 0}.
#'
#' \strong{Admissible set} (paper Assumption 3 and Section 5.3). A partition is
#' admissible when (i) \eqn{\Gamma_\alpha \ne 0} (the object of interest
#' depends on something that is estimated; a column counts as zero when it is
#' within 10x of its own finite-difference error), (ii) \eqn{W^{1/2}G_\alpha}
#' has full column rank \eqn{|S|} by the shared D1/D20 equilibrated,
#' finite-difference-aware rank test (\code{.ident_equilibrated_rank}), and
#' (iii) when \code{n_obs} is given, every singular value of the (raw)
#' \eqn{W^{1/2}G_\alpha} exceeds the Forneron (2024) cutoff
#' \eqn{\tau_n = \sqrt{\log n / n}} (paper eq. 18; not applied when
#' \code{n_obs} is \code{NULL}, since then \eqn{W} carries no sample-size
#' scale). \code{must_estimate}, \code{must_calibrate}, \code{min_estimated}
#' and \code{max_estimated} are the paper's researcher restrictions
#' \eqn{\mathcal R}.
#'
#' \strong{Selection.} \eqn{S^* = \arg\min K_S} over the admissible set. The
#' paper assumes a unique minimiser (Assumption 4); partitions within a
#' relative \code{1e-6} of the minimum are reported in \code{ties}, and the
#' tie is broken towards more estimated parameters (the paper prefers full
#' estimation whenever it is admissible) and then the better-identified
#' partition (larger smallest equilibrated singular value).
#'
#' \strong{Contributions} (paper Section 5.4). For each partition the squared
#' entries of the leading right singular vector of \eqn{A_S = D_{\gamma\beta}
#' \Sigma_S} give each fixed parameter's share of the worst-case perturbation
#' (paper Table 5); estimated parameters are \code{NA}.
#'
#' @param jacobian n_moment x d moment Jacobian with respect to ALL parameters
#'   (column names = parameter names), unweighted.
#' @param jacobian_2h The same Jacobian at step 2h (enables the
#'   finite-difference-aware rank tolerance), or \code{NULL}.
#' @param target_jacobian d_gamma x d Jacobian of the object of interest, with
#'   the same column names.
#' @param target_jacobian_2h Its step-2h counterpart for the zero-column test,
#'   or \code{NULL} (a column is then zero only when exactly zero).
#' @param weight n_moment x n_moment PSD weight W, or \code{NULL} (identity).
#' @param ranges Named positive normalisation widths \eqn{\Delta_j} for every
#'   parameter (paper: \eqn{\eta_{j,\max} - \eta_{j,\min}}).
#' @param current Character vector: the estimated set of the partition in use,
#'   reported for comparison (the paper's "Original" row), or \code{NULL}.
#' @param must_estimate,must_calibrate Parameter names forced into / out of
#'   the estimated set.
#' @param min_estimated,max_estimated Bounds on the number of estimated
#'   parameters.
#' @param n_obs Sample size for the \eqn{\tau_n} weak-identification cutoff,
#'   or \code{NULL}.
#' @param epsilon Normalised miscalibration size for the reported worst-case
#'   bias \eqn{\epsilon K_S} (default 0.05: 5\% of each range).
#' @param tol_rank Passed to the shared rank helper.
#' @param max_partitions Enumeration cap; above it nothing is enumerated and
#'   \code{status = "skipped"} (add restrictions to shrink the search).
#' @return list: \code{status} ("selected", "none_admissible", "skipped"),
#'   \code{selected} (estimated, calibrated, K, worst_bias, contributions),
#'   \code{ties}, \code{current}, \code{partitions} (data frame, one row per
#'   candidate, sorted by K), \code{contributions} (candidate x parameter),
#'   \code{ranges}, \code{epsilon}, \code{tau_n}, \code{n_candidates},
#'   \code{message}.
#' @references
#'   Alegre Canton, J. (2026). Choosing what to calibrate and what to estimate
#'   in structural models. arXiv:2606.25688 [econ.EM], June 2026.
#'
#'   Forneron, J.-J. (2024). Detecting identification failure in moment
#'   condition models. \emph{Journal of Econometrics}, 238(1), 105552.
#' @noRd
d26_calibration_choice <- function(jacobian,
                                   jacobian_2h = NULL,
                                   target_jacobian,
                                   target_jacobian_2h = NULL,
                                   weight = NULL,
                                   ranges,
                                   current = NULL,
                                   must_estimate = character(0),
                                   must_calibrate = character(0),
                                   min_estimated = 1L,
                                   max_estimated = Inf,
                                   n_obs = NULL,
                                   epsilon = 0.05,
                                   tol_rank = NULL,
                                   max_partitions = 4096L) {
  G <- as.matrix(jacobian)
  pn <- colnames(G)
  d <- ncol(G)
  if (is.null(pn) || anyDuplicated(pn) || any(!nzchar(pn)))
    .dynhr_abort("calibration choice: `jacobian` needs unique column names.")
  if (!all(is.finite(G)))
    .dynhr_abort("calibration choice: `jacobian` must be finite.")
  mom_names <- rownames(G) %||% paste0("m_", seq_len(nrow(G)))
  G2 <- if (!is.null(jacobian_2h)) as.matrix(jacobian_2h)[, pn, drop = FALSE]
  Gam <- as.matrix(target_jacobian)
  if (is.null(colnames(Gam)) || !setequal(colnames(Gam), pn) || ncol(Gam) != d)
    .dynhr_abort("calibration choice: columns of `target_jacobian` must be the ",
                 "parameters of `jacobian`.")
  Gam <- Gam[, pn, drop = FALSE]
  if (!all(is.finite(Gam)))
    .dynhr_abort("calibration choice: `target_jacobian` must be finite.")
  target_names <- rownames(Gam) %||% paste0("gamma_", seq_len(nrow(Gam)))
  if (!is.numeric(ranges) || is.null(names(ranges)) || !all(pn %in% names(ranges)))
    .dynhr_abort("calibration choice: `ranges` must be named and cover every parameter.")
  ranges <- ranges[pn]
  if (!all(is.finite(ranges)) || any(ranges <= 0))
    .dynhr_abort("calibration choice: `ranges` must be finite and > 0.")
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0)
    .dynhr_abort("calibration choice: `epsilon` must be a positive number.")
  tau_n <- NA_real_
  if (!is.null(n_obs)) {
    if (!is.numeric(n_obs) || length(n_obs) != 1L || !is.finite(n_obs) || n_obs <= 1)
      .dynhr_abort("calibration choice: `n_obs` must be a number > 1.")
    tau_n <- sqrt(log(n_obs) / n_obs)
  }
  must_estimate <- as.character(must_estimate)
  must_calibrate <- as.character(must_calibrate)
  bad <- setdiff(c(must_estimate, must_calibrate, current), pn)
  if (length(bad))
    .dynhr_abort("calibration choice: unknown parameter(s): ",
                 paste(bad, collapse = ", "), ".")
  both <- intersect(must_estimate, must_calibrate)
  if (length(both))
    .dynhr_abort("calibration choice: both must-estimate and must-calibrate: ",
                 paste(both, collapse = ", "), ".")

  Wh <- .d26_weight_sqrt(weight, mom_names)
  weigh <- function(J) if (is.null(Wh)) J else Wh %*% J
  Gw <- weigh(G)
  colnames(Gw) <- pn
  Gw2 <- if (!is.null(G2) && all(is.finite(G2))) {
    x <- weigh(G2)
    colnames(x) <- pn
    x
  }

  ## Assumption 3.i: a target column is zero when it is within 10x of its own
  ## finite-difference error (exactly zero when no 2h Jacobian is given).
  gam_noise <- if (!is.null(target_jacobian_2h)) {
    Gam2 <- as.matrix(target_jacobian_2h)[, pn, drop = FALSE]
    apply(abs(Gam - Gam2), 2, max)
  } else rep(0, d)
  gam_live <- apply(abs(Gam), 2, max) > 10 * gam_noise &
              apply(abs(Gam), 2, max) > 0

  ## ---- evaluate one estimated set --------------------------------------------
  eval_S <- function(S) {
    Sc <- setdiff(pn, S)
    out <- list(estimated = S, calibrated = Sc, rank = NA_integer_,
                rank_ok = FALSE, sv_min = NA_real_, sv_min_equilibrated = NA_real_,
                weak = NA, nontrivial = any(gam_live[S]), K = NA_real_,
                contributions = stats::setNames(rep(NA_real_, d), pn))
    Ja <- Gw[, S, drop = FALSE]
    rk <- .ident_equilibrated_rank(Ja, if (!is.null(Gw2)) Gw2[, S, drop = FALSE],
                                   tol_rank = tol_rank)
    out$rank <- rk$rank
    out$rank_ok <- rk$rank == length(S)
    out$sv_min_equilibrated <- min(rk$singular_values)
    sv_raw <- svd(Ja, nu = 0, nv = 0)$d
    out$sv_min <- if (length(sv_raw) < length(S)) 0 else min(sv_raw)
    out$weak <- if (is.na(tau_n)) NA else out$sv_min <= tau_n
    if (!out$rank_ok) return(out)
    if (!length(Sc)) {
      out$K <- 0
      return(out)
    }
    cs <- sqrt(colSums(Ja^2))
    D_ab <- -qr.coef(qr(sweep(Ja, 2, cs, "/")), Gw[, Sc, drop = FALSE]) / cs
    D_gb <- Gam[, S, drop = FALSE] %*% D_ab + Gam[, Sc, drop = FALSE]
    A <- sweep(D_gb, 2, ranges[Sc], "*")
    sv <- svd(A, nu = 0, nv = 1L)
    out$K <- sqrt(length(Sc)) * sv$d[1L]
    v <- if (sv$d[1L] > 0) sv$v[, 1L]^2 else rep(NA_real_, length(Sc))
    out$contributions[Sc] <- v
    out
  }

  ## ---- enumerate the restricted candidate set (paper: S in script-S n R) ------
  free <- setdiff(pn, c(must_estimate, must_calibrate))
  n_cand <- 2^length(free)
  base <- list(status = "skipped", selected = NULL, ties = list(),
               current = NULL, partitions = NULL, contributions = NULL,
               ranges = ranges, epsilon = epsilon, tau_n = tau_n,
               n_obs = n_obs, target = target_names, n_candidates = n_cand,
               message = "")
  if (n_cand > max_partitions) {
    base$message <- sprintf(
      "%d free parameters give %.0f candidate partitions (> max_partitions = %d); restrict the search with must_estimate / must_calibrate / min_estimated / max_estimated.",
      length(free), n_cand, as.integer(max_partitions))
    return(base)
  }
  cand <- lapply(seq_len(n_cand) - 1L, function(mask) {
    pick <- if (length(free)) bitwAnd(mask, as.integer(2^(seq_along(free) - 1L))) > 0 else logical(0)
    pn[pn %in% c(must_estimate, free[pick])]
  })
  sz <- lengths(cand)
  cand <- cand[sz >= max(1L, min_estimated) & sz <= max_estimated]
  cur_key <- if (!is.null(current)) paste(pn[pn %in% current], collapse = "\r")
  keys <- vapply(cand, paste, "", collapse = "\r")
  evals <- lapply(cand, eval_S)
  cur_in <- !is.null(cur_key) && cur_key %in% keys
  if (!is.null(current) && !cur_in) {
    evals <- c(evals, list(eval_S(pn[pn %in% current])))
    keys <- c(keys, cur_key)
  }
  restricted <- c(rep(TRUE, length(cand)), rep(FALSE, length(evals) - length(cand)))

  pt <- data.frame(
    estimated = vapply(evals, function(e) paste(e$estimated, collapse = ", "), ""),
    calibrated = vapply(evals, function(e) paste(e$calibrated, collapse = ", "), ""),
    n_estimated = vapply(evals, function(e) length(e$estimated), 1L),
    rank = vapply(evals, function(e) as.integer(e$rank), 1L),
    rank_ok = vapply(evals, function(e) e$rank_ok, TRUE),
    sv_min = vapply(evals, function(e) e$sv_min, 1),
    weak = vapply(evals, function(e) as.logical(e$weak), NA),
    nontrivial = vapply(evals, function(e) e$nontrivial, TRUE),
    K = vapply(evals, function(e) e$K, 1),
    stringsAsFactors = FALSE)
  pt$worst_bias <- epsilon * pt$K
  pt$admissible <- restricted & pt$rank_ok & pt$nontrivial &
                   (is.na(pt$weak) | !pt$weak)
  pt$current <- if (is.null(cur_key)) rep(FALSE, length(keys)) else keys == cur_key
  contrib <- do.call(rbind, lapply(evals, `[[`, "contributions"))
  if (is.null(contrib)) contrib <- matrix(NA_real_, 0L, d, dimnames = list(NULL, pn))
  sv_eq <- vapply(evals, function(e) e$sv_min_equilibrated, 1)

  ## ---- order: admissible by K, then the rest --------------------------------
  ord <- order(!pt$admissible, pt$K, -pt$n_estimated, -sv_eq, na.last = TRUE)
  pt <- pt[ord, , drop = FALSE]
  contrib <- contrib[ord, , drop = FALSE]
  evals <- evals[ord]
  sv_eq <- sv_eq[ord]
  pt$position <- ifelse(pt$admissible, cumsum(pt$admissible), NA_integer_)
  rownames(pt) <- NULL
  rownames(contrib) <- pt$estimated
  base$partitions <- pt
  base$contributions <- contrib
  if (!is.null(current)) {
    ic <- which(pt$current)[1L]
    base$current <- list(estimated = evals[[ic]]$estimated,
                         calibrated = evals[[ic]]$calibrated,
                         K = pt$K[ic], worst_bias = pt$worst_bias[ic],
                         admissible = pt$admissible[ic],
                         position = pt$position[ic],
                         contributions = contrib[ic, evals[[ic]]$calibrated])
  }
  adm <- which(pt$admissible)
  if (!length(adm)) {
    base$status <- "none_admissible"
    base$message <- "No admissible partition: no candidate estimated set is locally identified (and, with n_obs, clear of the weak-identification cutoff)."
    return(base)
  }
  K_adm <- pt$K[adm]
  tie_tol <- 1e-6 * max(K_adm)
  tie_rows <- adm[K_adm <= min(K_adm) + tie_tol]
  tie_rows <- tie_rows[order(-pt$n_estimated[tie_rows], -sv_eq[tie_rows])]
  sel <- tie_rows[1L]
  base$status <- "selected"
  base$selected <- list(estimated = evals[[sel]]$estimated,
                        calibrated = evals[[sel]]$calibrated,
                        K = pt$K[sel], worst_bias = pt$worst_bias[sel],
                        contributions = contrib[sel, evals[[sel]]$calibrated])
  base$ties <- lapply(tie_rows[-1L], function(i) evals[[i]]$estimated)
  base$message <- sprintf(
    "least-sensitive partition (Alegre Canton 2026) estimates {%s}, calibrates {%s}: K = %.3g (worst-case bias %.3g at epsilon = %g), %d admissible partition(s)%s.",
    pt$estimated[sel],
    if (nzchar(pt$calibrated[sel])) pt$calibrated[sel] else "nothing: full estimation",
    pt$K[sel], pt$worst_bias[sel], epsilon, length(adm),
    if (!is.null(base$current))
      sprintf("; current partition K = %s%s", format(signif(base$current$K, 3)),
              if (isTRUE(base$current$admissible))
                sprintf(" (ranked %d of %d)", base$current$position, length(adm))
              else " (not admissible)")
    else "")
  base
}


## D26 -> d26_calibration_choice(): resolve the `calibration_choice` options
## (target, ranges, restrictions) against D26's parameters and hand over the
## Jacobians D26 already computed. Returns NULL when the option is switched off.
.d26_choice_from_d26 <- function(opts, theta_c, theta_e, J_e, J_e2, J_c, J_c2,
                                 weight, tol_rank, eps, all_inert) {
  if (is.null(opts) || isFALSE(opts)) return(NULL)
  if (isTRUE(opts)) opts <- list()
  known <- c("target", "ranges", "epsilon", "n_obs", "must_estimate",
             "must_calibrate", "min_estimated", "max_estimated", "max_partitions")
  if (!is.list(opts) || (length(opts) && is.null(names(opts))) ||
      length(setdiff(names(opts), known)))
    .dynhr_abort("D26: `calibration_choice` must be NULL/FALSE or a named list ",
                 "with elements from: ", paste(known, collapse = ", "), ".")
  eta <- c(theta_c, theta_e)
  pn <- names(eta)
  if (all_inert)
    return(list(status = "inert", selected = NULL, message = paste0(
      "not assessable: no calibrated parameter moves any moment (see the D26 ",
      "status), so every calibration would look harmless.")))

  ## Normalisation widths Delta_j (paper Section 5.4). Default: |eta_j| (1 at
  ## zero), i.e. a calibration error is measured relative to the value -- the
  ## same convention as D26's elasticity and perturbation grid.
  delta <- stats::setNames(ifelse(eta == 0, 1, abs(eta)), pn)
  src <- stats::setNames(rep("relative (|value|)", length(pn)), pn)
  rg <- opts$ranges
  if (!is.null(rg)) {
    if (is.matrix(rg) || is.data.frame(rg)) {
      rg <- as.matrix(rg)
      if (ncol(rg) != 2L || is.null(rownames(rg)) || !is.numeric(rg))
        .dynhr_abort("D26: a matrix `calibration_choice$ranges` needs 2 numeric ",
                     "columns (min, max) and parameter row names.")
      rg <- stats::setNames(rg[, 2] - rg[, 1], rownames(rg))
    }
    if (!is.numeric(rg) || is.null(names(rg)))
      .dynhr_abort("D26: `calibration_choice$ranges` must be named widths or a ",
                   "(min, max) matrix with parameter row names.")
    bad <- setdiff(names(rg), pn)
    if (length(bad))
      .dynhr_abort("D26: `calibration_choice$ranges` names unknown parameter(s): ",
                   paste(bad, collapse = ", "), ".")
    if (!all(is.finite(rg)) || any(rg <= 0))
      .dynhr_abort("D26: `calibration_choice$ranges` must be finite and > 0.")
    delta[names(rg)] <- rg
    src[names(rg)] <- "user"
  }

  ## Object of interest gamma = Gamma(eta). Default: the currently estimated
  ## parameters; parameter targets are measured relative to their value, so
  ## with the default ranges K for D26's own partition is sqrt(n_cal) times
  ## the spectral norm of D26's elasticity matrix. A parameter that is the
  ## object of interest is always estimated (paper, remark on Assumption 3.i).
  tg <- opts$target %||% names(theta_e)
  forced <- character(0)
  Gam2 <- NULL
  if (is.function(tg)) {
    g0 <- tg(eta)
    if (!is.numeric(g0) || !length(g0) || !all(is.finite(g0)))
      .dynhr_abort("D26: `calibration_choice$target` returned no / non-finite ",
                   "values at the baseline.")
    tn <- names(g0) %||% paste0("gamma_", seq_along(g0))
    tfn <- function(x) tg(stats::setNames(x, pn))
    Gam <- .numerical_jacobian(tfn, eta, eps = eps)
    Gam2 <- .numerical_jacobian(tfn, eta, eps = 2 * eps)
    if (!all(is.finite(Gam)))
      .dynhr_abort("D26: the Jacobian of `calibration_choice$target` is non-finite.")
    dimnames(Gam) <- list(tn, pn)
    if (all(is.finite(Gam2))) dimnames(Gam2) <- list(tn, pn) else Gam2 <- NULL
  } else if (is.character(tg) && length(tg) && all(tg %in% pn)) {
    tg <- unique(tg)
    Gam <- matrix(0, length(tg), length(pn), dimnames = list(tg, pn))
    Gam[cbind(seq_along(tg), match(tg, pn))] <-
      1 / ifelse(eta[tg] == 0, 1, abs(eta[tg]))
    forced <- tg
  } else {
    .dynhr_abort("D26: `calibration_choice$target` must be parameter names or a ",
                 "function(eta) returning the object(s) of interest.")
  }

  G <- cbind(J_c, J_e)[, pn, drop = FALSE]
  G2 <- if (all(is.finite(J_c2)) && all(is.finite(J_e2)))
    cbind(J_c2, J_e2)[, pn, drop = FALSE]
  out <- d26_calibration_choice(
    jacobian = G, jacobian_2h = G2,
    target_jacobian = Gam, target_jacobian_2h = Gam2,
    weight = weight, ranges = delta, current = names(theta_e),
    must_estimate = union(forced, as.character(opts$must_estimate)),
    must_calibrate = as.character(opts$must_calibrate),
    min_estimated = opts$min_estimated %||% 1L,
    max_estimated = opts$max_estimated %||% Inf,
    n_obs = opts$n_obs, epsilon = opts$epsilon %||% 0.05,
    tol_rank = tol_rank, max_partitions = opts$max_partitions %||% 4096L)
  out$range_source <- src
  out
}
