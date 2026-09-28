## R/diag-pre-d30-arbitrary-precision.R
## --------------------------------------------------------------------------
## Phase H: D30 -- Arbitrary-Precision Rank Checks
##
## Re-computes the local identification rank (D1) with the WHOLE pipeline in
## multiple precision: the moment function is evaluated on mpfr parameters,
## the Jacobian is a central difference at an mpfr step, and the SVD is a
## one-sided Jacobi SVD in mpfr.  A double-precision finite-difference
## Jacobian carries truncation + rounding error of ~1e-10 relative, so a
## genuine singular value of ~1e-9 is indistinguishable from a structural
## zero in double; at 128 bits the finite-difference error is ~1e-26 and the
## two cases separate by many orders of magnitude.
##
## Converting a double Jacobian to mpfr gains NOTHING (the error is already
## baked in), so the high-precision path requires a moment function that is
## generic over Rmpfr numbers (`hp_solve_fn`).  Its precision is verified
## with a sub-double-epsilon probe before it is trusted.  Without such a
## function, or without Rmpfr, D30 reports the double baseline only and
## returns pass = NA -- it never claims a high-precision verdict it did not
## compute.
##
## Both precisions use the same rank rule as D1
## (`.ident_equilibrated_rank()`: row max-abs / column 2-norm
## equilibration, tolerance = max(machine, 10 * ||Je(h) - Je(2h)||_2)),
## with the machine epsilon of the respective precision.
## --------------------------------------------------------------------------

#' D30. Arbitrary-Precision Rank Checks
#'
#' Re-computes the D1 local-identification rank in multiple precision to
#' separate structural rank deficiency from finite-difference / rounding
#' noise.
#'
#' \enumerate{
#'   \item Double baseline: central-difference Jacobian of
#'     \code{model_solve_fn} at steps \code{eps} and \code{2 * eps}, ranked
#'     with \code{.ident_equilibrated_rank()} (identical rule to D1).
#'   \item High precision (Rmpfr): \code{hp_solve_fn} is evaluated on
#'     \code{Rmpfr::mpfr(theta, prec_bits)}; the Jacobian is a central
#'     difference at relative step \eqn{2^{-\lfloor prec/3 \rfloor}} (and
#'     twice that, for the noise estimate), equilibrated the same way, and
#'     its singular values come from a one-sided Jacobi SVD carried out in
#'     mpfr.  The tolerance is
#'     \eqn{\max(\max(m,n)\,\sigma_1 2^{-prec},\; 10\,\|J_e(h)-J_e(2h)\|_2)}.
#'   \item Verdict: directions below the high-precision tolerance are
#'     structural; directions the double check calls unidentified but the
#'     high-precision check resolves are numerical artefacts.
#' }
#'
#' \code{hp_solve_fn} must propagate mpfr numbers end to end (plain R
#' arithmetic and \code{exp}/\code{log}/\code{sqrt}/... do; \code{solve()},
#' \code{eigen()}, compiled code and \code{as.numeric()} do not).  Before
#' use it is checked to (a) return an mpfr vector of the right length,
#' (b) agree with \code{model_solve_fn} in double, and (c) respond to a
#' parameter perturbation of \eqn{2^{-3 prec/4}} (far below double epsilon)
#' with the slope of its own Jacobian -- a function that routes through
#' double fails (c).  If any check fails, or \code{hp_solve_fn} is
#' \code{NULL}, or Rmpfr is not installed, or \code{backend = "base"}, only
#' the double baseline is reported and \code{pass = NA}.
#'
#' \strong{This diagnostic is deliberately opt-in, and it is not the
#' package's answer to finite-difference noise.} That answer is
#' \code{.ident_equilibrated_rank()} -- the FD-noise-aware equilibrated rank
#' test shared by D1, D20 and D28, which compares the singular values against
#' a noise floor estimated from a second Jacobian at step \code{2 * eps}
#' instead of trying to out-precision the noise. That is what the literature
#' does too: Iskrev (2010) builds his rank test on \emph{analytic} derivatives
#' of the moment mapping, sidestepping FD error entirely, and Dynare's
#' \code{identification} command applies a relative singular-value cutoff
#' (\code{tol_rank}) to an analytically differentiated Jacobian. No
#' identification paper (Iskrev 2010, Ratto & Iskrev 2011, Mutschler 2015,
#' Qu & Tkachenko 2012/2017, Komunjer & Ng 2011, Kociecki 2018) and no
#' toolbox surveyed (Dynare, IRIS, MacroModelling.jl, RISE) proposes or ships
#' a multiple-precision Lyapunov/QZ path for identification, and mpfr
#' arithmetic buys nothing unless the \emph{input} derivatives are already
#' high precision -- which is exactly why a genuine mpfr-generic
#' \code{hp_solve_fn} is required rather than an internal double-to-mpfr
#' conversion. Use D30 when you have such a function and want to confirm that
#' a borderline direction is structural; otherwise read D1/D20.
#'
#' @param model_solve_fn  Function: theta -> named numeric vector of moments.
#' @param theta           Named numeric vector of parameter values.
#' @param param_names     Character vector of parameter names.
#' @param moment_names    Character vector of moment names.
#' @param hp_solve_fn     Optional mpfr-generic version of
#'   \code{model_solve_fn} (often the same function). \code{NULL} disables
#'   the high-precision check.
#' @param prec_bits       MPFR precision in bits (integer >= 64, default 128).
#' @param backend         \code{"auto"} / \code{"Rmpfr"} (high precision via
#'   Rmpfr; \code{"Rmpfr"} warns if it is not installed) or \code{"base"}
#'   (double baseline only, \code{pass = NA}).
#' @param eps             Double finite-difference step (default 1e-5).
#' @param weak_rel        Relative singular-value level below which an
#'   identified direction is reported as weak (default 1e-3, as D1).
#' @param verbose         Print progress messages.
#' @param meta            Optional plot metadata (\code{diag_meta()}).
#'
#' @return A \code{dynhr_diagnostic}; \code{result} holds
#'   \code{double_rank}, \code{double_tol}, \code{double_fd_noise},
#'   \code{high_prec_rank}, \code{high_prec_tol}, \code{high_prec_fd_noise},
#'   \code{sv_comparison} (equilibrated singular values at both
#'   precisions), \code{structural_dirs}, \code{artefact_dirs},
#'   \code{weak_dirs}, \code{null_space} (high-precision right singular
#'   vectors of the structural directions, double), \code{rank_consistent},
#'   \code{backend_used} (\code{"Rmpfr"} or \code{"none"}),
#'   \code{hp_status} (why the high-precision check did or did not run),
#'   \code{prec_bits}, \code{param_names}.  \code{pass} is \code{TRUE} when
#'   the high-precision rank is full, \code{FALSE} when it is deficient, and
#'   \code{NA} when no high-precision check was performed.
#' @noRd
d30_arbitrary_precision_rank <- function(model_solve_fn,
                                          theta,
                                          param_names = NULL,
                                          moment_names = NULL,
                                          hp_solve_fn = NULL,
                                          prec_bits = 128L,
                                          backend = c("auto", "Rmpfr", "base"),
                                          eps = 1e-5,
                                          weak_rel = 1e-3,
                                          verbose = FALSE,
                                          meta = NULL) {
  backend <- match.arg(backend)

  if (is.null(theta) || !is.function(model_solve_fn)) {
    return(.make_result(
      pass    = NA,
      summary = "D30 Arbitrary-Precision Rank: theta or model_solve_fn is NULL."
    ))
  }
  if (!is.numeric(prec_bits) || length(prec_bits) != 1L || !is.finite(prec_bits) ||
      prec_bits < 64) {
    .dynhr_abort("d30: `prec_bits` must be a single number >= 64.")
  }
  prec_bits <- as.integer(prec_bits)
  theta <- unlist(theta)
  n_par <- length(theta)
  if (is.null(param_names) || length(param_names) != n_par) {
    param_names <- names(theta) %||% paste0("theta_", seq_len(n_par))
  }
  names(theta) <- param_names

  f0 <- model_solve_fn(theta)
  if (is.null(f0) || length(f0) == 0L || !all(is.finite(f0))) {
    return(.make_result(
      pass    = NA,
      summary = "D30 Arbitrary-Precision Rank: model_solve_fn returned NULL, empty or non-finite moments."
    ))
  }
  if (is.null(moment_names) || length(moment_names) != length(f0)) {
    moment_names <- names(f0) %||% paste0("m_", seq_along(f0))
  }

  # ---- 1. Double baseline (same rule as D1) ----
  if (verbose) .dynhr_cat("[d30] Double-precision baseline...\n")
  J  <- .numerical_jacobian(model_solve_fn, theta, eps = eps)
  J2 <- .numerical_jacobian(model_solve_fn, theta, eps = 2 * eps)
  if (!all(is.finite(J)) || !all(is.finite(J2))) {
    return(.make_result(
      pass    = NA,
      summary = "D30 Arbitrary-Precision Rank: double Jacobian is non-finite; rank undetermined."
    ))
  }
  dimnames(J) <- list(moment_names, param_names)
  dimnames(J2) <- dimnames(J)
  rk_d <- .ident_equilibrated_rank(J, J2, weak_rel = weak_rel)
  sv_d <- unname(rk_d$singular_values)
  if (verbose) .dynhr_cat(sprintf("[d30] double: rank %d / %d (tol %.2e, %s)\n",
                                  rk_d$rank, n_par, rk_d$tol, rk_d$tol_source))

  # ---- 2. High-precision path: availability ----
  hp <- NULL
  hp_status <- NULL
  if (backend == "base") {
    hp_status <- "backend = 'base' requested: double baseline only"
  } else if (!requireNamespace("Rmpfr", quietly = TRUE)) {
    hp_status <- "Rmpfr is not installed (install.packages('Rmpfr'))"
    if (backend == "Rmpfr") .dynhr_warn(paste0("d30: backend = 'Rmpfr' requested but ", hp_status, "."))
  } else if (!is.function(hp_solve_fn)) {
    hp_status <- paste0("no mpfr-generic hp_solve_fn supplied (a Jacobian computed ",
                        "in double gains nothing from mpfr)")
  } else {
    if (verbose) .dynhr_cat(sprintf("[d30] Rmpfr path at %d bits...\n", prec_bits))
    hp <- .d30_hp_rank(hp_solve_fn, theta, f0, J, prec_bits, weak_rel)
    hp_status <- hp$status
    if (!isTRUE(hp$ok)) {
      .dynhr_warn(paste0("d30: high-precision check not performed: ", hp_status, "."))
      hp <- NULL
    }
  }
  hp_ran <- !is.null(hp)

  # ---- 3. Comparison ----
  rank_high <- if (hp_ran) hp$rank else NA_integer_
  sv_h <- if (hp_ran) hp$singular_values else rep(NA_real_, n_par)
  sv_comparison <- data.frame(
    index      = seq_len(n_par),
    double     = sv_d,
    high_prec  = sv_h,
    double_class = ifelse(sv_d <= rk_d$tol, "Unidentified",
                   ifelse(sv_d < rk_d$weak_threshold, "Weak", "Identified")),
    high_prec_class = if (hp_ran) hp$sv_class else NA_character_,
    stringsAsFactors = FALSE
  )
  structural <- if (hp_ran) which(hp$sv_class == "Unidentified") else integer(0)
  weak_dirs  <- if (hp_ran) which(hp$sv_class == "Weak") else integer(0)
  artefact   <- if (hp_ran) which(sv_comparison$double_class == "Unidentified" &
                                  hp$sv_class != "Unidentified") else integer(0)
  rank_consistent <- if (hp_ran) rk_d$rank == rank_high else NA
  pass_val <- if (hp_ran) rank_high == n_par else NA
  null_space <- if (hp_ran) hp$V[, structural, drop = FALSE] else NULL
  structural_params <- if (hp_ran && length(structural))
    unique(unlist(lapply(structural, function(k) {
      v <- abs(hp$V[, k]); param_names[v >= 0.5 * max(v)] }))) else character(0)

  # ---- 4. Plot ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots$sv_comparison <- .apply_meta(
      .d30_plot(sv_comparison, rk_d, hp, n_par, prec_bits, hp_status), meta)
  }

  # ---- 5. Result + text ----
  result <- list(
    double_rank        = rk_d$rank,
    double_tol         = rk_d$tol,
    double_fd_noise    = rk_d$fd_noise,
    high_prec_rank     = rank_high,
    high_prec_tol      = if (hp_ran) hp$tol else NA_real_,
    high_prec_fd_noise = if (hp_ran) hp$fd_noise else NA_real_,
    rank_consistent    = rank_consistent,
    sv_comparison      = sv_comparison,
    structural_dirs    = structural,
    structural_params  = structural_params,
    artefact_dirs      = artefact,
    weak_dirs          = weak_dirs,
    null_space         = null_space,
    backend_used       = if (hp_ran) "Rmpfr" else "none",
    hp_status          = hp_status,
    prec_bits          = prec_bits,
    param_names        = param_names
  )

  if (!hp_ran) {
    summary_str <- sprintf(paste0(
      "D30 Arbitrary-Precision Rank: double rank=%d/%d (tol=%.2e). ",
      "High-precision check NOT performed: %s."), rk_d$rank, n_par, rk_d$tol, hp_status)
    llm_str <- sprintf(paste0(
      "[INFO] D30 Arbitrary-Precision Rank double_rank=%d/%d high_prec=not_run\n",
      "  reason: %s\n",
      "  action: pass an mpfr-generic hp_solve_fn (d30_hp_solve_fn) with Rmpfr installed"),
      rk_d$rank, n_par, hp_status)
  } else {
    verdict <- if (isTRUE(pass_val)) {
      if (length(artefact)) sprintf(
        "Identified at %d bits: %d double-precision deficient direction(s) are finite-difference/rounding artefacts.",
        prec_bits, length(artefact))
      else sprintf("Full rank confirmed at %d bits.", prec_bits)
    } else {
      sprintf("STRUCTURAL rank deficiency confirmed at %d bits (%d direction(s); involving %s).",
              prec_bits, length(structural), paste(structural_params, collapse = ", "))
    }
    if (length(weak_dirs)) verdict <- paste0(verdict, sprintf(
      " %d weak but real direction(s) (sv < %.0e * sv_max).", length(weak_dirs), weak_rel))
    summary_str <- sprintf(
      "D30 Arbitrary-Precision Rank: double rank=%d/%d (tol=%.2e) vs %d-bit rank=%d/%d (tol=%.2e). %s",
      rk_d$rank, n_par, rk_d$tol, prec_bits, rank_high, n_par, hp$tol, verdict)
    llm_str <- sprintf(paste0(
      "[%s] D30 Arbitrary-Precision Rank double_rank=%d high_prec_rank=%d n_par=%d prec_bits=%d ",
      "structural=%d artefacts=%d weak=%d min_sv_hp=%.3e hp_tol=%.3e\n  action: %s"),
      if (isTRUE(pass_val)) "PASS" else "FAIL",
      rk_d$rank, rank_high, n_par, prec_bits, length(structural), length(artefact),
      length(weak_dirs), min(sv_h), hp$tol,
      if (isTRUE(pass_val)) "none (identification holds in high precision)"
      else paste0("reparameterise or calibrate one of: ", paste(structural_params, collapse = ", ")))
  }

  .make_result(result = result, pass = pass_val, plots = plots,
               summary = summary_str, llm_summary = llm_str)
}


# ==========================================================================
# High-precision helpers
# ==========================================================================

#' Rank of the identification Jacobian computed end-to-end in mpfr
#'
#' @param hp_solve_fn mpfr-generic moment function.
#' @param theta named double parameters; @param f0 double moments;
#' @param J double Jacobian (only used by the precision probe).
#' @return list(ok, status, rank, tol, fd_noise, singular_values, sv_class, V)
#' @noRd
.d30_hp_rank <- function(hp_solve_fn, theta, f0, J, prec_bits, weak_rel) {
  n_par <- length(theta); n_mom <- length(f0)
  fail <- function(msg) list(ok = FALSE, status = msg)
  mp <- function(x) Rmpfr::mpfr(x, precBits = prec_bits)
  th <- mp(unname(theta))
  eval_hp <- function(x) {
    out <- hp_solve_fn(x)
    if (!methods::is(out, "mpfr") || length(out) != n_mom) return(NULL)
    if (!all(is.finite(Rmpfr::asNumeric(out)))) return(NULL)
    out
  }

  g0 <- eval_hp(th)
  if (is.null(g0)) return(fail(sprintf(
    "hp_solve_fn did not return a finite mpfr vector of length %d", n_mom)))
  g0d <- Rmpfr::asNumeric(g0)
  if (max(abs(g0d - f0) / pmax(1, abs(f0))) > 1e-8)
    return(fail("hp_solve_fn disagrees with model_solve_fn in double (> 1e-8 relative)"))

  bump <- function(j, h) { x <- th; x[j] <- x[j] + h; x }
  scale_j <- pmax(1, abs(unname(theta)))
  h_rel <- mp(2)^(-(prec_bits %/% 3L))
  cols1 <- vector("list", n_par); cols2 <- vector("list", n_par)
  for (j in seq_len(n_par)) {
    h <- h_rel * scale_j[j]
    fp <- eval_hp(bump(j, h)); fm <- eval_hp(bump(j, -h))
    fp2 <- eval_hp(bump(j, 2 * h)); fm2 <- eval_hp(bump(j, -2 * h))
    if (is.null(fp) || is.null(fm) || is.null(fp2) || is.null(fm2))
      return(fail(sprintf("hp_solve_fn failed near theta[%d]", j)))
    cols1[[j]] <- (fp - fm) / (2 * h)
    cols2[[j]] <- (fp2 - fm2) / (4 * h)

    # Precision probe: a step far below double epsilon must reproduce the
    # column slope.  A function that routes through double returns 0 here.
    d <- mp(2)^(-((3L * prec_bits) %/% 4L)) * scale_j[j]
    fd <- eval_hp(bump(j, d))
    if (is.null(fd)) return(fail(sprintf("hp_solve_fn failed near theta[%d]", j)))
    slope <- Rmpfr::asNumeric((fd - g0) / d)
    colj <- Rmpfr::asNumeric(cols1[[j]])
    if (max(abs(colj)) > 0 && max(abs(slope - colj)) > 1e-6 * max(abs(colj)))
      return(fail(sprintf(paste0(
        "hp_solve_fn is not precise below double epsilon (probe on theta[%d]: ",
        "slope %.3g vs Jacobian %.3g) -- it probably converts to double internally"),
        j, max(abs(slope)), max(abs(colj)))))
  }

  # Equilibration identical to .ident_equilibrated_rank(), in mpfr.
  row_scale <- abs(cols1[[1]]); row_noise <- abs(cols1[[1]] - cols2[[1]])
  for (j in seq_len(n_par)[-1]) {
    row_scale <- Rmpfr::pmax(row_scale, abs(cols1[[j]]))
    row_noise <- Rmpfr::pmax(row_noise, abs(cols1[[j]] - cols2[[j]]))
  }
  keep <- as.logical(row_scale > 10 * row_noise & row_scale > 0)
  if (!any(keep)) return(list(ok = TRUE, status = "ok", rank = 0L, tol = 0,
                              fd_noise = NA_real_, singular_values = rep(0, n_par),
                              sv_class = rep("Unidentified", n_par), V = diag(n_par)))
  rs <- row_scale[keep]
  e1 <- lapply(cols1, function(cl) cl[keep] / rs)
  e2 <- lapply(cols2, function(cl) cl[keep] / rs)
  for (j in seq_len(n_par)) {
    cs <- sqrt(sum(e1[[j]]^2))
    if (Rmpfr::asNumeric(cs) > 0) { e1[[j]] <- e1[[j]] / cs; e2[[j]] <- e2[[j]] / cs }
  }
  s1 <- .d30_jacobi_svd(e1, prec_bits)
  noise <- .d30_jacobi_svd(Map(`-`, e1, e2), prec_bits)$d[1]
  sv_max <- s1$d[1]
  tol_machine <- max(sum(keep), n_par) * sv_max * 2^(-prec_bits)
  tol <- max(tol_machine, 10 * noise)
  weak_thr <- max(weak_rel * sv_max, tol)
  sv_class <- ifelse(s1$d <= tol, "Unidentified",
              ifelse(s1$d < weak_thr, "Weak", "Identified"))
  list(ok = TRUE, status = sprintf("performed with Rmpfr at %d bits", prec_bits),
       rank = sum(s1$d > tol), tol = tol, fd_noise = noise,
       singular_values = s1$d, sv_class = sv_class,
       V = matrix(s1$v, n_par, n_par, dimnames = list(names(theta), paste0("sv_", seq_len(n_par)))))
}


#' One-sided (Hestenes) Jacobi SVD in mpfr
#'
#' @param cols list of n mpfr column vectors (each length m).
#' @return list(d = double singular values (length n, decreasing),
#'   v = double n x n right singular vectors)
#' @noRd
.d30_jacobi_svd <- function(cols, prec_bits, max_sweeps = 60L) {
  n <- length(cols)
  one <- Rmpfr::mpfr(1, precBits = prec_bits)
  zero <- Rmpfr::mpfr(0, precBits = prec_bits)
  V <- lapply(seq_len(n), function(j) { v <- rep(zero, n); v[j] <- one; v })
  tolj <- Rmpfr::mpfr(2, precBits = prec_bits)^(-(prec_bits - 8L))
  norm2 <- lapply(cols, function(cl) sum(cl^2))
  big <- norm2[[1]]
  for (j in seq_len(n)[-1]) if (as.logical(norm2[[j]] > big)) big <- norm2[[j]]
  tiny <- big * tolj^2
  sweep_i <- 0L
  rotated <- n > 1L
  while (rotated && sweep_i < max_sweeps) {
    rotated <- FALSE
    sweep_i <- sweep_i + 1L
    for (p in seq_len(n - 1L)) for (q in (p + 1L):n) {
      alpha <- sum(cols[[p]]^2); beta <- sum(cols[[q]]^2)
      gamma <- sum(cols[[p]] * cols[[q]])
      # Columns at or below working precision relative to the largest one
      # are numerically zero: rotating them only stirs rounding noise.
      if (as.logical(alpha <= tiny | beta <= tiny)) next
      if (!as.logical(abs(gamma) > tolj * sqrt(alpha * beta))) next
      zeta <- (beta - alpha) / (2 * gamma)
      sgn <- if (as.logical(zeta >= 0)) 1 else -1
      t <- sgn / (abs(zeta) + sqrt(one + zeta^2))
      cc <- one / sqrt(one + t^2); ss <- cc * t
      ap <- cols[[p]]; aq <- cols[[q]]
      cols[[p]] <- cc * ap - ss * aq; cols[[q]] <- ss * ap + cc * aq
      vp <- V[[p]]; vq <- V[[q]]
      V[[p]] <- cc * vp - ss * vq; V[[q]] <- ss * vp + cc * vq
      rotated <- TRUE
    }
  }
  if (rotated) .dynhr_warn("d30: Jacobi SVD did not converge; singular values may be inaccurate.")
  d <- vapply(cols, function(cl) Rmpfr::asNumeric(sqrt(sum(cl^2))), numeric(1))
  ord <- order(d, decreasing = TRUE)
  v <- matrix(vapply(V, Rmpfr::asNumeric, numeric(n)), n, n)
  list(d = d[ord], v = v[, ord, drop = FALSE])
}


#' D30 singular-value comparison plot
#' @noRd
.d30_plot <- function(svc, rk_d, hp, n_par, prec_bits, hp_status) {
  floor_v <- 1e-300
  hp_lab <- sprintf("%d-bit (Rmpfr)", prec_bits)
  df <- data.frame(index = svc$index, value = pmax(svc$double, floor_v),
                   precision = "double", stringsAsFactors = FALSE)
  thr <- data.frame(y = rk_d$tol, precision = "double", stringsAsFactors = FALSE)
  if (!is.null(hp)) {
    df <- rbind(df, data.frame(index = svc$index, value = pmax(svc$high_prec, floor_v),
                               precision = hp_lab, stringsAsFactors = FALSE))
    thr <- rbind(thr, data.frame(y = max(hp$tol, floor_v), precision = hp_lab,
                                 stringsAsFactors = FALSE))
  }
  lv <- c("double", hp_lab)
  df$precision <- factor(df$precision, levels = lv)
  thr$precision <- factor(thr$precision, levels = lv)
  cols <- stats::setNames(unname(tol_vibrant[c("blue", "orange")]), lv)
  subtitle <- if (is.null(hp)) {
    sprintf("double rank %d / %d. High-precision check not performed: %s",
            rk_d$rank, n_par, hp_status)
  } else {
    sprintf("rank: double %d / %d, %d-bit %d / %d. Dashed lines = rank tolerance of each precision",
            rk_d$rank, n_par, prec_bits, hp$rank, n_par)
  }
  ggplot2::ggplot(df, ggplot2::aes(x = index, y = value, colour = precision)) +
    ggplot2::geom_hline(data = thr, ggplot2::aes(yintercept = y, colour = precision),
                        linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(shape = precision), size = 2.5) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(breaks = svc$index) +
    ggplot2::scale_colour_manual(values = cols, name = "Precision", drop = TRUE) +
    ggplot2::scale_shape_manual(values = c(16, 17), name = "Precision", drop = TRUE) +
    theme_dynhr_diagnostic() +
    ggplot2::labs(
      title = "D30: identification singular values, double vs high precision",
      subtitle = paste(strwrap(subtitle, 95), collapse = "\n"),
      x = "Singular value index (equilibrated Jacobian)",
      y = "Singular value (log scale)"
    )
}
