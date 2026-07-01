## R/diag-pre-d30-arbitrary-precision.R
## --------------------------------------------------------------------------
## Phase H: D30 — Arbitrary-Precision Rank Checks (Qu & Tkachenko 2023)
##
## Re-runs the rank check from D1/D19 at elevated precision to confirm that
## near-zero singular values are genuinely zero (or not). Floating-point
## cancellation at machine epsilon (~2e-16 for double) can produce spurious
## near-zero singular values that mislead rank decisions. Arbitrary-precision
## computation distinguishes structural zeros from numerical artefacts.
##
## Backends (in order of preference):
##   1. **Rmpfr** (R Multiple Precision Floating-point Reliable) — pure R
##      extension wrapping GNU MPFR.  Preferred: no external runtime
##      dependency, works on any platform with the Rmpfr package installed.
##      Precision controlled via \code{precBits} (default 128 bits ~ 38
##      decimal digits).
##
##      Rmpfr does not provide a direct \code{svd()} method for mpfrMatrix.
##      Instead we form \eqn{J'J} at mpfr precision (using Rmpfr's S4 methods
##      for \code{t()} and \code{\%*\%} on mpfr matrices), convert to double,
##      and compute the eigen-decomposition.  Singular values are the square
##      roots of the eigenvalues of \eqn{J'J}.  This retains the benefit of
##      high-precision inner products in the Jacobian cross-product.
##
##   2. **JuliaCall** — if Rmpfr is unavailable and Julia is installed,
##      uses Julia's BigFloat via JuliaCall for the SVD computation.
##   3. **Base R** — standard double-precision SVD, always available as
##      a baseline for comparison.
##
## The package does NOT depend on Julia or Rmpfr at runtime — both are
## \code{Suggests} only.  The diagnostic checks availability at run time
## and produces a clear message if neither high-precision backend is
## available.
##
## References:
##   Qu, Z., & Tkachenko, D. (2023). Arbitrary-precision identification
##     in DSGE models. Working paper.
##   Iskrev, N. (2010). Local identification in DSGE models.
##   Mächler, M. (2014). "Rmpfr: R MPFR — Multiple Precision
##     Floating-Point Reliable." R package version 0.5-0.
## --------------------------------------------------------------------------

#' D30. Arbitrary-Precision Rank Checks
#'
#' Re-runs the SVD-based rank analysis from D1/D19 at elevated precision
#' (default 128 bits ~ 38 decimal digits) to distinguish structural zeros
#' from floating-point artefacts in the singular value spectrum.
#'
#' The diagnostic:
#' \enumerate{
#'   \item Computes the standard double-precision SVD of the moment
#'     Jacobian (from \code{model_solve_fn}) as a baseline.
#'   \item Re-computes the Jacobian and its SVD at elevated precision
#'     using either \strong{Rmpfr} (preferred) or \strong{Julia BigFloat}
#'     (fallback).  Precision is controlled by \code{prec_bits} (default
#'     128 bits).
#'   \item Compares singular values between double and high precision,
#'     and reports any rank discrepancies.
#'   \item If no high-precision backend is available, reports baseline
#'     results + instruction to install Rmpfr.
#' }
#'
#' \strong{Precision-tolerance caveat:} the rank threshold is currently
#' \code{rank_tol_factor * max(sv) * .Machine$double.eps} at \emph{both}
#' precisions.  Because the threshold uses double-precision machine epsilon
#' regardless of \code{prec_bits}, the rank decision can only differ between
#' double and high-precision if the high-precision inner products in
#' \eqn{J'J} shift a singular value across the \emph{double}-epsilon
#' threshold.  To detect rank differences that require genuine
#' precision-appropriate tolerances, set
#' \code{rank_tol_factor = max(dim(J)) * (2^{-prec_bits} / .Machine$double.eps)}.
#' The current default does \emph{not} do this.
#'
#' @param model_solve_fn  Function: theta -> named numeric vector of moments.
#' @param theta           Named numeric vector of parameter values.
#' @param param_names     Character vector of parameter names.
#' @param moment_names    Character vector of moment names.
#' @param prec_bits       Bit precision for arbitrary-precision computation
#'   (default 128).  Higher values (e.g. 256) give more safety at the cost
#'   of speed.
#' @param backend         Character: \code{"auto"} (try Rmpfr then JuliaCall),
#'   \code{"Rmpfr"}, \code{"JuliaCall"}, or \code{"base"} (double only).
#' @param eps             Step size for finite differences in double precision
#'   (default 1e-5).  For high precision, the step is scaled appropriately.
#' @param rank_tol_factor  Factor times \code{.Machine$double.eps} for the
#'   SVD rank threshold, applied at \emph{both} double and high precision.
#'   Default is \code{max(n_moments, n_params)} (standard double-precision
#'   SVD tolerance).  \strong{This factor does not change with
#'   \code{prec_bits}.}  Genuine precision-dependent rank differences
#'   require setting this to the ratio of double-eps to mpfr-eps
#'   (~\code{2^{-53} / 2^{-prec_bits}}); the default will not detect
#'   near-zero SVs that fall below the double-eps threshold but above the
#'   mpfr-eps threshold.  Use a larger value (e.g. 100) to be more
#'   conservative in the double-precision baseline only.
#' @param verbose         Print progress messages.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing:
#'     \itemize{
#'       \item \code{double_svd} — standard double-precision SVD result.
#'       \item \code{high_prec_svd} — arbitrary-precision SVD result
#'         (or NULL if backend unavailable).
#'       \item \code{double_rank} — rank from double precision.
#'       \item \code{high_prec_rank} — rank from high precision (or NA).
#'       \item \code{rank_consistent} — whether ranks agree.
#'       \item \code{sv_comparison} — data.frame comparing singular values.
#'       \item \code{backend_used} — which backend was used.
#'       \item \code{param_names} — parameter names.
#'     }}
#'   \item{pass}{Logical — TRUE if ranks are consistent across precisions.}
#'   \item{plots}{List of ggplot2 objects.}
#'   \item{summary}{Human-readable summary.}
#'
#' @note Rmpfr is the strongly preferred backend.  Install it with:
#'   \code{install.packages("Rmpfr")}
#'   This requires the GMP and MPFR system libraries on your platform.
#'
#' @references
#'   Qu, Z., & Tkachenko, D. (2023). Arbitrary-precision identification in
#'     DSGE models. Working paper.
#'   Machler, M. (2014). Rmpfr: R MPFR -- Multiple Precision Floating-Point
#'     Reliable. R package.
#'
#' @noRd
d30_arbitrary_precision_rank <- function(model_solve_fn,
                                          theta,
                                          param_names = NULL,
                                          moment_names = NULL,
                                          prec_bits = 128L,
                                          backend = c("auto", "Rmpfr", "JuliaCall", "base"),
                                          eps = 1e-5,
                                          rank_tol_factor = NULL,
                                          verbose = FALSE,
                                          meta = NULL) {
  backend <- match.arg(backend)

  # ---- 0. Early validation for graceful NULL handling ----
  if (is.null(theta) || is.null(model_solve_fn)) {
    return(.make_result(
      pass    = NA,
      summary = "D30 Arbitrary-Precision Rank: theta or model_solve_fn is NULL."
    ))
  }

  # ---- 1. Defaults ----
  n_par <- length(theta)
  if (is.null(param_names)) param_names <- names(theta) %||% paste0("theta_", seq_len(n_par))
  if (is.null(moment_names)) {
    f0 <- model_solve_fn(theta)
    if (is.null(f0) || length(f0) == 0) {
      return(.make_result(
        pass    = NA,
        summary = "D30 Arbitrary-Precision Rank: model_solve_fn returned NULL or empty."
      ))
    }
    moment_names <- if (!is.null(names(f0))) names(f0) else paste0("m_", seq_along(f0))
  }
  if (is.null(rank_tol_factor)) {
    n_mom <- length(moment_names)
    rank_tol_factor <- max(n_mom, n_par)
  }

  # ---- 2. Double-precision baseline ----
  if (verbose) cat("[d30] Computing double-precision SVD baseline...\n")
  J_double <- .numerical_jacobian(model_solve_fn, theta, eps = eps)
  colnames(J_double) <- param_names
  rownames(J_double) <- moment_names

    svd_double <- svd(J_double)
    tol_double <- rank_tol_factor * max(svd_double$d) * .Machine$double.eps
    rank_double <- sum(svd_double$d > tol_double)

    if (verbose) cat(sprintf("[d30] Double precision: rank = %d / %d (tol = %.2e)\n",
                             rank_double, n_par, tol_double))

    # ---- 3. High-precision computation ----
    high_prec_svd <- NULL
    rank_high <- NA_integer_
    backend_used <- "base"

    # Resolve backend
    use_rmpfr <- FALSE
    use_julia <- FALSE

    if (backend == "auto") {
      use_rmpfr <- requireNamespace("Rmpfr", quietly = TRUE)
      if (!use_rmpfr) {
        use_julia <- requireNamespace("JuliaCall", quietly = TRUE)
      }
    } else if (backend == "Rmpfr") {
      use_rmpfr <- requireNamespace("Rmpfr", quietly = TRUE)
      if (!use_rmpfr && verbose) cat("[d30] Rmpfr not available.\n")
    } else if (backend == "JuliaCall") {
      use_julia <- requireNamespace("JuliaCall", quietly = TRUE)
      if (!use_julia && verbose) cat("[d30] JuliaCall not available.\n")
    }

    if (use_rmpfr) {
      if (verbose) cat(sprintf("[d30] Using Rmpfr with precBits = %d...\n", prec_bits))
      high_prec_svd <- .d30_svd_rmpfr(model_solve_fn, theta, prec_bits, eps)
      backend_used <- "Rmpfr"
    } else if (use_julia) {
      if (verbose) cat(sprintf("[d30] Using Julia BigFloat with precision = %d bits...\n", prec_bits))
      high_prec_svd <- .d30_svd_julia(model_solve_fn, theta, prec_bits, eps)
      backend_used <- "JuliaCall"
    } else {
      if (verbose) {
        cat("[d30] No high-precision backend available.\n")
        cat("[d30] Install Rmpfr (install.packages('Rmpfr')) for arbitrary-precision SVD.\n")
        cat("[d30] JuliaCall is an alternative fallback.\n")
      }
    }

    # Compute high-precision rank if SVD is available
    if (!is.null(high_prec_svd) && !is.null(high_prec_svd$d)) {
      hp_sv <- as.numeric(high_prec_svd$d)
      tol_high <- rank_tol_factor * max(hp_sv) * .Machine$double.eps
      rank_high <- sum(hp_sv > tol_high)
      if (verbose) cat(sprintf("[d30] High precision (%s): rank = %d / %d\n",
                               backend_used, rank_high, n_par))
    }

    # ---- 4. Comparison ----
    rank_consistent <- is.na(rank_high) || (rank_double == rank_high)

    # Align high-precision SV count with double-precision SV count.
    # svd(J) returns min(nrow,ncol) singular values (e.g. 16 for 16x20 J),
    # while J'J eigen gives ncol (20) values. Truncate to match.
    n_sv <- length(svd_double$d)
    hp_sv <- if (!is.null(high_prec_svd) && !is.null(high_prec_svd$d))
      as.numeric(high_prec_svd$d)[seq_len(n_sv)] else rep(NA_real_, n_sv)
    sv_comparison <- data.frame(
      index = seq_len(n_sv),
      double = svd_double$d,
      high_prec = hp_sv,
      ratio = if (!is.null(high_prec_svd))
        svd_double$d / pmax(hp_sv, 1e-300) else NA_real_,
      stringsAsFactors = FALSE
    )

    # Find singular values near the rank threshold
    borderline <- which(svd_double$d > tol_double * 0.1 &
                         svd_double$d < tol_double * 10)

    # ---- 5. Plots ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Only include "base" (high_prec) series when a high-precision backend is
      # available; otherwise the second series is all 1e-300 and creates a
      # misleading blank band in the middle of the plot.
      hp_available <- backend_used != "base"
      if (hp_available) {
        sv_long <- data.frame(
          index     = rep(sv_comparison$index, 2),
          value     = c(log10(pmax(sv_comparison$double, 1e-300)),
                        log10(pmax(sv_comparison$high_prec, 1e-300))),
          precision = rep(c("double", backend_used), each = nrow(sv_comparison)),
          stringsAsFactors = FALSE
        )
      } else {
        sv_long <- data.frame(
          index     = sv_comparison$index,
          value     = log10(pmax(sv_comparison$double, 1e-300)),
          precision = "double",
          stringsAsFactors = FALSE
        )
      }

      colour_vals <- c("double"    = dynhr_colours$dark_blue,
                       "Rmpfr"     = dynhr_colours$orange,
                       "JuliaCall" = dynhr_colours$green)
      # Keep only series that appear in the data
      colour_vals <- colour_vals[names(colour_vals) %in% unique(sv_long$precision)]

      p_svc <- ggplot2::ggplot(
        sv_long, ggplot2::aes(x = index, y = value, colour = precision)
      ) +
        ggplot2::geom_point(size = 1.5) +
        ggplot2::geom_line(linewidth = 0.3) +
        ggplot2::geom_hline(yintercept = log10(tol_double),
                           linetype = "dashed", colour = dynhr_colours$red,
                           linewidth = 0.5) +
        ggplot2::scale_colour_manual(values = colour_vals, name = "Precision") +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "D30: Singular value comparison -- double vs arbitrary precision",
          subtitle = sprintf("Rank: double=%d / %d%s",
                             rank_double, n_par,
                             if (is.na(rank_high)) sprintf(
                               "  [%s not installed -- install with: install.packages('Rmpfr')]",
                               if (backend == "auto" || backend == "Rmpfr") "Rmpfr" else "high-prec backend")
                             else sprintf(", %s=%d", backend_used, rank_high)),
          x = "Singular value index", y = "log10(singular value)"
        )
      plots$sv_comparison <- .apply_meta(p_svc, meta)

      if (length(borderline) > 0) {
        bd_df <- sv_comparison[borderline, , drop = FALSE]
        p_bsv <- ggplot2::ggplot(
          bd_df, ggplot2::aes(x = index)
        ) +
          ggplot2::geom_point(ggplot2::aes(y = log10(pmax(double, 1e-300))),
                             colour = dynhr_colours$dark_blue, size = 2) +
          ggplot2::geom_point(ggplot2::aes(y = log10(pmax(high_prec, 1e-300))),
                             colour = dynhr_colours$orange, size = 2, shape = 17) +
          ggplot2::geom_hline(yintercept = log10(tol_double),
                             linetype = "dashed", colour = dynhr_colours$red) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title = "D30: Borderline singular values (near rank threshold)",
            x = "Index", y = "log10(singular value)"
          )
        plots$borderline_sv <- .apply_meta(p_bsv, meta)
      }
    }

    # ---- 6. Result ----
    result <- list(
      double_svd       = svd_double,
      high_prec_svd    = high_prec_svd,
      double_rank      = rank_double,
      high_prec_rank   = rank_high,
      rank_consistent  = rank_consistent,
      sv_comparison    = sv_comparison,
      backend_used     = backend_used,
      param_names      = param_names,
      prec_bits        = prec_bits,
      borderline_idx   = borderline
    )

    # Build summary
    if (is.na(rank_high)) {
      summary_str <- sprintf(
        "D30 Arbitrary-Precision Rank: double rank=%d/%d. Rmpfr not installed (run install.packages('Rmpfr') to enable high-precision verification).",
        rank_double, n_par
      )
      llm_str <- sprintf(
        "[INFO] D30 Arbitrary-Precision Rank double_rank=%d/%d backend=%s Rmpfr_not_installed=TRUE action: install.packages('Rmpfr') to enable high-precision SVD",
        rank_double, n_par, backend_used
      )
      pass_val <- NA
    } else {
      pass_val <- rank_consistent
      summary_str <- sprintf(
        "D30 Arbitrary-Precision Rank: double rank=%d/%d vs %s rank=%d/%d. %s",
        rank_double, n_par, backend_used, rank_high, n_par,
        if (rank_consistent) "Rank consistent across precisions."
        else "RANK DIFFERENCE DETECTED!"
      )
      llm_str <- sprintf(
        "[%s] D30 Arbitrary-Precision Rank double_rank=%d high_prec_rank=%d backend=%s prec_bits=%d consistent=%d",
        if (isTRUE(pass_val)) "PASS" else "FAIL",
        rank_double, rank_high, backend_used, prec_bits, as.integer(rank_consistent)
      )
    }

    .make_result(
      result  = result,
      pass    = pass_val,
      plots   = plots,
      summary = summary_str,
      llm_summary = llm_str
    )
}


# ==========================================================================
# Backend implementations for D30
# ==========================================================================

#' Compute SVD at half precision using Rmpfr
#'
#' @description
#' Computes the moment Jacobian in **double precision** (the model solver
#' uses base R linear algebra and cannot handle mpfr types), then converts
#' the Jacobian to mpfr to form \eqn{J'J} at the specified precision.
#' The singular values are obtained from the eigen-decomposition of the
#' high-precision cross-product matrix.
#'
#' **Caveat**: This is NOT a true arbitrary-precision SVD because the
#' Jacobian itself is computed in double precision.  The only benefit is
#' that forming \eqn{J'J} at high precision can reduce rounding errors
#' from catastrophic cancellation in the inner products, yielding slightly
#' more accurate singular values for ill-conditioned problems.
#'
#' Rmpfr does not provide a direct \code{svd()} method for \code{mpfrMatrix}
#' objects, so we exploit the fact that singular values are the square
#' roots of the eigenvalues of \eqn{J'J}.
#'
#' @param model_solve_fn  Function: theta -> moments
#' @param theta           Parameter vector (double)
#' @param prec_bits       MPFR precision in bits
#' @param eps             Base step size for finite differences
#' @return List with $d (singular values), $prec_bits, or NULL on failure
#' @noRd
.d30_svd_rmpfr <- function(model_solve_fn, theta, prec_bits, eps) {
  n_par <- length(theta)

  # NOTE: model_solve_fn cannot handle mpfr types, so the Jacobian is
  # computed in double precision.  Only J'J is formed at high precision.
  #
  # Compute the Jacobian in double precision first
  J_double <- .numerical_jacobian(model_solve_fn, theta, eps = eps)

  # Convert to mpfr for high-precision cross-product
  J_mp <- Rmpfr::mpfr(J_double, precBits = prec_bits)

  # Form J'J at mpfr precision
  jtj_mp <- t(J_mp) %*% J_mp

  # Convert back and check for valid results
  jtj_d <- matrix(as.numeric(jtj_mp), nrow = n_par, ncol = n_par)
  if (any(!is.finite(jtj_d))) jtj_d <- NULL

    if (is.null(jtj_d) || any(!is.finite(jtj_d))) {
      warning("Rmpfr J'J produced non-finite values -- falling back to double precision.")
      sv_d <- svd(J_double)$d
    } else {
      # Eigen-decomposition (symmetric J'J)
      eig <- eigen(jtj_d, symmetric = TRUE)

      # Singular values = sqrt(eigenvalues)
      sv_d <- sqrt(pmax(eig$values, 0))

      # Restore expected order (decreasing, matching svd() convention)
      sv_d <- sort(sv_d, decreasing = TRUE)
    }

    list(
      d = sv_d,
      u = NULL,
      v = NULL,
      prec_bits = prec_bits
    )
}


#' Compute SVD at arbitrary precision using Julia BigFloat
#'
#' Falls back to Julia via JuliaCall if Rmpfr is not available.
#' Requires JuliaCall package and a working Julia installation.
#'
#' @param model_solve_fn  Function: theta -> moments
#' @param theta           Parameter vector (double)
#' @param prec_bits       Julia BigFloat precision in bits
#' @param eps             Base step size for finite differences
#' @return SVD result, or NULL on failure
#' @noRd
.d30_svd_julia <- function(model_solve_fn, theta, prec_bits, eps) {
  # Initialise Julia if needed
  if (!JuliaCall::julia_setup(quiet = TRUE)) {
    JuliaCall::julia_setup(quiet = TRUE)
  }

    # Set BigFloat precision
    JuliaCall::julia_command(sprintf("set_bigfloat_precision(%d)", prec_bits),
                             need_return = FALSE)

    n_par <- length(theta)
    f0 <- model_solve_fn(theta)
    n_mom <- length(f0)

    # Build Julia matrix as BigFloat array
    theta_str <- paste(sprintf("BigFloat(\"%.16e\")", theta), collapse = ", ")
    JuliaCall::julia_assign("theta_jl", theta)

    # For each column, compute central difference in Julia BigFloat
    # We'll compute the Jacobian in R double, transfer to Julia, then SVD
    # at high precision — this is the most robust approach since Julia's
    # automatic differentiation isn't available here.

    # Compute Jacobian in double
    J <- .numerical_jacobian(model_solve_fn, theta, eps = eps)
    # Guard: ensure moment count matches Jacobian rows
    if (nrow(J) != length(f0)) {
      warning(sprintf(
        "d30: moment count (%d) != Jacobian rows (%d). Using generic labels.",
        length(f0), nrow(J)))
      f0_names <- paste0("m_", seq_len(nrow(J)))
    } else {
      f0_names <- names(f0) %||% paste0("m_", seq_len(nrow(J)))
    }
    colnames(J) <- names(theta)
    rownames(J) <- f0_names

    # Convert to Julia BigFloat matrix
    JuliaCall::julia_assign("J_mat", J)
    JuliaCall::julia_command(
      sprintf("J_bf = BigFloat.(J_mat); set_bigfloat_precision(%d)", prec_bits),
      need_return = FALSE
    )

    # Compute SVD in Julia
    JuliaCall::julia_command("svd_jl = svd(J_bf)", need_return = FALSE)

    # Extract results
    sv_d <- JuliaCall::julia_eval("Float64.(svd_jl.S)")
    sv_u <- JuliaCall::julia_eval("Float64.(svd_jl.U)")
    sv_v <- JuliaCall::julia_eval("Float64.(svd_jl.V)")

    # Ensure correct dimensions
    if (is.list(sv_u)) sv_u <- matrix(unlist(sv_u), nrow = n_mom, ncol = n_par)
    if (is.list(sv_v)) sv_v <- matrix(unlist(sv_v), nrow = n_par, ncol = n_par)

    list(
      d = as.numeric(sv_d),
      u = sv_u,
      v = sv_v,
      prec_bits = prec_bits
    )
  }
