## R/solve-qz.R
## --------------------------------------------------------------------------
## Ordered generalized Schur (QZ) backend. ONE implementation -- the compiled
## LAPACK dgges with a stable-first SORT callback (src/ordered_qz.cpp), via
## .ordered_qz() -- serves .solve_from_system() and the HANK Reiter state
## spaces. The R-level backends (QZ package, geigen, pure R) are reachable only
## with options(dynhr.use_rcpp = FALSE), as a debugging reference; they are
## never a fallback for a failed dgges.
##
## Phase-1 split from perturbation-monolith.R (no logic changes).
## --------------------------------------------------------------------------

# =====================================================================
# QZ decomposition dispatcher
# =====================================================================

#' Generalized Schur (QZ) decomposition with multiple backends
#'
#' Compiled ordered dgges (see \code{.ordered_qz()}); a LAPACK failure aborts
#' with class \code{dynhr_error_qz_failed}. With
#' \code{options(dynhr.use_rcpp = FALSE)} only: QZ package, then geigen, then a
#' pure-R eigendecomposition.
#' Returns S, T, Q, Z matrices with stable eigenvalues ordered first.
#'
#' @param D       Left matrix of pencil D * z_t = E * z_{t-1}
#' @param E       Right matrix of pencil
#' @param n_minus Number of stable eigenvalues expected (for info only)
#' @param verbose Print progress
#' @param abort   \code{FALSE}: return \code{NULL} on a LAPACK failure instead
#'   of aborting (the per-draw estimation solve).
#' @return List with T, S, Q, Z, eigenvalues and the eigenvalue counts
#' @noRd
## Is the compiled ordered-dgges backend available (and not disabled)?
.HAS_RCPP_QZ <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("ordered_qz_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

## Ordered generalized Schur of the pencil (A, B): A = Q S Z', B = Q T Z', with
## the eigenvalues |alpha/beta| < critmod ("stable"; beta = 0 is unstable)
## ordered first. Returns S, T, Q, Z, ALPHAR, ALPHAI, BETA and sdim (the number
## of stable eigenvalues). The single entry point for every ordered QZ in the
## package.
##
## Compiled path: ordered_qz_cpp() (dgges with SORT). It reports LAPACK's INFO
## rather than throwing; any non-zero INFO means the factors are not a valid
## stable-first ordering, and aborts here (classed, INFO in the message) -- callers
## that treat an unsolvable draw as infeasible catch the class.
## options(dynhr.use_rcpp = FALSE): the QZ package's dgges + dtgsen, kept as an
## independent reference for debugging.
## @noRd
.ordered_qz <- function(A, B, critmod = 1 + 1e-6, what = "ordered_qz",
                        abort = TRUE) {
  N <- nrow(A)
  qz_fail <- function(info, backend) {
    ## abort = FALSE: hand the failure back (ok = FALSE) for a caller that
    ## treats it as an infeasible point (the per-draw estimation solve).
    if (!isTRUE(abort)) return(list(ok = FALSE, info = info, backend = backend))
    .dynhr_abort(sprintf(paste0(
      "%s: the ordered QZ decomposition failed (%s, LAPACK INFO = %s) on a ",
      "%d x %d pencil%s."), what, backend, format(info), N, N,
      if (is.finite(info) && info >= 1 && info <= N) " -- the QZ iteration did not converge"
      else if (isTRUE(info == N + 2L)) " -- after reordering, roundoff moved an eigenvalue across the stability threshold"
      else if (isTRUE(info == N + 3L)) " -- the eigenvalue reordering failed (ill-conditioned swap)"
      else ""),
      class = "dynhr_error_qz_failed")
  }
  if (.HAS_RCPP_QZ()) {
    oq <- ordered_qz_cpp(A, B, critmod)
    if (!isTRUE(oq$ok)) return(qz_fail(oq$info %||% NA_integer_, "dgges"))
    return(oq)
  }
  if (!requireNamespace("QZ", quietly = TRUE))
    .dynhr_abort(what, ": options(dynhr.use_rcpp = FALSE) needs the 'QZ' ",
                 "package for the ordered QZ decomposition.",
                 class = "dynhr_error_qz_failed")
  qzd <- QZ::qz.dgges(A, B)
  if (!isTRUE(qzd$INFO == 0L)) return(qz_fail(qzd$INFO, "QZ::qz.dgges"))
  al2 <- qzd$ALPHAR^2 + qzd$ALPHAI^2
  sel <- al2 < critmod^2 * qzd$BETA^2
  o <- QZ::qz.dtgsen(qzd$S, qzd$T, qzd$Q, qzd$Z, select = sel)
  if (!isTRUE(o$INFO == 0L)) return(qz_fail(o$INFO, "QZ::qz.dtgsen"))
  list(ok = TRUE, info = 0L, sdim = sum(sel), S = o$S, T = o$T, Q = o$Q,
       Z = o$Z, ALPHAR = o$ALPHAR, ALPHAI = o$ALPHAI, BETA = o$BETA)
}

## Classify the (already stable-ordered) generalized eigenvalues from a Schur
## decomposition's ALPHAR/ALPHAI/BETA and assemble the .solve_qz return list.
## Shared by the C++ (ordered dgges) and R (QZ package) paths so both classify
## identically.
##
## Unit-circle convention: |λ| < 1-tol stable, > 1+tol unstable, |λ-1| <= tol
## borderline (counted separately so the caller can reclassify for Ramsey/
## unit-root models). tol = 1e-6 matches the dgges SORT threshold (1+1e-6) and
## the convention used elsewhere in dynhr (backend-monolith.R: 1 - 1e-6).
.qz_classify <- function(S, T, Q, Z, ALPHAR, ALPHAI, BETA, tol = 1e-6) {
  eig_vals_ord <- ifelse(abs(BETA) > 1e-14,
    complex(real = ALPHAR, imaginary = ALPHAI) / BETA, Inf + 0i)
  eig_abs_ord <- ifelse(is.finite(Mod(eig_vals_ord)), Mod(eig_vals_ord), Inf)
  list(T = T, S = S, Q = Q, Z = Z,
       eigenvalues = eig_vals_ord,
       n_unstable_finite   = sum(is.finite(eig_abs_ord) & eig_abs_ord > 1 + tol),
       n_unstable_infinite = sum(!is.finite(eig_abs_ord)),
       n_on_unit_circle    = sum(is.finite(eig_abs_ord) &
                                 abs(eig_abs_ord - 1) <= tol))
}

.solve_qz <- function(D, E, n_minus, verbose = TRUE, abort = TRUE) {
  N <- nrow(D)

  if (N == 0) {
    return(list(T = matrix(0, 0, 0), S = matrix(0, 0, 0),
                Q = matrix(0, 0, 0), Z = matrix(0, 0, 0),
                eigenvalues = complex(0), n_stable = 0L))
  }

  # Ordered generalized Schur via LAPACK dgges + SORT callback (C++): stable
  # eigenvalues (|lambda| < 1+1e-6) are ordered first DURING the decomposition
  # (Dynare's mjdgges). A LAPACK failure aborts inside .ordered_qz(): there is
  # no second-library fallback, because QZ::qz() + qz.dtgsen() run the same
  # LAPACK routines and would meet the same condition.
  if (.HAS_RCPP_QZ()) {
    oq <- .ordered_qz(E, D, what = "solve_perturbation", abort = abort)
    if (!isTRUE(oq$ok)) return(NULL)        # abort = FALSE only
    if (verbose) .dynhr_cat("  Using ordered dgges (C++)...\n")
    return(.qz_classify(oq$S, oq$T, oq$Q, oq$Z, oq$ALPHAR, oq$ALPHAI, oq$BETA))
  }

  # options(dynhr.use_rcpp = FALSE) only -- R-level reference backends.
  if (requireNamespace("QZ", quietly = TRUE)) {
    if (verbose) .dynhr_cat("  Using QZ package with Eigenvalue update...\n")
    qz <- QZ::qz(E, D)
    eig_vals <- ifelse(abs(qz$BETA) > 1e-14,
      complex(real = qz$ALPHAR, imaginary = qz$ALPHAI) / qz$BETA, Inf + 0i)
    eig_abs <- ifelse(is.finite(Mod(eig_vals)), Mod(eig_vals), Inf)
    tol <- 1e-6
    qz_stable <- eig_abs <= 1 + tol
    qz_ord <- QZ::qz.dtgsen(qz$S, qz$T, qz$Q, qz$Z, select = qz_stable)
    return(.qz_classify(qz_ord$S, qz_ord$T, qz_ord$Q, qz_ord$Z,
                        qz_ord$ALPHAR, qz_ord$ALPHAI, qz_ord$BETA))
  }

  # Use geigen package if available
  if (requireNamespace("geigen", quietly = TRUE)) {
    if (verbose) .dynhr_cat("  Using geigen package...\n")
    ge <- geigen::geigen(E, D, symmetric = FALSE)
    eig_vals <- ge$values
    ord <- order(Mod(eig_vals))
    Z <- ge$vectors[, ord, drop = FALSE]
    eig_vals_ord <- eig_vals[ord]
    # geigen can return Inf/NaN eigenvalues when beta ~ 0 (infinite generalized
    # eigenvalues); treat any non-finite |λ| as infinite, matching .qz_classify.
    eig_abs_ord <- ifelse(is.finite(Mod(eig_vals_ord)), Mod(eig_vals_ord), Inf)
    tol <- 1e-6
    return(list(T = diag(eig_vals_ord, nrow = N), S = diag(N),
                Q = diag(N), Z = Z, eigenvalues = eig_vals_ord,
                n_stable = sum(Mod(eig_vals) <= 1 + 1e-6),
                n_unstable_finite   = sum(is.finite(eig_abs_ord) & eig_abs_ord > 1 + tol),
                n_unstable_infinite = sum(!is.finite(eig_abs_ord)),
                n_on_unit_circle    = sum(is.finite(eig_abs_ord) &
                                          abs(eig_abs_ord - 1) <= tol)))
  }

  # Pure R fallback. eigen() on D^{-1} %*% E is only valid when D is
  # nonsingular (i.e. the pencil has no infinite generalized eigenvalues);
  # regularising a near-singular D would silently produce a wrong decision
  # rule, so refuse instead.
  if (verbose) .dynhr_cat("  Using pure-R eigendecomposition fallback...\n")
  if (1 / kappa(D, exact = FALSE) < 1e-12) {
    stop("Pure-R QZ fallback cannot handle a numerically singular D matrix ",
         "(pencil has infinite generalized eigenvalues). Install the 'QZ' ",
         "or 'geigen' package for proper handling of this model.")
  }
  M <- solve(D) %*% E
  eig <- eigen(M)
  eig_vals <- eig$values
  ord <- order(Mod(eig_vals))
  eig_vals_ord <- eig_vals[ord]
  tol <- 1e-6
  list(T = diag(eig_vals_ord, nrow = N), S = diag(N),
       Q = diag(N), Z = eig$vectors[, ord, drop = FALSE],
       eigenvalues = eig_vals_ord,
       n_stable = sum(Mod(eig_vals) <= 1 + 1e-6),
       n_unstable_finite   = sum(Mod(eig_vals_ord) > 1 + tol),
       n_unstable_infinite = 0L,
       n_on_unit_circle    = sum(abs(Mod(eig_vals_ord) - 1) <= tol))
}
