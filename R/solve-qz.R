## R/solve-qz.R
## --------------------------------------------------------------------------
## QZ-decomposition backend with graceful fallback: tries QZ package first,
## then geigen, then a pure-R eigendecomposition + regularisation. All three
## return the same shape so .solve_from_system() does not branch.
##
## Phase-1 split from perturbation-monolith.R (no logic changes).
## --------------------------------------------------------------------------

# =====================================================================
# QZ decomposition dispatcher
# =====================================================================

#' Generalized Schur (QZ) decomposition with multiple backends
#'
#' Tries QZ package, then geigen, then pure-R eigendecomposition fallback.
#' Returns S, T, Q, Z matrices with stable eigenvalues ordered first.
#'
#' @param D       Left matrix of pencil D * z_t = E * z_{t-1}
#' @param E       Right matrix of pencil
#' @param n_minus Number of stable eigenvalues expected (for info only)
#' @param verbose Print progress
#' @return List with T, S, Q, Z, eigenvalues, n_stable, or NULL on failure
#' @noRd
## Is the compiled ordered-dgges backend available (and not disabled)?
.HAS_RCPP_QZ <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("ordered_qz_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
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

.solve_qz <- function(D, E, n_minus, verbose = TRUE) {
  N <- nrow(D)

  if (N == 0) {
    return(list(T = matrix(0, 0, 0), S = matrix(0, 0, 0),
                Q = matrix(0, 0, 0), Z = matrix(0, 0, 0),
                eigenvalues = complex(0), n_stable = 0L))
  }

  # Primary: ordered generalized Schur via LAPACK dgges + SORT callback (C++).
  # Orders stable eigenvalues (|λ| < 1+1e-6) first DURING the decomposition,
  # avoiding the QZ::qz.dtgsen post-hoc reorder, which has an input-independent
  # memory-corruption bug (see src/ordered_qz.cpp). Returns the same Schur form
  # and ALPHAR/ALPHAI/BETA as QZ::qz, so classification is unchanged.
  if (.HAS_RCPP_QZ()) {
    oq <- tryCatch(ordered_qz_cpp(E, D), error = function(e) NULL)
    if (!is.null(oq) && isTRUE(oq$ok)) {
      if (verbose) cat("  Using ordered dgges (C++)...\n")
      return(.qz_classify(oq$S, oq$T, oq$Q, oq$Z,
                          oq$ALPHAR, oq$ALPHAI, oq$BETA))
    }
    if (verbose)
      cat("  C++ ordered QZ unavailable/failed; trying QZ package...\n")
  }

  # Fallback: QZ package (QZ::qz + qz.dtgsen). Retained for builds without the
  # compiled DLL and for options(dynhr.use_rcpp = FALSE). NOTE: qz.dtgsen here
  # is the routine with the memory bug; this path is intentionally secondary.
  if (requireNamespace("QZ", quietly = TRUE)) {
    if (verbose) cat("  Using QZ package with Eigenvalue update...\n")
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
    if (verbose) cat("  Using geigen package...\n")
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
  if (verbose) cat("  Using pure-R eigendecomposition fallback...\n")
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
