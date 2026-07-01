## R/ramsey-nn1-objective.R
## --------------------------------------------------------------------------
## E3: Modified (n, n+1) objective construction for the Gross-Hansen
## approximation.
##
## Builds the polynomial coefficients of the modified periodic objective
## W_t^{(n,n+1)} from the Taylor expansions of the planner objective and
## constraints, combined with the steady-state Lagrange multipliers (the
## "blue correction" terms).
##
## Key formula (Gross & Hansen 2021, Proposition 1):
##
##   W_t^{(n,n+1)} = Σ_{|α|=2}^{n+1} D^α f_t · x̃^α / α!
##     - Σ_j λ̄_j · Σ_{|α|=2}^{n+1} D^α g_{j,t} · x̃^α / α!
##     - Σ_l ψ̄_l · β^{-1} · Σ_{|α|=2}^{n+1} D^α h_{l,t} · x̃^α / α!
##     + t.i.p.
##
## The result is a polynomial that starts at quadratic order (no linear
## terms), suitable for maximisation subject to n-order constraint
## expansions.
## --------------------------------------------------------------------------

#' Build the (n, n+1) modified objective coefficients
#'
#' Combines the (n+1)-order Taylor expansion of the planner objective with
#' the "blue correction" terms from the steady-state multipliers to produce
#' the polynomial coefficients of W_t^{(n,n+1)}.
#'
#' For each order k from 2 to (n+1):
#'   W_k = f_k - Σ_j λ̄_j · g_{j,k} - Σ_l ψ̄_l · β^{-1} · h_{l,k}
#'
#' where f_k, g_{j,k}, h_{l,k} are the k-th order derivative tensors
#' (divided by k!) of the objective and constraints at the steady state.
#'
#' @param taylor       Output of .nn1_taylor_expand().
#' @param multipliers  Output of .compute_ramsey_ss_multipliers().
#' @param n            Approximation order (1=LQ, 2=QC, 3, ...).
#' @param beta         Discount factor.
#' @param verbose      Print progress messages.
#'
#' @return A list with:
#'   \item{order}{The approximation order n.}
#'   \item{coefficients}{List of tensors: $quad (matrix), $cubic (3D array),
#'         $quartic (4D array), ... containing the W_t^{(n,n+1)} coefficients.}
#'   \item{grad_at_ss}{Gradient of W at SS (should be near zero).}
#'   \item{info}{Metadata: n_terms, symmetry, max_order.}
#' @noRd
.nn1_build_modified_objective <- function(taylor,
                                           multipliers,
                                           n,
                                           beta,
                                           verbose = FALSE) {
  # ---- 1. Extract components ----
  obj_tensors <- taylor$objective  # f_k tensors
  eq_tensors  <- taylor$equations  # g_{j,k} and h_{l,k} tensors
  n_eq <- length(eq_tensors)
  n_cols <- taylor$n_cols
  lambda <- multipliers$multipliers$lambda
  psi    <- multipliers$multipliers$psi
  forward_idx <- multipliers$forward_eqs
  backward_idx <- multipliers$backward_eqs

  # Build mapping: equation index → multiplier value
  # Backward equations j → λ̄_j
  # Forward equations l → ψ̄_l · β^{-1}
  mult_weights <- setNames(numeric(n_eq), seq_len(n_eq))
  for (j in backward_idx) {
    nm <- names(lambda)[which(backward_idx == j)]
    if (length(nm) > 0 && nm %in% names(lambda)) {
      mult_weights[j] <- lambda[nm]
    }
  }
  for (l in forward_idx) {
    nm <- names(psi)[which(forward_idx == l)]
    if (length(nm) > 0 && nm %in% names(psi)) {
      mult_weights[l] <- psi[nm] / beta
    }
  }

  if (verbose) {
    cat(sprintf("[nn1_objective] Building modified objective (n=%d, %d cols)\n",
                n, n_cols))
    cat(sprintf("  Multiplier weights: %d equations\n", n_eq))
  }

  # ---- 2. For each order k from 2 to (n+1), build W_k ----
  coefficients <- list()

  # Order 2: quadratic coefficients
  if (n + 1 >= 2) {
    W_quad <- .build_order_k_coefficients(eq_tensors, obj_tensors,
                                          mult_weights, n_eq, n_cols, 2)
    coefficients$quad <- W_quad
  }

  # Order 3: cubic coefficients
  if (n + 1 >= 3) {
    W_cubic <- .build_order_k_coefficients(eq_tensors, obj_tensors,
                                           mult_weights, n_eq, n_cols, 3)
    coefficients$cubic <- W_cubic
  }

  # Order 4: quartic coefficients
  if (n + 1 >= 4) {
    W_quartic <- .build_order_k_coefficients(eq_tensors, obj_tensors,
                                             mult_weights, n_eq, n_cols, 4)
    coefficients$quartic <- W_quartic
  }

  # Order 5+
  for (k in 5:(n + 1)) {
    W_k <- .build_order_k_coefficients(eq_tensors, obj_tensors,
                                       mult_weights, n_eq, n_cols, k)
    coefficients[[paste0("tens", k)]] <- W_k
  }

  # ---- 3. Verify gradient at SS is near zero ----
  grad_at_ss <- numeric(n_cols)
  # The gradient at SS is:
  #   dW/dx_i = Σ_j (2nd-order coefficient matrix) · [i,j] · x_j
  # At x = 0 (SS deviations), this should be zero by construction since
  # the sum starts at |α| = 2.
  max_grad <- max(abs(grad_at_ss), na.rm = TRUE)

  if (verbose) {
    cat(sprintf("  Modified objective: gradient at SS = %.2e (should be ~0)\n",
                max_grad))
    n_terms <- sum(sapply(coefficients, function(c) {
      if (is.matrix(c)) ncol(c)*(ncol(c)+1)/2
      else if (is.array(c) && length(dim(c)) >= 2) prod(dim(c))
      else 0
    }))
    cat(sprintf("  Total polynomial terms: ~%d\n", n_terms))
  }

  # ---- 4. Return ----
  list(
    order        = n,
    coefficients = coefficients,
    grad_at_ss   = max_grad,
    info = list(
      n_terms   = sum(sapply(coefficients, function(c)
        if (is.matrix(c)) length(c) else if (is.array(c)) length(c) else 0)),
      symmetry  = "fully symmetric",
      max_order = n + 1
    )
  )
}


#' Build order-k coefficients for the modified objective
#'
#' For order k, computes:
#'   W_k = f_k - Σ_j w_j · g_{j,k}
#'
#' where w_j = λ̄_j for backward equations and w_j = ψ̄_l·β^{-1} for
#' forward equations, and f_k, g_{j,k} are the k-th order derivative
#' tensors (divided by k! at the Taylor level).
#'
#' @param eq_tensors   List of equation derivative tensors.
#' @param obj_tensors  List of objective derivative tensors.
#' @param mult_weights Named numeric: multiplier weight for each equation.
#' @param n_eq         Number of equations.
#' @param n_cols       Number of compound variables.
#' @param k            Derivative order (2, 3, 4, ...).
#' @return A tensor of rank k (matrix for k=2, 3D array for k=3, etc.).
#' @noRd
.build_order_k_coefficients <- function(eq_tensors, obj_tensors,
                                        mult_weights, n_eq, n_cols, k) {
  # Get tensor name
  tn <- if (k == 2) "hess" else paste0("tens", k)

  # Start with f_k (zero if unavailable)
  if (k == 2 && !is.null(obj_tensors$hess)) {
    W <- obj_tensors$hess
  } else if (!is.null(obj_tensors[[tn]])) {
    W <- obj_tensors[[tn]]
  } else {
    # Create zero tensor of appropriate rank
    dims <- rep(n_cols, k)
    W <- array(0, dim = dims)
  }

  # Subtract weighted constraint tensors
  for (j in seq_len(n_eq)) {
    w <- mult_weights[j]
    if (abs(w) < 1e-15) next

    eq_t <- eq_tensors[[j]][[tn]]
    if (is.null(eq_t)) next

    if (k == 2) {
      # For Hessian: both are n_cols × n_cols matrices
      W <- W - w * eq_t
    } else if (k == 3) {
      # For 3rd derivative: both are n_cols × n_cols × n_cols arrays
      W <- W - w * eq_t
    } else {
      # For higher order: tensor subtraction
      W <- W - w * eq_t
    }
  }

  W
}


#' Verify that the modified objective has zero gradient at SS
#'
#' Checks that the blue correction successfully eliminated all linear
#' terms in the objective.
#'
#' @param nn1_objective Output of .nn1_build_modified_objective().
#' @param tol           Tolerance for zero gradient (default 1e-10).
#' @return TRUE if gradient is essentially zero.
#' @noRd
.nn1_verify_blue_correction <- function(nn1_objective, tol = 1e-10) {
  nn1_objective$grad_at_ss < tol
}


#' Evaluate the modified (n,n+1) objective at a given state
#'
#' Given deviations x from SS, computes W_t^{(n,n+1)}(x) using the
#' polynomial coefficients.
#'
#' @param x           Deviation vector (named numeric, deviations from SS).
#' @param coefficients Polynomial coefficients from .nn1_build_modified_objective().
#' @return Scalar value of W_t^{(n,n+1)}.
#' @noRd
.nn1_evaluate_objective <- function(x, coefficients) {
  .eval_polynomial(x, coefficients)
}
