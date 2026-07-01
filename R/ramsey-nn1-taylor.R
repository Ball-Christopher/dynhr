## R/ramsey-nn1-taylor.R
## --------------------------------------------------------------------------
## E2: Taylor expansion engine for the (n, n+1) approximation.
##
## Computes n-order Taylor expansions of the planner objective and each
## model constraint at the deterministic steady state, returning polynomial
## coefficients (tensors) suitable for constructing the modified objective
## W_t^{(n,n+1)} and the constraint approximations.
##
## Derivative sources:
##   Order 1 (Jacobian):  compiled$dynamic$jacobian_fn
##   Order 2 (Hessian):   compiled$dynamic$hessian2_fn + hess2_triplets
##   Order 3 (3rd deriv): compiled$dynamic$hessian3_fn + hess3_triplets
##   Order ≥ 4:           Recursive numerical finite differences (fallback)
##
## All tensors are stored in symmetric (hyper-triangular) form.
##
## References:
##   Gross, I. & Hansen, J. (2021). EER 140, 103918, Proposition 1.
##   Schmitt-Grohé & Uribe (2004). Solving DGGE Models Using a
##     Second-Order Approximation. JEDC 28(4), 755-775.
## --------------------------------------------------------------------------

#' Compute Taylor expansion coefficients for the (n, n+1) approximation
#'
#' Evaluates the derivatives of the planner objective and each model
#' equation up to the requested order at the deterministic steady state.
#'
#' @param compiled  dynhr_compiled with Jacobian, Hessian, etc.
#' @param ss        Named numeric steady state vector.
#' @param params    Named numeric parameter vector.
#' @param order     Maximum derivative order (1, 2, 3, ...).
#' @param obj_ast   Optional AST of the planner objective (from parse_expression).
#'                  If NULL, only equation derivatives are computed.
#' @param method    "symbolic" (preferred) or "numerical" (fallback).
#' @param h         Step size for numerical differentiation (default 1e-4).
#' @param verbose   Print progress messages.
#'
#' @return A list with:
#'   \item{objective}{List of derivative tensors for the planner objective:
#'         $grad (vector), $hess (matrix), $tens3 (3D array), $tens4 (4D array), ...}
#'   \item{equations}{List of n_eq elements, each with the same structure as $objective.}
#'   \item{n_cols}{Number of compound-variable columns (2m + p).}
#'   \item{col_names}{Character vector of column names (variable__leadlag).}
#'   \item{endo_names}{Endogenous variable names.}
#'   \item{exo_names}{Exogenous variable names.}
#'   \item{ss}{The steady state used.}
#' @noRd
.nn1_taylor_expand <- function(compiled,
                                ss,
                                params,
                                order,
                                obj_ast = NULL,
                                method = c("symbolic", "numerical"),
                                h = 1e-4,
                                verbose = FALSE) {
  method <- match.arg(method)
  dyn <- compiled$dynamic
  model <- compiled$model
  endo_names <- model$var_names
  exo_names <- model$varexo_names
  n_endo <- length(endo_names)
  n_exo <- length(exo_names)
  n_eq <- dyn$n_eq
  total_cols <- dyn$total_cols

  if (verbose) {
    cat(sprintf("[nn1_taylor] Expanding to order %d (%d cols, %d eqs)\n",
                order, total_cols, n_eq))
  }

  # ---- Build steady-state compound vector dy_ss ----
  dy_ss <- .build_dy_ss_nn1(compiled, ss)

  # ---- Column names for reference ----
  col_names <- .build_col_names(compiled)

  # ---- 1. Compute equation derivatives ----
  eq_tensors <- vector("list", n_eq)

  for (j in seq_len(n_eq)) {
    eq_tensors[[j]] <- list()
  }

  # Order 1: Jacobian (shared across all equations)
  J <- .compute_jacobian(compiled, dy_ss, params, ss)
  for (j in seq_len(n_eq)) {
    eq_tensors[[j]]$grad <- J[j, , drop = TRUE]
  }

  # Order 2: Hessian
  if (order >= 2) {
    H <- .compute_hessian_nn1(compiled, dy_ss, params, ss, h, method)
    if (!is.null(H)) {
      for (j in seq_len(n_eq)) {
        eq_tensors[[j]]$hess <- H[j, , , drop = TRUE]
      }
    }
  }

  # Order 3: 3rd derivative
  if (order >= 3) {
    T3 <- .compute_third_deriv_nn1(compiled, dy_ss, params, ss, h, method)
    if (!is.null(T3)) {
      for (j in seq_len(n_eq)) {
        eq_tensors[[j]]$tens3 <- T3[j, , , , drop = TRUE]
      }
    }
  }

  # Order 4+: Numerical FD
  if (order >= 4) {
    for (k in 4:order) {
      Tk <- .compute_kth_deriv_nn1(compiled, dy_ss, params, ss, k, h)
      if (!is.null(Tk)) {
        for (j in seq_len(n_eq)) {
          eq_tensors[[j]][[paste0("tens", k)]] <- Tk[j, , , , , drop = TRUE]
        }
      }
    }
  }

  # ---- 2. Compute objective derivatives (if obj_ast provided) ----
  obj_tensors <- list()

  if (!is.null(obj_ast)) {
    # Build a function that evaluates the objective at a given dy point
    obj_fn <- .make_objective_fn(obj_ast, endo_names, exo_names, params, ss)

    # Order 0: value at SS
    obj_tensors$value <- obj_fn(dy_ss)

    # Order 1: gradient (numerical)
    if (order >= 1) {
      obj_tensors$grad <- .numerical_gradient(obj_fn, dy_ss, h)
    }

    # Order 2: Hessian
    if (order >= 2) {
      obj_tensors$hess <- .numerical_hessian(obj_fn, dy_ss, h)
    }

    # Order 3+
    if (order >= 3) {
      for (k in 3:order) {
        tens_k <- .numerical_kth_tensor(
          obj_fn, dy_ss, n_out = 1L, k = k, h = h)
        # Squeeze out the first dimension (n_out=1 for scalar function)
        obj_tensors[[paste0("tens", k)]] <- drop(tens_k)
      }
    }
  }

  # ---- 3. Return ----
  list(
    objective  = obj_tensors,
    equations  = eq_tensors,
    n_cols     = total_cols,
    col_names  = col_names,
    endo_names = endo_names,
    exo_names  = exo_names,
    ss         = ss
  )
}


#' Build the compound variable vector at steady state
#'
#' Creates the named numeric vector used as input to Jacobian/Hessian
#' functions, mapping each (variable, lead_lag) pair to its SS value.
#'
#' @param compiled dynhr_compiled.
#' @param ss       Named numeric steady state.
#' @return Named numeric vector.
#' @noRd
.build_dy_ss_nn1 <- function(compiled, ss) {
  dyn <- compiled$dynamic
  model <- compiled$model
  endo <- model$var_names
  exo <- model$varexo_names

  dy <- numeric(nrow(dyn$dyn_col_map))
  keys <- character(nrow(dyn$dyn_col_map))

  for (k in seq_len(nrow(dyn$dyn_col_map))) {
    nm <- dyn$dyn_col_map$name[k]
    ll <- dyn$dyn_col_map$lead_lag[k]
    sfx <- if (ll == 0L) "__0"
    else if (ll > 0L) paste0("__p", ll)
    else paste0("__m", abs(ll))
    keys[k] <- paste0(nm, sfx)
    dy[k] <- if (nm %in% names(ss)) ss[[nm]] else 0
  }

  names(dy) <- keys

  # Append exogenous variables (they appear as additional columns in Jacobian)
  for (ex in exo) {
    dy[paste0(ex, "__0")] <- 0
  }

  dy
}


#' Build column names for the compound variable vector
#'
#' @param compiled dynhr_compiled.
#' @return Character vector of column names.
#' @noRd
.build_col_names <- function(compiled) {
  dyn <- compiled$dynamic
  model <- compiled$model
  endo <- model$var_names
  exo <- model$varexo_names

  names <- character(nrow(dyn$dyn_col_map))
  for (k in seq_len(nrow(dyn$dyn_col_map))) {
    nm <- dyn$dyn_col_map$name[k]
    ll <- dyn$dyn_col_map$lead_lag[k]
    sfx <- if (ll == 0L) "__0"
    else if (ll > 0L) paste0("__p", ll)
    else paste0("__m", abs(ll))
    names[k] <- paste0(nm, sfx)
  }

  for (ex in exo) {
    names <- c(names, paste0(ex, "__0"))
  }

  names
}


#' Compute the model Jacobian at steady state
#'
#' @param compiled dynhr_compiled.
#' @param dy_ss    Named compound vector at SS.
#' @param params   Named parameter vector.
#' @param ss       Named numeric steady state.
#' @return Matrix of size n_eq × total_cols.
#' @noRd
.compute_jacobian <- function(compiled, dy_ss, params, ss) {
  dyn <- compiled$dynamic

  J <- dyn$jacobian_fn(dy_ss, params, ss)

  if (is.null(J)) {
    # Fallback: numerical Jacobian
    fn <- function(x) {
      dyn$residuals_fn(x, rep(0, length(dy_ss)), params, ss)
    }
    if (requireNamespace("numDeriv", quietly = TRUE)) {
      J <- numDeriv::jacobian(fn, dy_ss)
    } else {
      # Simple finite-difference Jacobian
      n <- length(dy_ss)
      f0 <- fn(dy_ss)
      J <- matrix(0, nrow = length(f0), ncol = n)
      for (i in seq_len(n)) {
        xp <- dy_ss; xp[i] <- dy_ss[i] + h
        xm <- dy_ss; xm[i] <- dy_ss[i] - h
        J[, i] <- (fn(xp) - fn(xm)) / (2 * h)
      }
    }
  }

  # Handle any non-finite entries
  J[!is.finite(J)] <- 0

  J
}


#' Compute the model Hessian (2nd derivative tensor)
#'
#' Tries symbolic Hessian first, falls back to numerical.
#'
#' @param compiled dynhr_compiled.
#' @param dy_ss    Named compound vector at SS.
#' @param params   Named parameter vector.
#' @param ss       Named numeric steady state.
#' @param h        Step size for numerical FD.
#' @param method   "symbolic" or "numerical".
#' @return Array of size n_eq × total_cols × total_cols, or NULL.
#' @noRd
.compute_hessian_nn1 <- function(compiled, dy_ss, params, ss, h = 1e-4,
                                  method = "symbolic") {
  dyn <- compiled$dynamic

  # Try symbolic first
  if (method == "symbolic" && !is.null(dyn$hessian2_fn) &&
      !is.null(dyn$hess2_triplets)) {
    H <- .eval_symbolic_hessian(dyn, dy_ss, params, ss)
    if (!is.null(H)) return(H)
  }

  # Numerical fallback
  fn <- function(x) {
    dyn$residuals_fn(x, rep(0, length(dy_ss)), params, ss)
  }
  .numerical_hessian_tensor(fn, dy_ss, dyn$n_eq, h)
}


#' Evaluate the symbolic Hessian from the compiled model
#'
#' @param dyn     compiled$dynamic.
#' @param dy_ss   Named compound vector at SS.
#' @param params  Named parameter vector.
#' @param ss      Named numeric steady state.
#' @return Array of size n_eq × total_cols × total_cols, or NULL.
#' @noRd
.eval_symbolic_hessian <- function(dyn, dy_ss, params, ss) {
  n_eq <- dyn$n_eq
  total_cols <- dyn$total_cols
  n_hess <- dyn$n_hess %||% length(dyn$hess2_triplets)

  H <- array(0, dim = c(n_eq, total_cols, total_cols))
  if (n_hess == 0L) return(H)

  values <- dyn$hessian2_fn(dy_ss, params, ss)

  for (k in seq_len(n_hess)) {
    t <- dyn$hess2_triplets[[k]]
    v <- values[k]
    H[t$eq, t$col1, t$col2] <- v
    if (t$col1 != t$col2) H[t$eq, t$col2, t$col1] <- v
  }

  H
}


#' Compute the model 3rd derivative tensor
#'
#' Tries symbolic 3rd derivative first, falls back to numerical.
#'
#' @param compiled dynhr_compiled.
#' @param dy_ss    Named compound vector at SS.
#' @param params   Named parameter vector.
#' @param ss       Named numeric steady state.
#' @param h        Step size for numerical FD.
#' @param method   "symbolic" or "numerical".
#' @return Array of size n_eq × total_cols × total_cols × total_cols, or NULL.
#' @noRd
.compute_third_deriv_nn1 <- function(compiled, dy_ss, params, ss, h = 1e-4,
                                      method = "symbolic") {
  dyn <- compiled$dynamic

  # Try symbolic first
  if (method == "symbolic" && !is.null(dyn$hessian3_fn) &&
      !is.null(dyn$hess3_triplets)) {
    T3 <- .eval_symbolic_third_deriv(dyn, dy_ss, params, ss)
    if (!is.null(T3)) return(T3)
  }

  # Numerical fallback: differentiate the Hessian
  fn <- function(x) {
    dyn$residuals_fn(x, rep(0, length(dy_ss)), params, ss)
  }
  .numerical_third_tensor(fn, dy_ss, dyn$n_eq, h)
}


#' Evaluate the symbolic 3rd derivative from the compiled model
#'
#' @param dyn     compiled$dynamic.
#' @param dy_ss   Named compound vector at SS.
#' @param params  Named parameter vector.
#' @param ss      Named numeric steady state.
#' @return Array of size n_eq × total_cols × total_cols × total_cols, or NULL.
#' @noRd
.eval_symbolic_third_deriv <- function(dyn, dy_ss, params, ss) {
  n_eq <- dyn$n_eq
  total_cols <- dyn$total_cols
  n_hess3 <- dyn$n_hess3 %||% length(dyn$hess3_triplets)

  T3 <- array(0, dim = c(n_eq, total_cols, total_cols, total_cols))
  if (n_hess3 == 0L) return(T3)

  values <- dyn$hessian3_fn(dy_ss, params, ss)

  for (k in seq_len(n_hess3)) {
    t <- dyn$hess3_triplets[[k]]
    v <- values[k]
    T3[t$eq, t$col1, t$col2, t$col3] <- v
    # Fill symmetric permutations
    T3[t$eq, t$col1, t$col3, t$col2] <- v
    T3[t$eq, t$col2, t$col1, t$col3] <- v
    T3[t$eq, t$col2, t$col3, t$col1] <- v
    T3[t$eq, t$col3, t$col1, t$col2] <- v
    T3[t$eq, t$col3, t$col2, t$col1] <- v
  }

  T3
}


#' Compute k-th order derivative tensor via numerical FD
#'
#' Recursive finite-differences for k ≥ 4. Expensive (O(n_cols^k)).
#'
#' @param compiled dynhr_compiled.
#' @param dy_ss    Named compound vector at SS.
#' @param params   Named parameter vector.
#' @param ss       Named numeric steady state.
#' @param k        Derivative order (4, 5, ...).
#' @param h        Step size.
#' @return Array of size n_eq × total_cols^k, or NULL.
#' @noRd
.compute_kth_deriv_nn1 <- function(compiled, dy_ss, params, ss, k, h = 1e-4) {
  dyn <- compiled$dynamic
  fn <- function(x) {
    dyn$residuals_fn(x, rep(0, length(dy_ss)), params, ss)
  }
  .numerical_kth_tensor(fn, dy_ss, dyn$n_eq, k, h)
}


#' Numerical gradient of a scalar function
#'
#' @param fn  Function f(x) returning a scalar.
#' @param x0  Point at which to evaluate.
#' @param h   Step size.
#' @return Numeric gradient vector.
#' @noRd
.numerical_gradient <- function(fn, x0, h = 1e-4) {
  n <- length(x0)
  grad <- numeric(n)
  f0 <- fn(x0)
  for (i in seq_len(n)) {
    xp <- x0; xp[i] <- x0[i] + h
    xm <- x0; xm[i] <- x0[i] - h
    grad[i] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  grad
}


#' Numerical Hessian of a scalar function
#'
#' @param fn  Function f(x) returning a scalar.
#' @param x0  Point at which to evaluate.
#' @param h   Step size.
#' @return Numeric Hessian matrix.
#' @noRd
.numerical_hessian <- function(fn, x0, h = 1e-4) {
  n <- length(x0)
  H <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in i:n) {
      xpp <- x0; xpp[i] <- x0[i] + h; xpp[j] <- xpp[j] + h
      xpm <- x0; xpm[i] <- x0[i] + h; xpm[j] <- xpm[j] - h
      xmp <- x0; xmp[i] <- x0[i] - h; xmp[j] <- xmp[j] + h
      xmm <- x0; xmm[i] <- x0[i] - h; xmm[j] <- xmm[j] - h
      H[i, j] <- (fn(xpp) - fn(xpm) - fn(xmp) + fn(xmm)) / (4 * h^2)
      if (i != j) H[j, i] <- H[i, j]
    }
  }
  H
}


#' Numerical Hessian tensor for a vector-valued function
#'
#' @param fn  Function f(x) returning a vector.
#' @param x0  Point at which to evaluate.
#' @param n_out  Number of output dimensions.
#' @param h   Step size.
#' @return Array of size n_out × n × n.
#' @noRd
.numerical_hessian_tensor <- function(fn, x0, n_out, h = 1e-4) {
  n <- length(x0)
  H <- array(0, dim = c(n_out, n, n))
  f0 <- fn(x0)

  for (i in seq_len(n)) {
    for (j in i:n) {
      xpp <- x0; xpp[i] <- x0[i] + h; xpp[j] <- xpp[j] + h
      xpm <- x0; xpm[i] <- x0[i] + h; xpm[j] <- xpm[j] - h
      xmp <- x0; xmp[i] <- x0[i] - h; xmp[j] <- xmp[j] + h
      xmm <- x0; xmm[i] <- x0[i] - h; xmm[j] <- xmm[j] - h
      H[, i, j] <- (fn(xpp) - fn(xpm) - fn(xmp) + fn(xmm)) / (4 * h^2)
      if (i != j) H[, j, i] <- H[, i, j]
    }
  }

  H
}


#' Numerical 3rd derivative tensor for a vector-valued function
#'
#' @param fn     Function f(x) returning a vector.
#' @param x0     Point at which to evaluate.
#' @param n_out  Number of output dimensions.
#' @param h      Step size.
#' @return Array of size n_out × n × n × n.
#' @noRd
.numerical_third_tensor <- function(fn, x0, n_out, h = 1e-4) {
  n <- length(x0)
  T3 <- array(0, dim = c(n_out, n, n, n))

  for (i in seq_len(n)) {
    for (j in i:n) {
      for (k in j:n) {
        xppp <- x0; xppp[c(i,j,k)] <- x0[c(i,j,k)] + h
        xppm <- x0; xppm[c(i,j)] <- x0[c(i,j)] + h; xppm[k] <- x0[k] - h
        xpmp <- x0; xpmp[i] <- x0[i] + h; xpmp[j] <- x0[j] - h; xpmp[k] <- x0[k] + h
        xpmm <- x0; xpmm[i] <- x0[i] + h; xpmm[j] <- x0[j] - h; xpmm[k] <- x0[k] - h
        xmpp <- x0; xmpp[i] <- x0[i] - h; xmpp[j] <- x0[j] + h; xmpp[k] <- x0[k] + h
        xmpm <- x0; xmpm[i] <- x0[i] - h; xmpm[j] <- x0[j] + h; xmpm[k] <- x0[k] - h
        xmmp <- x0; xmmp[i] <- x0[i] - h; xmmp[j] <- x0[j] - h; xmmp[k] <- x0[k] + h
        xmmm <- x0; xmmm[c(i,j,k)] <- x0[c(i,j,k)] - h

        val <- (fn(xppp) - fn(xppm) - fn(xpmp) + fn(xpmm)
               - fn(xmpp) + fn(xmpm) + fn(xmmp) - fn(xmmm)) / (8 * h^3)

        # Fill all unique permutations of (i,j,k)
        if (i == j && j == k) {
          T3[, i, j, k] <- val
        } else if (i == j) {
          T3[, i, i, k] <- val; T3[, i, k, i] <- val; T3[, k, i, i] <- val
        } else if (j == k) {
          T3[, i, j, j] <- val; T3[, j, i, j] <- val; T3[, j, j, i] <- val
        } else {
          T3[, i, j, k] <- val; T3[, i, k, j] <- val; T3[, j, i, k] <- val
          T3[, j, k, i] <- val; T3[, k, i, j] <- val; T3[, k, j, i] <- val
        }
      }
    }
  }

  T3
}


#' Numerical k-th order derivative tensor for a vector-valued function
#'
#' General recursive FD for arbitrary order k. Uses central differences
#' on the (k-1)-order derivative. Very expensive for large n and k.
#'
#' @param fn     Function f(x) returning a vector.
#' @param x0     Point at which to evaluate.
#' @param n_out  Number of output dimensions.
#' @param k      Derivative order.
#' @param h      Step size.
#' @return Array of size n_out × n^k.
#' @noRd
.numerical_kth_tensor <- function(fn, x0, n_out, k, h = 1e-4) {
  if (k < 1) stop("k must be >= 1")
  if (k == 1) {
    # Jacobian
    n <- length(x0)
    J <- matrix(0, n_out, n)
    for (i in seq_len(n)) {
      xp <- x0; xp[i] <- x0[i] + h
      xm <- x0; xm[i] <- x0[i] - h
      J[, i] <- (fn(xp) - fn(xm)) / (2 * h)
    }
    return(J)
  }

  if (k == 2) {
    return(.numerical_hessian_tensor(fn, x0, n_out, h))
  }

  if (k == 3) {
    return(.numerical_third_tensor(fn, x0, n_out, h))
  }

  # For k >= 4, recursively difference the (k-1)-order tensor
  # This computes ∂^k f / (∂x_{i1} ... ∂x_{ik}) via:
  #   D^k f ≈ (D^{k-1}f(x + h·e_{ik}) - D^{k-1}f(x - h·e_{ik})) / (2h)
  n <- length(x0)
  T <- array(0, dim = c(n_out, rep(n, k)))

  for (ik in seq_len(n)) {
    xp <- x0; xp[ik] <- x0[ik] + h
    xm <- x0; xm[ik] <- x0[ik] - h

    # Compute D^{k-1}f at xp and xm (recursive)
    Tp <- .numerical_kth_tensor(fn, xp, n_out, k - 1, h)
    Tm <- .numerical_kth_tensor(fn, xm, n_out, k - 1, h)

    idx <- rep(list(TRUE), k)
    idx[[k]] <- ik
    T[,,, drop = FALSE] <- (Tp - Tm) / (2 * h)
  }

  T
}


#' Build a function that evaluates the planner objective at a dy point
#'
#' Creates a closure that maps the compound vector dy → scalar objective value.
#' The dy vector uses "__0", "__m1", "__p1" suffixed keys as expected by ast_eval.
#'
#' @param obj_ast    Parsed AST of the planner objective.
#' @param endo_names Endogenous variable names.
#' @param exo_names  Exogenous variable names.
#' @param params     Named parameter vector.
#' @param ss         Named numeric steady state.
#' @return Function f(dy) returning a scalar.
#' @noRd
.make_objective_fn <- function(obj_ast, endo_names, exo_names, params, ss) {
  function(dy) {
    # ast_eval expects var_values with __0, __m1, __p1 suffixes
    # Extract current-period values from dy and build full var_values
    var_values <- numeric(0)
    for (nm in endo_names) {
      sv <- if (!is.null(dy[[paste0(nm, "__0")]])) dy[[paste0(nm, "__0")]] else 0
      var_values[[paste0(nm, "__m1")]] <- sv
      var_values[[paste0(nm, "__0")]]  <- sv
      var_values[[paste0(nm, "__p1")]] <- sv
    }
    for (nm in exo_names) {
      var_values[[paste0(nm, "__0")]] <- if (!is.null(dy[[paste0(nm, "__0")]])) dy[[paste0(nm, "__0")]] else 0
    }
    ast_eval(obj_ast, var_values = var_values, param_values = params, ss_values = ss)
  }
}


#' Create an objective function from a character expression
#'
#' @param objective_text Character string: the planner objective.
#' @param model          dynhr_mod.
#' @param params         Named parameter vector.
#' @param ss             Named numeric steady state.
#' @return Function f(x) returning a scalar for a given compound vector x.
#' @noRd
.make_objective_from_text <- function(objective_text, model, params, ss) {
  all_var_names <- c(model$var_names, model$varexo_names, model$varexo_det_names)
  obj_ast <- parse_expression(objective_text,
                              var_names = all_var_names,
                              param_names = model$param_names)
  .make_objective_fn(obj_ast, model$var_names, model$varexo_names, params, ss)
}


#' Check if higher-order derivatives are available in the compiled model
#'
#' @param compiled dynhr_compiled.
#' @param k         Derivative order to check (4, 5, ...).
#' @return TRUE if symbolic k-th derivatives are available.
#' @noRd
.has_higher_order_derivatives <- function(compiled, k) {
  fn_name <- paste0("hessian", k, "_fn")
  !is.null(compiled$dynamic[[fn_name]])
}


#' Check if higher-order perturbation solver is available
#'
#' @param k Perturbation order (4, 5, ...).
#' @return TRUE if solve_perturbation_order{k} exists.
#' @noRd
.has_higher_order_solver <- function(k) {
  fn_name <- paste0("solve_perturbation_order", k)
  exists(fn_name, mode = "function")
}


#' Evaluate a polynomial given tensor coefficients and deviation vector
#'
#' Computes Σ_{|α|=2}^{order} c_α · x^{α} / α! where c_α are the tensor
#' coefficients stored by order in the `coefficients` list.
#'
#' @param x            Deviations from SS (named numeric vector).
#' @param coefficients List with elements $quad (matrix), $cubic (3D array),
#'                     $quartic (4D array), etc. May be NULL for missing orders.
#' @return Scalar polynomial value.
#' @noRd
.eval_polynomial <- function(x, coefficients) {
  val <- 0

  # Quadratic terms: x' · Q · x / 2!
  if (!is.null(coefficients$quad)) {
    Q <- coefficients$quad
    val <- val + 0.5 * as.numeric(t(x) %*% Q %*% x)
  }

  # Cubic terms: Σ_{ijk} T3[i,j,k] · x[i] · x[j] · x[k] / 3!
  if (!is.null(coefficients$cubic)) {
    T3 <- coefficients$cubic
    n <- length(x)
    cubic_val <- 0
    for (i in seq_len(n)) {
      for (j in i:n) {
        for (k in j:n) {
          w <- if (i == j && j == k) 1
          else if (i == j || j == k || i == k) 3
          else 6
          cubic_val <- cubic_val + w * T3[i, j, k] * x[i] * x[j] * x[k]
        }
      }
    }
    val <- val + cubic_val / 6  # divide by 3!
  }

  # Quartic terms
  if (!is.null(coefficients$quartic)) {
    T4 <- coefficients$quartic
    n <- length(x)
    quartic_val <- 0
    for (i in seq_len(n)) {
      for (j in i:n) {
        for (k in j:n) {
          for (l in k:n) {
            # Count distinct permutations
            counts <- table(c(i, j, k, l))
            w <- factorial(4) / prod(factorial(counts))
            quartic_val <- quartic_val + w * T4[i, j, k, l] * x[i] * x[j] * x[k] * x[l]
          }
        }
      }
    }
    val <- val + quartic_val / 24  # divide by 4!
  }

  val
}
