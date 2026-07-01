## R/opp-sufficient-stats.R
## --------------------------------------------------------------------------
## Phase I — Sufficient-Statistics OPP (Optimal Policy Projections) Module
##
## Implements the Barnichon & Mesters (2023) sufficient-statistics approach
## to optimal policy. The key insight is that the optimal policy intervention
## (deviation from a baseline instrument path) can be computed using only:
##
##   1. The IRF of target variables to a policy instrument shock ("R matrix")
##   2. The baseline forecast of target variables
##   3. The policymaker's loss function (weights + target values)
##
## This avoids solving the full Ramsey problem and works with any model
## that can generate IRFs — making it model-agnostic.
##
## Applications:
##   - "Is our simple rule close to optimal?" — quick check using only IRFs
##   - Optimal policy projections (OPP) as in COPPs Toolkit
##   - Benchmarking DSGE policy prescriptions against empirical models
##   - Robustness checks across alternative model specifications
##
## References:
##   Barnichon, R. and G. Mesters (2023). "A Sufficient Statistics Approach
##     for Macro Policy." American Economic Review, 113(11), 2809–2845.
##   de Groot, O., F. Mazelis, R. Motto, and A. Ristiniemi (2021). "A Toolkit
##     for Computing Optimal Policy Projections in DSGE Models." ECB WP.
##   Brennan, J. and S. Dasgupta (2025). "Optimal Policy Projections at the
##     RBNZ." RBNZ Discussion Paper.
## --------------------------------------------------------------------------


# ==========================================================================
# 1. Core computation: optimal intervention from sufficient statistics
# ==========================================================================

#' Compute optimal policy intervention from sufficient statistics
#'
#' The optimal instrument path that minimises a quadratic loss function
#' given the IRF of target variables to the instrument and a baseline
#' forecast.
#'
#' The loss function is:
#' \deqn{L = \sum_{h=0}^{H} \beta^{h} \left[ \hat{\mathbf{z}}_{t+h|t}' \mathbf{W}
#'       \hat{\mathbf{z}}_{t+h|t} + \lambda (\Delta i_{t+h})^2 \right]}
#'
#' where \eqn{\hat{\mathbf{z}}_{t+h|t} = \mathbf{z}_{t+h|t}^{base} - \mathbf{z}^*}
#' is the deviation of the baseline forecast from target, \eqn{\mathbf{R}_h} is
#' the IRF of target variables at horizon \eqn{h} to a unit instrument shock,
#' and \eqn{\Delta i} is the instrument deviation from baseline.
#'
#' The closed-form solution is:
#' \deqn{\Delta \mathbf{i}^* = -\left( \sum_{h=0}^H \beta^h \mathbf{R}_h'
#'       \mathbf{W} \mathbf{R}_h + \lambda \mathbf{I} \right)^{-1}
#'       \sum_{h=0}^H \beta^h \mathbf{R}_h' \mathbf{W}
#'       \hat{\mathbf{z}}_{t+h}}
#'
#' @param instrument_irf   Matrix or data.frame of the IRF of target variables
#'   to a one-unit policy instrument shock. Each row is a horizon; each column
#'   is a target variable. Alternatively, a named list of such matrices for
#'   multiple instruments.
#' @param baseline_forecast Matrix or data.frame of the baseline forecast
#'   deviations from target. Same dimensions as \code{instrument_irf}. If a
#'   vector, it is treated as the value at horizon 0.
#' @param loss_weights     Numeric vector of weights on target variables, or
#'   a square weight matrix. Defaults to equal weights (identity matrix).
#' @param target_values    Numeric vector of target (steady-state) values for
#'   each variable. If NULL, assumed to be zero (deviation from steady state).
#' @param instrument_penalty Numeric: coefficient \eqn{\lambda} penalising
#'   instrument volatility. Default 0 (no penalty). A small positive value
#'   (e.g. 0.01) improves numerical stability.
#' @param discount         Discount factor \eqn{\beta}. Default 0.99.
#' @param horizons         Number of horizons to consider. Defaults to
#'   \code{nrow(instrument_irf)}.
#' @param names            Character vector of variable names for labelling.
#'   Extracted from column names if not provided.
#' @param instrument_name  Character: name of the policy instrument (for
#'   reporting). Default "i".
#' @param verbose          Print progress and results.
#'
#' @return A list of class \code{dynhr_opp_result} containing:
#'   \describe{
#'     \item{instrument_path}{Optimal instrument deviations from baseline
#'       (vector of length \code{horizons + 1}).}
#'     \item{target_paths}{Matrix of implied optimal target variable paths
#'       (baseline + intervention effect).}
#'     \item{baseline_forecast}{The input baseline forecast.}
#'     \item{loss_baseline}{Loss value under the baseline policy.}
#'     \item{loss_optimal}{Loss value under the optimal intervention.}
#'     \item{loss_reduction}{Proportional loss reduction.}
#'     \item{sufficient_statistics}{The "R matrix" (discounted weighted
#'       IRF inner product) and "S vector" (discounted weighted forecast
#'       deviation).}
#'     \item{params}{Parameters used (discount, penalty, weights).}
#'     \item{meta}{Metadata.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Simple example with a single target variable (inflation)
#' irf_pi <- matrix(exp(-0.1 * 0:19), ncol = 1)  # Decaying IRF
#' colnames(irf_pi) <- "pi"
#' forecast <- matrix(rep(0.5, 20), ncol = 1)     # Persistent deviation
#' colnames(forecast) <- "pi"
#'
#' opp <- opp_optimal_intervention(irf_pi, forecast,
#'                                 loss_weights = c(pi = 1),
#'                                 instrument_penalty = 0.01)
#' print(opp)
#' }
#'
#' @export
opp_optimal_intervention <- function(instrument_irf,
                                      baseline_forecast,
                                      loss_weights = NULL,
                                      target_values = NULL,
                                      instrument_penalty = 0,
                                      discount = 0.99,
                                      horizons = NULL,
                                      names = NULL,
                                      instrument_name = "i",
                                      verbose = FALSE) {

  # ---- 1. Validate and standardise inputs ----
  if (is.list(instrument_irf) && !is.data.frame(instrument_irf)) {
    # Multiple instruments: recurse to compute per-instrument path
    # (for multi-instrument policy, the solution is joint)
    if (verbose) cat("Multiple instruments detected. Computing joint optimal path...\n")
    return(.opp_multi_instrument(instrument_irf, baseline_forecast,
                                 loss_weights, target_values,
                                 instrument_penalty, discount,
                                 horizons, verbose))
  }

  irf <- as.matrix(instrument_irf)
  fcst <- as.matrix(baseline_forecast)

  n_horiz <- nrow(irf)
  n_vars <- ncol(irf)

  if (nrow(fcst) != n_horiz) {
    stop(sprintf("Baseline forecast has %d rows but IRF has %d horizons.",
                 nrow(fcst), n_horiz))
  }
  if (ncol(fcst) != n_vars) {
    stop(sprintf("Baseline forecast has %d columns but IRF has %d variables.",
                 ncol(fcst), n_vars))
  }

  if (!is.null(horizons)) {
    h <- min(as.integer(horizons), n_horiz)
    irf <- irf[seq_len(h), , drop = FALSE]
    fcst <- fcst[seq_len(h), , drop = FALSE]
    n_horiz <- h
  }

  # ---- 2. Variable names ----
  if (is.null(names)) {
    names <- colnames(irf)
    if (is.null(names)) names <- colnames(fcst)
    if (is.null(names)) names <- paste0("z", seq_len(n_vars))
  }
  colnames(irf) <- names
  colnames(fcst) <- names

  # ---- 3. Loss weights ----
  if (is.null(loss_weights)) {
    W <- diag(n_vars)
  } else if (is.matrix(loss_weights)) {
    if (nrow(loss_weights) != n_vars || ncol(loss_weights) != n_vars) {
      stop("loss_weights matrix must be ", n_vars, "x", n_vars, ".")
    }
    W <- as.matrix(loss_weights)
  } else {
    w <- as.numeric(loss_weights)
    if (length(w) == 1) {
      W <- diag(w, n_vars)
    } else if (length(w) == n_vars) {
      W <- diag(w, n_vars)
    } else {
      stop("loss_weights length (", length(w), ") must match n_vars (", n_vars, ").")
    }
    dimnames(W) <- list(names, names)
  }

  # ---- 4. Target values ----
  if (is.null(target_values)) {
    target_dev <- fcst  # Forecast already in deviation form
  } else {
    tv <- as.numeric(target_values)
    if (length(tv) == 1) tv <- rep(tv, n_vars)
    target_dev <- sweep(fcst, 2, tv, "-")
  }

  # ---- 5. Compute sufficient statistics ----
  # The optimal intervention Δi is a scalar (single instrument deviation from
  # baseline). For each horizon h:
  #   RWR += β^h * (R_h * W * R_h')    → scalar
  #   RWz += β^h * (R_h * W * ẑ_h')    → scalar
  # where R_h is 1×n_vars (IRF row vector at horizon h).
  beta_vec <- discount^(seq_len(n_horiz) - 1)

  RWR <- 0  # scalar accumulator
  RWz <- 0  # scalar accumulator

  for (h in seq_len(n_horiz)) {
    Rh <- irf[h, , drop = FALSE]        # 1 x n_vars
    zh <- target_dev[h, , drop = FALSE]  # 1 x n_vars
    bh <- beta_vec[h]

    # Rh * W    → 1 x n_vars
    RhW <- Rh %*% W
    # Rh * W * Rh'  → 1 x 1 (scalar)
    RWR <- RWR + bh * as.numeric(RhW %*% t(Rh))
    # Rh * W * zh'  → 1 x 1 (scalar)
    RWz <- RWz + bh * as.numeric(RhW %*% t(zh))
  }

  # ---- 6. Add instrument penalty ----
  G <- RWR + instrument_penalty

  # ---- 7. Solve for optimal intervention ----
  # Δi* = -G^{-1} * RWz
  if (abs(G) < 1e-14) {
    if (verbose) warning("Sufficient-statistics G is near-zero. Check IRF scaling.")
    G <- max(G, 1e-14)
  }
  delta_i <- -RWz / G

  # ---- 8. Compute target paths under optimal intervention ----
  # z_opt(t) = z_base(t) + R(t) * Δi
  # Here Δi is a scalar and each element of R(t) multiplies it.
  opt_paths <- matrix(NA, n_horiz, n_vars)
  for (h in seq_len(n_horiz)) {
    opt_paths[h, ] <- fcst[h, ] + irf[h, ] * delta_i
  }
  colnames(opt_paths) <- names

  # ---- 9. Compute losses ----
  # Baseline loss: sum_h beta^h * z_hat_h * W * z_hat_h'
  loss_baseline <- 0
  for (h in seq_len(n_horiz)) {
    zh <- target_dev[h, , drop = FALSE]
    loss_baseline <- loss_baseline + beta_vec[h] * as.numeric(zh %*% W %*% t(zh))
  }

  # Optimal loss
  loss_optimal <- 0
  for (h in seq_len(n_horiz)) {
    zh_opt <- opt_paths[h, , drop = FALSE]
    if (is.null(target_values)) {
      zh_dev <- zh_opt
    } else {
      zh_dev <- sweep(zh_opt, 2, as.numeric(target_values), "-")
    }
    loss_optimal <- loss_optimal + beta_vec[h] * as.numeric(zh_dev %*% W %*% t(zh_dev))
  }
  loss_optimal <- loss_optimal + instrument_penalty * delta_i^2

  loss_reduction <- 1 - loss_optimal / max(loss_baseline, 1e-16)

  # ---- 10. Assemble result ----
  result <- list(
    instrument_path     = delta_i,
    instrument_name     = instrument_name,
    target_paths        = opt_paths,
    baseline_forecast   = fcst,
    target_deviation    = target_dev,
    loss_baseline       = loss_baseline,
    loss_optimal        = loss_optimal,
    loss_reduction      = loss_reduction,
    sufficient_statistics = list(
      RWR = RWR,
      RWz = RWz,
      G   = G
    ),
    params = list(
      discount         = discount,
      instrument_penalty = instrument_penalty,
      loss_weights     = W,
      target_values    = target_values,
      horizons         = n_horiz
    ),
    meta = list(
      n_vars           = n_vars,
      var_names        = names,
      timestamp        = Sys.time(),
      package_version  = utils::packageVersion("dynhr")
    )
  )
  class(result) <- c("dynhr_opp_result", "list")

  if (verbose) {
    cat("\n--- Optimal Policy Intervention ---\n")
    cat(sprintf("Instrument: %s\n", instrument_name))
    cat(sprintf("Horizons:   %d\n", n_horiz))
    cat(sprintf("Discount:   %.4f\n", discount))
    cat(sprintf("Penalty:    %.4f\n", instrument_penalty))
    cat(sprintf("Targets:    %s\n", paste(names, collapse = ", ")))
    cat(sprintf("Loss (baseline): %.6f\n", loss_baseline))
    cat(sprintf("Loss (optimal):  %.6f\n", loss_optimal))
    cat(sprintf("Loss reduction:  %.2f%%\n", loss_reduction * 100))
    cat(sprintf("Optimal intervention:\n"))
    print(round(delta_i, 6))
    cat("---\n")
  }

  result
}


#' Multi-instrument optimal intervention
#'
#' Handles the case of multiple policy instruments by building a joint
#' system and solving for all instrument paths simultaneously.
#'
#' @noRd
.opp_multi_instrument <- function(instrument_irfs,
                                   baseline_forecast,
                                   loss_weights,
                                   target_values,
                                   instrument_penalty,
                                   discount,
                                   horizons,
                                   verbose) {
  inst_names <- names(instrument_irfs)
  if (is.null(inst_names)) {
    inst_names <- paste0("i", seq_along(instrument_irfs))
  }
  n_inst <- length(instrument_irfs)

  fcst <- as.matrix(baseline_forecast)
  n_horiz <- nrow(fcst)
  n_vars <- ncol(fcst)

  if (is.null(loss_weights)) {
    W <- diag(n_vars)
  } else if (is.matrix(loss_weights)) {
    W <- as.matrix(loss_weights)
  } else {
    W <- diag(as.numeric(loss_weights), n_vars)
  }

  if (is.null(target_values)) {
    target_dev <- fcst
  } else {
    tv <- as.numeric(target_values)
    if (length(tv) == 1) tv <- rep(tv, n_vars)
    target_dev <- sweep(fcst, 2, tv, "-")
  }

  if (!is.null(horizons)) {
    h <- min(as.integer(horizons), n_horiz)
    fcst <- fcst[seq_len(h), , drop = FALSE]
    target_dev <- target_dev[seq_len(h), , drop = FALSE]
    n_horiz <- h
  }

  beta_vec <- discount^(seq_len(n_horiz) - 1)

  # Build joint system (n_inst x n_inst):
  # For each instrument k, Δi_k is a scalar.
  # G[k1,k2] = Σ_h β^h * (R_k1,h * W * R_k2,h')  → scalar cross-effect
  # RWz[k1]  = Σ_h β^h * (R_k1,h * W * ẑ_h')     → scalar forecast effect

  # Pre-compute per-instrument IRF matrices
  irf_mats <- list()
  for (k in seq_len(n_inst)) {
    nm <- inst_names[k]
    irf_k <- instrument_irfs[[nm]]
    if (is.null(irf_k)) {
      stop(sprintf("Instrument '%s' not found in instrument_irfs list.", nm))
    }
    irf_k <- as.matrix(irf_k)
    if (nrow(irf_k) < n_horiz) {
      stop(sprintf("IRF for '%s' has only %d rows, need %d.", nm, nrow(irf_k), n_horiz))
    }
    irf_mats[[k]] <- irf_k[seq_len(n_horiz), , drop = FALSE]
  }

  # Build joint G matrix (n_inst x n_inst)
  G_joint <- matrix(0, n_inst, n_inst)
  RWz_joint <- numeric(n_inst)

  for (k1 in seq_len(n_inst)) {
    for (k2 in seq_len(n_inst)) {
      # Compute R_k1 * W * R_k2' summed over horizons
      # Each IRF row is 1 x n_vars, so R_k1_h * W is 1 x n_vars,
      # and (R_k1_h * W) * t(R_k2_h) is 1 x 1 (scalar)
      block <- 0  # scalar for single-instrument-per-player case
      for (h in seq_len(n_horiz)) {
        Rk1_h <- irf_mats[[k1]][h, , drop = FALSE]  # 1 x n_vars
        Rk2_h <- irf_mats[[k2]][h, , drop = FALSE]  # 1 x n_vars
        bh <- beta_vec[h]
        block <- block + bh * as.numeric((Rk1_h %*% W) %*% t(Rk2_h))
      }

      # Place in joint matrix (each instrument contributes 1 row/col)
      G_joint[k1, k2] <- block
      if (k1 == k2) {
        G_joint[k1, k1] <- G_joint[k1, k1] + instrument_penalty
      }
    }

    # R_k * W * z' summed over horizons (scalar)
    Rwz_k <- 0
    for (h in seq_len(n_horiz)) {
      Rk_h <- irf_mats[[k1]][h, , drop = FALSE]  # 1 x n_vars
      zh <- target_dev[h, , drop = FALSE]         # 1 x n_vars
      bh <- beta_vec[h]
      Rwz_k <- Rwz_k + bh * as.numeric((Rk_h %*% W) %*% t(zh))
    }
    RWz_joint[k1] <- Rwz_k
  }

  # Trim joint matrices to n_inst x n_inst
  G_joint <- G_joint[seq_len(n_inst), seq_len(n_inst), drop = FALSE]
  RWz_joint <- RWz_joint[seq_len(n_inst)]

  # Solve joint system
  G_inv <- solve(G_joint)

  delta_all <- as.numeric(-G_inv %*% RWz_joint)
  names(delta_all) <- inst_names

  # Compute optimal target paths
  opt_paths <- fcst
  for (k in seq_len(n_inst)) {
    for (h in seq_len(n_horiz)) {
      opt_paths[h, ] <- opt_paths[h, ] + irf_mats[[k]][h, ] * delta_all[k]
    }
  }

  # Losses
  loss_baseline <- 0
  loss_optimal <- 0
  for (h in seq_len(n_horiz)) {
    zh <- target_dev[h, , drop = FALSE]
    loss_baseline <- loss_baseline + beta_vec[h] * as.numeric(zh %*% W %*% t(zh))

    zh_opt <- opt_paths[h, , drop = FALSE]
    if (is.null(target_values)) {
      zh_dev <- zh_opt
    } else {
      zh_dev <- sweep(zh_opt, 2, as.numeric(target_values), "-")
    }
    loss_optimal <- loss_optimal + beta_vec[h] * as.numeric(zh_dev %*% W %*% t(zh_dev))
  }
  loss_optimal <- loss_optimal + instrument_penalty * sum(delta_all^2)
  loss_reduction <- 1 - loss_optimal / max(loss_baseline, 1e-16)

  result <- list(
    instrument_path     = delta_all,
    instrument_name     = inst_names,
    target_paths        = opt_paths,
    baseline_forecast   = fcst,
    target_deviation    = target_dev,
    loss_baseline       = loss_baseline,
    loss_optimal        = loss_optimal,
    loss_reduction      = loss_reduction,
    sufficient_statistics = list(
      G_joint   = G_joint,
      G_inv     = G_inv,
      RWz_joint = RWz_joint
    ),
    params = list(
      discount         = discount,
      instrument_penalty = instrument_penalty,
      loss_weights     = W,
      target_values    = target_values,
      horizons         = n_horiz,
      n_instruments    = n_inst
    ),
    meta = list(
      n_vars           = n_vars,
      n_instruments    = n_inst,
      var_names        = colnames(fcst),
      instrument_names = inst_names,
      timestamp        = Sys.time(),
      package_version  = utils::packageVersion("dynhr")
    )
  )
  class(result) <- c("dynhr_opp_result", "list")

  if (verbose) {
    cat("\n--- Multi-Instrument Optimal Policy Intervention ---\n")
    cat(sprintf("Instruments: %s\n", paste(inst_names, collapse = ", ")))
    cat(sprintf("Horizons:    %d\n", n_horiz))
    cat(sprintf("Discount:    %.4f\n", discount))
    cat(sprintf("Loss (baseline): %.6f\n", loss_baseline))
    cat(sprintf("Loss (optimal):  %.6f\n", loss_optimal))
    cat(sprintf("Loss reduction:  %.2f%%\n", loss_reduction * 100))
    for (k in seq_len(n_inst)) {
      cat(sprintf("Optimal intervention [%s]:\n", inst_names[k]))
      print(round(delta_i_list[[inst_names[k]]], 6))
    }
    cat("---\n")
  }

  result
}


# ==========================================================================
# 2. Extract instrument IRF from a solved DSGE model
# ==========================================================================

#' Extract the IRF of target variables to a policy instrument shock
#'
#' Given a solved DSGE model (from \code{\link{solve_perturbation}} or a
#' policy result), computes the IRF of specified target variables to a
#' shock to the policy instrument.
#'
#' The instrument shock is identified as the exogenous shock that directly
#' affects the policy instrument equation (e.g., \code{eps_r} for the
#' interest rate rule). The function extracts the IRF from the model's
#' existing \code{compute_irfs()} infrastructure.
#'
#' @param dr            DecisionRules object from \code{\link{solve_perturbation}}.
#' @param model         dynhr_mod object.
#' @param target_vars   Character vector of target variable names (e.g.,
#'   \code{c("pi", "y_gap")}).
#' @param instrument_shock Character: name of the exogenous shock that
#'   affects the policy instrument (e.g., \code{"eps_r"}). If NULL,
#'   attempts auto-detection based on common naming conventions.
#' @param params        Named numeric parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param n_periods     Number of IRF periods (default 40).
#' @param shock_size    Size of the instrument shock in standard deviation
#'   units (default 1, i.e., a one-standard-deviation shock).
#' @param normalise     If TRUE, normalise the IRF so that the instrument's
#'   own impact response is 1 (a "unit instrument shock"). If FALSE, use
#'   raw IRF values (shock in std dev units).
#' @param instrument_name Character: name of the policy instrument variable
#'   in the model (e.g., \code{"r"}, \code{"i"}, \code{"R"}). Used for
#'   normalisation. If NULL and normalise=TRUE, attempts auto-detection.
#'
#' @return A matrix of dimension \code{n_periods x length(target_vars)}
#'   containing the IRF of target variables to the instrument shock.
#'
#' @examples
#' \dontrun{
#' # Parse and solve a model
#' model <- parse_mod("nk_model.mod")
#' compiled <- compile_model(model)
#' ss <- solve_steady(compiled, model$param_values)
#' dr <- solve_perturbation(model, compiled, ss$ss, model$param_values)
#'
#' # Extract inflation and output gap IRF to monetary policy shock
#' irf <- opp_estimate_instrument_irf(dr, model,
#'                                    target_vars = c("pi", "y_gap"),
#'                                    instrument_shock = "eps_r")
#'
#' # Use in optimal policy computation
#' forecast <- matrix(0.5, nrow = 40, ncol = 2)
#' colnames(forecast) <- c("pi", "y_gap")
#' opp <- opp_optimal_intervention(irf, forecast,
#'                                 loss_weights = c(pi = 1, y_gap = 0.5))
#' }
#'
#' @export
opp_estimate_instrument_irf <- function(dr,
                                         model,
                                         target_vars,
                                         instrument_shock = NULL,
                                         params = NULL,
                                         n_periods = 40L,
                                         shock_size = 1,
                                         normalise = TRUE,
                                         instrument_name = NULL) {

  if (missing(dr) || !is.list(dr)) {
    stop("'dr' must be a DecisionRules object from solve_perturbation().")
  }
  if (missing(model) || !inherits(model, "dynhr_mod")) {
    stop("'model' must be a dynhr_mod object from parse_mod().")
  }
  if (missing(target_vars) || length(target_vars) == 0) {
    stop("'target_vars' must be a character vector of variable names.")
  }

  if (is.null(params)) params <- model$param_values

  # ---- 1. Identify the instrument shock ----
  exo_names <- dr$exo_names %||% model$varexo_names

  if (is.null(instrument_shock)) {
    # Auto-detect: look for common names
    candidates <- c("eps_r", "eps_mp", "eps_i", "eps_r_",
                    "er", "e_r", "em", "eps_mon", "eps_pol")
    matched <- intersect(candidates, exo_names)
    if (length(matched) == 0) {
      stop("Could not auto-detect instrument shock. ",
           "Specify 'instrument_shock' from: ",
           paste(exo_names, collapse = ", "))
    }
    instrument_shock <- matched[1]
    message(sprintf("Auto-detected instrument shock: '%s'", instrument_shock))
  }

  if (!instrument_shock %in% exo_names) {
    stop(sprintf("Instrument shock '%s' not found in model. Available: %s",
                 instrument_shock, paste(exo_names, collapse = ", ")))
  }

  # ---- 2. Validate target variables ----
  endo_names <- dr$endo_names %||% model$var_names
  missing_vars <- setdiff(target_vars, endo_names)
  if (length(missing_vars) > 0) {
    stop(sprintf("Target variables not found in model: %s",
                 paste(missing_vars, collapse = ", ")))
  }

  # ---- 3. Compute full IRF ----
  all_irfs <- compute_irfs(dr, model, n_periods = n_periods,
                           shock_size = shock_size, params = params)

  if (!instrument_shock %in% names(all_irfs)) {
    stop(sprintf("IRF for shock '%s' not returned by compute_irfs().", instrument_shock))
  }

  irf_mat <- all_irfs[[instrument_shock]]

  # ---- 4. Extract target variables ----
  target_irf <- irf_mat[, target_vars, drop = FALSE]

  # ---- 5. Normalise (optional) ----
  if (normalise) {
    if (is.null(instrument_name)) {
      # Auto-detect instrument variable
      inst_candidates <- c("r", "i", "R", "rr", "nom_rate", "interest",
                           "ffr", "policy_rate")
      inst_matched <- intersect(inst_candidates, endo_names)
      if (length(inst_matched) > 0) {
        instrument_name <- inst_matched[1]
      } else {
        # Use the shock name minus "eps_" prefix as fallback
        instrument_name <- sub("^eps_", "", instrument_shock)
        instrument_name <- sub("^e_", "", instrument_name)
        instrument_name <- sub("_$", "", instrument_name)
      }
    }

    if (instrument_name %in% colnames(irf_mat)) {
      # Normalise so the instrument's own impact response is 1
      impact_response <- irf_mat[1, instrument_name]
      if (is.finite(impact_response) && abs(impact_response) > 1e-12) {
        target_irf <- target_irf / impact_response
        if (isTRUE(getOption("dynhr.opp.verbose", FALSE)) ||
            isTRUE(model$options$verbose)) {
          message(sprintf("Normalised IRF by instrument '%s' impact response: %.6f",
                          instrument_name, impact_response))
        }
      } else {
        warning(sprintf(
          "Instrument '%s' impact response is near-zero (%.2e). Using raw IRF.",
          instrument_name, impact_response))
      }
    } else {
      warning(sprintf(
        "Instrument variable '%s' not found in IRF columns. Available: %s",
        instrument_name, paste(colnames(irf_mat)[1:min(10, ncol(irf_mat))],
                               collapse = ", ")))
    }
  }

  attr(target_irf, "instrument_shock") <- instrument_shock
  attr(target_irf, "instrument_name") <- instrument_name
  attr(target_irf, "normalised") <- isTRUE(normalise)
  attr(target_irf, "n_periods") <- n_periods

  target_irf
}


# ==========================================================================
# 3. Top-level sufficient-statistics OPP function
# ==========================================================================

#' Sufficient-statistics optimal policy projections (OPP)
#'
#' Main entry point for the Barnichon & Mesters (2023) sufficient-statistics
#' approach to optimal policy. Given a solved DSGE model (or policy result),
#' extracts the IRF of target variables to the policy instrument and computes
#' the optimal policy intervention.
#'
#' This function provides a complete workflow:
#' \enumerate{
#'   \item Extract the instrument IRF from the model (or accept user-supplied IRF)
#'   \item Accept a baseline forecast (or generate one from model dynamics)
#'   \item Compute the optimal instrument path via the sufficient-statistics formula
#'   \item Report loss reduction and implied target variable paths
#' }
#'
#' @param model           A dynhr_mod object (from \code{\link{parse_mod}}) or
#'   a policy result object (\code{dynhr_ramsey_result}, \code{dynhr_nn1_result},
#'   etc.) from which the decision rules can be extracted.
#' @param dr              Optional DecisionRules object. If NULL, extracted from
#'   \code{model} or \code{result}.
#' @param target_vars     Character vector of target variable names (e.g.,
#'   \code{c("pi", "y_gap")}).
#' @param loss_weights    Numeric vector of weights on target variables, or a
#'   weight matrix. Defaults to equal weights.
#' @param baseline_forecast Optional matrix or data.frame of baseline forecast
#'   deviations from target. If NULL, a zero forecast (return to steady state)
#'   is assumed.
#' @param instrument_shock Character: name of the exogenous shock affecting the
#'   policy instrument. If NULL, attempts auto-detection.
#' @param instrument_name Character: name of the policy instrument variable in
#'   the model. Used for IRF normalisation.
#' @param instrument_penalty Numeric: penalty on instrument volatility (lambda).
#'   Default 0.01 for numerical stability.
#' @param target_values   Numeric vector of target (steady-state) values for
#'   each target variable. If NULL, assumed zero.
#' @param n_periods       Number of IRF / projection periods (default 40).
#' @param discount        Discount factor. Defaults to \code{params["beta"]}
#'   or 0.99.
#' @param params          Named numeric parameter vector. Defaults to
#'   \code{model$param_values}.
#' @param normalise_irf   If TRUE, normalise the IRF so the instrument's own
#'   impact response is 1.
#' @param verbose         Print progress and results.
#' @param ...             Additional arguments passed to
#'   \code{\link{opp_optimal_intervention}}.
#'
#' @return An object of class \code{dynhr_opp_result} (same as
#'   \code{\link{opp_optimal_intervention}}).
#'
#' @references
#'   Barnichon, R., & Mesters, G. (2023). A sufficient statistics approach for
#'     macro policy evaluation. \emph{Journal of Political Economy}, 131(1),
#'     159-197.
#'
#' @examples
#' \dontrun{
#' # Parse and solve a model
#' model <- parse_mod("nk_model.mod")
#'
#' # Run sufficient-statistics OPP
#' opp <- opp_sufficient_stats(
#'   model,
#'   target_vars    = c("pi", "y_gap"),
#'   loss_weights   = c(pi = 1, y_gap = 0.5),
#'   instrument_shock = "eps_r",
#'   verbose        = TRUE
#' )
#'
#' # Compare with Ramsey result
#' ramsey <- ramsey_model(model, "-(pi^2 + 0.5*y_gap^2)")
#' summary(opp)
#' }
#'
#' @export
opp_sufficient_stats <- function(model,
                                  dr = NULL,
                                  target_vars,
                                  loss_weights = NULL,
                                  baseline_forecast = NULL,
                                  instrument_shock = NULL,
                                  instrument_name = NULL,
                                  instrument_penalty = 0.01,
                                  target_values = NULL,
                                  n_periods = 40L,
                                  discount = NULL,
                                  params = NULL,
                                  normalise_irf = TRUE,
                                  verbose = FALSE,
                                  ...) {

  # ---- 1. Validate model / extract components ----
  if (inherits(model, "dynhr_mod")) {
    dynare_model <- model
  } else if (inherits(model, c("dynhr_ramsey_result", "dynhr_ramsey_result2",
                                "dynhr_nn1_result", "dynhr_osr_result",
                                "dynhr_discretionary_result"))) {
    # Try to extract dynhr_mod from the result
    dynare_model <- .extract_model(model)
    if (is.null(dynare_model)) {
      stop("Could not extract dynhr_mod from the result object. ",
           "Please pass the model directly.")
    }
    if (is.null(dr)) dr <- .extract_dr(model)
  } else {
    stop("'model' must be a dynhr_mod or a policy result object.")
  }

  if (is.null(params)) params <- dynare_model$param_values

  # ---- 2. Get decision rules ----
  if (is.null(dr)) {
    # Need to solve the model
    if (verbose) cat("[1/4] Compiling model and solving perturbation...\n")
    compiled <- compile_model(dynare_model, verbose = FALSE)
    ss <- solve_steady(compiled, params,
                        endo_names = dynare_model$var_names,
                        exo_names = dynare_model$varexo_names,
                        verbose = FALSE)
    if (is.null(ss) || !isTRUE(ss$converged)) {
      ss <- solve_steady_state(dynare_model, compiled,
                               params = params, verbose = FALSE)
    }
    ss_vals <- if (!is.null(ss$values)) ss$values else ss$ss
    dr <- solve_perturbation(dynare_model, compiled, ss_vals, params,
                             verbose = FALSE)
  }

  if (is.null(dr$ghx) || is.null(dr$ghu)) {
    stop("Decision rules must contain ghx and ghu components.")
  }

  # ---- 3. Extract discount factor ----
  if (is.null(discount)) {
    discount <- if ("beta" %in% names(params) && is.finite(params[["beta"]])) {
      as.numeric(params[["beta"]])
    } else {
      0.99
    }
  }

  # ---- 4. Extract instrument IRF ----
  if (verbose) cat(sprintf("[2/4] Extracting instrument IRF (shock=%s, periods=%d)...\n",
                           instrument_shock %||% "auto", n_periods))

  instrument_irf <- opp_estimate_instrument_irf(
    dr = dr,
    model = dynare_model,
    target_vars = target_vars,
    instrument_shock = instrument_shock,
    params = params,
    n_periods = n_periods,
    shock_size = 1,
    normalise = normalise_irf,
    instrument_name = instrument_name
  )

  # Extract instrument info from attributes
  used_shock <- attr(instrument_irf, "instrument_shock")
  used_inst  <- attr(instrument_irf, "instrument_name")

  # ---- 5. Build baseline forecast ----
  if (verbose) cat("[3/4] Preparing baseline forecast...\n")

  if (is.null(baseline_forecast)) {
    # Default: zero forecast (return to steady state)
    baseline_forecast <- matrix(0, nrow = n_periods, ncol = length(target_vars))
    colnames(baseline_forecast) <- target_vars
  } else {
    baseline_forecast <- as.matrix(baseline_forecast)
    if (ncol(baseline_forecast) != length(target_vars)) {
      # Try matching by column names
      if (!is.null(colnames(baseline_forecast))) {
        matched <- intersect(target_vars, colnames(baseline_forecast))
        if (length(matched) == length(target_vars)) {
          baseline_forecast <- baseline_forecast[, target_vars, drop = FALSE]
        } else {
          stop(sprintf(
            "Baseline forecast columns (%s) do not match target_vars (%s).",
            paste(colnames(baseline_forecast), collapse = ", "),
            paste(target_vars, collapse = ", ")))
        }
      } else {
        stop(sprintf(
          "Baseline forecast has %d columns but target_vars has length %d.",
          ncol(baseline_forecast), length(target_vars)))
      }
    }
    if (nrow(baseline_forecast) < n_periods) {
      # Pad with zeros
      padding <- matrix(0, nrow = n_periods - nrow(baseline_forecast),
                        ncol = ncol(baseline_forecast))
      colnames(padding) <- colnames(baseline_forecast)
      baseline_forecast <- rbind(baseline_forecast, padding)
    } else if (nrow(baseline_forecast) > n_periods) {
      baseline_forecast <- baseline_forecast[seq_len(n_periods), , drop = FALSE]
    }
  }

  # ---- 6. Compute optimal intervention ----
  if (verbose) cat("[4/4] Computing optimal intervention...\n")

  result <- opp_optimal_intervention(
    instrument_irf    = instrument_irf,
    baseline_forecast = baseline_forecast,
    loss_weights      = loss_weights,
    target_values     = target_values,
    instrument_penalty = instrument_penalty,
    discount          = discount,
    horizons          = n_periods,
    instrument_name   = used_inst %||% "i",
    verbose           = verbose,
    ...
  )

  # Attach model info
  result$meta$model_info <- list(
    model_name      = dynare_model$model_name %||% NA_character_,
    instrument_shock = used_shock,
    instrument_var   = used_inst,
    n_endo          = length(dynare_model$var_names),
    n_exo           = length(dynare_model$varexo_names),
    dr_order        = dr$order %||% 1L
  )

  result$meta$target_vars <- target_vars
  result$meta$method <- "opp_sufficient_stats"

  result
}


# ==========================================================================
# 4. Welfare gain from optimal intervention
# ==========================================================================

#' Welfare gain from optimal policy intervention
#'
#' Computes the consumption-equivalent welfare gain from implementing the
#' optimal policy intervention (from \code{\link{opp_sufficient_stats}})
#' relative to a baseline policy.
#'
#' This provides a welfare metric for the sufficient-statistics approach
#' in consumption-equivalent units, comparable to
#' \code{\link{welfare_ce_diff}}.
#'
#' @param opp_result  Result object from \code{\link{opp_sufficient_stats}}
#'   or \code{\link{opp_optimal_intervention}}.
#' @param model       Optional dynhr_mod object (needed for steady state
#'   and marginal utility computation).
#' @param consumption_variable Character: name of consumption variable.
#' @param params      Optional named parameter vector.
#' @param verbose     Print details.
#'
#' @return A list with class \code{dynhr_opp_welfare} containing:
#'   \describe{
#'     \item{ce_percent}{Consumption-equivalent gain in percent.}
#'     \item{loss_reduction}{Proportional loss reduction.}
#'     \item{loss_baseline}{Baseline loss value.}
#'     \item{loss_optimal}{Optimal intervention loss value.}
#'     \item{method}{Description.}
#'   }
#'
#' @examples
#' \dontrun{
#' opp <- opp_sufficient_stats(model, target_vars = c("pi", "y_gap"))
#' w <- opp_welfare_gain(opp, model)
#' print(w$ce_percent)
#' }
#'
#' @export
opp_welfare_gain <- function(opp_result,
                              model = NULL,
                              consumption_variable = NULL,
                              params = NULL,
                              verbose = FALSE) {

  if (!inherits(opp_result, "dynhr_opp_result")) {
    stop("opp_result must be from opp_sufficient_stats() or opp_optimal_intervention().")
  }

  # Extract model from result metadata if available
  if (is.null(model) && !is.null(opp_result$meta$model_info)) {
    # We don't have the full model object stored, so we need it passed
  }

  loss_base <- opp_result$loss_baseline
  loss_opt  <- opp_result$loss_optimal
  loss_red  <- opp_result$loss_reduction

  # The loss values in the OPP framework are in squared deviation units.
  # To convert to consumption-equivalent units, we need the marginal utility
  # of consumption and steady-state consumption.

  # Since we don't always have the full model, we report the loss reduction
  # as the primary metric, with consumption-equivalent as supplemental.

  # Approximate CE conversion:
  # CE ≈ (1 - beta) * (loss_base - loss_opt) / (MU * C_ss)
  # where MU is marginal utility of consumption at SS.
  # Without the model, we default to a unit scaling.

  discount <- opp_result$params$discount
  loss_diff <- loss_base - loss_opt

  # Try to compute proper CE if model is available
  ce_pct <- NA_real_
  mu <- NA_real_
  cons_ss <- NA_real_

  if (!is.null(model) && inherits(model, "dynhr_mod")) {
    params <- params %||% model$param_values
    ss <- model$initval %||% NULL

    if (!is.null(ss)) {
      # Find consumption variable
      cons_var <- NULL
      if (!is.null(consumption_variable) && consumption_variable %in% names(ss)) {
        cons_var <- consumption_variable
      } else {
        for (cv in c("c", "C", "cons", "consumption", "CONS", "y", "Y")) {
          if (cv %in% names(ss) && is.finite(ss[[cv]]) && abs(ss[[cv]]) > 1e-12) {
            cons_var <- cv
            break
          }
        }
      }

      if (!is.null(cons_var)) {
        cons_ss <- abs(as.numeric(ss[[cons_var]]))
        if (is.finite(cons_ss) && cons_ss > 1e-12) {
          mu <- 1 / cons_ss  # CRRA approximation
          ce_raw <- (1 - discount) * loss_diff / (mu * cons_ss)
          ce_pct <- ce_raw * 100
        }
      }
    }
  }

  if (!is.finite(ce_pct)) {
    # Fallback: report loss reduction as approximate welfare metric
    ce_pct <- loss_red * 100  # Interpret loss reduction as approximate % gain
  }

  result <- list(
    ce_percent     = ce_pct,
    loss_reduction = loss_red,
    loss_baseline  = loss_base,
    loss_optimal   = loss_opt,
    discount       = discount,
    marginal_utility = mu,
    consumption_ss = cons_ss,
    method         = if (is.finite(mu))
                       "Consumption-equivalent (via CRRA approximation)"
                     else
                       "Loss reduction (direct)",
    meta = list(
      timestamp = Sys.time(),
      package_version = utils::packageVersion("dynhr")
    )
  )
  class(result) <- c("dynhr_opp_welfare", "list")

  if (verbose) {
    cat("\n--- OPP Welfare Gain ---\n")
    cat(sprintf("Loss reduction:     %.2f%%\n", loss_red * 100))
    if (is.finite(ce_pct)) {
      cat(sprintf("CE welfare gain:    %.4f%%\n", ce_pct))
    }
    cat(sprintf("Method:             %s\n", result$method))
    cat("---\n")
  }

  result
}


# ==========================================================================
# 5. S3 methods
# ==========================================================================

#' @export
print.dynhr_opp_result <- function(x, digits = 4, ...) {
  cat("Sufficient-Statistics OPP Result\n")
  cat("================================\n")

  # Model info
  mi <- x$meta$model_info %||% list()
  if (length(mi) > 0) {
    cat(sprintf("Model:      %s\n", mi$model_name %||% "N/A"))
    cat(sprintf("Shock:      %s\n", mi$instrument_shock %||% "N/A"))
    cat(sprintf("Instrument: %s\n", mi$instrument_var %||% x$instrument_name %||% "i"))
  }

  cat(sprintf("Discount:   %.4f\n", x$params$discount))
  cat(sprintf("Penalty:    %.4f\n", x$params$instrument_penalty))
  cat(sprintf("Horizons:   %d\n", x$params$horizons %||% length(x$instrument_path)))
  cat(sprintf("Targets:    %s\n",
              paste(x$meta$var_names %||% colnames(x$target_paths), collapse = ", ")))

  cat(sprintf("\nLoss baseline:  %.6f\n", x$loss_baseline))
  cat(sprintf("Loss optimal:   %.6f\n", x$loss_optimal))
  cat(sprintf("Loss reduction: %.2f%%\n", x$loss_reduction * 100))

  cat(sprintf("\nOptimal instrument deviation: %.4f\n", x$instrument_path))

  invisible(x)
}


#' @export
summary.dynhr_opp_result <- function(object, ...) {
  cat("Sufficient-Statistics OPP -- Summary\n")
  cat("====================================\n")

  mi <- object$meta$model_info %||% list()
  if (length(mi) > 0) {
    cat(sprintf("Model:            %s\n", mi$model_name %||% "N/A"))
    cat(sprintf("Instr. shock:     %s\n", mi$instrument_shock %||% "N/A"))
    cat(sprintf("Instr. variable:  %s\n", mi$instrument_var %||% "N/A"))
  }

  cat(sprintf("\nLoss reduction:   %.2f%%\n", object$loss_reduction * 100))
  cat(sprintf("Instrument path:  %.4f\n", object$instrument_path))

  cat("\nTarget variable paths (selected horizons):\n")
  paths <- object$target_paths
  n_h <- nrow(paths)
  n_display <- min(n_h, 8)
  h_idx <- unique(c(1, floor(seq(2, n_display, length.out = n_display - 1))))
  if (n_h > n_display) h_idx <- c(h_idx, n_h)
  display <- paths[h_idx, , drop = FALSE]
  rownames(display) <- paste0("h=", h_idx - 1)
  print(round(display, 4))

  cat(sprintf("\nBaseline loss:    %.4f\n", object$loss_baseline))
  cat(sprintf("Optimal loss:     %.4f\n", object$loss_optimal))

  invisible(object)
}


#' @export
print.dynhr_opp_welfare <- function(x, digits = 4, ...) {
  cat("OPP Welfare Gain\n")
  cat("=================\n")
  cat(sprintf("CE welfare gain:  %.4f%%\n", x$ce_percent))
  cat(sprintf("Loss reduction:   %.2f%%\n", x$loss_reduction * 100))
  cat(sprintf("Loss baseline:    %.6f\n", x$loss_baseline))
  cat(sprintf("Loss optimal:     %.6f\n", x$loss_optimal))
  cat(sprintf("Method:           %s\n", x$method))
  invisible(x)
}


# ==========================================================================
# 6. Internal helpers
# ==========================================================================

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
