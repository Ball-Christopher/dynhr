## R/irf-matching.R
## --------------------------------------------------------------------------
## match_irfs() -- turnkey impulse-response-matching estimation driver.
##
## Minimum-distance estimation of a subset of structural parameters by
## matching model-implied impulse responses to a set of target IRFs (e.g.
## empirical IRFs from a VAR, or "true" model IRFs in a Monte-Carlo recovery
## exercise).  Wraps the existing solve + compute_irfs() pipeline inside an
## optim() call so users no longer have to assemble the loop by hand.
## --------------------------------------------------------------------------


## --------------------------------------------------------------------------
## Internal: extract the target IRF for `shock` as a (horizon x variable)
## matrix, regardless of whether the caller passed an IRFCollection, a plain
## named list of matrices, or a single matrix.
## --------------------------------------------------------------------------
.match_irfs_extract_target <- function(target_irfs, shock) {
  ## IRFCollection or named list keyed by shock name.
  if (is.list(target_irfs) && !is.data.frame(target_irfs)) {
    if (is.null(shock))
      stop("`shock` must be supplied when `target_irfs` is keyed by shock name.")
    if (!shock %in% names(target_irfs))
      stop(sprintf("shock '%s' not found in target_irfs (have: %s).",
                   shock, paste(names(target_irfs), collapse = ", ")))
    tgt <- target_irfs[[shock]]
  } else {
    tgt <- target_irfs
  }
  tgt <- as.matrix(tgt)
  if (is.null(colnames(tgt)))
    stop("`target_irfs` must have column names matching model variables.")
  tgt
}


#' Estimate structural parameters by impulse-response matching
#'
#' Minimum-distance estimation of a subset of a model's structural parameters
#' by matching model-implied impulse responses to a set of target IRFs.  This
#' is a turnkey driver around the existing solve + \code{\link{compute_irfs}}
#' pipeline: it builds the objective, calls \code{\link[stats]{optim}}, and
#' returns the estimated parameters together with the fitted IRFs.
#'
#' For a single shock the objective is the weighted sum of squared deviations
#' \deqn{Q(\theta) = \sum_{h,j} w_{h,j}\,\bigl(g_{h,j}(\theta) - \hat g_{h,j}\bigr)^2,}
#' where \eqn{g_{h,j}(\theta)} is the model-implied response of variable
#' \eqn{j} at horizon \eqn{h} (re-solving the steady state and the
#' perturbation at the candidate \eqn{\theta}), \eqn{\hat g_{h,j}} is the
#' target response, and \eqn{w_{h,j}} is the (optional) weight.  With the
#' default identity weight this is the unweighted distance \eqn{\lVert g(\theta) - \hat g\rVert^2}.
#'
#' Candidate \eqn{\theta} values that fail to solve (steady state does not
#' converge, Blanchard--Kahn violated, non-finite IRFs) are penalised with a
#' large finite objective so the optimiser steers away from them rather than
#' crashing.
#'
#' @param model A \code{dynhr_mod} (e.g. from \code{\link{parse_mod}}), or a
#'   path to a \code{.mod} file.
#' @param params_init Named numeric vector of starting values.  Must supply a
#'   value for every name in \code{free_params}.  Names not in
#'   \code{free_params} are treated as fixed calibration overrides for the
#'   duration of the fit (the model's own \code{param_values} provide any
#'   remaining defaults).
#' @param target_irfs The IRFs to match.  Either (a) a single
#'   \code{horizon x variable} numeric matrix (column names = model variable
#'   names) for the chosen \code{shock}, or (b) an \code{IRFCollection} /
#'   named list keyed by shock name (the \code{shock} entry is extracted).
#'   Only the variables present as columns are matched; extra model variables
#'   are ignored.
#' @param free_params Character vector of parameter names to estimate (a
#'   subset of \code{model$param_names}).
#' @param shock Name of the exogenous shock whose IRF is matched.  Required
#'   when \code{target_irfs} is keyed by shock; for a bare matrix it selects
#'   which model shock generated the target.
#' @param horizons Integer vector of horizons (1-based periods) to match.
#'   \code{NULL} (default) uses all horizons (rows) present in
#'   \code{target_irfs}.
#' @param weight Weighting of the squared deviations.  \code{NULL} (default)
#'   is the identity (unweighted).  May be a single scalar, a per-horizon
#'   numeric vector (length = number of matched horizons), or a full
#'   \code{horizon x variable} matrix conformable with the matched target.
#' @param compiled Optional pre-compiled model (from
#'   \code{\link{compile_model}}) to reuse across the optimisation; compiled
#'   once internally when \code{NULL}.
#' @param order Perturbation order used to generate model IRFs.  Default
#'   \code{1L}; \code{2L} uses the pruned second-order IRFs
#'   (\code{\link{compute_irfs_order2}}).
#' @param method Optimiser.  Either \code{"csminwel"} (Christopher Sims's
#'   DSGE-estimation BFGS minimiser; see \code{\link{csminwel}}) or any method
#'   passed to \code{\link[stats]{optim}}: \code{"Nelder-Mead"} (default),
#'   \code{"BFGS"}, \code{"L-BFGS-B"}, etc.  Nelder-Mead is the robust default;
#'   gradient methods (\code{"BFGS"}) can overshoot into the non-stationary
#'   region on poorly-scaled single-parameter problems, so prefer
#'   \code{"L-BFGS-B"} with \code{lower}/\code{upper} bounds, or
#'   \code{"csminwel"}, if you want a quasi-Newton method.  With
#'   \code{"csminwel"}, \code{control$reltol} and \code{control$maxit} map to
#'   csminwel's \code{crit}/\code{nit} and \code{lower}/\code{upper} are
#'   ignored (csminwel is unconstrained).
#' @param lower,upper Optional parameter bounds (named or positional, length
#'   \code{free_params}) forwarded to \code{optim} for box-constrained methods
#'   such as \code{"L-BFGS-B"}.
#' @param penalty Finite objective value returned for non-solvable candidate
#'   parameters.  Default \code{1e10}.
#' @param control List passed to \code{optim}'s \code{control} argument.
#' @param verbose Logical; if \code{TRUE} print the candidate objective at
#'   each evaluation.  Default \code{FALSE}.
#' @param ... Additional arguments forwarded to \code{\link[stats]{optim}}.
#'
#' @return A list of class \code{irf_match_result} with elements:
#'   \describe{
#'     \item{par}{Named numeric vector of estimated free parameters.}
#'     \item{params}{Full named parameter vector at the optimum (free +
#'       fixed).}
#'     \item{objective}{Objective (weighted squared distance) at the optimum.}
#'     \item{convergence}{\code{optim} convergence code (\code{0} = success).}
#'     \item{counts, message}{\code{optim}'s function/gradient counts and
#'       message.}
#'     \item{fitted_irf}{The model-implied \code{horizon x variable} IRF at
#'       the optimum (matched variables/horizons).}
#'     \item{target_irf}{The target IRF that was matched.}
#'     \item{optim}{The raw \code{optim} return value.}
#'   }
#'
#' @section Empirical targets and optimiser:
#' Use \code{\link{estimate_var}} + \code{\link{var_irf}} /
#' \code{\link{var_irf_bootstrap}} to construct empirical IRF targets (and
#' bootstrap bands) from data; the \code{horizon x variable} matrices they
#' return are directly consumable as \code{target_irfs}.  Pass
#' \code{method = "csminwel"} to estimate with Sims's \code{\link{csminwel}}
#' optimiser instead of \code{\link[stats]{optim}}.  Supply \code{weight}
#' directly if you have a precision/optimal-weighting matrix.
#'
#' @examples
#' \dontrun{
#' rbc <- parse_mod(system.file("extdata/models/rbc.mod", package = "dynhr"))
#' ## Target IRFs at the true calibration:
#' truth  <- stoch_simul(rbc, verbose = FALSE)
#' target <- truth$irfs[["eps_a"]]
#' ## Recover rho_a from a perturbed start:
#' fit <- match_irfs(rbc,
#'                   params_init = c(rho_a = 0.80),
#'                   target_irfs = target,
#'                   free_params = "rho_a",
#'                   shock       = "eps_a")
#' fit$par           # ~ 0.95
#' }
#'
#' @seealso \code{\link{compute_irfs}}, \code{\link{stoch_simul}}
#' @export
match_irfs <- function(model,
                       params_init,
                       target_irfs,
                       free_params,
                       shock,
                       horizons = NULL,
                       weight   = NULL,
                       compiled = NULL,
                       order    = 1L,
                       method   = "Nelder-Mead",
                       lower    = -Inf,
                       upper    =  Inf,
                       penalty  = 1e10,
                       control  = list(),
                       verbose  = FALSE,
                       ...) {

  ## ---- Parse / validate the model -------------------------------------
  if (is.character(model) && length(model) == 1L && file.exists(model)) {
    model <- parse_mod(model, verbose = FALSE)
  }
  if (!inherits(model, "dynhr_mod"))
    stop("`model` must be a dynhr_mod or a path to a .mod file.")

  order <- as.integer(order)
  if (!order %in% c(1L, 2L))
    stop("match_irfs() supports order = 1 or 2.")

  ## ---- Validate free_params -------------------------------------------
  if (!is.character(free_params) || length(free_params) == 0L)
    stop("`free_params` must be a non-empty character vector of parameter names.")
  unknown <- setdiff(free_params, model$param_names)
  if (length(unknown))
    stop(sprintf("free_params not in model$param_names: %s",
                 paste(unknown, collapse = ", ")))
  if (anyDuplicated(free_params))
    stop("`free_params` contains duplicate names.")

  ## ---- Starting values ------------------------------------------------
  if (is.null(names(params_init)))
    stop("`params_init` must be a *named* numeric vector.")
  missing_init <- setdiff(free_params, names(params_init))
  if (length(missing_init))
    stop(sprintf("params_init lacks starting values for: %s",
                 paste(missing_init, collapse = ", ")))
  theta0 <- as.numeric(params_init[free_params])
  names(theta0) <- free_params

  ## Base parameter vector: model defaults, then any fixed overrides from
  ## params_init (names that are NOT free are treated as calibration).
  base_params <- model$param_values
  fixed_over  <- setdiff(names(params_init), free_params)
  for (nm in fixed_over)
    if (nm %in% names(base_params)) base_params[[nm]] <- params_init[[nm]]

  ## ---- Target IRF -----------------------------------------------------
  if (missing(shock)) shock <- NULL
  tgt_full <- .match_irfs_extract_target(target_irfs, shock)

  ## Match only model variables that are present as target columns.
  match_vars <- intersect(colnames(tgt_full), model$var_names)
  if (length(match_vars) == 0L)
    stop("No target_irfs column names match model variables.")

  ## Horizons.
  n_h_target <- nrow(tgt_full)
  if (is.null(horizons)) {
    horizons <- seq_len(n_h_target)
  } else {
    horizons <- as.integer(horizons)
    if (any(horizons < 1L) || any(horizons > n_h_target))
      stop(sprintf("`horizons` out of range 1..%d.", n_h_target))
  }
  n_periods <- max(horizons)

  tgt <- tgt_full[horizons, match_vars, drop = FALSE]

  ## ---- Weight matrix --------------------------------------------------
  W <- .match_irfs_build_weight(weight, dim(tgt))

  ## ---- Compile once and reuse -----------------------------------------
  max_order <- max(1L, order)
  if (is.null(compiled))
    compiled <- compile_model(model, verbose = FALSE, max_order = max_order)

  ## ---- Build the model IRF for a candidate theta ----------------------
  .model_irf <- function(theta) {
    params <- base_params
    params[free_params] <- theta
    ## Re-solve steady state (handles SS-block params) then perturbation.
    ss <- tryCatch(
      solve_steady_state(model, compiled, params, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(ss) || !isTRUE(ss$converged)) return(NULL)

    dr <- tryCatch(
      solve_perturbation(model, compiled, ss, params,
                         order = order, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(dr) || !isTRUE(dr$bk_satisfied)) return(NULL)

    ## Stationarity guard: Blanchard-Kahn only flags explosive *forward-
    ## looking* roots; an explosive *exogenous* state (e.g. an AR(1) shock
    ## with rho >= 1) passes BK but yields a diverging, non-stationary IRF
    ## with a finite-but-meaningless SSE that can trap a gradient optimiser.
    ## Reject candidates whose state-transition block is not stable.
    Tmat <- dr$ghx[dr$state_idx, , drop = FALSE]
    if (length(Tmat) > 0L) {
      sr <- tryCatch(max(Mod(eigen(Tmat, only.values = TRUE)$values)),
                     error = function(e) Inf)
      if (!is.finite(sr) || sr >= 1 - 1e-8) return(NULL)
    }

    irfs <- tryCatch({
      if (order == 1L) {
        compute_irfs(dr, model, n_periods = n_periods, params = params)
      } else {
        compute_irfs_order2(dr, model, n_periods = n_periods, params = params)
      }
    }, error = function(e) NULL)
    if (is.null(irfs)) return(NULL)

    sh <- shock %||% names(irfs)[1L]
    if (!sh %in% names(irfs)) return(NULL)
    g <- irfs[[sh]]
    if (!all(match_vars %in% colnames(g))) return(NULL)
    g[horizons, match_vars, drop = FALSE]
  }

  ## ---- Objective ------------------------------------------------------
  obj <- function(theta) {
    g <- .model_irf(theta)
    if (is.null(g) || anyNA(g) || any(!is.finite(g))) {
      if (verbose) cat(sprintf("  obj = %.6g  [penalty]\n", penalty))
      return(penalty)
    }
    d <- g - tgt
    val <- sum(W * d * d)
    if (!is.finite(val)) val <- penalty
    if (verbose) cat(sprintf("  obj = %.6g\n", val))
    val
  }

  ## ---- Optimise -------------------------------------------------------
  ## `method = "csminwel"` routes to Sims's csminwel (see csminwel()); any
  ## other value is passed straight through to stats::optim() as before.
  if (identical(method, "csminwel")) {
    cs <- csminwel(obj, theta0,
                   crit    = control$reltol %||% 1e-7,
                   nit     = control$maxit  %||% 1000L,
                   verbose = verbose, ...)
    ## Adapt csminwel's return shape to the optim-style fields used below.
    opt <- list(par = cs$xh, value = cs$fh, convergence = cs$convergence,
                counts = c("function" = cs$fcount, gradient = NA_integer_),
                message = sprintf("csminwel retcode %d", cs$retcode),
                csminwel = cs)
  } else {
    bounded <- method %in% c("L-BFGS-B", "Brent")
    optim_args <- list(par = theta0, fn = obj, method = method,
                       control = control, ...)
    if (bounded) {
      optim_args$lower <- lower
      optim_args$upper <- upper
    }
    ## A single free parameter with a derivative-free method needs Brent or a
    ## 1-D-aware setup; BFGS/Nelder-Mead handle it, but optim() warns for
    ## Nelder-Mead with one parameter. Leave that to the caller's method choice.
    opt <- do.call(stats::optim, optim_args)
  }

  ## ---- Assemble result ------------------------------------------------
  par_hat <- opt$par
  names(par_hat) <- free_params
  params_hat <- base_params
  params_hat[free_params] <- par_hat
  fitted <- .model_irf(par_hat)

  result <- list(
    par         = par_hat,
    params      = params_hat,
    objective   = opt$value,
    convergence = opt$convergence,
    counts      = opt$counts,
    message     = opt$message,
    fitted_irf  = fitted,
    target_irf  = tgt,
    shock       = shock,
    free_params = free_params,
    horizons    = horizons,
    optim       = opt
  )
  class(result) <- "irf_match_result"
  result
}


## --------------------------------------------------------------------------
## Internal: coerce the `weight` argument into a matrix conformable with the
## matched target (dim = c(n_horizon, n_var)).
## --------------------------------------------------------------------------
.match_irfs_build_weight <- function(weight, dims) {
  n_h <- dims[1L]; n_v <- dims[2L]
  if (is.null(weight))
    return(matrix(1, n_h, n_v))
  if (is.matrix(weight)) {
    if (!all(dim(weight) == dims))
      stop(sprintf("`weight` matrix must be %d x %d (horizons x variables).",
                   n_h, n_v))
    if (any(weight < 0))
      stop("`weight` must be non-negative.")
    return(weight)
  }
  if (length(weight) == 1L) {
    if (weight < 0) stop("`weight` must be non-negative.")
    return(matrix(weight, n_h, n_v))
  }
  if (length(weight) == n_h) {
    if (any(weight < 0)) stop("`weight` must be non-negative.")
    return(matrix(weight, n_h, n_v))  # per-horizon, recycled across variables
  }
  stop(sprintf(paste0("`weight` must be NULL, a scalar, a length-%d per-horizon ",
                      "vector, or a %d x %d matrix."), n_h, n_h, n_v))
}


#' @export
print.irf_match_result <- function(x, ...) {
  cat("<irf_match_result>\n")
  cat(sprintf("  shock         : %s\n", x$shock %||% "(default)"))
  cat(sprintf("  free params   : %s\n", paste(x$free_params, collapse = ", ")))
  cat(sprintf("  objective     : %.6g\n", x$objective))
  cat(sprintf("  convergence   : %d%s\n", x$convergence,
              if (x$convergence == 0L) " (converged)" else ""))
  cat("  estimates:\n")
  for (nm in names(x$par))
    cat(sprintf("    %-12s = %.6g\n", nm, x$par[[nm]]))
  invisible(x)
}
