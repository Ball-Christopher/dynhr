## R/diag-pre-d39-stability-map.R
## --------------------------------------------------------------------------
## D39: Blanchard-Kahn stability mapping.
##
##   * prior mode (default): Ratto (2008) Monte-Carlo filtering -- classify
##     prior draws and rank parameters by the KS distance between the
##     determinate and non-determinate draws.
##   * grid mode (`grid = list(px = ..., py = ...)`): classify every point of a
##     2-D parameter grid (all other parameters held at `params`) and draw a
##     tile map of the stability regions with the calibration point marked.
##
## Every point is RE-SOLVED (steady state + system matrices + QZ) at its own
## parameter vector and assigned exactly one of the classes in
## .d39_classes(); solver failure, missing steady state, indeterminacy and
## explosiveness are never merged.
##
## Standalone by design: run_all_diagnostics() does not call D39 (it needs
## only a model + parameter vector; see ?run_diagnostics "Coverage").
## --------------------------------------------------------------------------

#' D39. Blanchard-Kahn stability map (Ratto 2008 / Monte-Carlo filtering)
#'
#' Classifies parameter vectors by the Blanchard-Kahn (BK) properties of the
#' first-order solution, re-solving the steady state and the linearised system
#' at every point. Each point gets exactly one status:
#' \describe{
#'   \item{\code{determinate}}{unique stable solution (BK satisfied).}
#'   \item{\code{indeterminate}}{fewer explosive roots than forward-looking
#'     variables (multiple stable solutions).}
#'   \item{\code{explosive}}{more explosive roots than forward-looking
#'     variables (no stable solution).}
#'   \item{\code{unit_root}}{a finite generalised eigenvalue lies on the unit
#'     circle, \eqn{||\lambda| - 1| \le} \code{unit_tol}: the point is on
#'     (numerically) the determinacy boundary, where the root count is
#'     ambiguous. Models with a structural unit root are \code{unit_root}
#'     everywhere.}
#'   \item{\code{invalid_rule}}{the root count matches, but the solver's
#'     decision rule is not a valid unique solution: either its realised state
#'     transition is explosive (the solver's post-solve guard rejected it), or
#'     the QZ block \eqn{Z_{11}} is singular -- the BK \emph{rank} condition
#'     fails (e.g. an explosive root belongs to a predetermined variable while
#'     a jump variable has a stable root), which since 0.9.4
#'     \code{solve_perturbation()} reports as \code{bk_satisfied = FALSE}
#'     with \code{bk_rank_deficient = TRUE} (it used to report
#'     \code{bk_satisfied = TRUE}).}
#'   \item{\code{no_steady_state}}{the steady-state solver did not converge.}
#'   \item{\code{solver_error}}{the steady-state or perturbation solver
#'     raised an error (message kept in \code{$errors}).}
#' }
#' Forward-looking variables are counted as in \code{solve_perturbation()}:
#' forward-only plus mixed variables, less mixed variables declared
#' predetermined.
#'
#' \strong{Prior mode} (\code{grid = NULL}): samples \code{n_draws} vectors from
#' the prior (the same sampler SMC uses, so draws match the density
#' \code{log_prior} scores) and ranks parameters by the two-sample
#' Kolmogorov-Smirnov statistic between determinate and all other draws
#' (\code{$drivers}), and between determinate draws and each failing class
#' separately (\code{$drivers_by_class}). A large KS statistic means the
#' parameter's value strongly separates the regions (Ratto 2008).
#'
#' \strong{Grid mode} (\code{grid} a named list of two numeric vectors):
#' classifies the Cartesian grid, the FIRST element on the x axis and the
#' SECOND on the y axis, with every other parameter at \code{params}. Returns
#' the grid table and a tile map (\code{$plots$map}) with the calibration point
#' (\code{params} at the two grid parameters) marked.
#'
#' @param model      Parsed model object (from \code{parse_mod()}).
#' @param compiled   Compiled model (from \code{compile_model()}). If
#'                   \code{NULL}, compiled internally.
#' @param prior_spec Prior-mode only: \code{data.frame} as returned by
#'                   \code{extract_prior_spec()}. If \code{NULL}, extracted
#'                   from \code{model}.
#' @param n_draws    Prior-mode only: number of prior draws (default
#'                   \code{2000L}).
#' @param grid       \code{NULL} (prior mode) or a named list of two finite
#'                   numeric vectors, \code{list(x_param = ..., y_param = ...)}.
#' @param params     Named base parameter vector (default
#'                   \code{model$param_values}); parameters not being varied
#'                   are held here, and in grid mode it is the calibration
#'                   point marked on the map.
#' @param unit_tol   Unit-circle tolerance for the \code{unit_root} class
#'                   (default \code{1e-6}, the solver's QZ threshold).
#' @param seed       Optional integer random seed. The result is reproducible
#'   for a fixed seed; the caller's global RNG stream is restored on exit.
#' @param verbose    Logical; report progress (default \code{FALSE}).
#' @param ...        Passed to \code{compile_model()} if compiling internally.
#'
#' @return A \code{dynhr_stability_map} list with elements:
#' \describe{
#'   \item{\code{mode}}{\code{"prior"} or \code{"grid"}.}
#'   \item{\code{status}}{factor (levels as above), one per draw / grid point.}
#'   \item{\code{class_counts}}{named integer counts of every class.}
#'   \item{\code{feasible_fraction}}{fraction of points that are
#'     \code{determinate}.}
#'   \item{\code{n_feasible}, \code{n_total}}{determinate / total counts.}
#'   \item{\code{errors}}{character vector of solver-error messages
#'     (one per \code{solver_error} point).}
#'   \item{\code{drivers}}{prior mode: data.frame sorted by \code{ks_stat}
#'     (descending): \code{param}, \code{ks_stat}, \code{p_value},
#'     \code{feasible_lo}, \code{feasible_hi}; \code{NULL} in grid mode.}
#'   \item{\code{drivers_by_class}}{prior mode: named list of such tables,
#'     determinate vs each failing class with >= 2 draws.}
#'   \item{\code{draws}}{prior mode: \code{n_draws x n_par} matrix of draws.}
#'   \item{\code{grid}}{grid mode: data.frame with \code{x}, \code{y},
#'     \code{status}, \code{n_unstable}, \code{n_forward}.}
#'   \item{\code{axes}}{grid mode: \code{c(x = , y = )} parameter names.}
#'   \item{\code{calibration}}{grid mode: the marked point (named, x then y).}
#'   \item{\code{plots}}{grid mode: \code{list(map = <ggplot>)}; else empty.}
#' }
#'
#' @references
#'   Ratto, M. (2008). Analysing DSGE models with global sensitivity analysis.
#'   \emph{Computational Economics}, 31(2), 115--139.
#'
#'   Blanchard, O. J. and Kahn, C. M. (1980). The solution of linear difference
#'   models under rational expectations. \emph{Econometrica}, 48(5),
#'   1305--1311.
#'
#' @export
diag_stability_map <- function(model,
                               compiled   = NULL,
                               prior_spec = NULL,
                               n_draws    = 2000L,
                               grid       = NULL,
                               params     = NULL,
                               unit_tol   = 1e-6,
                               seed       = NULL,
                               verbose    = FALSE,
                               ...) {

  .local_seed(seed)  # seeded, but the caller's RNG stream is restored on exit

  base_params <- if (is.null(params)) model$param_values else params
  if (is.null(base_params)) base_params <- numeric(0)
  if (length(base_params) > 0L && is.null(names(base_params)))
    .dynhr_abort("diag_stability_map: `params` must be a named numeric vector.")

  if (!is.null(grid)) .d39_check_grid(grid)

  if (is.null(compiled)) {
    if (verbose) .dynhr_inform("diag_stability_map: compiling model ...")
    compiled <- compile_model(model, verbose = FALSE, max_order = 1L, ...)
  }
  n_fwd <- .d39_n_forward(model)

  if (!is.null(grid)) {
    return(.d39_run_grid(model, compiled, grid, base_params, n_fwd,
                         unit_tol, verbose))
  }

  ## ---------------------------- prior mode --------------------------------
  if (is.null(prior_spec)) {
    prior_spec <- extract_prior_spec(model, verbose = FALSE)
  }
  if (is.null(prior_spec) || nrow(prior_spec) == 0L) {
    .dynhr_abort("diag_stability_map: no prior_spec available. ",
                 "Supply prior_spec= or ensure model has an estimated_params block.")
  }
  param_names <- prior_spec$name
  n_par       <- length(param_names)
  n_draws     <- as.integer(n_draws)
  if (length(n_draws) != 1L || is.na(n_draws) || n_draws < 1L)
    .dynhr_abort("diag_stability_map: `n_draws` must be a positive integer.")

  ## Shared SMC prior sampler: uniform support via .uniform_ab() (Dynare mean/sd or p3/p4 bounds);
  ## draws inv_gamma as IG1, i.e. the distribution log_prior() scores.
  prior_draw_fn <- .smc_make_prior_sampler(prior_spec)

  draw_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_par,
                     dimnames = list(NULL, param_names))
  status   <- character(n_draws)
  errors   <- character(0)

  if (verbose) .dynhr_inform(sprintf("diag_stability_map: classifying %d prior draws ...", n_draws))

  for (i in seq_len(n_draws)) {
    theta_i <- prior_draw_fn()[param_names]
    draw_mat[i, ] <- theta_i
    params_i <- base_params
    params_i[param_names] <- theta_i
    cl <- .d39_classify_point(model, compiled, params_i, n_fwd, unit_tol)
    status[i] <- cl$status
    if (!is.null(cl$error)) errors <- c(errors, cl$error)
  }

  status <- factor(status, levels = .d39_classes())
  counts <- table(status)
  n_feasible <- as.integer(counts[["determinate"]])

  if (verbose) {
    .dynhr_inform(sprintf("diag_stability_map: determinate fraction = %.3f  (%d / %d)",
                          n_feasible / n_draws, n_feasible, n_draws))
  }

  feas_idx <- which(status == "determinate")
  drivers  <- .d39_ks_ranking(draw_mat, feas_idx,
                              which(status != "determinate"), param_names)
  by_class <- list()
  for (cls in setdiff(.d39_classes(), "determinate")) {
    idx <- which(status == cls)
    if (length(idx) >= 2L)
      by_class[[cls]] <- .d39_ks_ranking(draw_mat, feas_idx, idx, param_names)
  }

  structure(
    list(
      mode              = "prior",
      status            = status,
      class_counts      = stats::setNames(as.integer(counts), names(counts)),
      feasible_fraction = n_feasible / n_draws,
      n_feasible        = n_feasible,
      n_total           = n_draws,
      errors            = errors,
      drivers           = drivers,
      drivers_by_class  = by_class,
      draws             = draw_mat,
      plots             = list()
    ),
    class = c("dynhr_stability_map", "list")
  )
}


## Ordered set of point classes (factor levels, legend order).
.d39_classes <- function() {
  c("determinate", "indeterminate", "explosive", "unit_root",
    "invalid_rule", "no_steady_state", "solver_error")
}

## Tol light fills per class (fixed mapping so a class keeps its colour
## whichever subset of classes a map shows).
.d39_class_fills <- function() {
  c(determinate     = unname(tol_light["light_blue"]),
    indeterminate   = unname(tol_light["light_yellow"]),
    explosive       = unname(tol_light["orange"]),
    unit_root       = unname(tol_light["mint"]),
    invalid_rule    = unname(tol_light["pink"]),
    no_steady_state = unname(tol_light["pear"]),
    solver_error    = unname(tol_light["pale_grey"]))
}

## Number of forward-looking (jump) variables, counted as solve_perturbation()
## does for its BK comparison (n_plus_bk_model).
.d39_n_forward <- function(model) {
  n_fwd <- model$n_forward %||% 0L
  n_mix <- model$n_mixed %||% 0L
  n_pred_mix <- sum(model$variable_classification$mixed %in%
                    (model$predetermined_vars %||% character(0)))
  as.integer(n_fwd + n_mix - n_pred_mix)
}


## Classify one parameter vector. Re-solves the steady state and the
## first-order system at `params` (nothing is reused across points).
## Returns list(status, n_unstable, error).
.d39_classify_point <- function(model, compiled, params, n_fwd, unit_tol) {
  ## One sweep point may make a solver throw (singular Jacobian, non-finite
  ## residual at an extreme draw). That is an OUTCOME to record as
  ## `solver_error`, not a reason to abort the whole map, so the error is
  ## caught here and kept (message returned) -- never folded into another class.
  err <- NULL
  out <- tryCatch({
    ss <- solve_steady(compiled, params, verbose = FALSE)
    if (is.null(ss) || !isTRUE(ss$converged) || !all(is.finite(ss$values))) {
      list(status = "no_steady_state", n_unstable = NA_integer_)
    } else {
      ## BK violations are reported by a warning AND in the return value; the
      ## return value is what we classify on, so muffle the per-point
      ## warnings. Since 0.9.4 a singular QZ block Z11 (the BK RANK condition,
      ## i.e. the unstable roots are not attached to the jump variables) sets
      ## `bk_satisfied = FALSE` and flags `bk_rank_deficient` in the returned
      ## rule, so this no longer has to match warning TEXT.
      dr <- suppressWarnings(
        solve_perturbation(model, compiled, ss$values, params,
                           order = 1L, verbose = FALSE))
      list(status = .d39_status_from_dr(dr, n_fwd, unit_tol,
                                        isTRUE(dr$bk_rank_deficient)),
           n_unstable = as.integer(dr$n_unstable %||% NA_integer_))
    }
  }, error = function(e) {
    list(status = "solver_error", n_unstable = NA_integer_,
         error = conditionMessage(e))
  })
  out
}

.d39_status_from_dr <- function(dr, n_fwd, unit_tol, rank_fail = FALSE) {
  ev  <- dr$eigenvalues
  mod <- Mod(ev)
  mod <- mod[is.finite(mod)]
  if (length(mod) > 0L && any(abs(mod - 1) <= unit_tol)) return("unit_root")
  if (isTRUE(dr$bk_satisfied)) return(if (rank_fail) "invalid_rule" else "determinate")
  nu <- dr$n_unstable
  if (is.null(nu) || !is.finite(nu)) return("solver_error")
  if (nu < n_fwd) return("indeterminate")
  if (nu > n_fwd) return("explosive")
  "invalid_rule"
}


.d39_check_grid <- function(grid) {
  ok <- is.list(grid) && length(grid) == 2L && !is.null(names(grid)) &&
    all(nzchar(names(grid))) && !anyDuplicated(names(grid))
  if (!ok)
    .dynhr_abort("diag_stability_map: `grid` must be a named list of two ",
                 "numeric vectors, e.g. list(phi_pi = ..., phi_y = ...).")
  for (nm in names(grid)) {
    v <- grid[[nm]]
    if (!is.numeric(v) || length(v) < 1L || !all(is.finite(v)))
      .dynhr_abort(sprintf("diag_stability_map: grid$%s must be a finite numeric vector.", nm))
  }
  invisible(TRUE)
}


.d39_run_grid <- function(model, compiled, grid, base_params, n_fwd,
                          unit_tol, verbose) {
  axes <- names(grid)
  missing_p <- setdiff(axes, names(base_params))
  if (length(missing_p) > 0L)
    .dynhr_abort("diag_stability_map: grid parameter(s) not in the model ",
                 "parameters: ", paste(missing_p, collapse = ", "))
  xs <- sort(unique(as.numeric(grid[[1L]])))
  ys <- sort(unique(as.numeric(grid[[2L]])))
  tab <- expand.grid(x = xs, y = ys, KEEP.OUT.ATTRS = FALSE)
  n <- nrow(tab)
  status <- character(n)
  n_unst <- integer(n)
  errors <- character(0)
  if (verbose)
    .dynhr_inform(sprintf("diag_stability_map: classifying %d x %d grid (%s, %s) ...",
                          length(xs), length(ys), axes[1L], axes[2L]))
  for (k in seq_len(n)) {
    p <- base_params
    p[[axes[1L]]] <- tab$x[k]
    p[[axes[2L]]] <- tab$y[k]
    cl <- .d39_classify_point(model, compiled, p, n_fwd, unit_tol)
    status[k] <- cl$status
    n_unst[k] <- cl$n_unstable
    if (!is.null(cl$error)) errors <- c(errors, cl$error)
  }
  tab$status     <- factor(status, levels = .d39_classes())
  tab$n_unstable <- n_unst
  tab$n_forward  <- n_fwd
  counts <- table(tab$status)
  n_feasible <- as.integer(counts[["determinate"]])
  calib <- stats::setNames(as.numeric(base_params[axes]), axes)

  res <- structure(
    list(
      mode              = "grid",
      status            = tab$status,
      class_counts      = stats::setNames(as.integer(counts), names(counts)),
      feasible_fraction = n_feasible / n,
      n_feasible        = n_feasible,
      n_total           = n,
      errors            = errors,
      drivers           = NULL,
      drivers_by_class  = list(),
      grid              = tab,
      axes              = c(x = axes[1L], y = axes[2L]),
      calibration       = calib,
      plots             = list()
    ),
    class = c("dynhr_stability_map", "list")
  )
  if (requireNamespace("ggplot2", quietly = TRUE))
    res$plots$map <- .d39_plot_map(res)
  res
}


.d39_plot_map <- function(x) {
  gg  <- .ensure_ggplot2()
  tab <- x$grid
  present <- levels(droplevels(tab$status))
  fills <- .d39_class_fills()[present]
  cal <- data.frame(x = x$calibration[[1L]], y = x$calibration[[2L]])
  counts <- x$class_counts[x$class_counts > 0L]
  sub <- paste0(x$n_total, " grid points: ",
                paste(sprintf("%s %d", names(counts), counts), collapse = ", "),
                ".  Cross = calibration point.")
  gg$ggplot(tab, gg$aes(x = .data$x, y = .data$y, fill = .data$status)) +
    gg$geom_tile(colour = NA) +
    gg$geom_point(data = cal, gg$aes(x = .data$x, y = .data$y),
                  inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.6,
                  colour = "black") +
    gg$scale_fill_manual(values = fills, breaks = present, drop = TRUE,
                         name = "BK class") +
    gg$scale_x_continuous(expand = c(0, 0)) +
    gg$scale_y_continuous(expand = c(0, 0)) +
    gg$labs(x = x$axes[["x"]], y = x$axes[["y"]],
            title = sprintf("D39 stability map: %s (x) vs %s (y)",
                            x$axes[["x"]], x$axes[["y"]]),
            subtitle = sub) +
    theme_dynhr()
}


## --------------------------------------------------------------------------
## .d39_ks_ranking: two-sample KS test per parameter
## --------------------------------------------------------------------------

.d39_ks_ranking <- function(draw_mat, feas_idx, infeas_idx, param_names) {

  n_par <- length(param_names)

  if (length(feas_idx) == 0L || length(infeas_idx) == 0L) {
    return(data.frame(
      param        = param_names,
      ks_stat      = rep(NA_real_, n_par),
      p_value      = rep(NA_real_, n_par),
      feasible_lo  = rep(NA_real_, n_par),
      feasible_hi  = rep(NA_real_, n_par),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(param_names, function(nm) {
    x_feas   <- draw_mat[feas_idx,   nm]
    x_infeas <- draw_mat[infeas_idx, nm]
    x_feas   <- x_feas[is.finite(x_feas)]
    x_infeas <- x_infeas[is.finite(x_infeas)]

    if (length(x_feas) < 2L || length(x_infeas) < 2L) {
      return(data.frame(param = nm, ks_stat = NA_real_, p_value = NA_real_,
                        feasible_lo = NA_real_, feasible_hi = NA_real_,
                        stringsAsFactors = FALSE))
    }

    ## ks.test warns about ties (point masses at clamped prior bounds); the
    ## statistic itself is still exact.
    kt <- suppressWarnings(stats::ks.test(x_feas, x_infeas))
    data.frame(
      param        = nm,
      ks_stat      = unname(kt$statistic),
      p_value      = kt$p.value,
      feasible_lo  = min(x_feas),
      feasible_hi  = max(x_feas),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[order(out$ks_stat, decreasing = TRUE, na.last = TRUE), ]
  rownames(out) <- NULL
  out
}


## --------------------------------------------------------------------------
## print method
## --------------------------------------------------------------------------

#' @export
print.dynhr_stability_map <- function(x, n_top = 10L, ...) {
  what <- if (identical(x$mode, "grid")) {
    sprintf("%s x %s grid", x$axes[["x"]], x$axes[["y"]])
  } else "prior draws"
  cat(sprintf(
    "D39 BK stability map (%s)\n  determinate: %d / %d (%.1f%%)\n",
    what, x$n_feasible, x$n_total, 100 * x$feasible_fraction
  ))
  cc <- x$class_counts
  if (length(cc) > 0L) {
    cc <- cc[cc > 0L]
    cat("  classes: ", paste(sprintf("%s=%d", names(cc), cc), collapse = ", "),
        "\n", sep = "")
  }
  if (length(x$errors) > 0L)
    cat(sprintf("  first solver error: %s\n", x$errors[1L]))
  if (identical(x$mode, "grid")) {
    cat(sprintf("  calibration: %s = %.4g, %s = %.4g\n",
                names(x$calibration)[1L], x$calibration[[1L]],
                names(x$calibration)[2L], x$calibration[[2L]]))
  }
  if (!is.null(x$drivers) && nrow(x$drivers) > 0L) {
    cat(sprintf("  top drivers, determinate vs rest (KS statistic, up to %d):\n", n_top))
    top <- utils::head(x$drivers, n_top)
    for (i in seq_len(nrow(top))) {
      r <- top[i, ]
      cat(sprintf("    %-20s  KS=%.3f  p=%.3g  determinate range=[%.4g, %.4g]\n",
                  r$param, r$ks_stat, r$p_value, r$feasible_lo, r$feasible_hi))
    }
  }
  invisible(x)
}
