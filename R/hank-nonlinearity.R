## R/hank-nonlinearity.R
## --------------------------------------------------------------------------
## How small is small? Two accuracy checks for the LINEAR sequence-space
## solution of a general hank_model:
##
##   hank_nonlinearity_scan()  -- for a grid of shock sizes, the exact
##       nonlinear perfect-foresight transition (hank_model_nonlinear_irf)
##       against the linear IRF (hank_model_irf) scaled by the size, per
##       output at impact / peak / cumulated. Bianchi & Kaplan (2026, NBER
##       w35311) show the linear solution overstating the consumption
##       response to a transfer by double digits at empirically relevant
##       sizes, because the constrained households' consumption function has
##       a KINK where the constraint stops binding. For the same reason the
##       comparator here is the nonlinear transition itself and never a
##       higher-order perturbation: around size 0 every higher-order term of a
##       constrained household is exactly zero (and with linearly interpolated
##       EGM policies, zero for the unconstrained too), so a local expansion
##       of any order cannot see the kink.
##
##   hank_residual_audit()  -- the model's EXACT equilibrium residuals (the
##       targets, e.g. market clearing, re-evaluated through every block's
##       nonlinear map) along any candidate path: the linear IRF at a large
##       shock, a stressed or counterfactual path. This is the off-path
##       certificate of Scheidegger & Schaab ("Equilibrium World Models",
##       arXiv 2606.23463): an approximate solution is certified only where
##       its exact residual has been checked, so check it where the answer
##       will be read. Each audit also reports the frozen-Jacobian Newton
##       step -H_U^{-1} R, the first-order distance from the candidate
##       unknowns to the exact transition (the step
##       hank_model_nonlinear_irf would take from there).
##
## Deviations in the scan are measured from the ZERO-SHOCK nonlinear solution
## (one extra solve), not from model$ss: that removes the steady-state solve
## tolerance from every gap, so the gap at size 0 is exactly 0 and the gap is
## the pure nonlinearity.
##
## The audit's paths are the aggregate unknowns (and exogenous) only: the
## nonlinear evaluator starts every household block from its steady-state
## distribution, so a counterfactual INITIAL distribution is not an auditable
## state here.
## --------------------------------------------------------------------------


#' Scan the linear-vs-nonlinear gap of a HANK model over shock sizes
#'
#' For each shock size \eqn{s}, solves the fully nonlinear perfect-foresight
#' transition to the exogenous path \eqn{Z_{ss} + s\,\mathrm{shock}} with
#' \code{\link{hank_model_nonlinear_irf}} and compares it, output by output,
#' with the linear sequence-space IRF (\code{\link{hank_model_irf}}) scaled by
#' \eqn{s}. The comparison is reported at impact, at the peak and cumulated,
#' and a size is flagged where the linear solution misstates the nonlinear one
#' by more than \code{tol} in relative terms.
#'
#' The comparator is the exact nonlinear transition, not a higher-order
#' perturbation: in a model with a borrowing constraint the constrained
#' households' consumption function is exactly linear up to a kink where the
#' constraint stops binding, so every higher-order term around \eqn{s = 0} is
#' zero for them and a local expansion of any order misses the kink (Bianchi
#' and Kaplan 2026). Expect the gap to be \eqn{O(s^2)} for small \eqn{s} and
#' asymmetric in the sign of \eqn{s}.
#'
#' Deviations of the nonlinear transition are measured from the zero-shock
#' nonlinear solution (one extra solve), so the steady-state solve tolerance
#' does not contaminate the gaps.
#'
#' @param model A \code{\link{hank_model}}.
#' @param shock Named list of deviation SHAPES, one length-\code{T_h} path per
#'   shocked exogenous (a subset of \code{model$exogenous}; the rest stay at
#'   steady state). A bare numeric vector is accepted when the model has
#'   exactly one exogenous. The shock at size \eqn{s} is \eqn{s} times this.
#' @param sizes Numeric vector of non-zero, finite shock sizes. Negative sizes
#'   are allowed (and informative: the gap is generally sign-asymmetric).
#' @param outputs Character: variables to report. Default: every unknown and
#'   every block output except the targets.
#' @param tol Relative misstatement above which a (size, output, measure) is
#'   flagged: \code{|linear - nonlinear| / |nonlinear| > tol}.
#' @param cum_horizon Integer: number of leading periods summed for the
#'   \code{"cumulative"} measure. Default \code{T_h}.
#' @param nl_tol,maxit Passed to \code{\link{hank_model_nonlinear_irf}} as
#'   \code{tol} and \code{maxit}.
#'
#' @return An object of class \code{hank_nonlinearity_scan}, a list with
#'   \describe{
#'     \item{\code{table}}{Data frame, one row per (size, output, measure),
#'       with columns \code{size}, \code{output}, \code{measure}
#'       (\code{"impact"}, \code{"peak"}, \code{"cumulative"}),
#'       \code{linear}, \code{nonlinear} (deviations), \code{gap}
#'       (\code{nonlinear - linear}), \code{ratio}
#'       (\code{nonlinear / linear}, \code{NA} when the linear value is 0),
#'       \code{rel_err} (\code{|gap| / |nonlinear|}, \code{NA} when the
#'       nonlinear value is 0), \code{flag} (\code{rel_err > tol}), and
#'       \code{t_linear}, \code{t_nonlinear} (the date of each value; the
#'       peaks may fall on different dates; \code{NA} for cumulative).}
#'     \item{\code{by_size}}{Data frame, one row per size: \code{converged},
#'       \code{iterations} and \code{max_resid} of the nonlinear solve,
#'       \code{linear_resid} (the largest exact target residual along the
#'       scaled linear path, see \code{\link{hank_residual_audit}}),
#'       \code{max_rel_err} over outputs and measures, and \code{flagged}.}
#'     \item{\code{paths}}{List with \code{linear} and \code{nonlinear}, each a
#'       named list (one per output) of \code{T_h x length(sizes)} deviation
#'       matrices.}
#'     \item{\code{sizes}, \code{tol}, \code{cum_horizon}}{The settings.}
#'   }
#'
#' @references Bianchi, J. and G. Kaplan (2026), "How Small is Small?
#'   Non-linearities in Heterogeneous Agent Models", NBER Working Paper 35311.
#' @seealso \code{\link{hank_residual_audit}}, \code{\link{hank_model_irf}},
#'   \code{\link{hank_model_nonlinear_irf}}
#' @examples
#' \donttest{
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ks  <- hank_ks_steady(hank_asset_grid(40, 30, 0), inc$Pi, inc$e,
#'                       beta = 0.95, eis = 0.5, alpha = 0.11, delta = 0.02)
#' m   <- hank_ks_model(ks, T_h = 40L)
#' sc  <- hank_nonlinearity_scan(m, list(Z = 0.7^(0:39)),
#'                               sizes = c(-0.05, 0.01, 0.05))
#' sc
#' }
#' @export
hank_nonlinearity_scan <- function(model, shock, sizes, outputs = NULL,
                                   tol = 0.05, cum_horizon = NULL,
                                   nl_tol = 1e-9, maxit = 50L) {
  if (!inherits(model, "hank_model"))
    .dynhr_abort("hank_nonlinearity_scan(): 'model' must be a hank_model.")
  T_h <- model$T_h
  dZ <- .hank_check_dZ(model, shock)
  if (!length(dZ) || all(vapply(dZ, function(v) all(v == 0), logical(1))))
    .dynhr_abort("hank_nonlinearity_scan(): 'shock' is empty or identically ",
                 "zero; supply at least one non-zero exogenous shape.")
  if (!is.numeric(sizes) || !length(sizes) || any(!is.finite(sizes)) ||
      any(sizes == 0))
    .dynhr_abort("hank_nonlinearity_scan(): 'sizes' must be a non-empty ",
                 "numeric vector of finite, non-zero shock sizes.")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0)
    .dynhr_abort("hank_nonlinearity_scan(): 'tol' must be one positive number.")
  if (is.null(cum_horizon)) cum_horizon <- T_h
  cum_horizon <- as.integer(cum_horizon)
  if (length(cum_horizon) != 1L || is.na(cum_horizon) || cum_horizon < 1L ||
      cum_horizon > T_h)
    .dynhr_abort("hank_nonlinearity_scan(): 'cum_horizon' must be an integer ",
                 "in 1..T_h (", T_h, ").")
  for (z in names(dZ))
    if (is.null(model$ss[[z]]))
      .dynhr_abort("hank_nonlinearity_scan(): model$ss has no steady-state ",
                   "level for the shocked exogenous '", z, "'.")

  exo_targets <- c(model$exogenous, model$targets)
  if (is.null(outputs))
    outputs <- setdiff(names(model$G), exo_targets)
  bad <- setdiff(outputs, names(model$G))
  if (length(bad))
    .dynhr_abort("hank_nonlinearity_scan(): unknown output(s) ",
                 paste0("'", bad, "'", collapse = ", "), "; available: ",
                 paste0("'", names(model$G), "'", collapse = ", "), ".")

  ## Linear: one unit-size solve, scaled per size (it is linear).
  lin1 <- hank_model_irf(model, dZ)
  ## Nonlinear baseline: the zero-shock transition. Deviations are measured
  ## from it so the steady-state solve tolerance cancels out of every gap.
  base <- hank_model_nonlinear_irf(model, list(), tol = nl_tol, maxit = maxit)
  if (!isTRUE(base$converged))
    .dynhr_abort("hank_nonlinearity_scan(): the zero-shock nonlinear ",
                 "transition did not converge (max residual ",
                 format(base$max_resid, digits = 3), "); the steady state in ",
                 "model$ss is not a fixed point of the model's blocks.")
  U_base <- base[model$unknowns]
  Z_base <- base[model$exogenous]

  n_s <- length(sizes)
  lin_paths <- setNames(lapply(outputs, function(v) matrix(NA_real_, T_h, n_s)),
                        outputs)
  nl_paths  <- lin_paths
  by_size <- data.frame(size = sizes, converged = NA, iterations = NA_integer_,
                        max_resid = NA_real_, linear_resid = NA_real_,
                        max_rel_err = NA_real_, flagged = NA)
  rows <- vector("list", n_s)

  measure <- function(x, what) {
    if (what == "impact") return(c(x[1L], 1))
    if (what == "peak") { k <- which.max(abs(x)); return(c(x[k], k)) }
    c(sum(x[seq_len(cum_horizon)]), NA_real_)
  }

  for (k in seq_len(n_s)) {
    s <- sizes[k]
    Z_paths <- lapply(names(dZ), function(z) model$ss[[z]] + s * dZ[[z]])
    names(Z_paths) <- names(dZ)
    nl <- hank_model_nonlinear_irf(model, Z_paths, tol = nl_tol, maxit = maxit)
    if (!isTRUE(nl$converged))
      .dynhr_warn("hank_nonlinearity_scan(): the nonlinear transition at ",
                  "size ", format(s), " did not converge in ", maxit,
                  " iterations (max residual ", format(nl$max_resid, digits = 3),
                  "); its row is reported but is not the exact transition.")
    by_size$converged[k]  <- isTRUE(nl$converged)
    by_size$iterations[k] <- nl$iterations
    by_size$max_resid[k]  <- nl$max_resid

    ## Exact residual of the scaled linear path, around the same baseline.
    U_lin <- lapply(model$unknowns, function(u) U_base[[u]] + s * lin1[[u]])
    names(U_lin) <- model$unknowns
    Z_lin <- Z_base
    for (z in names(dZ)) Z_lin[[z]] <- Z_paths[[z]]
    aud <- .hank_audit_levels(model, U_lin, Z_lin)
    by_size$linear_resid[k] <- max(vapply(aud$resid, function(r) max(abs(r)),
                                          numeric(1)))

    rk <- vector("list", length(outputs))
    for (j in seq_along(outputs)) {
      v <- outputs[j]
      lv <- s * lin1[[v]]
      nv <- nl[[v]] - base[[v]]
      lin_paths[[v]][, k] <- lv
      nl_paths[[v]][, k]  <- nv
      mm <- lapply(c("impact", "peak", "cumulative"), function(w) {
        a <- measure(lv, w); b <- measure(nv, w)
        data.frame(size = s, output = v, measure = w,
                   linear = a[1L], nonlinear = b[1L],
                   t_linear = a[2L], t_nonlinear = b[2L])
      })
      rk[[j]] <- do.call(rbind, mm)
    }
    rows[[k]] <- do.call(rbind, rk)
  }

  tab <- do.call(rbind, rows)
  tab$gap     <- tab$nonlinear - tab$linear
  tab$ratio   <- ifelse(tab$linear != 0, tab$nonlinear / tab$linear, NA_real_)
  tab$rel_err <- ifelse(tab$nonlinear != 0, abs(tab$gap) / abs(tab$nonlinear),
                        NA_real_)
  tab$flag    <- !is.na(tab$rel_err) & tab$rel_err > tol
  tab$t_linear    <- as.integer(tab$t_linear)
  tab$t_nonlinear <- as.integer(tab$t_nonlinear)
  tab <- tab[, c("size", "output", "measure", "linear", "nonlinear", "gap",
                 "ratio", "rel_err", "flag", "t_linear", "t_nonlinear")]
  rownames(tab) <- NULL
  for (k in seq_len(n_s)) {
    re <- tab$rel_err[tab$size == sizes[k]]
    by_size$max_rel_err[k] <- if (all(is.na(re))) NA_real_ else max(re, na.rm = TRUE)
    by_size$flagged[k] <- any(tab$flag[tab$size == sizes[k]])
  }

  structure(list(table = tab, by_size = by_size,
                 paths = list(linear = lin_paths, nonlinear = nl_paths),
                 sizes = sizes, tol = tol, cum_horizon = cum_horizon,
                 shocked = names(dZ)),
            class = "hank_nonlinearity_scan")
}


#' Print a HANK nonlinearity scan
#' @param x A \code{hank_nonlinearity_scan} from
#'   \code{\link{hank_nonlinearity_scan}}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.hank_nonlinearity_scan <- function(x, ...) {
  cat("=== hank_nonlinearity_scan ===\n")
  cat(sprintf("Shocked     : %s\n", paste(x$shocked, collapse = ", ")))
  cat(sprintf("Outputs     : %s\n",
              paste(unique(x$table$output), collapse = ", ")))
  cat(sprintf("Flag        : |linear - nonlinear| / |nonlinear| > %g\n", x$tol))
  bs <- x$by_size
  cat("Per size (nonlinear solve, exact residual of the linear path):\n")
  for (k in seq_len(nrow(bs)))
    cat(sprintf("  size %+10.4g  %s  max rel err %8.3g  linear-path resid %9.3g%s\n",
                bs$size[k],
                if (isTRUE(bs$converged[k])) "converged    " else "NOT converged",
                bs$max_rel_err[k], bs$linear_resid[k],
                if (isTRUE(bs$flagged[k])) "  FLAGGED" else ""))
  invisible(x)
}


#' Audit the exact equilibrium residuals along a candidate HANK path
#'
#' Re-evaluates every block of a \code{\link{hank_model}} NONLINEARLY along
#' candidate paths for the unknowns (and exogenous) and reports the
#' equilibrium residuals -- the model's targets, e.g. asset-market clearing --
#' that the candidate leaves. The exact nonlinear transition
#' (\code{\link{hank_model_nonlinear_irf}}) has residuals at the solver
#' tolerance; any approximate solution (the linear IRF at a large shock, a
#' hand-built stressed or counterfactual path) is certified only as far as
#' this residual is small (the off-path certificate of Scheidegger and Schaab
#' 2026).
#'
#' Each audit also reports the frozen-Jacobian Newton step
#' \eqn{-H_U^{-1} R}: the first-order correction from the candidate unknowns
#' to the exact transition, i.e. an estimate of the candidate's error in
#' units of the unknowns. For the linear IRF at size \eqn{s} the residual is
#' \eqn{O(s^2)} and the step matches the true error to relative \eqn{O(s)}.
#'
#' The household blocks start from their steady-state distribution; the
#' audited paths are aggregate paths only.
#'
#' @param model A \code{\link{hank_model}}.
#' @param paths Named list of length-\code{T_h} paths containing every
#'   \code{model$unknowns}; exogenous paths are optional (absent ones are held
#'   at steady state) and any other entries are ignored, so the whole output
#'   of \code{\link{hank_model_irf}} or \code{\link{hank_model_nonlinear_irf}}
#'   can be passed as is. Alternatively a NAMED list of such lists, one per
#'   scenario.
#' @param type \code{"deviation"} (paths are deviations from \code{model$ss},
#'   as returned by \code{\link{hank_model_irf}}) or \code{"level"} (as
#'   returned by \code{\link{hank_model_nonlinear_irf}}).
#'
#' @return An object of class \code{hank_residual_audit}, a list with
#'   \describe{
#'     \item{\code{summary}}{Data frame, one row per (scenario, variable):
#'       \code{scenario}, \code{variable}, \code{kind}
#'       (\code{"residual"} for a target, \code{"newton_step"} for an
#'       unknown), \code{max_abs} and the date \code{t_max} of the maximum.}
#'     \item{\code{scenarios}}{Named list, per scenario: \code{resid} (named
#'       list of target residual paths) and \code{newton_step} (named list of
#'       per-unknown correction paths).}
#'     \item{\code{max_resid}}{Largest absolute target residual over all
#'       scenarios.}
#'   }
#'
#' @references Scheidegger, S. and A. Schaab (2026), "Equilibrium World
#'   Models", arXiv:2606.23463.
#' @seealso \code{\link{hank_nonlinearity_scan}}
#' @examples
#' \donttest{
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ks  <- hank_ks_steady(hank_asset_grid(40, 30, 0), inc$Pi, inc$e,
#'                       beta = 0.95, eis = 0.5, alpha = 0.11, delta = 0.02)
#' m   <- hank_ks_model(ks, T_h = 40L)
#' dZ  <- list(Z = 0.05 * 0.7^(0:39))
#' hank_residual_audit(m, list(linear = hank_model_irf(m, dZ)))
#' }
#' @export
hank_residual_audit <- function(model, paths, type = c("deviation", "level")) {
  if (!inherits(model, "hank_model"))
    .dynhr_abort("hank_residual_audit(): 'model' must be a hank_model.")
  type <- match.arg(type)
  if (!is.list(paths) || !length(paths))
    .dynhr_abort("hank_residual_audit(): 'paths' must be a non-empty list.")
  multi <- all(vapply(paths, is.list, logical(1)))
  if (multi) {
    if (is.null(names(paths)) || any(!nzchar(names(paths))))
      .dynhr_abort("hank_residual_audit(): a list of scenarios must be NAMED.")
    scen <- paths
  } else {
    scen <- list(path = paths)
  }

  T_h <- model$T_h; ss <- model$ss
  out <- list(); rows <- list()
  for (nm in names(scen)) {
    p <- scen[[nm]]
    miss <- setdiff(model$unknowns, names(p))
    if (length(miss))
      .dynhr_abort("hank_residual_audit(): scenario '", nm, "' has no path ",
                   "for unknown(s) ", paste0("'", miss, "'", collapse = ", "),
                   ".")
    src <- intersect(c(model$unknowns, model$exogenous), names(p))
    short <- src[vapply(p[src], length, integer(1)) != T_h]
    if (length(short))
      .dynhr_abort("hank_residual_audit(): scenario '", nm, "': path(s) ",
                   paste0("'", short, "'", collapse = ", "),
                   " must have length T_h = ", T_h, ".")
    lev <- function(v) {
      x <- p[[v]]
      if (is.null(x)) return(rep(ss[[v]], T_h))
      if (type == "deviation") ss[[v]] + x else x
    }
    U <- setNames(lapply(model$unknowns, lev), model$unknowns)
    Z <- setNames(lapply(model$exogenous, lev), model$exogenous)
    a <- .hank_audit_levels(model, U, Z)
    out[[nm]] <- a
    for (t in names(a$resid)) {
      r <- abs(a$resid[[t]]); k <- which.max(r)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = nm, variable = t, kind = "residual",
        max_abs = r[k], t_max = k)
    }
    for (u in names(a$newton_step)) {
      d <- abs(a$newton_step[[u]]); k <- which.max(d)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = nm, variable = u, kind = "newton_step",
        max_abs = d[k], t_max = k)
    }
  }
  summ <- do.call(rbind, rows); rownames(summ) <- NULL
  structure(list(summary = summ, scenarios = out,
                 max_resid = max(summ$max_abs[summ$kind == "residual"])),
            class = "hank_residual_audit")
}


#' Print a HANK residual audit
#' @param x A \code{hank_residual_audit} from \code{\link{hank_residual_audit}}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.hank_residual_audit <- function(x, ...) {
  cat("=== hank_residual_audit ===\n")
  s <- x$summary
  for (k in seq_len(nrow(s)))
    cat(sprintf("  %-12s %-14s %-11s max|.| %10.3e at t = %d\n",
                s$scenario[k], s$variable[k], s$kind[k], s$max_abs[k],
                s$t_max[k]))
  invisible(x)
}


#' Exact target residuals and the frozen-Jacobian Newton step at LEVEL paths
#'
#' @param model A \code{\link{hank_model}}.
#' @param U,Z Named lists of level paths for every unknown and exogenous.
#' @return List with \code{resid} (named per target) and \code{newton_step}
#'   (named per unknown; \eqn{-H_U^{-1} R}, stacked as in
#'   \code{\link{hank_model_nonlinear_irf}}).
#' @keywords internal
.hank_audit_levels <- function(model, U, Z) {
  T_h <- model$T_h
  vals <- .hank_model_eval(model, c(U, Z))
  resid <- vals[model$targets]
  fac <- .hank_ge_factor(model$H_U)
  step <- -as.numeric(.hank_ge_solve(fac, model$H_U,
                                     do.call(c, resid)))
  newton <- lapply(seq_along(model$unknowns), function(k)
    step[((k - 1L) * T_h + 1L):(k * T_h)])
  names(newton) <- model$unknowns
  list(resid = resid, newton_step = newton)
}
