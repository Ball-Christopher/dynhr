################################################################################
# dynhr_jacobian.R  --  v0.2
# Phase 2 -- Symbolic differentiation, Jacobian generation, model compilation
#
# Part of the dynhr project: a native R implementation of core Dynare
# functionality (stoch_simul, estimation) without external dependencies.
#
# Depends on: dynhr_parser.R (must be sourced first or available)
#
# Author:  dynhr project
# Date:    2026-05-02
# License: MIT
#
# Changelog v0.2:
#   - Fix: row_labels parser handles "t" -> 0 correctly (was producing NA)
################################################################################

## R/compile-model.R
## --------------------------------------------------------------------------
## Top-level compile_model() orchestrator: invokes build_static_model and
## build_dynamic_model and returns a dynhr_compiled S3 object.
## Includes the model-validation check_model() and S3 print method.
##
## Phase-1 split from jacobian-monolith.R (no logic changes).
## --------------------------------------------------------------------------

#' Check a dynhr_mod for basic consistency
#'
#' @param model A dynhr_mod object.
#' @return Invisible TRUE if no fatal issues; prints diagnostics.
#' @noRd
check_model <- function(model) {
    issues <- character(0)

    n_endo <- length(model$var_names)
    n_eq   <- length(model$equations)

    # Check equation count.
    # For OccBin models each constraint contributes one extra (bind) equation
    # alongside the matching relax equation, so n_eq = n_endo + n_constraints is
    # EXPECTED.  Demote to a message in that case; the static/dynamic row-select
    # wrappers in solve_steady_state and extract_system_matrices handle it.
    if (n_eq != n_endo) {
        n_occbin <- length(model$occbin_constraints)
        if (n_occbin > 0L && n_eq == n_endo + n_occbin) {
            message(sprintf(
                "OccBin model: %d equations, %d endogenous variables (%d constraint pair(s)); relax-regime row selection will be applied.",
                n_eq, n_endo, n_occbin))
        } else {
            issues <- c(issues,
                sprintf("Equation/variable mismatch: %d equations but %d endogenous variables",
                        n_eq, n_endo))
        }
    }

    # Blanchard-Kahn count
    nfwd <- model$n_forward + model$n_mixed
    cat("Blanchard-Kahn: ", nfwd,
        " forward-looking variable(s) (need ", nfwd,
        " explosive eigenvalue(s))\n", sep = "")

    # Check for undeclared variables in equations
    all_vars <- c(model$var_names, model$varexo_names, model$varexo_det_names)
    all_refs <- do.call(rbind, lapply(model$equations, function(eq) {
        rbind(ast_collect_variables(eq$lhs, model$local_variables),
              ast_collect_variables(eq$rhs, model$local_variables))
    }))
    if (nrow(all_refs) > 0) {
        undeclared <- setdiff(unique(all_refs$name), all_vars)
        if (length(undeclared) > 0) {
            issues <- c(issues,
                paste("Undeclared variables used in equations:",
                      paste(undeclared, collapse = ", ")))
        }
    }

    # Report
    if (length(issues) > 0) {
        cat("ISSUES:\n")
        for (iss in issues) cat("  [!] ", iss, "\n")
    } else {
        cat("Model checks passed.\n")
    }

    invisible(length(issues) == 0)
}


#' Compile a dynhr_mod into evaluable functions
#'
#' @param model       A dynhr_mod object.
#' @param verbose     Logical; if TRUE, print progress.
#' @param max_order   Maximum derivative order to compute symbolically.
#'                    1 = Jacobian only, 2 = +Hessian, 3 = +3rd deriv, etc.
#'                    Higher orders are expensive and may cause node stack
#'                    overflow on models with complex equations. Default 5.
#' @param param_deriv Controls whether the symbolic parameter-Jacobian and
#'   second-order model Hessian are compiled to support the analytic
#'   (finite-difference-free) first-order gradient.  One of:
#'   \describe{
#'     \item{\code{"auto"} (default)}{Build the parameter-derivative
#'       machinery when the model is analytically parameter-differentiable
#'       (no \code{STEADY_STATE()} references, no \code{max}/\code{min}).
#'       This is equivalent to \code{"on"} for eligible models.}
#'     \item{\code{"on"}}{Same as \code{"auto"} -- always attempt to build
#'       the analytic-gradient support.}
#'     \item{\code{"off"}}{Skip the symbolic parameter-Jacobian and the
#'       hess2 build (when not already needed by \code{max_order >= 2}).
#'       \code{dynamic$hessian2_built} and \code{.can_use_analytic_primitive_deriv()}
#'       will be \code{FALSE}; the gradient layer then falls back cleanly to
#'       central finite differences.  Use this for large models where only
#'       first-order Kalman-filter likelihood evaluation is needed and the
#'       analytic-gradient compile time is unwanted.}
#'     \item{\code{"second"}}{As \code{"on"} PLUS the second-order parameter
#'       codegen (the var-param-param tensor \eqn{\partial^3F/\partial w\,
#'       \partial\theta\,\partial\theta} and the static second-order tensors)
#'       needed for the finite-difference-free EXACT posterior Hessian
#'       (\code{solution_derivatives_2} / \code{posterior_hessian}). More
#'       expensive to compile; opt in only when the analytic exact Hessian is
#'       wanted. Forces the third-order model Hessian even at \code{max_order = 1}.}
#'   }
#' @param cache Optional on-disk compilation cache. \code{NULL} (default) or
#'   \code{FALSE} disables it. \code{TRUE} caches under a per-session temp dir;
#'   a character path caches under that directory. Symbolic compilation depends
#'   only on the model structure (equations, names, lead/lag, \code{max_order},
#'   \code{param_deriv}) and NOT on parameter values, so a cached object is
#'   reused across runs/calibrations. The cache key embeds the package version,
#'   so upgrading dynhr invalidates stale entries. Useful for batch/replication
#'   workflows that recompile the same large models repeatedly.
#' @return An object of class "dynhr_compiled".
#' @seealso \code{\link{parse_mod}}, \code{\link{solve_steady}},
#'   \code{\link{solve_perturbation}}
#' @examples
#' model    <- parse_mod(system.file("extdata/models/rbc.mod",
#'                                   package = "dynhr"), verbose = FALSE)
#' compiled <- compile_model(model, verbose = FALSE)
#' compiled
#'
#' ## Raise max_order to unlock higher-order perturbation later
#' compiled2 <- compile_model(model, max_order = 2L, verbose = FALSE)
#' compiled2$max_order
#' @export
compile_model <- function(model, verbose = FALSE, max_order = 1L,
                          param_deriv = c("auto", "on", "off", "second"),
                          cache = NULL) {
    param_deriv <- match.arg(param_deriv)
    want_param_deriv  <- (param_deriv != "off")
    want_param_deriv2 <- (param_deriv == "second")

    # ---- Optional on-disk compilation cache --------------------------------
    # Symbolic compilation is value-independent (it depends on the equation
    # ASTs / names / lead-lag structure / max_order / param_deriv, NOT on the
    # parameter VALUES), so a model can be recompiled across runs from a cached
    # object. Opt-in: cache = TRUE uses a tempdir, cache = "<path>" a directory.
    # The key includes the package version so a package update invalidates it.
    cache_file <- NULL
    if (!is.null(cache) && !isFALSE(cache)) {
        cache_dir <- if (isTRUE(cache)) file.path(tempdir(), "dynhr-compile-cache")
                     else as.character(cache)
        ver <- tryCatch(as.character(utils::packageVersion("dynhr")),
                        error = function(e) "dev")
        key <- digest::digest(list(
            equations        = model$equations,
            var_names        = model$var_names,
            varexo_names     = model$varexo_names,
            varexo_det_names = model$varexo_det_names,
            param_names      = model$param_names,
            local_variables  = model$local_variables,
            lead_lag_incidence = model$lead_lag_incidence,
            predetermined    = model$predetermined_vars,
            max_order        = as.integer(max_order),
            param_deriv      = param_deriv,
            version          = ver))
        cache_file <- file.path(cache_dir, paste0("dynhr-cm-", key, ".rds"))
        if (file.exists(cache_file)) {
            if (verbose) cat("compile_model: loaded from cache (", cache_file, ")\n")
            cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
            if (inherits(cached, "dynhr_compiled")) return(cached)
        }
    }

    if (verbose) cat("Compiling model...\n")

    if (verbose) cat("  Building static model...\n")
    static <- build_static_model(model, want_param_deriv = want_param_deriv,
                                 want_param_deriv2 = want_param_deriv2)

    if (verbose) cat("  Building dynamic model...\n")
    dynamic <- build_dynamic_model(model, max_order = max_order,
                                   want_param_deriv = want_param_deriv,
                                   want_param_deriv2 = want_param_deriv2)

    result <- list(
        model     = model,
        static    = static,
        dynamic   = dynamic,
        max_order = as.integer(max_order)
    )
    class(result) <- "dynhr_compiled"

    # Attach OccBin metadata when the model has occbin_constraints.
    # This is used by solve_steady_state and extract_system_matrices* to select
    # the relax-regime equations (regime 0 = all constraints slack).
    if (length(model$occbin_constraints) > 0L) {
        pr <- occbin_parse_bind_relax(model)
        rm <- occbin_build_regime_map(pr)
        result$occbin <- list(
            parse_result   = pr,
            regime_map     = rm,
            n_constraints  = pr$n_constraints
        )
    }

    if (!is.null(cache_file)) {
        dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
        tryCatch({
            saveRDS(result, cache_file)
            if (verbose) cat("compile_model: cached to", cache_file, "\n")
        }, error = function(e)
            warning("compile_model: failed to write cache: ", conditionMessage(e)))
    }

    if (verbose) cat("  Done.\n")
    result
}


#' Print a compiled dynhr model
#'
#' Prints a one-screen summary of a `dynhr_compiled` object: variable and
#' parameter counts, and the shapes of the static and dynamic residual /
#' Jacobian blocks.
#'
#' @param x A `dynhr_compiled` object, as returned by [compile_model()].
#' @param ... Ignored; present for S3 generic compatibility.
#'
#' @return `x`, invisibly. Called for the side effect of printing.
#'
#' @export
print.dynhr_compiled <- function(x, ...) {
    cat("=== dynhr_compiled ===\n")
    cat("Endogenous vars: ", length(x$model$var_names), "\n")
    cat("Exogenous vars:  ", length(x$model$varexo_names), "\n")
    cat("Parameters:      ", length(x$model$param_names), "\n")
    cat("Equations:       ", x$static$n_eq, "\n")
    cat("\nStatic model:\n")
    cat("  Residuals fn:  ", x$static$n_eq, " equations\n")
    cat("  Jacobian:      ", x$static$n_eq, " x ", x$static$n_endo, "\n")
    cat("\nDynamic model:\n")
    cat("  Dynamic cols:  ", x$dynamic$total_cols,
        " (", x$dynamic$n_dyn_cols, " endo + ",
        x$dynamic$n_exo, " exo)\n")
    cat("  Jacobian:      ", x$dynamic$n_eq, " x ", x$dynamic$total_cols, "\n")
    n_nonzero <- length(x$dynamic$jac_triplets)
    total <- x$dynamic$n_eq * x$dynamic$total_cols
    cat("  Non-zero Jac entries: ", n_nonzero, " / ", total,
        " (", round(100 * n_nonzero / max(total, 1), 1), "%)\n")
    invisible(x)
}
