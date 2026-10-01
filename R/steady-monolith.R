################################################################################
# dynhr_steady.R  v0.2
# Phase 3 -- Steady state solvers (NA-robust)
# Depends on: dynhr_parser.R, dynhr_jacobian.R
################################################################################

# Dependencies: dynhr_parser.R, dynhr_jacobian.R (sourced via dynhr_load_all.R)

# ============================================================================
# 1. ANALYTICAL STEADY STATE SOLVER
# ============================================================================

## parse() is a pure function of its text, but it is comparatively expensive.
## solve_steady_state_analytical runs once per MCMC draw and re-parsed every
## steady_state_model assignment each time (~18 lines x thousands of draws for
## fs2000 -> ~5% of per-draw). Memoise the parsed expressions keyed on the text
## so each distinct line is parsed once per session. eval() (the actual
## computation) is unchanged, so results are bit-identical.
##
## A-SEC: the two parse caches (`.ssm_expr_cache` here and
## `.dynhr_expr_cache` in stochsimul-monolith.R) are the security CHOKE POINT
## for .mod text that is re-evaluated at solve time.  An expression is
## AST-checked against the .mod allowlist when it is FIRST inserted into a
## cache -- a disallowed call aborts with `dynhr_error_unsafe_mod_expression`
## before anything is evaluated -- so the check costs nothing per MCMC draw.
## As defence in depth the evaluation environment is a child of an
## emptyenv()-rooted allowlist sandbox (built once per session), never of
## baseenv(): a programmatically built model that skipped parse-time
## validation still cannot reach system(), file.create(), ...
.ssm_expr_cache <- new.env(parent = emptyenv())
.dynhr_mod_sandbox_cache <- new.env(parent = emptyenv())

#' Session-wide allowlist sandbox used as the PARENT of solve-time .mod
#' evaluation environments (steady_state_model, shocks-block *_expr).
#' Matrix allowlist plus `=` / `<-`, which only a top-level
#' `name = rhs` steady_state_model statement may use (the call check rejects
#' them anywhere else).
#' @noRd
.dynhr_mod_sandbox_parent <- function() {
  env <- .dynhr_mod_sandbox_cache$env
  if (is.null(env)) {
    env <- .dynhr_sandbox_env(NULL, c(.dynhr_safe_matrix_fn_names, "=", "<-"))
    assign("env", env, envir = .dynhr_mod_sandbox_cache)
  }
  env
}

#' Abort unless every call in a parsed .mod expression is allowlisted.
#'
#' @param ex      A parsed `expression` vector (from `parse(text = )`).
#' @param text    The source text, for the error message.
#' @param context Short label naming the block.
#' @param assignment When TRUE a top-level `name = rhs` / `name <- rhs` with a
#'   SYMBOL left-hand side is permitted and only `rhs` is checked.
#' @noRd
.dynhr_check_mod_expr <- function(ex, text, context, assignment = FALSE) {
  allowed <- .dynhr_safe_matrix_fn_names
  for (e in as.list(ex)) {
    if (assignment && is.call(e) && length(e) == 3L &&
        (identical(e[[1L]], as.name("=")) || identical(e[[1L]], as.name("<-"))) &&
        is.symbol(e[[2L]]))
      e <- e[[3L]]
    bad <- setdiff(.dynhr_expr_calls(e), allowed)
    if (length(bad) > 0L)
      .dynhr_abort(
        "dynhr refuses to evaluate ", context, " `", text, "`: it calls `",
        bad[1L], "`, which is not on the .mod expression allowlist. ",
        "Evaluating a model must never execute arbitrary code; move any ",
        "computation that needs `", bad[1L], "` out of the model.",
        class = "dynhr_error_unsafe_mod_expression")
  }
  invisible(TRUE)
}

#' Error-handler body for call sites that map evaluation errors to a fallback
#' value: re-raise an unsafe-expression refusal, return `value` otherwise.
#'
#' Must be used INSIDE the `error =` handler.  A separate
#' `dynhr_error_unsafe_mod_expression = function(e) stop(e)` handler in the
#' same tryCatch() would NOT work: tryCatch nests its handlers so the later
#' `error =` handler encloses the earlier one and would swallow the re-raise.
#' @noRd
.dynhr_reraise_unsafe <- function(e, value) {
  if (inherits(e, "dynhr_error_unsafe_mod_expression")) stop(e)
  value
}

.cached_parse <- function(txt) {
  e <- .ssm_expr_cache[[txt]]
  if (is.null(e)) {
    e <- parse(text = txt)
    .dynhr_check_mod_expr(e, txt, "the steady_state_model statement",
                          assignment = TRUE)
    assign(txt, e, envir = .ssm_expr_cache)
  }
  e
}

#' Fresh steady_state_model evaluation env: `params` bound in a child of the
#' session allowlist sandbox (never baseenv()).
#' @noRd
.ssm_eval_env <- function(params) {
  env <- new.env(parent = .dynhr_mod_sandbox_parent())
  for (nm in names(params)) assign(nm, params[[nm]], envir = env)
  env
}

solve_steady_state_analytical <- function(model, params) {
  ssm <- model$steady_state_model
  if (is.null(ssm) || length(ssm) == 0) return(NULL)
  env <- .ssm_eval_env(params)
  for (assignment in ssm) {
    val <- eval(.cached_parse(assignment$text), envir = env)
    if (is.numeric(val) && length(val) == 1)
      assign(assignment$name, val, envir = env)
  }
  ss <- setNames(numeric(length(model$var_names)), model$var_names)
  for (nm in model$var_names) {
    if (exists(nm, envir = env, inherits = FALSE))
      ss[[nm]] <- get(nm, envir = env, inherits = FALSE)
  }
  # Auto-set AUX_LEAD_* variables to their parent variable's SS value.
  # AUX_LEAD vars are created for lead > 1 (e.g. c(+2) => AUX_LEAD_c_1).
  # At steady state, leads equal their SS, so AUX_LEAD_c_1 = c.
  aux_vars <- grep("^AUX_LEAD_", model$var_names, value = TRUE)
  for (aux_nm in aux_vars) {
    if (ss[[aux_nm]] == 0 || is.na(ss[[aux_nm]])) {
      # Extract parent var name: AUX_LEAD_<parent>_<k>
      parent <- sub("^AUX_LEAD_", "", aux_nm)
      parent <- sub("_[0-9]+$", "", parent)
      if (parent %in% model$var_names && is.finite(ss[[parent]])) {
        ss[[aux_nm]] <- ss[[parent]]
      }
    }
  }
  # Also overwrite parameter values computed in the steady_state_model
  # block (e.g. gammax, delta, beta, r_star, b_star).  These are computed
  # from calibrated parameters and reflect the correct steady state, even
  # if the model file has dummy initial values (like r_star = 0).
  # We also include ALL SSM-assigned names (even those not in param_names),
  # because model equations may reference SSM intermediate values like
  # RSS, YSS, PSS, CSS as "parameters" (the parser treats unknown
  # identifiers as parameters).  Without these, params["RSS"] is NA and
  # the residual function returns NaN.
  updated_params <- params
  ssm_names <- unique(vapply(ssm, `[[`, character(1), "name"))
  all_param_names <- union(names(params),
    c(model$param_names %||% character(0), ssm_names))
  for (nm in all_param_names) {
    if (exists(nm, envir = env, inherits = FALSE)) {
      val <- get(nm, envir = env, inherits = FALSE)
      if (is.numeric(val) && length(val) == 1 && is.finite(val))
        updated_params[[nm]] <- val
    }
  }
  list(ss = ss, params = updated_params)
}

# ============================================================================
# 2. NEWTON-RAPHSON WITH LINE SEARCH
# ============================================================================

#' @references
#'   Juillard, M. (1996). Dynare: A program for the resolution and simulation
#'     of dynamic models with forward variables through the use of a relaxation
#'     algorithm. \emph{CEPREMAP Working Paper}.
#'   Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007).
#'     \emph{Numerical Recipes: The Art of Scientific Computing} (3rd ed.).
#'     Cambridge University Press. Section 9.7 (Newton with backtracking).
#' @param exo_init Optional named numeric vector of exogenous variable values
#'   to hold fixed when evaluating the static residuals (i.e. the steady-state
#'   level of the exogenous variables, from an \code{initval} block).  Names
#'   must match the model's \code{varexo} declarations; unnamed entries are
#'   matched in declaration order.  When \code{NULL} (default) all exogenous
#'   variables are treated as zero, which is the Dynare convention for models
#'   whose steady state is defined at \eqn{\varepsilon = 0}.
#' @export
solve_steady <- function(compiled, params, y0 = NULL,
                            endo_names = NULL, exo_names = NULL,
                            exo_init = NULL,
                            max_iter = 1000L, tol = 1e-10,
                            verbose = FALSE) {
  ## Auto-derive endo_names/exo_names from compiled model when not provided.
  if (is.null(endo_names) && inherits(compiled, "dynhr_compiled"))
    endo_names <- compiled$model$var_names
  if (is.null(exo_names) && inherits(compiled, "dynhr_compiled"))
    exo_names <- compiled$model$varexo_names

  ## Delegate to solve_steady_state for consistency; this preserves the
  ## public signature while benefiting from analytical-SS and linear-model
  ## fast paths in the internal function.
  model <- if (inherits(compiled, "dynhr_compiled")) compiled$model else compiled
  ss <- solve_steady_state(model, compiled, params, y0 = y0,
                            exo_init = exo_init,
                            method = "auto", max_iter = max_iter,
                            tol = tol, verbose = verbose)
  ## Remap solve_steady_state output fields to solve_steady convention.
  ## Include updated_params if the analytical SS filled in missing values.
  result <- list(values = ss$values, residuals = ss$residuals,
                 converged = ss$converged, iterations = ss$iterations,
                 max_residual = ss$max_residual,
                 method = ss$method_used %||% "Newton")
  if (!is.null(ss$params))
    result$params <- ss$params
  structure(result, class = c("dynhr_steady", "list"))
}

# ============================================================================
# 3. NLEQSLV WRAPPER
# ============================================================================

## Build the exogenous vector used by static evaluators.
## When exo_init is NULL, returns a zero vector (Dynare convention).
## When exo_init is a named numeric, values are matched by name; unmapped
## names stay zero.  When exo_init is an unnamed numeric of length n_exo,
## values are matched in declaration order.
.resolve_exo_init <- function(exo_init, exo_names) {
  x <- setNames(rep(0, length(exo_names)), exo_names)
  if (is.null(exo_init) || length(exo_init) == 0) return(x)
  exo_init <- as.numeric(exo_init)
  nms <- names(exo_init)
  if (!is.null(nms) && length(nms) > 0 && any(nzchar(nms))) {
    # Named: match by name
    common <- intersect(nms, exo_names)
    x[common] <- exo_init[common]
  } else {
    # Unnamed: match in declaration order
    n_fill <- min(length(exo_init), length(exo_names))
    x[seq_len(n_fill)] <- exo_init[seq_len(n_fill)]
  }
  x
}

solve_ss_nleqslv <- function(compiled, params, y0 = NULL,
                             endo_names = NULL, exo_names = NULL,
                             exo_init = NULL,
                             tol = 1e-10, verbose = FALSE) {
  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    if (verbose) .dynhr_cat("  nleqslv not available\n"); return(NULL)
  }
  if (inherits(compiled, "dynhr_compiled")) {
    res_fn <- compiled$static$residuals_fn
    jac_fn <- compiled$static$jacobian_fn
    if (is.null(endo_names)) endo_names <- compiled$model$var_names
    if (is.null(exo_names))  exo_names  <- compiled$model$varexo_names
  } else { res_fn <- compiled$residuals_fn; jac_fn <- compiled$jacobian_fn }
  n <- length(endo_names)
  x <- .resolve_exo_init(exo_init, exo_names)
  if (is.null(y0)) { y <- setNames(rep(0.5, n), endo_names)
  } else { y <- setNames(as.numeric(y0[endo_names]), endo_names); y[is.na(y)] <- 0.5 }
  fn_w <- function(yy) { names(yy) <- endo_names; res_fn(yy, x, params, yy) }
  jac_w <- function(yy) { names(yy) <- endo_names; jac_fn(yy, x, params, yy) }

  ## nleqslv::nleqslv() errors outright if the residuals or Jacobian at the
  ## starting guess are non-finite (e.g. log() of a non-positive parameter
  ## value during mode-finding/optimisation). Detect that case up front and
  ## report it as non-convergence instead, so callers (solve_steady_state's
  ## method cascade, make_log_posterior's -Inf-on-infeasible-draw path) can
  ## handle it the normal way rather than the whole estimation crashing.
  f0 <- fn_w(y)
  if (!all(is.finite(f0)))
    return(list(values = y, residuals = f0, converged = FALSE,
                iterations = 0L, max_residual = NA_real_, method = "nleqslv"))
  j0 <- jac_w(y)
  if (!all(is.finite(j0)))
    return(list(values = y, residuals = f0, converged = FALSE,
                iterations = 0L, max_residual = NA_real_, method = "nleqslv"))

  sol <- nleqslv::nleqslv(y, fn_w, jac = jac_w,
                          control = list(ftol = tol, xtol = tol))
  r_f <- fn_w(sol$x); names(sol$x) <- endo_names
  list(values = sol$x, residuals = r_f, converged = sol$termcd %in% c(1,2),
       iterations = sol$iter, max_residual = max(abs(r_f)), method = "nleqslv")
}

# ============================================================================
# 4. OPTIM FALLBACK
# ============================================================================

solve_ss_optim <- function(compiled, params, y0 = NULL,
                           endo_names = NULL, exo_names = NULL,
                           exo_init = NULL,
                           max_iter = 5000L, tol = 1e-10, verbose = FALSE) {
  if (inherits(compiled, "dynhr_compiled")) {
    res_fn <- compiled$static$residuals_fn
    if (is.null(endo_names)) endo_names <- compiled$model$var_names
    if (is.null(exo_names))  exo_names  <- compiled$model$varexo_names
  } else { res_fn <- compiled$residuals_fn }
  n <- length(endo_names)
  x <- .resolve_exo_init(exo_init, exo_names)
  if (is.null(y0)) { y <- setNames(rep(0.5, n), endo_names)
  } else { y <- setNames(as.numeric(y0[endo_names]), endo_names); y[is.na(y)] <- 0.5 }
  obj <- function(yy) {
    names(yy) <- endo_names
    r <- res_fn(yy, x, params, yy)
    if (any(!is.finite(r))) 1e20 else sum(r^2)
  }
  sol <- optim(y, obj, method = "L-BFGS-B",
               control = list(maxit = max_iter, factr = tol/.Machine$double.eps))
  y_sol <- sol$par; names(y_sol) <- endo_names
  r_f <- res_fn(y_sol, x, params, y_sol)
  list(values = y_sol, residuals = r_f,
       converged = max(abs(r_f)) < tol * 100,
       iterations = sol$counts[["function"]],
       max_residual = max(abs(r_f)), method = "optim")
}

# ============================================================================
# 5. HOMOTOPY SOLVER
# ============================================================================


# ============================================================================
# 6. MAIN STEADY STATE INTERFACE
# ============================================================================

#' Solve a model's deterministic steady state
#'
#' Finds the vector of endogenous values at which every static equation
#' residual is zero, given a parameter vector. This is the input the
#' perturbation and state-space routines expect, so it is the usual first step
#' when rebuilding a solution at a posterior mode or at a swept parameter.
#'
#' By default the solver runs a CASCADE, taking the first method that
#' converges: an analytical solve (using the model's own steady-state block if
#' it has one), then Newton, then \code{nleqslv}, then \code{optim}. Setting
#' \code{method} to one of those names pins it to that method and fails rather
#' than falling through, which is what you want in a replication script where a
#' silent change of method between runs would be a change in the numbers.
#'
#' A model declared \code{model(linear)} is a special case handled up front: it
#' is expressed in deviations, so its steady state is zero BY CONVENTION and
#' the static residual there is the linearisation point rather than an error.
#'
#' @param model A parsed model (see \code{\link{parse_mod}}).
#' @param compiled Optional compiled model from \code{\link{compile_model}}.
#'   \code{NULL} (default) compiles one internally; pass an existing object to
#'   avoid recompiling in a loop over parameters.
#' @param params Named numeric parameter vector. \code{NULL} (default) uses
#'   \code{model$param_values}. Errors if neither supplies values.
#' @param y0 Optional named numeric starting guess for the endogenous
#'   variables. \code{NULL} (default) uses the model's \code{initval} block
#'   where present, filling the remainder with 0.5 (or 1 when there is no
#'   \code{initval} block at all).
#' @param method One of \code{"auto"} (default; the cascade described above),
#'   \code{"analytical"}, \code{"Newton"}, \code{"nleqslv"} or \code{"optim"}.
#'   Anything other than \code{"auto"} pins the solver to that method.
#' @param exo_init Optional starting values for the EXOGENOUS variables,
#'   defaulting to all zero. A NAMED vector is matched by name; an UNNAMED one
#'   is matched in declaration order.
#' @param max_iter Maximum iterations for the iterative methods.
#' @param tol Convergence tolerance on the static residuals.
#' @param verbose Logical: print the method cascade's progress.
#' @param max_attempts Integer guard on how many times the whole cascade may be
#'   retried; each attempt tries the remaining methods in order. Raising it does
#'   not make a genuinely infeasible calibration solvable.
#'
#' @return An object of class \code{dynhr_steady} (a list) with the named
#'   steady-state vector in BOTH \code{ss} and \code{values} (the same object
#'   under two names, kept for back-compatibility), the static \code{residuals},
#'   \code{converged}, \code{max_residual}, and \code{method_used} recording
#'   which method actually succeeded. **Check \code{converged} before using the
#'   result**: a non-converged solve is returned rather than raised, so that
#'   callers sweeping a parameter can inspect the failure instead of aborting
#'   the sweep.
#'
#' @seealso \code{\link{compile_model}}, \code{\link{solve_perturbation}}
#' @examples
#' \donttest{
#' mod <- system.file("extdata", "models", "rbc", "rbc.mod", package = "dynhr")
#' if (nzchar(mod)) {
#'   m  <- parse_mod(mod)
#'   ss <- solve_steady_state(m)
#'   ss$converged
#' }
#' }
#' @export
solve_steady_state <- function(model, compiled = NULL, params = NULL,
                               y0 = NULL, method = "auto",
                               exo_init = NULL,
                               max_iter = 1000L, tol = 1e-10,
                               verbose = FALSE, max_attempts = 3L) {
  if (is.null(params)) params <- model$param_values
  if (length(params) == 0)
    stop("No parameter values available. Supply params argument.")
  n_na <- sum(is.na(params))
  if (n_na > 0 && verbose) {
    .dynhr_cat("WARNING:", n_na, "of", length(params), "parameters have NA values.\n")
    .dynhr_cat("  NA params:", paste(head(names(params)[is.na(params)], 10),
                              collapse = ", "),
        if (n_na > 10) ", ..." else "", "\n")
  }
  if (is.null(compiled)) {
    compiled <- compile_model(model, verbose = verbose)
  }
  endo_names <- model$var_names
  exo_names  <- model$varexo_names
  n <- length(endo_names)

  # ---- OccBin relax-regime row selection ------------------------------------
  # When the model has OccBin constraints, n_eq > n_endo (each constraint adds
  # one bind-equation alongside the relax-equation).  The SS solve needs the
  # square relax-regime subsystem (regime 0 = all constraints slack).
  # We wrap the static residuals_fn and jacobian_fn to select those rows,
  # giving downstream solvers (Newton, nleqslv, optim) a square n_endo system.
  # This is a LOCAL rebinding only - the compiled object on disk is unchanged.
  if (inherits(compiled, "dynhr_compiled") &&
      !is.null(compiled$occbin) &&
      compiled$static$n_eq > length(endo_names)) {
    relax_rows <- compiled$occbin$regime_map[[1L]]$eq_indices  # regime 0
    if (verbose)
      .dynhr_inform(sprintf("OccBin SS: selecting %d relax-regime rows from %d equations",
                      length(relax_rows), compiled$static$n_eq))
    orig_static <- compiled$static
    wrapped_res_fn <- local({
      rfn <- orig_static$residuals_fn; rr <- relax_rows
      function(y, x, params, ss) rfn(y, x, params, ss)[rr]
    })
    wrapped_jac_fn <- local({
      jfn <- orig_static$jacobian_fn; rr <- relax_rows
      function(y, x, params, ss) jfn(y, x, params, ss)[rr, , drop = FALSE]
    })
    compiled_local <- compiled
    compiled_local$static <- orig_static
    compiled_local$static$residuals_fn <- wrapped_res_fn
    compiled_local$static$jacobian_fn  <- wrapped_jac_fn
    compiled <- compiled_local
  }
  
  # Linear model fast path.
  # Core dynamics of a linear model are zero at steady state, but observation /
  # measurement equations (dy = y - y(-1) + ctrend; pinfobs = pinf + constepinf;
  # etc.) may set observable variables to non-zero constants via a
  # steady_state_model block.  Run the analytical solver first when such a
  # block exists so those constants are respected.  The convergence flag is
  # honest: it requires max|residual| < tol * 100.
  is_linear <- isTRUE(model$model_options$linear)
  if (is_linear) {
    if (verbose) .dynhr_cat("Linear model detected: checking for steady_state_model block.\n")
    ss <- setNames(rep(0, n), endo_names)
    ss_params <- params  # may be updated by analytical path
    has_ssm <- length(model$steady_state_model) > 0
    if (has_ssm) {
      if (verbose) .dynhr_cat("  steady_state_model block found - running analytical SS for obs constants.\n")
      ss_a_result <- solve_steady_state_analytical(model, params)
      if (!is.null(ss_a_result) && all(is.finite(ss_a_result$ss))) {
        ss        <- ss_a_result$ss
        ss_params <- ss_a_result$params
      } else if (verbose) {
        .dynhr_cat("  Analytical SS did not return finite values; falling back to zero.\n")
      }
    }
    x0 <- .resolve_exo_init(exo_init, exo_names)
    r  <- compiled$static$residuals_fn(ss, x0, ss_params, ss)
    max_r <- max(abs(r))

    ## ---- Affine constants (0.9.4) -----------------------------
    ## `model(linear)` does NOT imply "written in deviations from its own
    ## steady state". `r = rbar + phi*x` with rbar = 0.03 is a perfectly legal
    ## linear model whose steady state is r = 0.03; the pre-0.9.4 shortcut
    ## zero-filled it and reported "converged", so every downstream level
    ## (OBC bounds, hist decomp, forecasts) was off by the constant.
    ##
    ## The static system of a linear model is AFFINE, F(y) = A y + c, with
    ## c = F(0) = `r` above and A the (constant) static Jacobian. So the
    ## steady state is the exact linear solve A y = -c. When c = 0 this
    ## returns y = 0 and the old deviation convention is reproduced bit for
    ## bit; the branch is only entered when c is non-zero.
    ##
    ## A SINGULAR static Jacobian (a unit root: `y = y(-1) + e` contributes a
    ## zero row and column) used to skip the solve, keep SS = 0 and report
    ## converged = TRUE with the constant still sitting in the residual --
    ## `dy = y - y(-1) + 0.4` came back dy = 0, residual -0.4. The singular
    ## case is now solved in the least-squares / minimum-norm sense through
    ## the SVD pseudo-inverse, which is exact whenever A y = -c is CONSISTENT
    ## (dy = 0.4, the indeterminate unit-root level at its minimum-norm 0).
    ## An INCONSISTENT system (a drift, `y = y(-1) + 0.1 + e`, has no steady
    ## state) leaves a residual, and is reported as converged = FALSE with a
    ## classed warning -- never as a converged SS with a non-zero residual.
    aff_ok <- TRUE
    if (!has_ssm && all(is.finite(r)) && is.finite(max_r) && max_r > tol) {
      A <- compiled$static$jacobian_fn(ss, x0, ss_params, ss)
      A_ok <- is.matrix(A) && all(is.finite(A)) && nrow(A) == ncol(A)
      rc <- if (A_ok) rcond(A) else 0
      y_aff <- NULL
      if (is.finite(rc) && rc > .Machine$double.eps) {
        y_aff <- as.numeric(solve(A, -r))
      } else if (A_ok) {
        y_aff <- as.numeric(.safe_inv(A) %*% (-r))
        if (verbose)
          .dynhr_cat(paste0("  Linear model has non-zero constants and a ",
                            "singular static Jacobian: minimum-norm ",
                            "least-squares solve.\n"))
      }
      r_aff <- if (!is.null(y_aff)) {
        y_aff <- setNames(y_aff, endo_names)
        compiled$static$residuals_fn(y_aff, x0, ss_params, y_aff)
      }
      ## Accept only if it really is a steady state (guards an inconsistent
      ## singular system, and a model that is declared linear but is not
      ## affine in the static variables).
      if (!is.null(r_aff) && all(is.finite(r_aff)) &&
          max(abs(r_aff)) <= max(tol, max_r * 1e-8)) {
        ss    <- y_aff
        r     <- r_aff
        max_r <- max(abs(r_aff))
        if (verbose)
          .dynhr_cat(sprintf(paste0("  Linear model has non-zero constants: ",
                                    "solved the affine steady state ",
                                    "A y = -c (max|r| = %.3e).\n"), max_r))
      } else {
        aff_ok <- FALSE
        if (!is.null(r_aff) && all(is.finite(r_aff)) &&
            max(abs(r_aff)) < max_r) {
          ## Report the least-squares point, the closest the system gets.
          ss    <- y_aff
          r     <- r_aff
          max_r <- max(abs(r_aff))
        }
        .dynhr_warn(sprintf(paste0(
          "solve_steady_state: the linear model's static system A y = -c has ",
          "no solution (max|residual| = %.3e at the least-squares point%s). ",
          "Returning converged = FALSE. A drift term on a unit-root variable ",
          "(e.g. `y = y(-1) + c + e`) has no steady state."),
          max_r, if (A_ok && !(is.finite(rc) && rc > .Machine$double.eps))
            "; the static Jacobian is singular" else ""),
          class = "dynhr_warning_steady_state_inconsistent")
      }
    }

    # With a steady_state_model block (e.g. SW2007) the analytical SS SHOULD
    # drive the residual to ~0, and a large residual means the obs-only
    # constants were not resolved (the I10 silent-wrong-SS case). Without one,
    # the affine solve above has either zeroed the residual (a zero-constant
    # model never enters it: its residual at 0 is already 0) or flagged the
    # system as inconsistent; `converged` is TRUE only in the former case.
    conv <- if (has_ssm) (max_r < tol * 100) else aff_ok
    if (!conv && verbose)
      .dynhr_cat(sprintf("  Linear SS: max|r| = %.3e > tol; unconverged.\n", max_r))
    result_lin <- list(ss = ss, values = ss, residuals = r,
                       converged = conv,
                       iterations = 0L, max_residual = max_r,
                       method_used = "linear")
    if (has_ssm) result_lin$params <- ss_params
    return(.as_dynhr_steady(result_lin))
  }
  
  # Initial guess
  if (is.null(y0)) {
    if (length(model$initval) > 0) {
      y0 <- setNames(rep(0.5, n), endo_names)
      for (nm in names(model$initval))
        if (nm %in% endo_names) y0[[nm]] <- model$initval[[nm]]
    } else {
      y0 <- setNames(rep(1, n), endo_names)
    }
  }
  
  # Method cascade with safety max_attempts guard.
  # The cascade runs sequentially once (analytical -> Newton -> nleqslv -> optim).
  # max_attempts protects against theoretical edge cases where method dispatch
  # could retry; each attempt tries all remaining methods in order.
  for (.attempt in seq_len(max_attempts)) {
  result <- NULL
  
  # Analytical
  if (method %in% c("auto", "analytical")) {
    if (verbose) .dynhr_cat("Trying analytical steady state...\n")
    ss_a_result <- solve_steady_state_analytical(model, params)
    if (!is.null(ss_a_result)) {
      ss_a <- ss_a_result$ss
      # Update params with any parameter values computed in the
      # steady_state_model block (e.g. gammax, delta, beta).
      params <- ss_a_result$params
    }
    if (!is.null(ss_a_result) && all(is.finite(ss_a))) {
      x0 <- .resolve_exo_init(exo_init, exo_names)
      r <- compiled$static$residuals_fn(ss_a, x0, params, ss_a)
      max_r <- max(abs(r))
      if (max_r < tol * 100) {
        if (verbose) .dynhr_cat(sprintf("Analytical SS verified: max|r| = %.3e\n", max_r))
        result <- list(ss = ss_a, values = ss_a, residuals = r,
                       converged = TRUE, iterations = 0L,
                       max_residual = max_r, method_used = "analytical",
                       params = ss_a_result$params)
      } else {
        if (verbose) .dynhr_cat(sprintf("Analytical failed verification: max|r| = %.3e\n", max_r))
        y0 <- ss_a
      }
    } else if (method == "analytical") {
      stop("No steady_state_model block or evaluation failed")
    }
  }
  
  # Newton (uses nleqslv as the numerical workhorse)
  if (is.null(result) && method %in% c("auto", "Newton")) {
    if (verbose) .dynhr_cat("Trying Newton solver (nleqslv)...\n")
    result <- solve_ss_nleqslv(compiled, params, y0 = y0,
                               endo_names = endo_names, exo_names = exo_names,
                               exo_init = exo_init,
                               tol = tol, verbose = verbose)
    if (!is.null(result) && result$converged) {
      result$method_used <- "Newton"; result$ss <- result$values
    } else {
      if (verbose) {
        .dynhr_cat("Newton did not converge (max|r| =",
            if (!is.null(result)) result$max_residual else "NULL", ")\n")
      }
      if (method == "Newton") {
        if (!is.null(result))
          return(.as_dynhr_steady(c(result, list(ss = result$values, method_used = "Newton"))))
        result <- NULL
      } else {
        result <- NULL
      }
    }
  }
  
  # nleqslv
  if (is.null(result) && method %in% c("auto", "nleqslv")) {
    if (verbose) .dynhr_cat("Trying nleqslv...\n")
    result <- solve_ss_nleqslv(compiled, params, y0 = y0,
                               endo_names = endo_names, exo_names = exo_names,
                               exo_init = exo_init,
                               tol = tol, verbose = verbose)
    if (!is.null(result) && result$converged) {
      result$method_used <- "nleqslv"; result$ss <- result$values
    } else {
      if (method == "nleqslv" && !is.null(result))
        return(.as_dynhr_steady(c(result, list(ss = result$values, method_used = "nleqslv"))))
      result <- NULL
    }
  }
  
  # optim
  if (is.null(result) && method %in% c("auto", "optim")) {
    if (verbose) .dynhr_cat("Trying optim...\n")
    result <- solve_ss_optim(compiled, params, y0 = y0,
                             endo_names = endo_names, exo_names = exo_names,
                             exo_init = exo_init,
                             max_iter = max_iter * 5, tol = tol, verbose = verbose)
    if (!is.null(result)) { result$method_used <- "optim"; result$ss <- result$values }
  }
  
  if (!is.null(result)) break
  if (verbose && .attempt < max_attempts)
    .dynhr_cat(sprintf("Attempt %d failed, retrying (%d remaining)...\n", .attempt, max_attempts - .attempt))
  } # end attempt loop
  
  if (is.null(result)) stop("All steady state methods failed after ", max_attempts, " attempts")
  .as_dynhr_steady(result)
}

.as_dynhr_steady <- function(x) {
  if (!inherits(x, "dynhr_steady")) class(x) <- c("dynhr_steady", "list")
  x
}

#' @export
print.dynhr_steady <- function(x, ...) {
  conv <- if (isTRUE(x$converged)) "CONVERGED" else "NOT CONVERGED"
  meth <- if (!is.null(x$method_used)) x$method_used else "?"
  cat(sprintf("<dynhr_steady>  [%s, method = %s]\n", conv, meth))
  if (!is.null(x$values)) {
    n <- length(x$values)
    cat(sprintf("  %d endogenous variables at SS\n", n))
    max_res <- if (!is.null(x$residuals)) max(abs(x$residuals), na.rm = TRUE) else NA
    if (is.finite(max_res)) cat(sprintf("  Max residual: %.3e\n", max_res))
  }
  invisible(x)
}

# ============================================================================
# 7. STEADY STATE VERIFICATION (NA-robust)
# ============================================================================

verify_steady_state <- function(ss, compiled, params, tol = 1e-8) {
  endo_names <- compiled$model$var_names
  exo_names  <- compiled$model$varexo_names
  x0 <- setNames(rep(0, length(exo_names)), exo_names)

  # OccBin: for verification use only the relax-regime residuals
  static_res_fn <- compiled$static$residuals_fn
  dynamic_res_fn <- compiled$dynamic$residuals_fn
  if (inherits(compiled, "dynhr_compiled") &&
      !is.null(compiled$occbin) &&
      compiled$static$n_eq > length(endo_names)) {
    relax_rows_v <- compiled$occbin$regime_map[[1L]]$eq_indices
    static_res_fn <- local({
      rfn <- compiled$static$residuals_fn; rr <- relax_rows_v
      function(y, x, params, ss) rfn(y, x, params, ss)[rr]
    })
    if (!is.null(dynamic_res_fn)) {
      dynamic_res_fn <- local({
        rfn <- compiled$dynamic$residuals_fn; rr <- relax_rows_v
        function(dy, params, ss) rfn(dy, params, ss)[rr]
      })
    }
  }

  # Parameter completeness check
  if (any(is.na(params))) {
    na_params <- names(params)[is.na(params)]
    cat("=== Steady State Verification ===\n")
    cat("  WARNING:", length(na_params), "parameter(s) have NA values:\n")
    cat("    ", paste(head(na_params, 10), collapse = ", "),
        if (length(na_params) > 10) ", ..." else "", "\n")
  } else {
    cat("=== Steady State Verification ===\n")
  }

  # Static residuals (NA-safe)
  r_static <- static_res_fn(ss, x0, params, ss)
  max_static <- max(abs(r_static), na.rm = TRUE)
  if (!is.finite(max_static)) {
    cat("  Static residuals:  max|r| = NA/Inf  [FAIL]\n")
    cat("  (likely cause: missing parameter values or bad initial guess)\n")
    na_eqs <- which(is.na(r_static))
    if (length(na_eqs) > 0 && length(na_eqs) <= 10)
      cat("    NA in equations:", paste(na_eqs, collapse = ", "), "\n")
    max_static <- Inf
  } else {
    cat(sprintf("  Static residuals:  max|r| = %.3e", max_static))
    if (max_static < tol) cat("  [OK]\n") else cat("  [FAIL]\n")
  }
  
  # Report equations with large residuals
  if (is.finite(max_static) && max_static >= tol) {
    bad_eqs <- which(abs(r_static) > tol)
    if (length(bad_eqs) > 0 && length(bad_eqs) <= 10) {
      cat("  Equations with |residual| >", tol, ":\n")
      for (idx in bad_eqs)
        cat(sprintf("    eq %d: residual = %.6e\n", idx, r_static[idx]))
    }
  }
  
  # Dynamic residuals (NA-safe)
  dy <- numeric(0)
  for (k in seq_len(nrow(compiled$dynamic$dyn_col_map))) {
    nm  <- compiled$dynamic$dyn_col_map$name[k]
    ll  <- compiled$dynamic$dyn_col_map$lead_lag[k]
    sfx <- if (ll == 0L) "__0"
    else if (ll > 0L) paste0("__p", ll)
    else paste0("__m", abs(ll))
    key <- paste0(nm, sfx)
    val <- if (nm %in% names(ss)) ss[[nm]]
    else if (nm %in% exo_names) 0
    else 0
    dy[key] <- val
  }
  r_dyn <- dynamic_res_fn(dy, params, ss)
  max_dyn <- max(abs(r_dyn), na.rm = TRUE)
  if (!is.finite(max_dyn)) {
    cat("  Dynamic residuals: max|r| = NA/Inf  [FAIL]\n")
    dyn_ok <- FALSE
  } else {
    cat(sprintf("  Dynamic residuals: max|r| = %.3e", max_dyn))
    if (max_dyn < tol) cat("  [OK]\n") else cat("  [FAIL]\n")
    dyn_ok <- max_dyn < tol
  }
  
  cat("================================\n")
  result <- is.finite(max_static) && (max_static < tol) && isTRUE(dyn_ok)
  invisible(result)
}

# ============================================================================
# 8. DISPLAY UTILITIES
# ============================================================================
