#' Safe discount-factor lookup: checks for 'beta' or 'betta' in param names.
#' Both sides of a %||% call are evaluated eagerly in R, so params[["beta"]]
#' would throw "subscript out of bounds" when the parameter is named 'betta'.
#' @param params Named numeric parameter vector.
#' @return The discount factor value, or NULL if not found.
#' @noRd
.get_discount <- function(params) {
  if (is.null(params)) return(NULL)
  if ("beta" %in% names(params)) return(as.numeric(params[["beta"]]))
  if ("betta" %in% names(params)) return(as.numeric(params[["betta"]]))
  NULL
}

#' Resolve the Ramsey planner's discount factor
#'
#' Precedence: the caller's \code{discount} (a parameter name, applied
#' symbolically, or a number); else the \code{.mod}'s
#' \code{ramsey_model(planner_discount = ...)} /
#' \code{ramsey_policy(planner_discount = ...)} option (a parameter name is
#' kept symbolic, any other expression is evaluated at \code{params}); else a
#' parameter named \code{beta} / \code{betta}.  Errors (class
#' \code{dynhr_error_no_discount}) when none of these gives a finite number --
#' there is no silent 0.99 default.
#' @param discount NULL, a parameter name, or a number.
#' @param model    dynhr_mod (its parsed commands and parameter names).
#' @param params   Named numeric parameter vector the value is taken from.
#' @param fn       Caller name for the error message.
#' @return list(arg = discount to pass on (NULL = auto-detect beta/betta
#'   symbolically; a parameter name; or a number), value = numeric).
#' @noRd
.ramsey_discount <- function(discount, model, params, fn) {
  if (is.null(discount)) {
    pd <- .mod_command_options(model, c("ramsey_model", "ramsey_policy"))$planner_discount
    if (!is.null(pd)) {
      discount <- if (is.character(pd) && length(pd) == 1L &&
                      pd %in% model$param_names && pd %in% names(params)) pd
                  else .eval_mod_scalar_option(pd, params, "planner_discount")
    }
  }
  value <- if (is.character(discount) && length(discount) == 1L) {
    if (discount %in% names(params)) as.numeric(params[[discount]])
    else suppressWarnings(as.numeric(discount))
  } else if (is.numeric(discount) && length(discount) == 1L) {
    as.numeric(discount)
  } else {
    .get_discount(params)
  }
  if (is.null(value) || length(value) != 1L || !is.finite(value))
    .dynhr_abort(
      fn, "(): no discount factor for the planner. Pass discount = <value> ",
      "(Dynare: planner_discount =), write ramsey_model(planner_discount = ",
      "...) in the .mod, or calibrate a parameter named beta / betta.",
      class = "dynhr_error_no_discount")
  list(arg = discount, value = value)
}

#' Ramsey optimal policy workflow (v0.2 core support)
#'
#' Solves the model at supplied parameters, then evaluates the planner objective
#' under a stochastic simulation to obtain a practical welfare metric.
#'
#' This initial implementation is intentionally compatible with existing dynhr
#' parse/compile/solve pipelines and object structures.
#'
#' @param model dynhr_mod
#' @param compiled Optional dynhr_compiled
#' @param params Named parameter vector; defaults to \code{model$param_values}
#' @param planner_objective Optional objective expression text. If NULL, uses
#'   parsed \code{planner_objective(...);} from the model when available.
#' @param order Perturbation order for policy solution (1 or 2)
#' @param n_periods Simulation length used for welfare evaluation.  Default
#'   raised from 400 to 2000 in 0.9.4: the unconditional welfare is a
#'   Monte-Carlo mean over a persistent path and 400 periods left it visibly
#'   noisy.  Its batch-means standard error is reported as
#'   \code{$welfare$welfare_se}; check it before reading a welfare gap.
#' @param burn_in Burn-in discarded before welfare evaluation
#' @param discount Optional discount factor.  When \code{NULL} it is taken from
#'   the \code{.mod}'s \code{ramsey_policy(planner_discount = ...)} or
#'   \code{ramsey_model(planner_discount = ...)} option, else from the
#'   parameter \code{beta} or \code{betta} (Dynare's spelling); if none exists
#'   the call ERRORS (class \code{dynhr_error_no_discount}) rather than
#'   silently assuming 0.99, as it did before 0.9.4.
#' @param planner_discount Dynare-compatible alias for \code{discount}; used when
#'   \code{discount} is \code{NULL} (NEW-W1).
#' @param verbose Print progress
#' @return Object of class \code{dynhr_ramsey_result}
#' @export
ramsey_policy <- function(model,
                          compiled = NULL,
                          params = NULL,
                          planner_objective = NULL,
                          order = 1L,
                          n_periods = 2000L,
                          burn_in = 100L,
                          discount = NULL,
                          planner_discount = NULL,
                          verbose = FALSE) {
  # NEW-W1: accept Dynare's `planner_discount=` spelling as an alias for
  # `discount=` so verbatim-translated Ramsey .mod calls work unchanged.
  if (is.null(discount) && !is.null(planner_discount)) discount <- planner_discount
  if (is.null(params)) params <- model$param_values

  ## Free-instrument models have fewer equations than endogenous variables.
  ## ramsey_policy() uses solve_perturbation() on the ORIGINAL model, which
  ## requires a square system.  For non-square models, the caller must use
  ## ramsey_model() (augmented-system approach) instead.
  n_endo <- length(model$var_names)
  n_eq   <- length(model$equations)
  if (n_eq < n_endo) {
    stop(
      "ramsey_policy() cannot solve a free-instrument model (", n_eq,
      " equations, ", n_endo, " endogenous variables).\n",
      "The instrument '", paste(setdiff(model$var_names,
        unlist(lapply(model$equations, function(e) {
          if (!is.null(e$lhs) && e$lhs$type == "variable") e$lhs$name
        }))), collapse = ", "),
      "' has no own equation, so the original system is under-determined.\n",
      "Use ramsey_model() instead, which handles this via an augmented ",
      "Lagrangian approach (see ?ramsey_model)."
    )
  }

  ## Compile at the requested order: a default max_order = 1L makes
  ## ramsey_policy(order = 2) fail in solve_perturbation (it requires order-2
  ## derivatives). Resolve `order` first so the compile matches.
  order <- as.integer(order)
  if (!order %in% c(1L, 2L)) stop("ramsey_policy currently supports order = 1 or 2.")
  if (is.null(compiled))
    compiled <- compile_model(model, verbose = FALSE,
                              max_order = if (order >= 2L) 2L else 1L)

  objective_text <- planner_objective %||% .spec_planner_objective(model)
  if (is.null(objective_text) || !nzchar(trimws(objective_text))) {
    stop("No planner objective provided. Add planner_objective(...) or pass planner_objective=.")
  }

  ss <- solve_steady_state(model, compiled, params = params, verbose = FALSE)
  if (!isTRUE(ss$converged)) stop("Ramsey workflow failed: steady-state did not converge.")
  dr <- solve_perturbation(model, compiled, ss$ss, params, order = order, verbose = FALSE)

  ## Simulate at the solved order: at order = 2 the welfare below is computed
  ## from this path, so an order-1 simulate_model() here would silently discard
  ## the second-order correction (the whole point of order = 2).
  sim <- .simulate_dr_any_order(
    dr, n_periods = as.integer(n_periods), model = model,
    burn_in = as.integer(burn_in)
  )
  sim_levels <- attr(sim, "levels")

  ## Use the shared .get_discount() lookup, which ALSO recognises
  ## the Dynare spelling `betta`.  This block previously looked only for `beta`
  ## and silently fell back to 0.99, so every `betta`-calibrated model was
  ## discounted at the wrong rate with no diagnostic.  There is now no silent
  ## default: a model with no discount parameter must say so explicitly via
  ## `discount=` / `planner_discount=`.  The .mod's own
  ## ramsey_policy(planner_discount = ...) / ramsey_model(planner_discount =
  ## ...) option comes before the beta/betta lookup (it used to be ignored).
  discount <- .ramsey_discount(discount, model, params, "ramsey_policy")$value

  objective_ast <- parse_expression(
    text = objective_text,
    var_names = model$var_names,
    param_names = model$param_names
  )

  obj_t <- apply(sim_levels, 1L, function(row) {
    vals <- setNames(as.numeric(row), colnames(sim_levels))
    .eval_planner_ast(objective_ast, vals, params, ss$ss)
  })

  obj_mean <- mean(obj_t, na.rm = TRUE)
  welfare_unconditional <- obj_mean / max(1e-8, (1 - discount))
  welfare_steady <- .eval_planner_ast(objective_ast, ss$ss, params, ss$ss) /
    max(1e-8, (1 - discount))

  ## The unconditional welfare is a Monte-Carlo average over a
  ## SERIALLY CORRELATED simulated path, so report its uncertainty instead of
  ## letting callers read the last digits as signal.  Non-overlapping batch
  ## means (batch length ~ sqrt(n)) is the standard consistent estimator of the
  ## asymptotic variance of a mean from a stationary correlated series; it
  ## needs no bandwidth choice and handles the persistence of DSGE paths.
  obj_ok  <- obj_t[is.finite(obj_t)]
  n_ok    <- length(obj_ok)
  obj_se  <- NA_real_
  n_batch <- NA_integer_
  if (n_ok >= 16L) {
    b       <- max(2L, as.integer(floor(sqrt(n_ok))))
    n_batch <- as.integer(floor(n_ok / b))
    bm      <- colMeans(matrix(obj_ok[seq_len(n_batch * b)], nrow = b))
    ## Var(mean) = Var(batch means) / n_batch (batches ~ independent for b >> tau)
    obj_se  <- sqrt(stats::var(bm) / n_batch)
  }
  welfare_se <- obj_se / max(1e-8, (1 - discount))

  out <- list(
    objective = list(text = objective_text, ast = objective_ast),
    steady = ss,
    dr = dr,
    welfare = list(
      objective_mean = obj_mean,
      objective_mean_se = obj_se,
      steady_value = welfare_steady,
      unconditional_value = welfare_unconditional,
      unconditional_se = welfare_se,
      welfare_se = welfare_se,
      gap_vs_steady = welfare_unconditional - welfare_steady,
      gap_se = welfare_se,           # welfare_steady is deterministic
      n_batches = n_batch,
      discount = discount
    ),
    meta = list(order = order, n_periods = n_periods, burn_in = burn_in)
  )
  class(out) <- c("dynhr_ramsey_result", "list")

  if (isTRUE(verbose)) {
    .dynhr_inform(sprintf(
      "[dynhr] Ramsey welfare evaluated (unconditional=%.6f, steady=%.6f)",
      out$welfare$unconditional_value, out$welfare$steady_value
    ))
  }
  out
}

#' Evaluate planner objective from a Ramsey run
#'
#' Accepts either a \code{dynhr_ramsey_result} (from \code{ramsey_policy()})
#' or a \code{dynhr_ramsey_result2} (from \code{ramsey_model()}).
#'
#' @param x \code{dynhr_ramsey_result} or \code{dynhr_ramsey_result2}
#' @return Numeric scalar unconditional welfare value
#' @export
evaluate_planner_objective <- function(x) {
  if (inherits(x, "dynhr_ramsey_result")) {
    # ramsey_policy() result: welfare$unconditional_value is direct
    return(x$welfare$unconditional_value)
  }
  if (inherits(x, "dynhr_ramsey_result2")) {
    # ramsey_model() result: welfare$unconditional is direct; fall back to
    # steady_value when unconditional wasn't computed (no simulation run)
    val <- x$welfare$unconditional %||% x$welfare$steady_value
    if (!is.null(val) && is.finite(val)) return(val)
    return(x$welfare$steady_value)
  }
  stop("evaluate_planner_objective expects a dynhr_ramsey_result or dynhr_ramsey_result2.")
}

#' @export
print.dynhr_ramsey_result <- function(x, ...) {
  cat("\n<dynhr_ramsey_result>\n")
  cat(sprintf("  Order                 : %d\n", x$meta$order))
  cat(sprintf("  Objective mean (sim)  : %.6f\n", x$welfare$objective_mean))
  cat(sprintf("  Welfare (steady)      : %.6f\n", x$welfare$steady_value))
  cat(sprintf("  Welfare (unconditional): %.6f\n", x$welfare$unconditional_value))
  ## The unconditional value is a simulation mean: show its batch-means MC SE
  ## so a "welfare gap" is not read past the noise floor.
  if (!is.null(x$welfare$welfare_se) && is.finite(x$welfare$welfare_se))
    cat(sprintf("  MC SE (batch means)   : %.6f  (%d batches, %d periods)\n",
                x$welfare$welfare_se,
                as.integer(x$welfare$n_batches %||% NA_integer_),
                as.integer(x$meta$n_periods)))
  cat(sprintf("  Gap vs steady         : %.6f\n", x$welfare$gap_vs_steady))
  invisible(x)
}

.eval_planner_ast <- function(node, vars, params, ss) {
  if (is.null(node)) return(NA_real_)
  switch(
    node$type,
    number = node$value,
    variable = {
      nm <- node$name
      if (!is.null(vars[[nm]])) vars[[nm]] else NA_real_
    },
    parameter = {
      nm <- node$name
      if (!is.null(params[[nm]])) params[[nm]] else NA_real_
    },
    local_variable = {
      NA_real_
    },
    unaryop = {
      val <- .eval_planner_ast(node$operand, vars, params, ss)
      if (node$op == "-") -val else val
    },
    binop = {
      l <- .eval_planner_ast(node$left, vars, params, ss)
      r <- .eval_planner_ast(node$right, vars, params, ss)
      switch(
        node$op,
        "+" = l + r,
        "-" = l - r,
        "*" = l * r,
        "/" = l / r,
        "^" = l ^ r,
        stop(".eval_planner_ast: unsupported binary operator \"", node$op,
             "\" in planner objective.", call. = FALSE)
      )
    },
    funcall = {
      fname <- tolower(node$name)
      args <- lapply(node$args, .eval_planner_ast, vars = vars, params = params, ss = ss)
      if (fname == "steady_state" && length(node$args) == 1L && node$args[[1]]$type == "variable") {
        return(ss[[node$args[[1]]$name]] %||% NA_real_)
      }
      do.call(
        switch(
          fname,
          ln = log,
          log = log,
          exp = exp,
          sqrt = sqrt,
          abs = abs,
          min = min,
          max = max,
          sin = sin,
          cos = cos,
          tan = tan,
          ## Unknown function name: fail loud. Silently mapping it to pnorm (the
          ## old default) computed a completely wrong planner objective.
          stop(".eval_planner_ast: unsupported function \"", node$name,
               "\" in planner objective. Supported: ln/log, exp, sqrt, abs, ",
               "min, max, sin, cos, tan.", call. = FALSE)
        ),
        args
      )
    },
    NA_real_
  )
}


# ============================================================================
# Phase B: Augmented-system Ramsey engine
# ============================================================================

#' Ramsey optimal policy via augmented-system approach
#'
#' The main entry point for Ramsey optimal policy in dynhr. This function
#' implements the augmented-system approach (Bodenstein & Guerrieri 2019):
#'   1. Generates an augmented .mod file with Lagrange multipliers and FOCs
#'   2. Parses and compiles the augmented model
#'   3. Solves the competitive-equilibrium steady state
#'   4. Solves the augmented (Ramsey) steady state
#'   5. Runs perturbation on the augmented system
#'   6. Returns decision rules for the Ramsey-optimal policy
#'
#' @param model             A dynhr_mod object.
#' @param compiled          Optional compiled original model.
#' @param params            Named parameter vector. Defaults to model$param_values.
#' @param planner_objective Optional objective expression. If NULL, uses the
#'   parsed \code{planner_objective(...)} from the model.
#' @param order             Perturbation order: 1 (default), 2, or 3.
#' @param method            Solution method: \code{"augmented"} (default, Phase B,
#'   augmented-system approach) or \code{"nn1"} (Phase E, Gross-Hansen (n,n+1)
#'   approximation). When \code{method = "nn1"}, delegates to
#'   \code{\link{ramsey_nn1}}.
#' @param prefix            Prefix for Lagrange multiplier variable names (default "MULT").
#' @param discount          Discount factor for the Lagrangian used to derive
#'   the FOCs.  May be a character string naming a model parameter (applied
#'   symbolically, e.g. \code{"beta"}), a numeric value, or \code{NULL}
#'   (default) to take the \code{.mod}'s \code{ramsey_model(planner_discount =
#'   ...)} (or \code{ramsey_policy(planner_discount = ...)}) option, else to
#'   auto-detect a parameter named \code{"beta"} or \code{"betta"}.  If none
#'   of these exists the call errors (class \code{dynhr_error_no_discount});
#'   there is no silent 0.99 default.
#'   This must match the model's structural discount factor; the augmented
#'   decision rules are wrong otherwise (see \code{\link{ramsey_augment_mod}}).
#' @param planner_discount  Dynare-compatible alias for \code{discount}; used when
#'   \code{discount} is \code{NULL} (NEW-W1).
#' @param solve_orig_ss     If TRUE (default), solve the competitive-equilibrium SS.
#'   Set to FALSE to supply a pre-computed SS via \code{orig_ss}.
#' @param orig_ss           Pre-computed competitive-equilibrium SS (optional).
#' @param refine_ss         If TRUE (default), run Newton refinement on augmented SS.
#' @param verbose           Print progress messages.
#' @param ...               Additional arguments passed to \code{solve_perturbation()}
#'   or \code{\link{ramsey_nn1}} (when \code{method = "nn1"}).
#'
#' @return An object of class \code{dynhr_ramsey_result2} containing:
#'   \item{augmented_model}{Augmented dynhr_mod object.}
#'   \item{augmented_text}{Augmented .mod file text.}
#'   \item{ramsey_steady}{Augmented steady state (class \code{dynhr_ramsey_steady}).}
#'   \item{ramsey_dr}{Ramsey-optimal decision rules (class \code{dynhr_ramsey_dr}).}
#'   \item{competitive_ss}{Competitive-equilibrium steady state.}
#'   \item{multiplier_map}{Mapping of multiplier names to original equations.}
#'   \item{welfare_steady}{Steady-state welfare under Ramsey policy.}
#'   \item{meta}{Metadata (order, variable counts, timing).}
#' @references
#'   Bodenstein, M., & Guerrieri, L. (2011). The welfare-optimal degree of
#'     central bank transparency. \emph{Journal of Money, Credit and Banking},
#'     43(5), 857-889.
#'   Schmitt-Grohé, S., & Uribe, M. (2004). Solving dynamic general equilibrium
#'     models using a second-order approximation to the policy function.
#'     \emph{Journal of Economic Dynamics and Control}, 28(4), 755-775.
#' @export
ramsey_model <- function(model,
                          compiled = NULL,
                          params = NULL,
                          planner_objective = NULL,
                          order = 1L,
                          method = c("augmented", "nn1"),
                          prefix = "MULT",
                          discount = NULL,
                          planner_discount = NULL,
                          solve_orig_ss = TRUE,
                          orig_ss = NULL,
                          refine_ss = TRUE,
                          verbose = FALSE,
                          ...) {
  # ---- 1. Validate ----
  # NEW-W1: accept Dynare's `planner_discount=` as an alias for `discount=`.
  if (is.null(discount) && !is.null(planner_discount)) discount <- planner_discount
  if (!inherits(model, "dynhr_mod")) {
    stop("model must be a dynhr_mod object.")
  }
  order <- as.integer(order)
  if (!order %in% c(1L, 2L, 3L)) {
    stop("order must be 1, 2, or 3.")
  }
  if (is.null(params)) params <- model$param_values

  # ---- 1b. Fill in steady_state_model computed parameters ----
  # Parameters like kappa, Omega, lambda, vartheta may be declared as parameters
  # but computed in the steady_state_model block rather than calibrated directly.
  # Evaluate the SS block to fill these in so downstream steps (augmented model
  # parsing, ramsey_steady, ast_eval) can resolve them.
  ssm <- model$steady_state_model
  if (length(ssm) > 0) {
    ## A-SEC: allowlist-sandboxed env + the checked parse cache; an unsafe
    ## statement aborts (re-raised past the NA fallback), never evaluates.
    env <- .ssm_eval_env(params)
    for (assignment in ssm) {
      val <- tryCatch(eval(.cached_parse(assignment$text), envir = env),
                      error = function(e) .dynhr_reraise_unsafe(e, NA))
      if (is.numeric(val) && length(val) == 1 && is.finite(val)) {
        assign(assignment$name, val, envir = env)
        if (!assignment$name %in% names(params)) {
          params[[assignment$name]] <- val
        }
      }
    }
  }

  # ---- 2. Method dispatch ----
  method <- match.arg(method)

  # Resolve the discount: the user's discount= (a parameter name, applied
  # symbolically, or a number), else the .mod's
  # ramsey_model(planner_discount = ...) option (used to be ignored in favour
  # of beta), else beta/betta; no silent 0.99 when none of them exists.
  disc           <- .ramsey_discount(discount, model, params, "ramsey_model")
  discount       <- disc$arg
  discount_value <- disc$value

  if (method == "nn1") {
    # Delegate to the (n,n+1) approximation (Phase E)
    if (verbose) .dynhr_cat("[ramsey_model] Delegating to ramsey_nn1 (method='nn1')...\n")

    # The 'n' parameter for (n,n+1) approximation: (order, order+1)
    nn1_result <- ramsey_nn1(
      model            = model,
      planner_objective = planner_objective,
      n                = order,
      ramsey_result    = NULL,
      return_model     = FALSE,
      beta             = discount_value,
      orig_ss          = orig_ss,
      compiled         = compiled,
      verbose          = verbose,
      ...
    )

    # Map nn1 result to a compatible dynhr_ramsey_result2-like structure
    result <- list(
      augmented_model  = nn1_result$modified_model %||% model,
      augmented_text   = NULL,
      augmented_result = NULL,
      ramsey_steady    = list(ss = nn1_result$ss, converged = TRUE),
      ramsey_dr        = list(ramsey_dr = nn1_result$dr, bk_ok = nn1_result$bk_ok),
      competitive_ss   = nn1_result$ss,
      multiplier_map   = list(),
      welfare          = list(
        steady_value = nn1_result$welfare$steady_state %||% NA_real_,
        unconditional = nn1_result$welfare$unconditional %||% NA_real_,
        ## Why the unconditional welfare is NA, when it is
        ## ("order1_certainty_equivalent" at n = 1) -- D18 reports it.
        unconditional_reason =
          nn1_result$welfare$unconditional_reason %||% NA_character_,
        discount     = nn1_result$welfare$discount %||% discount_value
      ),
      meta = list(
        method           = "nn1",
        n                = nn1_result$n,
        order            = order,
        n_orig_vars      = length(model$var_names),
        n_multipliers    = length(nn1_result$multipliers$lambda),
        n_total_aug_vars = length(model$var_names),
        bk_ok            = nn1_result$bk_ok,
        prefix           = prefix,
        timestamp        = Sys.time()
      )
    )
    class(result) <- c("dynhr_ramsey_result2", "list")
    return(result)
  }

  # ---- 3. Get planner objective ----
  obj_text <- planner_objective %||% model$planner_objective$text %||% NULL
  if (is.null(obj_text) || !nzchar(trimws(obj_text))) {
    stop("No planner objective provided. Add planner_objective(...) or pass planner_objective=.")
  }

  # ---- 4. Generate augmented .mod ----
  if (verbose) .dynhr_cat("[ramsey_model] Generating augmented .mod...\n")
  aug_result <- ramsey_augment_mod(
    model,
    planner_objective = obj_text,
    prefix = prefix,
    discount = discount,
    verbose = verbose
  )

  # ---- 4. Parse augmented model ----
  if (verbose) .dynhr_cat("[ramsey_model] Parsing augmented model...\n")
  aug_model <- ramsey_parse_augmented(aug_result, verbose = verbose)

  # ---- 5. Solve competitive-equilibrium steady state ----
  if (is.null(orig_ss) && solve_orig_ss) {
    if (verbose) .dynhr_cat("[ramsey_model] Solving competitive-equilibrium SS...\n")
    if (is.null(compiled)) {
      compiled <- compile_model(model, verbose = verbose)
    }
    # Ramsey models have fewer equations than endogenous vars (an instrument
    # is free). The standard solve_steady_state expects a square system, so
    # handle non-square models by using initval / zeros as the CE steady state.
    n_endo <- length(model$var_names)
    n_eq   <- length(model$equations)
    if (n_eq < n_endo) {
      if (verbose) {
        .dynhr_cat(sprintf("  Non-square model: %d eqs for %d vars (instrument model).\n",
                    n_eq, n_endo))
        .dynhr_cat("  Using initval (or zeros) as competitive-equilibrium SS.\n")
      }
      competitive_ss <- setNames(rep(0, n_endo), model$var_names)
      if (length(model$initval) > 0) {
        for (nm in names(model$initval)) {
          if (nm %in% model$var_names) {
            competitive_ss[nm] <- as.numeric(model$initval[[nm]])
          }
        }
      }
    } else {
      ss_result <- solve_steady_state(model, compiled, params = params,
                                       verbose = verbose)
      if (!isTRUE(ss_result$converged)) {
        stop("Competitive-equilibrium steady state did not converge.")
      }
      competitive_ss <- ss_result$values
    }
  } else if (!is.null(orig_ss)) {
    competitive_ss <- orig_ss
  } else {
    stop("Cannot solve SS: set solve_orig_ss=TRUE or supply orig_ss.")
  }

  # ---- 6. Compile augmented model ----
  ## Compile once at the order ramsey_solve() will need, then reuse for both the
  ## SS solve and the perturbation (L10 perf: avoids a redundant recompile).
  if (verbose) .dynhr_cat("[ramsey_model] Compiling augmented model...\n")
  aug_compiled <- compile_model(aug_model,
                                max_order = if (order >= 2L) 2L else 1L,
                                verbose = verbose)

  # ---- 7. Solve augmented steady state ----
  if (verbose) .dynhr_cat("[ramsey_model] Solving augmented (Ramsey) SS...\n")
  ramsey_ss <- ramsey_steady(
    model = model,
    compiled = compiled,
    aug_model = aug_model,
    aug_compiled = aug_compiled,
    orig_ss = competitive_ss,
    params = params,
    planner_objective = obj_text,
    prefix = prefix,
    refine = refine_ss,
    verbose = verbose
  )

  # ---- 8. Solve perturbation on augmented system ----
  if (verbose) .dynhr_cat("[ramsey_model] Solving augmented perturbation...\n")
  ramsey_dr <- ramsey_solve(
    aug_model = aug_model,
    aug_steady = ramsey_ss,
    params = params,
    order = order,
    verbose = verbose,
    aug_compiled = aug_compiled,
    ...
  )

  # ---- 9. Compute steady-state welfare ----
  # Evaluate planner objective at competitive-equilibrium steady state
  all_vars <- c(model$var_names, model$varexo_names)
  obj_ast <- parse_expression(obj_text,
    var_names = all_vars,
    param_names = model$param_names)

  ss_vals <- competitive_ss
  welfare_val <- .eval_planner_ast(obj_ast, ss_vals, params, ss_vals)
  discount <- discount_value
  welfare_steady <- welfare_val / max(1e-8, (1 - discount))

  # ---- 10. Build result ----
  result <- list(
    augmented_model  = aug_model,
    augmented_text   = aug_result$augmented_text,
    augmented_result = aug_result,
    ramsey_steady    = ramsey_ss,
    ramsey_dr        = ramsey_dr,
    competitive_ss   = competitive_ss,
    multiplier_map   = aug_result$multiplier_map,
    welfare          = list(
      steady_value = welfare_steady,
      discount     = discount_value
    ),
    meta = list(
      order            = order,
      n_orig_vars      = length(model$var_names),
      n_multipliers    = length(aug_result$multiplier_names),
      n_total_aug_vars = length(aug_model$var_names),
      bk_ok            = ramsey_dr$bk_ok,
      prefix           = prefix,
      timestamp        = Sys.time()
    )
  )
  class(result) <- c("dynhr_ramsey_result2", "list")
  result
}


#' @export
print.dynhr_ramsey_result2 <- function(x, ...) {
  cat(sprintf("\n<dynhr_ramsey_result2>  [Ramsey optimal policy]\n"))
  cat(sprintf("  Order:             %d\n", x$meta$order))
  cat(sprintf("  Original vars:     %d\n", x$meta$n_orig_vars))
  cat(sprintf("  Multipliers:       %d\n", x$meta$n_multipliers))
  cat(sprintf("  Augmented vars:    %d\n", x$meta$n_total_aug_vars))
  cat(sprintf("  BK condition:      %s\n",
              if (isTRUE(x$meta$bk_ok)) "OK" else if (isFALSE(x$meta$bk_ok)) "FAIL" else "N/A"))
  cat(sprintf("  Welfare (steady):  %.6f\n", x$welfare$steady_value))
  cat(sprintf("  Discount factor:   %.4f\n", x$welfare$discount))

  if (!is.null(x$ramsey_dr$ramsey_dr$ghx)) {
    cat(sprintf("\n  Ramsey ghx:  %d x %d\n",
                nrow(x$ramsey_dr$ramsey_dr$ghx),
                ncol(x$ramsey_dr$ramsey_dr$ghx)))
  }
  if (!is.null(x$ramsey_dr$ramsey_dr$ghu)) {
    cat(sprintf("  Ramsey ghu:  %d x %d\n",
                nrow(x$ramsey_dr$ramsey_dr$ghu),
                ncol(x$ramsey_dr$ramsey_dr$ghu)))
  }
  invisible(x)
}
