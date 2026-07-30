## R/bk-wall-transform.R
## --------------------------------------------------------------------------
## Wall-coordinate reparameterization (attempt #12 ask (b) / pathological-DSGE
## RESEARCH_AGENDA Thread 4b, PRACTICAL scope). Samplers currently see
## loglik = -Inf outside the Blanchard-Kahn determinate region, so boundary
## proposals are wasted divergences. A GLOBAL diffeomorphism from R^P onto
## the determinate region is expected PROHIBITIVE (bk_distance()'s |lambda_c|
## - 1 is non-smooth at repeated/defective crossings, R/bk-distance.R); this
## file ships the constructible PRACTICAL target instead: a smooth map for
## ONE designated boundary-relevant parameter theta_j (e.g. the Taylor
## coefficient), holding the OTHER free parameters theta_{-j} fixed at the
## sampler's current position:
##
##   theta_j = wall_j(theta_{-j}) + exp(eta_j)                          (*)
##
## where wall_j(theta_{-j}) is the determinacy wall location along
## coordinate j, found by bracketed uniroot() on bk_distance()'s signed
## distance-to-wall. This maps eta_j in R onto the determinate side
## SMOOTHLY wherever the crossing eigenvalue at the wall is SIMPLE (the
## generic case) -- see the fail-loud checks in `.bk_wall_root()` below for
## what happens when it is not.
##
## LOG-JACOBIAN IS EXACT AND TRIVIAL. Because wall_j depends only on
## theta_{-j} -- never on theta_j itself -- the full Jacobian d theta / d eta
## (ordering eta_j last) is block lower-triangular:
##   * theta_{-j} = g(eta_{-j})            (ordinary per-coordinate rule,
##                                           independent of eta_j)
##   * theta_j    = wall_j(theta_{-j}(eta_{-j})) + exp(eta_j)
## so d theta_{-j} / d eta_j = 0 identically, and
##   det(d theta / d eta) = exp(eta_j) * prod_{k != j} d theta_k / d eta_k,
##   log|J| = eta_j + sum_{k != j} log|d theta_k / d eta_k|
## with NO wall_j term at all (it cancels: wall_j enters theta_j only as an
## additive shift, so d theta_j / d eta_j = exp(eta_j) regardless of
## wall_j's value). Exactly the "trivial log-Jacobian eta_j" the design
## note promised, and verified by full finite-difference determinant in
## tests/testthat/test-bk-wall-transform.R.
##
## DOCUMENTED LIMITATION: the GRADIENT chain rule is NOT triangular. For
## k != j, d theta_j / d eta_k = (d wall_j / d theta_{-j}) . (d theta_{-j} /
## d eta_k) is generally NONZERO, but `dtheta_deta()` here -- like
## `build_param_transform()`'s -- returns only the DIAGONAL (the interface
## `make_transformed_grad()` consumes assumes a diagonal Jacobian). An
## analytic `grad_fn` composed through `make_transformed_grad()` therefore
## DROPS this cross term and is WRONG for the non-wall coordinates. This is
## harmless with `dynhr_nuts()`/`dynhr_hmc()`'s DEFAULT `grad_fn = NULL`
## path: it differentiates the wrapped eta-space target NUMERICALLY (the
## Jacobian-included `make_transformed_logpost()` output, which IS exact),
## so the leapfrog gradient is correct up to FD error regardless of the
## triangular-vs-full Jacobian question. USE `grad_fn = NULL` WITH THIS
## TRANSFORM; an analytic grad_fn is not supported (silently wrong if
## supplied). Either way the Metropolis accept/reject step -- which is what
## actually guarantees the correct stationary measure -- always uses the
## EXACT `target_fn`/`log_jacobian()`, so even a wrong gradient would only
## cost sampling efficiency, not correctness; the numerical-gradient default
## sidesteps the question rather than relying on that argument.
##
## CACHING: wall_j root-finding costs a full BK pencil build + generalized
## eigendecomposition per call, so an EXACT-match cache (`.bk_wall_root()`)
## memoizes on the theta_{-j} value -- repeat evaluations at an identical
## point (leapfrog re-evaluating the same base point, energy checks, the
## eta_j-only coordinate of a numerical gradient which does not move
## theta_{-j} at all) are free; every genuinely new theta_{-j} still pays
## the full root-find. See the overhead-factor measurement (oracle-gated) in
## the certification test for the honest per-step cost.
## --------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Internal: cheap BK distance (value only, no gradient) for uniroot()
## ---------------------------------------------------------------------------

#' Signed distance-to-wall, VALUE ONLY (no gradient / no eigenvector solve)
#'
#' Cheaper than \code{bk_distance()} for root-finding: only the crossing
#' eigenvalue's modulus is needed, so this skips the left-eigenvector solve
#' and the finite-difference gradient loop entirely.
#' @keywords internal
.bk_wall_distance_value <- function(model, compiled, params, finite_tol = 1e-9) {
  pen <- .bk_pencil(model, compiled, params, ss = NULL)
  if (is.null(pen))
    stop("bk_wall_transform(): model has no dynamic (forward/mixed) block -- ",
         "the Blanchard-Kahn pencil is empty.")
  ge  <- geigen::geigen(pen$E, pen$D, symmetric = FALSE)
  lam <- ge$alpha / ge$beta
  fin <- is.finite(lam) & (abs(ge$beta) > finite_tol * max(abs(ge$alpha), 1))
  if (!any(fin))
    stop("bk_wall_transform(): no finite generalized eigenvalues in the pencil.")
  mods <- Mod(lam[fin])
  mods[which.min(abs(mods - 1))] - 1
}

#' Count of finite generalized eigenvalues within `near_tol` of |lambda| = 1
#'
#' Used to fail loud on a NON-SIMPLE crossing: if two or more distinct
#' eigenvalues sit at (or near) the unit circle simultaneously at the wall
#' root (e.g. a complex-conjugate pair crossing together), the wall location
#' is a KINK, not a smooth root, and (*) is not a valid smooth map there.
#' @keywords internal
.bk_wall_near_wall_count <- function(D, E, near_tol = 1e-4, finite_tol = 1e-9) {
  ge  <- geigen::geigen(E, D, symmetric = FALSE)
  lam <- ge$alpha / ge$beta
  fin <- is.finite(lam) & (abs(ge$beta) > finite_tol * max(abs(ge$alpha), 1))
  if (!any(fin)) return(0L)
  mods <- Mod(lam[fin])
  sum(abs(mods - 1) < near_tol)
}

## ---------------------------------------------------------------------------
## Internal: cached wall root-finder
## ---------------------------------------------------------------------------

#' Build a memoized wall_j(\eqn{\theta_{-j}}) root-finder
#'
#' @param model,compiled Parsed model + compiled form.
#' @param wall_param Name of the boundary-relevant parameter theta_j.
#' @param other_names Names of the OTHER free parameters \eqn{\theta_{-j}} the wall
#'   depends on (order matches the vector passed to the returned function).
#' @param base_params Full named parameter vector (all model parameters,
#'   including any not being sampled) used as the template; `wall_param` and
#'   `other_names` entries are overridden per call.
#' @param bracket length-2 numeric, the theta_j search bracket for uniroot().
#' @param tol uniroot() tolerance.
#' @param cache_size Max number of memoized (\eqn{\theta_{-j}}, root) pairs (simple
#'   FIFO eviction).
#' @return function(theta_other) -> wall_j (scalar).
#' @keywords internal
.bk_wall_maker <- function(model, compiled, wall_param, other_names,
                            base_params, bracket, tol = 1e-10,
                            near_tol = 1e-4, cache_size = 200L) {
  cache_keys <- character(0)
  cache_vals <- numeric(0)

  key_of <- function(theta_other)
    paste(formatC(theta_other, digits = 15, format = "g"), collapse = "|")

  function(theta_other) {
    if (is.null(names(theta_other))) names(theta_other) <- other_names
    key <- key_of(theta_other[other_names])
    hit <- match(key, cache_keys)
    if (!is.na(hit)) return(cache_vals[[hit]])

    f <- function(tj) {
      params <- base_params
      params[other_names] <- theta_other[other_names]
      params[[wall_param]] <- tj
      .bk_wall_distance_value(model, compiled, params)
    }
    f_lo <- f(bracket[1]); f_hi <- f(bracket[2])
    if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo * f_hi > 0) {
      stop(sprintf(paste0(
        "bk_wall_transform(): bracket [%.6g, %.6g] does not straddle the BK ",
        "wall for '%s' at theta_{-j} = %s (f(lo) = %.4g, f(hi) = %.4g). ",
        "Widen `bracket`."),
        bracket[1], bracket[2], wall_param,
        paste(sprintf("%s=%.6g", other_names, theta_other[other_names]),
              collapse = ", "), f_lo, f_hi))
    }
    ur   <- stats::uniroot(f, lower = bracket[1], upper = bracket[2], tol = tol)
    root <- ur$root

    ## fail-loud: root-find actually converged to the wall
    params_at_root <- base_params
    params_at_root[other_names] <- theta_other[other_names]
    params_at_root[[wall_param]] <- root
    d_at_root <- .bk_wall_distance_value(model, compiled, params_at_root)
    if (abs(d_at_root) > 1e-6)
      stop(sprintf(paste0(
        "bk_wall_transform(): root-find for '%s' converged to distance %.3e ",
        "(expected ~0) at theta_{-j} = %s; bracket or tol may be too loose."),
        wall_param, d_at_root,
        paste(sprintf("%s=%.6g", other_names, theta_other[other_names]),
              collapse = ", ")))

    ## fail-loud: non-simple crossing (repeated / near-repeated |lambda| = 1)
    pen <- .bk_pencil(model, compiled, params_at_root, ss = NULL)
    if (is.null(pen))
      stop("bk_wall_transform(): model has no dynamic block at the wall root.")
    n_near <- .bk_wall_near_wall_count(pen$D, pen$E, near_tol = near_tol)
    if (n_near > 1L)
      stop(sprintf(paste0(
        "bk_wall_transform(): NON-SIMPLE crossing at the '%s' wall -- %d ",
        "eigenvalues sit within %.1e of |lambda| = 1 simultaneously at ",
        "theta_{-j} = %s. The wall-coordinate map (*) is not smooth here ",
        "(a repeated/complex-conjugate-pair crossing); this is a genuine ",
        "obstruction, not a bug -- narrow the prior or pick a different ",
        "`wall_param` / `other_names` combination away from this region."),
        wall_param, n_near, near_tol,
        paste(sprintf("%s=%.6g", other_names, theta_other[other_names]),
              collapse = ", ")))

    ## fail-loud: defective eigenvector overlap at the (simple) crossing --
    ## reuse bk_distance()'s own y^H D x ~ 0 guard for the analytic-gradient
    ## degeneracy case (belt-and-braces on top of the multiplicity count).
    bk_at_root <- bk_distance(model, compiled, params = params_at_root,
                               param_names = wall_param, resolve_ss = FALSE)

    if (length(cache_keys) >= cache_size) {
      cache_keys <<- cache_keys[-1]; cache_vals <<- cache_vals[-1]
    }
    cache_keys <<- c(cache_keys, key)
    cache_vals <<- c(cache_vals, root)
    root
  }
}

## ---------------------------------------------------------------------------
## bk_wall_transform()
## ---------------------------------------------------------------------------

#' Determinacy-respecting wall-coordinate reparameterization
#'
#' Builds a \code{"dynhr_param_transform"}-compatible object (same public
#' interface as \code{\link{build_param_transform}}: \code{to_unconstrained},
#' \code{to_constrained}, \code{log_jacobian}, \code{dlog_jacobian},
#' \code{dtheta_deta}) that maps the full free-parameter vector \code{theta}
#' (\code{par_names}) onto an unconstrained \code{eta} such that
#' \code{wall_param}'s constrained value is ALWAYS on the Blanchard-Kahn
#' determinate side of the wall along that coordinate:
#' \deqn{\theta_j = \mathrm{wall}_j(\\eqn{\theta_{-j}}) + \exp(\eta_j),}
#' with \eqn{\mathrm{wall}_j(\\eqn{\theta_{-j}})} the root of
#' \code{\link{bk_distance}}\code{()$distance == 0} along coordinate
#' \code{wall_param}, holding the other free parameters \code{other_names =
#' setdiff(par_names, wall_param)} fixed at their current value (bracketed
#' \code{uniroot()}, memoized -- see the file header for the caching
#' policy). The other coordinates use the ordinary
#' \code{\link{build_param_transform}} rule from \code{prior_spec}.
#'
#' Drop-in for \code{dynhr_nuts()}/\code{dynhr_hmc()}/\code{rwmh()}'s
#' \code{transform} argument, exactly like a \code{build_param_transform()}
#' result -- \strong{use it with \code{grad_fn = NULL}} (the sampler's
#' default numerical-gradient path); an analytically-supplied \code{grad_fn}
#' silently drops the wall-normal cross term for the non-wall coordinates
#' (see the file header). The log-Jacobian is exact regardless (proven
#' separable/triangular; see the file header and the FD-determinant test),
#' so the sampler always targets the correct measure even when the gradient
#' used to steer the leapfrog trajectory is only approximate.
#'
#' \strong{Scope / limits (read before use):} this is the PRACTICAL
#' single-coordinate construction, not a global diffeomorphism (expected
#' prohibitive per the pathological-DSGE RESEARCH_AGENDA -- \code{bk_distance}
#' is non-smooth at repeated/defective crossings). It is valid only where
#' the crossing eigenvalue at the wall is SIMPLE; \code{bk_wall_transform()}
#' fails loud (via the returned transform's \code{to_unconstrained}/
#' \code{to_constrained}, called lazily at the wall root, NOT at
#' construction time) if a root-find lands on a non-simple (repeated /
#' complex-conjugate-pair) crossing or the search \code{bracket} does not
#' straddle the wall. Requires at least one non-wall free parameter
#' (\eqn{\theta_{-j}}); with a single free parameter the wall location is
#' fixed and an ordinary \code{build_param_transform()} log-transform at
#' that fixed threshold suffices.
#'
#' @param model,compiled Parsed \code{dynhr_mod} and its compiled form.
#' @param prior_spec data.frame (as consumed by
#'   \code{\link{build_param_transform}}) covering `other_names`'s support.
#' @param par_names Character vector, the full free-parameter vector order
#'   (same convention as \code{build_param_transform}).
#' @param wall_param Single parameter name in \code{par_names}: the
#'   boundary-relevant coordinate theta_j.
#' @param bracket length-2 numeric, the theta_j search range for uniroot()
#'   (must straddle the wall for every \eqn{\theta_{-j}} the sampler will visit;
#'   widen if `to_unconstrained`/`to_constrained` fail loud on a bracket
#'   error).
#' @param params Optional named numeric vector of ALL model parameters
#'   (defaults to \code{model$param_values}) used as the template for
#'   parameters NOT in \code{par_names} (fixed structural parameters).
#' @param tol uniroot() tolerance (default 1e-10).
#' @param near_tol Modulus tolerance for the non-simple-crossing guard
#'   (default 1e-4).
#' @param cache_size Max memoized (\eqn{\theta_{-j}}, root) pairs (default 200).
#' @return An object of class \code{c("dynhr_bk_wall_transform",
#'   "dynhr_param_transform")} with the fields listed above plus
#'   \code{wall_param}, \code{other_names}, \code{bracket}, and
#'   \code{wall_fn(theta_other)} (the memoized root-finder itself, exposed
#'   for direct inspection/testing).
#' @seealso \code{\link{bk_distance}}, \code{\link{build_param_transform}}
#' @export
bk_wall_transform <- function(model, compiled, prior_spec, par_names,
                               wall_param, bracket, params = NULL,
                               tol = 1e-10, near_tol = 1e-4,
                               cache_size = 200L) {
  stopifnot(is.character(par_names), length(par_names) >= 2L,
            is.character(wall_param), length(wall_param) == 1L,
            wall_param %in% par_names,
            is.numeric(bracket), length(bracket) == 2L, bracket[1] < bracket[2])

  other_names <- setdiff(par_names, wall_param)
  if (!length(other_names))
    stop("bk_wall_transform(): needs at least one non-wall free parameter ",
         "(theta_{-j}) besides `wall_param` -- with a single free parameter ",
         "the wall location is fixed and build_param_transform()'s ordinary ",
         "log-transform at that threshold already does the job.")

  base_params <- if (is.null(params)) model$param_values else params
  if (!all(par_names %in% names(base_params)))
    stop("bk_wall_transform(): `par_names` not all present in `params`/",
         "model$param_values: ",
         paste(setdiff(par_names, names(base_params)), collapse = ", "))

  other_tr <- build_param_transform(prior_spec, other_names)
  wall_fn  <- .bk_wall_maker(model, compiled, wall_param, other_names,
                              base_params, bracket, tol = tol,
                              near_tol = near_tol, cache_size = cache_size)

  clamp_eps <- 1e-10

  to_unconstrained <- function(theta) {
    if (is.null(names(theta))) names(theta) <- par_names
    eta <- setNames(numeric(length(par_names)), par_names)
    theta_other <- theta[other_names]
    eta[other_names] <- other_tr$to_unconstrained(theta_other)
    w <- wall_fn(theta_other)
    d <- theta[[wall_param]] - w
    if (d <= 0) {
      scale <- max(abs(w), abs(theta[[wall_param]]), 1)
      d <- clamp_eps * scale
    }
    eta[[wall_param]] <- log(d)
    eta
  }

  to_constrained <- function(eta) {
    if (is.null(names(eta))) names(eta) <- par_names
    theta <- setNames(numeric(length(par_names)), par_names)
    theta_other <- other_tr$to_constrained(eta[other_names])
    theta[other_names] <- theta_other
    w <- wall_fn(theta_other)
    theta[[wall_param]] <- w + exp(eta[[wall_param]])
    theta
  }

  log_jacobian <- function(eta) {
    if (is.null(names(eta))) names(eta) <- par_names
    other_tr$log_jacobian(eta[other_names]) + eta[[wall_param]]
  }

  dlog_jacobian <- function(eta) {
    if (is.null(names(eta))) names(eta) <- par_names
    out <- setNames(numeric(length(par_names)), par_names)
    out[other_names] <- other_tr$dlog_jacobian(eta[other_names])
    out[[wall_param]] <- 1
    out
  }

  ## DIAGONAL ONLY -- see file header. Correct for the wall coordinate
  ## itself; drops the wall-normal cross term for the other coordinates.
  ## Not used by dynhr_nuts()/dynhr_hmc() when grad_fn = NULL (recommended).
  dtheta_deta <- function(eta) {
    if (is.null(names(eta))) names(eta) <- par_names
    out <- setNames(numeric(length(par_names)), par_names)
    out[other_names] <- other_tr$dtheta_deta(eta[other_names])
    out[[wall_param]] <- exp(eta[[wall_param]])
    out
  }

  structure(
    list(
      par_names        = par_names,
      wall_param       = wall_param,
      other_names      = other_names,
      bracket          = bracket,
      to_unconstrained = to_unconstrained,
      to_constrained   = to_constrained,
      log_jacobian     = log_jacobian,
      dlog_jacobian    = dlog_jacobian,
      dtheta_deta      = dtheta_deta,
      wall_fn          = wall_fn
    ),
    class = c("dynhr_bk_wall_transform", "dynhr_param_transform")
  )
}
