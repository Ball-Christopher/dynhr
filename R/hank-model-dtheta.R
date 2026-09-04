## R/hank-model-dtheta.R
## --------------------------------------------------------------------------
## Parameter derivatives propagated through the sequence-space DAG.
##
## WHY.  hank_loglik_ar_structural_grad() factors the structural score as
## dl/dtheta_k = sum_z <dl/dTheta_z, dTheta_z/dtheta_k>. The first factor is
## exact and cheap (the autocovariance-kernel adjoint). The second was only
## available by central-differencing the whole model builder -- two model
## rebuilds per parameter -- and that is what dominates: measured end to end,
## the FD route is worth ~1.1x over plain finite differences of the likelihood,
## while an EXACT dTheta_z/dtheta_k gives 5.2x-7.7x (briefs/21 7). So the
## exact route is where the speedup lives, and this file is how a caller gets
## it WITHOUT hand-constructing T_h*n x T_h*n dH_U/dtheta matrices (a packing
## error there yields a plausible-but-wrong gradient, which is exactly the
## failure mode the package must not invite).
##
## HOW.  hank_model() accumulates G[[var]][[source]] = d(var)/d(source) along a
## topological order by the chain rule, G[o][s] = sum_i J_b[o][i] G[i][s].
## Differentiating in a parameter is THE SAME LOOP with one extra term,
##
##   dG[o][s] = sum_i ( dJ_b[o][i] G[i][s] + J_b[o][i] dG[i][s] ),
##
## seeded by dG[s][s] = 0 (a source's identity block carries no parameter).
## A block therefore only has to know the derivative of ITS OWN Jacobian --
## small, local, and verifiable block by block -- which is what
## hank_simple_block(djac_dtheta = ) declares. Packing the target rows of dG
## gives dH_U/dtheta and dH_Z/dtheta; differentiating the GE solve
## H_U dU + H_Z dZ = 0 gives
##
##   d(dU)/dtheta = -H_U^{-1} ( dH_U dU + dH_Z dZ ),
##
## one solve against the ALREADY CACHED LU (.hank_ge_factor), after which
## dTheta_z/dtheta_k is a propagation identical in shape to hank_model_irf().
##
## LIMIT, and it is the reason verification exists.  This machinery is exact
## for whatever the blocks declare, and BLIND to anything they do not: a
## parameter that also enters a het block (through the steady state and the
## fake-news Jacobian, which has no analytic parameter derivative here) will
## produce a confidently wrong derivative. A parameter declared by NO block
## fails loudly below; a parameter declared by some blocks but silently active
## in others cannot be detected here, which is why
## hank_loglik_ar_structural_grad() verifies a supplied dtheta_fn against one
## central difference per parameter by default.
## --------------------------------------------------------------------------


## Derivative of a block's own Jacobian w.r.t. one named parameter, or NULL
## when the block declares none (treated as an exact zero).
.hank_block_djac <- function(blk, ss, T_h, param) {
  dj <- blk$djac_dtheta
  if (is.null(dj) || is.null(dj[[param]])) return(NULL)
  out <- dj[[param]](ss, T_h)
  if (!is.list(out))
    stop("hank_model_dtheta(): block '", blk$name, "' djac_dtheta[[\"", param,
         "\"]] must return a nested list [[output]][[input]] of T_h x T_h ",
         "matrices, as `jac` does.")
  out
}


## Pack a G-shaped object's target rows into the (n_targets*T_h) x
## (n_cols*T_h) matrix layout hank_model() uses for H_U / H_Z.
.hank_pack_target_blocks <- function(G, targets, cols, T_h) {
  M <- matrix(0, T_h * length(targets), T_h * length(cols))
  for (ti in seq_along(targets)) for (cj in seq_along(cols)) {
    blk <- G[[targets[ti]]][[cols[cj]]]
    if (is.null(blk)) next
    M[((ti - 1) * T_h + 1):(ti * T_h),
      ((cj - 1) * T_h + 1):(cj * T_h)] <- blk
  }
  M
}


#' Propagate a parameter derivative through a sequence-space DAG
#'
#' Differentiates \code{\link{hank_model}}'s own chain-rule accumulation with
#' respect to one named structural parameter, using the per-block derivatives
#' declared via \code{\link{hank_simple_block}}'s \code{djac_dtheta}. Returns
#' the derivative of the accumulated Jacobians and of the packed GE matrices.
#'
#' @param model A \code{\link{hank_model}}.
#' @param param Character: the parameter name, as used in the blocks'
#'   \code{djac_dtheta} lists.
#'
#' @return A list with \code{param}, \code{dG} (same shape as
#'   \code{model$G}; \code{NULL} entries mean an exact zero), and
#'   \code{dH_U}, \code{dH_Z}.
#' @section Blindness: the result is exact for what the blocks declare and
#'   silently zero for what they do not. A parameter no block declares is an
#'   error; a parameter that ALSO acts through a channel with no declared
#'   derivative (typically a heterogeneous-agent block, whose fake-news
#'   Jacobian has no analytic parameter derivative here) cannot be detected,
#'   and will yield a confidently wrong result. Verify against a finite
#'   difference of the model builder -- which
#'   \code{\link{hank_loglik_ar_structural_grad}} does by default.
#' @seealso \code{\link{hank_dtheta_theta_list}},
#'   \code{\link{hank_loglik_ar_structural_grad}}, \code{\link{hank_model}}
#' @export
hank_model_dtheta <- function(model, param) {
  if (!inherits(model, "hank_model"))
    stop("hank_model_dtheta(): `model` must be a hank_model().")
  if (!is.character(param) || length(param) != 1L)
    stop("hank_model_dtheta(): `param` must be a single parameter name.")
  T_h <- model$T_h
  sources <- c(model$unknowns, model$exogenous)
  ss <- model$ss
  G <- model$G

  dG <- list()
  for (s in sources) dG[[s]] <- list()      # d(I)/dtheta = 0, i.e. all NULL
  declared <- FALSE

  for (b in model$block_order) {
    blk <- model$blocks[[b]]
    Jb  <- .hank_block_jacobian(blk, ss, T_h)          # cached
    dJb <- .hank_block_djac(blk, ss, T_h, param)
    if (!is.null(dJb)) declared <- TRUE
    for (o in blk$outputs) {
      Go <- stats::setNames(vector("list", length(sources)), sources)
      for (s in sources) {
        acc <- NULL
        for (i in blk$inputs) {
          Gis <- G[[i]][[s]]
          if (!is.null(dJb) && !is.null(Gis)) {
            dJ <- dJb[[o]][[i]]
            if (!is.null(dJ)) {
              contrib <- dJ %*% Gis
              acc <- if (is.null(acc)) contrib else acc + contrib
            }
          }
          dGis <- dG[[i]][[s]]
          if (!is.null(dGis)) {
            contrib <- Jb[[o]][[i]] %*% dGis
            acc <- if (is.null(acc)) contrib else acc + contrib
          }
        }
        if (!is.null(acc)) Go[[s]] <- acc
      }
      dG[[o]] <- Go
    }
  }

  if (!declared)
    stop("hank_model_dtheta(): no block declares a derivative for parameter '",
         param, "', so the propagated derivative would be identically zero -- ",
         "a silently wrong gradient. Declare it via ",
         "hank_simple_block(djac_dtheta = list(", param, " = function(ss, ",
         "T_h) ...)), or use the finite-difference route ",
         "(hank_loglik_ar_structural_grad with dtheta_fn = NULL).")

  list(param = param, dG = dG,
       dH_U = .hank_pack_target_blocks(dG, model$targets, model$unknowns, T_h),
       dH_Z = .hank_pack_target_blocks(dG, model$targets, model$exogenous, T_h))
}


#' Parameter derivative of a sequence-space impulse response
#'
#' Differentiates \code{\link{hank_model_irf}} with respect to a structural
#' parameter, holding the driving paths \code{dZ} fixed. The GE channel
#' contributes \eqn{d(dU)/d\theta = -H_U^{-1}(dH_U\,dU + dH_Z\,dZ)}, one solve
#' against the factorization \code{\link{hank_model_irf}} already cached.
#'
#' @param model A \code{\link{hank_model}}.
#' @param dZ Named list of exogenous deviation paths (as for
#'   \code{\link{hank_model_irf}}); parameter-independent.
#' @param dtheta Either a parameter name or a \code{\link{hank_model_dtheta}}
#'   result (pass the latter to reuse one propagation across several
#'   \code{dZ}).
#'
#' @return A named list of length-\code{T_h} derivative paths, one per model
#'   variable.
#' @seealso \code{\link{hank_model_irf}}, \code{\link{hank_model_dtheta}}
#' @export
hank_model_dtheta_irf <- function(model, dZ, dtheta) {
  T_h <- model$T_h
  dt <- if (is.character(dtheta)) hank_model_dtheta(model, dtheta) else dtheta
  if (!is.list(dt) || is.null(dt$dH_U))
    stop("hank_model_dtheta_irf(): `dtheta` must be a parameter name or a ",
         "hank_model_dtheta() result.")

  ## level paths: dsrc holds dU (the GE solution) and dZ itself
  dsrc <- .hank_irf_dsrc(model, dZ)
  stack <- function(nms) do.call(c, lapply(nms, function(s) dsrc[[s]]))
  z_stack <- stack(model$exogenous)
  u_stack <- stack(model$unknowns)

  ## differentiate H_U dU + H_Z dZ = 0 at fixed dZ
  fac <- .hank_ge_factor(model$H_U)
  rhs <- dt$dH_U %*% u_stack + dt$dH_Z %*% z_stack
  dU_d <- as.numeric(-.hank_ge_solve(fac, model$H_U, rhs))

  ddsrc <- stats::setNames(vector("list", length(dsrc)), names(dsrc))
  for (s in names(dsrc)) ddsrc[[s]] <- rep(0, T_h)
  for (k in seq_along(model$unknowns))
    ddsrc[[model$unknowns[k]]] <- dU_d[((k - 1) * T_h + 1):(k * T_h)]

  sources <- c(model$unknowns, model$exogenous)
  out <- list()
  for (v in names(model$G)) {
    dv <- rep(0, T_h)
    for (s in sources) {
      Gs <- model$G[[v]][[s]]
      if (!is.null(Gs)) dv <- dv + as.numeric(Gs %*% ddsrc[[s]])
      dGs <- dt$dG[[v]][[s]]
      if (!is.null(dGs)) dv <- dv + as.numeric(dGs %*% dsrc[[s]])
    }
    out[[v]] <- dv
  }
  out
}


#' Parameter derivatives of the per-shock MA coefficients
#'
#' The \code{dTheta_z/dtheta_k} that
#' \code{\link{hank_loglik_ar_structural_grad}} contracts against
#' \code{dl/dTheta_z}: one \code{\link{hank_model_dtheta}} propagation shared
#' by every shock, then one GE solve per shock.
#'
#' @param model A \code{\link{hank_model}}.
#' @param shock_specs Named list, one entry per exogenous shock, each with a
#'   \code{rho} element (the driving path is \code{rho^t}, as in
#'   \code{\link{hank_state_space}}).
#' @param obs_vars Character vector of observable names, in the column order
#'   of the \code{Theta} matrices.
#' @param param Character parameter name, or a \code{\link{hank_model_dtheta}}
#'   result.
#'
#' @return A named list of \code{T_h x n_obs} matrices, one per shock.
#' @seealso \code{\link{hank_dtheta_fn}},
#'   \code{\link{hank_loglik_ar_structural_grad}}
#' @export
hank_dtheta_theta_list <- function(model, shock_specs, obs_vars, param) {
  dt <- if (is.character(param)) hank_model_dtheta(model, param) else param
  T_h <- model$T_h
  exo <- model$exogenous
  stats::setNames(lapply(exo, function(z) {
    rho_z <- shock_specs[[z]]$rho
    if (is.null(rho_z))
      stop("hank_dtheta_theta_list(): shock_specs[['", z, "']] needs a `rho`.")
    dZ <- stats::setNames(list(rho_z^(seq_len(T_h) - 1L)), z)
    d <- hank_model_dtheta_irf(model, dZ, dt)
    matrix(vapply(obs_vars, function(o) d[[o]], numeric(T_h)),
           T_h, length(obs_vars), dimnames = list(NULL, obs_vars))
  }), exo)
}


#' Build an exact \code{dtheta_fn} for the structural score
#'
#' Wraps \code{\link{hank_dtheta_theta_list}} into the
#' \code{function(theta, param)} callback
#' \code{\link{hank_loglik_ar_structural_grad}} expects, so the exact
#' (declared-derivative) route needs one line at the call site. The model is
#' rebuilt at most once per distinct \code{theta} and each parameter's DAG
#' propagation is memoized, so a k-parameter gradient costs one rebuild plus k
#' propagations -- not the 2k rebuilds the finite-difference route pays.
#'
#' Seed the memo with \code{model} (plus the \code{theta} it was built at) when
#' the caller already holds the model -- a sampler or a block-coordinate
#' optimizer always does -- and the first gradient at that \code{theta} costs
#' NO rebuild at all. Combined with
#' \code{\link{hank_loglik_ar_structural_grad}}'s own \code{model}/
#' \code{Theta_list} arguments this takes a structural gradient from two
#' rebuilds to zero, which is the dominant term at scale.
#'
#' @param model_fn \code{function(theta)} returning a \code{\link{hank_model}}.
#' @param shock_specs,obs_vars As for \code{\link{hank_dtheta_theta_list}}.
#' @param model Optional prebuilt \code{\link{hank_model}} to seed the memo
#'   with; requires \code{theta}, and is used only while the requested
#'   \code{theta} is \code{identical()} to it (any other \code{theta} rebuilds
#'   through \code{model_fn} as usual). It is NOT checked against
#'   \code{model_fn(theta)} -- that check would cost the rebuild being saved.
#' @param theta The parameter vector \code{model} was built at.
#'
#' @section Moving shock processes: the returned closure takes an OPTIONAL
#'   third argument overriding the construction-time \code{shock_specs}. This
#'   matters in an estimation loop, where the persistences move with every
#'   draw: the expensive half of this memo (the model rebuild and each
#'   parameter's DAG propagation, \code{\link{hank_model_dtheta}}) does not
#'   depend on \code{shock_specs} at all, while the cheap half
#'   (\code{\link{hank_dtheta_theta_list}}, an application of the propagated
#'   derivative to the driving paths) does. Rebuilding the whole closure when
#'   \code{rho} moves would throw the expensive half away; passing the new
#'   \code{shock_specs} per call keeps it.
#'
#' @return A function \code{(theta, param, shock_specs = NULL)} returning the
#'   per-shock \code{dTheta/dparam} matrices; \code{shock_specs = NULL} uses
#'   the ones this was constructed with.
#' @seealso \code{\link{hank_loglik_ar_structural_grad}} (whose
#'   \code{verify = TRUE} default checks the result of this against a central
#'   difference -- worth keeping for the first call, since a declared-derivative
#'   set that omits a channel is exactly what this cannot self-detect)
#' @export
hank_dtheta_fn <- function(model_fn, shock_specs, obs_vars, model = NULL,
                           theta = NULL) {
  memo <- new.env(parent = emptyenv())
  memo$theta <- NULL; memo$model <- NULL; memo$dt <- list()
  if (!is.null(model)) {
    if (!inherits(model, "hank_model"))
      stop("hank_dtheta_fn(): `model` must be a hank_model() object or NULL.")
    if (is.null(theta))
      stop("hank_dtheta_fn(): seeding with `model` requires the `theta` it ",
           "was built at (the memo is keyed on theta, so an unkeyed model ",
           "would be reused at the wrong parameter vector).")
    memo$theta <- theta; memo$model <- model
  }
  function(theta, param, shock_specs_now = NULL) {
    if (is.null(memo$theta) || !identical(memo$theta, theta)) {
      memo$theta <- theta
      memo$model <- model_fn(theta)
      memo$dt <- list()
      if (!inherits(memo$model, "hank_model"))
        stop("hank_dtheta_fn(): `model_fn` must return a hank_model().")
    }
    if (is.null(memo$dt[[param]]))
      memo$dt[[param]] <- hank_model_dtheta(memo$model, param)
    ## The memo above is shock_specs-INDEPENDENT; only this last application
    ## is not, so a moving rho costs the application and not the rebuild.
    hank_dtheta_theta_list(memo$model,
                           shock_specs_now %||% shock_specs,
                           obs_vars, memo$dt[[param]])
  }
}
