## R/hank-model.R
## --------------------------------------------------------------------------
## General sequence-space model composition: assemble an arbitrary directed
## acyclic graph (DAG) of blocks into a general-equilibrium system and solve it,
## the sequence-jacobian `create_model` analog.  Generalizes the
## Krusell-Smith-specific solve in R/hank-ge.R.
##
## A model is a set of BLOCKS, each mapping named aggregate INPUT paths to named
## aggregate OUTPUT paths, plus:
##   - unknowns  U : endogenous aggregate paths solved for (free inputs),
##   - targets   H : block outputs that must be zero in equilibrium,
##   - exogenous Z : driving shock paths (free inputs).
##
## Each block exposes its sequence-space Jacobian J[out][in] (T x T). Simple
## blocks get theirs analytically (if supplied) or by finite differences around
## the steady state; het blocks via the fake-news algorithm
## (hank_het_jacobian).  The total Jacobian of every variable w.r.t. every
## source (unknown or exogenous) is accumulated along a topological order by the
## chain rule, giving H_U = d(targets)/d(unknowns) and H_Z = d(targets)/d(exog).
## Equilibrium: dU = -H_U^{-1} H_Z dZ  (ABRS 2021, eq. 30).
## --------------------------------------------------------------------------


#' Define a simple (representative-agent) sequence-space block
#'
#' A simple block maps aggregate input LEVEL paths to output LEVEL paths via a
#' deterministic function \code{fn(paths, ss)}, where \code{paths} is a named
#' list of length-\code{T} input paths and \code{ss} is a named list of
#' steady-state levels for all model variables (used to pad lags/leads, e.g.
#' \code{c(ss$K, K[-T])} for \code{K(-1)}).
#'
#' @param name Character block name.
#' @param inputs,outputs Character vectors of input/output variable names.
#' @param fn Function \code{fn(paths, ss)} returning a named list of length-T
#'   output paths.
#' @param jac Optional function \code{jac(ss, T_h)} returning the analytic block
#'   Jacobian as a nested list \code{[[output]][[input]]} of \code{T x T}
#'   matrices.  If \code{NULL}, the Jacobian is computed by central finite
#'   differences around the steady state.
#'
#' @return An object of class \code{hank_block} (kind \code{"simple"}).
#' @export
hank_simple_block <- function(name, inputs, outputs, fn, jac = NULL) {
  structure(list(name = name, kind = "simple", inputs = inputs,
                 outputs = outputs, fn = fn, jac = jac),
            class = "hank_block")
}


#' Wrap a heterogeneous-agent household as a sequence-space block
#'
#' @param name Character block name.
#' @param block A \code{\link{hank_het_block}} solved at steady state.
#' @param inputs,outputs Character vectors; must be supported by
#'   \code{\link{hank_het_jacobian}} (inputs a subset of \code{c("r","w")},
#'   outputs a subset of \code{c("A","C")}).
#'
#' @return An object of class \code{hank_block} (kind \code{"het"}).
#' @export
hank_het_block_spec <- function(name, block, inputs = c("r", "w"),
                                outputs = c("A", "C")) {
  structure(list(name = name, kind = "het", inputs = inputs,
                 outputs = outputs, block = block),
            class = "hank_block")
}


#' Wrap a K-type discount-factor MIXTURE household as a sequence-space block
#'
#' Generalizes \code{\link{hank_het_block_spec}} to a household block whose
#' aggregate response is a \code{K}-type preference mixture (see
#' \code{\link{hank_mixture_blocks}}): every type shares the asset grid and
#' prices, and MAY differ in its income process (\code{Pi}, \code{e}) -- e.g.
#' a per-type Rouwenhorst calibration for an income-risk heterogeneity axis --
#' as well as in \code{beta}. Because types interact solely through the common
#' aggregate prices \code{(r, w)}, the block's sequence-space Jacobian and
#' nonlinear map for the SCALAR aggregates (A, C) are EXACT omega-weighted
#' sums of the per-type objects (see \code{\link{hank_mixture_jacobian}},
#' \code{\link{hank_mixture_dist_jacobian}}) -- no cross term between types --
#' regardless of whether the income processes match. Pooling the per-type
#' stationary DISTRIBUTIONS into one cell-space object (\code{D} in
#' \code{\link{hank_mixture_dist}} etc.) is a separate question that DOES
#' require a common \code{(e, a)} cell space; see \code{same_income} below and
#' the guards in \code{\link{hank_mixture_dist}},
#' \code{\link{hank_mixture_dist_jacobian}}, and
#' \code{\link{hank_mixture_td_nonlinear}}.
#'
#' @param name Character block name.
#' @param blocks List of \code{\link{hank_het_block}} objects sharing an
#'   asset grid (as returned by \code{\link{hank_mixture_blocks}}, or built
#'   per-type with distinct income processes on a shared \code{a_grid}).
#' @param omega Numeric length-\code{length(blocks)} mixture weights,
#'   non-negative, summing to 1.
#' @param inputs,outputs Character vectors; must be supported by
#'   \code{\link{hank_mixture_jacobian}} (inputs a subset of
#'   \code{c("r","w")}, outputs a subset of \code{c("A","C")}).
#'
#' @return An object of class \code{hank_block} (kind \code{"het_mixture"}),
#'   with an additional logical field \code{same_income}: \code{TRUE} iff
#'   every block shares the identical \code{(Pi, e)} income process (so
#'   distribution-pooling functions may safely return a single summed cell
#'   vector), \code{FALSE} otherwise (downstream distribution-pooling
#'   functions then return per-type representations instead).
#' @export
hank_mixture_block_spec <- function(name, blocks, omega, inputs = c("r", "w"),
                                    outputs = c("A", "C")) {
  if (!is.list(blocks) || length(blocks) < 1L ||
      !all(vapply(blocks, inherits, logical(1), what = "hank_het_block")))
    stop("hank_mixture_block_spec(): 'blocks' must be a non-empty list of ",
         "hank_het_block objects.")
  grid_match <- vapply(blocks[-1L], function(b)
    isTRUE(all.equal(b$a_grid, blocks[[1L]]$a_grid)), logical(1))
  if (length(blocks) > 1L && !all(grid_match))
    stop("hank_mixture_block_spec(): all 'blocks' must share the same ",
         "asset grid ('a_grid'); the income process ('Pi', 'e') may differ ",
         "across types.")
  income_match <- vapply(blocks[-1L], function(b)
    isTRUE(all.equal(b$Pi, blocks[[1L]]$Pi)) &&
      isTRUE(all.equal(b$e, blocks[[1L]]$e)),
    logical(1))
  same_income <- length(blocks) <= 1L || all(income_match)
  .hank_mixture_check_omega(blocks, omega)
  structure(list(name = name, kind = "het_mixture", inputs = inputs,
                 outputs = outputs, blocks = blocks, omega = omega,
                 same_income = same_income),
            class = "hank_block")
}


#' Finite-difference sequence-space Jacobian of a simple block
#' @keywords internal
.hank_simple_jac_fd <- function(blk, ss, T_h, delta = 1e-6) {
  base <- lapply(blk$inputs, function(i) rep(ss[[i]], T_h))
  names(base) <- blk$inputs
  J <- setNames(lapply(blk$outputs, function(o)
    setNames(lapply(blk$inputs, function(i) matrix(0, T_h, T_h)), blk$inputs)),
    blk$outputs)
  for (i in blk$inputs) for (s in seq_len(T_h)) {
    pp <- base; pp[[i]][s] <- pp[[i]][s] + delta
    pm <- base; pm[[i]][s] <- pm[[i]][s] - delta
    yp <- blk$fn(pp, ss); ym <- blk$fn(pm, ss)
    for (o in blk$outputs)
      J[[o]][[i]][, s] <- (yp[[o]] - ym[[o]]) / (2 * delta)
  }
  J
}


#' Block Jacobian dispatch (simple: analytic or FD; het: fake-news)
#' @keywords internal
.hank_block_jacobian <- function(blk, ss, T_h) {
  if (blk$kind == "het")
    return(hank_het_jacobian(blk$block, T_h, inputs = blk$inputs,
                             outputs = blk$outputs))
  if (blk$kind == "het_mixture")
    return(hank_mixture_jacobian(blk$blocks, blk$omega, T_h,
                                 inputs = blk$inputs, outputs = blk$outputs))
  if (!is.null(blk$jac)) return(blk$jac(ss, T_h))
  .hank_simple_jac_fd(blk, ss, T_h)
}


#' Topologically order blocks by input/output dependencies
#' @keywords internal
.hank_topo_order <- function(blocks, sources) {
  produced_by <- list()
  for (b in seq_along(blocks))
    for (o in blocks[[b]]$outputs) produced_by[[o]] <- b
  available <- as.list(setNames(rep(TRUE, length(sources)), sources))
  done <- rep(FALSE, length(blocks)); order <- integer(0)
  repeat {
    progressed <- FALSE
    for (b in seq_along(blocks)) {
      if (done[b]) next
      ins <- blocks[[b]]$inputs
      if (all(vapply(ins, function(i) isTRUE(available[[i]]), logical(1)))) {
        order <- c(order, b); done[b] <- TRUE; progressed <- TRUE
        for (o in blocks[[b]]$outputs) available[[o]] <- TRUE
      }
    }
    if (all(done)) break
    if (!progressed)
      stop("hank_model: block DAG has a cycle or an undefined input ",
           "(an input is neither a source nor produced by any block).")
  }
  order
}


#' Assemble and solve a general sequence-space GE model
#'
#' @param blocks List of \code{\link{hank_simple_block}} /
#'   \code{\link{hank_het_block_spec}} objects.
#' @param unknowns Character: endogenous aggregate paths to solve for.
#' @param targets Character: block outputs that must be zero in equilibrium
#'   (same length as \code{unknowns}).
#' @param exogenous Character: driving shock paths.
#' @param ss Named list/vector of steady-state levels for all model variables
#'   (used by simple-block FD Jacobians to pad lags).
#' @param T_h Integer horizon.
#'
#' @return An object of class \code{hank_model} with the packed GE Jacobians
#'   \code{H_U}, \code{H_Z} (each \code{T*n x T*n}), the accumulated per-variable
#'   Jacobians \code{G} (\code{G[[var]][[source]]}), and metadata. Pass it to
#'   \code{\link{hank_model_irf}}.
#' @export
hank_model <- function(blocks, unknowns, targets, exogenous, ss, T_h) {
  if (length(unknowns) != length(targets))
    stop("hank_model: need one target per unknown (square GE system).")
  sources <- c(unknowns, exogenous)
  ord <- .hank_topo_order(blocks, sources)

  I <- diag(T_h)
  ## G[[var]][[source]] : d(var path)/d(source path), a T x T block.
  G <- list()
  for (s in sources) { G[[s]] <- list(); G[[s]][[s]] <- I }

  for (b in ord) {
    blk <- blocks[[b]]
    Jb  <- .hank_block_jacobian(blk, ss, T_h)
    for (o in blk$outputs) {
      Go <- setNames(vector("list", length(sources)), sources)
      for (s in sources) {
        acc <- NULL
        for (i in blk$inputs) {
          Gis <- G[[i]][[s]]
          if (is.null(Gis)) next            # input does not depend on source s
          contrib <- Jb[[o]][[i]] %*% Gis
          acc <- if (is.null(acc)) contrib else acc + contrib
        }
        if (!is.null(acc)) Go[[s]] <- acc
      }
      G[[o]] <- Go
    }
  }

  ## Pack H_U = d(targets)/d(unknowns), H_Z = d(targets)/d(exogenous).
  pack <- function(cols) {
    M <- matrix(0, T_h * length(targets), T_h * length(cols))
    for (ti in seq_along(targets)) for (cj in seq_along(cols)) {
      blk <- G[[targets[ti]]][[cols[cj]]]
      if (is.null(blk)) next
      M[((ti - 1) * T_h + 1):(ti * T_h),
        ((cj - 1) * T_h + 1):(cj * T_h)] <- blk
    }
    M
  }
  H_U <- pack(unknowns); H_Z <- pack(exogenous)

  structure(list(H_U = H_U, H_Z = H_Z, G = G,
                 unknowns = unknowns, targets = targets,
                 exogenous = exogenous, T_h = T_h, block_order = ord,
                 blocks = blocks, ss = ss),
            class = "hank_model")
}


#' Evaluate all model variables nonlinearly along candidate paths
#'
#' Runs each block's nonlinear map in topological order given LEVEL paths for
#' the sources (unknowns + exogenous), returning level paths for every variable
#' (simple blocks via their \code{fn}; het blocks via
#' \code{\link{hank_td_nonlinear}}).
#' @keywords internal
.hank_model_eval <- function(model, src_paths) {
  T_h <- model$T_h; ss <- model$ss
  vals <- src_paths                       # sources (unknowns + exogenous)
  for (b in model$block_order) {
    blk <- model$blocks[[b]]
    ins <- lapply(blk$inputs, function(i) vals[[i]]); names(ins) <- blk$inputs
    if (blk$kind == "het") {
      r_path <- ins[["r"]]; w_path <- ins[["w"]]
      td <- hank_td_nonlinear(blk$block, r_path = r_path, w_path = w_path,
                              T_h = T_h)
      for (o in blk$outputs) vals[[o]] <- td[[o]]
    } else if (blk$kind == "het_mixture") {
      r_path <- ins[["r"]]; w_path <- ins[["w"]]
      omega <- blk$omega
      td_k <- lapply(blk$blocks, function(bk)
        hank_td_nonlinear(bk, r_path = r_path, w_path = w_path, T_h = T_h))
      for (o in blk$outputs) {
        acc <- omega[1L] * td_k[[1L]][[o]]
        if (length(blk$blocks) > 1L)
          for (k in 2L:length(blk$blocks)) acc <- acc + omega[k] * td_k[[k]][[o]]
        vals[[o]] <- acc
      }
    } else {
      outs <- blk$fn(ins, ss)
      for (o in blk$outputs) vals[[o]] <- outs[[o]]
    }
  }
  vals
}


#' Nonlinear GE transition of a general sequence-space model (quasi-Newton)
#'
#' Solves the fully nonlinear perfect-foresight transition \eqn{H(U, Z) = 0} for
#' the unknown paths by Newton's method with the frozen steady-state Jacobian
#' \code{H_U} (ABRS eq. 38); each residual re-evaluates every block nonlinearly
#' via \code{\link{hank_td_nonlinear}} (het) and the block functions (simple).
#'
#' @param model A \code{\link{hank_model}}.
#' @param Z_paths Named list of exogenous LEVEL paths (length \code{T}), one per
#'   \code{model$exogenous}.
#' @param tol,maxit Newton tolerance (max abs target residual) and iteration cap.
#'
#' @return A list with level paths for every variable, plus \code{converged},
#'   \code{iterations}, \code{max_resid}.
#' @export
hank_model_nonlinear_irf <- function(model, Z_paths, tol = 1e-9, maxit = 50L) {
  T_h <- model$T_h; ss <- model$ss
  ## unknown level paths, initialized at steady state
  U <- setNames(lapply(model$unknowns, function(u) rep(ss[[u]], T_h)),
                model$unknowns)
  stack <- function(lst, nm) do.call(c, lst[nm])

  converged <- FALSE; it <- 0L; max_resid <- Inf
  for (it in seq_len(maxit)) {
    src <- c(U, Z_paths)
    vals <- .hank_model_eval(model, src)
    resid <- stack(vals, model$targets)          # targets must be 0
    max_resid <- max(abs(resid))
    if (max_resid < tol) { converged <- TRUE; break }
    dU <- as.numeric(solve(model$H_U, resid))
    U_stack <- stack(U, model$unknowns) - dU
    for (k in seq_along(model$unknowns))
      U[[model$unknowns[k]]] <- U_stack[((k - 1) * T_h + 1):(k * T_h)]
  }
  vals <- .hank_model_eval(model, c(U, Z_paths))
  vals$converged <- converged; vals$iterations <- it; vals$max_resid <- max_resid
  vals
}


#' Shared GE-solve step: per-source deviation paths from the unknowns/exogenous
#'
#' Solves \eqn{dU = -H_U^{-1} H_Z\, dZ} and returns the per-source (unknowns +
#' exogenous) deviation paths \code{dsrc}. Factored out of
#' \code{\link{hank_model_irf}} so that \code{\link{hank_model_dist_irf}} can
#' reuse the identical GE solve (including the ill-conditioning guard) without
#' duplicating it; \code{hank_model_irf} calls this and then only adds the
#' propagation to scalar variables via \code{model$G}.
#'
#' @param model A \code{\link{hank_model}}.
#' @param dZ Named list of length-\code{T} exogenous shock paths (deviations),
#'   one per \code{model$exogenous} (missing entries treated as zero).
#'
#' @return Named list \code{dsrc}: one length-\code{T_h} deviation path per
#'   entry of \code{model$unknowns} and \code{model$exogenous}.
#' @keywords internal
.hank_irf_dsrc <- function(model, dZ) {
  T_h <- model$T_h
  z_stack <- do.call(c, lapply(model$exogenous, function(z) {
    v <- dZ[[z]]; if (is.null(v)) rep(0, T_h) else v
  }))
  ## Well-posedness guard: a near-singular H_U silently yields a garbage IRF.
  rc <- rcond(model$H_U)
  if (!is.finite(rc) || rc < 1e-10)
    warning(sprintf(paste0("hank_model_irf(): H_U is ill-conditioned ",
                           "(rcond = %.2e); the GE solution may be unreliable ",
                           "(near-singular / indeterminate). See hank_determinacy()."),
                    rc))
  dU_stack <- as.numeric(-solve(model$H_U, model$H_Z %*% z_stack))

  ## per-source deviation paths
  dsrc <- list()
  for (k in seq_along(model$unknowns))
    dsrc[[model$unknowns[k]]] <- dU_stack[((k - 1) * T_h + 1):(k * T_h)]
  for (z in model$exogenous)
    dsrc[[z]] <- { v <- dZ[[z]]; if (is.null(v)) rep(0, T_h) else v }
  dsrc
}


#' Linear GE impulse response from a general sequence-space model
#'
#' Solves \eqn{dU = -H_U^{-1} H_Z\, dZ} and propagates to every model variable.
#'
#' @param model A \code{\link{hank_model}}.
#' @param dZ Named list of length-\code{T} exogenous shock paths (deviations),
#'   one per \code{model$exogenous} (missing entries treated as zero).
#'
#' @return A named list of deviation paths: the unknowns, plus every produced
#'   variable, plus the exogenous inputs. Also carries attribute
#'   \code{"target_resid"} (should be ~0).
#' @export
hank_model_irf <- function(model, dZ) {
  T_h <- model$T_h
  dsrc <- .hank_irf_dsrc(model, dZ)

  ## propagate to all variables via the accumulated G
  out <- list()
  for (v in names(model$G)) {
    dv <- rep(0, T_h)
    for (s in model$exogenous) {
      Gs <- model$G[[v]][[s]]
      if (!is.null(Gs)) dv <- dv + as.numeric(Gs %*% dsrc[[s]])
    }
    for (s in model$unknowns) {
      Gs <- model$G[[v]][[s]]
      if (!is.null(Gs)) dv <- dv + as.numeric(Gs %*% dsrc[[s]])
    }
    out[[v]] <- dv
  }
  ## target residual (should be ~0 by construction)
  resid <- vapply(model$targets, function(t) max(abs(out[[t]])), numeric(1))
  attr(out, "target_resid") <- resid
  out
}


#' Linear GE DISTRIBUTION impulse response from a general sequence-space model
#'
#' Companion to \code{\link{hank_model_irf}}: instead of propagating the GE
#' solution to scalar aggregate variables via \code{model$G}, propagates it to
#' the full cross-sectional distribution response of every \code{"het"} (and
#' \code{"het_mixture"}) block in the model, via that block's distribution
#' Jacobian (\code{\link{hank_het_dist_jacobian}}, or for a mixture block
#' \code{\link{hank_mixture_dist_jacobian}} -- the omega-weighted sum of the
#' per-type distribution Jacobians on the types' shared grid).
#'
#' Reuses the identical \code{dU = -H_U^{-1} H_Z\, dZ} GE solve as
#' \code{hank_model_irf} (via the shared internal \code{.hank_irf_dsrc}), so
#' the two functions always agree on the per-source deviation paths
#' \code{dsrc}; only the propagation step differs.
#'
#' For each het block \code{b}, the relevant deviation paths are those of
#' \code{b$inputs} itself (e.g. \code{r}, \code{w}) -- which are typically
#' intermediate variables produced by some OTHER block (e.g. a firm block),
#' not themselves unknowns or exogenous sources, so they do not appear
#' directly in \code{dsrc}. Exactly as \code{\link{hank_model_irf}} does for
#' every variable \code{v}, each input path is reconstructed by propagating
#' \code{dsrc} through the accumulated chain-rule Jacobian
#' \code{model$G[[input]][[s]]} over every source \code{s} (unknowns and
#' exogenous). Given that reconstructed deviation path \code{dI[i]} for each
#' input \code{i} in \code{b$inputs}, the distribution response at date
#' \code{t} is
#' \deqn{dD_t = \Sigma_i \Sigma_s J^D_i[t, s, ]\, dI[i][s]}
#' i.e. a sum over inputs and shock dates of the (per-cell) distribution
#' Jacobian rows, contracted against the realized input deviation path.
#' Recall (see \code{\link{hank_het_dist_jacobian}}) that row \code{t = 1} is
#' identically zero for every input/date, since the initial distribution is a
#' predetermined state that no anticipated shock can move.
#'
#' @param model A \code{\link{hank_model}}.
#' @param dZ Named list of length-\code{T} exogenous shock paths (deviations),
#'   one per \code{model$exogenous} (missing entries treated as zero).
#'
#' @return A list with:
#'   \item{dD}{Named list, one entry per \code{"het"}-kind block in
#'     \code{model$blocks} (keyed by block name). Each entry is an
#'     \code{(n_e*n_a) x T_h} matrix; column \code{t} is the distribution
#'     deviation \code{dD_t} for that block.}
#'   \item{dsrc}{The per-source (unknowns + exogenous) deviation paths, as
#'     returned internally by \code{\link{hank_model_irf}} -- exposed here so
#'     callers can confirm the two functions solved the identical GE system.}
#' @export
hank_model_dist_irf <- function(model, dZ) {
  T_h <- model$T_h
  dsrc <- .hank_irf_dsrc(model, dZ)

  ## Reconstruct the deviation path of an arbitrary model variable v (not
  ## necessarily a source) by propagating dsrc through the accumulated
  ## chain-rule Jacobian model$G[[v]][[s]] -- identical to the propagation
  ## step inside hank_model_irf(), just restricted to the one variable v
  ## instead of looped over every produced variable.
  dev_path <- function(v) {
    dv <- rep(0, T_h)
    for (s in c(model$exogenous, model$unknowns)) {
      Gs <- model$G[[v]][[s]]
      if (!is.null(Gs)) dv <- dv + as.numeric(Gs %*% dsrc[[s]])
    }
    dv
  }

  ## Contract a distribution Jacobian JD[[input]] (T_h x T_h x n_cell) against
  ## the realized input deviation paths: dD[,t] = Sum_i Sum_s JD[i][t,s,] *
  ## dsrc[i][s] -- the SAME contraction for a single het block and for a
  ## het_mixture block (whose JD is already the omega-weighted sum; see
  ## hank_mixture_dist_jacobian), just parameterized by n_cell and the input
  ## deviation-path lookup.
  contract_dist <- function(JD, inputs, n_cell) {
    dD_b <- matrix(0, n_cell, T_h)
    for (i in inputs) {
      ## blk$inputs (e.g. "r", "w") are generally intermediate variables
      ## produced by another block, NOT sources -- so look them up via
      ## dsrc[[i]] when i happens to be a source (unknown/exogenous itself),
      ## else reconstruct via dev_path(). Checking dsrc first avoids a
      ## redundant recomputation in the (rare) case a het block's input is
      ## itself a raw source.
      dz <- if (!is.null(dsrc[[i]])) dsrc[[i]] else dev_path(i)
      JDi <- JD[[i]]                        # T_h x T_h x n_cell
      for (tt in seq_len(T_h))
        dD_b[, tt] <- dD_b[, tt] +
          as.numeric(crossprod(matrix(JDi[tt, , ], nrow = T_h, ncol = n_cell), dz))
    }
    dD_b
  }

  dD <- list()
  for (blk in model$blocks) {
    if (identical(blk$kind, "het")) {
      JD <- hank_het_dist_jacobian(blk$block, T_h, inputs = blk$inputs)
      n_cell <- blk$block$n_e * blk$block$n_a
      dD[[blk$name]] <- contract_dist(JD, blk$inputs, n_cell)
    } else if (identical(blk$kind, "het_mixture")) {
      ## The pooled mixture distribution is only well-defined on a SHARED
      ## (e, a) cell space (all blk$blocks sharing a_grid/Pi/e -- see
      ## hank_mixture_block_spec's validation and same_income flag). When the
      ## types differ in their income process, skip the pooled dD entry here
      ## (rather than silently summing incomparable per-type cell vectors);
      ## callers needing per-type distribution paths for a hetinc mixture
      ## should call hank_mixture_dist_jacobian() per type directly.
      if (isTRUE(blk$same_income)) {
        JD_mix <- hank_mixture_dist_jacobian(blk$blocks, blk$omega, T_h,
                                             inputs = blk$inputs)
        n_cell <- blk$blocks[[1L]]$n_e * blk$blocks[[1L]]$n_a
        dD[[blk$name]] <- contract_dist(JD_mix, blk$inputs, n_cell)
      }
    }
  }

  list(dD = dD, dsrc = dsrc)
}


#' Build the Krusell-Smith model through the general block-DAG engine
#'
#' Convenience constructor assembling the KS economy (firm + household + asset
#' market) as a \code{\link{hank_model}}, for use and as a regression check that
#' the general engine reproduces the KS-specific solve in \code{R/hank-ge.R}.
#'
#' @param ks A \code{\link{hank_ks_steady}} steady state.
#' @param T_h Integer horizon.
#' @return A \code{\link{hank_model}} with unknown \code{K}, target
#'   \code{asset_mkt}, exogenous \code{Z}.
#' @export
hank_ks_model <- function(ks, T_h) {
  alpha <- ks$alpha; delta <- ks$delta

  firm <- hank_simple_block(
    "firm", inputs = c("K", "Z"), outputs = c("r", "w"),
    fn = function(paths, ss) {
      Klag <- c(ss$K, paths$K[-length(paths$K)])
      list(r = alpha * paths$Z * Klag^(alpha - 1) - delta,
           w = (1 - alpha) * paths$Z * Klag^(alpha))
    },
    jac = function(ss, T_h) {
      K <- ss$K; Z <- ss$Z
      lag <- rbind(0, cbind(diag(T_h - 1L), 0))
      list(r = list(K = alpha * (alpha - 1) * Z * K^(alpha - 2) * lag,
                    Z = alpha * K^(alpha - 1) * diag(T_h)),
           w = list(K = (1 - alpha) * alpha * Z * K^(alpha - 1) * lag,
                    Z = (1 - alpha) * K^(alpha) * diag(T_h)))
    })

  household <- hank_het_block_spec("household", ks$block,
                                   inputs = c("r", "w"), outputs = c("A", "C"))

  market <- hank_simple_block(
    "market", inputs = c("A", "K"), outputs = "asset_mkt",
    fn = function(paths, ss) list(asset_mkt = paths$A - paths$K))

  ss <- list(K = ks$K, Z = ks$Z, r = ks$r, w = ks$w,
             A = ks$block$A, C = ks$block$C, asset_mkt = 0)
  hank_model(list(firm, household, market),
             unknowns = "K", targets = "asset_mkt", exogenous = "Z",
             ss = ss, T_h = T_h)
}


#' Build the Krusell-Smith MIXTURE model through the general block-DAG engine
#'
#' Mirrors \code{\link{hank_ks_model}} exactly, replacing the single-type
#' household with a \code{\link{hank_mixture_block_spec}} household whose
#' aggregate response is the \code{K}-type discount-factor mixture's
#' omega-weighted sum (see \code{\link{hank_mixture_ks_steady}}). Firm and
#' market blocks are IDENTICAL to \code{hank_ks_model} (the firm's simple
#' block and market-clearing target depend only on aggregate \code{K}, \code{A}
#' regardless of what generates \code{A}).
#'
#' @param mks A \code{\link{hank_mixture_ks_steady}} steady state.
#' @param T_h Integer horizon.
#' @return A \code{\link{hank_model}} with unknown \code{K}, target
#'   \code{asset_mkt}, exogenous \code{Z}.
#' @export
hank_mixture_ks_model <- function(mks, T_h) {
  alpha <- mks$alpha; delta <- mks$delta

  firm <- hank_simple_block(
    "firm", inputs = c("K", "Z"), outputs = c("r", "w"),
    fn = function(paths, ss) {
      Klag <- c(ss$K, paths$K[-length(paths$K)])
      list(r = alpha * paths$Z * Klag^(alpha - 1) - delta,
           w = (1 - alpha) * paths$Z * Klag^(alpha))
    },
    jac = function(ss, T_h) {
      K <- ss$K; Z <- ss$Z
      lag <- rbind(0, cbind(diag(T_h - 1L), 0))
      list(r = list(K = alpha * (alpha - 1) * Z * K^(alpha - 2) * lag,
                    Z = alpha * K^(alpha - 1) * diag(T_h)),
           w = list(K = (1 - alpha) * alpha * Z * K^(alpha - 1) * lag,
                    Z = (1 - alpha) * K^(alpha) * diag(T_h)))
    })

  household <- hank_mixture_block_spec("household", mks$blocks, mks$omega,
                                       inputs = c("r", "w"), outputs = c("A", "C"))

  market <- hank_simple_block(
    "market", inputs = c("A", "K"), outputs = "asset_mkt",
    fn = function(paths, ss) list(asset_mkt = paths$A - paths$K))

  omega <- mks$omega
  A_ss <- sum(vapply(seq_along(mks$blocks), function(k) omega[k] * mks$blocks[[k]]$A,
                      numeric(1)))
  C_ss <- sum(vapply(seq_along(mks$blocks), function(k) omega[k] * mks$blocks[[k]]$C,
                      numeric(1)))
  ss <- list(K = mks$K, Z = mks$Z, r = mks$r, w = mks$w,
             A = A_ss, C = C_ss, asset_mkt = 0)
  hank_model(list(firm, household, market),
             unknowns = "K", targets = "asset_mkt", exogenous = "Z",
             ss = ss, T_h = T_h)
}
