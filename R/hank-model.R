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
#' @param djac_dtheta Optional named list of functions
#'   \code{function(ss, T_h)}, one per structural parameter this block's
#'   Jacobian depends on, each returning \eqn{\partial J/\partial\theta_k} in
#'   the same nested \code{[[output]][[input]]} shape as \code{jac} (omit an
#'   \code{[[output]][[input]]} entry that is identically zero).
#'   \code{\link{hank_model_dtheta}} propagates these through the DAG chain
#'   rule to \eqn{dH_U/d\theta_k}, \eqn{dH_Z/d\theta_k} and ultimately
#'   \eqn{d\Theta_z/d\theta_k} -- so a block only ever declares the derivative
#'   of its OWN small Jacobian, never a packed GE matrix. Used by the exact
#'   route of \code{\link{hank_loglik_ar_structural_grad}}; see
#'   \code{\link{hank_dtheta_fn}}.
#'
#' @return An object of class \code{hank_block} (kind \code{"simple"}).
#' @export
hank_simple_block <- function(name, inputs, outputs, fn, jac = NULL,
                              djac_dtheta = NULL) {
  if (!is.null(djac_dtheta)) {
    if (!is.list(djac_dtheta) || is.null(names(djac_dtheta)) ||
        any(!nzchar(names(djac_dtheta))))
      stop("hank_simple_block(): `djac_dtheta` must be a NAMED list ",
           "(one entry per structural parameter).")
    if (!all(vapply(djac_dtheta, is.function, logical(1))))
      stop("hank_simple_block(): every `djac_dtheta` entry must be a ",
           "function(ss, T_h) returning the same nested ",
           "[[output]][[input]] shape as `jac`.")
  }
  structure(list(name = name, kind = "simple", inputs = inputs,
                 outputs = outputs, fn = fn, jac = jac,
                 djac_dtheta = djac_dtheta),
            class = "hank_block")
}


#' Wrap a heterogeneous-agent household as a sequence-space block
#'
#' @param name Character block name.
#' @param block A \code{\link{hank_het_block}} solved at steady state.
#' @param inputs,outputs Character vectors; must be supported by
#'   \code{\link{hank_het_jacobian}} (inputs a subset of \code{c("r","w","Tr")}
#'   plus, for a block built with \code{Pi_fn}/\code{Pi_inputs}, the block's
#'   named transition-probability inputs -- e.g. \code{c("r","w","f","s")}
#'   for a \code{\link{hank_employment_income}} household whose job-finding
#'   and separation rates are produced by an upstream matching block;
#'   outputs a subset of \code{c("A","C")}). Validated here so a bad wiring
#'   fails at spec time, not inside the Jacobian dispatch.
#'
#' @return An object of class \code{hank_block} (kind \code{"het"}).
#' @export
hank_het_block_spec <- function(name, block, inputs = c("r", "w"),
                                outputs = c("A", "C")) {
  ## A two-asset block would otherwise slip through: the only validation below
  ## is on input NAMES, and a het2 block's defaults happen to satisfy it -- so
  ## it would be tagged kind = "het" and its n_e x n_b x n_a policies fed to
  ## the one-asset fake-news Jacobian, which reads them as n_e x n_a. Fail here
  ## instead, at spec time.
  if (inherits(block, "hank_het2_block"))
    stop("hank_het_block_spec(): this is a two-asset block ",
         "(hank_het2_block); use hank_het2_block_spec(), whose inputs are ",
         "('rb', 'ra', 'w') and outputs ('B', 'A', 'C').")
  if (!inherits(block, "hank_het_block"))
    stop("hank_het_block_spec(): 'block' must be a hank_het_block.")
  .hank_het_check_inputs(block, inputs)
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
#' as well as in \code{beta}, \code{eis}, and the borrowing constraint
#' \code{amin} (the WEALTH heterogeneity axis: per-type \code{amin} on the
#' shared \code{a_grid}; see \code{\link{hank_het_block}}). Because a
#' per-type \code{amin} lives on the shared grid, it keeps the
#' \code{(e, a)} cell space common across types -- so, unlike per-type
#' \code{Pi}/\code{e}, it does NOT flip \code{same_income} and pooled
#' distribution objects remain valid. Because types interact solely through
#' the common aggregate prices \code{(r, w)}, the block's sequence-space
#' Jacobian and
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


## Package-private block-Jacobian cache, same contract as .hank_ge_factor's LU
## cache: keyed by identical() on everything the Jacobian depends on, and a hit
## returns the BIT-IDENTICAL object the uncached path would have built.
##
## WHY. .hank_block_jacobian() is called once per block per hank_model() call,
## and for a het block that call IS the fake-news algorithm -- measured 0.465 s
## of a 0.558 s model build at T_h = 400 (83%), 0.239 s of 0.275 s at T_h = 200.
## Any workload that rebuilds the model at a new STRUCTURAL parameter pays it
## again, even when the parameter cannot touch the household at all (a Taylor
## rule coefficient, an NKPC slope): a frozen-Jacobian rebuild reproduces H_U
## to max|diff| = 0, so that work was provably redundant. Skipping it makes
## such a rebuild 16x cheaper at T_h = 400 (0.558 -> 0.035 s), which is what
## makes both affinity detection and FD-based dH_U/dtheta affordable -- see
## briefs/21-structural-score-api-scope.md.
##
## The key is list(blk, ss, T_h). Keying on the whole block is deliberate and
## conservative: for a het kind the Jacobian is a function of the solved block
## (its policies and Lambda), and for a simple kind of the closure in `fn`/`jac`
## -- and a closure rebuilt with new captured parameters is a DIFFERENT object
## to identical() (environments compare by reference), so a parameter change
## that only a closure can see still misses. `ss` is in the key because simple
## blocks difference around it. False misses cost a recomputation; there is no
## key under which a stale Jacobian can be returned.
##
## Several entries are retained (unlike the GE cache's single slot) because the
## caller pattern is a LOOP over blocks: one slot would thrash on every model
## with more than one cached block.
.hank_block_jac_cache <- new.env(parent = emptyenv())
.hank_block_jac_cache$entries <- list()
.hank_block_jac_cache$n_build <- 0L
.hank_block_jac_cache$n_hit   <- 0L
.hank_block_jac_max <- 16L

#' Reset the block-Jacobian cache (testing / memory reclamation)
#'
#' A cached het-block Jacobian is \code{n_inputs * n_outputs} dense
#' \code{T_h x T_h} matrices, so at production \code{T_h} the cache is large.
#' Call this to drop it, or set \code{options(dynhr.block_jac_cache = FALSE)}
#' to disable caching entirely.
#'
#' @return Invisibly, the telemetry counters as they stood before the reset.
#' @keywords internal
hank_block_jac_cache_reset <- function() {
  old <- list(n_build = .hank_block_jac_cache$n_build,
              n_hit   = .hank_block_jac_cache$n_hit,
              n_entries = length(.hank_block_jac_cache$entries))
  .hank_block_jac_cache$entries <- list()
  .hank_block_jac_cache$n_build <- 0L
  .hank_block_jac_cache$n_hit   <- 0L
  invisible(old)
}


#' Block Jacobian dispatch (simple: analytic or FD; het: fake-news)
#' @keywords internal
## WALL-CLOCK MUST NOT REACH THE CACHE KEY.
##
## The key below is compared with identical() on the whole block, which is
## deliberately conservative -- a false MISS costs a recomputation, and there
## is no key under which a stale Jacobian can be returned. But that same
## conservatism makes the key brittle to any run-VARYING field stored on the
## block: when hank_het_block()/hank_het2_block() started recording
## elapsed_solve/elapsed_dist (0.9.0.0035-36), two blocks built from an
## IDENTICAL calibration stopped being identical(), so the cache stopped
## hitting entirely -- silently, since a pure miss is still correct, just 16x
## slower on the rebuild this cache exists to make cheap.
## test-hank-block-jac-cache.R caught it.
##
## Only the timings are stripped. iterations/converged/gaps are deterministic
## given the calibration, and backend/threads are left IN the key on purpose:
## the R and compiled paths can differ at round-off, so keeping them cannot
## return a cross-backend result under the wrong key. Strip the minimum that
## restores the cache, not everything that looks like metadata.
.hank_block_cache_key <- function(blk) {
  drop <- c("elapsed_solve", "elapsed_dist", "elapsed")
  strip <- function(x) {
    if (!is.list(x)) return(x)
    x[intersect(names(x), drop)] <- NULL
    x
  }
  blk <- strip(blk)
  if (is.list(blk) && !is.null(blk$block)) blk$block <- strip(blk$block)
  blk
}


.hank_block_jacobian <- function(blk, ss, T_h) {
  if (!isTRUE(getOption("dynhr.block_jac_cache", TRUE)))
    return(.hank_block_jacobian_uncached(blk, ss, T_h))
  key <- list(blk = .hank_block_cache_key(blk), ss = ss, T_h = T_h)
  ents <- .hank_block_jac_cache$entries
  for (i in seq_along(ents)) {
    if (identical(ents[[i]]$key, key)) {
      .hank_block_jac_cache$n_hit <- .hank_block_jac_cache$n_hit + 1L
      ## move to front: the loop over blocks revisits the same few keys
      .hank_block_jac_cache$entries <- c(ents[i], ents[-i])
      return(ents[[i]]$J)
    }
  }
  J <- .hank_block_jacobian_uncached(blk, ss, T_h)
  .hank_block_jac_cache$n_build <- .hank_block_jac_cache$n_build + 1L
  ents <- c(list(list(key = key, J = J)), ents)
  if (length(ents) > .hank_block_jac_max)
    ents <- ents[seq_len(.hank_block_jac_max)]
  .hank_block_jac_cache$entries <- ents
  J
}


#' Block Jacobian dispatch, uncached (the reference path)
#' @keywords internal
.hank_block_jacobian_uncached <- function(blk, ss, T_h) {
  if (blk$kind == "het")
    return(hank_het_jacobian(blk$block, T_h, inputs = blk$inputs,
                             outputs = blk$outputs))
  if (blk$kind == "het2")
    return(hank_het2_jacobian(blk$block, T_h, inputs = blk$inputs,
                              outputs = blk$outputs))
  if (blk$kind == "het2d")
    return(hank_het2d_jacobian(blk$block, T_h, inputs = blk$inputs,
                               outputs = blk$outputs))
  if (blk$kind == "het3")
    return(hank_het3_jacobian(blk$block, T_h, inputs = blk$inputs,
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

  ## Fail loud on structural singularity: an unknown whose H_U column is
  ## identically zero (moves no target) or a target whose H_U row is
  ## identically zero (responds to no unknown) makes the GE Newton solve
  ## (.hank_ge_factor/.hank_ge_solve) singular. Check the raw per-(target,
  ## unknown) blocks in G BEFORE packing, so the offending name is known;
  ## exact-zero test only (no tolerance: these are structural, not numerical).
  for (u in unknowns) {
    moves_a_target <- FALSE
    for (t in targets) {
      blk <- G[[t]][[u]]
      if (!is.null(blk) && any(blk != 0)) { moves_a_target <- TRUE; break }
    }
    if (!moves_a_target)
      stop("hank_model: unknown '", u, "' moves no target (its Jacobian ",
           "column is identically zero); the GE system is singular. Check ",
           "the block wiring: does '", u, "' actually feed, directly or ",
           "through the DAG, into a block that produces one of the ",
           "targets (", paste(targets, collapse = ", "), ")?")
  }
  for (t in targets) {
    responds_to_an_unknown <- FALSE
    for (u in unknowns) {
      blk <- G[[t]][[u]]
      if (!is.null(blk) && any(blk != 0)) { responds_to_an_unknown <- TRUE; break }
    }
    if (!responds_to_an_unknown)
      stop("hank_model: target '", t, "' responds to no unknown (its ",
           "Jacobian row is identically zero); the GE system is singular. ",
           "Check the block wiring: does the block producing '", t, "' ",
           "actually depend, directly or through the DAG, on one of the ",
           "unknowns (", paste(unknowns, collapse = ", "), ")?")
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
            class = c("hank_model", "hank_block"))
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
      ## transition-probability inputs (HANK+SAM: e.g. f/s from an upstream
      ## matching block) ride along by name; NULL when the block is (r, w)-only.
      ## "Tr" is an AGGREGATE input (the lump-sum transfer), not a Pi input --
      ## without this exclusion it would be misrouted into pi_input_paths and
      ## die inside .hank_pi_path.
      pi_nm <- setdiff(blk$inputs, c("r", "w", "Tr", "r_minus"))
      pip   <- if (length(pi_nm)) ins[pi_nm] else NULL
      td <- hank_td_nonlinear(blk$block, r_path = r_path, w_path = w_path,
                              T_h = T_h, pi_input_paths = pip,
                              Tr_path = ins[["Tr"]],
                              r_minus_path = ins[["r_minus"]])
      for (o in blk$outputs) vals[[o]] <- td[[o]]
    } else if (blk$kind == "het2") {
      ## Two-asset household: three prices (rb, ra, w) rather than (r, w);
      ## transition-probability inputs ride along by name exactly as above.
      pi_nm <- setdiff(blk$inputs, c("rb", "ra", "w", "Tr", "theta_coll"))
      pip   <- if (length(pi_nm)) ins[pi_nm] else NULL
      td <- hank_td2_nonlinear(blk$block, rb_path = ins[["rb"]],
                               ra_path = ins[["ra"]], w_path = ins[["w"]],
                               T_h = T_h, pi_input_paths = pip,
                               Tr_path = ins[["Tr"]],
                               theta_path = ins[["theta_coll"]])
      for (o in blk$outputs) vals[[o]] <- td[[o]]
    } else if (blk$kind == "het3") {
      ## Three-asset household: four prices (rd, rf, ra, w) plus the optional
      ## foreign valuation px; transition-probability inputs ride along by
      ## name exactly as in the one- and two-asset arms above. px is absent
      ## from `ins` unless declared, and hank_td3_nonlinear reads NULL as
      ## "hold the block's steady-state px" -- so an undeclared valuation
      ## channel is held fixed rather than silently set to 1.
      pi_nm <- setdiff(blk$inputs, c("rd", "rf", "ra", "w", "px"))
      pip   <- if (length(pi_nm)) ins[pi_nm] else NULL
      td <- hank_td3_nonlinear(blk$block, rd_path = ins[["rd"]],
                               rf_path = ins[["rf"]], ra_path = ins[["ra"]],
                               w_path = ins[["w"]], px_path = ins[["px"]],
                               T_h = T_h, pi_input_paths = pip)
      for (o in blk$outputs) vals[[o]] <- td[[o]]
    } else if (blk$kind == "het2d") {
      ## Transition-probability inputs ride along by name, as in every other
      ## household arm.
      pi_nm <- setdiff(blk$inputs, c("rb", "ra", "w", "Tr"))
      pip   <- if (length(pi_nm)) ins[pi_nm] else NULL
      td <- hank_td2d_nonlinear(blk$block, rb_path = ins[["rb"]],
                                ra_path = ins[["ra"]], w_path = ins[["w"]],
                                Tr_path = ins[["Tr"]], T_h = T_h,
                                pi_input_paths = pip)
      for (o in blk$outputs) vals[[o]] <- td[[o]]
    } else if (blk$kind == "het_mixture") {
      ## The mixture path forwards ONLY (r, w) to the per-type transitions.
      ## Any other requested input (a Tr transfer, a Pi input) would be
      ## silently dropped here -- worse than an error, so refuse loudly.
      extra <- setdiff(blk$inputs, c("r", "w"))
      if (length(extra))
        stop("hank_model: mixture household block '", blk$name,
             "' requests input(s) ", paste0("'", extra, "'", collapse = ", "),
             ", but the nonlinear mixture transition forwards only ('r', ",
             "'w') to its type blocks; the extra input(s) would be silently ",
             "ignored. Use a single het block, or extend the mixture path.")
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
#' @param Z_paths Named list of exogenous LEVEL paths (length \code{T}) for a
#'   SUBSET of \code{model$exogenous}; any exogenous not named stays at its
#'   steady-state level for the whole horizon (the nonlinear analogue of
#'   \code{\link{hank_model_irf}}'s missing-shock-is-zero-deviation
#'   convention). Unknown names are an error.
#' @param tol,maxit Newton tolerance (max abs target residual) and iteration cap.
#'
#' @return A list with level paths for every variable, plus \code{converged},
#'   \code{iterations}, \code{max_resid}.
#' @export
hank_model_nonlinear_irf <- function(model, Z_paths, tol = 1e-9, maxit = 50L) {
  T_h <- model$T_h; ss <- model$ss
  ## Validate and complete the exogenous paths. Previously an exogenous the
  ## caller did not name silently had NO path at all, and the failure surfaced
  ## as a cryptic length error deep inside a downstream block (latent until a
  ## model with two exogenous was shocked in only one of them).
  bad <- setdiff(names(Z_paths), model$exogenous)
  if (length(bad))
    stop("hank_model_nonlinear_irf: Z_paths name(s) ",
         paste0("'", bad, "'", collapse = ", "),
         " are not exogenous in this model (exogenous: ",
         paste0("'", model$exogenous, "'", collapse = ", "), ").")
  short <- names(Z_paths)[vapply(Z_paths, length, integer(1)) != T_h]
  if (length(short))
    stop("hank_model_nonlinear_irf: Z_paths entr", if (length(short) > 1L)
         "ies " else "y ", paste0("'", short, "'", collapse = ", "),
         " must have length T_h = ", T_h, ".")
  for (z in setdiff(model$exogenous, names(Z_paths)))
    Z_paths[[z]] <- rep(ss[[z]], T_h)
  ## unknown level paths, initialized at steady state
  U <- setNames(lapply(model$unknowns, function(u) rep(ss[[u]], T_h)),
                model$unknowns)
  stack <- function(lst, nm) do.call(c, lst[nm])

  ## Quasi-Newton with the FROZEN steady-state H_U: every iteration solves
  ## against the same matrix, so one factorization serves the whole loop.
  fac <- .hank_ge_factor(model$H_U)

  converged <- FALSE; it <- 0L; max_resid <- Inf
  for (it in seq_len(maxit)) {
    src <- c(U, Z_paths)
    vals <- .hank_model_eval(model, src)
    resid <- stack(vals, model$targets)          # targets must be 0
    max_resid <- max(abs(resid))
    if (max_resid < tol) { converged <- TRUE; break }
    dU <- as.numeric(.hank_ge_solve(fac, model$H_U, resid))
    U_stack <- stack(U, model$unknowns) - dU
    for (k in seq_along(model$unknowns))
      U[[model$unknowns[k]]] <- U_stack[((k - 1) * T_h + 1):(k * T_h)]
  }
  vals <- .hank_model_eval(model, c(U, Z_paths))
  vals$converged <- converged; vals$iterations <- it; vals$max_resid <- max_resid
  vals
}


## Package-private GE factorization cache: the LU of H_U and its rcond, keyed
## on H_U ITSELF via identical(). Every sequence-space consumer -- one IRF per
## shock, one per finite-difference tap in the exact-AR gradient, one Newton
## iteration per nonlinear transition -- re-solved against the SAME H_U, and
## each of those calls paid a fresh O(n^3) LU *and* a fresh O(n^3) rcond.
## Measured on a dense n = 800 system: rcond 0.0060s + solve 0.0067s per call,
## against 0.0033s to factor once and 0.0003s per subsequent solve. Both
## discarded costs are cubic, so the waste grows with T_h x n_unknowns (the
## NZ HANK paper's production system is n = 8800).
##
## Exactness: reusing a dgetrf factorization through dgetrs is what base
## solve() does internally (dgesv = dgetrf + dgetrs), so a cache hit returns
## BIT-IDENTICAL floats, not merely close ones -- verified in
## test-hank-ge-factor.R at two sizes, and the standard the package's other
## caches already hold themselves to (see .hank_ar_autocov_slab).
##
## Keying on identical(H_U) rather than a fingerprint is deliberate: it cannot
## collide, and it is what makes in-place mutation (`bad$H_U[, 1] <- 0`, which
## the determinacy tests do) invalidate correctly, since R copies on modify.
## Two alternating model copies thrash the single slot; that costs an O(n^2)
## compare and a refactor, never a wrong answer.
.hank_ge_cache <- new.env(parent = emptyenv())
.hank_ge_cache$n_factor <- 0L    # telemetry: real factorizations performed
.hank_ge_cache$n_hit    <- 0L

#' Cached LU factorization and reciprocal condition number of \code{H_U}
#'
#' @param H_U The GE Jacobian of targets with respect to unknowns.
#' @return A list with \code{lu} (a \code{Matrix} LU, or \code{NULL} if the
#'   factorization failed and callers should fall back to \code{solve}) and
#'   \code{rcond}.
#' @keywords internal
.hank_ge_factor <- function(H_U) {
  if (!isTRUE(getOption("dynhr.hank_ge_cache", TRUE)))
    return(list(lu = NULL, rcond = rcond(H_U)))

  ent <- .hank_ge_cache$ent
  if (!is.null(ent) && identical(ent$H_U, H_U)) {
    .hank_ge_cache$n_hit <- .hank_ge_cache$n_hit + 1L
    return(ent)
  }
  ## A singular or otherwise pathological H_U must not become a hard error
  ## here: hank_model_irf() has always warned and carried on (returning
  ## whatever solve() gives), and hank_determinacy() relies on that. Degrade
  ## to the uncached path instead of propagating a factorization failure.
  lu <- tryCatch(Matrix::lu(H_U), error = function(e) NULL,
                 warning = function(w) NULL)
  ent <- list(H_U = H_U, lu = lu, rcond = rcond(H_U))
  .hank_ge_cache$ent      <- ent
  .hank_ge_cache$n_factor <- .hank_ge_cache$n_factor + 1L
  ent
}

#' Solve \code{H_U x = B} through the cached factorization
#'
#' @param fac A \code{\link{.hank_ge_factor}} entry.
#' @param H_U The matrix \code{fac} was built from (used only on fallback).
#' @param B Right-hand side.
#' @return The solution, identical to \code{solve(H_U, B)}.
#' @keywords internal
.hank_ge_solve <- function(fac, H_U, B) {
  if (is.null(fac$lu)) return(solve(H_U, B))
  out <- tryCatch(as.matrix(Matrix::solve(fac$lu, B)), error = function(e) NULL)
  if (is.null(out)) solve(H_U, B) else out
}

#' Reset the GE factorization cache (testing / memory reclamation)
#'
#' The cached LU roughly doubles the resident size of \code{H_U}, which
#' matters at production \code{T_h}. Call this to drop it, or set
#' \code{options(dynhr.hank_ge_cache = FALSE)} to disable caching entirely.
#'
#' @return Invisibly, the telemetry counters as they stood before the reset.
#' @keywords internal
.hank_ge_cache_clear <- function() {
  old <- list(n_factor = .hank_ge_cache$n_factor,
              n_hit    = .hank_ge_cache$n_hit)
  rm(list = ls(.hank_ge_cache), envir = .hank_ge_cache)
  .hank_ge_cache$n_factor <- 0L
  .hank_ge_cache$n_hit    <- 0L
  invisible(old)
}


#' Validate (and, for a single-exogenous model, coerce) a \code{dZ} argument
#'
#' Every sequence-space GE entry point takes the exogenous impulse as a NAMED
#' list, one length-\code{T_h} path per \code{model$exogenous}. Handed a bare
#' numeric vector, the old code reached \code{dZ[[z]]} with a character index
#' and died with "subscript out of bounds" -- a message that names neither the
#' argument nor the convention. That is not hypothetical: the package's own
#' HANK vignette made exactly this call, and it was one of the vignette
#' execution failures in every \code{R CMD check} log from 2026-07-15 on.
#'
#' The trap is real because two irf producers in this package disagree by
#' design: \code{\link{hank_ks_linear_irf}} takes a BARE vector (its model has
#' one shock and it returns \code{d}-prefixed names), while the general
#' \code{\link{hank_model_irf}} takes a named list keyed on the model's own
#' exogenous names and returns bare variable names.
#'
#' A bare numeric vector is therefore accepted when -- and only when -- the
#' model has exactly one exogenous, where the intent is unambiguous. Every
#' other malformed input gets a message that names the expected keys.
#'
#' @param model A \code{\link{hank_model}}.
#' @param dZ The user's \code{dZ} argument.
#'
#' @return A named list of exogenous paths (a subset of \code{model$exogenous};
#'   absent entries are treated as zero downstream).
#' @keywords internal
.hank_check_dZ <- function(model, dZ) {
  exo <- model$exogenous
  T_h <- model$T_h
  if (is.null(dZ)) return(list())
  if (is.numeric(dZ) && is.null(names(dZ))) {
    if (length(exo) != 1L)
      stop("dZ: a bare numeric vector is only accepted when the model has ",
           "exactly one exogenous; this model has ", length(exo), " (",
           paste0("'", exo, "'", collapse = ", "), "). Supply a named list, ",
           "e.g. dZ = list(", exo[1L], " = <length-", T_h, " path>).")
    dZ <- stats::setNames(list(as.numeric(dZ)), exo)
  }
  if (!is.list(dZ))
    stop("dZ must be a named list of exogenous paths (one per: ",
         paste0("'", exo, "'", collapse = ", "), "), or -- for a ",
         "single-exogenous model -- a bare numeric vector. Got ", class(dZ)[1L], ".")
  if (length(dZ) && is.null(names(dZ)))
    stop("dZ is an unnamed list; name its entries after the model's ",
         "exogenous (", paste0("'", exo, "'", collapse = ", "), ").")
  bad <- setdiff(names(dZ), exo)
  if (length(bad))
    stop("dZ name(s) ", paste0("'", bad, "'", collapse = ", "),
         " are not exogenous in this model (exogenous: ",
         paste0("'", exo, "'", collapse = ", "), ").")
  short <- names(dZ)[vapply(dZ, length, integer(1)) != T_h]
  if (length(short))
    stop("dZ entr", if (length(short) > 1L) "ies " else "y ",
         paste0("'", short, "'", collapse = ", "),
         " must have length T_h = ", T_h, ".")
  dZ
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
#'   one per \code{model$exogenous} (missing entries treated as zero). A bare
#'   numeric vector is accepted only when the model has exactly one exogenous.
#'
#' @return Named list \code{dsrc}: one length-\code{T_h} deviation path per
#'   entry of \code{model$unknowns} and \code{model$exogenous}.
#' @keywords internal
.hank_irf_dsrc <- function(model, dZ) {
  T_h <- model$T_h
  dZ <- .hank_check_dZ(model, dZ)
  z_stack <- do.call(c, lapply(model$exogenous, function(z) {
    v <- dZ[[z]]; if (is.null(v)) rep(0, T_h) else v
  }))
  ## Well-posedness guard: a near-singular H_U silently yields a garbage IRF.
  ## Both this rcond and the solve below come from one cached factorization
  ## (see .hank_ge_factor) -- the same H_U is hit once per shock, once per
  ## gradient tap, and once per Newton iteration.
  fac <- .hank_ge_factor(model$H_U)
  rc  <- fac$rcond
  if (!is.finite(rc) || rc < 1e-10)
    warning(sprintf(paste0("hank_model_irf(): H_U is ill-conditioned ",
                           "(rcond = %.2e); the GE solution may be unreliable ",
                           "(near-singular / indeterminate). See hank_determinacy()."),
                    rc))
  dU_stack <- as.numeric(-.hank_ge_solve(fac, model$H_U,
                                         model$H_Z %*% z_stack))

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
#'   one per \code{model$exogenous} (missing entries treated as zero). For a
#'   model with exactly one exogenous a bare numeric vector is also accepted;
#'   with more than one it is an error, since the intent would be ambiguous.
#'   Note the contrast with \code{\link{hank_ks_linear_irf}}, which takes a
#'   bare vector and returns \code{d}-prefixed names.
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
#' the full cross-sectional distribution response of every \code{"het"},
#' \code{"het2"}, \code{"het3"} (and \code{"het_mixture"}) block in the model,
#' via that block's distribution Jacobian (\code{\link{hank_het_dist_jacobian}}
#' for \code{"het"}, \code{\link{hank_het2_dist_jacobian}} for \code{"het2"},
#' \code{\link{hank_het3_dist_jacobian}} for \code{"het3"}, or for a mixture
#' block \code{\link{hank_mixture_dist_jacobian}} -- the omega-weighted sum of
#' the per-type distribution Jacobians on the types' shared grid).
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
#'   one per \code{model$exogenous} (missing entries treated as zero). For a
#'   model with exactly one exogenous a bare numeric vector is also accepted;
#'   with more than one it is an error, since the intent would be ambiguous.
#'   Note the contrast with \code{\link{hank_ks_linear_irf}}, which takes a
#'   bare vector and returns \code{d}-prefixed names.
#' @param ... Forwarded to \code{\link{hank_het2_dist_jacobian}} for any
#'   \code{"het2"}-kind block (e.g. \code{delta_in}, \code{delta_va},
#'   \code{delta_d}, \code{backend}, \code{threads}). Without this, the het2
#'   arm always used the factory's default finite-difference deltas, with no
#'   way to tighten them -- e.g. a \code{theta_coll} dist-IRF was stuck at
#'   only ~1e-3 accuracy (adversarial review, one-liner #2). Ignored by the
#'   \code{"het"}/\code{"het3"}/\code{"het_mixture"} arms, which have no
#'   tunable FD deltas exposed here.
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
hank_model_dist_irf <- function(model, dZ, ...) {
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
    } else if (identical(blk$kind, "het2")) {
      ## Like "het3" (single block, no pooling-across-types ambiguity): every
      ## het2 block has its own well-defined joint (e, b, a) cell space.
      ## `...` forwards FD-delta/backend/threads tuning (one-liner #2,
      ## adversarial review): without it this arm was pinned to the
      ## factory's default deltas with no way to tighten accuracy.
      JD <- hank_het2_dist_jacobian(blk$block, T_h, inputs = blk$inputs, ...)
      n_cell <- length(blk$block$D)
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
    } else if (identical(blk$kind, "het3")) {
      ## Unlike het_mixture, a "het3" block is a SINGLE block (no pooling
      ## across heterogeneous types sharing a grid), so there is no
      ## same_income-style ambiguity to guard against here -- every het3
      ## block has its own well-defined joint (e, d, f, a) cell space.
      JD <- hank_het3_dist_jacobian(blk$block, T_h, inputs = blk$inputs)
      n_cell <- length(blk$block$D)
      dD[[blk$name]] <- contract_dist(JD, blk$inputs, n_cell)
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
