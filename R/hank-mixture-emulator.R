## R/hank-mixture-emulator.R
## --------------------------------------------------------------------------
## D-agnostic emulator of the mixture-HANK GE/Jacobian solve, plus a
## simulation-based calibration (SBC) harness that uses it.
##
## This promotes the Wave-3 het-preferences-in-HANK research scratch
## (`.claude/orchestration/wave3-sbc/`) into tested package API. The bespoke
## scratch hard-coded a 2-D (centre, spread) tensor grid with bilinear
## interpolation (`wave3_lib.R`'s `.bw`/`emu_slice`/`build_proj_grid`/
## `emu_sigma`). Here the same idea is generalized to an arbitrary set of
## GE-coupled coordinates (`d = length(box)`, currently 2 for (centre, spread)
## but written so a 3rd coordinate like `omega1` can be added without a
## rewrite): a MULTILINEAR interpolant on a tensor grid (the direct
## `d`-dimensional generalization of `.bw`), and a THIN-PLATE-SPLINE RBF
## interpolant on ANY design (grid or scattered Latin-hypercube), for the
## higher-dimensional regime where a full tensor grid is infeasible.
##
## Cost structure this emulator exploits: each design node requires one
## mixture GE steady-state solve (`hank_mixture_ks_steady`, ~0.4s at the small
## test config) plus a distribution Jacobian (`hank_mixture_dist_jacobian`).
## Both are DATA-INDEPENDENT -- they depend only on theta, not on any observed
## sample -- so they are precomputed ONCE per design node and reused across
## every likelihood evaluation of every downstream SBC replication.
##
## The raw per-cell distribution response `dD_snap` does NOT interpolate
## cleanly (near-constraint structure shifts between theta-nodes make it
## spiky; see R/hank-reweighting.R's file section comment on the grid-
## invariant functional basis). But its PROJECTION onto a fixed reference
## basis `B_ref` (built once, at the design centre) is a low-order smooth
## functional of theta and interpolates well -- exactly the
## `hank_reweight_functional_basis`/`build_proj_grid` insight from the
## research scratch, reproduced here as package code.
## --------------------------------------------------------------------------


## ---- small internal helpers -----------------------------------------------

#' Bracket-and-weight a single coordinate against a sorted 1-D grid
#'
#' Direct generalization of the research scratch's `.bw`: clamps queries
#' outside the grid to the boundary node (weight 1 on that node), otherwise
#' returns the bracketing pair of grid indices and the linear weight on the
#' LOWER index.
#' @keywords internal
.hank_emu_bw <- function(grid, x) {
  n <- length(grid)
  if (x <= grid[1L]) return(list(i0 = 1L, i1 = 1L, w = 1))
  if (x >= grid[n]) return(list(i0 = n, i1 = n, w = 1))
  i1 <- findInterval(x, grid)
  i0 <- i1
  i1 <- i0 + 1L
  w <- (grid[i1] - x) / (grid[i1] - grid[i0])    # weight on i0
  list(i0 = i0, i1 = i1, w = w)
}


#' `d`-dimensional multilinear interpolation of a stored per-node array
#'
#' Generalizes `.bw` + `emu_slice` (bilinear, `d = 2`) to an arbitrary number
#' of GE-coupled coordinates: sums over the `2^d` corners of the bracketing
#' hyper-rectangle with the product of per-axis weights.
#'
#' @param axis_grids Named list, one sorted numeric grid vector per coordinate
#'   (in `names(box)` order).
#' @param A An array whose first `d = length(axis_grids)` dimensions index the
#'   tensor-grid nodes (in the SAME axis order as `axis_grids`) and whose
#'   remaining ("tail") dimensions are the stored field's own shape (e.g. a
#'   `T_h x n_obs` `Theta` matrix, or a length-`K` vector).
#' @param theta_vec Named numeric length-`d` query point (in `axis_grids`
#'   order; `names(theta_vec)` are not required, only position is used, but
#'   callers pass named vectors for clarity).
#'
#' @return The interpolated tail-shaped value: a numeric vector if `A` has one
#'   tail dimension, else an array of the tail dimensions.
#' @keywords internal
.hank_emu_multilinear <- function(axis_grids, A, theta_vec) {
  d <- length(axis_grids)
  dm <- dim(A)
  tail_dim <- dm[-seq_len(d)]
  m <- prod(tail_dim)
  Af <- A
  dim(Af) <- c(dm[seq_len(d)], m)

  bw <- lapply(seq_len(d), function(k) .hank_emu_bw(axis_grids[[k]], theta_vec[k]))

  ## Sum over the 2^d corners: corner b in 0..(2^d - 1), bit k selects i0/i1
  ## on axis k, weight = prod_k (bw[[k]]$w if bit==0 else 1 - bw[[k]]$w).
  acc <- NULL
  n_corner <- 2L^d
  for (b in 0:(n_corner - 1L)) {
    bits <- bitwAnd(bitwShiftR(b, seq_len(d) - 1L), 1L)   # 0 = i0, 1 = i1
    idx  <- vector("list", d)
    wgt  <- 1
    for (k in seq_len(d)) {
      idx[[k]] <- if (bits[k] == 0L) bw[[k]]$i0 else bw[[k]]$i1
      wgt <- wgt * (if (bits[k] == 0L) bw[[k]]$w else (1 - bw[[k]]$w))
    }
    if (wgt == 0) next
    corner_val <- do.call(`[`, c(list(Af), idx, list(TRUE)))
    acc <- if (is.null(acc)) wgt * corner_val else acc + wgt * corner_val
  }
  if (length(tail_dim) > 1L) array(acc, tail_dim) else as.numeric(acc)
}


#' Simple base-R Latin-hypercube design in `[0, 1]^d`
#'
#' A classic stratified LHS: each coordinate's `[0, 1]` range is split into
#' `n` equal strata, one point drawn uniformly within each stratum, and the
#' per-coordinate strata orderings independently permuted (so every 1-D
#' marginal is a jittered equal-stratification and no package dependency is
#' needed for this one design step).
#'
#' @param n Integer number of design points.
#' @param d Integer number of coordinates.
#' @return An `n x d` matrix in `[0, 1]^d`.
#' @keywords internal
.hank_emu_lhs <- function(n, d) {
  n <- as.integer(n)
  X <- matrix(0, n, d)
  for (k in seq_len(d)) {
    perm <- sample.int(n)
    X[, k] <- (perm - 1 + stats::runif(n)) / n
  }
  X
}


#' Thin-plate-spline RBF radial basis function `phi(r) = r^2 log(r)`, `phi(0) = 0`
#' @keywords internal
.hank_emu_tps <- function(r) {
  out <- numeric(length(r))
  pos <- r > 0
  out[pos] <- r[pos]^2 * log(r[pos])
  out
}


#' Fit a thin-plate-spline RBF interpolant (with linear tail) on normalized nodes
#'
#' Solves the saddle-point system once at BUILD time; the fitted `(W, C)`
#' weight matrices are data-independent of any later query.
#'
#' @param Xn `n x d` matrix of NORMALIZED (box-mapped to `[0,1]`) node
#'   coordinates.
#' @param Y `n x m` matrix of stacked node values (`m` = flattened field
#'   length).
#' @return A list with `W` (`n x m`), `C` (`(d+1) x m`), `Xn` (echoed back for
#'   prediction), and `ridge` (the diagonal ridge actually added, `0` if the
#'   system solved without one).
#' @keywords internal
.hank_emu_rbf_fit <- function(Xn, Y) {
  n <- nrow(Xn); d <- ncol(Xn)
  Dmat <- as.matrix(stats::dist(Xn))
  Phi <- .hank_emu_tps(Dmat)
  P <- cbind(1, Xn)                          # n x (d+1)
  Z0 <- matrix(0, d + 1L, d + 1L)
  Yz <- matrix(0, d + 1L, ncol(Y))

  solve_once <- function(ridge) {
    M <- rbind(cbind(Phi + ridge * diag(n), P), cbind(t(P), Z0))
    RHS <- rbind(Y, Yz)
    solve(M, RHS)
  }
  ridge <- 0
  WC <- tryCatch(solve_once(0), error = function(e) NULL)
  if (is.null(WC)) {
    ridge <- 1e-10
    WC <- solve_once(ridge)
  }
  list(W = WC[seq_len(n), , drop = FALSE],
       C = WC[(n + 1L):(n + d + 1L), , drop = FALSE],
       Xn = Xn, ridge = ridge)
}


#' Predict a fitted thin-plate-spline RBF interpolant at normalized query points
#' @param fit A `.hank_emu_rbf_fit()` return value.
#' @param xq Numeric length-`d` normalized query point.
#' @return Numeric length-`m` predicted (flattened) field value.
#' @keywords internal
.hank_emu_rbf_predict <- function(fit, xq) {
  dvec <- sqrt(rowSums(sweep(fit$Xn, 2, xq, `-`)^2))
  phi_q <- .hank_emu_tps(dvec)
  as.numeric(crossprod(phi_q, fit$W) + crossprod(c(1, xq), fit$C))
}


## ---- rebuild a KF state space from an interpolated Theta (single shock "Z")

#' Rebuild a `dsge_ss` state space from an interpolated `Theta` (`q x n_obs`)
#'
#' Reproduces `hank_state_space`'s single-shock ("Z"), lagged-timing,
#' block-shift-register construction EXACTLY, copied from the research
#' scratch's `mk_ss` (`wave3_lib.R`) -- the interpolated `Theta` plays the role
#' the direct GE solve's `Theta_list[["Z"]]` would.
#' @keywords internal
.hank_emu_mk_ss <- function(emu, Theta_mat) {
  cfg <- emu$config
  q <- nrow(Theta_mat); n_obs <- ncol(Theta_mat)
  sigma <- cfg$shock_specs$Z$sigma
  TT <- matrix(0, q, q)
  if (q >= 2L) TT[cbind(2:q, 1:(q - 1L))] <- 1
  RR <- matrix(0, q, 1L); RR[1L, 1L] <- 1
  Z_lag <- matrix(0, n_obs, q)
  if (q >= 2L) for (s in 1:(q - 1L)) Z_lag[, s] <- Theta_mat[s + 1L, ]
  D_lag <- matrix(Theta_mat[1L, ], n_obs, 1L)
  new_dsge_ss(T_mat = TT, R_mat = RR, Z_mat = Z_lag, D_mat = D_lag,
                      Sigma_e = matrix(sigma^2, 1L, 1L),
                      state_names = paste0("Z_lag", seq_len(q) - 1L),
                      obs_names = cfg$observables, shock_names = "Z",
                      timing = "lagged",
                      Theta_list = list(Z = Theta_mat), q = q)
}


## --------------------------------------------------------------------------
## PART 1 -- the emulator
## --------------------------------------------------------------------------


#' Build a D-agnostic emulator of the mixture-HANK GE/Jacobian solve
#'
#' Precomputes the data-independent general-equilibrium/Jacobian pieces of a
#' 2-type discount-factor mixture (steady state, macro state-space MA
#' coefficients, distribution-reweighting response) on a DESIGN of theta
#' points covering `box`, and returns an interpolator object. The design is
#' \strong{D-agnostic}: `box` is a named list of length-2 numeric ranges, one
#' per GE-coupled coordinate, and every piece of interpolation/design code
#' below is written for a general `d = length(box)` (currently the natural
#' choices are a subset of `c("centre", "spread", "omega1")`, but nothing here
#' hard-codes `d = 2`).
#'
#' \strong{theta -> economics map}: `betas = c(centre - spread, centre +
#' spread)`; `omega = c(theta$omega1, 1 - theta$omega1)` if `"omega1" %in%
#' names(box)`, else the fixed `config$omega` (default `c(0.5, 0.5)`).
#'
#' \strong{Per-node solve} (mirrors the research scratch's `snap_precompute.R`
#' `one()`): steady state (\code{\link{hank_mixture_ks_steady}}) -> GE model
#' (\code{\link{hank_mixture_ks_model}}) -> truncated-MA state space
#' (\code{\link{hank_state_space}}, whose `Theta_list[["Z"]]` is stored) ->
#' per-type blocks (\code{\link{hank_mixture_blocks}}) -> stationary
#' distribution `D0` (\code{\link{hank_mixture_dist}}) -> distribution
#' Jacobian `JD` (\code{\link{hank_mixture_dist_jacobian}}) -> single-date
#' snapshot reweighting response `dD_snap`
#' (\code{\link{hank_dist_response_snapshot}}).
#'
#' \strong{Projected reweighting pieces}: the raw per-cell `dD_snap` does not
#' interpolate cleanly (see the file section comment in R/hank-reweighting.R),
#' so only its projection onto a FIXED reference basis `B_ref`
#' (\code{\link{hank_reweight_functional_basis}}, built once at the node
#' nearest the box centre) is stored and interpolated: `gref = t(B_ref) \%*\%
#' dD_snap` (length `K`), `m1 = t(B_ref) \%*\% D0` (length `K`), `M2D0 =
#' t(B_ref) \%*\% (D0 * B_ref)` (`K x K`), `M2dD = t(B_ref) \%*\% (dD_snap *
#' B_ref)` (`K x K`) -- exactly the research scratch's `build_proj_grid`.
#'
#' @param box Named list of length-2 numeric `c(lo, hi)` ranges, one per
#'   GE-coupled coordinate, e.g. \code{list(centre = c(0.920, 0.940), spread =
#'   c(0.008, 0.016))}. `d = length(box)` is arbitrary.
#' @param config Fixed HANK configuration list with: `n_e, n_a, amax, alpha,
#'   delta, Z, eis, T_h` (household/GE sizing and calibration); `shock_specs`
#'   (named list, one `list(rho, sigma)` per exogenous macro shock, e.g.
#'   \code{list(Z = list(rho = 0.9, sigma = 0.007))}); `observables` (character
#'   vector, e.g. \code{c("C", "K")}); `rw_shock` (named list, one
#'   length-`T_h` reweighting shock path per input, e.g. \code{list(r = 0.01 *
#'   0.8^(seq_len(T_h) - 1))}); `t_star` (integer output date for the snapshot
#'   reweighting); `N_ref` (reference survey size, used only to document the
#'   design's intended scale -- the SBC harness passes its own `N`); `degree`
#'   (basis degree, default 8); `tail_trim` (default 0); `omega` (default
#'   `c(0.5, 0.5)`, used whenever `"omega1"` is not a `box` coordinate).
#' @param design `"grid"` (default surrogate `"multilinear"`): a tensor
#'   product of `n_grid` equally-spaced points per coordinate (`n_grid^d`
#'   nodes total). `"scatter"` (default surrogate `"rbf"`): `n_scatter`
#'   Latin-hypercube points in the box.
#' @param n_grid Integer points per axis for `design = "grid"`.
#' @param n_scatter Integer total design points for `design = "scatter"`
#'   (required in that case).
#' @param surrogate `"multilinear"` or `"rbf"`; `NULL` (default) picks
#'   `"multilinear"` for `design = "grid"` and `"rbf"` for `design =
#'   "scatter"`. `"multilinear"` REQUIRES `design = "grid"`.
#' @param seed RNG seed for the `"scatter"` Latin-hypercube design.
#'
#' @return An object of class `"hank_mixture_emulator"`: a list with `box`,
#'   `config`, `design` (`"grid"`/`"scatter"`), `surrogate`, `nodes` (`n x d`
#'   matrix of node coordinates in natural units, `colnames = names(box)`),
#'   `axis_grids` (named list of per-axis grids, only set for `design =
#'   "grid"`), `B_ref` (`n_cell x degree`), `a_cell`, `K` (`= ncol(B_ref)`),
#'   `Theta` (array: leading node dimension(s) + `T_h x n_obs` tail),
#'   `gref`, `m1` (each: leading node dimension(s) + length-`K` tail), `M2D0`,
#'   `M2dD` (each: leading node dimension(s) + `K x K` tail), `r`, `w` (each:
#'   an array over nodes), and (for `surrogate = "rbf"`) `rbf_fits` (named
#'   list of `.hank_emu_rbf_fit()` results, one per stored field).
#' @export
hank_mixture_emulator <- function(box, config, design = c("grid", "scatter"),
                                   n_grid = 5L, n_scatter = NULL,
                                   surrogate = NULL, seed = 1L) {
  design <- match.arg(design)
  nm <- names(box)
  d  <- length(box)
  if (d < 1L) stop("hank_mixture_emulator(): 'box' must have at least one coordinate.")
  if (is.null(surrogate)) surrogate <- if (design == "grid") "multilinear" else "rbf"
  if (surrogate == "multilinear" && design != "grid")
    stop("hank_mixture_emulator(): surrogate = 'multilinear' requires design = 'grid'.")
  if (design == "scatter" && is.null(n_scatter))
    stop("hank_mixture_emulator(): design = 'scatter' requires 'n_scatter'.")

  if (is.null(config$omega)) config$omega <- c(0.5, 0.5)
  if (is.null(config$degree)) config$degree <- 8L
  if (is.null(config$tail_trim)) config$tail_trim <- 0

  ## ---- design nodes (natural units) --------------------------------------
  axis_grids <- NULL
  if (design == "grid") {
    axis_grids <- setNames(lapply(nm, function(n) seq(box[[n]][1L], box[[n]][2L],
                                                        length.out = n_grid)), nm)
    grid_idx <- do.call(expand.grid, c(lapply(axis_grids, seq_along), KEEP.OUT.ATTRS = FALSE))
    nodes <- do.call(cbind, lapply(seq_len(d), function(k) axis_grids[[nm[k]]][grid_idx[[k]]]))
    colnames(nodes) <- nm
    node_dim <- unname(vapply(axis_grids, length, integer(1)))
  } else {
    set.seed(seed)
    unit <- .hank_emu_lhs(n_scatter, d)
    nodes <- do.call(cbind, lapply(seq_len(d), function(k) {
      rng <- box[[nm[k]]]; rng[1L] + unit[, k] * (rng[2L] - rng[1L])
    }))
    colnames(nodes) <- nm
    node_dim <- n_scatter
  }
  n_nodes <- nrow(nodes)

  ## ---- shared household primitives (fixed across nodes) ------------------
  inc <- hank_income_rouwenhorst(rho = 0.9, sigma = 0.7, n = config$n_e)
  ag  <- hank_asset_grid(amax = config$amax, n = config$n_a, amin = 0)
  a_cell <- rep(ag, times = config$n_e)             # asset-fast cell order

  ## ---- per-node GE solve (mirrors snap_precompute.R's one()) --------------
  one <- function(theta_row) {
    centre <- theta_row[["centre"]]; spread <- theta_row[["spread"]]
    betas <- c(centre - spread, centre + spread)
    omega <- if ("omega1" %in% nm) c(theta_row[["omega1"]], 1 - theta_row[["omega1"]]) else config$omega

    mks <- hank_mixture_ks_steady(ag, inc$Pi, inc$e, betas, omega, eis = config$eis,
                                  alpha = config$alpha, delta = config$delta, Z = config$Z)
    model  <- hank_mixture_ks_model(mks, config$T_h)
    sspace <- hank_state_space(model, config$shock_specs, config$observables, q = NULL)
    b   <- hank_mixture_blocks(ag, inc$Pi, inc$e, betas = betas, eis = config$eis,
                               r = mks$r, w = mks$w)
    D0  <- hank_mixture_dist(b, omega)$D
    JD  <- hank_mixture_dist_jacobian(b, omega, config$T_h, inputs = "r")
    dD_snap <- hank_dist_response_snapshot(JD, config$rw_shock, config$t_star)
    list(Theta = sspace$Theta_list[["Z"]], r = mks$r, w = mks$w, D0 = D0, dD_snap = dD_snap)
  }

  results <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) results[[i]] <- one(nodes[i, ])

  ## ---- reference basis B_ref at the node nearest the box centre -----------
  centre_target <- vapply(nm, function(n) mean(box[[n]]), numeric(1))
  d2 <- rowSums(sweep(nodes, 2, centre_target[nm], `-`)^2)
  i_centre <- which.min(d2)
  D0_centre <- results[[i_centre]]$D0
  B_ref <- hank_reweight_functional_basis(a_cell, D0_centre, degree = config$degree,
                                          tail_trim = config$tail_trim)
  K <- ncol(B_ref)

  ## ---- assemble per-node arrays (Theta, r, w, gref, m1, M2D0, M2dD) -------
  T_h_reg <- nrow(results[[1L]]$Theta); n_obs <- ncol(results[[1L]]$Theta)
  Theta_arr <- array(NA_real_, c(node_dim, T_h_reg, n_obs))
  gref_arr  <- array(NA_real_, c(node_dim, K))
  m1_arr    <- array(NA_real_, c(node_dim, K))
  M2D0_arr  <- array(NA_real_, c(node_dim, K, K))
  M2dD_arr  <- array(NA_real_, c(node_dim, K, K))
  r_arr <- array(NA_real_, node_dim); w_arr <- array(NA_real_, node_dim)

  ## flat-index helper: node i maps to the leading multi-index implied by
  ## node_dim (grid: the d-dim tensor index; scatter: a plain 1-D index).
  is_grid <- design == "grid"
  for (i in seq_len(n_nodes)) {
    res <- results[[i]]
    leading <- if (is_grid) as.list(grid_idx[i, ]) else list(i)
    Theta_arr <- do.call(`[<-`, c(list(Theta_arr), leading, list(TRUE, TRUE), list(value = res$Theta)))
    gref_val <- as.numeric(crossprod(B_ref, res$dD_snap))
    m1_val   <- as.numeric(crossprod(B_ref, res$D0))
    M2D0_val <- crossprod(B_ref, res$D0 * B_ref)
    M2dD_val <- crossprod(B_ref, res$dD_snap * B_ref)
    gref_arr <- do.call(`[<-`, c(list(gref_arr), leading, list(TRUE), list(value = gref_val)))
    m1_arr   <- do.call(`[<-`, c(list(m1_arr), leading, list(TRUE), list(value = m1_val)))
    M2D0_arr <- do.call(`[<-`, c(list(M2D0_arr), leading, list(TRUE, TRUE), list(value = M2D0_val)))
    M2dD_arr <- do.call(`[<-`, c(list(M2dD_arr), leading, list(TRUE, TRUE), list(value = M2dD_val)))
    r_arr <- do.call(`[<-`, c(list(r_arr), leading, list(value = res$r)))
    w_arr <- do.call(`[<-`, c(list(w_arr), leading, list(value = res$w)))
  }

  emu <- list(box = box, config = config, design = design, surrogate = surrogate,
              nodes = nodes, axis_grids = axis_grids, node_dim = node_dim,
              B_ref = B_ref, a_cell = a_cell, K = K,
              Theta = Theta_arr, gref = gref_arr, m1 = m1_arr,
              M2D0 = M2D0_arr, M2dD = M2dD_arr, r = r_arr, w = w_arr)

  if (surrogate == "rbf") {
    ## Normalize nodes to [0,1]^d and fit one RBF interpolant per stored field
    ## (each field flattened to n_nodes x m for .hank_emu_rbf_fit).
    Xn <- do.call(cbind, lapply(nm, function(n) {
      rng <- box[[n]]; (nodes[, n] - rng[1L]) / (rng[2L] - rng[1L])
    }))
    flatten_field <- function(A) {
      dm <- dim(A); tail_m <- prod(dm[-1L])
      matrix(A, n_nodes, tail_m)
    }
    emu$rbf_fits <- list(
      Theta = .hank_emu_rbf_fit(Xn, flatten_field(Theta_arr)),
      gref  = .hank_emu_rbf_fit(Xn, flatten_field(gref_arr)),
      m1    = .hank_emu_rbf_fit(Xn, flatten_field(m1_arr)),
      M2D0  = .hank_emu_rbf_fit(Xn, flatten_field(M2D0_arr)),
      M2dD  = .hank_emu_rbf_fit(Xn, flatten_field(M2dD_arr)))
    emu$rbf_field_dim <- list(Theta = dim(Theta_arr)[-1L], gref = K, m1 = K,
                              M2D0 = c(K, K), M2dD = c(K, K))
  }

  structure(emu, class = "hank_mixture_emulator")
}


#' Interpolate a stored per-node field of a `hank_mixture_emulator` at a query theta
#'
#' Internal dispatcher between the two surrogate backends: `"multilinear"`
#' (tensor-grid `d`-dimensional linear interpolation, generalizing the
#' research scratch's `.bw`/`emu_slice`) and `"rbf"` (thin-plate-spline with a
#' linear tail, on box-normalized coordinates).
#'
#' @param emu A \code{\link{hank_mixture_emulator}}.
#' @param theta_vec Named numeric length-`d` query point (or a list with
#'   `names(emu$box)` entries), in `names(emu$box)` order.
#' @param field One of `"Theta"`, `"gref"`, `"m1"`, `"M2D0"`, `"M2dD"`.
#' @return The interpolated field value, same shape as one node's slice of
#'   `emu[[field]]` (dropping the leading node dimension(s)).
#' @keywords internal
.emu_interp <- function(emu, theta_vec, field) {
  nm <- names(emu$box)
  theta_vec <- if (is.list(theta_vec)) vapply(nm, function(n) theta_vec[[n]], numeric(1)) else theta_vec[nm]

  if (emu$surrogate == "multilinear") {
    return(.hank_emu_multilinear(emu$axis_grids, emu[[field]], theta_vec))
  }

  ## rbf: normalize, predict the flattened field, reshape to its natural shape.
  x01 <- vapply(nm, function(n) {
    rng <- emu$box[[n]]; (theta_vec[[n]] - rng[1L]) / (rng[2L] - rng[1L])
  }, numeric(1))
  flat <- .hank_emu_rbf_predict(emu$rbf_fits[[field]], x01)
  fd <- emu$rbf_field_dim[[field]]
  if (length(fd) > 1L) matrix(flat, fd[1L], fd[2L]) else as.numeric(flat)
}


#' Emulated macro state space at a query theta
#'
#' Interpolates the emulator's stored `Theta` field at `theta` and rebuilds a
#' `dsge_ss` state space from it (\code{\link{.hank_emu_mk_ss}}, reproducing
#' \code{\link{hank_state_space}}'s single-shock shift-register construction
#' EXACTLY).
#'
#' @param emu A \code{\link{hank_mixture_emulator}}.
#' @param theta Named list or length-`d` named numeric in `names(emu$box)`
#'   order.
#' @return A `dsge_ss` object (see \code{\link{hank_kalman_loglik}}).
#' @export
hank_emulator_state_space <- function(emu, theta) {
  Theta_mat <- .emu_interp(emu, theta, "Theta")
  .hank_emu_mk_ss(emu, Theta_mat)
}


#' Emulated B_ref-projected reweighting mean at a query theta
#'
#' `rw_scale * ` the interpolated `gref` field: the B_ref-projected
#' model-implied reweighting mean at `theta`, for reweighting shock size
#' `rw_scale` (the emulator itself stores `gref` at the UNIT `config$rw_shock`
#' path; `rw_scale` rescales it, matching the research scratch's
#' `RW_SCALE * emu_slice(G, GP$gref, cc, ss)` usage).
#'
#' @inheritParams hank_emulator_state_space
#' @param rw_scale Scalar multiplier on the reweighting shock size.
#' @return Length-`emu$K` numeric vector.
#' @export
hank_emulator_reweight_mean <- function(emu, theta, rw_scale) {
  rw_scale * .emu_interp(emu, theta, "gref")
}


#' Emulated `K x K` covariance of the projected net reweighting at a query theta
#'
#' Rebuilds the covariance of `t(B_ref) (Dhat1 - Dhat0)` (two independent
#' size-`N` multinomial cross-sections, pre- and post-shock) from the
#' interpolated projected pieces `m1, M2D0, M2dD`, EXACTLY as the research
#' scratch's `emu_sigma`: with `D1 = D0 + rw_scale * dD_snap`,
#' \code{covD0 = M2D0 - m1 m1^T}, \code{covD1 = (M2D0 + rw_scale * M2dD) -
#' m1_D1 m1_D1^T} (\code{m1_D1 = m1 + rw_scale * gref}), \code{Sigma = (covD0 +
#' covD1) / N}.
#'
#' @inheritParams hank_emulator_reweight_mean
#' @param N Sample size of EACH of the two cross-sections.
#' @return A `emu$K x emu$K` numeric covariance matrix.
#' @export
hank_emulator_metric <- function(emu, theta, rw_scale, N) {
  m1   <- .emu_interp(emu, theta, "m1")
  gref <- .emu_interp(emu, theta, "gref")
  M2D0 <- .emu_interp(emu, theta, "M2D0")
  M2dD <- .emu_interp(emu, theta, "M2dD")
  m1_D1 <- m1 + rw_scale * gref
  covD0 <- M2D0 - tcrossprod(m1)
  covD1 <- (M2D0 + rw_scale * M2dD) - tcrossprod(m1_D1)
  (covD0 + covD1) / N
}


#' Emulated B_ref-projected stationary-distribution (LEVEL) mean at a query theta
#'
#' The interpolated `m1` field: the B_ref-projection `t(B_ref) D0(theta)` of the
#' candidate mixture's STATIONARY wealth distribution. This is the model mean of
#' the LEVEL cross-section observable -- a single stationary wealth survey, as
#' opposed to the price-shock reweighting RESPONSE
#' (\code{\link{hank_emulator_reweight_mean}}). The stationary distribution
#' SHAPE is far more sensitive to the discount-factor spread than the small
#' price-shock response is (the response differences most of that level signal
#' away), so this level observable is the sharp identifier of the mixture spread
#' -- see \code{\link{hank_mixture_sbc}}'s \code{channels} argument and Details.
#'
#' @inheritParams hank_emulator_state_space
#' @return Length-`emu$K` numeric vector.
#' @seealso \code{\link{hank_emulator_level_metric}}, \code{\link{hank_emulator_reweight_mean}}
#' @export
hank_emulator_level_mean <- function(emu, theta) {
  .emu_interp(emu, theta, "m1")
}


#' Emulated `K x K` covariance of a projected stationary (LEVEL) wealth survey
#'
#' The covariance of `t(B_ref) Dhat` for a SINGLE size-`N` multinomial
#' cross-section `Dhat` drawn from the stationary distribution `D0(theta)`,
#' rebuilt from the interpolated projected pieces `m1, M2D0`:
#' \deqn{\Sigma_L = cov_{D0}(B\_ref) / N = (M2D0 - m1\, m1^T) / N.}
#' This is the LEVEL-observable analogue of \code{\link{hank_emulator_metric}}
#' (which is the RESPONSE observable's covariance of a DIFFERENCE of two
#' cross-sections); the level uses only one survey, hence the single
#' `cov_D0` term.
#'
#' @inheritParams hank_emulator_state_space
#' @param N Sample size of the stationary cross-section.
#' @return A `emu$K x emu$K` numeric covariance matrix.
#' @seealso \code{\link{hank_emulator_level_mean}}, \code{\link{hank_emulator_metric}}
#' @export
hank_emulator_level_metric <- function(emu, theta, N) {
  m1   <- .emu_interp(emu, theta, "m1")
  M2D0 <- .emu_interp(emu, theta, "M2D0")
  (M2D0 - tcrossprod(m1)) / N
}


## --------------------------------------------------------------------------
## PART 2 -- the SBC harness
## --------------------------------------------------------------------------


#' Direct (non-emulated) GE-solve DGP truth for the SBC harness
#'
#' Mirrors the research scratch's `snap_sbc.R` `direct_truth()`: a fresh
#' mixture GE steady-state solve at `theta*` (NOT read off the emulator, so
#' the emulator does not mark its own homework), the truth stationary
#' distribution `D0`, and the GENUINE NONLINEAR post-shock snapshot
#' cross-section `D1` at `t_star` (\code{\link{hank_mixture_td_nonlinear}}, no
#' linear-response clipping).
#'
#' @param theta_star Named list with `centre`, `spread` (and `omega1` if it is
#'   a `box` coordinate).
#' @param emu A \code{\link{hank_mixture_emulator}} (used only for its `box`
#'   coordinate names, `config`, and shared household primitives via
#'   `emu$a_cell`/`emu$config`).
#' @param inc,ag Shared income process / asset grid (as built inside
#'   \code{\link{hank_mixture_emulator}}; passed in so the SBC harness builds
#'   them once, not once per replication).
#' @param rw_scale Reweighting shock-size multiplier: the direct-truth price
#'   path is `r_path = r + rw_scale * config$rw_shock$r` (passed explicitly
#'   rather than read off `emu$config`, so a caller's `rw_scale` argument to
#'   \code{\link{hank_mixture_sbc}} cannot silently diverge from the DGP's own
#'   shock size).
#' @param need_D1 Logical: whether to compute the nonlinear post-shock snapshot
#'   `D1` (needed only when the `"response"` channel is active). When `FALSE`
#'   (e.g. a macro + level-only SBC) the extra \code{\link{hank_mixture_td_nonlinear}}
#'   solve is skipped and `D1` is `NULL`.
#' @return A list with `ss_dir` (`dsge_ss`), `D0`, `D1` (length-`n_cell` or
#'   `NULL` when `need_D1 = FALSE`), `r`, `w`.
#' @keywords internal
.hank_sbc_direct_truth <- function(theta_star, emu, inc, ag, rw_scale, need_D1 = TRUE) {
  cfg <- emu$config
  nm  <- names(emu$box)
  centre <- theta_star$centre; spread <- theta_star$spread
  betas <- c(centre - spread, centre + spread)
  omega <- if ("omega1" %in% nm) c(theta_star$omega1, 1 - theta_star$omega1) else cfg$omega

  mks <- hank_mixture_ks_steady(ag, inc$Pi, inc$e, betas, omega, eis = cfg$eis,
                                alpha = cfg$alpha, delta = cfg$delta, Z = cfg$Z)
  model <- hank_mixture_ks_model(mks, cfg$T_h)
  ss_dir <- hank_state_space(model, cfg$shock_specs, cfg$observables, q = NULL)
  b  <- hank_mixture_blocks(ag, inc$Pi, inc$e, betas = betas, eis = cfg$eis,
                            r = mks$r, w = mks$w)
  D0 <- hank_mixture_dist(b, omega)$D
  D1 <- NULL
  if (need_D1) {
    td <- hank_mixture_td_nonlinear(b, omega,
                                    r_path = mks$r + rw_scale * cfg$rw_shock$r,
                                    w_path = rep(mks$w, cfg$T_h))
    D1 <- td$Dpath[, cfg$t_star]
  }
  list(ss_dir = ss_dir, D0 = D0, D1 = D1, r = mks$r, w = mks$w)
}


#' Simulate a macro series from a `dsge_ss` state space (single-shock MA form)
#'
#' Mirrors the research scratch's `sim_Y` (`wave3_lib.R`): draws iid `N(0,
#' sigma_Z^2)` innovations, forms the truncated-MA observation path, adds iid
#' measurement error, and demeans.
#' @keywords internal
.hank_sbc_sim_Y <- function(ss_obj, T_data, me_var, sigma_Z, observables) {
  Theta <- ss_obj$Theta_list[["Z"]]; q <- nrow(Theta); n_obs <- ncol(Theta)
  eps <- stats::rnorm(T_data, 0, sigma_Z)
  Y <- matrix(0, T_data, n_obs)
  for (t in seq_len(T_data)) {
    s <- 0:min(t - 1L, q - 1L)
    Y[t, ] <- colSums(Theta[s + 1L, , drop = FALSE] * eps[t - s])
  }
  if (me_var > 0) Y <- Y + matrix(stats::rnorm(T_data * n_obs, 0, sqrt(me_var)), T_data, n_obs)
  Y <- scale(Y, center = TRUE, scale = FALSE); attr(Y, "scaled:center") <- NULL
  matrix(Y, T_data, n_obs, dimnames = list(NULL, observables))
}


#' Exact-randomized weighted PIT of a test quantity against a grid pmf
#'
#' The rank/PIT of \code{g_star = g(theta_star)} in the posterior distribution
#' of the test quantity \code{g}, where the posterior over grid nodes is the
#' weight vector \code{w} (a normalized pmf). This is the test-quantity
#' generalization of the per-coordinate marginal PIT: ranking a JOINT functional
#' \code{g(theta)} (log-posterior, or a centre-by-spread interaction) catches
#' calibration errors in the dependence structure that a per-coordinate marginal
#' rank is blind to (Modrak et al. 2023). Ties are split exact-randomly,
#' matching the coordinate PIT.
#' @param g_star Scalar test-quantity value at the truth node.
#' @param g_all Numeric vector of the test quantity over all grid nodes.
#' @param w Numeric posterior weights over the grid nodes (need not be
#'   normalized).
#' @return A single PIT value in \code{[0, 1]}.
#' @keywords internal
.sbc_weighted_pit <- function(g_star, g_all, w) {
  w  <- w / sum(w)
  lt <- sum(w[g_all < g_star])
  eq <- sum(w[g_all == g_star])
  lt + stats::runif(1) * eq
}


#' Simulation-based calibration of the joint emulated estimator
#'
#' Runs simulation-based calibration (posterior rank-uniformity) of the joint
#' emulated mixture-HANK estimator, using `emu` (a
#' \code{\link{hank_mixture_emulator}}) for every POSTERIOR evaluation but the
#' DIRECT (non-emulated) GE solve for every DGP truth (\code{
#' \link{.hank_sbc_direct_truth}}) -- so the emulator is never asked to mark
#' its own homework. Mirrors the research scratch's `snap_sbc.R`.
#'
#' \strong{Observable channels}: `channels` selects which data sources enter
#' the likelihood, any non-empty subset of
#' \describe{
#'   \item{`"macro"`}{the aggregate time series `Y` filtered through the
#'     structural Kalman filter (\code{\link{hank_kalman_loglik}}) -- informative
#'     mainly about the mixture CENTRE.}
#'   \item{`"level"`}{a single stationary wealth cross-section
#'     (\code{\link{hank_emulator_level_mean}}) -- the discount-factor SPREAD
#'     reshapes the stationary distribution strongly, so this is the SHARP
#'     spread identifier (the Krusell-Smith / \dQuote{beta-heterogeneity <->
#'     wealth inequality} channel).}
#'   \item{`"response"`}{the price-shock net reweighting between two dated
#'     cross-sections (\code{\link{hank_emulator_reweight_mean}}) -- ROBUST (it
#'     differences out fixed cross-sectional heterogeneity) but, at a realistic
#'     survey size and shock, it also differences most of the level's spread
#'     signal away and is nearly uninformative about the spread; kept for
#'     completeness and as the original brief-17 observable.}
#' }
#' The default `c("macro", "response")` reproduces the joint estimator of
#' \code{\link{hank_mixture_joint_logpost}}; `c("macro", "level")` is the
#' recommended combination for sharply identifying BOTH the centre (macro) and
#' the spread (level). When `"response"` is not among `channels`, the
#' per-replication nonlinear post-shock snapshot solve is skipped.
#'
#' \strong{Prior / truth draws}: `theta*` is drawn per replication from a
#' DISCRETE uniform over the posterior evaluation grid's own nodes (rather
#' than a continuous uniform over `box`), so that the exact-randomized
#' probability-integral-transform (PIT) used for the rank statistic (step 4
#' below) is EXACT with no additional between-node apportionment -- matching
#' the research scratch's convention of drawing truth on the same fine grid
#' the posterior is evaluated over.
#'
#' \strong{Fixed reference measurement-error variance}: one reference macro
#' series is simulated at the box centre with `set.seed(T_data_seed)`, and
#' `me_var = (me_frac * min(per-column sd))^2` is fixed from it and shared
#' (not truth-dependent) across every replication.
#'
#' \strong{Per-replication steps} (`set.seed(seed + rep)`): (1) direct-GE DGP
#' truth at `theta*` (nonlinear snapshot `D1`, no clipping); (2) simulate a
#' macro series `Y` from the truth state space, and draw two independent
#' size-`N` multinomial cross-sections `Dhat0, Dhat1` from `D0, D1`, forming
#' the observed projected reweighting `m_hat = t(B_ref) (Dhat1 - Dhat0)`; (3)
#' build ONE `B_ref`-basis reweighting covariance for this replication from the
#' TRUTH `(D0, D1)` the DGP already computed (`rw_metric = "fixed_truth"`, the
#' default and orchestrator-certified choice): `Sigma = (cov_D0(B_ref) +
#' cov_D1(B_ref)) / N`, matching the shipped
#' \code{\link{hank_mixture_joint_logpost}}/\code{\link{hank_reweight_functional_metric}}
#' convention of a metric fixed across candidate thetas within one evaluation.
#' Because this covariance does NOT vary across eval nodes within a
#' replication, its log-determinant is an additive constant that cancels
#' exactly under the softmax normalization in step (5) -- so `ll_rw` is the
#' plain quadratic form with no log-det term, `-0.5 * t(d_k) Sigma^-1 d_k`,
#' `d_k = m_hat - rw_scale * gref_emu(theta_k)`. The alternative `rw_metric =
#' "emulated_varying"` (a THETA-VARYING emulator-interpolated covariance,
#' \code{\link{hank_emulator_metric}}, WITH its log-determinant term, since
#' there it does vary across eval nodes and does not cancel) is kept as an
#' option but is not the certified default; (4) evaluate the joint
#' log-posterior over every eval-grid node: `loglik_macro + loglik_rw`; (5)
#' normalize to a joint posterior pmf and compute, per coordinate, the
#' exact-randomized PIT and its floor-`L_ranks` integer rank; (6) record
#' posterior means and concentration.
#'
#' @param emu A \code{\link{hank_mixture_emulator}}.
#' @param n_rep Integer number of SBC replications.
#' @param prior Currently only the DEFAULT is implemented: a discrete uniform
#'   over the posterior evaluation grid's nodes (see Details). Reserved for a
#'   future `list(rsample = , in_support = )` user-supplied prior.
#' @param sampler Currently only `"grid"` is implemented: a `d`-dimensional
#'   tensor grid of `n_eval` points per axis over `emu$box`. This is
#'   deliberate, not a stub: the mixture posterior is a tilted, strongly
#'   anisotropic ridge (see \code{\link{hank_mixture_joint_logpost}}) on which
#'   a diagonal-metric random-walk sampler fails, so the exact grid is the
#'   recommended route for this posterior class. For a non-grid posterior draw,
#'   seed a curvature-aware sampler from \code{\link{hank_mixture_laplace}}
#'   (mode + Laplace covariance) rather than a diagonal RWM.
#' @param channels Character subset of `c("macro", "level", "response")`
#'   selecting the observable channels combined in the likelihood (see the
#'   Observable channels section). Default `c("macro", "response")`.
#' @param n_eval Integer points per axis for the posterior evaluation grid
#'   (`n_eval^d` total nodes).
#' @param N Survey cross-section sample size (each of `Dhat0`, `Dhat1`).
#' @param T_data Integer macro time-series length.
#' @param me_frac Fraction of the reference series' own per-column standard
#'   deviation used to set the shared, truth-independent `me_var` (see
#'   Details).
#' @param rw_scale Reweighting shock-size multiplier (as in
#'   \code{\link{hank_emulator_reweight_mean}}); also used to build the
#'   direct-truth nonlinear snapshot's price path (`r_path = r + rw_scale *
#'   rw_shock$r`). Default `0.005` (50bps): orchestrator-verified to calibrate
#'   (a smaller `0.003` leaves the reweighting channel too weak).
#' @param T_data_seed Seed for the fixed reference measurement-error series.
#' @param seed Base seed; replication `rep` uses `seed + rep`.
#' @param L_ranks Integer rank resolution: ranks are integers in
#'   `0..L_ranks-1`.
#' @param rw_metric `"fixed_truth"` (default; orchestrator-certified: one
#'   per-replication covariance built from the DGP's own truth `(D0, D1)`, no
#'   log-det term) or `"emulated_varying"` (theta-varying emulated covariance
#'   WITH its log-determinant term) -- see Details.
#'
#' @return An object of class `c("hank_mixture_sbc", "dynhr_sbc")`: a list
#'   with `pit` (`n_rep x d`, continuous), `ranks` (`n_rep x d` integer,
#'   `colnames = names(emu$box)`), `truth` (`n_rep x d`), `post_mean` (`n_rep x
#'   d`), `post_conc` (length-`n_rep`, `1 / sum(post^2)`), `uniformity`
#'   (\code{sbc_uniformity_test} on `ranks`), and `n_rep, box, sampler,
#'   channels, N, T_data, rw_scale, me_var, rw_metric`.
#' @seealso \code{\link{hank_mixture_emulator}}, \code{\link{hank_emulator_level_mean}},
#'   \code{sbc_uniformity_test}
#' @export
hank_mixture_sbc <- function(emu, n_rep, prior = NULL, sampler = "grid",
                              channels = c("macro", "response"),
                              n_eval = 25L, N = 5000L, T_data = 200L,
                              me_frac = 0.10, rw_scale = 0.005,
                              T_data_seed = 7L, seed = 1L, L_ranks = 1000L,
                              rw_metric = c("fixed_truth", "emulated_varying")) {
  if (!inherits(emu, "hank_mixture_emulator"))
    stop("hank_mixture_sbc(): 'emu' must be a hank_mixture_emulator.")
  sampler <- match.arg(sampler, "grid")
  rw_metric <- match.arg(rw_metric)
  channels <- match.arg(channels, c("macro", "level", "response"), several.ok = TRUE)
  use_macro <- "macro" %in% channels
  use_level <- "level" %in% channels
  use_resp  <- "response" %in% channels
  if (!is.null(prior))
    stop("hank_mixture_sbc(): only the default prior (uniform over the eval grid) is implemented.")

  box <- emu$box; nm <- names(box); d <- length(box)
  cfg <- emu$config

  inc <- hank_income_rouwenhorst(rho = 0.9, sigma = 0.7, n = cfg$n_e)
  ag  <- hank_asset_grid(amax = cfg$amax, n = cfg$n_a, amin = 0)

  ## ---- posterior evaluation grid + PRECOMPUTED emulated pieces (data-free) --
  ## Every stored emulator field an active channel needs is interpolated ONCE
  ## per eval node here, outside the replication loop (they are data-free).
  eval_axes <- setNames(lapply(nm, function(n) seq(box[[n]][1L], box[[n]][2L], length.out = n_eval)), nm)
  eval_idx  <- do.call(expand.grid, c(lapply(eval_axes, seq_along), KEEP.OUT.ATTRS = FALSE))
  n_eval_nodes <- nrow(eval_idx)
  eval_nodes <- do.call(cbind, lapply(seq_len(d), function(k) eval_axes[[nm[k]]][eval_idx[[k]]]))
  colnames(eval_nodes) <- nm

  ss_list <- if (use_macro) vector("list", n_eval_nodes) else NULL
  g_mat   <- if (use_resp)  matrix(NA_real_, n_eval_nodes, emu$K) else NULL   # response mean
  gL_mat  <- if (use_level) matrix(NA_real_, n_eval_nodes, emu$K) else NULL   # LEVEL mean
  for (k in seq_len(n_eval_nodes)) {
    th_k <- setNames(as.list(eval_nodes[k, ]), nm)
    if (use_macro) ss_list[[k]] <- hank_emulator_state_space(emu, th_k)
    if (use_resp)  g_mat[k, ]   <- hank_emulator_reweight_mean(emu, th_k, rw_scale)
    if (use_level) gL_mat[k, ]  <- hank_emulator_level_mean(emu, th_k)
  }

  ## ---- fixed reference measurement-error variance (truth-independent) -----
  ME_VAR <- 0
  if (use_macro) {
    set.seed(T_data_seed)
    centre_target <- vapply(nm, function(n) mean(box[[n]]), numeric(1))
    ss_ref <- hank_emulator_state_space(emu, setNames(as.list(centre_target), nm))
    Y_ref <- .hank_sbc_sim_Y(ss_ref, T_data, me_var = 0, sigma_Z = cfg$shock_specs$Z$sigma,
                            observables = cfg$observables)
    ME_VAR <- (me_frac * min(apply(Y_ref, 2, stats::sd)))^2
  }

  ## ---- outputs -------------------------------------------------------------
  pit_mat  <- matrix(NA_real_, n_rep, d, dimnames = list(NULL, nm))
  rank_mat <- matrix(NA_integer_, n_rep, d, dimnames = list(NULL, nm))
  truth_mat <- matrix(NA_real_, n_rep, d, dimnames = list(NULL, nm))
  postmean_mat <- matrix(NA_real_, n_rep, d, dimnames = list(NULL, nm))
  post_conc <- numeric(n_rep)

  ## JOINT test quantities (close the marginal-SBC dependency blind spot).
  ## The per-coordinate PIT below ranks each MARGINAL; a channel-likelihood /
  ## emulator error in the (centre, spread) DEPENDENCE with correct marginals
  ## is invisible to it. Since the full joint pmf `post` is computed per rep, we
  ## additionally rank joint functionals g(theta) against it (`.sbc_weighted_pit`):
  ##   tq_logpost -- the joint log-posterior at the truth node (most sensitive to
  ##                 overall posterior SHAPE / normalization);
  ##   tq_inter   -- the centred product of the first two coordinates (targets the
  ##                 cross-term / correlation specifically); only when d >= 2.
  tq_nm  <- if (d >= 2L) c("logpost", "inter") else "logpost"
  pit_tq  <- matrix(NA_real_, n_rep, length(tq_nm), dimnames = list(NULL, tq_nm))
  rank_tq <- matrix(NA_integer_, n_rep, length(tq_nm), dimnames = list(NULL, tq_nm))

  ## B_ref-basis functional covariance of a single projected size-N survey.
  cov_D <- function(D) crossprod(emu$B_ref, D * emu$B_ref) - tcrossprod(crossprod(emu$B_ref, D))

  for (rep in seq_len(n_rep)) {
    set.seed(seed + rep)

    ## ---- draw truth theta* from the (discrete) eval grid nodes -----------
    k_star <- sample.int(n_eval_nodes, 1L)
    theta_star <- setNames(as.list(eval_nodes[k_star, ]), nm)
    truth_mat[rep, ] <- eval_nodes[k_star, ]

    ## ---- (1) direct-GE DGP truth (nonlinear snapshot only if response used)
    tr <- .hank_sbc_direct_truth(theta_star, emu, inc, ag, rw_scale = rw_scale,
                                 need_D1 = use_resp)

    ## ---- (2) per-channel observed data + fixed_truth metrics --------------
    ## "fixed_truth" (default; orchestrator-certified): ONE covariance per rep
    ## and channel, built in the B_ref basis from the TRUTH distribution(s) the
    ## DGP already computed -- matching the shipped hank_mixture_joint_logpost /
    ## hank_reweight_functional_metric convention (a metric fixed across
    ## candidate thetas within an evaluation, -0.5*quadratic-form with NO
    ## log-det term, since a CONSTANT covariance's log-det is an additive
    ## constant across every eval node and cancels exactly in the softmax
    ## normalization below). "emulated_varying" is the alternative, theta-
    ## VARYING interpolated covariance INCLUDING its log-det term.
    if (use_macro)
      Y <- .hank_sbc_sim_Y(tr$ss_dir, T_data, me_var = ME_VAR, sigma_Z = cfg$shock_specs$Z$sigma,
                          observables = cfg$observables)
    if (use_resp) {
      Dhat0 <- as.numeric(stats::rmultinom(1, N, tr$D0)) / N
      Dhat1 <- as.numeric(stats::rmultinom(1, N, tr$D1)) / N
      m_hat <- as.numeric(crossprod(emu$B_ref, Dhat1 - Dhat0))
      if (rw_metric == "fixed_truth") Fc_R <- chol((cov_D(tr$D0) + cov_D(tr$D1)) / N)
    }
    if (use_level) {
      ## an INDEPENDENT stationary wealth cross-section (its own survey)
      DhatL <- as.numeric(stats::rmultinom(1, N, tr$D0)) / N
      m_lev <- as.numeric(crossprod(emu$B_ref, DhatL))
      if (rw_metric == "fixed_truth") Fc_L <- chol(cov_D(tr$D0) / N)
    }

    ## ---- (3) loglik over every eval node (sum of active channels) ---------
    loglik <- numeric(n_eval_nodes)
    for (k in seq_len(n_eval_nodes)) {
      ll <- 0
      if (use_macro) ll <- ll + hank_kalman_loglik(Y, ss_list[[k]], me_var = ME_VAR)
      if (use_resp) {
        d_k <- m_hat - g_mat[k, ]
        if (rw_metric == "fixed_truth") {
          ll <- ll - 0.5 * sum(backsolve(Fc_R, d_k, transpose = TRUE)^2)
        } else {
          Fc <- chol(hank_emulator_metric(emu, setNames(as.list(eval_nodes[k, ]), nm), rw_scale, N))
          ll <- ll - 0.5 * (sum(backsolve(Fc, d_k, transpose = TRUE)^2) + 2 * sum(log(diag(Fc))))
        }
      }
      if (use_level) {
        dL_k <- m_lev - gL_mat[k, ]
        if (rw_metric == "fixed_truth") {
          ll <- ll - 0.5 * sum(backsolve(Fc_L, dL_k, transpose = TRUE)^2)
        } else {
          FcL <- chol(hank_emulator_level_metric(emu, setNames(as.list(eval_nodes[k, ]), nm), N))
          ll <- ll - 0.5 * (sum(backsolve(FcL, dL_k, transpose = TRUE)^2) + 2 * sum(log(diag(FcL))))
        }
      }
      loglik[k] <- ll
    }

    ## ---- (4) normalize -> joint pmf; per-coordinate exact-randomized PIT --
    post <- exp(loglik - max(loglik)); post <- post / sum(post)
    post_conc[rep] <- 1 / sum(post^2)

    for (j in seq_len(d)) {
      marg <- tapply(post, eval_idx[[j]], sum)
      marg <- as.numeric(marg[order(as.integer(names(marg)))])
      idx_star <- eval_idx[[j]][k_star]
      u <- sum(marg[seq_len(idx_star - 1L)]) + stats::runif(1) * marg[idx_star]
      pit_mat[rep, j] <- u
      rank_mat[rep, j] <- min(floor(u * L_ranks), L_ranks - 1L)
      postmean_mat[rep, j] <- sum(marg * eval_axes[[nm[j]]])
    }

    ## ---- (4b) JOINT test-quantity PITs against the same joint pmf ----------
    u_lp <- .sbc_weighted_pit(loglik[k_star], loglik, post)
    pit_tq[rep, "logpost"]  <- u_lp
    rank_tq[rep, "logpost"] <- min(floor(u_lp * L_ranks), L_ranks - 1L)
    if (d >= 2L) {
      c1 <- eval_nodes[, 1L]; c2 <- eval_nodes[, 2L]
      g  <- (c1 - sum(post * c1)) * (c2 - sum(post * c2))   # centred cross-term
      u_in <- .sbc_weighted_pit(g[k_star], g, post)
      pit_tq[rep, "inter"]  <- u_in
      rank_tq[rep, "inter"] <- min(floor(u_in * L_ranks), L_ranks - 1L)
    }
  }

  structure(
    list(pit = pit_mat, ranks = rank_mat, truth = truth_mat, post_mean = postmean_mat,
         post_conc = post_conc, uniformity = sbc_uniformity_test(rank_mat, n_bins = NULL),
         ## JOINT test-quantity PITs/ranks + their own uniformity verdict. These
         ## bite on (centre, spread) DEPENDENCE errors the marginal `uniformity`
         ## is blind to; treat them as an ADDITIONAL calibration gate.
         pit_tq = pit_tq, ranks_tq = rank_tq,
         uniformity_tq = sbc_uniformity_test(rank_tq, n_bins = NULL),
         n_rep = n_rep, box = box, sampler = sampler, channels = channels, N = N, T_data = T_data,
         rw_scale = rw_scale, me_var = ME_VAR, rw_metric = rw_metric),
    class = c("hank_mixture_sbc", "dynhr_sbc"))
}


#' Print method for `hank_mixture_sbc`
#'
#' Summarizes recovery correlation per coordinate and the
#' \code{sbc_uniformity_test} verdict.
#' @param x A \code{\link{hank_mixture_sbc}} object.
#' @param ... Unused, present for S3 consistency.
#' @export
print.hank_mixture_sbc <- function(x, ...) {
  cat(sprintf("<hank_mixture_sbc> n_rep = %d, channels = {%s}, rw_metric = %s\n",
              x$n_rep, paste(x$channels, collapse = ", "), x$rw_metric))
  for (nmj in colnames(x$truth)) {
    cor_j <- stats::cor(x$truth[, nmj], x$post_mean[, nmj])
    cat(sprintf("  %-8s recovery cor = %.3f\n", nmj, cor_j))
  }
  cat(sprintf("  uniformity verdict: %s\n", x$uniformity$verdict))
  invisible(x)
}
