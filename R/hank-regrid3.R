## R/hank-regrid3.R
## --------------------------------------------------------------------------
## W3 (calibration-safe warm starts): grid continuation for the three-asset
## household. The paper's calibration ladder moves from coarse identification
## grids to a fine production grid; hank_egm3_regrid_values() interpolates a
## solved household's Vd/Vf/Va marginal values from one (d,f,a) grid onto
## another, preserving the (e,d,f,a) array convention, so the NEXT rung's
## solve can warm-start (via hank_het3_block's/hank_egm3_solve's Vd_init/
## Vf_init/Va_init) instead of starting cold on the refined grid.
## --------------------------------------------------------------------------


#' Regrid three-asset marginal values onto a new \code{(d, f, a)} grid
#'
#' Trilinear interpolation of a solved household's \code{Vd}/\code{Vf}/
#' \code{Va} marginal-value arrays from their current \code{(d, f, a)} grids
#' onto NEW \code{(d, f, a)} grids, preserving the \code{(e, d, f, a)} array
#' convention every other three-asset function uses. Built for grid
#' continuation: solve coarse, regrid the marginal values onto a refined grid,
#' and warm-start the refined solve from them
#' (\code{\link{hank_het3_block}}'s / \code{\link{hank_egm3_solve}}'s
#' \code{Vd_init}/\code{Vf_init}/\code{Va_init}) instead of starting cold.
#'
#' The income axis is NOT interpolated -- it is a discrete state, and this
#' function has no argument through which a caller could ask it to touch that
#' axis. Every \code{e}-slice is regridded independently over \code{(d, f,
#' a)} only.
#'
#' Extrapolation outside the SOURCE grid range is CLAMPED (constant), not
#' linear: \code{Vd}/\code{Vf}/\code{Va} are convex, decaying marginal values,
#' and linear extrapolation off the top of a grid can go NEGATIVE, which would
#' poison a warm start with an invalid \eqn{u'(c)}. This is implemented by
#' clamping the query points into \code{[min(source grid), max(source grid)]}
#' before calling the package's existing three-axis interpolator
## NOT \link{}: .hank3_interp3 is an unexported internal with no Rd page, and
## R CMD check flags a cross-reference to a topic that does not exist.
#' (\code{.hank3_interp3}, \code{R/hank-diagnostics3.R}) rather than by
#' writing a fourth independent interpolator: at a clamped boundary query, that
#' interpolator's weight is exactly \code{0} or \code{1} on the boundary knot
#' (see \code{.hank3_axis_interp}), so the result is the boundary VALUE, not a
#' linear projection past it. Interpolated output is asserted finite and
#' STRICTLY POSITIVE (a violation means the source marginal values were
#' already invalid, or a bug in this reconstruction) -- never returned
#' silently.
#'
#' @param hh A solved three-asset household: a \code{\link{hank_egm3_solve}}
#'   result or a \code{\link{hank_het3_block}}, i.e. anything carrying
#'   \code{Vd}, \code{Vf}, \code{Va} (each \code{n_e x n_d x n_f x n_a}) and
#'   its own \code{d_grid}, \code{f_grid}, \code{a_grid}.
#' @param d_grid,f_grid,a_grid New increasing grids to interpolate onto, same
#'   validity contract as \code{\link{hank_egm3_solve}}: \code{d_grid} and
#'   \code{a_grid} need length \code{>= 3}; \code{f_grid} needs length
#'   \code{>= 3} OR the singleton \code{0} (the two-asset reduction). All
#'   entries must be finite and strictly increasing.
#'
#' @param method How to build the target-grid marginal values.
#'   \code{"values"} (default) interpolates the three arrays independently.
#'   \code{"envelope"} instead interpolates the POLICIES and rebuilds
#'   \code{Vd}/\code{Vf}/\code{Va} from the solver's own envelope at the
#'   TARGET grid, keeping the three mutually consistent by construction.
#'
#'   \strong{\code{"envelope"} is better founded but is NOT a cure, and the
#'   default stays \code{"values"} on that evidence.} On the reference
#'   coarsening case it moved the failing cell from
#'   \code{(e=5, f=2, a=41)} to \code{(e=5, f=2, a=42)} and nothing else, so
#'   it is offered as an option rather than promoted on theory alone.
#'
#'   The distinction is not cosmetic. The three arrays are tied together
#'   by \eqn{V_d = u'(c)(1+r_d)}, \eqn{V_f = u'(c)(p_x(1+r_f) - \Psi_2)}
#'   and \eqn{V_a = u'(c)((1+r_a) - \Psi_2)}, whose \eqn{\Psi_2} terms are
#'   evaluated at the CURRENT asset holding. Interpolating them separately
#'   therefore carries the SOURCE grid's adjustment-cost geometry onto the
#'   target abscissae, producing a triple that is the marginal-value
#'   function of no policy at all. It looks healthy -- finite, positive,
#'   monotone -- and early EGM steps consume it, which is what makes the
#'   resulting failure hard to attribute: on a measured
#'   \code{(14,24,8,24)} to \code{(14,72,4,80)} continuation it ran 196
#'   clean iterations before the implied liquid grid went non-monotone,
#'   while a COLD start on the same target grid converged in 1,033.
#' @return A list with \code{Vd}, \code{Vf}, \code{Va} (each
#'   \code{n_e x length(d_grid) x length(f_grid) x length(a_grid)}, finite and
#'   strictly positive) and the echoed \code{d_grid}, \code{f_grid},
#'   \code{a_grid}, ready to pass straight through as
#'   \code{Vd_init}/\code{Vf_init}/\code{Va_init}. When the requested grids are
#'   IDENTICAL to \code{hh}'s own (\code{identical()}, entrywise), the input
#'   \code{Vd}/\code{Vf}/\code{Va} arrays are returned unchanged (no
#'   interpolation, no floating-point drift).
#' @seealso \code{\link{hank_het3_block}}, \code{\link{hank_egm3_solve}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' hh <- hank_egm3_solve(dg, fg, ag, e, Pi, rd = .01, rf = .015, ra = .02,
#'                       beta = .97, eis = .5, chi1 = .2, phi1 = .1,
#'                       tol = 1e-5, maxit = 250)
#' dg2 <- seq(-.1, .5, length.out = 5); fg2 <- seq(0, .5, length.out = 4)
#' ag2 <- seq(0, .6, length.out = 5)
#' rg <- hank_egm3_regrid_values(hh, dg2, fg2, ag2)
#' hh2 <- hank_egm3_solve(dg2, fg2, ag2, e, Pi, rd = .01, rf = .015, ra = .02,
#'                        beta = .97, eis = .5, chi1 = .2, phi1 = .1,
#'                        tol = 1e-5, maxit = 250,
#'                        Vd_init = rg$Vd, Vf_init = rg$Vf, Va_init = rg$Va)
#' hh2$converged
#' @export
hank_egm3_regrid_values <- function(hh, d_grid, f_grid, a_grid,
                                   method = c("values", "envelope")) {
  method <- match.arg(method)
  need <- c("Vd", "Vf", "Va", "d_grid", "f_grid", "a_grid")
  if (!is.list(hh) || !all(need %in% names(hh)))
    stop("hank_egm3_regrid_values: 'hh' must carry Vd/Vf/Va and its own ",
         "d_grid/f_grid/a_grid (a hank_egm3_solve result or a ",
         "hank_het3_block).")

  old_d <- hh$d_grid; old_f <- hh$f_grid; old_a <- hh$a_grid
  if (!is.array(hh$Vd) || length(dim(hh$Vd)) != 4L)
    stop("hank_egm3_regrid_values: 'hh$Vd' must be an e x d x f x a array.")
  ne <- dim(hh$Vd)[1L]
  want_dim <- c(ne, length(old_d), length(old_f), length(old_a))
  if (!identical(dim(hh$Vd), want_dim) || !identical(dim(hh$Vf), want_dim) ||
      !identical(dim(hh$Va), want_dim))
    stop("hank_egm3_regrid_values: 'hh$Vd'/'Vf'/'Va' dimensions are ",
         "inconsistent with 'hh's own d_grid/f_grid/a_grid.")

  ## Same validity contract as hank_egm3_solve(): d/a need length >= 3;
  ## f may additionally be the singleton-0 two-asset reduction. Checked on
  ## BOTH the new grids and hh's own (an hh built by something other than
  ## this package's solvers could carry a malformed source grid).
  .check_grid <- function(g, nm) {
    if (!is.numeric(g) || length(g) < 1L || !all(is.finite(g)))
      stop("hank_egm3_regrid_values: '", nm, "' must be finite numeric.")
    reduction <- length(g) == 1L && isTRUE(g[1L] == 0)
    if (nm %in% c("f_grid", "hh$f_grid")) {
      if (!(length(g) >= 3L || reduction))
        stop("hank_egm3_regrid_values: '", nm, "' must have length >= 3, ",
             "or be the singleton 0 (two-asset reduction).")
    } else if (length(g) < 3L) {
      stop("hank_egm3_regrid_values: '", nm, "' must have length >= 3.")
    }
    if (length(g) > 1L && any(diff(g) <= 0))
      stop("hank_egm3_regrid_values: '", nm, "' must be strictly increasing.")
  }
  .check_grid(d_grid, "d_grid"); .check_grid(f_grid, "f_grid"); .check_grid(a_grid, "a_grid")
  .check_grid(old_d, "hh$d_grid"); .check_grid(old_f, "hh$f_grid"); .check_grid(old_a, "hh$a_grid")

  ## COARSENING WARNING. Measured failure mode, not a guess. Continuing a solve
  ## onto a grid that REFINES (or holds) every axis works and is worth doing:
  ## (14,24,8,24) -> (14,72,8,80) converged in 509 iterations against 1,033 for
  ## a cold start on the same target. Continuing onto a grid that COARSENS an
  ## axis does not: (14,24,8,24) -> (14,72,4,80) tripped the endogenous-grid
  ## safeguard at iteration 197, and so did the same move onto four knots taken
  ## verbatim FROM the source grid -- so it is the loss of resolution, not the
  ## knot placement. A cold start on that same coarse target converged fine
  ## (1,033 iterations, 256 s), which is the recommended remedy.
  ##
  ## Warn rather than error: the regrid itself is well defined, the failure
  ## surfaces later in the EGM iteration, and a caller may legitimately want
  ## the coarse values for something other than a warm start.
  coarsened <- c(d = length(d_grid) < length(old_d),
                 f = length(f_grid) < length(old_f),
                 a = length(a_grid) < length(old_a))
  if (any(coarsened))
    warning("hank_egm3_regrid_values: this regrid COARSENS the ",
            paste(names(coarsened)[coarsened], collapse = "/"),
            " axis. Warm-starting hank_egm3_solve() from a coarsened regrid ",
            "is a MEASURED failure mode: the endogenous-grid safeguard trips ",
            "part-way through the iteration (at 197 of 2000 on the reference ",
            "case), and damping only postpones it. Refining or holding every ",
            "axis is safe and is worth ~2x against a cold start. If you must ",
            "coarsen, COLD-START the target grid instead -- that converges. ",
            "See ?hank_egm3_regrid_values.", call. = FALSE)

  ## Identity grids: return the input arrays unchanged, bit for bit -- no
  ## interpolation pass, and therefore no floating-point drift versus hh's own
  ## Vd/Vf/Va, even though mathematically an interpolation AT every source
  ## knot would reproduce them (up to rounding).
  if (identical(d_grid, old_d) && identical(f_grid, old_f) && identical(a_grid, old_a))
    return(list(Vd = hh$Vd, Vf = hh$Vf, Va = hh$Va,
                d_grid = d_grid, f_grid = f_grid, a_grid = a_grid,
                method = "identity"))

  nd <- length(d_grid); nf <- length(f_grid); na <- length(a_grid)

  ## CLAMPED (constant) extrapolation: see Details above for why linear
  ## extrapolation is rejected. expand.grid()'s default column order (first
  ## argument fastest) matches array(., c(nd, nf, na))'s fill order (first
  ## dimension fastest), so grid_q's rows line up 1:1 with the target array
  ## in the SAME order this reshape uses -- no separate index bookkeeping.
  qd <- pmin(pmax(d_grid, min(old_d)), max(old_d))
  qf <- pmin(pmax(f_grid, min(old_f)), max(old_f))
  qa <- pmin(pmax(a_grid, min(old_a)), max(old_a))
  grid_q <- expand.grid(d = qd, f = qf, a = qa, KEEP.OUT.ATTRS = FALSE)

  ## Policies may sit AT a grid floor of 0, so they get the same interpolation
  ## without the strict-positivity check that guards marginal values.
  regrid_one_signed <- function(V, nm) {
    out <- array(NA_real_, c(ne, nd, nf, na))
    for (ei in seq_len(ne)) {
      Z <- array(V[ei, , , ], c(length(old_d), length(old_f), length(old_a)))
      out[ei, , , ] <- array(.hank3_interp3(Z, old_d, old_f, old_a,
                                            grid_q$d, grid_q$f, grid_q$a),
                            c(nd, nf, na))
    }
    if (any(!is.finite(out)))
      stop("hank_egm3_regrid_values: interpolated '", nm, "' has non-finite ",
           "entries.")
    out
  }

  regrid_one <- function(V, nm) {
    out <- array(NA_real_, c(ne, nd, nf, na))
    for (ei in seq_len(ne)) {
      Z <- array(V[ei, , , ], c(length(old_d), length(old_f), length(old_a)))
      out[ei, , , ] <- array(.hank3_interp3(Z, old_d, old_f, old_a,
                                            grid_q$d, grid_q$f, grid_q$a),
                            c(nd, nf, na))
    }
    if (any(!is.finite(out)))
      stop("hank_egm3_regrid_values: interpolated '", nm, "' has non-finite ",
           "entries.")
    if (any(out <= 0))
      stop("hank_egm3_regrid_values: interpolated '", nm, "' has a ",
           "non-positive entry -- marginal values must stay strictly ",
           "positive; this would poison a warm start with an invalid u'(c).")
    out
  }
  if (identical(method, "values")) {
    Vd <- regrid_one(hh$Vd, "Vd"); Vf <- regrid_one(hh$Vf, "Vf")
    Va <- regrid_one(hh$Va, "Va")
    return(list(Vd = Vd, Vf = Vf, Va = Va, d_grid = d_grid,
                f_grid = f_grid, a_grid = a_grid, method = "values"))
  }

  ## ---- method = "envelope" (default) --------------------------------------
  ## Interpolate the POLICIES and rebuild the marginal values from the
  ## kernel's own envelope on the NEW grid, rather than interpolating the
  ## three value arrays independently.
  ##
  ## Why: the three are not free of one another. The kernel sets
  ##   Vd = u'(c) (1 + rd)
  ##   Vf = u'(c) (px (1 + rf) - Psi2(f', f))
  ##   Va = u'(c) ((1 + ra) - Psi2(a', a))
  ## so their RATIOS are pinned by adjustment-cost derivatives evaluated at the
  ## CURRENT asset holding. Interpolating each array separately carries the
  ## SOURCE grid's Psi2 geometry onto the target abscissae, and the result is
  ## a triple that is no marginal-value function of any policy. It looks
  ## perfectly healthy -- finite, positive, monotone -- and the first EGM steps
  ## consume it happily, which is what makes the failure so confusing: on the
  ## paper's (14,24,8,24) -> (14,72,4,80) continuation it ran 196 clean
  ## iterations before the implied liquid grid went non-monotone, while a COLD
  ## start on the SAME target grid converged in 1,033 iterations. Damping only
  ## postponed it (197 -> 689 -> 1368 as relax fell 1 -> 0.5 -> 0.25).
  ##
  ## Rebuilding from the envelope makes the triple mutually consistent BY
  ## CONSTRUCTION at the target grid, because Psi2 is re-evaluated there.
  need_pol <- c("c", "f", "a")
  if (!all(need_pol %in% names(hh)) ||
      any(vapply(hh[need_pol], is.null, TRUE)))
    stop("hank_egm3_regrid_values: method = \"envelope\" needs the policies ",
         "'c', 'f' and 'a' on 'hh' (any hank_egm3_solve or hank_het3_block ",
         "carries them). Use method = \"values\" for an object that does not.")
  px <- if (is.null(hh$px)) 1 else hh$px
  c_new <- regrid_one(hh$c, "c")
  f_new <- regrid_one_signed(hh$f, "f")
  a_new <- regrid_one_signed(hh$a, "a")
  FF <- .hank3_bcast_f(f_grid, ne, nd, nf, na)
  AA <- .hank3_bcast_a(a_grid, ne, nd, nf, na)
  p2f <- .hank_psi(f_new, FF, hh$rf, hh$phi0, hh$phi1, hh$phi2)$Psi2
  p2a <- .hank_psi(a_new, AA, hh$ra, hh$chi0, hh$chi1, hh$chi2)$Psi2
  uc <- c_new^(-1 / hh$eis)
  Vd <- uc * (1 + hh$rd)
  Vf <- uc * (px * (1 + hh$rf) - p2f)
  Va <- uc * ((1 + hh$ra) - p2a)
  for (nm in c("Vd", "Vf", "Va")) {
    V <- get(nm)
    if (any(!is.finite(V)) || any(V <= 0))
      stop("hank_egm3_regrid_values: envelope-rebuilt '", nm, "' has a ",
           "non-finite or non-positive entry. That means the regridded ",
           "policies imply an inadmissible marginal value at some cell; ",
           "method = \"values\" reproduces the pre-0.9.0.0017 behaviour.")
  }
  list(Vd = Vd, Vf = Vf, Va = Va, d_grid = d_grid, f_grid = f_grid,
       a_grid = a_grid, method = "envelope")
}
