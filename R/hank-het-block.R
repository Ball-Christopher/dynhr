## R/hank-het-block.R
## --------------------------------------------------------------------------
## Heterogeneous-agent "het block": packages a steady-state household solution
## (EGM policies + Young's-method distribution) with the operations needed to
## compute its sequence-space Jacobian (see R/hank-jacobian.R) and nonlinear
## perfect-foresight transitions.
##
## The reference household is a one-asset Krusell-Smith / one-asset-HANK
## household: CRRA utility, a single asset with return r, labour income
## y(e) = w * e over the idiosyncratic productivity states e.  Aggregate block
## OUTPUTS are:
##   A_t = sum_x D_t(x) a'(x)   (aggregate end-of-period assets / savings)
##   C_t = sum_x D_t(x) c(x)    (aggregate consumption)
## aggregated with the BEGINNING-of-period-t distribution D_t.  Block INPUTS are
## the aggregate paths {r_t, w_t}.
## --------------------------------------------------------------------------


#' Construct a one-asset heterogeneous-agent household block at steady state
#'
#' Solves the household EGM problem and its stationary distribution at fixed
#' aggregate prices \code{(r, w)}, and stores everything the sequence-space
#' Jacobian and nonlinear-transition routines need.
#'
#' @param a_grid Numeric: asset grid (see \code{\link{hank_asset_grid}}).
#' @param Pi Numeric \code{n_e x n_e}: income transition matrix.
#' @param e Numeric length-\code{n_e}: income levels (see
#'   \code{\link{hank_income_rouwenhorst}}).
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param r,w Steady-state real return and wage.
#' @param tol,maxit Passed to \code{\link{hank_egm_solve}}.
#' @param amin Numeric: borrowing constraint (minimum end-of-period assets).
#'   Must satisfy \code{amin >= a_grid[1]}. Defaults to \code{NULL}, which
#'   resolves to \code{a_grid[1L]} (the historical hardcoded behavior, so
#'   existing calls are byte-identical). Lets different household types on a
#'   shared \code{a_grid} face distinct borrowing limits (the wealth axis).
#'
#' @return An object of class \code{hank_het_block} with the steady-state
#'   policies (\code{a}, \code{c}), marginal value \code{Va}, distribution
#'   \code{D} (vector) and forward operator \code{Lambda}, aggregate
#'   steady-state outputs \code{A}, \code{C}, the borrowing constraint
#'   \code{amin}, and the calibration.
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 5)
#' ag  <- hank_asset_grid(50, 100, 0)
#' blk <- hank_het_block(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                       r = 0.01, w = 1.0)
#' blk$A   # aggregate assets
#' @export
hank_het_block <- function(a_grid, Pi, e, beta, eis, r, w,
                           tol = 1e-11, maxit = 5000L, amin = NULL) {
  if (is.null(amin)) amin <- a_grid[1L]
  y  <- w * e
  hh <- hank_egm_solve(a_grid, y = y, r = r, beta = beta, eis = eis, Pi = Pi,
                       tol = tol, maxit = maxit, amin = amin)
  if (!hh$converged)
    warning("hank_het_block: household EGM did not converge at steady state")
  Lam <- hank_forward_operator(hh$a, a_grid, Pi)
  sd  <- hank_stationary_dist(Lam)
  D   <- sd$d
  structure(
    list(a_grid = a_grid, Pi = Pi, e = e, beta = beta, eis = eis,
         r = r, w = w, amin = amin,
         a = hh$a, c = hh$c, Va = hh$Va,
         Lambda = Lam, D = D,
         A = hank_aggregate(D, hh$a),
         C = hank_aggregate(D, hh$c),
         n_e = length(e), n_a = length(a_grid),
         dist_converged = sd$converged),
    class = "hank_het_block")
}


#' One nonlinear EGM backward step at given aggregate prices
#'
#' Thin wrapper mapping block inputs \code{(r, w)} to the household EGM step.
#' @param block A \code{\link{hank_het_block}}.
#' @param Va_p Next-period marginal value (\code{n_e x n_a}).
#' @param r,w Aggregate return and wage this period.
#' @return List with \code{Va}, \code{a}, \code{c} (see \code{.hank_egm_step}).
#' @keywords internal
.hank_block_step <- function(block, Va_p, r, w) {
  amin <- if (!is.null(block$amin)) block$amin else block$a_grid[1L]
  .hank_egm_step(Va_p, block$a_grid, y = w * block$e, r = r,
                 beta = block$beta, eis = block$eis, Pi = block$Pi,
                 amin = amin)
}


#' Nonlinear perfect-foresight transition of a het block
#'
#' Given aggregate input PATHS \code{r_path}, \code{w_path} over horizon
#' \code{T} (with terminal conditions returning to steady state), computes the
#' aggregate output paths by a backward household solve followed by a forward
#' distribution simulation.  This is the block's nonlinear map from input paths
#' to output paths, and is the brute-force reference used to validate the
#' sequence-space Jacobian.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param r_path,w_path Numeric length-\code{T} input paths (levels).  Missing
#'   entries default to the steady-state value.
#' @param T_h Integer horizon (default \code{length(r_path)}).
#' @param keep_policies logical (default \code{FALSE}). When \code{TRUE},
#'   also return the per-period objects the solve builds and otherwise
#'   discards: \code{c_pol}/\code{a_pol} (length-\code{T_h} lists of
#'   \code{n_e x n_a} policy matrices from the backward pass) and
#'   \code{Lambda} (length-\code{T_h} list of sparse per-period forward
#'   operators from the forward pass; same no-transpose convention as
#'   \code{block$Lambda} -- distributions push forward via
#'   \code{t(Lambda[[t]]) \%*\% d}). Inputs for date-indexed distributional/
#'   welfare analysis along the transition. Off by default (\code{Lambda}
#'   costs \code{O(T_h)} sparse \code{n_cell x n_cell} matrices of memory).
#'
#' @return A list with numeric length-\code{T} paths \code{A} and \code{C}, and
#'   \code{Dpath} (\code{(n_e*n_a) x T} matrix): column \code{t} is the
#'   beginning-of-period-\code{t} distribution over which \code{A[t]},
#'   \code{C[t]} are aggregated. When \code{keep_policies = TRUE}, also
#'   \code{c_pol}, \code{a_pol}, and \code{Lambda} (see above); the exact
#'   aggregation identity \code{C[t] == hank_aggregate(Dpath[, t],
#'   c_pol[[t]])} holds by construction.
#' @export
hank_td_nonlinear <- function(block, r_path = NULL, w_path = NULL, T_h = NULL,
                              keep_policies = FALSE) {
  if (is.null(T_h))
    T_h <- max(length(r_path), length(w_path),
               if (is.null(r_path) && is.null(w_path)) 1L else 0L)
  if (is.null(r_path)) r_path <- rep(block$r, T_h)
  if (is.null(w_path)) w_path <- rep(block$w, T_h)
  stopifnot(length(r_path) == T_h, length(w_path) == T_h)

  ## Backward: terminal Va_{T+1} = Va_ss.
  a_pol <- vector("list", T_h)
  c_pol <- vector("list", T_h)
  Va <- block$Va
  for (t in T_h:1L) {
    step <- .hank_block_step(block, Va, r_path[t], w_path[t])
    a_pol[[t]] <- step$a
    c_pol[[t]] <- step$c
    Va <- step$Va
  }

  ## Forward: beginning distribution D_1 = D_ss.
  D <- block$D
  A <- numeric(T_h); C <- numeric(T_h)
  Dpath <- matrix(0, length(D), T_h)
  Lam_keep <- if (keep_policies) vector("list", T_h) else NULL
  for (t in seq_len(T_h)) {
    Dpath[, t] <- D
    A[t] <- hank_aggregate(D, a_pol[[t]])
    C[t] <- hank_aggregate(D, c_pol[[t]])
    Lam  <- hank_forward_operator(a_pol[[t]], block$a_grid, block$Pi)
    if (keep_policies) Lam_keep[[t]] <- Lam
    D    <- as.numeric(Matrix::t(Lam) %*% D)
  }
  out <- list(A = A, C = C, Dpath = Dpath)
  if (keep_policies) {
    out$c_pol  <- c_pol
    out$a_pol  <- a_pol
    out$Lambda <- Lam_keep
  }
  out
}
