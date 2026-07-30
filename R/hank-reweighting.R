## R/hank-reweighting.R
## --------------------------------------------------------------------------
## Phase-0 machinery for the het-preferences-in-HANK research line (brief-17):
## a K-type preference MIXTURE (types share the grid/income process; only the
## discount factor differs) plus the reweighting-loss diagnostic that compares
## an OBSERVED net survey reweighting against the model-implied one.
##
## Because types interact only through aggregate prices (r, w) -- never
## directly -- every mixture object (distribution, aggregate Jacobian,
## distribution Jacobian) is an EXACT omega-weighted sum of the per-type
## objects; there is no cross term.  This file is pure sequence-space
## PLUMBING: it builds the mixture objects, evaluates a quadratic loss between
## an observed and a model-implied net reweighting, and (hank_phase0_reweight_gate)
## assembles a synthetic-survey experiment that computes that loss over a grid
## of candidate (centre, spread) mixture calibrations. It does NOT judge
## whether the loss is small enough to "identify" the spread -- that is the
## orchestrator's call once it has the raw numbers.
##
## DATE CONVENTION (see hank_phase0_reweight_gate below): the synthetic survey
## step draws exactly TWO cross-sections, one from the pre-shock stationary
## mixture distribution D0 and one from the (first-order) post-shock
## distribution D1 = D0 + dD_star.  The observed net reweighting dD_hat is
## therefore a SINGLE realized quantity (Dhat1 - Dhat0), not a T_h-length
## path; it is compared against the model-implied net reweighting SUMMED over
## the shock horizon (dD_mod = Sum_t dD_model[, t], i.e. the total distribution
## displacement accumulated by the end of the impulse), matching what two
## snapshot cross-sections straddling the whole shock episode can identify.
## --------------------------------------------------------------------------


#' Build K heterogeneous-preference het blocks sharing one grid
#'
#' A K-type discount-factor mixture: every type shares the same asset grid,
#' income process and EIS -- only \code{beta} differs across types (the
#' object of interest for the het-preferences research line). Because types
#' never interact directly (only via aggregate prices), each block is an
#' ordinary \code{\link{hank_het_block}} solved independently at the SAME
#' \code{(r, w)}.
#'
#' @param a_grid Numeric asset grid (see \code{\link{hank_asset_grid}}),
#'   shared by every type.
#' @param Pi Numeric \code{n_e x n_e} income transition matrix, shared.
#' @param e Numeric length-\code{n_e} income levels, shared.
#' @param betas Numeric length-K vector of discount factors, one per type.
#' @param eis Elasticity of intertemporal substitution, shared (scalar; a
#'   single type-invariant value, since only \code{beta} varies here).
#' @param r,w Steady-state real return and wage (types interact only through
#'   these prices, so both are common across types).
#'
#' @return A list of length-K \code{\link{hank_het_block}} objects, one per
#'   entry of \code{betas}, in the same order.
#' @export
hank_mixture_blocks <- function(a_grid, Pi, e, betas, eis = 1, r, w) {
  if (!is.numeric(betas) || length(betas) < 1L)
    stop("hank_mixture_blocks(): 'betas' must be a non-empty numeric vector.")
  lapply(betas, function(b)
    hank_het_block(a_grid, Pi, e, beta = b, eis = eis, r = r, w = w))
}


#' Validate a mixture-weight vector against a list of blocks
#' @keywords internal
.hank_mixture_check_omega <- function(blocks, omega) {
  K <- length(blocks)
  if (length(omega) != K)
    stop(sprintf(
      "hank_mixture_*: length(omega) (%d) must equal length(blocks) (%d).",
      length(omega), K))
  if (any(!is.finite(omega)) || any(omega < 0))
    stop("hank_mixture_*: 'omega' must be finite and non-negative.")
  if (abs(sum(omega) - 1) > 1e-8)
    stop(sprintf(
      "hank_mixture_*: 'omega' must sum to 1 (sum = %.6f).", sum(omega)))
  invisible(TRUE)
}


#' Mixture stationary distribution and steady-state aggregates
#'
#' The omega-mixture steady-state aggregates \code{A, C = Sum_k omega_k
#' (A_k, C_k)} are ALWAYS valid, regardless of whether the K types share an
#' income process: each is a per-type SCALAR, so the omega-weighted sum never
#' requires the types to live on a common cell space.
#'
#' The pooled stationary distribution \code{D = Sum_k omega_k D_k}, by
#' contrast, is only well-defined when every type's distribution lives on
#' the IDENTICAL \code{(n_e*n_a)}-cell state space, i.e. all blocks share
#' \code{a_grid}, \code{Pi}, and \code{e} (see \code{\link{hank_mixture_blocks}}):
#' only then is cell \code{j} the same income-asset pair for every type, so
#' the pointwise weighted sum is a meaningful pooled distribution. A
#' per-type borrowing constraint \code{amin} (the wealth heterogeneity
#' axis) does NOT break this: \code{amin} lives on the shared
#' \code{a_grid}, so the cell space stays common and the pooled \code{D}
#' remains valid (a type simply carries zero mass below its own
#' \code{amin}). When the
#' blocks differ in their income process (distinct \code{Pi}/\code{e} -- an
#' income-risk heterogeneity axis), pooling cell-by-cell is semantically
#' invalid (cell \code{j} means a different income level per type), so this
#' function does NOT silently sum the per-type \code{D} vectors: it sets
#' \code{D = NULL} and instead returns \code{D_by_type}, a length-K list of
#' the per-type stationary distributions (each still on that type's own
#' \code{n_e*n_a} cells), plus \code{same_income = FALSE}.
#'
#' @param blocks List of K \code{\link{hank_het_block}} objects sharing an
#'   asset grid (as returned by \code{\link{hank_mixture_blocks}}, or built
#'   per-type with distinct income processes on a shared \code{a_grid}).
#' @param omega Numeric length-K mixture weights, non-negative, summing to 1.
#'
#' @return A list with \code{A}, \code{C} (mixture steady-state aggregates,
#'   always valid), \code{same_income} (logical: whether every block shares
#'   the identical \code{(Pi, e)}), \code{D} (length-\code{n_e*n_a} pooled
#'   mixture stationary distribution when \code{same_income} is \code{TRUE},
#'   else \code{NULL}), and \code{D_by_type} (length-K list of per-type
#'   stationary distributions; always populated, so callers can fall back to
#'   it when \code{D} is \code{NULL}).
#' @export
hank_mixture_dist <- function(blocks, omega) {
  .hank_mixture_check_omega(blocks, omega)
  A <- sum(vapply(seq_along(blocks), function(k) omega[k] * blocks[[k]]$A, numeric(1)))
  C <- sum(vapply(seq_along(blocks), function(k) omega[k] * blocks[[k]]$C, numeric(1)))
  D_by_type <- lapply(blocks, function(b) b$D)
  same_income <- length(blocks) <= 1L || all(vapply(blocks[-1L], function(b)
    isTRUE(all.equal(b$Pi, blocks[[1L]]$Pi)) &&
      isTRUE(all.equal(b$e, blocks[[1L]]$e)), logical(1)))
  D <- if (same_income) Reduce(`+`, Map(function(b, wk) wk * b$D, blocks, omega)) else NULL
  if (!same_income)
    message("hank_mixture_dist(): blocks differ in their income process ",
            "(Pi/e); a pooled cell-space D is not defined across differing ",
            "income spaces, so D = NULL. Use 'D_by_type' (per-type ",
            "distributions) instead.")
  list(D = D, A = A, C = C, D_by_type = D_by_type, same_income = same_income)
}


#' Mixture DISTRIBUTION sequence-space Jacobian
#'
#' The omega-weighted sum of the per-type distribution Jacobians
#' (\code{\link{hank_het_dist_jacobian}}). Exact by linearity: because types
#' interact only through \code{(r, w)}, a shock to input \code{i} moves EACH
#' type's own distribution by \code{omega_k dD_k}. When every type shares the
#' identical \code{(a_grid, Pi, e)} cell space, those per-type moves live on
#' the same cells and can be pooled into one array,
#' \code{JD_mix[[i]] = Sum_k omega_k JD_k[[i]]}, with NO cross term between
#' types.
#'
#' When the blocks differ in their income process (distinct \code{Pi}/\code{e}),
#' pooling cell-by-cell is semantically invalid (cell \code{j} indexes a
#' different income level per type), so this function does NOT silently sum
#' the per-type arrays: it returns \code{by_type} (a length-K list of the raw
#' per-type \code{JD_k} objects, each as returned by
#' \code{\link{hank_het_dist_jacobian}}) with the pooled \code{JD_mix}-style
#' top-level entries omitted, plus \code{same_income = FALSE}. Callers that
#' need a pooled array despite differing income spaces get a loud error
#' rather than a wrong number if they try to use the (absent) pooled fields.
#'
#' @param blocks List of K \code{\link{hank_het_block}} objects sharing an
#'   asset grid (income process \code{Pi}/\code{e} may differ across types).
#' @param omega Numeric length-K mixture weights, non-negative, summing to 1.
#' @param T_h Integer horizon.
#' @param inputs Character subset of \code{c("r", "w")}.
#'
#' @return When \code{same_income} is \code{TRUE}: named list
#'   \code{JD_mix[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_a)} (identical layout to
#'   \code{\link{hank_het_dist_jacobian}}'s return), PLUS a \code{same_income}
#'   attribute (\code{TRUE}) and a \code{by_type} attribute (the per-type
#'   Jacobians, for callers who want both views). When \code{same_income} is
#'   \code{FALSE}: a list with \code{same_income = FALSE} and \code{by_type}
#'   (length-K list of per-type \code{JD_k} objects); no pooled array is
#'   returned.
#' @export
hank_mixture_dist_jacobian <- function(blocks, omega, T_h, inputs = c("r", "w")) {
  .hank_mixture_check_omega(blocks, omega)
  per_type <- lapply(blocks, function(b)
    hank_het_dist_jacobian(b, T_h, inputs = inputs))
  same_income <- length(blocks) <= 1L || all(vapply(blocks[-1L], function(b)
    isTRUE(all.equal(b$Pi, blocks[[1L]]$Pi)) &&
      isTRUE(all.equal(b$e, blocks[[1L]]$e)), logical(1)))
  if (!same_income) {
    message("hank_mixture_dist_jacobian(): blocks differ in their income ",
            "process (Pi/e); a pooled cell-space distribution Jacobian is ",
            "not defined across differing income spaces. Returning ",
            "'by_type' (per-type Jacobians) with same_income = FALSE and no ",
            "pooled array.")
    return(list(same_income = FALSE, by_type = per_type))
  }
  JD_mix <- setNames(lapply(inputs, function(i) {
    acc <- omega[1L] * per_type[[1L]][[i]]
    if (length(blocks) > 1L)
      for (k in 2L:length(blocks)) acc <- acc + omega[k] * per_type[[k]][[i]]
    acc
  }), inputs)
  attr(JD_mix, "same_income") <- TRUE
  attr(JD_mix, "by_type") <- per_type
  JD_mix
}


#' Nonlinear perfect-foresight transition of a discount-factor MIXTURE
#'
#' The \code{omega}-weighted aggregate of the per-type nonlinear transitions
#' (\code{\link{hank_td_nonlinear}}). The types face the SAME aggregate price
#' paths \code{(r_path, w_path)}, interacting only through those prices, so
#' the mixture SCALAR aggregate paths are always exactly the
#' \code{omega}-weighted sums of the per-type paths,
#' \code{A^{mix}_t = sum_k omega_k A_{k,t}} (likewise \code{C}) --
#' regardless of whether the types share an income process.
#'
#' The mixture distribution PATH \code{Dpath^{mix}[, t] = sum_k omega_k
#' Dpath_k[, t]} (matching \code{\link{hank_mixture_dist}} at \code{t = 1})
#' additionally requires every type to live on the identical
#' \code{(a_grid, Pi, e)} cell space -- pooling cell-by-cell is otherwise
#' semantically invalid (cell \code{j} indexes a different income level per
#' type). When the blocks differ in their income process, this function does
#' NOT silently sum the per-type \code{Dpath} arrays: \code{Dpath = NULL} and
#' \code{Dpath_by_type} (length-K list of per-type \code{(n_e*n_a) x T_h}
#' paths) is returned instead, with \code{same_income = FALSE}. This is the
#' genuine nonlinear distribution truth against which the linear snapshot
#' reweighting \code{\link{hank_dist_response_snapshot}} is validated (when
#' \code{same_income} holds).
#'
#' @param blocks List of K \code{\link{hank_het_block}} objects sharing an
#'   asset grid (income process may differ across types).
#' @param omega Numeric length-K mixture weights (non-negative, summing to 1).
#' @param r_path,w_path Numeric length-\code{T_h} aggregate input paths (levels),
#'   applied in common to every type; missing entries default (per type) to that
#'   block's steady-state value. Intended use passes an explicit common path.
#' @param T_h Integer horizon (default from the path lengths).
#' @return A list with \code{A}, \code{C} (\code{omega}-weighted aggregate
#'   paths, always valid), \code{same_income} (logical), \code{Dpath}
#'   (\code{(n_e*n_a) x T_h} mixture distribution path when
#'   \code{same_income} is \code{TRUE}, else \code{NULL}), and
#'   \code{Dpath_by_type} (length-K list of per-type distribution paths;
#'   always populated).
#' @export
hank_mixture_td_nonlinear <- function(blocks, omega, r_path = NULL,
                                      w_path = NULL, T_h = NULL) {
  .hank_mixture_check_omega(blocks, omega)
  per <- lapply(blocks, hank_td_nonlinear, r_path = r_path, w_path = w_path,
                T_h = T_h)
  A <- Reduce(`+`, Map(function(p, wk) wk * p$A, per, omega))
  C <- Reduce(`+`, Map(function(p, wk) wk * p$C, per, omega))
  Dpath_by_type <- lapply(per, function(p) p$Dpath)
  same_income <- length(blocks) <= 1L || all(vapply(blocks[-1L], function(b)
    isTRUE(all.equal(b$Pi, blocks[[1L]]$Pi)) &&
      isTRUE(all.equal(b$e, blocks[[1L]]$e)), logical(1)))
  Dpath <- if (same_income)
    Reduce(`+`, Map(function(p, wk) wk * p$Dpath, per, omega)) else NULL
  if (!same_income)
    message("hank_mixture_td_nonlinear(): blocks differ in their income ",
            "process (Pi/e); a pooled cell-space Dpath is not defined ",
            "across differing income spaces, so Dpath = NULL. Use ",
            "'Dpath_by_type' (per-type paths) instead.")
  list(A = A, C = C, Dpath = Dpath, Dpath_by_type = Dpath_by_type,
       same_income = same_income)
}


#' Mixture AGGREGATE (scalar A/C) sequence-space Jacobian
#'
#' The omega-weighted sum of the per-type aggregate Jacobians
#' (\code{\link{hank_het_jacobian}}); companion to
#' \code{\link{hank_mixture_dist_jacobian}} for the scalar outputs. Exact by
#' the same linearity argument (shared grid, price-only interaction).
#'
#' @inheritParams hank_mixture_dist_jacobian
#' @param outputs Character subset of \code{c("A", "C")}.
#'
#' @return Nested list \code{J_mix[[output]][[input]]}, each a \code{T_h x
#'   T_h} matrix (identical layout to \code{\link{hank_het_jacobian}}'s
#'   return).
#' @export
hank_mixture_jacobian <- function(blocks, omega, T_h, inputs = c("r", "w"),
                                   outputs = c("A", "C")) {
  .hank_mixture_check_omega(blocks, omega)
  per_type <- lapply(blocks, function(b)
    hank_het_jacobian(b, T_h, inputs = inputs, outputs = outputs))
  setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) {
      acc <- omega[1L] * per_type[[1L]][[o]][[i]]
      if (length(blocks) > 1L)
        for (k in 2L:length(blocks)) acc <- acc + omega[k] * per_type[[k]][[o]][[i]]
      acc
    }), inputs)), outputs)
}


#' Quadratic reweighting loss: observed vs. model-implied net distribution change
#'
#' Given an OBSERVED net reweighting \code{dD_hat} (e.g. the difference of two
#' GREG-calibrated survey cross-sections; see
#' \code{\link{hank_phase0_reweight_gate}}), forms the model-implied net
#' reweighting from a candidate mixture distribution Jacobian \code{JD_mix}
#' and shock \code{dZ},
#' \deqn{dD^{mod}_t = \Sigma_i \Sigma_s JD_{mix,i}[t, s, ]\, dZ_i[s],}
#' and returns the (possibly date-weighted) quadratic loss between the two.
#'
#' \code{dD_hat} may be either a single length-\code{(n_e*n_a)} vector (the
#' single-snapshot-pair convention used by \code{hank_phase0_reweight_gate}:
#' compared against \code{dD_mod} SUMMED over all \code{T_h} dates, i.e. the
#' total accumulated displacement) or a full \code{(n_e*n_a) x T_h} matrix
#' (compared date-by-date against \code{dD_mod[, t]}), so the same loss
#' function serves both the Phase-0 harness (two snapshots) and any future
#' harness with a genuine per-date panel of reweightings.
#'
#' @param dD_hat Observed net reweighting: length-\code{(n_e*n_a)} vector (a
#'   single accumulated net change) or an \code{(n_e*n_a) x T_h} matrix (a
#'   per-date panel).
#' @param JD_mix Named list \code{JD_mix[[input]]} (as returned by
#'   \code{\link{hank_mixture_dist_jacobian}}), each a
#'   \code{T_h x T_h x (n_e*n_a)} array.
#' @param dZ Named list, one length-\code{T_h} shock path per entry of
#'   \code{names(JD_mix)} (missing entries treated as all-zero).
#' @param W Metric for the quadratic form: \code{NULL} (identity at every
#'   date), a single \code{(n_e*n_a) x (n_e*n_a)} matrix (used at every date /
#'   for the single-vector case), or a list of \code{T_h} such matrices (one
#'   per date, only meaningful when \code{dD_hat} is a matrix).
#'
#' @return A list with:
#'   \item{rho}{The scalar quadratic loss
#'     \code{Sum_t t(dD_hat_t - dD_mod_t) W_t (dD_hat_t - dD_mod_t)}
#'     (a single term when \code{dD_hat} is a vector).}
#'   \item{dD_mod}{The model-implied net reweighting, same shape as
#'     \code{dD_hat} (vector: summed over dates; matrix: per date).}
#'   \item{resid}{\code{dD_hat - dD_mod}, same shape.}
#' @export
hank_reweighting_loss <- function(dD_hat, JD_mix, dZ, W = NULL) {
  inputs <- names(JD_mix)
  T_h    <- dim(JD_mix[[1L]])[1L]
  n_cell <- dim(JD_mix[[1L]])[3L]

  ## Model-implied net reweighting at every date t: sum over inputs/shock dates.
  dD_mod_mat <- matrix(0, n_cell, T_h)
  for (i in inputs) {
    dz <- dZ[[i]]
    if (is.null(dz)) next
    JDi <- JD_mix[[i]]                       # T_h x T_h x n_cell
    for (tt in seq_len(T_h))
      dD_mod_mat[, tt] <- dD_mod_mat[, tt] +
        as.numeric(crossprod(matrix(JDi[tt, , ], nrow = T_h, ncol = n_cell), dz))
  }

  is_vec <- is.null(dim(dD_hat))
  if (is_vec) {
    if (length(dD_hat) != n_cell)
      stop(sprintf(
        "hank_reweighting_loss(): length(dD_hat) (%d) must equal n_cell (%d).",
        length(dD_hat), n_cell))
    dD_mod <- rowSums(dD_mod_mat)             # total accumulated displacement
    resid  <- dD_hat - dD_mod
    Wm <- if (is.null(W)) diag(n_cell) else W
    rho <- as.numeric(crossprod(resid, Wm %*% resid))
  } else {
    if (!identical(dim(dD_hat), c(n_cell, T_h)))
      stop(sprintf(
        "hank_reweighting_loss(): dim(dD_hat) (%s) must equal c(n_cell, T_h) = c(%d, %d).",
        paste(dim(dD_hat), collapse = ", "), n_cell, T_h))
    dD_mod <- dD_mod_mat
    resid  <- dD_hat - dD_mod
    rho <- 0
    for (tt in seq_len(T_h)) {
      Wt <- if (is.null(W)) diag(n_cell) else if (is.list(W)) W[[tt]] else W
      rho <- rho + as.numeric(crossprod(resid[, tt], Wt %*% resid[, tt]))
    }
  }

  list(rho = rho, dD_mod = dD_mod, resid = resid)
}


#' Phase-0 reweighting-gate harness: net-reweighting loss over a (centre,
#' spread) grid of candidate 2-type discount-factor mixtures
#'
#' Assembles the full Phase-0 synthetic-survey experiment for the
#' het-preferences-in-HANK research line (brief-17): (1) a TRUTH 2-type
#' mixture with discount factors \code{centre* -+ spread*}; (2) its
#' first-order net distribution reweighting under a small aggregate shock;
#' (3) two synthetic household surveys (pre- and post-shock) drawn from the
#' truth distributions with a deliberately mis-specified design weight, then
#' GREG-calibrated to the TRUE mixture macro totals; (4) a grid of candidate
#' \code{(centre, spread)} 2-type mixtures, each scored against the observed
#' survey-implied net reweighting by \code{\link{hank_reweighting_loss}}, plus
#' an aggregate-only comparator using \code{A}, \code{C} macro totals alone.
#'
#' This function performs the MECHANICS only: it returns the raw grids of
#' loss values (\code{rho}, \code{rho_agg}, \code{rho_identityW}), the truth
#' calibration, and each grid's argmin. It does NOT compute or assert a
#' verdict on whether reweighting "identifies" the spread -- that judgment is
#' left to the caller.
#'
#' \strong{Truth parameterization}: with \code{centre* = mean(betas)} and
#' \code{spread* = diff(betas) / 2} (so \code{betas = c(centre*-spread*,
#' centre*+spread*)}), the truth 2-type discount factors are
#' \code{beta_1* = centre* - spread*}, \code{beta_2* = centre* + spread*}.
#'
#' \strong{Date convention} (see the file header comment): the synthetic
#' survey draws exactly two cross-sections -- one from the pre-shock
#' stationary mixture distribution \code{D0} and one from the (first-order)
#' post-shock distribution \code{D1 = D0 + dD_star}, where \code{dD_star} is
#' the TOTAL accumulated distribution displacement over the horizon
#' \code{T_h} (i.e. \code{Sum_t dD_star[, t]} from the truth distribution
#' Jacobian). The observed net reweighting \code{dD_hat} is therefore a
#' SINGLE realized \code{(n_e*n_a)}-vector (calibrated \code{Dhat1 - Dhat0}),
#' compared by \code{\link{hank_reweighting_loss}} against each candidate
#' mixture's model-implied net reweighting, likewise summed over \code{T_h}.
#'
#' @param betas Length-2 truth discount factors; \code{centre* = mean(betas)},
#'   \code{spread* = diff(betas)/2} (see Details).
#' @param omega Length-2 truth mixture weights (non-negative, sum to 1).
#' @param n_e,n_a,amax Household grid size/extent (shared truth + candidate
#'   grid; small by default for a fast Phase-0 run).
#' @param r,w Steady-state prices (shared, exogenous to the mixture -- the
#'   harness treats \code{r}, \code{w} as fixed calibration targets, not a
#'   solved GE; a full GE closure is out of scope for this plumbing pass).
#' @param T_h Integer horizon for the Jacobians / shock path.
#' @param shock Named list of length-\code{T_h} shock paths (subset of
#'   \code{c("r","w")}); default a small AR(1)-decaying interest-rate shock.
#' @param N Sample size drawn for EACH synthetic survey (pre- and post-shock).
#' @param centre_grid,spread_grid Numeric vectors: the grid of candidate
#'   \code{centre}/\code{spread} mixture calibrations to score.
#' @param seed RNG seed for the synthetic sampling (design weights + draws).
#'
#' @return A list with:
#'   \item{centre_grid,spread_grid}{The grids as passed in.}
#'   \item{rho}{\code{length(centre_grid) x length(spread_grid)} matrix: the
#'     GREG-variance-weighted reweighting loss at every candidate.}
#'   \item{rho_identityW}{Same shape, identity-metric version (\code{W = I}).}
#'   \item{rho_agg}{Same shape: squared error of the candidate's
#'     \code{(A, C)} net response against the calibrated aggregate net
#'     change (the "aggregates-only" comparator).}
#'   \item{centre_star,spread_star}{The truth calibration
#'     (\code{mean(betas)}, \code{diff(betas)/2}).}
#'   \item{argmin_rho,argmin_rho_identityW,argmin_rho_agg}{Each a list
#'     \code{list(centre, spread)} at the grid argmin of the corresponding
#'     loss surface.}
#'   \item{dD_hat,dD_star}{The observed (survey-calibrated) and truth
#'     (model, noiseless) net reweighting vectors, for direct inspection.}
#' @export
hank_phase0_reweight_gate <- function(betas = c(0.95, 0.98),
                                       omega = c(0.5, 0.5),
                                       n_e = 3L, n_a = 50L, amax = 60,
                                       r = 0.01, w = 1.0,
                                       T_h = 20L,
                                       shock = list(r = 0.01 * 0.8^(seq_len(T_h) - 1L)),
                                       N = 5000L,
                                       centre_grid = seq(0.93, 0.99, length.out = 7L),
                                       spread_grid = seq(0.005, 0.035, length.out = 7L),
                                       seed = 1L) {
  if (length(betas) != 2L) stop("hank_phase0_reweight_gate(): 'betas' must have length 2.")
  set.seed(seed)

  inc <- hank_income_rouwenhorst(rho = 0.9, sigma = 0.7, n = n_e)
  ag  <- hank_asset_grid(amax = amax, n = n_a, amin = 0)
  shock_inputs <- names(shock)

  ## ---- 1-3. Truth mixture + synthetic surveys + GREG calibration -------
  ## Factored into .hank_phase0_survey() (verbatim extraction, no numerical
  ## change -- see that function's header) so hank_phase1_joint_gate() can
  ## reuse the identical truth/survey/dD_hat/W construction below.
  sv <- .hank_phase0_survey(ag, inc, betas, omega, eis = 1, r = r, w = w,
                             T_h = T_h, shock = shock, N = N)

  ## ---- 4. Grid of candidate (centre, spread) mixtures ------------------
  n_c <- length(centre_grid); n_s <- length(spread_grid)
  rho_mat      <- matrix(NA_real_, n_c, n_s)
  rho_id_mat   <- matrix(NA_real_, n_c, n_s)
  rho_agg_mat  <- matrix(NA_real_, n_c, n_s)

  for (ci in seq_len(n_c)) for (si in seq_len(n_s)) {
    cc <- centre_grid[ci]; ss <- spread_grid[si]
    betas_cand <- c(cc - ss, cc + ss)
    blocks_cand <- hank_mixture_blocks(ag, inc$Pi, inc$e, betas = betas_cand,
                                        eis = 1, r = r, w = w)
    JD_cand <- hank_mixture_dist_jacobian(blocks_cand, omega, T_h, inputs = shock_inputs)
    J_cand  <- hank_mixture_jacobian(blocks_cand, omega, T_h, inputs = shock_inputs,
                                      outputs = c("A", "C"))

    loss_greg <- hank_reweighting_loss(sv$dD_hat, JD_cand, shock, W = sv$W_greg)
    loss_id   <- hank_reweighting_loss(sv$dD_hat, JD_cand, shock, W = NULL)
    rho_mat[ci, si]    <- loss_greg$rho
    rho_id_mat[ci, si] <- loss_id$rho

    dAC_cand <- vapply(c("A", "C"), function(o) {
      acc <- 0
      for (i in shock_inputs) acc <- acc + sum(J_cand[[o]][[i]] %*% shock[[i]])
      acc
    }, numeric(1))
    rho_agg_mat[ci, si] <- sum((sv$dAC_hat - dAC_cand)^2)
  }

  argmin_idx <- function(M) {
    idx <- which(M == min(M), arr.ind = TRUE)[1L, ]
    list(centre = centre_grid[idx[1L]], spread = spread_grid[idx[2L]])
  }

  list(centre_grid = centre_grid, spread_grid = spread_grid,
       rho = rho_mat, rho_identityW = rho_id_mat, rho_agg = rho_agg_mat,
       centre_star = sv$centre_star, spread_star = sv$spread_star,
       argmin_rho = argmin_idx(rho_mat),
       argmin_rho_identityW = argmin_idx(rho_id_mat),
       argmin_rho_agg = argmin_idx(rho_agg_mat),
       dD_hat = sv$dD_hat, dD_star = sv$dD_star,
       dAC_hat = sv$dAC_hat, dAC_star = sv$dAC_star)
}


#' Phase-0 truth mixture + synthetic-survey + GREG-calibration construction
#' (internal, shared by \code{\link{hank_phase0_reweight_gate}} and
#' \code{\link{hank_phase1_joint_gate}})
#'
#' VERBATIM extraction of steps 1-3 of \code{hank_phase0_reweight_gate}'s body
#' (truth 2-type mixture -> first-order net distribution reweighting ->
#' synthetic pre-/post-shock surveys -> GREG calibration -> calibrated net
#' reweighting \code{dD_hat} and its GREG-variance metric \code{W_greg}): NO
#' numerical logic was changed when factoring this out, only the surrounding
#' function boundary. Callers must have already called \code{set.seed()} (the
#' caller owns the RNG seed, as before).
#'
#' @param ag Numeric asset grid (\code{\link{hank_asset_grid}} output).
#' @param inc List with \code{Pi} (income transition matrix) and \code{e}
#'   (income levels), as returned by \code{\link{hank_income_rouwenhorst}}.
#' @param betas,omega,eis,r,w,T_h,shock As in
#'   \code{\link{hank_phase0_reweight_gate}}.
#' @param N Sample size drawn for EACH synthetic survey (pre- and post-shock).
#'
#' @return A list with \code{centre_star, spread_star, truth_blocks, D0, D1,
#'   dD_star, dAC_star, JD_star, J_star, dD_hat, dAC_hat, W_greg,
#'   a_grid_cell, c_cell_mat, n_cell} -- every quantity
#'   \code{hank_phase0_reweight_gate}'s grid-search step (and
#'   \code{hank_phase1_joint_gate}'s) reads off this helper.
#' @keywords internal
.hank_phase0_survey <- function(ag, inc, betas, omega, eis = 1, r, w,
                                 T_h, shock, N = 5000L) {
  if (length(betas) != 2L) stop(".hank_phase0_survey(): 'betas' must have length 2.")
  n_e <- length(inc$e); n_a <- length(ag)
  n_cell <- n_e * n_a
  shock_inputs <- names(shock)

  ## ---- 1. Truth mixture -----------------------------------------------
  centre_star <- mean(betas)
  spread_star <- diff(betas) / 2
  truth_blocks <- hank_mixture_blocks(ag, inc$Pi, inc$e, betas = betas,
                                       eis = eis, r = r, w = w)
  truth_dist   <- hank_mixture_dist(truth_blocks, omega)
  D0 <- truth_dist$D
  JD_star <- hank_mixture_dist_jacobian(truth_blocks, omega, T_h, inputs = shock_inputs)
  J_star  <- hank_mixture_jacobian(truth_blocks, omega, T_h, inputs = shock_inputs,
                                    outputs = c("A", "C"))

  ## ---- 2. Model-implied net reweighting at truth (first-order) --------
  ## dD_star[, t] = Sum_i Sum_s JD_star[[i]][t, s, ] * shock[[i]][s]; summed
  ## over t = 1..T_h to get the TOTAL accumulated displacement (see Date
  ## convention in hank_phase0_reweight_gate's roxygen above).
  dD_star_mat <- matrix(0, n_cell, T_h)
  for (i in shock_inputs) {
    dz  <- shock[[i]]
    JDi <- JD_star[[i]]
    for (tt in seq_len(T_h))
      dD_star_mat[, tt] <- dD_star_mat[, tt] +
        as.numeric(crossprod(matrix(JDi[tt, , ], nrow = T_h, ncol = n_cell), dz))
  }
  dD_star <- rowSums(dD_star_mat)
  D1 <- pmax(D0 + dD_star, 0)
  D1 <- D1 / sum(D1)                      # renormalize for sampling only

  ## Truth aggregate (A, C) net response, summed over the horizon (same
  ## accumulation convention as dD_star, for consistency with rho_agg below).
  dAC_star <- vapply(c("A", "C"), function(o) {
    acc <- 0
    for (i in shock_inputs) acc <- acc + sum(J_star[[o]][[i]] %*% shock[[i]])
    acc
  }, numeric(1))

  ## ---- 3. Synthetic surveys: draw N units each from D0 and D1 ----------
  ## Distribution cell order is ASSET-FAST, income-slow: the het-block policy
  ## matrices are n_e x n_a and .hank_mat_to_vec() row-major-flattens them, so
  ## cell c has asset ((c-1) %% n_a)+1 and income floor((c-1)/n_a)+1 -- i.e.
  ## a_cell = rep(ag, times = n_e), e_cell = rep(e, each = n_a). (A stopifnot
  ## below LOCKS this against hank_mixture_dist()'s own A, C.)
  a_grid_cell <- rep(ag, times = n_e)         # cell -> asset level (distribution order)
  e_cell      <- rep(inc$e, each = n_a)       # cell -> income level (distribution order)

  draw_survey <- function(Dprob, n) {
    cells <- sample.int(n_cell, size = n, replace = TRUE, prob = Dprob)
    ## Deliberately mis-set design weight: base 1, times a unit-level
    ## multiplicative noise factor -- something GREG calibration must correct.
    d_raw <- rep(1, n) * (1 + 0.5 * stats::runif(n))
    list(cells = cells, d = d_raw)
  }
  survey0 <- draw_survey(D0, N)
  survey1 <- draw_survey(D1, N)

  ## Auxiliary matrix Z: per-cell indicator design (N x n_cell), so control
  ## totals X are exactly the TRUE mixture macro cell-mass*[asset,income,1]
  ## moments -- here we calibrate to (aggregate assets A, aggregate
  ## consumption C, total mass), the 3 truth macro totals named in the brief.
  ## Build per-cell (a, e*a-weighted-consumption-proxy) contributions from the
  ## TRUTH mixture policy (consumption at steady state r,w; first-order
  ## accurate for a small shock) so X is expressed in the same per-unit basis
  ## as Z's columns.
  c_cell_mat <- Reduce(`+`, Map(function(b, wk) wk * .hank_mat_to_vec(b$c),
                                 truth_blocks, omega))    # mixture ss consumption per cell

  ## LOCK the asset-cell-label order: assets are a COMMON state coordinate, so
  ## the mixture aggregate satisfies exactly A = Sum_k omega_k A_k =
  ## sum(D0 * a_grid_cell). Reproducing hank_mixture_dist()'s independently-
  ## computed A therefore catches any future asset/income transpose of
  ## a_grid_cell (the actual cell order must be asset-fast; see the comment
  ## above). (No analogous C lock: consumption is a per-TYPE policy, so
  ## sum(D0 * c_cell_mat) is the survey plug-in E_cell[mixture-avg c], which by
  ## design differs from the true mixture C = Sum_k omega_k C_k by the
  ## cross-type terms omega_j omega_k Sum_cell D_j c_k.)
  stopifnot(abs(sum(D0 * a_grid_cell) - truth_dist$A) < 1e-8 * (abs(truth_dist$A) + 1))

  build_ZX <- function(cells, target_D) {
    n <- length(cells)
    Z <- cbind(asset = a_grid_cell[cells],
               cons  = c_cell_mat[cells],
               mass  = 1)
    ## Control totals under the TRUE mixture distribution at this date:
    ## Sum_cell D_cell * (a_cell, c_cell, 1) = (A, C, 1).
    X <- c(asset = sum(target_D * a_grid_cell),
           cons  = sum(target_D * c_cell_mat),
           mass  = 1)
    list(Z = Z, X = X)
  }
  zx0 <- build_ZX(survey0$cells, D0)
  zx1 <- build_ZX(survey1$cells, D1)

  cal0 <- hank_calibrate_weights(survey0$d, zx0$Z, zx0$X, method = "greg")
  cal1 <- hank_calibrate_weights(survey1$d, zx1$Z, zx1$X, method = "greg")

  ## Calibrated weighted CELL distributions: sum of calibrated weight mass
  ## falling in each cell, normalized to sum to 1.
  cell_dist_from_cal <- function(cells, w_cal) {
    d <- numeric(n_cell)
    agg <- tapply(w_cal, cells, sum)
    d[as.integer(names(agg))] <- as.numeric(agg)
    d / sum(d)
  }
  Dhat0 <- cell_dist_from_cal(survey0$cells, cal0$w)
  Dhat1 <- cell_dist_from_cal(survey1$cells, cal1$w)
  dD_hat <- Dhat1 - Dhat0

  ## Calibrated aggregate (A, C) net change, from the same surveys' totals
  ## (the "aggregates-only" observed comparator for rho_agg).
  A_hat0 <- sum(Dhat0 * a_grid_cell); A_hat1 <- sum(Dhat1 * a_grid_cell)
  C_hat0 <- sum(Dhat0 * c_cell_mat);  C_hat1 <- sum(Dhat1 * c_cell_mat)
  dAC_hat <- c(A = A_hat1 - A_hat0, C = C_hat1 - C_hat0)

  ## GREG-variance metric: per-cell inverse-variance diagonal built from the
  ## post-shock calibration's vhat is NOT directly on the cell basis (vhat is
  ## p x p over the (asset, cons, mass) calibration variables, not n_cell x
  ## n_cell) -- project it onto the cell basis via the per-cell Z row so the
  ## metric is usable as an (n_cell x n_cell) diagonal weight, W_cell =
  ## diag(1 / (t(z_cell) vhat z_cell)), a per-cell scalar precision derived from
  ## the calibration-variable variance at that cell's own (asset, cons, 1)
  ## coordinates (a plug-in delta-method-style projection of the p x p
  ## variance handle down to the cell it directly informs).
  Zcell <- cbind(asset = a_grid_cell, cons = c_cell_mat, mass = 1)
  vhat_avg <- (cal0$vhat + cal1$vhat) / 2
  cell_var <- vapply(seq_len(n_cell), function(cc) {
    z <- Zcell[cc, ]
    as.numeric(z %*% vhat_avg %*% z)
  }, numeric(1))
  cell_var <- pmax(cell_var, 1e-12 * max(cell_var))   # floor against exact zeros
  W_greg <- diag(1 / cell_var)

  list(centre_star = centre_star, spread_star = spread_star,
       truth_blocks = truth_blocks, D0 = D0, D1 = D1,
       dD_star = dD_star, dAC_star = dAC_star,
       JD_star = JD_star, J_star = J_star,
       dD_hat = dD_hat, dAC_hat = dAC_hat, W_greg = W_greg,
       a_grid_cell = a_grid_cell, c_cell_mat = c_cell_mat, n_cell = n_cell)
}


## --------------------------------------------------------------------------
## Phase-1 Wave-1: JOINT macro + reweighting recovery.
##
## Phase-0 (above) showed that the model-implied net cross-section
## reweighting identifies the beta-mixture SPREAD that aggregate (A, C)
## moments alone confound. Its caveat: reweighting ALONE is biased in the
## beta-CENTRE, because the first-order aggregate IRF (which the net
## reweighting is built from, via the shock path) confounds centre vs.
## spread just as the aggregates do. Phase-1 closes that gap with a JOINT
## objective: a noisy macro time-path (A_t, C_t) pins the centre (the
## aggregate IRF is centre-sensitive, spread-insensitive at leading order),
## while the net reweighting pins the spread (as in Phase-0), and the sum of
## the two losses recovers BOTH.
## --------------------------------------------------------------------------


#' Proper multinomial-difference GLS precision for the net reweighting
#'
#' The observed net reweighting \code{dD_hat = Dhat1 - Dhat0} is a difference
#' of two INDEPENDENT size-\code{N} multinomial empirical cell-distributions
#' (one drawn from the pre-shock cross-section, one from the post-shock
#' cross-section). At leading order (each \code{Dhat_j[cell]} is itself a
#' sample proportion with the usual multinomial variance
#' \code{D_j[cell](1-D_j[cell])/N}, and the two cross-sections are
#' independent draws so the variance of their difference adds),
#' \deqn{Var(dD_hat[cell]) = (D0[cell](1 - D0[cell]) + D1[cell](1 -
#' D1[cell])) / N.}
#' This function returns the diagonal GLS precision built from that leading
#' variance, \code{W = diag(N / v)}, \code{v = D0(1-D0) + D1(1-D1)}.
#'
#' Cells with negligible support under EITHER truth cross-section (v below
#' \code{floor_frac * max(v)}) get weight \strong{0} rather than an
#' ill-conditioned huge weight from a near-zero variance -- i.e. unsupported
#' or degenerate cells are DROPPED from the quadratic form, not given a
#' \code{1/0} blow-up.
#'
#' This is the LEADING-ORDER multinomial-difference variance only. The
#' further variance reduction from GREG calibration to 3 auxiliary control
#' totals (asset/consumption/mass; see \code{\link{hank_phase0_reweight_gate}}'s
#' \code{W_greg}) is a SECOND-ORDER refinement on top of this leading term
#' (calibration can only reduce variance relative to the uncalibrated
#' multinomial baseline) and is \strong{not} included here -- this function
#' gives the plain-vanilla two-independent-multinomials metric, useful as a
#' baseline/sanity metric independent of any particular calibration design.
#'
#' @param D0,D1 Numeric length-\code{n_cell} pre-/post-shock TRUTH cell
#'   probability distributions (each non-negative, summing to 1; e.g. the
#'   \code{D0}, \code{D1} objects computed inside
#'   \code{\link{hank_phase0_reweight_gate}} / \code{\link{.hank_phase0_survey}}).
#' @param N Sample size of EACH of the two cross-sections (as in
#'   \code{\link{hank_phase0_reweight_gate}}'s \code{N} argument; both surveys
#'   are assumed the same size here).
#' @param floor_frac Cells with \code{v[cell] < floor_frac * max(v)} are
#'   dropped (weight set to 0) rather than inverted. Default \code{1e-6}.
#'
#' @return A list with:
#'   \item{W}{The \code{n_cell x n_cell} diagonal precision matrix
#'     \code{diag(w_diag)}.}
#'   \item{w_diag}{The length-\code{n_cell} diagonal itself.}
#'   \item{dropped}{Integer indices of cells set to weight 0 by the floor.}
#' @export
hank_reweight_metric <- function(D0, D1, N, floor_frac = 1e-6) {
  if (length(D0) != length(D1))
    stop(sprintf(
      "hank_reweight_metric(): length(D0) (%d) must equal length(D1) (%d).",
      length(D0), length(D1)))
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0)
    stop("hank_reweight_metric(): 'N' must be a single finite positive number.")

  v <- D0 * (1 - D0) + D1 * (1 - D1)
  dropped <- which(v < floor_frac * max(v))
  w_diag <- N / v
  w_diag[dropped] <- 0

  list(W = diag(w_diag), w_diag = w_diag, dropped = dropped)
}


#' Mixture aggregate IRF PATH: (A_t, C_t) response to a shock
#'
#' Convenience wrapper around \code{\link{hank_mixture_jacobian}}: contracts
#' the mixture aggregate Jacobian \code{J_mix[[o]][[i]]} (a \code{T_h x T_h}
#' matrix) against a shock path to return the actual length-\code{T_h}
#' aggregate response \code{path[[o]]}, rather than the Jacobian itself,
#' \deqn{path[[o]][t] = \Sigma_i \Sigma_s J_{mix}[[o]][[i]][t, s]\, shock[[i]][s].}
#'
#' @inheritParams hank_mixture_jacobian
#' @param shock Named list, one length-\code{T_h} shock path per entry of
#'   \code{inputs} (missing entries treated as all-zero; names not in
#'   \code{inputs} are ignored, mirroring \code{inputs = names(shock)} usage
#'   elsewhere in this file).
#'
#' @return Named list \code{path[[output]]}, each a length-\code{T_h} numeric
#'   vector: the aggregate response path.
#' @export
hank_mixture_agg_irf <- function(blocks, omega, T_h, shock,
                                  outputs = c("A", "C")) {
  inputs <- names(shock)
  J_mix <- hank_mixture_jacobian(blocks, omega, T_h, inputs = inputs, outputs = outputs)
  setNames(lapply(outputs, function(o) {
    acc <- numeric(T_h)
    for (i in inputs) {
      dz <- shock[[i]]
      if (is.null(dz)) next
      acc <- acc + as.numeric(J_mix[[o]][[i]] %*% dz)
    }
    acc
  }), outputs)
}


#' Aggregate IRF-path moment loss: observed macro path vs. candidate mixture
#'
#' Given an OBSERVED aggregate path \code{AC_hat} (one length-\code{T_h}
#' series per output, e.g. a noisy \code{(A_t, C_t)} macro time series) and a
#' candidate mixture, computes the candidate's implied path via
#' \code{\link{hank_mixture_agg_irf}} and returns the (possibly
#' per-output-scaled) sum-of-squares moment loss
#' \deqn{\rho_{macro} = \Sigma_o \Sigma_t \left( \frac{AC\_hat[[o]][t] -
#' AC\_cand[[o]][t]}{s_o} \right)^2.}
#'
#' @param AC_hat Named list, one length-\code{T_h} observed aggregate path per
#'   output (names determine \code{outputs} passed to
#'   \code{\link{hank_mixture_agg_irf}}).
#' @param blocks_cand List of candidate \code{\link{hank_het_block}} objects
#'   sharing a grid (as returned by \code{\link{hank_mixture_blocks}}).
#' @param omega,T_h,shock As in \code{\link{hank_mixture_agg_irf}}.
#' @param sigma_macro Named list/vector of per-output macro noise standard
#'   deviations \code{s_o} (matched by name to \code{AC_hat}); \code{NULL}
#'   (default) uses \code{s_o = 1} for every output, i.e. an unscaled sum of
#'   squares.
#'
#' @return A list with:
#'   \item{rho_macro}{The scalar aggregate-path moment loss.}
#'   \item{AC_cand}{The candidate's implied path (same shape as
#'     \code{AC_hat}), for direct inspection.}
#' @export
hank_macro_loss <- function(AC_hat, blocks_cand, omega, T_h, shock,
                             sigma_macro = NULL) {
  outputs <- names(AC_hat)
  AC_cand <- hank_mixture_agg_irf(blocks_cand, omega, T_h, shock, outputs = outputs)

  rho_macro <- 0
  for (o in outputs) {
    s_o <- if (is.null(sigma_macro)) 1 else sigma_macro[[o]]
    rho_macro <- rho_macro + sum(((AC_hat[[o]] - AC_cand[[o]]) / s_o)^2)
  }

  list(rho_macro = rho_macro, AC_cand = AC_cand)
}


#' Phase-1 joint macro + reweighting recovery demo: a JOINT (centre, spread)
#' grid search recovers BOTH mixture parameters
#'
#' Assembles the Phase-1 Wave-1 synthetic experiment: the SAME truth 2-type
#' mixture, synthetic-survey net-reweighting construction and GREG
#' calibration as \code{\link{hank_phase0_reweight_gate}} (reused via the
#' shared internal \code{.hank_phase0_survey()} helper, so the reweighting
#' side of this gate is numerically identical to Phase-0's), PLUS a synthetic
#' noisy MACRO aggregate time-path \code{(A_t, C_t)}. Over a
#' \code{(centre, spread)} grid of candidate 2-type mixtures it computes THREE
#' loss surfaces:
#' \itemize{
#'   \item \code{rho_macro} -- \code{\link{hank_macro_loss}} against the noisy
#'     macro path alone (expected to pin the CENTRE; the aggregate IRF is far
#'     more sensitive to the mixture's centre than to its spread at this
#'     order).
#'   \item \code{rho_rw} -- \code{\link{hank_reweighting_loss}} against the
#'     survey-calibrated net reweighting \code{dD_hat}, using the PROPER
#'     multinomial-difference metric \code{\link{hank_reweight_metric}}
#'     (expected to pin the SPREAD, as in Phase-0, but now with the
#'     statistically appropriate GLS weight rather than the GREG
#'     delta-method projection).
#'   \item \code{joint = rho_macro + lambda * rho_rw} -- the combined
#'     objective, expected to recover BOTH parameters where either loss alone
#'     is degenerate along one direction of the grid.
#' }
#' This function performs the MECHANICS only -- it returns the three raw grid
#' surfaces and their argmins; it does NOT compute or assert a verdict on
#' whether the joint objective "identifies" the truth (that judgment is left
#' to the caller, per the convention of \code{\link{hank_phase0_reweight_gate}}).
#'
#' \strong{Truth parameterization} (as in
#' \code{\link{hank_phase0_reweight_gate}}): with \code{centre* = mean(betas)}
#' and \code{spread* = diff(betas)/2}, the truth 2-type discount factors are
#' \code{beta_1* = centre* - spread*}, \code{beta_2* = centre* + spread*}.
#'
#' \strong{Macro observation}: \code{AC_hat_macro[[o]][t] = AC_star[[o]][t] +
#' rnorm(T_h, 0, s_o)}, with \code{AC_star} the TRUTH aggregate path
#' (\code{\link{hank_mixture_agg_irf}} at the truth mixture) and \code{s_o =
#' sigma_macro_frac * sd(AC_star[[o]])} (floored, see \code{Details}) -- i.e.
#' a fixed FRACTION of the truth path's own variation, a simple
#' scale-invariant way to inject "the econometrician observes this macro
#' series with some noise" without hand-picking absolute noise units per
#' output.
#'
#' @param betas Length-2 truth discount factors; \code{centre* = mean(betas)},
#'   \code{spread* = diff(betas)/2} (see Details).
#' @param omega Length-2 truth mixture weights (non-negative, sum to 1).
#' @param n_e,n_a,amax Household grid size/extent (shared truth + candidate
#'   grid; small by default for a fast run).
#' @param r,w Steady-state prices (shared, exogenous to the mixture; see
#'   \code{\link{hank_phase0_reweight_gate}}'s Details for the same caveat).
#' @param T_h Integer horizon for the Jacobians / shock / macro path.
#' @param shock Named list of length-\code{T_h} shock paths (subset of
#'   \code{c("r","w")}); default a small AR(1)-decaying interest-rate shock.
#' @param N Sample size drawn for EACH synthetic survey (pre- and post-shock;
#'   passed through to \code{.hank_phase0_survey()}).
#' @param sigma_macro_frac Macro noise scale as a FRACTION of each truth
#'   output path's own standard deviation (see Details); \code{0} gives a
#'   noiseless macro observation (\code{AC_hat_macro == AC_star} exactly),
#'   used by the noiseless-recovery test.
#' @param centre_grid,spread_grid Numeric vectors: the grid of candidate
#'   \code{centre}/\code{spread} mixture calibrations to score.
#' @param lambda Non-negative scalar weight on the reweighting loss in the
#'   joint objective \code{joint = rho_macro + lambda * rho_rw}.
#' @param seed RNG seed for the synthetic sampling (design weights, survey
#'   draws, and macro noise).
#'
#' @return A list with:
#'   \item{centre_grid,spread_grid}{The grids as passed in.}
#'   \item{rho_macro,rho_rw,joint}{Each a
#'     \code{length(centre_grid) x length(spread_grid)} matrix.}
#'   \item{truth}{\code{list(centre, spread)} at the truth calibration.}
#'   \item{argmin_macro,argmin_rw,argmin_joint}{Each a list
#'     \code{list(centre, spread)} at the grid argmin of the corresponding
#'     loss surface.}
#' @export
hank_phase1_joint_gate <- function(betas = c(0.95, 0.98),
                                    omega = c(0.5, 0.5),
                                    n_e = 3L, n_a = 50L, amax = 60,
                                    r = 0.01, w = 1.0,
                                    T_h = 20L,
                                    shock = list(r = 0.01 * 0.8^(seq_len(T_h) - 1L)),
                                    N = 5000L,
                                    sigma_macro_frac = 0.05,
                                    centre_grid, spread_grid,
                                    lambda = 1,
                                    seed = 1L) {
  if (length(betas) != 2L) stop("hank_phase1_joint_gate(): 'betas' must have length 2.")
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
    stop("hank_phase1_joint_gate(): 'lambda' must be a single non-negative finite number.")
  set.seed(seed)

  inc <- hank_income_rouwenhorst(rho = 0.9, sigma = 0.7, n = n_e)
  ag  <- hank_asset_grid(amax = amax, n = n_a, amin = 0)
  shock_inputs <- names(shock)
  outputs <- c("A", "C")

  ## ---- a. Truth mixture + reweighting side (identical construction to
  ## hank_phase0_reweight_gate, via the shared internal helper) -----------
  sv <- .hank_phase0_survey(ag, inc, betas, omega, eis = 1, r = r, w = w,
                             T_h = T_h, shock = shock, N = N)

  ## Truth aggregate path AC_star, and the proper multinomial-difference
  ## metric W_proper built from the two TRUTH cross-sections D0, D1 (not the
  ## GREG delta-method W_greg used by Phase-0 -- this gate uses the leading
  ## multinomial-variance metric per hank_reweight_metric()'s docs).
  AC_star  <- hank_mixture_agg_irf(sv$truth_blocks, omega, T_h, shock, outputs = outputs)
  W_proper <- hank_reweight_metric(sv$D0, sv$D1, N = N)$W

  ## ---- b. Synthetic MACRO observation: noisy (A_t, C_t) -----------------
  s_macro <- setNames(vapply(outputs, function(o) {
    s <- sigma_macro_frac * stats::sd(AC_star[[o]])
    max(s, 1e-10)                      # floor against a numerically-flat path
  }, numeric(1)), outputs)

  AC_hat_macro <- setNames(lapply(outputs, function(o) {
    if (sigma_macro_frac <= 0) return(AC_star[[o]])       # noiseless observation
    AC_star[[o]] + stats::rnorm(T_h, mean = 0, sd = s_macro[[o]])
  }), outputs)

  ## ---- c. Grid of candidate (centre, spread) mixtures -------------------
  n_c <- length(centre_grid); n_s <- length(spread_grid)
  rho_macro_mat <- matrix(NA_real_, n_c, n_s)
  rho_rw_mat    <- matrix(NA_real_, n_c, n_s)
  joint_mat     <- matrix(NA_real_, n_c, n_s)

  for (ci in seq_len(n_c)) for (si in seq_len(n_s)) {
    cc <- centre_grid[ci]; ss <- spread_grid[si]
    betas_cand <- c(cc - ss, cc + ss)
    blocks_cand <- hank_mixture_blocks(ag, inc$Pi, inc$e, betas = betas_cand,
                                        eis = 1, r = r, w = w)
    JD_cand <- hank_mixture_dist_jacobian(blocks_cand, omega, T_h, inputs = shock_inputs)

    loss_macro <- hank_macro_loss(AC_hat_macro, blocks_cand, omega, T_h, shock,
                                   sigma_macro = s_macro)
    loss_rw    <- hank_reweighting_loss(sv$dD_hat, JD_cand, shock, W = W_proper)

    rho_macro_mat[ci, si] <- loss_macro$rho_macro
    rho_rw_mat[ci, si]    <- loss_rw$rho
    joint_mat[ci, si]     <- loss_macro$rho_macro + lambda * loss_rw$rho
  }

  argmin_idx <- function(M) {
    idx <- which(M == min(M), arr.ind = TRUE)[1L, ]
    list(centre = centre_grid[idx[1L]], spread = spread_grid[idx[2L]])
  }

  list(centre_grid = centre_grid, spread_grid = spread_grid,
       rho_macro = rho_macro_mat, rho_rw = rho_rw_mat, joint = joint_mat,
       truth = list(centre = sv$centre_star, spread = sv$spread_star),
       argmin_macro = argmin_idx(rho_macro_mat),
       argmin_rw    = argmin_idx(rho_rw_mat),
       argmin_joint = argmin_idx(joint_mat))
}


## --------------------------------------------------------------------------
## GRID-INVARIANT functional reweighting loss.
##
## The cell-wise reweighting metric hank_reweight_metric() (diagonal
## N / (D0(1-D0) + D1(1-D1))) is GRID-PATHOLOGICAL: as the asset grid refines,
## each cell's mass shrinks, the diagonal treats every shrinking cell as an
## independent observation, and the implied Fisher information DIVERGES with
## n_a (empirically min-eig 9.2e7 -> 2.4e8 -> 9.3e8 as n_a = 50 -> 100 -> 200).
## It ignores the multinomial NEGATIVE correlations between neighbouring cells,
## which dominate as cells shrink -- so any information / conditioning /
## point-estimate quantity read off that diagonal metric is grid-dependent and
## not comparable across resolutions.
##
## The fix is to compare a MODERATE, grid-invariant set of smooth distribution
## FUNCTIONALS -- low-order orthonormal polynomials in a bounded wealth
## coordinate -- under their PROPER (full K x K) multinomial-difference
## covariance. Projecting the net reweighting onto K smooth test functions is a
## low-dimensional regularisation whose information converges as the grid
## refines (each functional averages over many cells), and whose covariance
## carries the multinomial cell-correlations the diagonal metric drops.
##
## Numerical detail that MATTERS: the polynomial basis must be built in a
## BOUNDED coordinate. Standardised-wealth powers z^p (z = (a - mean)/sd) are
## catastrophically ill-conditioned because z ranges to ~ +-20 on a curved grid
## with a long thin tail, so z^11 overflows the useful dynamic range and the
## Gram condition number reaches ~1e8, which DESTROYS the spread signal. Mapping
## wealth through its own D0-weighted CDF to u = 2 F(a) - 1 in [-1, 1] first (a
## classical bounded-domain orthogonalisation, cf. Gautschi's Stieltjes
## procedure) yields a basis orthonormal under D0 to machine precision, and a
## conditional-spread advantage over aggregates that is ~3 orders of magnitude
## (vs the understated ~80x from the ill-conditioned z-power basis).
## --------------------------------------------------------------------------


#' Contract a mixture distribution Jacobian against a shock into a net reweighting
#'
#' Convenience contraction shared by the reweighting losses: given a mixture
#' distribution Jacobian \code{JD_mix} (see
#' \code{\link{hank_mixture_dist_jacobian}}) and a shock path \code{dZ},
#' returns the model-implied net distribution change SUMMED over the horizon,
#' \deqn{dD[cell] = \Sigma_i \Sigma_t \Sigma_s JD_{mix,i}[t, s, cell]\, dZ_i[s],}
#' i.e. the total accumulated cell-mass displacement over all \code{T_h} dates
#' (the single-vector convention of \code{\link{hank_reweighting_loss}}).
#'
#' @param JD_mix Named list \code{JD_mix[[input]]}, each a
#'   \code{T_h x T_h x n_cell} array (as returned by
#'   \code{\link{hank_mixture_dist_jacobian}}).
#' @param dZ Named list, one length-\code{T_h} shock path per entry of
#'   \code{names(JD_mix)} (missing entries treated as all-zero).
#'
#' @return A length-\code{n_cell} numeric vector, the accumulated net
#'   reweighting.
#' @export
hank_dist_response <- function(JD_mix, dZ) {
  inputs <- names(JD_mix)
  T_h    <- dim(JD_mix[[1L]])[1L]
  n_cell <- dim(JD_mix[[1L]])[3L]
  dD <- numeric(n_cell)
  for (i in inputs) {
    dz <- dZ[[i]]
    if (is.null(dz)) next
    JDi <- JD_mix[[i]]
    for (tt in seq_len(T_h))
      dD <- dD + as.numeric(crossprod(matrix(JDi[tt, , ], nrow = T_h, ncol = n_cell), dz))
  }
  dD
}


#' Date-t* snapshot of a mixture distribution Jacobian contracted against a shock
#'
#' Like \code{\link{hank_dist_response}} but returns the distribution response at
#' a SINGLE output date \code{t_star} (a real survey compares two dated cross-
#' sections), \deqn{dD[cell] = \Sigma_i \Sigma_s JD_{mix,i}[t\_star, s, cell]\,
#' dZ_i[s],} i.e. the \code{t_star} output-date ROW of the Jacobian contracted
#' against the shock path, as opposed to \code{\link{hank_dist_response}}'s SUM
#' over all output dates. This is exactly column \code{t_star} of the per-date
#' model reweighting that \code{\link{hank_reweighting_loss}} forms internally.
#'
#' @param JD_mix Named list \code{JD_mix[[input]]}, each a \code{T_h x T_h x
#'   n_cell} array (from \code{\link{hank_mixture_dist_jacobian}}).
#' @param dZ Named list, one length-\code{T_h} shock path per entry of
#'   \code{names(JD_mix)} (missing entries treated as all-zero).
#' @param t_star Integer output date in \code{1..T_h} at which to snapshot.
#' @return Length-\code{n_cell} numeric vector: the date-\code{t_star} net
#'   reweighting.
#' @export
hank_dist_response_snapshot <- function(JD_mix, dZ, t_star) {
  inputs <- names(JD_mix)
  T_h    <- dim(JD_mix[[1L]])[1L]
  n_cell <- dim(JD_mix[[1L]])[3L]
  if (length(t_star) != 1L || t_star < 1L || t_star > T_h)
    stop(sprintf(
      "hank_dist_response_snapshot(): t_star must be a single integer in 1..%d.",
      T_h))
  dD <- numeric(n_cell)
  for (i in inputs) {
    dz <- dZ[[i]]
    if (is.null(dz)) next
    JDi <- JD_mix[[i]]
    dD <- dD + as.numeric(crossprod(matrix(JDi[t_star, , ], nrow = T_h,
                                           ncol = n_cell), dz))
  }
  dD
}


#' Grid-invariant orthonormal wealth-functional basis for the reweighting loss
#'
#' Builds a moderate basis of smooth distribution functionals -- low-order
#' monomials in a BOUNDED wealth coordinate, orthonormalised under the
#' reference distribution \code{D0} -- for projecting a net cross-section
#' reweighting onto a grid-invariant, well-conditioned subspace (see the file
#' section comment above for why the naive cell-diagonal metric and a raw
#' standardised-wealth power basis both fail).
#'
#' The bounded coordinate is the \code{D0}-weighted (midpoint) empirical CDF of
#' wealth mapped to \code{[-1, 1]}: with cells sorted by \code{a_cell},
#' \code{Fc = cumsum(D0) - 0.5 D0} and \code{u = 2 Fc - 1}. Monomials
#' \code{u^1, ..., u^degree} are then orthonormalised under the \code{D0} inner
#' product \code{<f, g> = Sum_c D0_c f_c g_c} (via a QR of
#' \code{sqrt(D0) [u^0..u^degree]}, dropping the constant column), so the
#' returned columns satisfy \code{t(B) diag(D0) B = I} to machine precision.
#' Cells carrying no reference mass (\code{D0 == 0}) contribute nothing to the
#' inner product and are assigned basis value 0.
#'
#' @param a_cell Length-\code{n_cell} numeric: the wealth (ordinating)
#'   coordinate of each distribution cell, in the SAME cell order as \code{D0}
#'   and the distribution Jacobian (asset-fast; see
#'   \code{\link{hank_mixture_dist}}). Income enters only through the tie
#'   ordering of equal-wealth cells -- the basis is a deliberate low-
#'   dimensional projection onto smooth functions of WEALTH, where the
#'   discount-factor-spread signal lives.
#' @param D0 Length-\code{n_cell} reference (e.g. pre-shock stationary)
#'   distribution, non-negative, summing to 1.
#' @param degree Integer polynomial degree = number of returned basis columns;
#'   default 11.
#' @param tail_trim Single number in \code{[0, 1)}, the cumulative-\code{D0}-
#'   mass fraction trimmed from EACH tail before building the basis; the
#'   orthonormal basis is then built using only the "core" cells whose
#'   midpoint CDF lies in \code{[tail_trim, 1 - tail_trim]}, with all other
#'   (extreme-tail) cells forced to basis value 0. Default 0 reproduces the
#'   untrimmed basis exactly (the tail amplification described above is left
#'   in place).
#'
#' @return An \code{n_cell x degree} matrix \code{B}, orthonormal under
#'   \code{D0} (\code{t(B) diag(D0) B = I}).
#' @export
hank_reweight_functional_basis <- function(a_cell, D0, degree = 11L, tail_trim = 0) {
  n_cell <- length(D0)
  if (length(a_cell) != n_cell)
    stop(sprintf("hank_reweight_functional_basis(): length(a_cell) (%d) must equal length(D0) (%d).",
                 length(a_cell), n_cell))
  degree <- as.integer(degree)
  if (is.na(degree) || degree < 1L)
    stop("hank_reweight_functional_basis(): 'degree' must be an integer >= 1.")
  if (!is.numeric(tail_trim) || length(tail_trim) != 1L || is.na(tail_trim) ||
      tail_trim < 0 || tail_trim >= 0.5)
    stop("hank_reweight_functional_basis(): 'tail_trim' must be a single number in [0, 0.5).")

  ord    <- order(a_cell)
  Fc     <- cumsum(D0[ord]) - 0.5 * D0[ord]    # D0-weighted midpoint CDF in [0, 1]
  u      <- numeric(n_cell)
  u[ord] <- 2 * Fc - 1                          # bounded coordinate in [-1, 1]

  if (tail_trim > 0) {
    core        <- logical(n_cell)
    core[ord]   <- Fc >= tail_trim & Fc <= (1 - tail_trim)   # CORE cells only
  } else {
    core <- rep(TRUE, n_cell)
  }

  B <- matrix(0, n_cell, degree)
  raw_core <- vapply(0:degree, function(p) u[core]^p, numeric(sum(core)))  # core x (degree+1), incl constant
  sq_core  <- sqrt(D0[core])
  Q_core   <- qr.Q(qr(sq_core * raw_core))
  B_core   <- Q_core / sq_core
  B_core[!is.finite(B_core)] <- 0               # D0 == 0 cells: no mass, zero basis value
  B[core, ] <- B_core[, -1L, drop = FALSE]      # drop the constant column; trimmed cells stay 0
  B
}


#' Proper K x K covariance of the projected net reweighting
#' @keywords internal
.hank_reweight_functional_cov <- function(B, D0, D1, N) {
  cvD <- function(D) {
    mu <- as.numeric(crossprod(B, D))            # t(B) D  (length K)
    crossprod(B, D * B) - outer(mu, mu)           # t(B) diag(D) B - mu mu^T
  }
  (cvD(D0) + cvD(D1)) / N
}


#' Grid-invariant functional reweighting metric (basis + covariance)
#'
#' Bundles the grid-invariant wealth-functional basis
#' (\code{\link{hank_reweight_functional_basis}}) with the PROPER \code{K x K}
#' multinomial-difference covariance of the projected net reweighting, for use
#' by \code{\link{hank_reweight_functional_loss}}. The covariance is that of
#' \code{t(B) (Dhat1 - Dhat0)} where \code{Dhat0}, \code{Dhat1} are independent
#' size-\code{N} multinomial samples from the reference and post-shock truth
#' distributions:
#' \deqn{\Sigma = (cov_{D0}(B) + cov_{D1}(B)) / N, \quad cov_D(B) = t(B)
#' diag(D) B - (t(B) D)(t(B) D)^T.}
#' Building the metric ONCE from a fixed reference \code{(D0, D1)} and reusing
#' it across candidate mixtures (each supplying only its own model-implied
#' \code{dD_mod}) is the intended grid-search usage.
#'
#' @param a_cell,D0,degree,tail_trim As in
#'   \code{\link{hank_reweight_functional_basis}}.
#' @param D1 Length-\code{n_cell} post-shock reference distribution (the second
#'   independent cross-section), non-negative, summing to 1.
#' @param N Sample size of EACH of the two cross-sections.
#'
#' @return A list with:
#'   \item{B}{The \code{n_cell x degree} \code{D0}-orthonormal basis.}
#'   \item{Sigma}{The \code{degree x degree} projected covariance.}
#'   \item{Sigma_inv}{Its inverse (the GLS precision used by the loss).}
#'   \item{orthonormality_err}{\code{max|t(B) diag(D0) B - I|} (a machine-
#'     precision diagnostic).}
#'   \item{cond}{Condition number of \code{Sigma} (basis-degree-dependent;
#'     reported so callers can reduce \code{degree} if it is too large).}
#' @export
hank_reweight_functional_metric <- function(a_cell, D0, D1, N, degree = 11L, tail_trim = 0) {
  if (length(D0) != length(D1))
    stop(sprintf("hank_reweight_functional_metric(): length(D0) (%d) must equal length(D1) (%d).",
                 length(D0), length(D1)))
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0)
    stop("hank_reweight_functional_metric(): 'N' must be a single finite positive number.")
  B   <- hank_reweight_functional_basis(a_cell, D0, degree = degree, tail_trim = tail_trim)
  Sig <- .hank_reweight_functional_cov(B, D0, D1, N)
  list(B = B, Sigma = Sig, Sigma_inv = solve(Sig),
       orthonormality_err = max(abs(crossprod(B, D0 * B) - diag(ncol(B)))),
       cond = kappa(Sig))
}


#' Grid-invariant functional LEVEL metric (single-survey basis + covariance)
#'
#' The single-survey LEVEL analogue of
#' \code{\link{hank_reweight_functional_metric}}: bundles the grid-invariant
#' wealth-functional basis (\code{\link{hank_reweight_functional_basis}}) with
#' the PROPER \code{K x K} multinomial covariance of a SINGLE projected
#' stationary cross-section, for use by \code{\link{hank_reweight_functional_loss}}.
#' Where the RESPONSE metric's covariance is that of a DIFFERENCE of two
#' independent cross-sections (\code{(cov_D0(B) + cov_D1(B)) / N}), the LEVEL
#' metric uses only ONE survey drawn from the stationary distribution
#' \code{D0} itself, so its covariance carries a SINGLE \code{cov_D0} term:
#' \deqn{\Sigma = cov_{D0}(B) / N, \quad cov_D(B) = t(B) diag(D) B - (t(B)
#' D)(t(B) D)^T.}
#' This is the covariance of \code{t(B) Dhat} for \code{Dhat} a single size-
#' \code{N} multinomial sample from \code{D0} (e.g. a stationary household
#' wealth survey), as opposed to \code{t(B) (Dhat1 - Dhat0)} of two dated
#' cross-sections straddling a shock.
#'
#' @param a_cell,D0,degree,tail_trim As in
#'   \code{\link{hank_reweight_functional_basis}}: \code{D0} is the reference
#'   (stationary) distribution the single survey is drawn from.
#' @param N Sample size of the single cross-section.
#'
#' @return A list with:
#'   \item{B}{The \code{n_cell x degree} \code{D0}-orthonormal basis.}
#'   \item{Sigma}{The \code{degree x degree} projected covariance,
#'     \code{cov_D0(B) / N} (a SINGLE \code{cov_D0} term; contrast
#'     \code{\link{hank_reweight_functional_metric}}'s two-survey sum).}
#'   \item{Sigma_inv}{Its inverse (the GLS precision used by the loss).}
#'   \item{orthonormality_err}{\code{max|t(B) diag(D0) B - I|} (a machine-
#'     precision diagnostic).}
#'   \item{cond}{Condition number of \code{Sigma} (basis-degree-dependent;
#'     reported so callers can reduce \code{degree} if it is too large).}
#' @seealso \code{\link{hank_reweight_functional_metric}},
#'   \code{\link{hank_reweight_functional_loss}}
#' @export
hank_reweight_level_metric <- function(a_cell, D0, N, degree = 11L, tail_trim = 0) {
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0)
    stop("hank_reweight_level_metric(): 'N' must be a single finite positive number.")
  B   <- hank_reweight_functional_basis(a_cell, D0, degree = degree, tail_trim = tail_trim)
  Sig <- .hank_reweight_functional_cov(B, D0, D0, 2 * N) # (cov_D0+cov_D0)/(2N) = cov_D0/N
  list(B = B, Sigma = Sig, Sigma_inv = solve(Sig),
       orthonormality_err = max(abs(crossprod(B, D0 * B) - diag(ncol(B)))),
       cond = kappa(Sig))
}


#' Grid-invariant functional reweighting loss
#'
#' The proper Mahalanobis (generalised-least-squares) discrepancy between an
#' OBSERVED net reweighting \code{dD_hat} and a candidate model-implied net
#' reweighting \code{dD_mod}, projected onto the grid-invariant wealth-
#' functional basis of \code{metric} and weighted by its \code{K x K}
#' covariance:
#' \deqn{\rho = t(d) \Sigma^{-1} d, \quad d = t(B) (dD\_hat - dD\_mod).}
#' Under the truth (\code{dD_mod} at the true mixture, \code{dD_hat} a genuine
#' pair of size-\code{N} cross-sections) \code{rho} is asymptotically
#' chi-square on \code{ncol(B)} degrees of freedom, so it is directly
#' comparable across grids -- unlike the cell-diagonal
#' \code{\link{hank_reweight_metric}} loss (see file section comment).
#'
#' @param dD_hat Observed net reweighting, a length-\code{n_cell} vector (e.g.
#'   \code{Dhat1 - Dhat0} of two calibrated survey cross-sections, or
#'   \code{\link{hank_dist_response}} of a truth Jacobian for a noiseless run).
#' @param dD_mod Candidate model-implied net reweighting, a length-\code{n_cell}
#'   vector (e.g. \code{\link{hank_dist_response}} of the candidate mixture
#'   Jacobian).
#' @param metric A metric object from
#'   \code{\link{hank_reweight_functional_metric}} (its \code{B} and
#'   \code{Sigma_inv} are used).
#'
#' @return A list with:
#'   \item{rho}{The scalar functional loss.}
#'   \item{d}{The length-\code{ncol(B)} projected discrepancy \code{t(B)
#'     (dD_hat - dD_mod)}.}
#'   \item{resid}{\code{dD_hat - dD_mod} on the cell basis.}
#' @export
hank_reweight_functional_loss <- function(dD_hat, dD_mod, metric) {
  if (length(dD_hat) != length(dD_mod))
    stop(sprintf("hank_reweight_functional_loss(): length(dD_hat) (%d) must equal length(dD_mod) (%d).",
                 length(dD_hat), length(dD_mod)))
  resid <- dD_hat - dD_mod
  d     <- as.numeric(crossprod(metric$B, resid))       # t(B) (dD_hat - dD_mod)
  rho   <- as.numeric(crossprod(d, metric$Sigma_inv %*% d))
  list(rho = rho, d = d, resid = resid)
}
