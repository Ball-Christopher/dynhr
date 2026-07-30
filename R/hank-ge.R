## R/hank-ge.R
## --------------------------------------------------------------------------
## General-equilibrium sequence-space solve for a heterogeneous-agent model,
## demonstrated on the canonical Krusell-Smith economy (het household + a
## representative firm + asset-market clearing).  Composes the het-block
## sequence-space Jacobian (R/hank-jacobian.R) with the firm's simple-block
## Jacobians into the GE system and solves
##
##   dU = - H_U^{-1} H_Z dZ                         (ABRS 2021, eq. 30)
##
## for the unknown aggregate path U (capital K) given an exogenous shock path Z
## (TFP), and the fully nonlinear transition by quasi-Newton with the frozen
## steady-state Jacobian (ABRS eq. 38).
##
## Krusell-Smith closure (Z_t aggregate TFP, L = E[e] = 1, capital predetermined):
##   firm:   r_t = alpha Z_t K_{t-1}^{alpha-1} - delta
##           w_t = (1-alpha) Z_t K_{t-1}^{alpha}
##   market: A_t(r, w) = K_t         (household assets = capital)
## Unknown U = {K_t}, target = asset-market clearing, exogenous Z = {Z_t}.
##
## The GE-solve/assembly here is written for this closure; generalising to an
## arbitrary block DAG (per the sequence-jacobian `create_model` API) is a
## follow-up increment.
## --------------------------------------------------------------------------


#' Krusell-Smith general-equilibrium steady state
#'
#' Solves for the equilibrium interest rate \eqn{r^*} at which household
#' aggregate assets equal the firm's capital demand, then returns the full
#' steady-state het block and firm quantities.  The market-clearing solve
#' doubles as an Aiyagari-style external check on the household block.
#'
#' @param a_grid,Pi,e Household grid, income transition, income levels.
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param alpha,delta Capital share and depreciation.
#' @param Z Steady-state TFP (default 1).
#' @param r_bracket Optional length-2 search bracket for \eqn{r^*}.
#'
#' @return A list of class \code{hank_ks} with \code{r}, \code{w}, \code{K},
#'   \code{Z}, the calibration, and the steady-state \code{block}
#'   (\code{\link{hank_het_block}}).
#' @export
hank_ks_steady <- function(a_grid, Pi, e, beta, eis, alpha, delta, Z = 1,
                           r_bracket = NULL) {
  K_of_r <- function(r) ((r + delta) / (alpha * Z))^(1 / (alpha - 1))
  w_of_r <- function(r) (1 - alpha) * Z * K_of_r(r)^alpha
  A_of_r <- function(r) {
    blk <- hank_het_block(a_grid, Pi, e, beta = beta, eis = eis,
                          r = r, w = w_of_r(r))
    blk$A
  }
  if (is.null(r_bracket))
    r_bracket <- c(-delta + 1e-4, 1 / beta - 1 - 1e-4)
  f <- function(r) A_of_r(r) - K_of_r(r)
  sol <- stats::uniroot(f, interval = r_bracket, tol = 1e-10)
  r <- sol$root; w <- w_of_r(r); K <- K_of_r(r)
  blk <- hank_het_block(a_grid, Pi, e, beta = beta, eis = eis, r = r, w = w)
  structure(list(r = r, w = w, K = K, Z = Z, alpha = alpha, delta = delta,
                 beta = beta, eis = eis, block = blk,
                 mkt_residual = blk$A - K),
            class = "hank_ks")
}


#' Krusell-Smith general-equilibrium steady state with a MIXTURE household
#'
#' Generalizes \code{\link{hank_ks_steady}} to a household that is a
#' \code{K}-type discount-factor MIXTURE (see \code{\link{hank_mixture_blocks}}):
#' every type shares the asset grid, income process and EIS, and interacts
#' with the rest of the economy only through the common prices \code{(r, w)},
#' so aggregate household assets at a candidate \code{r} are the EXACT
#' omega-weighted sum of the per-type steady-state assets,
#' \code{A_of_r(r) = Sum_k omega_k * hank_het_block(..., beta = betas[k], r =
#' r, w = w_of_r(r))$A}. The market-clearing solve is otherwise identical to
#' \code{hank_ks_steady}: \code{uniroot} on \code{A_of_r(r) - K_of_r(r)}.
#'
#' The upper end of the default search bracket, \code{1 / max(betas) - 1},
#' is set by the MOST-PATIENT type (the type with the highest discount
#' factor is the one whose Euler equation binds the equilibrium \eqn{r} from
#' above; any less-patient type's bound would be looser and could leave the
#' true root outside a bracket sized off it).
#'
#' @param a_grid,Pi,e Household grid, income transition, income levels
#'   (shared by every type).
#' @param betas Numeric length-K vector of discount factors, one per type.
#' @param omega Numeric length-K mixture weights, non-negative, summing to 1.
#' @param eis Elasticity of intertemporal substitution, shared (scalar).
#' @param alpha,delta Capital share and depreciation.
#' @param Z Steady-state TFP (default 1).
#' @param r_bracket Optional length-2 search bracket for \eqn{r^*}; defaults
#'   to \code{c(-delta + 1e-4, 1/max(betas) - 1 - 1e-4)} (see Details).
#'
#' @return A list of class \code{hank_mixture_ks} with \code{r}, \code{w},
#'   \code{K}, \code{Z}, the calibration (\code{betas}, \code{omega},
#'   \code{eis}, \code{alpha}, \code{delta}), the per-type steady-state
#'   \code{blocks} (list of \code{K} \code{\link{hank_het_block}} objects, all
#'   solved at the clearing \code{(r, w)}), and \code{mkt_residual}
#'   (\code{Sum_k omega_k*blocks[[k]]$A - K}, the mixture analog of
#'   \code{hank_ks}'s \code{mkt_residual}).
#' @export
hank_mixture_ks_steady <- function(a_grid, Pi, e, betas, omega, eis = 1,
                                   alpha, delta, Z = 1, r_bracket = NULL) {
  if (length(omega) != length(betas))
    stop(sprintf(
      "hank_mixture_ks_steady(): length(omega) (%d) must equal length(betas) (%d).",
      length(omega), length(betas)))
  if (any(!is.finite(omega)) || any(omega < 0))
    stop("hank_mixture_ks_steady(): 'omega' must be finite and non-negative.")
  if (abs(sum(omega) - 1) > 1e-8)
    stop(sprintf(
      "hank_mixture_ks_steady(): 'omega' must sum to 1 (sum = %.6f).", sum(omega)))
  K_of_r <- function(r) ((r + delta) / (alpha * Z))^(1 / (alpha - 1))
  w_of_r <- function(r) (1 - alpha) * Z * K_of_r(r)^alpha
  A_of_r <- function(r) {
    w <- w_of_r(r)
    acc <- 0
    for (k in seq_along(betas)) {
      blk <- hank_het_block(a_grid, Pi, e, beta = betas[k], eis = eis,
                            r = r, w = w)
      acc <- acc + omega[k] * blk$A
    }
    acc
  }
  ceil_r <- 1 / max(betas) - 1
  if (is.null(r_bracket))
    r_bracket <- c(-delta + 1e-4, ceil_r - 1e-4)
  f <- function(r) A_of_r(r) - K_of_r(r)

  ## Adaptive widen-toward-the-ceiling retry: when the most-patient type has
  ## small mixture weight, the pooled A_of_r(r) - K_of_r(r) can still be
  ## negative at the default upper bracket on a FINITE asset grid (the
  ## patient type's true divergence as r -> ceil_r is only resolved
  ## arbitrarily close to the ceiling, and/or needs a wider a_grid to show up
  ## before hitting the grid's own truncation) even though the bound is
  ## analytically correct in the continuum. Retry with the upper endpoint
  ## pushed successively closer to ceil_r before giving up -- this is a
  ## strictly WIDER (not different) bracket, so if it does find a sign
  ## change the root is still the genuine market-clearing rate.
  used_default_bracket <- identical(r_bracket, c(-delta + 1e-4, ceil_r - 1e-4))
  sol <- tryCatch(
    stats::uniroot(f, interval = r_bracket, tol = 1e-10),
    error = function(e) {
      ## Only retry-widen when we get here via the DEFAULT bracket (a
      ## caller-supplied r_bracket is respected as-is and fails loudly).
      if (used_default_bracket) {
        for (gap in c(1e-6, 1e-8, 1e-10, 1e-12)) {
          rb2 <- c(r_bracket[1L], ceil_r - gap)
          retry <- tryCatch(stats::uniroot(f, interval = rb2, tol = 1e-10),
                             error = function(e2) NULL)
          if (!is.null(retry)) return(retry)
        }
      }
      stop(sprintf(paste0(
        "hank_mixture_ks_steady(): uniroot() failed to bracket a root of ",
        "the mixture asset-market clearing condition on r_bracket = ",
        "[%.6f, %.6f] (f(lo) = %.6g, f(hi) = %.6g), even after retrying with ",
        "the upper endpoint pushed toward the most-patient type's ceiling ",
        "1/max(betas)-1 = %.6f. This usually means the asset grid's 'amax' ",
        "is too small to resolve the most-patient type's divergence near ",
        "its own Euler bound at this mixture weight -- try a wider a_grid, ",
        "or pass r_bracket explicitly. Original error: %s"),
        r_bracket[1L], r_bracket[2L],
        tryCatch(f(r_bracket[1L]), error = function(e2) NA_real_),
        tryCatch(f(r_bracket[2L]), error = function(e2) NA_real_),
        ceil_r, conditionMessage(e)))
    })
  r <- sol$root; w <- w_of_r(r); K <- K_of_r(r)
  blocks <- lapply(seq_along(betas), function(k)
    hank_het_block(a_grid, Pi, e, beta = betas[k], eis = eis, r = r, w = w))
  A_mix <- sum(vapply(seq_along(blocks), function(k) omega[k] * blocks[[k]]$A,
                       numeric(1)))
  structure(list(r = r, w = w, K = K, Z = Z, alpha = alpha, delta = delta,
                 betas = betas, omega = omega, eis = eis, blocks = blocks,
                 mkt_residual = A_mix - K),
            class = "hank_mixture_ks")
}


#' Krusell-Smith general-equilibrium steady state with a MIXTURE household
#' differing in its INCOME PROCESS and/or BORROWING CONSTRAINT
#' (income-risk and wealth heterogeneity axes)
#'
#' Sibling of \code{\link{hank_mixture_ks_steady}} for a \code{K}-type mixture
#' whose types differ in their income process (\code{Pi}, \code{e} -- e.g. a
#' per-type Rouwenhorst calibration), their borrowing constraint
#' (\code{amin}, the WEALTH heterogeneity axis: per-type \code{amin} on the
#' shared \code{a_grid}), and/or their EIS -- not (only) in \code{beta}.
#' Every type still shares the asset grid and prices, and interacts with the
#' rest of the economy only through the common \code{(r, w)}, so the
#' household side of market clearing is again the EXACT omega-weighted sum
#' of the per-type steady-state assets,
#' \code{A_of_r(r) = Sum_k omega_k * hank_het_block(a_grid, Pi_k, e_k, beta =
#' beta_k, eis = eis_k, amin = amin_k, r = r, w = w_of_r(r))$A}; the firm
#' FOCs and the \code{uniroot} market-clearing solve are otherwise IDENTICAL
#' to \code{\link{hank_mixture_ks_steady}} -- only the per-type block
#' construction changes. This is the market-clearing sibling of the
#' fixed-price \code{\link{hank_mixture_ks_assemble}} (same \code{types}
#' schema).
#'
#' @param a_grid Numeric asset grid (see \code{\link{hank_asset_grid}}),
#'   shared by every type.
#' @param types List of length K, one entry per type, each a list with
#'   elements \code{beta}, \code{Pi} (\code{n_e_k x n_e_k} income transition
#'   matrix), and \code{e} (length-\code{n_e_k} income levels), plus
#'   optionally \code{eis} (falls back to the \code{eis} argument) and
#'   \code{amin} (per-type borrowing constraint, falls back to
#'   \code{a_grid[1]}; must satisfy \code{amin >= a_grid[1]} -- see
#'   \code{\link{hank_het_block}}). \code{Pi}/\code{e} may differ in size and
#'   value across types (that is the point of this function); only
#'   \code{a_grid} is shared. Because a per-type \code{amin} lives on the
#'   SHARED grid, a pure wealth-axis mixture (identical \code{Pi}/\code{e},
#'   distinct \code{amin}) keeps a common cell space, so pooled
#'   distribution objects downstream remain valid.
#' @param omega Numeric length-K mixture weights, non-negative, summing to 1.
#' @param eis Default elasticity of intertemporal substitution for types
#'   that don't set their own.
#' @param alpha,delta Capital share and depreciation.
#' @param Z Steady-state TFP (default 1).
#' @param r_bracket Optional length-2 search bracket for \eqn{r^*}; defaults
#'   to \code{c(-delta + 1e-4, 1/max(betas) - 1 - 1e-4)}, exactly as in
#'   \code{\link{hank_mixture_ks_steady}}, with \code{betas} read off
#'   \code{types}.
#'
#' @return A list of class \code{hank_mixture_ks} (SAME shape as
#'   \code{\link{hank_mixture_ks_steady}}'s return, so it is a drop-in
#'   argument to \code{\link{hank_mixture_ks_model}} and every other
#'   \code{hank_mixture_ks} consumer) with \code{r}, \code{w}, \code{K},
#'   \code{Z}, the calibration (\code{betas}, \code{omega}, \code{eis},
#'   \code{alpha}, \code{delta}), the per-type steady-state \code{blocks}
#'   (list of \code{K} \code{\link{hank_het_block}} objects, each built from
#'   its OWN \code{Pi}/\code{e}, all solved at the clearing \code{(r, w)}),
#'   and \code{mkt_residual}.
#' @export
hank_mixture_ks_steady_hetinc <- function(a_grid, types, omega, eis = 1,
                                          alpha, delta, Z = 1,
                                          r_bracket = NULL) {
  if (!is.list(types) || length(types) < 1L)
    stop("hank_mixture_ks_steady_hetinc(): 'types' must be a non-empty list.")
  if (!all(vapply(types, function(tt)
    is.list(tt) && all(c("beta", "Pi", "e") %in% names(tt)), logical(1))))
    stop("hank_mixture_ks_steady_hetinc(): every entry of 'types' must be a ",
         "list with elements 'beta', 'Pi', 'e' (optionally 'eis', 'amin').")
  betas <- vapply(types, function(tt) tt$beta, numeric(1))
  if (length(omega) != length(types))
    stop(sprintf(
      "hank_mixture_ks_steady_hetinc(): length(omega) (%d) must equal length(types) (%d).",
      length(omega), length(types)))
  if (any(!is.finite(omega)) || any(omega < 0))
    stop("hank_mixture_ks_steady_hetinc(): 'omega' must be finite and non-negative.")
  if (abs(sum(omega) - 1) > 1e-8)
    stop(sprintf(
      "hank_mixture_ks_steady_hetinc(): 'omega' must sum to 1 (sum = %.6f).",
      sum(omega)))

  K_of_r <- function(r) ((r + delta) / (alpha * Z))^(1 / (alpha - 1))
  w_of_r <- function(r) (1 - alpha) * Z * K_of_r(r)^alpha
  type_block <- function(tt, r, w)
    hank_het_block(a_grid, tt$Pi, tt$e, beta = tt$beta,
                   eis = if (!is.null(tt$eis)) tt$eis else eis,
                   r = r, w = w,
                   amin = if (!is.null(tt$amin)) tt$amin else a_grid[1L])
  A_of_r <- function(r) {
    w <- w_of_r(r)
    acc <- 0
    for (k in seq_along(types))
      acc <- acc + omega[k] * type_block(types[[k]], r, w)$A
    acc
  }
  ceil_r <- 1 / max(betas) - 1
  if (is.null(r_bracket))
    r_bracket <- c(-delta + 1e-4, ceil_r - 1e-4)
  f <- function(r) A_of_r(r) - K_of_r(r)

  ## Identical adaptive widen-toward-the-ceiling retry as
  ## hank_mixture_ks_steady() -- see that function's comment for the
  ## rationale (finite-grid truncation can hide the most-patient type's
  ## divergence near its own Euler bound).
  used_default_bracket <- identical(r_bracket, c(-delta + 1e-4, ceil_r - 1e-4))
  sol <- tryCatch(
    stats::uniroot(f, interval = r_bracket, tol = 1e-10),
    error = function(e) {
      if (used_default_bracket) {
        for (gap in c(1e-6, 1e-8, 1e-10, 1e-12)) {
          rb2 <- c(r_bracket[1L], ceil_r - gap)
          retry <- tryCatch(stats::uniroot(f, interval = rb2, tol = 1e-10),
                             error = function(e2) NULL)
          if (!is.null(retry)) return(retry)
        }
      }
      stop(sprintf(paste0(
        "hank_mixture_ks_steady_hetinc(): uniroot() failed to bracket a root ",
        "of the mixture asset-market clearing condition on r_bracket = ",
        "[%.6f, %.6f] (f(lo) = %.6g, f(hi) = %.6g), even after retrying with ",
        "the upper endpoint pushed toward the most-patient type's ceiling ",
        "1/max(betas)-1 = %.6f. This usually means the asset grid's 'amax' ",
        "is too small to resolve the most-patient type's divergence near ",
        "its own Euler bound at this mixture weight -- try a wider a_grid, ",
        "or pass r_bracket explicitly. Original error: %s"),
        r_bracket[1L], r_bracket[2L],
        tryCatch(f(r_bracket[1L]), error = function(e2) NA_real_),
        tryCatch(f(r_bracket[2L]), error = function(e2) NA_real_),
        ceil_r, conditionMessage(e)))
    })
  r <- sol$root; w <- w_of_r(r); K <- K_of_r(r)
  blocks <- lapply(types, type_block, r = r, w = w)
  A_mix <- sum(vapply(seq_along(blocks), function(k) omega[k] * blocks[[k]]$A,
                       numeric(1)))
  structure(list(r = r, w = w, K = K, Z = Z, alpha = alpha, delta = delta,
                 betas = betas, omega = omega, eis = eis, blocks = blocks,
                 mkt_residual = A_mix - K),
            class = "hank_mixture_ks")
}


#' Assemble a fixed-price heterogeneous-mixture steady state (per-type amin +
#' income + EIS simultaneously)
#'
#' The fixed-price / K-matched analogue of \code{\link{hank_mixture_ks_steady}}
#' (which solves market clearing for \code{r}): builds a \code{hank_mixture_ks}
#' object from per-type household blocks at GIVEN prices \code{(r, w)}. Every
#' type may carry its OWN income process (\code{Pi}, \code{e}), borrowing limit
#' (\code{amin}), discount factor (\code{beta}) and EIS (\code{eis})
#' SIMULTANEOUSLY -- the combination the composition/interaction factorials need
#' and previously had to hand-splice onto a template mixture
#' (\code{mks$blocks <- ...}). No shipped constructor supported per-type
#' \code{amin} AND per-type income at once; this does.
#'
#' @param a_grid Shared asset grid (increasing).
#' @param types List of per-type specs; each a list with \code{beta}, \code{Pi},
#'   \code{e}, and optionally \code{eis} (falls back to the \code{eis} argument)
#'   and \code{amin} (falls back to \code{a_grid[1]}).
#' @param omega Type weights, length \code{length(types)}, non-negative, sum 1.
#' @param r,w Fixed prices at which every block is solved.
#' @param eis Default EIS for types that don't set their own.
#' @param alpha,delta,Z Firm/aggregate parameters carried on the object for
#'   \code{\link{hank_mixture_ks_model}}.
#' @param K Aggregate capital recorded on the object; \code{NULL} (default) sets
#'   it to the omega-weighted household asset supply \code{sum(omega * A_k)} at
#'   \code{(r, w)} (i.e. self-consistent), otherwise pins it to the supplied
#'   K-match target (the object then carries \code{mkt_residual = sum(omega*A) -
#'   K} so the caller can check the match).
#' @return A \code{hank_mixture_ks} object (same shape as
#'   \code{\link{hank_mixture_ks_steady}}), consumable by
#'   \code{\link{hank_mixture_ks_model}} / \code{\link{hank_model_irf}}.
#' @seealso \code{\link{hank_mixture_ks_steady}},
#'   \code{\link{hank_mixture_ks_steady_hetinc}}, \code{\link{hank_het_block}}
#' @export
hank_mixture_ks_assemble <- function(a_grid, types, omega, r, w,
                                     eis = 1, alpha, delta, Z = 1, K = NULL) {
  if (!is.list(types) || length(types) < 1L)
    stop("hank_mixture_ks_assemble(): 'types' must be a non-empty list.")
  if (!all(vapply(types, function(tt)
    is.list(tt) && all(c("beta", "Pi", "e") %in% names(tt)), logical(1))))
    stop("hank_mixture_ks_assemble(): every 'types' entry must be a list with ",
         "elements 'beta', 'Pi', 'e' (optionally 'eis', 'amin').")
  if (length(omega) != length(types))
    stop(sprintf(
      "hank_mixture_ks_assemble(): length(omega) (%d) must equal length(types) (%d).",
      length(omega), length(types)))
  if (any(!is.finite(omega)) || any(omega < 0) || abs(sum(omega) - 1) > 1e-8)
    stop("hank_mixture_ks_assemble(): 'omega' must be finite, non-negative, and sum to 1.")

  blocks <- lapply(types, function(tt)
    hank_het_block(a_grid, tt$Pi, tt$e, beta = tt$beta,
                   eis = if (!is.null(tt$eis)) tt$eis else eis,
                   r = r, w = w,
                   amin = if (!is.null(tt$amin)) tt$amin else a_grid[1L]))
  A_mix <- sum(vapply(seq_along(blocks),
                      function(k) omega[k] * blocks[[k]]$A, numeric(1)))
  if (is.null(K)) K <- A_mix
  betas <- vapply(types, function(tt) tt$beta, numeric(1))
  structure(list(r = r, w = w, K = K, Z = Z, alpha = alpha, delta = delta,
                 betas = betas, omega = omega, eis = eis, blocks = blocks,
                 a_grid = a_grid, mkt_residual = A_mix - K),
            class = "hank_mixture_ks")
}


#' GE sequence-space Jacobians for the Krusell-Smith model
#'
#' Assembles \code{H_K = dH/dK} and \code{H_Z = dH/dZ} of the asset-market
#' clearing target from the household Jacobians (fake-news) and the firm's
#' analytic simple-block Jacobians.
#'
#' @param ks A \code{\link{hank_ks_steady}} steady state.
#' @param T_h Integer horizon.
#'
#' @return A list with \code{H_K}, \code{H_Z} (\code{T x T}), the household
#'   Jacobians \code{J} (from \code{\link{hank_het_jacobian}}), and the firm
#'   Jacobian blocks \code{Jr_K, Jr_Z, Jw_K, Jw_Z}.
#' @export
hank_ks_ge_jacobian <- function(ks, T_h) {
  alpha <- ks$alpha; delta <- ks$delta; Z <- ks$Z; K <- ks$K
  ## Firm simple-block Jacobians (K_{t-1} predetermined -> lag = subdiagonal).
  dr_dKlag <- alpha * (alpha - 1) * Z * K^(alpha - 2)
  dr_dZ    <- alpha * K^(alpha - 1)
  dw_dKlag <- (1 - alpha) * alpha * Z * K^(alpha - 1)
  dw_dZ    <- (1 - alpha) * K^(alpha)
  I  <- diag(T_h)
  lag <- rbind(0, cbind(diag(T_h - 1L), 0))   # [t, t-1] = 1
  Jr_K <- dr_dKlag * lag
  Jr_Z <- dr_dZ * I
  Jw_K <- dw_dKlag * lag
  Jw_Z <- dw_dZ * I

  J <- hank_het_jacobian(ks$block, T_h, inputs = c("r", "w"),
                         outputs = c("A", "C"))
  J_Ar <- J[["A"]][["r"]]; J_Aw <- J[["A"]][["w"]]

  ## Chain rule through the DAG:  A depends on (r,w); (r,w) depend on (K,Z).
  H_K <- J_Ar %*% Jr_K + J_Aw %*% Jw_K - I
  H_Z <- J_Ar %*% Jr_Z + J_Aw %*% Jw_Z
  list(H_K = H_K, H_Z = H_Z, J = J,
       Jr_K = Jr_K, Jr_Z = Jr_Z, Jw_K = Jw_K, Jw_Z = Jw_Z,
       alpha = alpha, Z = Z, K = K)
}


#' Linear GE impulse response of the Krusell-Smith model
#'
#' Solves \eqn{dK = -H_K^{-1} H_Z\, dZ} and propagates to prices and aggregates.
#'
#' @param ks A \code{\link{hank_ks_steady}} steady state.
#' @param dZ Numeric length-\code{T} TFP shock path (deviations from \code{Z}).
#' @param ge Optional precomputed \code{\link{hank_ks_ge_jacobian}} result.
#'
#' @return A list of deviation paths \code{dK, dr, dw, dA} and the linear
#'   asset-market residual \code{mkt_resid} (should be ~0 by construction).
#' @export
hank_ks_linear_irf <- function(ks, dZ, ge = NULL) {
  T_h <- length(dZ)
  if (is.null(ge)) ge <- hank_ks_ge_jacobian(ks, T_h)
  dK <- as.numeric(-solve(ge$H_K, ge$H_Z %*% dZ))
  dr <- as.numeric(ge$Jr_K %*% dK + ge$Jr_Z %*% dZ)
  dw <- as.numeric(ge$Jw_K %*% dK + ge$Jw_Z %*% dZ)
  dA <- as.numeric(ge$J[["A"]][["r"]] %*% dr + ge$J[["A"]][["w"]] %*% dw)
  dC <- as.numeric(ge$J[["C"]][["r"]] %*% dr + ge$J[["C"]][["w"]] %*% dw)
  ## Output: Y_t = Z_t K_{t-1}^alpha (L=1). dK_lag has K_0 predetermined at ss.
  alpha <- ge$alpha; Zss <- ge$Z; Kss <- ge$K
  dK_lag <- c(0, dK[-length(dK)])
  dY <- alpha * Zss * Kss^(alpha - 1) * dK_lag + Kss^(alpha) * dZ
  list(dK = dK, dr = dr, dw = dw, dA = dA, dC = dC, dY = dY,
       mkt_resid = dA - dK)
}


#' Nonlinear GE transition of the Krusell-Smith model (quasi-Newton)
#'
#' Solves the fully nonlinear perfect-foresight transition for the capital path
#' given a TFP path, by Newton's method with the frozen steady-state Jacobian
#' \code{H_K} (ABRS eq. 38); each residual re-runs the nonlinear household
#' transition (\code{\link{hank_td_nonlinear}}).
#'
#' @param ks A \code{\link{hank_ks_steady}} steady state.
#' @param Z_path Numeric length-\code{T} TFP path (levels).
#' @param ge Optional precomputed \code{\link{hank_ks_ge_jacobian}} result.
#' @param tol,maxit Newton tolerance and iteration cap.
#'
#' @return A list with level paths \code{K, r, w, A}, plus \code{converged},
#'   \code{iterations}, and the final \code{max_resid}.
#' @export
hank_ks_nonlinear_irf <- function(ks, Z_path, ge = NULL, tol = 1e-9,
                                  maxit = 50L) {
  T_h <- length(Z_path)
  if (is.null(ge)) ge <- hank_ks_ge_jacobian(ks, T_h)
  alpha <- ks$alpha; delta <- ks$delta; Kss <- ks$K
  H_K_lu <- ge$H_K

  Kpath <- rep(Kss, T_h)
  converged <- FALSE; it <- 0L; max_resid <- Inf
  for (it in seq_len(maxit)) {
    Klag <- c(Kss, Kpath[-T_h])
    r_path <- alpha * Z_path * Klag^(alpha - 1) - delta
    w_path <- (1 - alpha) * Z_path * Klag^(alpha)
    td <- hank_td_nonlinear(ks$block, r_path = r_path, w_path = w_path,
                            T_h = T_h)
    resid <- td$A - Kpath
    max_resid <- max(abs(resid))
    if (max_resid < tol) { converged <- TRUE; break }
    Kpath <- Kpath - as.numeric(solve(H_K_lu, resid))
  }
  Klag <- c(Kss, Kpath[-T_h])
  r_path <- alpha * Z_path * Klag^(alpha - 1) - delta
  w_path <- (1 - alpha) * Z_path * Klag^(alpha)
  td <- hank_td_nonlinear(ks$block, r_path = r_path, w_path = w_path, T_h = T_h)
  list(K = Kpath, r = r_path, w = w_path, A = td$A,
       converged = converged, iterations = it, max_resid = max_resid)
}
