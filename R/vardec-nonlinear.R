## R/vardec-nonlinear.R
## --------------------------------------------------------------------------
## Variance decomposition for pruned order-2/3 solutions, with a Shapley
## (marginal-contribution) allocation of the cross-shock interaction.
##
## At order >= 2 the shocks interact (e1*e2 and e_i^2 terms, x1 (x) x2 state
## products), so the unconditional variance is NOT additive over shocks: the
## "one shock at a time" shares do not sum to one.  The coalition value
##
##   v(S) = Var(y) when only the (Cholesky-orthogonalised) shocks in S are
##          active,
##
## is computed exactly from the pruned augmented state space (the same
## .order2_aug_system / .order3_aug_system machinery as
## compute_moments_order2() / pruned_ss_moments3()), and the Shapley value
##
##   phi_i = sum_{S subset N\{i}} |S|! (n-|S|-1)! / n! * (v(S u {i}) - v(S))
##
## allocates v(N) - v(0) = Var(y) across shocks (efficiency), treats
## symmetric shocks symmetrically, and gives a shock that never matters zero.
##
## Relation to MacroModelling.jl (v0.2.0, get_variance_decomposition with
## marginal_contribution = true; docs/src/shapley_decompositions.md):
##   * same allocation rule (Shapley value of the coalition variance);
##   * same default "per-shock + Cross_shock_interaction" table for the
##     non-Shapley method (here method = "sequential");
##   * DIFFERENT coalition value: MacroModelling.jl projects the inner shock
##     cumulant block while keeping the STATE cumulants at their full-shock
##     values ("not a counterfactual recomputation"); here v(S) is the exact
##     variance of the pruned system re-solved with the inactive shocks'
##     innovations set to zero (every state moment recomputed), i.e. the
##     variance a pruned simulation with those shock columns zeroed would
##     produce.  v(empty) = 0 exactly.
##   * the decision rule itself (including the risk/sigma-correction terms
##     ghss, ghxss, ghuss of an order-3 solution) is held at the full-Sigma_e
##     solution: agents still price the full risk, only the REALISED
##     innovations of the inactive shocks are zero.  At order 2 the risk terms
##     only shift the mean, so this choice does not matter there.
##   * MacroModelling.jl evaluates the Shapley value by an Aumann-Shapley
##     path integral on the multilinear extension; here it is exact: full
##     enumeration (n_shocks <= max_exact), or the polynomial-coefficient
##     (Moebius / Owen) formula above that, which needs only the coalitions of
##     size <= order (see .vdnl_shapley_polynomial).  Permutation-sampling
##     Monte Carlo with a reported standard error remains available.
##   * correlated shocks: "active" is defined on the lower-Cholesky
##     orthogonalised innovations in the DECLARED shock order, the convention
##     of compute_moments() / conditional_variance_decomposition().
## --------------------------------------------------------------------------


#' Variance decomposition for pruned second/third-order solutions
#'
#' Decomposes the unconditional (and, optionally, the h-step conditional)
#' variance of every endogenous variable into shock contributions under the
#' pruned order-1, -2 or -3 perturbation solution.  Beyond first order the
#' shocks interact, so the variance is not additive over shocks; the default
#' \code{method = "shapley"} allocates the interaction with the Shapley value,
#' so the shares sum to one exactly.
#'
#' @details
#' \strong{Coalition value.}  Write \eqn{\Sigma_e = L L'} with \eqn{L} the
#' lower Cholesky factor in the DECLARED shock order and
#' \eqn{\varepsilon_t = L u_t}, \eqn{u_t} orthonormal.  For a set \eqn{S} of
#' active orthogonalised shocks, \eqn{v(S)} is the exact variance of each
#' variable under the pruned solution when the innovation covariance is
#' \eqn{L D_S L'} (\eqn{D_S} the 0/1 diagonal selector of \eqn{S}), with the
#' decision rule held fixed at the full-\eqn{\Sigma_e} solution (including
#' the risk-correction terms \code{ghss}, \code{ghxss}, \code{ghuss}).  All
#' state moments are recomputed for each coalition, so \eqn{v(S)} is the
#' variance a pruned simulation with the inactive shocks set to zero would
#' produce, and \eqn{v(\emptyset) = 0}.  With a diagonal \eqn{\Sigma_e},
#' "shock \eqn{k} active" simply means \eqn{\varepsilon_k} is switched on;
#' with correlated shocks the split of the shared variance depends on the
#' \code{varexo} order, exactly as in
#' \code{\link{conditional_variance_decomposition}}.
#'
#' \strong{Methods.}
#' \describe{
#'   \item{\code{"shapley"}}{\eqn{\phi_i = \sum_{S \subseteq N\setminus\{i\}}
#'     \frac{|S|!(n-|S|-1)!}{n!}\,(v(S\cup\{i\}) - v(S))}.  Efficiency gives
#'     \eqn{\sum_i \phi_i = v(N)}, so the shares sum to one.  Individual
#'     shares can fall outside \eqn{[0,1]} because \eqn{v} need not be
#'     monotone.  With at most \code{max_exact} shocks all \eqn{2^n}
#'     coalitions are enumerated; above that the polynomial-coefficient
#'     algorithm below gives the same exact value from far fewer coalitions.
#'     With \code{shapley_algorithm = "permutation"} \eqn{\phi} is instead
#'     estimated by averaging marginal contributions over \code{n_perm}
#'     random shock orderings (seeded, the caller's RNG stream is restored)
#'     and the Monte Carlo standard error of each share is returned in
#'     \code{mc_se}.  Every sampled ordering satisfies efficiency, so the
#'     estimated shares still sum to one exactly.}
#'   \item{\code{"sequential"}}{Each shock switched on on its own, one at a
#'     time: share \eqn{v(\{i\})/v(N)}, plus an \code{interaction} column
#'     \eqn{1 - \sum_i v(\{i\})/v(N)} holding the part of the variance that
#'     only exists when shocks are active together (the default table of
#'     MacroModelling.jl's \code{get_variance_decomposition}).  At first
#'     order the interaction is zero.}
#' }
#' At first order \eqn{v} is additive, so both methods reproduce the linear
#' decomposition of \code{\link{conditional_variance_decomposition}}.
#'
#' \strong{Polynomial-coefficient algorithm.}  Scale each orthogonalised
#' innovation \eqn{u_i} by an indicator \eqn{x_i \in \{0,1\}}.  At pruned
#' order \eqn{k} every variable is a polynomial of degree \eqn{\le k} in the
#' innovations, so each term of its variance is the expectation of a product
#' of at most \eqn{2k} innovations.  The innovations are independent with
#' mean zero, so a term survives only if every innovation index in it
#' appears at least twice: it involves at most \eqn{k} distinct shocks, and
#' because \eqn{x_i^m = x_i} on \eqn{\{0,1\}},
#' \deqn{v(S) = \sum_{T \subseteq N,\ |T| \le k} c_T \prod_{i \in T} x_i}
#' (Moebius / Harsanyi coefficients \eqn{c_T}).  The Shapley value of such a
#' game is \eqn{\phi_i = \sum_{T \ni i} c_T / |T|} (Owen 1972), and
#' \eqn{c_T} for \eqn{|T| \le k} needs \eqn{v} only on coalitions of size
#' \eqn{\le k}.  Substituting the Moebius inversion gives
#' \eqn{\phi_i = \sum_{|U| \le k} \omega(|U|, i \in U)\,(v(U) - v(\emptyset))}
#' with closed-form weights \eqn{\omega}, so the exact Shapley value costs
#' \eqn{\sum_{j \le k} \binom{n}{j} + 1} Lyapunov solves instead of
#' \eqn{2^n} (212 instead of \eqn{2^{20} \approx 10^6} at order 2 with 20
#' shocks).  Every coalition is solved on the reduced system driven by its
#' own \eqn{|S|} orthonormal innovations (\eqn{\varepsilon = L_{\cdot S}
#' u_S}), so the moment blocks of a small coalition are small too.  The
#' extra solve is \eqn{v(N)}: efficiency
#' \eqn{\sum_i \phi_i = v(N)} is not imposed but CHECKED, and its relative
#' violation is returned as \code{closure_error} (a warning of class
#' \code{dynhr_shapley_closure} is raised above \code{1e-8}).
#' MacroModelling.jl quotes the looser bound \eqn{2k} on the coalition size;
#' the zero-mean argument above tightens it to \eqn{k}, and
#' \code{closure_error} would expose a violation.
#'
#' \strong{Horizons.}  \code{Inf} is the unconditional variance (requires a
#' stable \code{hx}).  A finite \eqn{h} is the variance of \eqn{y_{t+h}}
#' conditional on starting at the deterministic steady state
#' (\eqn{x^{(1)} = x^{(2)} = x^{(3)} = 0}), i.e. the h-step forecast-error
#' variance from the steady state; at first order this is the Dynare FEVD of
#' \code{conditional_variance_decomposition}, at higher order it depends on
#' the starting point and the steady state is the one used.
#'
#' Relation to MacroModelling.jl (\code{get_variance_decomposition(...,
#' marginal_contribution = true)}): same Shapley allocation rule and the same
#' per-shock-plus-interaction table for the non-Shapley method, but its
#' coalition value keeps the state cumulants at their full-shock values,
#' whereas here every state moment is recomputed under the coalition.
#'
#' @param dr A \code{DecisionRules} object from
#'   \code{\link{solve_perturbation}} (order 1, 2 or 3; higher-order rules are
#'   used up to their order-3 terms).
#' @param model The parsed model (source of the shock covariance).
#' @param params Named parameter vector (default \code{model$param_values}).
#' @param order Pruning order to decompose at: 1, 2 or 3, at most the order of
#'   \code{dr}.  Default: the order of \code{dr}, capped at 3.
#' @param method \code{"shapley"} (default) or \code{"sequential"}; see
#'   Details.
#' @param horizons Numeric vector of horizons; \code{Inf} for the
#'   unconditional variance.  Default \code{Inf}.
#' @param shapley_algorithm \code{"auto"} (full enumeration when the number
#'   of shocks is at most \code{max_exact}, the polynomial-coefficient
#'   algorithm otherwise), \code{"exact"}, \code{"polynomial"} or
#'   \code{"permutation"}.  \code{"exact"} and \code{"polynomial"} return the
#'   same value (to round-off); see Details.
#' @param n_perm Number of random shock orderings for the permutation
#'   estimator.
#' @param seed Seed for the permutation estimator (\code{NULL}: use the
#'   current RNG stream).
#' @param max_exact Largest number of shocks for full enumeration under
#'   \code{shapley_algorithm = "auto"} (above it, \code{"polynomial"}).
#'
#' @return A list of class \code{"dynhr_vardec_nonlinear"} with
#'   \describe{
#'     \item{shares}{Array [variable x column x horizon] of variance shares;
#'       the columns are the shocks, plus \code{interaction} for
#'       \code{method = "sequential"}.  Each variable's shares sum to one
#'       (rows with zero variance are all zero).}
#'     \item{contributions}{Same array in variance units.}
#'     \item{variance}{Matrix [variable x horizon] of total variances
#'       \eqn{v(N)}.}
#'     \item{mc_se}{Monte Carlo standard errors of \code{shares} for the
#'       permutation estimator, else \code{NULL}.}
#'     \item{order, method, algorithm, n_perm, horizons}{What was computed;
#'       \code{algorithm} is \code{"exact"}, \code{"polynomial"},
#'       \code{"permutation"} or \code{"standalone"} (sequential).}
#'     \item{n_evaluations}{Number of distinct coalitions whose variance was
#'       computed.}
#'     \item{closure_error}{For \code{algorithm = "polynomial"}: the largest
#'       relative violation of efficiency, \eqn{|\sum_i \phi_i - v(N)|}
#'       divided by the largest \eqn{|v(N)|} of the same horizon (round-off
#'       when the degree bound holds); \code{NULL} otherwise.}
#'     \item{Sigma_e, chol_factor}{Shock covariance and its lower Cholesky
#'       factor (the orthogonalisation used).}
#'     \item{var_decomp_pct, var_decomp, var_cov}{When \code{Inf} is among the
#'       horizons: the unconditional shares in percent, the contributions and
#'       the full-coalition covariance matrix, in the layout of
#'       \code{compute_moments()} (so the object can be passed to the D10
#'       variance-decomposition diagnostic).}
#'   }
#'
#' @references
#'   Shapley, L. S. (1953). A value for n-person games. In \emph{Contributions
#'     to the Theory of Games II}, 307-317.
#'   Owen, G. (1972). Multilinear extensions of games. \emph{Management
#'     Science}, 18(5), 64-79.
#'   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F.
#'     (2018). The pruned state-space system for non-linear DSGE models.
#'     \emph{Review of Economic Studies}, 85(1), 1-49.
#'   Kockerols, T. MacroModelling.jl v0.2.0 (2026), "Shapley decompositions
#'     for pruned higher-order perturbation".
#'
#' @seealso \code{\link{conditional_variance_decomposition}} (first order),
#'   \code{\link{compute_moments_order2}}, \code{\link{pruned_ss_moments3}}.
#' @export
variance_decomposition_nonlinear <- function(dr, model, params = NULL,
                                             order = NULL,
                                             method = c("shapley", "sequential"),
                                             horizons = NULL,
                                             shapley_algorithm = c("auto", "exact",
                                                                   "polynomial",
                                                                   "permutation"),
                                             n_perm = 2000L, seed = 1L,
                                             max_exact = 12L) {
  method <- match.arg(method)
  shapley_algorithm <- match.arg(shapley_algorithm)
  if (!inherits(dr, "DecisionRules"))
    .dynhr_abort("variance_decomposition_nonlinear: `dr` must be a ",
                 "DecisionRules object from solve_perturbation().",
                 class = "dynhr_input_error")
  if (is.null(params)) params <- model$param_values

  dr_order <- if (inherits(dr, "DecisionRules3")) 3L
              else if (inherits(dr, "DecisionRules2")) 2L else 1L
  if (is.null(order)) order <- dr_order
  if (length(order) != 1L || !is.numeric(order) || !(order %in% 1:3))
    .dynhr_abort("variance_decomposition_nonlinear: `order` must be 1, 2 or 3.",
                 class = "dynhr_input_error")
  order <- as.integer(order)
  if (order > dr_order)
    .dynhr_abort(sprintf(paste0(
      "variance_decomposition_nonlinear: order = %d requested but `dr` is an ",
      "order-%d solution; re-solve with solve_perturbation(order = %d)."),
      order, dr_order, order), class = "dynhr_input_error")

  if (is.null(horizons)) horizons <- Inf
  if (!is.numeric(horizons) || length(horizons) == 0L || anyNA(horizons) ||
      any(horizons < 1) ||
      any(is.finite(horizons) & horizons != round(horizons)))
    .dynhr_abort("variance_decomposition_nonlinear: `horizons` must be ",
                 "positive integers or Inf.", class = "dynhr_input_error")
  n_perm <- as.integer(n_perm)
  if (length(n_perm) != 1L || is.na(n_perm) || n_perm < 2L)
    .dynhr_abort("variance_decomposition_nonlinear: `n_perm` must be >= 2.",
                 class = "dynhr_input_error")

  endo <- dr$endo_names
  exo  <- dr$exo_names
  n    <- length(exo)
  if (n == 0L)
    .dynhr_abort("variance_decomposition_nonlinear: the model has no shocks.",
                 class = "dynhr_input_error")

  Sigma_e <- .get_shock_cov(model, exo, params)
  L_e     <- .sigma_e_chol_lower(Sigma_e)
  dimnames(L_e) <- list(exo, exo)

  if (any(!is.finite(horizons)) && length(dr$state_idx) > 0L) {
    hx <- dr$ghx[dr$state_idx, , drop = FALSE]
    ev <- Mod(eigen(hx, only.values = TRUE)$values)
    if (any(ev >= 1 - getOption("dynhr.unit_root_tol", 1e-6)))
      .dynhr_abort(sprintf(paste0(
        "variance_decomposition_nonlinear: the unconditional variance does ",
        "not exist (max |eig(hx)| = %.6g). Use finite `horizons`, or ",
        "conditional_variance_decomposition() at first order."), max(ev)),
        class = "dynhr_vardec_nonstationary")
  }

  vfun <- .vdnl_value_fun(dr, order, horizons)
  ## Coalition value: v(active) for a logical vector over the orthogonalised
  ## shocks.  Returns an n_endo x n_h matrix of variances.  A non-empty
  ## coalition is solved on the REDUCED system whose |S| innovations are the
  ## active orthonormal u's (e = L[, S] u_S, Sigma_u = I): the same variance
  ## as masking Sigma_e to L D_S L' (checked in
  ## test-fix-0925-shapley-followups.R), but the order-2/3 moment blocks scale
  ## with |S| instead of n_exo, which is what makes the size-<= order
  ## coalitions of the polynomial algorithm cheap.
  value <- function(active) {
    if (!any(active)) return(vfun(0 * Sigma_e))
    vfun(diag(sum(active)),
         .vdnl_reduce_dr(dr, L_e[, active, drop = FALSE], order))
  }

  hz_names <- as.character(horizons)
  n_h <- length(horizons)
  mc_se  <- NULL
  se_raw <- NULL
  closure <- NULL
  if (method == "sequential") {
    algorithm <- "standalone"
    v_empty <- value(rep(FALSE, n))
    v_full  <- value(rep(TRUE, n))
    contrib <- array(0, c(length(endo), n + 1L, n_h))
    for (i in seq_len(n)) {
      act <- rep(FALSE, n); act[i] <- TRUE
      contrib[, i, ] <- value(act) - v_empty
    }
    contrib[, n + 1L, ] <- (v_full - v_empty) -
      apply(contrib[, seq_len(n), , drop = FALSE], c(1L, 3L), sum)
    col_names <- c(exo, "interaction")
    n_eval <- n + 2L
  } else {
    algorithm <- if (shapley_algorithm == "auto")
      (if (n <= max_exact) "exact" else "polynomial") else shapley_algorithm
    res <- switch(algorithm,
                  exact       = .vdnl_shapley_exact(value, n),
                  polynomial  = .vdnl_shapley_polynomial(value, n, order),
                  permutation = .vdnl_shapley_permutation(value, n, n_perm, seed))
    contrib <- res$phi
    v_full  <- res$v_full
    v_empty <- res$v_empty
    n_eval  <- res$n_eval
    se_raw  <- res$se
    closure <- res$closure_error
    col_names <- exo
    if (!is.null(closure) && closure > 1e-8)
      .dynhr_warn(sprintf(paste0(
        "variance_decomposition_nonlinear: polynomial Shapley closure error ",
        "%.3g (sum of contributions vs v(N)); the degree-%d bound on the ",
        "coalition value does not hold to round-off here. Use ",
        "shapley_algorithm = \"exact\" or \"permutation\"."), closure, order),
        class = "dynhr_shapley_closure")
  }
  dimnames(contrib) <- list(endo, col_names, hz_names)

  total <- v_full - v_empty
  dimnames(total) <- list(endo, hz_names)
  ## Rows with (numerically) zero variance have no decomposition: report zero
  ## shares, as compute_moments() does for nonstationary rows.
  denom <- total
  zero_row <- !(abs(total) > .Machine$double.eps)
  denom[zero_row] <- 1
  shares <- contrib
  for (k in seq_len(dim(contrib)[2L]))
    shares[, k, ] <- ifelse(zero_row, 0, contrib[, k, ] / denom)
  if (!is.null(se_raw)) {
    mc_se <- se_raw
    for (k in seq_len(n))
      mc_se[, k, ] <- ifelse(zero_row, 0, se_raw[, k, ] / denom)
    dimnames(mc_se) <- list(endo, exo, hz_names)
  }

  out <- list(
    shares = shares, contributions = contrib, variance = total,
    mc_se = mc_se, order = order, method = method, algorithm = algorithm,
    n_perm = if (algorithm == "permutation") n_perm else NULL,
    horizons = horizons, n_evaluations = n_eval, closure_error = closure,
    Sigma_e = Sigma_e, chol_factor = L_e
  )
  h_inf <- which(!is.finite(horizons))
  if (length(h_inf) > 0L) {
    h_inf <- h_inf[1L]
    out$var_decomp_pct <- 100 * shares[, , h_inf, drop = TRUE]
    out$var_decomp     <- contrib[, , h_inf, drop = TRUE]
    if (!is.matrix(out$var_decomp_pct)) {
      out$var_decomp_pct <- matrix(out$var_decomp_pct, length(endo),
                                   dimnames = list(endo, col_names))
      out$var_decomp <- matrix(out$var_decomp, length(endo),
                               dimnames = list(endo, col_names))
    }
    vc <- .vdnl_full_cov(dr, order, Sigma_e)
    dimnames(vc) <- list(endo, endo)
    out$var_cov <- vc
  }
  class(out) <- c("dynhr_vardec_nonlinear", "list")
  out
}


## ---- coalition value: variance under a given innovation covariance ----------

## Returns function(Sig, d = dr) -> n_endo x length(horizons) matrix of
## variances of every endogenous variable when the innovations have covariance
## `Sig`, with the decision rule `d` (truncated to `order`) held fixed.  `d`
## may be a .vdnl_reduce_dr() rule with fewer innovations (then `Sig` is its
## innovation covariance).
## @noRd
.vdnl_value_fun <- function(dr, order, horizons) {
  fin   <- is.finite(horizons)
  max_h <- if (any(fin)) as.integer(max(horizons[fin])) else 0L
  n_endo <- length(dr$endo_names)
  collect <- function(v_inf, v_fin) {
    out <- matrix(0, n_endo, length(horizons))
    if (any(!fin)) out[, !fin] <- v_inf
    if (any(fin))  out[, fin]  <- v_fin[, as.integer(horizons[fin]), drop = FALSE]
    out
  }

  if (order == 1L) {
    sidx <- dr$state_idx
    ghx <- dr$ghx
    hx <- ghx[sidx, , drop = FALSE]
    return(function(Sig, d = dr) {
      ghu <- d$ghu
      hu  <- ghu[sidx, , drop = FALSE]
      v_inf <- NULL
      if (any(!fin)) {
        V <- ghu %*% Sig %*% t(ghu)
        if (length(sidx) > 0L)
          V <- V + ghx %*% solve_lyapunov(hx, hu %*% Sig %*% t(hu)) %*% t(ghx)
        v_inf <- diag(V)
      }
      v_fin <- NULL
      if (max_h > 0L) {
        v_fin <- matrix(0, n_endo, max_h)
        P <- matrix(0, length(sidx), length(sidx))
        for (t in seq_len(max_h)) {
          v_fin[, t] <- diag(ghx %*% P %*% t(ghx) + ghu %*% Sig %*% t(ghu))
          P <- hx %*% P %*% t(hx) + hu %*% Sig %*% t(hu)
        }
      }
      collect(v_inf, v_fin)
    })
  }

  if (order == 2L) {
    return(function(Sig, d = dr) {
      sys <- .order2_aug_system(d, Sig)
      v_inf <- if (any(!fin)) diag(.order2_stationary_moments(sys)$var_cov)
      v_fin <- NULL
      if (max_h > 0L) {
        cm <- .order2_conditional_moments(sys, numeric(sys$n_s), max_h)
        v_fin <- vapply(cm$cov, diag, numeric(n_endo))
        v_fin <- matrix(v_fin, n_endo, max_h)
      }
      collect(v_inf, v_fin)
    })
  }

  function(Sig, d = dr) {
    sys <- .order3_aug_system(d, Sig)
    v_inf <- if (any(!fin)) diag(.order3_stationary_moments(sys)$var_cov)
    v_fin <- if (max_h > 0L) .vdnl_order3_transient_var(sys, max_h)
    collect(v_inf, v_fin)
  }
}

## Decision rule re-expressed in k = ncol(Lk) orthonormal innovations u with
## e = Lk u (Sigma_u = I_k).  Substituting into the pruned recursion (same
## Kronecker conventions as simulate_model_order3: e SLOW, states FAST):
##   ghu   -> ghu Lk            ghxu  -> ghxu  (Lk (x) I_s)
##   ghuu  -> ghuu (Lk (x) Lk)  ghxxu -> ghxxu (Lk (x) I_s (x) I_s)
##   ghxuu -> ghxuu (Lk (x) Lk (x) I_s)   ghuuu -> ghuuu (Lk (x) Lk (x) Lk)
##   ghuss -> ghuss Lk
## while the state-only and risk terms (ghx, ghxx, ghxxx, ghss, ghxss, ghs3)
## are unchanged -- the rule stays the full-Sigma_e solution.  Only the fields
## read by .order2_aug_system() / .order3_aug_system() / the order-1 branch
## of .vdnl_value_fun() are carried.
## @noRd
.vdnl_reduce_dr <- function(dr, Lk, order) {
  I_s <- diag(length(dr$state_idx))
  out <- list(ghx = dr$ghx, ghu = dr$ghu %*% Lk, ys = dr$ys,
              state_idx = dr$state_idx, endo_names = dr$endo_names)
  if (order >= 2L) {
    out$ghxx <- dr$ghxx
    out$ghss <- dr$ghss
    out$ghxu <- dr$ghxu %*% kronecker(Lk, I_s)
    out$ghuu <- dr$ghuu %*% kronecker(Lk, Lk)
  }
  if (order >= 3L) {
    LL <- kronecker(Lk, Lk)
    out$ghxxx <- dr$ghxxx
    out$ghxxu <- dr$ghxxu %*% kronecker(Lk, kronecker(I_s, I_s))
    out$ghxuu <- dr$ghxuu %*% kronecker(LL, I_s)
    out$ghuuu <- dr$ghuuu %*% kronecker(LL, Lk)
    out$ghxss <- dr$ghxss
    out$ghs3  <- dr$ghs3
    if (!is.null(dr$ghuss)) out$ghuss <- dr$ghuss %*% Lk
  }
  out
}

## Full-coalition covariance matrix (for the D10-compatible var_cov slot).
## @noRd
.vdnl_full_cov <- function(dr, order, Sigma_e) {
  if (order == 1L) {
    sidx <- dr$state_idx
    V <- dr$ghu %*% Sigma_e %*% t(dr$ghu)
    if (length(sidx) > 0L) {
      hx <- dr$ghx[sidx, , drop = FALSE]; hu <- dr$ghu[sidx, , drop = FALSE]
      V <- V + dr$ghx %*% solve_lyapunov(hx, hu %*% Sigma_e %*% t(hu)) %*%
        t(dr$ghx)
    }
    return(V)
  }
  if (order == 2L)
    return(.order2_stationary_moments(.order2_aug_system(dr, Sigma_e))$var_cov)
  .order3_stationary_moments(.order3_aug_system(dr, Sigma_e))$var_cov
}

## h-step variances (h = 1..max_h) of the order-3 pruned system started at the
## deterministic steady state (xi_0 = 0).  The transient analogue of
## .order3_stationary_moments(): with x1_0 = 0, E[x1_t] = 0 for all t, the raw
## innovation mean is the constant G[, j4] vec(Sigma_e) already in c_u, and
##   Var(y_t)   = Dxi S_t Dxi' + Gv Cr_t Gv' + Dxi Cxr_t Gv' + (.)'
##   S_{t+1}    = Tlin S_t Tlin' + G Cr_t G' + Tlin Cxr_t G' + (.)'
##   mu_{t+1}   = Tlin mu_t + cc + c_u
## with Cr_t = Cov(r_t) and Cxr_t = Cov(xi_t, r_t) evaluated at the CURRENT
## moments of (x1, x2) by the same .order3_cov_r / .order3_cov_xi_r used for
## the stationary point.  As t -> Inf this converges to the stationary
## variance (checked in tests/testthat/test-fix-0925-vardec-nonlinear.R).
## @noRd
.vdnl_order3_transient_var <- function(sys, max_h) {
  d <- sys$d
  mu <- numeric(d)
  S  <- matrix(0, d, d)
  out <- matrix(0, sys$n_endo, max_h)
  for (t in seq_len(max_h)) {
    Cr <- .order3_cov_r(mu[sys$ix1], S[sys$ix1, sys$ix1, drop = FALSE],
                        mu[sys$ix2], S[sys$ix2, sys$ix2, drop = FALSE],
                        sys$Sigma_e,
                        Cov_x2_x11 = S[sys$ix2, sys$ik2, drop = FALSE])
    Cxr <- .order3_cov_xi_r(S[, sys$ix1, drop = FALSE], sys$Sigma_e,
                            sys$jn, d, sys$Dr)
    cy <- sys$Dxi %*% Cxr %*% t(sys$Gv)
    out[, t] <- diag(sys$Dxi %*% S %*% t(sys$Dxi) +
                     sys$Gv %*% Cr %*% t(sys$Gv) + cy + t(cy))
    cs <- sys$Tlin %*% Cxr %*% t(sys$G)
    S  <- sys$Tlin %*% S %*% t(sys$Tlin) + sys$G %*% Cr %*% t(sys$G) + cs + t(cs)
    S  <- (S + t(S)) * 0.5
    mu <- as.numeric(sys$Tlin %*% mu + sys$cc + sys$c_u)
  }
  out
}


## ---- Shapley estimators -------------------------------------------------------

## Exact Shapley value by enumerating all 2^n coalitions.
## `value(active)` -> n_endo x n_h matrix.  Returns phi (n_endo x n x n_h).
## @noRd
.vdnl_shapley_exact <- function(value, n) {
  n_sub <- 2^n
  bits  <- 2^(seq_len(n) - 1L)
  vals  <- vector("list", n_sub)
  for (mask in 0:(n_sub - 1)) {
    active <- bitwAnd(mask, bits) > 0
    vals[[mask + 1L]] <- value(active)
  }
  v_empty <- vals[[1L]]
  v_full  <- vals[[n_sub]]
  ## w(s) = s! (n - s - 1)! / n!, s = |S| = 0..n-1 (lgamma: no overflow).
  w <- exp(lgamma(0:(n - 1L) + 1) + lgamma(n - 0:(n - 1L)) - lgamma(n + 1))
  phi <- array(0, c(dim(v_full)[1L], n, dim(v_full)[2L]))
  for (i in seq_len(n)) {
    acc <- 0 * v_full
    for (mask in 0:(n_sub - 1)) {
      if (bitwAnd(mask, bits[i]) > 0) next
      s <- sum(bitwAnd(mask, bits) > 0)
      acc <- acc + w[s + 1L] * (vals[[mask + bits[i] + 1L]] - vals[[mask + 1L]])
    }
    phi[, i, ] <- acc
  }
  list(phi = phi, v_full = v_full, v_empty = v_empty, n_eval = n_sub)
}

## Weights of the polynomial-coefficient (Moebius / Owen) Shapley formula.
##
## If the game is a multilinear polynomial of degree <= d in the membership
## indicators, v(S) = sum_{|T| <= d, T subset S} c_T with the Moebius
## coefficients c_T = sum_{U subset T} (-1)^{|T|-|U|} v(U), and the Shapley
## value is phi_i = sum_{T contains i} c_T / |T| (Owen 1972).  Swapping the
## two sums, every coalition U with |U| = u <= d enters phi_i with weight
##
##   w_in[u]  = sum_{t=u}^{d}   (-1)^{t-u} C(n-u,   t-u)   / t   (i in U)
##   w_out[u] = sum_{t=u+1}^{d} (-1)^{t-u} C(n-u-1, t-u-1) / t   (i not in U)
##
## (count the T of size t with U u {i} subset T).  A constant added to v has
## zero Moebius coefficient outside the empty set, so v(U) - v(empty) can be
## used in place of v(U).  With d = n these weights reproduce the classical
## Shapley formula for ANY game.  Returns list(w_in, w_out), each length d
## (index u = 1..d).
## @noRd
.shapley_poly_weights <- function(n, d) {
  w_in <- w_out <- numeric(d)
  for (u in seq_len(d)) {
    t_in <- u:d
    w_in[u] <- sum((-1)^(t_in - u) * choose(n - u, t_in - u) / t_in)
    if (u < d) {
      t_out <- (u + 1L):d
      w_out[u] <- sum((-1)^(t_out - u) * choose(n - u - 1, t_out - u - 1) /
                      t_out)
    }
  }
  list(w_in = w_in, w_out = w_out)
}

## All coalitions of size 1..d as an n x M logical matrix (columns ordered by
## size, then lexicographically), plus their sizes.
## @noRd
.shapley_small_coalitions <- function(n, d) {
  blocks <- lapply(seq_len(d), function(u) {
    cmb <- utils::combn(n, u)
    if (!is.matrix(cmb)) cmb <- matrix(cmb, nrow = u)
    act <- matrix(FALSE, n, ncol(cmb))
    act[cbind(as.vector(cmb), rep(seq_len(ncol(cmb)), each = u))] <- TRUE
    act
  })
  act <- do.call(cbind, blocks)
  list(active = act, size = colSums(act))
}

## Exact Shapley value of a game whose coalition value is a polynomial of
## degree <= `degree` in the membership indicators (pruned order-k variances:
## degree k, see the roxygen Details).  Needs v on the sum_{j<=d} C(n, j)
## coalitions of size <= d, plus v(N) for the efficiency (closure) check.
## Same return shape as .vdnl_shapley_exact(), plus `closure_error`.
## @noRd
.vdnl_shapley_polynomial <- function(value, n, degree) {
  d <- min(as.integer(degree), n)
  v_empty <- value(rep(FALSE, n))
  w <- .shapley_poly_weights(n, d)
  co <- .shapley_small_coalitions(n, d)
  phi <- array(0, c(dim(v_empty)[1L], n, dim(v_empty)[2L]))
  v_full <- NULL
  for (k in seq_len(ncol(co$active))) {
    act <- co$active[, k]
    u   <- co$size[k]
    vk  <- value(act)
    if (u == n) v_full <- vk
    dv  <- vk - v_empty
    for (i in seq_len(n))
      phi[, i, ] <- phi[, i, ] + (if (act[i]) w$w_in[u] else w$w_out[u]) * dv
  }
  n_eval <- ncol(co$active) + 1L
  if (is.null(v_full)) {
    v_full <- value(rep(TRUE, n))
    n_eval <- n_eval + 1L
  }
  list(phi = phi, v_full = v_full, v_empty = v_empty, n_eval = n_eval,
       closure_error = .shapley_closure_error(phi, v_full - v_empty))
}

## Relative efficiency violation max |sum_i phi_i - total| / max |total|,
## per horizon (slice 3 of phi), maximised over horizons.  `phi` is
## rows x players x horizons, `total` rows x horizons.
## @noRd
.shapley_closure_error <- function(phi, total) {
  total <- matrix(total, dim(phi)[1L], dim(phi)[3L])
  gap <- abs(apply(phi, c(1L, 3L), sum) - total)
  scale <- apply(abs(total), 2L, max)
  scale[!(scale > 0)] <- 1
  max(sweep(gap, 2L, scale, "/"))
}

## Permutation-sampling Shapley estimator: average marginal contributions over
## `n_perm` uniformly random orderings.  Coalition values are memoised (keyed
## by the 0/1 membership string), so no coalition is solved twice.  Each
## ordering's contributions telescope to v(N) - v(empty), so the estimate is
## efficient exactly; `se` is the Monte Carlo standard error per entry.
## @noRd
.vdnl_shapley_permutation <- function(value, n, n_perm, seed) {
  .local_seed(seed)
  memo <- new.env(parent = emptyenv())
  get_v <- function(active) {
    key <- paste(as.integer(active), collapse = "")
    if (is.null(memo[[key]])) assign(key, value(active), envir = memo)
    memo[[key]]
  }
  v_empty <- get_v(rep(FALSE, n))
  v_full  <- get_v(rep(TRUE, n))
  dims <- c(dim(v_full)[1L], n, dim(v_full)[2L])
  s1 <- array(0, dims); s2 <- array(0, dims)
  for (p in seq_len(n_perm)) {
    perm <- sample.int(n)
    active <- rep(FALSE, n)
    prev <- v_empty
    for (k in perm) {
      active[k] <- TRUE
      cur <- get_v(active)
      dlt <- cur - prev
      s1[, k, ] <- s1[, k, ] + dlt
      s2[, k, ] <- s2[, k, ] + dlt^2
      prev <- cur
    }
  }
  phi <- s1 / n_perm
  vr  <- pmax(s2 / n_perm - phi^2, 0) * n_perm / (n_perm - 1)
  list(phi = phi, se = sqrt(vr / n_perm), v_full = v_full, v_empty = v_empty,
       n_eval = length(ls(memo, all.names = TRUE)))
}
