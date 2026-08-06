## R/hank-jacobian.R
## --------------------------------------------------------------------------
## Sequence-space Jacobian of a heterogeneous-agent household block via the
## fake-news algorithm (Auclert, Bardoczy, Rognlie & Straub 2021, Econometrica).
##
## The block Jacobian J^{o,i}[t, s] = dO_t / dI_s is the response of aggregate
## output o at date t to a (perfect-foresight, anticipated) shock to aggregate
## input i at date s.  Computed naively it needs O(T) full backward+forward
## solves; the fake-news algorithm reduces it to a single backward sweep plus
## cheap bookkeeping (ABRS Prop. 1 / eq. 25):
##
##   1. Backward sweep -> curly-Y (date-0 outcome response to a shock s periods
##      ahead) and curly-D (the induced one-period-ahead distribution change).
##   2. Expectation vectors E_s = Lambda_ss^s y^o  (Lemma 3).
##   3. Fake-news matrix  F[0,s] = curlyY[s];  F[t,s] = <E_{t-1}, curlyD[s]>.
##   4. Jacobian          J[0,s] = F[0,s];  J[t,s] = J[t-1,s-1] + F[t,s].
##
## Per-step derivatives here are taken by central finite differences (the ABRS
## reference uses analytic within-step derivatives; econpizza uses autodiff --
## all three validate against the same brute-force numerical-differentiation
## reference, which is the mandatory oracle in test-hank-jacobian.R).
##
## IF YOU ARE HERE TO MAKE THIS FASTER, OR TO MAKE IT EXACT, READ FIRST:
## ROADMAP.md Tier 19.1 ("Exact het-block structural derivatives: the TANGENT
## SWEEP"), with the measurements in briefs/21-structural-score-api-scope.md
## sec. 11. Short version, so the wrong thing does not get built again:
##   * the het STEADY STATE is 0.4% of a structural FD tap -- a differentiable
##     EGM/steady-state adjoint buys nothing here; this sweep is 85-93% of it;
##   * the directional / IRF-form route to an exact derivative is MEASURED-DEAD
##     (it loses to a plain FD tap at any shock count >= 1);
##   * the tangent sweep is the only live route to exactness and its ceiling is
##     ~2x, because this sweep already central-differences the policy step, so
##     a tangent in a structural parameter needs MIXED SECOND derivatives of
##     .hank_block_step. Its prerequisite -- analytic within-step derivatives,
##     as in the paragraph above -- is worth doing on its own and flips that
##     cost accounting; Tier 19.1 lists the trigger conditions.
## Since 0.9.0.0025 the distributional half of this sweep is matrix-free
## (.hank_forward_push, R/hank-distribution.R); .hank_block_step is now the
## largest single piece of an iteration (57%, was 12%).
##
## Lambda_ss is the row-stochastic forward operator (hank_forward_operator), so
## E_s = Lambda_ss %*% E_{s-1} with NO transpose (E_s(x) = expected future
## outcome from state x), while distributions push forward as t(Lambda) %*% D.
##
## ENDOGENOUS TRANSITION PROBABILITIES (HANK+SAM): when the block carries a
## Pi_fn/Pi_inputs pair (see hank_het_block / hank_employment_income), the
## named transition-probability inputs (e.g. job-finding rate f, separation
## rate s) are perturbable alongside (r, w).  A date-s perturbation of such an
## input x_s moves Pi_s (applied between periods s and s+1) and therefore
## enters the fake-news sweep in TWO places at the shock date:
##   (a) the backward step -- date-s expectations use Pi_s
##       (Wa = beta * Pi_s %*% Va_{s+1}), so policies react at all t <= s via
##       the propagated value-function derivative exactly as for (r, w);
##   (b) the distribution update -- Pi_s enters Lambda_s DIRECTLY, so curly-D
##       at the shock date is the JOINT derivative of t(Lambda(a', Pi(x))) D_ss
##       in (policy, Pi), not the policy-only derivative used for (r, w).
## Steps 2-4 (expectation vectors, fake-news cumulation) are unchanged: they
## are steady-state objects.
## --------------------------------------------------------------------------


#' Validate requested het-block Jacobian inputs
#'
#' The admissible aggregate inputs are the prices \code{c("r", "w")} plus, for
#' a block built with \code{Pi_fn}/\code{Pi_inputs}, the block's named
#' transition-probability inputs (\code{names(block$Pi_inputs)}).
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param inputs Character vector of requested inputs.
#' @return \code{inputs}, validated.
#' @keywords internal
.hank_het_check_inputs <- function(block, inputs) {
  ## Shared by hank_het_jacobian, hank_het_jacobian_nd, hank_het_dist_jacobian
  ## and hank_het_block_spec, so one call covers every one-asset Jacobian route.
  .hank_reject_het2(block, "hank_het_jacobian", use = "hank_het2_jacobian")
  allowed <- c("r", "w", "Tr", "r_minus", names(block$Pi_inputs))
  bad <- setdiff(inputs, allowed)
  if (length(bad))
    stop("unsupported het-block input(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; this block supports ",
         paste0("'", allowed, "'", collapse = ", "),
         if (is.null(block$Pi_inputs))
           " (build the block with Pi_fn/Pi_inputs to add transition-probability inputs)"
         else "", ".")
  inputs
}


## =============================================================================
## E/U/N state and gross-flow OUTPUTS (three-state households, e.g.
## hank_employment_income3): indicator/flow aggregates over the employment
## margin, selectable in hank_het_jacobian()/hank_het_jacobian_nd() alongside
## "A"/"C". A block opts in by carrying idx_E/idx_U/idx_N index vectors (the
## combined-state row ranges hank_employment_income3() returns) -- attach them
## post-hoc if the constructor did not (hank_het_block() itself is untouched
## by this feature).
##
## MATH: for a state output ("E"/"U"/"N") the aggregate O_t = <D_t, y> uses a
## FIXED indicator vector y (1 on the state's rows, tiled across the asset
## grid, 0 elsewhere) that does not depend on the household POLICY at all, so
## curlyY (the date-0/F[0,s] row) is identically zero for every input -- the
## whole response runs through the standard curlyD/expectation-vector
## machinery already in .hank_curly_sweep()/hank_het_jacobian(), unchanged.
##
## For a flow output ("F_xy") O_t = <D_t, g_xy(Pi_t)> with g_xy(Pi)[i] =
## 1{i in X} * rowSum_{j in Y} Pi[i, j] -- gross mass moving X -> Y in period
## t. This has the SAME distribution-propagation part as a state output
## (using g_xy(Pi_ss), a fixed vector, as the E-vector seed), PLUS a
## CONTEMPORANEOUS product-rule term whenever the differentiated input is one
## of the block's Pi-rate inputs: Pi_t is perturbed only at the shock date
## t = s, which changes g_xy(Pi_s) directly (not through the distribution),
## adding a DIAGONAL correction <D_ss, dg_xy/drate> to every J[s, s] (the same
## scalar at every date, since to first order the distribution entering any
## date s is still D_ss along an isolated single-date anticipated shock).
## =============================================================================

#' Fixed indicator vector for a het-block employment-margin STATE output
#' @keywords internal
.hank_het_state_indicator <- function(idx, n_e, n_a) {
  m <- matrix(0, n_e, n_a)
  m[idx, ] <- 1
  .hank_mat_to_vec(m)
}


#' Gross-flow vector g_XY(Pi): mass-per-source-state moving X -> Y under Pi
#' @keywords internal
.hank_het_flow_vector <- function(idx_X, idx_Y, Pi, n_e, n_a) {
  m <- matrix(0, n_e, n_a)
  m[idx_X, ] <- rowSums(Pi[idx_X, idx_Y, drop = FALSE])
  .hank_mat_to_vec(m)
}


#' Named (fromIdxField, toIdxField) pairs for the six E/U/N gross-flow outputs
#' @keywords internal
.hank_het_flow_pairs <- list(
  F_eu = c("idx_E", "idx_U"), F_ue = c("idx_U", "idx_E"),
  F_un = c("idx_U", "idx_N"), F_nu = c("idx_N", "idx_U"),
  F_ne = c("idx_N", "idx_E"), F_en = c("idx_E", "idx_N")
)
.hank_het_state_output_names <- c("E", "U", "N")


#' Is a het-block output name a flow ("F_xy") output?
#' @keywords internal
.hank_het_is_flow_output <- function(o) o %in% names(.hank_het_flow_pairs)


#' Validate requested het-block Jacobian OUTPUTS
#'
#' Extends \code{c("A", "C")} with the state outputs \code{"E"}/\code{"U"}/
#' \code{"N"} and the gross-flow outputs \code{"F_eu"}/\code{"F_ue"}/
#' \code{"F_un"}/\code{"F_nu"}/\code{"F_ne"}/\code{"F_en"}, each available
#' only when \code{block} carries the corresponding \code{idx_E}/\code{idx_U}/
#' \code{idx_N} index vectors (see \code{\link{hank_employment_income3}}), and
#' with \code{"Omega"} -- the distribution-weighted aggregate of the block's
#' OWN transfer-incidence weight, \code{hank_aggregate(D_t, Tr_incidence)} --
#' available whenever \code{block} carries a \code{Tr_incidence} field (every
#' block built by \code{\link{hank_het_block}} does, Tier 1 or Tier 2; see
#' that constructor's \code{Tr_incidence} doc). \code{"Omega"} is what the
#' fiscal outlay \code{Tr_t * Omega_t} is written on when the incidence is
#' \code{(e, a)}-varying and so the outlay is not identically \code{Tr_t}
#' (Tier 2's \code{Omega_ss} is this same quantity's steady-state value).
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param outputs Character vector of requested outputs.
#' @return \code{outputs}, validated.
#' @keywords internal
.hank_het_check_outputs <- function(block, outputs) {
  have <- function(f) !is.null(block[[f]])
  allowed <- c("A", "C")
  if (have("Tr_incidence")) allowed <- c(allowed, "Omega")
  for (s in .hank_het_state_output_names)
    if (have(paste0("idx_", s))) allowed <- c(allowed, s)
  for (fo in names(.hank_het_flow_pairs)) {
    idxs <- .hank_het_flow_pairs[[fo]]
    if (have(idxs[1L]) && have(idxs[2L])) allowed <- c(allowed, fo)
  }
  bad <- setdiff(outputs, allowed)
  if (length(bad))
    stop("unsupported het-block output(s) ",
         paste0("'", bad, "'", collapse = ", "),
         "; this block supports ",
         paste0("'", allowed, "'", collapse = ", "),
         if (!have("idx_E") || !have("idx_U") || !have("idx_N"))
           paste0(" (attach 'idx_E'/'idx_U'/'idx_N' combined-state index ",
                  "vectors, e.g. from hank_employment_income3(), to enable ",
                  "state/flow outputs)")
         else "", ".")
  outputs
}


#' Fixed \code{n_e x n_a} incidence weight, flattened to distribution order
#'
#' The "y" vector for the \code{"Omega"} het-block output: the block's OWN
#' \code{Tr_incidence}, expanded to a full \code{n_e x n_a} matrix if it is
#' still the Tier-1 length-\code{n_e} vector form (constant across the asset
#' dimension, matching the broadcast \code{.hank_income_extra} relies on),
#' and flattened row-major to line up with the distribution.
#' @keywords internal
.hank_het_omega_vector <- function(block) {
  om <- block$Tr_incidence
  om_full <- if (is.matrix(om)) om else matrix(om, block$n_e, block$n_a)
  .hank_mat_to_vec(om_full)
}


#' Fixed "y" vector for a het-block state/flow/Omega OUTPUT at a given Pi
#'
#' \code{o} must be one of the non-"A"/"C" names \code{\link{.hank_het_check_outputs}}
#' allows for \code{block} ("Omega", a state name, or a flow name); dispatches
#' to \code{\link{.hank_het_omega_vector}}, \code{\link{.hank_het_state_indicator}}
#' or \code{\link{.hank_het_flow_vector}}. Like the state outputs, "Omega"'s
#' vector is FIXED (does not depend on the household policy), so it slots into
#' the same curlyY = 0 machinery -- see the file-header math note.
#' @keywords internal
.hank_het_output_vector <- function(block, o, Pi = block$Pi) {
  if (o == "Omega") return(.hank_het_omega_vector(block))
  if (o %in% .hank_het_state_output_names)
    return(.hank_het_state_indicator(block[[paste0("idx_", o)]],
                                     block$n_e, block$n_a))
  pr <- .hank_het_flow_pairs[[o]]
  .hank_het_flow_vector(block[[pr[1L]]], block[[pr[2L]]], Pi,
                        block$n_e, block$n_a)
}


#' Diagonal (contemporaneous) flow correction for a Pi-rate input
#'
#' The product-rule term \code{<D_ss, dg_XY/d(rate)>} (central FD), added to
#' every \code{J[s, s]} entry of a flow output's Jacobian w.r.t. a
#' transition-probability input -- see the file-header math note.
#' @keywords internal
.hank_het_flow_diag_correction <- function(block, o, i, delta_in) {
  pr <- .hank_het_flow_pairs[[o]]
  idx_X <- block[[pr[1L]]]; idx_Y <- block[[pr[2L]]]
  g_p <- .hank_het_flow_vector(idx_X, idx_Y, .hank_pi_perturb(block, i, +delta_in),
                               block$n_e, block$n_a)
  g_m <- .hank_het_flow_vector(idx_X, idx_Y, .hank_pi_perturb(block, i, -delta_in),
                               block$n_e, block$n_a)
  hank_aggregate(block$D, (g_p - g_m) / (2 * delta_in))
}


#' Transition matrix at a perturbed transition-probability input
#'
#' Evaluates \code{block$Pi_fn} with input \code{i} displaced by \code{delta}
#' from its steady-state value (all other transition inputs at steady state).
#' @keywords internal
.hank_pi_perturb <- function(block, i, delta) {
  args <- block$Pi_inputs
  args[[i]] <- args[[i]] + delta
  do.call(block$Pi_fn, args)
}


#' Brute-force numerical-differentiation Jacobian of a het block
#'
#' Reference sequence-space Jacobian: perturbs each aggregate input at each date
#' and re-runs the full nonlinear transition (\code{\link{hank_td_nonlinear}}),
#' central-differencing the aggregate outputs.  Correct but \eqn{O(T)} solves
#' per (input, date); used to validate \code{\link{hank_het_jacobian}}.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param T_h Integer horizon.
#' @param inputs Character subset of \code{c("r", "w", "Tr", "r_minus")} plus,
#'   for a block
#'   built with \code{Pi_fn}/\code{Pi_inputs}, the block's named
#'   transition-probability inputs (\code{names(block$Pi_inputs)}, e.g.
#'   \code{"f"}, \code{"s"} from \code{\link{hank_employment_income}}).
#' @param outputs Character subset of \code{c("A", "C")}, plus, for a block
#'   carrying \code{idx_E}/\code{idx_U}/\code{idx_N} index vectors (see
#'   \code{\link{hank_employment_income3}}), the state outputs \code{"E"} /
#'   \code{"U"} / \code{"N"} and the gross-flow outputs \code{"F_eu"} /
#'   \code{"F_ue"} / \code{"F_un"} / \code{"F_nu"} / \code{"F_ne"} /
#'   \code{"F_en"}.
#' @param delta Numeric FD step for the input perturbation.
#'
#' @return Nested list \code{J[[output]][[input]]}, each a \code{T x T} matrix
#'   with \code{[t, s] = dO_t/dI_s}.
#' @export
hank_het_jacobian_nd <- function(block, T_h,
                                 inputs = c("r", "w"),
                                 outputs = c("A", "C"),
                                 delta = 1e-5) {
  inputs  <- .hank_het_check_inputs(block, inputs)
  outputs <- .hank_het_check_outputs(block, outputs)
  r0 <- rep(block$r, T_h); w0 <- rep(block$w, T_h)
  extra <- setdiff(outputs, c("A", "C"))

  ## State/flow output PATH from a hank_td_nonlinear() result: Dpath is
  ## returned unconditionally; the per-period Pi_t needed for flow outputs is
  ## rebuilt with the SAME .hank_pi_path() helper hank_td_nonlinear() used
  ## internally, so this stays consistent with whichever pi_input_paths (if
  ## any) drove that particular run.
  .extra_paths <- function(out, pip) {
    if (!length(extra)) return(list())
    Pi_path <- .hank_pi_path(block, pip, T_h)
    setNames(lapply(extra, function(o) vapply(seq_len(T_h), function(t) {
      Pi_t <- if (is.null(Pi_path)) block$Pi else Pi_path[[t]]
      hank_aggregate(out$Dpath[, t], .hank_het_output_vector(block, o, Pi_t))
    }, numeric(1))), extra)
  }

  J <- setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) matrix(0, T_h, T_h)), inputs)), outputs)

  for (i in inputs) {
    for (s in seq_len(T_h)) {
      rp <- r0; wp <- w0; rm <- r0; wm <- w0
      pip_p <- NULL; pip_m <- NULL
      Trp <- NULL; Trm <- NULL; rmp <- NULL; rmm <- NULL
      if (i == "r")      { rp[s] <- rp[s] + delta; rm[s] <- rm[s] - delta }
      else if (i == "w") { wp[s] <- wp[s] + delta; wm[s] <- wm[s] - delta }
      else if (i == "Tr") {
        Tr0 <- rep(.hank_block_tr(block), T_h)
        Trp <- Tr0; Trp[s] <- Trp[s] + delta
        Trm <- Tr0; Trm[s] <- Trm[s] - delta
      }
      else if (i == "r_minus") {
        rm_ss <- if (is.null(block$r_minus)) block$r else block$r_minus
        rm0 <- rep(rm_ss, T_h)
        rmp <- rm0; rmp[s] <- rmp[s] + delta
        rmm <- rm0; rmm[s] <- rmm[s] - delta
      }
      else {
        x0 <- rep(block$Pi_inputs[[i]], T_h)
        xp <- x0; xp[s] <- xp[s] + delta
        xm <- x0; xm[s] <- xm[s] - delta
        pip_p <- setNames(list(xp), i)
        pip_m <- setNames(list(xm), i)
      }
      out_p <- hank_td_nonlinear(block, r_path = rp, w_path = wp, T_h = T_h,
                                 pi_input_paths = pip_p, Tr_path = Trp,
                                 r_minus_path = rmp)
      out_m <- hank_td_nonlinear(block, r_path = rm, w_path = wm, T_h = T_h,
                                 pi_input_paths = pip_m, Tr_path = Trm,
                                 r_minus_path = rmm)
      ext_p <- .extra_paths(out_p, pip_p)
      ext_m <- .extra_paths(out_m, pip_m)
      for (o in outputs) {
        vp <- if (o %in% c("A", "C")) out_p[[o]] else ext_p[[o]]
        vm <- if (o %in% c("A", "C")) out_m[[o]] else ext_m[[o]]
        J[[o]][[i]][, s] <- (vp - vm) / (2 * delta)
      }
    }
  }
  J
}


#' Brute-force numerical-differentiation DISTRIBUTION Jacobian of a het block
#'
#' Reference sequence-space distribution Jacobian \eqn{J^D[t,s,] = dD_t/dI_s}
#' (the change in the full \code{(n_e*n_a)}-vector cross-sectional distribution
#' at date \code{t} induced by a shock to aggregate input \code{i} at date
#' \code{s}), computed exactly like \code{\link{hank_het_jacobian_nd}} but
#' keeping the full \code{Dpath} instead of aggregating it into \code{A}/\code{C}.
#' Used to validate \code{\link{hank_het_dist_jacobian}}.
#'
#' @inheritParams hank_het_jacobian_nd
#'
#' @return Named list \code{JD_nd[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_a)} with \code{JD_nd[[i]][t, s, ] =
#'   dD_t/dI_s} (central difference).
#' @export
hank_het_dist_jacobian_nd <- function(block, T_h,
                                      inputs = c("r", "w"),
                                      delta = 1e-5) {
  inputs <- .hank_het_check_inputs(block, inputs)
  r0 <- rep(block$r, T_h); w0 <- rep(block$w, T_h)
  n_cell <- block$n_e * block$n_a

  JD_nd <- setNames(lapply(inputs, function(i) array(0, c(T_h, T_h, n_cell))),
                    inputs)

  for (i in inputs) {
    for (s in seq_len(T_h)) {
      rp <- r0; wp <- w0; rm <- r0; wm <- w0
      pip_p <- NULL; pip_m <- NULL
      Trp <- NULL; Trm <- NULL; rmp <- NULL; rmm <- NULL
      if (i == "r")      { rp[s] <- rp[s] + delta; rm[s] <- rm[s] - delta }
      else if (i == "w") { wp[s] <- wp[s] + delta; wm[s] <- wm[s] - delta }
      else if (i == "Tr") {
        Tr0 <- rep(.hank_block_tr(block), T_h)
        Trp <- Tr0; Trp[s] <- Trp[s] + delta
        Trm <- Tr0; Trm[s] <- Trm[s] - delta
      }
      else if (i == "r_minus") {
        rm_ss <- if (is.null(block$r_minus)) block$r else block$r_minus
        rm0 <- rep(rm_ss, T_h)
        rmp <- rm0; rmp[s] <- rmp[s] + delta
        rmm <- rm0; rmm[s] <- rmm[s] - delta
      }
      else {
        x0 <- rep(block$Pi_inputs[[i]], T_h)
        xp <- x0; xp[s] <- xp[s] + delta
        xm <- x0; xm[s] <- xm[s] - delta
        pip_p <- setNames(list(xp), i)
        pip_m <- setNames(list(xm), i)
      }
      out_p <- hank_td_nonlinear(block, r_path = rp, w_path = wp, T_h = T_h,
                                 pi_input_paths = pip_p, Tr_path = Trp,
                                 r_minus_path = rmp)
      out_m <- hank_td_nonlinear(block, r_path = rm, w_path = wm, T_h = T_h,
                                 pi_input_paths = pip_m, Tr_path = Trm,
                                 r_minus_path = rmm)
      dD <- (out_p$Dpath - out_m$Dpath) / (2 * delta)   # (n_cell x T_h), col t
      for (tt in seq_len(T_h)) JD_nd[[i]][tt, s, ] <- dD[, tt]
    }
  }
  JD_nd
}


#' Backward sweep shared by \code{hank_het_jacobian} and
#' \code{hank_het_dist_jacobian}: curly-Y (per-output date-0 outcome response)
#' and curly-D (induced one-period-ahead distribution change) to an anticipated
#' shock to input \code{i} at horizon \code{s = 1 .. T_h} (s=1 is the direct
#' current-period shock; s>=2 propagate via the value-function derivative).
#'
#' For a transition-probability input (\code{i} in
#' \code{names(block$Pi_inputs)}), the s=1 term perturbs Pi in BOTH places it
#' enters the shock date (see the file header): the backward step's
#' expectation uses \code{Pi_fn(x +/- delta_in)}, and curly-D is the JOINT
#' central difference of \eqn{t(\Lambda(a', \Pi))\,D_{ss}} along the direction
#' (policy change \code{dA}, unit input change) -- the policy-only
#' \code{curlyD_from_dA} would miss Pi's direct entry in the forward law of
#' motion.  The s>=2 anticipation terms are unchanged: at dates before the
#' shock both the step's Pi and Lambda's Pi sit at steady state.
#'
#' @return List with \code{curlyY} (named list over \code{outputs}, each a
#'   length-\code{T_h} vector) and \code{curlyD} (\code{(n_e*n_a) x T_h}
#'   matrix, column \code{s}).
#' @keywords internal
.hank_curly_sweep <- function(block, T_h, i, outputs,
                              delta_in, delta_va, delta_d) {
  a_grid <- block$a_grid; Pi <- block$Pi; D_ss <- block$D
  Va_ss <- block$Va; a_ss <- block$a
  is_pi_input <- !(i %in% c("r", "w", "Tr", "r_minus"))

  ## Helper: distributional response (curly-D) to a savings-policy change dA.
  ## MATRIX-FREE: this ran twice per date and was 78-83% of the whole sweep
  ## when it built two sparse Lambdas for a product that needs none (measured
  ## in .hank_forward_push's header). Same algebra, contracted rather than
  ## materialized -- parity with the sparse path at round-off.
  curlyD_from_dA <- function(dA) {
    (.hank_forward_push(a_ss + delta_d * dA, a_grid, Pi, D_ss) -
       .hank_forward_push(a_ss - delta_d * dA, a_grid, Pi, D_ss)) / (2 * delta_d)
  }

  curlyY <- setNames(lapply(outputs, function(o) numeric(T_h)), outputs)
  curlyD <- matrix(0, block$n_e * block$n_a, T_h)

  ## s = 1: direct input shock at the current date.
  if (i == "r") {
    sp <- .hank_block_step(block, Va_ss, block$r + delta_in, block$w)
    sm <- .hank_block_step(block, Va_ss, block$r - delta_in, block$w)
  } else if (i == "w") {
    sp <- .hank_block_step(block, Va_ss, block$r, block$w + delta_in)
    sm <- .hank_block_step(block, Va_ss, block$r, block$w - delta_in)
  } else if (i == "r_minus") {
    ## The borrowing-rate input. Around a SYMMETRIC block this perturbs
    ## r_minus away from r in both directions (the step dispatches to the
    ## wedge solver as soon as r_minus is non-NULL), which is exactly the
    ## directional derivative the DAG composition rb -> repricing -> r_minus
    ## needs at a zero-wedge steady state.
    rm0 <- if (is.null(block$r_minus)) block$r else block$r_minus
    sp <- .hank_block_step(block, Va_ss, block$r, block$w,
                           r_minus = rm0 + delta_in)
    sm <- .hank_block_step(block, Va_ss, block$r, block$w,
                           r_minus = rm0 - delta_in)
  } else if (i == "Tr") {
    ## Lump-sum transfer: enters the budget additively (y = w*e + Tr*omega),
    ## so its s = 1 column is the iMPC out of a date-1 transfer -- see
    ## hank_impc(), an independent implementation this is cross-validated
    ## against in test-hank-transfer.R.
    ##
    ## The incidence weight omega needs NO branch of its own here: it is block
    ## state, not a perturbable input, so .hank_block_step applies it and this
    ## stays a plain central difference in the scalar Tr. Nor does it reach
    ## curlyD below -- omega shifts income within a period but not the state
    ## coordinate, so the policy-only curlyD_from_dA branch remains correct
    ## (contrast theta_coll in .hank_curly_sweep2, which DOES move the
    ## coordinate and needs the date-1/date-2 rebasing corrections).
    Tr0 <- .hank_block_tr(block)
    sp <- .hank_block_step(block, Va_ss, block$r, block$w, Tr = Tr0 + delta_in)
    sm <- .hank_block_step(block, Va_ss, block$r, block$w, Tr = Tr0 - delta_in)
  } else {
    ## transition-probability input: the date-1 expectation uses perturbed Pi
    sp <- .hank_block_step(block, Va_ss, block$r, block$w,
                           Pi = .hank_pi_perturb(block, i, +delta_in))
    sm <- .hank_block_step(block, Va_ss, block$r, block$w,
                           Pi = .hank_pi_perturb(block, i, -delta_in))
  }
  dA  <- (sp$a  - sm$a)  / (2 * delta_in)
  dC  <- (sp$c  - sm$c)  / (2 * delta_in)
  dVa <- (sp$Va - sm$Va) / (2 * delta_in)
  if ("A" %in% outputs) curlyY[["A"]][1L] <- hank_aggregate(D_ss, dA)
  if ("C" %in% outputs) curlyY[["C"]][1L] <- hank_aggregate(D_ss, dC)
  if (is_pi_input) {
    ## Joint (policy, Pi) directional derivative of the forward update: Pi
    ## enters Lambda directly on the shock date, so perturb both together
    ## with the SAME step (unit direction in the input, dA in the policy).
    curlyD[, 1L] <-
      (.hank_forward_push(a_ss + delta_d * dA, a_grid,
                          .hank_pi_perturb(block, i, +delta_d), D_ss) -
         .hank_forward_push(a_ss - delta_d * dA, a_grid,
                            .hank_pi_perturb(block, i, -delta_d), D_ss)) /
      (2 * delta_d)
  } else {
    curlyD[, 1L] <- curlyD_from_dA(dA)
  }

  ## s >= 2: propagate the anticipation via the value-function derivative.
  dVa_prev <- dVa
  for (s in seq_len(T_h - 1L) + 1L) {              # s = 2 .. T_h, empty if T_h < 2
    h  <- delta_va / max(1, max(abs(dVa_prev)))
    sp <- .hank_block_step(block, Va_ss + h * dVa_prev, block$r, block$w)
    sm <- .hank_block_step(block, Va_ss - h * dVa_prev, block$r, block$w)
    dA  <- (sp$a  - sm$a)  / (2 * h)
    dC  <- (sp$c  - sm$c)  / (2 * h)
    dVa <- (sp$Va - sm$Va) / (2 * h)
    if ("A" %in% outputs) curlyY[["A"]][s] <- hank_aggregate(D_ss, dA)
    if ("C" %in% outputs) curlyY[["C"]][s] <- hank_aggregate(D_ss, dC)
    curlyD[, s] <- curlyD_from_dA(dA)
    dVa_prev <- dVa
  }

  list(curlyY = curlyY, curlyD = curlyD)
}


#' Sequence-space Jacobian of a het block via the fake-news algorithm
#'
#' For a block built with \code{Pi_fn}/\code{Pi_inputs} (e.g. via
#' \code{\link{hank_employment_income}}), the transition-probability inputs
#' are perturbable alongside \code{(r, w)}: pass their names in \code{inputs}
#' to get the Jacobian columns w.r.t. anticipated job-finding/separation-rate
#' paths (see the file header for how the Pi perturbation enters the sweep).
#'
#' @inheritParams hank_het_jacobian_nd
#' @param delta_in FD step for the input (r/w/transition-probability)
#'   perturbation in the s=0 term.
#' @param delta_va Relative FD step for the backward value-function propagation.
#' @param delta_d FD step for the distributional (curly-D) response.
#'
#' @return Nested list \code{J[[output]][[input]]}, each a \code{T x T} matrix
#'   with \code{[t, s] = dO_t/dI_s}.  (The date-0 response vectors are the first
#'   Jacobian row, \code{J[[o]][[i]][1, ]}.)
#' @export
hank_het_jacobian <- function(block, T_h,
                              inputs = c("r", "w"),
                              outputs = c("A", "C"),
                              delta_in = 1e-5, delta_va = 1e-6,
                              delta_d = 1e-6) {
  inputs  <- .hank_het_check_inputs(block, inputs)
  outputs <- .hank_het_check_outputs(block, outputs)
  Lam <- block$Lambda

  ## --- Step 2: expectation vectors E_s = Lambda^s y^o, s = 0 .. T-1 ---
  ## "A"/"C" seed on the per-agent policy (unchanged); state ("E"/"U"/"N")
  ## and flow ("F_xy") outputs seed on a FIXED steady-state indicator/flow
  ## vector -- see the file-header math note above .hank_het_state_indicator.
  y_out <- setNames(lapply(outputs, function(o) {
    if (o == "A") .hank_mat_to_vec(block$a)        # per-agent savings
    else if (o == "C") .hank_mat_to_vec(block$c)   # per-agent consumption
    else .hank_het_output_vector(block, o)         # state/flow indicator
  }), outputs)
  Elist <- setNames(vector("list", length(outputs)), outputs)
  for (o in outputs) {
    E <- vector("list", T_h)
    E[[1L]] <- y_out[[o]]                        # E_0
    ## s = 2 .. T_h, empty if T_h < 2
    for (s in seq_len(T_h - 1L) + 1L) E[[s]] <- as.numeric(Lam %*% E[[s - 1L]])
    Elist[[o]] <- E
  }

  J <- setNames(lapply(outputs, function(o)
    setNames(lapply(inputs, function(i) matrix(0, T_h, T_h)), inputs)), outputs)

  for (i in inputs) {
    ## --- Step 1: backward sweep -> curlyY[[o]][s], curlyD[, s] ---
    sweep  <- .hank_curly_sweep(block, T_h, i, outputs,
                                delta_in, delta_va, delta_d)
    curlyY <- sweep$curlyY
    curlyD <- sweep$curlyD

    ## --- Steps 3-4: fake-news matrix F then Jacobian J ---
    for (o in outputs) {
      Fm <- matrix(0, T_h, T_h)
      Fm[1L, ] <- curlyY[[o]]                       # F[0, s] = curlyY[s]
      E <- Elist[[o]]
      for (tt in seq_len(T_h - 1L) + 1L) {           # tt = 2 .. T_h, empty if T_h < 2
        Etm1 <- E[[tt - 1L]]                         # E_{t-1}
        Fm[tt, ] <- as.numeric(crossprod(curlyD, Etm1))  # <E_{t-1}, curlyD[,s]>
      }
      ## Diagonal cumulative sum.
      Jm <- matrix(0, T_h, T_h)
      Jm[1L, ] <- Fm[1L, ]
      for (tt in seq_len(T_h - 1L) + 1L) {           # tt = 2 .. T_h, empty if T_h < 2
        Jm[tt, 1L] <- Fm[tt, 1L]
        ## Body only runs when T_h >= 2, so the 2L:T_h slices are in range here.
        Jm[tt, 2L:T_h] <- Jm[tt - 1L, 1L:(T_h - 1L)] + Fm[tt, 2L:T_h]
      }
      ## Flow outputs get an extra CONTEMPORANEOUS product-rule term for
      ## Pi-rate inputs: Pi_t is perturbed only at the shock date t = s,
      ## which moves g_XY(Pi_s) directly -- see the file-header math note.
      if (.hank_het_is_flow_output(o) && i %in% names(block$Pi_inputs)) {
        corr <- .hank_het_flow_diag_correction(block, o, i, delta_in)
        diag(Jm) <- diag(Jm) + corr
      }
      J[[o]][[i]] <- Jm
    }
  }
  J
}


#' Sequence-space DISTRIBUTION Jacobian of a het block via the fake-news
#' algorithm
#'
#' Extends the fake-news algorithm (\code{\link{hank_het_jacobian}}) to expose
#' the full distributional response \eqn{J^D[t, s, ] = dD_t/dI_s} (the change
#' in the \code{(n_e*n_a)}-vector cross-sectional distribution at date \code{t}
#' from an anticipated shock to input \code{i} at date \code{s}), instead of
#' aggregating it into scalar outputs \code{A}/\code{C}.
#'
#' Reuses the same backward sweep (\code{curlyD}) as \code{hank_het_jacobian}.
#' Distributions push forward under the TRANSPOSE of the steady-state forward
#' operator (\code{t(Lambda) \%*\% d}; see the file header), so the
#' distribution fake-news matrix cumulates by repeatedly applying
#' \code{t(Lambda)} to \code{curlyD[,s]} rather than by dotting against an
#' expectation vector (the aggregate-output equivalent of that projection).
#'
#' TIMING: \code{D_t} is the distribution ENTERING period \code{t} (a
#' predetermined state), so \code{D_1 = D_ss} always and row \code{t = 1} of
#' \code{J^D} is identically zero for every shock date \code{s} -- no
#' anticipated shock can move the very first, already-fixed, initial
#' distribution.  \code{curlyD[, s]} is the response of the policy used IN
#' period \code{s} (see \code{\link{.hank_curly_sweep}}: its \code{s = 1} term
#' is built from a policy step at the current period, exactly as
#' \code{curlyY[1]} is in \code{hank_het_jacobian}), which the forward operator
#' then turns into a distribution change one calendar period later, at
#' \code{t = s + 1}.  So the fake-news distribution response first appears at
#' \code{t = s + 1}, not \code{t = s} as for the aggregate outputs -- i.e. the
#' whole cumulation is the aggregate-Jacobian recursion shifted down by one row
#' (see the ND-vs-fake-news diagnostic in test-hank-dist-jacobian.R, which
#' pinned this down to a clean off-by-one before this comment was written).
#'
#' @inheritParams hank_het_jacobian
#'
#' @return Named list \code{JD[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_a)} with \code{JD[[i]][t, s, ] = dD_t/dI_s}.
#' @export
hank_het_dist_jacobian <- function(block, T_h,
                                   inputs = c("r", "w"),
                                   delta_in = 1e-5, delta_va = 1e-6,
                                   delta_d = 1e-6) {
  ## own internal sweep, not yet wedge-aware (it would price the debt side at
  ## the saving rate) -- refuse loudly rather than return plausible nonsense
  .hank_reject_wedge(block, "hank_het_dist_jacobian")
  inputs <- .hank_het_check_inputs(block, inputs)
  if ("r_minus" %in% inputs)
    stop("hank_het_dist_jacobian: the 'r_minus' column is not wired for the ",
         "distribution Jacobian (its internal sweep is wedge-unaware).")
  Lam <- block$Lambda
  n_cell <- block$n_e * block$n_a
  P <- function(x) as.numeric(Matrix::t(Lam) %*% x)   ## distribution push

  JD <- setNames(lapply(inputs, function(i) array(0, c(T_h, T_h, n_cell))),
                inputs)
  if (T_h < 2L) return(JD)   ## row t=1 (D_ss, fixed) is the only row; all-zero

  for (i in inputs) {
    ## --- Step 1: backward sweep -> curlyD[, s] (curlyY unused here) ---
    sweep  <- .hank_curly_sweep(block, T_h, i, outputs = "A",
                                delta_in, delta_va, delta_d)
    curlyD <- sweep$curlyD

    ## --- Step 3: distribution fake-news F^D, indexed by CALENDAR date t ---
    ## FD[[2]][, s] = curlyD[, s]   (first possible response: t = s + 1 = 2
    ##                                relative offset 1, hit when s = 1)
    ## FD[[t]][, s] = P(FD[[t-1]][, s])   for t >= 3
    ## (FD[[1]] would be t=1, always zero -- omitted; loop starts at t=2.)
    FD <- vector("list", T_h)
    FD[[2L]] <- curlyD
    for (tt in seq_len(T_h - 2L) + 2L) {              # tt = 3 .. T_h, empty if T_h < 3
      prev <- FD[[tt - 1L]]
      cur  <- matrix(0, n_cell, T_h)
      for (s in seq_len(T_h)) cur[, s] <- P(prev[, s])
      FD[[tt]] <- cur
    }

    ## --- Step 4: diagonal cumulation, vector-valued per (t, s) ---
    ## Mirrors hank_het_jacobian's J[t,s] = J[t-1,s-1] + F[t,s] exactly, just
    ## shifted down one row (t=1 row is the fixed, unresponsive D_ss and stays
    ## all-zero; the recursion proper starts at t=2).
    JD[[i]][2L, 1L, ] <- FD[[2L]][, 1L]               # JD[2, 1, ] = FD[2][, 1]
    for (s in seq_len(T_h - 1L) + 1L)                 # s = 2 .. T_h
      JD[[i]][2L, s, ] <- FD[[2L]][, s]               # JD[2,s,]=FD[2][,s] (t-1=1 row is 0)
    for (tt in seq_len(T_h - 2L) + 2L) {               # tt = 3 .. T_h, empty if T_h < 3
      JD[[i]][tt, 1L, ] <- FD[[tt]][, 1L]             # JD[t, 1, ] = FD[t][, 1]
      for (s in seq_len(T_h - 1L) + 1L) {              # s = 2 .. T_h
        JD[[i]][tt, s, ] <- JD[[i]][tt - 1L, s - 1L, ] + FD[[tt]][, s]
      }
    }
  }
  JD
}
