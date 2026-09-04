## R/ms-smoother.R
## --------------------------------------------------------------------------
## ms_kim_smoother() -- Kim (1994) backward smoother for Markov-switching DSGE,
## with a DURBIN-KOOPMAN state pass.
##
## Companion to ms_kim_filter() (R/ms-filter.R).  Produces
##   * smoothed regime probabilities  Pr[s_t = j | y_{1:T}]
##   * smoothed states                E[s_t | y_{1:T}]  (and their covariances)
##   * smoothed structural shocks     E[eps_t | y_{1:T}]
##
## References: Kim, C.-J. (1994), "Dynamic linear models with Markov-switching",
##   Journal of Econometrics 60(1-2), 1-22 (the regime-probability pass and the
##   h^2 -> h collapse); Durbin, J. & Koopman, S. J. (2012), "Time Series
##   Analysis by State Space Methods", 2nd ed., sections 4.4-4.5 and 6.2 (the
##   state pass).
##
## THE FORWARD PASS IS THE FILTER.  This file never re-derives a Kim forward
## recursion, and never re-derives an innovation: it calls
## ms_kim_filter(return_state_path = TRUE) and smooths exactly the moments and
## the (v, F^{-1}, K) blocks that filter produced.  (Forward-consistency trap:
## a hand-copied forward pass silently drifts from the likelihood it is
## supposed to decompose the moment the filter is edited.)
##
## WHY DURBIN-KOOPMAN AND NOT RAUCH-TUNG-STRIEBEL (2026-09-03).
## RTS assumes y_{t+1} is conditionally independent of s_t given s_{t+1}.
## Under dynhr's lag-1 observation timing that is FALSE:
##     y_{t+1} = Z s_t + D eps_{t+1}
## loads s_t DIRECTLY, so an RTS state pass is not exact here.  This was first
## established for the single-regime smoother (E4-A, R/smoother-monolith.R,
## tests/testthat/test-kalman-smoother-exact.R, which pins it against a
## brute-force joint-Gaussian projection oracle); the same argument applies
## verbatim to every regime path of the Kim smoother, since each path is an
## ordinary Kalman system with the same timing.  This smoother's first version
## used the RTS form and missed kalman_smoother() by 4.7e-6 relative on the
## states and 5.1e-8 on the covariances under P = I -- small, but a genuine
## specification error, not round-off.  With the DK pass both agree to ~1e-13.
##
## THE REGIME PASS: WHY "joint" AND NOT KIM'S (F2-A, 2026-09-03).
##
## Kim (1994, eq. 10) builds the smoothed joint by
##     Pr[s_t=i, s_{t+1}=j | y_{1:T}]
##       ~= Pr[s_{t+1}=j | y_{1:T}] * Pr[s_t=i | s_{t+1}=j, y_{1:t}]        (K)
##       =  Pr[s_{t+1}=j | y_{1:T}] * Pr[s_t=i|y_{1:t}] P[i,j]
##                                     / Pr[s_{t+1}=j | y_{1:t}] ,
## i.e. it replaces the exact Pr[s_t=i | s_{t+1}=j, y_{1:T}] by the
## FILTERED-through-t conditional.  The approximation Kim needs is
##     y_{t+1:T}  _||_  s_t  |  s_{t+1}, y_{1:t} .                          (KA)
## In a standard state-space model (y_t = Z s_t + ...) that is only the usual
## GPB(2) collapse error.  Under dynhr's LAG-1 observation timing it is much
## worse, because
##     y_{t+1} = d_i + ZZ_i s_t + DD_j eps_{t+1}
## loads the regime i = s_t DIRECTLY (through ZZ_i and the intercept d_i), so
## y_{t+1} is a first-order-informative signal about s_t that (KA) discards.
## An orchestrator brute-force oracle (all 2^{T+1} regime paths enumerated on
## T = 8 with an exact joint-Gaussian projection per path) measured the shipped
## pass at up to 0.31 away from the exact smoothed regime probability under
## genuine switching -- while agreeing with Kim's own formula to 1e-16, i.e.
## the code was right and the FORMULA was the error.
##
## THE JOINT PASS uses the filter's own posterior joint one period later:
##     Pr[s_t=i, s_{t+1}=j | y_{1:T}]
##       ~= Pr[s_{t+1}=j | y_{1:T}] * Pr[s_t=i | s_{t+1}=j, y_{1:t+1}]      (J)
##       =  Pr[s_{t+1}=j | y_{1:T}]
##            * Pr[s_t=i, s_{t+1}=j | y_{1:t+1}] / Pr[s_{t+1}=j | y_{1:t+1}]
## and its numerator/denominator are BOTH exactly the Hamilton filter's
## `prob_joint` at period t+1 and its column sum -- nothing new is computed,
## the filter simply stops throwing them away (`joint_filt` in R/ms-filter.R).
## The conditional independence still assumed is
##     y_{t+2:T}  _||_  s_t  |  s_{t+1}, y_{1:t+1} .                        (JA)
## (JA) is STRICTLY WEAKER than (KA): it conditions on the same regime plus a
## strictly larger information set (y_{1:t+1} rather than y_{1:t}), and the one
## observation it adds, y_{t+1}, is exactly the one that loads s_t directly.
## Equivalently: (KA) must pretend y_{t+1} says nothing about s_t; (JA) need
## only pretend that y_{t+2:T} says nothing MORE about s_t once s_{t+1} and
## y_{1:t+1} are known -- which is the ordinary GPB(2) collapse error, the same
## error the FILTER already makes.  Neither pass is exact (both remain h^2 -> h
## collapses); the joint pass is exact in every case where Kim's is, and the
## degenerate case P = I makes the two identical, because both conditionals
## then reduce to delta_{ij}.
##
## The STATE pass consumes the SAME joint weights (w_from = column-normalised
## joint, w_to = row-normalised joint), so choosing the regime pass moves the
## smoothed states and shocks with the probabilities rather than leaving them
## on a different weighting.
##
## The FORWARD FILTER IS UNTOUCHED, so `loglik` is bit-identical across
## `regime_pass` (asserted in test-ms-smoother-joint.R).
##
## THE RECURSION.  Write alpha_t := s_{t-1}, so that
##     alpha_{t+1} = T alpha_t + R eps_t ,  y_t = Z alpha_t + D eps_t
## is DK's correlated-noise form.  DK's PREDICTION pair for alpha_t is this
## filter's UPDATED pair (s_{t-1|t-1}, P_{t-1|t-1}).  Per regime j, carry the
## adjoint (r_t^{(j)}, N_t^{(j)}) summarising y_{t+1:T} about s_t given
## s_t = j, with r_T = 0 and N_T = 0.  For a transition j (at t) -> k (at t+1),
## every block below is a stored by-product of the filter's own period-(t+1)
## measurement update on path (from = j, to = k):
##
##   r^{(j,k)} = a_{t+1}^{(j,k)} + L_{t+1}^{(j,k)'} r_{t+1}^{(k)}
##   N^{(j,k)} = M_{t+1}^{(j,k)} + L_{t+1}^{(j,k)'} N_{t+1}^{(k)} L_{t+1}^{(j,k)}
##   b^{(j,k)} = b_{t|t}^{(j)} + P_{t|t}^{(j)} r^{(j,k)}
##   P^{(j,k)} = P_{t|t}^{(j)} - P_{t|t}^{(j)} N^{(j,k)} P_{t|t}^{(j)}
##
## with a = Z' F^{-1} v, M = Z' F^{-1} Z, L = T - K Z.  These are then
## collapsed over k with Kim's smoothed transition weights
## w_jk = Pr[s_t=j, s_{t+1}=k | y_T] / Pr[s_t=j | y_T], the covariance carrying
## the usual across-path spread term.  Nothing observes s_T, so the smoothed
## pair at t = T is the FILTERED pair -- which is also what r_T = N_T = 0 gives.
##
## Smoothed shocks ride the same adjoint (DK section 4.4, correlated-noise form):
##   eps_{t|T}^{(i,j)} = Sigma_e^{(j)} ( G_t^{(i,j)} r_t^{(j)} + du_t^{(i,j)} )
## with G = (R - K D)' and du = D' F^{-1} v, again both stored by the filter.
## (Derivation: eps = Q(R' r + D' u), u = F^{-1} v - K' r, so the r-loading is
## R' - D'K' = (R - K D)'.)
##
## --------------------------------------------------------------------------
## THE GPB(3) BACKWARD PASS (F4-D, 2026-09-04): PAIR-INDEXED, TRIPLE JOINT.
##
## ms_kim_filter(collapse = "gpb3") keeps h^2 Gaussian components indexed by
## the PAIR (s_{t-1}, s_t) and collapses h^3 -> h^2 each period.  Its backward
## analogue is not a re-indexing of the pass above: it carries one
## Durbin-Koopman adjoint PER PAIR and collapses with the h^3 joint
## Pr[s_{t-2}, s_{t-1}, s_t | y_{1:T}].  Written out (h = the number of
## regimes, all indices 1..h):
##
##   COMPONENTS.  At period t the filter's surviving component g = (j, k)
##   [stored at j + (k-1)h] is the law of s_t given {s_{t-1}=j, s_t=k, y_{1:t}}
##   with moments (beta_filt[, g, t], P_filt[, , g, t]).  The period-t update
##   CELL is a TRIPLE (i, j, k) = (s_{t-2}, s_{t-1}, s_t), stored column-major
##   at m + (k-1)M with m = i + (j-1)h the incoming component; its post-update
##   group is g = (j, k), i.e. i is the coordinate collapsed away.  (At t = 1
##   there is no s_{-1}: M = h, the cell is the pair (i, k) and nothing is
##   collapsed.)
##
##   REGIME PASS ON TRIPLES.  Exactly,
##     Pr[s_{t-2}=i, s_{t-1}=j, s_t=k | y_{1:T}]
##       = Pr[s_{t-1}=j, s_t=k | y_{1:T}] * Pr[s_{t-2}=i | (j,k), y_{1:T}] ,
##   and the pass replaces the last conditional by the FILTERED one at the SAME
##   period t,
##     Pr[s_{t-2}=i | s_{t-1}=j, s_t=k, y_{1:t}]
##       = Pr[i, j, k | y_{1:t}] / Pr[j, k | y_{1:t}] ,                    (J3)
##   whose numerator is the filter's own per-cell posterior (`cell_filt`, the
##   per-triple filtered joint) and whose denominator is its group mass
##   (`joint_filt`).  Nothing new is computed -- the filter simply stops
##   throwing the uncollapsed cells away, exactly as F2-A did one level down.
##   Marginalising k gives the pair distribution one period back,
##     Pr[s_{t-2}=i, s_{t-1}=j | y_{1:T}] = sum_k Pr[i, j, k | y_{1:T}],
##   which is the recursion; it starts from Pr[s_{T-1}, s_T | y_{1:T}] =
##   `joint_filt[, , T]` (nothing observes past T), and the smoothed regime
##   marginal is Pr[s_t=k | y_T] = sum_j Pr[s_{t-1}=j, s_t=k | y_T].
##
##   THE SURVIVING CONDITIONAL-INDEPENDENCE ASSUMPTION is
##     y_{t+1:T}  _||_  s_{t-2}  |  s_{t-1}, s_t, y_{1:t} .                (JA3)
##   Compare (JA) of the GPB(2) joint pass, y_{t+1:T} _||_ s_{t-1} | s_t,
##   y_{1:t}: (JA3) conditions on the SAME data plus one extra regime, and the
##   variable it must render irrelevant is one period FURTHER back.  It is
##   therefore strictly weaker again, and it is exactly the backward twin of
##   the error the gpb3 FILTER still makes -- as (JA) is of the gpb2 filter's.
##   Both passes are exact whenever every regime path is self-contained
##   (P = I), where (J3) is a point mass.
##
##   STATE PASS.  One adjoint per pair, (r_t^{(j,k)}, N_t^{(j,k)}), zero at
##   t = T.  For the triple (i, j, k) at period t, whose DK blocks are the
##   filter's own by-products for that cell,
##     r^{(i,j),k} = a^{(ijk)} + L^{(ijk)'} r_t^{(j,k)}
##     N^{(i,j),k} = M^{(ijk)} + L^{(ijk)'} N_t^{(j,k)} L^{(ijk)}
##   and these are collapsed over k with the ROW-normalised smoothed triple
##   w[(i,j), k] = Pr[i,j,k | y_T] / Pr[i,j | y_T] to give the adjoint of the
##   period-(t-1) component (i, j) -- which is precisely the index of the
##   filtered pair the adjoint is applied to:
##     s_{t-1|T}^{(i,j)} = beta_filt[, m, t-1] + P_filt[, , m, t-1] r^{(i,j)} ,
##     V_{t-1|T}^{(i,j)} = P - P N^{(i,j)} P .
##   Smoothed shocks ride the same adjoint per cell,
##     eps^{(m,k)} = Sigma_e^{(k)} (G^{(mk)} r_t^{(g)} + du^{(mk)}),
##   and are collapsed WITHIN the group g = (j, k) over the merged coordinate i
##   with the COLUMN-normalised smoothed triple.  The adjoint the shock reads
##   is conditioned on the PAIR, not just on s_t, which is the whole point.
##
##   OUTPUT.  The reported per-REGIME moments collapse the h^2 pair components
##   onto h regimes with Pr[s_{t-1}=j | s_t=k, y_T], the covariance carrying
##   the across-pair spread term, so the object has exactly the shape and
##   meaning the gpb2 pass returns.  The prior conditionals stand in wherever a
##   smoothed mass underflows to 0 (and at t = 1), as in the gpb2 pass.
##
## The gpb2 pass below is UNTOUCHED by all of this (separate function, pinned
## bit-identical in test-ms-smoother-gpb3.R): the two recursions differ in what
## they index, not in a parameter, and merging them would put the release path
## at the mercy of an edit meant for the other.
## --------------------------------------------------------------------------
##
## ORACLES (test-ms-smoother.R): with P = I and IDENTICAL regimes this is the
## single-regime DK smoother, so states, covariances AND shocks must reproduce
## kalman_smoother() to ~1e-12; and with P = I but DIFFERENT regimes each regime
## path is a self-contained Kalman system, so the transition identity
## s_{t|T}^{(j)} = T s_{t-1|T}^{(j)} + R eps_{t|T}^{(j)} must hold exactly.
## --------------------------------------------------------------------------


#' Kim (1994) smoother for Markov-switching DSGE
#'
#' Backward smoother for the Kim-Nelson filter: returns smoothed regime
#' probabilities \eqn{\Pr[s_t = j \mid y_{1:T}]}, smoothed states
#' \eqn{E[s_t \mid y_{1:T}]} with their covariances, and smoothed structural
#' shocks \eqn{E[\varepsilon_t \mid y_{1:T}]}.  The forward pass is
#' \code{\link{ms_kim_filter}} itself (called with
#' \code{return_state_path = TRUE}), so the smoother and the likelihood can
#' never disagree about the filtered moments or the innovations.
#'
#' The state pass is the Durbin-Koopman \eqn{(r_t, N_t)} adjoint recursion,
#' applied per regime path and collapsed by Kim's weights.  A
#' Rauch-Tung-Striebel pass would NOT be exact here: under this package's lag-1
#' observation timing \eqn{y_{t+1} = Z s_t + D \varepsilon_{t+1}} loads
#' \eqn{s_t} directly, so \eqn{y_{t+1}} is not conditionally independent of
#' \eqn{s_t} given \eqn{s_{t+1}} (see \code{\link{kalman_smoother}}).
#'
#' @param data  Observation matrix (\code{n_obs x T}).  May contain \code{NA}s.
#' @param dr  Decision rule (output of \code{\link{solve_perturbation}}).
#' @param model  Compiled model object.
#' @param params  Named numeric vector of parameter values.
#' @param obs_vars  Character vector of observed variable names.
#' @param ms_spec  An \code{\link{ms_dsge_spec}} object.
#' @param me_variance  Scalar variance of TRUE i.i.d. measurement error
#'   (default \code{0}), forwarded to \code{\link{ms_kim_filter}}: it enters
#'   both the innovation covariance and the Joseph covariance update
#'   (\eqn{K (\code{me\_variance} I) K'}), so the smoothed moments are the
#'   exact conditional moments of the model WITH observation noise.  See the
#'   \emph{Measurement error} section of \code{\link{ms_kim_filter}}; note
#'   this is NOT the \eqn{F}-only regulariser convention of the multivariate
#'   \code{\link{kalman_filter}} methods.
#' @param lik_init  State covariance initialisation, forwarded to
#'   \code{\link{ms_kim_filter}}: \code{"auto"}, \code{"stationary"} or
#'   \code{"kappa"}.
#' @param regime_pass  Backward pass for the regime probabilities:
#'   \code{"joint"} (default) or \code{"kim"}.  See the
#'   \emph{Regime pass} section.
#'
#' @param collapse  Character; GPB collapse depth of the FORWARD pass, and
#'   with it the backward recursion.  \code{"gpb2"} (the default) is Kim's:
#'   the backward core carries one Durbin-Koopman adjoint and one filtered
#'   pair per REGIME and collapses them with the \eqn{h \times h} smoothed
#'   joint.  \code{"gpb3"} runs the forward pass at
#'   \code{collapse = "gpb3"} --- \eqn{h^2} components indexed by the PAIR
#'   \eqn{(s_{t-1}, s_t)} --- and smooths it with the matching PAIR-INDEXED
#'   backward pass: one adjoint per pair, collapsed with the
#'   \eqn{h \times h \times h} smoothed joint
#'   \eqn{\Pr[s_{t-2}, s_{t-1}, s_t \mid y_{1:T}]}, whose conditional comes
#'   from the filter's own uncollapsed per-cell posterior.  The returned
#'   object has the same shape either way (per-regime moments), plus the
#'   pair-indexed internals \code{smoothed_pair_probs},
#'   \code{smoothed_states_by_pair} and \code{smoothed_shocks_by_pair} under
#'   \code{"gpb3"}.  Use \code{"gpb3"} when the \code{collapse_max}
#'   diagnostic of a \code{"gpb2"} run is O(1) nats or more (see the
#'   \emph{Collapse quality} section of \code{\link{ms_kim_filter}}); it
#'   costs about \eqn{h} times as much.  \code{regime_pass = "kim"} is not
#'   defined for it (its backward step conditions on \eqn{y_{1:t}}, which
#'   discards exactly what the extra components are kept for) and is an
#'   error.
#'
#' @section Regime pass:
#' Both passes build the smoothed joint as
#' \eqn{\Pr[s_t = i, s_{t+1} = j \mid y_{1:T}] \approx
#'   \Pr[s_{t+1} = j \mid y_{1:T}] \Pr[s_t = i \mid s_{t+1} = j, \cdot]} and
#' differ only in the information set of the conditional.  Kim (1994, eq. 10)
#' conditions on \eqn{y_{1:t}}, which requires
#' \eqn{y_{t+1:T} \perp s_t \mid s_{t+1}, y_{1:t}}.  Under this package's lag-1
#' observation timing \eqn{y_{t+1} = d_i + Z_i s_t + D_j \varepsilon_{t+1}}
#' loads the regime \eqn{i = s_t} DIRECTLY, so that assumption throws away a
#' first-order-informative signal: a brute-force enumeration oracle (all regime
#' paths, exact joint-Gaussian projection per path) puts the Kim pass up to
#' 0.31 away from the exact smoothed regime probability under genuine
#' switching.  \code{regime_pass = "joint"} instead conditions on
#' \eqn{y_{1:t+1}}, reusing the Hamilton filter's own posterior joint
#' \eqn{\Pr[s_t = i, s_{t+1} = j \mid y_{1:t+1}]} (nothing extra is computed),
#' and assumes only
#' \eqn{y_{t+2:T} \perp s_t \mid s_{t+1}, y_{1:t+1}} -- strictly weaker,
#' because it conditions on a strictly larger information set that includes
#' exactly the informative observation.  The state and shock passes reuse the
#' same joint weights.  With \eqn{P = I} the two passes are identical; the
#' forward filter is untouched, so \code{loglik} does not depend on this
#' argument.  \code{"kim"} is retained for comparison and for reproducing
#' published Kim-smoother output.
#'
#' @return A list of class \code{"ms_kim_smoother"} with:
#'   \describe{
#'     \item{\code{loglik}}{Total log-likelihood from the forward pass.}
#'     \item{\code{filtered_probs}}{\code{h x T} matrix
#'       \eqn{\Pr[s_t = j \mid y_{1:t}]}.}
#'     \item{\code{smoothed_probs}}{\code{h x T} matrix
#'       \eqn{\Pr[s_t = j \mid y_{1:T}]}; columns sum to 1.}
#'     \item{\code{smoothed_states}}{\code{n_state x T} matrix
#'       \eqn{E[s_t \mid y_{1:T}]} (regime-averaged).}
#'     \item{\code{smoothed_states_by_regime}}{\code{n_state x h x T} array
#'       \eqn{E[s_t \mid s_t = j, y_{1:T}]}.}
#'     \item{\code{smoothed_cov}}{\code{n_state x n_state x T} array of
#'       regime-averaged smoothed state covariances (including the
#'       across-regime spread term).}
#'     \item{\code{smoothed_cov_by_regime}}{\code{n_state x n_state x h x T}.}
#'     \item{\code{smoothed_shocks}}{\code{n_shk x T} matrix
#'       \eqn{E[\varepsilon_t \mid y_{1:T}]} (regime-averaged).}
#'     \item{\code{smoothed_shocks_by_regime}}{\code{n_shk x h x T} array,
#'       conditional on \eqn{s_t = j}.}
#'     \item{\code{smoothed_initial}, \code{smoothed_initial_by_regime}}{
#'       \eqn{E[s_0 \mid y_{1:T}]}, free from the same backward recursion;
#'       needed to close the transition identity at \eqn{t = 1}.}
#'     \item{\code{filtered_states}}{\code{n_state x T} matrix
#'       \eqn{E[s_t \mid y_{1:t}]} (regime-averaged).}
#'     \item{\code{state_names}, \code{shock_names}, \code{regime_names}}{
#'       Labels.}
#'     \item{\code{regime_pass}, \code{collapse}}{The backward pass and the
#'       GPB depth actually used.}
#'     \item{\code{smoothed_pair_probs}, \code{smoothed_states_by_pair},
#'       \code{smoothed_shocks_by_pair}, \code{smoothed_initial_probs}}{
#'       \code{collapse = "gpb3"} only: the \eqn{h^2} PAIR-indexed
#'       components the GPB(3) backward pass actually carries, before their
#'       collapse onto regimes.  \code{smoothed_pair_probs} is
#'       \code{h^2 x T} with \eqn{\Pr[s_{t-1} = j, s_t = k \mid y_{1:T}]}
#'       at row \eqn{j + (k-1)h}; \code{smoothed_initial_probs} is
#'       \eqn{\Pr[s_0 \mid y_{1:T}]}.}
#'   }
#' @references Kim, C.-J. (1994). Dynamic linear models with Markov-switching.
#'   \emph{Journal of Econometrics} 60(1-2), 1-22.
#'
#'   Durbin, J. and Koopman, S. J. (2012). \emph{Time Series Analysis by State
#'   Space Methods}, 2nd ed. Oxford University Press, sections 4.4-4.5, 6.2.
#' @seealso \code{\link{ms_kim_filter}}, \code{\link{kalman_smoother}}
#' @export
ms_kim_smoother <- function(data, dr, model, params, obs_vars, ms_spec,
                             me_variance = 0,
                             lik_init = c("auto", "stationary", "kappa"),
                             regime_pass = c("joint", "kim"),
                             collapse = c("gpb2", "gpb3")) {

  lik_init    <- match.arg(lik_init)
  regime_pass <- match.arg(regime_pass)
  collapse    <- match.arg(collapse)

  ## GPB(3): a DIFFERENT backward recursion, not a re-indexing (F4-D).
  ## The gpb2 core carries one Durbin-Koopman adjoint per REGIME and collapses
  ## with the h x h smoothed joint; the gpb3 forward pass returns h^2
  ## PAIR-indexed components, so its backward twin carries an adjoint per PAIR
  ## and collapses with the h x h x h joint Pr[s_{t-1}, s_t, s_{t+1} | y_1:T].
  ## Running the gpb2 backward pass on a gpb3 forward pass would report
  ## smoothed moments that do not belong to the likelihood the caller asked
  ## for, so the two are separate functions and the choice is made HERE.
  gpb3 <- identical(collapse, "gpb3")
  if (gpb3 && identical(regime_pass, "kim"))
    stop("ms_kim_smoother: regime_pass = \"kim\" is not defined for ",
         "collapse = \"gpb3\". Kim (1994, eq. 10) conditions the backward ",
         "step on y_{1:t}, which discards exactly the observation the GPB(3) ",
         "components are kept to exploit; the GPB(3) pass is the triple-joint ",
         "one. Use regime_pass = \"joint\" (the default) with ",
         "collapse = \"gpb3\".", call. = FALSE)

  if (!inherits(ms_spec, "ms_dsge_spec"))
    stop("ms_kim_smoother: ms_spec must be an ms_dsge_spec object.",
         call. = FALSE)

  ## ---- forward pass = the filter, verbatim --------------------------------
  fwd <- ms_kim_filter(data, dr, model, params, obs_vars, ms_spec,
                       me_variance       = me_variance,
                       return_state_path = TRUE,
                       lik_init          = lik_init,
                       collapse          = collapse)

  if (!is.finite(fwd$loglik))
    stop("ms_kim_smoother: the forward filter returned a non-finite ",
         "log-likelihood (", fwd$loglik, "); there is nothing to smooth. ",
         "Check the model solution, the data scale, and lik_init.",
         call. = FALSE)

  if (gpb3)
    .ms_smoother_core_gpb3(fwd, ms_spec$transition, ms_spec$pi0,
                           ms_spec$regime_names)
  else
    .ms_smoother_core(fwd, ms_spec$transition, ms_spec$pi0,
                      ms_spec$regime_names, regime_pass = regime_pass)
}


## Internal: is a covariance materially indefinite?
##
## Returns 0 when the matrix is PSD up to round-off, and otherwise the
## RELATIVE minimum eigenvalue min(ev) / max|ev| (a negative number).  The
## cheap Cholesky screen with a 1e-10 relative jitter comes first so the
## O(n^3) eigen-decomposition only runs on the handful of periods that are
## actually suspect: a -1e-18 relative eigenvalue is float noise, not a
## covariance that a caller could be misled by.
## @noRd
.ms_cov_min_eig_rel <- function(V) {
  s <- max(abs(V))
  if (!is.finite(s) || s <= 0) return(0)
  n  <- nrow(V)
  Vs <- (V + t(V)) * 0.5
  ok <- tryCatch({ chol(Vs + diag(1e-10 * s, n)); TRUE },
                 error = function(e) FALSE)
  if (ok) return(0)
  ev <- eigen(Vs, symmetric = TRUE, only.values = TRUE)$values
  mx <- max(abs(ev))
  if (!is.finite(mx) || mx <= 0) return(0)
  min(ev) / mx
}


## Internal: scan the smoothed covariances and WARN if any is indefinite.
##
## Shared by both backward passes so the message is written once; `what`/`how`
## name the collapse that produced the covariance (the GPB(2) h^2 -> h one, or
## the GPB(3) h^3 -> h^2 one), and the wording is otherwise identical to the
## pre-F4-D warning.  See the call site in .ms_smoother_core() for the measured
## magnitudes and why this is reported rather than swallowed.
##
## @param cov_avg  n_state x n_state x T regime-averaged smoothed covariances.
## @param cov_by_j n_state x n_state x h x T per-regime smoothed covariances.
## @param what     Noun phrase naming the collapse (used in the message).
## @param how      Noun phrase naming what was collapsed onto what.
## @return Invisibly, the worst relative min-eigenvalue seen (0 if all PSD).
## @noRd
.ms_cov_psd_scan <- function(cov_avg, cov_by_j, what, how) {
  n_T <- dim(cov_avg)[3L]
  h   <- dim(cov_by_j)[3L]
  bad_n <- 0L; bad_rel <- 0
  for (t in seq_len(n_T)) {
    rel <- .ms_cov_min_eig_rel(cov_avg[, , t])
    if (rel < 0) { bad_n <- bad_n + 1L; bad_rel <- min(bad_rel, rel) }
  }
  bad_j <- 0L
  for (j in seq_len(h)) for (t in seq_len(n_T)) {
    rel <- .ms_cov_min_eig_rel(cov_by_j[, , j, t])
    if (rel < 0) { bad_j <- bad_j + 1L; bad_rel <- min(bad_rel, rel) }
  }
  if (bad_n > 0L || bad_j > 0L)
    warning(sprintf(paste0(
      "ms_kim_smoother: ", what, " produced a non-PSD smoothed ",
      "covariance in %d of %d periods (regime-averaged) and %d of %d ",
      "regime-period cells; worst relative min-eigenvalue %.2e. This is the ",
      "known approximation error of ", how, " ",
      "when the regimes have DIFFERENT state dynamics; the smoothed means, ",
      "shocks and regime probabilities are unaffected, but the smoothed ",
      "covariances should not be read as uncertainty bands for those periods."),
      bad_n, n_T, bad_j, h * n_T, bad_rel), call. = FALSE)
  invisible(bad_rel)
}


## Internal: the shared Kim (1994) + Durbin-Koopman BACKWARD pass.
##
## Both smoothers -- reduced-form (ms_kim_smoother, one decision rule with
## switching shock scales) and structural (ms_kim_smoother_struct, per-regime
## decision rules) -- run EXACTLY this recursion.  Everything regime-specific
## about the state space (TT_j, RR_j, ZZ_i, DD_i, Sigma_e^{(j)}) has already
## been absorbed by the forward pass into the per-path DK blocks
## (a, M, L, G, du), so the backward pass never touches a state-space matrix
## and there is only ONE copy of it to get right.
##
## @param fwd  Output of ms_kim_filter(return_state_path = TRUE) or
##   ms_kim_filter_struct(return_state_path = TRUE).  Consumed fields:
##   loglik, n_T, regime_probs, beta_filt, P_filt, dk_path, Sigma_e_list,
##   P0_list, state_names, shock_names.
## @param P   h x h transition matrix.
## @param pi0 Length-h initial regime distribution.
## @param regime_names Length-h character vector.
## @param regime_pass "joint" (default) or "kim"; see the file header.
## @return The `ms_kim_smoother` object (callers may append fields).
## @noRd
.ms_smoother_core <- function(fwd, P, pi0, regime_names,
                              regime_pass = c("joint", "kim")) {

  regime_pass <- match.arg(regime_pass)

  h    <- nrow(P)
  n_T  <- fwd$n_T
  ## Dimensions come from the stored moments, not from TT/RR: in the
  ## structural path those are per-regime LISTS.
  n_state <- dim(fwd$beta_filt)[1L]
  n_shk   <- nrow(fwd$Sigma_e_list[[1L]])

  filt_prob <- fwd$regime_probs        # h x T,  Pr[s_t = j | y_{1:t}]
  beta_f    <- fwd$beta_filt           # n_state x h x T
  P_f       <- fwd$P_filt              # n_state x n_state x h x T
  dk        <- fwd$dk_path             # [[t]][[i + (j-1)*h]] = list(a,M,L,G,du)
  Se_l      <- fwd$Sigma_e_list        # per-regime shock covariance

  ## ---- pass 1: smoothed regime probabilities ------------------------------
  ## Both passes have the same shape,
  ##   Pr[s_t=j, s_{t+1}=k | y_T] = Pr[s_{t+1}=k | y_T] * c_t[j, k]
  ## and differ ONLY in the conditional c_t[j, k] = Pr[s_t=j | s_{t+1}=k, .]:
  ##   "kim"   conditions on y_{1:t}    -> Pr[s_t=j|y_t] P[j,k] / Pr[s_{t+1}=k|y_t]
  ##   "joint" conditions on y_{1:t+1}  -> filter joint at t+1, column-normalised
  ## See the derivation in the file header: (JA) is strictly weaker than (KA),
  ## and under this package's lag-1 timing y_{t+1} loads s_t directly, so the
  ## extra observation is exactly the informative one.
  ## The state pass below needs the SAME joint weights, so they are computed
  ## once, stored, and reused rather than recomputed inside the state loop.
  use_joint <- identical(regime_pass, "joint")
  if (use_joint && is.null(fwd$joint_filt))
    stop("ms_kim_smoother: regime_pass = \"joint\" needs the filter's ",
         "posterior joint regime probabilities (`joint_filt`), which this ",
         "forward pass did not return. Re-run the filter with ",
         "return_state_path = TRUE, or use regime_pass = \"kim\".",
         call. = FALSE)

  sm_prob <- matrix(0, h, n_T)
  joint_l <- vector("list", n_T)       # joint_l[[t]] = Pr[s_t=., s_{t+1}=. |y_T]
  sm_prob[, n_T] <- filt_prob[, n_T]   # terminal condition

  if (n_T >= 2L) for (t in (n_T - 1L):1L) {
    if (use_joint) {
      ## num[j, k] = Pr[s_t=j, s_{t+1}=k | y_{1:t+1}] (the filter's own joint);
      ## den[k]    = Pr[s_{t+1}=k | y_{1:t+1}] = its column sum.
      num <- matrix(fwd$joint_filt[, , t + 1L], h, h)
      den <- colSums(num)
    } else {
      ## num[j, k] = Pr[s_t=j | y_{1:t}] P[j, k];  den[k] = Pr[s_{t+1}=k|y_{1:t}]
      num <- filt_prob[, t] * P
      den <- as.numeric(crossprod(P, filt_prob[, t]))
    }
    joint <- matrix(0, h, h)
    for (k in seq_len(h)) {
      if (den[k] <= 0 || sm_prob[k, t + 1L] <= 0) next
      joint[, k] <- sm_prob[k, t + 1L] * num[, k] / den[k]
    }
    p_sm_t <- rowSums(joint)
    ## Renormalise against accumulated round-off (the sum is 1 analytically);
    ## this keeps the reported probabilities a proper distribution at every t.
    tot <- sum(p_sm_t)
    if (is.finite(tot) && tot > 0) {
      joint  <- joint / tot
      p_sm_t <- p_sm_t / tot
    }
    sm_prob[, t] <- p_sm_t
    joint_l[[t]] <- joint
  }

  ## ---- pass 2: Durbin-Koopman state / shock adjoint ------------------------
  beta_s  <- array(0, c(n_state, h, n_T))
  P_s     <- array(0, c(n_state, n_state, h, n_T))
  eps_s   <- array(0, c(n_shk, h, n_T))
  ## r_cur[, j] = r_t^{(j)} ; N_cur[, , j] = N_t^{(j)}. Terminal: both zero,
  ## which is what makes the smoothed pair at t = T the FILTERED pair (no
  ## observation loads s_T).
  r_cur <- matrix(0, n_state, h)
  N_cur <- array(0, c(n_state, n_state, h))

  beta_s[, , n_T] <- beta_f[, , n_T]
  P_s[, , , n_T]  <- P_f[, , , n_T]

  ## Per-regime smoothed initial state s_{0|T}^{(j)} (filled at t = 1).
  b0_s <- matrix(0, n_state, h)

  ## PRIOR conditional regime weights, used at t = 1 (there is no smoothed
  ## joint at t = 0) and as the fallback whenever a regime's smoothed mass
  ## UNDERFLOWS to exactly 0 -- which routinely happens when one regime fits
  ## far better than another (Pr[s_t = 1 | y_T] reached 0 on the P = I /
  ## scale-3 fixture).  Falling back to zero weights instead would silently
  ## null that regime's adjoint, leaving its smoothed states equal to its
  ## filtered states and breaking the transition identity by ~1e-1.  The prior
  ## conditional is the only defensible answer when the data identify nothing,
  ## and it is EXACT in the degenerate case: with P = I both matrices below are
  ## the identity, so each regime path stays a self-contained Kalman system.
  ##   w_from_prior[i, j] = Pr[s_{t-1} = i | s_t = j]  (pi0-weighted Bayes)
  ##   w_to_prior[i, j]   = Pr[s_t = j | s_{t-1} = i]  = P[i, j]
  w_from_prior <- matrix(0, h, h)
  for (j in seq_len(h)) {
    col_j <- pi0 * P[, j]
    s_j   <- sum(col_j)
    w_from_prior[, j] <- if (s_j > 0) col_j / s_j else rep(1 / h, h)
  }
  w_to_prior <- P

  for (t in n_T:1L) {
    ## -- smoothed shocks at t, on each (from i -> to j) path ---------------
    ## eps^{(i,j)} = Sigma_e^{(j)} ( G^{(i,j)} r_t^{(j)} + du^{(i,j)} ), with
    ## r_t^{(j)} the adjoint BEFORE this period's backward step -- exactly the
    ## timing kalman_smoother() uses for its single regime.
    dk_t <- dk[[t]]
    eps_ij <- array(0, c(n_shk, h, h))
    for (i in seq_len(h)) for (j in seq_len(h)) {
      blk <- dk_t[[i + (j - 1L) * h]]
      eps_ij[, i, j] <- as.numeric(
        Se_l[[j]] %*% (blk$G %*% r_cur[, j] + blk$du))
    }
    ## Collapse over the FROM regime i: conditional on s_t = j the smoothed
    ## shock is the i-mixture with weights Pr[s_{t-1}=i | s_t=j, y_T], i.e. the
    ## COLUMN-normalised smoothed joint at t-1, with the prior conditional
    ## standing in for any column whose mass underflowed (and for t = 1).
    w_from <- if (t >= 2L && !is.null(joint_l[[t - 1L]])) {
      jm <- joint_l[[t - 1L]]                       # Pr[s_{t-1}=i, s_t=j | y_T]
      cs <- colSums(jm)
      out <- w_from_prior
      ok  <- cs > 0
      if (any(ok))
        out[, ok] <- sweep(jm[, ok, drop = FALSE], 2L, cs[ok], "/")
      out
    } else {
      w_from_prior
    }
    for (j in seq_len(h)) {
      acc <- numeric(n_shk)
      for (i in seq_len(h)) acc <- acc + w_from[i, j] * eps_ij[, i, j]
      eps_s[, j, t] <- acc
    }

    ## -- backward step: (r_t, N_t) -> (r_{t-1}, N_{t-1}), per FROM regime ---
    ## Path (i at t-1) -> (j at t) uses THIS period's blocks; the incoming
    ## adjoint is the one for regime j at time t.
    r_new <- matrix(0, n_state, h)
    N_new <- array(0, c(n_state, n_state, h))
    r_ij  <- array(0, c(n_state, h, h))
    N_ij  <- array(0, c(n_state, n_state, h, h))
    for (i in seq_len(h)) for (j in seq_len(h)) {
      blk <- dk_t[[i + (j - 1L) * h]]
      Lt  <- blk$L
      r_ij[, i, j]   <- blk$a + as.numeric(crossprod(Lt, r_cur[, j]))
      N_ij[, , i, j] <- blk$M + crossprod(Lt, N_cur[, , j]) %*% Lt
    }

    ## Weights for collapsing the (i -> j) adjoints back onto regime i at
    ## t-1: Pr[s_t = j | s_{t-1} = i, y_T], i.e. the ROW-normalised smoothed
    ## joint at t-1, with P[i, ] standing in for any row whose mass
    ## underflowed (and for t = 1, where the recursion still yields s_{0|T}).
    w_to <- if (t >= 2L && !is.null(joint_l[[t - 1L]])) {
      jm  <- joint_l[[t - 1L]]
      rs  <- rowSums(jm)
      out <- w_to_prior
      ok  <- rs > 0
      ## M / v recycles v down the ROWS, so this divides row i by rs[i].
      if (any(ok)) out[ok, ] <- jm[ok, , drop = FALSE] / rs[ok]
      out
    } else {
      w_to_prior
    }

    for (i in seq_len(h)) {
      rr <- numeric(n_state); NN <- matrix(0, n_state, n_state)
      for (j in seq_len(h)) {
        w <- w_to[i, j]
        if (w <= 0) next
        rr <- rr + w * r_ij[, i, j]
        NN <- NN + w * N_ij[, , i, j]
      }
      r_new[, i]   <- rr
      N_new[, , i] <- (NN + t(NN)) * 0.5
    }
    r_cur <- r_new
    N_cur <- N_new

    ## -- smoothed state / covariance for period t-1 ------------------------
    ##   s_{t-1|T}^{(i)} = s_{t-1|t-1}^{(i)} + P_{t-1|t-1}^{(i)} r_{t-1}^{(i)}
    ##   V_{t-1|T}^{(i)} = P^{(i)} - P^{(i)} N_{t-1}^{(i)} P^{(i)}
    ## (t = 1 uses the pre-sample pair (0, P_0^{(i)}) and gives s_{0|T}.)
    for (i in seq_len(h)) {
      if (t == 1L) {
        s_in <- numeric(n_state)
        P_in <- fwd$P0_list[[i]]
        b0_s[, i] <- s_in + as.numeric(P_in %*% r_cur[, i])
      } else {
        s_in <- beta_f[, i, t - 1L]
        P_in <- matrix(P_f[, , i, t - 1L], n_state, n_state)
        beta_s[, i, t - 1L] <- s_in + as.numeric(P_in %*% r_cur[, i])
        V_i <- P_in - P_in %*% N_cur[, , i] %*% P_in
        P_s[, , i, t - 1L] <- (V_i + t(V_i)) * 0.5
      }
    }
  }

  ## ---- regime-averaged summaries ------------------------------------------
  smoothed_states <- matrix(0, n_state, n_T)
  filtered_states <- matrix(0, n_state, n_T)
  smoothed_shocks <- matrix(0, n_shk, n_T)
  smoothed_cov    <- array(0, c(n_state, n_state, n_T))

  for (t in seq_len(n_T)) {
    w_s <- sm_prob[, t]
    w_f <- filt_prob[, t]
    bs  <- numeric(n_state); bf <- numeric(n_state); es <- numeric(n_shk)
    for (j in seq_len(h)) {
      bs <- bs + w_s[j] * beta_s[, j, t]
      bf <- bf + w_f[j] * beta_f[, j, t]
      es <- es + w_s[j] * eps_s[, j, t]
    }
    smoothed_states[, t] <- bs
    filtered_states[, t] <- bf
    smoothed_shocks[, t] <- es

    Vc <- matrix(0, n_state, n_state)
    for (j in seq_len(h)) {
      d  <- beta_s[, j, t] - bs
      Vc <- Vc + w_s[j] * (P_s[, , j, t] + tcrossprod(d))
    }
    smoothed_cov[, , t] <- (Vc + t(Vc)) * 0.5
  }

  ## ---- PSD scan on the smoothed covariances (fail LOUD) -------------------
  ## The Durbin-Koopman form V = P - P N P is PSD for an EXACT Kalman system,
  ## and it is exact here whenever each regime path is self-contained (P = I,
  ## or a single regime).  Under genuine switching the Kim/GPB(2) collapse
  ## replaces every path's own covariance by the h collapsed ones, so the
  ## incoming adjoint N_t^{(j)} can carry more information than the collapsed
  ## P^{(i)} it is subtracted from and V can come out INDEFINITE.  Measured
  ## 2026-09-03 on the structural rho_a = 0.95 / 0.10 fixture: min eigenvalue
  ## -1.0e-01 against a max of +2.4e-05 (regime 1), 50 of 300 periods affected
  ## in regime 2.  The reduced-form path (common TT/RR, switching shock scale
  ## only) shows nothing above round-off on the same kind of fixture
  ## (-3.3e-19, relative -2.1e-18) -- it is the switching DYNAMICS that break
  ## the inequality, not the switching variances.  The smoothed MEANS are not
  ## affected (they satisfy the transition identity exactly under P = I and
  ## reproduce the data to ~1e-2 relative under switching).
  ## Silently returning a non-PSD covariance as an uncertainty band would be
  ## the worst outcome, so it is reported.
  .ms_cov_psd_scan(smoothed_cov, P_s,
                   what = "the Kim (GPB(2)) collapse",
                   how  = "collapsing h^2 paths onto h covariances")

  ## s_{0|T}: weight the per-regime pre-sample states by pi0 (the only
  ## distribution over s_0 the model supplies).
  smoothed_initial <- as.numeric(b0_s %*% pi0)

  st_names <- fwd$state_names
  sh_names <- fwd$shock_names
  rownames(smoothed_states) <- st_names
  rownames(filtered_states) <- st_names
  rownames(smoothed_shocks) <- sh_names
  names(smoothed_initial)   <- st_names
  dimnames(smoothed_cov)    <- list(st_names, st_names, NULL)
  rownames(sm_prob)   <- regime_names
  rownames(filt_prob) <- regime_names
  dimnames(beta_s)    <- list(st_names, regime_names, NULL)
  dimnames(eps_s)     <- list(sh_names, regime_names, NULL)
  dimnames(b0_s)      <- list(st_names, regime_names)

  structure(
    list(
      loglik                     = fwd$loglik,
      filtered_probs             = filt_prob,
      smoothed_probs             = sm_prob,
      smoothed_states            = smoothed_states,
      smoothed_states_by_regime  = beta_s,
      smoothed_cov               = smoothed_cov,
      smoothed_cov_by_regime     = P_s,
      smoothed_shocks            = smoothed_shocks,
      smoothed_shocks_by_regime  = eps_s,
      smoothed_initial           = smoothed_initial,
      smoothed_initial_by_regime = b0_s,
      filtered_states            = filtered_states,
      state_names                = st_names,
      shock_names                = sh_names,
      regime_names               = regime_names,
      regime_pass                = regime_pass,
      collapse                   = "gpb2",
      n_regimes                  = h,
      n_T                        = n_T
    ),
    class = c("ms_kim_smoother", "list")
  )
}


## Internal: the GPB(3) PAIR-INDEXED backward pass.
##
## The gpb3 twin of .ms_smoother_core(): same inputs, same output object, but
## the recursion is carried on the h^2 PAIR components (s_{t-1}, s_t) the
## gpb3 forward pass produces, with the h^3 smoothed joint over triples
## (s_{t-2}, s_{t-1}, s_t) as the collapse weights.  The full derivation --
## the (J3) conditional, the surviving assumption (JA3), and the state pass --
## is in the file header.  Both smoothers route here when the forward pass was
## run with collapse = "gpb3"; everything regime-specific about the state space
## is already inside the filter's per-cell DK blocks, so this pass, like the
## gpb2 one, never touches a state-space matrix.
##
## WHY A SECOND FUNCTION and not a parameterised one: the two passes index
## different objects (per-regime vs per-pair adjoints, h x h vs h x h x h
## joints) and the gpb2 pass is on every release path.  A shared body would
## make each edit a risk to the other; the shared parts that ARE identical
## (the min-eigenvalue screen, the PSD scan, the regime-averaged summaries'
## arithmetic) are factored into helpers instead.
##
## @param fwd  Output of ms_kim_filter[_struct](return_state_path = TRUE,
##   collapse = "gpb3").  Consumed fields as in .ms_smoother_core(), plus
##   `cell_filt` (the per-triple filtered joint).
## @param P    h x h transition matrix.
## @param pi0  Length-h initial regime distribution.
## @param regime_names Length-h character vector.
## @return The `ms_kim_smoother` object (callers may append fields).
## @noRd
.ms_smoother_core_gpb3 <- function(fwd, P, pi0, regime_names) {

  h    <- nrow(P)
  hh   <- h * h
  n_T  <- fwd$n_T
  n_state <- dim(fwd$beta_filt)[1L]
  n_shk   <- nrow(fwd$Sigma_e_list[[1L]])

  if (is.null(fwd$cell_filt))
    stop("ms_kim_smoother: the GPB(3) backward pass needs the filter's ",
         "UNCOLLAPSED per-cell posteriors (`cell_filt`), which this forward ",
         "pass did not return. Re-run the filter with ",
         "return_state_path = TRUE and collapse = \"gpb3\".", call. = FALSE)
  if (!identical(dim(fwd$beta_filt)[2L], hh))
    stop("ms_kim_smoother: the GPB(3) backward pass expects h^2 = ", hh,
         " pair-indexed forward components, got ", dim(fwd$beta_filt)[2L],
         ".", call. = FALSE)

  beta_f <- fwd$beta_filt              # n_state x h^2 x T, pair (j, k)
  P_f    <- fwd$P_filt                 # n_state x n_state x h^2 x T
  dk     <- fwd$dk_path                # [[t]][[m + (k-1) M]] = list(a,M,L,G,du)
  Se_l   <- fwd$Sigma_e_list

  ## PRIOR conditionals, used at t = 1 and wherever a smoothed mass underflows
  ## to exactly 0 (see the gpb2 core for why zero weights would be worse):
  ##   w_from_prior[i, j] = Pr[s_{t-1} = i | s_t = j]  (pi0-weighted Bayes)
  ##   P[i, ]             = Pr[s_t = . | s_{t-1} = i]
  w_from_prior <- matrix(0, h, h)
  for (j in seq_len(h)) {
    col_j <- pi0 * P[, j]
    s_j   <- sum(col_j)
    w_from_prior[, j] <- if (s_j > 0) col_j / s_j else rep(1 / h, h)
  }

  ## ---- pass 1: the h^3 smoothed joint, backwards -------------------------
  ## sm_pair[, t] = Pr[s_{t-1} = j, s_t = k | y_{1:T}] at index j + (k-1)h;
  ## sm_cell[[t]] = Pr[cell | y_{1:T}] (M_t x h), the smoothed TRIPLE.
  ## Terminal: nothing observes past T, so the smoothed pair AT T is the
  ## filtered pair (`joint_filt`, the gpb3 group mass).
  sm_pair <- matrix(0, hh, n_T)
  sm_cell <- vector("list", n_T)
  sm_pair[, n_T] <- as.numeric(fwd$joint_filt[, , n_T])
  sm_s0   <- pi0                       # Pr[s_0 = i | y_{1:T}], filled at t = 1

  for (t in n_T:1L) {
    cf     <- fwd$cell_filt[[t]]       # M x h, Pr[cell | y_{1:t}]
    M      <- nrow(cf)
    n_coll <- M %/% h                  # cells merged per surviving pair
    pc     <- as.numeric(cf)
    ## Group (= surviving pair) of each column-major cell, and its filtered
    ## mass: the groups are CONTIGUOUS blocks of n_coll cells, exactly as in
    ## the forward collapse.
    gsum <- colSums(matrix(pc, n_coll, hh))
    smc  <- numeric(M * h)
    for (g in seq_len(hh)) {
      if (gsum[g] <= 0 || sm_pair[g, t] <= 0) next
      cells <- (g - 1L) * n_coll + seq_len(n_coll)
      smc[cells] <- sm_pair[g, t] * pc[cells] / gsum[g]
    }
    ## Renormalise against accumulated round-off (analytically the triple
    ## sums to 1), so every reported distribution is proper.
    tot <- sum(smc)
    if (is.finite(tot) && tot > 0) smc <- smc / tot
    sm_cell[[t]] <- matrix(smc, M, h)
    ## Marginalise s_t out of the triple -> the smoothed pair one period back
    ## (at t = 1 that is the smoothed distribution of s_0).
    prev <- rowSums(sm_cell[[t]])
    if (t >= 2L) sm_pair[, t - 1L] <- prev else sm_s0 <- prev
  }

  ## Smoothed regime marginal Pr[s_t = k | y_T] = sum_j Pr[(j, k) | y_T].
  sm_prob <- matrix(0, h, n_T)
  for (t in seq_len(n_T)) sm_prob[, t] <- colSums(matrix(sm_pair[, t], h, h))

  ## ---- pass 2: Durbin-Koopman adjoint, one per PAIR ----------------------
  beta_sp <- array(0, c(n_state, hh, n_T))          # per-pair smoothed mean
  P_sp    <- array(0, c(n_state, n_state, hh, n_T))
  eps_sp  <- array(0, c(n_shk, hh, n_T))
  r_cur   <- matrix(0, n_state, hh)
  N_cur   <- array(0, c(n_state, n_state, hh))
  beta_sp[, , n_T] <- beta_f[, , n_T]
  P_sp[, , , n_T]  <- P_f[, , , n_T]
  b0_s    <- matrix(0, n_state, h)

  for (t in n_T:1L) {
    dk_t   <- dk[[t]]
    M      <- length(dk_t) %/% h
    n_coll <- M %/% h
    ## FROM regime of each incoming component: the pair's SECOND index at
    ## t >= 2 (component m = i + (j-1)h carries s_{t-1} = j), the component
    ## itself at t = 1 (where components are the h values of s_0).
    from_m <- if (M == h) seq_len(h) else ((seq_len(M) - 1L) %/% h) + 1L
    smc    <- sm_cell[[t]]                          # M x h smoothed triple

    r_mk   <- array(0, c(n_state, M, h))
    N_mk   <- array(0, c(n_state, n_state, M, h))
    eps_mk <- array(0, c(n_shk, M, h))
    for (m in seq_len(M)) for (k in seq_len(h)) {
      cell <- m + (k - 1L) * M
      g    <- (cell - 1L) %/% n_coll + 1L           # post-update pair (j, k)
      blk  <- dk_t[[cell]]
      Lt   <- blk$L
      r_mk[, m, k]    <- blk$a + as.numeric(crossprod(Lt, r_cur[, g]))
      N_mk[, , m, k]  <- blk$M + crossprod(Lt, N_cur[, , g]) %*% Lt
      ## The shock adjoint reads the PAIR-conditioned r, not an s_t-only one.
      eps_mk[, m, k]  <- as.numeric(
        Se_l[[k]] %*% (blk$G %*% r_cur[, g] + blk$du))
    }

    ## -- smoothed shocks per surviving pair g = (j, k) ---------------------
    ## Collapse the n_coll cells of group g over the merged coordinate i with
    ## Pr[i | (j,k), y_T]; the prior conditional Pr[s_{t-2}=i | s_{t-1}=j]
    ## stands in when the group's smoothed mass underflowed (and at t = 1,
    ## where n_coll = 1 and the weight is 1 either way).
    eps_flat <- matrix(as.numeric(eps_mk), n_shk, M * h)
    smc_flat <- as.numeric(smc)
    for (g in seq_len(hh)) {
      cells <- (g - 1L) * n_coll + seq_len(n_coll)
      wv    <- smc_flat[cells]
      sw    <- sum(wv)
      if (!is.finite(sw) || sw <= 0) {
        j  <- (g - 1L) %% h + 1L                    # s_{t-1} of this pair
        wv <- if (n_coll == 1L) 1 else w_from_prior[, j]
      } else {
        wv <- wv / sw
      }
      acc <- numeric(n_shk)
      for (u in seq_len(n_coll)) acc <- acc + wv[u] * eps_flat[, cells[u]]
      eps_sp[, g, t] <- acc
    }

    ## -- backward step: adjoint of the period-(t-1) component m = (i, j) ---
    ## Collapse over the TO regime k with Pr[s_t = k | (i, j), y_T].
    r_new <- matrix(0, n_state, M)
    N_new <- array(0, c(n_state, n_state, M))
    for (m in seq_len(M)) {
      wv <- smc[m, ]
      sw <- sum(wv)
      wv <- if (is.finite(sw) && sw > 0) wv / sw else P[from_m[m], ]
      rr <- numeric(n_state); NN <- matrix(0, n_state, n_state)
      for (k in seq_len(h)) {
        if (wv[k] <= 0) next
        rr <- rr + wv[k] * r_mk[, m, k]
        NN <- NN + wv[k] * N_mk[, , m, k]
      }
      r_new[, m]   <- rr
      N_new[, , m] <- (NN + t(NN)) * 0.5
    }
    r_cur <- r_new
    N_cur <- N_new

    ## -- smoothed pair moments for period t-1 ------------------------------
    for (m in seq_len(M)) {
      if (t == 1L) {
        b0_s[, m] <- as.numeric(fwd$P0_list[[m]] %*% r_cur[, m])
      } else {
        s_in <- beta_f[, m, t - 1L]
        P_in <- matrix(P_f[, , m, t - 1L], n_state, n_state)
        beta_sp[, m, t - 1L] <- s_in + as.numeric(P_in %*% r_cur[, m])
        V_m <- P_in - P_in %*% N_cur[, , m] %*% P_in
        P_sp[, , m, t - 1L]  <- (V_m + t(V_m)) * 0.5
      }
    }
  }

  ## ---- collapse the h^2 pair components onto h regimes -------------------
  ## Conditional on s_t = k the smoothed law is the mixture over s_{t-1} = j
  ## with weights Pr[s_{t-1} = j | s_t = k, y_T]; the covariance carries the
  ## across-pair spread term, exactly as the forward collapse does.
  beta_s <- array(0, c(n_state, h, n_T))
  P_s    <- array(0, c(n_state, n_state, h, n_T))
  eps_s  <- array(0, c(n_shk, h, n_T))
  filt_prob <- fwd$regime_probs
  beta_fj   <- array(0, c(n_state, h, n_T))

  for (t in seq_len(n_T)) {
    pm <- matrix(sm_pair[, t], h, h)                # [j, k]
    fm <- matrix(fwd$joint_filt[, , t], h, h)
    for (k in seq_len(h)) {
      wj <- pm[, k]; sw <- sum(wj)
      wj <- if (is.finite(sw) && sw > 0) wj / sw else w_from_prior[, k]
      bs <- numeric(n_state); es <- numeric(n_shk)
      for (j in seq_len(h)) {
        g  <- j + (k - 1L) * h
        bs <- bs + wj[j] * beta_sp[, g, t]
        es <- es + wj[j] * eps_sp[, g, t]
      }
      Vk <- matrix(0, n_state, n_state)
      for (j in seq_len(h)) {
        g <- j + (k - 1L) * h
        d <- beta_sp[, g, t] - bs
        Vk <- Vk + wj[j] * (P_sp[, , g, t] + tcrossprod(d))
      }
      beta_s[, k, t] <- bs
      eps_s[, k, t]  <- es
      P_s[, , k, t]  <- (Vk + t(Vk)) * 0.5

      ## FILTERED per-regime mean, on the filter's own pair weights.
      wf <- fm[, k]; sf <- sum(wf)
      wf <- if (is.finite(sf) && sf > 0) wf / sf else w_from_prior[, k]
      bf <- numeric(n_state)
      for (j in seq_len(h)) bf <- bf + wf[j] * beta_f[, j + (k - 1L) * h, t]
      beta_fj[, k, t] <- bf
    }
  }

  ## ---- regime-averaged summaries -----------------------------------------
  smoothed_states <- matrix(0, n_state, n_T)
  filtered_states <- matrix(0, n_state, n_T)
  smoothed_shocks <- matrix(0, n_shk, n_T)
  smoothed_cov    <- array(0, c(n_state, n_state, n_T))

  for (t in seq_len(n_T)) {
    w_s <- sm_prob[, t]
    w_f <- filt_prob[, t]
    bs  <- numeric(n_state); bf <- numeric(n_state); es <- numeric(n_shk)
    for (j in seq_len(h)) {
      bs <- bs + w_s[j] * beta_s[, j, t]
      bf <- bf + w_f[j] * beta_fj[, j, t]
      es <- es + w_s[j] * eps_s[, j, t]
    }
    smoothed_states[, t] <- bs
    filtered_states[, t] <- bf
    smoothed_shocks[, t] <- es

    Vc <- matrix(0, n_state, n_state)
    for (j in seq_len(h)) {
      d  <- beta_s[, j, t] - bs
      Vc <- Vc + w_s[j] * (P_s[, , j, t] + tcrossprod(d))
    }
    smoothed_cov[, , t] <- (Vc + t(Vc)) * 0.5
  }

  ## The DK inequality V = P - P N P >= 0 can still fail under a collapse --
  ## GPB(3) merges h^3 cells onto h^2 pairs instead of h^2 onto h, which is a
  ## smaller error, not no error -- so the same scan runs here.
  .ms_cov_psd_scan(smoothed_cov, P_s,
                   what = "the GPB(3) pair collapse",
                   how  = "collapsing h^3 paths onto h^2 covariances")

  ## s_{0|T}: the gpb3 pass HAS a smoothed distribution over s_0 (the t = 1
  ## triple, marginalised), so it is used rather than the prior pi0 the gpb2
  ## pass has to fall back on.  With identical regimes, or P = I and a
  ## symmetric fit, the two coincide.
  smoothed_initial <- as.numeric(b0_s %*% sm_s0)

  st_names <- fwd$state_names
  sh_names <- fwd$shock_names
  rownames(smoothed_states) <- st_names
  rownames(filtered_states) <- st_names
  rownames(smoothed_shocks) <- sh_names
  names(smoothed_initial)   <- st_names
  dimnames(smoothed_cov)    <- list(st_names, st_names, NULL)
  rownames(sm_prob)   <- regime_names
  rownames(filt_prob) <- regime_names
  dimnames(beta_s)    <- list(st_names, regime_names, NULL)
  dimnames(eps_s)     <- list(sh_names, regime_names, NULL)
  dimnames(b0_s)      <- list(st_names, regime_names)

  structure(
    list(
      loglik                     = fwd$loglik,
      filtered_probs             = filt_prob,
      smoothed_probs             = sm_prob,
      smoothed_states            = smoothed_states,
      smoothed_states_by_regime  = beta_s,
      smoothed_cov               = smoothed_cov,
      smoothed_cov_by_regime     = P_s,
      smoothed_shocks            = smoothed_shocks,
      smoothed_shocks_by_regime  = eps_s,
      smoothed_initial           = smoothed_initial,
      smoothed_initial_by_regime = b0_s,
      filtered_states            = filtered_states,
      state_names                = st_names,
      shock_names                = sh_names,
      regime_names               = regime_names,
      regime_pass                = "joint",
      collapse                   = "gpb3",
      ## Pair-indexed internals, for callers that want the components the
      ## GPB(3) pass actually carries rather than their regime collapse.
      smoothed_pair_probs        = sm_pair,
      smoothed_states_by_pair    = beta_sp,
      smoothed_shocks_by_pair    = eps_sp,
      smoothed_initial_probs     = sm_s0,
      n_regimes                  = h,
      n_T                        = n_T
    ),
    class = c("ms_kim_smoother", "list")
  )
}


## --------------------------------------------------------------------------
## STRUCTURAL PATH
## --------------------------------------------------------------------------
##
## WHICH STEADY STATE ARE THE SMOOTHED STATES DEVIATIONS FROM? (2026-09-03)
##
## In the structural path every regime has its OWN steady state ys_s, so
## d_s = ys_s[obs_vars] differs across regimes and the question is real: Kim's
## collapse averages the h paths arriving in regime j, and averaging deviations
## taken around DIFFERENT means is a well-known way to manufacture a bias
## (the mixture mean would silently absorb the difference of the means).
##
## DECISION: a COMMON REFERENCE.  ms_kim_filter_struct()'s state recursion is
##     s_t = TT_j s_{t-1} + RR_j eps_t
## with NO regime-dependent intercept, while the per-regime steady state enters
## ONLY the observation equation, as the switching intercept d_i in
##     y_t = ZZ_i s_{t-1} + d_i + DD_i eps_t .
## So the filtered state vector is, by construction, one single coordinate
## system shared by all regimes -- there is no per-regime state mean in it at
## all, and Kim's collapse (in the filter) and the collapses in the backward
## pass are therefore mixtures of quantities measured from the SAME origin.
## That is what makes them legitimate; the trap is avoided by not having
## per-regime state means, not by correcting for them.
##
## The regime-specific means are then reported where they actually belong:
##   * observation level: y_hat_t = sum_j Pr[s_t=j|y_T] (ZZ_j s_{t-1|T}^{(j)}
##     + d_j + DD_j eps_{t|T}^{(j)}) -- computed by ms_smoothed_fit_struct()
##     below and pinned by the data-reproduction oracle;
##   * state LEVELS: `smoothed_states_level` mixes LEVELS, not deviations,
##       sum_j Pr[s_t=j|y_T] * (s_{t|T}^{(j)} + ys_j[state_vars]) ,
##     which is the only regime average of levels that is unbiased when the
##     ys_j differ.  `ss_states_by_regime` ships the ys_j[state_vars] columns
##     so a caller can redo the mapping against any other reference.
##
## KNOWN LIMITATION (inherited from the filter, NOT introduced here).
## solve_ms_perturbation() also returns the FRWZ (2016, Prop. 2) constant drift
## c_s = ghx_s sum_{s'} P[s,s'](ys_{s'} - ys_s), which the structural FILTER
## does not use; the smoother deliberately reproduces the filter's state space
## exactly rather than smoothing a different model from the one whose
## likelihood it decomposes.  Under P = I (every oracle below) c_s = 0
## identically, so the drift is not what any of the gates are absorbing.


#' Kim (1994) smoother for structural Markov-switching DSGE
#'
#' Backward smoother for \code{\link{ms_kim_filter_struct}}, i.e. for MS-DSGE
#' models solved by \code{\link{solve_ms_perturbation}}, where each regime has
#' its OWN decision rule (\code{TT_s, RR_s, ZZ_s, DD_s}) and steady state.
#' Returns smoothed regime probabilities, smoothed states with covariances and
#' smoothed structural shocks, in the same object class as
#' \code{\link{ms_kim_smoother}} (so \code{print} works on both) plus the
#' \code{ms_dr} it smoothed and the level-space fields described below.
#'
#' The forward pass is \code{ms_kim_filter_struct(return_state_path = TRUE)}
#' verbatim, and the backward pass is the SAME Durbin-Koopman code the
#' reduced-form smoother runs: every regime-specific state-space matrix has
#' already been folded into the per-path \eqn{(a, M, L, G, du)} blocks by the
#' forward pass, so there is exactly one backward recursion in the package.
#'
#' @section Deviation reference:
#' Smoothed states are deviations from the COMMON reference implied by the
#' structural filter's state equation \eqn{s_t = T_j s_{t-1} + R_j
#' \varepsilon_t}, which carries no regime-dependent intercept: the per-regime
#' steady state enters only as the switching observation intercept
#' \eqn{d_i = ys_i[obs\_vars]}.  All regimes' state vectors therefore live in
#' one coordinate system and Kim's collapse mixes like with like.  Regime
#' means are re-introduced only in level space, in
#' \code{smoothed_states_level}, which mixes LEVELS
#' \eqn{s_{t|T}^{(j)} + ys_j[\text{state}]} with the smoothed regime
#' probabilities.
#'
#' @section Smoothed covariances under switching dynamics:
#' The Durbin-Koopman covariance \eqn{V = P - P N P} is PSD for an exact
#' Kalman system, and this smoother is exact whenever each regime path is
#' self-contained (\eqn{P = I}, or one active regime).  Under genuine
#' switching of the state DYNAMICS the Kim / GPB(2) collapse replaces each
#' path's own covariance by the \eqn{h} collapsed ones, and the incoming
#' adjoint \eqn{N_t^{(j)}} can then carry more information than the collapsed
#' \eqn{P^{(i)}} it is subtracted from, so \eqn{V} can come out INDEFINITE
#' (measured on an RBC fixture whose regimes differ in \eqn{\rho_a}: minimum
#' eigenvalue \eqn{-1.0 \times 10^{-1}} against a positive maximum of
#' \eqn{+2.4 \times 10^{-5}}, 16 of 300 periods affected).  The smoothed means,
#' shocks and regime probabilities are unaffected.  Rather than hand back a
#' silent non-PSD uncertainty band, the smoother scans its own output and
#' \code{warning()}s with the count and the worst relative eigenvalue.  The
#' reduced-form path (\code{\link{ms_kim_smoother}}, common transition, only
#' shock variances switching) does not show the effect above round-off.
#'
#' @param data  Observation matrix (\code{n_obs x T}).  May contain \code{NA}s.
#' @param ms_dr  An \code{MsDecisionRules} object from
#'   \code{\link{solve_ms_perturbation}}.
#' @param model  Compiled model object.
#' @param params  Named numeric parameter vector (used only for
#'   \eqn{\Sigma_e} when \code{Sigma_e_by_regime} is \code{NULL}).
#' @param obs_vars  Character vector of observed variable names.
#' @param Sigma_e_by_regime  Optional length-h list of per-regime shock
#'   covariance matrices; \code{NULL} uses one common \eqn{\Sigma_e}.
#' @param me_variance  Scalar variance of TRUE i.i.d. measurement error
#'   (default \code{0}), forwarded to \code{\link{ms_kim_filter_struct}}: it
#'   enters both the innovation covariance and the Joseph covariance update
#'   (\eqn{K (\code{me\_variance} I) K'}).  See the \emph{Measurement error}
#'   section of \code{\link{ms_kim_filter}}; this is NOT the \eqn{F}-only
#'   regulariser convention of the multivariate \code{\link{kalman_filter}}
#'   methods.
#' @param lik_init  State covariance initialisation: \code{"auto"},
#'   \code{"stationary"} or \code{"kappa"}.
#' @param regime_pass  Backward pass for the regime probabilities:
#'   \code{"joint"} (default) or \code{"kim"}.  See the \emph{Regime pass}
#'   section of \code{\link{ms_kim_smoother}}: the default conditions the
#'   backward step on \eqn{y_{1:t+1}} rather than \eqn{y_{1:t}}, which under
#'   this package's lag-1 timing keeps the observation that loads \eqn{s_t}
#'   directly.  The forward filter is untouched, so \code{loglik} does not
#'   depend on it.
#'
#' @param collapse  Character; GPB collapse depth of the FORWARD pass, and
#'   with it the backward recursion.  \code{"gpb2"} (the default) is Kim's:
#'   the backward core carries one Durbin-Koopman adjoint and one filtered
#'   pair per REGIME and collapses them with the \eqn{h \times h} smoothed
#'   joint.  \code{"gpb3"} runs the forward pass at
#'   \code{collapse = "gpb3"} --- \eqn{h^2} components indexed by the PAIR
#'   \eqn{(s_{t-1}, s_t)} --- and smooths it with the matching PAIR-INDEXED
#'   backward pass: one adjoint per pair, collapsed with the
#'   \eqn{h \times h \times h} smoothed joint
#'   \eqn{\Pr[s_{t-2}, s_{t-1}, s_t \mid y_{1:T}]}, whose conditional comes
#'   from the filter's own uncollapsed per-cell posterior.  The returned
#'   object has the same shape either way (per-regime moments), plus the
#'   pair-indexed internals \code{smoothed_pair_probs},
#'   \code{smoothed_states_by_pair} and \code{smoothed_shocks_by_pair} under
#'   \code{"gpb3"}.  Use \code{"gpb3"} when the \code{collapse_max}
#'   diagnostic of a \code{"gpb2"} run is O(1) nats or more (see the
#'   \emph{Collapse quality} section of \code{\link{ms_kim_filter}}); it
#'   costs about \eqn{h} times as much.  \code{regime_pass = "kim"} is not
#'   defined for it (its backward step conditions on \eqn{y_{1:t}}, which
#'   discards exactly what the extra components are kept for) and is an
#'   error.
#'
#' @return A list of class \code{"ms_kim_smoother"} with all the fields
#'   documented in \code{\link{ms_kim_smoother}}, plus
#'   \describe{
#'     \item{\code{ms_dr}}{The \code{MsDecisionRules} object smoothed.}
#'     \item{\code{smoothed_states_level}}{\code{n_state x T} matrix of
#'       regime-averaged smoothed states in LEVELS (see the deviation-reference
#'       section).}
#'     \item{\code{ss_states_by_regime}}{\code{n_state x h} matrix of
#'       per-regime steady states of the state variables.}
#'   }
#' @seealso \code{\link{ms_kim_filter_struct}}, \code{\link{ms_kim_smoother}},
#'   \code{\link{solve_ms_perturbation}}
#' @export
ms_kim_smoother_struct <- function(data, ms_dr, model, params, obs_vars,
                                    Sigma_e_by_regime = NULL,
                                    me_variance = 0,
                                    lik_init = c("auto", "stationary",
                                                 "kappa"),
                                    regime_pass = c("joint", "kim"),
                                    collapse = c("gpb2", "gpb3")) {

  lik_init    <- match.arg(lik_init)
  regime_pass <- match.arg(regime_pass)
  collapse    <- match.arg(collapse)

  ## GPB(3): routed to the PAIR-INDEXED backward pass (F4-D); see the note in
  ## ms_kim_smoother() and the derivation in this file's header.
  gpb3 <- identical(collapse, "gpb3")
  if (gpb3 && identical(regime_pass, "kim"))
    stop("ms_kim_smoother_struct: regime_pass = \"kim\" is not defined for ",
         "collapse = \"gpb3\"; the GPB(3) pass is the triple-joint one. Use ",
         "regime_pass = \"joint\" (the default).", call. = FALSE)

  if (!inherits(ms_dr, "MsDecisionRules"))
    stop("ms_kim_smoother_struct: ms_dr must be an MsDecisionRules object ",
         "(output of solve_ms_perturbation()).", call. = FALSE)

  ## ---- forward pass = the structural filter, verbatim ---------------------
  fwd <- ms_kim_filter_struct(data, ms_dr, model, params, obs_vars,
                              Sigma_e_by_regime = Sigma_e_by_regime,
                              me_variance       = me_variance,
                              return_state_path = TRUE,
                              lik_init          = lik_init,
                              collapse          = collapse)

  if (!is.finite(fwd$loglik))
    stop("ms_kim_smoother_struct: the forward filter returned a non-finite ",
         "log-likelihood (", fwd$loglik, "); there is nothing to smooth. ",
         "Check the MS solution, the data scale, and lik_init.",
         call. = FALSE)

  h  <- length(ms_dr$dr)
  rn <- if (!is.null(names(ms_dr$dr))) names(ms_dr$dr)
        else paste0("regime", seq_len(h))

  out <- if (gpb3) .ms_smoother_core_gpb3(fwd, ms_dr$P, ms_dr$pi0, rn)
         else .ms_smoother_core(fwd, ms_dr$P, ms_dr$pi0, rn,
                                regime_pass = regime_pass)

  ## ---- level-space companion (see the deviation-reference note above) -----
  st_names <- fwd$state_names
  ss_by_j  <- vapply(ms_dr$dr, function(d) as.numeric(d$ys[st_names]),
                     numeric(length(st_names)))
  ss_by_j  <- matrix(ss_by_j, nrow = length(st_names), ncol = h,
                     dimnames = list(st_names, rn))

  lev <- matrix(0, length(st_names), fwd$n_T,
                dimnames = list(st_names, NULL))
  for (t in seq_len(fwd$n_T)) {
    w <- out$smoothed_probs[, t]
    for (j in seq_len(h))
      lev[, t] <- lev[, t] + w[j] * (out$smoothed_states_by_regime[, j, t] +
                                       ss_by_j[, j])
  }

  out$ms_dr                <- ms_dr
  out$smoothed_states_level <- lev
  out$ss_states_by_regime  <- ss_by_j
  out
}


#' Smoothed observation fit for a structural MS smoother
#'
#' Rebuilds the observables from a \code{\link{ms_kim_smoother_struct}} result
#' through the regime-weighted observation equation
#' \deqn{\hat y_t = \sum_i \Pr[s_{t-1} = i \mid y_{1:T}]
#'   (Z_i s_{t-1|T}^{(i)} + d_i)
#'   + \sum_j \Pr[s_t = j \mid y_{1:T}] D_j \varepsilon_{t|T}^{(j)},}
#' i.e. the check that the smoothed states and shocks actually reproduce the
#' data they were smoothed from.
#'
#' The regime index follows the structural filter's own timing: under dynhr's
#' lag-1 convention \eqn{y_t} loads \eqn{s_{t-1}} through the FROM-regime
#' matrices \eqn{(Z_i, d_i)}, while the contemporaneous shock loads the
#' TO-regime \eqn{D_j}; each expectation is therefore taken under the regime
#' distribution that actually governs that term (\eqn{\Pr[s_0 = i]} is
#' \code{pi0} at \eqn{t = 1}).  With \eqn{P = I} the two distributions
#' coincide and the identity is exact to round-off; with genuine switching a
#' residual remains, because the Kim filter/smoother is itself an
#' approximation (the \eqn{h^2 \to h} collapse keeps only \eqn{h} conditional
#' moments per period, so the state and the shock terms are averaged under
#' regime marginals that no longer share a single exact joint).
#'
#' A \code{collapse = "gpb3"} smoother uses the same formula (its pair index
#' at \eqn{t-1} is \eqn{(s_{t-2}, s_{t-1})}, so it is NOT the conditioning
#' pair of this equation, and summing the pair joint over the coordinate each
#' term does not condition on returns exactly this regime-marginal form), with
#' one refinement: at \eqn{t = 1} the FROM distribution is the smoothed
#' \eqn{\Pr[s_0 \mid y_{1:T}]} the GPB(3) pass supplies rather than the
#' prior \code{pi0}.
#'
#' @param sm  An object returned by \code{\link{ms_kim_smoother_struct}}.
#' @param data  The same observation matrix (\code{n_obs x T}) that was
#'   smoothed.
#' @return A list with \code{fitted} (\code{n_obs x T}), \code{residual}
#'   (\code{fitted - data}, \code{NA} columns kept as \code{NA}) and
#'   \code{max_abs_residual} over the observed entries.
#' @seealso \code{\link{ms_kim_smoother_struct}}
#' @export
ms_smoothed_fit_struct <- function(sm, data) {

  if (!inherits(sm, "ms_kim_smoother") || is.null(sm$ms_dr))
    stop("ms_smoothed_fit_struct: `sm` must come from ",
         "ms_kim_smoother_struct().", call. = FALSE)

  ms_dr   <- sm$ms_dr
  h       <- sm$n_regimes
  n_T     <- sm$n_T
  dr1     <- ms_dr$dr[[1L]]
  endo    <- dr1$endo_names
  obs_var <- rownames(data)
  if (is.null(obs_var))
    stop("ms_smoothed_fit_struct: `data` must have observable names as ",
         "rownames.", call. = FALSE)
  obs_idx <- match(obs_var, endo)
  if (anyNA(obs_idx))
    stop("ms_smoothed_fit_struct: observables not found in the model: ",
         paste(obs_var[is.na(obs_idx)], collapse = ", "), call. = FALSE)

  st_idx <- dr1$state_idx
  ZZ <- lapply(ms_dr$dr, function(d) d$ghx[obs_idx, , drop = FALSE])
  DD <- lapply(ms_dr$dr, function(d) d$ghu[obs_idx, , drop = FALSE])
  dd <- lapply(ms_dr$dr, function(d) as.numeric(d$ys[obs_var]))

  fitted <- matrix(0, length(obs_var), n_T,
                   dimnames = list(obs_var, colnames(data)))
  pi0 <- ms_dr$pi0

  ## THE PAIR COMPONENTS DO NOT SHORTCUT THIS IDENTITY (F4-D).  A GPB(3)
  ## smoother carries E[s_{t-1} | s_{t-2}, s_{t-1}, y_T] -- the pair one
  ## period back is (s_{t-2}, s_{t-1}), NOT the (s_{t-1}, s_t) this equation
  ## conditions on -- so pairing the period-(t-1) state component with the
  ## period-t shock component would mismatch the state's regime by one lag.
  ## (Measured: doing that lifts the relative residual on the T = 300
  ## switching fixture from 5.9e-3 to 1.6e-1.)  Summing the correct pair joint
  ## over the coordinate each term does not condition on returns exactly the
  ## regime-marginal form below, so gpb3 uses the same code; only the t = 1
  ## FROM distribution improves, because the pair pass HAS a smoothed
  ## Pr[s_0 | y_{1:T}] where the gpb2 pass can only offer the prior pi0.
  pri0 <- if (!is.null(sm$smoothed_initial_probs)) sm$smoothed_initial_probs
          else pi0

  for (t in seq_len(n_T)) {
    ## FROM-regime distribution (governs Z_i, d_i and s_{t-1}); at t = 1 the
    ## only distribution over s_0 the model supplies is pi0.
    w_from <- if (t == 1L) pri0 else sm$smoothed_probs[, t - 1L]
    ## TO-regime distribution (governs D_j and eps_t).
    w_to   <- sm$smoothed_probs[, t]
    acc <- numeric(length(obs_var))
    for (j in seq_len(h)) {
      prev <- if (t == 1L) sm$smoothed_initial_by_regime[, j]
              else         sm$smoothed_states_by_regime[, j, t - 1L]
      acc <- acc + w_from[j] * (as.numeric(ZZ[[j]] %*% prev) + dd[[j]]) +
        w_to[j] * as.numeric(DD[[j]] %*% sm$smoothed_shocks_by_regime[, j, t])
    }
    fitted[, t] <- acc
  }

  resid <- fitted - data
  list(fitted           = fitted,
       residual         = resid,
       max_abs_residual = max(abs(resid[is.finite(resid)])))
}


#' @export
#' @noRd
print.ms_kim_smoother <- function(x, ...) {
  cat(sprintf("<ms_kim_smoother>  %d regimes, %d states, T = %d%s\n",
              x$n_regimes, length(x$state_names), x$n_T,
              if (identical(x$collapse, "gpb3")) "  [GPB(3), pair-indexed]"
              else ""))
  cat(sprintf("  loglik: %.6f\n", x$loglik))
  cat("  Mean smoothed regime probabilities:\n")
  for (j in seq_len(x$n_regimes))
    cat(sprintf("    %s: %.4f\n", x$regime_names[j], mean(x$smoothed_probs[j, ])))
  invisible(x)
}
