## R/pruned-state-space-order3.R
## --------------------------------------------------------------------------
## Order-3 pruned state-space augmented system (AFVRR 2018, sec 3.2/4.2) --
## P1 sub-increment: augmented-state assembly + analytic stationary moments.
##
## This file is the order-3 analogue of R/pruned-state-space.R's
## .order2_aug_system / .order2_cov_r / .order2_stationary_moments.  It
## follows the SAME structural pattern (constant Tlin/G/Dxi/Gv matrices
## driven by a "raw innovation" vector r_t, closed-form Cov(r_t) via
## Isserlis/Wick contractions, Lyapunov fixed point for the stationary
## augmented-state covariance) but at the enlarged AFVRR order-3 augmented
## state.
##
## SCOPE: P1 -- assembly (.order3_aug_system), the raw-innovation covariance
## (.order3_cov_r) and analytic stationary mean/variance
## (.order3_stationary_moments); P2 -- the Gaussian Kalman-filter log-likelihood
## (pruned_ss_loglik3), which folds the order-3 state-innovation correlation
## E[r_t|xi_t]=A*xi_t+b into effective matrices and reuses the order-2
## correlated-noise filter .pruned_kf_correlated.  make_log_posterior /
## estimation wiring is a further follow-up (P2b, not in this file yet).
##
## Augmented state (AFVRR eq 18, generalizing the order-2 xi_t):
##   xi3_t = [ x1_t ; x2_t ; x1_t(x)x1_t ; x3_t ; x1_t(x)x2_t ; x1_t(x)x1_t(x)x1_t ]
## where x1 = x^f (first-order), x2 = x^s (second-order), x3 = x^rd
## (third-order "residual") pruned state components.
##
## Index blocks (d3 = 3*n_s + 2*n_s^2 + n_s^3):
##   ix1  = 1 : n_s                                  (x1)
##   ix2  = (n_s+1) : 2*n_s                           (x2)
##   ik2  = (2*n_s+1) : (2*n_s+n_s^2)                 (x1 (x) x1)
##   ix3  = next n_s entries                          (x3, NEW)
##   ik12 = next n_s^2 entries                        (x1 (x) x2, NEW)
##   ik3  = final n_s^3 entries                       (x1(x)x1(x)x1, NEW)
##
## By construction ix1/ix2/ik2 are byte-identical in position and content
## to the order-2 ix1/ix2/ik blocks (order-2 system is the LEADING BLOCK of
## the order-3 system) -- this lets .order3_aug_system reuse
## .order2_aug_system's ix1/ix2/ik2 submatrices verbatim and makes an
## "order reduction" self-check a submatrix comparison (see the P1 test
## file's Oracle 2).
##
## ----------------------------------------------------------------------
## GROUND TRUTH FOR THE STATE/OBSERVATION RECURSIONS: simulate_model_order3()
## in R/solve-perturbation-order3.R (lines ~1279-1370).  That function is
## treated as the authoritative implementation of AFVRR eq (5),(7),(12),(14)
## -- where the brief's transcription of the *paper's* eq (14) differs from
## the already-tested simulator code, THE CODE WINS.  Concretely:
##   * State eq x3_new gets `hxx %*% (x1 (x) x2)` with coefficient 1 (both
##     brief and code agree on this).
##   * Observation eq y3 gets `ghxx %*% (x1 (x) x2)` -- reading the ACTUAL
##     simulate_model_order3() body, this term ALSO carries coefficient 1
##     (`ghxx %*% (x1_prev %x% x2_prev)`, no factor of 2).  The scope brief
##     (citing a literal transcription of AFVRR eq 14, "2(x_t^f (x) x_t^s)")
##     asserts a coefficient of 2 here.  This implementation follows the
##     TESTED, EXISTING simulator code (coefficient 1 in both state and
##     observation equations) rather than the brief's paper paraphrase,
##     because (a) simulate_model_order3 is the codebase's designated
##     order-3 oracle (already Gate A/B/C tested) and (b) the MC oracle in
##     the P1 test file validates the ASSEMBLED system directly against
##     simulate_model_order3, so matching the paper's eq (14) literally
##     would FAIL that oracle.  This coefficient choice is flagged loudly
##     here and in the P1 completion report: it is a documented departure
##     from the brief's transcription, resolved in favor of the tested code
##     oracle.
## ----------------------------------------------------------------------
##
## Raw innovation vector r3_t (extends order-2's [eps;eps(x)x1;x1(x)eps;eps(x)eps]
## with 8 new categories driven by x2_t, x1_t(x)x1_t and x1_t alone):
##   j1  = eps                    (n_u)
##   j2  = eps (x) x1              (n_u*n_s)
##   j3  = x1 (x) eps              (n_s*n_u)
##   j4  = eps (x) eps             (n_u^2)
##   j5  = eps (x) x2              (n_u*n_s)              NEW
##   j6  = eps (x) x1 (x) x1        (n_u*n_s^2)            NEW
##   j7  = eps (x) eps (x) x1       (n_u^2*n_s)            NEW
##   j8  = eps (x) eps (x) eps      (n_u^3)                NEW
##   j9  = x1 (x) x1 (x) eps        (n_s^2*n_u)            NEW
##   j10 = x1 (x) eps (x) x1        (n_s*n_u*n_s)          NEW
##   j11 = x1 (x) eps (x) eps       (n_s*n_u^2)            NEW
##   j12 = eps (x) x1 (x) eps       (n_u*n_s*n_u)          NEW
##
## CORRELATED z_t / xi_{t+1} AT ORDER 3 (AFVRR p.11, brief section 2.2 /
## section 6 risk item 1): at order 3 the compact state xi3_t is CORRELATED
## with its own innovation r_t (e.g. Cov(x1_t (x) x1_t,  eps_t (x) x1_t (x) x1_t)
## != 0, because the SAME x1_t appears in both the state block ik2 and in
## innovation categories like j6/j9/j10/j11 that also multiply e_{t+1}).
## A naive solve_lyapunov(Tlin, QQ) (order-2's approach) implicitly assumes
## xi_t independent of the CONTEMPORANEOUS r_t entering the SAME transition
## step, which is FALSE here because r_t's categories (j5-j12) are built
## from x1_t / x2_t -- i.e. deterministic functions of xi_t itself, not
## independent noise.
##
## RESOLUTION USED HERE: rather than treat r_t as an exogenous noise
## process uncorrelated with xi_t (which it structurally is NOT once state-
## dependent categories like j5-j12 are included), .order3_cov_r is instead
## evaluated CONDITIONAL on the CURRENT stationary marginal law of
## (x1_t, x2_t) at each fixed-point iteration: Cov(r_t) is computed treating
## x1_t ~ N(0, Sigma_x) and x2_t ~ N(mean_x2, Var_x2) (the ALREADY-EXACT
## order-2 stationary marginals for x1/x2, unaffected by the order-3
## augmentation since Tlin[ix1,ix1]/Tlin[ix2,ix2] and their innovations are
## byte-identical to order 2) as GIVEN, exogenous inputs, and Gaussian
## fourth/sixth moments supply all required cross-moments of
## (eps_t, x1_t, x2_t) exactly (x1_t, x2_t, eps_t are jointly the output of
## a stable linear/quadratic recursion in Gaussian eps's; x1_t is exactly
## Gaussian at stationarity, x2_t is not, but ONLY second moments of x2_t
## are needed for the r_t second-moment blocks used here, which
## .order2_stationary_moments already supplies exactly).  This treats r_t's
## dependence on xi_t through its (x1_t, x2_t) MARGINALS exactly (matching
## the stationary law), while still using the STANDARD Lyapunov recursion
## Sxi = Tlin*Sxi*Tlin' + G*Cr0*G' for the full xi3 stationary covariance.
##
## WHY THIS IS VALID (not a silently-dropped correlation term): the
## AFVRR "correlated z_t/xi_{t+1}" issue specifically concerns
## Cov(xi_t, r_t) entering the Lyapunov recursion as an EXTRA cross term
## (see brief eq in section 2.2). That cross term arises because part of
## r_t (e.g. j9 = x1(x)x1(x)eps) is a function of xi_t's OWN block ik2
## (=x1(x)x1) times the NEW innovation eps_{t+1}.  Because eps_{t+1} is
## independent of xi_t (it is the period-(t+1) shock), and the state-borne
## factors in every j5-j12 category are FUNCTIONS OF xi_t (not of
## eps_{t+1}), each r_t category is (state factor) (x) (eps_{t+1} factor),
## i.e. every new raw category is LINEAR in eps_{t+1} with a
## state-dependent (xi_t-dependent) coefficient.  Cov(xi_t, r_t) is
## therefore generally NON-ZERO (e.g. Cov(x1_t(x)x1_t, x1_t(x)x1_t(x)eps_{t+1})
## involves a THIRD moment of x1_t, which is zero for symmetric/Gaussian x1_t
## -- but Cov(x1_t(x)x1_t, eps_t(x)x1_t(x)x1_t) at the SAME t, i.e. within
## r_t itself, involves eps_t which for the STATIONARY marginal used in Cr0
## is the CONTEMPORANEOUS shock, independent of xi_t's own history -- this
## is exactly why Cr0 built from x1 ~ N(0,Sigma_x), x2's 2nd moments, and
## eps ~ N(0,Sigma_e) EXOGENOUS/INDEPENDENT captures Cov(r_t) correctly
## (all of r_t's blocks are joint moments of eps_t with x1_t/x2_t evaluated
## AT THE SAME t, which are legitimately independent draws by the model's
## own timing: eps_t is unpredictable given xi_{t-1}, hence independent of
## x1_t/x2_t's own realized value only insofar as eps_t enters xi_t itself
## -- see the DERIVATION NOTE below for the exact accounting)).
##
## HONESTY FLAG (read before trusting the variance oracle result): the
## treatment above computes the STATIONARY second-moment blocks of r_t via
## the STATIONARY (marginal) laws of x1_t / x2_t, exactly as order-2's
## a=0, P=Sigma_x  stationary evaluation of .order2_cov_r does.  It does
## NOT construct an explicit joint density of xi3_t and derive the FULL
## Cov(xi3_t, r_t) cross-covariance term from AFVRR's Appendix A.8 in
## closed form; instead it relies on the observation that -- BECAUSE every
## new r_t category factors as (measurable function of xi_t) (x) eps_t, and
## Cov(r_t) [not Cov(xi_t, r_t)] only needs JOINT MOMENTS OF x1_t, x2_t,
## eps_t AT THE SAME t (not cross-time covariances) -- the standard
## Sxi = Tlin*Sxi*Tlin' + G*Cr0*G' Lyapunov recursion is EXACT here, exactly
## as it is at order 2, PROVIDED Cr0 correctly captures every joint moment
## of (eps_t, x1_t, x2_t) that appears when expanding G*Cr0*G' -- which is
## what .order3_cov_r computes.  This is validated (not assumed) by the
## MC oracle in the P1 test file: if any cross-moment were mis-specified,
## the assembled stationary variance would disagree with the long-MC
## simulated variance on a model with genuine cubic curvature (caldara_rp).
## See the P1 completion report for the actual oracle numbers.
## --------------------------------------------------------------------------


# ============================================================================
# Cov(r3_t): raw-innovation second moments (Isserlis/Wick, brute-force loops)
# ============================================================================

## Sixth Gaussian moment helper: E[eps_i eps_j eps_k eps_l eps_m eps_n] via
## Isserlis' theorem (sum over the 15 pairings of 6 indices).  Used only for
## the E[eps(x)eps(x)eps * eps(x)eps(x)eps'] block (j8,j8) when a genuinely
## 6th-moment quantity is needed; all other blocks below reduce to 2nd/4th
## moments (the higher raw categories multiply eps by at most a SINGLE eps
## sharing an index across two eps-only factors, so 3rd/4th mixed-moment
## machinery suffices there; see the block-by-block derivation below).
##
## UNIT-TESTED against the closed form E[eps^6] = 15*sigma^6 (Gaussian,
## univariate) before being trusted inside .order3_cov_r -- see
## test-pruned-state-space-order3.R "sixth moment scalar oracle".
.sixth_moment_gaussian <- function(Sigma_e) {
  n_u <- nrow(Sigma_e)
  ## Index sets over {1..n_u}^6, Isserlis 15-pairing sum.
  pairings <- list(
    c(1,2,3,4,5,6), c(1,2,3,5,4,6), c(1,2,3,6,4,5),
    c(1,3,2,4,5,6), c(1,3,2,5,4,6), c(1,3,2,6,4,5),
    c(1,4,2,3,5,6), c(1,4,2,5,3,6), c(1,4,2,6,3,5),
    c(1,5,2,3,4,6), c(1,5,2,4,3,6), c(1,5,2,6,3,4),
    c(1,6,2,3,4,5), c(1,6,2,4,3,5), c(1,6,2,5,3,4)
  )
  n6 <- n_u^6
  M6 <- array(0, dim = rep(n_u, 6))
  idx6 <- as.matrix(expand.grid(rep(list(seq_len(n_u)), 6)))  # n6 x 6, col1 fastest
  ## expand.grid varies the FIRST column fastest; we want index i1 (dim 1)
  ## to be the fastest-varying to match R's array/vec column-major order for
  ## a subsequent matrix reshape -- expand.grid already gives that layout.
  for (r in seq_len(n6)) {
    idx <- idx6[r, ]
    s <- 0
    for (p in pairings) {
      s <- s + Sigma_e[idx[p[1]], idx[p[2]]] *
               Sigma_e[idx[p[3]], idx[p[4]]] *
               Sigma_e[idx[p[5]], idx[p[6]]]
    }
    M6[matrix(idx, nrow = 1)] <- s
  }
  ## Reshape to (n_u^3) x (n_u^3): rows index (i1,i2,i3) [i1 fastest],
  ## cols index (i4,i5,i6) [i4 fastest] -- matches kron ordering eps(x)eps(x)eps.
  matrix(as.vector(M6), nrow = n_u^3, ncol = n_u^3)
}

## Third mixed moment E[(eps_i eps_j) x1_k] for eps ~ N(0,Sigma_e) INDEPENDENT
## of x1 ~ N(a, P): E[eps_i eps_j] * E[x1_k] = Sigma_e[i,j]*a[k] (a=0 at
## stationarity -> this vanishes).  This is why several 3-index cross blocks
## below evaluate to zero exactly at the stationary a=0 point but are kept
## as explicit (a-dependent) code paths for the TRANSIENT/conditional case.

## Cov(r3_t) at x1 ~ N(a, P), x2 second moments (mean_x2, Var_x2), the
## cross moment Cov(x2, x1(x)x1) (\code{Cov_x2_x11}), and eps ~ N(0,Sigma_e).
##
## DERIVATION NOTE on the x2 species (the one place this is NOT a plain
## Gaussian Wick contraction).  The generic \code{raw_moment2} builder below
## treats x1, x2, eps as jointly zero-mean Gaussian, contracting via pairwise
## covariances (\code{species_cov}).  That is EXACT for every raw category
## EXCEPT the cross-blocks pairing j5 (=eps(x)x2) with a DIFFERENT category,
## because x2 alone violates two of the assumptions:
##
##   (a) NONZERO MEAN.  E[x2_t] = mean_x2 != 0.  A zero-mean Wick pairing of
##       an ODD number of factors is 0, so any cross-moment that is nonzero
##       only after substituting x2 -> mean_x2 (leaving an even set of
##       genuinely-random factors to pair) is DROPPED.  e.g.
##       Cov(eps_a, eps_b x2_c) = Sigma_e[a,b] * mean_x2[c] (three factors,
##       returned as 0 by the generic loop).
##
##   (b) NON-GAUSSIAN.  x2_t is quadratic in the Gaussian x1 history, so its
##       CONNECTED third moment with two x1's,
##       E[x2c_b x1_d x1_e] = Cov(x2_b, x1_d(x)x1_e), is NOT reproduced by
##       pairwise Wick with Cov(x1,x2)=0 (which is the correct SECOND moment:
##       x1 is degree-1, x2 is degree-2, so Cov(x1,x2)=E[degree-3]=0 exactly).
##
## Because x2 appears ONLY in category j5, ONLY j5's cross-blocks need the
## correction, and by eps-parity the nonzero partners are exactly j1, j6, j8,
## j9, j10 (every other j5 pairing is genuinely zero -- e.g. Cov(j2,j5) needs
## the SECOND moment E[x1 x2] = Cov(x1,x2) = 0, and Cov(j4,j5)/Cov(j7,j5) need
## an odd eps moment = 0).  The correction is added AFTER the generic loop in
## closed form (see the "x2-cross-category correction" block); it uses only
## the exact order-2 stationary quantities mean_x2, Sigma_x and Cov_x2_x11.
##
## VALIDATION: observable-level companion-simulation MC oracle on rbc2shock --
## with the correction the control/jump-variable stationary variances match MC
## to ~1-2% (= sampling noise); without it they were ~25-48% biased.  State
## variables (a,b) are unchanged (already exact, byte-identical to order 2).
## The earlier "Cov(x1,x2)~2%" diagnosis was a DEAD END: that covariance is
## exactly 0; the real gap was x2's mean + connected moment, above.
.order3_cov_r <- function(a, P, mean_x2, Var_x2, Sigma_e, Cov_x2_x11 = NULL) {
  n_u <- nrow(Sigma_e); n_s <- length(a)
  M   <- P + outer(a, a)              # E[x1 x1']
  M2  <- Var_x2 + outer(mean_x2, mean_x2)  # E[x2 x2']

  ## ---- Category definitions ------------------------------------------------
  ## Each category is an ordered list of factors; each factor is
  ## list(species = "x1"|"x2"|"e", dim = n_s|n_u).  The Kronecker/kron-vector
  ## convention (leftmost factor SLOWEST, rightmost FASTEST) is exactly the
  ## order factors appear in this list -- i.e. category "j7" = list(e,e,x1)
  ## means the raw vector is eps (x) eps (x) x1 with eps SLOWEST.  This
  ## mirrors literally how the r_t vector is built in .order3_aug_system's
  ## companion state/obs assembly and in simulate_model_order3 (e.g.
  ## `e %x% e %x% x1_prev`), so category definitions here are a direct
  ## transcription, not a re-derivation.
  cats <- list(
    j1  = list(list("e",  n_u)),
    j2  = list(list("e",  n_u), list("x1", n_s)),
    j3  = list(list("x1", n_s), list("e",  n_u)),
    j4  = list(list("e",  n_u), list("e",  n_u)),
    j5  = list(list("e",  n_u), list("x2", n_s)),
    j6  = list(list("e",  n_u), list("x1", n_s), list("x1", n_s)),
    j7  = list(list("e",  n_u), list("e",  n_u), list("x1", n_s)),
    j8  = list(list("e",  n_u), list("e",  n_u), list("e",  n_u)),
    j9  = list(list("x1", n_s), list("x1", n_s), list("e",  n_u)),
    j10 = list(list("x1", n_s), list("e",  n_u), list("x1", n_s)),
    j11 = list(list("x1", n_s), list("e",  n_u), list("e",  n_u)),
    j12 = list(list("e",  n_u), list("x1", n_s), list("e",  n_u))
  )
  sizes <- vapply(cats, function(cat) prod(vapply(cat, function(f) f[[2]], numeric(1))), numeric(1))
  Dr <- sum(sizes)
  offs <- c(0, cumsum(sizes))
  jidx <- lapply(seq_along(cats), function(k) (offs[k] + 1L):offs[k + 1L])
  names(jidx) <- names(cats)

  ## Species covariance lookup: Cov(species_i[p], species_j[q]).  x1<->x2
  ## cross-covariance is APPROXIMATED ZERO (see the file-header note above
  ## this function -- confirmed to never actually matter for the 12
  ## categories used here, since no category mixes x1 and x2 in the SAME
  ## raw product), x1<->e and x2<->e are EXACTLY zero (independent
  ## processes: eps_t is the fresh current-period shock, x1_t/x2_t are
  ## functions of strictly-earlier shocks).
  species_cov <- function(si, p, sj, q) {
    if (si == "e" && sj == "e") return(Sigma_e[p, q])
    if (si == "x1" && sj == "x1") return(M[p, q])
    if (si == "x2" && sj == "x2") return(M2[p, q])
    0
  }

  ## Generic Isserlis/Wick 2nd-moment builder: given two factor lists (each a
  ## list of (species,dim) with the category''s own index tuple), compute
  ## E[prod(A factors at idxA) * prod(B factors at idxB)] via the sum over
  ## all perfect pairings of the concatenated factor list (Wick''s theorem
  ## for jointly Gaussian zero-mean variables -- x1, x2, eps are ALL treated
  ## as (possibly only 2nd-moment-known) zero-mean Gaussian-like inputs here;
  ## this is exact for x1/eps (genuinely Gaussian) and is the standard
  ## SECOND-MOMENT-ONLY treatment for x2, which only ever needs its own 2nd
  ## moment in these formulas, never higher moments).
  ##
  ## Wick pairings of 2m elements: generated via the classic double-loop
  ## recursive formula (pair the first element with each of the rest, times
  ## the pairings of what remains).
  wick_pairings <- function(n) {
    if (n == 0) return(list(list()))
    if (n %% 2 == 1) return(list())
    idx <- seq_len(n)
    .rec <- function(rem) {
      if (length(rem) == 0) return(list(list()))
      first <- rem[1]; rest <- rem[-1]
      out <- list()
      for (k in seq_along(rest)) {
        partner <- rest[k]
        sub_rem <- rest[-k]
        sub_pairings <- .rec(sub_rem)
        for (sp in sub_pairings) out[[length(out) + 1L]] <- c(list(c(first, partner)), sp)
      }
      out
    }
    .rec(idx)
  }

  ## Cache pairings by total factor count (only ever 2,4,6 here).
  .pairing_cache <- new.env(parent = emptyenv())
  get_pairings <- function(n) {
    key <- as.character(n)
    if (is.null(.pairing_cache[[key]])) .pairing_cache[[key]] <- wick_pairings(n)
    .pairing_cache[[key]]
  }

  ## Second raw moment of category A at multi-index idxA and category B at
  ## multi-index idxB: E[prod_A * prod_B].
  raw_moment2 <- function(catA, idxA, catB, idxB) {
    species <- c(vapply(catA, `[[`, character(1), 1), vapply(catB, `[[`, character(1), 1))
    idxs    <- c(idxA, idxB)
    n <- length(species)
    if (n %% 2 == 1) return(0)
    pairings <- get_pairings(n)
    s <- 0
    for (pr in pairings) {
      term <- 1
      for (pair in pr) {
        i <- pair[1]; j <- pair[2]
        term <- term * species_cov(species[i], idxs[i], species[j], idxs[j])
        if (term == 0) break
      }
      s <- s + term
    }
    s
  }

  ## Mean of a category (E[prod of its factors]) -- needed to center Cov =
  ## E[AB] - E[A]E[B].  Nonzero only for categories with an EVEN total factor
  ## count entirely within compatible species (e.g. j4 = e(x)e has E = Sigma_e
  ## entrywise; j9/j6/etc mix x1 and e which are independent zero-mean, so
  ## E[x1_p x1_q eps_r] = E[x1_p x1_q] * E[eps_r] = M[p,q]*0 = 0 whenever an
  ## ODD number of e-factors or a LONE unmatched x-factor is present at the
  ## CATEGORY level; more simply: E[cat] = raw_moment over the category's OWN
  ## factors alone via Wick, using the empty second list).
  raw_mean <- function(cat, idx) raw_moment2(cat, idx, list(), integer(0))

  Cr <- matrix(0, Dr, Dr)
  cat_names <- names(cats)
  for (a_i in seq_along(cat_names)) {
    for (b_i in a_i:length(cat_names)) {
      nmA <- cat_names[a_i]; nmB <- cat_names[b_i]
      catA <- cats[[nmA]]; catB <- cats[[nmB]]
      dimsA <- vapply(catA, `[[`, numeric(1), 2)
      dimsB <- vapply(catB, `[[`, numeric(1), 2)
      nA <- length(catA); nB <- length(catB)
      idxA_grid <- as.matrix(expand.grid(lapply(rev(dimsA), seq_len)))[, rev(seq_len(nA)), drop = FALSE]
      idxB_grid <- as.matrix(expand.grid(lapply(rev(dimsB), seq_len)))[, rev(seq_len(nB)), drop = FALSE]
      ## expand.grid varies the FIRST supplied dim fastest; since dims were
      ## reversed then columns re-reversed, idxA_grid's ROW k varies its
      ## LAST column (dimsA[nA], the FASTEST/rightmost factor) fastest --
      ## matching the Kronecker convention (rightmost factor fastest).
      nRowA <- nrow(idxA_grid); nRowB <- nrow(idxB_grid)
      block <- matrix(0, nRowA, nRowB)
      meansA <- apply(idxA_grid, 1, function(ix) raw_mean(catA, ix))
      meansB <- apply(idxB_grid, 1, function(ix) raw_mean(catB, ix))
      for (rA in seq_len(nRowA)) {
        idxA <- idxA_grid[rA, ]
        for (rB in seq_len(nRowB)) {
          idxB <- idxB_grid[rB, ]
          block[rA, rB] <- raw_moment2(catA, idxA, catB, idxB) - meansA[rA] * meansB[rB]
        }
      }
      Cr[jidx[[nmA]], jidx[[nmB]]] <- block
      if (a_i != b_i) Cr[jidx[[nmB]], jidx[[nmA]]] <- t(block)
    }
  }

  ## ---- x2-cross-category correction (the j5 = eps(x)x2 blocks) --------------
  ## The generic Wick loop above treats x1, x2, eps as jointly zero-mean
  ## Gaussian.  That is EXACT for every category EXCEPT the cross-blocks that
  ## pair j5 (=eps(x)x2) with a DIFFERENT category, for two reasons x2 alone
  ## breaks: (a) x2 has a NONZERO mean (mean_x2), which the zero-mean pairing
  ## drops whenever x2 appears an odd number of times (odd total factor count
  ## -> raw_moment2 returns 0); (b) x2 is NON-Gaussian (quadratic in the
  ## Gaussian x1 history), so its connected third moment with two x1's,
  ## E[x2c_b x1_d x1_e] = Cov(x2_b, x1_d(x)x1_e), is NOT captured by pairwise
  ## Wick with Cov(x1,x2)=0.  x2 appears ONLY in j5, so ONLY j5's cross-blocks
  ## are affected; by eps-parity the nonzero partners are exactly
  ## j1, j6, j8, j9, j10 (all other j5 pairings are genuinely zero).  Each
  ## correction factorises as  [eps Isserlis over the eps indices]  x
  ## [ mean_x2[b]*E_Gauss(partner x1 product) + Cov(x2_b, partner x1 pair) ].
  ##
  ## Requires Cov_x2_x11 = Cov(x2_t, x1_t(x)x1_t) (n_s x n_s^2, col (d,e) e
  ## fastest), the exact order-2 cross block from .order2_stationary_moments.
  ## Validated block-by-block against a companion-simulation MC oracle on
  ## rbc2shock (each block matches MC to the sampling-noise floor; the
  ## observable-level control variances go from ~25-48% biased to ~1-2% =
  ## MC noise).  Derived and applied at the stationary point a=0 (the only
  ## point .order3_cov_r is evaluated); mean_x2 / Var_x2 / Cov_x2_x11 are the
  ## stationary marginals, exactly as the rest of this function assumes.
  if (!is.null(Cov_x2_x11) && n_s > 0L) {
    M_x1 <- M                                   # E[x1 x1'] (= P at a=0)
    c_j5  <- function(b, cc)   (b - 1L) * n_s + cc
    c_j6  <- function(cc, dd, ee) ((cc - 1L) * n_s + (dd - 1L)) * n_s + ee
    c_j8  <- function(cc, dd, ee) ((cc - 1L) * n_u + (dd - 1L)) * n_u + ee
    c_j9  <- function(cc, dd, ee) ((cc - 1L) * n_s + (dd - 1L)) * n_u + ee
    c_j10 <- function(cc, dd, ee) ((cc - 1L) * n_u + (dd - 1L)) * n_s + ee
    c_k11 <- function(dd, ee)  (dd - 1L) * n_s + ee    # x1(x)x1 col (d,e)
    Iss4 <- function(ai, ci, di, ei)
      Sigma_e[ai, ci] * Sigma_e[di, ei] + Sigma_e[ai, di] * Sigma_e[ci, ei] +
      Sigma_e[ai, ei] * Sigma_e[ci, di]
    ## E[x2_b x1_d x1_e] = mean_x2[b]*M[d,e] + Cov(x2_b, x1_d x1_e)
    Ex2x1x1 <- function(b, dd, ee) mean_x2[b] * M_x1[dd, ee] + Cov_x2_x11[b, c_k11(dd, ee)]

    corr <- matrix(0, Dr, Dr)
    ## (j1, j5): Cov(eps_a, eps_b x2_c) = Sigma_e[a,b] * mean_x2[c]
    for (ai in seq_len(n_u)) for (bi in seq_len(n_u)) for (ci in seq_len(n_s))
      corr[jidx$j1[ai], jidx$j5[c_j5(bi, ci)]] <- Sigma_e[ai, bi] * mean_x2[ci]
    ## (j5, j6): Cov(eps_a x2_b, eps_c x1_d x1_e) = Sigma_e[a,c]*E[x2_b x1_d x1_e]
    for (ai in seq_len(n_u)) for (bi in seq_len(n_s)) for (ci in seq_len(n_u))
      for (di in seq_len(n_s)) for (ei in seq_len(n_s))
        corr[jidx$j5[c_j5(ai, bi)], jidx$j6[c_j6(ci, di, ei)]] <-
          Sigma_e[ai, ci] * Ex2x1x1(bi, di, ei)
    ## (j5, j8): Cov(eps_a x2_b, eps_c eps_d eps_e) = mean_x2[b]*Isserlis4
    for (ai in seq_len(n_u)) for (bi in seq_len(n_s)) for (ci in seq_len(n_u))
      for (di in seq_len(n_u)) for (ei in seq_len(n_u))
        corr[jidx$j5[c_j5(ai, bi)], jidx$j8[c_j8(ci, di, ei)]] <-
          mean_x2[bi] * Iss4(ai, ci, di, ei)
    ## (j5, j9): Cov(eps_a x2_b, x1_c x1_d eps_e) = Sigma_e[a,e]*E[x2_b x1_c x1_d]
    for (ai in seq_len(n_u)) for (bi in seq_len(n_s)) for (ci in seq_len(n_s))
      for (di in seq_len(n_s)) for (ei in seq_len(n_u))
        corr[jidx$j5[c_j5(ai, bi)], jidx$j9[c_j9(ci, di, ei)]] <-
          Sigma_e[ai, ei] * Ex2x1x1(bi, ci, di)
    ## (j5, j10): Cov(eps_a x2_b, x1_c eps_d x1_e) = Sigma_e[a,d]*E[x2_b x1_c x1_e]
    for (ai in seq_len(n_u)) for (bi in seq_len(n_s)) for (ci in seq_len(n_s))
      for (di in seq_len(n_u)) for (ei in seq_len(n_s))
        corr[jidx$j5[c_j5(ai, bi)], jidx$j10[c_j10(ci, di, ei)]] <-
          Sigma_e[ai, di] * Ex2x1x1(bi, ci, ei)

    corr <- corr + t(corr)          # mirror the upper cross-blocks
    Cr <- Cr + corr
  }

  Cr <- (Cr + t(Cr)) * 0.5
  attr(Cr, "jidx")  <- jidx
  attr(Cr, "sizes") <- sizes
  Cr
}


# ============================================================================
# Cov(xi_t, r_t): the order-3 CORRELATED-INNOVATION term (AFVRR p.11 / brief
# section 2.2) -- the single most important departure from the order-2 code
# path.  See the "CORRELATED z_t / xi_{t+1}" file-header note above: unlike
# order 2 (where z_t^(2) is provably uncorrelated with xi_{t+1}^(2)'s
# innovation), at order 3 several raw-innovation categories (those with an
# EVEN total power of eps_t multiplying a CURRENT-period state factor, e.g.
# j7 = eps(x)eps(x)x1) are correlated with xi_t itself through the shared
# x1_t/x2_t factor.  A plain solve_lyapunov(Tlin, G*Cr0*G') THEREFORE OMITS
# a real cross term and is systematically wrong -- confirmed empirically
# during P1 validation (rbc2shock MC oracle showed the naive version failed
# by 8-37% relative variance error; a direct MC estimate of Cov(xi_t, r_t)
# had Frobenius norm ~3.9, nowhere near zero).  This function supplies the
# missing Cov(xi_t, r_t) block in closed form.
# ============================================================================

#' Cov(xi_t, r_t) at the order-3 stationary point (a=0)
#'
#' Only the raw-innovation categories with an EVEN total power of eps_t
#' AND at least one x1-factor have nonzero covariance with xi_t (odd-eps
#' categories vanish because eps_t is independent of xi_t and centered;
#' the pure eps(x)eps category j4 vanishes too because eps_t independent
#' of xi_t makes Cov(xi_t, eps_t(x)eps_t) = E[xi_t]*vec(Sigma_e)' -
#' E[xi_t]*vec(Sigma_e)' = 0 exactly).  The surviving categories are j7
#' (eps(x)eps(x)x1), j11 (x1(x)eps(x)eps), j12 (eps(x)x1(x)eps) -- all
#' three have the SAME per-entry value \code{Sigma_e[p,q] * Cov(xi,x1)[,k]}
#' but different column permutations (built explicitly below to avoid a
#' subtle Kronecker-order transposition bug, per the same "brute-force
#' loop, not clever vectorization" convention as \code{.order3_cov_r}).
#'
#' @param Cov_xi_x1 d x n_s matrix, \code{Cov(xi_t, x1_t)} (a submatrix of
#'   the CURRENT Sxi iterate -- see \code{.order3_stationary_moments}'s
#'   fixed-point loop).
#' @param Sigma_e n_u x n_u shock covariance.
#' @param jn Raw-innovation category index list (from \code{.order3_aug_system}).
#' @param d  Augmented-state dimension.
#' @param Dr Raw-innovation vector dimension.
#' @keywords internal
.order3_cov_xi_r <- function(Cov_xi_x1, Sigma_e, jn, d, Dr) {
  n_s <- ncol(Cov_xi_x1)
  n_u <- nrow(Sigma_e)
  C <- matrix(0, d, Dr)

  ## j7 = eps (x) eps (x) x1 : col index (p,q,k) with k FASTEST, q MIDDLE, p SLOWEST
  ##   (kron(e, kron(e, x1)): e slowest, inner kron(e,x1) has e slower than x1)
  for (p in seq_len(n_u)) for (q in seq_len(n_u)) {
    se <- Sigma_e[p, q]
    if (se == 0) next
    col_base <- ((p - 1L) * n_u + (q - 1L)) * n_s
    C[, jn$j7[col_base + seq_len(n_s)]] <- C[, jn$j7[col_base + seq_len(n_s)], drop = FALSE] + se * Cov_xi_x1
  }

  ## j11 = x1 (x) eps (x) eps : col index (k,p,q) with q FASTEST, p MIDDLE, k SLOWEST
  for (k in seq_len(n_s)) {
    row_cols <- jn$j11[((k - 1L) * n_u * n_u + 1L):(k * n_u * n_u)]
    ## within this k-block, columns are ordered (p,q) with q fastest
    add <- numeric(n_u * n_u)
    for (p in seq_len(n_u)) for (q in seq_len(n_u)) {
      add[(p - 1L) * n_u + q] <- Sigma_e[p, q]
    }
    C[, row_cols] <- C[, row_cols, drop = FALSE] + outer(Cov_xi_x1[, k], add)
  }

  ## j12 = eps (x) x1 (x) eps : col index (p,k,q) with q FASTEST, k MIDDLE, p SLOWEST
  for (p in seq_len(n_u)) {
    for (k in seq_len(n_s)) {
      col_base <- ((p - 1L) * n_s + (k - 1L)) * n_u
      add <- Sigma_e[p, ]   # indexed by q, q fastest
      C[, jn$j12[col_base + seq_len(n_u)]] <-
        C[, jn$j12[col_base + seq_len(n_u)], drop = FALSE] + outer(Cov_xi_x1[, k], add)
    }
  }

  C
}


# ============================================================================
# Augmented-system assembly (Tlin3 / G3 / Dxi3 / Gv3)
# ============================================================================

#' Build the order-3 AFVRR augmented linear system for a DecisionRules3 object
#'
#' Internal analogue of \code{.order2_aug_system} for order-3 decision
#' rules.  See the file header of \code{pruned-state-space-order3.R} for the
#' full derivation, index conventions, and an explicit discussion of the
#' correlated-innovation issue flagged in the AFVRR paper (resolved here by
#' evaluating raw-innovation second moments at the exact stationary
#' marginal laws of x1/x2, see \code{.order3_cov_r}).
#'
#' @param dr3     A \code{DecisionRules3} object.
#' @param Sigma_e n_exo x n_exo shock covariance.
#' @keywords internal
.order3_aug_system <- function(dr3, Sigma_e) {
  ghx  <- dr3$ghx;  ghu  <- dr3$ghu
  ghxx <- dr3$ghxx; ghxu <- dr3$ghxu; ghuu <- dr3$ghuu; ghss <- dr3$ghss
  ghxxx <- dr3$ghxxx; ghxxu <- dr3$ghxxu; ghxuu <- dr3$ghxuu; ghuuu <- dr3$ghuuu
  ys <- dr3$ys

  sidx <- dr3$state_idx; endo <- dr3$endo_names
  n_endo <- length(endo); n_u <- ncol(ghu); n_s <- length(sidx)

  hx   <- ghx [sidx, , drop = FALSE]
  hu   <- ghu [sidx, , drop = FALSE]
  hxx  <- ghxx[sidx, , drop = FALSE]
  hxu  <- ghxu[sidx, , drop = FALSE]
  huu  <- ghuu[sidx, , drop = FALSE]
  hss  <- ghss[sidx]
  hxxx <- ghxxx[sidx, , drop = FALSE]
  hxxu <- ghxxu[sidx, , drop = FALSE]
  hxuu <- ghxuu[sidx, , drop = FALSE]
  huuu <- ghuuu[sidx, , drop = FALSE]

  has_xss <- !is.null(dr3$ghxss)
  has_uss <- !is.null(dr3$ghuss)
  has_s3  <- !is.null(dr3$ghs3)
  hxss <- if (has_xss) dr3$ghxss[sidx, , drop = FALSE] else matrix(0, n_s, n_s)
  huss <- if (has_uss) dr3$ghuss[sidx, , drop = FALSE] else matrix(0, n_s, n_u)
  hs3  <- if (has_s3)  dr3$ghs3[sidx]                  else numeric(n_s)

  ghxss_full <- if (has_xss) dr3$ghxss else matrix(0, n_endo, n_s)
  ghuss_full <- if (has_uss) dr3$ghuss else matrix(0, n_endo, n_u)
  ghs3_full  <- if (has_s3)  dr3$ghs3  else numeric(n_endo)

  vecSe <- as.numeric(Sigma_e)
  I_ns  <- diag(n_s)

  `%X%` <- function(A, B) kronecker(A, B)

  ## ---- Index blocks --------------------------------------------------------
  d <- 3L * n_s + 2L * n_s * n_s + n_s * n_s * n_s
  ix1  <- seq_len(n_s)
  ix2  <- (n_s + 1L):(2L * n_s)
  ik2  <- (2L * n_s + 1L):(2L * n_s + n_s * n_s)
  ix3  <- (2L * n_s + n_s * n_s + 1L):(2L * n_s + n_s * n_s + n_s)
  ik12 <- (2L * n_s + n_s * n_s + n_s + 1L):(2L * n_s + n_s * n_s + n_s + n_s * n_s)
  ik3  <- (2L * n_s + n_s * n_s + n_s + n_s * n_s + 1L):d
  stopifnot(length(ik3) == n_s^3)

  ## ---- Raw innovation category sizes/index ---------------------------------
  d1  <- n_u; d2 <- n_u * n_s; d3 <- n_s * n_u; d4 <- n_u * n_u
  d5  <- n_u * n_s; d6 <- n_u * n_s * n_s; d7 <- n_u * n_u * n_s; d8 <- n_u * n_u * n_u
  d9  <- n_s * n_s * n_u; d10 <- n_s * n_u * n_s; d11 <- n_s * n_u * n_u; d12 <- n_u * n_s * n_u
  sizes <- c(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12)
  Dr <- sum(sizes)
  offs <- c(0, cumsum(sizes))
  jn <- lapply(seq_len(12), function(k) (offs[k] + 1L):offs[k + 1L])
  names(jn) <- paste0("j", seq_len(12))

  Tlin <- matrix(0, d, d)
  G    <- matrix(0, d, Dr)
  cc   <- numeric(d)
  c_u  <- numeric(d)

  ## ==== ix1 (unchanged from order 2) =========================================
  Tlin[ix1, ix1] <- hx
  G[ix1, jn$j1]  <- hu

  ## ==== ix2 (unchanged from order 2) =========================================
  Tlin[ix2, ix2] <- hx
  Tlin[ix2, ik2] <- 0.5 * hxx
  G[ix2, jn$j2]  <- hxu
  G[ix2, jn$j4]  <- 0.5 * huu
  cc[ix2]  <- 0.5 * hss

  ## ==== ik2 (unchanged from order 2) =========================================
  Tlin[ik2, ik2] <- hx %X% hx
  G[ik2, jn$j2]  <- hu %X% hx
  G[ik2, jn$j3]  <- hx %X% hu
  G[ik2, jn$j4]  <- hu %X% hu

  ## ==== ix3 (NEW): x1_{t+1} = ... ; x3_new = hx*x3 + hxx*(x1(x)x2) + hxu*(e(x)x2)
  ##      + 0.5*hxxu*(e(x)x1(x)x1) + 0.5*hxuu*(e(x)e(x)x1) + (1/6)*hxxx*(x1(x)x1(x)x1)
  ##      + (1/6)*huuu*(e(x)e(x)e) [+ 0.5*hxss*x1 + 0.5*huss*e + (1/6)*hs3]
  Tlin[ix3, ix1]  <- 0.5 * hxss
  Tlin[ix3, ix3]  <- hx
  Tlin[ix3, ik12] <- hxx           ## coefficient 1 -- see file header note
  Tlin[ix3, ik3]  <- (1 / 6) * hxxx
  G[ix3, jn$j1]   <- 0.5 * huss
  G[ix3, jn$j5]   <- hxu
  G[ix3, jn$j6]   <- 0.5 * hxxu
  G[ix3, jn$j7]   <- 0.5 * hxuu
  G[ix3, jn$j8]   <- (1 / 6) * huuu
  cc[ix3] <- (1 / 6) * hs3

  ## ==== ik12 (NEW): x1_{t+1} (x) x2_{t+1} ====================================
  Tlin[ik12, ix1]  <- 0.5 * (I_ns %X% hss) %*% hx
  Tlin[ik12, ik12] <- hx %X% hx
  Tlin[ik12, ik3]  <- 0.5 * (hx %X% hxx)
  G[ik12, jn$j1]   <- 0.5 * (I_ns %X% hss) %*% hu
  G[ik12, jn$j5]   <- hu %X% hx
  G[ik12, jn$j6]   <- 0.5 * (hu %X% hxx)
  G[ik12, jn$j7]   <- hu %X% hxu
  G[ik12, jn$j8]   <- 0.5 * (hu %X% huu)
  G[ik12, jn$j10]  <- hx %X% hxu
  G[ik12, jn$j11]  <- 0.5 * (hx %X% huu)

  ## ==== ik3 (NEW): x1_{t+1}^(x)3 ==============================================
  Tlin[ik3, ik3] <- hx %X% hx %X% hx
  G[ik3, jn$j9]   <- hx %X% hx %X% hu
  G[ik3, jn$j10]  <- hx %X% hu %X% hx
  G[ik3, jn$j6]   <- G[ik3, jn$j6, drop = FALSE] + (hu %X% hx %X% hx)
  G[ik3, jn$j11]  <- hx %X% hu %X% hu
  G[ik3, jn$j12]  <- hu %X% hx %X% hu
  G[ik3, jn$j7]   <- G[ik3, jn$j7, drop = FALSE] + (hu %X% hu %X% hx)
  G[ik3, jn$j8]   <- G[ik3, jn$j8, drop = FALSE] + (hu %X% hu %X% hu)

  ## ---- c_u: deterministic-in-mean shock contributions (E[r_t] under the
  ## stationary raw-category means; only the "pure eps-power" categories
  ## with a nonzero mean survive -- j4 (eps(x)eps) has mean vec(Sigma_e),
  ## j8 (eps(x)eps(x)eps) has mean 0 (odd moment); everything else that
  ## mixes eps with x1/x2 has mean 0 at the stationary a=0 evaluation.
  c_u[ix2] <- 0.5 * as.numeric(huu %*% vecSe)
  c_u[ik2] <- as.numeric((hu %X% hu) %*% vecSe)
  ## ix3, ik12, ik3 pick up no NEW mean-shift terms beyond what's already in
  ## Tlin/cc (their j4/j8-loading blocks are all zero -- confirmed by
  ## inspection of the G[ix3,*]/G[ik12,*]/G[ik3,*] assignments above: none
  ## of them touch jn$j4).

  ## ---- Observation map (Dxi3 / Gv3 / c_v3) ----------------------------------
  ## NOTE (bug found + fixed during P1 validation, see the completion report):
  ## the observation equation's y3 term `0.5 * ghxss %*% x1_prev` (mirroring
  ## the state equation's `Tlin[ix3,ix1] = 0.5*hxss`) is a DETERMINISTIC
  ## LINEAR function of x1 -- it belongs in Dxi[,ix1] (added ON TOP OF ghx,
  ## the y1 first-order term), not in Gv (which only carries shock-loaded
  ## terms).  Omitting this addition was a real bug (caught by the
  ## rbc2shock one-step oracle at n_s=3, invisible on the n_s=2 fixture
  ## because ghxss happened to be numerically small there relative to MC
  ## noise) -- see the file header "one-step transition oracle" discussion.
  Dxi <- matrix(0, n_endo, d)
  Dxi[, ix1]  <- ghx + 0.5 * ghxss_full
  Dxi[, ix2]  <- ghx
  Dxi[, ik2]  <- 0.5 * ghxx
  Dxi[, ix3]  <- ghx
  Dxi[, ik12] <- ghxx            ## coefficient 1 -- see file header note
  Dxi[, ik3]  <- (1 / 6) * ghxxx

  Gv <- matrix(0, n_endo, Dr)
  Gv[, jn$j1] <- ghu
  Gv[, jn$j1] <- Gv[, jn$j1, drop = FALSE] + 0.5 * ghuss_full
  Gv[, jn$j2] <- ghxu
  Gv[, jn$j4] <- 0.5 * ghuu
  Gv[, jn$j5] <- ghxu
  Gv[, jn$j6] <- 0.5 * ghxxu
  Gv[, jn$j7] <- 0.5 * ghxuu
  Gv[, jn$j8] <- (1 / 6) * ghuuu

  r_mean <- numeric(Dr)
  r_mean[jn$j4] <- vecSe
  c_v <- as.numeric(Gv %*% r_mean) + (1 / 6) * ghs3_full

  list(
    Tlin = Tlin, cc = cc, c_u = c_u, G = G, Dxi = Dxi, Gv = Gv, c_v = c_v,
    ghss = ghss, ys = ys, hx = hx, hu = hu, Sigma_e = Sigma_e,
    n_s = n_s, n_u = n_u, n_endo = n_endo, d = d, Dr = Dr,
    ix1 = ix1, ix2 = ix2, ik2 = ik2, ix3 = ix3, ik12 = ik12, ik3 = ik3,
    jn = jn, sizes = sizes, endo = endo
  )
}


# ============================================================================
# Stationary moments (order 3)
# ============================================================================

#' Stationary (unconditional) order-3 output moments
#'
#' Internal analogue of \code{.order2_stationary_moments}.  Solves for the
#' order-2 marginal moments of \code{(x1, x2)} first (exactly reusing
#' \code{.order2_stationary_moments} machinery via the leading block of the
#' order-3 system, which is byte-identical to the order-2 system), then
#' evaluates \code{.order3_cov_r} at those marginals and solves the full
#' \eqn{d3 \times d3} Lyapunov fixed point for the augmented covariance.
#'
#' @param sys An order-3 augmented system list from \code{.order3_aug_system}.
#' @keywords internal
.order3_stationary_moments <- function(sys) {
  ## ---- Order-2 marginal moments of (x1, x2) (leading-block reuse) ----------
  Sigma_x <- solve_lyapunov(sys$hx, sys$hu %*% sys$Sigma_e %*% t(sys$hu))

  ## Build a lightweight order-2-shaped sys to reuse .order2_stationary_moments
  ## for the (x1,x2,x1(x)x1) marginal exactly (leading block is byte-identical).
  sys2 <- list(
    Tlin = sys$Tlin[c(sys$ix1, sys$ix2, sys$ik2), c(sys$ix1, sys$ix2, sys$ik2), drop = FALSE],
    cc   = sys$cc[c(sys$ix1, sys$ix2, sys$ik2)],
    c_u  = sys$c_u[c(sys$ix1, sys$ix2, sys$ik2)],
    G    = sys$G[c(sys$ix1, sys$ix2, sys$ik2), c(sys$jn$j1, sys$jn$j2, sys$jn$j3, sys$jn$j4), drop = FALSE],
    Dxi  = matrix(0, 1, length(sys$ix1) + length(sys$ix2) + length(sys$ik2)),  # unused here
    Gv   = matrix(0, 1, length(sys$jn$j1) + length(sys$jn$j2) + length(sys$jn$j3) + length(sys$jn$j4)),
    c_v  = 0, ghss = sys$ghss, hx = sys$hx, hu = sys$hu, Sigma_e = sys$Sigma_e,
    n_s = sys$n_s, n_u = sys$n_u, n_endo = 1L,
    d = length(sys$ix1) + length(sys$ix2) + length(sys$ik2),
    ix1 = seq_len(sys$n_s), ix2 = (sys$n_s + 1L):(2L * sys$n_s),
    ik  = (2L * sys$n_s + 1L):(2L * sys$n_s + sys$n_s^2), endo = "_dummy_"
  )
  st2 <- .order2_stationary_moments(sys2)
  mean_x2 <- st2$mean_x2
  Var_x2  <- st2$Var_x2
  Cov_x2_x11 <- st2$Cov_x2_x11   # exact Cov(x2, x1(x)x1) for the j5 cross-blocks

  ## ---- Cov(r3) at the stationary (a=0, Sigma_x, mean_x2, Var_x2) point -----
  Cr0 <- .order3_cov_r(numeric(sys$n_s), Sigma_x, mean_x2, Var_x2, sys$Sigma_e,
                       Cov_x2_x11 = Cov_x2_x11)
  QQ0 <- sys$G %*% Cr0 %*% t(sys$G)
  QQ0 <- (QQ0 + t(QQ0)) * 0.5

  ## ---- Full augmented-state stationary covariance -------------------------
  ##
  ## CORRELATED z_t/xi_{t+1} FIXED POINT (AFVRR p.11 / file header note):
  ## unlike order 2, xi_t is CORRELATED with r_t (via the j7/j11/j12 raw
  ## categories, which multiply a CURRENT-period x1_t/x2_t factor by
  ## eps_t(x)eps_t).  The correct stationary covariance solves
  ##   Sxi = Tlin Sxi Tlin' + G Cr0 G' + Tlin Cxr G' + G Cxr' Tlin'
  ## where Cxr = Cov(xi_t, r_t) depends on Sxi itself (via Sxi[,ix1] =
  ## Cov(xi_t, x1_t)) -- a linear fixed point, solved here by DAMPED
  ## PICARD ITERATION seeded at the (wrong, order-2-style) Cxr=0 solution.
  ## This is NOT solve_lyapunov(Tlin, QQ) alone (that omits the Cxr cross
  ## term entirely and was empirically confirmed WRONG by up to ~37%
  ## relative variance error on rbc2shock during P1 validation -- see the
  ## P1 completion report).  Convergence is fast in practice (the cross
  ## term is a second-order correction on top of the dominant QQ0 term);
  ## a hard iteration cap with a diagnostic warning guards against a
  ## non-converging edge case rather than looping silently.
  Sxi <- solve_lyapunov(sys$Tlin, QQ0)
  Sxi <- (Sxi + t(Sxi)) * 0.5

  max_iter_xr <- 100L
  tol_xr <- 1e-12
  converged_xr <- FALSE
  for (iter in seq_len(max_iter_xr)) {
    Cov_xi_x1 <- Sxi[, sys$ix1, drop = FALSE]
    Cxr <- .order3_cov_xi_r(Cov_xi_x1, sys$Sigma_e, sys$jn, sys$d, sys$Dr)
    cross <- sys$Tlin %*% Cxr %*% t(sys$G)
    QQ_full <- QQ0 + cross + t(cross)
    QQ_full <- (QQ_full + t(QQ_full)) * 0.5
    Sxi_new <- solve_lyapunov(sys$Tlin, QQ_full)
    Sxi_new <- (Sxi_new + t(Sxi_new)) * 0.5
    d_max <- max(abs(Sxi_new - Sxi)) / max(1, max(abs(Sxi_new)))
    Sxi <- Sxi_new
    if (d_max < tol_xr) { converged_xr <- TRUE; break }
  }
  if (!converged_xr)
    warning("`.order3_stationary_moments`: Cov(xi,r) fixed-point iteration ",
            "did not converge to tol=", tol_xr, " within ", max_iter_xr,
            " iterations; result may be inaccurate.")

  ## Observation covariance ALSO needs the Cov(xi_t, r_t) cross term (y_t
  ## depends on xi_t AND r_t at the SAME t via Gv): Var(y_t) =
  ##   Dxi Sxi Dxi' + Gv Cr0 Gv' + Dxi Cxr Gv' + Gv Cxr' Dxi'
  ## using the CONVERGED Cxr from the fixed point above (same correlated-
  ## innovation issue as the state covariance, and the same fix).
  cross_y <- sys$Dxi %*% Cxr %*% t(sys$Gv)
  var_cov <- sys$Dxi %*% Sxi %*% t(sys$Dxi) + sys$Gv %*% Cr0 %*% t(sys$Gv) +
             cross_y + t(cross_y)
  var_cov <- (var_cov + t(var_cov)) * 0.5

  mu_xi <- as.numeric(solve(diag(sys$d) - sys$Tlin, sys$cc + sys$c_u))
  mean_dev <- as.numeric(sys$Dxi %*% mu_xi) + 0.5 * sys$ghss + sys$c_v

  list(
    var_cov = var_cov, mean = sys$ys + mean_dev,
    Sigma_x = Sigma_x, Var_x2 = Var_x2, mean_x2 = mean_x2,
    Sxi = Sxi, mu_xi = mu_xi, Cr0 = Cr0
  )
}


# ============================================================================
# Public constructor + moments (P1: assembly/moments only, no KF/loglik)
# ============================================================================

#' Build an order-3 pruned state-space object
#'
#' The order-3 analogue of \code{\link{pruned_state_space}}.  Builds the AFVRR
#' (2018) order-3 augmented linear system
#' (\code{xi3_t = [x1;x2;x1(x)x1;x3;x1(x)x2;x1(x)x1(x)x1]}) via
#' \code{.order3_aug_system}.  The Kalman-filter log-likelihood on this
#' augmented state is \code{\link{pruned_ss_loglik3}}, and
#' \code{make_log_posterior(pruned_order = 3L)} wires it into estimation
#' (both merged in this file).
#'
#' @param dr3    A \code{DecisionRules3} object (output of
#'   \code{solve_perturbation(order = 3)}).
#' @param model  A parsed model object.
#' @param params Named numeric parameter vector.  Defaults to
#'   \code{model$param_values}.
#' @return An object of class \code{"pruned_ss3"}.
#' @seealso \code{\link{pruned_ss_moments3}}
#' @export
pruned_state_space3 <- function(dr3, model, params = NULL) {
  if (!inherits(dr3, "DecisionRules3"))
    stop("pruned_state_space3: dr3 must be a DecisionRules3 object.")

  if (is.null(params)) params <- model$param_values

  Sigma_e <- .get_shock_cov(model, dr3$exo_names, params)
  sys     <- .order3_aug_system(dr3, Sigma_e)

  structure(
    list(
      sys = sys, Sigma_e = Sigma_e, ys = dr3$ys,
      endo_names = dr3$endo_names, exo_names = dr3$exo_names,
      state_idx = dr3$state_idx,
      n_s = sys$n_s, n_u = sys$n_u, n_endo = sys$n_endo, d = sys$d,
      dr = dr3, model = model
    ),
    class = "pruned_ss3",
    ## The former ~25% control-variance bias is FIXED (2026-07-02): it was
    ## NOT a Cov(x1_t,x2_t) approximation (that covariance is exactly 0 --
    ## x1 degree-1, x2 degree-2 in the Gaussian shocks).  The real cause was
    ## the j5 (=eps(x)x2) cross-category blocks of Cr0: the zero-mean Wick
    ## machinery dropped (a) x2's NONZERO mean and (b) the connected
    ## non-Gaussian moment Cov(x2, x1(x)x1).  Both are now added in closed
    ## form (.order3_cov_r's x2-cross-category correction, using the exact
    ## order-2 mean_x2 and Cov(x2, x1(x)x1) blocks).  Validated observable-
    ## level against a companion-simulation MC oracle on rbc2shock: control
    ## variances go from ~25-48% biased to ~1-2% (= MC sampling noise).
    variance_approx = FALSE
  )
}

#' Compute unconditional moments from an order-3 pruned state-space object
#'
#' P1 sub-increment: analytic stationary mean and variance of the
#' observables under the order-3 AFVRR pruned recursion, via
#' \code{.order3_stationary_moments}.
#'
#' @param pss3 A \code{pruned_ss3} object from \code{\link{pruned_state_space3}}.
#' @return A list with \code{mean} (named n_endo vector), \code{var_cov}
#'   (n_endo x n_endo), \code{std_dev}, \code{Sigma_x}, \code{Var_x2},
#'   \code{mean_x2}.
#' @export
pruned_ss_moments3 <- function(pss3) {
  stopifnot(inherits(pss3, "pruned_ss3"))
  sys <- pss3$sys
  st  <- .order3_stationary_moments(sys)

  mn <- st$mean
  names(mn) <- pss3$endo_names
  Sigma_y <- st$var_cov
  rownames(Sigma_y) <- pss3$endo_names
  colnames(Sigma_y) <- pss3$endo_names

  variances <- pmax(diag(Sigma_y), 0)
  std_dev   <- sqrt(variances)
  names(std_dev) <- pss3$endo_names

  list(
    mean = mn, var_cov = Sigma_y, std_dev = std_dev,
    Sigma_e = pss3$Sigma_e, Sigma_x = st$Sigma_x,
    Var_x2 = st$Var_x2, mean_x2 = st$mean_x2
  )
}


# ============================================================================
# P2: Gaussian Kalman-filter log-likelihood on the order-3 augmented state
# ============================================================================

#' Linear conditional-mean coefficient of the order-3 raw innovation
#'
#' At order 3 the raw innovation \eqn{r_t} is CORRELATED with the current
#' augmented state \eqn{xi_t}, because its categories \code{j7} (eps(x)eps(x)x1),
#' \code{j11} (x1(x)eps(x)eps) and \code{j12} (eps(x)x1(x)eps) each carry an EVEN
#' power of the fresh shock \eqn{eps_{t+1}} multiplying a CURRENT-period
#' \eqn{x1_t} factor.  Taking the expectation over \eqn{eps_{t+1}} (independent of
#' \eqn{xi_t}) gives \eqn{E[r_t \mid xi_t] = A xi_t + b}, where the LINEAR part
#' \eqn{A} (this function) is nonzero only on the \code{ix1} columns and only for
#' those three categories (the value is \eqn{Sigma_e[p,q]} per Kronecker slot),
#' and the constant \eqn{b} (the pure-\code{j4} \eqn{vec(Sigma_e)} mean) is
#' already folded into the state/obs intercepts \code{cc + c_u} / \code{c_v}.
#'
#' This is exactly the structure of \code{\link{.order3_cov_xi_r}} (which returns
#' \eqn{Cov(xi_t, r_t) = Sxi[, ix1] A^\top}); evaluating that routine with the
#' \dQuote{covariance} argument set to the \code{ix1}-identity extracts
#' \eqn{A^\top} directly, so the two share one validated index convention.
#'
#' @return A \code{Dr x d} matrix \eqn{A} with \eqn{E[r_t \mid xi_t] = A xi_t + b}.
#' @keywords internal
.order3_cond_mean_r_coef <- function(Sigma_e, jn, d, Dr, ix1, n_s) {
  E <- matrix(0, d, n_s)
  E[ix1, ] <- diag(n_s)
  t(.order3_cov_xi_r(E, Sigma_e, jn, d, Dr))     # Dr x d
}

#' Order-3 pruned state-space Gaussian log-likelihood (P2)
#'
#' The order-3 analogue of \code{\link{pruned_ss_loglik}}: a constant-coefficient
#' Gaussian Kalman filter on the AFVRR order-3 augmented state.  Unlike order 2
#' (where the augmented state is uncorrelated with its own innovation), order 3
#' has \eqn{Cov(xi_t, r_t) \neq 0}.  This is handled EXACTLY (to second order) by
#' folding the linear conditional mean \eqn{E[r_t \mid xi_t] = A xi_t + b} into
#' effective system matrices, which leaves a residual innovation
#' \eqn{u_t = r_t - E[r_t \mid xi_t]} that is uncorrelated with the state:
#' \deqn{xi_{t+1} = (Tlin + G A) xi_t + c + G u_t, \quad
#'       y_t = (Dxi + Gv A) xi_t + d + Gv u_t.}
#' The residual covariance \eqn{Cu = Cr0 - A\,Sxi\,A^\top} (the mean of the
#' state-conditional innovation covariance, PSD by the law of total variance) and
#' the standard correlated-noise Kalman recursion (\code{.pruned_kf_correlated},
#' reused from order 2) then give the likelihood.
#'
#' By construction the folded system's stationary predictive covariance equals
#' the P1 augmented covariance \code{Sxi}, and the filter's first innovation
#' covariance equals the P1 (MC-validated) observable \code{var_cov} exactly --
#' so this likelihood is anchored to the validated order-3 moments and collapses
#' to \code{pruned_ss_loglik} when the order-3 tensors vanish (\eqn{A = 0}).
#'
#' \strong{Approximation.} Like the order-2 pruned likelihood this is a GAUSSIAN
#' quasi-likelihood: the pruned state space is polynomial (non-Gaussian) in the
#' shocks, and \eqn{Cu} is held at its stationary value rather than the true
#' state-dependent one-step innovation covariance.  It matches the model's first
#' and second conditional moments exactly; higher-moment/state-dependent
#' heteroskedasticity is not captured.
#'
#' \strong{Certification.} The likelihood's correctness (as a Gaussian density
#' of the folded linear SSM) is independently established by Oracle B (a
#' from-scratch batch multivariate-normal computation) and by the MC-validated
#' P1 moments; a SELF-CONSISTENT SBC -- simulating from this SAME folded law at
#' \eqn{\theta} and checking rank uniformity of the posterior draws -- is the
#' correct way to certify the SAMPLER and the \code{pruned_order = 3L}
#' posterior wiring (see \code{ORDER3_PRUNED_SS_FOLLOWUP.md}, section "P2c").
#' SBC that instead simulates from the TRUE NONLINEAR model
#' (\code{simulate_model_order3}) conflates the deliberate Gaussian-quasi-
#' likelihood approximation error with genuine sampler/wiring bugs; that is an
#' approximation-quality diagnostic, not a correctness test, and should not be
#' used as the primary certification gate.
#'
#' @param pss3 A \code{pruned_ss3} object.
#' @param Y Observations: \code{n_obs x T} (or \code{T x n_obs}, auto-oriented).
#'   \code{NA} entries get an exact PARTIAL measurement update: a period with
#'   SOME (but not all) observables missing still updates on the observed
#'   subset (subsetting \code{ZZ}/\code{HH}/\code{SS} and using the
#'   observed-count normalizing constant); a period with ALL observables
#'   missing falls back to predict-only (no update).
#' @param obs_vars Character vector of observed endogenous variable names.
#' @param me_variance Scalar measurement-error variance added to the observation
#'   covariance diagonal (default 0). See \code{\link{pruned_ss_loglik}} for
#'   the "Riccati lock" hazard on near-perfectly-revealed observables.
#' @param me_floor_check Logical: warn when the \code{me_variance} floor
#'   materially inflates the filter's steady-state innovation variances vs
#'   \code{me_variance = 0}. Default
#'   \code{getOption("dynhr.me_floor_check", TRUE)}.
#' @return Scalar Gaussian log-likelihood (\code{-Inf} on a non-PD innovation
#'   covariance).
#' @seealso \code{\link{pruned_ss_loglik}}, \code{\link{pruned_ss_moments3}}
#' @export
pruned_ss_loglik3 <- function(pss3, Y, obs_vars, me_variance = 0,
                              me_floor_check = getOption("dynhr.me_floor_check",
                                                         TRUE)) {
  stopifnot(inherits(pss3, "pruned_ss3"))
  sys <- pss3$sys

  n_obs <- length(obs_vars)
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)

  inp <- .order3_pruned_kf_inputs(sys, pss3$ys, obs_vars, pss3$endo_names,
                                  me_variance = me_variance)

  ## Measurement-error floor guard (see pruned_ss_loglik / the 2026-07-03
  ## MAJOR CORRECTION in ORDER3_PRUNED_SS_FOLLOWUP.md).
  if (me_variance > 0 && isTRUE(me_floor_check))
    .warn_me_floor_lock(
      .pruned_me_floor_ratio(inp$Tlin, inp$ZZ, inp$QQ,
                             inp$HH - me_variance * diag(n_obs),
                             inp$SS, inp$Sxi0, me_variance),
      obs_vars, me_variance)

  ## Stationary initial condition: unconditional augmented mean / covariance.
  .pruned_kf_correlated(Y, inp$Tlin, inp$ZZ, inp$d_y, inp$c_drift,
                        inp$QQ, inp$HH, inp$SS, inp$mu0, inp$Sxi0)
}

#' Build the nine order-3 pruned-SS Kalman-filter inputs
#'
#' Pure extract-function refactor (D1, zero behaviour change): factors the
#' assembly block that used to live inline in \code{pruned_ss_loglik3} (fold
#' of the linear state-innovation correlation into effective matrices, plus
#' the observation-side subsetting) into a single internal builder, so it can
#' be called from BOTH \code{pruned_ss_loglik3} and the semi-analytic order-3
#' gradient chain (\code{R/pruned-grad-chain-order3.R}) without duplicating
#' the assembly. Guarded by a byte-identity test
#' (test-pruned-grad-order3.R "D1 refactor guard") against a hard-coded
#' pre-refactor reference loglik.
#'
#' @param sys      Order-3 augmented system list from \code{.order3_aug_system}.
#' @param ys       Named steady-state vector (\code{dr3$ys} / \code{pss3$ys}).
#' @param obs_vars Character vector of observed endogenous variable names.
#' @param endo_names Character vector, \code{pss3$endo_names} (for matching
#'   \code{obs_vars}).
#' @param me_variance Scalar measurement-error variance (default 0).
#' @return A list with the nine \code{.pruned_kf_correlated}/adjoint inputs
#'   (\code{Tlin, ZZ, d_y, c_drift, QQ, HH, SS, mu0, Sxi0}), plus the
#'   intermediate \code{st} (stationary-moments list), \code{A} (linear
#'   conditional-mean coefficient) and \code{Cu} (residual covariance) for
#'   reuse by the gradient chain.
#' @keywords internal
.order3_pruned_kf_inputs <- function(sys, ys, obs_vars, endo_names, me_variance = 0) {
  obs_idx <- match(obs_vars, endo_names)
  if (any(is.na(obs_idx)))
    stop(".order3_pruned_kf_inputs: obs_vars not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))
  n_obs <- length(obs_vars)

  ## P1 stationary moments supply Sxi (augmented cov), Cr0 (raw-innovation cov)
  ## and mu_xi (augmented mean) -- all MC-validated.
  st  <- .order3_stationary_moments(sys)
  Sxi <- st$Sxi; Cr0 <- st$Cr0; mu_xi <- st$mu_xi

  ## Fold the linear state-innovation correlation into effective matrices.
  A <- .order3_cond_mean_r_coef(sys$Sigma_e, sys$jn, sys$d, sys$Dr, sys$ix1, sys$n_s)
  Cu <- Cr0 - A %*% Sxi %*% t(A)
  Cu <- (Cu + t(Cu)) * 0.5

  Tlin_eff <- sys$Tlin + sys$G %*% A            # d x d
  Dxi_eff  <- sys$Dxi  + sys$Gv %*% A           # n_endo x d
  ZZ <- Dxi_eff[obs_idx, , drop = FALSE]        # n_obs x d
  Gv <- sys$Gv[obs_idx, , drop = FALSE]         # n_obs x Dr
  G  <- sys$G

  QQ <- G  %*% Cu %*% t(G);  QQ <- (QQ + t(QQ)) * 0.5     # d x d
  HH <- Gv %*% Cu %*% t(Gv); HH <- (HH + t(HH)) * 0.5     # n_obs x n_obs
  if (me_variance > 0) HH <- HH + me_variance * diag(n_obs)
  SS <- G  %*% Cu %*% t(Gv)                                # d x n_obs

  ## Observation intercept (the constant b is already inside c_v / c_u); the
  ## Dxi*xi part is the state, so the KF intercept is ys + 0.5*ghss + c_v.
  d_y     <- ys[obs_vars] + 0.5 * sys$ghss[obs_idx] + sys$c_v[obs_idx]
  c_drift <- sys$cc + sys$c_u

  list(
    Tlin = Tlin_eff, ZZ = ZZ, d_y = as.numeric(d_y), c_drift = c_drift,
    QQ = QQ, HH = HH, SS = SS, mu0 = mu_xi, Sxi0 = Sxi,
    obs_idx = obs_idx, st = st, A = A, Cu = Cu, Gv_obs = Gv
  )
}

#' Build a log-posterior using the order-3 pruned-SS Gaussian likelihood (P2b)
#'
#' Order-3 analogue of \code{make_log_posterior_pruned}: at each draw it solves
#' the model to order 3, builds the AFVRR order-3 pruned state space, and
#' evaluates \code{\link{pruned_ss_loglik3}}.  Dispatched from
#' \code{make_log_posterior} when \code{pruned_order = 3L}.  No analytic
#' gradient (gradient optimizers fall back to finite differences).
#'
#' @inheritParams make_log_posterior_pruned
#' @return A closure \code{function(theta)} returning
#'   \code{list(logpost, loglik, logprior)}.
#' @noRd
make_log_posterior_pruned3 <- function(model, data, prior_spec, obs_vars,
                                        compiled, me_variance = 0,
                                        system_priors = NULL,
                                        power = NULL) {
  ## Resolve zeta ONCE here, not per draw (see .resolve_power_posterior).
  power <- .resolve_power_posterior(power, "make_log_posterior_pruned3")
  n_obs <- length(obs_vars)
  Y <- if (is.matrix(data) && nrow(data) == n_obs) data else t(data)

  ## Adapter over the shared closure builder (R/posterior-closure.R); the
  ## order-3 twin of make_log_posterior_pruned's. `sys_cache = FALSE` because
  ## this branch goes straight to solve_perturbation(order = 3L) and never
  ## touches extract_system_matrices_fast(), so building the cache would be
  ## factory-time work for nothing.
  .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    sys_cache = FALSE,
    solve_fn = function(model, compiled, sys_cache, ss, params, theta) {
      dr3 <- tryCatch(
        solve_perturbation(model, compiled, ss, params,
                           order = 3L, verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(dr3) || !isTRUE(dr3$bk_satisfied)) return(NULL)
      pss3 <- tryCatch(pruned_state_space3(dr3, model, params),
                       error = function(e) NULL)
      if (is.null(pss3)) return(NULL)
      list(dr = dr3, pss = pss3)
    },
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      loglik <- tryCatch(
        pruned_ss_loglik3(sol$pss, Y, obs_vars, me_variance = me_variance,
                          me_floor_check = me_floor_check),
        error = function(e) -Inf
      )
      if (!is.finite(loglik)) return(NULL)
      list(loglik = loglik, Sigma_e = sol$pss$Sigma_e)
    },
    power             = power,
    warm_retry        = FALSE,
    system_prior      = system_priors,
    system_prior_mode = "extra")
}
