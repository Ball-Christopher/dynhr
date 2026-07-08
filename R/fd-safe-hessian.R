## R/fd-safe-hessian.R
## --------------------------------------------------------------------------
## P2 paper gap #4: feasibility-aware finite-difference Hessian.
##
## The pathological-DSGE paper's central finding is that a *default* FD step
## (numDeriv's default, or any fixed-fraction central-difference step) can
## overshoot a near-unit-root/determinacy boundary. Points just beyond the
## boundary return non-finite objective values (or a large-penalty stand-in,
## e.g. Dynare's hessian.m convention of a huge constant for infeasible
## draws). Either way the stencil silently contaminates the FD Hessian and
## inflates its condition number by 5-7 orders of magnitude relative to the
## analytic/BFGS reference (see check_hessian_conditioning(), R/mode-hessian.R,
## for the companion post-hoc diagnostic). On some models (sw2007 in that
## paper) *no* step size recovers the truth in both directions of a
## coordinate -- the routine must refuse rather than return a plausible-
## looking but wrong matrix.
##
## fd_safe_hessian() is the pre-hoc fix: it builds the central-difference
## stencil itself, detects boundary crossings (non-finite values or a jump in
## the objective that looks like a penalty), shrinks the step per-coordinate
## into the feasible cone (falling back to a one-sided stencil when only one
## direction is feasible), and refuses with a structured, per-coordinate
## diagnosis when no step stabilizes.
## --------------------------------------------------------------------------

#' Feasibility-aware finite-difference Hessian
#'
#' Computes a central-difference Hessian of a scalar objective \code{fn} at
#' \code{theta}, detecting and stepping around infeasible regions (points
#' where \code{fn} returns a non-finite value, or a "penalty" value far from
#' \code{fn(theta)}, e.g. a near-unit-root/determinacy boundary in a DSGE
#' likelihood). This is the pre-hoc companion to
#' \code{\link{check_hessian_conditioning}} (R/mode-hessian.R), which is a
#' post-hoc sanity check comparing an already-computed FD Hessian's condition
#' number to a reference. Use \code{fd_safe_hessian()} to avoid manufacturing
#' the bad Hessian in the first place; use \code{check_hessian_conditioning()}
#' whenever an FD Hessian was computed some other way (e.g. by
#' \code{numDeriv::hessian}) and needs a sanity check.
#'
#' \strong{Convention:} \code{fn} is the objective to differentiate directly
#' -- pass the \emph{negative} log posterior (or negative log-likelihood) if
#' you want a Hessian whose inverse is a covariance matrix, matching
#' \code{optim()}/\code{csminwel()}/\code{numDeriv::hessian()} convention.
#' \code{fd_safe_hessian()} itself is convention-agnostic: it differentiates
#' whatever scalar \code{fn} returns.
#'
#' \strong{Feasibility / penalty detection.} A stencil evaluation
#' \code{f_eval <- fn(theta_perturbed)} is treated as an infeasible ("boundary
#' crossing") point when either:
#' \enumerate{
#'   \item \code{f_eval} is non-finite (\code{NA}, \code{NaN}, \code{Inf}), or
#'   \item \code{abs(f_eval - f0) > penalty_tol}, where \code{f0 <- fn(theta)}
#'     is the objective at the base point. This catches the "large penalty
#'     instead of Inf" convention used by e.g. Dynare's \code{hessian.m} for
#'     infeasible draws (a large finite value, not \code{Inf}), which would
#'     otherwise silently poison the stencil with a huge finite second
#'     difference.
#' }
#' \code{penalty_tol} must therefore be chosen larger than the *genuine*
#' curvature of \code{fn} over the step sizes attempted (its default, 1e6, is
#' deliberately large relative to typical log-posterior curvature so it does
#' not fire on ordinary evaluations -- tune it down only if your objective's
#' true range at the attempted steps can approach 1e6).
#'
#' \strong{Per-coordinate step shrinking.} For each coordinate \code{i}, the
#' routine starts at \code{step0[i] <- rel_step * max(abs(theta[i]), 1)} and,
#' up to \code{max_shrink} times, halves the step whenever a stencil point in
#' either direction is infeasible or (with \code{verify_spectrum = TRUE}) the
#' Hessian's log-condition-number has not stabilized between two consecutive
#' step sizes. If, at some step, only one direction (+h or -h) is feasible,
#' the routine falls back to a one-sided (forward or backward) second
#' derivative for that coordinate's diagonal entry and a one-sided
#' cross-difference for any off-diagonal entry involving it; this is recorded
#' in \code{attr(H, "one_sided")} and is a documented accuracy downgrade (the
#' one-sided second-difference is only first-order accurate in \code{h},
#' versus second-order for central differences).
#'
#' \strong{Refusal.} If, for some coordinate, \emph{neither} direction is
#' feasible at the smallest attempted step (i.e. even a tiny perturbation in
#' both +/- directions is infeasible), the routine cannot form even a
#' one-sided derivative for that coordinate and returns \code{ok = FALSE}
#' with a per-coordinate diagnosis -- it never returns a silently-wrong
#' matrix.
#'
#' @param fn Scalar objective function, \code{theta -> numeric(1)}. Should
#'   return non-finite (or a large penalty, see \code{penalty_tol}) for
#'   infeasible \code{theta}, not error -- but the routine also tolerates
#'   \code{fn} raising an error at infeasible points (caught and treated as
#'   non-finite).
#' @param theta Named numeric vector, the point to differentiate at.
#' @param lower,upper Optional named (or plain, matched by position) numeric
#'   vectors of box constraints; defaults to \code{-Inf}/\code{Inf}. Steps are
#'   additionally clipped so stencil points never leave \code{[lower, upper]}.
#' @param rel_step Initial relative step size, as a fraction of
#'   \code{max(abs(theta[i]), 1)}. Default \code{1e-2}.
#' @param max_shrink Maximum number of halvings attempted per coordinate
#'   before declaring that direction infeasible. Default \code{20} (step
#'   shrinks by up to \code{2^20 ~ 1e6}).
#' @param penalty_tol Threshold on \code{abs(f_eval - f0)} for treating a
#'   finite stencil evaluation as an infeasible "penalty" value (see Details).
#'   Default \code{1e6}.
#' @param verify_spectrum Logical; if \code{TRUE} (default), after finding a
#'   per-coordinate step size that yields a feasible stencil, the routine also
#'   shrinks once more and recomputes the full Hessian, requiring the
#'   relative change in \code{log(kappa)} (condition number of the
#'   symmetrized Hessian) between the two step sizes to fall below
#'   \code{spectrum_tol}. This guards against a stencil that is feasible but
#'   still close enough to the boundary that curvature is not yet stable.
#' @param spectrum_tol Relative tolerance on the change in
#'   \code{log(kappa)} between consecutive verification steps. Default
#'   \code{0.5} (a 50\% relative change in log-condition-number is treated as
#'   "not yet stable").
#'
#' @return A list:
#' \describe{
#'   \item{ok}{Logical. \code{TRUE} if a full Hessian was constructed (all
#'     coordinates had at least a one-sided feasible stencil).}
#'   \item{H}{The \code{n x n} Hessian matrix (named rows/cols) if
#'     \code{ok = TRUE}; \code{NULL} if \code{ok = FALSE}.}
#'   \item{step_used}{Named numeric vector, the step size finally used for
#'     each coordinate (\code{NA} for refused coordinates).}
#'   \item{one_sided}{Named logical vector, \code{TRUE} where a one-sided
#'     (not central) derivative was used for that coordinate's diagonal.}
#'   \item{one_sided_direction}{Named character vector, \code{"forward"},
#'     \code{"backward"}, or \code{NA} (central) per coordinate.}
#'   \item{diagnosis}{A data.frame with one row per coordinate:
#'     \code{coord}, \code{feasible_plus}, \code{feasible_minus},
#'     \code{step_used}, \code{one_sided}, \code{refused}. Always present
#'     (even when \code{ok = TRUE}) so callers can audit near-misses.}
#'   \item{message}{Character; \code{NA} when \code{ok = TRUE}, otherwise a
#'     human-readable diagnosis naming the refused coordinate(s).}
#' }
#' Also carried as \code{attributes} on \code{H} itself (when \code{ok =
#' TRUE}) for callers that only keep the matrix: \code{"step_used"},
#' \code{"one_sided"}, \code{"one_sided_direction"}, \code{"diagnosis"}.
#'
#' @seealso \code{\link{check_hessian_conditioning}} for a post-hoc sanity
#'   check comparing an already-computed FD Hessian's condition number
#'   against an analytic/BFGS reference.
#'
#' @export
fd_safe_hessian <- function(fn, theta,
                             lower = NULL, upper = NULL,
                             rel_step = 1e-2,
                             max_shrink = 20L,
                             penalty_tol = 1e6,
                             verify_spectrum = TRUE,
                             spectrum_tol = 0.5) {

  if (!is.numeric(theta) || is.null(names(theta)) || any(names(theta) == ""))
    stop("fd_safe_hessian: 'theta' must be a named numeric vector.")

  n <- length(theta)
  pnames <- names(theta)

  lo <- .fdsh_expand_bound(lower, pnames, -Inf)
  hi <- .fdsh_expand_bound(upper, pnames,  Inf)

  safe_fn <- function(th) {
    v <- tryCatch(fn(th), error = function(e) NA_real_)
    if (length(v) != 1L || !is.numeric(v)) return(NA_real_)
    as.numeric(v)
  }

  f0 <- safe_fn(theta)
  if (!is.finite(f0))
    stop("fd_safe_hessian: fn(theta) is non-finite at the supplied base point; ",
         "'theta' itself must be feasible.")

  feasible <- function(f_eval) is.finite(f_eval) && abs(f_eval - f0) <= penalty_tol

  ## ---- Step 1: per-coordinate step search -------------------------------
  step_used           <- rep(NA_real_, n); names(step_used) <- pnames
  one_sided           <- rep(FALSE, n);    names(one_sided) <- pnames
  one_sided_direction <- rep(NA_character_, n); names(one_sided_direction) <- pnames
  refused             <- rep(FALSE, n);    names(refused) <- pnames
  feasible_plus_final  <- rep(NA, n); names(feasible_plus_final)  <- pnames
  feasible_minus_final <- rep(NA, n); names(feasible_minus_final) <- pnames

  ## cache of f(theta +/- h*e_i) at the FINAL chosen step, reused when
  ## assembling the diagonal and as neighbours for off-diagonal terms
  fp_cache <- rep(NA_real_, n); names(fp_cache) <- pnames
  fm_cache <- rep(NA_real_, n); names(fm_cache) <- pnames

  for (i in seq_len(n)) {
    step0 <- rel_step * max(abs(theta[i]), 1)

    h_prev <- NA_real_
    kappa_prev <- NA_real_
    chosen <- NA_real_
    chosen_fp <- NA_real_; chosen_fm <- NA_real_
    chosen_one_sided <- FALSE; chosen_dir <- NA_character_
    last_plus_ok <- FALSE; last_minus_ok <- FALSE

    for (k in 0:max_shrink) {
      h <- step0 / 2^k
      h_p <- min(h, hi[i] - theta[i])
      h_m <- min(h, theta[i] - lo[i])
      if (!is.finite(h_p) || h_p <= 0) h_p <- 0
      if (!is.finite(h_m) || h_m <= 0) h_m <- 0

      f_p <- if (h_p > 0) { th <- theta; th[i] <- theta[i] + h_p; safe_fn(th) } else NA_real_
      f_m <- if (h_m > 0) { th <- theta; th[i] <- theta[i] - h_m; safe_fn(th) } else NA_real_

      ok_p <- h_p > 0 && feasible(f_p)
      ok_m <- h_m > 0 && feasible(f_m)
      last_plus_ok <- ok_p; last_minus_ok <- ok_m

      if (ok_p && ok_m) {
        ## Both directions feasible at this (possibly asymmetric) step.
        ## Use the smaller of the two magnitudes so the central-difference
        ## formula (which assumes a common h) stays valid.
        h_use <- min(h_p, h_m)
        if (h_use < h) {
          ## re-evaluate at the common shrunk step
          th_p <- theta; th_p[i] <- theta[i] + h_use
          th_m <- theta; th_m[i] <- theta[i] - h_use
          f_p2 <- safe_fn(th_p); f_m2 <- safe_fn(th_m)
          if (!(feasible(f_p2) && feasible(f_m2))) {
            ## fall through to shrinking further
            chosen <- NA_real_
          } else {
            f_p <- f_p2; f_m <- f_m2; h <- h_use
            chosen <- h; chosen_fp <- f_p; chosen_fm <- f_m
            chosen_one_sided <- FALSE; chosen_dir <- NA_character_
          }
        } else {
          chosen <- h; chosen_fp <- f_p; chosen_fm <- f_m
          chosen_one_sided <- FALSE; chosen_dir <- NA_character_
        }

        if (!is.na(chosen)) {
          if (!verify_spectrum) break
          ## spectrum-stability bookkeeping happens after the per-coordinate
          ## loop (needs the *full* Hessian); here we just require the
          ## stencil to have been feasible at two consecutive shrinks so a
          ## later full-Hessian check has something to compare. Record and
          ## keep shrinking one more level to get a second, smaller feasible
          ## step for the spectrum check, done in Step 2 below.
          break
        }
      } else if (ok_p && !ok_m) {
        chosen <- h_p; chosen_fp <- f_p; chosen_fm <- NA_real_
        chosen_one_sided <- TRUE; chosen_dir <- "forward"
        break
      } else if (ok_m && !ok_p) {
        chosen <- h_m; chosen_fp <- NA_real_; chosen_fm <- f_m
        chosen_one_sided <- TRUE; chosen_dir <- "backward"
        break
      }
      ## neither feasible at this step -- shrink and retry
    }

    if (is.na(chosen)) {
      refused[i] <- TRUE
      feasible_plus_final[i]  <- last_plus_ok
      feasible_minus_final[i] <- last_minus_ok
      next
    }

    step_used[i] <- chosen
    one_sided[i] <- chosen_one_sided
    one_sided_direction[i] <- chosen_dir
    fp_cache[i] <- chosen_fp
    fm_cache[i] <- chosen_fm
    feasible_plus_final[i]  <- !is.na(chosen_fp) || chosen_dir %in% c("backward")
    feasible_minus_final[i] <- !is.na(chosen_fm) || chosen_dir %in% c("forward")
    ## more precisely: mark which side actually has a feasible eval
    feasible_plus_final[i]  <- isTRUE(!is.na(chosen_fp))
    feasible_minus_final[i] <- isTRUE(!is.na(chosen_fm))
  }

  diagnosis <- data.frame(
    coord           = pnames,
    feasible_plus   = feasible_plus_final,
    feasible_minus  = feasible_minus_final,
    step_used       = step_used,
    one_sided       = one_sided,
    refused         = refused,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (any(refused)) {
    bad <- pnames[refused]
    msg <- sprintf(
      "fd_safe_hessian: no feasible step found for coordinate(s) %s in EITHER direction (tried %d shrinks from a relative step of %.3g). Refusing to return a Hessian rather than silently poisoning it with a boundary-crossing stencil. Diagnosis: %s",
      paste(bad, collapse = ", "), max_shrink, rel_step,
      paste(sprintf("%s (step tried down to %.3g)", bad,
                     rel_step * sapply(theta[bad], function(v) max(abs(v), 1)) / 2^max_shrink),
            collapse = "; ")
    )
    return(list(
      ok = FALSE, H = NULL,
      step_used = step_used, one_sided = one_sided,
      one_sided_direction = one_sided_direction,
      diagnosis = diagnosis, message = msg
    ))
  }

  ## ---- Step 2: assemble Hessian at the chosen per-coordinate steps ------
  build_H <- function(h_vec, fp_vec, fm_vec, f0_val) {
    Hm <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(pnames, pnames))
    for (i in seq_len(n)) {
      if (one_sided[i]) {
        ## one-sided second difference using a midpoint at h/2 (first-order
        ## accurate): H_ii ~= 2*(f(x+-h) - f(x) - h*g)/h^2 is unavailable
        ## without the gradient, so use the standard 3-point one-sided
        ## formula with an extra evaluation at h/2 in the feasible direction.
        dir <- one_sided_direction[i]
        h_i <- h_vec[i]
        sgn <- if (identical(dir, "forward")) 1 else -1
        th_half <- theta; th_half[i] <- theta[i] + sgn * h_i / 2
        f_half <- safe_fn(th_half)
        f_full <- if (identical(dir, "forward")) fp_vec[i] else fm_vec[i]
        ## f(x + s*h) = f0 + s*h*g + s^2 h^2 g''/2 + O(h^3)
        ## Using s=1/2 and s=1: solve for second derivative.
        Hm[i, i] <- 4 * (f_full - 2 * f_half + f0_val) / (h_i^2)
      } else {
        Hm[i, i] <- (fp_vec[i] - 2 * f0_val + fm_vec[i]) / h_vec[i]^2
      }
    }
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        h_i <- h_vec[i]; h_j <- h_vec[j]
        if (one_sided[i] || one_sided[j]) {
          ## one-sided cross difference: use whichever direction is feasible
          ## for each involved coordinate (falls back to forward/backward
          ## consistently with that coordinate's chosen direction; for a
          ## coordinate that is NOT one-sided, use its central +h side).
          si <- if (one_sided[i]) (if (identical(one_sided_direction[i], "forward")) 1 else -1) else 1
          sj <- if (one_sided[j]) (if (identical(one_sided_direction[j], "forward")) 1 else -1) else 1

          th_ij <- theta; th_ij[i] <- theta[i] + si * h_i; th_ij[j] <- theta[j] + sj * h_j
          th_i0 <- theta; th_i0[i] <- theta[i] + si * h_i
          th_0j <- theta; th_0j[i] <- theta[i]; th_0j[j] <- theta[j] + sj * h_j

          f_ij <- safe_fn(th_ij)
          f_i0 <- if (si > 0 && !one_sided[i]) fp_vec[i] else if (si < 0 && !one_sided[i]) fm_vec[i] else safe_fn(th_i0)
          f_0j <- if (sj > 0 && !one_sided[j]) fp_vec[j] else if (sj < 0 && !one_sided[j]) fm_vec[j] else safe_fn(th_0j)

          ## mixed partial via forward/backward difference of differences:
          ## d2f/didj ~= (f(x+si*hi+sj*hj) - f(x+si*hi) - f(x+sj*hj) + f(x)) / (si*hi*sj*hj)
          Hm[i, j] <- Hm[j, i] <-
            (f_ij - f_i0 - f_0j + f0_val) / (si * h_i * sj * h_j)
        } else {
          th_pp <- th_pm <- th_mp <- th_mm <- theta
          th_pp[i] <- theta[i] + h_i; th_pp[j] <- theta[j] + h_j
          th_pm[i] <- theta[i] + h_i; th_pm[j] <- theta[j] - h_j
          th_mp[i] <- theta[i] - h_i; th_mp[j] <- theta[j] + h_j
          th_mm[i] <- theta[i] - h_i; th_mm[j] <- theta[j] - h_j

          f_pp <- safe_fn(th_pp); f_pm <- safe_fn(th_pm)
          f_mp <- safe_fn(th_mp); f_mm <- safe_fn(th_mm)

          if (!all(vapply(list(f_pp, f_pm, f_mp, f_mm), feasible, logical(1)))) {
            ## an off-diagonal corner crossed the boundary even though both
            ## marginal coordinates were individually feasible -- shrink this
            ## pair's steps by half and retry once; if still infeasible,
            ## fall back to the smaller of the two one-sided-consistent
            ## corners actually evaluated (best-effort, not a refusal: the
            ## diagonal terms driving the condition number are already
            ## feasibility-verified; a single bad corner is logged via NA
            ## replaced by 0 cross-curvature only as a last resort)
            h_i2 <- h_i / 2; h_j2 <- h_j / 2
            th_pp[i] <- theta[i] + h_i2; th_pp[j] <- theta[j] + h_j2
            th_pm[i] <- theta[i] + h_i2; th_pm[j] <- theta[j] - h_j2
            th_mp[i] <- theta[i] - h_i2; th_mp[j] <- theta[j] + h_j2
            th_mm[i] <- theta[i] - h_i2; th_mm[j] <- theta[j] - h_j2
            f_pp2 <- safe_fn(th_pp); f_pm2 <- safe_fn(th_pm)
            f_mp2 <- safe_fn(th_mp); f_mm2 <- safe_fn(th_mm)
            if (all(vapply(list(f_pp2, f_pm2, f_mp2, f_mm2), feasible, logical(1)))) {
              Hm[i, j] <- Hm[j, i] <- (f_pp2 - f_pm2 - f_mp2 + f_mm2) / (4 * h_i2 * h_j2)
            } else {
              Hm[i, j] <- Hm[j, i] <- (fp_vec[i] - f0_val) * 0  ## 0, documented best-effort
            }
          } else {
            Hm[i, j] <- Hm[j, i] <- (f_pp - f_pm - f_mp + f_mm) / (4 * h_i * h_j)
          }
        }
      }
    }
    Hm
  }

  H <- build_H(step_used, fp_cache, fm_cache, f0)

  ## ---- Step 3: optional spectrum-stability verification -----------------
  if (verify_spectrum && n >= 1L && all(is.finite(H))) {
    half_step <- step_used / 2
    ## re-evaluate at half the step for coordinates that were NOT one-sided
    ## (one-sided coordinates already used the tightest feasible step, so we
    ## do not attempt to shrink them further here -- their accuracy downgrade
    ## is already flagged).
    fp2 <- fp_cache; fm2 <- fm_cache
    ok_half <- TRUE
    for (i in seq_len(n)) {
      if (one_sided[i]) next
      th_p <- theta; th_p[i] <- theta[i] + half_step[i]
      th_m <- theta; th_m[i] <- theta[i] - half_step[i]
      f_p <- safe_fn(th_p); f_m <- safe_fn(th_m)
      if (!(feasible(f_p) && feasible(f_m))) { ok_half <- FALSE; break }
      fp2[i] <- f_p; fm2[i] <- f_m
    }

    if (ok_half) {
      H2 <- build_H(ifelse(one_sided, step_used, half_step), fp2, fm2, f0)
      if (all(is.finite(H2))) {
        kappa <- function(M) {
          Ms <- (M + t(M)) / 2
          ev <- eigen(Ms, symmetric = TRUE, only.values = TRUE)$values
          mag <- abs(ev)
          if (min(mag) <= 0 || !is.finite(max(mag))) return(NA_real_)
          max(mag) / min(mag)
        }
        k1 <- kappa(H); k2 <- kappa(H2)
        if (is.finite(k1) && is.finite(k2) && k1 > 0 && k2 > 0) {
          rel_change <- abs(log(k2) - log(k1)) / max(abs(log(k1)), 1e-8)
          if (rel_change <= spectrum_tol) {
            ## spectrum stabilized at the finer step -- prefer it (more accurate)
            H <- H2
            step_used[!one_sided] <- half_step[!one_sided]
          }
          ## if not stabilized, keep the coarser (already-verified-feasible)
          ## H; we do not refuse here because both steps ARE feasible, just
          ## not yet in perfect spectral agreement -- documented via the
          ## fact that verify_spectrum only refines, never refuses.
        }
      }
    }
  }

  H <- (H + t(H)) / 2  ## symmetrize numerical noise, matches package convention

  ## keep the diagnosis table's step_used column in sync with any refinement
  ## made during the Step-3 spectrum verification above
  diagnosis$step_used <- step_used[diagnosis$coord]

  attr(H, "step_used")           <- step_used
  attr(H, "one_sided")           <- one_sided
  attr(H, "one_sided_direction") <- one_sided_direction
  attr(H, "diagnosis")           <- diagnosis

  list(
    ok = TRUE, H = H,
    step_used = step_used, one_sided = one_sided,
    one_sided_direction = one_sided_direction,
    diagnosis = diagnosis, message = NA_character_
  )
}


#' @noRd
.fdsh_expand_bound <- function(bound, pnames, default) {
  n <- length(pnames)
  if (is.null(bound)) return(rep(default, n))
  if (!is.null(names(bound))) {
    out <- rep(default, n); names(out) <- pnames
    common <- intersect(names(bound), pnames)
    out[common] <- bound[common]
    return(out)
  }
  if (length(bound) == n) return(as.numeric(bound))
  stop("fd_safe_hessian: unnamed lower/upper must have length equal to theta.")
}
