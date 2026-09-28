## R/diag-pre-d2-spectral.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R; updated Phase A.
##
## D23 identification check (Qu-Tkachenko 2012): entry point only.
## The implementation is d23_spectral_identification_impl() in
## diag-pre-d23-spectral.R (the two files are wrapper + body, not duplicates).
## --------------------------------------------------------------------------

## NOTE (2026-05-31): the `d2_spectral_identification` stub was removed. It was
## an admitted re-label of D23 ("D2 is a re-labelled D23") that only returned an
## "identical to D23 ... Skipped" placeholder, duplicating the diagnostic.

#' D23. Spectral identification check (Qu & Tkachenko 2012)
#'
#' Entry point used by \code{run_all_diagnostics()}; see
#' \code{d23_spectral_identification_impl()} for the method.
#'
#' @inheritParams d23_spectral_identification_impl
#' @param ... Ignored.
#' @return dynhr_diagnostic list
#' @noRd
d23_spectral_identification <- function(dr             = NULL,
                                        model_solve_fn = NULL,
                                        theta          = NULL,
                                        param_names    = NULL,
                                        Sigma_e        = NULL,
                                        n_freq         = 256L,
                                        eps            = 1e-5,
                                        tol_rank       = NULL,
                                        weak_rel       = 1e-3,
                                        meta           = NULL,
                                        ...) {
  if (is.null(dr)) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = "D23 Spectral identification: provide dr, theta, and model_solve_fn."
    ))
  }
  d23_spectral_identification_impl(
    dr             = dr,
    model_solve_fn = model_solve_fn,
    theta          = theta,
    param_names    = param_names,
    Sigma_e        = Sigma_e,
    n_freq         = n_freq,
    eps            = eps,
    tol_rank       = tol_rank,
    weak_rel       = weak_rel,
    meta           = meta
  )
}
