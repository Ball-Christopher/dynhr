## R/diag-pre-d2-spectral.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R; updated Phase A.
##
## D2 spectral / D23 identification check (Qu-Tkachenko 2012)
## --------------------------------------------------------------------------

## NOTE (2026-05-31): the `d2_spectral_identification` stub was removed. It was
## an admitted re-label of D23 ("D2 is a re-labelled D23") that only returned an
## "identical to D23 ... Skipped" placeholder, duplicating the diagnostic.
## `d23_spectral_identification` below is the real Qu-Tkachenko spectral check.

#' D23. Spectral identification check (Qu & Tkachenko 2012)
#'
#' Full implementation. Computes the Gram matrix of spectral density
#' derivatives from the first-order state-space representation.
#'
#' @inheritParams d23_spectral_identification
#' @return dynhr_diagnostic list
#' @noRd
d23_spectral_identification <- function(dr             = NULL,
                                        model_solve_fn = NULL,
                                        theta          = NULL,
                                        param_names    = NULL,
                                        Sigma_e        = NULL,
                                        n_freq         = 256L,
                                        eps            = 1e-5,
                                        ...) {
  if (is.null(dr) || is.null(theta)) {
    # Called without enough context — return informational
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = "D23 Spectral identification: provide dr, theta, and model_solve_fn."
    ))
  }
  # Delegate to the implementation in diag-pre-d23-spectral.R
  d23_spectral_identification_impl(
    dr             = dr,
    model_solve_fn = model_solve_fn,
    theta          = theta,
    param_names    = param_names,
    Sigma_e        = Sigma_e,
    n_freq         = n_freq,
    eps            = eps
  )
}
