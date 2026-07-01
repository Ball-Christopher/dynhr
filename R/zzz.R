## zzz.R
## --------------------------------------------------------------------------
## Package load / attach hooks.
## --------------------------------------------------------------------------

.onAttach <- function(libname, pkgname) {
  n_exports <- length(getNamespaceExports(pkgname))
  ver <- tryCatch(as.character(utils::packageVersion(pkgname)),
                  error = function(e) "")
  packageStartupMessage(
    "dynhr ", ver, " (", n_exports, " exported functions)\n",
    "  - Solver: parse_mod, compile_model, solve_steady, solve_perturbation,\n",
    "    stoch_simul, compute_irfs, compute_moments, simulate_model\n",
    "  - Filtering: kalman_filter, kalman_smoother, build_dsge_state_space,\n",
    "    historical_decomposition\n",
    "  - Optimal policy: ramsey_model, ramsey_nn1, osr, discretionary_policy\n",
    "  - Diagnostics: run_diagnostics, diag_expectations, write_report"
  )
}

.onLoad <- function(libname, pkgname) {
  ## Reserved for future use (registering S3 methods, options defaults, etc).
  ## Phase 0: no-op.
  invisible(NULL)
}
