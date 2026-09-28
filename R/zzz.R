## zzz.R
## --------------------------------------------------------------------------
## Package load / attach hooks.
## --------------------------------------------------------------------------

.onAttach <- function(libname, pkgname) {
  # Suppress startup banner during testthat runs
  # (testthat sets TESTTHAT=true; no need to load testthat on every attach)
  if (identical(Sys.getenv("TESTTHAT"), "true")) {
    return(invisible(NULL))
  }
  ## Keep this banner minimal: a hardcoded feature list drifts out of sync
  ## with the namespace and reads as noise in downstream (e.g. paper) sessions.
  n_exports <- length(getNamespaceExports(pkgname))
  ver <- tryCatch(as.character(utils::packageVersion(pkgname)),
                  error = function(e) "")
  packageStartupMessage(
    "dynhr ", ver, " (", n_exports, " exported functions)"
  )
}

.onLoad <- function(libname, pkgname) {
  ## Conditional S3 registration for the Suggests packages posterior and coda
  ## (as_draws*(), as.mcmc*(); R/interop.R): registered when -- or, if
  ## already loaded, as soon as -- their namespaces load, so dynhr need not
  ## Import them.
  .register_interop_methods()
  invisible(NULL)
}
