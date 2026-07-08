## zzz.R
## --------------------------------------------------------------------------
## Package load / attach hooks.
## --------------------------------------------------------------------------

.onAttach <- function(libname, pkgname) {
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
  ## Reserved for future use (registering S3 methods, options defaults, etc).
  ## Phase 0: no-op.
  invisible(NULL)
}
