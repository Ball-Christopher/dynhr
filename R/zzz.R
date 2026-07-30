## zzz.R
## --------------------------------------------------------------------------
## Package load / attach hooks.
## --------------------------------------------------------------------------

## RELEASE-ONLY DIVERGENCE (re-applied at every cut; see NEWS/README notes).
## Two things present on master are stripped here:
##
##  * the testthat is_testing() banner guard. testthat is not in this package's
##    Suggests -- the test harness is not shipped -- so referencing it is an
##    undeclared dependency and R CMD check raises it as a WARNING, not a NOTE.
##  * the optional plot-tools sourcing. inst/extdata/dynhr_plot_tools_v3.R is
##    not part of the curated release, so the block is dead code here.
##
## A wholesale sync from master reintroduces both. If R CMD check on the release
## tarball reports "'::' or ':::' import not declared from: 'testthat'", this
## edit has been lost.
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
