## R/native-build.R
## --------------------------------------------------------------------------
## Opt-in native (CPU-tuned) build of dynhr into a user library, without
## touching package files, ~/.R/Makevars or R's own Makeconf.
##
## The mechanism is R's R_MAKEVARS_USER: the install process reads the Makevars
## file that variable names IN PLACE of ~/.R/Makevars, so a temporary file
## carries the flags for exactly one `R CMD INSTALL` and disappears with it.
##
## WHICH VARIABLES. dynhr's src/Makevars declares CXX_STD = CXX20, and for a
## CXX20 package R compiles with $(CXX20) $(CXX20STD) ... $(CXX20FLAGS), not
## CXXFLAGS. Setting only CXXFLAGS would therefore leave the C++ at R's default
## -O2 and report a "native" build that is not one. The temporary Makevars sets
## the flag variable for every C++ standard R knows about (so it keeps working
## if the package standard changes), plus CFLAGS and FFLAGS/FCFLAGS for the C
## and Fortran sources. This REPLACES R's default flag set (-g -Wall -O2 ...)
## rather than appending to it, so that a later -O level wins unambiguously.
## --------------------------------------------------------------------------


## Flags that change floating-point semantics rather than speed alone. Fast-math
## reassociates sums and assumes finite values; the Kalman filter, the
## Lyapunov/QZ solvers and the NaN-aware missing-data paths rely on IEEE
## behaviour, so a build with these is refused rather than warned about.
.native_forbidden_flags <- c("-ffast-math", "-Ofast", "-funsafe-math-optimizations",
                             "-ffinite-math-only", "-fassociative-math")

## Forbidden flags present in `flags` (a character vector, split on whitespace).
.native_find_forbidden <- function(flags) {
  tok <- unlist(strsplit(paste(flags, collapse = " "), "[[:space:]]+"))
  tok <- tok[nzchar(tok)]
  unique(tok[tok %in% .native_forbidden_flags])
}

## Makevars text carrying `flags` for every compiler variable a dynhr
## installation can use (see the header for why CXX20FLAGS matters).
.native_makevars <- function(flags) {
  flags <- paste(flags, collapse = " ")
  vars <- c("CFLAGS", "CXXFLAGS", "CXX11FLAGS", "CXX14FLAGS", "CXX17FLAGS",
            "CXX20FLAGS", "CXX23FLAGS", "FFLAGS", "FCFLAGS")
  paste0(c("## Temporary Makevars written by dynhr_install_native().",
           paste0(vars, " = ", flags)), collapse = "\n")
}

## Stage a source directory as a tarball-free install tree in a temp dir:
## top-level dot entries (.git and similar), scratch output and every compiled
## object are left out, so the install compiles from source with the new flags
## (a stale -O0 object from a development load would otherwise be reused) and
## the caller's tree is not modified.
.native_stage_source <- function(dir) {
  stage <- tempfile("dynhr-native-src-")
  dir.create(stage)
  top <- list.files(dir, all.files = FALSE, no.. = TRUE)
  top <- setdiff(top, c("scratch", "tests_out"))
  ok <- file.copy(file.path(dir, top), stage, recursive = TRUE)
  if (!all(ok))
    .dynhr_abort("dynhr_install_native: could not copy the source directory ",
                 "to a temporary location.",
                 class = "dynhr_error_native_build")
  src <- file.path(stage, "src")
  if (dir.exists(src))
    unlink(list.files(src, pattern = "[.](o|so|dll|dylib)$", full.names = TRUE))
  stage
}

## What `dynhr_install_native()` runs: R executable, arguments, environment.
.native_install_command <- function(pkg, lib, makevars_file) {
  rbin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
  list(command = rbin,
       args = c("CMD", "INSTALL", "--preclean", "--no-test-load",
                paste0("--library=", shQuote(lib)), shQuote(pkg)),
       env = paste0("R_MAKEVARS_USER=", makevars_file))
}

.native_fingerprint_logpost <- -874.6253133934

## Script run in a fresh R process against the freshly installed library.
.native_verify_script <- function() {
  c("a <- commandArgs(trailingOnly = TRUE)",
    "library(dynhr, lib.loc = a[1L])",
    "p <- dynhr:::.bench_problem()",
    "saveRDS(list(path = find.package('dynhr'),",
    "             info = dynhr_system_info(),",
    "             logpost = p$logpost0), a[2L])")
}


#' Install dynhr with CPU-specific compiler flags, without administrator rights
#'
#' Builds dynhr from source with your own compiler flags (by default
#' \code{-O3 -march=native}) into a user library and checks the result in a
#' fresh R session. It needs a working C++20 toolchain (Rtools on Windows,
#' Xcode command-line tools on macOS) but no administrator rights and no
#' edits to R's configuration: the flags reach the compiler through a
#' temporary Makevars file named by the \code{R_MAKEVARS_USER} environment
#' variable for that one installation. \file{~/.R/Makevars} and R's own
#' \file{Makeconf} are never read or changed, and the temporary file is
#' deleted afterwards.
#'
#' The flags replace R's default flag set (\code{-g -Wall -O2} ...) for the C,
#' C++ (every standard, including the C++20 dynhr requires) and Fortran
#' compilers. Compiled dynhr code then reports what was actually built in
#' \code{\link{dynhr_system_info}} (columns \code{dynhr_simd},
#' \code{dynhr_fma}, \code{dynhr_fast_math}, ...).
#'
#' Flags that change floating-point semantics (\code{-ffast-math},
#' \code{-Ofast}, \code{-funsafe-math-optimizations},
#' \code{-ffinite-math-only}, \code{-fassociative-math}) are refused with an
#' error of class \code{dynhr_error_native_flags}: the filters and solvers
#' rely on IEEE arithmetic and results from such a build are not supported.
#'
#' @section Portability and accuracy:
#' A \code{-march=native} library uses instructions of the CPU that built it
#' and is NOT portable: loading it on an older or different CPU can crash R
#' with an illegal-instruction error. Build on the machine that will run it.
#' Fused multiply-add contraction, which \code{-march=native} usually enables,
#' changes the rounding of intermediate results (as does the more aggressive
#' vectorisation of \code{-O3}), so likelihoods move at round-off level: the
#' benchmark fingerprint moved by a relative \code{9e-15} on Apple silicon
#' (where fused multiply-add is always on, so the move comes from \code{-O3}
#' alone), far below any estimation tolerance but not bit-identical to a
#' portable build; seeded chains can therefore diverge after many iterations.
#' Do not mix native and portable builds when comparing seeded runs.
#'
#' @section Going back:
#' Reinstall dynhr the usual way (\code{install.packages()} or
#' \code{R CMD INSTALL} without this function) into the same library.
#'
#' @section Verification:
#' Unless \code{verify = FALSE}, a fresh \code{Rscript} process loads the
#' installed library, and the function checks that the loaded copy lives in
#' \code{lib}, records \code{\link{dynhr_system_info}} and compares the
#' log posterior of the fixed Smets-Wouters benchmark problem (the
#' \code{\link{dynhr_benchmark}} fingerprint, -874.6253133934) with its
#' reference. A relative difference above \code{fingerprint_tol} gives a
#' warning of class \code{dynhr_warning_native_fingerprint}.
#'
#' @param pkg Path to a dynhr source tarball (\file{dynhr_*.tar.gz}) or to a
#'   dynhr source directory (it is copied to a temporary location first, so
#'   the directory is not modified).
#' @param lib Library to install into. Default: your first library path.
#' @param flags Compiler flags, one string. Default \code{"-O3 -march=native"}.
#' @param dry_run If \code{TRUE}, return the Makevars text and the command that
#'   would run, without building anything.
#' @param verify Run the fresh-process verification after installing.
#' @param fingerprint_tol Relative tolerance for the fingerprint check.
#'
#' @return An object of class \code{dynhr_native_build}: a list with
#'   \code{makevars} (text), \code{command}, \code{flags}, \code{dry_run},
#'   and, after a real install, \code{status}, \code{lib} and \code{log}
#'   (compiler output) plus, when verified, \code{info} (the one-row
#'   \code{\link{dynhr_system_info}} of the new build), \code{logpost},
#'   \code{fingerprint_ref}, \code{fingerprint_reldiff} and
#'   \code{fingerprint_ok}.
#'
#' @examples
#' \dontrun{
#' # from a source tarball, into your user library, never touching R itself
#' dynhr_install_native("dynhr_0.9.5.tar.gz")
#' }
#' # what would be run, without building
#' dynhr_install_native(tempdir(), flags = "-O3 -march=native", dry_run = TRUE)
#' @seealso \code{\link{dynhr_system_info}}, \code{\link{dynhr_benchmark}}
#' @export
dynhr_install_native <- function(pkg, lib = .libPaths()[1L],
                                 flags = "-O3 -march=native",
                                 dry_run = FALSE, verify = TRUE,
                                 fingerprint_tol = 1e-8) {
  if (!is.character(flags) || length(flags) < 1L || anyNA(flags) ||
      !nzchar(trimws(paste(flags, collapse = " "))))
    .dynhr_abort("dynhr_install_native: `flags` must be a non-empty string.",
                 class = "dynhr_error_native_flags")
  bad <- .native_find_forbidden(flags)
  if (length(bad))
    .dynhr_abort("dynhr_install_native: refusing flag(s) ",
                 paste(bad, collapse = ", "),
                 ". They change floating-point semantics (reassociation, ",
                 "no NaN/Inf), which the Kalman filter and solvers do not ",
                 "support; use -O3 and -march=native instead.",
                 class = "dynhr_error_native_flags")
  if (!is.character(pkg) || length(pkg) != 1L || is.na(pkg) || !nzchar(pkg))
    .dynhr_abort("dynhr_install_native: `pkg` must be the path of a dynhr ",
                 "source tarball or source directory.",
                 class = "dynhr_error_native_build")
  if (!is.character(lib) || length(lib) != 1L || is.na(lib) || !nzchar(lib))
    .dynhr_abort("dynhr_install_native: `lib` must be a single directory.",
                 class = "dynhr_error_native_build")
  flags <- paste(trimws(flags), collapse = " ")
  makevars <- .native_makevars(flags)

  mk_file <- file.path(tempdir(), "dynhr-native-Makevars")
  cmd <- .native_install_command(pkg, lib, mk_file)
  if (isTRUE(dry_run))
    return(structure(list(makevars = makevars, command = cmd, flags = flags,
                          dry_run = TRUE), class = "dynhr_native_build"))

  if (!dir.exists(lib))
    .dynhr_abort("dynhr_install_native: library '", lib, "' does not exist.",
                 class = "dynhr_error_native_build")
  if (file.access(lib, 2L) != 0L)
    .dynhr_abort("dynhr_install_native: library '", lib, "' is not writable; ",
                 "pass a user library via `lib`.",
                 class = "dynhr_error_native_build")
  if (!file.exists(pkg))
    .dynhr_abort("dynhr_install_native: '", pkg, "' does not exist.",
                 class = "dynhr_error_native_build")
  staged <- NULL
  if (dir.exists(pkg)) {
    desc <- file.path(pkg, "DESCRIPTION")
    if (!file.exists(desc) ||
        !identical(unname(read.dcf(desc, fields = "Package")[1L, 1L]), "dynhr"))
      .dynhr_abort("dynhr_install_native: '", pkg, "' is not a dynhr source ",
                   "directory (no DESCRIPTION with Package: dynhr).",
                   class = "dynhr_error_native_build")
    staged <- .native_stage_source(normalizePath(pkg))
    pkg <- staged
  } else if (!grepl("^dynhr_.*[.]tar[.]gz$", basename(pkg))) {
    .dynhr_abort("dynhr_install_native: '", basename(pkg), "' is not a ",
                 "dynhr source tarball (dynhr_*.tar.gz).",
                 class = "dynhr_error_native_build")
  }
  lib <- normalizePath(lib, winslash = "/")
  mk_file <- normalizePath(file.path(tempdir(), "dynhr-native-Makevars"),
                           winslash = "/", mustWork = FALSE)
  writeLines(makevars, mk_file)
  on.exit(unlink(c(mk_file, staged), recursive = TRUE), add = TRUE)
  cmd <- .native_install_command(normalizePath(pkg, winslash = "/"), lib, mk_file)

  .dynhr_inform("dynhr_install_native: building with flags '", flags,
                "' into ", lib, " ...")
  out <- suppressWarnings(system2(cmd$command, cmd$args, env = cmd$env,
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status") %||% 0L
  if (!identical(as.integer(status), 0L))
    .dynhr_abort("dynhr_install_native: R CMD INSTALL failed (status ",
                 status, "). Last lines of output:\n",
                 paste(utils::tail(out, 25L), collapse = "\n"),
                 class = "dynhr_error_native_build")
  res <- list(makevars = makevars, command = cmd, flags = flags,
              dry_run = FALSE, status = 0L, lib = lib, log = out)

  if (isTRUE(verify)) {
    script <- tempfile("dynhr-native-verify-", fileext = ".R")
    rds <- tempfile("dynhr-native-verify-", fileext = ".rds")
    writeLines(.native_verify_script(), script)
    on.exit(unlink(c(script, rds)), add = TRUE)
    rscript <- file.path(R.home("bin"),
                         if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    vout <- suppressWarnings(system2(rscript,
                                     c(shQuote(script), shQuote(lib), shQuote(rds)),
                                     stdout = TRUE, stderr = TRUE))
    vstatus <- attr(vout, "status") %||% 0L
    if (!identical(as.integer(vstatus), 0L) || !file.exists(rds))
      .dynhr_abort("dynhr_install_native: the installed build failed to load ",
                   "in a fresh R session (status ", vstatus, "). If the CPU ",
                   "running R differs from the CPU that built it, rebuild ",
                   "there or use flags without -march=native. Output:\n",
                   paste(utils::tail(vout, 15L), collapse = "\n"),
                   class = "dynhr_error_native_build")
    chk <- readRDS(rds)
    loaded <- normalizePath(chk$path, winslash = "/")
    if (!startsWith(loaded, lib))
      .dynhr_abort("dynhr_install_native: the fresh session loaded dynhr from '",
                   loaded, "', not from '", lib, "'.",
                   class = "dynhr_error_native_build")
    ref <- .native_fingerprint_logpost
    rel <- abs(chk$logpost - ref) / abs(ref)
    res$info <- chk$info
    res$logpost <- chk$logpost
    res$fingerprint_ref <- ref
    res$fingerprint_reldiff <- rel
    res$fingerprint_ok <- is.finite(rel) && rel <= fingerprint_tol
    if (!res$fingerprint_ok)
      .dynhr_warn("dynhr_install_native: the benchmark fingerprint of the new ",
                  "build is ", format(chk$logpost, digits = 13), " against the ",
                  "reference ", format(ref, digits = 13), " (relative ",
                  "difference ", format(rel, digits = 3), ", tolerance ",
                  format(fingerprint_tol), ").",
                  class = "dynhr_warning_native_fingerprint")
  }
  structure(res, class = "dynhr_native_build")
}


#' @param x A \code{dynhr_native_build} object.
#' @param ... Ignored.
#' @rdname dynhr_install_native
#' @export
print.dynhr_native_build <- function(x, ...) {
  if (isTRUE(x$dry_run)) {
    cat("dynhr native build (dry run, nothing installed)\n\n")
    cat("Makevars (via R_MAKEVARS_USER):\n", x$makevars, "\n\n", sep = "")
    cat("Command:\n  ", x$command$env, " ", x$command$command, " ",
        paste(x$command$args, collapse = " "), "\n", sep = "")
    return(invisible(x))
  }
  cat("dynhr native build installed into ", x$lib, "\n", sep = "")
  cat("  flags: ", x$flags, "\n", sep = "")
  if (!is.null(x$info)) {
    i <- x$info
    cat(sprintf("  compiler %s | simd %s | optimised %s | fast-math %s\n",
                i$dynhr_compiler, i$dynhr_simd, i$dynhr_optimize,
                i$dynhr_fast_math))
    cat(sprintf("  fingerprint logpost %.10f (reference %.10f, relative diff %.2e) %s\n",
                x$logpost, x$fingerprint_ref, x$fingerprint_reldiff,
                if (isTRUE(x$fingerprint_ok)) "OK" else "** MISMATCH **"))
  }
  cat("  To go back, reinstall dynhr normally into the same library.\n")
  invisible(x)
}
