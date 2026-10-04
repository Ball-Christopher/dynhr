## R/benchmark.R
## --------------------------------------------------------------------------
## Cross-MACHINE benchmark: run a fixed, realistic estimation workload
## (Smets-Wouters 2007, 36 estimated parameters, 7 observables, 160 quarters)
## through RWMH at a sweep of core counts, and record enough system detail that
## a result from one machine can be compared with a result from another.
##
## WHY THE MEASURES ARE THROUGHPUTS, NOT TIMES. The comparison this exists for
## ("how does the new laptop compare to the old one?") has to survive different
## draw counts: a machine that is too slow to finish 100k draws in the budgeted
## time still has to be comparable with one that is. So every reported quantity
## is per-draw or per-second -- `draws_per_sec`, `us_per_draw`, `speedup`,
## `efficiency` -- and none of them is an elapsed time. Two runs with different
## `n_draws` are directly comparable; two runs with different MODELS are not,
## which is why the model is fixed and fingerprinted.
##
## WHY CHAINS, NOT BLAS THREADS, ARE THE CORE AXIS. A single RWMH chain is
## inherently serial, and this workload is R-overhead-bound on small matrices
## (see the perf-bottleneck-profile note): adding BLAS threads to one chain buys
## essentially nothing. The parallelism that matters in practice is independent
## chains, so the sweep holds the work PER CHAIN fixed and adds chains. Perfect
## scaling therefore shows up as constant elapsed time and linear throughput.
## --------------------------------------------------------------------------


## Shell out for one system fact, returning NA_character_ rather than failing:
## a benchmark must never die because a machine lacks `sysctl`.
.bench_sh <- function(cmd) {
  out <- tryCatch(
    suppressWarnings(system(cmd, intern = TRUE, ignore.stderr = TRUE)),
    error = function(e) character(0))
  if (!length(out) || !nzchar(out[1L])) NA_character_ else trimws(out[1L])
}


.bench_cpu_model <- function() {
  switch(Sys.info()[["sysname"]],
    Darwin  = .bench_sh("sysctl -n machdep.cpu.brand_string"),
    Linux   = {
      v <- .bench_sh("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2")
      if (is.na(v)) .bench_sh("lscpu | grep -m1 'Model name' | cut -d: -f2") else v
    },
    Windows = .bench_sh("wmic cpu get name /value | findstr Name="),
    NA_character_)
}


.bench_ram_gb <- function() {
  b <- switch(Sys.info()[["sysname"]],
    Darwin  = suppressWarnings(as.numeric(.bench_sh("sysctl -n hw.memsize"))),
    Linux   = suppressWarnings(as.numeric(
      .bench_sh("awk '/MemTotal/ {print $2 * 1024}' /proc/meminfo"))),
    Windows = suppressWarnings(as.numeric(.bench_sh(
      "wmic computersystem get TotalPhysicalMemory /value | findstr ="))),
    NA_real_)
  if (length(b) != 1L || is.na(b)) NA_real_ else round(b / 1024^3, 1)
}


## The BLAS / LAPACK R is using, and whether they are R's bundled REFERENCE
## implementations (Rblas / Rlapack, the netlib code), which are much slower on
## dense matrices than OpenBLAS, MKL or Accelerate. `blas` / `lapack` are the
## library paths as sessionInfo() reports them; `matprod` is
## sessionInfo()$matprod. R on Windows ships only the reference libraries, and
## sessionInfo() leaves the BLAS path empty there, so an empty BLAS entry on
## Windows is read as the bundled Rblas.dll (a replacement Rblas.dll cannot be
## told apart from the path alone). The arguments are injectable so the
## classification is testable on any OS.
.blas_report <- function(blas = utils::sessionInfo()$BLAS,
                         lapack = utils::sessionInfo()$LAPACK,
                         matprod = utils::sessionInfo()$matprod,
                         sysname = Sys.info()[["sysname"]]) {
  one <- function(v) if (length(v) == 1L && !is.na(v) && nzchar(v)) as.character(v) else NA_character_
  blas <- one(blas); lapack <- one(lapack)
  ref <- function(path, lib)
    !is.na(path) && grepl(paste0("^(lib)?", lib, "[.]"), basename(path))
  blas_ref <- ref(blas, "Rblas") ||
    (is.na(blas) && identical(sysname, "Windows"))
  lapack_ref <- ref(lapack, "Rlapack")
  note <- NA_character_
  if (blas_ref || lapack_ref) {
    which_ref <- if (blas_ref && lapack_ref) "BLAS and LAPACK"
                 else if (blas_ref) "BLAS" else "LAPACK"
    note <- if (!blas_ref) paste0(
      "R is using its reference LAPACK (Rlapack) with an optimised BLAS: ",
      "matrix products are fast, but QZ and eigen decompositions use the ",
      "unoptimised LAPACK routines.") else paste0(
      "R is using its reference ", which_ref, " (",
      paste(c(if (blas_ref) "Rblas", if (lapack_ref) "Rlapack"), collapse = ", "),
      "). Dense linear algebra (QZ, eigen, Lyapunov, Kalman filtering with a ",
      "large state vector) is then much slower than with OpenBLAS, MKL or ",
      "Accelerate. To switch: Windows -- replace R's Rblas.dll with an ",
      "optimised build (see the R for Windows FAQ); macOS -- link R's BLAS ",
      "to Accelerate/vecLib or OpenBLAS; Debian/Ubuntu -- install OpenBLAS ",
      "and select it with update-alternatives; Fedora -- FlexiBLAS (see R ",
      "Installation and Administration, section 'BLAS').")
  }
  list(blas = blas, lapack = lapack,
       lapack_version = La_version(),
       matprod = one(matprod), blas_reference = blas_ref,
       lapack_reference = lapack_ref, note = note)
}


#' System information for a benchmark record
#'
#' Everything needed to interpret a \code{\link{dynhr_benchmark}} result on a
#' different machine: hardware, OS, R build, the numerical libraries R is
#' actually linked against, and the exact dynhr build.
#'
#' The BLAS/LAPACK entries matter more than they look. R shipped with its
#' reference BLAS and R linked against Accelerate or OpenBLAS are different
#' machines for this purpose, and the difference is invisible in
#' \code{R.version}. Two benchmark results whose \code{blas} differ are not
#' measuring the same software stack, whatever the hardware says.
#'
#' \code{blas_reference} and \code{lapack_reference} are \code{TRUE} when R is
#' using its bundled reference BLAS / LAPACK (\code{Rblas} / \code{Rlapack}),
#' and \code{blas_note} then says what that costs and how to switch to an
#' optimised library (otherwise \code{NA}). \code{matprod} is R's
#' matrix-product setting as \code{sessionInfo()} reports it.
#'
#' The \code{dynhr_*} build columns describe dynhr's OWN compiled library as
#' the compiler actually built it, read from the predefined compiler macros of
#' the package's C++ -- unlike \code{cxxflags}, which is R's default flag set
#' and says nothing about a package installed with a personal Makevars (see
#' \code{\link{dynhr_install_native}}). \code{dynhr_optimize} is whether the
#' library was compiled with optimisation (\code{FALSE} for a
#' \code{devtools::load_all()} build, which is several times slower);
#' \code{dynhr_simd} lists the instruction sets enabled at compile time
#' (e.g. \code{"sse2"}, \code{"avx2+fma"}, \code{"neon"}); \code{dynhr_fma}
#' says whether fused multiply-add may be contracted; \code{dynhr_fast_math}
#' is \code{TRUE} if the library was compiled with \code{-ffast-math}, a
#' configuration whose results are not supported, in which case
#' \code{dynhr_build_note} says so; \code{dynhr_compiler} is the compiler name
#' and version and \code{dynhr_cplusplus} the C++ standard macro.
#'
#' @return A one-row \code{data.frame} of character/numeric fields.
#' @examples
#' dynhr_system_info()
#' @export
dynhr_system_info <- function() {
  si <- tryCatch(utils::sessionInfo(), error = function(e) NULL)
  br <- .blas_report(si$BLAS, si$LAPACK, si$matprod)
  nm <- Sys.info()
  bi <- .build_info_report()
  gc_stamp <- tryCatch({
    p <- system.file("GIT_COMMIT", package = "dynhr")
    if (nzchar(p)) readLines(p, warn = FALSE)[1L] else NA_character_
  }, error = function(e) NA_character_)
  data.frame(
    stringsAsFactors = FALSE,
    timestamp      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    hostname       = unname(nm[["nodename"]]),
    os             = paste(nm[["sysname"]], nm[["release"]]),
    arch           = unname(nm[["machine"]]),
    cpu_model      = .bench_cpu_model(),
    cores_physical = tryCatch(parallel::detectCores(logical = FALSE),
                              error = function(e) NA_integer_),
    cores_logical  = tryCatch(parallel::detectCores(logical = TRUE),
                              error = function(e) NA_integer_),
    ram_gb         = .bench_ram_gb(),
    r_version      = paste(R.version$major, R.version$minor, sep = "."),
    r_platform     = R.version$platform,
    blas           = if (!is.null(si$BLAS)) si$BLAS else NA_character_,
    lapack         = if (!is.null(si$LAPACK)) si$LAPACK else NA_character_,
    lapack_version = br$lapack_version,
    matprod        = br$matprod,
    blas_reference   = br$blas_reference,
    lapack_reference = br$lapack_reference,
    blas_note      = br$note,
    cxx            = .bench_sh(paste(shQuote(file.path(R.home("bin"), "R")),
                                     "CMD config CXX")),
    cxxflags       = .bench_sh(paste(shQuote(file.path(R.home("bin"), "R")),
                                     "CMD config CXXFLAGS")),
    dynhr_version  = as.character(utils::packageVersion("dynhr")),
    dynhr_commit   = gc_stamp,
    dynhr_optimize   = bi$optimize,
    dynhr_simd       = bi$simd,
    dynhr_fma        = bi$fma,
    dynhr_fast_math  = bi$fast_math,
    dynhr_compiler   = bi$compiler,
    dynhr_cplusplus  = bi$cplusplus,
    dynhr_build_note = bi$note
  )
}


## Summarise the compiled library's own build facts (see src/build_info.cpp)
## into the fields dynhr_system_info() reports. The probe is injectable so the
## fast-math note and the SIMD label are testable without rebuilding the DLL.
.build_info_report <- function(bi = .dynhr_build_info()) {
  simd <- c(if (isTRUE(bi$avx512f)) "avx512f",
            if (isTRUE(bi$avx2)) "avx2" else if (isTRUE(bi$avx)) "avx",
            if (isTRUE(bi$sse2)) "sse2",
            if (isTRUE(bi$neon)) "neon")
  fma <- isTRUE(bi$fma) || isTRUE(bi$fp_fast_fma)
  simd_lab <- paste(simd, collapse = "+")
  if (!nzchar(simd_lab)) simd_lab <- "none"
  if (isTRUE(bi$fma)) simd_lab <- paste0(simd_lab, "+fma")
  note <- NA_character_
  if (isTRUE(bi$fast_math))
    note <- paste0("dynhr's compiled library was built with -ffast-math: ",
                   "results from this build are not supported (fast-math ",
                   "reassociates floating-point sums and assumes no NaN/Inf, ",
                   "which breaks the Kalman filter and the likelihood ",
                   "fingerprints). Reinstall without fast-math flags.")
  else if (!isTRUE(bi$optimize))
    note <- paste0("dynhr's compiled library was built WITHOUT optimisation ",
                   "(for example by devtools::load_all()); compiled code runs ",
                   "several times slower than an installed build.")
  list(optimize = isTRUE(bi$optimize), simd = simd_lab, fma = fma,
       fast_math = isTRUE(bi$fast_math),
       compiler = trimws(paste(bi$compiler, bi$compiler_version)),
       cplusplus = as.numeric(bi$cplusplus), note = note)
}


## Build the fixed benchmark problem. Separated out so the pilot, the sweep and
## the tests all provably use the SAME workload -- if this were inlined, a
## change to the pilot's setup could silently make it a different problem from
## the one being timed.
.bench_problem <- function(model = "sw2007") {
  model <- match.arg(model, "sw2007")
  need <- function(f) {
    p <- system.file("extdata", "models", f, package = "dynhr")
    if (!nzchar(p))
      stop("dynhr_benchmark: benchmark asset '", f, "' is not installed. ",
           "It ships in inst/extdata/models; reinstall dynhr.", call. = FALSE)
    p
  }
  obs <- c("dy", "dc", "dinve", "labobs", "pinfobs", "dw", "robs")
  m   <- suppressWarnings(parse_mod(need("sw2007.mod"), verbose = FALSE))
  cm  <- suppressWarnings(compile_model(m, verbose = FALSE))
  pr  <- extract_prior_spec(m, verbose = FALSE)
  raw <- utils::read.csv(need("sw2007_data.csv"))
  ## first_obs = 71, matching the .mod's own estimation command. presample and
  ## lik_init are deliberately NOT reproduced -- see sw2007_SOURCE.md; this is a
  ## fixed workload, not a Dynare-parity run.
  dat <- as.matrix(raw[71:nrow(raw), obs, drop = FALSE])
  mode_v <- utils::read.csv(need("sw2007_mode.csv"))$mode
  hh     <- as.matrix(utils::read.csv(need("sw2007_hessian.csv")))
  if (length(mode_v) != nrow(pr))
    stop("dynhr_benchmark: mode has ", length(mode_v), " entries but the ",
         "model declares ", nrow(pr), " estimated parameters.", call. = FALSE)
  dimnames(hh) <- list(pr$name, pr$name)
  theta0 <- stats::setNames(mode_v, pr$name)
  lp <- make_log_posterior(m, dat, pr, obs, cm)
  lp0 <- lp(theta0)
  list(model = m, compiled = cm, prior_spec = pr, obs = obs, data = dat,
       log_post_fn = lp, theta0 = theta0, hessian = hh,
       ## Proposal = inverse mode Hessian, the standard RWMH metric. Symmetrised
       ## because a MAT-file Hessian is only symmetric to round-off and chol()
       ## is entitled to object.
       Sigma_prop = { S <- solve(hh); (S + t(S)) / 2 },
       n_par = nrow(pr), n_obs = ncol(dat), n_periods = nrow(dat),
       ## unname(): the log-posterior inherits theta's first name, so the
       ## fingerprint would otherwise be a NAMED scalar and every comparison of
       ## it -- including across machines -- would have to remember to strip it.
       logpost0 = unname(if (is.list(lp0)) lp0$logpost else lp0))
}


## The speed-relevant configuration the workload actually runs with. The
## benchmark pins nothing itself -- it inherits the package defaults, which are
## the fastest settings -- so a session that changed one (e.g.
## options(dynhr.use_rcpp = FALSE): the R Kalman path, ~5x slower on sw2007;
## the setting is replayed into the mirai daemons) would silently produce a
## slower, non-comparable result. Recorded in the result and checked against
## the defaults.
.bench_expected_config <- list(use_rcpp = TRUE, kf_method = "standard")

.bench_config <- function(p) {
  m  <- p$model
  ss <- solve_steady_state(m, p$compiled, m$param_values, verbose = FALSE)$ss
  dr <- solve_perturbation(m, p$compiled, ss, m$param_values, verbose = FALSE)
  ## same call shape as the Gaussian posterior's filter (method / lik_init
  ## "auto"); the method choice depends on the state dimension, the data and
  ## the init, not on theta
  kf <- kalman_filter(p$data, dr, m, m$param_values, p$obs)
  list(use_rcpp  = isTRUE(getOption("dynhr.use_rcpp", TRUE)),
       kf_method = kf$diagnostics$method_used %||% kf$method %||% NA_character_,
       lik_init  = kf$diagnostics$lik_init_used %||% kf$lik_init %||%
                   NA_character_)
}

## Names of the config fields that differ from .bench_expected_config.
.bench_config_deviations <- function(cfg) {
  exp <- .bench_expected_config
  names(exp)[!vapply(names(exp), function(k) identical(cfg[[k]], exp[[k]]),
                     logical(1L))]
}


## Row label for one chain within a core setting: 4 chains -> "4a".."4d".
## Falls back to "32.27" past the alphabet rather than recycling letters, which
## would make two different chains print the same label.
.bench_chain_label <- function(k, chain) {
  ifelse(chain >= 1L & chain <= 26L,
         paste0(k, letters[pmax(1L, pmin(26L, chain))]),
         paste0(k, ".", chain))
}


## Default core ladder: 1, 2, then even counts to the logical core limit.
.bench_core_ladder <- function(max_cores = NULL) {
  mx <- if (is.null(max_cores))
    tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
  else max_cores
  if (!is.finite(mx) || mx < 1) mx <- 1L
  mx <- as.integer(mx)
  unique(as.integer(c(1L, if (mx >= 2L) seq(2L, mx, by = 2L), mx)))
}


#' Benchmark this machine on a fixed DSGE estimation workload
#'
#' Runs Smets & Wouters (2007) -- 36 estimated parameters, 7 observables, 160
#' quarters -- through random-walk Metropolis at a sweep of core counts, and
#' returns the throughputs together with the system information needed to
#' compare the result against another machine.
#'
#' @section What is held fixed:
#' The model, the data, the starting point (the published posterior mode), the
#' proposal (the inverse published mode Hessian) and the seed. Only the number
#' of parallel chains varies across the sweep. Each setting runs \code{n_draws}
#' draws PER CHAIN, so the work per chain is identical at every core count and
#' perfect scaling appears as constant elapsed time with linearly rising
#' throughput.
#'
#' @section Comparing across runs:
#' Every reported measure is normalised per draw or per second, so runs with
#' different \code{n_draws} -- including a slow machine that had to run fewer --
#' are directly comparable. Before comparing two results, check that
#' \code{logpost_check} agrees to the last digit and that \code{system$blas}
#' matches: an identical fingerprint means the two machines ran the same
#' arithmetic, and a differing BLAS means they did not.
#'
#' @param cores Integer vector of core counts to sweep. \code{NULL} (default)
#'   uses 1, 2, and even counts up to the logical core limit.
#' @param seconds_per_setting Target wall-clock seconds for each core setting.
#'   \code{n_draws} is calibrated from a short pilot to hit this. The default of
#'   120 is chosen so the run is SUSTAINED: a laptop's advantage often lies in
#'   how long it holds its clocks, which a few-second burst cannot see. Ignored
#'   when \code{n_draws} is given.
#' @param n_draws Draws per chain. \code{NULL} (default) calibrates from the
#'   pilot; supply a value to force an exact draw count.
#' @param model Benchmark workload. Only \code{"sw2007"} at present.
#' @param seed Base seed. Fixed by default so a rerun on the same machine is
#'   reproducible.
#' @param progress Passed to the sampler.
#'
#' @return An object of class \code{dynhr_benchmark}: a list with
#'   \code{$system} (one-row data.frame, see \code{\link{dynhr_system_info}}),
#'   \code{$problem} (workload description and \code{logpost_check}),
#'   \code{$settings} (including \code{config}: \code{use_rcpp}, the Kalman
#'   \code{kf_method} and \code{lik_init} the workload runs with, and
#'   \code{config_default}, \code{FALSE} -- with a warning of class
#'   \code{dynhr_warning_benchmark_config} -- when they are not the package
#'   defaults, making the throughputs non-comparable), and \code{$results}, a
#'   data.frame with one row per core
#'   setting and columns \code{cores}, \code{n_draws}, \code{total_draws},
#'   \code{elapsed_sec}, \code{overhead_sec}, \code{draws_per_sec},
#'   \code{draws_per_sec_adj}, \code{us_per_draw}, \code{speedup},
#'   \code{efficiency}, \code{accept_rate}.
#'
#' @examples
#' \donttest{
#' # Quick smoke run; the real thing wants the defaults.
#' b <- dynhr_benchmark(cores = c(1, 2), n_draws = 200)
#' print(b)
#' }
#' @seealso \code{\link{dynhr_system_info}}
#' @export
dynhr_benchmark <- function(cores = NULL,
                            seconds_per_setting = 120,
                            n_draws = NULL,
                            model = "sw2007",
                            seed = 20260730L,
                            progress = FALSE) {
  ## Own the message epoch for this run: repeat-suppressed warnings
  ## (`.dynhr_warn(once = TRUE)`) are keyed within it and re-arm for the
  ## next run, and the close reports what it suppressed. A nested call
  ## inherits this epoch rather than opening a second one.
  .dynhr_run_epoch <- .dynhr_epoch("dynhr_benchmark")
  on.exit(.dynhr_close_epoch(.dynhr_run_epoch), add = TRUE)
  if (!is.null(cores)) {
    if (!is.numeric(cores) || !length(cores) || any(!is.finite(cores)) ||
        any(cores < 1) || any(cores != as.integer(cores)))
      stop("dynhr_benchmark: `cores` must be positive whole numbers.")
    cores <- sort(unique(as.integer(cores)))
  } else {
    cores <- .bench_core_ladder()
  }
  if (!is.null(n_draws)) {
    if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
        n_draws < 1)
      stop("dynhr_benchmark: `n_draws` must be a single positive number.")
    n_draws <- as.integer(n_draws)
  }
  if (!is.numeric(seconds_per_setting) || length(seconds_per_setting) != 1L ||
      !is.finite(seconds_per_setting) || seconds_per_setting <= 0)
    stop("dynhr_benchmark: `seconds_per_setting` must be a positive number.")
  ## Validated HERE, not on first use inside .bench_problem(), so an unknown
  ## model errors before the run announces it is building one.
  model <- match.arg(model, "sw2007")

  sys <- dynhr_system_info()
  .dynhr_inform("dynhr_benchmark: building the ", model, " workload ...")
  p <- .bench_problem(model)
  .dynhr_inform(sprintf("  %d parameters, %d observables, %d periods; logpost = %.6f",
                  p$n_par, p$n_obs, p$n_periods, p$logpost0))
  cfg <- .bench_config(p)
  dev <- .bench_config_deviations(cfg)
  if (length(dev))
    .dynhr_warn("dynhr_benchmark: the workload is not running with the ",
                "default (fastest) configuration -- ",
                paste(sprintf("%s = %s (default %s)", dev,
                              vapply(dev, function(k) format(cfg[[k]]), ""),
                              vapply(dev, function(k)
                                format(.bench_expected_config[[k]]), "")),
                      collapse = "; "),
                ". Throughputs are not comparable with default-configuration ",
                "runs on other machines.",
                class = "dynhr_warning_benchmark_config")

  ## Pilot: a single serial chain, used ONLY to size n_draws. Timed separately
  ## from the sweep so its cost never enters a reported throughput.
  pilot_draws <- 300L
  .dynhr_inform("dynhr_benchmark: pilot (", pilot_draws, " serial draws) ...")
  t0 <- proc.time()[["elapsed"]]
  pilot <- rwmh(p$log_post_fn, p$theta0, p$Sigma_prop,
                n_draws = pilot_draws, n_burn = 0L, verbose = FALSE)
  pilot_sec <- proc.time()[["elapsed"]] - t0
  per_draw <- pilot_sec / pilot_draws
  if (is.null(n_draws))
    n_draws <- max(200L, as.integer(ceiling(seconds_per_setting / per_draw)))
  .dynhr_inform(sprintf("  %.3f ms/draw serial -> %d draws per chain (~%.0f s/setting)",
                  per_draw * 1000, n_draws, n_draws * per_draw))

  ## run_mcmc_mirai reports via cat(), not message(), so it needs capturing
  ## rather than suppressing -- otherwise the sampler's own chatter interleaves
  ## with the benchmark's and the result table is unreadable.
  run_one <- function(k, nd) {
    out <- NULL
    t0 <- proc.time()[["elapsed"]]
    utils::capture.output(suppressMessages(
      out <- run_mcmc_mirai(
        prior_spec  = p$prior_spec,
        theta_mode  = p$theta0,
        Sigma_prop  = p$Sigma_prop,
        log_post_fn = p$log_post_fn,
        n_chains    = k, n_cores = k,
        n_draws     = nd, n_burn = 0L,
        ## Dynare's own mh_jscale for this model, and NOT a cosmetic choice.
        ## run_mcmc_mirai defaults to 1.50, which on an inverse-Hessian proposal
        ## is far too wide here: 0.9% acceptance against 17.0% at 0.20. The
        ## tempting reading is that per-draw cost is unaffected because the
        ## posterior is evaluated on every proposal accepted or not. That is
        ## WRONG, and measurably so -- a proposal that wide mostly lands outside
        ## the determinacy region, where the Blanchard-Kahn check fails and the
        ## posterior returns -Inf via an early exit WITHOUT solving the model.
        ## Measured: 1142 draws/s at scale 1.50 against 555 at 0.20, so the
        ## benchmark ran ~2x fast by mostly not doing the work it claims to
        ## time. A cheap rejection is not a draw.
        mh_scale    = 0.20,
        seed_base   = seed, progress = progress)))
    list(sec = proc.time()[["elapsed"]] - t0, out = out)
  }

  res <- vector("list", length(cores))
  by_chain <- vector("list", length(cores))
  for (i in seq_along(cores)) {
    k <- cores[i]
    ## Fixed pool/dispatch cost measured at the SAME core count with a token
    ## draw count. Without it the sweep charges high-k settings for daemon
    ## startup and reports a scaling defect that is really a fixed cost.
    ov <- run_one(k, 25L)$sec
    .dynhr_inform(sprintf("dynhr_benchmark: %2d core(s) ...", k))
    r <- run_one(k, n_draws)
    total <- as.numeric(k) * n_draws
    ## chain_stats is a data.frame (one row per chain); the old
    ## vapply(chain_stats, function(s) s$accept_rate) iterated its COLUMNS, so
    ## every call failed on `$` of an atomic vector and a catch-all handler
    ## turned that into NA -- the acceptance rate was always reported as NA.
    cs  <- r$out$chain_stats
    acc <- if (is.data.frame(cs) && "accept_rate" %in% names(cs))
      mean(cs$accept_rate, na.rm = TRUE) else NA_real_
    ## Overhead-adjusted throughput is reported ONLY when the sampling actually
    ## dominates the fixed pool cost. Dividing by a near-zero (or negative)
    ## `sec - ov` is not a sharper measurement, it is a divide-by-zero wearing a
    ## number: an earlier revision returned 1.8e18 draws/s for the 1-core row of
    ## a 400-draw smoke run and silently zeroed every speedup that referenced
    ## it. Below the 2x margin the honest answer is NA -- raise n_draws.
    adj <- if (r$sec > 2 * ov) total / (r$sec - ov) else NA_real_

    ## PER-CHAIN rows. The aggregate above is barrier-bound -- run_mcmc_mirai
    ## returns only when every chain has, so one slow chain sets the elapsed
    ## time for all of them and the aggregate cannot show WHY. On a
    ## heterogeneous CPU (Apple P/E cores, Intel P/E, any machine with boost
    ## asymmetry) that is exactly the information wanted: k chains doing
    ## identical work should take identical time, and any spread is the cores
    ## differing. chain_stats$elapsed_min is measured inside the worker, so it
    ## excludes host-side pool setup and is a cleaner per-core figure than the
    ## host wall clock.
    cs <- r$out$chain_stats
    bc <- NULL
    if (is.data.frame(cs) && nrow(cs) && "elapsed_min" %in% names(cs)) {
      sec_ch <- cs$elapsed_min * 60
      bc <- data.frame(
        cores = k,
        chain = cs$chain,
        label = .bench_chain_label(k, cs$chain),
        elapsed_sec = round(sec_ch, 3),
        draws_per_sec = n_draws / sec_ch,
        us_per_draw = 1e6 * sec_ch / n_draws,
        accept_rate = cs$accept_rate %||% NA_real_,
        stringsAsFactors = FALSE)
      bc <- bc[order(bc$chain), , drop = FALSE]
    }
    by_chain[[i]] <- bc
    ## Slowest/fastest per-chain throughput. 1.0 means a homogeneous set of
    ## cores; a step change as `cores` grows past the performance-core count is
    ## the signature of work landing on efficiency cores.
    spread <- if (is.null(bc) || nrow(bc) < 2L) NA_real_ else
      max(bc$draws_per_sec) / min(bc$draws_per_sec)

    res[[i]] <- data.frame(
      cores = k, n_draws = n_draws, total_draws = total,
      elapsed_sec = round(r$sec, 3), overhead_sec = round(ov, 3),
      draws_per_sec = total / r$sec,
      draws_per_sec_adj = adj,
      us_per_draw = 1e6 * r$sec / total,
      chain_spread = spread,
      accept_rate = acc)
    .dynhr_inform(sprintf("    %.1f draws/s (%.1f adj), %.1f us/draw%s",
                    res[[i]]$draws_per_sec, res[[i]]$draws_per_sec_adj,
                    res[[i]]$us_per_draw,
                    if (is.na(spread)) "" else
                      sprintf(", chain spread %.2fx", spread)))
  }
  out <- do.call(rbind, res)
  ## Speedup uses the adjusted throughput only if EVERY row has one -- mixing
  ## adjusted and raw across rows of the same column would compare two different
  ## measurements and call the difference scaling.
  basis <- if (anyNA(out$draws_per_sec_adj)) "raw" else "adjusted"
  if (basis == "raw")
    .dynhr_warn("dynhr_benchmark: pool overhead is not negligible at every core ",
            "count, so speedup is computed from RAW throughput and is ",
            "pessimistic at high core counts. Raise `n_draws` or ",
            "`seconds_per_setting` for a clean scaling curve.", call. = FALSE)
  thr <- if (basis == "raw") out$draws_per_sec else out$draws_per_sec_adj
  base <- thr[out$cores == min(out$cores)][1L]
  out$speedup    <- thr / base
  out$efficiency <- out$speedup / (out$cores / min(out$cores))

  bych <- if (all(vapply(by_chain, is.null, logical(1L)))) NULL else
    do.call(rbind, by_chain[!vapply(by_chain, is.null, logical(1L))])

  structure(list(
    system   = sys,
    problem  = list(model = model, n_par = p$n_par, n_obs = p$n_obs,
                    n_periods = p$n_periods, logpost_check = p$logpost0),
    settings = list(cores = cores, n_draws = n_draws, seed = seed,
                    seconds_per_setting = seconds_per_setting,
                    serial_ms_per_draw = per_draw * 1000,
                    speedup_basis = basis,
                    config = cfg, config_default = !length(dev)),
    results  = out,
    by_chain = bych
  ), class = "dynhr_benchmark")
}


#' @param x A \code{dynhr_benchmark} object.
#' @param ... Ignored.
#' @rdname dynhr_benchmark
#' @export
print.dynhr_benchmark <- function(x, ...) {
  s <- x$system
  cat("dynhr benchmark\n")
  cat(sprintf("  %s | %s | %s cores (%s physical) | %s GB\n",
              s$cpu_model, s$os, s$cores_logical, s$cores_physical, s$ram_gb))
  cat(sprintf("  R %s (%s) | BLAS: %s\n", s$r_version, s$r_platform,
              basename(s$blas %||% "?")))
  if (isTRUE(s$blas_reference) || isTRUE(s$lapack_reference))
    cat("  ** reference BLAS/LAPACK in use -- see dynhr_system_info()$blas_note **\n")
  cat(sprintf("  dynhr %s (%s) | %s\n", s$dynhr_version,
              substr(s$dynhr_commit %||% "?", 1, 8), s$timestamp))
  cat(sprintf("\n  workload: %s, %d params, %d obs, %d periods\n",
              x$problem$model, x$problem$n_par, x$problem$n_obs,
              x$problem$n_periods))
  cat(sprintf("  fingerprint (logpost): %.10f\n", x$problem$logpost_check))
  cat(sprintf("  %d draws/chain, seed %d\n", x$settings$n_draws,
              x$settings$seed))
  cfg <- x$settings$config
  if (!is.null(cfg))
    cat(sprintf("  config: use_rcpp = %s, Kalman method = %s, init = %s%s\n",
                cfg$use_rcpp, cfg$kf_method, cfg$lik_init,
                if (isTRUE(x$settings$config_default)) " (default)"
                else "  ** NOT the default configuration -- not comparable **"))
  cat("\n")
  r <- x$results
  cat(sprintf("  %5s %12s %12s %10s %8s %10s\n",
              "cores", "draws/s", "draws/s adj", "us/draw", "speedup", "effic."))
  for (i in seq_len(nrow(r)))
    cat(sprintf("  %5d %12.1f %12s %10.1f %8.2f %9.0f%%\n",
                r$cores[i], r$draws_per_sec[i],
                if (is.na(r$draws_per_sec_adj[i])) "-" else
                  sprintf("%.1f", r$draws_per_sec_adj[i]),
                r$us_per_draw[i], r$speedup[i], 100 * r$efficiency[i]))
  if (identical(x$settings$speedup_basis, "raw"))
    cat("\n  NOTE: speedup is from RAW throughput -- pool overhead was not",
        "\n  negligible at every core count. Raise n_draws for a clean curve.\n")

  if (!is.null(x$by_chain) && nrow(x$by_chain)) {
    cat("\n  Per-chain (each chain does the same work; spread = slowest vs",
        "\n  fastest core the scheduler gave us):\n\n")
    cat(sprintf("  %8s %12s %12s %10s %12s\n",
                "chain", "elapsed s", "draws/s", "us/draw", "accept"))
    b <- x$by_chain
    for (k in unique(b$cores)) {
      bk <- b[b$cores == k, , drop = FALSE]
      for (j in seq_len(nrow(bk)))
        cat(sprintf("  %8s %12.2f %12.1f %10.1f %11.1f%%\n",
                    bk$label[j], bk$elapsed_sec[j], bk$draws_per_sec[j],
                    bk$us_per_draw[j], 100 * bk$accept_rate[j]))
      if (nrow(bk) > 1L)
        cat(sprintf("  %8s %12s %12s %10s %12s\n", "",
                    sprintf("spread %.2fx",
                            max(bk$draws_per_sec) / min(bk$draws_per_sec)),
                    "", "", ""))
    }
    cat("\n  A spread that jumps once `cores` passes the performance-core",
        "\n  count is work landing on efficiency cores. Chain-to-core",
        "\n  assignment is the OS scheduler's, not dynhr's, so a chain label",
        "\n  identifies a chain, NOT a specific physical core.\n")
  }
  cat("\n  Comparable across runs with different n_draws: every column above",
      "\n  is per-draw or per-second. Compare only against runs with the same",
      "\n  fingerprint and BLAS.\n")
  invisible(x)
}


## --------------------------------------------------------------------------
## Posterior benchmark: time ONE user-supplied log-posterior closure and a
## short seeded RWMH run, and record enough build / system / closure detail
## that two runs (two builds, two machines) can be compared line by line.
## Independent of any particular model: the closure is the user's.
## --------------------------------------------------------------------------

## Flatten the closure's stats-like attributes (anything but srcref and the
## structural ones) into scalar fields named attr_<attribute>[_<field>]. Read
## generically so a closure that carries filter-method or cache/fallback
## counters is recorded whatever those attributes are called; a closure that
## carries none contributes no columns.
.bench_closure_attrs <- function(fn) {
  at <- attributes(fn)
  at <- at[setdiff(names(at), c("srcref", "class", "names", "dim", "dimnames"))]
  out <- list()
  scalar <- function(v) {
    if (is.null(v) || !length(v)) return(NA_character_)
    if (is.atomic(v)) return(paste(as.character(v), collapse = ";"))
    NA_character_
  }
  for (nm in names(at)) {
    v <- at[[nm]]
    if (is.list(v) && length(v) && !is.null(names(v))) {
      for (k in names(v)) out[[paste0("attr_", nm, "_", k)]] <- scalar(v[[k]])
    } else if (is.environment(v) || is.function(v)) {
      next
    } else {
      out[[paste0("attr_", nm)]] <- scalar(v)
    }
  }
  out
}

## T, n_obs, missing-cell count and first / last period holding a missing
## value, from a data matrix (NULL -> all NA).
.bench_data_shape <- function(data) {
  if (is.null(data) || !(is.matrix(data) || is.data.frame(data)))
    return(list(n_periods = NA_integer_, n_obs = NA_integer_,
                n_missing = NA_integer_, first_missing = NA_integer_,
                last_missing = NA_integer_))
  m <- as.matrix(data)
  miss_rows <- which(rowSums(is.na(m)) > 0L)
  list(n_periods = nrow(m), n_obs = ncol(m), n_missing = sum(is.na(m)),
       first_missing = if (length(miss_rows)) min(miss_rows) else NA_integer_,
       last_missing = if (length(miss_rows)) max(miss_rows) else NA_integer_)
}

## md5 of a character vector, via a temp file (base R only).
.bench_md5 <- function(txt) {
  f <- tempfile("dynhr-bench-digest-")
  on.exit(unlink(f), add = TRUE)
  writeLines(txt, f)
  unname(tools::md5sum(f))
}

#' Benchmark a log-posterior closure and a short seeded RWMH run
#'
#' Times repeated evaluations of your own log-posterior function and runs a
#' short seeded random-walk Metropolis chain, returning one row that records
#' the system, the dynhr build, the closure and the results. Appending rows
#' (\code{file}) from several machines or builds gives a table that compares
#' them like for like. It is independent of any model: the closure is yours;
#' \code{dynhr:::.bench_problem()} builds the Smets-Wouters closure used in
#' the examples.
#'
#' @section Columns:
#' \itemize{
#'   \item the \code{\link{dynhr_system_info}} fields (version, commit, R,
#'     compiler and flags, the compiled-library build columns, BLAS/LAPACK,
#'     matrix-product setting, cores);
#'   \item \code{attr_*}: every list-like or atomic attribute the closure
#'     carries (for example filter-method or cache and fallback counters),
#'     flattened as \code{attr_<attribute>_<field>}; none if it carries none;
#'   \item \code{n_periods}, \code{n_obs}, \code{n_state}, \code{n_missing},
#'     \code{first_missing}, \code{last_missing}: the data shape.
#'     \code{data} is read from the closure's environment when it holds a
#'     matrix called \code{data} (as \code{\link{make_log_posterior}} closures
#'     do) unless you pass it; \code{n_state} is \code{NA} unless you pass
#'     \code{n_state};
#'   \item \code{logpost}, \code{loglik}: the value at \code{theta}
#'     (\code{loglik} is \code{NA} if the closure returns a plain number);
#'   \item \code{ms_min}, \code{ms_median}, \code{ms_mean}: milliseconds per
#'     evaluation over \code{n_eval} calls after \code{n_warmup} untimed ones;
#'   \item \code{rwmh_elapsed_sec}, \code{rwmh_accept}: wall time and
#'     acceptance rate of the seeded chain;
#'   \item \code{rwmh_accept_digest}: md5 of the accept/reject sequence, equal
#'     between two builds exactly when they made identical decisions;
#'   \item \code{rwmh_draws_digest}: md5 of the draws at full precision, which
#'     changes with any round-off difference.
#' }
#'
#' @param log_post_fn Function of \code{theta} (and \code{...}) returning a
#'   number or a list with \code{logpost} and optionally \code{loglik}.
#' @param theta Evaluation point and starting value of the chain.
#' @param ... Further arguments passed to \code{log_post_fn}.
#' @param n_eval,n_warmup Timed and untimed warm-up evaluations.
#' @param rwmh_draws Draws of the seeded RWMH chain (no burn-in, scale 1).
#' @param Sigma_prop Proposal covariance. Default: diagonal with standard
#'   deviation \code{0.02 * max(|theta|, 0.1)} per parameter. Pass the same
#'   matrix on every machine you compare.
#' @param seed Seed of the chain. The caller's RNG state is restored.
#' @param file Optional CSV path; the row is appended (header written when the
#'   file is new). A file whose columns differ from this row is refused.
#' @param data Optional data matrix for the shape columns.
#' @param n_state Optional state dimension to record.
#'
#' @return A one-row \code{data.frame} (invisibly if \code{file} is given).
#' @examples
#' \donttest{
#' p <- dynhr:::.bench_problem()
#' r <- dynhr_benchmark_posterior(p$log_post_fn, p$theta0,
#'                                Sigma_prop = p$Sigma_prop,
#'                                n_eval = 5, rwmh_draws = 20)
#' r[, c("ms_median", "rwmh_accept", "rwmh_accept_digest")]
#' }
#' @seealso \code{\link{dynhr_benchmark}}, \code{\link{dynhr_system_info}}
#' @export
dynhr_benchmark_posterior <- function(log_post_fn, theta, ...,
                                      n_eval = 100L, n_warmup = 10L,
                                      rwmh_draws = 300L, Sigma_prop = NULL,
                                      seed = 1L, file = NULL,
                                      data = NULL, n_state = NA_integer_) {
  if (!is.function(log_post_fn))
    .dynhr_abort("dynhr_benchmark_posterior: `log_post_fn` must be a function.",
                 class = "dynhr_error_benchmark_args")
  if (!is.numeric(theta) || !length(theta) || anyNA(theta))
    .dynhr_abort("dynhr_benchmark_posterior: `theta` must be a numeric vector ",
                 "without missing values.", class = "dynhr_error_benchmark_args")
  whole <- function(v, min) is.numeric(v) && length(v) == 1L && is.finite(v) &&
    v >= min && v == as.integer(v)
  if (!whole(n_eval, 1) || !whole(n_warmup, 0) || !whole(rwmh_draws, 2) ||
      !whole(seed, -.Machine$integer.max))
    .dynhr_abort("dynhr_benchmark_posterior: `n_eval` >= 1, `n_warmup` >= 0, ",
                 "`rwmh_draws` >= 2 and `seed` must be whole numbers.",
                 class = "dynhr_error_benchmark_args")
  n_eval <- as.integer(n_eval); n_warmup <- as.integer(n_warmup)
  rwmh_draws <- as.integer(rwmh_draws); seed <- as.integer(seed)
  n_par <- length(theta)
  if (is.null(Sigma_prop)) {
    Sigma_prop <- diag((0.02 * pmax(abs(theta), 0.1))^2, nrow = n_par)
  } else if (!is.matrix(Sigma_prop) || any(dim(Sigma_prop) != n_par)) {
    .dynhr_abort("dynhr_benchmark_posterior: `Sigma_prop` must be a ", n_par,
                 " x ", n_par, " matrix.", class = "dynhr_error_benchmark_args")
  }

  extra <- list(...)
  f <- function(th) do.call(log_post_fn, c(list(th), extra))
  val <- function(r) {
    if (is.list(r))
      list(logpost = as.numeric(r$logpost)[1L],
           loglik = if (is.null(r$loglik)) NA_real_ else as.numeric(r$loglik)[1L])
    else list(logpost = as.numeric(r)[1L], loglik = NA_real_)
  }
  v0 <- val(f(theta))
  for (i in seq_len(n_warmup)) f(theta)
  ms <- numeric(n_eval)
  for (i in seq_len(n_eval)) {
    t0 <- proc.time()[["elapsed"]]
    f(theta)
    ms[i] <- 1000 * (proc.time()[["elapsed"]] - t0)
  }

  ## the seeded chain: restore the caller's RNG state afterwards
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  ch <- rwmh(function(th) list(logpost = val(f(th))$logpost), theta,Sigma_prop, n_draws = rwmh_draws, n_burn = 0L,
             scale = 1, verbose = FALSE)
  rw_sec <- proc.time()[["elapsed"]] - t0
  draws <- ch$full_chain
  ## accept/reject sequence: a step is an acceptance iff the state moved
  acc <- rowSums(abs(diff(draws))) > 0
  acc_digest <- .bench_md5(paste(as.integer(acc), collapse = ""))
  draws_digest <- .bench_md5(sprintf("%.17g", as.vector(draws)))

  if (is.null(data)) {
    d <- get0("data", envir = environment(log_post_fn), inherits = FALSE)
    if (is.matrix(d) || is.data.frame(d)) data <- d
  }
  shape <- .bench_data_shape(data)

  ca <- .bench_closure_attrs(log_post_fn)
  ca <- if (length(ca)) as.data.frame(ca, stringsAsFactors = FALSE) else
    data.frame(row.names = 1L)
  row <- cbind(
    dynhr_system_info(),
    ca,
    data.frame(n_periods = shape$n_periods, n_obs = shape$n_obs,
               n_state = as.integer(n_state), n_missing = shape$n_missing,
               first_missing = shape$first_missing,
               last_missing = shape$last_missing,
               n_par = n_par, n_eval = n_eval,
               logpost = v0$logpost, loglik = v0$loglik,
               ms_min = min(ms), ms_median = stats::median(ms),
               ms_mean = mean(ms),
               rwmh_draws = rwmh_draws, rwmh_elapsed_sec = rw_sec,
               rwmh_accept = mean(acc),
               rwmh_accept_digest = acc_digest,
               rwmh_draws_digest = draws_digest,
               stringsAsFactors = FALSE))

  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file))
      .dynhr_abort("dynhr_benchmark_posterior: `file` must be a single path.",
                   class = "dynhr_error_benchmark_args")
    if (file.exists(file) && file.size(file) > 0L) {
      hdr <- names(utils::read.csv(file, nrows = 1L, check.names = FALSE))
      if (!identical(hdr, names(row)))
        .dynhr_abort("dynhr_benchmark_posterior: '", file, "' has different ",
                     "columns from this result (a different closure ",
                     "attribute set or dynhr version); use a new file.",
                     class = "dynhr_error_benchmark_args")
      utils::write.table(row, file, sep = ",", append = TRUE,
                         col.names = FALSE, row.names = FALSE, qmethod = "double")
    } else {
      utils::write.table(row, file, sep = ",", col.names = TRUE,
                         row.names = FALSE, qmethod = "double")
    }
    return(invisible(row))
  }
  row
}
