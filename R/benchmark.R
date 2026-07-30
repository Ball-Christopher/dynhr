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
#' @return A one-row \code{data.frame} of character/numeric fields.
#' @examples
#' dynhr_system_info()
#' @export
dynhr_system_info <- function() {
  si <- tryCatch(utils::sessionInfo(), error = function(e) NULL)
  nm <- Sys.info()
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
    lapack_version = tryCatch(La_version(), error = function(e) NA_character_),
    cxx            = .bench_sh(paste(shQuote(file.path(R.home("bin"), "R")),
                                     "CMD config CXX")),
    cxxflags       = .bench_sh(paste(shQuote(file.path(R.home("bin"), "R")),
                                     "CMD config CXXFLAGS")),
    dynhr_version  = as.character(utils::packageVersion("dynhr")),
    dynhr_commit   = gc_stamp
  )
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
#'   \code{$settings}, and \code{$results}, a data.frame with one row per core
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
  message("dynhr_benchmark: building the ", model, " workload ...")
  p <- .bench_problem(model)
  message(sprintf("  %d parameters, %d observables, %d periods; logpost = %.6f",
                  p$n_par, p$n_obs, p$n_periods, p$logpost0))

  ## Pilot: a single serial chain, used ONLY to size n_draws. Timed separately
  ## from the sweep so its cost never enters a reported throughput.
  pilot_draws <- 300L
  message("dynhr_benchmark: pilot (", pilot_draws, " serial draws) ...")
  t0 <- proc.time()[["elapsed"]]
  pilot <- rwmh(p$log_post_fn, p$theta0, p$Sigma_prop,
                n_draws = pilot_draws, n_burn = 0L, verbose = FALSE)
  pilot_sec <- proc.time()[["elapsed"]] - t0
  per_draw <- pilot_sec / pilot_draws
  if (is.null(n_draws))
    n_draws <- max(200L, as.integer(ceiling(seconds_per_setting / per_draw)))
  message(sprintf("  %.3f ms/draw serial -> %d draws per chain (~%.0f s/setting)",
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
    message(sprintf("dynhr_benchmark: %2d core(s) ...", k))
    r <- run_one(k, n_draws)
    total <- as.numeric(k) * n_draws
    acc <- tryCatch(mean(vapply(r$out$chain_stats,
                                function(s) s$accept_rate %||% NA_real_, 0),
                         na.rm = TRUE), error = function(e) NA_real_)
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
    message(sprintf("    %.1f draws/s (%.1f adj), %.1f us/draw%s",
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
    warning("dynhr_benchmark: pool overhead is not negligible at every core ",
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
                    speedup_basis = basis),
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
  cat(sprintf("  dynhr %s (%s) | %s\n", s$dynhr_version,
              substr(s$dynhr_commit %||% "?", 1, 8), s$timestamp))
  cat(sprintf("\n  workload: %s, %d params, %d obs, %d periods\n",
              x$problem$model, x$problem$n_par, x$problem$n_obs,
              x$problem$n_periods))
  cat(sprintf("  fingerprint (logpost): %.10f\n", x$problem$logpost_check))
  cat(sprintf("  %d draws/chain, seed %d\n\n", x$settings$n_draws,
              x$settings$seed))
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
