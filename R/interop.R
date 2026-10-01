## R/interop.R
## --------------------------------------------------------------------------
## Ecosystem interop for posterior output.
##
##   posterior (Suggests): as_draws(), as_draws_array(), as_draws_df()
##   coda      (Suggests): as.mcmc(), as.mcmc.list()
##   base/stats          : as.data.frame(), coef(), vcov(), nobs() on
##                         dynhr_chains; coef(), vcov(), logLik(), nobs() on
##                         dynhr_mode_result.
##
## The posterior / coda methods are registered CONDITIONALLY from .onLoad
## (.s3_register(), the vctrs pattern) so dynhr does not Import either
## package: the method is registered when the generic's namespace loads.
##
## Every method reads the draws through ONE normaliser, .interop_draws(), so
## the chain structure, the log-posterior column and the SMC weights are
## decided in one place.
## --------------------------------------------------------------------------


# ============================================================================
# Conditional S3 registration
# ============================================================================

#' Register an S3 method for a generic in a Suggests package
#'
#' The standard vctrs \code{s3_register()} pattern: registers
#' \code{method} for \code{generic} ("pkg::generic") now if \code{pkg} is
#' already loaded, and from a \code{packageEvent(pkg, "onLoad")} hook
#' otherwise, so the package never has to be Imported.
#'
#' @param generic "pkg::generic".
#' @param class   Class to register the method for.
#' @param method  The method function.
#' @return \code{NULL}, invisibly.
#' @noRd
.s3_register <- function(generic, class, method) {
  ## Force now: `register` runs later (at the generic package's onLoad), and
  ## a lazy `class` promise would then see the caller's loop variable at its
  ## LAST value -- every method would be registered for the last class only.
  force(class)
  pieces <- strsplit(generic, "::", fixed = TRUE)[[1L]]
  if (length(pieces) != 2L || !is.function(method))
    .dynhr_abort(".s3_register: `generic` must be \"pkg::generic\" and ",
                 "`method` a function.")
  pkg <- pieces[[1L]]
  gen <- pieces[[2L]]
  register <- function(...) {
    ns <- asNamespace(pkg)
    if (exists(gen, envir = ns, inherits = FALSE)) {
      registerS3method(gen, class, method, envir = ns)
    } else if (identical(Sys.getenv("NOT_CRAN"), "true")) {
      .dynhr_warn(sprintf(
        "dynhr: cannot find generic `%s` in package %s to register the %s method.",
        gen, pkg, class))
    }
    invisible(NULL)
  }
  setHook(packageEvent(pkg, "onLoad"), register)
  if (isNamespaceLoaded(pkg)) register()
  invisible(NULL)
}


#' Register the posterior / coda interop methods (called from .onLoad)
#' @noRd
.register_interop_methods <- function() {
  for (cls in c("dynhr_chains", "dynhr_posterior_result")) {
    .s3_register("posterior::as_draws",       cls, as_draws.dynhr_chains)
    .s3_register("posterior::as_draws_array", cls, as_draws_array.dynhr_chains)
    .s3_register("posterior::as_draws_df",    cls, as_draws_df.dynhr_chains)
    .s3_register("coda::as.mcmc",             cls, as.mcmc.dynhr_chains)
    .s3_register("coda::as.mcmc.list",        cls, as.mcmc.list.dynhr_chains)
  }
  invisible(NULL)
}


# ============================================================================
# The one draw normaliser
# ============================================================================

#' Normalise a posterior result to per-chain draw matrices
#'
#' \describe{
#'   \item{dynhr_chains, multi-chain}{When \code{$chain_list} holds the
#'     per-chain results whose row-bind reproduces \code{$chain} (what
#'     \code{.compute_convergence()} builds), each becomes one chain, in
#'     order.}
#'   \item{dynhr_chains, weighted SMC cloud}{When \code{$chain} IS the
#'     particle cloud (no \code{$resample_idx}) and \code{$smc_weights} is a
#'     non-uniform vector of matching length, the cloud is one chain and the
#'     normalised weights are returned. Resampled SMC results (and a final
#'     stage that resampled, i.e. uniform weights) are plain equal-weight
#'     draws.}
#'   \item{dynhr_chains, otherwise}{One chain (\code{$chain}). A DIME
#'     ensemble is one chain: its walkers are coupled, not independent.}
#'   \item{dynhr_posterior_result}{Every stored chain of every method, in
#'     \code{names($chains)} order, each once. Estimation drivers resample
#'     SMC to equal weights, so these are never weighted.}
#' }
#'
#' @param x A \code{dynhr_chains} or \code{dynhr_posterior_result}.
#' @return \code{list(chains, logpost, weights)}: a list of draw matrices
#'   (named columns), a parallel list of log-posterior vectors or
#'   \code{NULL} (only when EVERY chain carries one of matching length), and
#'   the normalised weight vector or \code{NULL}.
#' @noRd
.interop_draws <- function(x) {
  name_cols <- function(m) {
    if (!is.matrix(m)) m <- as.matrix(m)
    if (is.null(colnames(m)))
      colnames(m) <- paste0("theta_", seq_len(ncol(m)))
    m
  }
  chain_of <- function(r) if (is.list(r) && !is.null(r$chain)) r$chain else r
  usable   <- function(m) is.matrix(m) && nrow(m) > 0L && ncol(m) > 0L

  weights <- NULL
  if (inherits(x, "dynhr_posterior_result")) {
    res <- list()
    for (m in names(x$chains))
      for (r in x$chains[[m]]$chains)
        if (usable(chain_of(r))) res <- c(res, list(r))
    if (length(res) == 0L)
      .dynhr_abort("This dynhr_posterior_result holds no posterior draws.",
                   class = "dynhr_error_no_draws")
    mats <- lapply(res, function(r) name_cols(chain_of(r)))
    lps  <- lapply(res, function(r) if (is.list(r)) r$post_logpost else NULL)
  } else {
    ch <- x$chain
    if (is.null(ch) || !usable(as.matrix(ch)))
      .dynhr_abort("This dynhr_chains object holds no posterior draws ",
                   "($chain is empty).", class = "dynhr_error_no_draws")
    ch <- name_cols(ch)
    w  <- x$smc_weights
    weighted <- is.null(x$resample_idx) && length(w) == nrow(ch) &&
      length(w) >= 2L && !all(w == w[1L])
    mats <- list(ch)
    lps  <- list(x$post_logpost)
    if (weighted) {
      weights <- as.numeric(w) / sum(w)
    } else if (length(x$chain_list) >= 1L) {
      cl   <- x$chain_list
      keep <- vapply(cl, function(r) usable(chain_of(r)), logical(1))
      sub  <- lapply(cl[keep], function(r) name_cols(chain_of(r)))
      if (length(sub) >= 1L &&
          identical(dim(do.call(rbind, sub)), dim(ch)) &&
          all(do.call(rbind, sub) == ch)) {
        mats <- sub
        lps  <- lapply(cl[keep], function(r) if (is.list(r)) r$post_logpost else NULL)
      } else {
        .dynhr_warn("dynhr_chains: $chain_list does not reproduce $chain; ",
                    "exporting the pooled $chain as a single chain.")
      }
    }
  }

  lp_ok <- all(vapply(seq_along(mats), function(k)
    is.numeric(lps[[k]]) && length(lps[[k]]) == nrow(mats[[k]]), logical(1)))
  list(chains  = mats,
       logpost = if (lp_ok) lapply(lps, as.numeric) else NULL,
       weights = weights)
}


#' Stack the normalised draws into one long data.frame
#' @noRd
.interop_long <- function(d, lp_name) {
  n_k <- vapply(d$chains, nrow, integer(1))
  out <- as.data.frame(do.call(rbind, d$chains), stringsAsFactors = FALSE)
  if (!is.null(d$logpost)) out[[lp_name]] <- unlist(d$logpost, use.names = FALSE)
  out$.chain     <- rep(seq_along(n_k), n_k)
  out$.iteration <- unlist(lapply(n_k, seq_len), use.names = FALSE)
  out
}


#' @noRd
.interop_need_posterior <- function() {
  if (!requireNamespace("posterior", quietly = TRUE))
    .dynhr_abort("The 'posterior' package is required: install.packages(\"posterior\").",
                 class = "dynhr_error_missing_package")
}


# ============================================================================
# posterior: as_draws(), as_draws_array(), as_draws_df()
# ============================================================================

#' Convert dynhr posterior output to other packages' formats
#'
#' Methods that hand a \code{\link{dynhr_chains}} object (or a
#' \code{dynhr_posterior_result} from \code{\link{run_posterior_estimation}}
#' / \code{\link{pmmh}}) to the \pkg{posterior}, \pkg{coda} and base-R
#' ecosystems, so \code{posterior::summarise_draws()}, \pkg{bayesplot} and
#' \pkg{loo}-style tooling work on dynhr output directly.
#'
#' @section Chain structure:
#' A multi-chain \code{dynhr_chains} (\code{$chain_list} present, as
#' \code{run_full_estimation()} / \code{run_posterior_estimation()} build
#' for \code{n_chains > 1}) keeps its chains: iterations x chains x
#' variables, in chain order. A single sampler run is one chain; a DIME
#' ensemble is one chain (its walkers are coupled, not independent chains).
#' An SMC particle set is exported as one "chain" too, but its particles are
#' not a Markov chain: R-hat / ESS computed on it are meaningless.
#' A \code{dynhr_posterior_result} exports every stored chain of every method
#' once, in method order. When per-draw log-posteriors are stored for every
#' chain they are added as the variable \code{lp__} (Stan's name, which
#' \pkg{bayesplot} recognises).
#'
#' @section Weighted SMC clouds:
#' When \code{$chain} is SMC's weighted particle cloud (non-uniform
#' \code{$smc_weights}, not yet resampled), the draws are the particles and
#' the weights are attached with \code{posterior::weight_draws()} -- nothing
#' is resampled, so no Monte-Carlo noise is added. Read them with
#' \code{weights()} (posterior's method for \code{stats::weights}). Note that
#' \code{posterior::summarise_draws()} IGNORES weights; call \code{posterior::resample_draws()} first for
#' weighted summaries. SMC results already resampled to equal weights (what
#' the estimation drivers return) are exported as ordinary draws.
#' \code{coef()} and \code{vcov()} use the weights; \pkg{coda} has no weight
#' concept, so \code{as.mcmc()} / \code{as.mcmc.list()} refuse a weighted
#' cloud rather than silently dropping the weights.
#'
#' @section Registration:
#' \pkg{posterior} and \pkg{coda} are Suggests: the methods for their
#' generics are registered when those packages load, and
#' \code{as_draws()} returns a \code{draws_df} (which, unlike a
#' \code{draws_array}, also admits chains of unequal length).
#'
#' @param x,object A \code{dynhr_chains} object (the posterior and coda
#'   methods also accept a \code{dynhr_posterior_result}).
#' @param stat \code{"mean"} (default) or \code{"median"}: the posterior
#'   point summary \code{coef()} returns (weighted for a weighted cloud).
#' @param row.names,optional Ignored (base \code{as.data.frame()} signature).
#' @param ... Ignored.
#' @return \code{as_draws()} / \code{as_draws_df()}: a \code{draws_df};
#'   \code{as_draws_array()}: a \code{draws_array} (chains must have equal
#'   length); \code{as.mcmc.list()}: a \code{coda} \code{mcmc.list} with one
#'   \code{mcmc} per chain; \code{as.mcmc()}: an \code{mcmc} (single chain
#'   only); \code{as.data.frame()}: one row per draw with columns
#'   \code{.chain}, \code{.iteration}, the parameters, \code{logpost} (when
#'   stored) and \code{.weight} (weighted clouds only); \code{coef()}: named
#'   posterior mean/median; \code{vcov()}: posterior covariance matrix;
#'   \code{nobs()}: the number of observations when recorded, else
#'   \code{NA}.
#' @seealso \code{\link{dynhr_chains}}, \code{\link{summary.dynhr_chains}},
#'   \code{\link{coef.dynhr_mode_result}}
#' @examples
#' ## A tiny two-chain object (normally returned by dynhr_mcmc() / run_full_estimation())
#' set.seed(1)
#' ch1 <- cbind(a = rnorm(50, 1), b = rnorm(50, -1))
#' ch2 <- cbind(a = rnorm(50, 1), b = rnorm(50, -1))
#' x <- structure(list(chain = rbind(ch1, ch2), sampler = "rwmh", n_chains = 2L,
#'                     chain_list = list(list(chain = ch1), list(chain = ch2))),
#'                class = c("dynhr_chains", "list"))
#' coef(x)
#' vcov(x)
#' head(as.data.frame(x))
#' if (requireNamespace("posterior", quietly = TRUE))
#'   posterior::summarise_draws(posterior::as_draws(x))
#' @name dynhr_interop
NULL


#' @rdname dynhr_interop
#' @usage \method{as_draws}{dynhr_chains}(x, ...)
as_draws.dynhr_chains <- function(x, ...) {
  as_draws_df.dynhr_chains(x, ...)
}

#' @rdname dynhr_interop
#' @usage \method{as_draws_df}{dynhr_chains}(x, ...)
as_draws_df.dynhr_chains <- function(x, ...) {
  .interop_need_posterior()
  d   <- .interop_draws(x)
  out <- posterior::as_draws_df(.interop_long(d, "lp__"))
  if (!is.null(d$weights))
    out <- posterior::weight_draws(out, weights = d$weights)
  out
}

#' @rdname dynhr_interop
#' @usage \method{as_draws_array}{dynhr_chains}(x, ...)
as_draws_array.dynhr_chains <- function(x, ...) {
  .interop_need_posterior()
  d   <- .interop_draws(x)
  n_k <- vapply(d$chains, nrow, integer(1))
  if (length(unique(n_k)) != 1L)
    .dynhr_abort("as_draws_array(): the chains have unequal lengths (",
                 paste(n_k, collapse = ", "), "); a draws_array needs equal ",
                 "lengths. Use as_draws_df() instead.",
                 class = "dynhr_error_unequal_chains")
  mats <- d$chains
  if (!is.null(d$logpost))
    mats <- lapply(seq_along(mats), function(k)
      cbind(mats[[k]], lp__ = d$logpost[[k]]))
  vars <- colnames(mats[[1L]])
  arr  <- array(NA_real_, dim = c(n_k[1L], length(mats), length(vars)),
                dimnames = list(iteration = NULL, chain = NULL, variable = vars))
  for (k in seq_along(mats)) arr[, k, ] <- mats[[k]]
  out <- posterior::as_draws_array(arr)
  if (!is.null(d$weights))
    out <- posterior::weight_draws(out, weights = d$weights)
  out
}

#' @rdname dynhr_interop
#' @usage \method{as_draws}{dynhr_posterior_result}(x, ...)
as_draws.dynhr_posterior_result <- as_draws.dynhr_chains

#' @rdname dynhr_interop
#' @usage \method{as_draws_df}{dynhr_posterior_result}(x, ...)
as_draws_df.dynhr_posterior_result <- as_draws_df.dynhr_chains

#' @rdname dynhr_interop
#' @usage \method{as_draws_array}{dynhr_posterior_result}(x, ...)
as_draws_array.dynhr_posterior_result <- as_draws_array.dynhr_chains


# ============================================================================
# coda: as.mcmc(), as.mcmc.list()
# ============================================================================

#' @rdname dynhr_interop
#' @usage \method{as.mcmc.list}{dynhr_chains}(x, ...)
as.mcmc.list.dynhr_chains <- function(x, ...) {
  d <- .interop_draws(x)
  if (!is.null(d$weights))
    .dynhr_abort("coda has no representation for weighted draws: this SMC ",
                 "cloud carries non-uniform $smc_weights. Resample it first ",
                 "(or use posterior::as_draws(), which keeps the weights).",
                 class = "dynhr_error_weighted_draws")
  ## coda::mcmc() / coda::mcmc.list() build exactly these structures; they
  ## are constructed here so the coercion does not need coda attached.
  structure(lapply(d$chains, function(m) {
    storage.mode(m) <- "double"
    structure(m, mcpar = c(1, nrow(m), 1), class = "mcmc")
  }), class = "mcmc.list")
}

#' @rdname dynhr_interop
#' @usage \method{as.mcmc}{dynhr_chains}(x, ...)
as.mcmc.dynhr_chains <- function(x, ...) {
  ml <- as.mcmc.list.dynhr_chains(x)
  if (length(ml) != 1L)
    .dynhr_abort("as.mcmc(): this object has ", length(ml), " chains; ",
                 "use as.mcmc.list() to keep them apart.",
                 class = "dynhr_error_multiple_chains")
  ml[[1L]]
}

#' @rdname dynhr_interop
#' @usage \method{as.mcmc.list}{dynhr_posterior_result}(x, ...)
as.mcmc.list.dynhr_posterior_result <- as.mcmc.list.dynhr_chains

#' @rdname dynhr_interop
#' @usage \method{as.mcmc}{dynhr_posterior_result}(x, ...)
as.mcmc.dynhr_posterior_result <- as.mcmc.dynhr_chains


# ============================================================================
# base / stats generics on dynhr_chains
# ============================================================================

#' @rdname dynhr_interop
#' @export
as.data.frame.dynhr_chains <- function(x, row.names = NULL, optional = FALSE,
                                       ...) {
  d   <- .interop_draws(x)
  lng <- .interop_long(d, "logpost")
  par <- setdiff(names(lng), c(".chain", ".iteration"))
  out <- lng[, c(".chain", ".iteration", par), drop = FALSE]
  if (!is.null(d$weights)) out$.weight <- d$weights
  rownames(out) <- NULL
  out
}

#' @rdname dynhr_interop
#' @export
coef.dynhr_chains <- function(object, stat = c("mean", "median"), ...) {
  stat <- match.arg(stat)
  d <- .interop_draws(object)
  X <- do.call(rbind, d$chains)
  w <- d$weights %||% rep(1 / nrow(X), nrow(X))
  out <- if (stat == "mean") {
    colSums(X * w)
  } else {
    ## Unweighted: stats::median(). Weighted: the smallest draw whose
    ## cumulative normalised weight reaches 1/2.
    apply(X, 2L, function(v) {
      if (is.null(d$weights)) return(stats::median(v))
      o <- order(v)
      v[o][which(cumsum(w[o]) >= 0.5 - 1e-12)[1L]]
    })
  }
  setNames(as.numeric(out), colnames(X))
}

#' @rdname dynhr_interop
#' @importFrom stats vcov
#' @export
vcov.dynhr_chains <- function(object, ...) {
  d <- .interop_draws(object)
  X <- do.call(rbind, d$chains)
  V <- if (is.null(d$weights)) stats::cov(X) else
    stats::cov.wt(X, wt = d$weights, method = "unbiased")$cov
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

#' @rdname dynhr_interop
#' @importFrom stats nobs
#' @export
nobs.dynhr_chains <- function(object, ...) {
  n <- object$nobs %||% object$n_obs
  if (is.null(n)) {
    .dynhr_inform("nobs(): this dynhr_chains object does not record the ",
                  "number of observations; returning NA.")
    return(NA_integer_)
  }
  as.integer(n)
}


# ============================================================================
# stats generics on dynhr_mode_result
# ============================================================================

#' Point estimate, covariance and likelihood at the posterior mode
#'
#' Base-R accessors for a \code{dynhr_mode_result} from
#' \code{\link{run_mode_finding}}.
#'
#' \code{vcov()} returns the inverse-Hessian covariance exactly as
#' \code{run_mode_finding()} builds it for the RWMH proposal (the
#' \code{.make_pd} repair of an ill-conditioned Hessian and the eigen-direction
#' cap at the prior scale included), before the \eqn{2.38^2 / n_p} proposal
#' scaling: \code{$V_mode} when stored, otherwise derived from the exact
#' analytic Hessian (\code{run_mode_finding(use_exact_hessian = TRUE)}). The
#' finite-difference Hessian of the default path is not kept in the result,
#' and \code{$Sigma_prop} cannot be un-scaled unambiguously (the
#' \code{proposal_cov_method = "full"} path stores it unscaled), so
#' \code{vcov()} errors rather than guess in that case.
#'
#' \code{logLik()} evaluates the stored log-posterior closure once at the
#' mode and returns its \code{$loglik} (the likelihood, without the prior),
#' with \code{df} = the number of estimated parameters and \code{nobs} = the
#' number of time periods in \code{$data}. For a particle likelihood this is
#' one (noisy) filter evaluation.
#'
#' @param object A \code{dynhr_mode_result}.
#' @param ... Ignored.
#' @return \code{coef()}: the named mode vector; \code{vcov()}: an
#'   \eqn{n_p \times n_p} covariance matrix; \code{logLik()}: a
#'   \code{"logLik"} object; \code{nobs()}: the number of time periods
#'   (rows of \code{$data}).
#' @seealso \code{\link{run_mode_finding}}, \code{\link{dynhr_interop}}
#' @name coef.dynhr_mode_result
NULL

#' @rdname coef.dynhr_mode_result
#' @export
coef.dynhr_mode_result <- function(object, ...) {
  object$theta_mode
}

#' @rdname coef.dynhr_mode_result
#' @export
vcov.dynhr_mode_result <- function(object, ...) {
  V <- object$V_mode
  if (is.null(V) && !is.null(object$hessian_exact))
    V <- .proposal_cov_from_hessian(object$hessian_exact, object$prior_spec,
                                    object$theta_mode, verbose = FALSE)$V_mode
  if (is.null(V))
    .dynhr_abort("vcov(): this mode result stores no inverse-Hessian ",
                 "covariance ($V_mode is NULL and there is no $hessian_exact). ",
                 "Re-run run_mode_finding(use_exact_hessian = TRUE), or use ",
                 "posterior draws.", class = "dynhr_error_no_vcov")
  nm <- names(object$theta_mode)
  dimnames(V) <- list(nm, nm)
  V
}

#' @rdname coef.dynhr_mode_result
#' @export
nobs.dynhr_mode_result <- function(object, ...) {
  if (is.null(object$data)) return(NA_integer_)
  NROW(object$data)
}

#' @rdname coef.dynhr_mode_result
#' @importFrom stats logLik
#' @export
logLik.dynhr_mode_result <- function(object, ...) {
  if (!is.function(object$log_post_fn))
    .dynhr_abort("logLik(): this mode result has no log-posterior closure ",
                 "($log_post_fn) to evaluate at the mode.",
                 class = "dynhr_error_no_loglik")
  ll <- object$log_post_fn(object$theta_mode)$loglik
  structure(as.numeric(ll), df = length(object$theta_mode),
            nobs = nobs.dynhr_mode_result(object), class = "logLik")
}
