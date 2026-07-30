## R/stochastic-volatility.R
## --------------------------------------------------------------------------
## Call-level constructors and internal helpers for the stochastic_volatility
## block: latent per-shock AR(1) log-variance (Justiniano-Primiceri 2008;
## Kim-Shephard-Chib 1998), the estimable/latent analogue of the deterministic
## heteroskedastic_shocks schedule.
##
## Syntax in .mod files:
##   stochastic_volatility;
##     var e_a; mu = mu_a; rho = rho_a; sigma_eta = sig_eta_a;
##     var e_r; mu = mu_r; rho = rho_r; sigma_eta = sig_eta_r;
##   end;
##
## Semantics: shock i has time-varying standard deviation
##   stderr_t(e_i) = sigma_i * exp(h_{i,t} / 2),
##   h_{i,t} = (1 - rho_i) mu_i + rho_i h_{i,t-1} + sigma_eta_i eta_{i,t},
##   eta_{i,t} ~ N(0, 1)  (marginalized by the RB particle filter),
## where sigma_i is the BASELINE shock stderr from the shocks block. Unlike
## heteroskedastic_shocks (a literal exogenous `scales` schedule), the volatility
## path h is a LATENT state integrated out by the filter, and mu/rho/sigma_eta
## are (typically estimated) hyperparameters referenced by name.
##
## Structure mirrors R/heteroskedastic-shocks.R.
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 1.  Constructors: stochastic_volatility() and sv_entry()
## ---------------------------------------------------------------------------

#' Specify a single stochastic-volatility entry
#'
#' Helper for \code{\link{stochastic_volatility}()} defining the latent AR(1)
#' log-variance for one exogenous shock. Each hyperparameter (\code{mu},
#' \code{rho}, \code{sigma_eta}) is either a numeric literal (a fixed value) or
#' a character string naming a model parameter — usually one listed in
#' \code{estimated_params} so that the SV hyperparameter is itself estimated.
#'
#' @param shock     Name of the exogenous shock (character scalar).
#' @param mu        Log-level of the volatility process: the unconditional mean
#'   of \eqn{h_{i,t}}. Numeric or a parameter-name string.
#' @param rho       AR(1) persistence of \eqn{h_{i,t}}, in \eqn{(-1, 1)} for
#'   stationarity. Numeric or a parameter-name string.
#' @param sigma_eta Standard deviation of the log-variance innovation
#'   \eqn{\eta_{i,t}} (the "vol-of-vol"); must be positive. Numeric or a
#'   parameter-name string.
#' @return An \code{"sv_entry"} list for use inside
#'   \code{\link{stochastic_volatility}()}.
#' @seealso \code{\link{stochastic_volatility}}
#' @examples
#' \dontrun{
#' sv_entry("e_a", mu = "mu_a", rho = "rho_a", sigma_eta = "sig_eta_a")
#' sv_entry("e_r", mu = 0, rho = 0.9, sigma_eta = 0.2)   # fixed hyperparameters
#' }
#' @export
sv_entry <- function(shock, mu, rho, sigma_eta) {
  if (!is.character(shock) || length(shock) != 1L || nchar(shock) == 0L)
    stop("sv_entry: 'shock' must be a non-empty character scalar.", call. = FALSE)

  chk <- function(x, nm) {
    if (length(x) != 1L)
      stop(sprintf("sv_entry: '%s' for shock '%s' must be a scalar.", nm, shock),
           call. = FALSE)
    if (is.character(x)) {
      if (nchar(x) == 0L)
        stop(sprintf("sv_entry: '%s' for shock '%s' is an empty string.", nm, shock),
             call. = FALSE)
      return(x)
    }
    if (is.numeric(x) && is.finite(x)) return(as.numeric(x))
    stop(sprintf(
      "sv_entry: '%s' for shock '%s' must be a finite numeric or a parameter name.",
      nm, shock), call. = FALSE)
  }

  entry <- list(shock = shock,
                mu        = chk(mu,        "mu"),
                rho       = chk(rho,       "rho"),
                sigma_eta = chk(sigma_eta, "sigma_eta"))
  class(entry) <- "sv_entry"
  entry
}


#' Specify call-level stochastic volatility for a DSGE model
#'
#' Constructs an \code{sv_spec} object — the same structure produced by a
#' \code{stochastic_volatility ... end;} block in a \code{.mod} file — so that
#' latent per-shock stochastic volatility can be supplied or overridden at
#' estimation time without editing the model file.
#'
#' Pass the result to the \code{stochastic_volatility} argument of
#' \code{\link{run_full_estimation}()} / \code{\link{make_log_posterior_sv_rbpf}()}.
#' Pass \code{FALSE} to disable a mod-file block.
#'
#' @param ... One or more entries created by \code{\link{sv_entry}()}.
#' @return An object of class \code{"sv_spec"}: a list with one element
#'   \code{$sv}, a data.frame with columns \code{shock} (character) and list
#'   columns \code{mu}, \code{rho}, \code{sigma_eta} (each element a numeric
#'   literal or a parameter-name string).
#' @seealso \code{\link{sv_entry}}, \code{\link{make_log_posterior_sv_rbpf}}
#' @examples
#' \dontrun{
#' sv <- stochastic_volatility(
#'   sv_entry("e_a", mu = "mu_a", rho = "rho_a", sigma_eta = "sig_eta_a"),
#'   sv_entry("e_r", mu = "mu_r", rho = "rho_r", sigma_eta = "sig_eta_r"))
#' }
#' @export
stochastic_volatility <- function(...) {
  entries <- list(...)
  ## Flatten a single list-of-entries argument (same idiom as heteroskedastic_shocks).
  if (length(entries) == 1L && is.list(entries[[1L]]) &&
      !inherits(entries[[1L]], "sv_entry"))
    entries <- entries[[1L]]

  for (i in seq_along(entries)) {
    if (!inherits(entries[[i]], "sv_entry"))
      stop(sprintf(
        "stochastic_volatility: argument %d is not an sv_entry() (got %s).",
        i, class(entries[[i]])[1L]), call. = FALSE)
  }

  if (length(entries) == 0L) {
    empty <- data.frame(shock = character(0L), stringsAsFactors = FALSE)
    empty$mu <- list(); empty$rho <- list(); empty$sigma_eta <- list()
    out <- list(sv = empty)
    class(out) <- "sv_spec"
    return(out)
  }

  shocks <- vapply(entries, function(e) e$shock, character(1L))
  if (anyDuplicated(shocks))
    stop(sprintf("stochastic_volatility: duplicate shock(s): %s.",
                 paste(unique(shocks[duplicated(shocks)]), collapse = ", ")),
         call. = FALSE)

  df <- data.frame(shock = shocks, stringsAsFactors = FALSE)
  df$mu        <- lapply(entries, function(e) e$mu)
  df$rho       <- lapply(entries, function(e) e$rho)
  df$sigma_eta <- lapply(entries, function(e) e$sigma_eta)

  out <- list(sv = df)
  class(out) <- "sv_spec"
  out
}


## ---------------------------------------------------------------------------
## 2.  Internal: resolve call-level stochastic_volatility arg onto the model
## ---------------------------------------------------------------------------

#' Resolve call-level stochastic_volatility argument for estimation
#'
#' @param model                 Parsed \code{dynhr_mod}.
#' @param stochastic_volatility User-supplied argument (NULL = use mod-file
#'   block; FALSE = disable; an \code{sv_spec} = override).
#' @return Possibly-modified \code{model} list.
#' @noRd
.resolve_stochastic_volatility <- function(model, stochastic_volatility) {
  if (is.null(stochastic_volatility))
    return(model)
  if (isFALSE(stochastic_volatility)) {
    empty <- data.frame(shock = character(0L), stringsAsFactors = FALSE)
    empty$mu <- list(); empty$rho <- list(); empty$sigma_eta <- list()
    model$stochastic_volatility <- list(sv = empty)
    return(model)
  }
  if (inherits(stochastic_volatility, "sv_spec")) {
    model$stochastic_volatility <- stochastic_volatility
    return(model)
  }
  stop(paste0(
    "stochastic_volatility must be NULL (use mod-file block), FALSE (disable), ",
    "or an sv_spec object from stochastic_volatility()."), call. = FALSE)
}


## ---------------------------------------------------------------------------
## 3.  Internal: shock index + hyperparameter resolution
## ---------------------------------------------------------------------------

#' Map SV entries to positions in the shock (dr$exo_names) order
#'
#' @param sv_spec   An \code{sv_spec} (\code{model$stochastic_volatility}).
#' @param exo_names Character vector of shock names in dr$exo_names order
#'   (Landmine 9: NOT model$varexo_names).
#' @return \code{NULL} if the spec is empty; otherwise an integer vector, one
#'   position per SV entry, into \code{exo_names} (in the spec's row order).
#' @noRd
.sv_shock_index <- function(sv_spec, exo_names) {
  if (is.null(sv_spec)) return(NULL)
  df <- sv_spec$sv
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  idx <- match(df$shock, exo_names)
  if (any(is.na(idx)))
    stop(sprintf(
      ".sv_shock_index: SV shock(s) not in dr$exo_names (%s): %s.",
      paste(exo_names, collapse = ", "),
      paste(df$shock[is.na(idx)], collapse = ", ")), call. = FALSE)
  idx
}


#' Resolve SV hyperparameters (mu, rho, sigma_eta) to numeric per entry
#'
#' Each of \code{mu}/\code{rho}/\code{sigma_eta} is either a numeric literal or
#' a parameter-name string looked up in the named \code{params} vector.
#'
#' @param sv_spec An \code{sv_spec}.
#' @param params  Named numeric parameter vector (after
#'   \code{.apply_theta_to_params}); character references are resolved against
#'   its names.
#' @return A numeric matrix \code{n_sv x 3} with columns
#'   \code{c("mu","rho","sigma_eta")} and rownames = shock names, OR a character
#'   scalar error tag \code{"__sv_infeasible__"} attached as
#'   \code{attr(, "reason")} when a resolved value violates a domain constraint
#'   (\eqn{|rho| < 1}, \eqn{sigma_eta > 0}) — so per-draw callers can return
#'   \code{-Inf} rather than \code{stop()}. Missing parameter NAMES are a
#'   specification error and still \code{stop()}.
#' @noRd
.sv_resolve_hyperparams <- function(sv_spec, params) {
  df <- sv_spec$sv
  n  <- nrow(df)
  out <- matrix(NA_real_, nrow = n, ncol = 3L,
                dimnames = list(df$shock, c("mu", "rho", "sigma_eta")))

  resolve1 <- function(x, nm, shock) {
    if (is.character(x)) {
      if (!x %in% names(params))
        stop(sprintf(
          paste0("stochastic_volatility: hyperparameter '%s' for shock '%s' ",
                 "names a parameter '%s' that is not in the model's parameters."),
          nm, shock, x), call. = FALSE)
      return(as.numeric(params[[x]]))
    }
    as.numeric(x)
  }

  for (i in seq_len(n)) {
    out[i, "mu"]        <- resolve1(df$mu[[i]],        "mu",        df$shock[i])
    out[i, "rho"]       <- resolve1(df$rho[[i]],       "rho",       df$shock[i])
    out[i, "sigma_eta"] <- resolve1(df$sigma_eta[[i]], "sigma_eta", df$shock[i])
  }

  ## Domain constraints -> soft-infeasible (return -Inf upstream), not stop().
  if (any(!is.finite(out)) ||
      any(abs(out[, "rho"]) >= 1) ||
      any(out[, "sigma_eta"] <= 0)) {
    bad <- structure("__sv_infeasible__",
                     reason = "rho must be in (-1,1) and sigma_eta > 0")
    return(bad)
  }
  out
}
