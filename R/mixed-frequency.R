## R/mixed-frequency.R
## ---------------------------------------------------------------------------
## Mixed-frequency observation blocks: temporal aggregation of a
## higher-frequency latent into a lower-frequency observable.
##
## THE GAP THIS CLOSES.  Ragged edges were already handled by row-dropping in
## the filter's missing-data branch, but that branch keeps a CONSTANT `ZZ`:
## every observable is a contemporaneous function of the state.  A quarterly
## flow built out of three monthly latents is not of that form, so nowcasting
## from a monthly model was out of reach.
##
## THE FORM.  Declare, per observable,
##
##   obs_aggregation = list(
##     gdp_q = list(of = "gdp_m", type = "flow_sum", k = 3L)
##   )
##
## meaning "the observable `gdp_q` is the k-period aggregate of the model
## variable `gdp_m`".  The model itself stays at the HIGH frequency; the
## aggregate is observed only every k-th period and is `NA` in between (the
## existing missing-data path handles those periods, unchanged).
##
## THE IMPLEMENTATION: fixed-weight state augmentation (Harvey 1989 §6.3;
## Mariano & Murasawa 2003).  Every aggregator supported here is a fixed
## linear filter of the underlying latent,
##
##   y^agg_t = w_1 x_t + w_2 x_{t-1} + ... + w_m x_{t-m+1},
##
## so instead of a periodically-reset accumulator (which needs a TIME-VARYING
## transition the recursions cannot express) we carry the m-1 lags of x as
## extra states.  With
##
##   x_t = Zb s_{t-1} + Db eps_t            (dynhr's lagged-state convention)
##
## the augmented system z_t = [s_t; x_t; x_{t-1}; ...; x_{t-m+2}] is
##
##   z_t = [ TT  0 ] z_{t-1} + [ RR ] eps_t
##         [ Zb  0 ]           [ Db ]
##         [ 0 I 0 ]           [ 0  ]
##
##   y^agg_t = [ w_1 Zb | w_2 ... w_m ] z_{t-1} + w_1 Db eps_t.
##
## `ZZ` and `TT` are therefore still CONSTANT: the filter's hot path, its
## steady-state lock and both C++ kernels are untouched, and a model with no
## `obs_aggregation` produces byte-identical log-likelihoods (pinned in
## tests/testthat/test-mixed-frequency.R).  The price is m-1 extra states per
## aggregated observable, which is the standard Mariano-Murasawa cost.
##
## The m = 1 case (`k = 1`, or `stock_end` at any k) adds NO states and the
## observation row is bit-identical to the plain filter's.
##
## WEIGHTS (`type`), for aggregation length k:
##   flow_sum   w = (1, 1, ..., 1)                    -- k-period sum of a flow
##   flow_mean  w = (1/k, ..., 1/k)                   -- k-period average
##   stock_end  w = (1)                               -- end-of-period stock
##   triangle   w = (1, 2, ..., k, k-1, ..., 1) / k   -- Mariano-Murasawa
##
## `triangle` is the Mariano-Murasawa (2003) approximation used when the
## monthly latent is a GROWTH RATE and the quarterly observable is the growth
## rate of the quarterly average: it spans 2k-1 lags and its weights sum to k,
## like `flow_sum`.
##
## STEADY STATE.  The aggregate's mean is sum(w) times the underlying
## variable's steady state, so `kalman_filter()` scales the `dr$ys` offset by
## sum(w) for aggregated rows.  Data must be supplied on the aggregate's own
## scale (a quarterly SUM for `flow_sum`, not a monthly average).
##
## WHAT DOES NOT SUPPORT IT.  Only the Gaussian Kalman branch consumes the
## augmentation.  Every other likelihood in `make_log_posterior()` (the
## particle filters `tpf`/`ppf`/`copf`/`sv_rbpf`/`global_pf`, the pruned and
## PSKF higher-order filters, `whittle`, `cumulant`, `student_t`) and the
## Markov-switching Kim filter evaluate the observation equation themselves
## and would SILENTLY ignore the spec, so `make_log_posterior()` errors loudly
## on them rather than returning a quietly-wrong likelihood.
##
## The ANALYTIC score (`make_posterior_grad()` and the adjoint/tangent Kalman
## gradients) differentiates its own copy of the observation equation and has
## no aggregation awareness, so a mixed-frequency posterior must be sampled
## with a gradient-free sampler (RWMH / SMC / DA) or a finite-difference
## gradient. Wiring the augmentation through the adjoint is left to a
## follow-up.
## ---------------------------------------------------------------------------


.MF_TYPES <- c("flow_sum", "flow_mean", "stock_end", "triangle")


# ---------------------------------------------------------------------------
#' Temporal-aggregation weights
#'
#' The fixed linear filter behind each \code{obs_aggregation} \code{type}:
#' the aggregate observable equals \code{sum(w[i] * x[t - i + 1])} over
#' \code{i = 1..length(w)}, where \code{x} is the higher-frequency latent.
#'
#' @param type One of \code{"flow_sum"} (k-period sum), \code{"flow_mean"}
#'   (k-period average), \code{"stock_end"} (end-of-period stock -- a single
#'   unit weight, no lags) or \code{"triangle"} (the Mariano-Murasawa
#'   log-approximation weights \code{(1, 2, ..., k, k-1, ..., 1) / k},
#'   spanning \code{2k - 1} lags).
#' @param k Aggregation length (number of high-frequency periods per
#'   low-frequency observation). \code{k = 1} returns \code{1} for every type.
#' @return Numeric vector of weights, most recent period first.
#' @examples
#' mf_aggregation_weights("flow_sum", 3)
#' mf_aggregation_weights("triangle", 3)
#' @seealso \code{\link{kalman_filter}}, \code{\link{mf_augment_state_space}}
#' @export
# ---------------------------------------------------------------------------
mf_aggregation_weights <- function(type = .MF_TYPES, k = 1L) {
  type <- match.arg(type, .MF_TYPES)
  k    <- .mf_check_k(k)
  switch(type,
         flow_sum  = rep(1, k),
         flow_mean = rep(1 / k, k),
         stock_end = 1,
         triangle  = c(seq_len(k), rev(seq_len(k - 1L))) / k)
}


#' @noRd
.mf_check_k <- function(k, what = "k") {
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) ||
      k < 1 || k != as.integer(k))
    stop(sprintf("obs_aggregation: `%s` must be a single integer >= 1.", what),
         call. = FALSE)
  as.integer(k)
}


## Validate one obs_aggregation entry and return list(of, type, k, w).
#' @noRd
.mf_parse_entry <- function(entry, nm) {
  if (!is.list(entry))
    stop(sprintf(
      "obs_aggregation[['%s']] must be a list with `of`, `type` and `k`.", nm),
      call. = FALSE)
  if (is.null(entry$of) || !is.character(entry$of) || length(entry$of) != 1L)
    stop(sprintf("obs_aggregation[['%s']]$of must be a single variable name.",
                 nm), call. = FALSE)
  type <- entry$type %||% .MF_TYPES
  if (!is.character(type) || length(type) < 1L)
    stop(sprintf("obs_aggregation[['%s']]$type must be a character.", nm),
         call. = FALSE)
  type <- tryCatch(match.arg(type, .MF_TYPES),
                   error = function(e)
                     stop(sprintf(
                       "obs_aggregation[['%s']]$type must be one of %s; got '%s'.",
                       nm, paste(sQuote(.MF_TYPES), collapse = ", "), type[1]),
                       call. = FALSE))
  k <- .mf_check_k(entry$k %||% 1L, sprintf("obs_aggregation[['%s']]$k", nm))
  list(of = entry$of, type = type, k = k,
       w = mf_aggregation_weights(type, k))
}


## Resolve an obs_aggregation spec against a vector of observable names.
##
## `obs_names`   the observables the caller asked for (aggregate names).
## `known`       names the `of` fields must resolve into (model endo names,
##               or an existing state-space's obs_names). NULL disables the
##               membership check.
##
## Returns NULL when nothing is aggregated (no spec, or every weight vector is
## the scalar 1 AND every k is 1 -- i.e. the identity aggregator), otherwise a
## list with, aligned to `obs_names`:
##   base    character  underlying variable per observable
##   w_list  list       weight vector per observable (scalar 1 when plain)
##   k       integer    aggregation length per observable (1 when plain)
##   type    character  aggregator per observable (NA when plain)
##   scale   numeric    sum(w) per observable (steady-state multiplier)
##   agg     logical    which observables are aggregated
#' @noRd
.mf_resolve <- function(spec, obs_names, known = NULL, what = "obs_aggregation") {
  if (is.null(spec) || length(spec) == 0L) return(NULL)
  if (!is.list(spec) || is.null(names(spec)) || any(!nzchar(names(spec))))
    stop(what, " must be a NAMED list, one entry per aggregated observable.",
         call. = FALSE)
  if (anyDuplicated(names(spec)))
    stop(what, ": duplicated observable name(s): ",
         paste(unique(names(spec)[duplicated(names(spec))]), collapse = ", "),
         call. = FALSE)

  n <- length(obs_names)
  out <- list(base   = obs_names,
              w_list = rep(list(1), n),
              k      = rep(1L, n),
              type   = rep(NA_character_, n),
              scale  = rep(1, n),
              agg    = rep(FALSE, n))

  for (nm in names(spec)) {
    i <- match(nm, obs_names)
    if (is.na(i))
      stop(what, ": '", nm, "' is not among the observables (",
           paste(obs_names, collapse = ", "), ").", call. = FALSE)
    e <- .mf_parse_entry(spec[[nm]], nm)
    if (!is.null(known) && !(e$of %in% known))
      stop(what, "[['", nm, "']]$of = '", e$of, "' is not available (",
           paste(known, collapse = ", "), ").", call. = FALSE)
    out$base[i]      <- e$of
    out$w_list[[i]]  <- e$w
    out$k[i]         <- e$k
    out$type[i]      <- e$type
    out$scale[i]     <- sum(e$w)
    out$agg[i]       <- TRUE
  }

  ## Identity aggregator: nothing to augment and nothing to rescale, so hand
  ## back NULL and let every caller take its untouched constant-ZZ path.
  ## `base != obs_names` still counts as active (the observable is a RENAME of
  ## a model variable and the observation row must follow the rename).
  if (all(vapply(out$w_list, function(w) length(w) == 1L && w == 1, TRUE)) &&
      identical(out$base, as.character(obs_names)))
    return(NULL)
  out
}


## Build the augmented (TT, RR, ZZ, DD) from the plain state block and the
## per-observable base rows.
##
##   TT, RR   n_state x n_state / n_state x n_shock  (state block)
##   Zb, Db   n_obs x n_state / n_obs x n_shock      (rows of the UNDERLYING
##                                                    variable per observable)
##   w_list   per-observable weight vectors
##
## Returns list(TT, RR, ZZ, DD, aug_names, n_aug).
#' @noRd
.mf_augment_matrices <- function(TT, RR, Zb, Db, w_list, obs_names) {
  n_state <- nrow(TT); n_shock <- ncol(RR); n_obs <- nrow(Zb)
  m_vec   <- vapply(w_list, length, 1L)
  n_lag   <- m_vec - 1L                       # extra states per observable
  n_aug   <- sum(n_lag)

  TT_a <- matrix(0, n_state + n_aug, n_state + n_aug)
  RR_a <- matrix(0, n_state + n_aug, n_shock)
  ZZ_a <- matrix(0, n_obs, n_state + n_aug)
  DD_a <- matrix(0, n_obs, n_shock)

  TT_a[seq_len(n_state), seq_len(n_state)] <- TT
  RR_a[seq_len(n_state), ] <- RR

  aug_names <- character(0)
  off <- n_state
  for (i in seq_len(n_obs)) {
    w <- w_list[[i]]
    ## Contemporaneous part: w_1 * (Zb s_{t-1} + Db eps_t).
    ZZ_a[i, seq_len(n_state)] <- w[1L] * Zb[i, ]
    DD_a[i, ] <- w[1L] * Db[i, ]
    if (n_lag[i] == 0L) next
    idx <- off + seq_len(n_lag[i])
    ## Lag block: first row generates x_t from s_{t-1}/eps_t, the rest shift.
    TT_a[idx[1L], seq_len(n_state)] <- Zb[i, ]
    RR_a[idx[1L], ] <- Db[i, ]
    if (n_lag[i] > 1L)
      TT_a[cbind(idx[-1L], idx[-length(idx)])] <- 1
    ## z_{t-1} carries x_{t-1}, ..., x_{t-m+1}: weights w_2 .. w_m.
    ZZ_a[i, idx] <- w[-1L]
    aug_names <- c(aug_names,
                   paste0(obs_names[i], "_L", seq_len(n_lag[i]) - 1L))
    off <- off + n_lag[i]
  }

  list(TT = TT_a, RR = RR_a, ZZ = ZZ_a, DD = DD_a,
       aug_names = aug_names, n_aug = n_aug)
}


## Warn when an aggregated column's non-NA pattern is inconsistent with its
## aggregation length: a quarterly series stored at monthly frequency must
## carry observations k periods apart (ragged edges and whole missing
## quarters are fine -- every gap is still a MULTIPLE of k). A misaligned
## column is otherwise silently filtered as if it were high-frequency.
##
## `data` is n_obs x T (the filter's orientation).
#' @noRd
.mf_check_pattern <- function(data, mf, obs_names) {
  for (i in which(mf$agg)) {
    if (mf$k[i] <= 1L) next
    obs_t <- which(!is.na(data[i, ]))
    if (length(obs_t) < 2L) next
    if (any(diff(obs_t) %% mf$k[i] != 0L))
      warning(sprintf(
        paste0("kalman_filter: observable '%s' is declared as a %d-period ",
               "aggregate but its non-missing observations are not spaced in ",
               "multiples of %d. Low-frequency observations must sit at the ",
               "LAST high-frequency period of each aggregation window, with ",
               "NA elsewhere."),
        obs_names[i], mf$k[i], mf$k[i]), call. = FALSE)
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
#' Place low-frequency observations on a high-frequency grid
#'
#' Convenience helper for building the \code{NA}-padded column an aggregated
#' observable needs: the \code{j}-th low-frequency value is written at
#' high-frequency period \code{offset + j * k}, every other period is
#' \code{NA}.  With the default \code{offset = 0} the first observation sits
#' at period \code{k}, i.e. at the END of the first complete aggregation
#' window -- which is what the weights in
#' \code{\link{mf_aggregation_weights}} assume.
#'
#' @param x Numeric vector of low-frequency observations (\code{NA} allowed).
#' @param k Aggregation length (high-frequency periods per observation).
#' @param n_T Length of the high-frequency sample. Default: just long enough
#'   to hold every element of \code{x}.
#' @param offset Number of leading high-frequency periods before the first
#'   aggregation window closes (default \code{0}).
#' @return Numeric vector of length \code{n_T}, \code{NA} except at the
#'   aggregation dates.
#' @examples
#' mf_expand_observations(c(1.2, 0.8), k = 3)
#' @seealso \code{\link{mf_aggregation_weights}}
#' @export
# ---------------------------------------------------------------------------
mf_expand_observations <- function(x, k, n_T = NULL, offset = 0L) {
  k <- .mf_check_k(k)
  if (length(offset) != 1L || !is.finite(offset) || offset < 0)
    stop("mf_expand_observations: `offset` must be a non-negative integer.",
         call. = FALSE)
  offset <- as.integer(offset)
  pos    <- offset + seq_along(x) * k
  if (is.null(n_T)) n_T <- if (length(pos)) max(pos) else offset
  n_T <- .mf_check_k(n_T, "n_T")
  if (length(pos) && max(pos) > n_T)
    stop(sprintf(
      "mf_expand_observations: %d observations at k = %d (offset %d) need ",
      length(x), k, offset),
      sprintf("n_T >= %d; got %d.", max(pos), n_T), call. = FALSE)
  out      <- rep(NA_real_, n_T)
  out[pos] <- as.numeric(x)
  out
}


# ---------------------------------------------------------------------------
#' Augment a DSGE state-space with temporal-aggregation observables
#'
#' Turns a compact state-space from \code{\link{build_dsge_state_space}} into
#' one whose observables are \code{k}-period aggregates of the higher-frequency
#' model variables, by appending the lag states the aggregator needs (see the
#' file header of \code{R/mixed-frequency.R} for the algebra).  The result is
#' an ordinary \code{dsge_ss}, so \code{\link{kalman_smoother}} smooths it
#' unchanged -- which is how a monthly latent is recovered from quarterly
#' data.
#'
#' The state block, the shock block and every non-aggregated observation row
#' are carried over unchanged, so the identity spec is a no-op.
#'
#' @param ss A \code{dsge_ss} from \code{\link{build_dsge_state_space}}, built
#'   with \code{obs_vars} naming the HIGH-frequency model variables.
#' @param obs_aggregation Named list, one entry per aggregated observable:
#'   \code{list(gdp_q = list(of = "gdp_m", type = "flow_sum", k = 3L))}. The
#'   name is the aggregate observable's name in the result; \code{of} must be
#'   one of \code{ss$obs_names}. See \code{\link{mf_aggregation_weights}} for
#'   the \code{type} menu.
#' @return A \code{dsge_ss} with \code{T_mat}, \code{R_mat}, \code{Z_mat},
#'   \code{D_mat}, \code{state_names}, \code{obs_names}, \code{n_state} and
#'   the observation intercept \code{d} (rescaled by \code{sum(w)}, since an
#'   aggregate's steady state is \code{sum(w)} times the underlying
#'   variable's) updated, plus an \code{mf} field recording the resolved spec. The
#'   \code{ghx}/\code{ghu} fields are dropped: they describe the
#'   UNaugmented model, so \code{\link{historical_decomposition}} (which
#'   consumes them) is not available on an augmented state-space.
#' @seealso \code{\link{kalman_filter}} (which does the same augmentation
#'   internally from \code{model$obs_aggregation}),
#'   \code{\link{mf_expand_observations}}
#' @export
# ---------------------------------------------------------------------------
mf_augment_state_space <- function(ss, obs_aggregation) {
  if (!is.list(ss) || is.null(ss$T_mat))
    stop("mf_augment_state_space: `ss` must be a state-space list from ",
         "build_dsge_state_space().", call. = FALSE)
  if (!is.null(ss$timing) && !identical(ss$timing, "lagged"))
    ss <- ss_convert_timing(ss)

  ## The spec is keyed by the AGGREGATE name but `of` points at an existing
  ## obs row, so resolve positionally through `of` and then rename the row.
  base_names <- as.character(ss$obs_names)
  spec_by_of <- list()
  new_names  <- base_names
  for (nm in names(obs_aggregation)) {
    e <- obs_aggregation[[nm]]
    of <- if (is.list(e)) e$of else NULL
    if (is.null(of) || !is.character(of) || length(of) != 1L)
      stop("mf_augment_state_space: obs_aggregation[['", nm,
           "']]$of must be a single variable name.", call. = FALSE)
    i <- match(of, base_names)
    if (is.na(i))
      stop("mf_augment_state_space: obs_aggregation[['", nm, "']]$of = '", of,
           "' is not an observable of `ss` (",
           paste(base_names, collapse = ", "), ").", call. = FALSE)
    spec_by_of[[of]] <- e
    new_names[i]     <- nm
  }
  mf <- .mf_resolve(spec_by_of, base_names, known = base_names,
                    what = "mf_augment_state_space: obs_aggregation")
  if (is.null(mf)) {
    ss$obs_names <- new_names
    return(ss)
  }

  aug <- .mf_augment_matrices(ss$T_mat, ss$R_mat, ss$Z_mat, ss$D_mat,
                              mf$w_list, new_names)

  ss$T_mat       <- aug$TT
  ss$R_mat       <- aug$RR
  ss$Z_mat       <- aug$ZZ
  ss$D_mat       <- aug$DD
  ss$state_names <- c(ss$state_names, aug$aug_names)
  ss$obs_names   <- new_names
  ss$n_state     <- ss$n_state + aug$n_aug
  ss$ghx         <- NULL
  ss$ghu         <- NULL
  ## An aggregate's steady state is sum(w) times the underlying variable's, so
  ## the observation intercept the smoother subtracts has to be rescaled with
  ## it -- the same rescaling kalman_filter() applies to its own `d`.
  if (!is.null(ss$d)) {
    if (length(ss$d) != length(mf$scale))
      stop(sprintf(paste0("mf_augment_state_space: the state space's ",
                          "observation intercept has %d entries but there are ",
                          "%d observables to rescale."),
                   length(ss$d), length(mf$scale)), call. = FALSE)
    ss$d <- ss$d * mf$scale
  }
  ss$mf          <- mf
  ss$mf_base     <- base_names
  ss
}
