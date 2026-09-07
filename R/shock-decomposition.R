## R/shock-decomposition.R
## --------------------------------------------------------------------------
## Shock-decomposition drivers built on top of historical_decomposition()
## (R/smoother-monolith.R) and historical_decomposition_obc()
## (R/obc-decomposition.R).
##
## Provides:
##   smoothed_initial_state()  -- s_{0|T} from a kalman_smoother() result
##   realtime_decomposition()  -- expanding-window (vintage) decomposition
##   plot.dynhr_shock_decomposition() -- stacked contribution bars
##   .resolve_shock_groups() / .group_contributions()  (internal)
##
## WHY s_{0|T} MATTERS.  The per-shock recursion in historical_decomposition()
## starts every shock's state at zero, so the shock columns alone reconstruct
## the smoothed path only when the smoothed initial state happens to equal the
## steady state.  With a non-trivial s_{0|T} the missing piece is a
## deterministic, shock-free trajectory
##
##   s^{(0)}_t = T^t s_{0|T},   x^{(0)}_t = ghx s^{(0)}_{t-1},
##
## which is exactly the "initial" component.  Adding it makes
##
##   sum_j x^{(j)}_t + x^{(0)}_t  ==  ghx s_{t-1|T} + ghu eps_{t|T}
##
## an identity (to smoother round-off), because both sides then satisfy the
## same linear recursion from the same starting point.
## --------------------------------------------------------------------------


# ---------------------------------------------------------------------------
#' Smoothed initial state \eqn{s_{0|T}} from a Kalman-smoother result
#'
#' \code{kalman_smoother()} returns \eqn{s_{t|T}} for \eqn{t = 1, \dots, T};
#' the historical decomposition also needs the smoothed \emph{pre-sample}
#' state \eqn{s_{0|T}} for its initial-condition component.
#'
#' The smoother's Durbin-Koopman backward pass now produces \eqn{s_{0|T}}
#' directly -- it is the \eqn{t = 1} entry of the same
#' \eqn{s_{t-1|T} = s_{t-1|t-1} + P_{t-1|t-1} r_{t-1}} recursion -- so this
#' function returns \code{smoother$smoothed_initial} when that field is
#' present.
#'
#' \strong{Legacy fallback.}  For a smoother result predating that field this
#' falls back to the Rauch-Tung-Striebel step at \eqn{t = 0},
#'
#' \deqn{s_{0|T} = s_{0|0} + J_0 (s_{1|T} - s_{1|0}), \qquad
#'       J_0 = P_{0|0} T' P_{1|0}^{-1},}
#'
#' using the smoother's own initialisation \eqn{s_{0|0} = 0} (hence
#' \eqn{s_{1|0} = T s_{0|0} = 0}) and \eqn{P_{0|0}} = the unconditional state
#' covariance (or the diffuse fallback when the Lyapunov solve fails, exactly
#' as \code{kalman_smoother()} does).  That step is NOT exact under dynhr's
#' lag-1 observation timing (\eqn{y_1 = Z s_0 + D \varepsilon_1} loads
#' \eqn{s_0} directly, which is what RTS assumes away), so it is only a
#' back-compatibility path.
#'
#' @param smoother A list returned by \code{\link{kalman_smoother}} (uses
#'   \code{$smoothed_initial}; falls back to \code{$smoothed_states} and
#'   \code{$predicted_cov}).
#' @param ss       State-space list from \code{\link{build_dsge_state_space}}.
#' @param Q        Shock covariance; default \code{NULL} uses
#'   \code{ss$Sigma_e} (identity only for hand-built \code{ss} lists that
#'   carry no \code{Sigma_e}), matching \code{kalman_smoother()}'s default.
#'   Used by the legacy fallback only.
#' @return Numeric vector of length \code{ss$n_state}, in \code{ss$state_names}
#'   order.
#' @export
# ---------------------------------------------------------------------------
smoothed_initial_state <- function(smoother, ss, Q = NULL) {

  if (!is.list(smoother) || is.null(smoother$smoothed_states))
    stop("smoothed_initial_state: `smoother` must be a kalman_smoother() ",
         "result with a $smoothed_states matrix.", call. = FALSE)

  ## `pre_sample` FIRST: with a backfill the smoother returns its series
  ## trimmed to the caller's sample, but `smoothed_initial` is still s_{0|T}
  ## for the PADDED one -- the state k periods earlier. The anchor the
  ## decomposition needs is the period immediately before the returned rows,
  ## which is the LAST pre-sample row (they are chronological). Using
  ## smoothed_initial there put the initial-condition trajectory k periods out
  ## of phase and broke the adding-up (measured 3.8e-2 at k = 2 on the
  ## two-shock fixture), with every dimension still correct.
  if (!is.null(smoother$presample_states) &&
      nrow(smoother$presample_states) > 0L) {
    out <- as.numeric(smoother$presample_states[nrow(smoother$presample_states), ])
    names(out) <- ss$state_names
    return(out)
  }

  ## Exact route: the DK backward pass already computed s_{0|T}.
  if (!is.null(smoother$smoothed_initial)) {
    out <- as.numeric(smoother$smoothed_initial)
    names(out) <- ss$state_names
    return(out)
  }

  if (is.null(smoother$predicted_cov))
    stop("smoothed_initial_state: `smoother` carries no $smoothed_initial ",
         "and no $predicted_cov; s_{0|T} cannot be recovered. ",
         "Pass `s0` explicitly instead.", call. = FALSE)

  T_mat <- ss$T_mat
  R_mat <- ss$R_mat
  n_s   <- nrow(T_mat)

  if (is.null(Q)) Q <- ss$Sigma_e
  if (is.null(Q)) Q <- diag(ncol(R_mat))

  ## P_{0|0}: same initialisation (and same fallback) as kalman_smoother().
  RQR  <- R_mat %*% Q %*% t(R_mat)
  P_00 <- solve_lyapunov(T_mat, RQR)
  if (anyNA(P_00)) P_00 <- .DIFFUSE_SCALE * diag(n_s)

  P_10 <- smoother$predicted_cov[, , 1L]
  J_0  <- P_00 %*% t(T_mat) %*% MASS::ginv(P_10)

  ## s_{0|0} = 0 and s_{1|0} = T s_{0|0} = 0, so the RTS step collapses to
  ## s_{0|T} = J_0 s_{1|T}.
  out <- as.numeric(J_0 %*% smoother$smoothed_states[1L, ])
  names(out) <- ss$state_names
  out
}


# ---------------------------------------------------------------------------
# Shock groups
# ---------------------------------------------------------------------------

#' Normalise and validate a shock_groups specification.
#'
#' @param shock_groups Explicit argument (wins), or NULL.
#' @param model        Parsed model; \code{model$shock_groups} is the fallback.
#' @param shock_names  Character vector of the model's shock names.
#' @return A named list of character vectors, or NULL when no grouping applies.
#' @noRd
.resolve_shock_groups <- function(shock_groups, model, shock_names) {

  if (is.null(shock_groups) && !is.null(model)) shock_groups <- model$shock_groups
  ## An EMPTY grouping means "no grouping", not a malformed one.  parse_mod()
  ## always populates model$shock_groups now (list() when the .mod carries no
  ## shock_groups block), so length 0 has to behave exactly like NULL here --
  ## otherwise every model without the block would fail the named-list check
  ## below.  (The field is always present, never absent: `$` partial matching
  ## would otherwise resolve m$shock_groups to m$shock_groups_blocks.)
  if (is.null(shock_groups) || length(shock_groups) == 0L) return(NULL)

  if (!is.list(shock_groups) || is.null(names(shock_groups)) ||
      any(!nzchar(names(shock_groups))))
    stop("shock_groups must be a NAMED list, e.g. ",
         "list(supply = c(\"e_a\", \"e_z\"), demand = \"e_g\").", call. = FALSE)

  shock_groups <- lapply(shock_groups, as.character)

  bad_name <- intersect(names(shock_groups), c("other", "initial", "constraint"))
  if (length(bad_name))
    stop("shock_groups: the group name(s) ", paste(bad_name, collapse = ", "),
         " are reserved by the decomposition.", call. = FALSE)

  if (anyDuplicated(names(shock_groups)))
    stop("shock_groups: duplicated group name(s): ",
         paste(unique(names(shock_groups)[duplicated(names(shock_groups))]),
               collapse = ", "), call. = FALSE)

  members <- unlist(shock_groups, use.names = FALSE)
  unknown <- setdiff(members, shock_names)
  if (length(unknown))
    stop("shock_groups references unknown shock(s): ",
         paste(unknown, collapse = ", "), ". Model shocks: ",
         paste(shock_names, collapse = ", "), call. = FALSE)

  dup <- unique(members[duplicated(members)])
  if (length(dup))
    stop("shock_groups: shock(s) ", paste(dup, collapse = ", "),
         " appear in more than one group; groups must partition the shocks ",
         "(adding-up would otherwise double-count them).", call. = FALSE)

  shock_groups
}

#' Collapse per-shock contribution matrices into group sums.
#'
#' Ungrouped shocks land in a single \code{"other"} column, which is omitted
#' when every shock is grouped. Works for either orientation (the matrices are
#' only ever added elementwise).
#'
#' @noRd
.group_contributions <- function(contributions, groups) {

  shock_names <- names(contributions)
  zero <- contributions[[1L]] * 0

  out <- lapply(groups, function(g) {
    Reduce(`+`, contributions[g], accumulate = FALSE)
  })
  names(out) <- names(groups)

  ungrouped <- setdiff(shock_names, unlist(groups, use.names = FALSE))
  if (length(ungrouped))
    out$other <- Reduce(`+`, contributions[ungrouped])

  ## Preserve dimnames from the per-shock matrices.
  out <- lapply(out, function(M) { dimnames(M) <- dimnames(zero); M })
  out
}


# ---------------------------------------------------------------------------
#' Real-time (vintage) shock decomposition
#'
#' Re-runs the Kalman smoother on each expanding-window vintage
#' \code{data[1:v, ]} and records the decomposition of that vintage's LAST
#' observation, i.e. period \code{v}. This is Dynare's
#' \code{realtime_shock_decomposition} object: it answers "what did we think
#' drove period \code{v} when \code{v} was the end of the sample?", which is
#' NOT the same as the full-sample answer, because later data revise the
#' smoothed shocks.
#'
#' At \code{v = nrow(data)} the row is, by construction, the full-sample
#' decomposition of that period (same smoother, same data, bit for bit).
#'
#' @param model     Parsed model object (from \code{\link{parse_mod}}).
#' @param dr        Decision rules (from the perturbation solver).
#' @param data      \code{T x n_obs} matrix or data frame of observables in
#'   \strong{levels} (the model's steady state is subtracted internally, as in
#'   \code{\link{kalman_filter}()}). When it carries column names,
#'   \code{obs_vars} selects and orders the columns.
#' @param obs_vars  Character vector of observable variable names.
#' @param vintages  Integer vector of sample ends (each in \code{2:nrow(data)}).
#' @param shock_groups Optional named list grouping shocks, e.g.
#'   \code{list(supply = c("e_a", "e_z"), demand = "e_g")}. Defaults to
#'   \code{model$shock_groups} when the model carries one.
#' @param params    Named parameter vector for the shock covariance
#'   (default \code{model$param_values}).
#' @param Q         Optional shock covariance override passed to
#'   \code{\link{kalman_smoother}}.
#' @param verbose   Print state-space construction diagnostics.
#' @return Object of class \code{dynhr_realtime_decomposition}:
#'   \describe{
#'     \item{decomposition}{\code{n_vintage x n_component x n_endo} array.}
#'     \item{total}{\code{n_vintage x n_endo} matrix (sum over components).}
#'     \item{vintages}{The integer vintage ends.}
#'     \item{components}{Component names (shocks or groups, plus
#'       \code{"initial"}).}
#'   }
#' @export
# ---------------------------------------------------------------------------
realtime_decomposition <- function(model, dr, data, obs_vars, vintages,
                                   shock_groups = NULL,
                                   params = model$param_values,
                                   Q = NULL, verbose = FALSE) {

  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data)) data <- as.matrix(data)
  if (!is.null(colnames(data)) && all(obs_vars %in% colnames(data)))
    data <- data[, obs_vars, drop = FALSE]
  if (ncol(data) != length(obs_vars))
    stop(sprintf(
      "realtime_decomposition: data has %d column(s) but obs_vars names %d.",
      ncol(data), length(obs_vars)), call. = FALSE)

  n_T <- nrow(data)
  vintages <- as.integer(vintages)
  if (!length(vintages)) stop("realtime_decomposition: `vintages` is empty.",
                              call. = FALSE)
  if (anyNA(vintages) || any(vintages < 2L) || any(vintages > n_T))
    stop(sprintf(
      "realtime_decomposition: `vintages` must be integers in 2:%d; got %s.",
      n_T, paste(vintages, collapse = ", ")), call. = FALSE)
  vintages <- sort(unique(vintages))

  ss <- build_dsge_state_space(model, dr, obs_vars, verbose = verbose,
                               params = params)

  groups <- .resolve_shock_groups(shock_groups, model, ss$shock_names)

  per_vintage <- lapply(vintages, function(v) {
    sm <- .kalman_smoother_ss(data[seq_len(v), , drop = FALSE], ss, Q = Q)
    hd <- historical_decomposition(sm, ss, shock_groups = groups)
    ## Row v of every component: the decomposition of the vintage's LAST obs.
    list(row   = vapply(hd$contributions, function(C) C[v, ], numeric(ss$n_endo)),
         total = hd$total[v, ])
  })

  comp_names <- colnames(per_vintage[[1L]]$row)
  n_comp     <- length(comp_names)

  dec <- array(0, dim = c(length(vintages), n_comp, ss$n_endo),
               dimnames = list(as.character(vintages), comp_names,
                               ss$endo_names))
  tot <- matrix(0, length(vintages), ss$n_endo,
                dimnames = list(as.character(vintages), ss$endo_names))
  for (i in seq_along(vintages)) {
    ## per_vintage[[i]]$row is n_endo x n_comp; the array wants comp x endo.
    dec[i, , ] <- t(per_vintage[[i]]$row)
    tot[i, ]   <- per_vintage[[i]]$total
  }

  structure(
    list(decomposition = dec,
         total         = tot,
         vintages      = vintages,
         components    = comp_names,
         endo_names    = ss$endo_names,
         obs_vars      = obs_vars,
         shock_groups  = groups),
    class = c("dynhr_realtime_decomposition", "list")
  )
}


# ---------------------------------------------------------------------------
#' Plot a historical shock decomposition
#'
#' Stacked contribution bars (positive contributions above zero, negative
#' below) with the reconstructed total overlaid as a line. The
#' initial-condition component is drawn like any other component but is given
#' a neutral grey so it reads as "not a shock".
#'
#' Handles both orientations: the linear decomposition stores
#' \code{T x n_endo} matrices, the OBC one \code{n_endo x T}.
#'
#' @param x        A \code{dynhr_shock_decomposition} object.
#' @param variable Name (or index) of the endogenous variable to plot;
#'   default: the first one with non-zero variation.
#' @param dates    Optional vector of period labels (length T).
#' @param col      Optional colour vector, one per component.
#' @param legend   Draw a legend (default \code{TRUE}).
#' @param main     Plot title.
#' @param ...      Passed to \code{barplot()}.
#' @return Invisibly, the \code{n_component x T} contribution matrix plotted.
#' @export
# ---------------------------------------------------------------------------
plot.dynhr_shock_decomposition <- function(x, variable = NULL, dates = NULL,
                                           col = NULL, legend = TRUE,
                                           main = NULL, ...) {

  by_time_rows <- identical(x$orientation, "time_endo")
  vnames <- if (by_time_rows) colnames(x$total) else rownames(x$total)
  if (is.null(vnames))
    vnames <- paste0("var", seq_len(if (by_time_rows) ncol(x$total)
                                    else nrow(x$total)))

  if (is.null(variable)) {
    spread <- if (by_time_rows) apply(x$total, 2L, stats::var)
              else apply(x$total, 1L, stats::var)
    variable <- vnames[which.max(spread)]
  }
  vi <- if (is.character(variable)) match(variable[1L], vnames) else as.integer(variable)
  if (is.na(vi))
    stop("plot.dynhr_shock_decomposition: variable '", variable[1L],
         "' not found.", call. = FALSE)

  ## comp: n_component x T
  comp <- t(vapply(x$contributions,
                   function(C) if (by_time_rows) C[, vi] else C[vi, ],
                   numeric(if (by_time_rows) nrow(x$total) else ncol(x$total))))
  rownames(comp) <- names(x$contributions)
  total <- if (by_time_rows) x$total[, vi] else x$total[vi, ]
  n_t   <- ncol(comp)

  if (is.null(dates)) dates <- seq_len(n_t)
  if (is.null(col)) {
    n_shock_col <- nrow(comp) - 1L
    col <- c(grDevices::hcl.colors(max(n_shock_col, 1L), "Dark 3"), "grey70")
    col <- col[seq_len(nrow(comp))]
    ## "initial" always gets the neutral grey, wherever it sits.
    ini <- match("initial", rownames(comp))
    if (!is.na(ini)) col[ini] <- "grey70"
  }
  if (is.null(main))
    main <- sprintf("Historical decomposition: %s", vnames[vi])

  pos <- pmax(comp, 0); neg <- pmin(comp, 0)
  ylim <- range(c(0, colSums(pos), colSums(neg), total), finite = TRUE)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (isTRUE(legend)) graphics::par(mar = c(4.5, 4.5, 3, 8) + 0.1)

  bp <- graphics::barplot(pos, col = col, border = NA, space = 0, ylim = ylim,
                          names.arg = dates, main = main, ylab = vnames[vi],
                          xlab = "period", ...)
  graphics::barplot(neg, col = col, border = NA, space = 0, add = TRUE,
                    axes = FALSE, axisnames = FALSE)
  graphics::abline(h = 0, col = "black", lwd = 1)
  graphics::lines(bp, total, lwd = 2, col = "black")

  if (isTRUE(legend)) {
    graphics::par(xpd = TRUE)
    graphics::legend("topright", inset = c(-0.22, 0), legend = rownames(comp),
                     fill = col, border = NA, bty = "n", cex = 0.8)
  }

  invisible(comp)
}
