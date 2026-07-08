## R/obc-regime-path.R
## --------------------------------------------------------------------------
## Canonical regime-path object for occasionally-binding-constraint (OBC)
## solvers.  dynhr has several independent OBC/complementarity path solvers
## (Boehl spell-duration, LCP Lemke, MCP semi-smooth Newton, nonlinear OccBin,
## PF-Newton) whose returns encode "which constraints bind when" three different
## ways -- an integer bitfield (\code{regime_path}), the same bitfield under a
## different name (\code{active_set}), and a logical matrix (\code{regime}) --
## and store the solved paths in two orientations (n_endo x T vs T x n_endo).
##
## This module provides ONE canonical object and adapters onto it.  It is
## ADDITIVE: it does not modify any solver (item A1a).  Later increments (A1b)
## may migrate solvers to emit it natively; a shared cache (A1c) can key on its
## \code{path_hash}.
## --------------------------------------------------------------------------


#' Decode an integer regime-bitfield to a logical binding matrix
#'
#' @param codes Integer vector (length T): bit \code{b-1} set means constraint
#'   \code{b} binds in that period.
#' @param n_constr Number of constraints.
#' @return A \code{T x n_constr} logical matrix.
#' @keywords internal
.obc_bitfield_to_binding <- function(codes, n_constr) {
  codes <- as.integer(codes)
  B <- matrix(FALSE, length(codes), n_constr)
  for (b in seq_len(n_constr))
    B[, b] <- bitwAnd(codes, bitwShiftL(1L, b - 1L)) != 0L
  B
}

#' Encode a logical binding matrix as an integer regime-bitfield
#' @param binding A \code{T x n_constr} logical matrix.
#' @return Integer vector of length \code{T}.
#' @keywords internal
.obc_binding_to_bitfield <- function(binding) {
  binding <- matrix(as.logical(binding), nrow = nrow(binding))
  w <- bitwShiftL(1L, seq_len(ncol(binding)) - 1L)
  as.integer(binding %*% w)
}


#' Construct a canonical OBC regime-path object
#'
#' @param binding A \code{T x n_constr} logical matrix: which constraints bind
#'   in which periods.
#' @param paths Optional \code{T x n_endo} matrix of solved endogenous paths
#'   (canonical orientation: rows = periods).
#' @param converged Logical convergence flag (or \code{NA}).
#' @param shadow_values Optional \code{T x n_constr} matrix of constraint shadow
#'   values / multipliers (\code{NULL} if the solver does not report them).
#' @param terminal_state Optional terminal state vector.
#' @param constraints Optional character vector naming the constrained variables.
#'
#' @return An object of class \code{obc_regime_path} with fields \code{binding},
#'   \code{bitfield}, \code{active_periods}, \code{n_periods},
#'   \code{n_constraints}, \code{paths}, \code{shadow_values},
#'   \code{terminal_state}, \code{converged}, \code{constraints}, and a
#'   \code{path_hash} (a deterministic key for caching).
#' @export
new_obc_regime_path <- function(binding, paths = NULL, converged = NA,
                                shadow_values = NULL, terminal_state = NULL,
                                constraints = NULL) {
  binding <- matrix(as.logical(binding), nrow = nrow(binding))
  codes <- .obc_binding_to_bitfield(binding)
  structure(list(
    binding        = binding,
    bitfield       = codes,
    active_periods = which(rowSums(binding) > 0),
    n_periods      = nrow(binding),
    n_constraints  = ncol(binding),
    constraints    = constraints,
    paths          = paths,
    shadow_values  = shadow_values,
    terminal_state = terminal_state,
    converged      = converged,
    path_hash      = paste(codes, collapse = "-")
  ), class = "obc_regime_path")
}


#' Coerce an OBC-solver return to the canonical regime-path object
#'
#' Normalizes the fragmented solver returns: the binding sequence from
#' \code{regime_path} / \code{active_set} (integer bitfield) or \code{regime}
#' (logical matrix), and the solved paths to the canonical \code{T x n_endo}
#' orientation.
#'
#' @param res A list returned by an OBC solver (Boehl/LCP/MCP/OccBin/PF-Newton).
#' @param n_constr Number of constraints (required when the binding sequence is
#'   an integer bitfield).
#' @param paths_endo_by_T Optional logical: set \code{TRUE} if \code{res$paths}
#'   is stored as \code{n_endo x T} (as in the DR-based solvers), \code{FALSE}
#'   if \code{T x n_endo}.  If \code{NULL}, inferred from the period count.
#' @param constraints Optional character vector naming the constrained variables.
#'
#' @return An \code{\link{new_obc_regime_path}} object.
#' @export
as_obc_regime_path <- function(res, n_constr = NULL, paths_endo_by_T = NULL,
                               constraints = NULL) {
  if (!is.null(res$regime) && is.matrix(res$regime) && is.logical(res$regime)) {
    binding <- res$regime
  } else {
    codes <- res$regime_path %||% res$active_set
    if (is.null(codes))
      stop("as_obc_regime_path(): no 'regime'/'regime_path'/'active_set' field.")
    if (is.null(n_constr))
      stop("as_obc_regime_path(): n_constr is required to decode an integer ",
           "regime bitfield.")
    binding <- .obc_bitfield_to_binding(codes, n_constr)
  }
  paths <- res$paths %||% res$Y
  if (!is.null(paths)) {
    n_per <- nrow(binding)
    transpose <- if (!is.null(paths_endo_by_T)) isTRUE(paths_endo_by_T)
                 else (nrow(paths) != n_per && ncol(paths) == n_per)
    if (transpose) paths <- t(paths)
  }
  new_obc_regime_path(binding, paths = paths,
                      converged = res$converged %||% NA,
                      terminal_state = res$terminal_state,
                      constraints = constraints)
}


#' @export
print.obc_regime_path <- function(x, ...) {
  cat(sprintf("<obc_regime_path> %d periods, %d constraint(s); %d binding period(s)\n",
              x$n_periods, x$n_constraints, length(x$active_periods)))
  cat(sprintf("  converged: %s   hash: %s\n",
              as.character(x$converged),
              if (nchar(x$path_hash) > 40) paste0(substr(x$path_hash, 1, 40), "...")
              else x$path_hash))
  invisible(x)
}
