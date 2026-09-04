## R/hank-threeasset-tools.R
## --------------------------------------------------------------------------
## Production tooling around the three-asset block: provenance fingerprints,
## an arbitrary-path budget/Walras audit, checkpointed long-horizon Jacobian
## assembly, and cheap spot verification of individual fake-news columns.
##
## These exist because the three-asset block is used at a scale where the
## ordinary idioms stop working: 1.29M states, T_h = 400, seven inputs. At
## that size you cannot hash a block by serialising it, cannot re-run the full
## ND battery to check a suspicion, and cannot afford to lose a multi-hour
## Jacobian to a restart.
## --------------------------------------------------------------------------


#' Content fingerprint of a three-asset household block
#'
#' A short, stable hash of everything that makes a \code{\link{hank_het3_block}}
#' the block it is: the three grids, the income process (\code{Pi}, \code{e}),
#' preferences, the adjustment-cost parameters, the prices (including
#' \code{px} and any transfer), the transition inputs, the converged policies
#' and the stationary distribution.
#'
#' \strong{What this is for.} Keying a cache of expensive derived objects --
#' a long-horizon Jacobian, a GE solve -- on the block that produced them.
#' Hashing the serialised block itself is the obvious alternative and is a
#' trap at production scale: it is multi-gigabyte, and it also folds in fields
#' that have no bearing on the mathematics (timings, iteration counts, the
#' peak-RSS reading), so two identical households solved on different days
#' hash differently and the cache never hits.
#'
#' \strong{What it deliberately excludes}, for that reason: \code{elapsed_*},
#' \code{iterations}, \code{converged}, \code{backend}, \code{threads}. The
#' compiled and R backends are required to agree, and thread count is
#' bit-identical by contract, so a block is the same block either way. If you
#' need to distinguish HOW a block was produced, that is what
#' \code{\link{hank_het3_manifest}} is for -- fingerprint answers "is this the
#' same household?", the manifest answers "how did this run go?".
#'
#' The policies and distribution are included even though they are implied by
#' the calibration, because they are what downstream objects were actually
#' built from: a block re-solved at a looser \code{tol} is a genuinely
#' different object to differentiate around, and must not share a cache entry.
#'
#' @param block A \code{\link{hank_het3_block}}.
#' @param n Number of hex characters to return (1-16, default 16).
#' @return A length-1 character string of \code{n} lowercase hex digits.
#' @seealso \code{\link{hank_het3_manifest}} (run provenance, not identity)
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' hank_het3_fingerprint(blk)
#' @keywords internal
#' @export
hank_het3_fingerprint <- function(block, n = 16L) {
  if (!inherits(block, "hank_het3_block"))
    stop("hank_het3_fingerprint: 'block' must be a hank_het3_block.")
  ## Order is fixed and explicit, NOT names(block): a future field appended to
  ## the block would otherwise silently change every existing fingerprint.
  parts <- list(
    block$d_grid, block$f_grid, block$a_grid, block$Pi, block$e,
    block$beta, block$eis,
    block$rd, block$rf, block$ra, block$w,
    if (is.null(block$px)) 1 else block$px,
    .hank_block_tr(block), .hank_block_omega(block),
    block$chi0, block$chi1, block$chi2,
    block$phi0, block$phi1, block$phi2,
    block$d, block$f, block$a, block$c, block$D,
    ## Transition inputs by NAME as well as value: two blocks whose Pi_fn
    ## arguments differ only in name are different wirings of the DAG.
    names(block$Pi_inputs), unlist(block$Pi_inputs, use.names = FALSE)
  )
  .hank_fingerprint_hash(parts, n, "hank_het3_fingerprint")
}


#' Aggregate budget (Walras) audit along an arbitrary three-asset path
#'
#' Decomposes the economy-wide household budget date by date along any path
#' produced by \code{\link{hank_td3_nonlinear}}, returning each term SEPARATELY
#' rather than only the residual.
#'
#' \strong{Why the components and not just the residual.} A valuation error in
#' a three-asset economy is almost always a counterparty error: some quantity
#' is booked as a source on one side of the accounts and never as a use on the
#' other. The net residual can be small while two large components are each
#' wrong, and a residual-only check cannot localise anything. With the terms
#' separated, the side that fails to move is the side with the bug.
#'
#' The identity, at each date \eqn{t}, is
#' \deqn{\underbrace{w_t \bar e + Tr_t + (1+rd_t) D_{t-1} + p_{x,t}(1+rf_t) F_{t-1} + (1+ra_t) A_{t-1}}_{sources}
#'       = \underbrace{C_t + D_t + p_{x,t} F_t + A_t + \chi_t + \phi_t}_{uses},}
#' where the stocks carried IN are aggregates of the beginning-of-period
#' distribution over the STATE grids, and the stocks carried out are
#' aggregates of the policies. \eqn{\chi} and \eqn{\phi} are real goods
#' absorbed by rebalancing, which is exactly why the block reports them.
#'
#' @param block A \code{\link{hank_het3_block}}.
#' @param td The result of \code{\link{hank_td3_nonlinear}} on \code{block}
#'   (its \code{Dpath} is required, so call it without \code{D0} tricks that
#'   would desynchronise the two).
#' @param rd_path,rf_path,ra_path,w_path,px_path,Tr_path The paths \code{td}
#'   was computed at. \code{NULL} means the block's steady-state level, held
#'   constant -- the same convention as \code{\link{hank_td3_nonlinear}}, so
#'   passing the same arguments to both is correct by construction.
#' @return A \code{data.frame} with one row per date and columns
#'   \code{labour}, \code{transfer}, \code{d_payoff}, \code{f_payoff},
#'   \code{a_payoff}, \code{sources}, \code{consumption}, \code{d_purchase},
#'   \code{f_purchase}, \code{a_purchase}, \code{chi}, \code{phi},
#'   \code{uses}, \code{residual} (\code{sources - uses}), and
#'   \code{rel_residual} (scaled by \code{sources}).
#' @seealso \code{\link{hank_td3_nonlinear}}, \code{\link{hank_het3_block}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' td <- hank_td3_nonlinear(blk, T_h = 3)
#' hank_walras3(blk, td)$rel_residual
#' @export
hank_walras3 <- function(block, td,
                         rd_path = NULL, rf_path = NULL, ra_path = NULL,
                         w_path = NULL, px_path = NULL, Tr_path = NULL) {
  if (!inherits(block, "hank_het3_block"))
    stop("hank_walras3: 'block' must be a hank_het3_block.")
  if (is.null(td$Dpath))
    stop("hank_walras3: 'td' must be a hank_td3_nonlinear result (no Dpath ",
         "found).")
  T_h <- ncol(td$Dpath)
  fill <- function(p, v) if (is.null(p)) rep(v, T_h) else {
    if (length(p) != T_h)
      stop("hank_walras3: every supplied path must have length ", T_h,
           " (the horizon of 'td').")
    p
  }
  rd <- fill(rd_path, block$rd); rf <- fill(rf_path, block$rf)
  ra <- fill(ra_path, block$ra); w  <- fill(w_path,  block$w)
  px <- fill(px_path, if (is.null(block$px)) 1 else block$px)
  Tr <- fill(Tr_path, .hank_block_tr(block))

  ## State-grid aggregates over the BEGINNING-of-period distribution. These
  ## must be built from the grids, not from the policies: the policy is what
  ## the household chooses to carry OUT, and the sources side needs what it
  ## carried IN. Package cell order is e slowest, then d, f, with a fastest.
  ne <- block$n_e; nd <- block$n_d; nf <- block$n_f; na <- block$n_a
  e_st <- rep(block$e, each = nd * nf * na)
  d_st <- rep(rep(block$d_grid, each = nf * na), ne)
  f_st <- rep(rep(block$f_grid, each = na), ne * nd)
  a_st <- rep(block$a_grid, ne * nd * nf)
  omega <- .hank_block_omega(block)
  om_st <- rep(omega, each = nd * nf * na)

  agg <- function(x) as.numeric(crossprod(td$Dpath, x))
  e_bar <- agg(e_st); om_bar <- agg(om_st)
  D_in  <- agg(d_st); F_in <- agg(f_st); A_in <- agg(a_st)

  labour   <- w * e_bar
  transfer <- Tr * om_bar
  d_pay    <- (1 + rd) * D_in
  f_pay    <- px * (1 + rf) * F_in
  a_pay    <- (1 + ra) * A_in
  sources  <- labour + transfer + d_pay + f_pay + a_pay
  uses     <- td$C + td$D + px * td$F + td$A + td$CHI + td$PHI
  resid    <- sources - uses
  data.frame(labour = labour, transfer = transfer,
             d_payoff = d_pay, f_payoff = f_pay, a_payoff = a_pay,
             sources = sources,
             consumption = td$C, d_purchase = td$D,
             f_purchase = px * td$F, a_purchase = td$A,
             chi = td$CHI, phi = td$PHI, uses = uses,
             residual = resid,
             rel_residual = resid / pmax(abs(sources), .Machine$double.eps))
}


#' Checkpointed three-asset fake-news Jacobian
#'
#' \code{\link{hank_het3_jacobian}} one input at a time, writing each input's
#' columns to \code{dir} as they complete and skipping any already present.
#' An interrupted run resumes where it stopped instead of restarting.
#'
#' At the production scale this exists for -- 1.29M states, \code{T_h = 400},
#' seven inputs -- a single Jacobian is hours of work, and losing all of it to
#' a restart is the actual risk. Because the inputs are INDEPENDENT sweeps
#' (each does its own backward pass; only the steady-state expectation vectors
#' are shared, and those are recomputed cheaply), per-input is the natural
#' checkpoint granularity: no partial state has to be serialised, and every
#' file on disk is a complete, usable column block.
#'
#' \strong{Staleness is checked, not assumed.} Each file records the block's
#' \code{\link{hank_het3_fingerprint}} and the horizon; a cached file that
#' does not match the current block and \code{T_h} is REFUSED rather than
#' silently reused, because a Jacobian from a neighbouring calibration is the
#' most expensive possible thing to accept by mistake.
#'
#' @param block A \code{\link{hank_het3_block}}.
#' @param T_h Integer horizon.
#' @param dir Directory for the checkpoint files; created if absent.
#' @param inputs,outputs As in \code{\link{hank_het3_jacobian}}.
#' @param verbose Report each input as it is computed or reused.
#' @param ... Further arguments passed to \code{\link{hank_het3_jacobian}}
#'   (\code{delta_in}, \code{delta_v}, \code{delta_d}, \code{threads}).
#' @return The same nested list \code{J[[output]][[input]]} that
#'   \code{\link{hank_het3_jacobian}} returns, assembled from the checkpoints.
#' @seealso \code{\link{hank_het3_jacobian}},
#'   \code{\link{hank_het3_jacobian_spot}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' dir <- tempfile("jac3"); 
#' J <- hank_het3_jacobian_checkpoint(blk, T_h = 2, dir = dir,
#'                                    inputs = "rd", outputs = "C")
#' dim(J$C$rd)
#' @keywords internal
#' @export
hank_het3_jacobian_checkpoint <- function(block, T_h, dir, inputs = NULL,
                                          outputs = .hank3_jac_outputs,
                                          verbose = TRUE, ...) {
  if (!inherits(block, "hank_het3_block"))
    stop("hank_het3_jacobian_checkpoint: block must be hank_het3_block")
  if (is.null(inputs)) inputs <- .hank3_jac_inputs(block)
  inputs <- .hank3_check_inputs(block, inputs,
                                "hank_het3_jacobian_checkpoint")
  bad <- setdiff(outputs, .hank3_jac_outputs)
  if (length(bad))
    stop("hank_het3_jacobian_checkpoint: unsupported output(s) ",
         paste(bad, collapse = ", "))
  if (!dir.exists(dir))
    dir.create(dir, recursive = TRUE)

  fp <- hank_het3_fingerprint(block)
  J <- setNames(vector("list", length(outputs)), outputs)
  for (o in outputs) J[[o]] <- setNames(vector("list", length(inputs)), inputs)

  for (i in inputs) {
    path <- file.path(dir, paste0("jac3_", i, ".rds"))
    piece <- NULL
    if (file.exists(path)) {
      cached <- readRDS(path)
      ok <- identical(cached$fingerprint, fp) &&
        identical(as.integer(cached$T_h), as.integer(T_h)) &&
        all(outputs %in% names(cached$columns))
      if (ok) {
        piece <- cached$columns
        if (verbose) message("hank_het3_jacobian_checkpoint: reusing '", i, "'")
      } else if (!identical(cached$fingerprint, fp) ||
                 !identical(as.integer(cached$T_h), as.integer(T_h))) {
        ## Refuse rather than overwrite: a mismatch means this directory is
        ## being shared between calibrations, and the OTHER files in it are
        ## probably stale too. Silently recomputing one would leave the run
        ## mixing horizons or households.
        stop("hank_het3_jacobian_checkpoint: '", path, "' was written for a ",
             "different block or horizon (cached fingerprint ",
             cached$fingerprint, " at T_h = ", cached$T_h, "; this block is ",
             fp, " at T_h = ", T_h, "). Point 'dir' at a fresh directory, or ",
             "delete the stale checkpoints deliberately.")
      }
    }
    if (is.null(piece)) {
      if (verbose) message("hank_het3_jacobian_checkpoint: computing '", i, "'")
      one <- hank_het3_jacobian(block, T_h, inputs = i, outputs = outputs, ...)
      piece <- setNames(lapply(outputs, function(o) one[[o]][[i]]), outputs)
      ## Write to a temporary name and rename: a kill between "file exists"
      ## and "file is complete" would otherwise leave a truncated checkpoint
      ## that the resume path would happily read.
      tmp <- paste0(path, ".part")
      saveRDS(list(fingerprint = fp, T_h = as.integer(T_h), input = i,
                   columns = piece), tmp)
      file.rename(tmp, path)
    }
    for (o in outputs) J[[o]][[i]] <- piece[[o]]
  }
  J
}


#' Spot-verify individual fake-news Jacobian columns against numerical differentiation
#'
#' Runs the numerical-differentiation oracle at a FEW shock dates instead of
#' all of them, and compares against \code{\link{hank_het3_jacobian}}.
#'
#' \code{\link{hank_het3_jacobian_nd}} is the acceptance gate, and at
#' \code{T_h = 400} it is unaffordable: its cost is \eqn{O(T_h)} full
#' nonlinear transitions per (input, date), i.e. quadratic in the horizon.
#' But the fake-news recursion builds every column of a given input from ONE
#' backward sweep, so an error in that sweep shows up at every date -- which
#' makes a handful of dates strong evidence at a fraction of the cost. Use
#' this to check a production block; use the full oracle to accept a code
#' change.
#'
#' The dates default to \code{s = 1} and the largest requested date, which is
#' the informative pair: \code{s = 1} is the contemporaneous term (the only
#' one where a transition input's direct-\code{Pi} channel enters), and a far
#' date exercises the anticipation recursion after many steps of propagation.
#'
#' @param block A \code{\link{hank_het3_block}}.
#' @param T_h Integer horizon.
#' @param inputs,outputs As in \code{\link{hank_het3_jacobian}}.
#' @param dates Integer shock dates (columns) to verify; default
#'   \code{c(1, T_h)}.
#' @param delta FD step for the numerical oracle.
#' @param J Optional precomputed \code{\link{hank_het3_jacobian}} result to
#'   check (e.g. one read back from
#'   \code{\link{hank_het3_jacobian_checkpoint}}); computed here if absent.
#' @param ... Further arguments passed to \code{\link{hank_het3_jacobian}}
#'   when \code{J} is not supplied.
#' @return A \code{data.frame} with one row per (input, output, date) and
#'   columns \code{input}, \code{output}, \code{date}, \code{max_abs_dev},
#'   \code{scale} (the largest absolute ND entry in that column) and
#'   \code{rel_dev}. A column that is identically zero in BOTH routes reports
#'   \code{rel_dev = 0}.
#' @seealso \code{\link{hank_het3_jacobian_nd}} (the full oracle),
#'   \code{\link{hank_het3_jacobian_checkpoint}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' hank_het3_jacobian_spot(blk, T_h = 3, inputs = "rd", outputs = "C")
#' @keywords internal
#' @export
hank_het3_jacobian_spot <- function(block, T_h, inputs = NULL,
                                    outputs = .hank3_jac_outputs,
                                    dates = NULL, delta = 1e-5, J = NULL,
                                    ...) {
  if (!inherits(block, "hank_het3_block"))
    stop("hank_het3_jacobian_spot: block must be hank_het3_block")
  if (is.null(inputs)) inputs <- .hank3_jac_inputs(block)
  inputs <- .hank3_check_inputs(block, inputs, "hank_het3_jacobian_spot")
  bad <- setdiff(outputs, .hank3_jac_outputs)
  if (length(bad))
    stop("hank_het3_jacobian_spot: unsupported output(s) ",
         paste(bad, collapse = ", "))
  if (is.null(dates)) dates <- unique(c(1L, as.integer(T_h)))
  dates <- as.integer(dates)
  if (any(dates < 1L) || any(dates > T_h))
    stop("hank_het3_jacobian_spot: 'dates' must lie in 1:T_h.")
  if (is.null(J)) J <- hank_het3_jacobian(block, T_h, inputs = inputs,
                                          outputs = outputs, ...)

  rows <- list()
  for (i in inputs) {
    ## One ND column per requested date: this is the whole saving, and it is
    ## why the oracle is re-derived here rather than calling
    ## hank_het3_jacobian_nd, which has no way to sweep a subset of dates.
    for (s in dates) {
      nd_col <- .hank3_nd_column(block, T_h, i, outputs, s, delta)
      for (o in outputs) {
        dev <- max(abs(J[[o]][[i]][, s] - nd_col[[o]]))
        sc  <- max(abs(nd_col[[o]]))
        rows[[length(rows) + 1L]] <- data.frame(
          input = i, output = o, date = s,
          max_abs_dev = dev, scale = sc,
          rel_dev = if (sc > 0) dev / sc else 0,
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}


## One numerically-differentiated Jacobian COLUMN: the response of every
## output at every date to a perturbation of `input` at date `s`. This is the
## inner loop of hank_het3_jacobian_nd, factored out so the spot checker can
## sweep a subset of dates without recomputing the rest.
#' @keywords internal
.hank3_nd_column <- function(block, T_h, input, outputs, s, delta) {
  base <- list(rd = rep(block$rd, T_h), rf = rep(block$rf, T_h),
               ra = rep(block$ra, T_h), w = rep(block$w, T_h),
               px = rep(if (is.null(block$px)) 1 else block$px, T_h),
               Tr = rep(.hank_block_tr(block), T_h))
  as_paths <- function(p) setNames(p, paste0(names(p), "_path"))
  p <- m <- base; pip <- pim <- NULL
  if (.hank3_is_pi_input(block, input)) {
    x0 <- rep(block$Pi_inputs[[input]], T_h)
    xp <- x0; xp[s] <- xp[s] + delta
    xm <- x0; xm[s] <- xm[s] - delta
    pip <- setNames(list(xp), input); pim <- setNames(list(xm), input)
  } else if (input == "px" && .hank3_px_is_inert(block)) {
    return(setNames(lapply(outputs, function(o) rep(0, T_h)), outputs))
  } else {
    p[[input]][s] <- p[[input]][s] + delta
    m[[input]][s] <- m[[input]][s] - delta
  }
  run <- function(pp, pi_paths) do.call(
    hank_td3_nonlinear,
    c(list(block = block, T_h = T_h, pi_input_paths = pi_paths), as_paths(pp)))
  op <- run(p, pip); om <- run(m, pim)
  setNames(lapply(outputs, function(o) (op[[o]] - om[[o]]) / (2 * delta)),
           outputs)
}
