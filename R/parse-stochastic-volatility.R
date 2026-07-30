## R/parse-stochastic-volatility.R
## --------------------------------------------------------------------------
## Parser for the `stochastic_volatility ... end;` block:
##
##   stochastic_volatility;
##     var e_a; mu = mu_a; rho = rho_a; sigma_eta = sig_eta_a;
##     var e_r; mu = 0;    rho = 0.9;   sigma_eta = 0.2;
##   end;
##
## Structure mirrors parse-heteroskedastic-shocks.R. Each `var NAME; mu = ...;
## rho = ...; sigma_eta = ...;` group becomes one row of a data.frame with
## list-columns mu/rho/sigma_eta. Each hyperparameter value is either a numeric
## literal (fixed) or an identifier naming a model parameter (usually estimated).
## Output is identical in structure to the stochastic_volatility() constructor.
## --------------------------------------------------------------------------

#' Parse the stochastic_volatility block into a structured list
#'
#' @param body Body text of the stochastic_volatility block (between
#'   \code{stochastic_volatility;} and \code{end;}).
#' @return A list with one element \code{sv}: a data.frame with column
#'   \code{shock} (character) and list columns \code{mu}, \code{rho},
#'   \code{sigma_eta} (each element a length-1 numeric or character).
#' @noRd
parse_stochastic_volatility_block <- function(body) {
  entries <- list()

  stmts <- strsplit(body, ";")[[1]]
  stmts <- trimws(stmts)
  stmts <- stmts[nchar(stmts) > 0]

  cur_shock <- NULL
  cur_vals  <- list(mu = NULL, rho = NULL, sigma_eta = NULL)

  ## A hyperparameter RHS is a numeric literal or a parameter-name identifier.
  parse_val <- function(rhs, key, shock) {
    rhs <- trimws(rhs)
    num <- suppressWarnings(as.numeric(rhs))
    if (!is.na(num)) return(num)
    if (grepl("^[A-Za-z_][A-Za-z0-9_]*$", rhs)) return(rhs)
    stop(sprintf(
      paste0("stochastic_volatility: %s for shock '%s' must be a number or a ",
             "parameter name (got '%s')."), key, shock, rhs), call. = FALSE)
  }

  flush <- function() {
    if (is.null(cur_shock)) return(invisible())
    for (k in c("mu", "rho", "sigma_eta")) {
      if (is.null(cur_vals[[k]]))
        stop(sprintf("stochastic_volatility: var %s has no '%s' statement.",
                     cur_shock, k), call. = FALSE)
    }
    entries[[length(entries) + 1L]] <<- list(
      shock     = cur_shock,
      mu        = cur_vals$mu,
      rho       = cur_vals$rho,
      sigma_eta = cur_vals$sigma_eta
    )
    cur_shock <<- NULL
    cur_vals  <<- list(mu = NULL, rho = NULL, sigma_eta = NULL)
  }

  for (s in stmts) {
    if (grepl("^\\s*end\\s*$", s, ignore.case = TRUE)) next

    m_var <- regmatches(s, regexec(
      "^\\s*var\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*$", s, perl = TRUE))[[1]]
    if (length(m_var) > 0 && nchar(m_var[1]) > 0) {
      flush()
      cur_shock <- m_var[2]
      next
    }

    ## key = value  (key in {mu, rho, sigma_eta})
    m_kv <- regmatches(s, regexec(
      "^\\s*(mu|rho|sigma_eta)\\s*=\\s*(.+)$", s, perl = TRUE))[[1]]
    if (length(m_kv) > 0 && nchar(m_kv[1]) > 0) {
      if (is.null(cur_shock))
        stop(sprintf(
          "stochastic_volatility: '%s' statement without a preceding 'var'.",
          m_kv[2]), call. = FALSE)
      cur_vals[[m_kv[2]]] <- parse_val(m_kv[3], m_kv[2], cur_shock)
      next
    }

    stop(sprintf("stochastic_volatility: unrecognised statement '%s'.", s),
         call. = FALSE)
  }
  flush()

  if (length(entries) == 0L) {
    empty <- data.frame(shock = character(0), stringsAsFactors = FALSE)
    empty$mu <- list(); empty$rho <- list(); empty$sigma_eta <- list()
    return(list(sv = empty))
  }

  shocks <- vapply(entries, function(x) x$shock, character(1))
  if (anyDuplicated(shocks))
    stop(sprintf("stochastic_volatility: duplicate shock(s): %s.",
                 paste(unique(shocks[duplicated(shocks)]), collapse = ", ")),
         call. = FALSE)

  df <- data.frame(shock = shocks, stringsAsFactors = FALSE)
  df$mu        <- lapply(entries, function(x) x$mu)
  df$rho       <- lapply(entries, function(x) x$rho)
  df$sigma_eta <- lapply(entries, function(x) x$sigma_eta)
  list(sv = df)
}
