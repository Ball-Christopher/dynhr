# ===========================================================================
# dynhr_aux_expansion.R
# Auxiliary variable expansion for Dynare .mod files with leads/lags > 1
#
# Dynare's first-order perturbation requires max lead = 1 and max lag = 1.
# Variables with higher leads/lags (e.g., dp(4) in a Taylor rule) must be
# reduced via chain variables:
#
#   dp(4) in equation  =>  introduce:
#     AUX_LEAD_dp_1 = dp(+1)
#     AUX_LEAD_dp_2 = AUX_LEAD_dp_1(+1)
#     AUX_LEAD_dp_3 = AUX_LEAD_dp_2(+1)
#   and rewrite dp(4) => AUX_LEAD_dp_3(+1)
#
# Same logic for lags > 1:
#   x(-3) in equation  =>  introduce:
#     AUX_LAG_x_1 = x(-1)
#     AUX_LAG_x_2 = AUX_LAG_x_1(-1)
#   and rewrite x(-3) => AUX_LAG_x_2(-1)
#
# This function operates on the raw model block text, BEFORE parse_model_block
# is called, so all downstream code (variable classification, LLI construction,
# Jacobian compilation) works automatically.
# ===========================================================================

.expand_aux_timing <- function(model_body, var_names, exo_names = character(0),
                               verbose = TRUE) {

  # exo_names: names that are EXOGENOUS. These must be treated specially: an
  # exogenous variable at ANY nonzero lead/lag (incl. +-1) becomes an
  # AUX_EXO_LEAD_* / AUX_EXO_LAG_* ENDOGENOUS aux chain (Dynare convention),
  # because the perturbation state space cannot carry an exogenous variable at
  # a nonzero timing. Endogenous +-1 leads/lags are handled natively and get
  # no aux var (unchanged behaviour). See H8.
  exo_names <- setdiff(exo_names, var_names)  # an endo name is never exo

  # Collect all variable references with their timings
  # Matches patterns like: varname(+4), varname(-3), varname(4), varname( +4 )
  # Must be careful not to match parameter names or function calls.
  # Strategy: find all name(number) patterns, filter to known var_names.
  
  # Build regex that matches any endogenous variable name followed by (timing)
  # We escape names and join with | for alternation.
  # Sort by length descending so longer names match first (e.g., pxstar before pstar)
  all_names    <- c(var_names, exo_names)
  sorted_names <- all_names[order(nchar(all_names), decreasing = TRUE)]
  escaped_names <- gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", sorted_names)
  name_pattern <- paste(escaped_names, collapse = "|")
  
  # Pattern: variable name, then '(' with optional spaces, optional +/-, digits, ')'
  # Capture groups: (1) variable name, (2) full timing string including sign
  full_pattern <- paste0(
    "\\b(", name_pattern, ")\\s*\\(\\s*([+-]?\\s*\\d+)\\s*\\)"
  )
  
  # Find all matches
  matches <- gregexpr(full_pattern, model_body, perl = TRUE)
  match_strings <- regmatches(model_body, matches)[[1]]
  
  if (length(match_strings) == 0) {
    if (verbose) cat("  No leads/lags > 1 found; no auxiliary expansion needed.\n")
    return(list(
      model_body   = model_body,
      aux_var_names = character(0),
      aux_info     = data.frame(
        aux_name         = character(0),
        original_var     = character(0),
        timing_type      = character(0),
        chain_index      = integer(0),
        represents_timing = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # Parse each match to extract variable name and timing
  timing_records <- data.frame(
    var_name = character(0),
    timing   = integer(0),
    stringsAsFactors = FALSE
  )
  
  for (ms in match_strings) {
    m <- regmatches(ms, regexec(full_pattern, ms, perl = TRUE))[[1]]
    if (length(m) >= 3) {
      vname <- m[2]
      # Remove internal spaces from timing string (e.g., "+ 4" -> "+4")
      tstr <- gsub("\\s+", "", m[3])
      tval <- as.integer(tstr)
      if (!is.na(tval) && vname %in% all_names) {
        timing_records <- rbind(timing_records, data.frame(
          var_name = vname,
          timing   = tval,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (nrow(timing_records) == 0) {
    if (verbose) cat("  No endogenous leads/lags > 1 found.\n")
    return(list(
      model_body    = model_body,
      aux_var_names = character(0),
      aux_info      = data.frame(
        aux_name         = character(0),
        original_var     = character(0),
        timing_type      = character(0),
        chain_index      = integer(0),
        represents_timing = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # Split endogenous vs exogenous timing records.
  endo_records <- timing_records[timing_records$var_name %in% var_names, , drop = FALSE]
  exo_records  <- timing_records[timing_records$var_name %in% exo_names,  , drop = FALSE]

  # Endogenous: aux only needed for |timing| > 1 (the native LLI handles +-1).
  leads_gt1 <- endo_records[endo_records$timing > 1,  , drop = FALSE]
  lags_gt1  <- endo_records[endo_records$timing < -1, , drop = FALSE]

  # Exogenous: aux needed for ANY nonzero timing (incl. +-1).
  exo_leads <- exo_records[exo_records$timing >= 1, , drop = FALSE]
  exo_lags  <- exo_records[exo_records$timing <= -1, , drop = FALSE]

  if (nrow(leads_gt1) == 0 && nrow(lags_gt1) == 0 &&
      nrow(exo_leads) == 0 && nrow(exo_lags) == 0) {
    if (verbose) cat("  All leads/lags are within [-1, +1]; no auxiliary expansion needed.\n")
    return(list(
      model_body    = model_body,
      aux_var_names = character(0),
      aux_info      = data.frame(
        aux_name         = character(0),
        original_var     = character(0),
        timing_type      = character(0),
        chain_index      = integer(0),
        represents_timing = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # For each variable, find the maximum lead and maximum lag depth
  # We need chains up to max(lead)-1 for leads and max(|lag|)-1 for lags
  aux_var_names <- character(0)
  aux_equations <- character(0)
  aux_info_rows <- list()
  
  new_body <- model_body
  
  # ---- Process LEADS > 1 ----
  if (nrow(leads_gt1) > 0) {
    lead_vars <- unique(leads_gt1$var_name)
    
    for (vn in lead_vars) {
      max_lead <- max(leads_gt1$timing[leads_gt1$var_name == vn])
      
      if (verbose) {
        cat(sprintf("  Expanding %s with max lead %d: creating %d auxiliary variable(s)\n",
                    vn, max_lead, max_lead - 1L))
      }
      
      # Build the chain: AUX_LEAD_vn_1, AUX_LEAD_vn_2, ..., AUX_LEAD_vn_{max_lead-1}
      chain_names <- character(max_lead - 1L)
      for (k in seq_len(max_lead - 1L)) {
        aux_name <- paste0("AUX_LEAD_", vn, "_", k)
        chain_names[k] <- aux_name
        aux_var_names <- c(aux_var_names, aux_name)
        
        # Build the chain equation:
        #   k=1: AUX_LEAD_vn_1 = vn(1)          i.e., vn(+1)
        #   k=2: AUX_LEAD_vn_2 = AUX_LEAD_vn_1(1)
        #   k=j: AUX_LEAD_vn_j = AUX_LEAD_vn_{j-1}(1)
        if (k == 1L) {
          rhs <- paste0(vn, "(1)")
        } else {
          rhs <- paste0(chain_names[k - 1L], "(1)")
        }
        eq_str <- paste0(aux_name, " = ", rhs, ";")
        aux_equations <- c(aux_equations, eq_str)
        
        # Record info
        aux_info_rows[[length(aux_info_rows) + 1L]] <- data.frame(
          aux_name          = aux_name,
          original_var      = vn,
          timing_type       = "lead",
          chain_index       = k,
          represents_timing = k + 1L,
          stringsAsFactors  = FALSE
        )
      }
      
      # Rewrite references in the model body.
      # Process from highest timing downward to avoid partial replacement issues.
      # vn(N) for N > 1 => AUX_LEAD_vn_{N-1}(1)
      timings_to_fix <- sort(unique(leads_gt1$timing[leads_gt1$var_name == vn]),
                             decreasing = TRUE)
      
      for (tt in timings_to_fix) {
        # Build pattern to match vn(tt) or vn(+tt) with optional spaces
        # e.g., dp(4), dp(+4), dp( 4 ), dp( +4 )
        pat <- paste0(
          "\\b", gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", vn),
          "\\s*\\(\\s*\\+?\\s*", tt, "\\s*\\)"
        )
        replacement <- paste0(chain_names[tt - 1L], "(1)")
        
        new_body <- gsub(pat, replacement, new_body, perl = TRUE)
        
        if (verbose) {
          cat(sprintf("    Rewrote %s(%d) => %s(1)\n", vn, tt, chain_names[tt - 1L]))
        }
      }
    }
  }
  
  # ---- Process LAGS < -1 ----
  if (nrow(lags_gt1) > 0) {
    lag_vars <- unique(lags_gt1$var_name)
    
    for (vn in lag_vars) {
      max_lag_depth <- max(abs(lags_gt1$timing[lags_gt1$var_name == vn]))
      
      if (verbose) {
        cat(sprintf("  Expanding %s with max lag %d: creating %d auxiliary variable(s)\n",
                    vn, max_lag_depth, max_lag_depth - 1L))
      }
      
      # Build the chain: AUX_LAG_vn_1, AUX_LAG_vn_2, ..., AUX_LAG_vn_{max_lag_depth-1}
      chain_names <- character(max_lag_depth - 1L)
      for (k in seq_len(max_lag_depth - 1L)) {
        aux_name <- paste0("AUX_LAG_", vn, "_", k)
        chain_names[k] <- aux_name
        aux_var_names <- c(aux_var_names, aux_name)
        
        # Chain equation:
        #   k=1: AUX_LAG_vn_1 = vn(-1)
        #   k=j: AUX_LAG_vn_j = AUX_LAG_vn_{j-1}(-1)
        if (k == 1L) {
          rhs <- paste0(vn, "(-1)")
        } else {
          rhs <- paste0(chain_names[k - 1L], "(-1)")
        }
        eq_str <- paste0(aux_name, " = ", rhs, ";")
        aux_equations <- c(aux_equations, eq_str)
        
        aux_info_rows[[length(aux_info_rows) + 1L]] <- data.frame(
          aux_name          = aux_name,
          original_var      = vn,
          timing_type       = "lag",
          chain_index       = k,
          represents_timing = -(k + 1L),
          stringsAsFactors  = FALSE
        )
      }
      
      # Rewrite references from deepest lag to shallowest.
      # vn(-N) for N > 1 => AUX_LAG_vn_{N-1}(-1)
      timings_to_fix <- sort(unique(abs(lags_gt1$timing[lags_gt1$var_name == vn])),
                             decreasing = TRUE)
      
      for (tt in timings_to_fix) {
        pat <- paste0(
          "\\b", gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", vn),
          "\\s*\\(\\s*-\\s*", tt, "\\s*\\)"
        )
        replacement <- paste0(chain_names[tt - 1L], "(-1)")
        
        new_body <- gsub(pat, replacement, new_body, perl = TRUE)
        
        if (verbose) {
          cat(sprintf("    Rewrote %s(-%d) => %s(-1)\n", vn, tt, chain_names[tt - 1L]))
        }
      }
    }
  }
  
  # ---- Process EXOGENOUS LEADS >= +1 (H8) ----
  # An exo var x at lead +N gets N aux ENDOGENOUS vars:
  #   AUX_EXO_LEAD_x_1 = x(+1)
  #   AUX_EXO_LEAD_x_k = AUX_EXO_LEAD_x_{k-1}(+1)   (k = 2..N)
  # and every reference x(+t) is rewritten to AUX_EXO_LEAD_x_t (at time t).
  if (nrow(exo_leads) > 0) {
    lead_vars <- unique(exo_leads$var_name)

    for (vn in lead_vars) {
      max_lead <- max(exo_leads$timing[exo_leads$var_name == vn])

      if (verbose) {
        cat(sprintf("  Expanding exo %s with max lead %d: creating %d AUX_EXO_LEAD variable(s)\n",
                    vn, max_lead, max_lead))
      }

      chain_names <- character(max_lead)
      for (k in seq_len(max_lead)) {
        aux_name <- paste0("AUX_EXO_LEAD_", vn, "_", k)
        chain_names[k] <- aux_name
        aux_var_names <- c(aux_var_names, aux_name)

        if (k == 1L) {
          rhs <- paste0(vn, "(1)")
        } else {
          rhs <- paste0(chain_names[k - 1L], "(1)")
        }
        aux_equations <- c(aux_equations, paste0(aux_name, " = ", rhs, ";"))

        aux_info_rows[[length(aux_info_rows) + 1L]] <- data.frame(
          aux_name          = aux_name,
          original_var      = vn,
          timing_type       = "exo_lead",
          chain_index       = k,
          represents_timing = k,
          stringsAsFactors  = FALSE
        )
      }

      # Rewrite x(+t) => AUX_EXO_LEAD_x_t  (no explicit timing: aux is at t).
      # Highest timing first to avoid clobbering shorter timings.
      timings_to_fix <- sort(unique(exo_leads$timing[exo_leads$var_name == vn]),
                             decreasing = TRUE)
      for (tt in timings_to_fix) {
        pat <- paste0(
          "\\b", gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", vn),
          "\\s*\\(\\s*\\+?\\s*", tt, "\\s*\\)"
        )
        new_body <- gsub(pat, chain_names[tt], new_body, perl = TRUE)
        if (verbose) {
          cat(sprintf("    Rewrote %s(+%d) => %s\n", vn, tt, chain_names[tt]))
        }
      }
    }
  }

  # ---- Process EXOGENOUS LAGS <= -1 (H8 / NEW-1) ----
  # An exo var x at lag -N gets N aux ENDOGENOUS vars with the ENTRY node at t:
  #   AUX_EXO_LAG_x_1 = x          (entry node -> gives ghu[AUX_1, x] = 1)
  #   AUX_EXO_LAG_x_k = AUX_EXO_LAG_x_{k-1}(-1)   (k = 2..N)
  # and x(-N) is rewritten to AUX_EXO_LAG_x_N(-1).  (Endo lag-N uses N-1 because
  # the endo var already carries its own t-1 column; an exo has no base column,
  # hence the extra entry node.)
  if (nrow(exo_lags) > 0) {
    lag_vars <- unique(exo_lags$var_name)

    for (vn in lag_vars) {
      max_lag_depth <- max(abs(exo_lags$timing[exo_lags$var_name == vn]))

      if (verbose) {
        cat(sprintf("  Expanding exo %s with max lag %d: creating %d AUX_EXO_LAG variable(s)\n",
                    vn, max_lag_depth, max_lag_depth))
      }

      chain_names <- character(max_lag_depth)
      for (k in seq_len(max_lag_depth)) {
        aux_name <- paste0("AUX_EXO_LAG_", vn, "_", k)
        chain_names[k] <- aux_name
        aux_var_names <- c(aux_var_names, aux_name)

        if (k == 1L) {
          rhs <- vn                       # entry node at time t
        } else {
          rhs <- paste0(chain_names[k - 1L], "(-1)")
        }
        aux_equations <- c(aux_equations, paste0(aux_name, " = ", rhs, ";"))

        aux_info_rows[[length(aux_info_rows) + 1L]] <- data.frame(
          aux_name          = aux_name,
          original_var      = vn,
          timing_type       = "exo_lag",
          chain_index       = k,
          represents_timing = -(k - 1L),
          stringsAsFactors  = FALSE
        )
      }

      # Rewrite x(-t) => AUX_EXO_LAG_x_t(-1).  Deepest lag first.
      timings_to_fix <- sort(unique(abs(exo_lags$timing[exo_lags$var_name == vn])),
                             decreasing = TRUE)
      for (tt in timings_to_fix) {
        pat <- paste0(
          "\\b", gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", vn),
          "\\s*\\(\\s*-\\s*", tt, "\\s*\\)"
        )
        replacement <- paste0(chain_names[tt], "(-1)")
        new_body <- gsub(pat, replacement, new_body, perl = TRUE)
        if (verbose) {
          cat(sprintf("    Rewrote %s(-%d) => %s(-1)\n", vn, tt, chain_names[tt]))
        }
      }
    }
  }

  # ---- Append auxiliary equations to model body ----
  if (length(aux_equations) > 0) {
    aux_block <- paste0(
      "\n",
      paste(aux_equations, collapse = "\n"),
      "\n"
    )
    new_body <- paste0(new_body, aux_block)
  }
  
  # ---- Build aux_info data.frame ----
  if (length(aux_info_rows) > 0) {
    aux_info <- do.call(rbind, aux_info_rows)
    rownames(aux_info) <- NULL
  } else {
    aux_info <- data.frame(
      aux_name          = character(0),
      original_var      = character(0),
      timing_type       = character(0),
      chain_index       = integer(0),
      represents_timing = integer(0),
      stringsAsFactors  = FALSE
    )
  }
  
  if (verbose && length(aux_var_names) > 0) {
    cat(sprintf("  Total auxiliary variables added: %d\n", length(aux_var_names)))
    cat(sprintf("  Total auxiliary equations added: %d\n", length(aux_equations)))
    cat("  Auxiliary variables:", paste(aux_var_names, collapse = ", "), "\n")
  }
  
  list(
    model_body    = new_body,
    aux_var_names = aux_var_names,
    aux_info      = aux_info
  )
}