# ==============================================================================
# dynhr_log_parser.R
# Parse Dynare/Julia log files to extract expected results for cross-checking
# ==============================================================================

# ---- Internal helpers --------------------------------------------------------

.trim <- function(x) gsub("^\\s+|\\s+$", "", x)

# Core table parser for Dynare.jl Unicode box-drawing tables.
#
# Layout:
#   section header text (e.g. "  Steady state")
#   ----------???-----           <-- top border (contains ???)
#   header | col1  col2 ...    <-- (optional) column header row(s)
#   ----------------------     <-- (optional) all-dash separator line
#   row1   | val1  val2 ...    <-- data rows
#   ----------???-----           <-- bottom border (contains ???)
#
# The ??? and ??? delimit the content area.  Inside that area:
#   - Lines containing | are content lines (header or data).
#   - An all-dash line (only - / spaces / ???) is the separator between
#     the HEADER row(s) and the DATA rows.
#   - If there is NO separator, all |-lines are data (e.g. the steady-
#     state table has no column header because it is just name | value).
#
# Returns: list(col_names, row_names, values_matrix, end_idx)

.parse_unicode_table <- function(lines, start_idx) {

    n <- length(lines)

    # 1. Find ??? (top) and ??? (bottom) ----------------------------------------
    top_idx <- NA
    bot_idx <- NA
    for (i in start_idx:min(n, start_idx + 500)) {
        if (grepl("\u252c", lines[i])) {
            if (is.na(top_idx)) top_idx <- i
        }
        if (grepl("\u2534", lines[i])) {
            bot_idx <- i
            break
        }
    }

    if (is.na(top_idx) || is.na(bot_idx) || bot_idx <= top_idx + 1) {
        return(list(col_names = NULL, row_names = NULL,
                    values = NULL, end_idx = ifelse(is.na(bot_idx), start_idx, bot_idx)))
    }

    # 2. Collect content lines between ??? and ??? --------------------------------
    content_lines <- lines[(top_idx + 1):(bot_idx - 1)]

    # Classify each content line
    is_pipe <- grepl("\u2502", content_lines)             # has |
    is_sep  <- grepl("^\u2500", gsub("\\s", "", content_lines)) &
               !is_pipe                                    # all-dash, no |

    # 3. Find the separator that divides header from data ----------------------
    sep_positions <- which(is_sep)
    if (length(sep_positions) > 0) {
        first_sep <- sep_positions[1]
        header_lines <- content_lines[seq_len(first_sep - 1)]
        data_lines   <- content_lines[(first_sep + 1):length(content_lines)]
        # Keep only pipe lines
        header_lines <- header_lines[grepl("\u2502", header_lines)]
        data_lines   <- data_lines[grepl("\u2502", data_lines)]
    } else {
        # No separator -> no header, all are data
        header_lines <- character(0)
        data_lines   <- content_lines[is_pipe]
    }

    # 4. Parse header -> column names ------------------------------------------
    col_names <- NULL
    if (length(header_lines) > 0) {
        # Use the LAST header line (in case there are spanning headers)
        hline <- header_lines[length(header_lines)]
        parts <- strsplit(hline, "\u2502")[[1]]
        # First part is the row-label column header (often blank); rest are
        # value column names (may be packed into one string, whitespace-separated)
        if (length(parts) >= 2) {
            val_part <- paste(parts[-1], collapse = " ")
            col_names <- strsplit(.trim(val_part), "\\s{2,}|\\s+")[[1]]
            col_names <- col_names[nchar(col_names) > 0]
        }
    }

    # 5. Parse data rows -------------------------------------------------------
    row_names <- c()
    value_list <- list()

    for (dline in data_lines) {
        parts <- strsplit(dline, "\u2502")[[1]]
        if (length(parts) < 2) next
        rn <- .trim(parts[1])
        val_str <- paste(parts[-1], collapse = " ")
        vals <- strsplit(.trim(val_str), "\\s+")[[1]]
        vals <- vals[nchar(vals) > 0]
        vals_num <- as.numeric(vals)

        if (nchar(rn) == 0 && all(is.na(vals_num))) next

        row_names <- c(row_names, rn)
        value_list[[length(value_list) + 1]] <- vals_num
    }

    # 6. Build matrix ----------------------------------------------------------
    mat <- NULL
    if (length(value_list) > 0) {
        # Pad to equal lengths
        max_len <- max(sapply(value_list, length))
        value_list <- lapply(value_list, function(v) {
            if (length(v) < max_len) c(v, rep(NA_real_, max_len - length(v)))
            else v
        })
        mat <- do.call(rbind, value_list)
        rownames(mat) <- row_names
        if (!is.null(col_names) && length(col_names) == ncol(mat)) {
            colnames(mat) <- col_names
        }
    }

    list(col_names = col_names, row_names = row_names,
         values = mat, end_idx = bot_idx)
}


# ---- Main parser function ----------------------------------------------------

parse_dynare_log <- function(log_file) {

    lines <- readLines(log_file, warn = FALSE, encoding = "UTF-8")
    n <- length(lines)

    result <- list(
        file             = log_file,
        solved           = FALSE,
        bk_failed        = FALSE,
        bk_message       = NULL,
        steady_state     = NULL,
        decision_rules   = NULL,
        moments          = NULL,
        variance_decomp  = NULL,
        correlation      = NULL,
        autocorrelation  = NULL,
        is_estimation    = FALSE,
        command          = NA_character_,
        preprocessing_ok = FALSE,
        error_message    = NULL
    )

    # Check preprocessing completion
    if (any(grepl("End of preprocessing|End parser", lines))) {
        result$preprocessing_ok <- TRUE
    }

    # Detect BK failure
    bk_fail_patterns <- c(
        "Blanchard.*Kahn.*not\\s+(satisfied|met)",
        "BK.*not\\s+satisfied",
        "no\\s+stable\\s+equilibrium",
        "indeterminacy",
        "rank\\s+condition.*NOT\\s+satisfied",
        "too\\s+many\\s+unstable",
        "could\\s+not\\s+solve",
        "There\\s+are\\s+\\d+\\s+eigenvalue.*but.*forward"
    )
    for (pat in bk_fail_patterns) {
        m <- grep(pat, lines, ignore.case = TRUE, value = TRUE)
        if (length(m) > 0) {
            result$bk_failed  <- TRUE
            result$bk_message <- m[1]
            break
        }
    }

    # Section presence flags
    has_ss   <- any(grepl("Steady state", lines))
    has_dr   <- any(grepl("Coefficients of approximate solution", lines))
    has_mom  <- any(grepl("THEORETICAL MOMENTS", lines))
    has_vdec <- any(grepl("VARIANCE DECOMPOSITION", lines))
    has_corr <- any(grepl("CORRELATION MATRIX", lines))
    has_ac   <- any(grepl("AUTOCORRELATION COEFFICIENTS", lines))

    if (has_mom || has_dr) result$solved <- TRUE

    # If preprocessing ok but no output and no explicit BK message
    if (result$preprocessing_ok && !result$solved && !result$bk_failed) {
        fail_pats <- c(
            "steady\\s+state.*not\\s+found",
            "impossible\\s+to\\s+find\\s+the\\s+steady\\s+state",
            "Error", "FAILED"
        )
        for (pat in fail_pats) {
            m <- grep(pat, lines, ignore.case = FALSE, value = TRUE)
            if (length(m) > 0) {
                result$bk_failed  <- TRUE
                result$bk_message <- m[1]
                break
            }
        }
    }

    # ---- Parse Steady State --------------------------------------------------
    if (has_ss) {
        idx <- grep("^\\s*Steady state\\s*$", lines)
        if (length(idx) > 0) {
            tbl <- .parse_unicode_table(lines, idx[1])
            if (!is.null(tbl$values) && nrow(tbl$values) > 0) {
                ss_vals <- tbl$values[, 1]
                names(ss_vals) <- tbl$row_names
                result$steady_state <- ss_vals
            }
        }
    }

    # ---- Parse Decision Rules ------------------------------------------------
    if (has_dr) {
        idx <- grep("Coefficients of approximate solution", lines)
        if (length(idx) > 0) {
            tbl <- .parse_unicode_table(lines, idx[1])
            if (!is.null(tbl$values)) {
                result$decision_rules <- list(
                    row_names = tbl$row_names,
                    col_names = tbl$col_names,
                    values    = tbl$values
                )
            }
        }
    }

    # ---- Parse Theoretical Moments -------------------------------------------
    if (has_mom) {
        idx <- grep("THEORETICAL MOMENTS", lines)
        if (length(idx) > 0) {
            tbl <- .parse_unicode_table(lines, idx[1])
            if (!is.null(tbl$values) && nrow(tbl$values) > 0) {
                nc <- ncol(tbl$values)
                result$moments <- data.frame(
                    variable = tbl$row_names,
                    mean     = if (nc >= 1) tbl$values[, 1] else NA,
                    std_dev  = if (nc >= 2) tbl$values[, 2] else NA,
                    variance = if (nc >= 3) tbl$values[, 3] else NA,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    # ---- Parse Variance Decomposition ----------------------------------------
    if (has_vdec) {
        idx <- grep("VARIANCE DECOMPOSITION", lines)
        if (length(idx) > 0) {
            tbl <- .parse_unicode_table(lines, idx[1])
            if (!is.null(tbl$values)) {
                result$variance_decomp <- list(
                    row_names = tbl$row_names,
                    col_names = tbl$col_names,
                    values    = tbl$values
                )
            }
        }
    }

    # ---- Parse Correlation Matrix --------------------------------------------
    if (has_corr) {
        idx <- grep("CORRELATION MATRIX", lines)
        if (length(idx) > 0) {
            tbl <- .parse_unicode_table(lines, idx[1])
            if (!is.null(tbl$values)) {
                result$correlation <- list(
                    row_names = tbl$row_names,
                    col_names = tbl$col_names,
                    values    = tbl$values
                )
            }
        }
    }

    # ---- Parse Autocorrelation Coefficients ----------------------------------
    if (has_ac) {
        idx <- grep("AUTOCORRELATION COEFFICIENTS", lines)
        if (length(idx) > 0) {
            tbl <- .parse_unicode_table(lines, idx[1])
            if (!is.null(tbl$values)) {
                result$autocorrelation <- list(
                    row_names = tbl$row_names,
                    col_names = tbl$col_names,
                    values    = tbl$values
                )
            }
        }
    }

    # ---- Detect estimation ---------------------------------------------------
    if (any(grepl("estimation|mode_compute|mh_replic|posterior",
                  lines, ignore.case = TRUE))) {
        result$is_estimation <- TRUE
    }

    result
}


# ---- Parse .mod file to detect command ----------------------------------------

detect_mod_command <- function(mod_file) {
    lines <- readLines(mod_file, warn = FALSE, encoding = "UTF-8")
    full  <- paste(lines, collapse = "\n")

    if (grepl("estimation\\s*[\\(;]", full))
        return("estimation")
    if (grepl("stoch_simul\\s*[\\(;]", full))
        return("stoch_simul")
    if (grepl("perfect_foresight_solver|simul\\s*\\(", full))
        return("perfect_foresight")

    NA_character_
}


# ---- DR row/column name helpers ----------------------------------------------

.clean_dr_row_name <- function(name) {
    name <- .trim(name)
    # ??(x) -> x   (Unicode phi U+03D5 or U+03C6)
    m <- regmatches(name, regexec("[\u03d5\u03c6]\\((.+)\\)", name))[[1]]
    if (length(m) == 2) return(m[2])
    # ASCII phi(x)
    m <- regmatches(name, regexec("phi\\((.+)\\)", name, ignore.case = TRUE))[[1]]
    if (length(m) == 2) return(m[2])
    # Strip trailing _t from shock rows (eps_a_t -> eps_a)
    sub("_t$", "", name)
}

.clean_dr_col_name <- function(name) {
    sub("_t$", "", .trim(name))
}
