## R/parse-equations.R
## --------------------------------------------------------------------------
## AST node constructors, tokenizer, expression parser (Pratt-style), and
## the model{...} block parser. The AST nodes produced here are consumed by
## the symbolic Jacobian builder (see jacobian-monolith.R / phase-1d split).
##
## Phase-1 split from parser-monolith.R (no logic changes).
## --------------------------------------------------------------------------

#' Create an AST node for a numeric literal
#' @noRd
ast_number <- function(value) {
  list(type = "number", value = as.numeric(value))
}

#' Create an AST node for a variable reference (endogenous or exogenous)
#'
#' @param name  Variable name.
#' @param lead_lag  Integer: 0 = current period, +1 = lead, -1 = lag, etc.
#' @noRd
ast_variable <- function(name, lead_lag = 0L) {
  list(type = "variable", name = name, lead_lag = as.integer(lead_lag))
}

#' Create an AST node for a parameter reference
#' @noRd
ast_parameter <- function(name) {
  list(type = "parameter", name = name)
}

#' Create an AST node for a binary operation (+, -, *, /, ^)
#' @noRd
ast_binop <- function(op, left, right) {
  list(type = "binop", op = op, left = left, right = right)
}

#' Create an AST node for a unary operation (unary minus or plus)
#' @noRd
ast_unaryop <- function(op, operand) {
  list(type = "unaryop", op = op, operand = operand)
}

#' Create an AST node for a function call
#'
#' @param name  Function name (e.g. "exp", "log", "STEADY_STATE").
#' @param args  List of AST nodes (the arguments).
#' @noRd
ast_funcall <- function(name, args) {
  list(type = "funcall", name = name, args = args)
}

#' Create an AST node for a model-local variable (#-defined)
#' @noRd
ast_local_variable <- function(name) {
  list(type = "local_variable", name = name)
}


#' Convert an AST node to a human-readable mathematical string
#'
#' @param node An AST node (list with $type).
#' @return Character string.
#' @noRd
ast_to_string <- function(node) {
  if (is.null(node)) return("0")
  switch(node$type,
         "number" = {
           format(node$value, scientific = FALSE)
         },
         "variable" = {
           if (node$lead_lag == 0L) {
             node$name
           } else if (node$lead_lag > 0L) {
             paste0(node$name, "(+", node$lead_lag, ")")
           } else {
             paste0(node$name, "(", node$lead_lag, ")")
           }
         },
         "parameter" = {
           node$name
         },
         "local_variable" = {
           paste0("#", node$name)
         },
         "binop" = {
           left_str  <- ast_to_string(node$left)
           right_str <- ast_to_string(node$right)
           if (node$op %in% c("*", "/", "^")) {
             if (node$left$type == "binop" && node$left$op %in% c("+", "-"))
               left_str <- paste0("(", left_str, ")")
             if (node$right$type == "binop" && node$right$op %in% c("+", "-"))
               right_str <- paste0("(", right_str, ")")
           }
           if (node$op == "^" && node$right$type == "binop")
             right_str <- paste0("(", right_str, ")")
           paste0(left_str, " ", node$op, " ", right_str)
         },
         "unaryop" = {
           operand_str <- ast_to_string(node$operand)
           if (node$operand$type %in% c("binop", "unaryop"))
             operand_str <- paste0("(", operand_str, ")")
           paste0(node$op, operand_str)
         },
         "funcall" = {
           args_str <- paste(vapply(node$args, ast_to_string, character(1)),
                             collapse = ", ")
           paste0(node$name, "(", args_str, ")")
         },
         stop("ast_to_string: unknown node type '", node$type, "'")
  )
}


#' Convert an AST to an evaluable R expression string
#'
#' Variable naming convention in output:
#'   - Current period:  y__0   (or simply the endogenous vector index)
#'   - Lead +k:         y__pk  (e.g. y__p1)
#'   - Lag  -k:         y__mk  (e.g. y__m1)
#' Parameters use their original name.
#'
#' @param node        AST node.
#' @param var_names   Character vector of endogenous + exogenous variable names.
#' @param param_names Character vector of parameter names.
#' @return Character string parseable by R's parse().
#' @noRd
ast_to_r_expr <- function(node, var_names = character(0),
                          param_names = character(0)) {
  if (is.null(node)) return("0")
  switch(node$type,
         "number" = {
           format(node$value, scientific = TRUE)
         },
         "variable" = {
           suffix <- if (node$lead_lag == 0L) "__0"
           else if (node$lead_lag > 0L) paste0("__p", node$lead_lag)
           else paste0("__m", abs(node$lead_lag))
           paste0(node$name, suffix)
         },
         "parameter" = {
           node$name
         },
         "local_variable" = {
           paste0("LOCAL_", node$name)
         },
         "binop" = {
           l <- ast_to_r_expr(node$left, var_names, param_names)
           r <- ast_to_r_expr(node$right, var_names, param_names)
           paste0("(", l, " ", node$op, " ", r, ")")
         },
         "unaryop" = {
           operand <- ast_to_r_expr(node$operand, var_names, param_names)
           paste0("(", node$op, operand, ")")
         },
         "funcall" = {
           rname <- switch(node$name,
                           "ln"           = "log",
                           "normcdf"      = "pnorm",
                           "normpdf"      = "dnorm",
                           "cbrt"         = "function(x) x^(1/3)",
                           "STEADY_STATE" = "STEADY_STATE",
                           node$name
           )
           args_str <- paste(
             vapply(node$args,
                    function(a) ast_to_r_expr(a, var_names, param_names),
                    character(1)),
             collapse = ", "
           )
           paste0(rname, "(", args_str, ")")
         },
         stop("ast_to_r_expr: unknown node type '", node$type, "'")
  )
}


#' Collect all variable references from an AST (name + lead_lag pairs)
#'
#' @param node AST node.
#' @return data.frame with columns: name (character), lead_lag (integer).
#' @noRd
ast_collect_variables <- function(node, local_vars = list()) {
  if (is.null(node)) return(data.frame(name = character(0),
                                       lead_lag = integer(0),
                                       stringsAsFactors = FALSE))
  switch(node$type,
         "number"    = data.frame(name = character(0), lead_lag = integer(0),
                                  stringsAsFactors = FALSE),
         "parameter" = data.frame(name = character(0), lead_lag = integer(0),
                                  stringsAsFactors = FALSE),
         "local_variable" = {
           # Resolve the local variable: look up its definition and collect
           # variable references from there, so timing information (lead/lag)
           # on underlying endogenous vars is not lost.
           def <- local_vars[[node$name]]
           if (!is.null(def)) {
             ast_collect_variables(def, local_vars)
           } else {
             data.frame(name = character(0), lead_lag = integer(0),
                        stringsAsFactors = FALSE)
           }
         },
         "variable"  = data.frame(name = node$name, lead_lag = node$lead_lag,
                                  stringsAsFactors = FALSE),
         "binop"     = rbind(ast_collect_variables(node$left, local_vars),
                             ast_collect_variables(node$right, local_vars)),
         "unaryop"   = ast_collect_variables(node$operand, local_vars),
         "funcall"   = do.call(rbind,
                               lapply(node$args, ast_collect_variables,
                                      local_vars = local_vars)),
         stop("ast_collect_variables: unknown node type '", node$type, "'")
  )
}


#' Tokenize an expression string into a list of tokens
#'
#' Each token is a list with fields:
#'   - type:  "NUMBER", "IDENT", "OP", "LPAREN", "RPAREN",
#'            "COMMA", "SEMICOLON", "EQUALS", "HASH", "LBRACKET", "RBRACKET"
#'   - value: the matched text
#'   - pos:   character position in input (for error messages)
#'
#' @param text Character string of a mathematical expression.
#' @return List of token lists.
#' @noRd
tokenize_expr <- function(text) {
  tokens <- list()
  i <- 1L
  n <- nchar(text)

  while (i <= n) {
    ch <- substr(text, i, i)

    # Skip whitespace
    if (grepl("\\s", ch)) { i <- i + 1L; next }

    pos <- i

    # Number: digits, optional decimal, optional exponent
    if (grepl("[0-9]", ch) ||
        (ch == "." && i + 1L <= n && grepl("[0-9]", substr(text, i+1, i+1)))) {
      j <- i
      while (j <= n && grepl("[0-9]", substr(text, j, j))) j <- j + 1L
      if (j <= n && substr(text, j, j) == ".") {
        j <- j + 1L
        while (j <= n && grepl("[0-9]", substr(text, j, j))) j <- j + 1L
      }
      # Exponent
      if (j <= n && substr(text, j, j) %in% c("e", "E")) {
        j <- j + 1L
        if (j <= n && substr(text, j, j) %in% c("+", "-")) j <- j + 1L
        while (j <= n && grepl("[0-9]", substr(text, j, j))) j <- j + 1L
      }
      tokens <- c(tokens, list(list(type = "NUMBER",
                                    value = substr(text, i, j - 1L),
                                    pos = pos)))
      i <- j
      next
    }

    # Identifier: letter or underscore, then alphanumeric or underscore
    if (grepl("[A-Za-z_]", ch)) {
      j <- i
      while (j <= n && grepl("[A-Za-z0-9_]", substr(text, j, j)))
        j <- j + 1L
      tokens <- c(tokens, list(list(type = "IDENT",
                                    value = substr(text, i, j - 1L),
                                    pos = pos)))
      i <- j
      next
    }

    # Two-character relational operators: maximal munch BEFORE single-char.
    # '<=' and '>=' must be checked before '<' and '>' so they are not split.
    if (ch %in% c("<", ">") && i + 1L <= n) {
      ch2 <- substr(text, i + 1L, i + 1L)
      if (ch2 == "=") {
        tok <- list(type = "OP", value = paste0(ch, "="), pos = pos)
        tokens <- c(tokens, list(tok))
        i <- i + 2L
        next
      }
    }

    # Operators and punctuation
    tok <- switch(ch,
                  "+" = list(type = "OP", value = "+"),
                  "-" = list(type = "OP", value = "-"),
                  "*" = list(type = "OP", value = "*"),
                  "/" = list(type = "OP", value = "/"),
                  "^" = list(type = "OP", value = "^"),
                  "<" = list(type = "OP", value = "<"),
                  ">" = list(type = "OP", value = ">"),
                  "(" = list(type = "LPAREN", value = "("),
                  ")" = list(type = "RPAREN", value = ")"),
                  "," = list(type = "COMMA", value = ","),
                  ";" = list(type = "SEMICOLON", value = ";"),
                  "=" = list(type = "EQUALS", value = "="),
                  "#" = list(type = "HASH", value = "#"),
                  "[" = list(type = "LBRACKET", value = "["),
                  "]" = list(type = "RBRACKET", value = "]"),
                  "'" = list(type = "APOSTROPHE", value = "'"),
                  NULL
    )

    if (!is.null(tok)) {
      tok$pos <- pos
      tokens <- c(tokens, list(tok))
      i <- i + 1L
      next
    }

    # Skip unrecognised characters (e.g. $ from LaTeX names)
    i <- i + 1L
  }

  tokens
}


#' Create a new expression parser environment
#'
#' The parser is closure-based: call new_expr_parser() to obtain an
#' environment, then call env$parse_expression() to parse.
#'
#' @param tokens       List of tokens from tokenize_expr().
#' @param var_names    Character vector of declared variable names
#'                     (endogenous + exogenous).
#' @param param_names  Character vector of declared parameter names.
#' @param local_names  Character vector of model-local variable names
#'                     (defined with #).
#' @return An environment with $parse_expression() method.
#' @noRd
new_expr_parser <- function(tokens, var_names = character(0),
                            param_names = character(0),
                            local_names = character(0)) {
  env <- new.env(parent = emptyenv())
  env$tokens <- tokens
  env$pos    <- 1L
  env$var_names   <- var_names
  env$param_names <- param_names
  env$local_names <- local_names

  # -- Helper functions --------------------------------------------------

  env$peek <- function() {
    if (env$pos <= length(env$tokens)) env$tokens[[env$pos]]
    else NULL
  }

  env$advance <- function() {
    tok <- env$peek()
    env$pos <- env$pos + 1L
    tok
  }

  env$expect <- function(type, value = NULL) {
    tok <- env$peek()
    if (is.null(tok))
      stop("Unexpected end of expression; expected ", type,
           if (!is.null(value)) paste0(" '", value, "'"))
    if (tok$type != type || (!is.null(value) && tok$value != value))
      stop("Expected ", type,
           if (!is.null(value)) paste0(" '", value, "'"),
           " but got ", tok$type, " '", tok$value,
           "' at position ", tok$pos)
    env$advance()
  }

  env$at_end <- function() {
    env$pos > length(env$tokens)
  }

  # -- Grammar rules -----------------------------------------------------

  # expression -> relational
  env$parse_expression <- function() {
    env$parse_relational()
  }

  # relational -> additive (('<=' | '>=' | '<' | '>') additive)?
  # Relational operators have LOWER precedence than +/-.
  # Non-associative: a < b < c is a syntax error (unusual in DSGE models).
  env$parse_relational <- function() {
    left <- env$parse_additive()
    tok <- env$peek()
    if (!is.null(tok) && tok$type == "OP" &&
        tok$value %in% c("<=", ">=", "<", ">")) {
      op <- env$advance()$value
      right <- env$parse_additive()
      left <- ast_binop(op, left, right)
    }
    left
  }

  # additive -> multiplicative (('+' | '-') multiplicative)*
  env$parse_additive <- function() {
    left <- env$parse_multiplicative()
    while (!env$at_end()) {
      tok <- env$peek()
      if (!is.null(tok) && tok$type == "OP" && tok$value %in% c("+", "-")) {
        op <- env$advance()$value
        right <- env$parse_multiplicative()
        left <- ast_binop(op, left, right)
      } else {
        break
      }
    }
    left
  }

  # multiplicative -> unary (('*' | '/') unary)*
  env$parse_multiplicative <- function() {
    left <- env$parse_unary()
    while (!env$at_end()) {
      tok <- env$peek()
      if (!is.null(tok) && tok$type == "OP" && tok$value %in% c("*", "/")) {
        op <- env$advance()$value
        right <- env$parse_unary()
        left <- ast_binop(op, left, right)
      } else {
        break
      }
    }
    left
  }

  # unary -> ('-' | '+') unary | power
  #
  # Unary minus binds LOOSER than '^', matching MATLAB/Dynare semantics where
  # `-a^b` parses as `-(a^b)` (e.g. -2^2 == -4), NOT `(-a)^b`.  Placing unary
  # above power (rather than below it) is what makes the base of a power a bare
  # `primary`, so the sign attaches to the whole power expression.  Right
  # recursion here also permits chained/signed forms like `--x` and `2^-3`.
  env$parse_unary <- function() {
    tok <- env$peek()
    if (!is.null(tok) && tok$type == "OP" && tok$value %in% c("-", "+")) {
      op <- env$advance()$value
      operand <- env$parse_unary()
      if (op == "-") {
        # Optimise: -NUMBER -> negative number literal (only a bare literal,
        # never the base of a pending '^', which parse_power has consumed).
        if (operand$type == "number") {
          operand$value <- -operand$value
          return(operand)
        }
        return(ast_unaryop("-", operand))
      }
      return(operand)  # unary + is a no-op
    }
    env$parse_power()
  }

  # power -> primary ('^' unary)?    [right-associative; signed exponents]
  env$parse_power <- function() {
    base <- env$parse_primary()
    tok <- env$peek()
    if (!is.null(tok) && tok$type == "OP" && tok$value == "^") {
      env$advance()
      exponent <- env$parse_unary()  # unary -> power gives right-assoc + sign
      base <- ast_binop("^", base, exponent)
    }
    base
  }

  # primary -> NUMBER
  #          | IDENT '(' lead_lag ')'        -- variable with timing
  #          | IDENT '(' arg_list ')'        -- function call
  #          | IDENT                          -- bare variable, parameter, or local
  #          | '(' expression ')'
  env$parse_primary <- function() {
    tok <- env$peek()
    if (is.null(tok))
      stop("Unexpected end of expression in primary")

    # Number literal
    if (tok$type == "NUMBER") {
      env$advance()
      return(ast_number(tok$value))
    }

    # Identifier
    if (tok$type == "IDENT") {
      name <- env$advance()$value

      # Check if followed by '('
      nxt <- env$peek()
      if (!is.null(nxt) && nxt$type == "LPAREN") {
        is_var <- name %in% env$var_names
        is_param <- name %in% env$param_names

        if (is_var || is_param) {
          # Try to parse as lead/lag: IDENT '(' [+/-] INT ')'
          saved_pos <- env$pos
          result <- {
            env$expect("LPAREN")
            sign <- 1L
            tok2 <- env$peek()
            if (!is.null(tok2) && tok2$type == "OP" &&
                tok2$value %in% c("+", "-")) {
              if (env$advance()$value == "-") sign <- -1L
            }
            tok3 <- env$peek()
            if (!is.null(tok3) && tok3$type == "NUMBER") {
              ll <- as.integer(env$advance()$value) * sign
              env$expect("RPAREN")
              # Parameters with timing (e.g. tau_c(+1)) are constant;
              # return ast_parameter rather than ast_variable.
              if (is_param) ast_parameter(name) else ast_variable(name, ll)
            } else {
              NULL
            }
          }
          if (!is.null(result)) return(result)
          # Backtrack: not a lead/lag, parse as function call
          env$pos <- saved_pos
        }

        # Function call (or STEADY_STATE/EXPECTATION)
        env$expect("LPAREN")

        # EXPECTATION(k)(expr) has two pairs of parens
        if (name == "EXPECTATION") {
          k_tok <- env$peek()
          k_val <- 0L
          if (!is.null(k_tok) && k_tok$type == "NUMBER") {
            k_val <- as.integer(env$advance()$value)
          }
          env$expect("RPAREN")
          env$expect("LPAREN")
          arg <- env$parse_expression()
          env$expect("RPAREN")
          return(ast_funcall("EXPECTATION",
                             list(ast_number(k_val), arg)))
        }

        # Regular function call with comma-separated arguments
        args <- list()
        if (is.null(env$peek()) || env$peek()$type != "RPAREN") {
          args <- c(args, list(env$parse_expression()))
          while (!is.null(env$peek()) &&
                 env$peek()$type == "COMMA") {
            env$advance()  # consume comma
            args <- c(args, list(env$parse_expression()))
          }
        }
        env$expect("RPAREN")

        if (name == "STEADY_STATE") {
          return(ast_funcall("STEADY_STATE", args))
        }

        return(ast_funcall(name, args))
      }

      # Bare identifier: classify
      if (name %in% env$var_names)   return(ast_variable(name, 0L))
      if (name %in% env$param_names) return(ast_parameter(name))
      if (name %in% env$local_names) return(ast_local_variable(name))
      # Unknown identifier -- treat as parameter (may be defined later)
      return(ast_parameter(name))
    }

    # Parenthesised sub-expression
    if (tok$type == "LPAREN") {
      env$advance()
      expr <- env$parse_expression()
      env$expect("RPAREN")
      return(expr)
    }

    stop("Unexpected token '", tok$value, "' (", tok$type,
         ") at position ", tok$pos)
  }

  env
}


#' Parse an expression string into an AST
#'
#' Convenience wrapper around the tokenizer and parser.
#'
#' @param text        Expression string.
#' @param var_names   Declared variable names.
#' @param param_names Declared parameter names.
#' @param local_names Model-local variable names.
#' @return AST node.
#' @noRd
parse_expression <- function(text, var_names = character(0),
                             param_names = character(0),
                             local_names = character(0)) {
  tokens <- tokenize_expr(trimws(text))
  if (length(tokens) == 0) return(ast_number(0))
  parser <- new_expr_parser(tokens, var_names, param_names, local_names)
  result <- parser$parse_expression()
  result
}


#' Parse the model block body into a list of equation objects
#'
#' Each equation object is a list with:
#'   - lhs:    AST of left-hand side
#'   - rhs:    AST of right-hand side (ast_number(0) if no '=')
#'   - tag:    equation tag string (or NA)
#'   - text:   original equation text
#'
#' Also extracts model-local variable definitions (lines starting with #).
#'
#' @param body        Body text of the model block.
#' @param var_names   Declared variable names (endo + exo).
#' @param param_names Declared parameter names.
#' @return A list with:
#'   - equations:  list of equation objects
#'   - local_vars: named list of (name -> AST expression)
#' @noRd
parse_model_block <- function(body, var_names, param_names) {
  equations  <- list()
  local_vars <- list()
  local_names <- character(0)

  # Split on semicolons
  stmts <- strsplit(body, ";")[[1]]
  stmts <- trimws(stmts)
  stmts <- stmts[nchar(stmts) > 0]

  for (stmt in stmts) {
    # Skip lines that are just 'end'
    if (grepl("^\\s*end\\s*$", stmt, ignore.case = TRUE)) next

    # Extract equation tag: [name='...'] or [tag='...']
    # Also store tag_raw (full bracket content) for bind/relax attribute parsing
    tag <- NA_character_
    tag_raw <- NA_character_
    tag_match <- regmatches(stmt,
                            regexec("\\[([^\\]]+)\\]", stmt, perl = TRUE))[[1]]
    if (length(tag_match) > 0 && nchar(tag_match[1]) > 0) {
      tag_raw <- tag_match[2]
      name_match <- regmatches(tag_raw,
                               regexec("name\\s*=\\s*'([^']*)'",
                                       tag_raw, perl = TRUE))[[1]]
      if (length(name_match) > 1) tag <- name_match[2]
      else tag <- tag_raw
      stmt <- sub("\\s*\\[[^\\]]+\\]\\s*", " ", stmt, perl = TRUE)
      stmt <- trimws(stmt)
    }

    if (nchar(stmt) == 0) next

    # Check for model-local variable definition: # local_name = expr
    if (grepl("^\\s*#", stmt)) {
      stmt <- sub("^\\s*#\\s*", "", stmt)
      if (grepl("=", stmt)) {
        parts <- strsplit(stmt, "\\s*=\\s*", perl = TRUE)[[1]]
        lname <- trimws(parts[1])
        lexpr_text <- trimws(paste(parts[-1], collapse = "="))
        local_names <- c(local_names, lname)
        local_vars[[lname]] <- parse_expression(
          lexpr_text, var_names, param_names, local_names
        )
      }
      next
    }

    original_text <- stmt

    # Split on '=' to separate LHS and RHS
    if (grepl("=", stmt)) {
      eq_pos <- regexpr("=", stmt)
      lhs_text <- trimws(substr(stmt, 1, eq_pos - 1))
      rhs_text <- trimws(substr(stmt, eq_pos + 1, nchar(stmt)))
    } else {
      lhs_text <- stmt
      rhs_text <- "0"
    }

    lhs_ast <- parse_expression(lhs_text, var_names, param_names,
                                local_names)
    rhs_ast <- parse_expression(rhs_text, var_names, param_names,
                                local_names)

    equations <- c(equations, list(list(
      lhs     = lhs_ast,
      rhs     = rhs_ast,
      tag     = tag,
      tag_raw = tag_raw,
      text    = original_text
    )))
  }

  list(equations = equations, local_vars = local_vars)
}
