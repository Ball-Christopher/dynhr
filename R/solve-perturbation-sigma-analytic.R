## R/solve-perturbation-sigma-analytic.R
## --------------------------------------------------------------------------
## Analytic, moment-based forcing for the order-4/5 sigma-correction terms
## (ghxxss, ghss2, ghxxxss, ghxss2).  Replaces the slow finite-difference
## Gauss-Hermite path (.build_phi_fd_sigma): exact for Gaussian shocks and
## machine-precision accurate (the FD path was 4%-200%+ off on strongly
## nonlinear models -- see the exp closed-form golden, test-sigma-exp-golden.R).
##
## Architecture ("compose THEN fold").  Carry the future shock `e` as explicit
## perturbation modes  w = (x_1..x_ns, sigma, e_1..e_nu),  compose F o dy to
## order K = k_x + 2m, then fold the e-moments at the END via the binomial
## identity (h(sigma)=E_eps[R(x,sigma,sigma*eps)], e=sigma*eps linear in sigma):
##
##   Phi_{x^k_x sigma^2m} = sum_{i=0,2,..,2m} C(2m,i)
##                          [d_x^{k_x} d_sigma^{2m-i} d_e^i (F o dy)] : E[eps^i]
##
## with  E[eps^2] = vec(Sigma),  E[eps^4] = Isserlis(Sigma) (3 pairings).
## The resulting Phi feeds the SAME operator solve as the deterministic k-th
## state term:  A_L*g + fp*g*hx^{(x)k} = -Phi  (.solve_kron_direct, k>=1) or
## (A_L+fp)*g = -Phi  (k=0).
##
## Reuses the deterministic symbolic machinery: .build_combined_policy_derivs,
## .dense_derivs_to_triplets, .fdb_compose, .build_F_triplets, .build_dy_ss_o2,
## .multiset_perms.
## --------------------------------------------------------------------------

# ---- Gaussian moment vectors ----------------------------------------------

#' E[eps (x) eps] = vec(Sigma), length nu^2 (column-major).
#' @noRd
.M2_vec <- function(Sigma_e) as.numeric(Sigma_e)

#' E[eps (x) eps (x) eps (x) eps] = Isserlis/Wick sum of the 3 pairings,
#' length nu^4 (column-major, first index fastest).
#' @noRd
.M4_vec <- function(Sigma_e) {
  nu <- nrow(Sigma_e)
  M4 <- numeric(nu^4); idx <- 0L
  for (d in seq_len(nu)) for (cc in seq_len(nu)) for (b in seq_len(nu)) for (a in seq_len(nu)) {
    idx <- idx + 1L
    M4[idx] <- Sigma_e[a, b] * Sigma_e[cc, d] + Sigma_e[a, cc] * Sigma_e[b, d] +
               Sigma_e[a, d] * Sigma_e[b, cc]
  }
  M4
}

# ---- tensor mode helpers --------------------------------------------------

#' Embed a symmetric k-tensor over `p_sub` modes into the FIRST p_sub modes of
#' a k-tensor over `full` modes (remaining modes zero).  M_sub: n x p_sub^k.
#' @noRd
.embed_first <- function(M_sub, k, p_sub, full) {
  if (p_sub == full) return(M_sub)
  n <- nrow(M_sub); out <- matrix(0, n, full^k)
  idxg <- as.matrix(expand.grid(rep(list(seq_len(p_sub)), k)))   # col1 fastest
  for (a in seq_len(n)) {
    A <- array(0, dim = rep(full, k))
    A[idxg] <- M_sub[a, ]
    out[a, ] <- as.numeric(A)
  }
  out
}

#' Select the same subset `sel` of modes on every axis of a symmetric k-tensor.
#' M: n x full^k  ->  n x length(sel)^k  (sub-modes in the order of `sel`).
#' @noRd
.select_modes <- function(M, k, full, sel) {
  n <- nrow(M); p <- length(sel)
  if (k == 0L) return(M)
  out <- matrix(0, n, p^k)
  selL <- rep(list(sel), k)
  for (a in seq_len(n)) {
    A <- array(M[a, ], dim = rep(full, k))
    sub <- do.call(`[`, c(list(A), selL, list(drop = FALSE)))
    out[a, ] <- as.numeric(sub)
  }
  out
}

#' Order-5 deterministic combined-(state,shock) policy derivative block.
#' Returns n x nz^5.  Mirrors the order-4 layout logic in
#' .build_combined_policy_derivs (reversed-dims reshape, place at all perms);
#' .build_combined_policy_derivs itself only goes to order 4.
#' @noRd
.build_P5 <- function(dr, ns, nu, n) {
  nz <- ns + nu; xs <- seq_len(ns)
  mp <- .multiset_perms
  # The mixed blocks are symmetric in their state slots AND in their shock
  # slots, so the old all-combinations loops (ns^a * nu^b iterations) visited
  # every distinct canonical entry ~a!*b! times -- each redundant visit
  # recomputing mp() and re-assigning the same value to the same positions.
  # Iterate only the SORTED (canonical) state/shock tuples: one mp() call per
  # distinct entry, vectorized matrix-index assignment over its permutations.
  # Bit-identical output (same values scattered to the same array cells).
  st4 <- .sorted_multiindices(ns, 4L)
  st3 <- .sorted_multiindices(ns, 3L)
  st2 <- .sorted_multiindices(ns, 2L)
  ut4 <- if (nu > 0L) .sorted_multiindices(nu, 4L) else NULL
  ut3 <- if (nu > 0L) .sorted_multiindices(nu, 3L) else NULL
  ut2 <- if (nu > 0L) .sorted_multiindices(nu, 2L) else NULL
  P5 <- matrix(0, n, nz^5)
  for (j in seq_len(n)) {
    A <- array(0, dim = rep(nz, 5))
    A[xs, xs, xs, xs, xs] <- array(dr$ghxxxxx[j, ], dim = rep(ns, 5))
    if (nu > 0L) {
      Gxxxxu <- array(dr$ghxxxxu[j, ], dim = c(nu, ns, ns, ns, ns))   # [u,s4,s3,s2,s1]
      Gxxxuu <- array(dr$ghxxxuu[j, ], dim = c(nu, nu, ns, ns, ns))   # [u2,u1,s3,s2,s1]
      Gxxuuu <- array(dr$ghxxuuu[j, ], dim = c(nu, nu, nu, ns, ns))   # [u3,u2,u1,s2,s1]
      Gxuuuu <- array(dr$ghxuuuu[j, ], dim = c(nu, nu, nu, nu, ns))   # [u4,u3,u2,u1,s1]
      Guuuuu <- array(dr$ghuuuuu[j, ], dim = rep(nu, 5))
      us <- ns + seq_len(nu)
      # xxxxu: 4 states (sorted) + 1 shock  (G symmetric in the 4 state slots)
      for (si in seq_len(nrow(st4))) { s <- st4[si, ]
        for (k in seq_len(nu)) {
          v <- Gxxxxu[k, s[4L], s[3L], s[2L], s[1L]]; if (v == 0) next
          A[do.call(rbind, mp(c(s, ns + k)))] <- v
        }
      }
      # xxxuu: 3 states (sorted) + 2 shocks (sorted)
      for (si in seq_len(nrow(st3))) { s <- st3[si, ]
        for (ui in seq_len(nrow(ut2))) { u <- ut2[ui, ]
          v <- Gxxxuu[u[2L], u[1L], s[3L], s[2L], s[1L]]; if (v == 0) next
          A[do.call(rbind, mp(c(s, ns + u)))] <- v
        }
      }
      # xxuuu: 2 states (sorted) + 3 shocks (sorted)
      for (si in seq_len(nrow(st2))) { s <- st2[si, ]
        for (ui in seq_len(nrow(ut3))) { u <- ut3[ui, ]
          v <- Gxxuuu[u[3L], u[2L], u[1L], s[2L], s[1L]]; if (v == 0) next
          A[do.call(rbind, mp(c(s, ns + u)))] <- v
        }
      }
      # xuuuu: 1 state + 4 shocks (sorted)
      for (a in xs) {
        for (ui in seq_len(nrow(ut4))) { u <- ut4[ui, ]
          v <- Gxuuuu[u[4L], u[3L], u[2L], u[1L], a]; if (v == 0) next
          A[do.call(rbind, mp(c(a, ns + u)))] <- v
        }
      }
      A[us, us, us, us, us] <- Guuuuu
    }
    P5[j, ] <- as.numeric(A)
  }
  P5
}

#' Lead-policy derivative tensors over A = (state, shock, sigma),
#' qa = ns + nu + 1 modes (sigma last).  Deterministic part to order K, sigma
#' part to order K-1 (the order-K sigma term is the unknown, excluded).
#' Returns list of dense matrices n x qa^k for k = 1..K.
#' @noRd
.build_GFULL_dense <- function(dr, ns, nu, n, K) {
  qa  <- ns + nu + 1L
  sig <- qa
  Pdet <- .build_combined_policy_derivs(dr, ns, nu, min(K, 4L))
  if (K >= 5L) Pdet[[5L]] <- .build_P5(dr, ns, nu, n)
  G <- vector("list", K)
  for (k in seq_len(K)) G[[k]] <- .embed_first(Pdet[[k]], k, ns + nu, qa)

  add_perm <- function(Ak, idxvec, val) {                  # assign val at all perms
    for (pm in .multiset_perms(idxvec)) Ak[matrix(pm, 1L)] <- val
    Ak
  }
  # k=2: ghss at (sig,sig)
  if (K >= 2L && !is.null(dr$ghss)) {
    ghss <- as.numeric(dr$ghss)
    for (j in seq_len(n)) {
      A <- array(G[[2L]][j, ], dim = rep(qa, 2L))
      A[sig, sig] <- A[sig, sig] + ghss[j]
      G[[2L]][j, ] <- as.numeric(A)
    }
  }
  # k=3: ghxss (state,sig,sig), ghuss (shock,sig,sig)
  if (K >= 3L && K - 1L >= 3L) {
    ghxss <- if (!is.null(dr$ghxss)) dr$ghxss else matrix(0, n, ns)
    ghuss <- if (!is.null(dr$ghuss)) dr$ghuss else matrix(0, n, nu)
    for (j in seq_len(n)) {
      A <- array(G[[3L]][j, ], dim = rep(qa, 3L))
      for (s in seq_len(ns)) if (ghxss[j, s] != 0) A <- add_perm(A, c(s, sig, sig), ghxss[j, s])
      for (l in seq_len(nu)) if (ghuss[j, l] != 0) A <- add_perm(A, c(ns + l, sig, sig), ghuss[j, l])
      G[[3L]][j, ] <- as.numeric(A)
    }
  }
  # k=4: ghxxss (state,state,sig,sig), ghss2 (sig^4)   [order-5 build only]
  if (K >= 4L && K - 1L >= 4L) {
    ghxxss <- dr$ghxxss; ghss2 <- as.numeric(dr$ghss2)
    for (j in seq_len(n)) {
      A <- array(G[[4L]][j, ], dim = rep(qa, 4L))
      Mxx <- matrix(ghxxss[j, ], ns, ns)                  # [s2 fast, s1] symmetric
      for (s1 in seq_len(ns)) for (s2 in seq_len(ns)) {
        v <- Mxx[s2, s1]; if (v != 0) A <- add_perm(A, c(s1, s2, sig, sig), v)
      }
      A[matrix(rep(sig, 4L), 1L)] <- A[matrix(rep(sig, 4L), 1L)] + ghss2[j]
      G[[4L]][j, ] <- as.numeric(A)
    }
  }
  G
}

#' Separable lead-block composition  GLEAD[[k]] = D^k_w [ g(x_now(x,sigma), e, sigma) ].
#'
#' The future shock e enters the lead policy g ONLY as its shock-argument, and
#' that dependence is linear and identity (A_shock = e); nothing else depends on
#' e (x_now and sigma do not).  Hence in the Faa-di-Bruno expansion every e w-slot
#' is forced to be a singleton block mapping (via the identity) to one of g's
#' shock-arguments, while the (x,sigma) w-slots compose through the reduced inner
#' map Mr : (x,sigma) -> (x_now(x,sigma), sigma).  So for a canonical w-column with
#' `b` e-modes (shock multi-index S) and `a = k-b` reduced modes,
#'
#'   GLEAD[j, v] = [ D^a_{(x,sigma)} ( G^{(S)} o Mr ) ]_j  at reduced column v_reduced,
#'
#' where G^{(S)}[[p]][j; R] = g^{(p+b)}_j[ R (reduced args), S (shock args) ] is g's
#' (p+b)-th derivative with b of its arguments fixed to the shock modes S.  This
#' composes over ns+1 reduced modes instead of rp = ns+1+nu, the analog of the
#' deterministic .restrict_input_modes (nz -> ns) win.  Value-identical to the
#' dense GFULL o MINNER compose (the e-slots contribute factor 1).
#' @noRd
.build_glead_separated <- function(K, GFULL_d, state_idx, ns, nu, n, lead_rows) {
  qa    <- ns + nu + 1L
  nr    <- ns + 1L                          # reduced modes: x(1..ns), sigma(=nr)
  sig_A <- qa
  rp    <- ns + 1L + nu
  red2qa   <- c(seq_len(ns), sig_A)         # reduced mode r -> qa A-mode
  shock2qa <- ns + seq_len(nu)              # shock index s  -> qa A-mode

  # Reduced inner map Mr[[s]] : nr x nr^s  (state rows = x_now reduced derivs;
  # sigma row = identity at order 1, zero higher).
  Mr <- vector("list", K)
  for (s in seq_len(K)) {
    CURr_s <- .select_modes(GFULL_d[[s]], s, qa, red2qa)      # n x nr^s
    Ms <- matrix(0, nr, nr^s)
    Ms[seq_len(ns), ] <- CURr_s[state_idx, , drop = FALSE]
    if (s == 1L) Ms[nr, nr] <- 1
    Mr[[s]] <- Ms
  }

  GLEAD <- vector("list", K)
  for (k in seq_len(K)) {
    full_combos <- .sorted_multiindices(rp, k)               # Ncf x k (sorted asc)
    Ncf   <- nrow(full_combos)
    pw_rp <- rp^(seq_len(k) - 1L)
    flat_full   <- as.integer((full_combos - 1L) %*% pw_rp)
    col_of_flat <- integer(rp^k); col_of_flat[flat_full + 1L] <- seq_len(Ncf)
    Dfold <- matrix(0, n, Ncf)

    for (b in 0:k) {
      a    <- k - b
      Sset <- if (b == 0L) matrix(integer(0), 1L, 0L) else .sorted_multiindices(nu, b)
      Sw_pw <- if (b > 0L) (ns + 1L) else 0L
      for (si in seq_len(nrow(Sset))) {
        S   <- if (b > 0L) Sset[si, ] else integer(0)
        Sqa <- if (b > 0L) shock2qa[S] else integer(0)
        Sw  <- if (b > 0L) (ns + 1L) + S else integer(0)      # e w-modes (sorted)

        if (a == 0L) {
          # all e-modes: GLEAD value = g^{(k)}[j; S] directly (no reduced compose)
          mflat <- sum((Sqa - 1L) * qa^(seq_len(b) - 1L)) + 1L
          v     <- Sw
          cc    <- col_of_flat[sum((v - 1L) * pw_rp) + 1L]
          Dfold[lead_rows, cc] <- Dfold[lead_rows, cc] + GFULL_d[[k]][lead_rows, mflat]
          next
        }

        # outer map G^{(S)}[[p]] (p = 1..a) over reduced modes, as triplets
        G_S <- vector("list", a)
        for (p in seq_len(a)) {
          m       <- p + b
          rcombos <- .sorted_multiindices(nr, p)              # Nr_p x p (sorted)
          pw_qa   <- qa^(seq_len(m) - 1L)
          eqv <- integer(0); colmat <- vector("list", 0L); valv <- numeric(0); ri <- 0L
          for (rj in seq_len(nrow(rcombos))) {
            R    <- rcombos[rj, ]
            flat <- sum((c(red2qa[R], Sqa) - 1L) * pw_qa) + 1L
            col  <- GFULL_d[[m]][lead_rows, flat]
            nzl  <- which(col != 0)
            for (l in nzl) {
              ri <- ri + 1L
              eqv[ri] <- lead_rows[l]; colmat[[ri]] <- R; valv[ri] <- col[l]
            }
          }
          G_S[[p]] <- list(eq   = eqv,
                           cols = if (ri) do.call(rbind, colmat) else matrix(integer(0), 0L, p),
                           val  = valv)
        }

        red        <- .fdb_compose_folded(a, G_S, Mr, n, nr)  # n x C(nr+a-1,a)
        red_combos <- .sorted_multiindices(nr, a)
        for (rc in seq_len(nrow(red_combos))) {
          v  <- c(red_combos[rc, ], Sw)                       # already sorted (reduced < e)
          cc <- col_of_flat[sum((v - 1L) * pw_rp) + 1L]
          Dfold[, cc] <- Dfold[, cc] + red[, rc]
        }
      }
    }
    GLEAD[[k]] <- .expand_folded(Dfold, k, rp)
  }
  GLEAD
}


#' Full symmetric K-tensor Phi over w = (x, sigma, e), declaration rows.
#' This is the un-folded forcing; .extract_sigma_block2 folds the e-moments.
#' @noRd
.build_phi_sigma_full <- function(dyn, dr, ss, params, state_idx,
                                  endo_names, exo_names, ns, nu, n,
                                  K, res_perm) {
  qa <- ns + nu + 1L; sig_A <- qa
  rp <- ns + 1L + nu                       # w-modes
  sig_w <- ns + 1L                          # w: sigma
  e_w   <- if (nu > 0L) (ns + 2L):(ns + 1L + nu) else integer(0)  # w: shocks

  GFULL_d  <- .build_GFULL_dense(dr, ns, nu, n, K)

  # GLEAD[[k]][j, ] is only ever read back below for FORWARD-looking
  # (lead_lag == +1) endogenous rows j; the other rows scatter into no DY column.
  dcm0 <- dyn$dyn_col_map
  lead_rows <- integer(0)
  for (kc in seq_len(nrow(dcm0))) {
    if (dcm0$lead_lag[kc] == 1L && !(dcm0$name[kc] %in% exo_names)) {
      jj <- which(endo_names == dcm0$name[kc])
      if (length(jj) == 1L) lead_rows <- c(lead_rows, jj)
    }
  }
  lead_rows <- sort(unique(lead_rows))
  if (length(lead_rows) == 0L) lead_rows <- seq_len(n)

  # CUR: current shock-free policy g(x,0,sigma) over w  (n x rp^k)
  CUR <- vector("list", K)
  for (k in seq_len(K)) {
    sub <- .select_modes(GFULL_d[[k]], k, qa, c(seq_len(ns), sig_A))   # (ns+1) sub-modes
    CUR[[k]] <- .embed_first(sub, k, ns + 1L, rp)                      # first ns+1 w-modes
  }

  # lead block g(x_now(x,sigma), e, sigma) = GFULL o M  (n x rp^k).  The future
  # shock e enters g ONLY as its (linear) shock-argument and nothing else depends
  # on it, so the e-modes are clean pass-throughs: a separable Faa-di-Bruno over
  # the reduced (x,sigma) modes (ns+1, not rp) with e fixed as spectator
  # shock-arguments.  See .build_glead_separated.
  GLEAD <- .build_glead_separated(K, GFULL_d, state_idx, ns, nu, n, lead_rows)

  # DY (total_cols x rp^k)
  total_cols <- dyn$total_cols
  dcm <- dyn$dyn_col_map
  DY <- vector("list", K)
  for (k in seq_len(K)) DY[[k]] <- matrix(0, total_cols, rp^k)
  for (kc in seq_len(nrow(dcm))) {
    c  <- dcm$col[kc]; nm <- dcm$name[kc]; ll <- dcm$lead_lag[kc]
    if (nm %in% exo_names) {
      l <- which(exo_names == nm)
      if (length(l) == 1L) DY[[1L]][c, e_w[l]] <- 1
      next
    }
    j <- which(endo_names == nm); if (length(j) != 1L) next
    if (ll == -1L) {
      s <- which(state_idx == j); if (length(s) == 1L) DY[[1L]][c, s] <- 1
    } else if (ll == 0L) {
      for (k in seq_len(K)) DY[[k]][c, ] <- CUR[[k]][j, ]
    } else if (ll == 1L) {
      for (k in seq_len(K)) DY[[k]][c, ] <- GLEAD[[k]][j, ]
    }
  }

  # outer Faa di Bruno F o DY
  dy_ss <- .build_dy_ss_o2(list(dynamic = dyn,
              model = list(varexo_names = exo_names)), ss)
  Flist <- .build_F_triplets(dyn, dy_ss, params, ss, K)
  # Fold the outer compose too (DY stays dense, indexed by dense col-major flat);
  # expand once to the dense layout the e-moment extractor consumes.
  Phi <- .expand_folded(.fdb_compose_folded(K, Flist, DY, n, rp), K, rp)
  Phi[res_perm, , drop = FALSE]                      # -> declaration order
}

#' Fold the e-moments out of the full Phi tensor for target g_{x^k_x sigma^2m}.
#' Column order over x: last-state-fastest (FD odometer) via aperm(.,k_x:1), so
#' the result matches the FD path's Phi column layout exactly.
#' @noRd
.extract_sigma_block2 <- function(Phi_full, K, ns, nu, k_x, twom, n, M2, M4) {
  rp    <- ns + 1L + nu
  sig_w <- ns + 1L
  e_w   <- if (nu > 0L) (ns + 2L):(ns + 1L + nu) else integer(0)
  x_w   <- seq_len(ns)
  ncols <- if (k_x == 0L) 1L else ns^k_x
  out   <- matrix(0, n, ncols)
  i_seq <- seq(0L, twom, by = 2L)
  for (i in i_seq) {
    j_sig <- twom - i
    cf <- choose(twom, i)
    Mom <- if (i == 0L) 1 else if (i == 2L) M2 else M4
    sel <- c(rep(list(x_w), k_x), rep(list(sig_w), j_sig), rep(list(e_w), i))
    for (e in seq_len(n)) {
      A   <- array(Phi_full[e, ], dim = rep(rp, K))
      sub <- do.call(`[`, c(list(A), sel, list(drop = FALSE)))   # ns(xk_x),1(xj_sig),nu(xi)
      adim <- c(rep(ns, k_x), rep(nu, i))
      sub  <- if (length(adim)) array(sub, dim = adim) else as.numeric(sub)
      if (i == 0L) {
        valx <- as.numeric(sub)
      } else {
        mat  <- matrix(sub, nrow = ncols, ncol = nu^i)           # rows=x (col-major), cols=e
        valx <- as.numeric(mat %*% Mom)
      }
      if (k_x == 0L) {
        out[e, ] <- out[e, ] + cf * valx
      } else {
        ax <- array(valx, dim = rep(ns, k_x))                    # x_1 fastest
        out[e, ] <- out[e, ] + cf * as.numeric(aperm(ax, k_x:1))
      }
    }
  }
  out
}
