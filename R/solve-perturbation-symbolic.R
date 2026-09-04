## R/solve-perturbation-symbolic.R
## --------------------------------------------------------------------------
## Analytic (symbolic) forcing-term assembler for the order-4 / order-5
## perturbation solvers.
##
## Replaces the finite-difference "FD-forcing" builder (removed 2026-09) with an
## exact Faa-di-Bruno assembly of the K-th total derivative of the composed
## map  F( dy(z) ),  where
##
##   z   = (x, u)            combined (state, shock) perturbation, dim nz = ns+nu
##   dy  : R^nz -> R^tc      the dynamic compound vector built from the
##                           order-(K-1) policy (current + lead + lag + shock)
##   F   : R^tc -> R^n       the model residual (rows in compiled-eq order)
##
## Because the policy used to build dy only carries terms up to order K-1, the
## K-th derivative of dy vanishes on the current block, so D^K(F o dy) equals
## exactly the part of the perturbation forcing that does NOT involve the
## unknown order-K policy term.  The generalized Sylvester solve
## (.solve_kron_direct) then reconstructs the order-K decision-rule tensor.
## This is the same quantity the FD builder computes by finite differences,
## but assembled analytically (no residual re-evaluation), ~100-1000x faster.
##
## The central routine .fdb_compose_folded() is fully generic in the order K:
## it is used both to compose the lead-block policy (P o Q) and to assemble
## the outer Phi (F o dy).
## --------------------------------------------------------------------------


# Module-level cache for data-independent fold/expand plans (keyed by K, r).
.fdb_cache <- new.env(parent = emptyenv())


# =====================================================================
# Combinatorial helpers: set partitions and multiset permutations
# =====================================================================

#' All set partitions of {1, ..., K}, each as a list of integer blocks.
#' @noRd
.set_partitions <- function(K) {
  if (K == 0L) return(list(list()))
  prev <- .set_partitions(K - 1L)
  out  <- list()
  for (p in prev) {
    # place element K into each existing block, or into a new singleton block
    for (b in seq_along(p)) {
      q <- p
      q[[b]] <- c(q[[b]], K)
      out[[length(out) + 1L]] <- q
    }
    q <- p
    q[[length(q) + 1L]] <- K
    out[[length(out) + 1L]] <- q
  }
  out
}


#' All distinct permutations of a (possibly repeated) integer vector.
#' Returns a list of integer vectors.
#'
#' MEMOIZED (module-level `.mp_cache`): the same base multi-index recurs heavily
#' -- the policy-tensor scatter calls mp() once per endogenous row j for a base
#' that depends only on the mode-tuple (not j), and .build_combined_policy_derivs
#' is invoked several times per solve (order-4/5 x det/sigma) with identical
#' inputs; the compose engine also re-permutes recurring triplet b-vectors.  The
#' recursion calls the cached entry point so sub-multisets are memoized too,
#' collapsing the exponential recursion to near-linear with reuse.  Returned
#' lists are treated read-only by all callers, so sharing cached objects is safe.
#' @noRd
.mp_cache <- new.env(parent = emptyenv())

.multiset_perms <- function(v) {
  n <- length(v)
  if (n <= 1L) return(list(v))
  key    <- paste0(v, collapse = ",")
  cached <- .mp_cache[[key]]
  if (!is.null(cached)) return(cached)
  # recursive: pick each distinct first element
  out  <- list()
  uniq <- unique(v)
  for (e in uniq) {
    pos <- match(e, v)            # remove one occurrence of e
    rest <- v[-pos]
    for (sub in .multiset_perms(rest)) {
      out[[length(out) + 1L]] <- c(e, sub)
    }
  }
  .mp_cache[[key]] <- out
  out
}


# =====================================================================
# Generic Faa-di-Bruno composition derivative
# =====================================================================


# =====================================================================
# Folded (canonical-column) Faa-di-Bruno composition
# =====================================================================
#
# A dense Faa-di-Bruno composition materializes the full symmetric n_out x r^K
# tensor.  Because that tensor is symmetric, only its C(r+K-1, K) canonical
# (sorted-multi-index) columns are distinct -- at K=5, r=12 that is 4368 vs
# 248832 (~57x fewer).  .fdb_compose_folded() computes exactly those canonical
# columns and never allocates an r^K vector or runs the per-instance aperm.
#
# Math (value-identical to the dense routine, up to floating-point summation
# order): for a canonical column c with sorted mode-tuple v, the symmetric
# tensor value at v is
#
#   D[a, v] = sum_shapes sum_{triplet (a,b,val)} val
#               sum_{partition instances pi of shape} sum_{bp in perms(b)}
#                 prod_i  Hlist[[s_i]][ bp_i, denseflat( v[pos of block i] ) ]
#
# i.e. each partition instance assigns the K mode-slots to blocks; block i's
# modes select a (symmetric) entry of the s_i-th inner derivative.  This is the
# same sum the dense code accumulates for the natural column v, evaluated
# directly via factorized Kronecker indexing.  Vectorized over all canonical
# columns at once (the per-(shape,instance,slot) dense flat-index vectors depend
# only on K, r and the partition geometry, so they are precomputed and cached).

#' Data-independent fold plan for (K, r): canonical columns + per-shape,
#' per-instance, per-slot dense column-index vectors (length N_c).
#' @noRd
.fdb_fold_plan <- function(K, r) {
  key <- paste0("plan_", K, "_", r)
  cached <- .fdb_cache[[key]]
  if (!is.null(cached)) return(cached)

  combos <- .sorted_multiindices(r, K)          # N_c x K, rows sorted ascending
  Nc     <- nrow(combos)
  parts  <- .set_partitions(K)
  sigs   <- vapply(parts, function(p)
    paste(sort(vapply(p, length, 1L), decreasing = TRUE), collapse = "-"), "")
  groups <- split(seq_along(parts), sigs)

  shapes <- lapply(groups, function(gi) {
    grp    <- parts[gi]
    m      <- length(grp[[1L]])
    csizes <- sort(vapply(grp[[1L]], length, 1L), decreasing = TRUE)
    insts  <- lapply(grp, function(part) {
      ord    <- order(vapply(part, length, 1L), decreasing = TRUE)
      blocks <- part[ord]                       # blocks in canonical size order
      lapply(seq_len(m), function(i) {
        pos <- sort(blocks[[i]])                # ascending positions
        s   <- csizes[i]
        # combos[, pos] picks ascending positions of an ascending row -> the
        # sub-multiset is already sorted; flatten col-major (mode 1 fastest).
        sub  <- combos[, pos, drop = FALSE]     # N_c x s
        flat <- sub[, 1L]
        if (s > 1L) for (j in 2:s) flat <- flat + (sub[, j] - 1L) * r^(j - 1L)
        as.integer(flat)
      })
    })
    list(m = m, csizes = csizes, insts = insts)
  })

  plan <- list(Nc = Nc, combos = combos, shapes = shapes)
  .fdb_cache[[key]] <- plan
  plan
}


#' Is the compiled Faa-di-Bruno backend available (and not disabled)?
#' @noRd
.HAS_RCPP_FDB <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("fdb_compose_folded_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

#' Faa-di-Bruno composition in folded form: returns n_out x C(r+K-1, K),
#' columns enumerated by .sorted_multiindices(r, K).  Hlist is in the dense
#' format (row b reshaped col-major to rep(r, k)); the folded output holds
#' exactly the canonical columns of the dense symmetric tensor.
#'
#' Dispatches to the compiled backend fdb_compose_folded_cpp (src/fdb_compose.cpp)
#' when available; the pure-R .fdb_compose_folded_R is the fallback (value-identical
#' up to float summation order -- see test-fdb-rcpp-parity.R).  The fold plan is
#' built/cached in R either way.
#' @noRd
.fdb_compose_folded <- function(K, Glist, Hlist, n_out, r) {
  plan <- .fdb_fold_plan(K, r)
  if (.HAS_RCPP_FDB())
    return(fdb_compose_folded_cpp(plan$Nc, n_out, plan$shapes, Glist, Hlist))
  .fdb_compose_folded_R(plan, Glist, Hlist, n_out)
}

#' Pure-R fallback for .fdb_compose_folded (takes the prebuilt fold plan).
#' @noRd
.fdb_compose_folded_R <- function(plan, Glist, Hlist, n_out) {
  Nc   <- plan$Nc
  D    <- matrix(0, n_out, Nc)

  for (sh in plan$shapes) {
    m <- sh$m
    G <- Glist[[m]]
    if (is.null(G) || length(G$val) == 0L) next
    csizes <- sh$csizes
    cols   <- G$cols
    if (is.null(dim(cols))) cols <- matrix(cols, ncol = m)
    Hsz <- lapply(csizes, function(s) Hlist[[s]])

    for (t in seq_along(G$val)) {
      val <- G$val[t]
      if (val == 0) next
      a     <- G$eq[t]
      perms <- .multiset_perms(cols[t, ])
      contrib <- numeric(Nc)
      for (inst in sh$insts) {
        for (bp in perms) {
          p <- Hsz[[1L]][bp[1L], inst[[1L]]]    # length N_c
          if (m > 1L) for (i in 2:m) p <- p * Hsz[[i]][bp[i], inst[[i]]]
          contrib <- contrib + p
        }
      }
      D[a, ] <- D[a, ] + val * contrib
    }
  }
  D
}


#' Dense-expansion index for (K, r): expand_idx[densecol] = canonical column of
#' the sorted tuple of densecol.  So  folded[, expand_idx]  is the full dense
#' n_out x r^K symmetric tensor.  Cached (O(r^K) build, data-independent).
#' @noRd
.fdb_expand_plan <- function(K, r) {
  key <- paste0("expand_", K, "_", r)
  cached <- .fdb_cache[[key]]
  if (!is.null(cached)) return(cached)

  combos <- .sorted_multiindices(r, K)
  Nc     <- nrow(combos)
  rK     <- r^K
  pw     <- r^(seq_len(K) - 1L)                 # col-major place values

  # Digit matrix of every dense column (col-major, mode 1 fastest), 1-based.
  M   <- matrix(0L, rK, K)
  tmp <- 0:(rK - 1L)
  for (j in seq_len(K)) { M[, j] <- tmp %% r + 1L; tmp <- tmp %/% r }
  # Sort each row ascending via a vectorized bubble compare-exchange network.
  if (K > 1L) for (i in seq_len(K - 1L)) for (j in seq_len(K - i)) {
    a <- M[, j]; b <- M[, j + 1L]
    sw <- a > b
    if (any(sw)) { M[sw, j] <- b[sw]; M[sw, j + 1L] <- a[sw] }
  }
  sorted_flat <- as.integer((M - 1L) %*% pw)             # 0-based flat of sorted tuple
  canon_flat  <- as.integer((combos - 1L) %*% pw)        # 0-based flat of each canonical row
  lookup <- integer(rK); lookup[canon_flat + 1L] <- seq_len(Nc)
  expand_idx <- lookup[sorted_flat + 1L]

  .fdb_cache[[key]] <- expand_idx
  expand_idx
}


#' Expand a folded n_out x C(r+K-1,K) matrix to the dense n_out x r^K layout.
#' @noRd
.expand_folded <- function(Dfold, K, r) {
  if (K == 1L) return(Dfold)                    # folded == dense at K = 1
  Dfold[, .fdb_expand_plan(K, r), drop = FALSE]
}


# =====================================================================
# Dense-array -> canonical-triplet conversion
# =====================================================================

#' Convert a dense derivative matrix (n x p^m, col-major over m p-modes) into
#' canonical sparse triplets (b-indices sorted ascending).
#' @noRd
.dense_derivs_to_triplets <- function(M, m, p, tol = 0) {
  n <- nrow(M)
  if (m == 1L) {
    eq <- integer(0); cols <- integer(0); val <- numeric(0)
    for (a in seq_len(n)) {
      row <- M[a, ]
      nz  <- which(abs(row) > tol)
      if (length(nz)) {
        eq   <- c(eq, rep(a, length(nz)))
        cols <- c(cols, nz)
        val  <- c(val, row[nz])
      }
    }
    return(list(eq = eq, cols = matrix(cols, ncol = 1L), val = val))
  }
  # enumerate canonical (sorted, with repetition) multi-indices over 1..p
  combos <- .sorted_multiindices(p, m)        # N x m matrix, each row sorted asc
  # flat col-major index of each combo (mode 1 fastest)
  flat <- combos[, 1L]
  if (m > 1L) for (k in 2:m) flat <- flat + (combos[, k] - 1L) * p^(k - 1L)
  eq <- integer(0); rows <- list(); val <- numeric(0)
  ri <- 0L
  for (a in seq_len(n)) {
    vals <- M[a, flat]
    nz   <- which(abs(vals) > tol)
    for (j in nz) {
      ri <- ri + 1L
      eq[ri]    <- a
      rows[[ri]] <- combos[j, ]
      val[ri]   <- vals[j]
    }
  }
  cols <- if (ri) do.call(rbind, rows) else matrix(integer(0), 0L, m)
  list(eq = eq, cols = cols, val = val)
}


#' Restrict a dense derivative matrix (n x nz^m, col-major over m nz-modes) to
#' its first `keep` input indices in EVERY mode, returning n x keep^m.  Used to
#' drop the shock input-modes of the lead-block outer map P, which multiply the
#' identically-zero shock rows of the inner map.
#' @noRd
.restrict_input_modes <- function(M, m, nz, keep) {
  if (keep == nz) return(M)
  # The kept columns are the col-major flat indices of the sub-block
  # [1:keep, ..., 1:keep] of the nz^m array -- data-independent in (m, nz, keep),
  # so cache the index vector and reduce the per-row array()/[ ] to one gather.
  key <- paste0("restrict_", m, "_", nz, "_", keep)
  idx <- .fdb_cache[[key]]
  if (is.null(idx)) {
    N   <- keep^m
    dig <- matrix(0L, N, m)                 # 0-based digits, mode 1 fastest
    tmp <- 0:(N - 1L)
    for (j in seq_len(m)) { dig[, j] <- tmp %% keep; tmp <- tmp %/% keep }
    pw  <- nz^(seq_len(m) - 1L)
    idx <- as.integer(dig %*% pw) + 1L
    .fdb_cache[[key]] <- idx
  }
  M[, idx, drop = FALSE]
}


#' Canonical (non-decreasing) multi-indices of length m over 1..p.
#' Returns an N x m integer matrix, each row sorted ascending.
#' @noRd
.sorted_multiindices <- function(p, m) {
  out <- list()
  rec <- function(start, prefix) {
    if (length(prefix) == m) { out[[length(out) + 1L]] <<- prefix; return(invisible()) }
    for (i in start:p) rec(i, c(prefix, i))
  }
  rec(1L, integer(0))
  do.call(rbind, out)
}


# =====================================================================
# Combined-z policy derivative tensors
# =====================================================================

#' Build combined-(state,shock) derivative tensors of the policy up to a
#' requested order.
#'
#' The policy y = ss + g(z) with z = (x, u) has symmetric derivative tensors
#'   P1 : n x nz          (= [ghx | ghu])
#'   P2 : n x nz^2        (from ghxx, ghxu, ghuu)
#'   P3 : n x nz^3        (from ghxxx, ghxxu, ghxuu, ghuuu)
#'   P4 : n x nz^4        (from ghxxxx, ghxxxu, ghxxuu, ghxuuu, ghuuuu)
#' returned in Hlist format (row j reshaped col-major to rep(nz,k) is the
#' symmetric k-th derivative of policy component j).
#'
#' @param dr        decision-rules object carrying gh* up to `pol_order`
#' @param ns,nu     #states, #shocks ; nz = ns + nu
#' @param pol_order highest policy order available/needed (2..4)
#' @return list Plist with Plist[[k]] for k = 1..pol_order
#' @noRd
.build_combined_policy_derivs <- function(dr, ns, nu, pol_order) {
  n  <- nrow(dr$ghx)
  nz <- ns + nu
  xs <- seq_len(ns)
  us <- if (nu > 0L) ns + seq_len(nu) else integer(0)
  Plist <- vector("list", pol_order)

  # ---- order 1 ----
  P1 <- matrix(0, n, nz)
  P1[, xs] <- dr$ghx
  if (nu > 0L) P1[, us] <- dr$ghu
  Plist[[1L]] <- P1

  # ---- order 2 ----
  if (pol_order >= 2L) {
    P2 <- matrix(0, n, nz^2)
    for (j in seq_len(n)) {
      A <- array(0, dim = c(nz, nz))
      Gxx <- matrix(dr$ghxx[j, ], ns, ns)
      A[xs, xs] <- Gxx
      if (nu > 0L) {
        # stored ghxu: standard Kron, col=(s-1)*nu+u (exo FAST) -> dims (nu,ns)=[u,s]
        Gxu <- array(dr$ghxu[j, ], dim = c(nu, ns))
        Guu <- matrix(dr$ghuu[j, ], nu, nu)
        A[xs, us] <- t(Gxu)
        A[us, xs] <- Gxu
        A[us, us] <- Guu
      }
      P2[j, ] <- as.numeric(A)
    }
    Plist[[2L]] <- P2
  }

  # ---- order 3 ----
  if (pol_order >= 3L) {
    P3 <- matrix(0, n, nz^3)
    for (j in seq_len(n)) {
      A <- array(0, dim = c(nz, nz, nz))
      Gxxx <- array(dr$ghxxx[j, ], dim = c(ns, ns, ns))
      A[xs, xs, xs] <- Gxxx
      if (nu > 0L) {
        # stored standard Kron (rightmost index FAST):
        #   ghxxu col=(s1-1)*ns*nu+(s2-1)*nu+u -> dims (nu,ns,ns)=[u,s2,s1]
        #   ghxuu col=(s1-1)*nu^2+(u1-1)*nu+u2 -> dims (nu,nu,ns)=[u2,u1,s1]
        Gxxu <- array(dr$ghxxu[j, ], dim = c(nu, ns, ns))
        Gxuu <- array(dr$ghxuu[j, ], dim = c(nu, nu, ns))
        Guuu <- array(dr$ghuuu[j, ], dim = c(nu, nu, nu))
        for (a in xs) for (b in xs) for (k in seq_len(nu)) {
          v <- Gxxu[k, b, a]
          A[a, b, ns + k] <- v; A[a, ns + k, b] <- v; A[ns + k, a, b] <- v
        }
        for (a in xs) for (k in seq_len(nu)) for (l in seq_len(nu)) {
          v <- Gxuu[l, k, a]
          A[a, ns + k, ns + l] <- v; A[ns + k, a, ns + l] <- v
          A[ns + k, ns + l, a] <- v
        }
        A[us, us, us] <- Guuu
      }
      P3[j, ] <- as.numeric(A)
    }
    Plist[[3L]] <- P3
  }

  # ---- order 4 ----
  if (pol_order >= 4L) {
    P4 <- matrix(0, n, nz^4)
    perms4 <- .multiset_perms                    # alias
    for (j in seq_len(n)) {
      A <- array(0, dim = rep(nz, 4))
      Gxxxx <- array(dr$ghxxxx[j, ], dim = rep(ns, 4))
      A[xs, xs, xs, xs] <- Gxxxx
      if (nu > 0L) {
        # stored standard Kron (rightmost index FAST): reshape with reversed dims
        #   ghxxxu (x,x,x,u) -> dims (nu,ns,ns,ns)=[u,s3,s2,s1]
        #   ghxxuu (x,x,u,u) -> dims (nu,nu,ns,ns)=[u2,u1,s2,s1]
        #   ghxuuu (x,u,u,u) -> dims (nu,nu,nu,ns)=[u3,u2,u1,s1]
        Gxxxu <- array(dr$ghxxxu[j, ], dim = c(nu, ns, ns, ns))
        Gxxuu <- array(dr$ghxxuu[j, ], dim = c(nu, nu, ns, ns))
        Gxuuu <- array(dr$ghxuuu[j, ], dim = c(nu, nu, nu, ns))
        Guuuu <- array(dr$ghuuuu[j, ], dim = rep(nu, 4))
        # place each canonical entry into all distinct slot orderings
        for (a in xs) for (b in xs) for (c in xs) for (k in seq_len(nu)) {
          v <- Gxxxu[k, c, b, a]
          if (v == 0) next
          base <- c(a, b, c, ns + k)
          A[do.call(rbind, .multiset_perms(base))] <- v
        }
        for (a in xs) for (b in xs) for (k in seq_len(nu)) for (l in seq_len(nu)) {
          v <- Gxxuu[l, k, b, a]
          if (v == 0) next
          base <- c(a, b, ns + k, ns + l)
          A[do.call(rbind, .multiset_perms(base))] <- v
        }
        for (a in xs) for (k in seq_len(nu)) for (l in seq_len(nu)) for (mm in seq_len(nu)) {
          v <- Gxuuu[mm, l, k, a]
          if (v == 0) next
          base <- c(a, ns + k, ns + l, ns + mm)
          A[do.call(rbind, .multiset_perms(base))] <- v
        }
        A[us, us, us, us] <- Guuuu
      }
      P4[j, ] <- as.numeric(A)
    }
    Plist[[4L]] <- P4
  }

  Plist
}


# =====================================================================
# Compound-vector derivative tensors  dy(z)
# =====================================================================

#' Build the compound-vector derivative tensors DY[[k]] (k = 1..K).
#'
#' DY[[k]] is a (total_cols x nz^k) matrix in Hlist format.  Rows correspond
#' to dynamic compound columns (dyn_col_map): current/lead policy blocks, lag
#' state block, and shock block.
#'
#'   lag state s  (ll=-1): dy = ss + x_s          -> DY1 unit, higher 0
#'   shock k      (ll= 0): dy = u_k               -> DY1 unit, higher 0
#'   current j    (ll= 0): dy = ss + P_j(z)       -> DY[[k]] = Plist[[k]]_j
#'   lead j       (ll=+1): dy = ss + P_j(Q(z))    -> Faa-di-Bruno composition,
#'                          Q(z) = P(z)[state rows]
#'
#' @param K      highest derivative order to build (4 or 5)
#' @param dyn    compiled$dynamic
#' @param Plist  combined policy derivatives (>= order K-1) from
#'               .build_combined_policy_derivs()
#' @param state_idx,endo_names,exo_names,ns,nu  bookkeeping
#' @return list DY with DY[[k]] for k=1..K  (each total_cols x nz^k)
#' @noRd
.build_DY <- function(K, dyn, Plist, state_idx, endo_names, exo_names, ns, nu) {
  total_cols <- dyn$total_cols
  nz  <- ns + nu
  dcm <- dyn$dyn_col_map
  pol_order <- length(Plist)

  # ---- lead-block policy composition  P o Q, derivatives 1..K ----
  # y_lead = P( z' ) with z' = ( x_next, 0 ),  x_next = P(z)[state rows].
  # Outer map g = P : R^nz -> R^n (b-indices over the full combined space nz).
  # Inner map h : R^nz -> R^nz, z |-> ( state-policy(z), 0 ): first ns output
  # rows are the state-policy derivatives, the nu shock rows are zero.
  # Policy is degree pol_order, so derivatives of order > pol_order vanish; pad
  # with explicit zero matrices / empty triplets so all 1..K orders are indexable.
  n  <- length(endo_names)
  # The inner map h(z) = (state-policy(z), 0) is ZERO outside its first ns
  # (state) output rows, so only the ns state input-modes of the outer map P
  # contribute to D^K(P o h).  Restricting the composition to inner-output
  # dimension p = ns (instead of nz) is exact and cuts the Faa-di-Bruno work by
  # (nz/ns)^K -- the difference between feasible and intractable at K = 5.
  # Only the policy rows for FORWARD-looking (lead-lag == +1) endogenous
  # variables are ever read back from `lead[[k]]` below.  Computing the other
  # rows of the (expensive) lead composition is pure waste -- restrict the
  # outer-map triplets to the forward rows, cutting the Faa-di-Bruno triplet
  # loop by n / n_forward (e.g. 5x when 2 of 10 vars are forward).  Exact: the
  # dropped rows are never indexed (they would scatter into no DY column).
  lead_rows <- integer(0)
  for (kc in seq_len(nrow(dcm))) {
    if (dcm$lead_lag[kc] == 1L && !(dcm$name[kc] %in% exo_names)) {
      jj <- which(endo_names == dcm$name[kc])
      if (length(jj) == 1L) lead_rows <- c(lead_rows, jj)
    }
  }
  lead_rows <- sort(unique(lead_rows))

  Qlist <- vector("list", K)               # each ns x nz^k  (nonzero rows of h)
  Glead <- vector("list", K)               # outer map P, b restricted to 1..ns
  for (k in seq_len(K)) {
    if (k <= pol_order) {
      Qlist[[k]] <- Plist[[k]][state_idx, , drop = FALSE]            # ns x nz^k
      Glead[[k]] <- .dense_derivs_to_triplets(
        .restrict_input_modes(Plist[[k]], k, nz, ns), k, ns)         # b over 1..ns
    } else {
      Qlist[[k]] <- matrix(0, ns, nz^k)
      Glead[[k]] <- list(eq = integer(0),
                         cols = matrix(integer(0), 0L, k),
                         val = numeric(0))
    }
  }
  # Drop triplets whose output row is not a forward variable (never read).
  if (length(lead_rows) < n) {
    for (k in seq_len(K)) {
      G <- Glead[[k]]
      if (length(G$val) == 0L) next
      sel <- G$eq %in% lead_rows
      Glead[[k]] <- list(eq   = G$eq[sel],
                         cols = G$cols[sel, , drop = FALSE],
                         val  = G$val[sel])
    }
  }
  lead <- vector("list", K)                 # lead[[k]] : n x nz^k (forward rows only)
  for (k in seq_len(K)) {
    # Fold the (dominant) lead composition to canonical columns, then expand
    # back to the dense nz^k layout the outer compose / DY contract expects.
    lead[[k]] <- .expand_folded(
      .fdb_compose_folded(k, Glead, Qlist, n, nz), k, nz)
  }

  DY <- vector("list", K)
  for (k in seq_len(K)) DY[[k]] <- matrix(0, total_cols, nz^k)

  xslot <- seq_len(ns)
  for (kc in seq_len(nrow(dcm))) {
    c  <- dcm$col[kc]
    nm <- dcm$name[kc]
    ll <- dcm$lead_lag[kc]
    is_exo <- nm %in% exo_names

    if (is_exo) {
      kk <- which(exo_names == nm)
      if (length(kk) == 1L) DY[[1L]][c, ns + kk] <- 1
      next
    }
    j <- which(endo_names == nm)
    if (length(j) != 1L) next

    if (ll == -1L) {
      s <- which(state_idx == j)
      if (length(s) == 1L) DY[[1L]][c, s] <- 1
    } else if (ll == 0L) {
      for (k in seq_len(min(K, pol_order))) DY[[k]][c, ] <- Plist[[k]][j, ]
      # k > pol_order: current block derivative is 0 (policy is degree pol_order)
    } else if (ll == 1L) {
      for (k in seq_len(K)) DY[[k]][c, ] <- lead[[k]][j, ]
    }
  }
  DY
}


# =====================================================================
# Model residual derivative triplets  F1..FK  (compiled-eq row order)
# =====================================================================

#' Assemble the outer-map (model residual) derivative triplets F1..FK.
#' F1 from the Jacobian (dense); F2..FK from the compiled hessian triplets.
#' b-indices range over 1..total_cols.  Rows in compiled-equation order.
#' @noRd
.build_F_triplets <- function(dyn, dy_ss, params, ss, K) {
  n  <- dyn$n_eq
  tc <- dyn$total_cols
  Flist <- vector("list", K)

  # F1: Jacobian (n x tc)
  J <- dyn$jacobian_fn(dy_ss, params, ss)
  Flist[[1L]] <- .dense_derivs_to_triplets(J, 1L, tc)

  trip_fields <- list(
    `2` = list(fn = "hessian2_fn", tr = "hess2_triplets"),
    `3` = list(fn = "hessian3_fn", tr = "hess3_triplets"),
    `4` = list(fn = "hessian4_fn", tr = "hess4_triplets"),
    `5` = list(fn = "hessian5_fn", tr = "hess5_triplets")
  )
  for (m in 2:K) {
    spec <- trip_fields[[as.character(m)]]
    fn   <- dyn[[spec$fn]]
    trip <- dyn[[spec$tr]]
    if (is.null(fn) || is.null(trip) || length(trip) == 0L) {
      Flist[[m]] <- list(eq = integer(0),
                         cols = matrix(integer(0), 0L, m),
                         val = numeric(0))
      next
    }
    vals <- fn(dy_ss, params, ss)
    eq   <- vapply(trip, function(t) t$eq, 1L)
    cols <- t(vapply(trip, function(t) {
      as.integer(c(t$col1, t$col2,
                   if (m >= 3L) t$col3 else NULL,
                   if (m >= 4L) t$col4 else NULL,
                   if (m >= 5L) t$col5 else NULL))
    }, integer(m)))
    Flist[[m]] <- list(eq = eq, cols = cols, val = as.numeric(vals))
  }
  Flist
}


# =====================================================================
# Slice the symmetric Phi (n x nz^K) into the decision-rule blocks
# =====================================================================

#' Slice a symmetric K-tensor forcing matrix (n x nz^K, col-major, mode 1
#' fastest) into the FD-compatible blocks.  Block layout matches
#' phi_type layout: within each block the LAST slot varies fastest
#' (i.e. as.numeric(aperm(sub, K:1))).
#'
#' @param Phi   n x nz^K matrix (rows already in declaration order)
#' @param K     order
#' @param ns,nu dims
#' @return named list of blocks; names depend on K
#' @noRd
.slice_phi_blocks <- function(Phi, K, ns, nu) {
  n  <- nrow(Phi)
  nz <- ns + nu
  xs <- seq_len(ns)
  us <- if (nu > 0L) ns + seq_len(nu) else integer(0)
  rev_perm <- K:1

  extract <- function(n_x) {
    # block with n_x state slots (positions 1..n_x) and (K-n_x) shock slots
    out <- matrix(0, n, ns^n_x * max(nu, 1L)^(K - n_x))
    sel <- c(rep(list(xs), n_x), rep(list(us), K - n_x))
    for (e in seq_len(n)) {
      A   <- array(Phi[e, ], dim = rep(nz, K))
      sub <- do.call(`[`, c(list(A), sel, list(drop = FALSE)))
      out[e, ] <- as.numeric(aperm(sub, rev_perm))
    }
    out
  }

  if (K == 4L) {
    list(xxxx = extract(4L), xxxu = extract(3L), xxuu = extract(2L),
         xuuu = extract(1L), uuuu = extract(0L))
  } else if (K == 5L) {
    list(xxxxx = extract(5L), xxxxu = extract(4L), xxxuu = extract(3L),
         xxuuu = extract(2L), xuuuu = extract(1L), uuuuu = extract(0L))
  } else {
    stop("Only K = 4 or 5 supported.")
  }
}


# =====================================================================
# Top-level analytic forcing assembler
# =====================================================================

#' Analytic replacement for the finite-difference Phi builder at order 4 / 5.
#'
#' Computes the order-K forcing blocks (the part of the perturbation forcing
#' that excludes the unknown order-K policy term) by exact Faa-di-Bruno
#' assembly of D^K( F o dy ).  Output blocks match the FD builder's layout and
#' sign convention (positive forcing; caller negates for the Sylvester solve).
#'
#' @param dyn        compiled$dynamic
#' @param dr         decision-rules object carrying policy up to order K-1
#' @param ss,params  steady state / parameters
#' @param state_idx,endo_names,exo_names,n_s,n_u,n  bookkeeping
#' @param order      K (4 or 5)
#' @param res_perm   permutation mapping compiled-eq rows -> declaration order
#' @return named list of forcing blocks (same names/layout as the FD builder)
#' @noRd
.build_phi_analytic <- function(dyn, dr, ss, params,
                                state_idx, endo_names, exo_names,
                                n_s, n_u, n, order, res_perm = NULL) {
  if (is.null(res_perm)) res_perm <- seq_len(n)
  K  <- order
  nz <- n_s + n_u
  pol_order <- K - 1L

  dy_ss <- .build_dy_ss_o2(list(dynamic = dyn,
                                model = list(varexo_names = exo_names)), ss)

  Plist <- .build_combined_policy_derivs(dr, n_s, n_u, pol_order)
  DY    <- .build_DY(K, dyn, Plist, state_idx, endo_names, exo_names, n_s, n_u)
  Flist <- .build_F_triplets(dyn, dy_ss, params, ss, K)

  # Fold the outer F o dy compose too (DY stays dense, indexed by dense flat),
  # then expand to the dense layout .slice_phi_blocks expects.
  Phi <- .expand_folded(.fdb_compose_folded(K, Flist, DY, n, nz), K, nz)
  Phi <- Phi[res_perm, , drop = FALSE]            # -> declaration order

  .slice_phi_blocks(Phi, K, n_s, n_u)
}
