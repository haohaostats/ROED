.pi_columns <- function(n_scenarios) {
  paste0("pi_", seq_len(n_scenarios))
}

.generate_candidates <- function(dose, scenarios, admissible, cache,
                                 tox_null, eff_null, alpha, target_power,
                                 n_min, n_max) {
  n_scenarios <- nrow(scenarios$toxicity)
  pi_names <- .pi_columns(n_scenarios)
  sole <- which(rowSums(admissible) == 1L & admissible[, dose])
  rows <- vector("list", 0L)
  counter <- 0L

  for (n in n_min:n_max) {
    m_t_values <- 0:(n - 1L)
    m_e_values <- 1:n
    a_t_values <- stats::pbinom(m_t_values, n, tox_null)
    a_e_values <- stats::pbinom(
      m_e_values - 1L, n, eff_null, lower.tail = FALSE
    )
    m_t_values <- m_t_values[a_t_values <= alpha + 1e-15]
    a_t_values <- a_t_values[a_t_values <= alpha + 1e-15]
    m_e_values <- m_e_values[a_e_values <= alpha + 1e-15]
    a_e_values <- a_e_values[a_e_values <= alpha + 1e-15]
    if (!length(m_t_values) || !length(m_e_values)) {
      next
    }

    for (it in seq_along(m_t_values)) {
      for (ie in seq_along(m_e_values)) {
        a_t <- a_t_values[it]
        a_e <- a_e_values[ie]
        a <- max(a_t, a_e)
        if (a > alpha + 1e-15) {
          next
        }
        pi <- numeric(n_scenarios)
        for (s in seq_len(n_scenarios)) {
          pi[s] <- .selection_probability(
            cache, n, m_t_values[it], m_e_values[ie],
            scenarios$toxicity[s, dose],
            scenarios$efficacy[s, dose],
            scenarios$joint[s, dose]
          )
        }
        if (length(sole) && min(pi[sole]) < target_power - 1e-13) {
          next
        }
        counter <- counter + 1L
        rows[[counter]] <- c(
          n = n,
          m_t = m_t_values[it],
          m_e = m_e_values[ie],
          a_t = a_t,
          a_e = a_e,
          a = a,
          r = -log1p(-a),
          stats::setNames(pi, pi_names)
        )
      }
    }
  }
  if (!length(rows)) {
    return(data.frame())
  }
  out <- as.data.frame(do.call(rbind, rows), check.names = FALSE)
  out$n <- as.integer(out$n)
  out$m_t <- as.integer(out$m_t)
  out$m_e <- as.integer(out$m_e)
  out
}

.dominance_matrix <- function(pool, dose, admissible) {
  pi_names <- .pi_columns(nrow(admissible))
  pi <- as.matrix(pool[, pi_names, drop = FALSE])
  direction <- ifelse(
    matrix(admissible[, dose], nrow(pi), ncol(pi), byrow = TRUE),
    -pi,
    pi
  )
  cbind(n = pool$n, r = pool$r, direction)
}

.prune_dominated <- function(pool, dose, admissible) {
  if (!nrow(pool)) {
    return(pool)
  }
  vec <- .dominance_matrix(pool, dose, admissible)
  order_idx <- order(pool$n, pool$r, pool$m_t, -pool$m_e)
  front <- integer()
  tol <- 2e-14

  for (idx in order_idx) {
    v <- vec[idx, ]
    if (length(front)) {
      f <- vec[front, , drop = FALSE]
      le <- rowSums(f <= matrix(v + tol, nrow(f), ncol(f), byrow = TRUE)) ==
        ncol(f)
      strict <- rowSums(f < matrix(v - tol, nrow(f), ncol(f), byrow = TRUE)) >
        0L
      if (any(le & strict)) {
        next
      }
      same_n <- abs(f[, "n"] - v[["n"]]) <= tol
      dominated <- same_n &
        rowSums(matrix(v, nrow(f), ncol(f), byrow = TRUE) <= f + tol) ==
          ncol(f) &
        rowSums(matrix(v, nrow(f), ncol(f), byrow = TRUE) < f - tol) > 0L
      if (any(dominated)) {
        front <- front[!dominated]
      }
    }
    front <- c(front, idx)
  }
  out <- pool[front, , drop = FALSE]
  out <- out[order(out$n, out$r, out$m_t, -out$m_e), , drop = FALSE]
  rownames(out) <- NULL
  out
}

