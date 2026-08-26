.state_key <- function(p_t, p_e, q) {
  paste(
    formatC(c(p_t, p_e, q), digits = 17L, format = "fg", flag = "#"),
    collapse = "|"
  )
}

.reverse_cumulative <- function(dist) {
  nr <- nrow(dist)
  nc <- ncol(dist)
  lower_t <- apply(dist, 2L, cumsum)
  if (is.null(dim(lower_t))) {
    lower_t <- matrix(lower_t, nrow = nr, ncol = nc)
  }
  reversed <- lower_t[, nc:1L, drop = FALSE]
  tail_e <- t(apply(reversed, 1L, cumsum))
  if (is.null(dim(tail_e))) {
    tail_e <- matrix(tail_e, nrow = nr, ncol = nc)
  }
  tail_e[, nc:1L, drop = FALSE]
}

.joint_surface_series <- function(p_t, p_e, q, n_max) {
  eta <- c(
    eta00 = 1 - p_t - p_e + q,
    eta10 = p_t - q,
    eta01 = p_e - q,
    eta11 = q
  )
  if (any(eta < -1e-12)) {
    .roed_stop("Invalid bivariate Bernoulli cell probabilities.")
  }
  eta <- pmax(eta, 0)

  surfaces <- vector("list", n_max + 1L)
  dist <- matrix(1, 1L, 1L)
  surfaces[[1L]] <- dist
  if (n_max == 0L) {
    return(surfaces)
  }

  for (n in seq_len(n_max)) {
    old_n <- nrow(dist)
    nxt <- matrix(0, old_n + 1L, old_n + 1L)
    idx <- seq_len(old_n)
    nxt[idx, idx] <- nxt[idx, idx] + eta[["eta00"]] * dist
    nxt[idx + 1L, idx] <- nxt[idx + 1L, idx] + eta[["eta10"]] * dist
    nxt[idx, idx + 1L] <- nxt[idx, idx + 1L] + eta[["eta01"]] * dist
    nxt[idx + 1L, idx + 1L] <-
      nxt[idx + 1L, idx + 1L] + eta[["eta11"]] * dist
    dist <- nxt
    surfaces[[n + 1L]] <- .reverse_cumulative(dist)
  }
  surfaces
}

.build_probability_cache <- function(scenarios, n_max, progress = FALSE) {
  cache <- new.env(parent = emptyenv(), hash = TRUE)
  nr <- nrow(scenarios$toxicity)
  nc <- ncol(scenarios$toxicity)
  states <- character()
  state_values <- list()
  for (s in seq_len(nr)) {
    for (j in seq_len(nc)) {
      values <- c(
        scenarios$toxicity[s, j],
        scenarios$efficacy[s, j],
        scenarios$joint[s, j]
      )
      key <- .state_key(values[1L], values[2L], values[3L])
      if (!key %in% states) {
        states <- c(states, key)
        state_values[[key]] <- values
      }
    }
  }

  pb <- NULL
  if (isTRUE(progress) && interactive()) {
    pb <- utils::txtProgressBar(min = 0, max = length(states), style = 3)
    on.exit(close(pb), add = TRUE)
  }
  for (i in seq_along(states)) {
    key <- states[[i]]
    values <- state_values[[key]]
    assign(
      key,
      .joint_surface_series(values[1L], values[2L], values[3L], n_max),
      envir = cache
    )
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  cache
}

.selection_probability <- function(cache, n, m_t, m_e, p_t, p_e, q) {
  key <- .state_key(p_t, p_e, q)
  surfaces <- get(key, envir = cache, inherits = FALSE)
  surfaces[[n + 1L]][m_t + 1L, m_e + 1L]
}

.local_error <- function(n, m_t, m_e, tox_null, eff_null) {
  a_t <- stats::pbinom(m_t, n, tox_null)
  a_e <- stats::pbinom(m_e - 1L, n, eff_null, lower.tail = FALSE)
  c(a_t = a_t, a_e = a_e, a = max(a_t, a_e))
}

.classify_scenarios <- function(scenarios, tox_null, eff_null,
                                require_coverage = FALSE) {
  admissible <- scenarios$toxicity < tox_null &
    scenarios$efficacy > eff_null
  dimnames(admissible) <- dimnames(scenarios$toxicity)
  if (require_coverage) {
    empty <- rowSums(admissible) == 0L
    if (any(empty)) {
      .roed_stop(
        "Every planning scenario must contain at least one admissible dose; failed for: %s.",
        paste(scenarios$scenario_names[empty], collapse = ", ")
      )
    }
    uncovered <- colSums(admissible) == 0L
    if (any(uncovered)) {
      .roed_stop(
        "Every candidate dose must be admissible in at least one planning scenario; failed for: %s.",
        paste(scenarios$dose_names[uncovered], collapse = ", ")
      )
    }
  }
  admissible
}

.evaluate_rules <- function(rules, scenarios, tox_null, eff_null,
                            cache = NULL, progress = FALSE) {
  n_max <- max(rules$n)
  if (is.null(cache)) {
    cache <- .build_probability_cache(scenarios, n_max, progress = progress)
  }
  admissible <- .classify_scenarios(
    scenarios, tox_null, eff_null, require_coverage = FALSE
  )
  nr <- nrow(scenarios$toxicity)
  nc <- ncol(scenarios$toxicity)
  pi <- matrix(NA_real_, nr, nc,
               dimnames = list(scenarios$scenario_names, scenarios$dose_names))
  for (s in seq_len(nr)) {
    for (j in seq_len(nc)) {
      pi[s, j] <- .selection_probability(
        cache,
        rules$n[j], rules$m_t[j], rules$m_e[j],
        scenarios$toxicity[s, j],
        scenarios$efficacy[s, j],
        scenarios$joint[s, j]
      )
    }
  }

  g1 <- g2 <- numeric(nr)
  for (s in seq_len(nr)) {
    adm <- admissible[s, ]
    no_false <- if (any(!adm)) prod(1 - pi[s, !adm]) else 1
    any_true <- if (any(adm)) 1 - prod(1 - pi[s, adm]) else 0
    g1[s] <- no_false * any_true
    g2[s] <- any_true
  }

  local <- vapply(
    seq_len(nc),
    function(j) .local_error(
      rules$n[j], rules$m_t[j], rules$m_e[j], tox_null, eff_null
    ),
    numeric(3L)
  )
  local <- t(local)
  fwer <- 1 - prod(1 - local[, "a"])

  scenario_table <- data.frame(
    scenario = scenarios$scenario_names,
    g1 = g1,
    g2 = g2,
    stringsAsFactors = FALSE
  )
  dose_table <- expand.grid(
    scenario = scenarios$scenario_names,
    dose = scenarios$dose_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  idx_s <- match(dose_table$scenario, scenarios$scenario_names)
  idx_d <- match(dose_table$dose, scenarios$dose_names)
  dose_table$admissible <- admissible[cbind(idx_s, idx_d)]
  dose_table$selection_probability <- pi[cbind(idx_s, idx_d)]
  dose_table$p_t <- scenarios$toxicity[cbind(idx_s, idx_d)]
  dose_table$p_e <- scenarios$efficacy[cbind(idx_s, idx_d)]
  dose_table$rho <- scenarios$correlation[cbind(idx_s, idx_d)]

  list(
    fwer = fwer,
    local = local,
    min_g1 = min(g1),
    mean_g1 = mean(g1),
    min_g2 = min(g2),
    scenario = scenario_table,
    dose = dose_table,
    selection = pi,
    admissible = admissible
  )
}

