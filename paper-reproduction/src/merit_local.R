# Local MERIT design search and operating-characteristic engine.
#
# This is an independent implementation of Yang et al. (2024), Statistics in
# Medicine 43:2972--2986, DOI 10.1002/sim.10093.  It follows the paper's
# simulation search, latent bivariate-normal binary outcomes, common decision
# boundaries, generalized power I, and dose-order PAVA transformation.

suppressPackageStartupMessages(library(jsonlite))

row_cumsum <- function(x) {
  if (ncol(x) > 1L) {
    for (k in 2:ncol(x)) x[, k] <- x[, k] + x[, k - 1L]
  }
  storage.mode(x) <- "integer"
  x
}

simulate_bank <- function(J, n_max, nsim, rho, t1, t0, e0, e1) {
  stopifnot(abs(rho) < 1, t1 < t0, e0 < e1)
  qt <- qnorm(c(t1, t0))
  qe <- qnorm(c(e0, e1))
  bank <- vector("list", J)
  for (j in seq_len(J)) {
    zt <- matrix(rnorm(nsim * n_max), nsim, n_max)
    ze <- rho * zt + sqrt(1 - rho^2) * matrix(rnorm(nsim * n_max), nsim, n_max)
    bank[[j]] <- list(
      T1 = row_cumsum(zt <= qt[1]), T0 = row_cumsum(zt <= qt[2]),
      E0 = row_cumsum(ze <= qe[1]), E1 = row_cumsum(ze <= qe[2])
    )
  }
  bank
}

# Vectorized equal-weight PAVA.  For increasing isotonic regression,
# fitted_i = max_{a <= i} min_{b >= i} mean(y_a,...,y_b).
pava_rows <- function(x) {
  x <- as.matrix(x)
  J <- ncol(x)
  out <- matrix(0, nrow(x), J)
  for (i in seq_len(J)) {
    outer <- rep(-Inf, nrow(x))
    for (a in seq_len(i)) {
      inner <- rep(Inf, nrow(x))
      for (b in i:J) {
        avg <- if (a == b) x[, a] else rowMeans(x[, a:b, drop = FALSE])
        inner <- pmin(inner, avg)
      }
      outer <- pmax(outer, inner)
    }
    out[, i] <- outer
  }
  out
}

scenario_counts <- function(bank, n, states) {
  J <- length(states)
  nsim <- nrow(bank[[1]]$T1)
  tox <- eff <- matrix(0, nsim, J)
  for (j in seq_len(J)) {
    tox[, j] <- bank[[j]][[if (states[[j]]$t == 1L) "T1" else "T0"]][, n]
    eff[, j] <- bank[[j]][[if (states[[j]]$e == 1L) "E1" else "E0"]][, n]
  }
  list(T = pava_rows(tox), E = pava_rows(eff))
}

null_states <- function(J) {
  rows <- list()
  for (s in 0:J) for (k in s:J) {
    states <- c(
      replicate(s, list(t = 1L, e = 0L), simplify = FALSE),
      replicate(k - s, list(t = 0L, e = 0L), simplify = FALSE),
      replicate(J - k, list(t = 0L, e = 1L), simplify = FALSE)
    )
    rows[[length(rows) + 1L]] <- list(name = sprintf("H0_s%d_k%d", s, k), states = states)
  }
  rows
}

lfc_states <- function(J) {
  lapply(seq_len(J), function(j) list(
    name = sprintf("H1_%d", j), admissible = j,
    states = c(
      replicate(j - 1L, list(t = 1L, e = 0L), simplify = FALSE),
      list(list(t = 1L, e = 1L)),
      replicate(J - j, list(t = 0L, e = 1L), simplify = FALSE)
    )
  ))
}

null_surface <- function(T, E, n) {
  ans <- matrix(0, n + 1L, n + 1L)
  for (mt in 0:n) {
    eligible_e <- E
    eligible_e[T > mt] <- -Inf
    emax <- apply(eligible_e, 1L, max)
    ans[mt + 1L, ] <- vapply(0:n, function(me) mean(emax >= me), numeric(1))
  }
  ans
}

power1_surface <- function(T, E, n, j) {
  J <- ncol(T)
  low_max <- if (j == 1L) rep(-Inf, nrow(T)) else apply(E[, seq_len(j - 1L), drop = FALSE], 1L, max)
  high_min <- if (j == J) rep(Inf, nrow(T)) else apply(T[, (j + 1L):J, drop = FALSE], 1L, min)
  ans <- matrix(0, n + 1L, n + 1L)
  for (mt in 0:n) {
    tox_ok <- T[, j] <= mt & high_min > mt
    ans[mt + 1L, ] <- vapply(
      0:n, function(me) mean(tox_ok & E[, j] >= me & low_max < me), numeric(1)
    )
  }
  ans
}

search_design <- function(spec, cfg) {
  J <- as.integer(spec$J)
  nsim <- as.integer(spec$simulations)
  n_min <- as.integer(spec$n_min)
  n_max <- as.integer(spec$n_max)
  set.seed(as.integer(spec$seed))
  bank <- simulate_bank(
    J, n_max, nsim, as.numeric(spec$merit_correlation),
    as.numeric(spec$phi_T1), as.numeric(spec$phi_T0),
    as.numeric(spec$phi_E0), as.numeric(spec$phi_E1)
  )
  nulls <- null_states(J)
  alternatives <- lfc_states(J)
  trace <- list()

  for (n in n_min:n_max) {
    alpha_surface <- matrix(0, n + 1L, n + 1L)
    for (sc in nulls) {
      counts <- scenario_counts(bank, n, sc$states)
      alpha_surface <- pmax(alpha_surface, null_surface(counts$T, counts$E, n))
    }
    power_surface <- matrix(1, n + 1L, n + 1L)
    for (sc in alternatives) {
      counts <- scenario_counts(bank, n, sc$states)
      power_surface <- pmin(
        power_surface, power1_surface(counts$T, counts$E, n, sc$admissible)
      )
    }
    feasible <- which(
      alpha_surface <= as.numeric(spec$alpha) + 1e-12 &
        power_surface >= as.numeric(spec$target) - 1e-12,
      arr.ind = TRUE
    )
    trace[[length(trace) + 1L]] <- list(
      n = n, best_power = max(power_surface[alpha_surface <= as.numeric(spec$alpha)], na.rm = TRUE),
      feasible = nrow(feasible)
    )
    if (nrow(feasible)) {
      values <- data.frame(
        row = feasible[, 1], col = feasible[, 2],
        power = power_surface[feasible], alpha = alpha_surface[feasible]
      )
      # n is the primary objective.  At the minimal n, prefer the strongest
      # generalized power I, then the greatest valid use of the error budget.
      values <- values[order(-values$power, -values$alpha, values$col, values$row), ]
      chosen <- values[1, ]
      return(list(
        n_per_arm = n,
        mT = as.integer(chosen$row - 1L),
        mE = as.integer(chosen$col - 1L),
        total_n = J * n,
        estimated_global_type1 = as.numeric(chosen$alpha),
        estimated_global_power1 = as.numeric(chosen$power),
        design_simulations = nsim,
        design_seed = as.integer(spec$seed),
        correlation = as.numeric(spec$merit_correlation),
        isotonic_transformation = TRUE,
        search_trace = trace
      ))
    }
  }
  stop(sprintf("No feasible MERIT design found through n=%d", n_max))
}

oc_scenarios <- function(J) {
  rows <- null_states(J)
  for (k in 0:(J - 1L)) {
    states <- lapply(seq_len(J), function(j) {
      if (j <= k) list(t = 1L, e = 0L) else if (j == k + 1L) list(t = 1L, e = 1L) else list(t = 0L, e = 1L)
    })
    rows[[length(rows) + 1L]] <- list(name = sprintf("W%d", k + 1L), states = states)
  }
  if (J > 1L) for (k in 0:(J - 2L)) {
    states <- lapply(seq_len(J), function(j) {
      if (j <= k) list(t = 1L, e = 0L) else if (j %in% c(k + 1L, k + 2L)) list(t = 1L, e = 1L) else list(t = 0L, e = 1L)
    })
    rows[[length(rows) + 1L]] <- list(name = sprintf("P%d", k + 1L), states = states)
  }
  rows
}

state_parameters <- function(states, spec) {
  unlist(lapply(states, function(st) c(
    if (st$t == 1L) spec$phi_T1 else spec$phi_T0,
    if (st$e == 1L) spec$phi_E1 else spec$phi_E0
  )), use.names = FALSE)
}

run_oc <- function(spec, cfg, design) {
  J <- as.integer(spec$J)
  n <- as.integer(design$n_per_arm)
  nsim <- as.integer(spec$oc_simulations)
  scenarios <- oc_scenarios(J)
  results <- list()
  correlations <- as.numeric(unlist(cfg$oc_correlations))
  for (rho_index in seq_along(correlations)) {
    rho <- correlations[rho_index]
    oc_seed <- as.integer((as.numeric(spec$seed) + 104729 * rho_index) %% 2147483646 + 1)
    set.seed(oc_seed)
    bank <- simulate_bank(
      J, n, nsim, rho, spec$phi_T1, spec$phi_T0, spec$phi_E0, spec$phi_E1
    )
    rows <- vector("list", length(scenarios))
    for (i in seq_along(scenarios)) {
      sc <- scenarios[[i]]
      counts <- scenario_counts(bank, n, sc$states)
      selected <- counts$T <= design$mT & counts$E >= design$mE
      admissible <- vapply(sc$states, function(st) st$t == 1L && st$e == 1L, logical(1))
      if (!any(admissible)) {
        value <- mean(rowSums(selected) > 0L)
        kind <- "null"; metric <- "Type I error"
      } else {
        value <- mean(rowSums(selected[, admissible, drop = FALSE]) > 0L &
                        rowSums(selected[, !admissible, drop = FALSE]) == 0L)
        kind <- "alternative"; metric <- "Power"
      }
      rows[[i]] <- list(
        scenario_index = i, scenario = sc$name, kind = kind,
        parameters = state_parameters(sc$states, spec), metric = metric,
        value = value, average_sample_size = n
      )
    }
    null_values <- vapply(rows[vapply(rows, function(x) x$kind == "null", logical(1))], `[[`, numeric(1), "value")
    power_values <- vapply(rows[vapply(rows, function(x) x$kind == "alternative", logical(1))], `[[`, numeric(1), "value")
    key <- sprintf("rho_%.2f", rho)
    results[[key]] <- list(
      rho = rho, simulations = nsim, seed = oc_seed,
      max_finite_null_type1 = max(null_values), min_power = min(power_values),
      mean_power = mean(power_values), rows = rows
    )
  }
  results
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript merit_local.R input.json output.json")
request <- fromJSON(args[1], simplifyVector = FALSE)
spec <- request$spec
cfg <- request$config
started <- proc.time()[[3]]
design <- search_design(spec, cfg)
oc <- if (isTRUE(spec$run_oc)) run_oc(spec, cfg, design) else NULL
answer <- list(
  engine = list(name = "local MERIT R", version = "1.0.0", language = R.version.string,
                reference_doi = "10.1002/sim.10093"),
  design = design, operating_characteristics = oc,
  elapsed_seconds = unname(proc.time()[[3]] - started)
)
write_json(answer, args[2], auto_unbox = TRUE, pretty = TRUE, digits = 15, null = "null")
