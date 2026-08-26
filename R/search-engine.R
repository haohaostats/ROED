.candidate_metrics <- function(rows, admissible) {
  pi_names <- .pi_columns(nrow(admissible))
  pi <- t(vapply(rows, function(x) {
    as.numeric(x[pi_names])
  }, numeric(length(pi_names))))
  if (nrow(pi) != ncol(admissible)) {
    pi <- matrix(pi, nrow = ncol(admissible), byrow = TRUE)
  }
  pi <- t(pi)

  g1 <- g2 <- numeric(nrow(admissible))
  for (s in seq_len(nrow(admissible))) {
    adm <- admissible[s, ]
    no_false <- if (any(!adm)) prod(1 - pi[s, !adm]) else 1
    any_true <- if (any(adm)) 1 - prod(1 - pi[s, adm]) else 0
    g1[s] <- no_false * any_true
    g2[s] <- any_true
  }
  a <- vapply(rows, function(x) as.numeric(x[["a"]]), numeric(1L))
  n <- vapply(rows, function(x) as.integer(x[["n"]]), integer(1L))
  list(
    total_n = sum(n),
    fwer = 1 - prod(1 - a),
    min_g1 = min(g1),
    mean_g1 = mean(g1),
    min_g2 = min(g2),
    g1 = g1,
    g2 = g2
  )
}

.candidate_objective <- function(rows, metrics) {
  c(
    metrics$total_n,
    -metrics$min_g1,
    -metrics$mean_g1,
    metrics$fwer,
    vapply(rows, function(x) as.numeric(x[["n"]]), numeric(1L)),
    vapply(rows, function(x) as.numeric(x[["m_t"]]), numeric(1L)),
    -vapply(rows, function(x) as.numeric(x[["m_e"]]), numeric(1L))
  )
}

.lexicographically_less <- function(x, y, tol = 1e-13) {
  if (is.null(y)) {
    return(TRUE)
  }
  for (i in seq_along(x)) {
    if (x[i] < y[i] - tol) {
      return(TRUE)
    }
    if (x[i] > y[i] + tol) {
      return(FALSE)
    }
  }
  FALSE
}

.optimistic_g1 <- function(chosen, depth, pools, admissible) {
  n_scenarios <- nrow(admissible)
  pi_names <- .pi_columns(n_scenarios)
  bound <- 1
  for (s in seq_len(n_scenarios)) {
    no_false <- 1
    no_true <- 1
    if (length(chosen)) {
      for (j in seq_along(chosen)) {
        p <- as.numeric(chosen[[j]][[pi_names[s]]])
        if (admissible[s, j]) {
          no_true <- no_true * (1 - p)
        } else {
          no_false <- no_false * (1 - p)
        }
      }
    }
    if (depth < length(pools)) {
      for (j in (depth + 1L):length(pools)) {
        values <- pools[[j]][[pi_names[s]]]
        if (admissible[s, j]) {
          no_true <- no_true * (1 - max(values))
        } else {
          no_false <- no_false * (1 - min(values))
        }
      }
    }
    bound <- min(bound, no_false * (1 - no_true))
  }
  bound
}

.available_allocations <- function(by_n, total_n) {
  j_count <- length(by_n)
  values <- lapply(by_n, function(x) sort(as.integer(names(x))))
  suffix_min <- suffix_max <- numeric(j_count + 1L)
  if (j_count) {
    for (j in j_count:1L) {
      suffix_min[j] <- suffix_min[j + 1L] + min(values[[j]])
      suffix_max[j] <- suffix_max[j + 1L] + max(values[[j]])
    }
  }
  out <- list()
  counter <- 0L
  recurse <- function(j, remaining, prefix) {
    if (j > j_count) {
      if (remaining == 0L) {
        counter <<- counter + 1L
        out[[counter]] <<- prefix
      }
      return(invisible(NULL))
    }
    for (n in values[[j]]) {
      rem <- remaining - n
      if (rem < suffix_min[j + 1L]) {
        break
      }
      if (rem > suffix_max[j + 1L]) {
        next
      }
      recurse(j + 1L, rem, c(prefix, n))
    }
    invisible(NULL)
  }
  recurse(1L, total_n, integer())
  out
}

.checkpoint_write <- function(path, state) {
  if (is.null(path)) {
    return(invisible(NULL))
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  saveRDS(state, temporary, version = 3)
  if (file.exists(path)) {
    unlink(path)
  }
  if (!file.rename(temporary, path)) {
    .roed_stop("Could not write checkpoint file '%s'.", path)
  }
  invisible(path)
}

.checkpoint_read <- function(path, arguments) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  state <- readRDS(path)
  if (!is.list(state) || is.null(state$status) || is.null(state$arguments)) {
    .roed_stop("The checkpoint file is not a valid ROED checkpoint.")
  }
  comparable <- arguments
  comparable$progress <- state$arguments$progress
  comparable$checkpoint <- state$arguments$checkpoint
  if (!identical(comparable, state$arguments)) {
    .roed_stop(
      "The checkpoint was created from different design arguments."
    )
  }
  state
}

.search_design <- function(pools, admissible, alpha, target_power,
                           start_total = NULL, checkpoint = NULL,
                           checkpoint_arguments = NULL, progress = FALSE) {
  j_count <- length(pools)
  by_n <- lapply(pools, function(pool) split(pool, pool$n, drop = TRUE))
  min_total <- sum(vapply(pools, function(x) min(x$n), numeric(1L)))
  max_total <- sum(vapply(pools, function(x) max(x$n), numeric(1L)))
  if (!is.null(start_total)) {
    min_total <- max(min_total, start_total)
  }
  error_budget <- -log1p(-alpha)
  best_rows <- NULL
  best_metrics <- NULL
  best_objective <- NULL

  totals <- min_total:max_total
  pb <- NULL
  if (isTRUE(progress) && interactive()) {
    pb <- utils::txtProgressBar(min = 0, max = length(totals), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (total_index in seq_along(totals)) {
    total_n <- totals[total_index]
    allocations <- .available_allocations(by_n, total_n)

    for (allocation in allocations) {
      local_pools <- lapply(seq_len(j_count), function(j) {
        by_n[[j]][[as.character(allocation[j])]]
      })
      order_idx <- order(vapply(local_pools, nrow, integer(1L)))
      ordered_pools <- local_pools[order_idx]
      ordered_admissible <- admissible[, order_idx, drop = FALSE]

      recurse <- function(depth, chosen, cost) {
        if (cost > error_budget + 1e-14) {
          return(invisible(NULL))
        }
        if (depth < j_count &&
            .optimistic_g1(chosen, depth, ordered_pools,
                           ordered_admissible) < target_power - 1e-12) {
          return(invisible(NULL))
        }
        if (depth == j_count) {
          restored <- vector("list", j_count)
          for (position in seq_len(j_count)) {
            restored[[order_idx[position]]] <- chosen[[position]]
          }
          metrics <- .candidate_metrics(restored, admissible)
          if (metrics$fwer <= alpha + 1e-12 &&
              metrics$min_g1 >= target_power - 1e-12) {
            objective <- .candidate_objective(restored, metrics)
            if (.lexicographically_less(objective, best_objective)) {
              best_rows <<- restored
              best_metrics <<- metrics
              best_objective <<- objective
            }
          }
          return(invisible(NULL))
        }

        pool <- ordered_pools[[depth + 1L]]
        for (row in seq_len(nrow(pool))) {
          candidate <- as.list(pool[row, , drop = FALSE])
          recurse(
            depth + 1L,
            c(chosen, list(candidate)),
            cost + as.numeric(candidate$r)
          )
        }
        invisible(NULL)
      }
      recurse(0L, list(), 0)
    }

    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, total_index)
    }
    if (!is.null(best_rows)) {
      return(list(rows = best_rows, metrics = best_metrics,
                  total_n = total_n))
    }
    .checkpoint_write(
      checkpoint,
      list(
        status = "searching",
        arguments = checkpoint_arguments,
        last_completed_total = total_n
      )
    )
  }
  .roed_stop(
    paste0(
      "No feasible design was found for the requested sample-size range. ",
      "Increase 'n_max' or revise the planning scenarios."
    )
  )
}

