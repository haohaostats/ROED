# Construct a robustly optimal exact design.
roed_design <- function(scenarios, tox_null, eff_null, alpha = 0.05,
                        target_power = 0.80, n_min = 10L, n_max = 150L,
                        checkpoint = NULL, progress = interactive()) {
  started <- proc.time()[["elapsed"]]
  .validate_roed_scenarios(scenarios)
  .assert_scalar_probability(tox_null, "tox_null", open = TRUE)
  .assert_scalar_probability(eff_null, "eff_null", open = TRUE)
  .assert_scalar_probability(alpha, "alpha", open = TRUE)
  .assert_scalar_probability(target_power, "target_power", open = TRUE)
  .assert_whole_scalar(n_min, "n_min", minimum = 1L)
  .assert_whole_scalar(n_max, "n_max", minimum = 1L)
  n_min <- as.integer(n_min)
  n_max <- as.integer(n_max)
  if (n_min > n_max) {
    .roed_stop("'n_min' cannot exceed 'n_max'.")
  }
  if (!is.null(checkpoint) &&
      (!is.character(checkpoint) || length(checkpoint) != 1L ||
       !nzchar(checkpoint))) {
    .roed_stop("'checkpoint' must be NULL or one non-empty file path.")
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    .roed_stop("'progress' must be TRUE or FALSE.")
  }

  arguments <- list(
    scenarios = scenarios,
    tox_null = tox_null,
    eff_null = eff_null,
    alpha = alpha,
    target_power = target_power,
    n_min = n_min,
    n_max = n_max,
    checkpoint = checkpoint,
    progress = progress
  )
  state <- .checkpoint_read(checkpoint, arguments)
  if (!is.null(state) && identical(state$status, "complete")) {
    return(state$result)
  }
  start_total <- if (!is.null(state)) state$last_completed_total + 1L else NULL

  admissible <- .classify_scenarios(
    scenarios, tox_null, eff_null, require_coverage = TRUE
  )
  if (isTRUE(progress)) {
    message("ROED: computing exact probability surfaces")
  }
  cache <- .build_probability_cache(scenarios, n_max, progress = progress)

  j_count <- ncol(scenarios$toxicity)
  raw_counts <- pruned_counts <- integer(j_count)
  pools <- vector("list", j_count)
  for (j in seq_len(j_count)) {
    if (isTRUE(progress)) {
      message(sprintf(
        "ROED: generating candidate rules for %s",
        scenarios$dose_names[j]
      ))
    }
    raw <- .generate_candidates(
      j, scenarios, admissible, cache,
      tox_null, eff_null, alpha, target_power, n_min, n_max
    )
    raw_counts[j] <- nrow(raw)
    if (!nrow(raw)) {
      .roed_stop(
        paste0(
          "No viable dose-level rule was found for '%s'. Increase 'n_max' ",
          "or revise the planning assumptions."
        ),
        scenarios$dose_names[j]
      )
    }
    pools[[j]] <- .prune_dominated(raw, j, admissible)
    pruned_counts[j] <- nrow(pools[[j]])
    if (!pruned_counts[j]) {
      .roed_stop("All candidate rules were removed for '%s'.",
                 scenarios$dose_names[j])
    }
  }

  if (isTRUE(progress)) {
    message(sprintf(
      "ROED: retained %s candidate rules; starting exact search",
      paste(pruned_counts, collapse = ", ")
    ))
  }
  search <- .search_design(
    pools = pools,
    admissible = admissible,
    alpha = alpha,
    target_power = target_power,
    start_total = start_total,
    checkpoint = checkpoint,
    checkpoint_arguments = arguments,
    progress = progress
  )

  rows <- search$rows
  rules <- data.frame(
    dose = scenarios$dose_names,
    n = vapply(rows, function(x) as.integer(x$n), integer(1L)),
    m_t = vapply(rows, function(x) as.integer(x$m_t), integer(1L)),
    m_e = vapply(rows, function(x) as.integer(x$m_e), integer(1L)),
    a_t = vapply(rows, function(x) as.numeric(x$a_t), numeric(1L)),
    a_e = vapply(rows, function(x) as.numeric(x$a_e), numeric(1L)),
    local_error = vapply(rows, function(x) as.numeric(x$a), numeric(1L)),
    stringsAsFactors = FALSE
  )
  oc <- .evaluate_rules(
    rules, scenarios, tox_null, eff_null, cache = cache
  )
  elapsed <- proc.time()[["elapsed"]] - started

  result <- structure(
    list(
      call = match.call(),
      version = "1.0.0",
      inputs = list(
        tox_null = tox_null,
        eff_null = eff_null,
        alpha = alpha,
        target_power = target_power,
        n_min = n_min,
        n_max = n_max
      ),
      scenarios = scenarios,
      admissible = admissible,
      design = rules,
      operating_characteristics = oc,
      search = list(
        raw_candidate_counts = stats::setNames(
          raw_counts, scenarios$dose_names
        ),
        retained_candidate_counts = stats::setNames(
          pruned_counts, scenarios$dose_names
        ),
        elapsed_seconds = unname(elapsed),
        deterministic = TRUE,
        single_threaded = TRUE
      )
    ),
    class = "roed_design"
  )
  .checkpoint_write(
    checkpoint,
    list(status = "complete", arguments = arguments, result = result)
  )
  result
}

# Resume an interrupted design search.
roed_resume <- function(checkpoint) {
  if (!is.character(checkpoint) || length(checkpoint) != 1L ||
      !nzchar(checkpoint) || !file.exists(checkpoint)) {
    .roed_stop("'checkpoint' must identify an existing ROED checkpoint file.")
  }
  state <- readRDS(checkpoint)
  if (!is.list(state) || is.null(state$status) || is.null(state$arguments)) {
    .roed_stop("The checkpoint file is not a valid ROED checkpoint.")
  }
  if (identical(state$status, "complete")) {
    return(state$result)
  }
  state$arguments$checkpoint <- checkpoint
  do.call(roed_design, state$arguments)
}
