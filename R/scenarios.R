# Define planning or evaluation scenarios for ROED.
roed_scenarios <- function(toxicity, efficacy, correlation = 0,
                           scenario_names = NULL, dose_names = NULL) {
  toxicity <- .as_probability_matrix(toxicity, "toxicity")
  efficacy <- .as_probability_matrix(efficacy, "efficacy")
  if (!identical(dim(toxicity), dim(efficacy))) {
    .roed_stop("'toxicity' and 'efficacy' must have identical dimensions.")
  }

  nr <- nrow(toxicity)
  nc <- ncol(toxicity)
  correlation <- .expand_correlation(correlation, nr, nc)
  scenario_names <- .validate_names(
    scenario_names, nr, "Scenario ", "scenario_names"
  )
  dose_names <- .validate_names(dose_names, nc, "Dose ", "dose_names")

  scale <- sqrt(toxicity * (1 - toxicity) * efficacy * (1 - efficacy))
  joint <- toxicity * efficacy + correlation * scale
  lower <- pmax(0, toxicity + efficacy - 1)
  upper <- pmin(toxicity, efficacy)
  bad <- joint < lower - 1e-12 | joint > upper + 1e-12
  if (any(bad)) {
    location <- which(bad, arr.ind = TRUE)[1L, ]
    .roed_stop(
      paste0(
        "The requested correlation is infeasible for scenario '%s', dose ",
        "'%s'. Choose a value compatible with the Bernoulli margins."
      ),
      scenario_names[location[1L]], dose_names[location[2L]]
    )
  }
  joint <- pmin(pmax(joint, lower), upper)
  dimnames(toxicity) <- list(scenario_names, dose_names)
  dimnames(efficacy) <- list(scenario_names, dose_names)
  dimnames(correlation) <- list(scenario_names, dose_names)
  dimnames(joint) <- list(scenario_names, dose_names)

  structure(
    list(
      toxicity = toxicity,
      efficacy = efficacy,
      correlation = correlation,
      joint = joint,
      scenario_names = scenario_names,
      dose_names = dose_names
    ),
    class = "roed_scenarios"
  )
}

print.roed_scenarios <- function(x, ...) {
  cat("ROED scenario set\n")
  cat(sprintf("  Scenarios: %d\n", nrow(x$toxicity)))
  cat(sprintf("  Candidate doses: %d\n", ncol(x$toxicity)))
  cat("  Scenario names:", paste(x$scenario_names, collapse = ", "), "\n")
  invisible(x)
}

as.data.frame.roed_scenarios <- function(x, row.names = NULL, optional = FALSE,
                                         ...) {
  grid <- expand.grid(
    scenario = x$scenario_names,
    dose = x$dose_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  idx_s <- match(grid$scenario, x$scenario_names)
  idx_d <- match(grid$dose, x$dose_names)
  grid$p_t <- x$toxicity[cbind(idx_s, idx_d)]
  grid$p_e <- x$efficacy[cbind(idx_s, idx_d)]
  grid$rho <- x$correlation[cbind(idx_s, idx_d)]
  grid$q <- x$joint[cbind(idx_s, idx_d)]
  grid
}

