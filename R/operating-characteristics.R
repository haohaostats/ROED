# Evaluate exact operating characteristics.
roed_oc <- function(design, scenarios = NULL) {
  .validate_design(design)
  if (is.null(scenarios)) {
    scenarios <- design$scenarios
    oc <- design$operating_characteristics
  } else {
    .validate_roed_scenarios(scenarios)
    if (ncol(scenarios$toxicity) != nrow(design$design)) {
      .roed_stop(
        "The evaluation scenarios must contain the same number of doses as the design."
      )
    }
    if (!identical(scenarios$dose_names, design$scenarios$dose_names)) {
      .roed_stop(
        "Evaluation dose names and order must match the fitted design."
      )
    }
    oc <- .evaluate_rules(
      design$design,
      scenarios,
      design$inputs$tox_null,
      design$inputs$eff_null
    )
  }
  structure(
    list(
      scenario = oc$scenario,
      dose = oc$dose,
      fwer = oc$fwer,
      min_g1 = oc$min_g1,
      mean_g1 = oc$mean_g1,
      min_g2 = oc$min_g2,
      scenarios = scenarios
    ),
    class = "roed_oc"
  )
}

print.roed_oc <- function(x, ...) {
  cat("Exact ROED operating characteristics\n")
  cat(sprintf("  Strong FWER: %.6f\n", x$fwer))
  cat(sprintf("  Minimum G1: %.6f\n", x$min_g1))
  cat(sprintf("  Mean G1: %.6f\n", x$mean_g1))
  cat(sprintf("  Minimum G2: %.6f\n", x$min_g2))
  invisible(x)
}

as.data.frame.roed_oc <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$dose
}

