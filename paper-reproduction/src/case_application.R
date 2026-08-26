#!/usr/bin/env Rscript

# Reproduce the manuscript's sotorasib application and write CSV only.
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else file.path("results", "case-application")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ROED)

scenarios <- roed_scenarios(
  toxicity = rbind(
    c(.20, .40), c(.20, .20), c(.20, .20),
    c(.20, .40), c(.20, .20), c(.20, .20)
  ),
  efficacy = rbind(
    c(.30, .30), c(.10, .30), c(.30, .30),
    c(.30, .30), c(.10, .30), c(.30, .30)
  ),
  correlation = c(0, 0, 0, .30, .30, .30),
  scenario_names = c(
    "S_L_rho_0", "S_H_rho_0", "S_P_rho_0",
    "S_L_rho_0.30", "S_H_rho_0.30", "S_P_rho_0.30"
  ),
  dose_names = c("240 mg", "960 mg")
)

write_outputs <- function(fit, prefix, observed_t = NULL, observed_e = NULL) {
  design <- data.frame(
    dose = fit$design$dose, n = fit$design$n,
    max_toxicities = fit$design$m_t,
    min_responses = fit$design$m_e,
    local_worst_error = fit$design$local_error
  )
  write.csv(design, file.path(out_dir, paste0(prefix, "_design.csv")), row.names = FALSE)
  oc <- fit$operating_characteristics
  metrics <- data.frame(
    quantity = c("total_n", "strong_fwer", "minimum_g1", "mean_g1", "minimum_g2"),
    target = c(NA, .10, .80, NA, NA),
    attained = c(sum(fit$design$n), oc$fwer, oc$min_g1, oc$mean_g1, oc$min_g2)
  )
  write.csv(metrics, file.path(out_dir, paste0(prefix, "_operating_characteristics.csv")), row.names = FALSE)
  if (!is.null(observed_t)) {
    selected <- as.data.frame(roed_select(fit, toxicity = observed_t, efficacy = observed_e))
    write.csv(selected, file.path(out_dir, paste0(prefix, "_selection.csv")), row.names = FALSE)
  }
}

prospective <- roed_design(
  scenarios, tox_null = .40, eff_null = .10, alpha = .10,
  target_power = .80, n_min = 10, n_max = 150, progress = FALSE
)
write_outputs(prospective, "prospective")

fixed_n <- roed_design(
  scenarios, tox_null = .40, eff_null = .10, alpha = .10,
  target_power = .80, n_min = 104, n_max = 104, progress = FALSE
)
write_outputs(fixed_n, "fixed_n", c(20, 37), c(26, 34))
