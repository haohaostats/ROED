small_scenarios <- function() {
  roed_scenarios(
    toxicity = rbind(c(0.10, 0.40), c(0.40, 0.10)),
    efficacy = rbind(c(0.60, 0.10), c(0.10, 0.60)),
    correlation = 0,
    scenario_names = c("S1", "S2"),
    dose_names = c("D1", "D2")
  )
}

small_design <- function(checkpoint = NULL) {
  roed_design(
    small_scenarios(),
    tox_null = 0.30,
    eff_null = 0.30,
    alpha = 0.20,
    target_power = 0.50,
    n_min = 3,
    n_max = 20,
    checkpoint = checkpoint,
    progress = FALSE
  )
}

