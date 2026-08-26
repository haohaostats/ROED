print.roed_design <- function(x, ...) {
  cat("Robustly Optimal Exact Design (ROED)\n")
  cat(sprintf("  Candidate doses: %d\n", nrow(x$design)))
  cat(sprintf("  Total sample size: %d\n", sum(x$design$n)))
  cat(sprintf(
    "  Strong FWER: %.6f (target %.3f)\n",
    x$operating_characteristics$fwer, x$inputs$alpha
  ))
  cat(sprintf(
    "  Minimum G1: %.6f (target %.3f)\n",
    x$operating_characteristics$min_g1, x$inputs$target_power
  ))
  cat("\nDose-specific rules:\n")
  display <- x$design[, c("dose", "n", "m_t", "m_e", "local_error")]
  names(display) <- c("Dose", "n", "Max toxicity", "Min efficacy",
                      "Local error")
  print(display, row.names = FALSE, digits = 6)
  invisible(x)
}

summary.roed_design <- function(object, ...) {
  structure(
    list(
      version = object$version,
      n_doses = nrow(object$design),
      n_scenarios = nrow(object$scenarios$toxicity),
      total_n = sum(object$design$n),
      target_fwer = object$inputs$alpha,
      attained_fwer = object$operating_characteristics$fwer,
      target_power = object$inputs$target_power,
      min_g1 = object$operating_characteristics$min_g1,
      mean_g1 = object$operating_characteristics$mean_g1,
      min_g2 = object$operating_characteristics$min_g2,
      design = object$design,
      scenario = object$operating_characteristics$scenario,
      search = object$search
    ),
    class = "summary.roed_design"
  )
}

print.summary.roed_design <- function(x, ...) {
  cat(sprintf("ROED %s design summary\n", x$version))
  cat(sprintf("  Candidate doses:              %d\n", x$n_doses))
  cat(sprintf("  Planning scenarios:           %d\n", x$n_scenarios))
  cat(sprintf("  Total sample size:            %d\n", x$total_n))
  cat(sprintf("  Target strong FWER:           %.3f\n", x$target_fwer))
  cat(sprintf("  Attained strong FWER:         %.6f\n", x$attained_fwer))
  cat(sprintf("  Target minimum G1:            %.3f\n", x$target_power))
  cat(sprintf("  Attained minimum G1:          %.6f\n", x$min_g1))
  cat(sprintf("  Mean G1:                      %.6f\n", x$mean_g1))
  cat(sprintf("  Minimum G2:                   %.6f\n", x$min_g2))
  cat(sprintf("  Search time (seconds):        %.2f\n",
              x$search$elapsed_seconds))
  cat("  Search engine:                deterministic, single-threaded\n")
  cat("\nDose-specific decision rules:\n")
  display <- x$design[, c("dose", "n", "m_t", "m_e", "local_error")]
  names(display) <- c("Dose", "n", "Max toxicity", "Min efficacy",
                      "Local error")
  print(display, row.names = FALSE, digits = 6)
  cat("\nScenario-specific generalized powers:\n")
  print(x$scenario, row.names = FALSE, digits = 6)
  invisible(x)
}

as.data.frame.roed_design <- function(x, row.names = NULL, optional = FALSE,
                                      ...) {
  x$design
}

coef.roed_design <- function(object, ...) {
  object$design[, c("dose", "n", "m_t", "m_e")]
}

