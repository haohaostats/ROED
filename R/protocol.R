# Produce protocol-ready ROED decision rules.
roed_protocol <- function(design) {
  .validate_design(design)
  rules <- design$design
  out <- data.frame(
    dose = rules$dose,
    target_enrollment = rules$n,
    maximum_toxicities = rules$m_t,
    minimum_responses = rules$m_e,
    toxicity_rule = sprintf("X_T <= %d", rules$m_t),
    efficacy_rule = sprintf("X_E >= %d", rules$m_e),
    selection_rule = sprintf(
      "Select if X_T <= %d and X_E >= %d", rules$m_t, rules$m_e
    ),
    local_worst_error = rules$local_error,
    stringsAsFactors = FALSE
  )
  attr(out, "total_n") <- sum(rules$n)
  attr(out, "fwer") <- design$operating_characteristics$fwer
  attr(out, "min_g1") <- design$operating_characteristics$min_g1
  class(out) <- c("roed_protocol", "data.frame")
  out
}

print.roed_protocol <- function(x, ...) {
  cat("Protocol-ready ROED decision table\n")
  cat(sprintf("  Total target enrollment: %d\n", attr(x, "total_n")))
  cat(sprintf("  Attained strong FWER: %.6f\n", attr(x, "fwer")))
  cat(sprintf("  Attained minimum G1: %.6f\n\n", attr(x, "min_g1")))
  NextMethod("print", x, row.names = FALSE, ...)
  invisible(x)
}

as.data.frame.roed_protocol <- function(x, row.names = NULL, optional = FALSE,
                                        ...) {
  class(x) <- "data.frame"
  x
}

