# Apply a completed ROED design to final trial counts.
roed_select <- function(design, toxicity, efficacy, enrolled = NULL) {
  .validate_design(design)
  n_doses <- nrow(design$design)
  toxicity <- .check_count_vector(toxicity, "toxicity", n_doses)
  efficacy <- .check_count_vector(efficacy, "efficacy", n_doses)
  if (is.null(enrolled)) {
    enrolled <- design$design$n
  } else {
    enrolled <- .check_count_vector(enrolled, "enrolled", n_doses)
  }
  if (any(enrolled != design$design$n)) {
    .roed_stop(
      paste0(
        "Final enrollment must equal the fixed dose-specific quotas in the ",
        "design. Use roed_oc() for sensitivity analyses."
      )
    )
  }
  if (any(toxicity > enrolled)) {
    .roed_stop("A toxicity count cannot exceed enrollment.")
  }
  if (any(efficacy > enrolled)) {
    .roed_stop("An efficacy count cannot exceed enrollment.")
  }
  toxicity_pass <- toxicity <= design$design$m_t
  efficacy_pass <- efficacy >= design$design$m_e
  out <- data.frame(
    dose = design$design$dose,
    enrolled = enrolled,
    toxicity = toxicity,
    efficacy = efficacy,
    maximum_toxicities = design$design$m_t,
    minimum_responses = design$design$m_e,
    toxicity_pass = toxicity_pass,
    efficacy_pass = efficacy_pass,
    selected = toxicity_pass & efficacy_pass,
    stringsAsFactors = FALSE
  )
  class(out) <- c("roed_selection", "data.frame")
  out
}

print.roed_selection <- function(x, ...) {
  cat("ROED final dose-selection result\n")
  selected <- x$dose[x$selected]
  if (length(selected)) {
    cat("  Selected:", paste(selected, collapse = ", "), "\n\n")
  } else {
    cat("  No candidate dose satisfied both decision criteria.\n\n")
  }
  NextMethod("print", x, row.names = FALSE, ...)
  invisible(x)
}

as.data.frame.roed_selection <- function(x, row.names = NULL,
                                         optional = FALSE, ...) {
  class(x) <- "data.frame"
  x
}

