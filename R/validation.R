.roed_stop <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

.assert_scalar_probability <- function(x, name, open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (open) {
    ok <- ok && x > 0 && x < 1
  } else {
    ok <- ok && x >= 0 && x <= 1
  }
  if (!ok) {
    interval <- if (open) "(0, 1)" else "[0, 1]"
    .roed_stop("'%s' must be one finite number in %s.", name, interval)
  }
  invisible(TRUE)
}

.assert_whole_scalar <- function(x, name, minimum = 0L) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x == floor(x) && x >= minimum
  if (!ok) {
    .roed_stop("'%s' must be one integer greater than or equal to %d.",
               name, minimum)
  }
  invisible(TRUE)
}

.as_probability_matrix <- function(x, name) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.vector(x) && is.numeric(x)) {
    x <- matrix(x, nrow = 1L)
  }
  if (!is.matrix(x) || !is.numeric(x) || length(x) == 0L) {
    .roed_stop("'%s' must be a non-empty numeric matrix or vector.", name)
  }
  storage.mode(x) <- "double"
  if (any(!is.finite(x)) || any(x < 0 | x > 1)) {
    .roed_stop("Every value in '%s' must be finite and lie in [0, 1].", name)
  }
  x
}

.expand_correlation <- function(x, nr, nc) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x) || any(!is.finite(x))) {
    .roed_stop("'correlation' must contain only finite numeric values.")
  }
  if (length(x) == 1L) {
    out <- matrix(x, nr, nc)
  } else if (is.null(dim(x)) && length(x) == nr) {
    out <- matrix(rep(x, each = nc), nr, nc, byrow = TRUE)
  } else if (is.matrix(x) && identical(dim(x), c(nr, nc))) {
    out <- x
  } else {
    .roed_stop(
      paste0(
        "'correlation' must be a scalar, a vector with one value per ",
        "scenario, or a matrix matching 'toxicity'."
      )
    )
  }
  storage.mode(out) <- "double"
  if (any(out < -1 | out > 1)) {
    .roed_stop("Every correlation must lie in [-1, 1].")
  }
  out
}

.validate_roed_scenarios <- function(x) {
  if (!inherits(x, "roed_scenarios")) {
    .roed_stop("'scenarios' must be created by roed_scenarios().")
  }
  required <- c("toxicity", "efficacy", "correlation", "joint",
                "scenario_names", "dose_names")
  if (!all(required %in% names(x))) {
    .roed_stop("The scenario object is incomplete or corrupted.")
  }
  invisible(TRUE)
}

.validate_names <- function(x, n, prefix, argument) {
  if (is.null(x)) {
    return(paste0(prefix, seq_len(n)))
  }
  if (!is.character(x) || length(x) != n || anyNA(x) ||
      any(!nzchar(x)) || anyDuplicated(x)) {
    .roed_stop(
      "'%s' must contain %d unique, non-empty character values.",
      argument, n
    )
  }
  x
}

.validate_design <- function(x) {
  if (!inherits(x, "roed_design")) {
    .roed_stop("'design' must be an object returned by roed_design().")
  }
  invisible(TRUE)
}

.check_count_vector <- function(x, name, n_doses) {
  if (!is.numeric(x) || length(x) != n_doses || any(!is.finite(x)) ||
      any(x < 0) || any(x != floor(x))) {
    .roed_stop(
      "'%s' must contain one non-negative integer for each of the %d doses.",
      name, n_doses
    )
  }
  as.integer(x)
}

