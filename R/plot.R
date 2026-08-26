plot.roed_design <- function(x,
                             type = c("design", "power", "selection", "error"),
                             ...) {
  type <- match.arg(type)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)

  if (type == "design") {
    values <- rbind(
      Sample.size = x$design$n,
      Max.toxicity = x$design$m_t,
      Min.efficacy = x$design$m_e
    )
    graphics::barplot(
      values,
      beside = TRUE,
      names.arg = x$design$dose,
      col = c("white", "grey75", "grey35"),
      border = "grey20",
      ylab = "Participants or event count",
      main = "ROED dose-specific design",
      ...
    )
    graphics::legend(
      "topleft",
      legend = c("Sample size", "Maximum toxicities", "Minimum responses"),
      fill = c("white", "grey75", "grey35"),
      border = "grey20",
      bty = "n",
      cex = 0.85
    )
  } else if (type == "power") {
    tab <- x$operating_characteristics$scenario
    graphics::matplot(
      seq_len(nrow(tab)),
      cbind(tab$g1, tab$g2),
      type = "b",
      lty = c(1, 2),
      pch = c(1, 15),
      col = c("grey20", "grey55"),
      xaxt = "n",
      ylim = c(0, 1),
      xlab = "Planning scenario",
      ylab = "Generalized power",
      main = "Scenario-specific generalized power",
      ...
    )
    graphics::axis(1, at = seq_len(nrow(tab)), labels = tab$scenario,
                   las = 2, cex.axis = 0.75)
    graphics::abline(h = x$inputs$target_power, lty = 3, col = "grey45")
    graphics::legend(
      "bottomright", legend = c("G1", "G2", "Target G1"),
      lty = c(1, 2, 3), pch = c(1, 15, NA),
      col = c("grey20", "grey55", "grey45"), bty = "n"
    )
  } else if (type == "selection") {
    pi <- x$operating_characteristics$selection
    graphics::matplot(
      seq_len(nrow(pi)), pi,
      type = "b",
      lty = seq_len(ncol(pi)),
      pch = seq_len(ncol(pi)),
      col = grDevices::gray(seq(0.2, 0.7, length.out = ncol(pi))),
      xaxt = "n",
      ylim = c(0, 1),
      xlab = "Planning scenario",
      ylab = "Selection probability",
      main = "Dose-selection probabilities",
      ...
    )
    graphics::axis(
      1, at = seq_len(nrow(pi)),
      labels = rownames(pi), las = 2, cex.axis = 0.75
    )
    graphics::legend(
      "right", legend = colnames(pi),
      lty = seq_len(ncol(pi)), pch = seq_len(ncol(pi)),
      col = grDevices::gray(seq(0.2, 0.7, length.out = ncol(pi))), bty = "n"
    )
  } else {
    errors <- x$design$local_error
    graphics::barplot(
      errors,
      names.arg = x$design$dose,
      col = "grey75",
      border = "grey20",
      ylim = c(0, max(c(errors, x$inputs$alpha)) * 1.1),
      xlab = "Candidate dose",
      ylab = "Worst-case local error",
      main = sprintf(
        "Local errors; attained strong FWER = %.4f",
        x$operating_characteristics$fwer
      ),
      ...
    )
  }
  invisible(x)
}
