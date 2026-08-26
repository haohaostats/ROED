\name{plot.roed_design}
\alias{plot.roed_design}
\title{Plot ROED Design Diagnostics}
\usage{
\method{plot}{roed_design}(
  x,
  type = c("design", "power", "selection", "error"),
  ...
)
}
\arguments{
  \item{x}{A fitted \code{roed_design} object.}
  \item{type}{Diagnostic plot type.}
  \item{...}{Additional graphical arguments.}
}
\value{The design object, invisibly.}
\description{
Produces base-R diagnostic plots for the integer design, generalized powers,
dose-selection probabilities, or exact local errors.
}

