\name{roed_select}
\alias{roed_select}
\alias{print.roed_selection}
\alias{as.data.frame.roed_selection}
\title{Apply the Final ROED Dose-Selection Rule}
\usage{
roed_select(design, toxicity, efficacy, enrolled = NULL)

\method{print}{roed_selection}(x, ...)

\method{as.data.frame}{roed_selection}(x, row.names = NULL,
  optional = FALSE, ...)
}
\arguments{
  \item{design}{A fitted \code{roed_design} object.}
  \item{toxicity}{One final toxicity count per candidate dose.}
  \item{efficacy}{One final efficacy count per candidate dose.}
  \item{enrolled}{Optional enrollment counts. These must equal the fixed
  dose-specific design quotas.}
  \item{x}{A \code{roed_selection} object.}
  \item{row.names,optional,...}{Arguments for standard S3 methods.}
}
\value{A \code{roed_selection} data frame identifying doses that satisfy both
integer decision criteria.}
\description{
Applies the protocol-specified final rules to completed trial counts.
}

