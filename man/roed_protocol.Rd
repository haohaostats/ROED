\name{roed_protocol}
\alias{roed_protocol}
\alias{print.roed_protocol}
\alias{as.data.frame.roed_protocol}
\title{Create a Protocol-Ready Decision Table}
\usage{
roed_protocol(design)

\method{print}{roed_protocol}(x, ...)

\method{as.data.frame}{roed_protocol}(x, row.names = NULL,
  optional = FALSE, ...)
}
\arguments{
  \item{design}{A fitted \code{roed_design} object.}
  \item{x}{A \code{roed_protocol} object.}
  \item{row.names,optional,...}{Arguments for standard S3 methods.}
}
\value{A data-frame subclass containing the fixed arm quotas and final
decision rules.}
\description{
Formats the optimized integer design as a table suitable for a trial protocol
or statistical analysis plan.
}

