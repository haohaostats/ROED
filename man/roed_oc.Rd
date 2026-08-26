\name{roed_oc}
\alias{roed_oc}
\alias{print.roed_oc}
\alias{as.data.frame.roed_oc}
\title{Evaluate Exact Operating Characteristics}
\usage{
roed_oc(design, scenarios = NULL)

\method{print}{roed_oc}(x, ...)

\method{as.data.frame}{roed_oc}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
  \item{design}{A fitted \code{roed_design} object.}
  \item{scenarios}{Optional new \code{roed_scenarios} object. When omitted,
  the planning scenarios are used.}
  \item{x}{A \code{roed_oc} object.}
  \item{row.names,optional,...}{Arguments for standard S3 methods.}
}
\value{
An object of class \code{roed_oc} containing exact FWER, generalized powers,
and dose-specific selection probabilities.
}
\description{
Evaluates a saved design under its planning scenarios or under new sensitivity
scenarios without re-optimizing the design.
}

