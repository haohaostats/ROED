\name{roed_scenarios}
\alias{roed_scenarios}
\alias{print.roed_scenarios}
\alias{as.data.frame.roed_scenarios}
\title{Define ROED Scenarios}
\usage{
roed_scenarios(
  toxicity,
  efficacy,
  correlation = 0,
  scenario_names = NULL,
  dose_names = NULL
)
}
\arguments{
  \item{toxicity}{Numeric matrix of toxicity probabilities. Rows represent
  scenarios and columns represent candidate doses. A vector defines one
  scenario.}
  \item{efficacy}{Numeric matrix of efficacy probabilities with dimensions
  matching \code{toxicity}.}
  \item{correlation}{A scalar, one value per scenario, or a matrix matching
  \code{toxicity}.}
  \item{scenario_names}{Optional unique scenario labels.}
  \item{dose_names}{Optional unique candidate-dose labels.}
}
\value{An object of class \code{roed_scenarios}.}
\description{
Creates and validates planning or evaluation scenarios, including the
within-participant bivariate Bernoulli association.
}
\examples{
s <- roed_scenarios(
  toxicity = rbind(c(0.10, 0.40), c(0.40, 0.10)),
  efficacy = rbind(c(0.60, 0.10), c(0.10, 0.60)),
  correlation = 0,
  dose_names = c("Low", "High")
)
print(s)
as.data.frame(s)
}

