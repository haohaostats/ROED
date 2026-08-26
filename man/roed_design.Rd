\name{roed_design}
\alias{roed_design}
\alias{print.roed_design}
\alias{summary.roed_design}
\alias{print.summary.roed_design}
\alias{as.data.frame.roed_design}
\alias{coef.roed_design}
\title{Construct a Robustly Optimal Exact Design}
\usage{
roed_design(
  scenarios,
  tox_null,
  eff_null,
  alpha = 0.05,
  target_power = 0.80,
  n_min = 10L,
  n_max = 150L,
  checkpoint = NULL,
  progress = interactive()
)

\method{print}{roed_design}(x, ...)

\method{summary}{roed_design}(object, ...)

\method{coef}{roed_design}(object, ...)

\method{as.data.frame}{roed_design}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
  \item{scenarios}{A \code{roed_scenarios} object.}
  \item{tox_null}{Maximum acceptable toxicity probability.}
  \item{eff_null}{Minimum clinically relevant efficacy probability.}
  \item{alpha}{Target strong familywise error rate.}
  \item{target_power}{Target minimum generalized power I.}
  \item{n_min,n_max}{Inclusive dose-specific sample-size search limits.}
  \item{checkpoint}{Optional path for an automatically updated search
  checkpoint.}
  \item{progress}{Logical; whether to report search progress.}
  \item{x,object}{A fitted \code{roed_design} object.}
  \item{row.names,optional,...}{Arguments for standard S3 methods.}
}
\value{
\code{roed_design()} returns an object of class \code{roed_design} containing
the dose-specific integer rules, exact operating characteristics, planning
inputs, and search diagnostics.
}
\description{
Jointly determines dose-specific sample sizes and toxicity and efficacy
decision thresholds using a deterministic, single-threaded exact search.
}
\details{
A dose is selected when its observed toxicity count is no larger than its
toxicity threshold and its observed efficacy count is no smaller than its
efficacy threshold. Internal exact-probability caching is automatic and is not
a user option.
}
\examples{
s <- roed_scenarios(
  toxicity = rbind(c(0.10, 0.40), c(0.40, 0.10)),
  efficacy = rbind(c(0.60, 0.10), c(0.10, 0.60)),
  dose_names = c("Low", "High")
)
fit <- roed_design(
  s, tox_null = 0.30, eff_null = 0.30,
  alpha = 0.20, target_power = 0.50,
  n_min = 3, n_max = 20, progress = FALSE
)
summary(fit)
}

