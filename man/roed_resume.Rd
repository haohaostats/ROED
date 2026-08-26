\name{roed_resume}
\alias{roed_resume}
\title{Resume an Interrupted ROED Search}
\usage{
roed_resume(checkpoint)
}
\arguments{
  \item{checkpoint}{Path to a checkpoint previously created by
  \code{roed_design()}.}
}
\value{A completed \code{roed_design} object, or the result of continuing the
saved deterministic search.}
\description{
Restarts an interrupted search from the first total sample size not already
exhaustively evaluated.
}

