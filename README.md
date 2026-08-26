# ROED

ROED implements Robustly Optimal Exact Designs for randomized
dose-optimization trials with co-primary binary toxicity and efficacy
endpoints.

The package is intended for designing and implementing a trial. It is not a
paper-reproduction package and does not include simulation-study comparators or
manuscript figure-generation code.

## Installation

~~~r
# install.packages("remotes")
# remotes::install_github("haohaostats/ROED")
~~~

## Basic workflow

~~~r
library(ROED)

scenarios <- roed_scenarios(
  toxicity = rbind(
    c(0.15, 0.15),
    c(0.35, 0.15)
  ),
  efficacy = rbind(
    c(0.45, 0.20),
    c(0.20, 0.45)
  ),
  correlation = c(0, 0.20),
  scenario_names = c("dose_1_active", "dose_2_active"),
  dose_names = c("Low", "High")
)

fit <- roed_design(
  scenarios = scenarios,
  tox_null = 0.30,
  eff_null = 0.25,
  alpha = 0.05,
  target_power = 0.80,
  n_min = 10,
  n_max = 100
)

summary(fit)
roed_protocol(fit)

roed_select(
  fit,
  toxicity = c(5, 8),
  efficacy = c(25, 31)
)
~~~

The design search is exact, deterministic, and single-threaded. Probability
surfaces are cached internally during a search; users do not need to configure
or manage the cache.

## User-facing functions

* roed_scenarios() defines clinically plausible planning or sensitivity
  scenarios.
* roed_design() constructs the protocol-ready exact design.
* roed_oc() evaluates exact operating characteristics without re-optimizing.
* roed_protocol() returns a decision table for the protocol or analysis plan.
* roed_select() applies the final integer rules to completed trial counts.
* roed_resume() continues an interrupted search from an explicit checkpoint.

The package does not expose internal cache settings and does not start parallel
workers or background processes.
