# ROED

ROED implements Robustly Optimal Exact Designs for randomized
dose-optimization trials with co-primary binary toxicity and efficacy
endpoints.

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
  efficacy = c(28, 31)
)
~~~

## Example result

The example produces a two-dose design with a total target enrollment of 143
participants. The reported operating characteristics are exact for the
specified planning scenarios.

### Design performance

| Quantity | Target | Attained |
|:--|--:|--:|
| Strong familywise error rate (FWER) | 0.050 | 0.047051 |
| Minimum generalized power I (G1) | 0.800 | 0.807474 |
| Mean generalized power I (G1) | -- | 0.810183 |
| Minimum generalized power II (G2) | -- | 0.807990 |

The attained FWER is below its prespecified limit, while the minimum G1 is
above its target across the planning scenarios.

### Protocol-ready decision rules

| Dose | Target enrollment | Maximum toxicities | Minimum responses | Local worst-case error |
|:--|--:|--:|--:|--:|
| Low | 69 | 13 | 25 | 0.025624 |
| High | 74 | 14 | 27 | 0.021991 |

A dose is selected only when both conditions are satisfied: its observed
toxicity count is no greater than the maximum toxicity threshold, and its
observed efficacy count is no smaller than the minimum response threshold.
For example, the Low dose is selected when `X_T <= 13` and `X_E >= 25`.

### Scenario-specific power

| Planning scenario | G1 | G2 |
|:--|--:|--:|
| `dose_1_active` | 0.807474 | 0.807990 |
| `dose_2_active` | 0.812893 | 0.812893 |

G1 is the probability of selecting at least one admissible dose while
selecting no inadmissible dose. G2 is the probability of selecting at least
one admissible dose, regardless of whether an inadmissible dose is also
selected.

### Applying the design to observed counts

For `toxicity = c(5, 8)` and `efficacy = c(28, 31)`, the final result is:

| Dose | Enrolled | Toxicities | Responses | Toxicity pass | Efficacy pass | Selected |
|:--|--:|--:|--:|:--:|:--:|:--:|
| Low | 69 | 5 | 28 | Yes | Yes | Yes |
| High | 74 | 8 | 31 | Yes | Yes | Yes |

Both doses are therefore selected. The values supplied to `roed_select()` are
observed participant counts, not rates or percentages. Search time is not
reported here because it depends on the user's hardware.

The design search is exact and deterministic.

## User-facing functions

* roed_scenarios() defines clinically plausible planning or sensitivity
  scenarios.
* roed_design() constructs the protocol-ready exact design.
* roed_oc() evaluates exact operating characteristics without re-optimizing.
* roed_protocol() returns a decision table for the protocol or analysis plan.
* roed_select() applies the final integer rules to completed trial counts.
* roed_resume() continues an interrupted search from an explicit checkpoint.
