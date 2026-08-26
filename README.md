# ROED

ROED implements Robustly Optimal Exact Designs for randomized
dose-optimization trials with co-primary binary toxicity and efficacy
endpoints.

## Installation

~~~r
# install.packages("remotes")
remotes::install_github("haohaostats/ROED")
~~~

## Worked example: randomized sotorasib dose comparison

This example reconstructs the application presented in the accompanying
manuscript. The published study compared sotorasib 240 mg and 960 mg in
previously treated advanced non-small-cell lung cancer with a
*KRAS G12C* mutation. The public aggregate data are described by
[Hochmair et al. (2024)](https://doi.org/10.1016/j.ejca.2024.114204).

The toxicity null boundary is 0.40, the efficacy null boundary is 0.10,
the target strong FWER is 0.10, and the target minimum G1 is 0.80. Three
clinical configurations are evaluated at correlations 0 and 0.30, producing
six planning scenarios.

### Define the planning scenarios

~~~r
library(ROED)

sotorasib_scenarios <- roed_scenarios(
  toxicity = rbind(
    c(0.20, 0.40), # S_L, rho = 0
    c(0.20, 0.20), # S_H, rho = 0
    c(0.20, 0.20), # S_P, rho = 0
    c(0.20, 0.40), # S_L, rho = 0.30
    c(0.20, 0.20), # S_H, rho = 0.30
    c(0.20, 0.20)  # S_P, rho = 0.30
  ),
  efficacy = rbind(
    c(0.30, 0.30),
    c(0.10, 0.30),
    c(0.30, 0.30),
    c(0.30, 0.30),
    c(0.10, 0.30),
    c(0.30, 0.30)
  ),
  correlation = c(0, 0, 0, 0.30, 0.30, 0.30),
  scenario_names = c(
    "S_L_rho_0", "S_H_rho_0", "S_P_rho_0",
    "S_L_rho_0.30", "S_H_rho_0.30", "S_P_rho_0.30"
  ),
  dose_names = c("240 mg", "960 mg")
)
~~~

Here `S_L` represents a configuration in which only 240 mg is admissible,
`S_H` one in which only 960 mg is admissible, and `S_P` a plateau
configuration in which both doses are admissible.

### Construct the prospective design

~~~r
fit <- roed_design(
  scenarios = sotorasib_scenarios,
  tox_null = 0.40,
  eff_null = 0.10,
  alpha = 0.10,
  target_power = 0.80,
  n_min = 10,
  n_max = 150
)

summary(fit)
roed_protocol(fit)
~~~

The resulting protocol-ready design is:

| Dose | Target enrollment | Maximum toxicities | Minimum responses | Local worst-case error |
|:--|--:|--:|--:|--:|
| 240 mg | 41 | 11 | 8 | 0.056659 |
| 960 mg | 45 | 12 | 9 | 0.044631 |

| Quantity | Target | Attained |
|:--|--:|--:|
| Total enrollment | -- | 86 |
| Strong familywise error rate (FWER) | 0.100 | 0.098761 |
| Minimum generalized power I (G1) | 0.800 | 0.820166 |
| Mean generalized power I (G1) | -- | 0.874144 |
| Minimum generalized power II (G2) | -- | 0.853046 |

The design therefore meets both prespecified requirements: its exact strong
FWER is below 0.10 and its minimum G1 is above 0.80 across all six planning
scenarios.

### Replay the published aggregate data

The published trial had 104 treated participants in each dose group, more
than the prospective ROED quotas. The correct replay therefore fixes both
sample sizes at 104 and reoptimizes the integer decision boundaries before
applying the observed counts.

~~~r
replay_fit <- roed_design(
  scenarios = sotorasib_scenarios,
  tox_null = 0.40,
  eff_null = 0.10,
  alpha = 0.10,
  target_power = 0.80,
  n_min = 104,
  n_max = 104
)

roed_protocol(replay_fit)

replay_result <- roed_select(
  replay_fit,
  toxicity = c(20, 37),
  efficacy = c(26, 34)
)

replay_result
~~~

The fixed-size rules and observed decisions are:

| Dose | Maximum toxicities | Minimum responses | Observed toxicities | Observed responses | Selected |
|:--|--:|--:|--:|--:|:--:|
| 240 mg | 34 | 20 | 20 | 26 | Yes |
| 960 mg | 30 | 19 | 37 | 34 | No |

The 240-mg dose passes both requirements. The 960-mg dose passes the efficacy
requirement but fails the toxicity requirement because `37 > 30`.
Consequently, **240 mg is the only dose selected in this illustration**.

This is an illustrative statistical reconstruction based on prespecified
binary endpoints and public aggregate counts. It is not a claim of clinical
superiority and does not replace a complete benefit-risk assessment.

## Reading the operating characteristics

* Strong FWER is the worst-case probability of selecting at least one
  inadmissible dose from the candidate set.
* G1 is the probability of selecting at least one admissible dose while
  selecting no inadmissible dose.
* G2 is the probability of selecting at least one admissible dose, regardless
  of whether an inadmissible dose is also selected.
* Values supplied to `roed_select()` are participant counts, not rates or
  percentages.

The design search is exact and deterministic. Runtime depends on the number
of doses, planning scenarios, search bounds, and the user's hardware.

## User-facing functions

* `roed_scenarios()` defines clinically plausible planning or sensitivity
  scenarios.
* `roed_design()` constructs the protocol-ready exact design.
* `roed_oc()` evaluates exact operating characteristics without re-optimizing.
* `roed_protocol()` returns a decision table for the protocol or analysis plan.
* `roed_select()` applies the final integer rules to completed trial counts.
* `roed_resume()` continues an interrupted search from an explicit checkpoint.
