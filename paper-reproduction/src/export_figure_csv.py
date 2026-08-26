"""Export the numerical inputs used by manuscript figures (CSV only)."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from statistics import mean, median

from scipy.stats import binom


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figure-data"
METHODS = ("MERIT", "Tabata", "ROED-EA", "ROED-B", "ROED")
METHOD_X = {method: index + 1 for index, method in enumerate(METHODS)}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def number(value: str) -> float:
    return float(value) if value not in ("", None) else math.nan


def study12() -> None:
    rows = read_csv(RESULTS / "study1" / "study1_all_designs.csv")
    centers = {2: -.13, 3: 0, 4: .13}
    settings = ((e, a, p) for e in (.30, .60) for a in (.05, .10) for p in (.80, .90))
    offsets = {key: -.035 + index * .01 for index, key in enumerate(settings)}
    tabata = {
        (int(r["J"]), number(r["alpha"]), number(r["target"]), number(r["phi_E1"])): number(r["total_n"])
        for r in rows if r["method"] == "Tabata"
    }
    first, second = [], []
    for row in rows:
        j, method = int(row["J"]), row["method"]
        alpha, target, e1 = number(row["alpha"]), number(row["target"]), number(row["phi_E1"])
        x = METHOD_X[method] + centers[j] + offsets[(e1, alpha, target)]
        first.append({"J": j, "method": method, "xplot": x,
                      "fwer_ratio": number(row["fwer"]) / alpha})
        second.append({
            "J": j, "method": method, "xplot": x,
            "sample_ratio": number(row["total_n"]) / tabata[(j, alpha, target, e1)],
            "power_excess": number(row["min_g1"]) - target,
            "mean_g1": number(row["mean_g1"]), "min_g2": number(row["min_g2"]),
        })
    for j in (2, 3, 4):
        write_csv(FIGURES / "study1" / f"J{j}.csv", [r for r in first if r["J"] == j])
        write_csv(FIGURES / "study2" / f"J{j}.csv", [r for r in second if r["J"] == j])
    for folder, data, metrics in (
        ("study1", first, ("fwer_ratio",)),
        ("study2", second, ("sample_ratio", "power_excess", "mean_g1", "min_g2")),
    ):
        for metric in metrics:
            summary = []
            for method in METHODS:
                values = [r[metric] for r in data if r["method"] == method and not math.isnan(r[metric])]
                y = median(values)
                summary.extend(({"x": METHOD_X[method] - .22, "y": y},
                                {"x": METHOD_X[method] + .22, "y": y},
                                {"x": "nan", "y": "nan"}))
            write_csv(FIGURES / folder / f"median_{metric}.csv", summary)

    common = read_csv(RESULTS / "study1" / "merit_common_bernoulli_oc.csv")
    grouped = {}
    for row in common:
        if row["kind"] != "null":
            continue
        key = (int(row["J"]), number(row["alpha"]), number(row["target"]),
               number(row["phi_E1"]), number(row["rho"]))
        grouped[key] = max(grouped.get(key, 0), number(row["fwer"]) / number(row["alpha"]))
    rho_rows = []
    for (j, alpha, target, e1, rho), ratio in grouped.items():
        x = (1 if rho == 0 else 2) + centers[j] + offsets[(e1, alpha, target)]
        rho_rows.append({"J": j, "rho": rho, "xplot": x, "fwer_ratio": ratio})
    for j in (2, 3, 4):
        write_csv(FIGURES / "study1" / f"merit_rho_J{j}.csv",
                  [r for r in rho_rows if r["J"] == j])


def study3() -> None:
    path = ROOT / "_checkpoints" / "study1" / "exact" / "J4_a0p05_p0p80_e0p60_0p40.json"
    item = json.loads(path.read_text(encoding="utf-8"))["payload"]
    methods = METHODS[1:]
    offsets = {"Tabata": -.27, "ROED-EA": -.09, "ROED-B": .09, "ROED": .27}
    bars = []
    for method in methods:
        rules = item["designs"][method]["rules"]
        write_csv(FIGURES / "study3" / f"{method.lower().replace('-', '_')}.csv",
                  [{"dose": i, "a": r["a"], "n": r["n"]} for i, r in enumerate(rules, 1)])
        bars.extend({"method": method, "dose": i, "xplot": i + offsets[method], "n": r["n"]}
                    for i, r in enumerate(rules, 1))
        powers = item["designs"][method]["g1"]
        write_csv(FIGURES / "study3" / f"power_{method.lower().replace('-', '_')}.csv",
                  [{"index": i, "scenario": name, "g1": powers[name]}
                   for i, name in enumerate(item["scenarios"], 1)])
    write_csv(FIGURES / "study3" / "sample_size_bars.csv", bars)
    setting, rules = item["setting"], item["designs"]["ROED"]["rules"]
    tails = []
    for group, rule in (("Dose 1", rules[0]), ("Doses 2--4", rules[1])):
        for h in (-1, 0, 1):
            tails.append({"group": group, "endpoint": "Toxicity", "h": h,
                          "tail": binom.cdf(rule["mT"] + h, rule["n"], setting["phi_T0"])})
            tails.append({"group": group, "endpoint": "Efficacy", "h": h,
                          "tail": binom.sf(rule["mE"] + h - 1, rule["n"], setting["phi_E0"])})
    write_csv(FIGURES / "study3" / "threshold_neighborhood.csv", tails)
    for group_id, group in (("d1", "Dose 1"), ("d234", "Doses 2--4")):
        for endpoint_id, endpoint in (("t", "Toxicity"), ("e", "Efficacy")):
            write_csv(FIGURES / "study3" / f"tail_{group_id}_{endpoint_id}.csv",
                      [r for r in tails if r["group"] == group and r["endpoint"] == endpoint])


def study4() -> None:
    effect = read_csv(RESULTS / "study4" / "effect_contraction.csv")
    corr = read_csv(RESULTS / "study4" / "correlation_sensitivity.csv")
    roed_s = read_csv(RESULTS / "study4" / "roed_s_comparison_summary.csv")
    bench = read_csv(RESULTS / "study4" / "exact_search_performance.csv")
    heat = []
    for dt in (0, .025, .05):
        for de in (0, .025, .05):
            subset = [r for r in effect if number(r["delta_T"]) == dt and number(r["delta_E"]) == de]
            heat.append({"delta_T": dt, "delta_E": de,
                         "worst_g1": min(number(r["min_g1"]) for r in subset),
                         "mean_g1": mean(number(r["min_g1"]) for r in subset)})
    write_csv(FIGURES / "study4" / "effect_grid.csv", heat)
    correlation = []
    for j in (2, 3, 4):
        for rho in (-.1, 0, .15, .3):
            subset = [r for r in corr if int(r["J"]) == j and number(r["rho"]) == rho]
            correlation.append({"J": j, "rho": rho,
                                "mean_min_g1": mean(number(r["min_g1"]) for r in subset),
                                "worst_min_g1": min(number(r["min_g1"]) for r in subset)})
    write_csv(FIGURES / "study4" / "correlation.csv", correlation)
    comparison = []
    for j in (2, 3, 4):
        subset = [r for r in roed_s if int(r["J"]) == j]
        comparison.append({"J": j,
            "roeds_n": mean(number(r["mean_roed_s_total_n"]) for r in subset),
            "roed_n": mean(number(r["mean_roed_total_n"]) for r in subset),
            "roeds_nominal": mean(number(r["mean_roed_s_nominal_g1"]) for r in subset),
            "roeds_full": mean(number(r["mean_roed_s_full_min_g1"]) for r in subset),
            "roed_full": mean(number(r["mean_roed_full_min_g1"]) for r in subset)})
    write_csv(FIGURES / "study4" / "roed_s.csv", comparison)
    timing = []
    for j in (2, 3, 4):
        subset = [r for r in bench if int(r["J"]) == j]
        times = sorted(number(r["elapsed_seconds"]) for r in subset)
        retained = [number(r["retained_candidates"]) for r in subset]
        timing.append({"J": j, "runs": len(subset), "median_retained": median(retained),
                       "median_time": median(times), "q1_time": median(times[:len(times)//2]),
                       "q3_time": median(times[(len(times)+1)//2:]), "max_time": max(times)})
    write_csv(FIGURES / "study4" / "exact_time.csv", timing)


def main() -> None:
    study12()
    study3()
    study4()


if __name__ == "__main__":
    main()
