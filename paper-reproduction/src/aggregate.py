"""Build analysis-ready outputs for Studies 1--4 from completed atomic tasks."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from statistics import mean, median

import numpy as np

from .exact_adapter import load_engine


METHODS = ("MERIT", "Tabata", "ROED-EA", "ROED-B", "ROED")


def _checkpoints(root: Path) -> Path:
    return root / "_checkpoints"


def _results(root: Path) -> Path:
    return root / "results"


def _read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("status") != "complete":
        raise RuntimeError(f"Incomplete result: {path}")
    return value


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _exact_files(root: Path) -> list[Path]:
    files = sorted((_checkpoints(root) / "study1" / "exact").glob("J*.json"))
    if len(files) != 24:
        raise RuntimeError(f"Expected 24 formal exact results; found {len(files)}")
    return files


def _merit_files(root: Path) -> list[Path]:
    files = sorted((_checkpoints(root) / "study1" / "merit").glob("J*.json"))
    if len(files) != 24:
        raise RuntimeError(f"Expected 24 local MERIT results; found {len(files)}")
    return files


def build_study1(root: Path) -> list[dict]:
    rows = []
    merit_oc_rows = []
    for path in _exact_files(root):
        envelope = _read(path)
        item = envelope["payload"]
        st = item["setting"]
        for method, design in item["designs"].items():
            rows.append({
                "J": st["J"], "alpha": st["alpha"], "target": st["target"],
                "phi_E1": st["phi_E1"], "phi_E0": st["phi_E0"],
                "method": method, "source": "local exact engine",
                "total_n": design["total_n"], "fwer": design["fwer"],
                "min_g1": design["min_g1"], "mean_g1": design["mean_g1"],
                "min_g2": design["min_g2"],
            })
    for path in _merit_files(root):
        item = _read(path)
        st, design = item["input"], item["design"]
        oc = item.get("operating_characteristics") or {}
        oc_values = list(oc.values())
        finite_fwer = max((entry["max_finite_null_type1"] for entry in oc_values), default="")
        min_g1 = min((entry["min_power"] for entry in oc_values), default="")
        all_power = [row["value"] for entry in oc_values for row in entry["rows"]
                     if row["kind"] == "alternative"]
        rows.append({
            "J": st["J"], "alpha": st["alpha"], "target": st["target"],
            "phi_E1": st["phi_E1"], "phi_E0": st["phi_E0"],
            "method": "MERIT",
            "source": f"{item['source']['application']} {item['source']['version']}",
            "total_n": design["total_n"], "fwer": finite_fwer,
            "min_g1": min_g1,
            "mean_g1": mean(all_power) if all_power else "", "min_g2": "",
        })
        for oc_item in oc_values:
            for result in oc_item["rows"]:
                merit_oc_rows.append({
                    "J": st["J"], "alpha": st["alpha"], "target": st["target"],
                    "phi_E1": st["phi_E1"], "phi_E0": st["phi_E0"],
                    "rho": oc_item["rho"], "simulations": oc_item["simulations"],
                    "seed": oc_item["seed"], "scenario": result["scenario"],
                    "kind": result["kind"], "parameters": ";".join(map(str, result["parameters"])),
                    "metric": result["metric"], "value": result["value"],
                    "average_sample_size": result["average_sample_size"],
                    "engine_version": item["source"]["version"],
                })

    # Replace native-model MERIT OC summaries with the common-model
    # evaluation of the same frozen design rules.  The detailed native runs
    # remain available as provenance in merit_operating_characteristics.csv.
    common_path = _checkpoints(root) / "study1" / "merit_common_bernoulli_oc.csv"
    if common_path.exists():
        with common_path.open(encoding="utf-8", newline="") as stream:
            common_rows = list(csv.DictReader(stream))
        grouped = {}
        for result in common_rows:
            key = (
                int(result["J"]), float(result["alpha"]), float(result["target"]),
                float(result["phi_E1"]), float(result["phi_E0"]),
            )
            grouped.setdefault(key, []).append(result)
        for row in rows:
            if row["method"] != "MERIT":
                continue
            key = (row["J"], row["alpha"], row["target"], row["phi_E1"], row["phi_E0"])
            subset = grouped[key]
            row["fwer"] = max(float(x["fwer"]) for x in subset if x["kind"] == "null")
            g1 = [float(x["g1"]) for x in subset if x["kind"] == "alternative"]
            g2 = [float(x["g2"]) for x in subset if x["kind"] == "alternative"]
            row["min_g1"], row["mean_g1"], row["min_g2"] = min(g1), mean(g1), min(g2)
            row["source"] += "; common Bernoulli OC, 100000 trials/scenario"
    rows.sort(key=lambda r: (r["phi_E1"], r["J"], r["alpha"], r["target"], METHODS.index(r["method"])))
    _write_csv(_results(root) / "study1" / "study1_all_designs.csv", rows)
    _write_csv(_results(root) / "study1" / "merit_operating_characteristics.csv", merit_oc_rows)
    if common_path.exists():
        target = _results(root) / "study1" / "merit_common_bernoulli_oc.csv"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(common_path.read_bytes())
    return rows


def build_study2(root: Path, study1_rows: list[dict]) -> None:
    lookups = {}
    for row in study1_rows:
        key = (row["J"], row["alpha"], row["target"], row["phi_E1"], row["phi_E0"])
        lookups.setdefault(key, {})[row["method"]] = row
    rows = []
    for key, values in sorted(lookups.items(), key=lambda x: x[0]):
        tabata_n = float(values["Tabata"]["total_n"])
        for method in METHODS:
            row = values[method]
            rows.append({
                "J": row["J"], "alpha": row["alpha"], "target": row["target"],
                "phi_E1": row["phi_E1"], "phi_E0": row["phi_E0"],
                "method": method, "total_n": row["total_n"],
                "sample_ratio_to_tabata": float(row["total_n"]) / tabata_n,
                "power_excess": "" if row["min_g1"] == "" else float(row["min_g1"]) - float(row["target"]),
                "min_g1": row["min_g1"], "mean_g1": row["mean_g1"], "min_g2": row["min_g2"],
                "source": row["source"],
            })
    _write_csv(_results(root) / "study2" / "study2_power_sample_metrics.csv", rows)


def build_study3(root: Path) -> None:
    detail = []
    for path in _exact_files(root):
        item = _read(path)["payload"]
        st, d = item["setting"], item["designs"]
        nt, nea, nb, nr = (d[m]["total_n"] for m in ("Tabata", "ROED-EA", "ROED-B", "ROED"))
        dp, da, df = nt - nea, nt - nb, nt - nr
        ns = np.asarray([r["n"] for r in d["ROED"]["rules"]], dtype=float)
        detail.append({
            "J": st["J"], "alpha": st["alpha"], "target": st["target"],
            "phi_E1": st["phi_E1"], "phi_E0": st["phi_E0"],
            "delta_prod": dp, "delta_alloc": da,
            "delta_int": df - dp - da, "delta_full": df,
            "u_alpha_tabata": d["Tabata"]["fwer"] / st["alpha"],
            "u_alpha_roed": d["ROED"]["fwer"] / st["alpha"],
            "cv_n_roed": float(ns.std(ddof=0) / ns.mean()),
        })
    detail.sort(key=lambda r: (r["phi_E1"], r["J"], r["alpha"], r["target"]))
    _write_csv(_results(root) / "study3" / "study3_ablation_by_configuration.csv", detail)
    aggregate = []
    for e1, e0 in ((0.30, 0.10), (0.60, 0.40)):
        for j in (2, 3, 4):
            subset = [r for r in detail if r["phi_E1"] == e1 and r["J"] == j]
            aggregate.append({"phi_E1": e1, "phi_E0": e0, "J": j, **{
                key: mean(r[key] for r in subset) for key in (
                    "delta_prod", "delta_alloc", "delta_int", "delta_full",
                    "u_alpha_tabata", "u_alpha_roed", "cv_n_roed",
                )
            }})
    _write_csv(_results(root) / "study3" / "study3_ablation_aggregated.csv", aggregate)


def _fixed_design_metrics(item: dict, rules: list[dict], phi_t1: float, phi_e1: float,
                          rhos: tuple[float, ...]) -> dict:
    eng = load_engine()
    st = item["setting"]
    scenarios = eng.make_scenarios(st["J"], phi_t1, st["phi_T0"], phi_e1, st["phi_E0"], rhos=rhos)
    nmax = max(r["n"] for r in rules)
    cache = eng.ProbabilityCache(scenarios, n_max=nmax)
    candidates = []
    for j, rule in enumerate(rules):
        pi = np.asarray([
            cache.selection(rule["n"], rule["mT"], rule["mE"], sc.p_t[j], sc.p_e[j], sc.rho)
            for sc in scenarios
        ])
        a = float(rule["a"])
        candidates.append(eng.Candidate(rule["n"], rule["mT"], rule["mE"],
                                        rule["aT"], rule["aE"], a, -math.log1p(-a), pi))
    return eng.design_metrics(candidates, scenarios)


def build_study4(root: Path, config: dict) -> None:
    effects, correlations, performance = [], [], []
    contractions = config["study4"]["effect_contractions"]
    evaluation_rho = config["study4"]["evaluation_rho"]
    for path in _exact_files(root):
        envelope = _read(path)
        item = envelope["payload"]
        st, roed = item["setting"], item["designs"]["ROED"]
        base = {"J": st["J"], "alpha": st["alpha"], "target": st["target"],
                "phi_E1": st["phi_E1"], "phi_E0": st["phi_E0"]}
        for dt in contractions:
            for de in contractions:
                met = _fixed_design_metrics(item, roed["rules"], st["phi_T1"] + dt,
                                            st["phi_E1"] - de, (0.0, 0.30))
                effects.append({**base, "delta_T": dt, "delta_E": de,
                                "min_g1": met["min_g1"], "mean_g1": met["mean_g1"],
                                "min_g2": met["min_g2"], "fwer": roed["fwer"]})
        for rho in evaluation_rho:
            met = _fixed_design_metrics(item, roed["rules"], st["phi_T1"], st["phi_E1"], (rho,))
            correlations.append({**base, "rho": rho, "min_g1": met["min_g1"],
                                 "mean_g1": met["mean_g1"], "min_g2": met["min_g2"],
                                 "fwer": roed["fwer"]})
        performance.append({**base, "repetition": 0,
                            "elapsed_seconds": envelope["elapsed_seconds"],
                            "retained_candidates": sum(roed["candidate_counts"])})

    benchmark_files = sorted((_checkpoints(root) / "study4" / "benchmark").glob("benchmark_*.json"))
    if benchmark_files:
        performance = []
        for path in benchmark_files:
            envelope = _read(path)
            item = envelope["payload"]
            st, roed = item["setting"], item["designs"]["ROED"]
            repetition = int(path.stem.rsplit("_r", 1)[1])
            performance.append({
                "J": st["J"], "alpha": st["alpha"], "target": st["target"],
                "phi_E1": st["phi_E1"], "phi_E0": st["phi_E0"],
                "repetition": repetition, "elapsed_seconds": envelope["elapsed_seconds"],
                "retained_candidates": sum(roed["candidate_counts"]),
            })
    _write_csv(_results(root) / "study4" / "effect_contraction.csv", effects)
    _write_csv(_results(root) / "study4" / "correlation_sensitivity.csv", correlations)
    _write_csv(_results(root) / "study4" / "exact_search_performance.csv", performance)

    roed_s_files = sorted((_checkpoints(root) / "study4" / "roed_s").glob("J*.json"))
    if len(roed_s_files) != 24:
        raise RuntimeError(f"Expected 24 ROED-S results; found {len(roed_s_files)}")
    roed_s_rows = []
    for path in roed_s_files:
        item = _read(path)["payload"]
        st = item["setting"]
        exact_name = path.name
        exact_item = _read(_checkpoints(root) / "study1" / "exact" / exact_name)["payload"]
        roed = exact_item["designs"]["ROED"]
        roed_s_rows.append({
            "J": st["J"], "alpha": st["alpha"], "target": st["target"],
            "phi_E1": st["phi_E1"], "phi_E0": st["phi_E0"],
            "nominal_scenario": item["nominal_scenario"],
            "roed_s_total_n": item["total_n"], "roed_total_n": roed["total_n"],
            "roed_s_nominal_g1": item["nominal_min_g1"],
            "roed_s_full_min_g1": item["full_min_g1"],
            "roed_full_min_g1": roed["min_g1"],
            "roed_s_power_loss": item["nominal_min_g1"] - item["full_min_g1"],
            "roed_s_fwer": item["fwer"], "roed_fwer": roed["fwer"],
            "roed_s_elapsed_seconds": item["elapsed_seconds"],
        })
    roed_s_rows.sort(key=lambda r: (r["phi_E1"], r["J"], r["alpha"], r["target"]))
    _write_csv(_results(root) / "study4" / "roed_s_comparison_by_configuration.csv", roed_s_rows)
    roed_s_summary = []
    for e1, e0 in ((0.30, 0.10), (0.60, 0.40)):
        for j in (2, 3, 4):
            subset = [r for r in roed_s_rows if r["phi_E1"] == e1 and r["J"] == j]
            roed_s_summary.append({
                "phi_E1": e1, "phi_E0": e0, "J": j,
                "mean_roed_s_total_n": mean(r["roed_s_total_n"] for r in subset),
                "mean_roed_total_n": mean(r["roed_total_n"] for r in subset),
                "mean_roed_s_nominal_g1": mean(r["roed_s_nominal_g1"] for r in subset),
                "mean_roed_s_full_min_g1": mean(r["roed_s_full_min_g1"] for r in subset),
                "mean_roed_full_min_g1": mean(r["roed_full_min_g1"] for r in subset),
                "mean_roed_s_power_loss": mean(r["roed_s_power_loss"] for r in subset),
            })
    _write_csv(_results(root) / "study4" / "roed_s_comparison_summary.csv", roed_s_summary)


def build_all(root: Path, config: dict) -> None:
    rows = build_study1(root)
    build_study2(root, rows)
    build_study3(root)
    build_study4(root, config)
