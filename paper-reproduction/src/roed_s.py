"""Nominal-scenario ROED-S comparator for Numerical Study 4."""

from __future__ import annotations

import math
import time

import numpy as np

from .exact_adapter import load_engine


def _nominal_window(engine, j_count: int, phi_t1: float, phi_t0: float,
                    phi_e1: float, phi_e0: float, rho: float):
    # Manuscript convention: k0 = ceil(J/2), with dose labels starting at 1.
    k0 = math.ceil(j_count / 2)
    pt, pe = [], []
    for dose in range(1, j_count + 1):
        if dose < k0:
            pt.append(phi_t1); pe.append(phi_e0)
        elif dose == k0:
            pt.append(phi_t1); pe.append(phi_e1)
        else:
            pt.append(phi_t0); pe.append(phi_e1)
    admissible = tuple(t < phi_t0 and e > phi_e0 for t, e in zip(pt, pe))
    return k0, engine.Scenario(
        name=f"W{k0}_rho{rho:.2f}", rho=float(rho),
        p_t=tuple(pt), p_e=tuple(pe), admissible=admissible,
    )


def _evaluate(engine, rules: list[dict], scenarios, phi_t0: float, phi_e0: float) -> dict:
    cache = engine.ProbabilityCache(scenarios, n_max=max(rule["n"] for rule in rules))
    candidates = []
    for dose, rule in enumerate(rules):
        pi = np.asarray([
            cache.selection(rule["n"], rule["mT"], rule["mE"],
                            scenario.p_t[dose], scenario.p_e[dose], scenario.rho)
            for scenario in scenarios
        ])
        a = float(rule["a"])
        candidates.append(engine.Candidate(
            int(rule["n"]), int(rule["mT"]), int(rule["mE"]),
            float(rule["aT"]), float(rule["aE"]), a, -math.log1p(-a), pi,
        ))
    return engine.design_metrics(candidates, scenarios)


def solve_roed_s(spec: dict, nominal_rho: float = 0.15,
                 evaluation_rhos=(0.0, 0.30)) -> dict:
    started = time.perf_counter()
    engine = load_engine()
    j_count = int(spec["J"])
    k0, nominal = _nominal_window(
        engine, j_count, spec["phi_T1"], spec["phi_T0"],
        spec["phi_E1"], spec["phi_E0"], nominal_rho,
    )
    nominal_scenarios = [nominal]
    cache = engine.ProbabilityCache(nominal_scenarios, n_max=int(spec["n_max"]))
    raw_counts, pools = [], []
    for dose in range(j_count):
        raw = engine.generate_candidates(
            dose, nominal_scenarios, cache, spec["phi_T0"], spec["phi_E0"],
            spec["alpha"], spec["target"], spec["n_min"], spec["n_max"], None,
        )
        raw_counts.append(len(raw))
        pools.append(engine.prune_dominated(raw, dose, nominal_scenarios))

    design, nominal_metrics = engine.search_design(
        pools, nominal_scenarios, spec["alpha"], spec["target"],
        equal_n=False, common_thresholds=False, error_mode="product",
    )
    rules = [{
        "n": candidate.n, "mT": candidate.m_t, "mE": candidate.m_e,
        "aT": candidate.a_t, "aE": candidate.a_e, "a": candidate.a,
    } for candidate in design]

    full_scenarios = engine.make_scenarios(
        j_count, spec["phi_T1"], spec["phi_T0"], spec["phi_E1"], spec["phi_E0"],
        rhos=tuple(evaluation_rhos),
    )
    full_metrics = _evaluate(engine, rules, full_scenarios, spec["phi_T0"], spec["phi_E0"])
    return {
        "setting": {k: spec[k] for k in (
            "J", "alpha", "target", "phi_T1", "phi_T0", "phi_E1", "phi_E0",
            "n_min", "n_max",
        )},
        "method": "ROED-S",
        "nominal_scenario": nominal.name,
        "nominal_rho": nominal_rho,
        "evaluation_rhos": list(evaluation_rhos),
        "rules": rules,
        "total_n": nominal_metrics["total_n"],
        "fwer": nominal_metrics["fwer"],
        "nominal_min_g1": nominal_metrics["min_g1"],
        "nominal_mean_g1": nominal_metrics["mean_g1"],
        "nominal_min_g2": nominal_metrics["min_g2"],
        "full_min_g1": full_metrics["min_g1"],
        "full_mean_g1": full_metrics["mean_g1"],
        "full_min_g2": full_metrics["min_g2"],
        "full_g1": dict(zip((scenario.name for scenario in full_scenarios), full_metrics["g1"])),
        "raw_candidate_counts": raw_counts,
        "retained_candidate_counts": [len(pool) for pool in pools],
        "elapsed_seconds": time.perf_counter() - started,
    }
