"""Exact design engine for ROED numerical study 1.

This module implements the finite probability calculations and the four
exact comparator designs described in manuscript/manuscript.tex:

* Tabata: equal n, common thresholds, Bonferroni local errors;
* ROED-EA: equal n, dose-specific thresholds, exact product FWER;
* ROED-B: unequal n and thresholds, Bonferroni local errors;
* ROED: unequal n and thresholds, exact product FWER.

The implementation is deliberately self-contained and deterministic.  It
uses scipy only for stable binomial tails; all bivariate Bernoulli count
probabilities are obtained by the finite recursion in the manuscript.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix
from scipy.stats import binom


@dataclass(frozen=True)
class Scenario:
    name: str
    rho: float
    p_t: tuple[float, ...]
    p_e: tuple[float, ...]
    admissible: tuple[bool, ...]


@dataclass
class Candidate:
    n: int
    m_t: int
    m_e: int
    a_t: float
    a_e: float
    a: float
    r: float
    pi: np.ndarray

    def key(self) -> tuple[int, int, int]:
        return (self.n, self.m_t, self.m_e)


def make_scenarios(
    j_count: int,
    phi_t1: float,
    phi_t0: float,
    phi_e1: float,
    phi_e0: float,
    rhos: Sequence[float] = (0.0, 0.30),
) -> list[Scenario]:
    patterns: list[tuple[str, list[float], list[float]]] = []
    for k in range(j_count):
        pt, pe = [], []
        for j in range(j_count):
            if j < k:
                pt.append(phi_t1); pe.append(phi_e0)
            elif j == k:
                pt.append(phi_t1); pe.append(phi_e1)
            else:
                pt.append(phi_t0); pe.append(phi_e1)
        patterns.append((f"W{k + 1}", pt, pe))
    for k in range(j_count - 1):
        pt, pe = [], []
        for j in range(j_count):
            if j < k:
                pt.append(phi_t1); pe.append(phi_e0)
            elif j in (k, k + 1):
                pt.append(phi_t1); pe.append(phi_e1)
            else:
                pt.append(phi_t0); pe.append(phi_e1)
        patterns.append((f"P{k + 1}", pt, pe))

    out: list[Scenario] = []
    for name, pt, pe in patterns:
        adm = tuple(t < phi_t0 and e > phi_e0 for t, e in zip(pt, pe))
        for rho in rhos:
            out.append(
                Scenario(
                    name=f"{name}_rho{rho:.2f}",
                    rho=float(rho),
                    p_t=tuple(pt),
                    p_e=tuple(pe),
                    admissible=adm,
                )
            )
    return out


def bernoulli_cells(p_t: float, p_e: float, rho: float) -> tuple[float, ...]:
    q = p_t * p_e + rho * math.sqrt(p_t * (1.0 - p_t) * p_e * (1.0 - p_e))
    lower = max(0.0, p_t + p_e - 1.0)
    upper = min(p_t, p_e)
    if q < lower - 1e-12 or q > upper + 1e-12:
        raise ValueError(f"Infeasible Bernoulli correlation: pT={p_t}, pE={p_e}, rho={rho}")
    q = min(max(q, lower), upper)
    return (1.0 - p_t - p_e + q, p_t - q, p_e - q, q)


def reverse_cumulative_surfaces(
    p_t: float, p_e: float, rho: float, n_max: int
) -> list[np.ndarray | None]:
    """Return Q[n][u,v] = P(X_T <= u, X_E >= v), n=0,...,n_max."""
    eta00, eta10, eta01, eta11 = bernoulli_cells(p_t, p_e, rho)
    dist = np.ones((1, 1), dtype=float)
    out: list[np.ndarray | None] = [np.ones((1, 1), dtype=float)] + [None] * n_max
    for n in range(1, n_max + 1):
        nxt = np.zeros((n + 1, n + 1), dtype=float)
        nxt[:n, :n] += eta00 * dist
        nxt[1:, :n] += eta10 * dist
        nxt[:n, 1:] += eta01 * dist
        nxt[1:, 1:] += eta11 * dist
        dist = nxt
        lower_t = np.cumsum(dist, axis=0)
        out[n] = np.flip(np.cumsum(np.flip(lower_t, axis=1), axis=1), axis=1)
    return out


class ProbabilityCache:
    def __init__(self, scenarios: Sequence[Scenario], n_max: int):
        states = sorted(
            {(sc.p_t[j], sc.p_e[j], sc.rho) for sc in scenarios for j in range(len(sc.p_t))}
        )
        self.surfaces = {
            state: reverse_cumulative_surfaces(*state, n_max=n_max) for state in states
        }

    def selection(self, n: int, m_t: int, m_e: int, p_t: float, p_e: float, rho: float) -> float:
        return float(self.surfaces[(p_t, p_e, rho)][n][m_t, m_e])


def design_metrics(candidates: Sequence[Candidate], scenarios: Sequence[Scenario]) -> dict:
    a = np.array([c.a for c in candidates])
    fwer = 1.0 - float(np.prod(1.0 - a))
    g1, g2 = [], []
    for s, sc in enumerate(scenarios):
        p = np.array([c.pi[s] for c in candidates])
        adm = np.array(sc.admissible, dtype=bool)
        no_false = float(np.prod(1.0 - p[~adm])) if np.any(~adm) else 1.0
        any_true = 1.0 - float(np.prod(1.0 - p[adm]))
        g1.append(no_false * any_true)
        g2.append(any_true)
    return {
        "total_n": int(sum(c.n for c in candidates)),
        "fwer": fwer,
        "min_g1": min(g1),
        "mean_g1": float(np.mean(g1)),
        "min_g2": min(g2),
        "g1": g1,
        "g2": g2,
    }


def lex_objective(candidates: Sequence[Candidate], scenarios: Sequence[Scenario]) -> tuple:
    met = design_metrics(candidates, scenarios)
    return (
        met["total_n"],
        -met["min_g1"],
        -met["mean_g1"],
        met["fwer"],
        tuple(c.n for c in candidates),
        tuple(c.m_t for c in candidates),
        tuple(-c.m_e for c in candidates),
    )


def generate_candidates(
    dose_index: int,
    scenarios: Sequence[Scenario],
    cache: ProbabilityCache,
    phi_t0: float,
    phi_e0: float,
    alpha_star: float,
    target: float,
    n_min: int,
    n_max: int,
    local_cap: float | None,
) -> list[Candidate]:
    """Generate scientifically viable dose-level rules.

    The sole-admissible W_j scenarios imply pi_j >= target for every
    association variant.  Enforcing that necessary condition removes many
    candidates without changing the feasible design set.
    """
    relevant_sole = [
        s for s, sc in enumerate(scenarios)
        if sum(sc.admissible) == 1 and sc.admissible[dose_index]
    ]
    cap = alpha_star if local_cap is None else local_cap
    out: list[Candidate] = []
    for n in range(n_min, n_max + 1):
        a_t = binom.cdf(np.arange(0, n), n, phi_t0)
        a_e = binom.sf(np.arange(1, n + 1) - 1, n, phi_e0)
        mt_ok = np.flatnonzero(a_t <= cap + 1e-15)
        me_ok = np.flatnonzero(a_e <= cap + 1e-15) + 1
        if len(mt_ok) == 0 or len(me_ok) == 0:
            continue
        for m_t in mt_ok:
            for m_e in me_ok:
                at = float(a_t[m_t]); ae = float(a_e[m_e - 1]); aa = max(at, ae)
                if aa > cap + 1e-15:
                    continue
                pi = np.array([
                    cache.selection(n, int(m_t), int(m_e), sc.p_t[dose_index], sc.p_e[dose_index], sc.rho)
                    for sc in scenarios
                ])
                if relevant_sole and min(pi[relevant_sole]) + 1e-13 < target:
                    continue
                out.append(Candidate(n, int(m_t), int(m_e), at, ae, aa, -math.log1p(-aa), pi))
    return out


def _dominance_vectors(cands: Sequence[Candidate], dose_index: int, scenarios: Sequence[Scenario]) -> np.ndarray:
    rows = []
    for c in cands:
        direction = [(-c.pi[s] if sc.admissible[dose_index] else c.pi[s]) for s, sc in enumerate(scenarios)]
        rows.append([float(c.n), c.r, *direction])
    return np.asarray(rows)


def prune_dominated(cands: Sequence[Candidate], dose_index: int, scenarios: Sequence[Scenario]) -> list[Candidate]:
    """Exact Pareto pruning using vectorized blocks.

    Candidate counts after the necessary single-window power filter are
    modest.  Sorting by n and error cost lets the skyline be maintained
    without quadratic storage.
    """
    if not cands:
        return []
    vec = _dominance_vectors(cands, dose_index, scenarios)
    order = np.lexsort((vec[:, 1], vec[:, 0]))
    front_idx: list[int] = []
    tol = 2e-14
    for idx in order:
        v = vec[idx]
        if front_idx:
            f = vec[front_idx]
            le = np.all(f <= v + tol, axis=1)
            strict = np.any(f < v - tol, axis=1)
            if np.any(le & strict):
                continue
            # A current rule can only dominate rules with the same n because
            # processing order guarantees all prior rules have n <= current n.
            same_n = np.abs(f[:, 0] - v[0]) <= tol
            dominated = same_n & np.all(v <= f + tol, axis=1) & np.any(v < f - tol, axis=1)
            if np.any(dominated):
                front_idx = [k for k, drop in zip(front_idx, dominated) if not drop]
        front_idx.append(int(idx))
    kept = [cands[i] for i in front_idx]
    kept.sort(key=lambda c: (c.n, c.r, c.m_t, -c.m_e))
    return kept


def _optimistic_g1(
    chosen: Sequence[Candidate],
    next_dose: int,
    pools: Sequence[Sequence[Candidate]],
    scenarios: Sequence[Scenario],
) -> float:
    min_bound = 1.0
    for s, sc in enumerate(scenarios):
        no_false = 1.0
        no_true = 1.0
        for j, c in enumerate(chosen):
            if sc.admissible[j]:
                no_true *= 1.0 - c.pi[s]
            else:
                no_false *= 1.0 - c.pi[s]
        for j in range(next_dose, len(pools)):
            vals = [c.pi[s] for c in pools[j]]
            if sc.admissible[j]:
                no_true *= 1.0 - max(vals)
            else:
                no_false *= 1.0 - min(vals)
        min_bound = min(min_bound, no_false * (1.0 - no_true))
    return min_bound


def milp_minimum_total(
    pools: Sequence[Sequence[Candidate]],
    scenarios: Sequence[Scenario],
    alpha_star: float,
    target: float,
    equal_n: bool,
    error_mode: str,
) -> int:
    """Obtain an exact lower bound (and normally the optimum) for total N.

    In every single-window scenario W_j there is exactly one admissible dose,
    so log(G1) is additive across selected dose-level rules.  These necessary
    constraints, the multiplicity constraint, and equal-n restrictions are a
    multiple-choice MILP.  Plateau scenarios are checked by the exhaustive
    search that follows; consequently this routine can never skip a feasible
    smaller total.
    """
    offsets = np.cumsum([0] + [len(p) for p in pools])
    nvar = int(offsets[-1])
    objective = np.zeros(nvar)
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    lb: list[float] = []
    ub: list[float] = []

    def add_row(coeffs: Iterable[tuple[int, float]], lower: float, upper: float) -> None:
        row = len(lb)
        for col, value in coeffs:
            rows.append(row); cols.append(col); vals.append(float(value))
        lb.append(lower); ub.append(upper)

    for j, pool in enumerate(pools):
        start = int(offsets[j])
        for k, cand in enumerate(pool):
            objective[start + k] = cand.n
        add_row(((start + k, 1.0) for k in range(len(pool))), 1.0, 1.0)

    budget = -math.log1p(-alpha_star) if error_mode == "product" else alpha_star
    add_row(
        (
            (int(offsets[j]) + k, cand.r if error_mode == "product" else cand.a)
            for j, pool in enumerate(pools) for k, cand in enumerate(pool)
        ),
        -np.inf,
        budget,
    )

    for s, sc in enumerate(scenarios):
        if sum(sc.admissible) != 1:
            continue
        coeffs = []
        for j, pool in enumerate(pools):
            for k, cand in enumerate(pool):
                p = min(max(float(cand.pi[s]), 1e-300), 1.0 - 1e-16)
                value = math.log(p) if sc.admissible[j] else math.log1p(-p)
                coeffs.append((int(offsets[j]) + k, value))
        add_row(coeffs, math.log(target), np.inf)

    if equal_n:
        for j in range(1, len(pools)):
            coeffs = []
            for k, cand in enumerate(pools[0]):
                coeffs.append((k, float(cand.n)))
            for k, cand in enumerate(pools[j]):
                coeffs.append((int(offsets[j]) + k, -float(cand.n)))
            add_row(coeffs, 0.0, 0.0)

    matrix = coo_matrix((vals, (rows, cols)), shape=(len(lb), nvar)).tocsr()
    ans = milp(
        c=objective,
        integrality=np.ones(nvar, dtype=int),
        bounds=Bounds(np.zeros(nvar), np.ones(nvar)),
        constraints=LinearConstraint(matrix, np.asarray(lb), np.asarray(ub)),
        options={"mip_rel_gap": 0.0, "presolve": True},
    )
    if not ans.success:
        raise RuntimeError(f"MILP lower-bound search failed: {ans.message}")
    return int(round(float(ans.fun)))


def search_design(
    pools: Sequence[Sequence[Candidate]],
    scenarios: Sequence[Scenario],
    alpha_star: float,
    target: float,
    equal_n: bool,
    common_thresholds: bool,
    error_mode: str = "product",
) -> tuple[list[Candidate], dict]:
    j_count = len(pools)
    if error_mode not in {"product", "bonferroni"}:
        raise ValueError("error_mode must be 'product' or 'bonferroni'")
    budget = -math.log1p(-alpha_star) if error_mode == "product" else alpha_star

    if common_thresholds:
        lookup = [{c.key(): c for c in pool} for pool in pools]
        keys = sorted(set.intersection(*(set(x) for x in lookup)))
        best = None
        for key in keys:
            design = [lookup[j][key] for j in range(j_count)]
            met = design_metrics(design, scenarios)
            if met["fwer"] <= alpha_star + 1e-12 and met["min_g1"] >= target - 1e-12:
                if best is None or lex_objective(design, scenarios) < lex_objective(best, scenarios):
                    best = design
        if best is None:
            raise RuntimeError("No feasible common-threshold design")
        return best, design_metrics(best, scenarios)

    best_design: list[Candidate] | None = None
    best_obj = None
    min_total = milp_minimum_total(
        pools, scenarios, alpha_star, target, equal_n=equal_n, error_mode=error_mode
    )
    max_total = sum(max(c.n for c in pool) for pool in pools)

    for total_n in range(min_total, max_total + 1):
        by_n = []
        feasible_total = True
        for pool in pools:
            d: dict[int, list[Candidate]] = {}
            for c in pool:
                d.setdefault(c.n, []).append(c)
            by_n.append(d)
            if not d:
                feasible_total = False
        if not feasible_total:
            continue

        allocations = []
        if equal_n:
            if total_n % j_count:
                continue
            n = total_n // j_count
            if all(n in d for d in by_n):
                allocations = [(n,) * j_count]
        else:
            ranges = [sorted(d) for d in by_n]
            suffix_min = [0] * (j_count + 1)
            suffix_max = [0] * (j_count + 1)
            for j in range(j_count - 1, -1, -1):
                suffix_min[j] = suffix_min[j + 1] + ranges[j][0]
                suffix_max[j] = suffix_max[j + 1] + ranges[j][-1]

            def build_alloc(j: int, remaining: int, prefix: tuple[int, ...]) -> None:
                if j == j_count:
                    if remaining == 0:
                        allocations.append(prefix)
                    return
                low_rest = suffix_min[j + 1]
                high_rest = suffix_max[j + 1]
                for n in ranges[j]:
                    rem = remaining - n
                    if rem < low_rest:
                        break
                    if rem > high_rest:
                        continue
                    build_alloc(j + 1, rem, prefix + (n,))

            build_alloc(0, total_n, ())
        if not allocations:
            continue

        for ns in allocations:
            local_pools = [by_n[j][ns[j]] for j in range(j_count)]
            # Search difficult doses first, then restore the scientific order.
            order = sorted(range(j_count), key=lambda j: len(local_pools[j]))
            ordered_pools = [local_pools[j] for j in order]
            ordered_scenarios = [
                Scenario(sc.name, sc.rho, tuple(sc.p_t[j] for j in order), tuple(sc.p_e[j] for j in order), tuple(sc.admissible[j] for j in order))
                for sc in scenarios
            ]

            def rec(depth: int, chosen: list[Candidate], cost: float) -> None:
                nonlocal best_design, best_obj
                if cost > budget + 1e-14:
                    return
                if depth < j_count and _optimistic_g1(chosen, depth, ordered_pools, ordered_scenarios) < target - 1e-12:
                    return
                if depth == j_count:
                    restored = [None] * j_count
                    for pos, original in enumerate(order):
                        restored[original] = chosen[pos]
                    design = list(restored)
                    met = design_metrics(design, scenarios)
                    if met["fwer"] <= alpha_star + 1e-12 and met["min_g1"] >= target - 1e-12:
                        obj = lex_objective(design, scenarios)
                        if best_obj is None or obj < best_obj:
                            best_obj = obj
                            best_design = design
                    return
                for cand in ordered_pools[depth]:
                    increment = cand.r if error_mode == "product" else cand.a
                    rec(depth + 1, chosen + [cand], cost + increment)

            rec(0, [], 0.0)

        if best_design is not None:
            return best_design, design_metrics(best_design, scenarios)

    raise RuntimeError("No feasible design in the requested search range")


def solve_configuration(
    j_count: int,
    alpha_star: float,
    target: float,
    phi_e1: float,
    phi_e0: float,
    phi_t1: float = 0.20,
    phi_t0: float = 0.40,
    n_min: int = 10,
    n_max: int = 150,
) -> dict:
    started = time.time()
    scenarios = make_scenarios(j_count, phi_t1, phi_t0, phi_e1, phi_e0)
    cache = ProbabilityCache(scenarios, n_max=n_max)
    specs = {
        "Tabata": (True, True, alpha_star / j_count, "product"),
        "ROED-EA": (True, False, None, "product"),
        # The manuscript's sotorasib table and its footnote define ROED-B by
        # the Bonferroni sum constraint sum_j a_j <= alpha*, not by equal
        # per-dose caps.  This is also the definition that reproduces N=87.
        "ROED-B": (False, False, None, "bonferroni"),
        "ROED": (False, False, None, "product"),
    }
    result = {
        "setting": {
            "J": j_count,
            "alpha": alpha_star,
            "target": target,
            "phi_T1": phi_t1,
            "phi_T0": phi_t0,
            "phi_E1": phi_e1,
            "phi_E0": phi_e0,
            "n_min": n_min,
            "n_max": n_max,
        },
        "scenarios": [sc.name for sc in scenarios],
        "designs": {},
    }
    pool_cache: dict[float | None, list[list[Candidate]]] = {}
    for method, (equal_n, common, cap, error_mode) in specs.items():
        method_started = time.time()
        print(f"[{method}] preparing candidates", flush=True)
        if cap not in pool_cache:
            raw_pools, pruned_pools = [], []
            for j in range(j_count):
                raw = generate_candidates(j, scenarios, cache, phi_t0, phi_e0, alpha_star, target, n_min, n_max, cap)
                raw_pools.append(raw)
                pruned_pools.append(prune_dominated(raw, j, scenarios))
            pool_cache[cap] = pruned_pools
        pools = pool_cache[cap]
        print(f"[{method}] retained candidates {[len(x) for x in pools]}; searching", flush=True)
        design, met = search_design(
            pools, scenarios, alpha_star, target, equal_n, common, error_mode=error_mode
        )
        result["designs"][method] = {
            "rules": [
                {"n": c.n, "mT": c.m_t, "mE": c.m_e, "aT": c.a_t, "aE": c.a_e, "a": c.a}
                for c in design
            ],
            **{k: v for k, v in met.items() if k not in ("g1", "g2")},
            "g1": dict(zip(result["scenarios"], met["g1"])),
            "g2": dict(zip(result["scenarios"], met["g2"])),
            "candidate_counts": [len(x) for x in pools],
        }
        print(f"[{method}] N={met['total_n']} FWER={met['fwer']:.8f} minG1={met['min_g1']:.8f} elapsed={time.time()-method_started:.2f}s", flush=True)
    result["elapsed_seconds"] = time.time() - started
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--J", type=int, required=True)
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--target", type=float, required=True)
    parser.add_argument("--phi-e1", type=float, required=True)
    parser.add_argument("--phi-e0", type=float, required=True)
    parser.add_argument("--n-max", type=int, default=150)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ans = solve_configuration(args.J, args.alpha, args.target, args.phi_e1, args.phi_e0, n_max=args.n_max)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(ans, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(args.output), "elapsed_seconds": ans["elapsed_seconds"], "designs": {k: {"total_n": v["total_n"], "fwer": v["fwer"], "min_g1": v["min_g1"], "counts": v["candidate_counts"]} for k, v in ans["designs"].items()}}, indent=2))


if __name__ == "__main__":
    main()
