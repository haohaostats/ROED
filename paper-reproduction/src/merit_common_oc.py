"""Re-evaluate frozen MERIT rules under the manuscript's Bernoulli model.

The MERIT design search is not repeated.  For each saved (n, mT, mE) rule,
this script simulates the common bivariate-Bernoulli data-generating model,
applies MERIT's dose-order PAVA transformation, and records FWER, G1, and G2.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "_checkpoints" / "study1" / "merit"
OUTPUT = ROOT / "_checkpoints" / "study1" / "merit_common_bernoulli_oc.csv"


def pava_rows(x: np.ndarray) -> np.ndarray:
    """Equal-weight increasing isotonic regression, row by row."""
    rows, j_count = x.shape
    out = np.zeros((rows, j_count), dtype=float)
    for i in range(j_count):
        outer = np.full(rows, -np.inf)
        for a in range(i + 1):
            inner = np.full(rows, np.inf)
            for b in range(i, j_count):
                inner = np.minimum(inner, x[:, a:b + 1].mean(axis=1))
            outer = np.maximum(outer, inner)
        out[:, i] = outer
    return out


def cells(p_t: float, p_e: float, rho: float) -> tuple[float, ...]:
    q = p_t * p_e + rho * math.sqrt(p_t * (1.0 - p_t) * p_e * (1.0 - p_e))
    lower, upper = max(0.0, p_t + p_e - 1.0), min(p_t, p_e)
    if q < lower - 1e-12 or q > upper + 1e-12:
        raise ValueError((p_t, p_e, rho, q, lower, upper))
    q = min(max(q, lower), upper)
    return q, p_t - q, p_e - q, 1.0 - p_t - p_e + q


def evaluate_one(args: tuple[str, int]) -> list[dict]:
    path_text, nsim = args
    path = Path(path_text)
    item = json.loads(path.read_text(encoding="utf-8"))
    st, design = item["input"], item["design"]
    base_seed = int(st["seed"]) + 700_001
    result: list[dict] = []
    for rho_index, oc in enumerate(item["operating_characteristics"].values()):
        rho = float(oc["rho"])
        for scenario_index, row in enumerate(oc["rows"]):
            pars = [float(x) for x in row["parameters"]]
            # A separate reproducible stream per configuration/rho/scenario.
            seed = (base_seed + 104_729 * (rho_index + 1) + 1_009 * (scenario_index + 1)) % (2**32 - 1)
            rng = np.random.default_rng(seed)
            tox = np.empty((nsim, int(st["J"])), dtype=float)
            eff = np.empty_like(tox)
            admissible = np.zeros(int(st["J"]), dtype=bool)
            for j in range(int(st["J"])):
                p_t, p_e = pars[2*j], pars[2*j + 1]
                draws = rng.multinomial(int(design["n_per_arm"]), cells(p_t, p_e, rho), size=nsim)
                tox[:, j] = draws[:, 0] + draws[:, 1]
                eff[:, j] = draws[:, 0] + draws[:, 2]
                admissible[j] = abs(p_t - float(st["phi_T1"])) < 1e-12 and abs(p_e - float(st["phi_E1"])) < 1e-12
            selected = (pava_rows(tox) <= int(design["mT"])) & (pava_rows(eff) >= int(design["mE"]))
            any_selected = selected.any(axis=1)
            if admissible.any():
                any_true = selected[:, admissible].any(axis=1)
                no_false = ~selected[:, ~admissible].any(axis=1) if (~admissible).any() else np.ones(nsim, dtype=bool)
                g1 = float(np.mean(any_true & no_false))
                g2 = float(np.mean(any_true))
                fwer = math.nan
                metric = "alternative"
            else:
                fwer = float(np.mean(any_selected))
                g1 = g2 = math.nan
                metric = "null"
            result.append({
                "configuration": path.stem, "J": int(st["J"]),
                "alpha": float(st["alpha"]), "target": float(st["target"]),
                "phi_E1": float(st["phi_E1"]), "phi_E0": float(st["phi_E0"]),
                "rho": rho, "scenario": row["scenario"], "kind": metric,
                "simulations": nsim, "seed": int(seed),
                "fwer": fwer, "g1": g1, "g2": g2,
            })
    return result


def recompute(simulations: int = 100_000, workers: int = 4) -> None:
    files = sorted(INPUT.glob("J*.json"))
    if len(files) != 24:
        raise RuntimeError(f"Expected 24 frozen MERIT designs; found {len(files)}")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        chunks = list(pool.map(evaluate_one, [(str(p), simulations) for p in files]))
    rows = [row for chunk in chunks for row in chunk]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(OUTPUT)
    print(f"wrote {len(rows)} rows to {OUTPUT}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulations", type=int, default=100_000)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    recompute(args.simulations, args.workers)


if __name__ == "__main__":
    main()
