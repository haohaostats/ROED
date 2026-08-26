"""Generate exact surface data for the representative study-1 figure."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPRESENTATIVE = ROOT / "_checkpoints" / "study1" / "exact" / "J2_a0p10_p0p80_e0p30_0p10.json"
MERIT_REPRESENTATIVE = ROOT / "_checkpoints" / "study1" / "merit" / "J2_a0p10_p0p80_e0p30_0p10.json"
OUT = ROOT / "results" / "figure-data" / "study1"


def select_probability(n: int, m_t: int, m_e: int, p_t: float, p_e: float) -> float:
    """Worst association at fixed marginals: q at the lower Frechet bound."""
    q = max(0.0, p_t + p_e - 1.0)
    eta00 = 1.0 - p_t - p_e + q
    eta10 = p_t - q
    eta01 = p_e - q
    eta11 = q
    dist = np.ones((1, 1), dtype=float)
    for r in range(1, n + 1):
        nxt = np.zeros((r + 1, r + 1), dtype=float)
        nxt[:r, :r] += eta00 * dist
        nxt[1:, :r] += eta10 * dist
        nxt[:r, 1:] += eta01 * dist
        nxt[1:, 1:] += eta11 * dist
        dist = nxt
    return float(dist[: m_t + 1, m_e :].sum())


def joint_count_pmf(n: int, p_t: float, p_e: float, rho: float = 0.0) -> np.ndarray:
    """Exact joint PMF of the toxicity and efficacy counts."""
    scale = np.sqrt(p_t * (1.0 - p_t) * p_e * (1.0 - p_e))
    q = p_t * p_e + rho * scale
    q = min(max(q, max(0.0, p_t + p_e - 1.0)), min(p_t, p_e))
    eta00 = 1.0 - p_t - p_e + q
    eta10 = p_t - q
    eta01 = p_e - q
    eta11 = q
    dist = np.ones((1, 1), dtype=float)
    for r in range(1, n + 1):
        nxt = np.zeros((r + 1, r + 1), dtype=float)
        nxt[:r, :r] += eta00 * dist
        nxt[1:, :r] += eta10 * dist
        nxt[:r, 1:] += eta01 * dist
        nxt[1:, 1:] += eta11 * dist
        dist = nxt
    return dist


def merit_surface(grid: np.ndarray) -> list[dict[str, str]]:
    """Exact rho=0 MERIT surface with the second dose at (0.40, 0.30)."""
    saved = json.loads(MERIT_REPRESENTATIVE.read_text(encoding="utf-8"))
    design = saved["design"]
    n = int(design["n_per_arm"])
    m_t, m_e = int(design["mT"]), int(design["mE"])

    counts = np.arange(n + 1)
    t1, e1 = np.meshgrid(counts, counts, indexing="ij")
    t1, e1 = t1.ravel(), e1.ravel()
    t2, e2 = t1.copy(), e1.copy()

    # Equal-weight two-dose PAVA, evaluated for every pair of arm-level
    # toxicity and efficacy counts.  The decision mask is data-independent.
    tox_violation = t1[:, None] > t2[None, :]
    eff_violation = e1[:, None] > e2[None, :]
    pav_t1 = np.where(tox_violation, (t1[:, None] + t2[None, :]) / 2.0, t1[:, None])
    pav_t2 = np.where(tox_violation, (t1[:, None] + t2[None, :]) / 2.0, t2[None, :])
    pav_e1 = np.where(eff_violation, (e1[:, None] + e2[None, :]) / 2.0, e1[:, None])
    pav_e2 = np.where(eff_violation, (e1[:, None] + e2[None, :]) / 2.0, e2[None, :])
    any_selected = ((pav_t1 <= m_t) & (pav_e1 >= m_e)) | ((pav_t2 <= m_t) & (pav_e2 >= m_e))

    fixed = joint_count_pmf(n, 0.40, 0.30, rho=0.0).ravel()
    conditional = any_selected @ fixed
    records: list[dict[str, str]] = []
    for p_e in grid:
        for p_t in grid:
            is_null = p_t >= 0.40 - 1e-12 or p_e <= 0.10 + 1e-12
            if is_null:
                probability = float(joint_count_pmf(n, p_t, p_e, rho=0.0).ravel() @ conditional)
                z = f"{probability:.10f}"
            else:
                z = "nan"
            records.append({"pT": f"{p_t:.6f}", "pE": f"{p_e:.6f}", "fwer": z})
    return records


def main() -> None:
    result = json.loads(REPRESENTATIVE.read_text(encoding="utf-8"))["payload"]
    OUT.mkdir(parents=True, exist_ok=True)
    grid = np.linspace(0.0, 1.0, 31)
    merit_records = merit_surface(grid)
    with (OUT / "merit.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["pT", "pE", "fwer"])
        writer.writeheader(); writer.writerows(merit_records)
    print(OUT / "merit.csv")
    for method in ("Tabata", "ROED-EA", "ROED-B", "ROED"):
        rules = result["designs"][method]["rules"]
        first, second = rules
        a2 = second["a"]
        records = []
        for p_e in grid:
            for p_t in grid:
                is_null = p_t >= 0.40 - 1e-12 or p_e <= 0.10 + 1e-12
                if is_null:
                    p1 = select_probability(first["n"], first["mT"], first["mE"], p_t, p_e)
                    fwer = 1.0 - (1.0 - p1) * (1.0 - a2)
                    z = f"{fwer:.10f}"
                else:
                    z = "nan"
                records.append({"pT": f"{p_t:.6f}", "pE": f"{p_e:.6f}", "fwer": z})
        path = OUT / f"{method.lower().replace('-', '_')}.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["pT", "pE", "fwer"])
            writer.writeheader(); writer.writerows(records)
        print(path)


if __name__ == "__main__":
    main()
