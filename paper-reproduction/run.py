"""One-command, resumable and parallel runner for the formal simulations."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from functools import partial
from pathlib import Path


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.exact_adapter import ENGINE_PATH, engine_digest, load_engine
from src.aggregate import build_all
from src.export_figure_csv import main as export_figure_csv
from src.export_surface_csv import main as export_surface_csv
from src.io_utils import atomic_write_json
from src.merit_local import SCRIPT as MERIT_R_SCRIPT, run_local_merit
from src.roed_s import solve_roed_s
from src.merit_common_oc import recompute as recompute_merit_common_oc


def load_config() -> dict:
    return json.loads((ROOT / "config.json").read_text(encoding="utf-8"))


def label(spec: dict) -> str:
    raw = (f"J{spec['J']}_a{spec['alpha']:.2f}_p{spec['target']:.2f}_"
           f"e{spec['phi_E1']:.2f}_{spec['phi_E0']:.2f}")
    return raw.replace(".", "p")


def spec_digest(spec: dict) -> str:
    return hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def primary_specs(config: dict) -> list[dict]:
    grid = config["design_grid"]
    rows = []
    for e1, e0 in grid["efficacy_pairs"]:
        for j, alpha, target in itertools.product(grid["J"], grid["alpha"], grid["target_power"]):
            row = {
                "J": j, "alpha": alpha, "target": target,
                "phi_T1": grid["phi_T1"], "phi_T0": grid["phi_T0"],
                "phi_E1": e1, "phi_E0": e0,
                "n_min": grid["n_min"], "n_max": grid["n_max"],
            }
            row["label"] = label(row)
            rows.append(row)
    return rows


def valid_result(path: Path, spec: dict) -> bool:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        if value.get("status") != "complete" or value.get("task_id") != spec["task_id"]:
            return False
        recorded = value.get("spec_sha256")
        if recorded is not None:
            return recorded == spec_digest(spec)
        # Backward-compatible validation for the initial smoke-test artifacts.
        if "input" in value:
            return all(value["input"].get(k) == spec.get(k) for k in value["input"])
        setting = value.get("payload", {}).get("setting", {})
        mapping = {"alpha": "alpha", "target": "target", "J": "J",
                   "phi_E1": "phi_E1", "phi_E0": "phi_E0",
                   "phi_T1": "phi_T1", "phi_T0": "phi_T0"}
        return all(setting.get(dst) == spec.get(src) for src, dst in mapping.items())
    except (OSError, json.JSONDecodeError):
        return False


def exact_worker(spec: dict, output_text: str) -> str:
    output = Path(output_text)
    task_id = spec["task_id"]
    started = time.perf_counter()
    eng = load_engine()
    try:
        payload = eng.solve_configuration(
            spec["J"], spec["alpha"], spec["target"], spec["phi_E1"], spec["phi_E0"],
            phi_t1=spec["phi_T1"], phi_t0=spec["phi_T0"],
            n_min=spec["n_min"], n_max=spec["n_max"],
        )
        envelope = {
            "schema_version": 1, "task_id": task_id, "status": "complete",
            "spec_sha256": spec_digest(spec),
            "engine": {"path": str(ENGINE_PATH), "sha256": engine_digest()},
            "payload": payload, "elapsed_seconds": time.perf_counter() - started,
            "finished_at": datetime.now(timezone.utc).isoformat(),
        }
        atomic_write_json(output, envelope)
        return str(output)
    except Exception:
        raise RuntimeError(f"{task_id}\n{traceback.format_exc()}")


def merit_worker(spec: dict, output_text: str, merit_config: dict) -> str:
    output = Path(output_text)
    try:
        result = run_local_merit(spec, output.parent, merit_config)
        result["spec_sha256"] = spec_digest(spec)
        atomic_write_json(output, result)
        return str(output)
    except Exception:
        raise RuntimeError(f"{spec['task_id']}\n{traceback.format_exc()}")


def roed_s_worker(spec: dict, output_text: str, nominal_rho: float) -> str:
    output = Path(output_text)
    started = time.perf_counter()
    try:
        payload = solve_roed_s(spec, nominal_rho=nominal_rho)
        envelope = {
            "schema_version": 1, "task_id": spec["task_id"], "status": "complete",
            "spec_sha256": spec_digest(spec),
            "engine": {"path": str(ENGINE_PATH), "sha256": engine_digest()},
            "payload": payload, "elapsed_seconds": time.perf_counter() - started,
            "finished_at": datetime.now(timezone.utc).isoformat(),
        }
        atomic_write_json(output, envelope)
        return str(output)
    except Exception:
        raise RuntimeError(f"{spec['task_id']}\n{traceback.format_exc()}")


def run_pool(tasks, worker, workers: int, retries: int, dry_run: bool = False,
             retry_backoff_seconds: float = 0) -> None:
    task_ids = {spec["task_id"] for spec, _ in tasks}
    pending = [(spec, path) for spec, path in tasks if not valid_result(path, spec)]
    skipped = len(tasks) - len(pending)
    print(f"tasks={len(tasks)} pending={len(pending)} resumed/skipped={skipped}", flush=True)
    if dry_run or not pending:
        return
    attempts = {spec["task_id"]: 0 for spec, _ in pending}
    while pending:
        failed = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(worker, spec, str(path)): (spec, path) for spec, path in pending}
            for future in as_completed(futures):
                spec, path = futures[future]
                try:
                    future.result()
                    prior_failure = ROOT / "_checkpoints" / "state" / "failures" / f"{spec['task_id']}.log"
                    if prior_failure.exists():
                        prior_failure.unlink()
                    print(f"complete {spec['task_id']}", flush=True)
                except Exception as exc:
                    attempts[spec["task_id"]] += 1
                    error_path = ROOT / "_checkpoints" / "state" / "failures" / f"{spec['task_id']}.log"
                    error_path.parent.mkdir(parents=True, exist_ok=True)
                    error_path.write_text(str(exc), encoding="utf-8")
                    print(f"failed {spec['task_id']} attempt={attempts[spec['task_id']]}", flush=True)
                    if attempts[spec["task_id"]] <= retries:
                        failed.append((spec, path))
        pending = failed
        if pending and retry_backoff_seconds > 0:
            print(
                f"retrying {len(pending)} task(s) after "
                f"{retry_backoff_seconds:g}s backoff",
                flush=True,
            )
            time.sleep(retry_backoff_seconds)
    failure_dir = ROOT / "_checkpoints" / "state" / "failures"
    failure_logs = list(failure_dir.glob("*.log")) if failure_dir.exists() else []
    unresolved = [
        p for p in failure_logs
        if p.stem in task_ids and not any(
            valid_result(path, spec)
            for spec, path in tasks if spec["task_id"] == p.stem
        )
    ]
    if unresolved:
        raise RuntimeError(f"Unresolved tasks: {', '.join(p.stem for p in unresolved)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=(
        "all", "exact", "merit", "merit-oc", "roed-s", "benchmark", "aggregate", "case"
    ), default="all")
    parser.add_argument("--workers", type=int)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--smoke", action="store_true", help="Run one configuration with 100 MERIT simulations")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--match", help="Run only configurations whose label or task id contains this text"
    )
    args = parser.parse_args()
    config = load_config()
    merit_code_sha256 = hashlib.sha256(MERIT_R_SCRIPT.read_bytes()).hexdigest()
    merit_config_sha256 = hashlib.sha256(
        json.dumps(config["merit"], sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    specs = primary_specs(config)
    if args.match:
        specs = [spec for spec in specs if args.match in spec["label"]]
        if not specs:
            parser.error(f"--match did not select a primary configuration: {args.match}")
    if args.smoke:
        specs = specs[:1]
    workers = args.workers or config["default_workers"]
    retries = config["max_retries"]
    output_root = ROOT / "_checkpoints"
    if args.smoke:
        output_root = output_root / "smoke"

    manifest = {
        "schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
        "config_sha256": hashlib.sha256((ROOT / "config.json").read_bytes()).hexdigest(),
        "exact_engine_sha256": engine_digest(), "exact_engine_path": str(ENGINE_PATH),
        "python": sys.version, "command": sys.argv,
    }
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    atomic_write_json(ROOT / "_checkpoints" / "state" / "runs" / f"{run_id}.json", manifest)

    if args.stage in ("all", "exact"):
        tasks = []
        for base in specs:
            spec = dict(base, task_id=f"exact_{base['label']}")
            path = output_root / "study1" / "exact" / f"{base['label']}.json"
            if args.force and path.exists():
                path.unlink()
            tasks.append((spec, path))
        run_pool(tasks, exact_worker, workers, retries, args.dry_run)

    if args.stage in ("all", "merit"):
        tasks = []
        for base in specs:
            spec = dict(
                base, task_id=f"merit_{base['label']}",
                merit_correlation=config["merit"]["design_correlation"],
                simulations=100 if args.smoke else config["merit"]["design_simulations"],
                oc_simulations=100 if args.smoke else config["merit"]["oc_simulations"],
                run_oc=True,
                merit_engine=config["merit"]["primary_search"],
                merit_code_sha256=merit_code_sha256,
                merit_config_sha256=merit_config_sha256,
                seed=config["master_seed"],
            )
            path = output_root / "study1" / "merit" / f"{base['label']}.json"
            if args.force and path.exists():
                path.unlink()
            tasks.append((spec, path))
        merit_cfg = config["merit"]
        run_pool(
            tasks,
            partial(merit_worker, merit_config=merit_cfg),
            min(workers, merit_cfg["max_parallel_sessions"]),
            merit_cfg.get("max_retries", retries), args.dry_run,
            merit_cfg.get("retry_backoff_seconds", 0),
        )

    if args.stage in ("all", "roed-s"):
        roed_s_tasks = []
        for base in specs:
            spec = dict(base, task_id=f"roed_s_{base['label']}")
            path = output_root / "study4" / "roed_s" / f"{base['label']}.json"
            if args.force and path.exists():
                path.unlink()
            roed_s_tasks.append((spec, path))
        run_pool(
            roed_s_tasks,
            partial(roed_s_worker, nominal_rho=config["merit"]["design_correlation"]),
            workers, retries, args.dry_run,
        )

    if args.stage in ("all", "benchmark"):
        benchmark_tasks = []
        repetitions = 1 if args.smoke else config["study4"]["timing_repetitions"]
        bench_specs = specs if not args.smoke else specs[:1]
        for base in bench_specs:
            for repetition in range(1, repetitions + 1):
                task_id = f"benchmark_{base['label']}_r{repetition:02d}"
                spec = dict(base, task_id=task_id, benchmark_repetition=repetition)
                path = output_root / "study4" / "benchmark" / f"{task_id}.json"
                if args.force and path.exists():
                    path.unlink()
                benchmark_tasks.append((spec, path))
        run_pool(benchmark_tasks, exact_worker, workers, retries, args.dry_run)

    if args.stage in ("all", "merit-oc", "aggregate") and not args.dry_run and not args.smoke:
        print("evaluating frozen MERIT rules under the common Bernoulli model (no design search)", flush=True)
        recompute_merit_common_oc(simulations=100_000, workers=workers)

    if args.stage in ("all", "aggregate") and not args.dry_run and not args.smoke:
        build_all(ROOT, config)
        export_figure_csv()
        export_surface_csv()
        print("complete aggregate outputs for Studies 1--4", flush=True)

    if args.stage in ("all", "case") and not args.dry_run and not args.smoke:
        rscript = os.environ.get("ROED_RSCRIPT") or shutil.which("Rscript") or "Rscript"
        subprocess.run(
            [rscript, str(ROOT / "src" / "case_application.R"),
             str(ROOT / "results" / "case-application")], check=True
        )
        print("complete case-application CSV outputs", flush=True)


if __name__ == "__main__":
    main()
