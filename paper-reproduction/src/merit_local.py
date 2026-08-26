"""Python adapter for the local R implementation of MERIT."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


SCRIPT = Path(__file__).with_name("merit_local.R")
R_CANDIDATES = (
    Path(r"D:\Soft\R-4.5.2\bin\x64\Rscript.exe"),
    Path(r"C:\Program Files\R\R-4.5.2\bin\Rscript.exe"),
)


def _rscript() -> str:
    override = os.environ.get("ROED_RSCRIPT")
    if override:
        return override
    discovered = shutil.which("Rscript")
    if discovered:
        return discovered
    for path in R_CANDIDATES:
        if path.exists():
            return str(path)
    return "Rscript"


def run_local_merit(spec: dict, destination: Path, merit_config: dict) -> dict:
    destination.mkdir(parents=True, exist_ok=True)
    started_at = datetime.now(timezone.utc).isoformat()
    request = {"spec": spec, "config": merit_config}
    with tempfile.TemporaryDirectory(prefix="merit_", dir=destination) as work:
        work_path = Path(work)
        input_path = work_path / "input.json"
        output_path = work_path / "output.json"
        input_path.write_text(json.dumps(request, ensure_ascii=False), encoding="utf-8")
        completed = subprocess.run(
            [_rscript(), str(SCRIPT), str(input_path), str(output_path)],
            capture_output=True, text=True, encoding="utf-8", errors="replace",
        )
        if completed.returncode:
            raise RuntimeError(
                "Local MERIT R engine failed\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        payload = json.loads(output_path.read_text(encoding="utf-8"))

    design = payload["design"]
    design.pop("search_trace", None)
    return {
        "schema_version": 1,
        "task_id": spec["task_id"],
        "status": "complete",
        "source": {
            "application": "local MERIT R implementation",
            "version": payload["engine"]["version"],
            "reference_doi": payload["engine"]["reference_doi"],
            "script": str(SCRIPT),
            "script_sha256": hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
            "r_version": payload["engine"]["language"],
        },
        "input": {key: spec[key] for key in (
            "J", "alpha", "target", "phi_T0", "phi_T1", "phi_E0", "phi_E1",
            "merit_correlation", "simulations", "seed",
        )},
        "design": design,
        "operating_characteristics": payload.get("operating_characteristics"),
        "started_at": started_at,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": payload["elapsed_seconds"],
    }
