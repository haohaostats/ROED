# ROED paper reproduction

This directory contains the self-contained computation code and frozen numeric
results for Numerical Studies 1--4 and the sotorasib application. 

## Obtain the manuscript results without recomputation

The exact CSV files used for the current manuscript are stored under
`reference-results/`. To copy them into `results/` immediately:

```bash
python install_reference_results.py
python verify.py
```

This operation performs no simulation or design search.

## Recompute from source

Install Python dependencies with `pip install -r requirements.txt`, install the
ROED R package, and ensure `Rscript` is on `PATH` (or set `ROED_RSCRIPT`). Then
run one of:

```text
run_all.bat       Windows Command Prompt
run_all.ps1       Windows PowerShell
./run_all.sh      macOS or Linux
```

The workflow is checkpointed under `_checkpoints/`. Re-running the same command
skips completed tasks and continues after interruption. `--workers N` controls
the number of worker processes; use `--workers 1` on a limited computer.

The full run comprises 24 primary exact configurations, 24 local-R MERIT
configurations, 24 ROED-S configurations, the exact-search timing benchmark,
the common-model MERIT operating-characteristic evaluation, CSV aggregation,
figure-data CSV export, and the case application. It can require many hours.

Useful partial commands include:

```bash
python run.py --stage exact --workers 1
python run.py --stage merit --workers 1
python run.py --stage roed-s --workers 1
python run.py --stage aggregate
python run.py --stage case
```

After a full recomputation, `python verify.py` compares generated CSVs with the
frozen manuscript results. Runtime columns are excluded because they depend on
hardware. MERIT Monte Carlo quantities use a small numerical tolerance.

## Provenance

`src/exact_engine.py`, `src/merit_local.R`, and `src/roed_s.py` are copies of
the engines used for the current manuscript computations. `SHA256SUMS.csv`
records the archived source and frozen-result hashes. No LaTeX drawing or table
code is included.
