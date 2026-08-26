"""Create the SHA-256 inventory for bundled source and frozen CSV files."""

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main() -> None:
    files = sorted(
        [p for p in (ROOT / "src").rglob("*") if p.is_file() and p.suffix in {".py", ".R"}]
        + list((ROOT / "reference-results").rglob("*.csv"))
        + [ROOT / "config.json"]
    )
    with (ROOT / "SHA256SUMS.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("path", "sha256"))
        for path in files:
            writer.writerow((path.relative_to(ROOT).as_posix(), hashlib.sha256(path.read_bytes()).hexdigest()))


if __name__ == "__main__":
    main()
