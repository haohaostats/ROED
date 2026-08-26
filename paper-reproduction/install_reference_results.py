"""Copy the manuscript's frozen CSV results into results/ without recomputing."""

from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "reference-results"
TARGET = ROOT / "results"


def main() -> None:
    for source in SOURCE.rglob("*.csv"):
        target = TARGET / source.relative_to(SOURCE)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    print(f"Installed {len(list(SOURCE.rglob('*.csv')))} frozen CSV files in {TARGET}")


if __name__ == "__main__":
    main()
