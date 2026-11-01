#!/usr/bin/env python3
"""Display existing Faster Firm Dynamics timing measurements without running a model."""
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    data = json.loads((root / "results/benchmarks.json").read_text())
    print("Faster Firm Dynamics: recorded performance on the Real Credit Cycles models\n")
    print("| Comparison | Wall speedup | Aggregation |")
    print("| --- | ---: | --- |")
    for row in data["comparisons"]:
        print(f"| {row['comparison']} | {row['wall_speedup']:.4f}x | {row['aggregation']} |")
    print("\n" + data["aggregation_note"])
    print("See docs/benchmarks.md and docs/accuracy.md for the numerical evidence.")


if __name__ == "__main__":
    main()
