#!/usr/bin/env python3
"""Search explicit parameter bounds using the integrated profile evaluator."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

from build import digest

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--spec", type=Path, required=True,
                        help="JSON with model, bounds, sigma_pi_bounds, seed, max_evaluations")
    parser.add_argument("--inputs-dir", type=Path, default=ROOT / "data/fixed")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args()
    import numpy as np
    from scipy.optimize import differential_evolution

    spec = json.loads(args.spec.read_text())
    build_dir = args.build.resolve()
    build = json.loads((build_dir / "build.json").read_text())
    if not build["compiled"] or build["mode"] != "sigma-profile":
        parser.error("Search requires a compiled sigma-profile build")
    if spec["model"] != build["model"]:
        parser.error("Model does not match the build")
    n = {"RE": 8, "DE": 9}[spec["model"]]
    bounds = np.asarray(spec["bounds"], dtype=float)
    sigma_bounds = np.asarray(spec["sigma_pi_bounds"], dtype=float)
    if (bounds.shape != (n, 2) or not np.isfinite(bounds).all()
            or not np.all(bounds[:, 0] < bounds[:, 1])):
        parser.error(f"Provide {n} finite lower/upper pairs with lower < upper")
    if (sigma_bounds.shape != (2,) or not np.isfinite(sigma_bounds).all()
            or not 0 <= sigma_bounds[0] < sigma_bounds[1]):
        parser.error("Invalid sigma_pi_bounds")
    budget, seed = spec["max_evaluations"], spec["seed"]
    if type(budget) is not int or budget < 1 or type(seed) is not int or seed < 0:
        parser.error("Use a positive integer budget and nonnegative integer seed")
    if args.threads < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("Threads and timeout must be positive")
    output = args.output.resolve()
    if output.exists():
        parser.error("Choose a new output directory")
    output.mkdir(parents=True)
    record = {"spec": spec, "algorithm": "scipy.optimize.differential_evolution",
              "strategy": "best1bin", "population_multiplier": 5,
              "initialization": "latinhypercube", "polish": False,
              "workers": 1, "evaluations": 0, "best": None,
              "build_record_sha256": digest(build_dir / "build.json"),
              "status": "running"}
    destination = output / "search.json"

    def save():
        temporary = output / "search.tmp"
        temporary.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
        temporary.replace(destination)

    class BudgetReached(Exception):
        pass

    def objective(parameters):
        if record["evaluations"] >= budget:
            raise BudgetReached()
        i = record["evaluations"]
        point = output / f"point-{i:06d}.json"
        point.write_text(json.dumps({"model": spec["model"], "parameters": parameters.tolist()}) + "\n")
        runtime = output / f"evaluation-{i:06d}"
        record["evaluations"] += 1
        save()
        command = [sys.executable, str(ROOT / "run.py"), "--build", str(build_dir),
                   "--inputs-dir", str(args.inputs_dir.resolve()), "--parameters", str(point),
                   "--output", str(runtime), "--sigma-bounds", *map(str, sigma_bounds),
                   "--threads", str(args.threads), "--timeout", str(args.timeout)]
        with (output / f"evaluation-{i:06d}.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        result = json.loads((runtime / "result.json").read_text())
        if record["best"] is None or result["objective"] < record["best"]["objective"]:
            record["best"] = {"evaluation": i, "parameters": result["parameters"],
                              "objective": result["objective"],
                              "result": str(runtime.relative_to(output) / "result.json")}
        save()
        return result["objective"]

    save()
    try:
        result = differential_evolution(objective, bounds, seed=seed, popsize=5,
                                        maxiter=budget, polish=False, workers=1,
                                        updating="immediate", tol=1e-7, atol=0)
        record.update(status="finished", message=str(result.message), success=bool(result.success))
    except BudgetReached:
        record.update(status="evaluation_budget_reached", success=False)
    except (Exception, KeyboardInterrupt) as exc:
        record.update(status="failed", message=str(exc), success=False)
        save()
        raise
    save()
    print(destination)


if __name__ == "__main__":
    main()
