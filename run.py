#!/usr/bin/env python3
"""Run one integrated direct or conditional-profile evaluation."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parent


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True, help="Directory produced by build.py")
    parser.add_argument("--inputs-dir", type=Path, default=ROOT / "data/fixed")
    parser.add_argument("--sigma-bounds", type=float, nargs=2, metavar=("LOWER", "UPPER"))
    parser.add_argument("--parameters", type=Path, required=True, help="JSON with model and parameters")
    parser.add_argument("--output", type=Path, required=True, help="A new runtime directory")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--timeout", type=float, default=300)
    args = parser.parse_args()
    manifest = json.loads((ROOT / "results/integrated-manifest.json").read_text())
    build_dir = args.build.resolve()
    build = json.loads((build_dir / "build.json").read_text())
    point = json.loads(args.parameters.read_text())
    if not build["compiled"]:
        parser.error("This runner requires a compiled build")
    if point["model"] != build["model"]:
        parser.error("Parameter model does not match the build")
    params = list(point["parameters"])
    count = {"RE": 9, "DE": 10}[build["model"]]
    if not all(isinstance(x, (int, float)) and not isinstance(x, bool)
               and math.isfinite(x) for x in params):
        parser.error("Provide finite numerical parameters")
    if build["mode"] == "sigma-profile":
        if (args.sigma_bounds is None or not all(map(math.isfinite, args.sigma_bounds))
                or not 0 <= args.sigma_bounds[0] < args.sigma_bounds[1]):
            parser.error("Profile evaluation requires --sigma-bounds LOWER UPPER, 0 <= LOWER < UPPER")
        try:
            import numpy
        except ImportError:
            parser.error("Profile selection requires NumPy")
        if len(params) == count-1:
            params.append(0.0)
        if len(params) != count or params[-1] != 0:
            parser.error("Supply the policy coordinates, or a complete vector ending in zero")
    elif len(params) != count or args.sigma_bounds is not None:
        parser.error("Direct evaluation requires a complete parameter vector and no sigma bounds")
    if args.threads < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("Threads and timeout must be positive")
    if build["source_manifest_sha256"] != digest(ROOT / "results/integrated-manifest.json"):
        parser.error("Build belongs to a different source manifest")
    binary = build_dir / build["binary"]
    if digest(binary) != build["binary_sha256"]:
        parser.error("Binary differs from build.json")
    inputs = args.inputs_dir.resolve()
    for name, expected in manifest["runtime_input_sha256"].items():
        path = inputs / name
        if not path.is_file() or digest(path) != expected:
            parser.error(f"Missing or changed fixed input: {path}")
    output = args.output.resolve()
    if output.exists():
        parser.error(f"Output already exists; choose a new directory: {output}")
    output.mkdir(parents=True)
    for name in manifest["runtime_input_sha256"]:
        shutil.copyfile(inputs / name, output / name)
    # Suppress inherited solver overrides so a shell cannot silently alter fidelity.
    env = {k: v for k, v in os.environ.items() if not k.startswith("RCC_")}
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC="FALSE")
    if build["mode"] == "sigma-profile":
        env["RCC_SIGMAPI_PROFILE"] = "1"
    record = {"model": build["model"], "mode": build["mode"], "parameters": params,
              "sigma_pi_bounds": args.sigma_bounds,
              "parameter_hex": [float(x).hex() for x in params],
              "build_record_sha256": digest(build_dir / "build.json"),
              "binary_sha256": build["binary_sha256"],
              "input_sha256": manifest["runtime_input_sha256"], "threads": args.threads,
              "timeout_seconds": args.timeout, "timed_out": False}
    start = time.perf_counter()
    with (output / "stdout.txt").open("w") as stdout, (output / "stderr.txt").open("w") as stderr:
        process = subprocess.Popen([str(binary), *[repr(float(x)) for x in params]],
                                   cwd=output, env=env, stdout=stdout, stderr=stderr)
        record["pid"] = process.pid
        try:
            process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
            record["timed_out"] = True
        record["returncode"] = process.returncode
    record["wall_seconds"] = time.perf_counter() - start
    record["stdout_sha256"] = digest(output / "stdout.txt")
    record["stderr_sha256"] = digest(output / "stderr.txt")
    (output / "run.json").write_text(json.dumps(record, indent=2) + "\n")
    print(output / "run.json")
    if record["timed_out"] or record["returncode"]:
        raise SystemExit("Evaluation did not complete successfully; inspect stdout.txt and stderr.txt")
    from sigma_profile import tagged_values, from_stdout
    stdout_text = (output / "stdout.txt").read_text()
    try:
        if build["mode"] == "sigma-profile":
            result = from_stdout(stdout_text, inputs, args.sigma_bounds)
            result["parameters"] = params[:-1] + [result["sigma_pi"]]
        else:
            moments = tagged_values(stdout_text, "RCC_MOMENTS", 17)
            objective = tagged_values(stdout_text, "RCC_RAW_Q", 1)[0]
            if not math.isfinite(objective):
                raise ValueError("Nonfinite objective")
            result = {"parameters": params, "objective": objective,
                      "moments": [x if math.isfinite(x) else None for x in moments]}
        result.update(model=build["model"], mode=build["mode"],
                      parameter_hex=[float(x).hex() for x in result["parameters"]],
                      run_record_sha256=digest(output / "run.json"))
        (output / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    except (ValueError, OSError) as exc:
        record["result_error"] = str(exc)
        (output / "run.json").write_text(json.dumps(record, indent=2) + "\n")
        raise SystemExit(f"Could not form a result: {exc}") from exc
    print(output / "result.json")


if __name__ == "__main__":
    main()
