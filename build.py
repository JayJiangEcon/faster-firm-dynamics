#!/usr/bin/env python3
"""Prepare or compile the integrated RE or DE Fortran evaluator."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", choices=("RE", "DE"))
    parser.add_argument("--mode", choices=("direct", "sigma-profile"), default="direct")
    parser.add_argument("--output", type=Path, required=True, help="A new build directory")
    parser.add_argument("--prepare-only", action="store_true", help="Generate source without compiling")
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--flags", default="-O3 -funroll-loops -fopenmp -m64")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "results/integrated-manifest.json").read_text())
    output = args.output.resolve()
    source = ROOT / "src/fortran"

    for name, expected in manifest["models"][args.model]["files"].items():
        path = ROOT / name
        if not path.is_file() or digest(path) != expected:
            parser.error(f"Missing or changed packaged source: {path}")
    if output.exists():
        parser.error(f"Output already exists; choose a new directory: {output}")
    compiler = None if args.prepare_only else shutil.which(args.fc)
    if not args.prepare_only and not compiler:
        parser.error(f"Compiler not found: {args.fc}")
    output.mkdir(parents=True)
    shutil.copyfile(source / "base_lib.f90", output / "base_lib.f90")
    for name in ("diag_firms_solve.f90", "evaluate.f90"):
        shutil.copyfile(source / args.model / name, output / name)
    solver = output / "diag_firms_solve.f90"
    actual = digest(solver)
    record = {
        "model": args.model, "mode": args.mode,
        "source_manifest_sha256": digest(ROOT / "results/integrated-manifest.json"),
        "generated_solver_sha256": actual, "recorded_source_match": True,
        "packaged_source_sha256": manifest["models"][args.model]["files"],
        "driver_sha256": digest(output / "evaluate.f90"), "compiled": False,
    }
    if compiler:
        flags = shlex.split(args.flags)
        env = os.environ.copy()
        env["PATH"] = str(Path(compiler).parent) + os.pathsep + env.get("PATH", "")
        if sys.platform == "darwin" and shutil.which("xcrun"):
            env.setdefault("SDKROOT", subprocess.check_output(
                ["xcrun", "--show-sdk-path"], text=True).strip())
        for name in ("base_lib", "diag_firms_solve", "evaluate"):
            subprocess.run([compiler, *flags, "-I", ".", "-c", f"{name}.f90", "-o", f"{name}.o"],
                           cwd=output, env=env, check=True)
        binary = output / f"rcc-{args.model.lower()}-{args.mode}"
        subprocess.run([compiler, *flags, "base_lib.o", "diag_firms_solve.o", "evaluate.o", "-o", str(binary)],
                       cwd=output, env=env, check=True)
        record.update({"compiled": True, "compiler": compiler,
                       "compiler_sha256": digest(Path(compiler)),
                       "compiler_version": subprocess.check_output([compiler, "--version"], text=True, env=env),
                       "flags": flags, "binary": binary.name, "binary_sha256": digest(binary)})
    (output / "build.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Recorded source identity matched: {args.model} {args.mode}")
    print(output / "build.json")


if __name__ == "__main__":
    main()
