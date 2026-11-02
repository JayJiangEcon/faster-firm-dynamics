#!/usr/bin/env python3
"""Compare independently written arrays without importing either implementation."""
import argparse
import json
from pathlib import Path


def compare(cpu, gpu, item_size, indices=None):
    a, b = cpu.read_bytes(), gpu.read_bytes()
    if len(a) != len(b) or len(a) % item_size:
        raise ValueError(f"incompatible array lengths: {cpu}, {gpu}")
    if indices is None and a == b:
        return {"entries": len(a) // item_size, "bit_mismatches": 0}
    indices = range(len(a) // item_size) if indices is None else indices
    count = bad = 0
    for i in indices:
        count += 1
        j = i * item_size
        bad += a[j:j+item_size] != b[j:j+item_size]
    return {"entries": count, "bit_mismatches": bad}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu", type=Path, required=True)
    parser.add_argument("--gpu", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--final-moments", type=Path, help="Independently written final pipeline moments")
    args = parser.parse_args()
    mapping = {
        "policy": ("RCC_kbpolind.bin", "policy.bin", 4),
        "default_policy": ("RCC_dfpolind.bin", "policy_default.bin", 4),
        "value": ("RCC_Vtheta.bin", "value.bin", 8),
        "old_value": ("RCC_Vthetaold.bin", "old.bin", 8),
        "price": ("RCC_qmat.bin", "q.bin", 8),
        "distribution": ("RCC_dist.bin", "dist.bin", 8),
        "profit_table": ("RCC_pimat.bin", "profit_table.bin", 8),
        "macro_output": ("RCC_sim_Ysim.bin", "output.bin", 8),
        "macro_default": ("RCC_sim_dfsim.bin", "default.bin", 8),
        "macro_spread": ("RCC_sim_spreadsim.bin", "spread.bin", 8),
        "moments": ("RCC_GMMmodelmom.bin", "moments.bin", 8),
    }
    results = {
        label: compare(args.cpu / cpu, args.final_moments if label == "moments" and args.final_moments else args.gpu / gpu, size)
        for label, (cpu, gpu, size) in mapping.items()
    }
    for label, cpu, gpu, begin, end in [
        ("micro_investment", "RCC_sim_ifirm.bin", "investment.bin", 1, 58),
        ("micro_profit", "RCC_sim_pifirm.bin", "profit.bin", 1, 58),
        ("micro_forecast", "RCC_sim_pifcstfirm.bin", "forecast.bin", 2, 59),
    ]:
        indices = (firm * 302 + time for firm in range(1000) for time in range(begin, end))
        results[label] = compare(args.cpu / cpu, args.gpu / gpu, 8, indices)
    report = {
        "cpu_directory": str(args.cpu.resolve()),
        "gpu_directory": str(args.gpu.resolve()),
        "comparison": "exact bytes from independently written outputs",
        "moment_output": str(args.final_moments.resolve()) if args.final_moments else str((args.gpu / "moments.bin").resolve()),
        "arrays": results,
        "total_entries": sum(r["entries"] for r in results.values()),
        "total_bit_mismatches": sum(r["bit_mismatches"] for r in results.values()),
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["total_entries", "total_bit_mismatches"]}))
    return 2 if report["total_bit_mismatches"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
