"""Conditional measurement-noise minimization from a quadratic moment basis."""
from __future__ import annotations

import math
from pathlib import Path
import re


def fortran_float(value: str) -> float:
    value = value.replace("D", "E").replace("d", "e")
    if "e" not in value.lower():
        value = re.sub(r"(?<=\d)([+-]\d+)$", r"E\1", value)
    return float(value)


def tagged_values(stdout: str, tag: str, count: int) -> list[float]:
    matches = re.findall(rf"^{re.escape(tag)}=\s*(.+)$", stdout, re.MULTILINE)
    if len(matches) != 1:
        raise ValueError(f"Expected one {tag} output")
    values = [fortran_float(x) for x in matches[0].split()]
    if len(values) != count:
        raise ValueError(f"Expected {count} values for {tag}")
    return values


def minimize_basis(m0, c1, c2, target, weight, bounds, selected=range(1, 17)) -> dict:
    """Minimize Q(s) on the supplied interval; selected indices are zero based."""
    import numpy as np

    vectors = [np.asarray(x, dtype=float) for x in (m0, c1, c2, target)]
    m0, c1, c2, target = vectors
    weight = np.asarray(weight, dtype=float)
    selected = list(selected)
    if any(v.shape != (17,) for v in vectors) or weight.shape != (17, 17):
        raise ValueError("Expected 17 moments and a 17 by 17 weight matrix")
    if (not selected or len(set(selected)) != len(selected)
            or any(type(i) is not int or not 0 <= i < 17 for i in selected)):
        raise ValueError("Invalid selected moment indices")
    w = weight[np.ix_(selected, selected)]
    if not all(np.isfinite(v[selected]).all() for v in vectors) or not np.isfinite(w).all():
        raise ValueError("Objective inputs must be finite")
    lower, upper = map(float, bounds)
    if not (math.isfinite(lower) and math.isfinite(upper) and 0 <= lower < upper):
        raise ValueError("Use finite sigma bounds with 0 <= lower < upper")

    # Express sigma = lower + width*t to solve the polynomial on [0, 1].
    width = upper - lower
    residual = np.stack((
        target[selected] - (m0[selected] + lower*c1[selected] + lower*lower*c2[selected]),
        -width*(c1[selected] + 2*lower*c2[selected]),
        -width*width*c2[selected],
    ))
    coefficients = np.zeros(5)
    for i in range(3):
        for j in range(3):
            coefficients[i+j] += residual[i] @ w @ residual[j]
    if not np.isfinite(coefficients).all():
        raise ValueError("Nonfinite objective coefficients")
    derivative = np.arange(1, 5)*coefficients[1:]
    derivative = np.trim_zeros(derivative, "b")
    candidates = [lower, upper]
    if len(derivative) > 1:
        derivative = derivative / np.max(np.abs(derivative))
        for root in np.polynomial.polynomial.polyroots(derivative):
            if abs(root.imag) <= 1e-9 and -1e-12 <= root.real <= 1+1e-12:
                t = min(1.0, max(0.0, float(root.real)))
                candidates.append(min(upper, max(lower, lower + width*t)))
    rows = []
    for sigma in sorted(set(candidates)):
        moments = m0 + sigma*c1 + sigma*sigma*c2
        diff = target[selected] - moments[selected]
        q = float(diff @ w @ diff)
        if math.isfinite(q):
            rows.append({"sigma_pi": sigma, "objective": q,
                         "moments": [float(x) if math.isfinite(x) else None for x in moments]})
    if not rows:
        raise ValueError("No finite conditional objective value")
    best = min(rows, key=lambda row: (row["objective"], row["sigma_pi"]))
    return {**best, "sigma_pi_bounds": [lower, upper],
            "candidate_sigmas": [row["sigma_pi"] for row in rows],
            "candidate_objectives": [row["objective"] for row in rows],
            "selected_moments": [i+1 for i in selected],
            "selection": "interval endpoints and real stationary points"}


def from_stdout(stdout: str, inputs: Path, bounds) -> dict:
    import numpy as np

    basis = {key: tagged_values(stdout, tag, 17) for key, tag in (
        ("m0", "RCC_SIGMAPI_M0"), ("c1", "RCC_SIGMAPI_C1"), ("c2", "RCC_SIGMAPI_C2"))}
    target = np.loadtxt(inputs / "_GMMdatamom.txt", delimiter=",")
    weight = np.loadtxt(inputs / "_GMMdatawgt.txt", delimiter=",")
    return {**minimize_basis(**basis, target=target, weight=weight, bounds=bounds),
            "basis": {key: [x if math.isfinite(x) else None for x in values]
                      for key, values in basis.items()}}
