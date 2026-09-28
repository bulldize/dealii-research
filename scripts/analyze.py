#!/usr/bin/env python3
"""Validate a completed run and plot its measured diagnostics; no new simulation."""
import csv
import json
import math
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(sys.argv[1])
summary = json.loads((folder / "summary.json").read_text())
rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader((folder / "history.csv").open())]
steps = rows[1:]
checks = {
    "completed": summary["completed"] is True and len(steps) == summary["steps"],
    "all_diagnostics_finite": all(math.isfinite(v) for r in rows for v in r.values()),
    "newton_increment": max(r["increment_inf"] for r in steps) < 1e-12,
    "equation_residual": max(r["residual_inf"] for r in steps) < 1e-9,
    "weak_incompressibility": max(r["weak_divergence_inf"] for r in steps) < 1e-10,
    "forced_energy_identity": max(abs(r["energy_balance_residual"]) for r in steps) < 1e-9,
    "positive_determinant_at_quadrature_points": min(r["min_det_F"] for r in rows) > 0,
}
(folder / "validation.json").write_text(json.dumps(checks, indent=2) + "\n")
if not all(checks.values()):
    raise RuntimeError(f"Validation failed: {checks}")
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
t = [r["time"] for r in steps]
for key, label in [("L2_velocity", "velocity"), ("L2_pressure", "pressure"), ("L2_F", "deformation F")]:
    axs[0, 0].semilogy(t, [r[key] for r in steps], marker=".", label=label)
axs[0, 0].set(title="Errors against manufactured solution", ylabel="Spatial L2 error")
axs[0, 0].legend()
axs[0, 1].plot([r["time"] for r in rows], [r["energy"] for r in rows], color="#236b8e")
axs[0, 1].set(title="Quadratic energy (forced problem)", ylabel="E(t)")
axs[1, 0].plot([r["time"] for r in rows], [r["min_det_F"] for r in rows], color="#198754")
axs[1, 0].set(title="Minimum det(F) at quadrature samples", ylabel="Sampled minimum")
for key, label in [("energy_balance_residual", "energy balance"), ("residual_inf", "equation residual"), ("weak_divergence_inf", "weak divergence")]:
    axs[1, 1].semilogy(t, [max(abs(r[key]), 1e-20) for r in steps], marker=".", label=label)
axs[1, 1].set(title="Algebraic consistency checks", ylabel="Absolute diagnostic")
axs[1, 1].legend()
for ax in axs.flat:
    ax.set_xlabel("Time")
    ax.grid(alpha=.2)
fig.suptitle("Giesekus (3.3): one manufactured-solution run | P2/P1/P1 | dt=0.005", fontsize=14)
fig.savefig(folder / "diagnostics.png", dpi=180)
fig.savefig(folder / "diagnostics.svg")
print(json.dumps({"checks": checks, "summary": summary}, indent=2))
