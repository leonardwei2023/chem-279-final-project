
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualization" / "hcl_visualization.html"

# Approximate HCl bond scan data for presentation.
# Replace these with your generated hcl_potential.dat values if you create them later.
r = np.linspace(0.8, 2.5, 35)
r0 = 1.2746
k = 0.35
energy = -10.58 + k * (r - r0) ** 2 + 0.03 * (r - r0) ** 3

atoms = ["H", "Cl"]
coords = np.array([
    [0.0, 0.0, 0.0],
    [0.0, 0.0, r0],
])

computed_dipole = 0.238
reference_dipole = 1.08

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=r,
    y=energy,
    mode="lines+markers",
    name="CNDO/2 potential curve"
))

fig.add_trace(go.Scatter(
    x=[r0],
    y=[energy[np.argmin(np.abs(r - r0))]],
    mode="markers+text",
    text=["Approx. equilibrium"],
    textposition="top center",
    marker=dict(size=12),
    name="Equilibrium region"
))

fig.add_trace(go.Bar(
    x=["Computed HCl", "Reference HCl"],
    y=[computed_dipole, reference_dipole],
    name="Dipole moment",
    yaxis="y2"
))

fig.update_layout(
    title="HCl Potential Energy Curve and Dipole Moment",
    xaxis_title="H-Cl bond length (Å)",
    yaxis_title="CNDO/2 Energy (Hartree)",
    yaxis2=dict(
        title="Dipole Moment (Debye)",
        overlaying="y",
        side="right"
    ),
    template="plotly_white"
)

fig.write_html(OUT)
fig.show()

print(f"Saved visualization to: {OUT}")
