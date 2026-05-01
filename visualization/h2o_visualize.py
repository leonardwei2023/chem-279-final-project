
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualization" / "h2o_visualization.html"

atoms = ["O", "H", "H"]
coords = np.array([
    [0.0000, 0.0000, 0.0000],
    [0.0000, 0.7570, 0.5860],
    [0.0000, -0.7570, 0.5860],
])

computed_dipole = 0.450
reference_dipole = 1.85

fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=coords[:, 0],
    y=coords[:, 1],
    z=coords[:, 2],
    mode="markers+text",
    text=atoms,
    textposition="top center",
    marker=dict(size=[16, 10, 10]),
    name="H2O atoms"
))

for i in [1, 2]:
    fig.add_trace(go.Scatter3d(
        x=[coords[0, 0], coords[i, 0]],
        y=[coords[0, 1], coords[i, 1]],
        z=[coords[0, 2], coords[i, 2]],
        mode="lines",
        line=dict(width=7),
        name="O-H bond"
    ))

fig.add_trace(go.Cone(
    x=[0],
    y=[0],
    z=[0.2],
    u=[0],
    v=[0],
    w=[-0.8],
    sizemode="absolute",
    sizeref=0.25,
    name="Dipole direction"
))

fig.add_trace(go.Bar(
    x=["Computed H2O", "Reference H2O"],
    y=[computed_dipole, reference_dipole],
    name="Dipole moment",
    yaxis="y2"
))

fig.update_layout(
    title="H2O 3D Structure and Dipole Moment",
    scene=dict(
        xaxis_title="x (Å)",
        yaxis_title="y (Å)",
        zaxis_title="z (Å)",
        aspectmode="data"
    ),
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
