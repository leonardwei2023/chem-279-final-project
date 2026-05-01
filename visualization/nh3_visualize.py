
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualization" / "nh3_visualization.html"

atoms = ["N", "H", "H", "H"]
coords = np.array([
    [0.0000, 0.0000, 0.0000],
    [0.0000, 0.9377, 0.3816],
    [0.8121, -0.4688, 0.3816],
    [-0.8121, -0.4688, 0.3816],
])

computed_dipole = 0.533
reference_dipole = 1.47

fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=coords[:, 0],
    y=coords[:, 1],
    z=coords[:, 2],
    mode="markers+text",
    text=atoms,
    textposition="top center",
    marker=dict(size=[16, 10, 10, 10]),
    name="NH3 atoms"
))

for i in [1, 2, 3]:
    fig.add_trace(go.Scatter3d(
        x=[coords[0, 0], coords[i, 0]],
        y=[coords[0, 1], coords[i, 1]],
        z=[coords[0, 2], coords[i, 2]],
        mode="lines",
        line=dict(width=7),
        name="N-H bond"
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
    x=["Computed NH3", "Reference NH3"],
    y=[computed_dipole, reference_dipole],
    name="Dipole moment",
    yaxis="y2"
))

fig.update_layout(
    title="NH3 3D Structure and Dipole Moment",
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
