
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualization" / "h2_visualization.html"

atoms = ["H", "H"]
coords = np.array([
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.7414],
])

computed_freq = 4404.48
reference_freq = 4400.0

fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=coords[:, 0],
    y=coords[:, 1],
    z=coords[:, 2],
    mode="markers+text",
    text=atoms,
    textposition="top center",
    marker=dict(size=12),
    name="H2 atoms"
))

fig.add_trace(go.Scatter3d(
    x=coords[:, 0],
    y=coords[:, 1],
    z=coords[:, 2],
    mode="lines",
    line=dict(width=8),
    name="H-H bond"
))

fig.add_trace(go.Bar(
    x=["Computed H2", "Reference H2"],
    y=[computed_freq, reference_freq],
    name="Frequency",
    yaxis="y2"
))

fig.update_layout(
    title="H2 Vibrational Frequency and Molecular Structure",
    scene=dict(
        xaxis_title="x (Å)",
        yaxis_title="y (Å)",
        zaxis_title="z (Å)",
        aspectmode="data"
    ),
    yaxis2=dict(
        title="Frequency (cm⁻¹)",
        overlaying="y",
        side="right"
    ),
    template="plotly_white"
)

fig.write_html(OUT)
fig.show()

print(f"Saved visualization to: {OUT}")
