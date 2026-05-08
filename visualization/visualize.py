import json
import argparse
import numpy as np
import py3Dmol

def build_frames(mode, base, atoms, n_frames=30, amplitude=0.5):

    base = np.array(base).reshape(-1, 3)
    disp = np.array([[a["dx"], a["dy"], a["dz"]] for a in mode["displacements"]])

    frames = []

    for t in np.linspace(0, 2*np.pi, n_frames):

        coords = base + amplitude * np.sin(t) * disp

        xyz = f"{len(atoms)}\nmode\n"
        for i, a in enumerate(atoms):
            x, y, z = coords[i]
            xyz += f"{a} {x:.6f} {y:.6f} {z:.6f}\n"

        frames.append(xyz)

    return frames

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Visualize vibrational modes from JSON"
    )

    

    parser.add_argument(
        "input",
        help="Input modes JSON file"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="all_modes.html",
        help="Output HTML file"
    )

    args = parser.parse_args()

    with open(args.input) as f:
        data = json.load(f)

    atoms = data["atoms"]
    base = data["base_coords"]
    modes = data["modes"]

    animations = []

    for mode in modes:
        frames = build_frames(mode, base, atoms, amplitude=0.3)
        animations.append(frames)

    views = []

    for i, frames in enumerate(animations):

        view = py3Dmol.view(width=300, height=300)

        view.addModelsAsFrames("".join(frames), "xyz")

        view.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
        view.animate({"loop": "forward"})
        view.zoomTo()

        view.rotate(-90, "x")
        view.rotate(0, "y")
        view.rotate(0, "z")

        views.append(view)

    html = "<html><body style='display:flex;flex-wrap:wrap;'>"

    for i, view in enumerate(views):
        html += f"<div style='margin:10px'>"
        html += view._make_html()
        html += "</div>"

    html += "</body></html>"

    with open(args.output, "w") as f:
        f.write(html)

    print(f"Wrote {args.output}")