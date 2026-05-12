import json
import argparse
import numpy as np
import py3Dmol

def read_xyz(filename):
    atoms = []

    with open(filename, "r") as f:
        lines = f.readlines()

    num_atoms = int(lines[0].strip())

    for line in lines[2:2 + num_atoms]:
        parts = line.split()
        if len(parts) < 4:
            continue

        symbol = parts[0]
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])

        atoms.append({
            "elem": symbol,
            "x": x,
            "y": y,
            "z": z
        })

    return atoms


def show_dipole(dipole_vector, dipole_magnitude, xyz, scale=100, output_file="molecule_dipole.html"):
    # Create viewer
    viewer = py3Dmol.view(width=700, height=500)
    viewer.addModel(xyz, "xyz")

    viewer.setStyle({
        "stick": {"radius": 0.15},
        "sphere": {"scale": 0.3}
    })

    # Dipole visualization
    origin = np.array([0.0, 0.0, 0.0])
    # scale = 100.0  # visual scaling factor for dipole length
    # end = origin + scale * dipole_vector
    end = origin + np.linalg.norm(dipole_vector) * scale

    viewer.addArrow({
        "start": {
            "x": float(origin[0]),
            "y": float(origin[1]),
            "z": float(origin[2]),
        },
        "end": {
            "x": float(end[0]),
            "y": float(end[1]),
            "z": float(end[2]),
        },
        "radius": 0.08,
        "color": "red",
        "mid": 0.7
    })

    # Label showing magnitude
    viewer.addLabel(
        f"|μ| = {dipole_magnitude:.3f} D",
        {
            "position": {
                "x": float(end[0]+1*scale),
                "y": float(end[1]+1*scale),
                "z": float(end[2]+1*scale),
            },
            "backgroundColor": "white",
            "fontColor": "black"
        }
    )

    viewer.setBackgroundColor("0xeeeeee")
    viewer.zoomTo()
    viewer.spin(True)

    # output_file = f"visualization/output/{output_file}_dipole.html"
    html = viewer._make_html()

    with open(output_file, "w") as f:
        f.write(html)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Visualize dipoles modes from JSON"
    )

    parser.add_argument(
        "xyz",
        help="Input molecule xyz file"
    )

    parser.add_argument(
        "dipole",
        help="Input dipole JSON file"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="dipole.html",
        help="Output HTML file"
    )

    args = parser.parse_args()

    atoms = read_xyz(args.xyz)
    dipole_file = args.dipole

    xyz = f"{len(atoms)}\nfrom xyz file\n"
    for a in atoms:
        xyz += f"{a['elem']} {a['x']} {a['y']} {a['z']}\n"

    with open(dipole_file, "r") as f:
        data = json.load(f)

    mu = data["mu"]
    dipole_vector = np.array([mu["x"], mu["y"], mu["z"]])

    show_dipole(dipole_vector, data["magnitude"], xyz, scale=1, output_file=args.output)

    print(f"Wrote {args.output}")