import sys
import psi4

if len(sys.argv) != 3:
    print("Usage: python scripts/psi4_reference.py input.xyz output_freq.dat")
    sys.exit(1)

xyz_file = sys.argv[1]
output_file = sys.argv[2]

with open(xyz_file, "r") as f:
    lines = f.readlines()

num_atoms = int(lines[0])
atom_lines = lines[2:2 + num_atoms]

geometry = "0 1\n"
for line in atom_lines:
    geometry += line

mol = psi4.geometry(geometry)

psi4.set_options({
    "basis": "sto-3g",
    "scf_type": "pk"
})

energy, wfn = psi4.frequency("hf/sto-3g", return_wfn=True)

freqs = wfn.frequency_analysis["omega"].data

with open(output_file, "w") as f:
    for freq in freqs:
        if freq > 0:
            f.write(f"{freq}\n")

print("Psi4 frequencies written to", output_file)
