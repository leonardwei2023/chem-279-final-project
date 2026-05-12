"""
Script to run all visualizations
"""

from pathlib import Path
import subprocess

BASE_DIR = Path(__file__).resolve().parent

def list_files_no_ext(directory):
    return [p.stem for p in Path(directory).iterdir() if p.is_file()]


###
# Save all dipole vis to htmls
###
files = list_files_no_ext(f"{BASE_DIR}/dipoles/")
files = [f[:-7] for f in files]
print("Processing dipoles: ", files)

(BASE_DIR / "output/dipole").mkdir(parents=True, exist_ok=True)

for file in files:
    input_file = f"input/{file}.xyz"
    dipole_file = f"visualization/dipoles/{file}_dipole.json"
    output_file = f"{BASE_DIR}/output/dipole/{file}_dipole.html"

    subprocess.call(["python", f"{BASE_DIR}/dipole_vis.py", input_file, dipole_file, "-o", output_file])


###
# Save all freq vis to htmls
###
files = list_files_no_ext(f"{BASE_DIR}/modes/")
files = [f[:-6] for f in files]
print("Processing vibrational frequencies: ", files)

(BASE_DIR / "output/modes").mkdir(parents=True, exist_ok=True)

for file in files:
    input_file = f"visualization/modes/{file}_modes.json"
    output_file = f"{BASE_DIR}/output/modes/{file}_modes.html"

    subprocess.call(["python", f"{BASE_DIR}/freq_vis.py", input_file, "-o", output_file])