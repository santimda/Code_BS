#!/usr/bin/env python3
"""Regression test for the numerical summary printed by the BD+43 run (Tests/cases/BD43).

The code is copied to a temporary directory, with the frozen system_parameters.f90
and controles.f90 of Tests/cases/BD43, the grid of global.f90 fixed (GRID) and a
radiative FS (system_parameters.f90), and the values printed in the terminal are
compared with reference_terminal_output.json.
The sources in Extended/ are not modified.

  python3 Tests/test_terminal_output.py                     # compare
  python3 Tests/test_terminal_output.py --update-reference  # rewrite the reference
"""

import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]
REFERENCE = pathlib.Path(__file__).with_name("reference_terminal_output.json")
CONFIGURATION = pathlib.Path(__file__).parent / "cases" / "BD43"

# global.f90 also contains code: only the grid is fixed, so that the reference
# does not depend on the resolution used in Extended/
GRID = {"m_theta": "300", "mE_e": "300", "mE_p": "300", "mE_f": "500", "mE_ph": "500", "mphi": "30"}
# system_parameters.f90: the reference is for a radiative FS
RADIATIVE_FS = ".true."

NUMBER = r"([0-9.]+E[+-]?[0-9]+)"
PATTERNS = {
    "injected_power": (r"L_inj_w =\s+" + NUMBER + r".*?L_inj\(RS\) =\s+" + NUMBER
                       + r".*?L_inj\(FS\) =\s+" + NUMBER),
    "magnetic_field_apex": r"\n\s+B0 =\s*" + NUMBER,
    "temperature_apex": r"\n\s+T0 =\s*" + NUMBER,
    "density_apex": r"n_a0 =\s*" + NUMBER,
    "free_free_RS": r"f-f emission from the RS:.*?L_ion_tot =\s+" + NUMBER,
    "synchrotron_unabsorbed": r"L_syn\(unabsorbed\)=\s+" + NUMBER,
    "bremsstrahlung": r"L_Brem=\s+" + NUMBER,
    "ic_dust": r"L_IC_IR =\s+" + NUMBER,
    "ic_stellar": r"L_IC_st =\s+" + NUMBER + r"\s*\n\s*Luminosity between",
    "x_rays": r"L_X =\s+" + NUMBER,
    "proton_proton": r"L_pp=\s+" + NUMBER,
    "free_free_hot_layer": r"\n\s*L_HL_tot =\s+" + NUMBER,
    "free_free_isothermal_layer": r"isothermal layer:.*?L_ion_tot =\s+" + NUMBER,
    "maps_synchrotron": r"Pre-conv:  S_sync\(\w+\) =\s+" + NUMBER,
    "maps_convolved": r"Post-conv: S\(\w+\)\s+=\s+" + NUMBER,
}
# synchrotron fluxes of the RS at the 7 printed frequencies
PATTERNS["synchrotron_radio_RS"] = r"ELECTRONS:\s*L_syn\(unabsorbed\)=\s+[0-9.E+-]+" + \
    "".join(r"\s*S_\w+\[mJy\]=\s+" + NUMBER for _ in range(7))


def values(pattern, output):
    matches = re.findall(pattern, output, flags=re.DOTALL)
    return [float(value) for match in matches for value in (match if isinstance(match, tuple) else (match,))]


def run_checked(command, cwd):
    try:
        return subprocess.run(command, cwd=cwd, check=True, capture_output=True, text=True, errors="replace")
    except subprocess.CalledProcessError as error:
        print(f"Command failed with exit status {error.returncode}: {' '.join(error.cmd)}", file=sys.stderr)
        if error.stdout:
            print("--- stdout ---", file=sys.stderr)
            print(error.stdout[-5000:], file=sys.stderr, end="")
        if error.stderr:
            print("--- stderr ---", file=sys.stderr)
            print(error.stderr[-5000:], file=sys.stderr, end="")
        raise


def fix_grid(path):
    text = path.read_text()
    for name, value in GRID.items():
        pattern = re.compile(r"^(\s*integer\s*,\s*parameter\s*::\s*" + name + r"\s*=\s*)\d+", re.MULTILINE)
        text, n = pattern.subn(lambda m: m.group(1) + value, text)
        if n != 1:
            raise SystemExit(f"global.f90: parameter {name} found {n} times (expected 1)")
    path.write_text(text)


def fix_radiative_fs(path):
    text = path.read_text()
    pattern = re.compile(r"^(\s*logical\s*,\s*parameter\s*::\s*radiative_FS\s*=\s*)\.\w+\.", re.MULTILINE | re.IGNORECASE)
    text, n = pattern.subn(lambda m: m.group(1) + RADIATIVE_FS, text)
    if n != 1:
        raise SystemExit(f"system_parameters.f90: radiative_FS found {n} times (expected 1)")
    path.write_text(text)


def main():
    update = "--update-reference" in sys.argv[1:]
    with tempfile.TemporaryDirectory(prefix="code-bs-regression-") as temporary_directory:
        test_root = pathlib.Path(temporary_directory)
        for source in ROOT.glob("*.f90"):
            shutil.copy2(source, test_root / source.name)
        for filename in ("Makefile", "Makefile.depend"):
            shutil.copy2(ROOT / filename, test_root / filename)
        shutil.copy2(ROOT / "gauntff.dat", test_root / "gauntff.dat")
        for filename in ("controles.f90", "system_parameters.f90"):
            shutil.copy2(CONFIGURATION / filename, test_root / filename)
        fix_grid(test_root / "global.f90")
        fix_radiative_fs(test_root / "system_parameters.f90")

        jobs = str(max(1, min(8, os.cpu_count() or 1)))
        run_checked(["make", "-j", jobs, "LIBS=", "run.x"], test_root)   # lapack/blas are not used
        result = run_checked(["./run.x"], test_root)

    found = {name: values(pattern, result.stdout) for name, pattern in PATTERNS.items()}
    if update:
        lines = [f'  "{name}": {json.dumps(found[name])}' for name in PATTERNS]
        REFERENCE.write_text("{\n" + ",\n".join(lines) + "\n}\n")
        print(f"Reference rewritten: {REFERENCE}")
        return 0

    reference = json.loads(REFERENCE.read_text())
    failures = []
    for name, expected in reference.items():
        actual = found.get(name, [])
        if len(actual) != len(expected):
            failures.append(f"{name}: expected {len(expected)} values, found {len(actual)}")
            continue
        for index, (actual_value, expected_value) in enumerate(zip(actual, expected), start=1):
            if not math.isclose(actual_value, expected_value, rel_tol=2e-3, abs_tol=1e-30):
                failures.append(f"{name}[{index}]: expected {expected_value:g}, found {actual_value:g}")

    if failures:
        print("Terminal-output regression failed:")
        print("\n".join(f"- {failure}" for failure in failures))
        return 1

    print("Terminal-output regression passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())