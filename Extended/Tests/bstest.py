"""Helpers for the tests of the Extended bow-shock code (see Tests/README.md).

Each test case is built in its own temporary directory: the Fortran sources are
copied there, some files are replaced by frozen copies (Tests/cases/), some
parameters are changed (patch_parameters), and the code is compiled with make and
run. The sources in Extended/ are never modified.
The results are read from the terminal output and from the files in output/.
"""

import os
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np

EXTENDED = pathlib.Path(__file__).resolve().parents[1]
ONE_ZONE = EXTENDED.parent / "One-Zone"
TESTS = pathlib.Path(__file__).resolve().parent

eV = 1.602176e-12          # as in constants.f90
mpc2 = 1836.0 * 0.511e6 * eV

FLOAT = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?)"


class TestFailure(Exception):
    pass


# ---------------------------------------------------------------------------
# Build and run
# ---------------------------------------------------------------------------

def patch_parameters(path, values):
    """Changes the value of Fortran parameters: 'name = old' -> 'name = new'.

    Only declarations with '::' are changed, and each name must appear exactly
    once in the file, so that a renamed parameter makes the test fail instead
    of being silently ignored."""
    text = path.read_text(errors="surrogateescape")
    for name, value in values.items():
        pattern = re.compile(
            r"^(?P<head>[ \t]*(?:real|integer|logical|character)[^!\n]*::[ \t]*"
            + re.escape(name) + r"[ \t]*=[ \t]*)(?P<value>[^!\n]*?)(?P<tail>[ \t]*(?:!.*)?)$",
            re.IGNORECASE | re.MULTILINE)
        n = len(pattern.findall(text))
        if n != 1:
            raise TestFailure(f"{path.name}: parameter '{name}' found {n} times (expected 1)")
        text = pattern.sub(lambda m: m.group("head") + str(value) + m.group("tail"), text)
    path.write_text(text, errors="surrogateescape")


def replace_text(path, pairs):
    """Replaces exact pieces of code [(old, new), ...] in a copy of a source file
    (only in the test build). Each 'old' must appear exactly once."""
    text = path.read_text(errors="surrogateescape")
    for old, new in pairs:
        n = text.count(old)
        if n != 1:
            raise TestFailure(f"{path.name}: text to replace found {n} times (expected 1):\n{old}")
        text = text.replace(old, new)
    path.write_text(text, errors="surrogateescape")


def merge(*dicts):
    """Merges {file: {parameter: value}} dictionaries (the last one wins)."""
    out = {}
    for d in dicts:
        for f, values in d.items():
            out.setdefault(f, {}).update(values)
    return out


def build(workdir, source=EXTENDED, files=None, patches=None, make_vars=None, replace=None):
    """Copies the code to workdir, replaces 'files' {name: path}, applies the
    patches {file: {parameter: value}} and the text replacements
    {file: [(old, new), ...]} to the copies, and compiles run.x. The sources in
    Extended/ are never modified."""
    workdir = pathlib.Path(workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    for f in source.glob("*.f90"):
        shutil.copy2(f, workdir / f.name)
    for name in ("Makefile", "Makefile.depend"):
        shutil.copy2(source / name, workdir / name)
    shutil.copy2(EXTENDED / "gauntff.dat", workdir / "gauntff.dat")
    for name, path in (files or {}).items():
        shutil.copy2(path, workdir / name)
    for filename, values in (patches or {}).items():
        patch_parameters(workdir / filename, values)
    for filename, pairs in (replace or {}).items():
        replace_text(workdir / filename, pairs)

    # LIBS is emptied: lapack and blas are not used by the code
    jobs = str(max(1, min(8, os.cpu_count() or 1)))
    command = ["make", "-j", jobs, "LIBS="] + [f"{k}={v}" for k, v in (make_vars or {}).items()] + ["run.x"]
    result = subprocess.run(command, cwd=workdir, capture_output=True, text=True, errors="replace")
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-30:])
        raise TestFailure(f"compilation failed in {workdir}:\n{tail}")
    return workdir


def run(workdir, threads=None, timeout=7200):
    """Runs ./run.x (with OMP_NUM_THREADS = threads, if given) and returns
    the terminal output (also saved in run.log)."""
    env = dict(os.environ)
    if threads is not None:
        env["OMP_NUM_THREADS"] = str(threads)
    result = subprocess.run(["./run.x"], cwd=workdir, capture_output=True, text=True,
                            errors="replace", env=env, timeout=timeout)
    (pathlib.Path(workdir) / "run.log").write_text(result.stdout + result.stderr)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-30:])
        raise TestFailure(f"run.x failed in {workdir} (exit code {result.returncode}):\n{tail}")
    return result.stdout


# ---------------------------------------------------------------------------
# Read the results
# ---------------------------------------------------------------------------

def section(stdout, start, end=None):
    """Terminal output from 'start' to 'end' (e.g. 'Reverse Shock', 'Forward shock')."""
    i = stdout.find(start)
    if i < 0:
        raise TestFailure(f"'{start}' not found in the terminal output")
    j = stdout.find(end, i + len(start)) if end else -1
    return stdout[i:j] if j > 0 else stdout[i:]


def value(text, label, index=0):
    """Number printed after 'label' (index: which occurrence; -1 = the last one)."""
    found = re.findall(re.escape(label) + r"\s*" + FLOAT, text)
    if not found:
        raise TestFailure(f"'{label}' not found in the terminal output")
    return float(found[index].replace("D", "E").replace("d", "e"))


def load(workdir, name):
    """Numeric table output/<name> (2D array). Fortran numbers such as
    1.0-100 (exponent without 'E') are also read."""
    path = pathlib.Path(workdir) / "output" / name
    if not path.exists():
        raise TestFailure(f"output file {name} not found")
    text = re.sub(r"(\d)([-+]\d{3})(?=\s|$)", r"\1E\2", path.read_text())
    return np.loadtxt(text.splitlines(), ndmin=2)


def first(workdir, name):
    """First number of output/<name> (files with a single total)."""
    return float(load(workdir, name)[0, 0])


def sed_luminosity(workdir, name):
    """int L_E dE of a SED file with columns E [eV], E*L_E [erg/s] (log grid)."""
    d = load(workdir, name)
    return float(np.sum(d[:, 1]) * np.log(d[1, 0] / d[0, 0]))


def at_frequency(table, nu):
    """Row of a table whose first column (frequency) is the closest to nu (in log)."""
    return table[np.argmin(np.abs(np.log(table[:, 0] / nu)))]


def read_dist(workdir, name, m_theta):
    """N(l, E) from diste_*.dat or distp_*.dat: E [eV] (n_E) and N [1/erg] (m_theta, n_E)."""
    d = load(workdir, name)
    n_E = d.shape[0] // m_theta
    return d[:n_E, 1], d[:, 2].reshape(m_theta, n_E)


def read_times(workdir, name, m_theta):
    """CD_tcoole/CD_tcoolp files: array (column, l, E)."""
    d = load(workdir, name)
    n_E = d.shape[0] // m_theta
    return d.T.reshape(d.shape[1], m_theta, n_E)


def power(E, N, t, n_emit):
    """sum_l sum_E N*E*dE/t for the whole BS [erg/s] (N for one linear emitter)."""
    Ee = E * eV
    dE = Ee * (E[1] / E[0] - 1.0)
    return float(np.sum(N * Ee * dE / t)) * n_emit


def outflow(E, N, t_cell, n_emit):
    """Energy flux of the particles advected out of the last cell [erg/s]."""
    Ee = E * eV
    dE = Ee * (E[1] / E[0] - 1.0)
    return float(np.sum(N[-1] * Ee * dE) / t_cell[-1]) * n_emit


def bad_numbers(workdir):
    """Output files with NaN or Infinity."""
    bad = []
    for f in sorted((pathlib.Path(workdir) / "output").glob("*.dat")):
        if re.search(r"nan|inf", f.read_text(errors="replace"), re.IGNORECASE):
            bad.append(f.name)
    return bad


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

class Checker:
    """Prints and collects the results of the checks of one test case."""

    def __init__(self, case):
        self.case = case
        self.failures = []
        self.n = 0

    def report(self, ok, message):
        self.n += 1
        print(f"  [{'ok' if ok else 'FAIL'}] {message}")
        sys.stdout.flush()
        if not ok:
            self.failures.append(f"{self.case}: {message}")
        return ok

    def true(self, condition, message):
        return self.report(bool(condition), message)

    def at_most(self, a, b, message, factor=1.0):
        """a <= factor*b"""
        lim = f"{factor:g} x " if factor != 1.0 else ""
        return self.report(a <= factor * b, f"{message}: {a:.3e} <= {lim}{b:.3e} (ratio {a / b:.3g})")

    def between(self, x, lo, hi, message):
        return self.report(lo <= x <= hi, f"{message}: {x:.3g} in [{lo:g}, {hi:g}]")

    def close(self, a, b, message, rtol):
        err = abs(a - b) / max(abs(b), 1e-300)
        return self.report(err <= rtol, f"{message}: {a:.6e} vs {b:.6e} (rel. diff. {err:.1e}, tol. {rtol:g})")

    def same_order(self, a, b, message, factor=10.0):
        """a and b agree within a factor 'factor'."""
        ratio = a / b if b != 0 else np.inf
        return self.report(1.0 / factor <= ratio <= factor,
                           f"{message}: {a:.3e} vs {b:.3e} (ratio {ratio:.3g}, max. factor {factor:g})")
