"""
Command line options and run information shared by all the plotting scripts.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg") # files only, no windows
import matplotlib.pyplot as plt
import numpy as np

from . import style
from .fortran import read_controls, read_distance, read_parameter, read_string
from .io import Output

PLOTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_DIR = os.path.dirname(PLOTS_DIR) # Extended/


def base_parser(description):
    p = argparse.ArgumentParser(description=description,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", default=os.path.join(BASE_DIR, "output"), help="output files of the code")
    p.add_argument("--controls", default=os.path.join(BASE_DIR, "controles.f90"))
    p.add_argument("--figdir", default=None, help="figures folder (default: figures/<system>/)")
    p.add_argument("--system", default=None, help="block of plot_config.SYSTEMS (default: system)")
    p.add_argument("--label", default=None, help="system label in the plots")
    p.add_argument("--fontsize", type=float, default=style.FONTSIZE, help="[pt]")
    p.add_argument("--font", choices=["times", "cm"], default="times")
    p.add_argument("--no-usetex", action="store_true")
    p.add_argument("--fmt", nargs="+", default=["pdf"], help="output formats")
    p.add_argument("--jobs", nargs="+", default=None, help="only these plots (make_plots.py)")
    p.add_argument("--dry-run", action="store_true", help="only list the plots (make_plots.py)")
    return p


def default_figdir(system):
    # figures/<system>/, with only safe characters in the folder name
    name = "".join(ch if ch.isalnum() or ch in "+-_." else "_" for ch in system.strip())
    return os.path.join(BASE_DIR, "figures", name)


class Run:
    """Everything a plotting script needs to know about the run."""

    def __init__(self, parser, argv=None):
        a = self.args = parser.parse_args(argv)
        style.set_style(a.fontsize, a.font, not a.no_usetex)

        code_dir = os.path.dirname(os.path.abspath(a.controls))
        self.sources = [a.controls] + [os.path.join(code_dir, f)
                                       for f in ("system_parameters.f90", "global.f90")]
        self.controls = read_controls(a.controls)
        main_f90 = os.path.join(code_dir, "main.f90")
        self.radiative_fs = os.path.isfile(main_f90) and read_controls(main_f90).get("radiative_fs", False)

        # Plot configuration of the system (plot_config.py)
        sys.path.insert(0, PLOTS_DIR)
        import plot_config as cfg
        self.system = read_string(self.sources, "system", "")
        key = a.system or self.system
        self.conf = {**cfg.DEFAULTS, **cfg.SYSTEMS.get(key, {})}
        self.conf["distance"] = self.conf["distance"] or read_distance(self.sources)
        self.label = a.label if a.label is not None else self.conf["label"] or key
        self.obs_dir, self.sens_dir = cfg.OBS_DIR, cfg.SENS_DIR

        self.out = Output(a.data, self.system)
        self.figdir = a.figdir or default_figdir(self.system)
        os.makedirs(self.figdir, exist_ok=True)

    def flag(self, name):
        return self.controls.get(name.lower(), False)

    def shocks(self):
        return [s for s in ("RS", "FS") if self.flag(s)]

    def shocks_nt(self):
        # shocks with non-thermal particles (none in a radiative FS)
        return [s for s in self.shocks() if not (s == "FS" and self.radiative_fs)]

    def parameter(self, name, default=None):
        return read_parameter(self.sources, name, default)

    def save(self, fig, name):
        for fmt in self.args.fmt:
            path = os.path.join(self.figdir, f"{name}.{fmt}")
            fig.savefig(path)
            print(f"  -> {os.path.relpath(path)}")
        plt.close(fig)


def int_ticks(ax, axis, lo, hi, max_major=6):
    # integer ticks, major ones every 'step' so that there are at most max_major
    ticks = range(int(np.ceil(lo)), int(np.floor(hi)) + 1)
    step = max(1, int(np.ceil(len(ticks) / max_major)))
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_ticks([t for t in ticks if t % step == 0])
    a.set_ticks([t for t in ticks if t % step], minor=True)
