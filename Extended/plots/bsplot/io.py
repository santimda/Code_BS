"""
Output files of the code (output/). All file names are in FILES.
"""
import os

import numpy as np

FILES = {
    # CD_termo
    "termonorm": "termonorm_{shock}.dat",
    "rad_adi": "rad_adi_{sys}_{shock}.dat",
    "vs_RH": "new_vs_RH_{shock}.dat",
    # CD_diste, CD_distp: N of one linear emitter; totals for one quadrant (x n_emit/4)
    "diste": "diste_{sys}_{shock}.dat",
    "diste_total": "diste_{shock}_{sys}_total.dat",
    "distp": "distp_{sys}_{shock}.dat",
    "distp_total": "distp_{shock}_total.dat",
    "tcoole": "CD_tcoole_{sys}_{shock}.dat",
    "tcoolp": "CD_tcoolp_{sys}_{shock}.dat",
    # Whole shock SEDs: E [eV], E L_E [erg/s]
    "radsyn": "radsyn_{sys}_{shock}_tot.dat",
    "radic_st": "radic_st_{sys}_{shock}_tot.dat",
    "radic_ir": "radic_ir_{sys}_{shock}_tot.dat",
    "radbr": "radbr_{sys}_{shock}_tot.dat",
    "radpp": "radpp_{sys}_{shock}_tot.dat",
    "th_ion": "ionized_ff_{sys}_{shock}.dat", # RS: shocked wind; FS: isothermal layer
    "th_HL": "{shock}_ff_{sys}.dat", # hot layer (radiative FS)
    # Radio: nu [Hz], S [mJy]
    "S_syn": "S_syn_{sys}_{shock}.dat",
    "S_syn_abs": "S_syn_abs_{sys}_{shock}.dat", # nu, S_unabs, S_abs
    "S_ff": "S_ff_{sys}_{shock}.dat", # hot layer
    "S_ionized": "S_ionized_{sys}_{shock}.dat", # RS: shocked wind; FS: isothermal layer
    # Maps: x, y ["], value
    "map": "SW_conv_{label}abs_{sys}.dat",
    "R0_proj": "R0_proj_{sys}.dat",
    "halpha": "S_alpha_convmap_{sys}.dat",
}

# Columns of CD_tcoole and CD_tcoolp
TCOOLE_COLS = ["l", "E", "t_IC", "t_syn", "t_br", "t_diff", "t_acc", "t_cool", "t_adv",
               "t_cell", "t_IC_star", "t_IC_IR", "t_esc", "t_tot", "t_adi", "t_coul", "t_ion"]
TCOOLP_COLS = ["l", "E", "t_acc", "t_pp", "t_adv", "t_diff", "t_cool", "t_cell",
               "t_esc", "t_tot", "t_adi"]


class Output:

    def __init__(self, directory, system):
        self.dir = directory
        self.system = system

    def path(self, key, **fmt):
        return os.path.join(self.dir, FILES[key].format(sys=self.system, **fmt))

    def exists(self, key, **fmt):
        return os.path.isfile(self.path(key, **fmt))

    def read(self, key, **fmt):
        # numeric columns, or None if the file is missing
        if not self.exists(key, **fmt):
            print(f"  missing {self.path(key, **fmt)}")
            return None
        return np.atleast_2d(np.loadtxt(self.path(key, **fmt)))


def load_gridded(path):
    # x, y, value columns -> X, Y, Z on the regular grid (as gnuplot splot)
    x, y, z = np.loadtxt(path, unpack=True)
    xu, yu = np.unique(x), np.unique(y)
    Z = np.full((yu.size, xu.size), np.nan)
    Z[np.searchsorted(yu, y), np.searchsorted(xu, x)] = z
    X, Y = np.meshgrid(xu, yu)
    return X, Y, Z


def log10(x):
    # log10 with NaN where x <= 0 (the curve is cut there)
    x = np.asarray(x, dtype=float)
    return np.log10(np.where(x > 0, x, np.nan))
