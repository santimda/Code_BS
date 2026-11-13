#!/usr/bin/env python3
"""
Timescales of electrons or protons at a few positions of the shock.
One figure per shock and angle: t{e,p}_{RS,FS}_th{theta}. The cell of each
angle is taken from termonorm_{shock}.dat (row l).
"""
import numpy as np

from bsplot import style
from bsplot.common import Run, base_parser
from bsplot.io import TCOOLE_COLS, TCOOLP_COLS

cb = style.cb_colors

# column: label, line style, colour, width
CURVES = {
    "e": {"t_acc": (r"$t_{\rm acc}$", "--", cb["black"], 2.), "t_diff": (r"$t_{\rm diff}$", "-.", cb["blue"], 2.),
          "t_syn": (r"$t_{\rm syn}$", "-", cb["green"], 2.), "t_adv": (r"$t_{\rm adv}$", "-.", "cyan", 2.),
          "t_IC_star": (None, "-", "purple", 1.), "t_IC_IR": (None, "-", "orange", 1.),
          "t_IC": (r"$t_{\rm IC}$", "-", cb["red"], 2.), "t_tot": (r"$t_{\rm tot}$", "-", "magenta", 2.),
          "t_adi": (r"$t_{\rm adi}$", ":", "grey", 2.)},
    "p": {"t_acc": (r"$t_{\rm acc}$", "--", cb["black"], 2.), "t_diff": (r"$t_{\rm diff}$", "-.", cb["blue"], 2.),
          "t_pp": (r"$t_{pp}$", "-", "tomato", 2.), "t_adv": (r"$t_{\rm adv}$", "-.", "cyan", 2.),
          "t_adi": (r"$t_{\rm adi}$", ":", "grey", 2.), "t_tot": (r"$t_{\rm tot}$", "-", "magenta", 2.)},
}
COLS = {"e": TCOOLE_COLS, "p": TCOOLP_COLS}


def main(argv=None):
    p = base_parser(__doc__)
    p.add_argument("--species", choices=["e", "p"], default="e")
    p.add_argument("--thetas", type=float, nargs="+", default=[10., 45., 120.], help="[deg]")
    run = Run(p, argv)
    sp = run.args.species

    for shock in run.shocks_nt():
        d, termo = run.out.read(f"tcool{sp}", shock=shock), run.out.read("termonorm", shock=shock)
        if d is None or termo is None:
            continue
        col = {name: d[:, i] for i, name in enumerate(COLS[sp])}
        if sp == "e":
            col["t_IC"] = 1. / (1. / col["t_IC_star"] + 1. / col["t_IC_IR"])
        theta = np.degrees(termo[:, 0])

        for th0 in run.args.thetas:
            l = np.argmin(np.abs(theta - th0)) # row l of termonorm = cell l + 1
            m = col["l"] == l + 1
            E = col["E"][m]

            fig, ax = style.new_figure()
            logt = []
            for name, (label, ls, color, lw) in CURVES[sp].items():
                ax.plot(E, col[name][m], label=label, ls=ls, color=color, lw=lw, alpha=0.8)
                t = col[name][m]
                logt.append(np.log10(t[(t > 0) & np.isfinite(t)]))
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(E.min(), E.max())
            logt = np.concatenate(logt) # y range: 1-99 percentiles
            ax.set_ylim(10**np.floor(np.percentile(logt, 1)), 10**np.ceil(np.percentile(logt, 99)))
            style.finish_axes(ax)
            ax.set_xlabel(rf"$E_{sp}$ [eV]")
            ax.set_ylabel(r"$t$ [s]")
            ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.), ncol=4,
                      handlelength=1.5, columnspacing=1.)
            style.system_label(ax, f"{run.label}\n{shock}, " + rf"$\theta = {theta[l]:.0f}^\circ$",
                               box=True)
            run.save(fig, f"t{sp}_{shock}_th{theta[l]:.0f}")


if __name__ == "__main__":
    main()
