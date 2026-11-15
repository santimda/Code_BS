#!/usr/bin/env python3
"""
Broadband SED, log10(eps L_eps) vs log10(eps): non-thermal processes of each
shock (RS dotted, FS dashed), free-free of the shocks, total, instrument
sensitivities and upper limits (plot_config.py).
"""
import os

import numpy as np
from matplotlib.lines import Line2D

from bsplot import style
from bsplot.common import Run, base_parser, int_ticks
from bsplot.io import log10

EV, H, JY = 1.602176e-12, 6.6260755e-27, 1.e-23
TO_EV = {"eV": 1., "keV": 1.e3, "MeV": 1.e6, "GeV": 1.e9}

# file, control, particles, label, colour
PROCESSES = [
    ("radsyn", "radsyn_sw_calc", "electrons", "syn", "green"),
    ("radic_ir", "radic_IR_calc", "electrons", r"IC$_{\rm IR}$", "orange"),
    ("radic_st", "radic_st_calc", "electrons", r"IC$_*$", "purple"),
    ("radpp", "radpp_calc", "protons", r"$pp$", "tomato"),
    ("radbr", "radbr_sw_calc", "electrons", r"$br$", "royalblue"),
]
LS = {"RS": ":", "FS": "--"}


def sensitivity(run, s):
    # sensitivity curve -> log10 eps [eV], log10 eps L_eps [erg/s]
    path = os.path.join(run.sens_dir, s["file"])
    if not os.path.isfile(path):
        print(f"  missing {path}")
        return None
    x, y = np.loadtxt(path, skiprows=s.get("skiprows", 0), usecols=(0, 1), unpack=True)
    y = y * s.get("factor", 1.)
    if s["x"] == "GHz":
        nu = 1.e9 * x
        E = H * nu / EV
    else:
        E = x * TO_EV[s["x"]]
        nu = E * EV / H
    F = 1.e-6 * JY * y * nu if s["y"] == "uJy" else y # erg cm^-2 s^-1
    return log10(E), log10(4. * np.pi * run.conf["distance"]**2 * F)


def main(argv=None):
    run = Run(base_parser(__doc__), argv)
    fig, ax = style.new_figure()
    curves, handles = [], []

    def add(E, L, label, color, ls="-"):
        curves.append((E, L))
        ax.plot(log10(E), log10(L), color=color, ls=ls, lw=2, alpha=0.8)
        if label not in [h.get_label() for h in handles]:
            handles.append(Line2D([], [], color=color, lw=2, alpha=0.8, label=label))

    #-----------------------------------------------------------
    # Non-thermal processes (only one shock: solid lines)

    shocks = [s for s in run.shocks_nt() if any(run.out.exists(p[0], shock=s) for p in PROCESSES)]
    ls = LS if len(shocks) > 1 else {s: "-" for s in shocks}
    for key, flag, particles, label, color in PROCESSES:
        if run.flag(particles) and run.flag(flag):
            for s in shocks:
                if run.out.exists(key, shock=s):
                    d = run.out.read(key, shock=s)
                    add(d[:, 0], d[:, 1], label, color, ls[s])

    #-----------------------------------------------------------
    # Free-free of the shocks: RS, and hot and isothermal layers of the FS
    if run.flag("xrad_calc"):
        keys = {"RS": ("th_ion",),
                "FS": ("th_ion", "th_HL") if run.radiative_fs else ("th_ion",)}
        th = [run.out.read(k, shock=s) for s in run.shocks() for k in keys[s]
              if run.out.exists(k, shock=s)]
        if th:
            E = th[0][:, 0]
            add(E, sum(np.interp(E, d[:, 0], d[:, 1], left=0., right=0.) for d in th),
                "free-free", "saddlebrown")

    if not curves:
        print("  nothing to plot")
        return

    # Total, on a common grid
    E = np.geomspace(min(c[0][c[0] > 0].min() for c in curves), max(c[0].max() for c in curves), 400)
    L = sum(np.interp(E, c[0], c[1], left=0., right=0.) for c in curves)
    ax.plot(log10(E), log10(L), color="black", lw=2, alpha=0.8)
    handles.append(Line2D([], [], color="black", lw=2, alpha=0.8, label="total"))

    #-----------------------------------------------------------
    # Sensitivities and upper limits. The upper limits are band fluxes F; they
    # are drawn at <eps L_eps> = 4 pi d^2 F / ln(E2/E1) across the band.

    for s in run.conf["sensitivities"]:
        c = sensitivity(run, s)
        if c is not None:
            ax.plot(*c, color="black", lw=1., alpha=0.8)
            if s.get("text"):
                ax.text(*s["text"], s["label"], fontsize=0.85 * run.args.fontsize)

    ul_top = []
    for u in run.conf["sed_ul"]:
        E1, E2 = u["E"]
        y = np.log10(4. * np.pi * run.conf["distance"]**2 * u["F"] / np.log(E2 / E1))
        x0, x1 = np.log10(E1), np.log10(E2)
        ax.plot([x0, x1], [y, y], color=u["color"], lw=1.5, zorder=4)
        ax.annotate("", xy=((x0 + x1) / 2., y - 0.6), xytext=((x0 + x1) / 2., y), zorder=4,
                    arrowprops=dict(arrowstyle="-|>", color=u["color"], lw=1.5))
        if u.get("label"):
            handles.append(Line2D([], [], color=u["color"], lw=1.5, label=u["label"]))
        ul_top.append(y)

    #-----------------------------------------------------------
    # Axes: 8 decades below the maximum

    xlim = (np.floor(np.log10(E[0])), np.ceil(np.log10(E[-1])))
    top = np.ceil(max([np.nanmax(log10(L))] + ul_top) + 0.5)
    ax.set_xlim(*xlim)
    ax.set_ylim(top - 8., top)
    int_ticks(ax, "x", *xlim, max_major=8)
    int_ticks(ax, "y", top - 8., top, max_major=5)
    style.finish_axes(ax)
    ax.set_xlabel(r"$\log_{10}(\epsilon$ [eV]$)$")
    ax.set_ylabel(r"$\log_{10}(\epsilon L_\epsilon\,[{\rm erg\,s^{-1}}])$")

    # add_artist clips to the axes: the legend inside goes first, the one above last
    if len(shocks) > 1:
        styles = [Line2D([], [], color="grey", ls=LS[s], lw=1.5, label=s) for s in shocks]
        ax.add_artist(ax.legend(handles=styles, loc="best", handlelength=1.5,
                                title=run.label, alignment="left"))
    else:
        style.system_label(ax, run.label)
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.), borderaxespad=0.2,
              ncol=len(handles) if len(handles) <= 6 else (len(handles) + 1) // 2, columnspacing=1.)
    run.save(fig, "SED")


if __name__ == "__main__":
    main()
