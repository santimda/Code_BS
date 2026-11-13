#!/usr/bin/env python3
"""
Radio SED, S_nu [mJy] vs nu [GHz]: synchrotron (absorbed if abssyn_calc, the
unabsorbed one dotted), isothermal layer (IL) and total, with the observations
of plot_config.py. The hot layer of the FS only enters the total.
"""
import os

import numpy as np
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

from bsplot import style
from bsplot.common import Run, base_parser


def observations(run, o):
    # o["cols"] = (nu, dnu, S, dS); dS can be a relative error
    path = os.path.join(run.obs_dir, o["file"])
    if not os.path.isfile(path):
        print(f"  missing {path}")
        return None
    d = np.atleast_2d(np.loadtxt(path))
    i_nu, i_dnu, i_S, i_dS = o["cols"]
    dnu = d[:, i_dnu] / (2. if o.get("dnu_half") else 1.)
    S = d[:, i_S] - (d[:, o["S_minus"]] if "S_minus" in o else 0.)
    dS = i_dS * S if isinstance(i_dS, float) else d[:, i_dS]
    return d[:, i_nu], dnu, o.get("scale", 1.) * S, o.get("scale", 1.) * dS


def main(argv=None):
    run = Run(base_parser(__doc__), argv)
    conf = run.conf.get("radio_sed", {})
    fig, ax = style.new_figure()
    ax.set_xscale("log")
    ax.set_yscale("log")
    parts = [] # components of the total: (nu [GHz], S [mJy])

    def read(key, shock, col=1):
        d = run.out.read(key, shock=shock) if run.out.exists(key, shock=shock) else None
        return None if d is None else (d[:, 0] / 1.e9, d[:, col])

    #-----------------------------------------------------------
    # Synchrotron: absorbed (abs_maps), or unabsorbed (CD_radsyn)

    if run.flag("electrons") and run.flag("radsyn_sw_calc"):
        tags = ["BOTH"] if run.flag("both") else run.shocks_nt()
        absorbed = [read("S_syn_abs", t, 2) for t in tags] if run.flag("abssyn_calc") else []
        absorbed = [c for c in absorbed if c is not None]
        if absorbed:
            unabs = [read("S_syn_abs", t, 1) for t in tags if run.out.exists("S_syn_abs", shock=t)]
        else:
            unabs = [c for c in (read("S_syn", s) for s in run.shocks_nt()) if c is not None]
        if unabs:
            nu = unabs[0][0]
            ax.plot(nu, sum(np.interp(nu, *c) for c in unabs), color="royalblue",
                    lw=1 if absorbed else 2, ls=":" if absorbed else "-", alpha=0.8,
                    label=None if absorbed else "synchrotron")
        if absorbed:
            nu = absorbed[0][0]
            S = sum(np.interp(nu, *c) for c in absorbed)
            ax.plot(nu, S, color="royalblue", lw=2, alpha=0.8, label="synchrotron")
            parts.append((nu, S))
        elif unabs:
            parts.append((nu, sum(np.interp(nu, *c) for c in unabs)))

    #-----------------------------------------------------------
    # Free-free of the FS: isothermal layer (plotted) and hot layer

    if run.flag("xrad_calc") and "FS" in run.shocks():
        il, hl = read("S_ionized", "FS"), read("S_ff", "FS")
        if il is not None:
            ax.plot(*il, color="orange", lw=2, alpha=0.8, label="IL")
        parts += [c for c in (il, hl) if c is not None]

    if not parts:
        print("  nothing to plot")
        return

    nu = np.geomspace(min(c[0].min() for c in parts), max(c[0].max() for c in parts), 400)
    S = sum(np.interp(nu, *c, left=0., right=0.) for c in parts)
    ax.plot(nu, S, color="darkred", lw=3, alpha=0.8, label="total")
    model = ax.get_legend_handles_labels()

    #-----------------------------------------------------------
    # Observations

    obs = []
    for o in run.conf["radio_obs"]:
        if "vline" in o: # flux range at one frequency (e.g. WENSS)
            nu0, (s0, s1) = o["vline"]
            obs += ax.plot([nu0, nu0], [s0, s1], color=o["color"], lw=2, alpha=0.8,
                           zorder=2, label=o["label"])
            continue
        d = observations(run, o)
        if d is not None:
            obs.append(ax.errorbar(d[0], d[2], xerr=d[1], yerr=d[3], fmt="none",
                                   ecolor=o["color"], alpha=o.get("alpha", 0.8),
                                   lolims=o.get("lolims", False), zorder=2, label=o["label"]))

    #-----------------------------------------------------------
    # Axes

    ax.set_xlim(*conf.get("xlim", (nu[0], nu[-1])))
    ax.set_ylim(*conf.get("ylim", (S.max() / 30., 2. * S.max())))
    for axis, subs in ((ax.xaxis, (1.,)), (ax.yaxis, (1., 2.5, 5.))):
        axis.set_major_locator(LogLocator(subs=subs))
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_locator(LogLocator(subs=np.arange(1, 10)))
        axis.set_minor_formatter(NullFormatter())
    if "yticks" in conf:
        ax.set_yticks(conf["yticks"])
    style.finish_axes(ax)
    ax.set_xlabel(r"$\nu$ [GHz]")
    ax.set_ylabel(r"$S_\nu$ [mJy]")

    # add_artist clips to the axes: the legend inside goes first, the one above last
    if obs:
        ax.add_artist(ax.legend(handles=obs, loc="lower left"))
    ax.legend(*model, loc="lower center", bbox_to_anchor=(0.5, 1.), ncol=len(model[0]),
              columnspacing=1.2)
    style.system_label(ax, run.label)
    run.save(fig, "SED_radio")


if __name__ == "__main__":
    main()
