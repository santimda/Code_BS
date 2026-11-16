#!/usr/bin/env python3
"""
Thermodynamics along the shocks (CD_thermo):
  thermo_{RS,FS}: normalised quantities vs theta
  rad_adi: log10(t_cool/t_adv) vs theta
  thermo_vs_RH_{RS,FS}: ratio to the Rankine-Hugoniot values
"""
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

from bsplot import style
from bsplot.common import Run, base_parser
from bsplot.io import log10

# column, label, colour, line style
THERMO = [(1, r"$\rho$", "olive", "-"), (2, r"$P$", "royalblue", "-"),
         (4, r"$T$", "darkred", "-"), (5, r"$B$", "mediumorchid", "-"),
         (7, r"$c_s$", "grey", "--"), (3, r"$v_t$", "coral", "--"),
         (6, r"$v_\perp$", "springgreen", "--")]
RH = [(3, r"$\rho$", "olive", "-"), (2, r"$P$", "royalblue", "-"), (1, r"$v_t$", "coral", "-")]


def theta_axis(ax, theta):
    ax.set_xlim(0., 10. * np.ceil(theta.max() / 10.))
    ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.set_xlabel(r"$\theta\,[^\circ]$")


def ylim_without_apex(values, ref):
    # y range without l = 1, where t_adv -> infinity and v_t -> 0; includes ref
    v = np.concatenate(values)
    v = v[np.isfinite(v)]
    lo, hi = min(np.percentile(v, 0.5), ref), max(np.percentile(v, 99.5), ref)
    pad = 0.08 * (hi - lo) or 0.1 # all equal to ref (e.g. radiative FS vs RH)
    return lo - pad, hi + pad


def plot(run, name, theta, curves, ylabel, label, ylim=None, ref=None):
    # curves: (y, label, colour, line style); ref: dashed horizontal line
    fig, ax = style.new_figure()
    for y, lab, color, ls in curves:
        ax.plot(theta, y, color=color, ls=ls, lw=2, alpha=0.8, label=lab)
    if ref is not None: # as a curve, so that loc="best" avoids it
        ax.plot([0., 180.], [ref, ref], color="black", ls="--", lw=1, scalex=False)
        ylim = ylim_without_apex([c[0][1:] for c in curves], ref)
    theta_axis(ax, theta)
    ax.set_ylim(*ylim)
    if ref is None:
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    else:
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    style.finish_axes(ax)
    ax.set_ylabel(ylabel)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.), ncol=len(curves),
              columnspacing=0.8, handlelength=1.2)
    style.system_label(ax, label, box=True)
    run.save(fig, name)


def main(argv=None):
    run = Run(base_parser(__doc__), argv)
    adi = []
    for shock in run.shocks():
        label = f"{run.label}\n{shock}"
        d = run.out.read("thermonorm", shock=shock)
        if d is not None:
            plot(run, f"thermo_{shock}", np.degrees(d[:, 0]),
                 [(d[:, c], *rest) for c, *rest in THERMO], "normalized value", label, ylim=(0., 1.05))
        d = run.out.read("vs_RH", shock=shock)
        if d is not None:
            plot(run, f"thermo_vs_RH_{shock}", np.degrees(d[:, 0]),
                 [(d[:, c], *rest) for c, *rest in RH], r"$X / X_{\rm RH}$", label, ref=1.)
        d = run.out.read("rad_adi", shock=shock)
        if d is not None:
            adi.append((np.degrees(d[:, 0]), log10(d[:, 1]), shock))

    if adi:
        plot(run, "rad_adi", adi[0][0],
             [(np.interp(adi[0][0], th, y), s, {"RS": "royalblue", "FS": "tomato"}[s], "-")
              for th, y, s in adi],
             r"$\log_{10}(t_{\rm cool}/t_{\rm adv})$", run.label, ref=0.)


if __name__ == "__main__":
    main()
