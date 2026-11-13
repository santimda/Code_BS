#!/usr/bin/env python3
"""
Particle distributions, log10(E^2 N(E)) vs log10(E): one curve per position
of the shock (coloured by theta) and the total of the shock (black).
One figure per shock with non-thermal particles: dist{e,p}_{RS,FS}.
"""
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize

from bsplot import style
from bsplot.common import Run, base_parser, int_ticks
from bsplot.io import log10

eV = 1.602176e-12


def main(argv=None):
    p = base_parser(__doc__)
    p.add_argument("--species", choices=["e", "p"], default="e")
    run = Run(p, argv)
    sp = run.args.species
    n_emit = 4. * (run.parameter("mphi", 30.) - 1.) # linear emitters in a ring (global.f90)

    for shock in run.shocks_nt():
        d = run.out.read(f"dist{sp}", shock=shock)
        if d is None:
            continue
        theta, E = np.degrees(d[:, 0]), d[:, 1]
        E2N = (E * eV)**2 * d[:, 2] * n_emit # whole ring

        fig, ax = style.new_figure()
        norm = Normalize(0., theta.max())
        ys = []
        for th in np.unique(theta):
            i = np.where((theta == th) & (E2N > 0))[0][1:]
            if i.size:
                ax.plot(log10(E[i]), log10(E2N[i]), color=cm.plasma(norm(th)), lw=1.)
                ys.append(log10(E2N[i]))

        # Total of the shock (the file has one quadrant)
        t = run.out.read(f"dist{sp}_total", shock=shock)
        if t is not None:
            ys.append(log10(4. * (t[:, 0] * eV)**2 * t[:, 1]))
            ax.plot(log10(t[:, 0]), ys[-1], color="black", lw=3)

        # Axes: 11 decades below the maximum
        xlim = (np.floor(np.log10(E[E > 0].min())), np.ceil(np.log10(E.max())))
        top = np.ceil(np.nanmax(np.concatenate(ys)) + 0.5)
        ax.set_xlim(*xlim)
        ax.set_ylim(top - 11., top)
        int_ticks(ax, "x", *xlim, max_major=4)
        int_ticks(ax, "y", top - 11., top)
        style.finish_axes(ax)
        ax.set_xlabel(rf"$\log_{{10}}\,(E_{sp} \, [{{\rm eV}}])$")
        ax.set_ylabel(rf"$\log_{{10}}\,(E_{sp}^2 N_{sp}(E_{sp})\,[{{\rm erg}}])$")

        # Colour bar with the height of the axes
        fig.canvas.draw()
        pos = ax.get_position()
        cbar = fig.colorbar(cm.ScalarMappable(norm, cm.plasma),
                            cax=fig.add_axes([pos.x1 + 0.02, pos.y0, 0.03, pos.height]))
        cbar.set_ticks(np.arange(0., theta.max() + 1., 30.))
        cbar.set_label(r"$\theta\,[^\circ]$")

        style.system_label(ax, f"{run.label}\n{shock}", box=True)
        run.save(fig, f"dist{sp}_{shock}")


if __name__ == "__main__":
    main()
