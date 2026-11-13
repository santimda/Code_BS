#!/usr/bin/env python3
"""
Emission maps: one per frequency of map_catalogue with 'yes', plus map_avg,
from the convolved maps (convolution = yes), and the H-alpha map (rad_alpha
and convolution_alpha). Same scale in x and y; the beam is drawn in the
upper right corner. Geometry and contours: plot_config.py.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator

from bsplot import style
from bsplot.common import Run, base_parser
from bsplot.fortran import read_map_catalogue
from bsplot.io import load_gridded

KB, MP, AU, PC, MSUN_YR = 1.380649e-16, 1.6726219e-24, 1.495978707e13, 3.085677581e18, 1.989e33 / 3.15576e7


def rotate(x, y, ang, sign=1):
    # rotation of the gnuplot scripts, from arcsec to arcmin (sign = 1 also flips y)
    return ((x * np.cos(ang) + sign * y * np.sin(ang)) / 60.,
            (x * np.sin(ang) - sign * y * np.cos(ang)) / 60.)


def R0_proj(p, distance, incli):
    # analytical stagnation point distance, projected [arcmin]
    alpha = KB * p["T_ISM"] / (p["v_st"]**2 * MP * p["mu_ISM"])
    R0 = np.sqrt(p["v_w"] * p["Mdot"] * MSUN_YR / (4. * np.pi * (1. + alpha) * MP * p["n_ISM"] * p["mu_ISM"] * p["v_st"]**2))
    return (R0 / AU) / (distance / PC) / 60. * np.cos(incli)


def draw(run, path, name, text, cbar_label, mc, fc, beam=None):
    # mc: map settings of the system; fc: settings of this frequency
    X, Y, Z = load_gridded(path)
    rot = np.deg2rad(mc["rot_deg"])
    incli = np.pi / 2. - np.deg2rad(mc["inclination_deg"])
    X, Y = rotate(X, Y, rot)

    fig, ax = style.new_figure()
    cmap = plt.get_cmap(fc.get("cmap", "viridis"))
    mesh = ax.pcolormesh(X, Y, Z, cmap=cmap, shading="auto", rasterized=True)
    ax.set_facecolor(cmap(0.)) # outside the grid of the code: zero flux
    levels = fc.get("levels", list(np.nanmax(Z) * np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
    if levels:
        ax.contour(X, Y, Z, levels=sorted(levels), colors="black", linewidths=1.)
    if fc.get("levels_red"):
        ax.contour(X, Y, Z, levels=sorted(fc["levels_red"]), colors="red", linewidths=1.)

    #-----------------------------------------------------------
    # Star, stagnation point of the code (cyan) and analytical one (magenta)

    if run.out.exists("R0_proj"):
        x, y = rotate(*run.out.read("R0_proj")[0], rot, sign=-1)
        x, y = x * np.cos(incli), y * np.cos(incli)
        ax.scatter(x, y, s=200, marker="+", c="cyan", linewidths=1.2, zorder=5)
        ax.annotate("", xy=(x, y), xytext=(0., 0.), zorder=6,
                    arrowprops=dict(arrowstyle="-|>", color="#6495ED", lw=1.8))
    if mc.get("R0"):
        r0, ang = R0_proj(mc["R0"], run.conf["distance"], incli), np.deg2rad(mc["arrow_deg"])
        ax.annotate("", xy=(-r0 * np.sin(ang), r0 * np.cos(ang)), xytext=(0., 0.), zorder=6,
                    arrowprops=dict(arrowstyle="-|>", color="magenta", lw=1.8))
    ax.scatter([0.], [0.], s=260, marker="*", c="purple", edgecolors="none", zorder=6)

    #-----------------------------------------------------------
    # Limits: those of plot_config, widened (centred) to fill the axes box
    # with the same scale in x and y

    (x0, x1), (y0, y1) = fc.get("xlim", mc["xlim"]), fc.get("ylim", mc["ylim"])
    if (y1 - y0) < style.BOX_ASPECT * (x1 - x0):
        pad = (style.BOX_ASPECT * (x1 - x0) - (y1 - y0)) / 2.
        y0, y1 = y0 - pad, y1 + pad
    else:
        pad = ((y1 - y0) / style.BOX_ASPECT - (x1 - x0)) / 2.
        x0, x1 = x0 - pad, x1 + pad
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect(1.)

    if beam is not None: # sigma (x, y) of map_catalogue ["]
        ax.add_patch(Ellipse((x0 + 0.92 * (x1 - x0), y0 + 0.90 * (y1 - y0)),
                             width=beam[1] / 60. / mc["cosdec"], height=beam[0] / 60.,
                             angle=mc["rot_deg"], fc="grey", ec="black", lw=0.8, zorder=7))

    # ticks every 5' (every 1' in small maps)
    for axis, size in ((ax.xaxis, x1 - x0), (ax.yaxis, y1 - y0)):
        step = 5. if size > 8. else 1.
        axis.set_major_locator(MultipleLocator(step))
        axis.set_minor_locator(MultipleLocator(step / 5. if step == 5. else step / 4.))
    style.finish_axes(ax)
    ax.set_xlabel(r"$\Delta$RA [arcmin]")
    ax.set_ylabel(r"$\Delta$DEC [arcmin]")

    # Colour bar with the height of the map
    fig.canvas.draw()
    pos = ax.get_position()
    fig.colorbar(mesh, cax=fig.add_axes([pos.x1 + 0.03, pos.y0, 0.03, pos.height])).set_label(cbar_label)

    style.system_label(ax, f"{run.label}\n{text}", loc="upper left",
                       color=fc.get("label_color", "lavender"))
    run.save(fig, name)


def main(argv=None):
    run = Run(base_parser(__doc__), argv)
    mc = run.conf["maps"]
    if mc is None:
        print(f"  no map settings for {run.system} in plot_config.py")
        return

    if run.flag("maps_calc") and run.flag("convolution"):
        catalogue, avg = read_map_catalogue(run.args.controls)
        for f in [f for f in catalogue + [avg] if f is not None and f.on]:
            fc = mc["per_freq"].get(f.label, {})
            if run.out.exists("map", label=f.label):
                draw(run, run.out.path("map", label=f.label), f"map_{f.label}",
                     fc.get("text", f.label.replace("MHz", " MHz").replace("GHz", " GHz")),
                     r"Flux [mJy beam$^{-1}$]", mc, fc, beam=f.sigma)
            else:
                print(f"  missing {run.out.path('map', label=f.label)}")

    if run.flag("rad_alpha") and run.flag("convolution_alpha") and run.out.exists("halpha"):
        fc = {"cmap": "Reds", "label_color": "black", "levels": [], **mc.get("halpha", {})}
        draw(run, run.out.path("halpha"), "map_Halpha", r"H$\alpha$",
             "Surface brightness [R]", mc, fc)


if __name__ == "__main__":
    main()
