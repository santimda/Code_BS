"""
Plot style: A&A fonts and sizes, same axes box in every figure.
"""
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

COLUMN_WIDTH = 3.54399 # A&A column width [in]
FONTSIZE = 9. # A&A body text [pt]
BOX_ASPECT = 2. / 3. # height/width of the axes box, the same in all figures

# Colour-blind safe palette (Okabe & Ito)
cb_colors = {"black": "#000000", "orange": "#E69F00", "green": "#009E73",
             "blue": "#0072B2", "red": "#D55E00", "purple": "#CC79A7"}


def set_style(fontsize=FONTSIZE, font="times", usetex=True):
    # font: 'times' (txfonts, as A&A) or 'cm' (Computer Modern)
    if usetex and shutil.which("latex") is None:
        print("LaTeX not found: using mathtext")
        usetex = False
    mpl.rcParams.update({
        "text.usetex": usetex,
        "font.family": "serif",
        "font.size": fontsize, "axes.labelsize": fontsize, "legend.fontsize": fontsize,
        "xtick.labelsize": fontsize, "ytick.labelsize": fontsize,
        "axes.linewidth": 1.5,
        "xtick.major.size": 4, "xtick.minor.size": 2, "ytick.major.size": 4, "ytick.minor.size": 2,
        "xtick.major.width": 1.2, "xtick.minor.width": 1., "ytick.major.width": 1.2, "ytick.minor.width": 1.,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "legend.frameon": False, "legend.handlelength": 1.,
        "savefig.dpi": 500, "savefig.bbox": "tight", "savefig.pad_inches": 0.01,
    })
    if usetex:
        mpl.rcParams["text.latex.preamble"] = r"\usepackage{txfonts}" if font == "times" else ""
    elif font == "times":
        mpl.rcParams.update({"font.serif": ["STIXGeneral", "Times New Roman"], "mathtext.fontset": "stix"})
    else:
        mpl.rcParams["mathtext.fontset"] = "cm"


def new_figure():
    fig, ax = plt.subplots(figsize=(COLUMN_WIDTH, BOX_ASPECT * COLUMN_WIDTH))
    ax.set_box_aspect(BOX_ASPECT)
    return fig, ax


def finish_axes(ax):
    # ticks inside on the four sides, no labels on the minor ones
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())


def system_label(ax, text, loc="best", color="black", box=False):
    # An empty legend with a title, so that loc="best" places it where it covers
    # the fewest curves. It is kept as a separate artist, so other legends of
    # the axes neither replace it nor are replaced by it.
    old = ax.get_legend()
    leg = ax.legend(handles=[], title=text, loc=loc, frameon=box, framealpha=0.75,
                    edgecolor="none", borderpad=0.3, alignment="left",
                    title_fontsize=mpl.rcParams["font.size"])
    leg.get_title().set_color(color)
    leg.set_zorder(10)
    ax.add_artist(leg)
    ax.legend_ = old
