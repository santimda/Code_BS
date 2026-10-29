import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.ticker import LogLocator, LogFormatter
from matplotlib.ticker import LogLocator, NullFormatter

# Fonts
text_width = 441.01773/72.27
mpl.rcParams['figure.dpi'] = 500

# Esthetics for plots
nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": 'serif',
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 11,
        "font.size": 11,
        "axes.linewidth":1.5,
        "axes.titlesize":11,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "ytick.major.size": 4,     # major tick size in points
        "ytick.minor.size": 2,       # minor tick size in points
        "xtick.major.size": 4,     # major tick size in points
        "xtick.minor.size": 2,       # minor tick size in points
        "xtick.major.width": 1.2,     # major tick width in points
        "xtick.minor.width": 1.,     # minor tick width in points
}
mpl.rcParams.update(nice_fonts)

####################################################################
####################################################################

####################################################################
####################################################################

def load_gridded_data(path):
    data = np.loadtxt(path)
    x, y, flux = data[:, 0], data[:, 1], data[:, 2]
    
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    nx, ny = len(x_unique), len(y_unique)
    
    if nx * ny != len(flux):
        raise ValueError(
            f"Data doesn't match grid format ({nx}x{ny} "
        )
        
    # índice de cada punto dentro de la grilla
    xi = np.searchsorted(x_unique, x)
    yi = np.searchsorted(y_unique, y)
    
    Z = np.full((ny, nx), np.nan)
    Z[yi, xi] = flux
    
    X, Y = np.meshgrid(x_unique, y_unique)
    return X/60., Y/60., Z
 
####################################################################
####################################################################

def read_file(path):
    """
    Reads a file and creates a np array with its first two columns
    """
    datos = []

    with open(path, 'r') as f:
        for linea in f:
            linea = linea.strip()
            if not linea:
                continue

            valores = linea.split()
            if len(valores) >= 2:
                datos.append([float(valores[0]), float(valores[1])])

    return np.array(datos, dtype=np.float64)

####################################################################
####################################################################

def safe_log10(x, replace_inf=True):
    """Safe log10 that replaces -inf with 0 for zeros"""
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.log10(x)
        if replace_inf:
            result = np.where(np.isinf(result), 0, result)
            result = np.where(np.isnan(result), 0, result)
    return result

####################################################################
####################################################################
def rotate(x, y, ang, sign=+1):
    '''
    '''
    a = sign * ang + np.pi
    c, s = np.cos(a), np.sin(a)
    xr = x * c - y * s
    yr = x * s + y * c
    return xr, yr

####################################################################
####################################################################

def f_lineal(ax, x, ycurve, x0, dx, slope,
             dy=1.0, dtext=20.0,
             label='',
             **plot_kwargs):

    # Punto de anclaje sobre la curva
    y0 = np.interp(x0, x, ycurve)

    # Recta guía
    xx = np.linspace(x0, x0 + dx, 1000)
    yy = y0 + dy + slope*(xx - x0)

    ax.plot(xx, yy, '--', color='black', **plot_kwargs)

    # Punto medio de la recta
    i = len(xx)//2
    xmid = xx[i]
    ymid = yy[i]

    # Ángulo real en coordenadas de pantalla
    p0 = ax.transData.transform((xx[0], yy[0]))
    p1 = ax.transData.transform((xx[-1], yy[-1]))

    dxs = p1[0] - p0[0]
    dys = p1[1] - p0[1]

    angle = np.degrees(np.arctan2(dys, dxs))

    # Vector normal unitario (en pantalla)
    nx = -dys
    ny = dxs
    norm = np.hypot(nx, ny)

    nx /= norm
    ny /= norm

    # Punto medio en pantalla
    pmid = ax.transData.transform((xmid, ymid))

    # Posición del texto (desplazada dtext píxeles)
    ptext = (
        pmid[0] + dtext*nx/2.,
        pmid[1] + dtext*ny/2.
    )

    # Volver a coordenadas de datos
    xtext, ytext = ax.transData.inverted().transform(ptext)

    # Dibujar texto
    ax.text(
        xtext,
        ytext,
        label,
        rotation=angle,
        rotation_mode='anchor',
        ha='center',
        va='center'
    )

    return xx, yy


def f_vertical(ax, x,
           label='',
           y=40,
           dx=0.1,
           rotation=90,
           ymin=None,
           ymax=None,
           **line_kwargs):
	"""
	Dibuja una recta vertical y una etiqueta asociada.

	Parameters
	----------
	ax : matplotlib axis
	x : float
	    Posición de la recta.
	label : str
	    Texto de la etiqueta.
	y : float
	    Posición vertical de la etiqueta.
	dx : float
	    Desplazamiento horizontal respecto de la recta.
	rotation : float
	    Rotación del texto (por defecto 90°).
	line_kwargs :
	    Argumentos adicionales para ax.axvline().
	"""

	defaults = {
	    "color": "black",
	    "linestyle": ":"
	}

	defaults.update(line_kwargs)

	ax.plot([x, x], [ymin, ymax], **defaults)

	ax.text(
	    x + dx,
	    y,
	    label,
	    rotation=rotation,
	    va='bottom',
	    ha='left'
	)