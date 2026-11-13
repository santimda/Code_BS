#!/usr/bin/env python3
"""
Plots according to controles.f90.

  python3 plots/make_plots.py
  python3 plots/make_plots.py --dry-run # list the plots, without making them
  python3 plots/make_plots.py --jobs sed maps # only these plots
  python3 plots/make_plots.py --fontsize 10 --fmt pdf png

The remaining options (--fontsize, --font, --data, --figdir, --label, --system,
--fmt, ...) are passed to each plotting script.
"""
import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import plot_dist
import plot_maps
import plot_radio_sed
import plot_sed
import plot_termo
import plot_timescales
from bsplot.common import base_parser, default_figdir
from bsplot.fortran import read_controls, read_string


def any_of(f, *names):
    return any(f.get(n.lower(), False) for n in names)


def all_of(f, *names):
    return all(f.get(n.lower(), False) for n in names)


#===========================================================
# Plots: name, condition on the controls, script, extra arguments.
# A plot is made only if the module that writes its data was run (main.f90).

JOBS = [
    ("termo", lambda f: all_of(f, "termo_calc"),
     plot_termo, []),
    ("diste", lambda f: all_of(f, "electrons", "diste_calc"),
     plot_dist, ["--species", "e"]),
    ("distp", lambda f: all_of(f, "protons", "distp_calc"),
     plot_dist, ["--species", "p"]),
    ("tcoole", lambda f: all_of(f, "electrons", "diste_calc", "escribir_e"),
     plot_timescales, ["--species", "e"]),
    ("tcoolp", lambda f: all_of(f, "protons", "distp_calc", "escribir_p"),
     plot_timescales, ["--species", "p"]),
    ("sed", lambda f: (f.get("electrons") and any_of(f, "radsyn_sw_calc", "radbr_sw_calc",
                                                     "radic_st_calc", "radic_ir_calc"))
                      or all_of(f, "protons", "radpp_calc"),
     plot_sed, []),
    ("radio_sed", lambda f: all_of(f, "electrons", "radsyn_sw_calc") or all_of(f, "xrad_calc"),
     plot_radio_sed, []),
    ("maps", lambda f: all_of(f, "maps_calc", "convolution", "rotate")
                       or all_of(f, "rad_alpha", "convolution_alpha"),
     plot_maps, []),
]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    args = base_parser(__doc__).parse_args(argv)

    #-----------------------------------------------------------
    # Run information

    if not os.path.isfile(args.controls):
        raise SystemExit(f"cannot find {args.controls}")
    code_dir = os.path.dirname(os.path.abspath(args.controls))
    flags = read_controls(args.controls)
    main_f90 = os.path.join(code_dir, "main.f90")
    rad_fs = read_controls(main_f90).get("radiative_fs") if os.path.isfile(main_f90) else None
    system = read_string([os.path.join(code_dir, "system_parameters.f90")], "system", "")

    print(f"Controls: {args.controls}")
    if rad_fs is not None:
        print(f"FS: {'radiative' if rad_fs else 'adiabatic'} (main.f90)")
    print(f"System: {system or '(not found in system_parameters.f90)'}")
    print(f"Data: {args.data}")
    print(f"Figures: {args.figdir or default_figdir(system)}\n")

    #-----------------------------------------------------------
    # Plots to make

    to_run = [j for j in JOBS if j[1](flags) and (args.jobs is None or j[0] in args.jobs)]
    skipped = [j[0] for j in JOBS if not j[1](flags)]
    print("Plot list: " + (", ".join(j[0] for j in to_run) or "none"))
    if skipped:
        print("Not plotted: " + ", ".join(skipped))
    if args.dry_run:
        return 0

    #-----------------------------------------------------------
    # Make them (an error in one plot does not stop the others)

    failed = []
    for name, _, script, extra in to_run:
        print(f"\n[{name}]")
        try:
            script.main(extra + argv)
        except SystemExit as e:
            if e.code not in (0, None):
                failed.append(name)
        except Exception:
            traceback.print_exc()
            failed.append(name)
    print("\nDone." if not failed else f"\nFailed: {', '.join(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
