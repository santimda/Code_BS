# Code_BS

Code to compute the non-thermal, free-free and H-alpha emission of the bow shock (BS) produced by the wind of a runaway massive star. It has two versions: a one-zone model (`One-Zone/`), approximate and fast, and an extended (multi-zone) model (`Extended/`), more detailed and slower. The model is described in [del Palacio et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..13D/abstract).


## Quick start

Requirements: `gfortran` (with OpenMP), `make`, and Python 3 with `numpy` and `matplotlib`.

```
cd Extended
make test                    # ~4 min: checks that the code works and returns the expected outputs for BD+43 3654
./make_plots_runaway.sh     # run and plot a runaway star
./make_plots.sh             # run and plot a (hypothetical) hypervelocity star
```

The two `make_plots` scripts do the whole cycle: they ask for a system, copy its `system_parameters.f90` (from `Runaway-Stars/` or `Hypervelocity_stars/`), compile, run, and make the plots in `figures/<system>/`. Extra arguments are passed to the plotting script, e.g. `sh make_plots_runaway.sh --no-usetex`.


## Model summary

The stellar wind and the interstellar medium (ISM) form a BS with a reverse shock (RS) in the wind and a forward shock (FS) in the ISM, separated by a contact discontinuity (CD). In the extended model both shocks are thin surfaces on the CD, whose shape is that of [Christie et al. (2016)](https://scixplorer.org/abs/2016MNRAS.459.2420C/abstract) for a warm ISM ([Wilkin 1996](https://ui.adsabs.harvard.edu/abs/1996ApJ...459L..31W/abstract) for a cold one). The emitter is divided into cells along the CD, with analytical thermodynamics in each cell.

Electrons and protons are accelerated at the adiabatic shocks (broken power law with indices `alpha_1`, `alpha_2` and an exponential cutoff) and are convected along linear emitters, rotated around the symmetry axis to build the 3-D emitter. They cool and radiate by synchrotron, inverse Compton (IC) on the stellar and dust photons ([Khangulyan et al. 2014](https://ui.adsabs.harvard.edu/abs/2014ApJ...783..100K/abstract), anisotropic), relativistic Bremsstrahlung and proton-proton collisions. The synchrotron emission is free-free absorbed by the wind and the ionized shocked ISM.

The FS (`radiative_FS` in `system_parameters.f90`) can be:

- **adiabatic**, like the RS, with particle acceleration ([Martinez et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.102M/abstract)); relevant for hypervelocity stars;
- **radiative**, a hot layer followed by an isothermal layer (IL) ionized by the star, with only thermal free-free emission ([Martinez et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...680A..99M/abstract)).

The code also estimates the emission of the particles that escape into the HII region, and makes synthetic maps (synchrotron, free-free, H-alpha) convolved with the instrumental beam.

The main free parameters are the ratio of magnetic to thermal pressure, `eta_B` (`frac_sh = sqrt(eta_B)`), the fraction of the wind kinetic power that goes to relativistic particles, `f_NT` (`frac_NT`, of which `inj_e` goes to electrons), and, for a radiative FS, the reduction of the IL width, `eta_H` ≤ 1 (`frac_H`).

![Sketch (not to scale) of the BS of a runaway star. The FS and the RS are separated by the CD (black line). The post-shock ISM, the shocked medium ionised by the star and the shocked wind are shown in red, orange and blue; the arrows are the streamlines of the shocked gas.](docs/Radiative_BS.jpg)


## Usage

**Running a system.** Use `make_plots_runaway.sh` or `make_plots.sh` (see Quick start). By hand: copy a `system_parameters.f90` from `Runaway-Stars/<system>/` or `Hypervelocity_stars/<scenario>/` to `Extended/`, then `make`, `./run.x` and `python3 plots/make_plots.py`. The outputs go to `output/`; `make clean` removes them and the compiled files; `make DEBUG=Y` adds runtime checks. The number of OpenMP threads is `OMP_NUM_THREADS`, or `nth` (`global.f90`) if it is not set.

**Parameters.**

| file | contents | changed |
|---|---|---|
| `system_parameters.f90` | the system: distance, inclination, star, wind, ISM, `L_IR`, `radiative_FS`, `frac_sh`, injection (`alpha_1`, `alpha_2`, `frac_NT`, `inj_e`), `frac_H` | for each system |
| `controls.f90` | modules to run (`yes`/`no`), map frequencies, `verbose` | for each run |
| `global.f90` | resolution (`m_theta`, `mE_e`, `mphi`, ...), threads, energy ranges, output directory | rarely |
| `common_parameters.f90` | parameters common to all the systems and derived quantities | rarely |
| `constants.f90` | physical constants | never |

A new parameter of `system_parameters.f90` must be added to all its copies, including `Tests/cases/`.

**Plots.** `python3 plots/make_plots.py` plots the modules that were run (according to `controls.f90`): `--dry-run` lists them, `--jobs sed maps` makes only some, `--no-usetex` avoids LaTeX. The observations shown for each system are set in `plots/plot_config.py`; the data are in `data/`.

**Tests.** `make test`: unit tests, regression of the BD+43 3654 values (`Tests/reference_terminal_output.json`) and physical checks (energy budgets, comparison with One-Zone, OpenMP). If the results change on purpose, update the reference with `python3 Tests/test_terminal_output.py --update-reference`. See [`Extended/Tests/README.md`](Extended/Tests/README.md).


## Modules (Extended)

| file | contents |
|---|---|
| `main.f90` | calls the modules |
| `initialize.f90` | allocations, frequency grid |
| `CD.f90`, `CD_rot.f90` | shape of the CD and 3-D emitter |
| `CD_thermo.f90` | thermodynamics and injected power along the CD |
| `shock_thermal.f90` | radiative FS: hot and isothermal layers, free-free emission |
| `CD_diste.f90`, `CD_distp.f90` | electron and proton distributions, transport along the linear emitters |
| `CD_radsyn.f90`, `CD_radic.f90`, `CD_radbr.f90`, `CD_radpp.f90` | synchrotron, IC, relativistic Bremsstrahlung, p-p SEDs |
| `CD_abssyn.f90`, `opac_ff.f90` | free-free absorption (wind, IL, HII region), projected emission at the map frequencies |
| `radstr.f90` | particles escaped into the HII region |
| `rad_H_alpha.f90`, `SW_conv.f90` | H-alpha map, convolution of the maps with the beam |
| `BB.f90`, `diffusion.f90`, `gaunt_factor.f90`, `numerical_utils.f90`, `map_freqs.f90` | blackbody, diffusion, Gaunt factor (van Hoof et al. 2014), numerical routines, map frequencies |

**One-Zone** has the same parameter files, `diste.f90`/`distp.f90` for the particle distributions and `rad*.f90` for the SEDs. `sh make_plots.sh` runs a scenario (folder with its `system_parameters.f90`, e.g. `B2-60-d/`) and plots it with gnuplot.


## References

If you use this code, please cite [del Palacio et al. (2018), A&A, 617, A13](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..13D/abstract) (model), [Martinez et al. (2022), A&A, 661, A102](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.102M/abstract) (adiabatic FS) and [Martinez et al. (2023), A&A, 680, A99](https://ui.adsabs.harvard.edu/abs/2023A%26A...680A..99M/abstract) (radiative FS).
