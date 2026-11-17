# Code_BS

Code to compute the non-thermal, free-free and H-alpha emission from the bow shock (BS) produced by the wind of a runaway massive star. Two versions of the model are available: i) a one-zone model (`One-Zone/`), which is more approximate and runs faster, and ii) an extended (multi-zone) model (`Extended/`), which is more detailed and time-consuming. The model is detailed in [del Palacio et al. (2018), A&A, 617, A13](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..13D/abstract).


## Quick start

Requirements: `gfortran` (with OpenMP), `make`, and Python 3 with `numpy` (and `matplotlib` for the plots).

```
cd Extended
make test                    # ~7 min: consistency checks and checks that the code returns the expected outputs for BD+43 3654
sh make_plots_runaway.sh     # choose a runaway star: copies its parameters, compiles, runs and plots
sh make_plots.sh             # choose a hypothetical hyper-velocity star: copies its parameters, compiles, runs and plots
```

`make test` compiles and runs copies of the code in temporary directories (the sources are not modified). It prints the main luminosities and radio fluxes of BD+43 3654 next to the reference values, and ends with `All the tests passed.` See [`Extended/Tests/README.md`](Extended/Tests/README.md).

## Model summary

The supersonic stellar wind interacts with the interstellar medium (ISM) and forms a BS with two shocks separated by a contact discontinuity (CD): a reverse shock (RS) in the stellar wind and a forward shock (FS) in the ISM. In the extended model, both shocks are treated as thin, axisymmetric surfaces co-spatial with the CD, whose shape follows the analytical solution of [Christie et al. 2016](https://scixplorer.org/abs/2016MNRAS.459.2420C/abstract) for a star moving through a warm ISM (it reduces to [Wilkin 1996](https://ui.adsabs.harvard.edu/abs/1996ApJ...459L..31W/abstract) for a cold ISM). The emitter is divided into cells along the CD, labelled by the polar angle θ measured from the apex. The shocked gas flows along the CD away from the apex, and the thermodynamic quantities in each cell are calculated using analytical prescriptions. The calculation extends up to the angle θ_max where the shock vanishes.

Relativistic electrons and protons are accelerated via diffusive shock acceleration at the adiabatic shocks. They are injected with a broken power law with an exponential cutoff (indices `alpha_1` and `alpha_2` below and above `E_break`), with the maximum energy set by the balance between acceleration and energy losses. The particles are convected with the shocked flow along a set of one-dimensional linear emitters, and the full 3-D emitter is obtained by rotating them around the symmetry axis. Along the way, the particles cool and emit radiation via synchrotron, inverse-Compton (IC) scattering of the stellar and dust infrared photon fields, relativistic Bremsstrahlung, and proton-proton (pp) collisions. The IC emission is computed with the formulae of [Khangulyan et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...783..100K/abstract), which account for the anisotropy of the stellar photon field.

The FS can be modeled in two ways (`radiative_FS` in `system_parameters.f90`):

- **Adiabatic FS**: Relevant in the context of hypervelocity-stars. like the RS, the FS is a strong adiabatic shock where particles are accelerated. The thermodynamic quantities of both shocks are derived from Bernoulli's equation, distinguishing between the subsonic and supersonic regions of the shocked flow.
- **Radiative FS**: the shocked ISM is described by a double-layer model, with a hot layer (HL) followed by an isothermal layer (IL) kept warm (T ~ 8000 K) and ionized by the star. In this case, only the thermal free-free emission of the FS is computed.

The synchrotron emission can be absorbed via free-free absorption by the unshocked stellar wind and by the IL. The code also estimates the emission of the relativistic particles that escape the BS into the HII region (Strömgren sphere) ionized by the star, produces synthetic emission maps (synchrotron, free-free, and H-alpha), and convolves them with the instrumental beam to compare with observations.

Besides system-specific parameters, the model has two free parameters that have a strong impact on the non-thermal output: the ratio between magnetic and thermal pressure, `eta_B` (`frac_sh = sqrt(eta_B)` in `system_parameters.f90`), and the fraction of the available wind kinetic power that is transferred to relativistic particles, `f_NT` (`frac_NT`, of which `inj_e` goes to electrons). In the radiative FS model, the parameter `eta_H` ≤ 1 (`frac_H`) reduces the width of the IL to account for processes not included in the model, such as thermal conduction and instabilities.

For additional details of the model we refer to [del Palacio et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..13D/abstract), [Martinez et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.102M/abstract), and [Martinez et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023A%26A...680A..99M/abstract).

![Sketch (not to scale) of a BS of a runaway star. The FS and the RS are separated by the CD, represented by a solid black line. The immediate post-shock external medium, the shocked medium ionised by the star, and the shocked stellar wind are shown in red, orange, and blue, respectively. The solid lines with arrows indicate the shocked gas streamlines, from the injection position up to the tail of the BS.](docs/Radiative_BS.jpg)


## Code usage and configuration

The modular code is written in Fortran 90, articulated by a *Main* program, compiled with *Make* (`make`), and run as `./run.x`. The output files are written to `output/` and can be removed, together with the compiled files, with `make clean`. `make DEBUG=Y` compiles with runtime checks (array bounds, floating-point exceptions).

The extended model is parallelized with OpenMP. The number of threads is set once at the beginning of the run: the value of the `OMP_NUM_THREADS` environment variable is used if it is defined, and `nth` (in `global.f90`) otherwise.

### Parameters

The parameters of the extended model are divided according to how often they change:

| file | contents | when to change it |
|---|---|---|
| `system_parameters.f90` | physical parameters of the system: distance and inclination, star and wind, ISM, dust luminosity `L_IR`, FS model (`radiative_FS`), magnetic field (`frac_sh`), injection (`alpha_1`, `alpha_2`, `frac_NT`, `inj_e`), IL width (`frac_H`) | for each system |
| `controls.f90` | which modules are run (`yes`/`no`), frequencies of the radio maps (`map_catalogue`), `verbose` | for each run |
| `global.f90` | numerical resolution (`m_theta`, `mE_e`, `mE_p`, `mphi`, ...), number of threads, energy ranges and `E_break`, output directory; also the shared variables and physical functions | rarely |
| `common_parameters.f90` | parameters common to all the systems (abundances, `R0`, derived wind and ISM quantities, ...) | rarely |
| `constants.f90` | physical and numerical constants | never |

A copy of `system_parameters.f90` for each system is kept in `Runaway-Stars/<system>/` and `Hypervelocity_stars/<scenario>/`. Before a run, copy the desired file to `Extended/` (`make_plots_runaway.sh` and `make_plots.sh` do it, for runaway and hypervelocity stars, respectively). A new parameter must be added to all the copies, including the frozen test cases in `Tests/cases/`.

### Plots

`python3 plots/make_plots.py` makes the plots of the modules that were run (according to `controls.f90`) in `figures/<system>/`; `--dry-run` lists them, `--jobs sed maps` makes only some of them, and `--no-usetex` avoids LaTeX. The observational data and the instrument sensitivities are in `data/`, and the observations shown for each system are set in `plots/plot_config.py`.

### Tests

`make test` runs `Tests/test_set.sh`: unit tests of the numerical routines and physical functions, a regression of the values printed in the terminal for BD+43 3654 (`Tests/reference_terminal_output.json`), and physical checks (energy budgets, comparison with the One-Zone model, OpenMP). If the results change on purpose, the reference is rewritten with `python3 Tests/test_terminal_output.py --update-reference`. Details in [`Extended/Tests/README.md`](Extended/Tests/README.md).


## Extended

**main.f90**: main program that calls the modules.

**constants.f90**, **system_parameters.f90**, **common_parameters.f90**, **global.f90**, **controls.f90**: parameters (see above). `global.f90` imports the first three and is used by all the other modules.

**initialize.f90**: allocates the CD vectors and writes the frequency grid.

**CD.f90**: calculates the position of the CD. The extended emitter is made up of emission cells at these same positions (for both shocks, FS and RS).

**CD_rot.f90**: rotates the CD to build the 3-D emitter.

**CD_thermo.f90**: calculates the thermodynamics and the injected power at each position along the CD.

**shock_thermal.f90**: calculates the structure of the radiative FS (hot and isothermal layers) and its thermal free-free emission.

**CD_diste.f90**, **CD_distp.f90**: calculate the non-thermal electron and proton populations at each position of the CD (for each shock), including the transport along the linear emitters.

**CD_rad(XX).f90**: calculates the SED of each process (syn = synchrotron, ic = inverse Compton, br = relativistic Bremsstrahlung, pp = proton-proton) from the particle populations of **CD_dist(X).f90**.

**CD_abssyn.f90**: applies free-free absorption by the stellar wind and the IL (only for a radiative FS) to the synchrotron emission, and projects the emission at the map frequencies in the plane of the sky.

**opac_ff.f90**: free-free absorption and emission of the unshocked stellar wind and of the HII region.

**BB.f90**: blackbody emission of the star.

**radstr.f90**: estimates the emission (relativistic Bremsstrahlung and pp) of the relativistic particles that escape from the BS into the HII region (Strömgren sphere) ionized by the star.

**rad_H_alpha.f90**: calculates the H-alpha emission map of the IL and convolves it with the instrumental beam.

**SW_conv.f90**: convolves the emission maps with the instrumental (synthesized) beam to compare with the observed maps.

**map_freqs.f90**, **map_freq_type.f90**: frequencies of the radio maps (set in `controls.f90`).

**diffusion.f90**: diffusion coefficients. **gaunt_factor.f90**: free-free Gaunt factor (table of van Hoof et al. 2014, `gauntff.dat`). **numerical_utils.f90**: logarithmic grids and integration.


## One-Zone

**system_parameters.f90**, **global.f90**, **constants.f90**, **controls.f90**: parameters, as in the extended model. A copy of `system_parameters.f90` for each scenario is kept in its own folder (e.g. `B2-60-d/`).

**diste.f90**, **distp.f90**: calculate the non-thermal electron and proton populations.

**rad(XX).f90**: calculate the SED of each process from the particle populations of **dist(X).f90**; **abssyn.f90** and **opac_ff.f90**: free-free absorption.


## References and citation

If you use this code in your research, please cite the relevant papers:

- **Model description:** [del Palacio et al. (2018), A&A, 617, A13](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..13D/abstract)
- **Model update: adiabatic FS model:** [Martinez et al. (2022), A&A, 661, A102](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.102M/abstract)
- **Model update: Radiative double-layer FS model:** [Martinez et al. (2023), A&A, 680, A99](https://ui.adsabs.harvard.edu/abs/2023A%26A...680A..99M/abstract)
