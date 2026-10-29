# Code_BS

Software to compute the non-thermal and thermal (free-free) emission from the bow shock (BS) produced by the wind of a runaway massive star. Two versions of the model are available: i) a one-zone model, which is more approximate and runs faster, and ii) an extended (multi-zone) model, which is more detailed and time-consuming. The model is detailed in [del Palacio et al. (2018), A&A, 617, A13](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..13D/abstract).


## Model summary

The supersonic stellar wind interacts with the interstellar medium (ISM) and forms a BS with two shocks separated by a contact discontinuity (CD): a reverse shock (RS) in the stellar wind and a forward shock (FS) in the ISM. In the extended model, both shocks are treated as thin, axisymmetric surfaces co-spatial with the CD, whose shape follows [Wilkin (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...459L..31W/abstract). The emitter is divided into cells along the CD, labelled by the polar angle θ measured from the apex. The shocked gas flows along the CD away from the apex, and the thermodynamic quantities in each cell are calculated using analytical prescriptions. The calculation extends up to the angle θ_max where the shock vanishes.

Relativistic electrons and protons are accelerated via diffusive shock acceleration at the adiabatic shocks. They are injected with a power-law spectrum with an exponential cutoff (or a broken power law), with the maximum energy set by the balance between acceleration and energy losses. The particles are convected with the shocked flow along a set of one-dimensional linear emitters, and the full 3-D emitter is obtained by rotating them around the symmetry axis. Along the way, the particles cool and emit radiation via synchrotron, inverse-Compton (IC) scattering of the stellar and dust infrared photon fields, relativistic Bremsstrahlung, and proton-proton (pp) collisions. The IC emission is computed with the formulae of [Khangulyan et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...783..100K/abstract), which account for the anisotropy of the stellar photon field.

The FS can be modeled in two ways:

- **Adiabatic FS**: like the RS, the FS is a strong adiabatic shock where particles are accelerated. The thermodynamic quantities of both shocks are derived from Bernoulli's equation, distinguishing between the subsonic and supersonic regions of the shocked flow.
- **Radiative FS**: the shocked ISM is described by a double-layer model, with a hot layer (HL) followed by an isothermal layer (IL) kept warm (T ~ 8000 K) and ionized by the star. In this case, only the thermal free-free emission of the FS is computed.

The synchrotron emission can be absorbed via free-free absorption by the unshocked stellar wind and by the IL. The code also estimates the emission of the relativistic particles that escape the BS into the HII region (Strömgren sphere) ionized by the star, produces synthetic emission maps (synchrotron, free-free, and H-alpha), and convolves them with the instrumental beam to compare with observations.

Besides system-specific parameters, the model has two free parameters that have a strong impact on the non-thermal output: the ratio between magnetic and thermal pressure (`eta_B`) and the fraction of the available wind kinetic power that is transferred to relativistic particles (`f_NT`, split between electrons and protons). In the radiative FS model, the parameter `eta_H` (≤ 1) reduces the width of the IL to account for processes not included in the model, such as thermal conduction and instabilities.

For additional details of the model we refer to [del Palacio et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..13D/abstract), [Martinez et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.102M/abstract), and [Martinez et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023A%26A...680A..99M/abstract).


## Code usage

The modular code is written in Fortran 90, articulated by a *Main* program, compiled with *Make* (simply run `make` in the terminal), and run as `./run.x`. Runtime-generated data files are written to `output/`. Generated files can be safely removed using `make clean`.

The extended model is parallelized with OpenMP. The number of threads is set once at the beginning of the run: the value of the `OMP_NUM_THREADS` environment variable is used if it is defined, and `nth` (in `global.f90`) otherwise.

Some further details of each code are given below.


## One-Zone

**Makefile**: instructions for compiling using *make*.

**Makefile.depend**: details which modules depend on which other modules (needed for compilation).

**constants.f90**: useful physical and numerical constants. It is imported by **global.f90** to be used by all other modules as well.

**system_parameters.f90**: parameters of the specific system being studied. It is imported by **global.f90** to be used by all other modules as well.

**global.f90**: constants, variables and number of points in the simulation to be used by all other modules.

**controles.f90**: this sets which modules are going to run in the simulation by setting variables to *.True.* or *.False.*.

**diste.f90**: calculates the non-thermal electron population and saves it in a plain text file.

**distp.f90**: calculates the non-thermal proton population and saves it in a plain text file.

**rad(XX).f90**: calculates the SED by the respective process reading the corresponding non-thermal particle population from the file generated by **dist(X).f90**.


## Extended

**Makefile**: instructions for compiling using *make*.

**Makefile.depend**: details which modules depend on which other modules (needed for compilation).

**constants.f90**: useful physical and numerical constants. It is imported by **global.f90** to be used by all other modules as well.

**system_parameters.f90**: parameters of the specific system being studied. It is imported by **global.f90** to be used by all other modules as well.

**global.f90**: constants, variables, model parameters and the number of points in the simulation to be used by all other modules.

**controles.f90**: this sets which modules are going to run in the simulation by setting variables to *.True.* or *.False.*.

**main.f90**: main program that calls the modules. The FS model (radiative or adiabatic) is selected with the `radiative_FS` parameter at the beginning of this file.

**CD.f90**: calculates the position of the CD. The extended emitter is assumed to be made up of emission cells at these same positions (for both shocks, FS and RS).

**CD_rot.f90**: rotates the CD to simulate a 3-D emitter.

**CD_termo.f90**: calculates the thermodynamics at each position along the CD.

**shock_thermal.f90**: calculates the structure of the radiative FS (hot and isothermal layers) and its thermal free-free emission.

**CD_diste.f90**: calculates the non-thermal electron population at each position of the CD (for each shock) and saves it in a plain text file.

**CD_distp.f90**: calculates the non-thermal proton population at each position of the CD (for each shock) and saves it in a plain text file.

**CD_rad(XX).f90**: calculates the SED by the respective process (syn = synchrotron, ic = inverse Compton, br = relativistic Bremsstrahlung, pp = proton-proton) reading the corresponding non-thermal particle population from the file generated by **CD_dist(X).f90**.

**CD_abssyn.f90**: applies free-free absorption by the stellar wind and the IL to the synchrotron emission, and projects the emission at some particular frequencies in the plane of the sky, generating the input for making synthetic maps. Takes as an input the intrinsic synchrotron SED obtained with **CD_radsyn.f90**. The absorption by the IL is only applied if the radiative FS was computed.

**radstr.f90**: estimates the emission (relativistic Bremsstrahlung and pp) of the relativistic particles that escape from the BS into the HII region (Strömgren sphere) ionized by the star.

**rad_H_alpha.f90**: calculates the H-alpha emission map of the IL and convolves it with the instrumental beam.

**SW_conv.f90**: convolves the emission maps with the instrumental (synthesized) beam to compare with the observed maps. Takes as an input the emission maps obtained with **CD_abssyn.f90**.

**numerical_utils.f90**: reusable numerical routines (logarithmic grids and integration).


## Configuration workflow

The run configuration is divided into three files:

- **controles.f90**: selects which calculation modules are executed.
- **global.f90**: defines global run variables, model parameters, and numerical resolutions, such as the number of spatial and energy points.
- **system_parameters.f90**: defines the physical parameters of the system being modeled.

When modeling a different system, normally only `system_parameters.f90` needs to be changed. System-specific copies can be kept in separate folders. Before a run, copy the desired system file to the project root as `system_parameters.f90`.

## References and citation

If you use this code in your research, please cite the relevant papers:

- **Model description:** [del Palacio et al. (2018), A&A, 617, A13](https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..13D/abstract)
- **Model update: adiabatic FS modelling:** [Martinez et al. (2022), A&A, 661, A102](https://ui.adsabs.harvard.edu/abs/2022A%26A...661A.102M/abstract)
- **Model update: Radiative double-layer FS model:** [Martinez et al. (2023), A&A, 680, A99](https://ui.adsabs.harvard.edu/abs/2023A%26A...680A..99M/abstract)
