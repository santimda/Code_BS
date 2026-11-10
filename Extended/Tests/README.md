# Tests of the Extended model

Run from `Extended/`:

```
make test                                     # everything (Tests/test_set.sh), ~7 min with 2 cores
sh Tests/test_set.sh --only BD43 HVS          # all the unit tests and the regression, some physical cases
python3 Tests/test_terminal_output.py         # only the regression
python3 Tests/test_physics.py --only HVS      # only some physical cases
```

`make test` needs this target in the `Makefile` of `Extended/`:

```
test:
	sh Tests/test_set.sh
```

Requirements: `gfortran`, `make`, Python 3 with `numpy`. **The sources in `Extended/` are
not modified**: each test copies the code to a temporary directory, changes the copy and
compiles it with the `Makefile` of `Extended/` (with `LIBS=`: lapack and blas are not
used). The number of OpenMP threads is `OMP_NUM_THREADS`, or `nth` if it is not set.

The structure follows the tests of the binaries code (`CWB_non-thermal_emission/tests`),
plus the physical checks.

## 1. Unit tests (Fortran)

- `test_numerical_utilities.f90`: `vector_log`, `integrate`, `unit_vector_diff`, `safe_acos_dot`.
- `test_global_functions.f90`: normalization of `nBB`, `F_1..F_4(0) = pi^2/6`, `sigma_pp`,
  continuity of `lambda_cool`, f-f power of the hot layer <= `L_inj`.

## 2. Regression (`test_terminal_output.py`)

Runs the BD+43 case and compares the values printed in the terminal (luminosities, apex
values, radio fluxes, maps) with `reference_terminal_output.json` (relative tolerance 2e-3).
The run uses a radiative FS (`radiative_FS = .true.` in the copy of `main.f90`).

`Tests/cases/BD43/` has frozen copies of `system_parameters.f90` and `controles.f90`; the
grid of `global.f90` is fixed in `GRID`. So the reference does not depend on the system or
the options being used in `Extended/`. If a new parameter is added to `system_parameters.f90`
or `controles.f90`, it must also be added to the copies in `Tests/cases/BD43/`.

If the results change on purpose (a bug fix, a new physical ingredient), rewrite the
reference with `python3 Tests/test_terminal_output.py --update-reference` and commit it
together with the change.

## 3. Physical checks (`test_physics.py`)

| case | configuration | checks |
|---|---|---|
| `BD43` | default (RS + radiative FS, maps, convolutions, H-alpha) | energy budgets; SED summed over the linear emitters = SED of `CD_radsyn`; `S_abs <= S_unabs`; sum of the maps = total flux; flux conserved in the convolutions; thermodynamics; no NaN |
| `termo_RH` | Rankine-Hugoniot conditions in the adiabatic shocks (in the copy of `CD_termo.f90`; same `L_inj`) | luminosities and radio fluxes within a factor 2 of Martinez et al. 2022 |
| `HVS` | B2-60-d (`Tests/cases/B2-60-d`), R0 ~ 0.5 AU, adiabatic FS | `0.1 <= L_IC,st/L_inj,e <= 1`; energy budgets of both shocks |
| `HVS_IR` | B2-60-d with `L_st x 1e10` (only the dust IR field changes) | radiated <= injected; electron transport <= 1.3 |
| `dense_ISM` | BD+43 with `n_ISM = 6e6`, adiabatic FS, protons only | p-p losses >= 0.3 `L_inj,p`; `L_pp <= L_inj,p/3`; energy budgets |
| `one_zone` | Extended and One-Zone with the BD+43 parameters and the same injection (power law, index 2); adiabatic FS | `L_syn`, `S(1.4 GHz)`, `L_br`, `L_pp`, `L_IC` (RS) and `B` within a factor 10 |
| `adiabatic_FS` | BD+43 with an adiabatic FS (`radiative_FS = .false.`), coarse grids, all the modules | the run (maps and convolutions) finishes; flux conserved in the convolutions; no NaN |
| `openmp` | coarse grids, all the modules | same results with 1 and 3 OpenMP threads; run compiled with `DEBUG=Y` |

Energy budgets (each shock with that process):

- free-free: `L_ff >= 0`; `L_ff <= L_inj` (adiabatic shock); `L_ff(hot layer) + L_ff(isothermal layer) <= L_inj` (radiative FS).
- electrons: `L_syn + L_br + L_IC,IR + L_IC,st <= L_inj,e`; transport in `CD_diste`:
  `(sum N E/t_cool + energy advected out of the last cell)/L_inj,e` in [0.9, 1.3];
  `L_syn` and `L_IC,IR` agree with `sum N E/t_sy` and `sum N E/t_IC,IR` within 25%
  (`L_IC,st` within a factor 2: the emission is anisotropic, the cooling isotropic).
- protons: `L_pp <= L_inj,p`; `L_pp <=` energy lost in p-p collisions; transport in `CD_distp`
  (lower bound 0.7 in `dense_ISM`: protons that cool to 1 GeV of total energy leave the grid
  with their rest-mass energy).
- injection: `L_inj(RS) <= L_w`; `L_inj(FS) <= 0.5 rho_ISM v^3 pi y_max^2`.

## Checks that fail with the current code

`CD_diste` and `CD_distp` multiply the density by `E_dot(index2)/E_dot(index1)` when the
particles change of energy bin, which overestimates the number of particles when they cool
(see the proposed changes):

- `HVS`: electron transport 1.34 (1.07 with the proposed changes).
- `adiabatic_FS`: the run stops when the maps are saved (`merge_map` in `main.f90`): with an
  adiabatic FS there is no hot layer, but `abs_maps` and `SW_conv_run` read its files.
- `dense_ISM`: proton transport ~100 and p-p losses ~6 `L_inj,p` (protons piled up below the
  p-p threshold, where `E_dot` drops by ~1e4).
