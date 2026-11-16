"""
Plot settings that are not in the Fortran code: label of each system,
observations, upper limits, sensitivities and map geometry.

The block of SYSTEMS is chosen with 'system' of system_parameters.f90
(or --system). Systems without a block use DEFAULTS. The distance is taken from
system_parameters.f90.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OBS_DIR = os.path.join(HERE, "..", "data", "radio")
SENS_DIR = os.path.join(HERE, "..", "sensitivities")

MeV = 1.602176e-6 # erg

# Sensitivity curves, scaled with distance as 4*pi*distance^2 assuming isotropy.
# x: units of column 1 (GHz, eV, keV, MeV); y: uJy or erg (erg cm^-2 s^-1)
# text: position of the name, (log10 eps, log10 eps L_eps)
RADIO_SENS = [
    dict(file="sens_VLA.dat", label="VLA", x="GHz", y="uJy", skiprows=1, text=(-7., 26.5)),
    dict(file="sens_ngVLA.dat", label="ngVLA", x="GHz", y="uJy", skiprows=1, factor=2**-0.5, text=(-3.5, 26.)),
    dict(file="sens_SKA.dat", label="SKA", x="GHz", y="uJy", skiprows=1, factor=2**-0.5, text=(-7.5, 25.25)),
]

# Hypervelocity stars (make_plots.sh) and any system without a block
DEFAULTS = dict(label=None, distance=None, sensitivities=RADIO_SENS, sed_ul=[],
                radio_obs=[], maps=None)

#Any system can be added here
SYSTEMS = {

    "BD+43": dict(
        label=r"BD+43$^\circ$ 3654",
        sensitivities=[],
        # Upper limits as band fluxes [erg cm^-2 s^-1], E: band [eV]
        sed_ul=[
            # XMM-Newton, Gamma = 1.5 (Toala+2016; del Palacio+2018, Sect. 3.3)
            dict(label=r"\textit{XMM}", E=(0.4e3, 4.e3), F=3.6e-14, color="darkcyan"),
            # Fermi-LAT, 95% c.l., Gamma = 2 (Schulz+2014, Table 1)
            dict(label=r"\textit{Fermi}", E=(0.1e9, 0.74e9), F=1.00e-6 * MeV, color="darkorange"),
            dict(E=(0.74e9, 5.5e9), F=0.33e-6 * MeV, color="darkorange"),
            dict(E=(5.5e9, 41.e9), F=1.05e-6 * MeV, color="darkorange"),
            dict(E=(41.e9, 300.e9), F=1.19e-6 * MeV, color="darkorange"),
        ],
        # Radio observations. cols: columns of (nu [GHz], dnu, S [mJy], dS); dS can be
        # a relative error. scale: factor for S and dS; S_minus: column subtracted from S
        radio_obs=[
            dict(vline=(0.325, (1000., 1500.)), label="WENSS", color="darkcyan", marker="s"),
            dict(file="Benaglia_2010.dat", cols=(0, 1, 2, 0.3), label="B2010",
                 color="magenta", marker="^", lolims=True),
            dict(file="Benaglia_2021.dat", cols=(0, 2, 1, 3), scale=1e3, label="B2021-A",
                 color="darkgray", marker="v"),
            dict(file="Benaglia_2021_B.dat", cols=(0, 2, 1, 3), scale=1e3, label="B2021-B",
                 color="k", marker="v"),
            dict(file="Moutzouri_2022.dat", cols=(0, 1, 2, 0.05), dnu_half=True, S_minus=3,
                 label="M2022", color="green", marker="D"),
        ],
        radio_sed=dict(xlim=(0.1, 11), ylim=(80, 2500), yticks=[100, 250, 500, 1000, 2000]),
        # Maps: rotation to match the observed map, angle of the analytical R0 arrow,
        # cos(dec) for the beam width, inclination, analytical R0 (cgs, Mdot in Msun/yr)
        maps=dict(
            rot_deg=3.3, arrow_deg=351.5, cosdec=0.7195, inclination_deg=60.,
            R0=dict(v_w=3000e5, Mdot=4e-6, v_st=50e5, n_ISM=6., mu_ISM=0.615, T_ISM=8e3),
            xlim=(-10, 10), ylim=(-5, 10),
            per_freq={
                "3GHz": dict(levels=[1.04, 1.5, 2, 3, 4, 5, 6.], text="3 GHz"),
                "325MHz": dict(levels=[44., 54., 64.], levels_red=[74., 84.], text="325 MHz"),
                "avg": dict(levels=[0.3, 0.6, 0.9], levels_red=[1.8], text="4--12 GHz"),
                "150MHz": dict(levels=[10.], levels_red=[20., 40., 80., 120.], text="150 MHz"),
            },
        ),
    ),

    "BD+60": dict(
        label=r"BD+60$^\circ$ 2522",
        sensitivities=[],
        radio_sed=dict(xlim=(0.1, 11)),
        maps=dict(
            rot_deg=28., arrow_deg=28., cosdec=0.4819, inclination_deg=56.,
            R0=dict(v_w=2000e5, Mdot=1.3e-6, v_st=42e5, n_ISM=27., mu_ISM=0.615, T_ISM=8e3),
            xlim=(-2, 2), ylim=(-1.5, 1.5),
            per_freq={
                "avg": dict(levels=[0.5, 1., 1.5, 3.], text="4--12 GHz"),
                "110MHz": dict(levels=[0.5, 1, 1.75], text="110 MHz"),
            },
        ),
    ),

    "G1": dict(
        radio_obs=[
            dict(file="G1.dat", cols=(0, 1, 2, 0.05), dnu_half=True, label="VdE2022",
                 color="blue", marker="s"),
        ],
    ),

    "G3": dict(
        label="HD 319881",
        radio_obs=[
            dict(file="HD-319881.dat", cols=(0, 1, 2, 0.05), dnu_half=True, label="VdE2022",
                 color="blue", marker="s"),
        ],
    ),

    "Vela_X-1": dict(
        label="Vela X-1",
        sensitivities=[],
        radio_sed=dict(xlim=(0.1, 11)),
        maps=dict(
            rot_deg=4.6, arrow_deg=4.6, cosdec=0.759786128, inclination_deg=65.,
            R0=dict(v_w=700e5, Mdot=1e-6, v_st=90e5, n_ISM=2., mu_ISM=0.615, T_ISM=8e3),
            xlim=(-3, 3), ylim=(-2.25, 2.25),
            per_freq={},
            halpha=dict(levels=[20., 43.]),
        ),
    ),
}