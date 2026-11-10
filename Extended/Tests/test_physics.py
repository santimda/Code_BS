#!/usr/bin/env python3
"""Physical checks of the Extended bow-shock code (called by test_set.sh, make test).

The sources in Extended/ are not modified: each case is copied, changed and
compiled in a temporary directory (by default <tmp>/code_bs_Tests/<case>).

Test cases:

  BD43        BD+43 3654 with the default configuration (Tests/cases/BD43): energy
              budgets, conservation of the flux in the sums over phi and the maps,
              absorption, thermodynamics.
  termo_RH    BD+43 with Rankine-Hugoniot conditions in the adiabatic shocks
              (instead of Martinez et al. 2022): the NT emission must not change much.
  HVS         Hypervelocity star B2-60-d (Tests/cases/B2-60-d, R0 ~ 0.5 AU): the
              electrons cool fast by IC with the stellar photons; energy budgets.
  HVS_IR      B2-60-d with L_st x 1e10 (only the dust IR field changes): extreme
              cooling, the electrons are not accelerated above ~1 MeV.
  dense_ISM   Very dense ISM (adiabatic FS): protons lose their energy by p-p.
  one_zone    BD+43 with the Extended and the One-Zone codes (same injection,
              alpha = 2): the luminosities must agree within an order of magnitude.
  adiabatic_FS  BD+43 with an adiabatic FS (radiative_FS = .false.), coarse grids, all
              the modules: the run (maps and convolutions included) must finish.
  openmp      Coarse grids, all the modules: same results with 1 and several
              OpenMP threads, and a run compiled with DEBUG=Y (bounds and FPE checks).

Usage:
  python3 Tests/test_physics.py                  # all the cases
  python3 Tests/test_physics.py --only BD43 HVS  # some of them
  python3 Tests/test_physics.py --build-dir DIR  # build in DIR instead of <tmp>
"""

import argparse
import pathlib
import re
import shutil
import sys
import tempfile
import time

import numpy as np

sys.dont_write_bytecode = True          # no __pycache__ in Tests/
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import bstest as bt  # noqa: E402

TESTS = bt.TESTS
BUILD = pathlib.Path(tempfile.gettempdir()) / "code_bs_tests"      # changed with --build-dir
CASE_BD43 = TESTS / "cases" / "BD43"
SYSTEM = "BD+43"                      # 'system' of the case being checked

# Frozen BD+43 case: system parameters and controls. global.f90 and main.f90
# (which also contain code) are taken from the repository, with the grid and
# the options below fixed, so that the results do not depend on local changes.
FILES_BD43 = {"system_parameters.f90": CASE_BD43 / "system_parameters.f90",
              "controles.f90": CASE_BD43 / "controles.f90"}
PATCH_BD43 = {
    "global.f90": dict(m_theta="300", mE_e="300", mE_p="300", mE_f="500", mE_ph="500", mphi="30"),
    "main.f90": dict(radiative_FS=".true."),
}

# Only the SEDs (no absorption, maps, convolutions or H-alpha)
SEDS_ONLY = {"controles.f90": dict(maps_calc="no", abssyn_calc="no", convolution="no",
                                   rad_alpha="no", convolution_alpha="no", calc_N_bin="no")}

def n_emit(patches):
    return 4 * (int(patches["global.f90"]["mphi"]) - 1)


def m_theta(patches):
    return int(patches["global.f90"]["m_theta"])


def build_and_run(name, patches, files=FILES_BD43, threads=None, make_vars=None, replace=None):
    t0 = time.time()
    w = bt.build(BUILD / name, files=files, patches=patches, make_vars=make_vars, replace=replace)
    out = bt.run(w, threads=threads)
    print(f"  (compiled and run in {time.time() - t0:.0f} s: {w})")
    return w, out


def map_frequencies(controls):
    """{label: nu} of the radio maps of controles.f90."""
    text = pathlib.Path(controls).read_text()
    found = re.findall(r"map_freq_t\(\s*([0-9.dDeE+-]+)\s*,\s*'(\w+)'", text)
    return {label: float(nu.replace("d", "e").replace("D", "e")) for nu, label in found if label != "avg"}


# ---------------------------------------------------------------------------
# Quantities and generic checks
# ---------------------------------------------------------------------------

def luminosities(w, shock, ne):
    """Injected and radiated luminosities of one shock (whole BS) [erg/s]."""
    L = {"inj": float(np.sum(bt.load(w, f"Linj_{shock}.dat")[:, 2])) * ne}
    ff = bt.load(w, f"ionized_ff_{SYSTEM}_{shock}.dat")
    L["ff"] = float(np.sum(ff[:, 1]) * np.log(ff[1, 0] / ff[0, 0]))
    if (w / "output" / f"{shock}_ff.dat").exists():                 # hot layer (radiative FS)
        L["ff_HL"] = bt.sed_luminosity(w, f"{shock}_ff.dat")
    if (w / "output" / f"L_inj_e_{shock}.dat").exists():
        L["inj_e"] = bt.first(w, f"L_inj_e_{shock}.dat") * ne
        L["syn"] = bt.first(w, f"radsyn_tot_{shock}.dat")
        L["br"] = bt.sed_luminosity(w, f"radbr_{SYSTEM}_{shock}_tot.dat")
        L["IC_IR"] = bt.sed_luminosity(w, f"radic_ir_{SYSTEM}_{shock}_tot.dat")
        L["IC_st"] = bt.first(w, f"rad_IC_st_tot_{shock}.dat")
    if (w / "output" / f"L_inj_p_{shock}.dat").exists():
        L["inj_p"] = bt.first(w, f"L_inj_p_{shock}.dat") * ne
        L["pp"] = bt.first(w, f"radpp_tot_{shock}.dat")
    return L


def check_free_free(c, L, shock, radiative):
    c.true(L["ff"] >= 0 and L.get("ff_HL", 0.0) >= 0,
           f"{shock}: f-f luminosities >= 0 (L_ff = {L['ff']:.3e}" + (f", L_HL = {L['ff_HL']:.3e})" if "ff_HL" in L else ")"))
    if radiative:
        c.at_most(L["ff_HL"], L["inj"], f"{shock}: f-f of the hot layer <= L_inj")
        c.at_most(L["ff_HL"] + L["ff"], L["inj"], f"{shock}: f-f of the hot + isothermal layers <= L_inj")
    else:
        c.at_most(L["ff"], L["inj"], f"{shock}: f-f luminosity <= L_inj")


def check_electrons(c, w, L, shock, ne, mt, lo=0.9, consistency=True):
    """Energy budget of the electrons and consistency between the cooling
    times (CD_diste) and the emission modules."""
    L_rad = L["syn"] + L["br"] + L["IC_IR"] + L["IC_st"]
    c.at_most(L_rad, L["inj_e"], f"{shock}: L_syn + L_br + L_IC,IR + L_IC,st <= L_inj,e")

    E, N = bt.read_dist(w, f"diste_{SYSTEM}_{shock}.dat", mt)
    t = bt.read_times(w, f"CD_tcoole_{SYSTEM}_{shock}.dat", mt)
    # columns: 2 t_ic, 3 t_sy, 4 t_br, 7 t_coole, 9 t_cell, 10 t_ic_st, 11 t_ic_IR
    P_cool = bt.power(E, N, t[7], ne)
    P_out = bt.outflow(E, N, t[9][:, 0], ne)
    c.between((P_cool + P_out) / L["inj_e"], lo, 1.3,
              f"{shock}: electron transport, (cooling + advected out)/L_inj,e")
    if not consistency:
        return
    c.same_order(L["syn"], bt.power(E, N, t[3], ne), f"{shock}: L_syn vs sum N*E/t_sy", factor=1.25)
    c.same_order(L["IC_IR"], bt.power(E, N, t[11], ne), f"{shock}: L_IC,IR vs sum N*E/t_IC,IR", factor=1.25)
    c.same_order(L["IC_st"], bt.power(E, N, t[10], ne),
                 f"{shock}: L_IC,st (anisotropic) vs sum N*E/t_IC,st (isotropic)", factor=2.0)


def check_protons(c, w, L, shock, ne, mt, lo=0.9):
    c.at_most(L["pp"], L["inj_p"], f"{shock}: L_pp <= L_inj,p")
    E, N = bt.read_dist(w, f"distp_{SYSTEM}_{shock}.dat", mt)
    t = bt.read_times(w, f"CD_tcoolp_{SYSTEM}_{shock}.dat", mt)
    # columns: 3 t_pp, 6 t_coolp, 7 t_cell
    P_pp = bt.power(E, N, t[3], ne)
    c.at_most(L["pp"], P_pp, f"{shock}: L_pp <= energy lost by the protons in p-p collisions")
    P_out = bt.outflow(E, N, t[7][:, 0], ne)
    c.between((bt.power(E, N, t[6], ne) + P_out) / L["inj_p"], lo, 1.3,
              f"{shock}: proton transport, (cooling + advected out)/L_inj,p")


def check_thermodynamics(c, w, shock):
    t = bt.load(w, f"termo_{shock}.dat")
    # columns: l rho P v_t T B v_perp theta Mach r rho_med rho_IL
    c.true(np.all(t[:, [1, 2, 4, 5, 11]] > 0), f"{shock}: rho, P, T, B and rho_IL > 0")
    H = bt.load(w, f"Hnorm_{shock}.dat")[:, 1]
    c.true(np.all(H > 0), f"{shock}: width of the shocked layer H > 0 ({np.sum(H <= 0)} cells with H <= 0)")
    c.true(np.all((t[:, 9] >= 1.0) & (t[:, 9] <= 4.0 + 1e-12)),
           f"{shock}: compression ratio of the jump in [1, 4] (max {t[:, 9].max():.4f})")


def check_injection(c, w, out, L):
    c.at_most(bt.value(out, "L_inj(RS)/L_w ="), 1.0, "RS: L_inj <= L_w (wind kinetic power)")
    # FS: at most the kinetic power of the ISM through the projected area of the BS
    H = bt.load(w, "Hnorm_FS.dat")
    t = bt.load(w, "termo_FS.dat")
    y_max = H[:, 2].max()                          # R_acel = max(y, R0) [cm]
    L_kin = 0.5 * t[0, 10] * t[0, 6] ** 3 * np.pi * y_max ** 2
    c.at_most(L["FS"]["inj"], L_kin, "FS: L_inj <= 0.5*rho_ISM*v**3*pi*y_max**2")


def check_numbers(c, w):
    bad = bt.bad_numbers(w)
    c.true(not bad, "no NaN/Infinity in the output files" + (f" ({', '.join(bad[:5])})" if bad else ""))


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------

def m22_values(w, ne):
    """Luminosities and radio fluxes of the RS, compared with the RH case."""
    v = {f"L_{k}": x for k, x in luminosities(w, "RS", ne).items()}
    syn = bt.load(w, f"S_syn_{SYSTEM}_RS.dat")
    for label, nu in map_frequencies(CASE_BD43 / "controles.f90").items():
        v[f"S_{label}"] = float(bt.at_frequency(syn, nu)[1])
    return v


def case_BD43(c, args, results):
    w, out = build_and_run("BD43", PATCH_BD43)
    ne, mt = n_emit(PATCH_BD43), m_theta(PATCH_BD43)
    results["BD43"] = m22_values(w, ne)

    print(" Energy budgets")
    L = {s: luminosities(w, s, ne) for s in ("RS", "FS")}
    check_free_free(c, L["RS"], "RS", radiative=False)
    check_free_free(c, L["FS"], "FS", radiative=True)
    check_injection(c, w, out, L)
    check_electrons(c, w, L["RS"], "RS", ne, mt)
    check_protons(c, w, L["RS"], "RS", ne, mt)
    c.at_most(bt.value(out, "L_NT_(RS) ="), bt.value(out, "L_inj(RS) ="), "RS: L_NT <= L_inj")

    print(" Sums over phi, absorption and maps")
    syn = bt.load(w, f"S_syn_{SYSTEM}_RS.dat")
    sabs = bt.load(w, f"S_syn_abs_{SYSTEM}_BOTH.dat")
    radio = sabs[:, 0] < 1e12
    err = np.max(np.abs(sabs[radio, 1] / syn[radio, 1] - 1.0))
    c.true(err < 1e-3, f"unabsorbed SED summed over the linear emitters (CD_abssyn) = SED of CD_radsyn "
                       f"(max. rel. diff. {err:.1e} below 1 THz)")
    c.true(np.all(sabs[:, 2] <= sabs[:, 1] * (1 + 1e-9)), "absorbed flux <= unabsorbed flux at all frequencies")
    c.between(float(bt.at_frequency(sabs, 1e10)[2] / bt.at_frequency(sabs, 1e10)[1]), 0.99, 1.0,
              "S_abs/S_unabs at 10 GHz")
    for label, nu in map_frequencies(CASE_BD43 / "controles.f90").items():
        if f"Pre-conv:  S_sync({label})" not in out:
            continue
        pre = bt.section(out, f"Pre-conv:  S_sync({label})", "Post-conv")
        s_syn = bt.value(pre, f"S_sync({label}) =")
        s_th = bt.value(pre, f"S_ff({label})   =")
        s_hot = bt.value(pre, f"S_hot({label})  =")
        c.close(s_syn, float(bt.at_frequency(sabs, nu)[2]), f"{label}: sum of the synchrotron map = S_abs", 1e-2)
        c.close(bt.value(out, f"Post-conv: S({label})      ="), s_syn + s_th + s_hot,
                f"{label}: flux conserved in the convolution", 1.5e-2)

    print(" Thermodynamics and numbers")
    check_thermodynamics(c, w, "RS")
    check_thermodynamics(c, w, "FS")
    check_numbers(c, w)


# Rankine-Hugoniot thermodynamics in the adiabatic shocks (copy of CD_termo.f90 in the
# test build only; L_inj is not changed)
TERMO_RH = {"CD_termo.f90": [("""         rho    = rho_new
         P      = P_new
         v_t    = v_t_new
         B      = B_new
         C_s_sh = C_s_sh_new
         V_A    = v_A_new
         T      = T_new
""", """         rho    = rho_RH
         P      = P_RH
         v_t    = v_t_RH
         B      = B_RH
         C_s_sh = sqrt(g_ad*P_RH/rho_RH)
         V_A    = v_ARH
         T      = T_RH
""")]}


def case_termo_RH(c, args, results):
    patches = bt.merge(PATCH_BD43, SEDS_ONLY)
    w, out = build_and_run("termo_RH", patches, replace=TERMO_RH)
    ne, mt = n_emit(patches), m_theta(patches)
    if "BD43" in results:
        ref, src = results["BD43"], "BD43 case"
    else:
        w0, _ = build_and_run("termo_M22", patches)
        ref, src = m22_values(w0, ne), "termo_M22"
    print(f" Martinez et al. 2022 ({src}) vs Rankine-Hugoniot: within a factor 2")
    L = luminosities(w, "RS", ne)
    for k in ("syn", "IC_st", "IC_IR", "br", "pp"):
        c.same_order(L[k], ref[f"L_{k}"], f"RS: L_{k} (RH vs M22)", factor=2.0)
    syn = bt.load(w, f"S_syn_{SYSTEM}_RS.dat")
    for label, nu in map_frequencies(CASE_BD43 / "controles.f90").items():
        c.same_order(float(bt.at_frequency(syn, nu)[1]), ref[f"S_{label}"], f"RS: S_syn({label}) (RH vs M22)", 2.0)
    print(" Energy budgets")
    check_free_free(c, L, "RS", radiative=False)
    check_electrons(c, w, L, "RS", ne, mt)
    check_protons(c, w, L, "RS", ne, mt)
    check_thermodynamics(c, w, "RS")
    check_numbers(c, w)


# RS only, SEDs only, 4 rotation angles
SATURATION = bt.merge(PATCH_BD43, SEDS_ONLY, {
    "controles.f90": dict(FS="no", BOTH="no", stromgren="no", opacff_calc="no"),
    "global.f90": dict(mphi="4")})


# Hypervelocity star B2-60-d (Tests/cases/B2-60-d): v_st = 60000 km/s, R0 ~ 0.5 AU. The
# stellar photon field is so intense that most electrons cool by IC before leaving the BS.
# The FS is adiabatic (t_th > t_adv everywhere: with radiative_FS = .true. the width of the
# isothermal layer, and its f-f emission, would be negative).
FILES_HVS = {"system_parameters.f90": TESTS / "cases" / "B2-60-d" / "system_parameters.f90",
             "controles.f90": CASE_BD43 / "controles.f90"}


def case_HVS(c, args, results):
    global SYSTEM
    patches = bt.merge(PATCH_BD43, SEDS_ONLY, {"global.f90": dict(mphi="10"), "controles.f90": dict(stromgren="no"),
                                               "main.f90": dict(radiative_FS=".false.")})
    SYSTEM = "B2-60-d"
    try:
        w, out = build_and_run("HVS", patches, files=FILES_HVS)
        ne, mt = n_emit(patches), m_theta(patches)
        L = luminosities(w, "RS", ne)
        c.between(L["IC_st"] / L["inj_e"], 0.1, 1.0, "RS: L_IC,st/L_inj,e (fast IC cooling)")
        check_free_free(c, L, "RS", radiative=False)
        LF = luminosities(w, "FS", ne)
        check_free_free(c, LF, "FS", radiative=False)
        check_electrons(c, w, L, "RS", ne, mt)       # inj_p = 0 in B2-60-d: no protons
        check_electrons(c, w, LF, "FS", ne, mt)
        check_thermodynamics(c, w, "RS")
        check_thermodynamics(c, w, "FS")
        check_numbers(c, w)
    finally:
        SYSTEM = "BD+43"


def case_HVS_IR(c, args, results):
    # L_st x 1e10 only changes the dust IR field (L_IR = L_st/350, T_IR ~ L_st**(1/6)): the IC
    # losses exceed the acceleration gains above ~1 MeV, so E_max ~ E_min. Most electrons
    # cool below E_min and thermalize, so the energy balance has no lower bound here. The
    # emission formulae (synchrotron, IC) are for ultra-relativistic electrons and the
    # synchrotron of ~MeV electrons falls below the frequency grid: the consistency between
    # emission and cooling is not checked
    global SYSTEM
    patches = bt.merge(SATURATION, {"system_parameters.f90": dict(L_st="1.d10*4.d0*pi*R_st**2*sigma*T_st**4")})
    SYSTEM = "B2-60-d"
    try:
        w, out = build_and_run("HVS_IR", patches, files=FILES_HVS)
        ne, mt = n_emit(patches), m_theta(patches)
        L = luminosities(w, "RS", ne)
        check_electrons(c, w, L, "RS", ne, mt, lo=0.0, consistency=False)
        check_numbers(c, w)
    finally:
        SYSTEM = "BD+43"


def case_dense_ISM(c, args, results):
    # n_ISM x 1e6, adiabatic FS with protons only: t_pp << t_adv, so most of the energy of
    # the protons above the p-p threshold is lost in p-p collisions. At most ~1/3 of it
    # goes to gamma-rays (pi0). Note: L_inj,p includes the rest-mass energy of the protons
    patches = bt.merge(PATCH_BD43, SEDS_ONLY, {
        "controles.f90": dict(RS="no", BOTH="no", stromgren="no", opacff_calc="no", electrons="no"),
        "global.f90": dict(mphi="4"), "main.f90": dict(radiative_FS=".false."),
        "system_parameters.f90": dict(n_ISM="6.d6")})
    w, out = build_and_run("dense_ISM", patches)
    ne, mt = n_emit(patches), m_theta(patches)
    L = luminosities(w, "FS", ne)
    E, N = bt.read_dist(w, f"distp_{SYSTEM}_FS.dat", mt)
    t = bt.read_times(w, f"CD_tcoolp_{SYSTEM}_FS.dat", mt)
    c.between(bt.power(E, N, t[3], ne) / L["inj_p"], 0.3, 1.0, "FS: energy lost in p-p/L_inj,p")
    c.at_most(L["pp"], L["inj_p"] / 3.0, "FS: L_pp <= L_inj,p/3")
    # the protons cool to the bottom of the grid (1 GeV of total energy ~ rest mass) and are
    # removed with their rest-mass energy (L_inj,p uses the total energy): lower bound 0.7
    check_protons(c, w, L, "FS", ne, mt, lo=0.7)
    # with this density the FS is radiative: the code must warn that the f-f luminosity
    # of the (assumed adiabatic) shock is larger than L_inj
    c.true(L["ff"] <= L["inj"] or "the shock is not adiabatic" in out,
           f"FS: L_ff/L_inj = {L['ff'] / L['inj']:.3g}; if > 1, warning 'the shock is not adiabatic' printed")
    check_numbers(c, w)


# BD+43 parameters in One-Zone/system_parameters.f90 (the other parameters of
# the One-Zone model, e.g. R_eff, are not changed)
ONE_ZONE_BD43 = {"system_parameters.f90": dict(
    distance="1.72d0*kpc", inclination="60.d0*pi/180.d0", v_st="50.d5", T_st="4.07d4",
    R_st="19.d0*R_sun", L_st="8.7d5*L_sun", Mdot="4.d-6*M_sun_yr", v_w="3000.d5",
    n_ISM="6.d0", T_ISM="8.d3", mu_ISM="4.d0/(3.d0 + 5.d0*0.73d0 - 0.02d0)",
    frac_B="sqrt(0.2d0)", frac_NT="0.13d0", inj_e="frac_NT*0.1d0", inj_p="frac_NT*0.9d0",
    alpha_NT="2.d0",
    T_IR="27.d0*(0.2d0)**(-1.d0/6.d0)*(L_st/1.d38)**(1.d0/6.d0)*(R0/pc)**(-1.d0/3.d0)")}


def case_one_zone(c, args, results):
    # Same injection in both codes: single power law with index 2 (alpha_1 = 2 and, in the
    # copy of CD_termo.f90 of the test build, alpha_sh = 2 above E_break)
    patches = bt.merge(PATCH_BD43, SEDS_ONLY, {
        "global.f90": dict(mphi="10"), "system_parameters.f90": dict(alpha_1="2.d0"),
        "main.f90": dict(radiative_FS=".false.")})
    w, out = build_and_run("one_zone_extended", patches,
                           replace={"CD_termo.f90": [("alpha_sh = 3._dp", "alpha_sh = 2._dp")]})
    ne, mt = n_emit(patches), m_theta(patches)
    t0 = time.time()
    w1 = bt.build(BUILD / "one_zone", source=bt.ONE_ZONE, patches=ONE_ZONE_BD43)
    out1 = bt.run(w1)
    print(f"  (One-Zone compiled and run in {time.time() - t0:.0f} s: {w1})")

    print(" Extended vs One-Zone: within a factor 10")
    # (the One-Zone header also prints 'Forward shock is: ...')
    for shock, start, end in (("RS", " Reverse Shock\n", " Forward shock\n"), ("FS", " Forward shock\n", " END\n")):
        L = luminosities(w, shock, ne)
        oz = bt.section(out1, start, end)
        syn = bt.load(w, f"S_syn_{SYSTEM}_{shock}.dat")
        c.same_order(L["syn"], bt.value(oz, "L_syn(unabs) ="), f"{shock}: L_syn")
        c.same_order(float(bt.at_frequency(syn, 1.4e9)[1]), bt.value(bt.section(oz, "Unabsorbed:"), "S(1.4GHz)="),
                     f"{shock}: S_syn(1.4 GHz)")
        c.same_order(L["br"], bt.value(oz, "L_Br(tot)="), f"{shock}: L_br")
        c.same_order(L["pp"], bt.value(oz, "L_pp ="), f"{shock}: L_pp")
        if shock == "RS":       # in the FS the One-Zone photon fields use R_eff_FS = 5 R0
            c.same_order(L["IC_IR"], bt.value(oz, "L_IC =", 0), f"{shock}: L_IC,IR")
            c.same_order(L["IC_st"], bt.value(oz, "L_IC =", 1), f"{shock}: L_IC,st")
        c.same_order(float(bt.load(w, f"termo_{shock}.dat")[0, 5]), bt.value(out1, f"B_{shock}[G]="),
                     f"{shock}: B at the apex")
        print(f" Energy budgets ({shock}, Extended)")
        check_free_free(c, L, shock, radiative=False)
        check_electrons(c, w, L, shock, ne, mt)
        check_protons(c, w, L, shock, ne, mt)
    check_numbers(c, w)


def case_adiabatic_FS(c, args, results):
    # With an adiabatic FS there is no hot layer: the maps and the convolutions must not
    # read its files (abs_maps in main.f90, SW_conv_run in SW_conv.f90)
    patches = bt.merge(PATCH_BD43, {"global.f90": dict(m_theta="100", mE_e="120", mE_p="120", mE_f="150",
                                                       mE_ph="150", mphi="4"),
                                    "main.f90": dict(radiative_FS=".false.")})
    try:
        w, out = build_and_run("adiabatic_FS", patches)
    except bt.TestFailure as e:
        msg = [line.strip() for line in str(e).splitlines() if "At line" in line or "runtime error" in line]
        c.true(False, "run with an adiabatic FS, maps and convolutions: " + " | ".join(msg[-2:] or ["error"]))
        return
    c.true(True, "run with an adiabatic FS, maps and convolutions finished")
    for label in map_frequencies(CASE_BD43 / "controles.f90"):
        if f"Pre-conv:  S_sync({label})" not in out:
            continue
        pre = bt.section(out, f"Pre-conv:  S_sync({label})", "Post-conv")
        total = bt.value(pre, f"S_sync({label}) =") + bt.value(pre, f"S_ff({label})   =")
        c.close(bt.value(out, f"Post-conv: S({label})      ="), total, f"{label}: flux conserved in the convolution", 1.5e-2)
    check_numbers(c, w)


def case_openmp(c, args, results):
    # All the modules with coarse grids
    patches = bt.merge(PATCH_BD43, {"global.f90": dict(m_theta="100", mE_e="120", mE_p="120", mE_f="150",
                                                       mE_ph="150", mphi="4")})
    w = bt.build(BUILD / "small_grid_1", files=FILES_BD43, patches=patches)
    w2 = BUILD / "small_grid_n"
    shutil.rmtree(w2, ignore_errors=True)
    shutil.copytree(w, w2)
    t0 = time.time()
    bt.run(w, threads=1)
    bt.run(w2, threads=3)
    print(f"  (runs with 1 and 3 threads: {time.time() - t0:.0f} s)")
    files = sorted(p.name for p in (w / "output").glob("*.dat"))
    differ, worst = [], (0.0, "")
    for name in files:
        a, b = bt.load(w, name), bt.load(w2, name)
        if a.shape != b.shape:
            differ.append(name)
            continue
        err = np.max(np.abs(a - b) / (np.abs(a) + np.abs(b) + 1e-300)) if a.size else 0.0
        worst = max(worst, (float(err), name))
        if err > 1e-9:
            differ.append(name)
    c.true(not differ, f"OpenMP: same results with 1 and 3 threads in the {len(files)} output files "
                       f"(largest rel. diff. {worst[0]:.1e}, {worst[1]})" + (f"; differ: {differ[:6]}" if differ else ""))
    check_numbers(c, w)

    t0 = time.time()
    wd = bt.build(BUILD / "small_grid_debug", files=FILES_BD43, patches=patches, make_vars={"DEBUG": "Y"})
    try:
        bt.run(wd)
        c.true(True, "DEBUG=Y (-fcheck=all, -ffpe-trap=invalid,zero,overflow): run without errors "
                     f"({time.time() - t0:.0f} s)")
        a, b = bt.first(w, "radsyn_tot_RS.dat"), bt.first(wd, "radsyn_tot_RS.dat")
        c.close(b, a, "DEBUG=Y vs optimized (-O3 -ffast-math): L_syn(RS)", 1e-3)
    except bt.TestFailure as e:
        c.true(False, f"DEBUG=Y run: {e}")


CASES = {"BD43": case_BD43, "termo_RH": case_termo_RH, "HVS": case_HVS, "HVS_IR": case_HVS_IR,
         "dense_ISM": case_dense_ISM, "one_zone": case_one_zone, "adiabatic_FS": case_adiabatic_FS, "openmp": case_openmp}


def main():
    global BUILD
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", nargs="+", choices=list(CASES), help="run only these cases")
    parser.add_argument("--build-dir", type=pathlib.Path, help=f"directory for the test builds (default {BUILD})")
    args = parser.parse_args()
    if args.build_dir:
        BUILD = args.build_dir.resolve()
    BUILD.mkdir(parents=True, exist_ok=True)
    print(f"Test builds in {BUILD}")
    names = args.only or list(CASES)
    results, checkers, t_start = {}, [], time.time()
    for name in names:
        print(f"\n=== {name} ===")
        sys.stdout.flush()
        c = bt.Checker(name)
        checkers.append(c)
        try:
            CASES[name](c, args, results)
        except bt.TestFailure as e:
            c.true(False, f"ERROR: {e}")

    n = sum(c.n for c in checkers)
    failures = [f for c in checkers for f in c.failures]
    print(f"\n{n - len(failures)}/{n} checks passed in {time.time() - t_start:.0f} s")
    if failures:
        print("FAILED:")
        print("\n".join(f" - {f}" for f in failures))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
