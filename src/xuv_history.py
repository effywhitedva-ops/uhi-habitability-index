r"""Stellar X-ray/EUV history, cumulative exposure, and atmosphere-retention helpers.

The main outputs are lifetime fluence in erg cm^-2, present luminosities in erg s^-1,
and dimensionless retention probabilities. Stellar luminosity enters in solar units,
orbital distance in AU, ages in Gyr, planetary radii/masses in Earth units, and escape
speeds in km s^-1. Formulas are intentionally compact proxies for population-level
prioritization; they are not a hydrodynamic atmosphere-evolution code.

Primary sources:
Pass, Charbonneau & Vanderburg (2025, ApJL 986, L3; doi:10.3847/2041-8213/adda39),
equations 6-11 and Table 1, define the fully convective saturation fluence, Earth's
reference XUV fluence (2.0e18 erg cm^-2) and atmospheric-retention metric (ARM).
Johnstone, Bartel & Gudel (2021, A&A 649, A96) provides Mors X-ray/EUV tracks
with Spada et al. (2013) stellar evolution for heavier stars. Jackson, Davis & Wheatley
(2012, MNRAS 422, 2024), Table 2, is an alternative age-activity relation. Sanz-Forcada
et al. (2011, A&A 532, A6) converts L_X to L_EUV. Owen & Wu (2017, ApJ 847, 29)
supplies the energy-limited efficiency scaling used for Hycean envelopes. Innes, Tsai &
Pierrehumbert (2023, ApJ 953, 168), Table 1, provides runaway limits for H2-He atmospheres.
Berta-Thompson et al. (2026, ApJ 1011, 123) motivates the optional bolometric shoreline.

All extrapolations and normalizations not stated by those papers are labelled as assumptions
here and in docs/assumptions_v6_es.md, sections D2, D4 and "fallback checks"."""
import os
import numpy as np
from scipy.special import ndtr as PHI

# Constants are cgs unless named otherwise. GYR_S converts Gyr to seconds.
LSUN_CGS = 3.828e33; AU_CM = 1.495978707e13; GYR_S = 3.15576e16
I_XUV_EARTH = 2.0e18                    # erg cm^-2 (Pass+2025)
LOG_RSAT_FC = -3.05                     # Wright+2018 via Pass+2025 (+0.05/-0.06)
PASS_M = np.array([0.10, 0.12, 0.15, 0.20, 0.25, 0.30])     # Pass+2025 Table 1
PASS_FL = np.array([1.19, 1.17, 1.13, 1.13, 1.16, 1.22])
STILL_SATURATED = {"TRAPPIST-1"}       # Pass+2025 Sec. 2.2: TRAPPIST-1 is still saturated today (P_rot 3.3 d) -> integrate up to its present age
# Jackson+2012 Table 2: (B-V)0 lo, hi, log(Lx/Lbol)_sat, log tau_sat [yr], alpha
# Each Jackson et al. (2012) row: B-V lower/upper colour edge, saturated log10(L_X/L_bol), log10 saturation time
# in years, and post-saturation power-law index.
JACKSON = np.array([[0.290, 0.450, -4.28, 7.87, 1.22], [0.450, 0.565, -4.24, 8.35, 1.24], [0.565, 0.675, -3.67, 7.84, 1.13],
                    [0.675, 0.790, -3.71, 8.03, 1.28], [0.790, 0.935, -3.36, 7.90, 1.40], [0.935, 1.275, -3.35, 8.28, 1.09],
                    [1.275, 1.410, -3.14, 8.21, 1.18]])
# Innes, Tsai & Pierrehumbert 2023, ApJ 953, 168, Table 1: runaway limit [S0] for 1 and 10 bar of H2-He
INNES_TEFF = np.array([3457.0, 5772.0])          # M star = K2-18 luminosity case (K2-18 Teff taken from the catalogue); G = Sun
INNES_S1 = np.array([0.320, 0.389]); INNES_S10 = np.array([0.0860, 0.0676])   # v6.1: M-star 10 bar = 0.0860 (published Table 1; v6.0 used the preprint value 0.0852)

_GRID = None


def mors_grid():
    r"""Load and cache data/mors_grid.npz, returning an empty dictionary if absent.

    The grid contains masses (M_sun), rotation percentiles, ages (Gyr), bolometric/X-ray/EUV
    luminosities (erg s^-1) and cumulative X-ray energy (erg). It is a precomputed product of
    build_mors_grid.py; an empty result deliberately triggers the documented Jackson fallback."""
    global _GRID
    if _GRID is None:
        p = os.path.join(os.path.dirname(__file__), "..", "data", "mors_grid.npz")
        _GRID = dict(np.load(p)) if os.path.exists(p) else {}
    return _GRID


def fl_pass(mstar):
    r"""Return the dimensionless pre-main-sequence luminosity factor f_L.

    mstar is in solar masses. Values are interpolated from Pass et al. (2025) Table 1 between
    0.10 and 0.30 M_sun and held constant outside that interval (assumption, D2)."""
    return np.interp(mstar, PASS_M, PASS_FL)   # held constant outside 0.10-0.30 M_sun (ASSUMPTION)


def tsat_pass(mstar, scale=1.0):
    r"""Return the X-ray saturation duration in Gyr for fully convective stars.

    mstar is solar mass and scale is an optional dimensionless uncertainty factor. Pass et al.
    (2025) equation 7 is calibrated over 0.1-0.3 M_sun and extrapolated by them to 0.35;
    the 0.1-M_sun lower clamp and 0.01-Gyr minimum are numerical safeguards (assumption)."""
    return np.maximum((5.9 - 15.4 * np.maximum(mstar, 0.10)) * scale, 0.01)


def fluence_fc(L, a_au, age, mstar, log_rsat, tsat, f_flare, saturated_now=False, post_drop=1e-2):
    r"""Return lifetime X-ray fluence and current X-ray luminosity for an M dwarf.

    L is solar luminosity, a_au AU, age/tsat Gyr, mstar solar mass, log_rsat log10(L_X/L_bol)
    during saturation, f_flare a dimensionless multiplier and post_drop a fraction of that
    saturated ratio after saturation ends. Return (fluence in erg cm^-2, L_X in erg s^-1).

    This implements the Pass et al. (2025) equation-10 saturated fluence. The low post-saturation
    emission is an assumption added for conservative continuity; Pass et al. treat the later
    fluence as small. Stars in STILL_SATURATED integrate saturated output until their present age."""
    Lc = L * LSUN_CGS; rs = 10 ** log_rsat; fl = fl_pass(mstar)
    t_on = np.where(saturated_now, age, np.minimum(tsat, age))
    E = f_flare * rs * fl * Lc * t_on * GYR_S + np.where(saturated_now, 0.0, post_drop * rs * Lc * np.maximum(age - tsat, 0.0) * GYR_S)
    Lx_now = np.where(saturated_now | (age < tsat), f_flare * rs * Lc, post_drop * rs * Lc)
    return E / (4 * np.pi * (a_au * AU_CM) ** 2), Lx_now


def fluence_mors(L, a_au, age, mstar, pct):
    r"""Interpolate Mors cumulative X-ray energy and current X-ray/EUV luminosity.

    L is observed bolometric luminosity (L_sun), a_au AU, age Gyr, mstar M_sun, and pct the
    rotation percentile. Return (fluence erg cm^-2, L_X erg s^-1, L_EUV erg s^-1), or None
    when the grid is missing. Values are trilinearly interpolated in mass, percentile and age,
    clipped to the grid boundaries, then rescaled by observed/track bolometric luminosity.
    That rescaling and the boundary clipping are assumptions (D2), chosen so catalogue
    luminosity is respected without extrapolating outside the tracks."""
    g = mors_grid()
    if not g:
        return None
    M, P, A = g["masses"], g["pcts"], g["ages"]
    m = np.clip(mstar, M[0], M[-1]); p = np.clip(pct, P[0], P[-1]); t = np.clip(age, A[0], A[-1])
    i = np.clip(np.searchsorted(M, m) - 1, 0, len(M) - 2); wm = (m - M[i]) / (M[i + 1] - M[i])
    j = np.clip(np.searchsorted(P, p) - 1, 0, len(P) - 2); wp = (p - P[j]) / (P[j + 1] - P[j])
    k = np.clip(np.searchsorted(A, t) - 1, 0, len(A) - 2); wt = (t - A[k]) / (A[k + 1] - A[k])

    def interp(Q):
        r"""Return one Mors grid quantity interpolated at the current draw coordinates.

        Q has axes (mass, percentile, age); the enclosing function provides the bracketing indices
        and fractional weights. Output units are those of the selected grid quantity."""
        out = 0.0
        for di, fm in ((0, 1 - wm), (1, wm)):
            for dj, fp in ((0, 1 - wp), (1, wp)):
                out = out + fm * fp * (Q[i + di, j + dj, k] * (1 - wt) + Q[i + di, j + dj, k + 1] * wt)
        return out
    lbol_tr = interp(g["Lbol"]); scale = L * LSUN_CGS / np.maximum(lbol_tr, 1e20)
    E = interp(g["CLx"]) * scale
    return E / (4 * np.pi * (a_au * AU_CM) ** 2), interp(g["Lx"]) * scale, interp(g["Leuv"]) * scale


def fluence_jackson(L, a_au, age, bv, u_scatter=None):
    r"""Return lifetime X-ray fluence and current X-ray luminosity from Jackson et al. (2012).

    L is solar luminosity, a_au AU, age Gyr, and bv intrinsic B-V colour. Table 2 gives
    saturation level, saturation time and power-law decline by colour bin. Return
    (erg cm^-2, erg s^-1). Out-of-table colours return NaN so the caller retains the previous
    model. Extinction/reddening is neglected (assumption), and u_scatter is currently unused."""
    bv = np.asarray(bv, float)
    idx = np.searchsorted(JACKSON[:, 1], bv, side="right"); ok = (bv >= 0.29) & (bv < 1.41)
    row = JACKSON[np.clip(idx, 0, 6)]
    rs = 10 ** row[..., 2]; tau = 10 ** row[..., 3] / 1e9; al = row[..., 4]
    Lc = L * LSUN_CGS
    sat = rs * Lc * np.minimum(age, tau)
    post = np.where(age > tau, rs * Lc * tau / (1 - al) * ((np.maximum(age, tau) / tau) ** (1 - al) - 1), 0.0)
    E = (sat + post) * GYR_S
    Lx_now = np.where(age < tau, rs * Lc, rs * Lc * (np.maximum(age, tau) / tau) ** -al)
    return np.where(ok, E / (4 * np.pi * (a_au * AU_CM) ** 2), np.nan), np.where(ok, Lx_now, np.nan)


_EEM = None


def bv_from_teff(teff):
    r"""Estimate intrinsic dwarf B-V colour from effective temperature in K.

    The table is E. Mamajek's mean dwarf colour/temperature sequence (Pecaut & Mamajek 2013;
    version 2022.04.16, md5 recorded in docs/assumptions_v6_es.md). It is used only by the Jackson
    variant because the catalogue snapshot lacks B and V. Main-sequence colour is assumed."""
    global _EEM
    if _EEM is None:
        rows = []
        p = os.path.join(os.path.dirname(__file__), "..", "data", "EEM_dwarf_UBVIJHK_colors_Teff.txt")
        sec = 0
        for line in open(p):
            if line.startswith("#SpT"):
                sec += 1; continue
            if sec != 1 or line.startswith("#") or not line.strip():
                continue
            f = line.split()
            try:
                rows.append((float(f[1]), float(f[8])))
            except (ValueError, IndexError):
                pass
        a = np.array(sorted(rows)); _EEM = a
    return np.interp(teff, _EEM[:, 0], _EEM[:, 1])


def leuv_from_lx(Lx):
    r"""Return L_EUV in erg s^-1 from L_X in erg s^-1 using Sanz-Forcada et al. (2011).

    The empirical cgs relation log L_EUV = 4.80 + 0.860 log L_X is clamped at L_X >= 1 erg s^-1
    only to avoid a logarithm of zero."""
    return 10 ** (4.80 + 0.860 * np.log10(np.maximum(Lx, 1.0)))   # Sanz-Forcada+2011


def arm(vesc_kms, I_rel):
    r"""Return the Pass et al. (2025) atmospheric-retention metric ARM, dimensionless.

    vesc_kms is planetary escape speed in km s^-1; I_rel is lifetime XUV fluence relative to
    Earth's 2.0e18 erg cm^-2. ARM is 4 log10(v_esc) - log10(I_rel) - 3.16. A tiny floor protects
    the logarithm; ARM itself is not a probability."""
    return 4 * np.log10(vesc_kms) - np.log10(np.maximum(I_rel, 1e-30)) - 3.16   # Pass+2025 Eq. 11


def p_ret_xuv(vesc_kms, I_rel, shore_logc=0.0, shore_sig=1.0):
    r"""Map ARM to a dimensionless probability of retaining an atmosphere.

    shore_logc shifts the cosmic shoreline in dex; shore_sig is the natural-log width.
    The normal CDF form and the 0/1.0 central values inherited from v5.2 are assumptions.
    Sobol analysis shows these two parameters dominate the rocky v6.1 ranking variance."""
    return PHI((arm(vesc_kms, I_rel) + shore_logc) * np.log(10.0) / shore_sig)


def p_ret_bt3d(vesc_kms, S_rel, L_rel, shore_sig=1.0, p=5.9, q=1.17, logL_hz=-2.23):
    r"""Return the optional bolometric-shoreline retention probability.

    vesc_kms is km s^-1, S_rel is Earth-normalized bolometric instellation, and L_rel is
    solar luminosity. p and q are the slopes adopted for log S_crit; the normalization places
    an Earth-like HZ planet at the reference late-M luminosity. This variant is used only for
    sensitivity testing and is not the production shoreline model."""
    const = np.log10(1.0) - p * np.log10(11.186) - q * logL_hz
    log_scrit = p * np.log10(vesc_kms) + q * np.log10(L_rel) + const
    return PHI((log_scrit - np.log10(np.maximum(S_rel, 1e-30))) * np.log(10.0) / shore_sig)


def s_crit_innes(teff, p_h2_bar):
    r"""Return the Hycean runaway instellation limit S_crit in Earth units.

    teff is K; p_h2_bar is H2-He envelope pressure in bar. Innes et al. (2023) Table 1
    supplies 1- and 10-bar limits for an M-star case and the Sun. log S is linearly
    interpolated in log pressure and in Teff, then held constant outside the anchors; pressure
    extrapolation beyond 1-10 bar is an explicit assumption."""
    s1 = np.interp(teff, INNES_TEFF, INNES_S1); s10 = np.interp(teff, INNES_TEFF, INNES_S10)
    return 10 ** (np.log10(s1) + (np.log10(s10) - np.log10(s1)) * np.log10(p_h2_bar))


def envelope_t_loss_gyr(rp_re, mp_me, p_h2_bar, fxuv_cgs, eta_fac=1.0):
    r"""Return the energy-limited time to remove a thin H2 envelope, in Gyr.

    rp_re/mp_me are Earth radius/mass; p_h2_bar bar; fxuv_cgs present X+EUV flux in erg cm^-2
    s^-1; eta_fac is a dimensionless sensitivity multiplier. Envelope mass is 4 pi R^2 P/g.
    Mass loss uses Owen & Wu (2017) energy-limited scaling with eta = 0.1 (v_esc / 15 km s^-1)^-2.
    The result is a present-day stability timescale, not a full time-integrated envelope history."""
    G = 6.674e-11; R = rp_re * 6.371e6; M = mp_me * 5.9722e24
    g = G * M / R ** 2; menv = 4 * np.pi * R ** 2 * p_h2_bar * 1e5 / g
    vesc = np.sqrt(2 * G * M / R) / 1e3
    eta = 0.1 * (vesc / 15.0) ** -2 * eta_fac
    mdot = eta * np.pi * R ** 3 * (fxuv_cgs * 1e-3) / (G * M)      # erg/cm2/s -> W/m2 (x1e-3)
    return menv / np.maximum(mdot, 1e-30) / GYR_S
