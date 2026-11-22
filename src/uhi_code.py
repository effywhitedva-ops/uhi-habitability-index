r"""UHI numerical core inherited from v5.2 and reused by the v6.1 model.

The Unified Habitability Index (UHI, concept by Effy Viridis) is a weighted geometric mean
of seven dimensionless sub-indices: S_fis (size/mass), S_clima (climate), S_atm (atmosphere
retention), S_estrella (stellar activity), S_geo (geology), S_mag (magnetic shielding)
and S_obs (observability). Weights are sampled from Dirichlet(2,2,2,2,1,1,2), and catalogue
uncertainties are propagated by Monte Carlo; medians and 68/95% intervals are reported.

Specification items taken directly from the UHI v2.0 documents: the geometric-mean form,
the seven components and Dirichlet prior, Kepler's third law, S = L/a^2, equilibrium temperature,
transit depth, fHZ from a simple energy-balance model with a Kopparapu HZ prior, TSM/ESM-based
atmospheric observability (Kempton et al. 2018), activity-age stellar stress, tidal heating,
a broad magnetic prior, logistic/trapezoid normalization and T0-T3 tiers. All other thresholds,
slopes, constants and tier rules are documented assumptions (lines labelled ASSUMPTION).

Version history: v2 adopted a stellar-Teff-dependent Bond albedo (Shields et al. 2013); v3 added
probabilistic detection and a cosmic shoreline; v4 added the ELT/PCS reflected-light model for
non-transiting planets (Kasper et al. 2021); v5 separated retention from instrumentation and added
an exact zero-score veto; v5.1 added the Goldblatt et al. (2013) 282 W m^-2 runaway-greenhouse
limit; v5.2 introduced a continuous Yang et al. (2013) albedo ramp and an ESI-exponent erratum.

Units: catalogue masses/radii in Earth units, stars in solar units, a in AU, P in days, ages in
Gyr, fluxes S in Earth units unless named otherwise, constants and internal physics in SI.
v6.1 calls simulate_planet only in compatibility tests; production uses uhi_code_v6.py, which
imports most helpers here and retains this module's random-draw order."""
import zlib
import numpy as np
import pandas as pd

# ---------------- physical constants (SI) ----------------
G = 6.674e-11; SIGMA = 5.670374e-8; KB = 1.380649e-23; MH = 1.6735575e-27
MSUN = 1.98847e30; RSUN = 6.957e8; LSUN = 3.828e26; AU = 1.495978707e11
MEARTH = 5.9722e24; REARTH = 6.371e6; DAY = 86400.0; YR = 365.25 * DAY
S0_FLUX = 1361.0                     # W/m^2 at 1 AU
F_XUV_EARTH = 4.64e-3                # W/m^2, present solar XUV at 1 AU (~4.64 erg/cm2/s, Ribas+2005)  # ASSUMPTION (literature value)
SUB = ["S_fis", "S_clima", "S_atm", "S_estrella", "S_geo", "S_mag", "S_obs"]
DIRICHLET_ALPHA = np.array([2, 2, 2, 2, 1, 1, 2], float)   # from the Toolkit, section 7
TIER_AVAIL = {"T0": 1.0, "T1": 0.75, "T2": 0.5, "T3": 0.25}  # ASSUMPTION
# Observability vector per sub-index (Annex A of the UHI v2.0 document -> assigned tier)  # ASSUMPTION in the exact mapping
TIER_VECTOR = "fis:T0|clima:T1|atm:T1-T2|estrella:T0|geo:T3*|mag:T3*|obs:T0"


# ---------------- utilities ----------------
def logistic_log(x, x0, k=4.0):
    r"""Map positive x to a dimensionless logistic score in log10 space.

    The output is 0.5 at x0 (same units as x); k=4 changes it from about 0.02 to 0.98 across +/-1 dex.
    The 1e-30 floor avoids log of zero. Shape and steepness are assumptions."""
    x = np.maximum(x, 1e-30)
    return 1.0 / (1.0 + np.exp(-k * (np.log10(x) - np.log10(x0))))


def trapezoid(x, a, b, c, d):
    r"""Return a dimensionless trapezoid score: 0 outside [a,d], 1 on [b,c], linear on both flanks.

    x and the breakpoints share units. It is used where the specification calls for bounded
    normalization; individual breakpoints are assumptions unless cited at the call site."""
    y = np.zeros_like(x, dtype=float)
    up = (x > a) & (x < b); y[up] = (x[up] - a) / (b - a)
    y[(x >= b) & (x <= c)] = 1.0
    dn = (x > c) & (x < d); y[dn] = (d - x[dn]) / (d - c)
    return y


def split_normal(rng, val, e1, e2, n, frac_default, lo=None, hi=None):
    r"""Draw n values from a split normal using catalogue asymmetric uncertainties.

    val, e1 (positive error) and e2 (negative error) are in the catalogue column's units. Missing or
    zero errors use frac_default * |val|; optional lo/hi clip the draw. rng order is relied upon by
    v6.1 for exact reproduction. Split normals are a pragmatic representation of asymmetric errors,
    not the published posteriors, and parameter correlations are not available in the snapshot."""
    val = float(val)
    s_up = abs(e1) if np.isfinite(e1) and e1 != 0 else abs(val) * frac_default
    s_dn = abs(e2) if np.isfinite(e2) and e2 != 0 else abs(val) * frac_default
    z = rng.standard_normal(n)
    x = np.where(z >= 0, val + z * s_up, val + z * s_dn)
    if lo is not None or hi is not None:
        x = np.clip(x, lo if lo is not None else -np.inf, hi if hi is not None else np.inf)
    return x


# Chen & Kipping (2017) 'Forecaster' M-R (deterministic form + intrinsic scatter)
_CK_T1, _CK_T2 = 2.04, 131.6                    # M_earth (2.04 Me; 0.414 MJ)
_CK_R1 = 1.008 * _CK_T1 ** 0.279               # radius at transition 1
_CK_R2 = _CK_R1 * (_CK_T2 / _CK_T1) ** 0.589   # radius at transition 2


def ck_radius(m, rng):
    r"""Draw radius in R_E from mass in M_E using Chen & Kipping (2017) Forecaster pieces.

    The deterministic broken power law is multiplied by fractional intrinsic scatter, clipped to
    0.5-1.5. rng supplies the normal scatter. Used only in v5.2 legacy and some prior tests."""
    r = np.where(m < _CK_T1, 1.008 * m ** 0.279,
         np.where(m < _CK_T2, _CK_R1 * (m / _CK_T1) ** 0.589, _CK_R2 * (m / _CK_T2) ** -0.044))
    sig = np.where(m < _CK_T1, 0.0403, np.where(m < _CK_T2, 0.146, 0.0737))
    return r * (1 + sig * rng.standard_normal(m.shape)).clip(0.5, 1.5)


def ck_mass(r, rng):
    r"""Draw mass in M_E from radius in R_E using inverted Chen & Kipping (2017) relations.

    Radius scatter is propagated through the local power-law exponent. For the degenerate Jovian
    branch, mass is drawn log-uniformly from the transition to 4000 M_E (assumption)."""
    m = np.where(r < _CK_R1, (r / 1.008) ** (1 / 0.279),
         np.where(r < _CK_R2, _CK_T1 * (r / _CK_R1) ** (1 / 0.589), np.nan))
    sig = np.where(r < _CK_R1, 0.0403 / 0.279, 0.146 / 0.589)   # radius scatter propagated to mass
    m = m * np.exp(sig * rng.standard_normal(r.shape))
    giant = ~np.isfinite(m)                                       # degenerate Jovian branch  # ASSUMPTION
    m[giant] = np.exp(rng.uniform(np.log(_CK_T2), np.log(4000.0), giant.sum()))
    return m


def kopparapu_seff(teff):
    r"""Return Kopparapu et al. (2014) 1-M_E HZ instellations in Earth flux units.

    teff is K, clipped to the 2600-7200 K fit range. Keys rv, rg, mg and em denote recent Venus,
    runaway greenhouse, maximum greenhouse and early Mars. These are boundary priors, not proof of habitability."""
    T = np.clip(teff, 2600, 7200) - 5780.0
    co = {"rv": (1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15),
          "rg": (1.107, 1.332e-4, 1.580e-8, -8.308e-12, -1.931e-15),
          "mg": (0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16),
          "em": (0.320, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16)}
    return {k: s + a * T + b * T ** 2 + c * T ** 3 + d * T ** 4 for k, (s, a, b, c, d) in co.items()}


def fhz_ebm(S, albedo, A=203.3, D=0.649):
    r"""Return the fraction of planetary area with 0-100 deg C surface temperature.

    S is Earth-normalized instellation, albedo dimensionless, A the OLR intercept in W m^-2, and D
    heat transport in W m^-2 C^-1. This is an annual-mean two-mode 1-D EBM (North 1975; North,
    Cahalan & Coakley 1981) with Earth obliquity s2=-0.482, B=2.09 W m^-2 C^-1 and constant albedo.
    It has no ice-albedo feedback, so it is an inexpensive climate proxy only."""
    B, s2 = 2.09, -0.482
    Q = S * S0_FLUX / 4.0
    absorbed = Q * (1 - albedo)
    T0 = (absorbed - A) / B
    T2 = s2 * absorbed / (B + 6 * D)
    # T(x) = T0 + T2*(3x^2-1)/2 is linear in u = x^2 (x = sin(latitude); area is uniform in x)
    c0 = T0 - T2 / 2.0; c1 = 1.5 * T2
    with np.errstate(divide="ignore", invalid="ignore"):
        u1 = (0.0 - c0) / c1; u2 = (100.0 - c0) / c1
    lo = np.clip(np.minimum(u1, u2), 0, 1); hi = np.clip(np.maximum(u1, u2), 0, 1)
    f = np.sqrt(hi) - np.sqrt(lo)
    flat = np.abs(c1) < 1e-9
    f = np.where(flat, ((c0 > 0) & (c0 < 100)).astype(float), f)
    return np.clip(f, 0, 1)


# ---------------- v2: Bond albedo dependent on stellar Teff ----------------
# Shields et al. 2013, Astrobiology 13:715-739, Table 2, column "Earth-like gases, clouds, Rayleigh
# scattering on" (64% cloud cover): broadband planetary albedo for ocean and land surfaces,
# computed with the SEDs of HD128167 (F2V), the Sun (G2V), HD22049 (K2V) and AD Leo (M3V).
SH13_OCEAN = {"M": 0.23372, "K": 0.30235, "G": 0.31948, "F": 0.32865}
SH13_LAND = {"M": 0.33165, "K": 0.40133, "G": 0.41484, "F": 0.41428}
# Teff values assigned to those stars (NOT taken from Shields+2013)  # ASSUMPTION
#   AD Leo 3477 K, HD22049 (eps Eri) 5020 K (NASA Exoplanet Archive), Sol 5772 K, HD128167 (sigma Boo) 6700 K
ALB_TEFF = np.array([3477.0, 5020.0, 5772.0, 6700.0])
ALB_OCEAN = np.array([SH13_OCEAN[k] for k in "MKGF"])
ALB_LAND = np.array([SH13_LAND[k] for k in "MKGF"])
F_LAND_MAX = 0.5   # land fraction ~ U(0, 0.5): from ocean world to ~1.7x Earth's (0.29)  # ASSUMPTION


def albedo_teff(teff, u, f_land_max=None, shift=0.0, late_m_scale=1.0):
    r"""Return Bond albedo as a function of stellar Teff in K and a U(0,1) land draw.

    Ocean and land albedos are from Shields et al. (2013), Table 2 (Earth-like gases, clouds and
    Rayleigh scattering). The temperatures assigned to their four spectral templates, linear
    interpolation, endpoint extension and land fraction up to 0.5 are assumptions. shift and
    late_m_scale are sensitivity-only perturbations. Output is clipped to 0-0.9."""
    fl = F_LAND_MAX if f_land_max is None else f_land_max
    a_oc = np.interp(teff, ALB_TEFF, ALB_OCEAN)
    a_la = np.interp(teff, ALB_TEFF, ALB_LAND)
    A = a_oc + fl * u * (a_la - a_oc)
    # v3 (sensitivity only): late_m_scale < 1 lowers albedo below 3477 K; shift = cloud-driven offset
    A = np.where(teff < ALB_TEFF[0], A * late_m_scale, A) + shift
    return np.clip(A, 0.0, 0.9)


def planck_ratio(T1, T2, lam=7.5e-6):
    r"""Return the 7.5-micron Planck ratio B(T1)/B(T2) used by the Kempton et al. (2018) ESM.

    T1 and T2 are temperatures in K; lam is wavelength in meters. The result is dimensionless."""
    h, c, k = 6.62607015e-34, 2.99792458e8, KB
    return np.expm1(h * c / (lam * k * T2)) / np.expm1(h * c / (lam * k * T1))


def tsm_scale(rp):
    r"""Return the Kempton et al. (2018) radius-binned TSM scale factor for radius rp in R_E."""
    return np.where(rp < 1.5, 0.190, np.where(rp < 2.75, 1.26, np.where(rp < 4.0, 1.28, 1.15)))


def tsm_threshold(rp):
    r"""Return Kempton et al. (2018) recommended TSM thresholds: 10 below 1.5 R_E, else 90."""
    return np.where(rp < 1.5, 10.0, 90.0)   # thresholds recommended by Kempton et al. 2018


ESM_THR = 7.5
ELT_IWA_AS = 2 * 1.0e-6 / 39.0 * 206265.0   # 2 lambda/D, lambda=1 um, D=39 m -> ~0.0106"  # (v1-v3; unused since v4)
CONTRAST_REF = 1e-8                          # ELT reference contrast in reflected light  # (v1-v3; unused since v4)



# ---------------- v3: S_atm and S_obs with a probabilistic detection model ----------------
from scipy.special import ndtr as PHI   # Phi = standard normal CDF

V3 = dict(
    snr_thr=5.0,      # SNR detection threshold (5-sigma convention; 3 sigma is the alternate scenario)  # ASSUMPTION
    sig_model=0.5,    # model uncertainty of the predicted SNR (ln units)                     # ASSUMPTION
    cloud_min=0.5,    # cloud/haze suppression factor ~ U(cloud_min, 1) in transmission      # ASSUMPTION
    t_prog=10.0,      # reference program hours for S_atm (Kempton: TSM ~ S/N in 10 h with NIRISS)
    n_ecl=2.0,        # secondary eclipses for S_atm (Kempton: "at least two eclipses" for S/N>10 on GJ 1132b)
    floor_ppm=20.0,   # systematic noise floor (Greene et al. 2016: 20 ppm NIRISS SOSS)
    n_H=2.0,          # number of scale heights spanned by spectral features              # ASSUMPTION
    n_bins=10.0,      # independent spectral channels over which the noise floor averages  # ASSUMPTION
    ntr_cap=10.0,     # maximum number of transits for S_obs                                # ASSUMPTION
    shore_logc=0.0,   # log10 factor on the Mars-calibrated cosmic shoreline               # ASSUMPTION
    shore_sig=1.0,    # width (ln) of the shoreline transition                             # ASSUMPTION
    # --- v4: non-transiting planets, reflected light with ELT/PCS (XAO + coronagraph + HDS) ---
    t_hci=10.0,       # observing hours (Kasper+2021: 10 h example; Snellen+2015: 1 night = 10 h)
    n_ph_j26=1.0e5,   # photons in 10 h from a J=26 object (Kasper+2021: ELT, J band 50 nm, 10% throughput)
    h_15=5.0e-5,      # residual XAO halo contrast at the smallest separation (Kasper+2021: ~5e-5)
    h_ratio=10.0,     # halo improvement from 15 to 100 mas (ASSUMPTION: equal to the factor 10 in the PCS requirement 1e-8 -> 1e-9)
    iwa_lo=2.0, iwa_hi=4.0,   # IWA ~ U(iwa_lo, iwa_hi) lambda/D per draw, lambda=0.75 um (ASSUMPTION; 15 mas ~ 4 lambda/D)
    ag_fac=2.0 / 3.0, # A_g = (2/3) A_Bond (Lambertian sphere)
)
LAM_HCI = 0.75e-6                                   # centre of 0.6-0.9 um (Snellen+2015); O2 A band ~0.76 um
LD_MAS = LAM_HCI / 39.0 * 206264.806 * 1e3          # ELT lambda/D (D = 39 m) in mas (~3.97)


def halo_contrast(sep_mas, P=V3):
    r"""Return residual stellar halo contrast for ELT/PCS observations.

    sep_mas is angular separation in milliarcseconds. Kasper et al. (2021) suggests ~5e-5 at the
    smallest separations and a tenfold improvement over 15-100 mas. The power-law interpolation and
    flat extension outside that range are assumptions; output is a dimensionless contrast."""
    s = np.clip(sep_mas, 15.0, 100.0)
    slope = np.log10(P["h_ratio"]) / np.log10(100.0 / 15.0)
    return P["h_15"] * (s / 15.0) ** (-slope)


def snr_reflected(rp, a_au, dist_pc, J, albedo_bond, u_iwa, P=V3):
    r"""Return (reflected-light SNR, planet contrast, separation in mas) for an ELT/PCS proxy.

    rp is R_E, a_au AU, dist_pc pc, J magnitude, albedo_bond dimensionless, u_iwa U(0,1). The model
    uses halo-photon-limited high-dispersion detection (Kasper et al. 2021), Lambert phase 1/pi at
    maximum elongation and geometric albedo (2/3)A_B. SNR is zero inside a sampled 2-4 lambda/D IWA."""
    sep_mas = a_au / dist_pc * 1e3
    contrast = P["ag_fac"] * albedo_bond / np.pi * (rp * REARTH / (a_au * AU)) ** 2
    nstar = P["n_ph_j26"] * 10 ** (-0.4 * (J - 26.0)) * P["t_hci"] / 10.0
    snr = contrast * np.sqrt(nstar / halo_contrast(sep_mas, P))
    iwa = (P["iwa_lo"] + (P["iwa_hi"] - P["iwa_lo"]) * u_iwa) * LD_MAS
    return np.where(sep_mas >= iwa, snr, 0.0), contrast, sep_mas


MARS_I, MARS_VESC, EARTH_VESC = 1 / 1.524 ** 2, 5.03, 11.186   # insolation (S_Earth) and v_esc (km/s)


def p_retain(S, mp, rp, P=V3):
    r"""Return v3/v5 cosmic-shoreline retention probability.

    S is Earth-normalized instellation; mp/rp are Earth units. Following Zahnle & Catling (2017),
    critical insolation scales as v_esc^4. The shoreline is normalized through Mars and shifted by
    shore_logc; P = Phi(ln(I_shore/I)/shore_sig). v6.1 uses lifetime XUV ARM instead in production."""
    vesc = EARTH_VESC * np.sqrt(np.maximum(mp, 1e-6) / rp)
    i_shore = 10 ** P["shore_logc"] * MARS_I * (vesc / MARS_VESC) ** 4
    return PHI(np.log(i_shore / np.maximum(S, 1e-9)) / P["shore_sig"])


def atm_obs_v3(transit, rp, mp, rs, teff, a_au, per_d, teq, g, J, K, S, avail, dist_pc, albedo, u_iwa, z1, z2, z3, u_cl, u_mu, P=V3):
    r"""Return (S_atm, S_obs, TSM, ESM, P_retain) for the legacy v3 observability model.

    Inputs follow atm_obs_v5 units. It combines retention with detection in S_atm, so it is retained
    only to document history and is not called by v6.1 production."""
    teq0 = teff * np.sqrt(rs * RSUN / (a_au * AU)) * 0.25 ** 0.25
    tsm = tsm_scale(rp) * rp ** 3 * teq0 / (np.maximum(mp, 0.05) * rs ** 2) * 10 ** (-J / 5)
    delta = (rp * REARTH / (rs * RSUN)) ** 2
    esm = 4.29e6 * planck_ratio(1.10 * teq0, teff) * delta * 10 ** (-K / 5)
    m1, m2, m3 = np.exp(P["sig_model"] * z1), np.exp(P["sig_model"] * z2), np.exp(P["sig_model"] * z3)
    cl = P["cloud_min"] + (1 - P["cloud_min"]) * u_cl
    # --- S_atm (transiting): P(retain) x P(detect an atmosphere with TSM or ESM)
    tsm_hab = np.where((rp >= 1.5) & (rp < 2.0), tsm * 2.3 / 18.0, tsm)   # Kempton Sec. 3.1 rescaling (1.5-2.0 R_E)
    snr_t = tsm_hab * np.sqrt(P["t_prog"] / 10.0) * m1 * cl
    snr_e = esm * np.sqrt(P["n_ecl"]) * m2
    p_det_tr = 1 - (1 - PHI(snr_t - P["snr_thr"])) * (1 - PHI(snr_e - P["snr_thr"]))
    # --- S_atm (non-transiting): ELT reflected light
    # v4: physical ELT/PCS model (see snr_reflected); non-finite distance -> 0
    d_ok = np.isfinite(dist_pc) & (dist_pc > 0)
    snr_r0, _, _ = snr_reflected(rp, a_au, np.where(d_ok, dist_pc, 1.0), J, albedo, u_iwa, P)
    snr_r = np.where(d_ok, snr_r0, 0.0) * m1
    p_det_nt = PHI(snr_r - P["snr_thr"])
    pret = p_retain(S, mp, rp, P)
    S_atm = pret * np.where(transit, p_det_tr, p_det_nt)
    # --- S_obs (transiting): secondary (N2/CO2) or H2 atmosphere, N_tr transits, noise floor
    mu = np.where(rp < 1.6, np.where(u_mu < 0.5, 28.0, 44.0), 2.3)
    mu_k = np.where(rp < 1.5, 18.0, 2.3)               # mu assumed by Kempton/Louie in each bin
    snr10 = tsm * mu_k / mu                              # S/N in 10 h (5 h in transit) for that mu
    H = KB * teq / (mu * MH * g)
    A_sig = P["n_H"] * 2 * H * rp * REARTH / (rs * RSUN) ** 2
    t14_h = per_d * 24 / np.pi * np.arcsin(np.clip(rs * RSUN / (a_au * AU), 0, 1))
    ntr = np.minimum(365.25 / per_d, P["ntr_cap"])
    # total photon noise implied by the Kempton S/N, scaled to the available in-transit time
    sig_ph = A_sig / np.maximum(snr10, 1e-12) * np.sqrt(5.0 / np.maximum(ntr * t14_h, 1e-6))
    # systematic floor per channel (Greene+2016; does not drop with more transits), averaged over n_bins channels
    sig_fl = P["floor_ppm"] * 1e-6 / np.sqrt(P["n_bins"])
    sig = np.sqrt(sig_ph ** 2 + sig_fl ** 2)
    snr_o = A_sig / sig * m3 * cl
    snr_o = np.where(transit, snr_o, snr_r)
    S_obs = PHI(snr_o - P["snr_thr"]) * avail
    return S_atm, S_obs, tsm, esm, pret


# ======================= v5 (Neon / R-Aetherios audit) =======================
EPS = 1e-12            # floor before logarithms; S_i <= EPS is treated as an exact veto (UHI = 0)
V3.update(
    # log-sigmoid mu(Rp) (radius-valley centre: Fulton et al. 2017, deficit at 1.5-2.0 R_E)
    R0=1.65, R0_sd=0.08,   # centre and scatter of the centre (Neon spec values; 1.65 lies inside 1.5-2.0)  # ASSUMPTION
    dR=0.12,               # transition width [R_E] (Neon spec)                                            # ASSUMPTION
    mu_prim=2.3,           # primordial H2/He envelope
    rho_ovr=5.5,           # override: rho > 5.5 g/cc with M/sigma_M > 3 -> w_H2 = 0 (Neon spec)           # ASSUMPTION
    # climate of synchronous M dwarfs (Teff < 3800 K)
    teff_md=3800.0,
    yang_eff_lo=0.5,       # efficacy of the substellar-cloud albedo increase ~ U(yang_eff_lo, 1)            # ASSUMPTION
                           # (Kopparapu+2016: faster rotation of late M dwarfs weakens the cloud effect)
    kin_lo=1.5, kin_hi=2.0,  # factor on the 1-D inner edge (Kopparapu+2014 rg, rv) ~ U(kin_lo, kin_hi):
                           # Yang+2013 "twice the stellar flux"; Kopparapu+2016 "up to ~25% lower" -> 2 x 0.75 = 1.5
    s_eye_max=0.8,         # the outer-edge floor applies for S_eff < 0.8 (Neon spec)
    lpco2_lo=-3.0, lpco2_hi=2.0,   # prior on available CO2 pressure: log10(P/bar) ~ U(lo, hi)                  # ASSUMPTION
    f_eye_c=0.30, f_eye_w=0.15,    # habitable fraction of the eyeball/dense-CO2 floor ~ U(c-w, c+w)           # ASSUMPTION
    # v5.1: outgoing thermal radiation limit of a moist atmosphere (Goldblatt, Robinson, Zahnle & Crisp 2013,
    # Nature Geosci 6, 661-667, doi:10.1038/ngeo1892; abstract: "thermal radiation limit of 282 W m-2").
    # If absorbed flux S*S0/4*(1-A) > F_rg -> runaway greenhouse -> fHZ = 0 (veto).
    # HYBRID variant (Effy's decision): only NON-synchronous planets or Teff >= teff_md; synchronous M planets
    # keep the 3-D inner edge (Yang+2013 x k_in), because the 1-D limit with /4 redistribution does not apply to them.
    F_rg=282.0, rg_veto=True,
    # v5.2: ramp for the Yang+2013 weight (Neon v5.2 spec) and Teff boundary of the veto (= teff_md by default)
    yang_ramp_lo=800.0, yang_ramp_hi=1000.0, teff_rg=3800.0,
)
# Yang, Cowan & Abbot 2013 (ApJL 771, L45), Table 1, group 3 (CAM3, 1:1 rotation, 3400 K star):
# planetary albedo vs S0 [W/m2]. v5.1: only S0 >= 1000 (step); v5.2: weight w_eff ramps over 800-1000 W/m2; above 2200
# 0.53 is kept (ASSUMPTION); below 1000 W/m2 the Shields albedo is kept (no extrapolation).
YANG_S0 = np.array([1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0])
YANG_ALB = np.array([0.42, 0.45, 0.47, 0.49, 0.50, 0.51, 0.53])
# Turbet et al. 2018 (A&A 612, A86), synchronous TRAPPIST-1 planets (S in S_E, from their Table 1):
#  - S >= 0.662 (TRAPPIST-1e; Proxima b ~890 W/m2): substellar liquid water "whatever its atmosphere" (if water-rich)
#  - TRAPPIST-1f (0.382): CO2 > ~1 bar (less -> collapse); TRAPPIST-1g (0.258): CO2 > ~5 bar
#  - outer edge with 70 bar of CO2: S_eff = 0.225; never below that.
TURBET_S_EYE = 0.662
TURBET_S = np.array([0.225, 0.258, 0.382])
TURBET_LOGP = np.log10(np.array([70.0, 5.0, 1.0]))


def p_co2_req(S):
    r"""Return CO2 pressure in bar needed for surface liquid water on a synchronous M-dwarf planet.

    S is Earth-normalized instellation. Anchors are from Turbet et al. (2018), Table 1: none at or
    above 0.662 S_E, about 1 bar at 0.382, 5 bar at 0.258, and 70 bar at 0.225. Interpolation choices
    are assumptions; below 0.225 S_E, the required pressure is infinite."""
    S = np.asarray(S, float)
    lin = 1.0 * (TURBET_S_EYE - S) / (TURBET_S_EYE - TURBET_S[-1])
    logp = 10 ** np.interp(S, TURBET_S, TURBET_LOGP)
    return np.where(S >= TURBET_S_EYE, 0.0, np.where(S >= TURBET_S[-1], lin, np.where(S >= TURBET_S[0], logp, np.inf)))


def climate_v5(S, teff, locked, alb_sh, A_olr, u_eff, u_kin, u_pco2, u_eye, P=V3, D=0.649, phz=(1.0, 0.5, 0.1)):
    r"""Return v5.2 climate quantities: (Bond albedo, fHZ, HZ prior, eyeball flag).

    S is Earth-normalized flux, teff K, locked Boolean, alb_sh dimensionless, A_olr W m^-2, and u_*
    U(0,1) draws. For synchronous M dwarfs, Yang et al. (2013) raises albedo through a continuous
    800-1000 W m^-2 ramp and k_in shifts Kopparapu inner edges. Turbet et al. (2018) supplies a
    thick-CO2/eyeball floor. Other planets receive the Goldblatt et al. (2013) runaway veto."""
    md = locked & (teff < P["teff_md"])
    s0 = S * S0_FLUX
    a_y = np.interp(s0, YANG_S0, YANG_ALB)
    eff = P["yang_eff_lo"] + (1.0 - P["yang_eff_lo"]) * u_eff
    # v5.2: continuous linear ramp w_eff(S0) = clip((S0 - 800)/200, 0, 1) instead of the step at S0 = 1000 W/m2
    # (v5.1: md & (s0 >= YANG_S0[0])). Below 1000 W/m2 a_y stays at the first tabulated value (0.42, np.interp
    # does not extrapolate) and its weight falls linearly to 0 at 800 W/m2: the Yang+2013 regime is not extrapolated.
    w_eff = np.clip((s0 - P["yang_ramp_lo"]) / (P["yang_ramp_hi"] - P["yang_ramp_lo"]), 0.0, 1.0)
    alb = np.where(md, alb_sh + w_eff * eff * (a_y - alb_sh), alb_sh)
    fhz = fhz_ebm(S, alb, A_olr, D)
    k = kopparapu_seff(teff)
    kin = np.where(md, P["kin_lo"] + (P["kin_hi"] - P["kin_lo"]) * u_kin, 1.0)
    cons = (S <= k["rg"] * kin) & (S >= k["mg"])
    opt = (S <= k["rv"] * kin) & (S >= k["em"])
    p_hz = np.where(cons, phz[0], np.where(opt, phz[1], phz[2]))
    pco2 = 10 ** (P["lpco2_lo"] + (P["lpco2_hi"] - P["lpco2_lo"]) * u_pco2)
    eye = md & (S < P["s_eye_max"]) & (pco2 >= p_co2_req(S))
    f_eye = np.clip(P["f_eye_c"] + P["f_eye_w"] * (2 * u_eye - 1), 0.0, 1.0)
    fhz = np.where(eye, np.maximum(fhz, f_eye), fhz)
    if P.get("rg_veto", True):   # v5.1: hybrid runaway-greenhouse veto (Goldblatt+2013)
        md_rg = locked & (teff < P.get("teff_rg", P["teff_md"]))   # v5.2: separable veto boundary (1-D sweep)
        fhz = np.where((~md_rg) & (S * S0_FLUX / 4.0 * (1 - alb) > P["F_rg"]), 0.0, fhz)
    return alb, fhz, p_hz, eye


def mu_v5(rp, rho_gcc, dens_ok, u_mu, z_R0, P=V3):
    r"""Return mean molecular weight and H2 weight for legacy v5 transmission calculations.

    rp is R_E, rho_gcc density in g cm^-3, dens_ok whether mass is precise enough, u_mu U(0,1),
    and z_R0 N(0,1). A log-sigmoid transition near the radius valley (Fulton et al. 2017) mixes
    secondary 28/44 and primordial 2.3 molecular weights. Dense precise planets are forced to w_H2=0."""
    R0 = P["R0"] + P["R0_sd"] * z_R0
    w = 1.0 / (1.0 + np.exp(-(rp - R0) / P["dR"]))
    w = np.where(dens_ok & (rho_gcc > P["rho_ovr"]), 0.0, w)
    mu_sec = np.where(u_mu < 0.5, 28.0, 44.0)
    return np.exp((1 - w) * np.log(mu_sec) + w * np.log(P["mu_prim"])), w


def atm_obs_v5(transit, rp, mp, rs, teff, a_au, per_d, teq, g, J, K, S, avail, dist_pc, albedo, u_iwa,
               z1, z2, z3, u_cl, u_mu, z_R0, dens_ok, P=V3):
    r"""Return (S_atm, S_obs, TSM, ESM, P_retain, w_H2) for the v5.2 observability model.

    Units: rp/mp Earth, rs solar radii, teff/teq K, a_au AU, per_d days, g m s^-2, J/K magnitudes,
    S Earth flux, avail dimensionless, dist_pc pc. S_atm is now retention only. S_obs uses JWST
    transmission or emission for transiting planets (either detection channel counts), and ELT/PCS
    reflected light otherwise. Model-error, cloud and noise-floor terms are assumptions from V3."""
    teq0 = teff * np.sqrt(rs * RSUN / (a_au * AU)) * 0.25 ** 0.25
    tsm = tsm_scale(rp) * rp ** 3 * teq0 / (np.maximum(mp, 0.05) * rs ** 2) * 10 ** (-J / 5)
    delta = (rp * REARTH / (rs * RSUN)) ** 2
    esm = 4.29e6 * planck_ratio(1.10 * teq0, teff) * delta * 10 ** (-K / 5)
    m1, m2, m3 = np.exp(P["sig_model"] * z1), np.exp(P["sig_model"] * z2), np.exp(P["sig_model"] * z3)
    cl = P["cloud_min"] + (1 - P["cloud_min"]) * u_cl
    pret = p_retain(S, mp, rp, P)
    S_atm = pret * 1.0                                   # f_flare = 1
    # transmission
    rho = 5.513 * mp / rp ** 3
    mu, w = mu_v5(rp, rho, dens_ok, u_mu, z_R0, P)
    mu_k = np.where(rp < 1.5, 18.0, 2.3)                 # mu assumed by Kempton/Louie in each bin (TSM normalization)
    snr10 = tsm * mu_k / mu
    H = KB * teq / (mu * MH * g)
    A_sig = P["n_H"] * 2 * H * rp * REARTH / (rs * RSUN) ** 2
    t14_h = per_d * 24 / np.pi * np.arcsin(np.clip(rs * RSUN / (a_au * AU), 0, 1))
    ntr = np.minimum(365.25 / per_d, P["ntr_cap"])
    sig_ph = A_sig / np.maximum(snr10, 1e-12) * np.sqrt(5.0 / np.maximum(ntr * t14_h, 1e-6))
    sig_fl = P["floor_ppm"] * 1e-6 / np.sqrt(P["n_bins"])
    snr_tr = A_sig / np.sqrt(sig_ph ** 2 + sig_fl ** 2) * m3 * cl
    snr_em = esm * np.sqrt(P["n_ecl"]) * m2
    p_tr = 1 - (1 - PHI(snr_tr - P["snr_thr"])) * (1 - PHI(snr_em - P["snr_thr"]))
    # reflected light (non-transiting)
    d_ok = np.isfinite(dist_pc) & (dist_pc > 0)
    snr_r0, _, _ = snr_reflected(rp, a_au, np.where(d_ok, dist_pc, 1.0), J, albedo, u_iwa, P)
    snr_r = np.where(d_ok, snr_r0, 0.0) * m1
    p_nt = PHI(snr_r - P["snr_thr"])
    S_obs = np.where(transit, p_tr, p_nt) * avail
    return S_atm, S_obs, tsm, esm, pret, w


def aggregate_v5(Smat, W):
    r"""Return per-draw UHI from a (7,n) component matrix and (7,n) weights.

    Scores above 1 are capped, values at or below EPS veto the draw exactly, and positive draws use
    exp(sum_i w_i log S_i). This veto is why P_alive is reported separately from conditional score."""
    alive = np.all(Smat > EPS, axis=0)
    uhi = np.exp((W * np.log(np.maximum(np.minimum(Smat, 1.0), EPS))).sum(0))
    return np.where(alive, uhi, 0.0)

# ---------------- catalogue preparation ----------------
def prepare(df):
    r"""Filter and annotate the NASA Exoplanet Archive snapshot for simulation.

    A row needs stellar Teff/radius/mass, an orbital period or semimajor axis, planetary radius or mass,
    and J or K magnitude. Missing J uses K (and vice versa), an explicit assumption. The function adds
    radius_measured, mass_mode and dens_ok. Calculated radii are not treated as measured; transiting
    Msini values are treated as true masses because sin i is near one. dens_ok requires M/sigma_M>3."""
    d = df.copy()
    d["J_use"] = d["sy_jmag"].fillna(d["sy_kmag"])       # ASSUMPTION: J ~ K when J is missing
    d["K_use"] = d["sy_kmag"].fillna(d["sy_jmag"])
    have_star = d.st_teff.notna() & d.st_rad.notna() & d.st_mass.notna()
    have_orbit = d.pl_orbper.notna() | d.pl_orbsmax.notna()
    have_size = d.pl_rade.notna() | d.pl_bmasse.notna()
    have_mag = d.J_use.notna()
    ok = have_star & have_orbit & have_size & have_mag
    d = d[ok].reset_index(drop=True)
    calc_r = d.pl_rade_reflink.fillna("").str.contains("CALCULATED", case=False)
    d["radius_measured"] = (d.tran_flag == 1) & d.pl_rade.notna() & ~calc_r
    d["mass_mode"] = np.where(d.pl_bmasse.isna(), "mr",
                       np.where(d.pl_bmassprov.isin(["Mass", "Msin(i)/sin(i)"]), "measured",
                       np.where(d.pl_bmassprov == "Msini", np.where(d.tran_flag == 1, "measured", "msini"), "mr")))
    # transiting with Msini: sin i ~ 1, treated as a measured mass  # ASSUMPTION
    sig_m = 0.5 * (d.pl_bmasseerr1.abs() + d.pl_bmasseerr2.abs())
    d["dens_ok"] = (d.mass_mode == "measured") & (d.pl_bmasse / sig_m > 3.0)   # v5: M/sigma_M > 3
    return d


def nominal_tier(row):
    r"""Assign a heuristic observing tier T0-T3 from nominal catalogue values.

    T0 requires a transiting planet with a precise mass and a Kempton-threshold TSM/ESM. T1 needs
    at least half the threshold. T2 covers larger transiting planets with a weaker metric or a
    non-transiting nominal ELT/PCS SNR of at least 10% of threshold. T3 is the remainder. The rules
    and availability values (1.0/0.75/0.5/0.25) are assumptions, used only as an S_obs multiplier."""
    rp, rs, teff, ms = row.rp_nom, row.st_rad, row.st_teff, row.st_mass
    a_au = row.a_nom
    teq0 = teff * np.sqrt(rs * RSUN / (a_au * AU)) * 0.25 ** 0.25
    if row.tran_flag == 1:
        mp = max(row.mp_nom, 0.05)
        tsm = tsm_scale(rp) * rp ** 3 * teq0 / (mp * rs ** 2) * 10 ** (-row.J_use / 5)
        esm = 4.29e6 * planck_ratio(1.10 * teq0, teff) * (rp * REARTH / (rs * RSUN)) ** 2 * 10 ** (-row.K_use / 5)
        x = max(tsm / tsm_threshold(rp), esm / ESM_THR)
        mrel = np.nan
        if row.mass_mode == "measured" and np.isfinite(row.pl_bmasse):
            e1 = abs(row.pl_bmasseerr1) if np.isfinite(row.pl_bmasseerr1) else np.inf
            e2 = abs(row.pl_bmasseerr2) if np.isfinite(row.pl_bmasseerr2) else np.inf
            mrel = 0.5 * (e1 + e2) / row.pl_bmasse
        if row.mass_mode == "measured" and np.isfinite(mrel) and mrel <= 0.5 and x >= 1:
            return "T0"
        if x >= 0.5:
            return "T1"
        if rp >= 1.6 and x >= 0.1:
            return "T2"
        return "T3"
    else:
        if not np.isfinite(row.sy_dist) or row.sy_dist <= 0:
            return "T3"
        # v4: nominal ELT/PCS SNR (A_B=0.3, IWA=2 lambda/D) >= 0.1 x threshold, analogous to the transiting "metric>=0.1"  # ASSUMPTION
        snr, _, _ = snr_reflected(rp, a_au, row.sy_dist, row.J_use, 0.3, 0.0)
        if snr >= 0.1 * V3["snr_thr"]:
            return "T2"
        return "T3"


# ---------------- per-planet Monte Carlo ----------------
def simulate_planet(r, n, rng, alpha=None, keep_draws=False):
    r"""Simulate one planet using the v5.2 model and return summary statistics.

    r is a prepared catalogue row, n the number of draws, rng the catalogue-order generator, alpha an
    optional Dirichlet concentration vector, and keep_draws stores positive score draws for ranking.
    The v6.1 compatibility test calls the v6 wrapper rather than this function, but uhi_code_v6.py
    reproduces its random-number order exactly. Numbered comments explain each model stage."""
    out = {}
    # --- star ---
    teff = split_normal(rng, r.st_teff, r.st_tefferr1, r.st_tefferr2, n, 0.02, 2300, 50000)
    rs = split_normal(rng, r.st_rad, r.st_raderr1, r.st_raderr2, n, 0.07, 0.05, None)
    ms = split_normal(rng, r.st_mass, r.st_masserr1, r.st_masserr2, n, 0.07, 0.05, None)
    L = rs ** 2 * (teff / 5772.0) ** 4                       # L* from R* and Teff (Toolkit table 0.1)
    if np.isfinite(r.st_age):
        age = split_normal(rng, r.st_age, r.st_ageerr1, r.st_ageerr2, n, 0.3, 0.01, 13.8)
    else:
        age = rng.uniform(1.0, 10.0, n)                       # ASSUMPTION: unknown age ~ U(1,10) Gyr
    # --- orbit ---
    if np.isfinite(r.pl_orbper):
        P = split_normal(rng, r.pl_orbper, r.pl_orbpererr1, r.pl_orbpererr2, n, 1e-4, 1e-3, None)
        a = ((G * ms * MSUN) / (4 * np.pi ** 2)) ** (1 / 3) * (P * DAY) ** (2 / 3) / AU   # Kepler (Toolkit 1)
    else:
        a = split_normal(rng, r.pl_orbsmax, r.pl_orbsmaxerr1, r.pl_orbsmaxerr2, n, 0.05, 1e-4, None)
        P = np.sqrt(a ** 3 / ms) * 365.25
    if np.isfinite(r.pl_orbeccen):
        e = split_normal(rng, r.pl_orbeccen, r.pl_orbeccenerr1, r.pl_orbeccenerr2, n, 0.0, 0.0, 0.95)
        if not (np.isfinite(r.pl_orbeccenerr1) or np.isfinite(r.pl_orbeccenerr2)):
            e = np.clip(r.pl_orbeccen + 0.05 * rng.standard_normal(n), 0, 0.95)   # ASSUMPTION
    else:
        e = np.zeros(n)                                       # ASSUMPTION: unknown e = 0
    # --- planet: mass and radius ---
    if r.radius_measured:
        rp = split_normal(rng, r.pl_rade, r.pl_radeerr1, r.pl_radeerr2, n, 0.10, 0.1, None)
    else:
        rp = None
    if r.mass_mode == "measured":
        mp = split_normal(rng, r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, n, 0.20, 0.01, None)
    elif r.mass_mode == "msini":
        msini = split_normal(rng, r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, n, 0.20, 0.01, None)
        cosi = rng.uniform(0, 1, n)                           # isotropic inclination (non-transiting)
        mp = msini / np.sqrt(np.clip(1 - cosi ** 2, 1e-4, 1))
        mp = np.minimum(mp, 26600.0)
    else:
        mp = None
    if rp is None and mp is None:
        rp = split_normal(rng, r.pl_rade, r.pl_radeerr1, r.pl_radeerr2, n, 0.15, 0.1, None)
    if mp is None:
        mp = ck_mass(rp, rng)
    if rp is None:
        rp = ck_radius(mp, rng)
    # --- flux and temperatures (Toolkit 1) ---
    S = L / a ** 2
    alb_sh = albedo_teff(teff, rng.uniform(0.0, 1.0, n))      # v2: Teff-dependent albedo (Shields+2013, Table 2); f=1
    # v5: extra draws from a separate generator (seed = crc32 of the name + 1) so the main sequence is unchanged
    rng5 = np.random.default_rng(zlib.crc32(str(r.pl_name).encode()) + 1)
    u5 = rng5.uniform(size=(4, n)); z_R0 = rng5.standard_normal(n)
    locked = (r.st_mass < 0.6) & (np.median(a) < 0.2)       # ASSUMPTION (synchronous-rotation rule of v1-v4)
    A_olr = rng.uniform(150.0, 203.3, n)                      # ASSUMPTION: uncertain greenhouse (Earth -> enhanced CO2)
    albedo, fhz, p_hz, eye = climate_v5(S, teff, locked, alb_sh, A_olr, u5[0], u5[1], u5[2], u5[3])
    f_red = 1.0
    teq = ((1 - albedo) * L * LSUN / (16 * np.pi * SIGMA * (a * AU) ** 2 * f_red)) ** 0.25
    delta = (rp * REARTH / (rs * RSUN)) ** 2
    g = G * mp * MEARTH / (rp * REARTH) ** 2

    # S_fis: rocky size/mass (trapezoids)  # ASSUMPTION thresholds
    f_r = trapezoid(rp, 0.3, 0.5, 1.6, 2.5)
    f_m = trapezoid(mp, 0.03, 0.1, 5.0, 10.0)
    S_fis = np.sqrt(f_r * f_m)

    # S_clima = fHZ(EBM, + v5 synchronous M-dwarf climate) * Kopparapu prior (v5: shifted inner edge)
    S_clima = fhz * p_hz

    # v3: probabilistic S_atm and S_obs (see atm_obs_v3)
    tr = (r.tran_flag == 1)
    z = rng.standard_normal((3, n)); u_cl = rng.uniform(size=n); u_mu = rng.uniform(size=n)
    # v4: per-draw IWA from a separate generator (seed = crc32 of the name) so the main sequence is unchanged
    u_iwa = np.random.default_rng(zlib.crc32(str(r.pl_name).encode())).uniform(size=n)
    S_atm, S_obs, tsm, esm, pret, w_h2 = atm_obs_v5(tr, rp, mp, rs, teff, a, P, teq, g, r.J_use, r.K_use, S,
                                               TIER_AVAIL[r.tier], r.sy_dist, albedo, u_iwa, z[0], z[1], z[2], u_cl, u_mu,
                                               z_R0, bool(r.dens_ok))

    # S_estrella: XUV (activity-age) + flare risk
    tsat = 0.1 if r.st_mass >= 0.6 else 1.0                    # ASSUMPTION
    rx = np.where(age < tsat, 10 ** -3.5, 10 ** -3.5 * (age / tsat) ** -1.23)   # Ribas+2005 slope -1.23
    rx_sun = 10 ** -3.5 * (4.6 / 0.1) ** -1.23
    fxuv_rel = rx * S / rx_sun                                 # F_XUV in units of present Earth
    stress = logistic_log(fxuv_rel, 10.0)                      # ASSUMPTION: 50% stress at 10x Earth
    flare = 0.5 if r.st_mass < 0.35 else (0.25 if r.st_mass < 0.6 else 0.0)   # ASSUMPTION
    S_estrella = (1 - stress) * (1 - flare)

    # S_geo: internal heat flux (radiogenic proxy + tides) -> regime; erosion as a sub-Neptune flag
    n_orb = 2 * np.pi / (P * DAY)
    k2Q = 0.3 / 100.0                                          # ASSUMPTION k2=0.3, Q=100
    Edot = 10.5 * k2Q * G * (ms * MSUN) ** 2 * (rp * REARTH) ** 5 * n_orb * e ** 2 / (a * AU) ** 6
    q_tide = Edot / (4 * np.pi * (rp * REARTH) ** 2)
    q_rad = 0.087 * (mp / rp ** 2) * 2.0 ** (-(age - 4.5) / 3.0)   # ASSUMPTION
    lq = np.log10(np.maximum(q_rad + q_tide, 1e-6))
    S_regime = trapezoid(lq, -2.0, np.log10(0.03), 0.0, 1.0)  # ASSUMPTION: optimum 0.03-1 W/m2
    fxuv_abs = fxuv_rel * F_XUV_EARTH
    mdot = 0.1 * np.pi * (rp * REARTH) ** 3 * fxuv_abs / (G * mp * MEARTH)   # eps=0.1, K=1  # ASSUMPTION
    t_loss = 0.01 * mp * MEARTH / np.maximum(mdot, 1e-30) / (1e9 * YR)       # envelope = 1% of M_p  # ASSUMPTION
    risk = 1 - logistic_log(t_loss, 1.0)
    subnep = (rp >= 1.6) & (rp <= 4.0)
    S_geo = S_regime * np.where(subnep, 1 - 0.5 * risk, 1.0)

    # S_mag: broad Beta prior with indirect updates  # ASSUMPTION
    aB, bB = (2.0, 3.0) if locked else (2.0, 2.0)
    aB = np.full(n, aB); bB = np.where(mp < 0.3, bB + 1.0, bB)
    S_mag = rng.beta(aB, bB)

    # UHI aggregation with per-draw Dirichlet weights (Toolkit 7)
    Smat = np.vstack([S_fis, S_clima, S_atm, S_estrella, S_geo, S_mag, S_obs])
    W = rng.dirichlet(DIRICHLET_ALPHA if alpha is None else alpha, n).T
    uhi = aggregate_v5(Smat, W)                                # v5: EPS floor; S_i <= EPS => UHI = 0 exactly (veto)

    # ESI (Schulze-Makuch et al. 2011)
    rho = mp / rp ** 3                      # density relative to Earth
    vesc = np.sqrt(mp / rp)                 # relative escape speed
    Ts = teq + 33.0                         # ASSUMPTION: T_surface = T_eq + 33 K (Earth's greenhouse effect)
    # v5.2 (erratum): Schulze-Makuch+2011: ESI_I = prod_{R,rho} term^(w/n) with n = 2 per group, ESI = sqrt(ESI_I*ESI_S)
    # => total exponent w/4 (same as the 4-property form with n = 4). v1-v5.1 used w/4 INSIDE each group
    # plus the square root (total exponent w/8), which compressed ESI toward 1. Does not affect UHI (comparison only).
    def esi_term(x, x0, w):
        r"""Return one ESI factor (Schulze-Makuch et al. 2011) for value x, Earth reference x0 and weight w."""
        return (1 - np.abs((x - x0) / (x + x0))) ** (w / 2.0)
    esi_i = esi_term(rp, 1.0, 0.57) * esi_term(rho, 1.0, 1.07)
    esi_s = esi_term(vesc, 1.0, 0.70) * esi_term(Ts, 288.0, 5.58)
    esi = np.sqrt(esi_i * esi_s)

    q = lambda x, p: float(np.nanpercentile(x, p))
    out.update(UHI_med=q(uhi, 50), UHI_p16=q(uhi, 16), UHI_p84=q(uhi, 84), UHI_p2_5=q(uhi, 2.5), UHI_p97_5=q(uhi, 97.5),
               UHI_mean=float(np.mean(uhi)), P_UHI_gt0=float(np.mean(uhi > 0)),
               P_alive=float(np.mean(uhi > 0)),
               UHI_cond_med=q(uhi[uhi > 0], 50) if (uhi > 0).any() else 0.0,
               UHI_cond_q25=q(uhi[uhi > 0], 25) if (uhi > 0).any() else 0.0,
               UHI_cond_q75=q(uhi[uhi > 0], 75) if (uhi > 0).any() else 0.0,
               albedo_med=q(albedo, 50), f_eye=float(np.mean(eye)), wH2_med=q(w_h2, 50),
               ESI_med=q(esi, 50), ESI_p16=q(esi, 16), ESI_p84=q(esi, 84))
    for name, arr in zip(SUB, Smat):
        out[name + "_med"] = q(arr, 50)
    out.update(fHZ_med=q(fhz, 50), Rp_med=q(rp, 50), Mp_med=q(mp, 50), S_med=q(S, 50), Teq_med=q(teq, 50),
               a_med=q(a, 50), FXUV_rel_med=q(fxuv_rel, 50), TSM_med=q(tsm, 50) if r.tran_flag == 1 else np.nan, P_retain_med=q(pret, 50),
               ESM_med=q(esm, 50) if r.tran_flag == 1 else np.nan)
    out["UHI_rank"] = out["P_alive"] * out["UHI_cond_med"]     # v5: hurdle, UHI_rank = P_alive x conditional median
    if keep_draws and out["P_alive"] > 0:
        out["_draws"] = uhi.astype(np.float32)
    return out


def run(df, n_draws=10000, seed=42, progress=None, alpha=None, joint=True):
    r"""Run the legacy v5.2 catalogue simulation and return a ranked DataFrame.

    df is the raw Archive table, n_draws draws per planet, seed the main RNG seed, progress an optional
    callback, alpha optional Dirichlet weights and joint enables per-draw ranks. Ranking uses the v5.2
    convention: dead draws are placed after survivors and count toward P_top10 if within ten places."""
    d = prepare(df)
    # nominal values for the tier
    a_nom = np.where(d.pl_orbper.notna(),
                     ((G * d.st_mass * MSUN) / (4 * np.pi ** 2)) ** (1 / 3) * (d.pl_orbper * DAY) ** (2 / 3) / AU,
                     d.pl_orbsmax)
    d["a_nom"] = a_nom
    rng0 = np.random.default_rng(0)
    mp_nom = d.pl_bmasse.copy()
    rp_nom = d.pl_rade.copy()
    need_m = mp_nom.isna()
    if need_m.any():
        mp_nom[need_m] = np.where(rp_nom[need_m] < _CK_R1, (rp_nom[need_m] / 1.008) ** (1 / 0.279),
                                  _CK_T1 * (rp_nom[need_m] / _CK_R1) ** (1 / 0.589))
    need_r = rp_nom.isna()
    if need_r.any():
        m = mp_nom[need_r].values
        rp_nom[need_r] = np.where(m < _CK_T1, 1.008 * m ** 0.279, np.where(m < _CK_T2, _CK_R1 * (m / _CK_T1) ** 0.589, _CK_R2 * (m / _CK_T2) ** -0.044))
    d["mp_nom"] = mp_nom; d["rp_nom"] = rp_nom
    d["tier"] = [nominal_tier(row) for row in d.itertuples()]
    rng = np.random.default_rng(seed)
    rows = []
    for i, row in enumerate(d.itertuples()):
        res = simulate_planet(row, n_draws, rng, alpha=alpha, keep_draws=joint)
        res.update(planet=row.pl_name, host=row.hostname, method=row.discoverymethod, transiting=int(row.tran_flag),
                   tier=row.tier, tier_vector=TIER_VECTOR, mass_source=row.mass_mode,
                   controversial=int(row.pl_controv_flag) if np.isfinite(row.pl_controv_flag) else 0,
                   st_teff=row.st_teff, st_mass=row.st_mass, sy_dist=row.sy_dist, age_known=int(np.isfinite(row.st_age)))
        rows.append(res)
        if progress and i % 500 == 0:
            progress(i, len(d))
    draws = {i: r.pop("_draws") for i, r in enumerate(rows) if "_draws" in r}
    res = pd.DataFrame(rows)
    res["P_top10"] = 0.0; res["rank_med"] = np.nan; res["rank_ci5"] = np.nan; res["rank_ci95"] = np.nan
    if joint and draws:
        idx = np.array(sorted(draws)); M = np.vstack([draws[i] for i in idx])          # (surviving planets, draws)
        npos = (M > 0).sum(0)                                                           # survivors per draw
        order = np.argsort(-M, axis=0, kind="stable"); rk = np.empty_like(order)
        np.put_along_axis(rk, order, np.arange(1, M.shape[0] + 1)[:, None].repeat(M.shape[1], 1), axis=0)
        rk = np.where(M > 0, rk, npos[None, :] + 1)                                     # UHI=0 -> rank = survivors+1
        res.loc[idx, "P_top10"] = (rk <= 10).mean(1)
        res.loc[idx, "rank_med"] = np.median(rk, 1)
        res.loc[idx, "rank_ci5"] = np.percentile(rk, 5, 1); res.loc[idx, "rank_ci95"] = np.percentile(rk, 95, 1)
    res["ESI_rank"] = res.ESI_med.rank(ascending=False, method="min").astype("Int64")
    # v5: rank by UHI_rank = P_alive x conditional median; ties broken by the posterior mean
    res = res.sort_values(["UHI_rank", "UHI_mean"], ascending=False).reset_index(drop=True)
    res["rank_pos"] = np.arange(1, len(res) + 1)
    first = ["rank_pos", "planet", "UHI_rank", "P_alive", "UHI_cond_med", "UHI_cond_q25", "UHI_cond_q75", "P_top10",
             "rank_med", "rank_ci5", "rank_ci95", "UHI_med", "UHI_p16", "UHI_p84", "UHI_p2_5", "UHI_p97_5", "UHI_mean",
             "P_UHI_gt0", "tier", "ESI_med", "ESI_rank"]
    return res[first + [c for c in res.columns if c not in first]].sort_values("rank_pos")


if __name__ == "__main__":
    import sys, time
    df = pd.read_csv(sys.argv[1]); n = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    t = time.time(); out = run(df, n, progress=lambda i, N: print(i, N, flush=True))
    out.to_csv(sys.argv[3] if len(sys.argv) > 3 else "uhi_results.csv", index=False)
    print("planets", len(out), "seconds", round(time.time() - t, 1))
