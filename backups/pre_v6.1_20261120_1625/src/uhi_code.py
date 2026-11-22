"""
Version de codigo: v5.2 RC (18 nov 2026): rampa lineal w_eff(S0)=clip((S0-800)/200,0,1) para el albedo de Yang+2013 (antes escalon en S0=1000) y frontera Teff del veto como parametro teff_rg; v5.1 = uhi_code_v51.py.
v5.1 (18 nov 2026): veto hibrido de invernadero desbocado (Goldblatt+2013, 282 W/m2) en climate_v5; v5 = uhi_code_v5.py.
v4 (18 nov 2026): observabilidad de no transitantes con modelo fisico ELT/PCS (luz reflejada, Lambert, halo XAO + HDS, Kasper+2021).\nv3 (S_atm/S_obs con Phi + linea de costa cosmica) = uhi_code_v3.py;\nv2 (albedo dependiente de Teff, decidido por Effy) = uhi_code_v2.py; v1 = uhi_code_v1.py
UHI v2.0 (Unified Habitability Index, Effy Viridis) -- implementacion probabilistica de referencia.

Fuentes de la especificacion (solo lectura):
  * "Indice Unificado de Habitabilidad (UHI) - Version 2.0"  (UHI_v2_integrado1)
  * "UHI Mathematical Toolkit - T0-T2 (2025-2035)"            (UHI_T0_Toolkit)

Lo que viene TEXTUALMENTE de los documentos:
  - UHI = prod_i S_i^{w_i}, sum w_i = 1, w_i >= 0 (media geometrica ponderada).
  - 7 subindices: S_fis, S_clima, S_atm, S_estrella, S_geo, S_mag, S_obs.
  - Pesos iniciales desde un prior de Dirichlet alpha = [2,2,2,2,1,1,2].
  - Incertidumbre por Monte Carlo sobre las incertidumbres del catalogo; reporte mediana + IC 68/95 %.
  - a por 3a ley de Kepler; S = L/a^2; T_eq = ((1-A) L / (16 pi sigma a^2 f))^(1/4); delta = (Rp/R*)^2.
  - S_clima basado en fHZ (fraccion de superficie 273-373 K con un EBM 0-D/1-D) + etiqueta HZ de Kopparapu como prior.
  - S_atm/S_obs con TSM/ESM (Kempton et al. 2018), H = kT/(mu g) con escenarios de mu (N2/CO2/H2),
    SNR_tran-spec ~ (H/Rp) * delta * sqrt(N_tr) / sigma_phot;  S_obs = logistic(SNR_eff) * Availability(instrument).
  - S_estrella: F_XUV ~ L_XUV / a^2 con relaciones actividad-edad; riesgo de flares.
  - S_geo: calentamiento de marea (21/2)(k2/Q)... e^2 ; fotoevaporacion energy-limited como bandera de erosion (sub-Neptunos).
  - S_mag: prior Beta amplio, actualizado solo con restricciones indirectas.
  - Normalizacion logistica o trapezoidal a [0,1]; tiers de observabilidad T0-T3.

TODO lo demas (umbrales, pendientes, constantes, reglas de tier) son SUPUESTOS documentados en
assumptions.md y marcados abajo como  # SUPUESTO.
"""
import zlib
import numpy as np
import pandas as pd

# ---------------- constantes fisicas ----------------
G = 6.674e-11; SIGMA = 5.670374e-8; KB = 1.380649e-23; MH = 1.6735575e-27
MSUN = 1.98847e30; RSUN = 6.957e8; LSUN = 3.828e26; AU = 1.495978707e11
MEARTH = 5.9722e24; REARTH = 6.371e6; DAY = 86400.0; YR = 365.25 * DAY
S0_FLUX = 1361.0                     # W/m^2 a 1 AU
F_XUV_EARTH = 4.64e-3                # W/m^2, XUV solar actual a 1 AU (~4.64 erg/cm2/s, Ribas+2005)  # SUPUESTO (valor de literatura)
SUB = ["S_fis", "S_clima", "S_atm", "S_estrella", "S_geo", "S_mag", "S_obs"]
DIRICHLET_ALPHA = np.array([2, 2, 2, 2, 1, 1, 2], float)   # del Toolkit, seccion 7
TIER_AVAIL = {"T0": 1.0, "T1": 0.75, "T2": 0.5, "T3": 0.25}  # SUPUESTO
# Vector de observabilidad por subindice (Anexo A del documento UHI v2.0 -> tier asignado)  # SUPUESTO en el mapeo exacto
TIER_VECTOR = "fis:T0|clima:T1|atm:T1-T2|estrella:T0|geo:T3*|mag:T3*|obs:T0"


# ---------------- utilidades ----------------
def logistic_log(x, x0, k=4.0):
    """logistica en log10: 0.5 en x = x0; k=4 -> pasa de ~0.02 a ~0.98 en +-1 dex.  # SUPUESTO"""
    x = np.maximum(x, 1e-30)
    return 1.0 / (1.0 + np.exp(-k * (np.log10(x) - np.log10(x0))))


def trapezoid(x, a, b, c, d):
    """0 fuera de [a,d], 1 en [b,c], lineal en los flancos."""
    y = np.zeros_like(x, dtype=float)
    up = (x > a) & (x < b); y[up] = (x[up] - a) / (b - a)
    y[(x >= b) & (x <= c)] = 1.0
    dn = (x > c) & (x < d); y[dn] = (d - x[dn]) / (d - c)
    return y


def split_normal(rng, val, e1, e2, n, frac_default, lo=None, hi=None):
    """Muestra una normal partida con err1 (+) y err2 (-) del catalogo; si faltan, usa frac_default*val."""
    val = float(val)
    s_up = abs(e1) if np.isfinite(e1) and e1 != 0 else abs(val) * frac_default
    s_dn = abs(e2) if np.isfinite(e2) and e2 != 0 else abs(val) * frac_default
    z = rng.standard_normal(n)
    x = np.where(z >= 0, val + z * s_up, val + z * s_dn)
    if lo is not None or hi is not None:
        x = np.clip(x, lo if lo is not None else -np.inf, hi if hi is not None else np.inf)
    return x


# Chen & Kipping (2017) 'Forecaster' M-R (forma deterministica + dispersion intrinseca)
_CK_T1, _CK_T2 = 2.04, 131.6                    # M_earth (2.04 Me; 0.414 MJ)
_CK_R1 = 1.008 * _CK_T1 ** 0.279               # radio en la transicion 1
_CK_R2 = _CK_R1 * (_CK_T2 / _CK_T1) ** 0.589   # radio en la transicion 2


def ck_radius(m, rng):
    r = np.where(m < _CK_T1, 1.008 * m ** 0.279,
         np.where(m < _CK_T2, _CK_R1 * (m / _CK_T1) ** 0.589, _CK_R2 * (m / _CK_T2) ** -0.044))
    sig = np.where(m < _CK_T1, 0.0403, np.where(m < _CK_T2, 0.146, 0.0737))
    return r * (1 + sig * rng.standard_normal(m.shape)).clip(0.5, 1.5)


def ck_mass(r, rng):
    m = np.where(r < _CK_R1, (r / 1.008) ** (1 / 0.279),
         np.where(r < _CK_R2, _CK_T1 * (r / _CK_R1) ** (1 / 0.589), np.nan))
    sig = np.where(r < _CK_R1, 0.0403 / 0.279, 0.146 / 0.589)   # dispersion en R propagada a M
    m = m * np.exp(sig * rng.standard_normal(r.shape))
    giant = ~np.isfinite(m)                                       # rama joviana degenerada  # SUPUESTO
    m[giant] = np.exp(rng.uniform(np.log(_CK_T2), np.log(4000.0), giant.sum()))
    return m


def kopparapu_seff(teff):
    """Kopparapu et al. (2014), 1 M_earth: recent Venus, runaway GH, max GH, early Mars."""
    T = np.clip(teff, 2600, 7200) - 5780.0
    co = {"rv": (1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15),
          "rg": (1.107, 1.332e-4, 1.580e-8, -8.308e-12, -1.931e-15),
          "mg": (0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16),
          "em": (0.320, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16)}
    return {k: s + a * T + b * T ** 2 + c * T ** 3 + d * T ** 4 for k, (s, a, b, c, d) in co.items()}


def fhz_ebm(S, albedo, A=203.3, D=0.649):
    """fHZ con un EBM 1-D de dos modos (North 1975; North, Cahalan & Coakley 1981), media anual,
    oblicuidad terrestre (s2=-0.482), OLR lineal A+B*T (A terrestre=203.3 W/m2, B=2.09 W/m2/C), D=0.649 W/m2/C,
    albedo constante (sin retroalimentacion hielo-albedo). Devuelve la fraccion de area con 0<T<100 C."""
    B, s2 = 2.09, -0.482
    Q = S * S0_FLUX / 4.0
    absorbed = Q * (1 - albedo)
    T0 = (absorbed - A) / B
    T2 = s2 * absorbed / (B + 6 * D)
    # T(x) = T0 + T2*(3x^2-1)/2  es lineal en u = x^2 (x = sen(latitud), area uniforme en x)
    c0 = T0 - T2 / 2.0; c1 = 1.5 * T2
    with np.errstate(divide="ignore", invalid="ignore"):
        u1 = (0.0 - c0) / c1; u2 = (100.0 - c0) / c1
    lo = np.clip(np.minimum(u1, u2), 0, 1); hi = np.clip(np.maximum(u1, u2), 0, 1)
    f = np.sqrt(hi) - np.sqrt(lo)
    flat = np.abs(c1) < 1e-9
    f = np.where(flat, ((c0 > 0) & (c0 < 100)).astype(float), f)
    return np.clip(f, 0, 1)


# ---------------- v2: albedo de Bond dependiente de la Teff de la estrella ----------------
# Shields et al. 2013, Astrobiology 13:715-739, Tabla 2, columna "Earth-like gases, clouds, Rayleigh
# scattering on" (64% de nubes): albedo planetario de banda ancha para superficie oceano y tierra,
# calculado con los SED de HD128167 (F2V), el Sol (G2V), HD22049 (K2V) y AD Leo (M3V).
SH13_OCEAN = {"M": 0.23372, "K": 0.30235, "G": 0.31948, "F": 0.32865}
SH13_LAND = {"M": 0.33165, "K": 0.40133, "G": 0.41484, "F": 0.41428}
# Teff asignadas a esas estrellas (NO vienen de Shields+2013)  # SUPUESTO
#   AD Leo 3477 K, HD22049 (eps Eri) 5020 K (NASA Exoplanet Archive), Sol 5772 K, HD128167 (sigma Boo) 6700 K
ALB_TEFF = np.array([3477.0, 5020.0, 5772.0, 6700.0])
ALB_OCEAN = np.array([SH13_OCEAN[k] for k in "MKGF"])
ALB_LAND = np.array([SH13_LAND[k] for k in "MKGF"])
F_LAND_MAX = 0.5   # fraccion de tierra ~ U(0, 0.5): de mundo oceano a ~1.7x la Tierra (0.29)  # SUPUESTO


def albedo_teff(teff, u, f_land_max=None, shift=0.0, late_m_scale=1.0):
    """Albedo de Bond = A_oceano(Teff) + f_tierra*(A_tierra(Teff) - A_oceano(Teff)), f_tierra = F_LAND_MAX*u.
    Interpolacion lineal en Teff entre las 4 estrellas ancla; fuera de [3477, 6700] K se usa el valor
    del extremo (constante).  # SUPUESTO (interpolacion, extrapolacion constante, fraccion de tierra)
    Superficie sin hielo: el EBM de este codigo no tiene retroalimentacion hielo-albedo."""
    fl = F_LAND_MAX if f_land_max is None else f_land_max
    a_oc = np.interp(teff, ALB_TEFF, ALB_OCEAN)
    a_la = np.interp(teff, ALB_TEFF, ALB_LAND)
    A = a_oc + fl * u * (a_la - a_oc)
    # v3 (solo para sensibilidad): late_m_scale <1 reduce el albedo por debajo de 3477 K; shift = desplazamiento por nubes
    A = np.where(teff < ALB_TEFF[0], A * late_m_scale, A) + shift
    return np.clip(A, 0.0, 0.9)


def planck_ratio(T1, T2, lam=7.5e-6):
    h, c, k = 6.62607015e-34, 2.99792458e8, KB
    return np.expm1(h * c / (lam * k * T2)) / np.expm1(h * c / (lam * k * T1))


def tsm_scale(rp):
    return np.where(rp < 1.5, 0.190, np.where(rp < 2.75, 1.26, np.where(rp < 4.0, 1.28, 1.15)))


def tsm_threshold(rp):
    return np.where(rp < 1.5, 10.0, 90.0)   # umbrales recomendados de Kempton et al. 2018


ESM_THR = 7.5
ELT_IWA_AS = 2 * 1.0e-6 / 39.0 * 206265.0   # 2 lambda/D, lambda=1 um, D=39 m -> ~0.0106"  # (v1-v3; ya no se usa en v4)
CONTRAST_REF = 1e-8                          # contraste de referencia ELT en luz reflejada  # (v1-v3; ya no se usa en v4)



# ---------------- v3: S_atm y S_obs con modelo probabilistico de deteccion ----------------
from scipy.special import ndtr as PHI   # Phi = CDF normal estandar

V3 = dict(
    snr_thr=5.0,      # umbral de deteccion en SNR (convencion 5 sigma)                      # SUPUESTO
    sig_model=0.5,    # incertidumbre de modelo de la SNR predicha (ln)                       # SUPUESTO
    cloud_min=0.5,    # factor de supresion por nubes/neblina ~ U(cloud_min, 1) en transmision  # SUPUESTO
    t_prog=10.0,      # horas de programa de referencia para S_atm (Kempton: TSM ~ S/N en 10 h NIRISS)
    n_ecl=2.0,        # eclipses secundarios para S_atm (Kempton: "at least two eclipses" para S/N>10 en GJ 1132b)
    floor_ppm=20.0,   # piso de ruido sistematico (Greene et al. 2016: 20 ppm NIRISS SOSS)
    n_H=2.0,          # escalas de altura que abarcan los rasgos espectrales                   # SUPUESTO
    n_bins=10.0,      # canales espectrales independientes en que se promedia el piso de ruido  # SUPUESTO
    ntr_cap=10.0,     # maximo de transitos para S_obs                                         # SUPUESTO
    shore_logc=0.0,   # log10 del factor sobre la linea de costa calibrada en Marte            # SUPUESTO
    shore_sig=1.0,    # ancho (ln) de la transicion de la linea de costa                       # SUPUESTO
    # --- v4: no transitantes, luz reflejada con ELT/PCS (XAO + coronografo + HDS) ---
    t_hci=10.0,       # horas de observacion (Kasper+2021: ejemplo de 10 h; Snellen+2015: 1 noche = 10 h)
    n_ph_j26=1.0e5,   # fotones en 10 h de un objeto J=26 (Kasper+2021: ELT, banda J 50 nm, 10% throughput)
    h_15=5.0e-5,      # contraste del halo XAO residual a la menor separacion (Kasper+2021: ~5e-5)
    h_ratio=10.0,     # mejora del halo de 15 a 100 mas (SUPUESTO: igual al factor 10 del requisito PCS 1e-8 -> 1e-9)
    iwa_lo=2.0, iwa_hi=4.0,   # IWA ~ U(iwa_lo, iwa_hi) lambda/D por draw, lambda=0.75 um (SUPUESTO; 15 mas ~ 4 lambda/D)
    ag_fac=2.0 / 3.0, # A_g = (2/3) A_Bond (esfera lambertiana)
)
LAM_HCI = 0.75e-6                                   # centro de 0.6-0.9 um (Snellen+2015); banda A de O2 ~0.76 um
LD_MAS = LAM_HCI / 39.0 * 206264.806 * 1e3          # lambda/D del ELT (D = 39 m) en mas (~3.97)


def halo_contrast(sep_mas, P=V3):
    """Contraste del halo estelar residual (XAO+coronografo) en la posicion del planeta.
    Kasper+2021: ~5e-5 a las menores separaciones; el requisito de PCS mejora x10 de 15 a 100 mas (1e-8 -> 1e-9).
    SUPUESTO: ley de potencias entre 15 y 100 mas, constante dentro de 15 mas y mas alla de 100 mas."""
    s = np.clip(sep_mas, 15.0, 100.0)
    slope = np.log10(P["h_ratio"]) / np.log10(100.0 / 15.0)
    return P["h_15"] * (s / 15.0) ** (-slope)


def snr_reflected(rp, a_au, dist_pc, J, albedo_bond, u_iwa, P=V3):
    """SNR de deteccion en luz reflejada con XAO+coronografo+HDS limitado por ruido fotonico del halo (Kasper+2021):
    SNR = C * sqrt(N_star / h), con N_star = n_ph_j26 * 10^(-0.4 (J-26)) * t/10h.
    Contraste C = A_g * Phi_L(90 deg) * (Rp/a)^2, Phi_L(90) = 1/pi (Lambert, maxima elongacion), A_g = ag_fac*A_B.
    Separacion en maxima elongacion = a/d; cero si queda dentro del IWA."""
    sep_mas = a_au / dist_pc * 1e3
    contrast = P["ag_fac"] * albedo_bond / np.pi * (rp * REARTH / (a_au * AU)) ** 2
    nstar = P["n_ph_j26"] * 10 ** (-0.4 * (J - 26.0)) * P["t_hci"] / 10.0
    snr = contrast * np.sqrt(nstar / halo_contrast(sep_mas, P))
    iwa = (P["iwa_lo"] + (P["iwa_hi"] - P["iwa_lo"]) * u_iwa) * LD_MAS
    return np.where(sep_mas >= iwa, snr, 0.0), contrast, sep_mas


MARS_I, MARS_VESC, EARTH_VESC = 1 / 1.524 ** 2, 5.03, 11.186   # insolacion (S_tierra) y v_esc (km/s)


def p_retain(S, mp, rp, P=V3):
    """Linea de costa cosmica (Zahnle & Catling 2017): frontera I proporcional a v_esc^4.
    Normalizacion (SUPUESTO): la linea pasa por Marte (I=0.431, v_esc=5.03 km/s), por 10**shore_logc.
    P(retener) = Phi( ln(I_shore / I) / shore_sig )."""
    vesc = EARTH_VESC * np.sqrt(np.maximum(mp, 1e-6) / rp)
    i_shore = 10 ** P["shore_logc"] * MARS_I * (vesc / MARS_VESC) ** 4
    return PHI(np.log(i_shore / np.maximum(S, 1e-9)) / P["shore_sig"])


def atm_obs_v3(transit, rp, mp, rs, teff, a_au, per_d, teq, g, J, K, S, avail, dist_pc, albedo, u_iwa, z1, z2, z3, u_cl, u_mu, P=V3):
    """Devuelve (S_atm, S_obs, tsm, esm). Todos arrays del mismo tamano (o broadcastables)."""
    teq0 = teff * np.sqrt(rs * RSUN / (a_au * AU)) * 0.25 ** 0.25
    tsm = tsm_scale(rp) * rp ** 3 * teq0 / (np.maximum(mp, 0.05) * rs ** 2) * 10 ** (-J / 5)
    delta = (rp * REARTH / (rs * RSUN)) ** 2
    esm = 4.29e6 * planck_ratio(1.10 * teq0, teff) * delta * 10 ** (-K / 5)
    m1, m2, m3 = np.exp(P["sig_model"] * z1), np.exp(P["sig_model"] * z2), np.exp(P["sig_model"] * z3)
    cl = P["cloud_min"] + (1 - P["cloud_min"]) * u_cl
    # --- S_atm (transitantes): P(retener) x P(detectar atmosfera con TSM o ESM)
    tsm_hab = np.where((rp >= 1.5) & (rp < 2.0), tsm * 2.3 / 18.0, tsm)   # re-escala de Kempton 3.1 (1.5-2.0 R_E)
    snr_t = tsm_hab * np.sqrt(P["t_prog"] / 10.0) * m1 * cl
    snr_e = esm * np.sqrt(P["n_ecl"]) * m2
    p_det_tr = 1 - (1 - PHI(snr_t - P["snr_thr"])) * (1 - PHI(snr_e - P["snr_thr"]))
    # --- S_atm (no transitantes): luz reflejada ELT
    # v4: modelo fisico ELT/PCS (ver snr_reflected); dist no finita -> 0
    d_ok = np.isfinite(dist_pc) & (dist_pc > 0)
    snr_r0, _, _ = snr_reflected(rp, a_au, np.where(d_ok, dist_pc, 1.0), J, albedo, u_iwa, P)
    snr_r = np.where(d_ok, snr_r0, 0.0) * m1
    p_det_nt = PHI(snr_r - P["snr_thr"])
    pret = p_retain(S, mp, rp, P)
    S_atm = pret * np.where(transit, p_det_tr, p_det_nt)
    # --- S_obs (transitantes): atmosfera secundaria (N2/CO2) o H2, N_tr transitos, piso de ruido
    mu = np.where(rp < 1.6, np.where(u_mu < 0.5, 28.0, 44.0), 2.3)
    mu_k = np.where(rp < 1.5, 18.0, 2.3)               # mu supuesto por Kempton/Louie en cada bin
    snr10 = tsm * mu_k / mu                              # S/N en 10 h (5 h en transito) para ese mu
    H = KB * teq / (mu * MH * g)
    A_sig = P["n_H"] * 2 * H * rp * REARTH / (rs * RSUN) ** 2
    t14_h = per_d * 24 / np.pi * np.arcsin(np.clip(rs * RSUN / (a_au * AU), 0, 1))
    ntr = np.minimum(365.25 / per_d, P["ntr_cap"])
    # ruido fotonico total implicito en la S/N de Kempton, escalado al tiempo en transito disponible
    sig_ph = A_sig / np.maximum(snr10, 1e-12) * np.sqrt(5.0 / np.maximum(ntr * t14_h, 1e-6))
    # piso sistematico por canal (Greene+2016, no baja con mas transitos); se promedia sobre n_bins canales
    sig_fl = P["floor_ppm"] * 1e-6 / np.sqrt(P["n_bins"])
    sig = np.sqrt(sig_ph ** 2 + sig_fl ** 2)
    snr_o = A_sig / sig * m3 * cl
    snr_o = np.where(transit, snr_o, snr_r)
    S_obs = PHI(snr_o - P["snr_thr"]) * avail
    return S_atm, S_obs, tsm, esm, pret


# ======================= v5 (auditoria Neon/R-Aetherios) =======================
EPS = 1e-12            # piso antes de logaritmos; S_i <= EPS se trata como veto exacto (UHI = 0)
V3.update(
    # mu(Rp) log-sigmoide (centro del valle de radios: Fulton et al. 2017, deficit en 1.5-2.0 R_E)
    R0=1.65, R0_sd=0.08,   # centro y dispersion del centro (valores del spec de Neon; 1.65 cae dentro de 1.5-2.0)  # SUPUESTO
    dR=0.12,               # ancho de la transicion [R_E] (spec de Neon)                                          # SUPUESTO
    mu_prim=2.3,           # envoltura primordial H2/He
    rho_ovr=5.5,           # override: rho > 5.5 g/cc con M/sigma_M > 3 -> w_H2 = 0 (spec de Neon)               # SUPUESTO
    # clima de enanas M sincronicas (Teff < 3800 K)
    teff_md=3800.0,
    yang_eff_lo=0.5,       # eficacia del aumento de albedo por nubes subestelares ~ U(yang_eff_lo, 1)            # SUPUESTO
                           # (Kopparapu+2016: rotacion mas rapida de M tardias reduce el efecto de nubes)
    kin_lo=1.5, kin_hi=2.0,  # factor sobre el borde interno 1-D (Kopparapu+2014 rg, rv) ~ U(kin_lo, kin_hi):
                           # Yang+2013 "twice the stellar flux"; Kopparapu+2016 "up to ~25% lower" -> 2 x 0.75 = 1.5
    s_eye_max=0.8,         # el piso de borde externo aplica para S_eff < 0.8 (spec de Neon)
    lpco2_lo=-3.0, lpco2_hi=2.0,   # prior de presion de CO2 disponible: log10(P/bar) ~ U(lo, hi)                # SUPUESTO
    f_eye_c=0.30, f_eye_w=0.15,    # fraccion habitable del piso "eyeball"/CO2 denso ~ U(c-w, c+w)              # SUPUESTO
    # v5.1: limite de radiacion termica saliente de una atmosfera humeda (Goldblatt, Robinson, Zahnle & Crisp 2013,
    # Nature Geosci 6, 661-667, doi:10.1038/ngeo1892; abstract: "thermal radiation limit of 282 W m-2").
    # Si el flujo absorbido S*S0/4*(1-A) > F_rg -> invernadero desbocado -> fHZ = 0 (veto).
    # Variante HIBRIDA (decision de Effy): solo planetas NO sincronicos o Teff >= teff_md; las M sincronicas
    # conservan el borde interno 3D (Yang+2013 x k_in), porque el limite 1-D con redistribucion /4 no les aplica.
    F_rg=282.0, rg_veto=True,
    # v5.2: rampa del peso de Yang+2013 (spec de Neon v5.2) y frontera Teff del veto (= teff_md por defecto)
    yang_ramp_lo=800.0, yang_ramp_hi=1000.0, teff_rg=3800.0,
)
# Yang, Cowan & Abbot 2013 (ApJL 771, L45), Tabla 1, grupo 3 (CAM3, rotacion 1:1, estrella 3400 K):
# albedo planetario vs S0 [W/m2]. v5.1: solo S0 >= 1000 (escalon); v5.2: peso w_eff en rampa 800-1000 W/m2; por encima de 2200 se
# mantiene 0.53 (SUPUESTO); por debajo de 1000 W/m2 se mantiene el albedo de Shields (no se extrapola).
YANG_S0 = np.array([1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0])
YANG_ALB = np.array([0.42, 0.45, 0.47, 0.49, 0.50, 0.51, 0.53])
# Turbet et al. 2018 (A&A 612, A86), planetas sincronicos de TRAPPIST-1 (S en S_E, de su Tabla 1):
#  - S >= 0.662 (TRAPPIST-1e; Proxima b ~890 W/m2): agua liquida subestelar "whatever its atmosphere" (si rico en agua)
#  - TRAPPIST-1f (0.382): CO2 > ~1 bar (menos -> colapso); TRAPPIST-1g (0.258): CO2 > ~5 bar
#  - borde externo con 70 bar de CO2: S_eff = 0.225; por debajo, nunca.
TURBET_S_EYE = 0.662
TURBET_S = np.array([0.225, 0.258, 0.382])
TURBET_LOGP = np.log10(np.array([70.0, 5.0, 1.0]))


def p_co2_req(S):
    """Presion de CO2 [bar] requerida para agua liquida superficial en un planeta sincronico de enana M.
    Anclas de Turbet+2018; interpolacion (lineal en P entre 0.382 y 0.662, log-lineal entre 0.225 y 0.382)  # SUPUESTO"""
    S = np.asarray(S, float)
    lin = 1.0 * (TURBET_S_EYE - S) / (TURBET_S_EYE - TURBET_S[-1])
    logp = 10 ** np.interp(S, TURBET_S, TURBET_LOGP)
    return np.where(S >= TURBET_S_EYE, 0.0, np.where(S >= TURBET_S[-1], lin, np.where(S >= TURBET_S[0], logp, np.inf)))


def climate_v5(S, teff, locked, alb_sh, A_olr, u_eff, u_kin, u_pco2, u_eye, P=V3, D=0.649, phz=(1.0, 0.5, 0.1)):
    """S_clima v5. Para planetas sincronicos de enanas M (locked & Teff<teff_md):
    (1) borde interno: albedo = A_Shields + eff*(A_Yang(S0) - A_Shields) dentro del rango de Yang+2013;
        bordes rg y rv de Kopparapu+2014 multiplicados por k_in ~ U(kin_lo, kin_hi) en el prior p_HZ;
    (2) borde externo: si S < s_eye_max y P_CO2 (prior) >= P_req(S) (Turbet+2018, incluye colapso),
        fHZ = max(fHZ_EBM, f_eye). Resto de planetas: igual que v4.
    Devuelve (albedo, fhz, p_hz, eye_flag)."""
    md = locked & (teff < P["teff_md"])
    s0 = S * S0_FLUX
    a_y = np.interp(s0, YANG_S0, YANG_ALB)
    eff = P["yang_eff_lo"] + (1.0 - P["yang_eff_lo"]) * u_eff
    # v5.2: rampa lineal continua w_eff(S0) = clip((S0 - 800)/200, 0, 1) en lugar del escalon en S0 = 1000 W/m2
    # (v5.1: md & (s0 >= YANG_S0[0])). Por debajo de 1000 W/m2 a_y queda en el primer valor tabulado (0.42, np.interp
    # no extrapola) y su peso cae linealmente a 0 en 800 W/m2: no se extrapola el regimen de Yang+2013.
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
    if P.get("rg_veto", True):   # v5.1: veto hibrido de invernadero desbocado (Goldblatt+2013)
        md_rg = locked & (teff < P.get("teff_rg", P["teff_md"]))   # v5.2: frontera del veto separable (barrido 1-D)
        fhz = np.where((~md_rg) & (S * S0_FLUX / 4.0 * (1 - alb) > P["F_rg"]), 0.0, fhz)
    return alb, fhz, p_hz, eye


def mu_v5(rp, rho_gcc, dens_ok, u_mu, z_R0, P=V3):
    """Peso molecular medio con transicion log-sigmoide en el valle de radios:
    w_H2 = 1/(1+exp(-(Rp-R0)/dR)), R0 ~ N(R0, R0_sd); ln mu = (1-w) ln mu_sec + w ln mu_prim (interpolacion geometrica).
    mu_sec = 28 (N2) o 44 (CO2) con prob. 1/2 (como v3-v4). Override: masa con M/sigma_M>3 y rho>rho_ovr -> w=0."""
    R0 = P["R0"] + P["R0_sd"] * z_R0
    w = 1.0 / (1.0 + np.exp(-(rp - R0) / P["dR"]))
    w = np.where(dens_ok & (rho_gcc > P["rho_ovr"]), 0.0, w)
    mu_sec = np.where(u_mu < 0.5, 28.0, 44.0)
    return np.exp((1 - w) * np.log(mu_sec) + w * np.log(P["mu_prim"])), w


def atm_obs_v5(transit, rp, mp, rs, teff, a_au, per_d, teq, g, J, K, S, avail, dist_pc, albedo, u_iwa,
               z1, z2, z3, u_cl, u_mu, z_R0, dens_ok, P=V3):
    """v5: desacople ortogonal.
    S_atm = P_ret(v_esc, I) * f_flare, f_flare = 1 (no hay trazador de actividad por planeta en el catalogo).
            P_ret = linea de costa cosmica (Zahnle & Catling 2017) con la misma (Mp, Rp) del draw. Sin terminos instrumentales.
    S_obs = P_obs(SNR_eff) * Availability(tier):
            transitantes: 1-(1-Phi(SNR_trans-u))(1-Phi(SNR_emis-u)) con SNR_trans (mu log-sigmoide, piso de ruido) y
            SNR_emis = ESM*sqrt(n_ecl)  [combinacion O, como v3 en S_atm]  # SUPUESTO
            no transitantes: Phi(SNR_reflejada-u) (ELT/PCS, igual que v4).
    Devuelve (S_atm, S_obs, tsm, esm, pret, w_H2)."""
    teq0 = teff * np.sqrt(rs * RSUN / (a_au * AU)) * 0.25 ** 0.25
    tsm = tsm_scale(rp) * rp ** 3 * teq0 / (np.maximum(mp, 0.05) * rs ** 2) * 10 ** (-J / 5)
    delta = (rp * REARTH / (rs * RSUN)) ** 2
    esm = 4.29e6 * planck_ratio(1.10 * teq0, teff) * delta * 10 ** (-K / 5)
    m1, m2, m3 = np.exp(P["sig_model"] * z1), np.exp(P["sig_model"] * z2), np.exp(P["sig_model"] * z3)
    cl = P["cloud_min"] + (1 - P["cloud_min"]) * u_cl
    pret = p_retain(S, mp, rp, P)
    S_atm = pret * 1.0                                   # f_flare = 1
    # transmision
    rho = 5.513 * mp / rp ** 3
    mu, w = mu_v5(rp, rho, dens_ok, u_mu, z_R0, P)
    mu_k = np.where(rp < 1.5, 18.0, 2.3)                 # mu supuesto por Kempton/Louie en cada bin (normalizacion del TSM)
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
    # reflejada (no transitantes)
    d_ok = np.isfinite(dist_pc) & (dist_pc > 0)
    snr_r0, _, _ = snr_reflected(rp, a_au, np.where(d_ok, dist_pc, 1.0), J, albedo, u_iwa, P)
    snr_r = np.where(d_ok, snr_r0, 0.0) * m1
    p_nt = PHI(snr_r - P["snr_thr"])
    S_obs = np.where(transit, p_tr, p_nt) * avail
    return S_atm, S_obs, tsm, esm, pret, w


def aggregate_v5(Smat, W):
    """UHI = prod S_i^w_i con piso EPS antes del log; cualquier S_i <= EPS es veto exacto (UHI = 0)."""
    alive = np.all(Smat > EPS, axis=0)
    uhi = np.exp((W * np.log(np.maximum(np.minimum(Smat, 1.0), EPS))).sum(0))
    return np.where(alive, uhi, 0.0)

# ---------------- preparacion del catalogo ----------------
def prepare(df):
    d = df.copy()
    d["J_use"] = d["sy_jmag"].fillna(d["sy_kmag"])       # SUPUESTO: J ~ K si falta J
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
    # transitante con Msini: sin i ~ 1, se trata como masa medida  # SUPUESTO
    sig_m = 0.5 * (d.pl_bmasseerr1.abs() + d.pl_bmasseerr2.abs())
    d["dens_ok"] = (d.mass_mode == "measured") & (d.pl_bmasse / sig_m > 3.0)   # v5: M/sigma_M > 3
    return d


def nominal_tier(row):
    """Tier planetario T0-T3 (cuando es caracterizable de forma util), sobre valores nominales.  # SUPUESTO
    T0: transita, masa medida (sigma/M<=0.5) y TSM>=umbral o ESM>=7.5 (JWST hoy; 'espectros de objetos grandes').
    T1: transita y max(TSM/umbral, ESM/7.5) >= 0.5 (objetos brillantes; masas por RV de alta precision/ELT).
    T2: transita con Rp>=1.6 (no-terrestre, ARIEL) y metrica>=0.1, o no transita con SNR nominal ELT/PCS
        (v4, snr_reflected con A_B=0.3 e IWA=2 lambda/D) >= 0.1 x snr_thr (caracterizacion selectiva ELT).
    T3: el resto (requiere observatorios post-2035 / HWO)."""
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
        # v4: SNR nominal ELT/PCS (A_B=0.3, IWA=2 lambda/D) >= 0.1 x umbral, analogo a "metrica>=0.1" de transitantes  # SUPUESTO
        snr, _, _ = snr_reflected(rp, a_au, row.sy_dist, row.J_use, 0.3, 0.0)
        if snr >= 0.1 * V3["snr_thr"]:
            return "T2"
        return "T3"


# ---------------- Monte Carlo por planeta ----------------
def simulate_planet(r, n, rng, alpha=None, keep_draws=False):
    out = {}
    # --- estrella ---
    teff = split_normal(rng, r.st_teff, r.st_tefferr1, r.st_tefferr2, n, 0.02, 2300, 50000)
    rs = split_normal(rng, r.st_rad, r.st_raderr1, r.st_raderr2, n, 0.07, 0.05, None)
    ms = split_normal(rng, r.st_mass, r.st_masserr1, r.st_masserr2, n, 0.07, 0.05, None)
    L = rs ** 2 * (teff / 5772.0) ** 4                       # L* desde R* y Teff (tabla 0.1 del Toolkit)
    if np.isfinite(r.st_age):
        age = split_normal(rng, r.st_age, r.st_ageerr1, r.st_ageerr2, n, 0.3, 0.01, 13.8)
    else:
        age = rng.uniform(1.0, 10.0, n)                       # SUPUESTO: edad desconocida ~ U(1,10) Gyr
    # --- orbita ---
    if np.isfinite(r.pl_orbper):
        P = split_normal(rng, r.pl_orbper, r.pl_orbpererr1, r.pl_orbpererr2, n, 1e-4, 1e-3, None)
        a = ((G * ms * MSUN) / (4 * np.pi ** 2)) ** (1 / 3) * (P * DAY) ** (2 / 3) / AU   # Kepler (Toolkit 1)
    else:
        a = split_normal(rng, r.pl_orbsmax, r.pl_orbsmaxerr1, r.pl_orbsmaxerr2, n, 0.05, 1e-4, None)
        P = np.sqrt(a ** 3 / ms) * 365.25
    if np.isfinite(r.pl_orbeccen):
        e = split_normal(rng, r.pl_orbeccen, r.pl_orbeccenerr1, r.pl_orbeccenerr2, n, 0.0, 0.0, 0.95)
        if not (np.isfinite(r.pl_orbeccenerr1) or np.isfinite(r.pl_orbeccenerr2)):
            e = np.clip(r.pl_orbeccen + 0.05 * rng.standard_normal(n), 0, 0.95)   # SUPUESTO
    else:
        e = np.zeros(n)                                       # SUPUESTO: e desconocida = 0
    # --- planeta: masa y radio ---
    if r.radius_measured:
        rp = split_normal(rng, r.pl_rade, r.pl_radeerr1, r.pl_radeerr2, n, 0.10, 0.1, None)
    else:
        rp = None
    if r.mass_mode == "measured":
        mp = split_normal(rng, r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, n, 0.20, 0.01, None)
    elif r.mass_mode == "msini":
        msini = split_normal(rng, r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, n, 0.20, 0.01, None)
        cosi = rng.uniform(0, 1, n)                           # inclinacion isotropica (no transita)
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
    # --- flujo y temperaturas (Toolkit 1) ---
    S = L / a ** 2
    alb_sh = albedo_teff(teff, rng.uniform(0.0, 1.0, n))      # v2: albedo dependiente de Teff (Shields+2013, Tabla 2); f=1
    # v5: draws adicionales con generador propio (semilla = crc32 del nombre + 1) para no alterar la secuencia principal
    rng5 = np.random.default_rng(zlib.crc32(str(r.pl_name).encode()) + 1)
    u5 = rng5.uniform(size=(4, n)); z_R0 = rng5.standard_normal(n)
    locked = (r.st_mass < 0.6) & (np.median(a) < 0.2)       # SUPUESTO (regla de rotacion sincronica de v1-v4)
    A_olr = rng.uniform(150.0, 203.3, n)                      # SUPUESTO: invernadero incierto (Tierra -> CO2 reforzado)
    albedo, fhz, p_hz, eye = climate_v5(S, teff, locked, alb_sh, A_olr, u5[0], u5[1], u5[2], u5[3])
    f_red = 1.0
    teq = ((1 - albedo) * L * LSUN / (16 * np.pi * SIGMA * (a * AU) ** 2 * f_red)) ** 0.25
    delta = (rp * REARTH / (rs * RSUN)) ** 2
    g = G * mp * MEARTH / (rp * REARTH) ** 2

    # S_fis: tamano/masa rocosos (trapecios)  # SUPUESTO umbrales
    f_r = trapezoid(rp, 0.3, 0.5, 1.6, 2.5)
    f_m = trapezoid(mp, 0.03, 0.1, 5.0, 10.0)
    S_fis = np.sqrt(f_r * f_m)

    # S_clima = fHZ(EBM, + v5 clima de enanas M sincronicas) * prior Kopparapu (v5: borde interno desplazado)
    S_clima = fhz * p_hz

    # v3: S_atm y S_obs probabilisticos (ver atm_obs_v3)
    tr = (r.tran_flag == 1)
    z = rng.standard_normal((3, n)); u_cl = rng.uniform(size=n); u_mu = rng.uniform(size=n)
    # v4: IWA por draw con generador propio (semilla = crc32 del nombre) para no alterar la secuencia principal
    u_iwa = np.random.default_rng(zlib.crc32(str(r.pl_name).encode())).uniform(size=n)
    S_atm, S_obs, tsm, esm, pret, w_h2 = atm_obs_v5(tr, rp, mp, rs, teff, a, P, teq, g, r.J_use, r.K_use, S,
                                               TIER_AVAIL[r.tier], r.sy_dist, albedo, u_iwa, z[0], z[1], z[2], u_cl, u_mu,
                                               z_R0, bool(r.dens_ok))

    # S_estrella: XUV (actividad-edad) + riesgo de flares
    tsat = 0.1 if r.st_mass >= 0.6 else 1.0                    # SUPUESTO
    rx = np.where(age < tsat, 10 ** -3.5, 10 ** -3.5 * (age / tsat) ** -1.23)   # Ribas+2005 pendiente -1.23
    rx_sun = 10 ** -3.5 * (4.6 / 0.1) ** -1.23
    fxuv_rel = rx * S / rx_sun                                 # F_XUV en unidades de la Tierra actual
    stress = logistic_log(fxuv_rel, 10.0)                      # SUPUESTO: 50% de estres a 10x Tierra
    flare = 0.5 if r.st_mass < 0.35 else (0.25 if r.st_mass < 0.6 else 0.0)   # SUPUESTO
    S_estrella = (1 - stress) * (1 - flare)

    # S_geo: flujo interno (radiogenico proxy + marea) -> regimen; erosion como bandera para sub-Neptunos
    n_orb = 2 * np.pi / (P * DAY)
    k2Q = 0.3 / 100.0                                          # SUPUESTO k2=0.3, Q=100
    Edot = 10.5 * k2Q * G * (ms * MSUN) ** 2 * (rp * REARTH) ** 5 * n_orb * e ** 2 / (a * AU) ** 6
    q_tide = Edot / (4 * np.pi * (rp * REARTH) ** 2)
    q_rad = 0.087 * (mp / rp ** 2) * 2.0 ** (-(age - 4.5) / 3.0)   # SUPUESTO
    lq = np.log10(np.maximum(q_rad + q_tide, 1e-6))
    S_regime = trapezoid(lq, -2.0, np.log10(0.03), 0.0, 1.0)  # SUPUESTO: optimo 0.03-1 W/m2
    fxuv_abs = fxuv_rel * F_XUV_EARTH
    mdot = 0.1 * np.pi * (rp * REARTH) ** 3 * fxuv_abs / (G * mp * MEARTH)   # eps=0.1, K=1  # SUPUESTO
    t_loss = 0.01 * mp * MEARTH / np.maximum(mdot, 1e-30) / (1e9 * YR)       # envolvente 1% M_p  # SUPUESTO
    risk = 1 - logistic_log(t_loss, 1.0)
    subnep = (rp >= 1.6) & (rp <= 4.0)
    S_geo = S_regime * np.where(subnep, 1 - 0.5 * risk, 1.0)

    # S_mag: prior Beta amplio con actualizaciones indirectas  # SUPUESTO
    aB, bB = (2.0, 3.0) if locked else (2.0, 2.0)
    aB = np.full(n, aB); bB = np.where(mp < 0.3, bB + 1.0, bB)
    S_mag = rng.beta(aB, bB)

    # Agregacion UHI con pesos Dirichlet por muestra (Toolkit 7)
    Smat = np.vstack([S_fis, S_clima, S_atm, S_estrella, S_geo, S_mag, S_obs])
    W = rng.dirichlet(DIRICHLET_ALPHA if alpha is None else alpha, n).T
    uhi = aggregate_v5(Smat, W)                                # v5: piso EPS; S_i <= EPS => UHI = 0 exacto (veto)

    # ESI (Schulze-Makuch et al. 2011)
    rho = mp / rp ** 3                      # densidad relativa a la Tierra
    vesc = np.sqrt(mp / rp)                 # velocidad de escape relativa
    Ts = teq + 33.0                         # SUPUESTO: T_superficie = T_eq + 33 K (efecto invernadero terrestre)
    # v5.2 (errata): Schulze-Makuch+2011: ESI_I = prod_{R,rho} term^(w/n) con n = 2 por grupo, ESI = sqrt(ESI_I*ESI_S)
    # => exponente total w/4 (equivale a la forma de 4 propiedades con n = 4). v1-v5.1 usaban w/4 DENTRO de cada grupo
    # y ademas la raiz (exponente total w/8), lo que comprimia el ESI hacia 1. No afecta al UHI (solo comparacion).
    def esi_term(x, x0, w):
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
    out["UHI_rank"] = out["P_alive"] * out["UHI_cond_med"]     # v5: hurdle, UHI_rank = P_alive x mediana condicional
    if keep_draws and out["P_alive"] > 0:
        out["_draws"] = uhi.astype(np.float32)
    return out


def run(df, n_draws=10000, seed=42, progress=None, alpha=None, joint=True):
    d = prepare(df)
    # nominales para el tier
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
        idx = np.array(sorted(draws)); M = np.vstack([draws[i] for i in idx])          # (planetas vivos, draws)
        npos = (M > 0).sum(0)                                                           # vivos por draw
        order = np.argsort(-M, axis=0, kind="stable"); rk = np.empty_like(order)
        np.put_along_axis(rk, order, np.arange(1, M.shape[0] + 1)[:, None].repeat(M.shape[1], 1), axis=0)
        rk = np.where(M > 0, rk, npos[None, :] + 1)                                     # UHI=0 -> rango = vivos+1
        res.loc[idx, "P_top10"] = (rk <= 10).mean(1)
        res.loc[idx, "rank_med"] = np.median(rk, 1)
        res.loc[idx, "rank_ci5"] = np.percentile(rk, 5, 1); res.loc[idx, "rank_ci95"] = np.percentile(rk, 95, 1)
    res["ESI_rank"] = res.ESI_med.rank(ascending=False, method="min").astype("Int64")
    # v5: ranking por UHI_rank = P_alive x mediana condicional; empates por la media del posterior
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
