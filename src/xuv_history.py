"""UHI v6 - historia XUV y retencion atmosferica (Fase 1, parte A).

Fuentes (ver Obsidian Referencias/Fase_1.md):
  * Pass, Charbonneau & Vanderburg 2025, ApJL 986, L3 (doi:10.3847/2041-8213/adda39): Ec. 6-11, Tabla 1.
      I_XUV = f_flare * 10^-3.05 * f_L * L_bol * t_sat / (4 pi a^2);  t_sat = 5.9 - 15.4 M*/Msun [Gyr] (0.1-0.3 Msun)
      I_XUV,Earth = 2.0e18 erg cm^-2 (integral de Lammer+2009 Ec.2 de 0 a 4.6 Gyr);  ARM = 4 log v_esc - log I_XUV - 3.16
  * Johnstone, Bartel & Gudel 2021, A&A 649, A96 (codigo Mors; pistas Spada+2013) para M* > 0.35 Msun.
  * Jackson, Davis & Wheatley 2012, MNRAS 422, 2024, Tabla 2 (variante/fallback con (B-V)_0 del catalogo).
  * Sanz-Forcada+2011, A&A 532, A6: log L_EUV = 4.80 + 0.860 log L_X (cgs).
  * Owen & Wu 2017, ApJ 847, 29, Ec. 19 y 31: Mdot = eta pi R^3 F_XUV / (G M), eta = 0.1 (v_esc / 15 km/s)^-2.
  * Berta-Thompson+2026, ApJ 1011, 123: costa 3D (variante; ver p_ret_bt3d).
Lo marcado SUPUESTO no tiene fuente directa (documentado en assumptions_v6.md).
"""
import os
import numpy as np
from scipy.special import ndtr as PHI

LSUN_CGS = 3.828e33; AU_CM = 1.495978707e13; GYR_S = 3.15576e16
I_XUV_EARTH = 2.0e18                    # erg cm^-2 (Pass+2025)
LOG_RSAT_FC = -3.05                     # Wright+2018 via Pass+2025 (+0.05/-0.06)
PASS_M = np.array([0.10, 0.12, 0.15, 0.20, 0.25, 0.30])     # Pass+2025 Tabla 1
PASS_FL = np.array([1.19, 1.17, 1.13, 1.13, 1.16, 1.22])
STILL_SATURATED = {"TRAPPIST-1"}       # Pass+2025 Sec. 2.2: TRAPPIST-1 sigue saturada hoy (Prot 3.3 d) -> se integra hasta la edad
# Jackson+2012 Tabla 2: (B-V)0 lo, hi, log(Lx/Lbol)sat, log tau_sat [yr], alpha
JACKSON = np.array([[0.290, 0.450, -4.28, 7.87, 1.22], [0.450, 0.565, -4.24, 8.35, 1.24], [0.565, 0.675, -3.67, 7.84, 1.13],
                    [0.675, 0.790, -3.71, 8.03, 1.28], [0.790, 0.935, -3.36, 7.90, 1.40], [0.935, 1.275, -3.35, 8.28, 1.09],
                    [1.275, 1.410, -3.14, 8.21, 1.18]])
# Innes, Tsai & Pierrehumbert 2023, ApJ 953, 168, Tabla 1: limite de desbocado [S0] para 1 y 10 bar de H2-He
INNES_TEFF = np.array([3457.0, 5772.0])          # estrella M = luminosidad de K2-18 (Teff de K2-18 del catalogo); G = Sol
INNES_S1 = np.array([0.320, 0.389]); INNES_S10 = np.array([0.0860, 0.0676])   # v6.1: M 10 bar = 0.0860 (Tabla 1 publicada; v6.0 usaba 0.0852 del preprint)

_GRID = None


def mors_grid():
    global _GRID
    if _GRID is None:
        p = os.path.join(os.path.dirname(__file__), "..", "data", "mors_grid.npz")
        _GRID = dict(np.load(p)) if os.path.exists(p) else {}
    return _GRID


def fl_pass(mstar):
    return np.interp(mstar, PASS_M, PASS_FL)   # extrapolacion constante fuera de 0.10-0.30 (SUPUESTO)


def tsat_pass(mstar, scale=1.0):
    """Pass+2025 Ec. 7 (ajuste 0.1-0.3 Msun; Pass la extrapolan a 0.35). M<0.1 se fija en 0.1 (SUPUESTO). Piso 0.01 Gyr."""
    return np.maximum((5.9 - 15.4 * np.maximum(mstar, 0.10)) * scale, 0.01)


def fluence_fc(L, a_au, age, mstar, log_rsat, tsat, f_flare, saturated_now=False, post_drop=1e-2):
    """Fluencia en rayos X (banda ROSAT, proxy XUV) [erg cm^-2] para M completamente convectivas (Pass+2025 Ec. 10).
    post_drop: L_X/L_bol tras t_sat = post_drop * R_sat (SUPUESTO conservador; Pass: caida > 2 dex, fluencia despreciable).
    Devuelve (fluencia, L_X hoy [erg/s])."""
    Lc = L * LSUN_CGS; rs = 10 ** log_rsat; fl = fl_pass(mstar)
    t_on = np.where(saturated_now, age, np.minimum(tsat, age))
    E = f_flare * rs * fl * Lc * t_on * GYR_S + np.where(saturated_now, 0.0, post_drop * rs * Lc * np.maximum(age - tsat, 0.0) * GYR_S)
    Lx_now = np.where(saturated_now | (age < tsat), f_flare * rs * Lc, post_drop * rs * Lc)
    return E / (4 * np.pi * (a_au * AU_CM) ** 2), Lx_now


def fluence_mors(L, a_au, age, mstar, pct):
    """Johnstone+2021 (Mors): interpolacion bilineal en (masa, percentil) de la integral acumulada de L_X.
    Se reescala por L_obs / L_bol,pista(edad) para respetar la luminosidad observada (SUPUESTO). Masa fuera de 0.10-1.20 -> extremo (Mors falla a 1.25 Msun: edad 5 Gyr fuera de la pista).
    Devuelve (fluencia [erg cm^-2], L_X hoy, L_EUV hoy) o None si no hay malla."""
    g = mors_grid()
    if not g:
        return None
    M, P, A = g["masses"], g["pcts"], g["ages"]
    m = np.clip(mstar, M[0], M[-1]); p = np.clip(pct, P[0], P[-1]); t = np.clip(age, A[0], A[-1])
    i = np.clip(np.searchsorted(M, m) - 1, 0, len(M) - 2); wm = (m - M[i]) / (M[i + 1] - M[i])
    j = np.clip(np.searchsorted(P, p) - 1, 0, len(P) - 2); wp = (p - P[j]) / (P[j + 1] - P[j])
    k = np.clip(np.searchsorted(A, t) - 1, 0, len(A) - 2); wt = (t - A[k]) / (A[k + 1] - A[k])

    def interp(Q):
        out = 0.0
        for di, fm in ((0, 1 - wm), (1, wm)):
            for dj, fp in ((0, 1 - wp), (1, wp)):
                out = out + fm * fp * (Q[i + di, j + dj, k] * (1 - wt) + Q[i + di, j + dj, k + 1] * wt)
        return out
    lbol_tr = interp(g["Lbol"]); scale = L * LSUN_CGS / np.maximum(lbol_tr, 1e20)
    E = interp(g["CLx"]) * scale
    return E / (4 * np.pi * (a_au * AU_CM) ** 2), interp(g["Lx"]) * scale, interp(g["Leuv"]) * scale


def fluence_jackson(L, a_au, age, bv, u_scatter=None):
    """Jackson+2012 Tabla 2 por (B-V)_0 del catalogo (enrojecimiento ignorado: SUPUESTO). Ley partida:
    L_X/L_bol = R_sat (t < tau), R_sat (t/tau)^-alpha (t > tau). bv fuera de [0.29, 1.41] -> NaN (usar otro modelo)."""
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
    """(B-V)_0 de enana de la tabla de Pecaut & Mamajek 2013 (version 2022.04.16, md5 afa3e0c8b249e8b9511b917c03aceb84),
    interpolada en Teff. Se usa solo en la variante Jackson porque el snapshot no trae B ni V. SUPUESTO: estrella de secuencia principal."""
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
    return 10 ** (4.80 + 0.860 * np.log10(np.maximum(Lx, 1.0)))   # Sanz-Forcada+2011


def arm(vesc_kms, I_rel):
    return 4 * np.log10(vesc_kms) - np.log10(np.maximum(I_rel, 1e-30)) - 3.16   # Pass+2025 Ec. 11


def p_ret_xuv(vesc_kms, I_rel, shore_logc=0.0, shore_sig=1.0):
    """P_ret = Phi((ARM + c) ln10 / sigma_ln). c y sigma_ln con los mismos valores base que v5.2 (0 y 1.0). SUPUESTO."""
    return PHI((arm(vesc_kms, I_rel) + shore_logc) * np.log(10.0) / shore_sig)


def p_ret_bt3d(vesc_kms, S_rel, L_rel, shore_sig=1.0, p=5.9, q=1.17, logL_hz=-2.23):
    """Variante Berta-Thompson+2026 (costa 3D bolometrica): log S_crit = p log v_esc + q log L* + const.
    NO VERIFICADO: la constante de normalizacion no se extrajo del paper; aqui se fija para que la Tierra en la HZ
    (S=1) quede sobre la costa con L* = 10^logL_hz y v_esc = 11.186 km/s, que es como el resumen describe el
    limite (~M4V). Solo se usa como variante de sensibilidad."""
    const = np.log10(1.0) - p * np.log10(11.186) - q * logL_hz
    log_scrit = p * np.log10(vesc_kms) + q * np.log10(L_rel) + const
    return PHI((log_scrit - np.log10(np.maximum(S_rel, 1e-30))) * np.log(10.0) / shore_sig)


def s_crit_innes(teff, p_h2_bar):
    """Limite de desbocado Hycean [S0]: log S lineal en log P entre las anclas de 1 y 10 bar (Innes+2023 Tabla 1),
    interpolado linealmente en Teff entre la estrella M (3457 K) y G (5772 K), constante fuera (SUPUESTO).
    Fuera de 1-10 bar se extrapola con la misma pendiente (Innes: OLR_max ~ 1/masa seca; SUPUESTO)."""
    s1 = np.interp(teff, INNES_TEFF, INNES_S1); s10 = np.interp(teff, INNES_TEFF, INNES_S10)
    return 10 ** (np.log10(s1) + (np.log10(s10) - np.log10(s1)) * np.log10(p_h2_bar))


def envelope_t_loss_gyr(rp_re, mp_me, p_h2_bar, fxuv_cgs, eta_fac=1.0):
    """t_loss = M_env / Mdot [Gyr]. M_env = 4 pi R^2 P / g (envoltura fina hidrostatica); Owen & Wu 2017 Ec. 19/31."""
    G = 6.674e-11; R = rp_re * 6.371e6; M = mp_me * 5.9722e24
    g = G * M / R ** 2; menv = 4 * np.pi * R ** 2 * p_h2_bar * 1e5 / g
    vesc = np.sqrt(2 * G * M / R) / 1e3
    eta = 0.1 * (vesc / 15.0) ** -2 * eta_fac
    mdot = eta * np.pi * R ** 3 * (fxuv_cgs * 1e-3) / (G * M)      # erg/cm2/s -> W/m2 (x1e-3)
    return menv / np.maximum(mdot, 1e-30) / GYR_S
