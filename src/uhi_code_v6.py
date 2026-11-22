"""UHI v6 - Fase 1 (Effy Viridis): transicion climatica continua por rotacion/acoplamiento de marea (D1),
historia XUV acumulada en S_atm (D2, D3) y division Rocky/Hycean (D4, D5).

Construido SOBRE uhi_code.py v5.2 (copia identica, md5 0288ce510341d50e9b3a769307ca2652) sin modificarlo.
Con CFG_COMPAT (todos los modulos apagados) reproduce v5.2 bit a bit (tests/test_v6.py).

Reglas de reproducibilidad:
  * la secuencia del generador principal (semilla 42, planetas en orden) es IDENTICA a v5.2: se consumen los mismos
    draws en el mismo orden (incluido ck_mass/ck_radius aunque luego se sustituyan) -> las diferencias son solo fisica;
  * todos los draws nuevos salen de rng6 = default_rng(crc32(nombre) + 2), con un bloque fijo de 20 uniformes y
    8 normales por planeta (agregar parametros no altera los existentes).
Lo marcado SUPUESTO no tiene fuente directa; ver assumptions_v6.md."""
import zlib
import numpy as np
import pandas as pd
import uhi_code as u
from uhi_code import (G, SIGMA, KB, MH, MSUN, RSUN, LSUN, AU, MEARTH, REARTH, DAY, YR, S0_FLUX, F_XUV_EARTH, SUB,
                      DIRICHLET_ALPHA, TIER_AVAIL, TIER_VECTOR, EPS, V3, PHI, YANG_S0, YANG_ALB, split_normal, trapezoid,
                      logistic_log, ck_mass, ck_radius, albedo_teff, fhz_ebm, kopparapu_seff, p_co2_req, climate_v5,
                      atm_obs_v5, aggregate_v5, prepare, nominal_tier, tsm_scale, planck_ratio, snr_reflected, p_retain,
                      _CK_R1, _CK_T1, _CK_T2, _CK_R2)
import xuv_history as X
import composition as C

VERSION = "v6.1-fase1 (2026-11-20)"

CFG_V6 = dict(climate_mode="rot", xuv=True, shore_model="xuv4", xuv_star="mors", d3=True, branches=True,
              hycean_climate="innes", class_prior=(1.0, 1.0, 1.0),
              msini="iso", mr="otegi", rock_line="h2o50", rs_floor=True)   # v6.1: piso fisico de R* (fix 2)   # pruebas de prior: msini iso|min|howard; mr otegi|ck; rock_line h2o50|h2o100|earthlike
CFG_COMPAT = dict(climate_mode="v52", xuv=False, shore_model="xuv4", xuv_star="mors", d3=False, branches=False,
                  hycean_climate="innes", class_prior=(1.0, 1.0, 1.0))

V6 = dict(
    # D1: probabilidad de rotacion sincronica (reemplaza locked = M<0.6 & a<0.2 de v1-v5.2, que queda como centro)
    lock_m_lo=0.55, lock_m_hi=0.65, lock_m_w=0.05,     # M_lock ~ U, ancho logistico [Msun]          # SUPUESTO
    lock_a_lo=0.15, lock_a_hi=0.25, lock_a_w=0.03,     # a_lock ~ U, ancho [AU]                       # SUPUESTO
    # eficacia 3D vs periodo (Kopparapu+2016 ~10 d; Haqq-Misra+2018: rapido <5 d, Rhines 5-20 d, lento >20 d)
    prot_lo_lo=4.0, prot_lo_hi=6.0, prot_hi_lo=15.0, prot_hi_hi=25.0,
    # variante 'teff': w3D = p_sync * sigmoide((T_c - Teff)/dT)
    tc_lo=3700.0, tc_hi=3900.0, dT_lo=50.0, dT_hi=150.0,
    # D2: XUV (Pass+2025 para completamente convectivas; Mors para M > 0.35)
    log_rsat=-3.05, log_rsat_sd=0.06, tsat_sd_ln=0.3, fflare_lo=1.0, fflare_hi=1.5, post_drop=1e-2,
    blend_m=0.35, blend_w=0.03, pct_lo=5.0, pct_hi=95.0,
    # D4: clases
    wv_lo=0.08, wv_hi=0.20,                            # ancho del valle de radios [R_E]               # SUPUESTO
    mv_lo=3.0, mv_hi=6.0, mv_s=0.3,                    # ruta de masa (sin radio)                      # SUPUESTO
    # Hycean
    lph2_lo=0.0, lph2_hi=2.0,                          # log10 P_H2 [bar] ~ U (Innes: 1-10 bar; extendido)  # SUPUESTO
    fcl_lo=1.0, fcl_hi=2.0,                            # factor de tolerancia sobre S_crit (Innes)    # SUPUESTO
    eta_sd=0.5, tloss_ref=1.0, tloss_sig=1.0,          # S_atm,H = Phi(ln(t_loss/1 Gyr)/1)             # SUPUESTO
    madhu_teq=500.0,                                   # variante Madhusudhan+2021 (borde interno T_eq ~ 500 K)
    hyc_alb=0.12,                                      # albedo de Bond de Innes+2023
    hf_a=1.0, hf_b=2.0, hf_c=10.0, hf_d=15.0,          # S_fis,H = trap(M; 1, 2, 10, 15) M_E           # SUPUESTO
)
PARAMS_DEFAULT = {**V3, **V6}


def r_ms(m):
    """Relacion masa-radio de secuencia principal, Demircan & Kahraman 1991 (Ap&SS 181, 313):
    R = 1.06 M^0.945 (M < 1.66 Msun); R = 1.33 M^0.555 (M > 1.66)."""
    m = np.maximum(m, 0.08)
    return np.where(m < 1.66, 1.06 * m ** 0.945, 1.33 * m ** 0.555)


def apply_rs_floor(rs, ms, rs_cat, ms_cat, f=0.5):
    """v6.1 (fix 2): piso fisico R* >= f * R_MS(M*_draw) (f = 0.5, SUPUESTO conservador; ver assumptions_v6.md).
    Solo si el radio de catalogo es compatible con secuencia principal (rs_cat >= f R_MS(ms_cat)); asi no se tocan
    enanas blancas ni anfitriones compactos raros. Se recorta (no se re-muestrea) para no alterar la secuencia del generador."""
    if not (np.isfinite(rs_cat) and np.isfinite(ms_cat)) or rs_cat < f * float(r_ms(ms_cat)):
        return rs
    return np.maximum(rs, f * r_ms(ms))


def sig(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -500, 500)))


def climate_v6(S, teff, sync, alb_sh, A_olr, eff, kin, u_pco2, u_eye, P, D=0.649, phz=(1.0, 0.5, 0.1)):
    """Igual que climate_v5 pero con el regimen 3D decidido por draw (sync) y eficacia/k_in fisicos (D1).
    Ramas sincronicas: albedo de Yang+2013 (rampa 800-1000 W/m2) x eff, borde interno x k_in, piso eyeball de Turbet+2018,
    SIN el veto 1-D. Ramas no sincronicas: EBM 1-D con veto de Goldblatt+2013 (282 W/m2)."""
    s0 = S * S0_FLUX
    a_y = np.interp(s0, YANG_S0, YANG_ALB)
    w_eff = np.clip((s0 - P["yang_ramp_lo"]) / (P["yang_ramp_hi"] - P["yang_ramp_lo"]), 0.0, 1.0)
    alb = np.where(sync, alb_sh + w_eff * eff * (a_y - alb_sh), alb_sh)
    fhz = fhz_ebm(S, alb, A_olr, D)
    k = kopparapu_seff(teff)
    kin = np.where(sync, kin, 1.0)
    cons = (S <= k["rg"] * kin) & (S >= k["mg"])
    opt = (S <= k["rv"] * kin) & (S >= k["em"])
    p_hz = np.where(cons, phz[0], np.where(opt, phz[1], phz[2]))
    pco2 = 10 ** (P["lpco2_lo"] + (P["lpco2_hi"] - P["lpco2_lo"]) * u_pco2)
    eye = sync & (S < P["s_eye_max"]) & (pco2 >= p_co2_req(S))
    f_eye = np.clip(P["f_eye_c"] + P["f_eye_w"] * (2 * u_eye - 1), 0.0, 1.0)
    fhz = np.where(eye, np.maximum(fhz, f_eye), fhz)
    if P.get("rg_veto", True):
        fhz = np.where((~sync) & (S * S0_FLUX / 4.0 * (1 - alb) > P["F_rg"]), 0.0, fhz)
    return alb, fhz, p_hz, eye


def obs_v6(transit, rp, mp, rs, teff, a_au, per_d, teq, g, J, K, avail, dist_pc, albedo, u_iwa, z1, z2, z3, u_cl, mu, P):
    """S_obs de atm_obs_v5 con mu impuesto por la clase (R: 28/44; H: 2.3). Formulas identicas a v5.2."""
    teq0 = teff * np.sqrt(rs * RSUN / (a_au * AU)) * 0.25 ** 0.25
    tsm = tsm_scale(rp) * rp ** 3 * teq0 / (np.maximum(mp, 0.05) * rs ** 2) * 10 ** (-J / 5)
    delta = (rp * REARTH / (rs * RSUN)) ** 2
    esm = 4.29e6 * planck_ratio(1.10 * teq0, teff) * delta * 10 ** (-K / 5)
    m1, m2, m3 = np.exp(P["sig_model"] * z1), np.exp(P["sig_model"] * z2), np.exp(P["sig_model"] * z3)
    cl = P["cloud_min"] + (1 - P["cloud_min"]) * u_cl
    mu_k = np.where(rp < 1.5, 18.0, 2.3)
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
    d_ok = np.isfinite(dist_pc) & (dist_pc > 0)
    snr_r0, _, _ = snr_reflected(rp, a_au, np.where(d_ok, dist_pc, 1.0), J, albedo, u_iwa, P)
    snr_r = np.where(d_ok, snr_r0, 0.0) * m1
    p_nt = PHI(snr_r - P["snr_thr"])
    return np.where(transit, p_tr, p_nt) * avail, tsm, esm


def xuv_rocky(r, L, a, age, ms, U, Z, P, cfg, teff=None):
    """Fluencia XUV acumulada relativa a la Tierra (I/I_E) y L_X, L_EUV actuales [erg/s]. Mezcla log en M*:
    w = sigmoide((M - 0.35)/0.03): Pass+2025 (w->0) / Mors o Jackson (w->1)."""
    host = str(getattr(r, "hostname", ""))
    log_rsat = P["log_rsat"] + P["log_rsat_sd"] * Z[0]
    tsat = X.tsat_pass(ms, P.get("tsat_fac", 1.0) * np.exp(P["tsat_sd_ln"] * Z[1]))
    ff = P["fflare_lo"] + (P["fflare_hi"] - P["fflare_lo"]) * U[7]
    sat_now = host in X.STILL_SATURATED
    I_fc, Lx_fc = X.fluence_fc(L, a, age, ms, log_rsat, tsat, ff, saturated_now=sat_now, post_drop=P["post_drop"])
    Leuv_fc = X.leuv_from_lx(Lx_fc)
    w = sig((ms - P["blend_m"]) / P["blend_w"])
    if sat_now:
        w = np.zeros_like(ms)                    # Pass+2025: TRAPPIST-1 aun saturada; Mors no se aplica
    pct = P["pct_lo"] + (P["pct_hi"] - P["pct_lo"]) * U[8]
    mo = X.fluence_mors(L, a, age, ms, pct)
    if mo is None:                                # sin malla Mors -> fallback Jackson+2012 (documentado)
        cfg = dict(cfg, xuv_star="jackson")
        I_m, Lx_m, Leuv_m = I_fc, Lx_fc, Leuv_fc
    else:
        I_m, Lx_m, Leuv_m = mo
    if cfg.get("xuv_star") == "jackson":
        bv = (r.sy_bmag - r.sy_vmag) if (np.isfinite(getattr(r, "sy_bmag", np.nan)) and np.isfinite(getattr(r, "sy_vmag", np.nan))) else np.nan
        bv = np.full_like(L, bv) if np.isfinite(bv) else X.bv_from_teff(teff)     # snapshot sin B/V -> Pecaut & Mamajek (Teff)
        I_j, Lx_j = X.fluence_jackson(L, a, age, bv)
        okj = np.isfinite(I_j)
        I_m = np.where(okj, I_j, I_m); Lx_m = np.where(okj, Lx_j, Lx_m); Leuv_m = np.where(okj, X.leuv_from_lx(Lx_j), Leuv_m)
    lg = lambda x: np.log(np.maximum(x, 1e-300))
    I = np.exp((1 - w) * lg(I_fc) + w * lg(I_m))
    Lx = np.exp((1 - w) * lg(Lx_fc) + w * lg(Lx_m)); Leuv = np.exp((1 - w) * lg(Leuv_fc) + w * lg(Leuv_m))
    return I / X.I_XUV_EARTH, Lx, Leuv


def simulate_planet_v6(r, n, rng, cfg=CFG_V6, P=None, alpha=None, keep_draws=False):
    P = PARAMS_DEFAULT if P is None else P
    out = {}
    # --- estrella / orbita / planeta: IDENTICO a v5.2 (misma secuencia del generador principal) ---
    teff = split_normal(rng, r.st_teff, r.st_tefferr1, r.st_tefferr2, n, 0.02, 2300, 50000)
    rs = split_normal(rng, r.st_rad, r.st_raderr1, r.st_raderr2, n, 0.07, 0.05, None)
    ms = split_normal(rng, r.st_mass, r.st_masserr1, r.st_masserr2, n, 0.07, 0.05, None)
    if cfg.get("rs_floor", False):
        rs = apply_rs_floor(rs, ms, r.st_rad, r.st_mass, P.get("rs_floor_f", 0.5))   # v6.1 (fix 2); no consume draws
    L = rs ** 2 * (teff / 5772.0) ** 4
    if np.isfinite(r.st_age):
        age = split_normal(rng, r.st_age, r.st_ageerr1, r.st_ageerr2, n, 0.3, 0.01, 13.8)
    else:
        age = rng.uniform(1.0, 10.0, n)
    if np.isfinite(r.pl_orbper):
        Porb = split_normal(rng, r.pl_orbper, r.pl_orbpererr1, r.pl_orbpererr2, n, 1e-4, 1e-3, None)
        a = ((G * ms * MSUN) / (4 * np.pi ** 2)) ** (1 / 3) * (Porb * DAY) ** (2 / 3) / AU
    else:
        a = split_normal(rng, r.pl_orbsmax, r.pl_orbsmaxerr1, r.pl_orbsmaxerr2, n, 0.05, 1e-4, None)
        Porb = np.sqrt(a ** 3 / ms) * 365.25
    if np.isfinite(r.pl_orbeccen):
        e = split_normal(rng, r.pl_orbeccen, r.pl_orbeccenerr1, r.pl_orbeccenerr2, n, 0.0, 0.0, 0.95)
        if not (np.isfinite(r.pl_orbeccenerr1) or np.isfinite(r.pl_orbeccenerr2)):
            e = np.clip(r.pl_orbeccen + 0.05 * rng.standard_normal(n), 0, 0.95)
    else:
        e = np.zeros(n)
    rp_obs = mp_obs = None
    if r.radius_measured:
        rp = split_normal(rng, r.pl_rade, r.pl_radeerr1, r.pl_radeerr2, n, 0.10, 0.1, None); rp_obs = rp
    else:
        rp = None
    if r.mass_mode == "measured":
        mp = split_normal(rng, r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, n, 0.20, 0.01, None); mp_obs = mp
    elif r.mass_mode == "msini":
        msini = split_normal(rng, r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, n, 0.20, 0.01, None)
        cosi = rng.uniform(0, 1, n)
        mp = np.minimum(msini / np.sqrt(np.clip(1 - cosi ** 2, 1e-4, 1)), 26600.0); mp_obs = mp
        _msini = msini
    else:
        mp = None
    if rp is None and mp is None:
        rp = split_normal(rng, r.pl_rade, r.pl_radeerr1, r.pl_radeerr2, n, 0.15, 0.1, None); rp_obs = rp
    if mp is None:
        mp = ck_mass(rp, rng)
    if rp is None:
        rp = ck_radius(mp, rng)
    # --- prueba de prior (b): masa verdadera desde M sin i (solo afecta a mp_obs; el generador principal no cambia) ---
    ms_mode = cfg.get("msini", "iso")
    if r.mass_mode == "msini" and ms_mode != "iso":
        if ms_mode == "min":            # M = M sin i (cota inferior; sin i = 1)
            mp_obs = np.minimum(_msini, 26600.0)
        elif ms_mode == "howard":
            # posterior de la masa verdadera con prior dN/dlogM ~ M^alpha (Howard+2010, alpha = -0.48): el muestreo isotropico
            # equivale a prior plano en M, asi que se repondera por (M/Msini)^(alpha-1) = (1/sin i)^(alpha-1). Solo se pondera la
            # proyeccion (no la incertidumbre de Msini, que ya es la del catalogo) para no divergir hacia M -> 0.
            w = (mp_obs / np.maximum(_msini, 1e-6)) ** (P.get("mf_alpha", -0.48) - 1.0); w = w / w.sum()
            idx = np.random.default_rng(zlib.crc32(str(r.pl_name).encode()) + 3).choice(n, n, p=w)
            mp_obs = mp_obs[idx]
    S = L / a ** 2
    alb_sh = albedo_teff(teff, rng.uniform(0.0, 1.0, n))
    nm = str(r.pl_name).encode()
    rng5 = np.random.default_rng(zlib.crc32(nm) + 1)
    u5 = rng5.uniform(size=(4, n)); z_R0 = rng5.standard_normal(n)
    rng6 = np.random.default_rng(zlib.crc32(nm) + 2)
    U = rng6.uniform(size=(20, n)); Z = rng6.standard_normal((8, n))
    locked = (r.st_mass < 0.6) & (np.median(a) < 0.2)          # v5.2 (se conserva para S_mag)
    A_olr = rng.uniform(P.get("olr_Amin", 150.0), 203.3, n)

    # --- D1: clima ---
    cm = cfg["climate_mode"]
    if cm == "v52":
        albedo, fhz, p_hz, eye = climate_v5(S, teff, locked, alb_sh, A_olr, u5[0], u5[1], u5[2], u5[3], P=P)
        p_sync = np.full(n, float(locked and np.median(teff) < P["teff_md"]))
        sync = np.full(n, bool(locked)) & (teff < P["teff_md"])
    else:
        m_lock = P["lock_m_lo"] + (P["lock_m_hi"] - P["lock_m_lo"]) * U[0]
        a_lock = P["lock_a_lo"] + (P["lock_a_hi"] - P["lock_a_lo"]) * U[1]
        p_sync = sig((m_lock - ms) / P["lock_m_w"]) * sig((a_lock - a) / P["lock_a_w"])
        if cm == "teff":
            t_c = P["tc_lo"] + (P["tc_hi"] - P["tc_lo"]) * U[5]; dT = P["dT_lo"] + (P["dT_hi"] - P["dT_lo"]) * U[6]
            p_sync = p_sync * sig((t_c - teff) / dT)
        sync = U[2] < p_sync
        p_lo = P["prot_lo_lo"] + (P["prot_lo_hi"] - P["prot_lo_lo"]) * U[3]
        p_hi = P["prot_hi_lo"] + (P["prot_hi_hi"] - P["prot_hi_lo"]) * U[4]
        e3d = np.clip((Porb - p_lo) / (p_hi - p_lo), 0.0, 1.0)
        if cm == "teff":     # variante: eficacia/k_in como v5.2 (aleatorios), solo cambia la frontera
            eff = P["yang_eff_lo"] + (1 - P["yang_eff_lo"]) * u5[0]; kin = P["kin_lo"] + (P["kin_hi"] - P["kin_lo"]) * u5[1]
        else:
            eff = P["yang_eff_lo"] + (1 - P["yang_eff_lo"]) * e3d
            kin = P["kin_hi"] * (1 - 0.25 * (1 - e3d))           # 2 x [1 - 0.25 (1 - e)]: 2.0 lento, 1.5 rapido
        albedo, fhz, p_hz, eye = climate_v6(S, teff, sync, alb_sh, A_olr, eff, kin, u5[2], u5[3], P)
    teq = ((1 - albedo) * L * LSUN / (16 * np.pi * SIGMA * (a * AU) ** 2)) ** 0.25
    S_clima = fhz * p_hz

    tr = (r.tran_flag == 1)
    z = rng.standard_normal((3, n)); u_cl = rng.uniform(size=n); u_mu = rng.uniform(size=n)
    u_iwa = np.random.default_rng(zlib.crc32(nm)).uniform(size=n)
    avail = TIER_AVAIL[r.tier]

    # --- S_estrella (D3: solo riesgo superficial ACTUAL; formula v5.2 sin cambios) ---
    tsat5 = 0.1 if r.st_mass >= 0.6 else P.get("tsat_M", 1.0)
    rx = np.where(age < tsat5, 10 ** -3.5, 10 ** -3.5 * (age / tsat5) ** -1.23)
    rx_sun = 10 ** -3.5 * (4.6 / 0.1) ** -1.23
    fxuv_rel = rx * S / rx_sun
    stress = logistic_log(fxuv_rel, P.get("xuv_mid", 10.0))
    flare = P.get("flare_late", 0.5) if r.st_mass < 0.35 else (P.get("flare_early", 0.25) if r.st_mass < 0.6 else 0.0)
    S_estrella = (1 - stress) * (1 - flare)

    def s_geo(mp_, rp_):
        n_orb = 2 * np.pi / (Porb * DAY)
        k2Q = 0.3 / 100.0
        Edot = 10.5 * k2Q * G * (ms * MSUN) ** 2 * (rp_ * REARTH) ** 5 * n_orb * e ** 2 / (a * AU) ** 6
        q_tide = Edot / (4 * np.pi * (rp_ * REARTH) ** 2)
        q_rad = 0.087 * (mp_ / rp_ ** 2) * 2.0 ** (-(age - 4.5) / 3.0)
        lq = np.log10(np.maximum(q_rad + q_tide, 1e-6))
        S_regime = trapezoid(lq, -2.0, np.log10(0.03), 0.0, 1.0)
        if cfg["d3"]:
            return S_regime                                         # D3: erosion fuera de S_geo
        fxuv_abs = fxuv_rel * F_XUV_EARTH
        mdot = 0.1 * np.pi * (rp_ * REARTH) ** 3 * fxuv_abs / (G * mp_ * MEARTH)
        t_loss = 0.01 * mp_ * MEARTH / np.maximum(mdot, 1e-30) / (1e9 * YR)
        risk = 1 - logistic_log(t_loss, 1.0)
        subnep = (rp_ >= 1.6) & (rp_ <= 4.0)
        return S_regime * np.where(subnep, 1 - 0.5 * risk, 1.0)

    aB, bB = (2.0, 3.0) if locked else (2.0, 2.0)
    aB = np.full(n, aB); bB = np.where(mp < 0.3, bB + 1.0, bB)
    S_mag = rng.beta(aB, bB)                                      # v5.2 (masa v5.2 para no alterar la secuencia)
    W = rng.dirichlet(DIRICHLET_ALPHA if alpha is None else alpha, n).T

    # --- D2: XUV acumulado ---
    if cfg["xuv"] or cfg["branches"]:
        I_rel, Lx_now, Leuv_now = xuv_rocky(r, L, a, age, ms, U, Z, P, cfg, teff)
    q = lambda x, p: float(np.nanpercentile(x, p)) if np.size(x) else np.nan

    def ret_rocky(mp_, rp_):
        vesc = u.EARTH_VESC * np.sqrt(np.maximum(mp_, 1e-6) / rp_)
        if not cfg["xuv"]:
            return p_retain(S, mp_, rp_, P), np.full(n, np.nan)
        if cfg["shore_model"] == "bt3d":
            return X.p_ret_bt3d(vesc, S, L, P["shore_sig"]), X.arm(vesc, I_rel)
        return X.p_ret_xuv(vesc, I_rel, P["shore_logc"], P["shore_sig"]), X.arm(vesc, I_rel)

    blocks = {}
    if not cfg["branches"]:
        # ---------- rama unica (v5.2 + modulos D1/D2/D3 opcionales) ----------
        g = G * mp * MEARTH / (rp * REARTH) ** 2
        S_atm5, S_obs, tsm, esm, pret5, w_h2 = atm_obs_v5(tr, rp, mp, rs, teff, a, Porb, teq, g, r.J_use, r.K_use, S, avail,
                                                         r.sy_dist, albedo, u_iwa, z[0], z[1], z[2], u_cl, u_mu, z_R0,
                                                         bool(r.dens_ok), P)
        pret, armv = ret_rocky(mp, rp)
        S_atm = pret if cfg["xuv"] else S_atm5
        S_fis = np.sqrt(trapezoid(rp, 0.3, 0.5, 1.6, 2.5) * trapezoid(mp, 0.03, 0.1, 5.0, 10.0))
        Smat = np.vstack([S_fis, S_clima, S_atm, S_estrella, s_geo(mp, rp), S_mag, S_obs])
        uhi = aggregate_v5(Smat, W)
        blocks["combined"] = (uhi, Smat, np.ones(n, bool), dict(P_retain_med=q(pret, 50), ARM_med=q(armv, 50), TSM_med=q(tsm, 50) if tr else np.nan,
                                                               ESM_med=q(esm, 50) if tr else np.nan, wH2_med=q(w_h2, 50)))
        cls = np.zeros(n, int); evid = "v52"
    else:
        # ---------- D4: clases por draw ----------
        pi = cfg["class_prior"]
        z1, z2 = Z[2], Z[3]
        if rp_obs is not None and mp_obs is not None and bool(r.dens_ok):
            evid = "densidad"; mpR = mpH = mp_obs; rpR = rpH = rp_obs
            cls = C.classify_density(mp_obs, rp_obs, teq, U[9], pi, cfg.get("rock_line", "h2o50"))
        elif rp_obs is not None:
            evid = "radio"; rpR = rpH = rp_obs
            if mp_obs is not None:
                mpR = mpH = mp_obs
            elif cfg.get("mr", "otegi") == "ck":       # prueba (c): Chen & Kipping 2017 (masa v5.2) para ambas clases
                mpR = mpH = mp
            else:
                mpR = C.otegi_mass(rp_obs, True, z1, z2); mpH = C.otegi_mass(rp_obs, False, z1, z2)
            w_v = P["wv_lo"] + (P["wv_hi"] - P["wv_lo"]) * U[11]
            cls, _ = C.classify_radius(rp_obs, ms, U[9], U[10], Z[4], w_v, mpH, pi)
        else:
            evid = "prior"; mpR = mpH = mp_obs
            if cfg.get("mr", "otegi") == "ck":
                rpR = rpH = ck_radius(mp_obs, np.random.default_rng(zlib.crc32(nm) + 4))   # prueba (c)
            else:
                rpR = C.otegi_radius(mp_obs, True, z1, z2); rpH = C.otegi_radius(mp_obs, False, z1, z2)
            mv = P["mv_lo"] + (P["mv_hi"] - P["mv_lo"]) * U[12]
            cls, _ = C.classify_mass(mp_obs, U[9], U[10], mv, P["mv_s"], rpH, pi)
        # ---------- rama Rocky ----------
        gR = G * mpR * MEARTH / (rpR * REARTH) ** 2
        muR = np.where(u_mu < 0.5, 28.0, 44.0)
        S_obsR, tsmR, esmR = obs_v6(tr, rpR, mpR, rs, teff, a, Porb, teq, gR, r.J_use, r.K_use, avail, r.sy_dist, albedo, u_iwa,
                                   z[0], z[1], z[2], u_cl, muR, P)
        pretR, armR = ret_rocky(mpR, rpR)
        S_fisR = np.sqrt(np.clip((rpR - 0.3) / 0.2, 0, 1) * trapezoid(mpR, 0.03, 0.1, 5.0, 10.0))   # sin borde superior en R
        SmR = np.vstack([S_fisR, S_clima, pretR, S_estrella, s_geo(mpR, rpR), S_mag, S_obsR])
        uR = aggregate_v5(SmR, W) * (cls == 0)
        # ---------- rama Hycean ----------
        lph2 = P["lph2_lo"] + (P["lph2_hi"] - P["lph2_lo"]) * U[13]
        if cfg["hycean_climate"] == "madhu":
            teq0 = (L * LSUN / (16 * np.pi * SIGMA * (a * AU) ** 2)) ** 0.25          # A_B = 0 (SUPUESTO)
            S_climH = (teq0 <= P["madhu_teq"]).astype(float)
        else:
            f_cl = P["fcl_lo"] + (P["fcl_hi"] - P["fcl_lo"]) * U[14]
            S_climH = (S <= f_cl * X.s_crit_innes(teff, 10 ** lph2) * (1 - 0.12) / (1 - P["hyc_alb"])).astype(float)   # limite en flujo absorbido (SUPUESTO)
        fx_now = (Lx_now + Leuv_now) / (4 * np.pi * (a * u.AU * 100) ** 2)           # erg cm^-2 s^-1
        t_loss = X.envelope_t_loss_gyr(rpH, mpH, 10 ** lph2, fx_now, P.get("eta_fac", 1.0) * np.exp(P["eta_sd"] * Z[5]))
        S_atmH = PHI(np.log(np.maximum(t_loss, 1e-30) / P["tloss_ref"]) / P["tloss_sig"])
        S_fisH = trapezoid(mpH, P["hf_a"], P["hf_b"], P["hf_c"], P["hf_d"])
        teqH = ((1 - P["hyc_alb"]) * L * LSUN / (16 * np.pi * SIGMA * (a * AU) ** 2)) ** 0.25
        gH = G * mpH * MEARTH / (rpH * REARTH) ** 2
        S_obsH, tsmH, esmH = obs_v6(tr, rpH, mpH, rs, teff, a, Porb, teqH, gH, r.J_use, r.K_use, avail, r.sy_dist,
                                   np.full(n, P["hyc_alb"]), u_iwa, z[0], z[1], z[2], u_cl, 2.3, P)
        SmH = np.vstack([S_fisH, S_climH, S_atmH, S_estrella, s_geo(mpH, rpH), S_mag, S_obsH])
        uH = aggregate_v5(SmH, W) * (cls == 1)
        mR, mH = cls == 0, cls == 1
        blocks["rocky"] = (uR, SmR, mR, dict(P_retain_med=q(pretR[mR], 50), ARM_med=q(armR[mR], 50), ARM_all_med=q(armR, 50),
                                             TSM_med=q(tsmR[mR], 50) if tr else np.nan, Rp_cls_med=q(rpR[mR], 50), Mp_cls_med=q(mpR[mR], 50)))
        blocks["hycean"] = (uH, SmH, mH, dict(tloss_med=q(t_loss[mH], 50), tloss_all_med=q(t_loss, 50), P_clima_H=float(np.mean(S_climH)),
                                              TSM_med=q(tsmH[mH], 50) if tr else np.nan, Rp_cls_med=q(rpH[mH], 50), Mp_cls_med=q(mpH[mH], 50)))
        blocks["combined"] = (np.where(mR, uR, np.where(mH, uH, 0.0)), np.where(mR[None], SmR, SmH), mR | mH, {})
    # --- ESI (como v5.2, con la masa/radio v5.2) ---
    rho = mp / rp ** 3; vrel = np.sqrt(mp / rp); Ts = teq + 33.0
    et = lambda x, x0, w: (1 - np.abs((x - x0) / (x + x0))) ** (w / 2.0)
    esi_i = et(rp, 1.0, 0.57) * et(rho, 1.0, 1.07); esi_s = et(vrel, 1.0, 0.70) * et(Ts, 288.0, 5.58)
    esi = np.sqrt(esi_i * esi_s)
    common = dict(ESI_med=q(esi, 50), ESI_p16=q(esi, 16), ESI_p84=q(esi, 84), albedo_med=q(albedo, 50), f_eye=float(np.mean(eye)),
                  P_sync=float(np.mean(p_sync)), f_sync=float(np.mean(sync)), fHZ_med=q(fhz, 50), Rp_med=q(rp, 50), Mp_med=q(mp, 50),
                  S_med=q(S, 50), Teq_med=q(teq, 50), a_med=q(a, 50), FXUV_rel_med=q(fxuv_rel, 50),
                  P_R=float(np.mean(cls == 0)), P_H=float(np.mean(cls == 1)), P_G=float(np.mean(cls == 2)), evidence=evid)
    if cfg["xuv"] or cfg["branches"]:
        common.update(IXUV_rel_med=q(I_rel, 50), IXUV_rel_p16=q(I_rel, 16), IXUV_rel_p84=q(I_rel, 84))
    for b, (uhi, Smat, mask, extra) in blocks.items():
        o = dict(common); o.update(extra)
        alive = uhi > 0
        o.update(UHI_med=q(uhi, 50), UHI_p16=q(uhi, 16), UHI_p84=q(uhi, 84), UHI_p2_5=q(uhi, 2.5), UHI_p97_5=q(uhi, 97.5),
                 UHI_mean=float(np.mean(uhi)), P_UHI_gt0=float(np.mean(alive)), P_alive=float(np.mean(alive)),
                 UHI_cond_med=q(uhi[alive], 50) if alive.any() else 0.0,
                 UHI_cond_q25=q(uhi[alive], 25) if alive.any() else 0.0,
                 UHI_cond_q75=q(uhi[alive], 75) if alive.any() else 0.0)
        for name, arr in zip(SUB, Smat):
            o[name + "_med"] = q(arr[mask], 50) if mask.any() else np.nan     # mediana condicionada a la clase
        o["UHI_rank"] = o["P_alive"] * o["UHI_cond_med"]
        if keep_draws and o["P_alive"] > 0:
            o["_draws"] = uhi.astype(np.float32)
        out[b] = o
    return out


def _tiers(d):
    a_nom = np.where(d.pl_orbper.notna(), ((G * d.st_mass * MSUN) / (4 * np.pi ** 2)) ** (1 / 3) * (d.pl_orbper * DAY) ** (2 / 3) / AU,
                     d.pl_orbsmax)
    d["a_nom"] = a_nom
    mp_nom = d.pl_bmasse.copy(); rp_nom = d.pl_rade.copy()
    need_m = mp_nom.isna()
    if need_m.any():
        mp_nom[need_m] = np.where(rp_nom[need_m] < _CK_R1, (rp_nom[need_m] / 1.008) ** (1 / 0.279), _CK_T1 * (rp_nom[need_m] / _CK_R1) ** (1 / 0.589))
    need_r = rp_nom.isna()
    if need_r.any():
        m = mp_nom[need_r].values
        rp_nom[need_r] = np.where(m < _CK_T1, 1.008 * m ** 0.279, np.where(m < _CK_T2, _CK_R1 * (m / _CK_T1) ** 0.589, _CK_R2 * (m / _CK_T2) ** -0.044))
    d["mp_nom"] = mp_nom; d["rp_nom"] = rp_nom
    d["tier"] = [nominal_tier(row) for row in d.itertuples()]      # tier nominal de v5.2 (usa V3['snr_thr'] global)
    return d


def finalize(rows, draws, joint=True, strict=True):
    res = pd.DataFrame(rows)
    res["P_top10"] = 0.0; res["rank_med"] = np.nan; res["rank_ci5"] = np.nan; res["rank_ci95"] = np.nan
    if joint and draws:
        idx = np.array(sorted(draws)); M = np.vstack([draws[i] for i in idx])
        npos = (M > 0).sum(0)
        order = np.argsort(-M, axis=0, kind="stable"); rk = np.empty_like(order)
        np.put_along_axis(rk, order, np.arange(1, M.shape[0] + 1)[:, None].repeat(M.shape[1], 1), axis=0)
        rk = np.where(M > 0, rk, npos[None, :] + 1)
        # v6 (strict): P_top10 cuenta solo draws vivos. v5.2 contaba tambien rango muerto = vivos+1 <= 10 cuando hay <10 vivos por draw
        # (artefacto irrelevante en v5.2 con ~500 vivos, pero dominante en la rama Hycean). CFG_COMPAT conserva la convencion v5.2.
        res.loc[idx, "P_top10"] = (((rk <= 10) & (M > 0)) if strict else (rk <= 10)).mean(1)
        res.loc[idx, "rank_med"] = np.median(rk, 1)
        res.loc[idx, "rank_ci5"] = np.percentile(rk, 5, 1); res.loc[idx, "rank_ci95"] = np.percentile(rk, 95, 1)
        del M, order, rk
    res["ESI_rank"] = res.ESI_med.rank(ascending=False, method="min").astype("Int64")
    res = res.sort_values(["UHI_rank", "UHI_mean"], ascending=False).reset_index(drop=True)
    res["rank_pos"] = np.arange(1, len(res) + 1)
    first = ["rank_pos", "planet", "UHI_rank", "P_alive", "UHI_cond_med", "UHI_cond_q25", "UHI_cond_q75", "P_top10",
             "rank_med", "rank_ci5", "rank_ci95", "UHI_med", "UHI_p16", "UHI_p84", "UHI_p2_5", "UHI_p97_5", "UHI_mean",
             "P_UHI_gt0", "tier", "ESI_med", "ESI_rank"]
    return res[first + [c for c in res.columns if c not in first]].sort_values("rank_pos")


def run_v6(df, n_draws=10000, seed=42, cfg=CFG_V6, P=None, progress=None, alpha=None, joint=True, names=None):
    """Devuelve dict rama -> DataFrame (rocky, hycean, combined; o solo combined si branches=False)."""
    d = _tiers(prepare(df))
    if names is not None:
        d = d[d.pl_name.isin(names)].reset_index(drop=True)
    rng = np.random.default_rng(seed)
    rows = {}; draws = {}
    for i, row in enumerate(d.itertuples()):
        res = simulate_planet_v6(row, n_draws, rng, cfg=cfg, P=P, alpha=alpha, keep_draws=joint)
        meta = dict(planet=row.pl_name, host=row.hostname, method=row.discoverymethod, transiting=int(row.tran_flag),
                    tier=row.tier, tier_vector=TIER_VECTOR, mass_source=row.mass_mode,
                    controversial=int(row.pl_controv_flag) if np.isfinite(row.pl_controv_flag) else 0,
                    st_teff=row.st_teff, st_mass=row.st_mass, sy_dist=row.sy_dist, age_known=int(np.isfinite(row.st_age)),
                    dens_ok=bool(row.dens_ok), radius_measured=bool(row.radius_measured))
        for b, o in res.items():
            o.update(meta)
            if "_draws" in o:
                draws.setdefault(b, {})[i] = o.pop("_draws")
            rows.setdefault(b, []).append(o)
        if progress and i % 500 == 0:
            progress(i, len(d))
    strict = dict(cfg) != dict(CFG_COMPAT)
    return {b: finalize(rows[b], draws.get(b, {}), joint, strict) for b in rows}
