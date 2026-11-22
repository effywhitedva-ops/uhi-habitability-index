"""Recalculo a mano INDEPENDIENTE (no importa uhi_code / uhi_code_v6 / xuv_history / composition) de GJ 273 b, Proxima Cen b
y TRAPPIST-1 d. Constantes y formulas tecleadas desde las fuentes (Pass+2025 Ec. 6-11 y Tabla 1; Kepler III; Zeng+2019; Otegi+2020).
Valor nominal + MC propio (semilla 7, 200k draws) para las medianas."""
import numpy as np, pandas as pd, json
from scipy.stats import norm
from scipy.interpolate import RegularGridInterpolator
G = 6.674e-11; MSUN = 1.98847e30; AU = 1.495978707e11; DAY = 86400.0; LSUN = 3.828e33; GYR = 3.15576e16
I_E = 2.0e18; VE = 11.186
rng = np.random.default_rng(7); N = 200_000
cat = pd.read_csv("data/pscomppars_2026-11-18.csv").set_index("pl_name")
sig = lambda x: 1 / (1 + np.exp(-x))
def sn(v, e1, e2, frac, lo=None, hi=None, n=N):
    su = abs(e1) if np.isfinite(e1) and e1 != 0 else abs(v) * frac; sd = abs(e2) if np.isfinite(e2) and e2 != 0 else abs(v) * frac
    z = rng.standard_normal(n); x = np.where(z >= 0, v + z * su, v + z * sd); return np.clip(x, lo if lo is not None else -np.inf, hi if hi is not None else np.inf)
# Pass+2025 Tabla 1 (f_L) y Ec. 7 (t_sat)
PM = [0.10, 0.12, 0.15, 0.20, 0.25, 0.30]; PFL = [1.19, 1.17, 1.13, 1.13, 1.16, 1.22]
g = np.load("data/mors_grid.npz")
ip = {k: RegularGridInterpolator((g["masses"], g["pcts"], np.log(g["ages"])), np.log(np.maximum(g[k], 1e-300)), bounds_error=False, fill_value=None) for k in ("CLx", "Lbol")}
def mors_I(L, a_au, age, m, pct):
    pts = np.c_[np.clip(m, .1, 1.2), np.clip(pct, 5, 95), np.log(np.clip(age, g["ages"][0], g["ages"][-1]))]
    E = np.exp(ip["CLx"](pts)) * L * LSUN / np.exp(ip["Lbol"](pts))
    return E / (4 * np.pi * (a_au * AU * 100) ** 2) / I_E
# Zeng+2019
Z = lambda f: np.loadtxt(f"data/zeng2019/{f}")
SIL, W50 = Z("massradiusmgsio3.txt"), Z("massradius_50percentH2O_300K_1mbar.txt")
rz = lambda tb, m: np.interp(np.log(m), np.log(tb[:, 0]), tb[:, 1])
rhmax = lambda m: 2.3 + 0.3 * np.log2(np.maximum(m, 1e-3) / 5.0)
out = {}
for p in ["GJ 273 b", "Proxima Cen b", "TRAPPIST-1 d"]:
    r = cat.loc[p]; o = {}
    # --- nominal ---
    L0 = r.st_rad ** 2 * (r.st_teff / 5772) ** 4
    a0 = (G * r.st_mass * MSUN * (r.pl_orbper * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / AU
    o["L_nom"] = L0; o["a_nom_AU"] = a0; o["S_nom"] = L0 / a0 ** 2
    o["Teq_A0_nom"] = 278.6 * o["S_nom"] ** 0.25
    # --- MC propio ---
    teff = sn(r.st_teff, r.st_tefferr1, r.st_tefferr2, .02, 2300, 50000); rs = sn(r.st_rad, r.st_raderr1, r.st_raderr2, .07, .05)
    ms = sn(r.st_mass, r.st_masserr1, r.st_masserr2, .07, .05); L = rs ** 2 * (teff / 5772) ** 4
    age = sn(r.st_age, r.st_ageerr1, r.st_ageerr2, .3, .01, 13.8) if np.isfinite(r.st_age) else rng.uniform(1, 10, N)
    Porb = sn(r.pl_orbper, r.pl_orbpererr1, r.pl_orbpererr2, 1e-4, 1e-3)
    a = (G * ms * MSUN * (Porb * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3) / AU; S = L / a ** 2
    o["S_MCmed"] = float(np.median(S))
    # --- D1 p_sync: E sobre M_lock~U(.55,.65), a_lock~U(.15,.25) ---
    mlk = rng.uniform(.55, .65, N); alk = rng.uniform(.15, .25, N)
    ps = sig((mlk - ms) / .05) * sig((alk - a) / .03); o["p_sync"] = float(ps.mean())
    o["p_sync_nominal_analytic"] = float(np.mean(sig((np.linspace(.55, .65, 401)[:, None] - r.st_mass) / .05) * sig((np.linspace(.15, .25, 401)[None, :] - a0) / .03)))
    plo = rng.uniform(4, 6, N); phi = rng.uniform(15, 25, N); e3 = np.clip((Porb - plo) / (phi - plo), 0, 1)
    o["e3D_med"] = float(np.median(e3)); o["k_in_med"] = float(np.median(2 * (1 - .25 * (1 - e3))))
    # --- D2: Pass+2025 FC ---
    lrs = -3.05 + .06 * rng.standard_normal(N); ts = np.maximum((5.9 - 15.4 * np.maximum(ms, .1)) * np.exp(.3 * rng.standard_normal(N)), .01)
    ff = rng.uniform(1, 1.5, N); fl = np.interp(ms, PM, PFL); Lc = L * LSUN; Rs = 10 ** lrs
    satnow = r.hostname == "TRAPPIST-1"
    ton = age if satnow else np.minimum(ts, age)
    E = ff * Rs * fl * Lc * ton * GYR + (0 if satnow else 1e-2 * Rs * Lc * np.maximum(age - ts, 0) * GYR)
    Ifc = E / (4 * np.pi * (a * AU * 100) ** 2) / I_E
    w = np.zeros(N) if satnow else sig((ms - .35) / .03)
    Im = mors_I(L, a, age, ms, rng.uniform(5, 95, N))
    I = np.exp((1 - w) * np.log(Ifc) + w * np.log(Im))
    o["I_rel_Pass_only_med"] = float(np.median(Ifc)); o["blend_w_med"] = float(np.median(w)); o["I_rel_med"] = float(np.median(I))
    o["I_rel_p16"], o["I_rel_p84"] = [float(x) for x in np.percentile(I, [16, 84])]
    # nominal analitico Pass (f_flare 1.25, R_sat 10^-3.05, t_sat nominal, edad nominal o 5.5 Gyr)
    agen = r.st_age if np.isfinite(r.st_age) else 5.5; tsn = 5.9 - 15.4 * max(r.st_mass, .1)
    En = 1.25 * 10 ** -3.05 * np.interp(r.st_mass, PM, PFL) * L0 * LSUN * GYR * (agen if satnow else min(tsn, agen))
    o["I_rel_Pass_nominal"] = float(En / (4 * np.pi * (a0 * AU * 100) ** 2) / I_E); o["t_sat_nominal_Gyr"] = tsn
    # --- masa/radio y clase ---
    if r.pl_bmassprov == "Msini" and r.tran_flag == 0:
        msini = sn(r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, .2, .01); ci = rng.uniform(0, 1, N)
        M = np.minimum(msini / np.sqrt(np.clip(1 - ci ** 2, 1e-4, 1)), 26600)
        z1, z2 = rng.standard_normal(N), rng.standard_normal(N)
        rR = (M / np.maximum(.9 + .06 * z1, .5)) ** (1 / (3.45 + .12 * z2)); rV = (M / np.maximum(1.74 + .38 * z1, .5)) ** (1 / (1.58 + .10 * z2))
        mv = rng.uniform(3, 6, N); PR = sig((np.log(mv) - np.log(M)) / .3)
        isR = rng.uniform(size=N) < PR; canH = rV <= rhmax(M); isH = ~isR & canH & (rng.uniform(size=N) < .5)
        o.update(evidence="prior (masa)", M_med=float(np.median(M)), P_R=float(isR.mean()), P_H=float(isH.mean()), P_G=float(1 - isR.mean() - isH.mean()),
                 P_R_analytic_mean=float(PR.mean()), frac_canH_given_volatile=float(canH.mean()))
        RpR = rR; MpR = M
    else:
        M = sn(r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, .2, .01); R = sn(r.pl_rade, r.pl_radeerr1, r.pl_radeerr2, .1, .1)
        rw, rsil, rh = rz(W50, M), rz(SIL, M), rhmax(M)      # T_eq < 300 K -> tabla 300 K
        cR, cH, cG = R <= rw, (R > rsil) & (R <= rh), R > rsil
        tot = cR.astype(float) + cH + cG; u = rng.uniform(size=N) * np.where(tot > 0, tot, 1)
        wR = np.where(tot > 0, cR, 1.0)
        cls = np.where(u < wR, 0, np.where(u < wR + cH, 1, 2))
        o.update(evidence="densidad", R_50H2O_nom=float(rz(W50, r.pl_bmasse)), R_MgSiO3_nom=float(rz(SIL, r.pl_bmasse)), R_Hmax_nom=float(rhmax(r.pl_bmasse)),
                 frac_R_compat=float(cR.mean()), frac_H_compat=float(cH.mean()), frac_G_compat=float(cG.mean()),
                 P_R=float((cls == 0).mean()), P_H=float((cls == 1).mean()), P_G=float((cls == 2).mean()))
        RpR, MpR = R, M
    vesc = VE * np.sqrt(MpR / RpR); ARM = 4 * np.log10(vesc) - np.log10(I) - 3.16
    o["vesc_med"] = float(np.median(vesc)); o["ARM_med"] = float(np.median(ARM)); o["P_ret_med"] = float(np.median(norm.cdf(ARM * np.log(10))))
    # Hycean clima (Innes): P(S <= f_cl * S_crit(P_H2, Teff)), log10 P ~ U(0,2), f_cl ~ U(1,2)
    s1 = np.interp(teff, [3457, 5772], [.320, .389]); s10 = np.interp(teff, [3457, 5772], [.0852, .0676]); lp = rng.uniform(0, 2, N)
    sc = 10 ** (np.log10(s1) + (np.log10(s10) - np.log10(s1)) * lp); o["P_clima_H_innes"] = float((S <= rng.uniform(1, 2, N) * sc).mean())
    out[p] = o
pipe = {b: pd.read_csv(f"results/uhi_v6_{b}_main_5s.csv").set_index("planet") for b in ("rocky", "hycean")}
comp = []
MAP = [("S_MCmed", "rocky", "S_med"), ("p_sync", "rocky", "P_sync"), ("I_rel_med", "rocky", "IXUV_rel_med"), ("I_rel_p16", "rocky", "IXUV_rel_p16"),
       ("I_rel_p84", "rocky", "IXUV_rel_p84"), ("ARM_med", "rocky", "ARM_all_med"), ("P_R", "rocky", "P_R"), ("P_H", "rocky", "P_H"), ("P_G", "rocky", "P_G"),
       ("P_clima_H_innes", "hycean", "P_clima_H")]
for p, o in out.items():
    for k, b, c in MAP:
        pv = float(pipe[b].loc[p, c]); hv = o[k]
        comp.append(dict(planet=p, quantity=c, hand=round(hv, 5), pipeline=round(pv, 5), abs_diff=round(hv - pv, 5),
                         rel_diff=round((hv - pv) / pv, 4) if pv != 0 else np.nan))
cp = pd.DataFrame(comp); cp.to_csv("results/analysis_phase1/hand_recompute_vs_pipeline.csv", index=False)
json.dump(out, open("results/analysis_phase1/hand_recompute_values.json", "w"), indent=1, default=float)
pd.set_option("display.width", 200); print(cp.to_string())
for p, o in out.items(): print(p, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in o.items()})
