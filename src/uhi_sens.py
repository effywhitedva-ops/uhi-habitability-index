"""
Analisis de sensibilidad global del UHI v5 (implementacion asistida por agente de IA; ver docs/SOFTWARE_ACKNOWLEDGMENTS.md).
Salida por planeta: UHI_rank = P_alive x mediana condicional (hurdle v5).
Version vectorizada (planetas x draws) con numeros aleatorios comunes: los draws de las
mediciones del catalogo se generan UNA vez; cada punto del espacio de supuestos solo
re-transforma esos draws. Asi la diferencia entre evaluaciones se debe solo a los supuestos.
La fisica es la misma de uhi_code.py; los numeros que alli son fijos aqui son parametros.
"""
import zlib
import numpy as np, pandas as pd
from scipy.special import betaincinv
import uhi_code as U

# nombre, valor base, min, max, descripcion
PARAMS = [
 ("fR_a", 0.3, 0.2, 0.4, "S_fis: borde inferior (0) del trapecio de radio [R_E]"),
 ("fR_b", 0.5, 0.4, 0.8, "S_fis: inicio del plateau de radio [R_E]"),
 ("fR_c", 1.6, 1.4, 2.0, "S_fis: fin del plateau de radio [R_E]"),
 ("fR_d", 2.5, 2.0, 3.0, "S_fis: borde superior (0) del trapecio de radio [R_E]"),
 ("fM_b", 0.1, 0.05, 0.3, "S_fis: inicio del plateau de masa [M_E] (borde 0 = 0.3*fM_b)"),
 ("fM_c", 5.0, 3.0, 8.0, "S_fis: fin del plateau de masa [M_E]"),
 ("fM_d", 10.0, 8.0, 15.0, "S_fis: borde superior (0) del trapecio de masa [M_E]"),
 ("pHZ_opt", 0.5, 0.2, 0.8, "S_clima: prior Kopparapu en HZ optimista (fuera de la conservadora)"),
 ("pHZ_out", 0.1, 0.0, 0.3, "S_clima: prior Kopparapu fuera de la HZ optimista"),
 ("olr_Amin", 150.0, 120.0, 190.0, "S_clima: limite inferior del intercepto OLR A (invernadero max) [W/m2]"),
 ("ebm_D", 0.649, 0.4, 0.9, "S_clima: difusion meridional del EBM [W/m2/C]"),
 ("f_land_max", 0.5, 0.2, 1.0, "Albedo v2: fraccion de tierra ~ U(0, f_land_max)"),
 ("alb_shift", 0.0, -0.05, 0.05, "Albedo v2: desplazamiento aditivo (incertidumbre de nubes)"),
 ("late_m_scale", 1.0, 0.8, 1.0, "Albedo v2: factor para Teff < 3477 K (1 = valor M3V constante)"),
 ("t_hci", 10.0, 5.0, 50.0, "v4 no transitantes: horas ELT/PCS (Kasper+2021 10 h; su simulacion de rendimiento usa 50 h)"),
 ("h_15", 5.0e-5, 1.0e-5, 1.0e-4, "v4 no transitantes: contraste del halo XAO a 15 mas (Kasper+2021 ~5e-5)"),
 ("h_ratio", 10.0, 3.0, 30.0, "v4 no transitantes: mejora del halo de 15 a 100 mas"),
 ("iwa_lo", 2.0, 1.5, 3.0, "v5 no transitantes: IWA minimo en lambda/D; IWA ~ U(iwa_lo, iwa_hi)"),
 ("iwa_hi", 4.0, 3.0, 6.0, "v4 no transitantes: IWA maximo en lambda/D (lambda=0.75 um); IWA ~ U(2, iwa_hi)"),
 ("ag_fac", 2.0/3.0, 0.5, 1.2, "v4 no transitantes: A_g / A_Bond (2/3 = Lambert)"),
 ("tsm_small", 10.0, 5.0, 20.0, "Tiers: umbral TSM para Rp<1.5 (Kempton recomienda 10)"),
 ("snr_thr", 5.0, 3.0, 10.0, "v3: umbral de deteccion SNR en P=Phi(SNR-umbral)"),
 ("sig_model", 0.5, 0.2, 1.0, "v3: incertidumbre (ln) de la SNR predicha"),
 ("cloud_min", 0.5, 0.1, 1.0, "v3: supresion por nubes en transmision ~ U(cloud_min,1)"),
 ("n_ecl", 2.0, 1.0, 10.0, "v3: eclipses secundarios para S_atm (ESM = S/N por eclipse)"),
 ("floor_ppm", 20.0, 10.0, 50.0, "v3: piso de ruido sistematico por canal [ppm] (Greene+2016: 20 NIRISS)"),
 ("n_H", 2.0, 1.0, 5.0, "v3: escalas de altura que abarcan los rasgos espectrales"),
 ("n_bins", 10.0, 1.0, 30.0, "v3: canales independientes sobre los que se promedia el piso"),
 ("shore_logc", 0.0, -0.5, 0.5, "v3: log10 factor de la linea de costa cosmica (0 = pasa por Marte)"),
 ("shore_sig", 1.0, 0.5, 2.0, "v3: ancho (ln) de la transicion de la linea de costa"),
 ("ntr_cap", 10.0, 5.0, 30.0, "S_obs transitantes: maximo de transitos disponibles"),
 ("avail_T2", 0.5, 0.3, 0.7, "S_obs: Availability del tier T2"),
 ("avail_T3", 0.25, 0.1, 0.5, "S_obs: Availability del tier T3 (T0=1, T1=0.75 fijos)"),
 ("tier_x", 1.0, 0.5, 2.0, "Tiers: factor sobre los umbrales TSM/ESM (T0>=1x, T1>=0.5x, T2>=0.1x)"),
 ("flare_late", 0.5, 0.2, 0.8, "S_estrella: penalizacion por flares, M*<0.35"),
 ("flare_early", 0.25, 0.0, 0.5, "S_estrella: penalizacion por flares, 0.35<=M*<0.6"),
 ("xuv_mid", 10.0, 3.0, 30.0, "S_estrella: F_XUV (x Tierra) con 50% de estres"),
 ("tsat_M", 1.0, 0.5, 3.0, "S_estrella: tiempo de saturacion XUV para M*<0.6 [Gyr]"),
 ("q_lo", 0.03, 0.015, 0.06, "S_geo: inicio del plateau de flujo interno [W/m2]"),
 ("q_hi", 1.0, 0.3, 3.0, "S_geo: fin del plateau de flujo interno [W/m2]"),
 ("log_env", -2.0, -3.0, -1.3, "S_geo: log10 fraccion de envolvente para erosion (sub-Neptunos)"),
 ("mag_c", 2.0, 1.0, 4.0, "S_mag: concentracion del Beta base (alfa=beta=c)"),
 ("mag_lock_b", 3.0, 2.0, 5.0, "S_mag: beta del Beta si rotacion sincronica probable (alfa=2)"),
 ("age_hi", 10.0, 6.0, 13.0, "Edad desconocida ~ U(1, age_hi) [Gyr]"),
 ("R0", 1.65, 1.50, 1.80, "v5 mu(Rp): centro de la log-sigmoide [R_E] (Fulton+2017: deficit en 1.5-2.0)"),
 ("dR", 0.12, 0.06, 0.24, "v5 mu(Rp): ancho de la transicion [R_E]"),
 ("yang_eff_lo", 0.5, 0.0, 1.0, "v5 clima M: eficacia minima del albedo de nubes subestelares de Yang+2013 ~ U(lo,1)"),
 ("kin_lo", 1.5, 1.0, 2.0, "v5 clima M: factor minimo del borde interno ~ U(kin_lo, 2) (Yang+2013 x2; Kopparapu+2016 -25%)"),
 ("lpco2_hi", 2.0, 0.0, 2.0, "v5 clima M: log10 de la presion maxima de CO2 del prior [bar] (lo = -3)"),
 ("f_eye_c", 0.30, 0.15, 0.50, "v5 clima M: centro de la fraccion habitable del piso eyeball/CO2 ~ U(c-0.15, c+0.15)"),
 ("s_eye_max", 0.8, 0.7, 0.9, "v5 clima M: S_eff maximo donde aplica el piso de borde externo"),
]
NAMES = [p[0] for p in PARAMS]
BASE = {p[0]: p[1] for p in PARAMS}
BOUNDS = [[p[2], p[3]] for p in PARAMS]
KEY = ["TOI-7166 b", "GJ 3378 b", "GJ 273 b", "Proxima Cen b", "Proxima Cen d", "TRAPPIST-1 d", "TRAPPIST-1 e",
       "TRAPPIST-1 f", "TOI-700 d", "TOI-700 e", "K2-18 b", "LHS 1140 b", "Kepler-186 f", "TOI-4336 A b"]


def lg(x, x0, k):
    return 1.0 / (1.0 + np.exp(-k * (np.log10(np.maximum(x, 1e-30)) - np.log10(x0))))


def trap(x, a, b, c, d):
    y = np.zeros_like(x)
    y = np.where((x > a) & (x < b), (x - a) / (b - a), y)
    y = np.where((x >= b) & (x <= c), 1.0, y)
    y = np.where((x > c) & (x < d), (d - x) / (d - c), y)
    return y


class Model:
    def __init__(self, df, planets, n=300, seed=7):
        d = U.prepare(df)
        d = d[d.pl_name.isin(planets)].reset_index(drop=True)
        self.names = d.pl_name.values
        rng = np.random.default_rng(seed)
        P = len(d); self.P, self.n = P, n
        A = {k: np.zeros((P, n)) for k in ["teff", "rs", "ms", "a", "per", "e", "rp", "mp", "age_k",
                                            "u_age", "u_alb", "u_A", "u_mu", "u_mag", "u_cl", "z1", "z2", "z3", "u_iwa", "u_eff", "u_kin", "u_pco2", "u_eye", "z_R0"]}
        W = rng.dirichlet(U.DIRICHLET_ALPHA, (P, n))          # pesos del documento, fijos
        for i, r in enumerate(d.itertuples()):
            teff = U.split_normal(rng, r.st_teff, r.st_tefferr1, r.st_tefferr2, n, 0.02, 2300, 50000)
            rs = U.split_normal(rng, r.st_rad, r.st_raderr1, r.st_raderr2, n, 0.07, 0.05, None)
            ms = U.split_normal(rng, r.st_mass, r.st_masserr1, r.st_masserr2, n, 0.07, 0.05, None)
            if np.isfinite(r.st_age):
                A["age_k"][i] = U.split_normal(rng, r.st_age, r.st_ageerr1, r.st_ageerr2, n, 0.3, 0.01, 13.8)
            else:
                A["age_k"][i] = np.nan
            if np.isfinite(r.pl_orbper):
                per = U.split_normal(rng, r.pl_orbper, r.pl_orbpererr1, r.pl_orbpererr2, n, 1e-4, 1e-3, None)
                a = ((U.G * ms * U.MSUN) / (4 * np.pi ** 2)) ** (1 / 3) * (per * U.DAY) ** (2 / 3) / U.AU
            else:
                a = U.split_normal(rng, r.pl_orbsmax, r.pl_orbsmaxerr1, r.pl_orbsmaxerr2, n, 0.05, 1e-4, None)
                per = np.sqrt(a ** 3 / ms) * 365.25
            if np.isfinite(r.pl_orbeccen):
                e = U.split_normal(rng, r.pl_orbeccen, r.pl_orbeccenerr1, r.pl_orbeccenerr2, n, 0.0, 0.0, 0.95)
                if not (np.isfinite(r.pl_orbeccenerr1) or np.isfinite(r.pl_orbeccenerr2)):
                    e = np.clip(r.pl_orbeccen + 0.05 * rng.standard_normal(n), 0, 0.95)
            else:
                e = np.zeros(n)
            rp = U.split_normal(rng, r.pl_rade, r.pl_radeerr1, r.pl_radeerr2, n, 0.10, 0.1, None) if r.radius_measured else None
            if r.mass_mode == "measured":
                mp = U.split_normal(rng, r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, n, 0.20, 0.01, None)
            elif r.mass_mode == "msini":
                msini = U.split_normal(rng, r.pl_bmasse, r.pl_bmasseerr1, r.pl_bmasseerr2, n, 0.20, 0.01, None)
                cosi = rng.uniform(0, 1, n)
                mp = np.minimum(msini / np.sqrt(np.clip(1 - cosi ** 2, 1e-4, 1)), 26600.0)
            else:
                mp = None
            if rp is None and mp is None:
                rp = U.split_normal(rng, r.pl_rade, r.pl_radeerr1, r.pl_radeerr2, n, 0.15, 0.1, None)
            if mp is None: mp = U.ck_mass(rp, rng)
            if rp is None: rp = U.ck_radius(mp, rng)
            for k, v in dict(teff=teff, rs=rs, ms=ms, a=a, per=per, e=e, rp=rp, mp=mp).items():
                A[k][i] = v
            for k in ["u_age", "u_alb", "u_A", "u_mu", "u_mag", "u_cl"]:
                A[k][i] = rng.uniform(size=n)
            for k in ["z1", "z2", "z3"]:
                A[k][i] = rng.standard_normal(n)
            A["u_iwa"][i] = np.random.default_rng(zlib.crc32(str(r.pl_name).encode())).uniform(size=n)
            r5 = np.random.default_rng(zlib.crc32(str(r.pl_name).encode()) + 1)    # igual que uhi_code v5
            u5 = r5.uniform(size=(4, n)); A["z_R0"][i] = r5.standard_normal(n)
            for j, k in enumerate(["u_eff", "u_kin", "u_pco2", "u_eye"]):
                A[k][i] = u5[j]
        self.A, self.W = A, W
        # cantidades que no dependen de los supuestos
        teff, rs, ms, a, rp, mp = A["teff"], A["rs"], A["ms"], A["a"], A["rp"], A["mp"]
        self.L = rs ** 2 * (teff / 5772.0) ** 4
        self.S = self.L / a ** 2
        k = U.kopparapu_seff(teff)
        self.cons = (self.S <= k["rg"]) & (self.S >= k["mg"])
        self.opt = (self.S <= k["rv"]) & (self.S >= k["em"])
        self.delta = (rp * U.REARTH / (rs * U.RSUN)) ** 2
        self.g = U.G * mp * U.MEARTH / (rp * U.REARTH) ** 2
        teq0 = teff * np.sqrt(rs * U.RSUN / (a * U.AU)) * 0.25 ** 0.25
        J = d.J_use.values[:, None]; K = d.K_use.values[:, None]
        self.J = J; self.Kmag = K
        self.tran = (d.tran_flag.values == 1)[:, None]
        self.tsm_raw = U.tsm_scale(rp) * rp ** 3 * teq0 / (np.maximum(mp, 0.05) * rs ** 2) * 10 ** (-J / 5)
        self.esm = 4.29e6 * U.planck_ratio(1.10 * teq0, teff) * self.delta * 10 ** (-K / 5)
        self.dist = d.sy_dist.values[:, None].astype(float)
        self.small = rp < 1.6
        n_orb = 2 * np.pi / (A["per"] * U.DAY)
        Edot = 10.5 * (0.3 / 100) * U.G * (ms * U.MSUN) ** 2 * (rp * U.REARTH) ** 5 * n_orb * A["e"] ** 2 / (a * U.AU) ** 6
        self.q_tide = Edot / (4 * np.pi * (rp * U.REARTH) ** 2)
        self.subnep = (rp >= 1.6) & (rp <= 4.0)
        self.mlow = mp < 0.3
        st_mass = d.st_mass.values[:, None]
        self.st_mass = st_mass
        self.locked = ((d.st_mass.values < 0.6) & (np.median(a, 1) < 0.2))[:, None]
        self.dens_ok = d.dens_ok.values.astype(bool)[:, None]
        # tiers: ingredientes nominales
        a_nom = np.where(d.pl_orbper.notna(), ((U.G * d.st_mass * U.MSUN) / (4 * np.pi ** 2)) ** (1 / 3) * (d.pl_orbper * U.DAY) ** (2 / 3) / U.AU, d.pl_orbsmax).astype(float)
        rp_nom = np.nanmedian(rp, 1); mp_nom = np.nanmedian(mp, 1)
        rp_nom = np.where(d.pl_rade.notna(), d.pl_rade, rp_nom); mp_nom = np.where(d.pl_bmasse.notna(), d.pl_bmasse, mp_nom)
        teq0n = d.st_teff.values * np.sqrt(d.st_rad.values * U.RSUN / (a_nom * U.AU)) * 0.25 ** 0.25
        self.tsm_n = U.tsm_scale(rp_nom) * rp_nom ** 3 * teq0n / (np.maximum(mp_nom, 0.05) * d.st_rad.values ** 2) * 10 ** (-d.J_use.values / 5)
        self.esm_n = 4.29e6 * U.planck_ratio(1.10 * teq0n, d.st_teff.values) * (rp_nom * U.REARTH / (d.st_rad.values * U.RSUN)) ** 2 * 10 ** (-d.K_use.values / 5)
        self.rp_nom = rp_nom
        e1 = np.abs(d.pl_bmasseerr1.fillna(np.inf).values); e2 = np.abs(d.pl_bmasseerr2.fillna(np.inf).values)
        self.mrel_ok = (d.mass_mode.values == "measured") & (0.5 * (e1 + e2) / d.pl_bmasse.values <= 0.5)
        self.tran_n = d.tran_flag.values == 1
        self.dist_n = d.sy_dist.values
        self.a_nom = a_nom; self.J_n = d.J_use.values

    def tiers(self, p):
        thr = np.where(self.rp_nom < 1.5, p["tsm_small"], 90.0)
        x = np.maximum(self.tsm_n / thr, self.esm_n / 7.5) / p["tier_x"]
        t = np.full(self.P, 3)
        tr = self.tran_n
        t = np.where(tr & (self.rp_nom >= 1.6) & (x >= 0.1), 2, t)
        t = np.where(tr & (x >= 0.5), 1, t)
        t = np.where(tr & self.mrel_ok & (x >= 1), 0, t)
        PP = dict(U.V3); PP.update({k: p[k] for k in U.V3 if k in p})
        dok = np.isfinite(self.dist_n) & (self.dist_n > 0)
        snr_n, _, _ = U.snr_reflected(self.rp_nom, self.a_nom, np.where(dok, self.dist_n, 1.0), self.J_n, 0.3, 0.0, PP)
        nt = ~tr & dok & (snr_n >= 0.1 * p["snr_thr"])
        t = np.where(nt, 2, t)
        return t

    def evaluate(self, p):
        A = self.A
        S = self.S
        # S_fis
        f_r = trap(A["rp"], p["fR_a"], p["fR_b"], p["fR_c"], p["fR_d"])
        f_m = trap(A["mp"], 0.3 * p["fM_b"], p["fM_b"], p["fM_c"], p["fM_d"])
        S_fis = np.sqrt(f_r * f_m)
        # clima
        PP = dict(U.V3); PP.update({k: p[k] for k in U.V3 if k in p})
        alb_sh = U.albedo_teff(A["teff"], A["u_alb"], p["f_land_max"], p["alb_shift"], p["late_m_scale"])
        Aolr = p["olr_Amin"] + A["u_A"] * (203.3 - p["olr_Amin"])
        albedo, fhz, p_hz, _ = U.climate_v5(S, A["teff"], self.locked, alb_sh, Aolr, A["u_eff"], A["u_kin"], A["u_pco2"],
                                            A["u_eye"], PP, p["ebm_D"], (1.0, p["pHZ_opt"], p["pHZ_out"]))
        S_clima = fhz * p_hz
        # atm + obs (v3, funcion compartida con uhi_code)
        tier = self.tiers(p)
        avail = np.array([1.0, 0.75, p["avail_T2"], p["avail_T3"]])[tier][:, None]
        teq = ((1 - albedo) * self.L * U.LSUN / (16 * np.pi * U.SIGMA * (A["a"] * U.AU) ** 2)) ** 0.25
        S_atm, S_obs, _, _, _, _ = U.atm_obs_v5(self.tran, A["rp"], A["mp"], A["rs"], A["teff"], A["a"], A["per"], teq, self.g,
                                              self.J, self.Kmag, S, avail, self.dist, albedo, A["u_iwa"],
                                              A["z1"], A["z2"], A["z3"], A["u_cl"], A["u_mu"], A["z_R0"], self.dens_ok, PP)
        # estrella
        age = np.where(np.isnan(A["age_k"]), 1.0 + A["u_age"] * (p["age_hi"] - 1.0), A["age_k"])
        tsat = np.where(self.st_mass >= 0.6, 0.1, p["tsat_M"])
        rx = np.where(age < tsat, 10 ** -3.5, 10 ** -3.5 * (age / tsat) ** -1.23)
        rx_sun = 10 ** -3.5 * (4.6 / 0.1) ** -1.23
        fxuv = rx * S / rx_sun
        stress = lg(fxuv, p["xuv_mid"], 4.0)
        flare = np.where(self.st_mass < 0.35, p["flare_late"], np.where(self.st_mass < 0.6, p["flare_early"], 0.0))
        S_est = (1 - stress) * (1 - flare)
        # geo
        q_rad = 0.087 * (A["mp"] / A["rp"] ** 2) * 2.0 ** (-(age - 4.5) / 3.0)
        lq = np.log10(np.maximum(q_rad + self.q_tide, 1e-6))
        S_reg = trap(lq, -2.0, np.log10(p["q_lo"]), np.log10(p["q_hi"]), 1.0)
        mdot = 0.1 * np.pi * (A["rp"] * U.REARTH) ** 3 * fxuv * U.F_XUV_EARTH / (U.G * A["mp"] * U.MEARTH)
        t_loss = 10 ** p["log_env"] * A["mp"] * U.MEARTH / np.maximum(mdot, 1e-30) / (1e9 * U.YR)
        risk = 1 - lg(t_loss, 1.0, 4.0)
        S_geo = S_reg * np.where(self.subnep, 1 - 0.5 * risk, 1.0)
        # mag
        aB = np.where(self.locked, 2.0, p["mag_c"])
        bB = np.where(self.locked, p["mag_lock_b"], p["mag_c"]) + np.where(self.mlow, 1.0, 0.0)
        aB = np.broadcast_to(aB, S.shape); bB = np.broadcast_to(bB, S.shape)
        S_mag = betaincinv(aB, bB, A["u_mag"])
        Sm = np.stack([S_fis, S_clima, S_atm, S_est, S_geo, S_mag, S_obs], -1)
        alive = np.all(Sm > U.EPS, -1)
        uhi = np.where(alive, np.exp((self.W * np.log(np.maximum(np.minimum(Sm, 1.0), U.EPS))).sum(-1)), 0.0)
        p_alive = alive.mean(1)
        with np.errstate(all="ignore"):
            cmed = np.nanmedian(np.where(alive, uhi, np.nan), 1)
        score = p_alive * np.nan_to_num(cmed)               # v5: UHI_rank = P_alive x mediana condicional
        mean = uhi.mean(1)
        return score, mean

    @staticmethod
    def score(med, mean):
        # orden: UHI_rank, empates por media (igual que uhi_code.run v5)
        return med + 1e-9 * mean
