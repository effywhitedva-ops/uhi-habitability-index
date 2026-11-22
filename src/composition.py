"""UHI v6 - clase de composicion por draw: R (terrestre, atmosfera secundaria), H (Hycean), G (mini-Neptuno rico en gas).
Fuentes: Zeng+2019 (tablas M-R descargadas de lweb.cfa.harvard.edu/~lzeng, md5 en assumptions_v6.md);
Otegi+2020 (M = 0.9 R^3.45 rocosos; M = 1.74 R^1.58 volatiles, con errores de coeficientes); Venturini+2024
(R_valley = 1.8197 M*^0.14 (+0.02/-0.01)); Madhusudhan+2021 (R_H,max = 2.3 R_E a 5 M_E, 2.6 R_E a 10 M_E)."""
import os
import numpy as np

_D = os.path.join(os.path.dirname(__file__), "..", "data", "zeng2019")
_ROCK = np.loadtxt(os.path.join(_D, "massradiusEarthlikeRocky.txt"))
_H2O = {T: np.loadtxt(os.path.join(_D, f"massradius_50percentH2O_{T}K_1mbar.txt")) for T in (300, 500, 700, 1000)}   # 50% agua
_H2O100 = {T: np.loadtxt(os.path.join(_D, f"massradius_100percentH2O_{T}K_1mbar.txt")) for T in (300, 500, 700, 1000)}
_SIL = np.loadtxt(os.path.join(_D, "massradiusmgsio3.txt"))                                                              # roca sin hierro
_TS = np.array([300.0, 500.0, 700.0, 1000.0])


def r_rock(m):
    return np.interp(np.log(m), np.log(_ROCK[:, 0]), _ROCK[:, 1])


def r_silicate(m):
    """Radio maximo de un cuerpo rocoso sin volatiles (MgSiO3 puro, Zeng+2019)."""
    return np.interp(np.log(m), np.log(_SIL[:, 0]), _SIL[:, 1])


def r_h2o(m, teq, tables=None):
    """Linea de 50% agua + 50% nucleo terrestre de Zeng+2019 (1 mbar; mundos de agua), interpolada linealmente en T entre 300-1000 K (constante fuera; SUPUESTO)."""
    tb = _H2O if tables is None else tables
    rs = np.vstack([np.interp(np.log(m), np.log(tb[T][:, 0]), tb[T][:, 1]) for T in (300, 500, 700, 1000)])
    t = np.clip(teq, 300.0, 1000.0); k = np.clip(np.searchsorted(_TS, t) - 1, 0, 2); w = (t - _TS[k]) / (_TS[k + 1] - _TS[k])
    cols = np.arange(rs.shape[1]) if rs.ndim == 2 else 0
    return rs[k, cols] * (1 - w) + rs[k + 1, cols] * w


def r_hycean_max(m):
    """Madhusudhan+2021: 2.3 R_E (5 M_E), 2.6 R_E (10 M_E); lineal en log M, extrapolado (SUPUESTO)."""
    return 2.3 + (2.6 - 2.3) * np.log(np.maximum(m, 1e-3) / 5.0) / np.log(2.0)


def otegi_mass(rp, rocky, z1, z2):
    """Masa por clase desde el radio (Otegi+2020), muestreando los coeficientes con sus errores (z1, z2 ~ N(0,1))."""
    c_r, e_r = 0.9 + 0.06 * z1, 3.45 + 0.12 * z2
    c_v, e_v = 1.74 + 0.38 * z1, 1.58 + 0.10 * z2
    return np.where(rocky, np.maximum(c_r, 0.5) * rp ** e_r, np.maximum(c_v, 0.5) * rp ** e_v)


def otegi_radius(mp, rocky, z1, z2):
    c_r, e_r = 0.9 + 0.06 * z1, 3.45 + 0.12 * z2
    c_v, e_v = 1.74 + 0.38 * z1, 1.58 + 0.10 * z2
    return np.where(rocky, (mp / np.maximum(c_r, 0.5)) ** (1 / e_r), (mp / np.maximum(c_v, 0.5)) ** (1 / e_v))


def r_valley(mstar, z):
    """Venturini+2024 con exponente normal partido (+0.02/-0.01)."""
    ex = 0.14 + np.where(z >= 0, 0.02, 0.01) * z
    return 1.8197 * np.maximum(mstar, 0.08) ** ex


def classify_density(mp, rp, teq, u, pi=(1.0, 1.0, 1.0), rock_line="h2o50"):
    """Ruta de densidad: hipotesis compatibles con (M, R) del draw, sorteo con priors pi renormalizados.
    R (sin envoltura de H2: rocoso o mundo de agua hasta 50% H2O): Rp <= R_50%H2O(M, Teq);
    H: R_MgSiO3(M) < Rp <= R_Hmax(M);  G: Rp > R_MgSiO3(M)  (por debajo de la linea de silicato puro no hace falta ningun volatil).
    Si ninguna es compatible -> R (planeta mas denso que todas las lineas)."""
    # variantes de frontera (prueba de sensibilidad (d)):
    #   h2o50 (principal): R <= linea 50% H2O, volatiles > MgSiO3;  h2o100: R <= 100% H2O, volatiles > MgSiO3;
    #   earthlike: frontera unica en la linea Earth-like (R <= R_Earthlike < H/G), sin solape (estricta, "rocoso seco")
    if rock_line == "earthlike":
        rr = r_rock(mp); rw = rr.copy()
    elif rock_line == "h2o100":
        rr = r_silicate(mp); rw = r_h2o(mp, teq, _H2O100)
    else:
        rr = r_silicate(mp); rw = r_h2o(mp, teq)
    rh = r_hycean_max(mp)
    cR = rp <= rw; cH = (rp > rr) & (rp <= rh); cG = rp > rr
    wR, wH, wG = pi[0] * cR, pi[1] * cH, pi[2] * cG
    tot = wR + wH + wG
    tot = np.where(tot > 0, tot, 1.0); wR = np.where((wR + wH + wG) > 0, wR, 1.0)   # sin hipotesis -> R (denso)
    x = u * tot
    return np.where(x < wR, 0, np.where(x < wR + wH, 1, 2))      # 0=R, 1=H, 2=G


def prior_vol(p_vol, pi):
    """v6.1 (fix 3): aplica el prior de clase a una probabilidad de 'volatil' (H o G) calculada con prior plano R:V = 1:1.
    pi_V = (pi_H + pi_G)/2 (SUPUESTO: asi 1:1:1 deja todo identico a v6.0, decision D4 aprobada):
    P'(V) = pi_V p / (pi_R (1 - p) + pi_V p). Si pi_R == pi_V devuelve p sin tocar (identidad exacta)."""
    piV = 0.5 * (pi[1] + pi[2])
    if pi[0] == piV:
        return p_vol
    return piV * p_vol / (pi[0] * (1.0 - p_vol) + piV * p_vol)


def classify_radius(rp, mstar, u_h2, u_hg, z_v, w, mp_vol, pi=(1.0, 1.0, 1.0)):
    """Ruta de radio: P(H2) = sigmoide((Rp - R_v(M*))/w); dentro de H2, H si Rp <= R_Hmax(M_vol) con pi_H/(pi_H+pi_G)."""
    ph2 = prior_vol(1.0 / (1.0 + np.exp(-(rp - r_valley(mstar, z_v)) / w)), pi)   # v6.1: con pi_R
    h2 = u_h2 < ph2
    canH = rp <= r_hycean_max(mp_vol)
    isH = canH & (u_hg < pi[1] / (pi[1] + pi[2]))
    return np.where(~h2, 0, np.where(isH, 1, 2)), ph2


def classify_mass(mp, u, u_hg, mv, s, rp_vol, pi=(1.0, 1.0, 1.0)):
    """No transitantes: P_R(M) = sigmoide((ln M_v - ln M)/s) (SUPUESTO, M_v ~ U(3,6) M_E, s = 0.3)."""
    pr = 1.0 / (1.0 + np.exp(-(np.log(mv) - np.log(np.maximum(mp, 1e-3))) / s))
    piV = 0.5 * (pi[1] + pi[2])
    if pi[0] != piV:                                  # v6.1 (fix 3): pi_R tambien en la ruta de masa
        pr = pi[0] * pr / (pi[0] * pr + piV * (1.0 - pr))
    rocky = u < pr
    isH = (rp_vol <= r_hycean_max(mp)) & (u_hg < pi[1] / (pi[1] + pi[2]))
    return np.where(rocky, 0, np.where(isH, 1, 2)), pr
