r"""Per-draw composition classes for UHI v6.1.

Classes are encoded as 0 = rocky or water-rich world without an H2 envelope (R),
1 = Hycean, and 2 = gas-rich mini-Neptune (G). Inputs are planetary mass/radius in
Earth units, equilibrium temperature in K, stellar mass in solar masses and independent
U(0,1) or N(0,1) draws. Functions return radii in R_E, masses in M_E, class codes, or
dimensionless probabilities.

Sources: Zeng et al. (2019, PNAS 116, 9723) mass-radius tables (downloaded copies and
md5 checksums in docs/assumptions_v6_es.md); Otegi et al. (2020, A&A 634, A43)
rocky/volatile M-R relations with coefficient errors; Venturini et al. (2024, A&A 686,
L9) stellar-mass-dependent radius valley; Madhusudhan et al. (2021, ApJ 918, 1)
maximum Hycean radii. Class priors are equal by default (assumption D4): when data
admit several hypotheses, the model represents ambiguity instead of choosing the most
habitable answer. Linear interpolations/extrapolations are documented assumptions."""
import os
import numpy as np

_D = os.path.join(os.path.dirname(__file__), "..", "data", "zeng2019")
_ROCK = np.loadtxt(os.path.join(_D, "massradiusEarthlikeRocky.txt"))
_H2O = {T: np.loadtxt(os.path.join(_D, f"massradius_50percentH2O_{T}K_1mbar.txt")) for T in (300, 500, 700, 1000)}   # 50% water
_H2O100 = {T: np.loadtxt(os.path.join(_D, f"massradius_100percentH2O_{T}K_1mbar.txt")) for T in (300, 500, 700, 1000)}
_SIL = np.loadtxt(os.path.join(_D, "massradiusmgsio3.txt"))                                                              # iron-free rock
_TS = np.array([300.0, 500.0, 700.0, 1000.0])


def r_rock(m):
    r"""Return the Earth-like rocky radius in R_E for mass m in M_E (Zeng et al. 2019).

    Interpolation is linear in log mass. Values outside the table use NumPy edge clamping."""
    return np.interp(np.log(m), np.log(_ROCK[:, 0]), _ROCK[:, 1])


def r_silicate(m):
    r"""Return the pure-MgSiO3 radius in R_E for mass m in M_E (Zeng et al. 2019).

    This iron-free silicate curve marks the largest radius that does not require volatiles; it is
    used as the lower bound for volatile-rich classes. Interpolation is linear in log mass."""
    return np.interp(np.log(m), np.log(_SIL[:, 0]), _SIL[:, 1])


def r_h2o(m, teq, tables=None):
    r"""Return the radius in R_E of a Zeng et al. (2019) 50% water / 50% Earth-like world.

    m is M_E and teq is equilibrium temperature in K; optional tables selects another water line
    (the 100% sensitivity case). Tables at 300/500/700/1000 K and 1 mbar are interpolated
    linearly in temperature, with temperatures clipped to 300-1000 K (assumption)."""
    tb = _H2O if tables is None else tables
    rs = np.vstack([np.interp(np.log(m), np.log(tb[T][:, 0]), tb[T][:, 1]) for T in (300, 500, 700, 1000)])
    t = np.clip(teq, 300.0, 1000.0); k = np.clip(np.searchsorted(_TS, t) - 1, 0, 2); w = (t - _TS[k]) / (_TS[k + 1] - _TS[k])
    cols = np.arange(rs.shape[1]) if rs.ndim == 2 else 0
    return rs[k, cols] * (1 - w) + rs[k + 1, cols] * w


def r_hycean_max(m):
    r"""Return the maximum Hycean radius in R_E for mass m in M_E.

    The 2.3 R_E at 5 M_E and 2.6 R_E at 10 M_E anchors follow Madhusudhan et al. (2021).
    Linear interpolation in log mass and extrapolation beyond those anchors are assumptions."""
    return 2.3 + (2.6 - 2.3) * np.log(np.maximum(m, 1e-3) / 5.0) / np.log(2.0)


def otegi_mass(rp, rocky, z1, z2):
    r"""Return class-dependent planetary mass in M_E from radius rp in R_E.

    rocky is a Boolean array; z1 and z2 are standard normals applied to Otegi et al. (2020)
    coefficient and exponent errors, using the same draws for both populations. The 0.5
    coefficient floor prevents nonphysical negative or near-zero normalizations from tails."""
    c_r, e_r = 0.9 + 0.06 * z1, 3.45 + 0.12 * z2
    c_v, e_v = 1.74 + 0.38 * z1, 1.58 + 0.10 * z2
    return np.where(rocky, np.maximum(c_r, 0.5) * rp ** e_r, np.maximum(c_v, 0.5) * rp ** e_v)


def otegi_radius(mp, rocky, z1, z2):
    r"""Invert the Otegi et al. (2020) relations: return radius in R_E from mass in M_E.

    rocky and z1/z2 have the same meaning as in otegi_mass. It is used for non-transiting draws
    where mass is the available constraint. The coefficient floor matches otegi_mass."""
    c_r, e_r = 0.9 + 0.06 * z1, 3.45 + 0.12 * z2
    c_v, e_v = 1.74 + 0.38 * z1, 1.58 + 0.10 * z2
    return np.where(rocky, (mp / np.maximum(c_r, 0.5)) ** (1 / e_r), (mp / np.maximum(c_v, 0.5)) ** (1 / e_v))


def r_valley(mstar, z):
    r"""Return the radius-valley location in R_E for stellar mass mstar in M_sun.

    R_valley = 1.8197 M_star^(0.14) follows Venturini et al. (2024); z is a standard normal for
    the asymmetric exponent error (+0.02/-0.01). The 0.08-M_sun lower clamp is numerical."""
    ex = 0.14 + np.where(z >= 0, 0.02, 0.01) * z
    return 1.8197 * np.maximum(mstar, 0.08) ** ex


def classify_density(mp, rp, teq, u, pi=(1.0, 1.0, 1.0), rock_line="h2o50"):
    r"""Assign classes from paired mass/radius draws when density is informative.

    mp in M_E, rp in R_E, teq in K, u U(0,1), pi three nonnegative class weights, and rock_line
    selects the boundary test. A draw is compatible with R at or below the 50% water line, H
    between pure silicate and the Madhusudhan maximum, and G above pure silicate. Weights are
    renormalized over compatible hypotheses. A draw denser than every line is assigned R.
    rock_line variants are prior test (d): h2o50 production, h2o100, or strict earthlike."""
    # boundary variants (sensitivity test (d)):
    #   h2o50 (production): R <= 50% H2O line, volatiles > MgSiO3;  h2o100: R <= 100% H2O, volatiles > MgSiO3;
    #   earthlike: one boundary at the Earth-like line (R <= R_Earthlike < H/G), no overlap (strict 'dry rocky')
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
    tot = np.where(tot > 0, tot, 1.0); wR = np.where((wR + wH + wG) > 0, wR, 1.0)   # no compatible hypothesis -> R (dense)
    x = u * tot
    return np.where(x < wR, 0, np.where(x < wR + wH, 1, 2))      # 0=R, 1=H, 2=G


def prior_vol(p_vol, pi):
    r"""Apply class priors to a probability of having volatiles.

    p_vol is computed under an implicit 1:1 rocky:volatile prior. pi_V is the average of
    Hycean and gas-rich weights so equal (1,1,1) priors reproduce v6.0 exactly. v6.1 fix 3
    uses Bayes' rule P'(V)=pi_V p/[pi_R(1-p)+pi_V p]. Returns a dimensionless probability."""
    piV = 0.5 * (pi[1] + pi[2])
    if pi[0] == piV:
        return p_vol
    return piV * p_vol / (pi[0] * (1.0 - p_vol) + piV * p_vol)


def classify_radius(rp, mstar, u_h2, u_hg, z_v, w, mp_vol, pi=(1.0, 1.0, 1.0)):
    r"""Assign classes from radius when no informative density measurement is available.

    rp in R_E, mstar in M_sun, u_h2/u_hg U(0,1), z_v N(0,1), w valley width in R_E, mp_vol the
    volatile-route mass in M_E, and pi class weights. P(H2) is a logistic function of distance
    from the valley, adjusted by prior_vol. H2 draws are Hycean only if compatible with the
    Hycean maximum radius and selected using pi_H/(pi_H+pi_G). Return (classes, P(H2))."""
    ph2 = prior_vol(1.0 / (1.0 + np.exp(-(rp - r_valley(mstar, z_v)) / w)), pi)   # v6.1: includes pi_R
    h2 = u_h2 < ph2
    canH = rp <= r_hycean_max(mp_vol)
    isH = canH & (u_hg < pi[1] / (pi[1] + pi[2]))
    return np.where(~h2, 0, np.where(isH, 1, 2)), ph2


def classify_mass(mp, u, u_hg, mv, s, rp_vol, pi=(1.0, 1.0, 1.0)):
    r"""Assign classes from mass alone, mainly for non-transiting planets.

    mp in M_E, u/u_hg U(0,1), mv transition mass in M_E, s logistic width in natural-log units,
    rp_vol volatile radius in R_E, and pi class weights. The logistic rocky probability with
    M_v ~ U(3,6) M_E and s=0.3 is an assumption (D4); v6.1 fix 3 applies pi_R in this route.
    Return (classes, P(rocky))."""
    pr = 1.0 / (1.0 + np.exp(-(np.log(mv) - np.log(np.maximum(mp, 1e-3))) / s))
    piV = 0.5 * (pi[1] + pi[2])
    if pi[0] != piV:                                  # v6.1 (fix 3): pi_R also applied in the mass route
        pr = pi[0] * pr / (pi[0] * pr + piV * (1.0 - pr))
    rocky = u < pr
    isH = (rp_vol <= r_hycean_max(mp)) & (u_hg < pi[1] / (pi[1] + pi[2]))
    return np.where(rocky, 0, np.where(isH, 1, 2)), pr
