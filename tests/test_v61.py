"""Tests unitarios de los fixes v6.1 (2026-11-20).  python3 tests/test_v61.py
 fix 2: piso fisico de R*; fix 3: pi_R en rutas de radio y masa; fix 4: Innes 0.0860; fix 5: reporte Sobol (Ishigami + degeneradas);
 regresion: con 1:1:1 las rutas de clase dan resultados identicos bit a bit a v6.0 (copia en backups/)."""
import os, sys, importlib.util, warnings; warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np, pandas as pd
import uhi_code as u, uhi_code_v6 as v6, composition as C, xuv_history as X, sobol_report_v6 as SR
BK = os.path.join(ROOT, "backups", "pre_v6.1_20261120_1625", "src")
def load_old(name):
    import types; m = types.ModuleType(name + "_v60"); m.__file__ = os.path.join(ROOT, "src", name + ".py")   # rutas de datos relativas a src/
    exec(compile(open(os.path.join(BK, name + ".py")).read(), name + "_v60", "exec"), m.__dict__); return m
OK = []
def check(name, cond, info=""):
    assert cond, f"FALLA {name}: {info}"; OK.append(name); print("ok", name, info)

# ---- fix 4 ----
check("fix4_innes_M_10bar", X.INNES_S10[0] == 0.0860 and X.INNES_S10[1] == 0.0676 and X.INNES_S1.tolist() == [0.320, 0.389])
check("fix4_s_crit_10bar_M", abs(X.s_crit_innes(3000.0, 10.0) - 0.0860) < 1e-12, X.s_crit_innes(3000.0, 10.0))

# ---- fix 2 ----
check("fix2_r_ms_sun", abs(float(v6.r_ms(1.0)) - 1.06) < 1e-12)
check("fix2_r_ms_DK91_hi", abs(float(v6.r_ms(2.0)) - 1.33 * 2 ** 0.555) < 1e-12)
rs = np.array([0.05, 0.2, 0.68, 1.2]); ms = np.full(4, 0.902)
f = v6.apply_rs_floor(rs, ms, 0.68, 0.902)
check("fix2_floor_values", np.allclose(f, np.maximum(rs, 0.5 * 1.06 * 0.902 ** 0.945)) and f.min() > 0.47, f.round(3))
wd = v6.apply_rs_floor(np.array([0.013, 0.012]), np.array([0.6, 0.6]), 0.0131, 0.6)
check("fix2_wd_exempt", np.array_equal(wd, [0.013, 0.012]))
df = pd.read_csv(os.path.join(ROOT, "data", "pscomppars_2026-11-18.csv")); d = v6._tiers(u.prepare(df))
row = next(r for r in d.itertuples() if r.pl_name == "Kepler-399 d")
P = {**u.V3, **v6.V6}
o_on = v6.simulate_planet_v6(row, 20000, np.random.default_rng(1), cfg=dict(v6.CFG_V6), P=P)
o_off = v6.simulate_planet_v6(row, 20000, np.random.default_rng(1), cfg=dict(v6.CFG_V6, rs_floor=False), P=P)
check("fix2_kepler399_hycean_tail_removed", o_on["hycean"]["P_alive"] < o_off["hycean"]["P_alive"],
      f"P_alive H {o_off['hycean']['P_alive']:.4f} -> {o_on['hycean']['P_alive']:.4f}")
row_t = next(r for r in d.itertuples() if r.pl_name == "TRAPPIST-1 d")
a = v6.simulate_planet_v6(row_t, 5000, np.random.default_rng(3), cfg=dict(v6.CFG_V6), P=P)["rocky"]
b = v6.simulate_planet_v6(row_t, 5000, np.random.default_rng(3), cfg=dict(v6.CFG_V6, rs_floor=False), P=P)["rocky"]
check("fix2_trappist_unchanged", a["UHI_rank"] == b["UHI_rank"] and a["S_med"] == b["S_med"])
check("fix2_compat_off", "rs_floor" not in v6.CFG_COMPAT or not v6.CFG_COMPAT["rs_floor"])

# ---- fix 3 ----
C0 = load_old("composition"); rng = np.random.default_rng(5); n = 200000
rp = rng.uniform(1.0, 3.0, n); mst = rng.uniform(0.1, 0.8, n); uu = rng.uniform(size=(3, n)); z = rng.standard_normal(n); mv = rng.uniform(3, 6, n)
mp = rng.lognormal(1.0, 0.8, n); rpv = C.otegi_radius(mp, False, z, z)
new = C.classify_radius(rp, mst, uu[0], uu[1], z, 0.14, mp, (1, 1, 1))[0]; old = C0.classify_radius(rp, mst, uu[0], uu[1], z, 0.14, mp, (1, 1, 1))[0]
check("fix3_radius_111_identical", np.array_equal(new, old))
new = C.classify_mass(mp, uu[0], uu[1], mv, 0.3, rpv, (1, 1, 1))[0]; old = C0.classify_mass(mp, uu[0], uu[1], mv, 0.3, rpv, (1, 1, 1))[0]
check("fix3_mass_111_identical", np.array_equal(new, old))
for pi in [(2.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 2.0, 1.0)]:
    piV = 0.5 * (pi[1] + pi[2])
    ph = 1 / (1 + np.exp(-(rp - C.r_valley(mst, z)) / 0.14)); exp_R = np.mean(pi[0] * (1 - ph) / (pi[0] * (1 - ph) + piV * ph))
    got = np.mean(C.classify_radius(rp, mst, uu[0], uu[1], z, 0.14, mp, pi)[0] == 0)
    check(f"fix3_radius_posterior_{pi}", abs(got - exp_R) < 4 * np.sqrt(exp_R * (1 - exp_R) / n), f"P_R {got:.4f} vs {exp_R:.4f}")
    pr = 1 / (1 + np.exp(-(np.log(mv) - np.log(mp)) / 0.3)); exp_R = np.mean(pi[0] * pr / (pi[0] * pr + piV * (1 - pr)))
    got = np.mean(C.classify_mass(mp, uu[0], uu[1], mv, 0.3, rpv, pi)[0] == 0)
    check(f"fix3_mass_posterior_{pi}", abs(got - exp_R) < 4 * np.sqrt(exp_R * (1 - exp_R) / n), f"P_R {got:.4f} vs {exp_R:.4f}")
pr2 = np.mean(C.classify_radius(rp, mst, uu[0], uu[1], z, 0.14, mp, (2, 1, 1))[0] == 0); pr1 = np.mean(C.classify_radius(rp, mst, uu[0], uu[1], z, 0.14, mp, (1, 1, 1))[0] == 0)
check("fix3_A1_raises_PR_radius", pr2 > pr1, f"{pr1:.3f} -> {pr2:.3f}")

# ---- fix 5 ----
# Ishigami (a=7, b=0.1): S_T exactos 0.5576, 0.4424, 0.2437 (Saltelli+2010)
from SALib.sample import sobol as ss
prob = dict(num_vars=3, names=["x1", "x2", "x3"], bounds=[[-np.pi, np.pi]] * 3)
Xs = ss.sample(prob, 1024, calc_second_order=False, seed=1)
f = lambda x: np.sin(x[:, 0]) + 7 * np.sin(x[:, 1]) ** 2 + 0.1 * x[:, 2] ** 4 * np.sin(x[:, 0])
Y = np.c_[f(Xs), np.where(np.arange(len(Xs)) == 7, 1.0, 0.0), np.full(len(Xs), 3.0)]
res, flags, N = SR.analyze(Y, ["rocky:ish", "rocky:casi0", "rocky:const"], 3, branches=("rocky",))
st = res["rocky"]["m"]
check("fix5_ishigami_ST", np.allclose(st, [0.5576, 0.4424, 0.2437], atol=0.04), st.round(3))
check("fix5_ishigami_CI_covers", all(res["rocky"]["lo"][i] <= v <= res["rocky"]["hi"][i] for i, v in enumerate([0.5576, 0.4424, 0.2437])))
check("fix5_degenerate_excluded", flags[0] == "" and flags[1].startswith("casi") and flags[2] == "constante", flags)
print(f"TODOS LOS TESTS v6.1 OK ({len(OK)})")
