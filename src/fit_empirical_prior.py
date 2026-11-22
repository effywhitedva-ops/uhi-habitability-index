"""Prior empirico de la ruta de masa: ajuste logistico P_R(M) = sigmoide((ln M_v - ln M)/s) a los planetas del snapshot con
radio medido y M/sigma_M > 3 (dens_ok), M < 30 M_E, etiquetados R si Rp <= linea 50% H2O de Zeng+2019 a T_eq (pl_eqt, 300-1000 K).
SESGO: muestra de transito+RV (favorece planetas densos/cortos); se usa solo como escenario de sensibilidad."""
import os, sys, warnings; warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, uhi_code as u, composition as C
from scipy.optimize import minimize
d = u.prepare(pd.read_csv(os.path.join(HERE, "..", "data", "pscomppars_2026-11-18.csv")))
s = d[d.dens_ok & d.radius_measured & (d.pl_bmasse < 30)].copy()
s["R"] = (s.pl_rade.values <= C.r_h2o(s.pl_bmasse.values, s.pl_eqt.fillna(300).clip(300, 1000).values)).astype(int)
def fit(x):
    lm = np.log(x.pl_bmasse.values); y = x.R.values
    p = lambda q: 1 / (1 + np.exp(-(q[0] - lm) / q[1]))
    nll = lambda q: -np.sum(y * np.log(p(q) + 1e-12) + (1 - y) * np.log(1 - p(q) + 1e-12))
    r = minimize(nll, [np.log(4), 0.3], bounds=[(np.log(0.5), np.log(30)), (0.05, 2)])
    return dict(M_v=round(float(np.exp(r.x[0])), 2), s=round(float(r.x[1]), 2), N=len(x), frac_R=round(float(y.mean()), 3))
for name, x in [("todos", s), ("M*<0.6", s[s.st_mass < 0.6]), ("M*<0.6 y P>5 d", s[(s.st_mass < 0.6) & (s.pl_orbper > 5)])]:
    print(name, fit(x))
