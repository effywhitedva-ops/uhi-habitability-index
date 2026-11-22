r"""Fit an empirical mass-route class prior for a sensitivity scenario.

Model: P_R(M)=sigmoid((ln M_v - ln M)/s). The sample comprises snapshot planets with measured radius,
M/sigma_M>3, and M<30 M_E. A planet is labelled R if its radius is at or below the Zeng et al. (2019)
50% water line at pl_eqt (missing values set to 300 K and clipped to 300-1000 K).

Caveat: transit+RV samples favour dense, short-period planets. The printed M_v [M_E] and s are used only
in run_uhi_v6.py variants pA4_empM and pA5_empM_P5, never in the production model."""
import os, sys, warnings; warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, uhi_code as u, composition as C
from scipy.optimize import minimize
d = u.prepare(pd.read_csv(os.path.join(HERE, "..", "data", "pscomppars_2026-11-18.csv")))
s = d[d.dens_ok & d.radius_measured & (d.pl_bmasse < 30)].copy()
s["R"] = (s.pl_rade.values <= C.r_h2o(s.pl_bmasse.values, s.pl_eqt.fillna(300).clip(300, 1000).values)).astype(int)
def fit(x):
    r"""Fit the logistic class prior by maximum likelihood and return a summary dictionary.

    x is a DataFrame with pl_bmasse in M_E and Boolean label R. Bounds keep M_v between 0.5 and 30 M_E
    and s between 0.05 and 2 natural-log units. Output: M_v, s, N and the fraction labelled R."""
    lm = np.log(x.pl_bmasse.values); y = x.R.values
    p = lambda q: 1 / (1 + np.exp(-(q[0] - lm) / q[1]))
    nll = lambda q: -np.sum(y * np.log(p(q) + 1e-12) + (1 - y) * np.log(1 - p(q) + 1e-12))
    r = minimize(nll, [np.log(4), 0.3], bounds=[(np.log(0.5), np.log(30)), (0.05, 2)])
    return dict(M_v=round(float(np.exp(r.x[0])), 2), s=round(float(r.x[1]), 2), N=len(x), frac_R=round(float(y.mean()), 3))
for name, x in [("todos", s), ("M*<0.6", s[s.st_mass < 0.6]), ("M*<0.6 y P>5 d", s[(s.st_mass < 0.6) & (s.pl_orbper > 5)])]:
    print(name, fit(x))
