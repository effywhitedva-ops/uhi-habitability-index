"""UHI v5.2 RC: TODAS las tablas y estadisticas desde los CSV, sin pasos manuales.
Entradas: uhi_results_v52.csv, uhi_results_v52_3sigma.csv, anchors_v52*.csv, sweep_v52/*.csv, phl/hwc.csv (opcional),
          uhi_results_v51.csv (opcional, para el efecto de la rampa), pscomppars.csv.
Salidas: tables_v52/*.md, tables_v52/*.csv, tables_v52/stats_v52.json.
Uso: python make_tables.py [directorio_base]"""
import os, sys, json, hashlib
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from scipy.special import ndtri
B = sys.argv[1] if len(sys.argv) > 1 else "."
P = lambda *a: os.path.join(B, *a)
OUT = P("tables_v52"); os.makedirs(OUT, exist_ok=True)
AVAIL = {"T0": 1.0, "T1": 0.75, "T2": 0.5, "T3": 0.25}
# ---- umbrales explicitos del benchmark ----
ESI_HI = 0.8          # "ESI alto" (PHL: ESI > 0.8 ~ "Earth-like"; Schulze-Makuch+2011)
PDET_UNOBS = 1e-3     # "SNR_eff ~ 0": P_det mediana = S_obs_med / Availability(tier) < 1e-3  (SNR_eff mediana < 5 - 3.09 = 1.91 a 5 sigma)
PDET_RESC = 0.5       # rescate por luz reflejada: no transita (sin TSM/ESM JWST) y P_det mediana >= 0.5 (SNR_eff mediana >= umbral) via ELT/PCS
S = {}


def kopp(teff):
    T = np.clip(teff, 2600, 7200) - 5780.0
    co = {"rv": (1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15), "rg": (1.107, 1.332e-4, 1.580e-8, -8.308e-12, -1.931e-15),
          "mg": (0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16), "em": (0.320, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16)}
    return {k: s + a * T + b * T ** 2 + c * T ** 3 + d * T ** 4 for k, (s, a, b, c, d) in co.items()}


def add_cols(r):
    k = kopp(r.st_teff.values)
    s = r.S_med.values
    r["HZ_Kopp"] = np.where((s <= k["rg"]) & (s >= k["mg"]), "conservadora", np.where((s <= k["rv"]) & (s >= k["em"]), "optimista", "fuera"))
    r["P_det_med"] = (r.S_obs_med / r.tier.map(AVAIL)).clip(0, 1)
    return r


def snr_eff(p, thr):
    p = np.clip(np.asarray(p, float), 1e-300, 1 - 1e-16)
    return thr + ndtri(p)


def md(df, cols, fmt=None):
    fmt = fmt or {}
    h = "| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n"
    for _, q in df.iterrows():
        h += "| " + " | ".join((fmt[c](q[c]) if c in fmt else str(q[c])) for c in cols) + " |\n"
    return h


f3 = lambda x: f"{x:.3f}" if pd.notna(x) else "—"
f2 = lambda x: f"{x:.2f}" if pd.notna(x) else "—"
fi = lambda x: f"{int(x)}" if pd.notna(x) else "—"

S["catalog_md5"] = hashlib.md5(open(P("pscomppars.csv"), "rb").read()).hexdigest()
r5 = add_cols(pd.read_csv(P("uhi_results_v52.csv")).sort_values("rank_pos").reset_index(drop=True))
r3 = add_cols(pd.read_csv(P("uhi_results_v52_3sigma.csv")).sort_values("rank_pos").reset_index(drop=True))
S["n_scored"] = len(r5); S["n_alive_5s"] = int((r5.P_alive > 0).sum()); S["n_alive_3s"] = int((r3.P_alive > 0).sum())
S["top10_5s"] = r5.planet.head(10).tolist(); S["top10_3s"] = r3.planet.head(10).tolist()
m53 = r5[["planet", "UHI_rank", "rank_pos", "P_alive"]].merge(r3[["planet", "UHI_rank", "rank_pos", "P_alive"]], on="planet", suffixes=("_5", "_3"))
live53 = (m53.P_alive_5 > 0) | (m53.P_alive_3 > 0)
S["spearman_5v3_live"] = float(spearmanr(m53.UHI_rank_5[live53], m53.UHI_rank_3[live53])[0])

# ---------- 1. anclas ----------
an_rows = []
for tag, f in [("5σ", "anchors_v52.csv"), ("3σ", "anchors_v52_3sigma.csv")]:
    if not os.path.exists(P(f)): continue
    a = pd.read_csv(P(f))
    for _, q in a.iterrows():
        nm = q.planet.split()[1]
        ok = (q.pct_vs_catalog > 99.0) if nm == "Earth" else (q.P_alive == 0 and q.UHI_rank == 0)
        an_rows.append(dict(run=tag, planet=q.planet, UHI_rank=q.UHI_rank, P_alive=q.P_alive, pct=q.pct_vs_catalog,
                            pos=int(q.pos_if_inserted), criterio="pct>99" if nm == "Earth" else "P_alive=0", pasa="SÍ" if ok else "NO"))
an = pd.DataFrame(an_rows); an.to_csv(os.path.join(OUT, "anchors_v52.csv"), index=False)
S["anchors_all_pass"] = bool((an.pasa == "SÍ").all())
S["anchors"] = an.to_dict("records")

# ---------- 2. top 20 ----------
def top_table(r, n=20):
    t = r.head(n).copy()
    return md(t, ["rank_pos", "planet", "UHI_rank", "P_alive", "UHI_cond_med", "P_top10", "rank_ci5", "rank_ci95", "tier", "transiting", "HZ_Kopp"],
              dict(UHI_rank=f3, P_alive=f2, UHI_cond_med=f3, P_top10=f2, rank_ci5=fi, rank_ci95=fi))
r5.head(20).to_csv(os.path.join(OUT, "top20_v52_5sigma.csv"), index=False); r3.head(20).to_csv(os.path.join(OUT, "top20_v52_3sigma.csv"), index=False)

# ---------- 3. efecto de la rampa (v5.1 -> v5.2) ----------
ramp_md = ""
if os.path.exists(P("uhi_results_v51.csv")):
    r51 = pd.read_csv(P("uhi_results_v51.csv"))
    m = r51[["planet", "UHI_rank", "rank_pos", "P_alive", "albedo_med"]].merge(r5[["planet", "UHI_rank", "rank_pos", "P_alive", "albedo_med", "S_med", "st_teff"]], on="planet", suffixes=("_51", "_52"))
    lv = (m.P_alive_51 > 0) | (m.P_alive_52 > 0)
    S["ramp_spearman_v51_v52_live"] = float(spearmanr(m.UHI_rank_51[lv], m.UHI_rank_52[lv])[0])
    S["ramp_top10_same_set"] = set(r51.sort_values("rank_pos").planet.head(10)) == set(S["top10_5s"])
    S["ramp_top10_same_order"] = r51.sort_values("rank_pos").planet.head(10).tolist() == S["top10_5s"]
    m["d_alb"] = m.albedo_med_52 - m.albedo_med_51
    ch = m[np.abs(m.d_alb) > 1e-9].copy(); S["ramp_n_albedo_changed"] = int(len(ch))
    ch["S0_Wm2"] = ch.S_med * 1361.0
    ch = ch.reindex(ch.UHI_rank_52.sort_values(ascending=False).index)
    ch.to_csv(os.path.join(OUT, "ramp_changes_v51_v52.csv"), index=False)
    ramp_md = md(ch.head(15), ["planet", "st_teff", "S0_Wm2", "albedo_med_51", "albedo_med_52", "rank_pos_51", "rank_pos_52", "P_alive_51", "P_alive_52", "UHI_rank_51", "UHI_rank_52"],
                 dict(st_teff=lambda x: f"{x:.0f}", S0_Wm2=lambda x: f"{x:.0f}", albedo_med_51=f3, albedo_med_52=f3, P_alive_51=f2, P_alive_52=f2, UHI_rank_51=f3, UHI_rank_52=f3))

# ---------- 4. barrido 1-D del veto ----------
sw_md = ""
SW = P("sweep_v52")
if os.path.exists(os.path.join(SW, "base.csv")):
    base = pd.read_csv(os.path.join(SW, "base.csv")).sort_values("rank_pos")
    b10 = base.planet.head(10).tolist()
    rows = []
    order = ["Frg_275", "Frg_280", "base", "Frg_286", "Frg_290", "Frg_295", "Trg_3700", "Trg_3750", "Trg_3850", "Trg_3900",
             "Tmd_3700", "Tmd_3750", "Tmd_3850", "Tmd_3900"]
    for tag in order:
        f = os.path.join(SW, tag + ".csv")
        if not os.path.exists(f): continue
        x = pd.read_csv(f).sort_values("rank_pos"); a = pd.read_csv(os.path.join(SW, tag + "_anchors.csv"))
        mm = base[["planet", "UHI_rank", "P_alive"]].merge(x[["planet", "UHI_rank", "P_alive"]], on="planet", suffixes=("_b", "_x"))
        lv = (mm.P_alive_b > 0) | (mm.P_alive_x > 0)
        rho_l = spearmanr(mm.UHI_rank_b[lv], mm.UHI_rank_x[lv])[0]; rho_a = spearmanr(mm.UHI_rank_b, mm.UHI_rank_x)[0]
        t10 = x.planet.head(10).tolist()
        k2 = x[x.planet == "K2-125 b"].iloc[0]
        ven = a[a.planet.str.contains("Venus")]; ear = a[a.planet.str.contains("Earth")]; mar = a[a.planet.str.contains("Mars")]
        par = {"base": "F_rg=282, T_rg=3800"}.get(tag, tag.replace("Frg_", "F_rg=").replace("Trg_", "T_rg=").replace("Tmd_", "T_md=T_rg="))
        ok = bool((ven.P_alive == 0).all() and set(t10) == set(b10) and rho_l > 0.99)
        rows.append(dict(config=par, tag=tag, n_alive=int((x.P_alive > 0).sum()), venus_P_alive_max=float(ven.P_alive.max()),
                         earth_pct_min=float(ear.pct_vs_catalog.min()), mars_P_alive_max=float(mar.P_alive.max()),
                         top10_overlap=len(set(t10) & set(b10)), top10_same_set=set(t10) == set(b10), top10_same_order=t10 == b10, spearman_live=rho_l, spearman_all=rho_a,
                         K2_125b_pos=int(k2.rank_pos), K2_125b_P_alive=k2.P_alive, K2_125b_UHI_rank=k2.UHI_rank,
                         K2_125b_in_top10=bool("K2-125 b" in t10), pasa=ok,
                         out_of_top10=";".join(sorted(set(b10) - set(t10))), into_top10=";".join(sorted(set(t10) - set(b10)))))
    sw = pd.DataFrame(rows); sw.to_csv(os.path.join(OUT, "sweep_veto_v52.csv"), index=False)
    S["sweep"] = sw.to_dict("records")
    prim = sw[sw.tag.str.startswith(("Frg", "Trg", "base"))]
    S["sweep_primary_all_pass"] = bool(prim.pasa.all()); S["sweep_secondary_all_pass"] = bool(sw[sw.tag.str.startswith("Tmd")].pasa.all())
    sw_md = md(sw, ["config", "n_alive", "venus_P_alive_max", "earth_pct_min", "top10_overlap", "top10_same_set", "top10_same_order", "spearman_live", "spearman_all",
                    "K2_125b_pos", "K2_125b_P_alive", "K2_125b_UHI_rank", "pasa", "out_of_top10", "into_top10"],
               dict(venus_P_alive_max=f3, earth_pct_min=f2, spearman_live=lambda x: f"{x:.4f}", spearman_all=lambda x: f"{x:.4f}",
                    K2_125b_P_alive=f3, K2_125b_UHI_rank=f3))

# ---------- 5. benchmark ESI / PHL / Kopparapu ----------
live = r5.P_alive > 0
S["spearman_UHI_ESI_all"] = float(spearmanr(r5.UHI_rank, r5.ESI_med)[0])
S["spearman_UHI_ESI_live"] = float(spearmanr(r5.UHI_rank[live], r5.ESI_med[live])[0])
S["spearman_UHI3_ESI_live"] = float(spearmanr(r3.UHI_rank[r3.P_alive > 0], r3.ESI_med[r3.P_alive > 0])[0])
S["HZ_Kopp_counts_all"] = r5.HZ_Kopp.value_counts().to_dict(); S["HZ_Kopp_counts_live"] = r5[live].HZ_Kopp.value_counts().to_dict()
S["HZ_Kopp_counts_top20"] = r5.head(20).HZ_Kopp.value_counts().to_dict()
phl_ok = os.path.exists(P("phl", "hwc.csv"))
r5["ESI_rank_SM11"] = r5.ESI_med.rank(ascending=False, method="min").astype(int)
if phl_ok:
    h = pd.read_csv(P("phl", "hwc.csv"))[["P_NAME", "P_ESI", "P_HABITABLE", "P_UPDATE"]].drop_duplicates("P_NAME")
    S["phl_rows"] = int(len(h)); S["phl_md5"] = hashlib.md5(open(P("phl", "hwc.csv"), "rb").read()).hexdigest(); S["phl_last_update"] = str(h.P_UPDATE.max())
    r5 = r5.merge(h, left_on="planet", right_on="P_NAME", how="left")
    r5["PHL_ESI_rank"] = r5.P_ESI.rank(ascending=False, method="min")
    mt = r5.P_ESI.notna(); S["phl_matched"] = int(mt.sum())
    S["spearman_UHI_PHLESI_matched"] = float(spearmanr(r5.UHI_rank[mt], r5.P_ESI[mt])[0])
    S["spearman_UHI_PHLESI_live_matched"] = float(spearmanr(r5.UHI_rank[mt & (r5.P_alive > 0)], r5.P_ESI[mt & (r5.P_alive > 0)])[0])
    S["spearman_ESI_SM11_vs_PHLESI"] = float(spearmanr(r5.ESI_med[mt], r5.P_ESI[mt])[0])
    ph = r5[r5.P_HABITABLE > 0]; S["phl_habitable_in_catalog"] = int(len(ph)); S["phl_habitable_alive"] = int((ph.P_alive > 0).sum())
    S["phl_habitable_in_UHI_top20"] = int((ph.rank_pos <= 20).sum()); S["phl_habitable_in_UHI_top50"] = int((ph.rank_pos <= 50).sum())
    S["UHI_top20_in_PHL"] = int(r5.head(20).P_ESI.notna().sum()); S["UHI_top20_PHL_habitable"] = int((r5.head(20).P_HABITABLE > 0).sum())
r5["UHI3_pos"] = r5.planet.map(r3.set_index("planet").rank_pos)
cols = ["rank_pos", "UHI3_pos", "planet", "UHI_rank", "ESI_med", "ESI_rank_SM11"] + (["P_ESI", "PHL_ESI_rank", "P_HABITABLE"] if phl_ok else []) + ["HZ_Kopp", "transiting", "P_det_med"]
r5.head(20)[cols].to_csv(os.path.join(OUT, "benchmark_top20_v52.csv"), index=False)
bench_md = md(r5.head(20), cols, dict(UHI_rank=f3, ESI_med=f3, P_ESI=f3, PHL_ESI_rank=fi, P_HABITABLE=fi, UHI3_pos=fi, P_det_med=f3))
# falsos positivos inobservables
for lab, col in [("SM11", "ESI_med")] + ([("PHL", "P_ESI")] if phl_ok else []):
    hi = r5[r5[col] >= ESI_HI]; un = hi[hi.P_det_med < PDET_UNOBS]
    S[f"highESI_{lab}_n"] = int(len(hi)); S[f"highESI_{lab}_unobs_n"] = int(len(un))
    S[f"highESI_{lab}_unobs_rate"] = float(len(un) / max(len(hi), 1))
    S[f"highESI_{lab}_unobs_UHIrank_max"] = float(un.UHI_rank.max()) if len(un) else 0.0
    S[f"highESI_{lab}_unobs_alive_n"] = int((un.P_alive > 0).sum())
    S[f"highESI_{lab}_unobs_best_pos"] = int(un.rank_pos.min()) if len(un) else None
    un = un.assign(SNR_eff_med=snr_eff(un.P_det_med, 5.0)).sort_values(col, ascending=False)
    un.to_csv(os.path.join(OUT, f"unobservable_highESI_{lab}_v52.csv"), index=False)
fp = pd.read_csv(os.path.join(OUT, "unobservable_highESI_SM11_v52.csv"))
fp_md = md(fp.head(15), ["planet", "ESI_med"] + (["P_ESI"] if phl_ok else []) + ["HZ_Kopp", "transiting", "tier", "sy_dist", "P_det_med", "rank_pos", "UHI_rank"],
           dict(ESI_med=f3, P_ESI=f3, sy_dist=lambda x: f"{x:.0f}", P_det_med=lambda x: f"{x:.1e}", UHI_rank=f3))
# rescatados por luz reflejada
res = r5[(r5.transiting == 0) & (r5.P_alive > 0) & (r5.P_det_med >= PDET_RESC)].copy()
res["SNR_eff_med"] = snr_eff(res.P_det_med, 5.0)
res.to_csv(os.path.join(OUT, "reflected_rescued_v52.csv"), index=False)
S["rescued_reflected_n_live"] = int(len(res)); S["rescued_reflected_in_top20"] = int((res.rank_pos <= 20).sum())
S["rescued_reflected_in_top50"] = int((res.rank_pos <= 50).sum()); S["rescued_reflected_top20_list"] = res[res.rank_pos <= 20].planet.tolist()
S["nontransiting_in_top20"] = int((r5.head(20).transiting == 0).sum())
rs_md = md(res.head(20), ["rank_pos", "planet", "UHI_rank", "P_alive", "P_det_med", "sy_dist", "a_med", "HZ_Kopp", "ESI_med"] + (["P_ESI", "P_HABITABLE"] if phl_ok else []),
           dict(UHI_rank=f3, P_alive=f2, P_det_med=f3, sy_dist=f2, a_med=lambda x: f"{x:.4f}", ESI_med=f3, P_ESI=f3, P_HABITABLE=fi))

# ---------- 6. escritura ----------
r5.drop(columns=[c for c in ["P_NAME"] if c in r5.columns]).to_csv(os.path.join(OUT, "uhi_results_v52_with_benchmark.csv"), index=False)
json.dump(S, open(os.path.join(OUT, "stats_v52.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o), ensure_ascii=False)
J = lambda k: S.get(k, "—")
txt = f"""# Tablas UHI v5.2-final (generadas por make_tables.py; no editar a mano)

Catálogo md5 {J('catalog_md5')} · {J('n_scored')} planetas puntuados · vivos 5σ {J('n_alive_5s')} · vivos 3σ {J('n_alive_3s')} · Spearman 5σ vs 3σ (vivos) {S['spearman_5v3_live']:.4f}

## 1. Anclas del Sistema Solar a 10 pc (criterios: Tierra > percentil 99; Venus y Marte P_alive = 0)
{md(an, ['run','planet','UHI_rank','P_alive','pct','pos','criterio','pasa'], dict(UHI_rank=f3, P_alive=f3, pct=f2))}
Todas pasan: **{J('anchors_all_pass')}**

## 2. Top 20 v5.2 (5σ)
{top_table(r5)}
## 3. Top 20 v5.2 (3σ)
{top_table(r3)}
## 4. Efecto de la rampa w_eff(S0) (v5.1 → v5.2)
Spearman v5.1 vs v5.2 (vivos en alguna): {S.get('ramp_spearman_v51_v52_live', float('nan')):.4f} · top 10 mismo conjunto: {J('ramp_top10_same_set')} · mismo orden: {J('ramp_top10_same_order')} · planetas con albedo mediano cambiado: {J('ramp_n_albedo_changed')}

{ramp_md}
## 5. Barrido 1-D del veto híbrido (uno a la vez; base F_rg = 282 W m⁻², T_rg = 3800 K)
Criterio: Venus sintético P_alive = 0 (tr y nt), top 10 con la misma composición que la base, Spearman (vivos) > 0.99.
T_rg mueve solo la frontera del veto; T_md=T_rg (secundario) mueve también todo el régimen M-sincrónico (Yang, k_in, eyeball).

{sw_md}
Primario (F_rg y T_rg) pasa: **{J('sweep_primary_all_pass')}** · secundario (T_md) pasa: **{J('sweep_secondary_all_pass')}**

## 6. Benchmark ESI / PHL / Kopparapu
- ESI_SM11: Schulze-Makuch et al. 2011, 4 propiedades (R, ρ, v_esc, T_s = T_eq + 33 K), exponente w_i/4, mediana Monte Carlo.
- PHL: Habitable Worlds Catalog (hwc.csv, última actualización {J('phl_last_update')}; md5 {J('phl_md5')}); su ESI es la versión radio + flujo estelar.
- HZ: Kopparapu et al. 2014 (1 M⊕; conservadora = [max greenhouse, runaway greenhouse], optimista = [early Mars, recent Venus]) con S mediana y T_eff nominal.
- Spearman UHI_rank vs ESI_SM11: todos {S['spearman_UHI_ESI_all']:.4f} (N={J('n_scored')}); vivos {S['spearman_UHI_ESI_live']:.4f} (N={J('n_alive_5s')}); vivos 3σ {S['spearman_UHI3_ESI_live']:.4f}
"""
if phl_ok:
    txt += f"""- Spearman UHI_rank vs ESI PHL: emparejados {S['spearman_UHI_PHLESI_matched']:.4f} (N={J('phl_matched')}); vivos emparejados {S['spearman_UHI_PHLESI_live_matched']:.4f}; ESI_SM11 vs ESI PHL {S['spearman_ESI_SM11_vs_PHLESI']:.4f}
- PHL potencialmente habitables en nuestro catálogo: {J('phl_habitable_in_catalog')} (vivos en UHI: {J('phl_habitable_alive')}; en top 20 UHI: {J('phl_habitable_in_UHI_top20')}; en top 50: {J('phl_habitable_in_UHI_top50')}). Top 20 UHI presentes en PHL: {J('UHI_top20_in_PHL')}, marcados habitables: {J('UHI_top20_PHL_habitable')}
"""
txt += f"""- HZ Kopparapu (todos): {J('HZ_Kopp_counts_all')} · (vivos): {J('HZ_Kopp_counts_live')} · (top 20): {J('HZ_Kopp_counts_top20')}

### 6.1 Top 20 UHI v5.2 vs ESI/PHL
{bench_md}
### 6.2 Falsos positivos inobservables (ESI ≥ {ESI_HI} y P_det mediana < {PDET_UNOBS:g}, es decir SNR_eff mediana < 1.91 a 5σ)
- ESI_SM11: {J('highESI_SM11_unobs_n')} de {J('highESI_SM11_n')} ({100*S['highESI_SM11_unobs_rate']:.1f} %); vivos entre ellos: {J('highESI_SM11_unobs_alive_n')}; UHI_rank máx {S['highESI_SM11_unobs_UHIrank_max']:.3f}; mejor posición {J('highESI_SM11_unobs_best_pos')}
"""
if phl_ok:
    txt += f"- ESI PHL: {J('highESI_PHL_unobs_n')} de {J('highESI_PHL_n')} ({100*S['highESI_PHL_unobs_rate']:.1f} %); vivos entre ellos: {J('highESI_PHL_unobs_alive_n')}; UHI_rank máx {S['highESI_PHL_unobs_UHIrank_max']:.3f}; mejor posición {J('highESI_PHL_unobs_best_pos')}\n"
txt += f"""
{fp_md}
### 6.3 Rescatados por luz reflejada (no transitan, vivos, P_det mediana ≥ {PDET_RESC} vía ELT/PCS)
Total vivos: {J('rescued_reflected_n_live')} · en top 20: {J('rescued_reflected_in_top20')} · en top 50: {J('rescued_reflected_in_top50')}

{rs_md}"""
txt += """
## 7. Limitación aceptada: frontera T_eff del veto (v5.2-final, sin ajuste)
El barrido 1-D de la sección 5 pasa para F_rg ∈ [275, 295] W m⁻² (incluye 282 W m⁻² de Goldblatt et al. 2013 y 286 W m⁻² del corrigendum 2014) y para T_rg = 3850 y 3900 K. **Falla** el criterio de composición invariante del top 10 en T_rg = 3700 y 3750 K: K2-125 b (T_eff = 3654 K) sale del #10 y entra Ross 128 b.

En la corrida certificada (T_rg = 3800 K) el margen es marginal: UHI_rank 0.333 (K2-125 b, #10) frente a 0.323 (Ross 128 b, #11). Los otros 9 del top 10 no cambian en el barrido primario de T_rg. Detalle numérico (columna K2_125b_* de la sección 5):

- T_rg = 3700 K: K2-125 b pasa a #18 (UHI_rank 0.265, P_alive 0.652).
- T_rg = 3750 K: K2-125 b pasa a #13 (UHI_rank 0.306, P_alive 0.766).
- En ambos valores, Venus sintético mantiene P_alive = 0 y Spearman (vivos) es 0.9973 y 0.9993, ambos > 0.99.

Esta sensibilidad se documenta como **limitación aceptada**. No se retocó T_rg (sigue en 3800 K = teff_md) para hacer pasar el criterio.

El análisis de Sobol post-v5 no incluye F_rg ni la rampa w_eff(S0). La errata del ESI (exponente w/8 → w/4) afecta solo a los valores de ESI de documentos anteriores, no al UHI.
"""
open(os.path.join(OUT, "tables_v52.md"), "w").write(txt)
print(json.dumps({k: v for k, v in S.items() if k not in ("sweep", "anchors")}, indent=1, default=str, ensure_ascii=False))
