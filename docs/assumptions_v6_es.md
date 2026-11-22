# UHI v6 Fase 1: supuestos y fuentes (2026-11-20, hora de El Salvador)

Base: `src/uhi_code.py` = copia idéntica de v5.2 (md5 `0288ce510341d50e9b3a769307ca2652`). Snapshot `data/pscomppars_2026-11-18.csv` (md5 `f2f70b24e67a9054940267e8081c81e8`), semilla 42, 10,000 draws.
**SUPUESTO** = elección propia sin fuente directa. **NO VERIFICADO** = la fuente no se pudo abrir o extraer; se usa el fallback indicado.

## Datos externos descargados (md5)
| Archivo | Fuente | md5 |
|---|---|---|
| `ext/fs255_grid.tar.gz` | Pistas de Spada+2013 (las de Mors), Zenodo 15729101 | f76987cf3d1da50435547f44a484a97f |
| `ext/Mors` | github ColinPhilipJohnstone/Mors, commit e4cfd253 (MIT). Parche para Python 3.13 en `patches/mors_py313_exec_fix.patch` (cambia `exec` por un closure; no cambia la física) | — |
| `data/mors_grid.npz` | Malla precalculada con Mors (`src/build_mors_grid.py`): M★ 0.10–1.20 (paso 0.05) × percentiles de rotación 5–95 × 300 edades | 2d735cdcc9df4216a017a563a139fcbc |
| `data/zeng2019/*.txt` | Tablas M–R de Zeng+2019 (lweb.cfa.harvard.edu/~lzeng): Earth-like, MgSiO3, 50% y 100% H2O a 300/500/700/1000 K | ver `md5sum data/zeng2019/*` |
| `data/EEM_dwarf_UBVIJHK_colors_Teff.txt` | Pecaut & Mamajek 2013, versión 2022.04.16 | afa3e0c8b249e8b9511b917c03aceb84 |

## D1: transición climática por rotación y acoplamiento de marea
- p_sync = σ((M_lock − M★)/0.05) · σ((a_lock − a)/0.03), con M_lock ~ U(0.55, 0.65) M☉ y a_lock ~ U(0.15, 0.25) AU. Es la regla de v5.2 (M < 0.6, a < 0.2) suavizada. **SUPUESTO.** La fórmula de Gladman+1996/Barnes 2017 ya está **VERIFICADA**, pero no se usa en Fase 1; el contraste está en «Verificación de los cinco fallbacks», punto 1.
- Eficacia 3D e(P) = clip((P − P_lo)/(P_hi − P_lo), 0, 1), con P_lo ~ U(4, 6) d y P_hi ~ U(15, 25) d. Fuentes: Kopparapu+2016 (transición hacia 10 d) y Haqq-Misra+2018 (rápido < 5 d, Rhines 5–20 d, lento > 20 d). La forma lineal es **SUPUESTO**.
- En los draws sincrónicos: k_in = 2[1 − 0.25(1 − e)] (Yang+2013 "el doble"; Kopparapu+2016 "hasta 25% menos") y eficacia de Yang = 0.5 + 0.5e. No se aplica el veto 1-D de Goldblatt. En los draws no sincrónicos: EBM 1-D con veto de 282 W/m².
- `locked` de v5.2 se conserva solo para el prior Beta de S_mag. **SUPUESTO**: así no se altera la secuencia del generador.
- Variante `teff`: w3D = p_sync · σ((T_c − Teff)/ΔT), con T_c ~ U(3700, 3900) K y ΔT ~ U(50, 150) K; eficacia y k_in aleatorios como en v5.2.

## D2: XUV acumulado en S_atm
- M★ completamente convectivas: Pass+2025, Ec. 6–11 y Tabla 1. log R_sat ~ N(−3.05, 0.06); t_sat = (5.9 − 15.4 M★)·lnN(0, 0.3) Gyr, con piso de 0.01 Gyr y M★ < 0.1 tratada como 0.1 (**SUPUESTO**). f_L se interpola de la Tabla 1, con valor constante fuera de 0.10–0.30 (**SUPUESTO**). f_flare ~ U(1, 1.5) (Pass usa 1.25). TRAPPIST-1 se integra saturada hasta su edad (Pass, Sec. 2.2).
- Tras la saturación, L_X/L_bol = 10⁻²·R_sat. **SUPUESTO**: Pass solo dice que la contribución es despreciable.
- M★ > 0.35: Johnstone+2021 (Mors). Se usa la integral acumulada de L_X (banda ROSAT, igual que Pass), reescalada por L_obs/L_pista(edad) (**SUPUESTO**). Percentil de rotación ~ U(5, 95) (**SUPUESTO**). Masas > 1.20 M☉ se fijan en 1.20, porque Mors falla a 1.25 M☉ (edad fuera de la pista).
- Mezcla logarítmica w = σ((M★ − 0.35)/0.03). A 5 Gyr y 0.35 M☉ los dos modelos difieren solo +0.16 dex (`figures/fig_v6_xuv_blend.png`).
- Normalización: I_E = 2.0×10¹⁸ erg cm⁻² (Pass+2025). Control: el Sol con Mors da 1.53 I_E y con Jackson+2012 1.35 I_E.
- Retención: P_ret = Φ((ARM + c)·ln10/σ), con ARM = 4 log v_esc − log(I/I_E) − 3.16 (Pass Ec. 11), c = 0 y σ = 1.0 ln. Son los mismos valores base que la costa de v5.2 (**SUPUESTO**).
- Variante `jackson`: Jackson+2012 Tabla 2 con (B−V)₀. El snapshot no trae B ni V, así que (B−V)₀ sale de Pecaut & Mamajek según Teff (**SUPUESTO**: enana de secuencia principal, sin enrojecimiento). Fuera de 0.29–1.41 se usa Mors.
- Variante `bt3d` (Berta-Thompson+2026): log S_crit = 5.9 log v_esc + 1.17 log L★ + const. La constante está verificada de forma indirecta, con la afirmación log L = −2.23 del paper (punto 5 de la verificación). Se fija para que v_esc terrestre y L★ = 10^−2.23 queden sobre la costa en S = 1. Es solo variante.

## D3: XUV contado una vez
- S_atm = P_ret acumulada (D2). S_estrella mantiene la fórmula de v5.2 (estrés XUV actual y flares). S_geo = S_régimen, sin la bandera de erosión de sub-Neptunos.

## D4: clases por draw (R = sin envoltura de H₂, rocoso o mundo de agua; H = Hycean; G = mini-Neptuno rico en gas)
- Ruta de **densidad** (radio medido y M/σ_M > 3):
  - R es compatible si Rp ≤ R_50%H2O(M, T_eq), usando Zeng+2019 interpolado en T entre 300 y 1000 K, constante fuera.
  - H es compatible si R_MgSiO3(M) < Rp ≤ R_Hmax(M).
  - G es compatible si Rp > R_MgSiO3(M).
  - Se asignan priors iguales (D4) entre las hipótesis compatibles; si ninguna es compatible, la clase es R.
  - Usar el 50% de agua como límite de R y el silicato puro como mínimo de los volátiles es **SUPUESTO**. Con la línea Earth-like y 100% agua, TRAPPIST-1 d salía 1/3–1/3–1/3, lo cual no es físico.
- R_Hmax: Madhusudhan+2021, 2.3 R⊕ a 5 M⊕ y 2.6 R⊕ a 10 M⊕. Lineal en log M y extrapolado (**SUPUESTO**).
- Ruta de **radio**: P(H₂) = σ((Rp − R_v)/w), con R_v = 1.8197 M★^0.14 (+0.02/−0.01) (Venturini+2024) y w ~ U(0.08, 0.20) R⊕ (**SUPUESTO**). La masa sin medir sale de Otegi+2020, rocosa o volátil según la clase, con coeficientes muestreados dentro de sus errores. Así v6 ya no usa la masa de Chen & Kipping para el score. Dentro de H₂: H si Rp ≤ R_Hmax, con prior 1/2.
- Ruta de **masa** (sin radio; evidencia "prior"): P_R = σ((ln M_v − ln M)/0.3), con M_v ~ U(3, 6) M⊕. **SUPUESTO completo.** El radio sale de Otegi según la clase.
- Rama Rocky: S_fis sin borde superior de radio (la clase ya decide), μ = 28 o 44 y S_atm = P_ret XUV.
- Rama Hycean:
  - S_fis,H = trap(M; 1, 2, 10, 15) (**SUPUESTO**).
  - S_clima,H = 1 si S ≤ f_cl · S_crit(P_H2, Teff) · 0.88/(1 − A_H), y 0 si no. S_crit sale de Innes+2023 Tabla 1 (1 y 10 bar; estrella M = 3457 K, G = 5772 K): log-lineal en P, lineal en Teff y constante fuera de ese rango. log P_H2 ~ U(0, 2) bar y f_cl ~ U(1, 2) (**SUPUESTO**). No hay borde externo: Madhusudhan+2021 admite Hyceans fríos.
  - S_atm,H = Φ(ln(t_loss / 1 Gyr)). Se calcula con M_env = 4πR²P/g y Ṁ de Owen & Wu 2017, η = 0.1(v_esc/15)⁻²·lnN(0, 0.5). Usa el XUV actual (L_X + EUV de Sanz-Forcada+2011, o L_EUV de Mors). Se usa K = 1 (**SUPUESTO**). La fórmula de Erkaev ya está verificada y su efecto es despreciable con Innes (punto 3 de la verificación).
  - S_obs con μ = 2.3 y A_B = 0.12.
- Variante `madhu` (optimista): veto si T_eq(A = 0) > 500 K. Cuál es el albedo que corresponde a ese T_eq queda **NO VERIFICADO** y se toma A = 0.
- Sin cambios respecto a v5.2: tier nominal (calculado con masas de Chen & Kipping), S_mag (masa de v5.2) y ESI.

## D5: rankings
- Tres salidas: Rocky = UHI_R·[clase = R], Hycean = UHI_H·[clase = H] y combinado (la rama de la clase del draw; G = 0).
- Ranking hurdle igual que v5.2. P_top10 ahora cuenta solo los draws vivos. En la rama Hycean, la convención de v5.2 (rango de muerto = vivos + 1) daba P_top10 = 1 artificial porque hay menos de 10 vivos por draw. CFG_COMPAT conserva la convención de v5.2.
- snr_thr es un escenario: corridas separadas a 5σ y 3σ, cada una con su propio Sobol.

## Generador aleatorio
- La secuencia principal es idéntica a v5.2: se consumen ck_mass, ck_radius y Beta de S_mag aunque se reemplacen.
- Los draws nuevos salen de default_rng(crc32(nombre) + 2), con un bloque fijo de 20 uniformes y 8 normales.

## Pruebas de sensibilidad al prior de clase (añadidas el 2026-11-20 a pedido de Effy)
Todas usan la configuración principal (Innes, Mors, rot) y cambian una sola cosa. Corren a 10,000 draws, 5σ y 3σ. Código: `run_uhi_v6.py --variant pXX`, `src/prior_sensitivity_v6.py`.
- **(a) Prior de clase**:
  - A1 favorece lo rocoso: π_R:π_H:π_G = 2:1:1 y M_v ~ U(5, 10).
  - A2 favorece lo volátil: 0.5:1:1 y M_v ~ U(1.5, 3).
  - A3 usa M_v = 2.04 M⊕, la transición Terran→Neptunian de Chen & Kipping 2017 (la misma constante que ya usa `ck_mass` de v5.2).
  - A4 usa un prior empírico: ajuste logístico a los planetas del snapshot con dens_ok, M★ < 0.6 (N = 101) y M < 30 M⊕. Da M_v = 4.10 M⊕, s = 0.40 (`src/fit_empirical_prior.py`). Tiene **SESGO** de detección por tránsito+RV, así que es solo un escenario.
  - A5 es igual que A4 con P > 5 d (N = 44): M_v = 2.83, s = 0.70.
  - El prior de ocurrencia en radio (Fulton/Petigura) **no** se implementó: convertirlo a P(R | M sin i) exige una M–R y una completitud que no verifiqué hoy.
- **(b) M sin i → M**: el código de v5.2 (y el de v6) **ya** muestrea cos i ~ U(0, 1), que es inclinación isotrópica. Eso equivale a un prior **plano** en M (lo comprobé analíticamente: la densidad de muestreo es proporcional a la verosimilitud por un factor constante igual a M sin i). Variantes:
  - B1: M = M sin i (cota inferior).
  - B2: prior de función de masa dN/dlogM ∝ M^−0.48 (Howard+2010, Science 330, 653, doi:10.1126/science.1194854, verificado en el resumen y en el texto). Se aplica repesando por (M / M sin i)^(α−1), es decir solo la proyección. Repesar la masa total divergía hacia M → 0 por la cola de error de M sin i.
  - Advertencia: ese ajuste viene de FGK con P < 50 d y aquí se aplica a enanas M (**SUPUESTO**).
  - Las 15 entradas "Msin(i)/sin(i)" y los planetas transitantes con M sin i se tratan como masa verdadera (igual que en v5.2).
  - La condición de "no transita" (i lejos de 90°) no se impone; su efecto es pequeño (**SUPUESTO**).
- **(c) M–R**: C1 usa Chen & Kipping 2017 (la masa y el radio de v5.2) en vez de Otegi+2020 para las magnitudes sin medir.
- **(d) Frontera rocosa** (solo afecta la ruta de densidad):
  - D1, línea Earth-like única: R si Rp ≤ R_Earthlike, volátil si no.
  - D2: R hasta 100% H2O.
  - El principal usa 50% H2O, con los volátiles por encima de MgSiO3.
- **Métricas**, siempre contra el principal:
  - Solapamiento del top 10 y del top 20.
  - Kendall τ y Spearman sobre la unión de los top 20; Spearman sobre los vivos.
  - Rango [min, max] de posición, UHI_rank, P_R y P_top10 por planeta.
- **Veredicto por planeta**:
  - "robusto": está en el top 10 en todos los escenarios.
  - "estable": está en el top 20 en todos y la posición varía ≤ 10.
  - "DEPENDE DEL PRIOR": cualquier otro caso.

## Verificación de los cinco fallbacks (2026-11-20, 16:00–17:00 hora de El Salvador)
Scripts: `src/analysis/a6_fallbacks.py`; salidas en `results/analysis_phase1/` (`fallback_checks.json`, `zeng_md5_check.txt`).

1. **Fórmula del bloqueo de marea (Gladman+1996; Barnes 2017)**: **fórmula VERIFICADA** como t_lock ≈ ω a⁶ I Q / (3 G M★² k₂ R⁵) (Gladman+1996, Icarus 122, 166; la misma forma de marea de equilibrio que usa Barnes 2017, CeMDA 129, 509). **La regla p_sync sigue siendo SUPUESTO** y ahora tiene un contraste cuantitativo. Con un planeta terrestre, Q = 100, k₂ = 0.3, P₀ = 24 h e I = 0.33 MR², t_lock < 1 Gyr hasta a ≈ 0.22 AU (0.1 M☉), 0.32 (0.3 M☉), 0.40 (0.6 M☉) y 0.48 AU (1 M☉). Con Q entre 10 y 1000, el radio va de 0.15–0.33 a 0.33–0.70 AU. Es decir, el corte a_lock ~ U(0.15, 0.25) AU y M★ < 0.6 M☉ es **conservador**: según Gladman, también se bloquean planetas de la HZ de K tardías y de M tempranas a 0.2–0.4 AU. No se cambia en Fase 1 porque las mareas térmicas atmosféricas (Leconte+2015) pueden impedir la sincronización, y porque Effy aprobó D1 con la regla suavizada. En el Sobol, lock_a_c tiene S_T ≈ 0.07 en Rocky, así que no es despreciable. Propuesta para Fase 2: p_sync = P(t_lock(Q, k₂, P₀) < edad), con log Q ~ U(1, 3).
2. **Tablas de Zeng+2019**: **VERIFICADAS**. Los 10 archivos de `data/zeng2019/` tienen el mismo md5 que los de lweb.cfa.harvard.edu/~lzeng/tables/, descargados de nuevo hoy. Como control de cordura, la línea Earth-like da 1.000 R⊕ a 1 M⊕. Lo que queda como SUPUESTO es el uso que hacemos de ellas: 50% H₂O como límite de R, MgSiO₃ como mínimo de los volátiles e interpolación en T.
3. **K de Erkaev+2007**: **fórmula VERIFICADA** como K(ξ) = 1 − 3/(2ξ) + 1/(2ξ³), con ξ = (M_p/3M★)^{1/3} a/R_p (Erkaev+2007, A&A 472, 329, Ec. 17; Salz+2016, Ec. 2). El código usa **K = 1 (SUPUESTO)**, que subestima Ṁ en 1/K. Calculada en los 242 Hycean vivos (Innes), K tiene mediana 0.958 y mínimo 0.861, y 10 planetas quedan por debajo de 0.9. En el top 10 de Innes, K va de 0.956 a 0.994. Equivale a desplazar ln t_loss en −0.04 (mediana), dentro del ancho σ = 1 de S_atm,H, así que el efecto en los rankings es despreciable. En la variante de Madhusudhan, el mínimo es 0.615 y 153 de 1213 planetas quedan por debajo de 0.9; ahí convendría incluir K en Fase 2.
4. **MIST en vez de BHAC15**: **VERIFICADO.** No hace falta ninguno de los dos modelos en el código. f_L sale directamente de la Tabla 1 de Pass+2025, y la nota de esa tabla dice "We estimate the PMS correction, f_L, from MIST models" (suppdata de ApJL 986 L3). Comprobé la tabla completa: f_L = 1.19, 1.17, 1.13, 1.13, 1.16 y 1.22 para 0.10–0.30 M☉, y t_sat = 5.9 − 15.4 M★ reproduce 4.4/4.1/3.6/2.8/2.1/1.3 Gyr. Además, recalculando la Ec. 10 con f_flare = 1.25, R_sat = 10^−3.05 e I_E = 2.0×10¹⁸ erg cm⁻² salen I_XUV,corr = 0.108/a² (Pass: 0.11/a²) a 0.10 M☉ y 0.379/a² (Pass: 0.38/a²) a 0.30 M☉, lo que **verifica también la normalización I_E**. Para M★ > 0.35, Mors usa las pistas de Spada+2013, no BHAC15. BHAC15 queda sin uso.
5. **Normalización de Berta-Thompson+2026 (variante `bt3d`)**: **VERIFICADA de forma indirecta.** El paper (arXiv 2507.02136v2; ApJ, doi:10.3847/1538-4357/ae77fd) dice: "If we take f_hz = f_⊕ as a crude approximation … f_shoreline < f_hz for stars less luminous than log10(L★/L☉) = −2.23 (+0.18/−0.21)" para planetas del tamaño de la Tierra. Nuestra constante impone exactamente S_crit(v_esc,⊕, L = 10^−2.23) = 1, así que coincide con la normalización implícita del paper. El intercepto numérico propio del ajuste sigue sin extraerse. La incertidumbre de ±0.2 dex en log L se traduce en ±0.23 dex en log S_crit. Sigue siendo solo una variante.

### Otros números revisados en el análisis
- **Innes+2023, Tabla 1** (versión publicada en ApJ 953, 168; suppdata): con 1 bar, 0.389 S₀ (G) y 0.320 S₀ (M); con 10 bar, 0.0676 S₀ (G) y **0.0860 S₀ (M)**. El código usa 0.0852 (116 W m⁻²/1361), que parece salir del preprint. La diferencia es del 0.9%, con efecto despreciable. Conviene corregirla en la próxima corrida. Innes no tiene caso húmedo a 100 bar, así que la extrapolación de log P_H₂ entre 1 y 2 es SUPUESTO, como ya estaba marcado. El espectro M de Innes es de unos 3500 K y el código usa 3457 K, otra diferencia menor.
- **Compatibilidad con v5.2 (9.2e-15)**: corregí la explicación. **No** viene del redondeo del CSV, porque el CSV de v5.2 guarda repr de 17 cifras y relee exacto. Viene del **entorno**: el propio `uhi_code.py` de v5.2 (md5 idéntico), reejecutado hoy en este box (numpy 2.5.4, Python 3.13.5), difiere del CSV publicado exactamente en 9.16e-15 (UHI_cond_med) y 1.4e-16 (UHI_rank) (`results/analysis_phase1/compat_origin_check.txt`). v6 con CFG_COMPAT es idéntico bit a bit a v5.2 corrido en el mismo entorno.
- **Colas no físicas del muestreo de la estrella (heredado de v5.2)**: `split_normal` recorta R★ en 0.05 R☉. Cuando el catálogo da σ(R★)/R★ ≳ 0.3 (sobre todo en Kepler), una fracción de draws tiene L ∝ R★² hasta unas 100 veces menor. Ese caso es una estrella de 5500 K con 0.05 R☉, que no es físico. Es **SUPUESTO heredado, ahora con impacto cuantificado**: 181 de los 242 Hycean vivos de Innes y 338 de los 495 Rocky vivos solo viven por esas colas (su S mediana está fuera del criterio climático). El top 20 Rocky no se ve afectado (el mejor caso es el #36). Sí aparecen en el top 10 Hycean de Innes a 5σ: Kepler-399 d, Kepler-230 c y Kepler-386 c, con UHI ≤ 0.006. Propuesta: muestrear R★ en log o truncar a ±3σ con un piso físico según T_eff, y reportar los vivos con P_alive ≥ 0.01. Con ese umbral quedan Rocky 140, Hycean 22 y combinado 166.
- **Prior de clase π_R en las variantes A1/A2**: `classify_radius` y `classify_mass` **no usan π_R**; solo usan π_H/(π_H+π_G). Por eso A1 (2:1:1) y A2 (0.5:1:1) cambian π_R únicamente en la ruta de densidad. En la ruta de masa actúan a través de M_v, y en la ruta de radio no actúan. Afecta solo a la interpretación de A1/A2; el principal (1:1:1) no cambia. Queda **pendiente de decisión** cómo aplicar π_R en la ruta de radio.

## v6.1: correcciones aprobadas por Effy (2026-11-20, 4:22 PM hora de El Salvador)
- **Fix 1. Compatibilidad con v5.2:** Effy acepta la diferencia de 9.2e-15 como «reproduce v5.2». Se debe al entorno numérico (ver la sección de verificación).
- **Fix 2. Piso físico de R★ (`apply_rs_floor` en `uhi_code_v6.py`; `CFG_V6["rs_floor"] = True`, apagado en CFG_COMPAT):**
  - Regla: R★ ≥ f·R_MS(M★,draw), con R_MS de Demircan & Kahraman 1991 (Ap&SS 181, 313): R = 1.06 M^0.945 si M < 1.66 M☉ y 1.33 M^0.555 si no (verificado hoy). Usamos f = 0.5.
  - Por qué un piso y no muestreo en log: el muestreo logarítmico con σ_ln ≈ σ/R ≈ 0.5 (Kepler) sigue dando colas de un factor 4–5 hacia abajo, sin ningún límite físico. El piso ataca el problema real: una estrella que quema H no puede ser mucho más compacta que la secuencia principal de su masa. Además se aplica recortando, sin consumir draws, así que la secuencia del generador es idéntica a v6.0.
  - Por qué f = 0.5 (**SUPUESTO conservador**): las subenanas pobres en metales y las estrellas ZAMS quedan a no menos de unos 0.7–0.9 R_MS. Con 0.5 solo se eliminan los draws imposibles y no se recorta la incertidumbre legítima.
  - El piso **solo** se aplica si el radio de catálogo es compatible con secuencia principal (R_cat ≥ 0.5 R_MS(M_cat)). Así no se tocan enanas blancas ni anfitriones compactos raros.
  - El radio del planeta (pl_rade) se sigue muestreando aparte, como en v5.2. Correlacionarlo con R★ queda para la Fase 2.
  - Los conjuntos de vivos se reportan además con **P_alive ≥ 0.01** (`summary_*.json["alive_P001"]`).
- **Fix 3. π_R también en las rutas de radio y masa (`composition.prior_vol`, `classify_mass`):**
  - La probabilidad de volátil p, calculada con el valle de radios o con la sigmoide de masa, se trata como posterior con prior plano R:V = 1:1 y se reponderan las odds: P′(V) = π_V p / (π_R(1−p) + π_V p).
  - π_V = (π_H + π_G)/2 (**SUPUESTO**). Se eligió así para que 1:1:1, la D4 aprobada, deje la corrida principal idéntica bit a bit (lo comprueba un test).
  - Ahora A1 (2:1:1) y A2 (0.5:1:1) actúan en las tres rutas.
  - En el Sobol, log2_piH también desplaza R frente a V, porque π_V cambia.
- **Fix 4.** Innes+2023, Tabla 1 publicada: 10 bar, estrella M = **0.0860 S₀** (antes 0.0852).
- **Fix 5.** `sobol_report_v6.py` excluye salidas constantes, casi siempre iguales (menos del 10% de evaluaciones distintas del valor modal) e inestables (S_T > 1.5). El IC95 de la media sale de un bootstrap Jansen (B = 500), y se reportan ancho completo y semiancho.
- **Fix 6.** El p_sync basado en el t_lock de Gladman+1996/Barnes 2017 queda para la **Fase 2**. En v6.1 sigue la regla suavizada (SUPUESTO conservador).
