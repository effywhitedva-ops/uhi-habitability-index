# UHI v2.0 — corrida sobre el NASA Exoplanet Archive: supuestos (código v5.2 RC)

> **Actualizado el 18 nov 2026 para el código v4** (`uhi_code.py`, resultados en `uhi_results_v4.csv`; variante 3σ en `uhi_results_v4_3sigma.csv`). v4 reemplaza el proxy de observabilidad de los planetas no transitantes (log_cref, snr_adv, A_g = 0.3 fijo, IWA = 2λ/D a 1 µm, d ≤ 20 pc) por un modelo físico ELT/PCS de luz reflejada (punto 10b). Todo lo demás es igual que en v3 (5σ principal). El código v3 está en `uhi_code_v3.py` y este documento en versión v3 en `assumptions_v3.md`.

> **Actualizado el 18 nov 2026 para el código v3** (`uhi_code.py`, resultados en `uhi_results_v3.csv`). v3 reemplaza las logísticas de S_atm y S_obs (pendiente k_log = 4) por un modelo probabilístico de detección (puntos 10 y 11). Umbral de detección **decidido por Effy: 5σ principal, 3σ variante** (`uhi_results_v3_3sigma.csv`). Se conservan el albedo dependiente de T_eff (v2) y fR_d = 2.5. Versiones anteriores: `uhi_code_v1.py` → `uhi_results.csv`; `uhi_code_v2.py` → `uhi_results_v2.csv` (este documento en versión v2: `assumptions_v2.md`).

Corrida del 18 nov 2026. Fuentes de la especificación, leídas sin modificarlas: *Índice Unificado de Habitabilidad (UHI) — Versión 2.0* (UHI_v2_integrado1) y *UHI Mathematical Toolkit — T0–T2* (UHI_T0_Toolkit).

## A. Lo que se tomó tal cual de los documentos
- UHI = Π S_i^{w_i}, con Σw_i = 1 y w_i ≥ 0 (media geométrica ponderada). Si un S_i vale 0, el UHI vale 0 exactamente, porque es un producto.
- Siete subíndices: S_fis, S_clima, S_atm, S_estrella, S_geo, S_mag y S_obs.
- Pesos tomados de un prior de Dirichlet con α = [2,2,2,2,1,1,2] (Toolkit §7). En cada draw de Monte Carlo se muestrea un vector de pesos nuevo, así la incertidumbre de los pesos también se propaga.
- Monte Carlo sobre las incertidumbres del catálogo. Se reportan la mediana y los intervalos 68 % (p16–p84) y 95 % (p2.5–p97.5).
- Fórmulas del Toolkit §1: a por la 3.ª ley de Kepler; S = L/a²; T_eq = ((1−A)L/(16πσa²f))^{1/4}; δ = (Rp/R*)². L* se calcula desde R* y T_eff (tabla 0.1).
- S_clima: fHZ (fracción de la superficie entre 273 y 373 K según un EBM), con la clase de HZ de Kopparapu como prior.
- S_atm y S_obs: TSM/ESM de Kempton et al. (2018); H = kT/(μg) con escenarios de μ (N₂/CO₂/H₂); SNR ∝ (H/Rp)·δ·√N_tr/σ_phot; S_obs = logistic(SNR_eff) × Availability. **En v3 la "logistic" se reemplaza por Φ(SNR − umbral)** (ver puntos 10 y 11).
- S_estrella: F_XUV ∝ L_XUV/a², con L_XUV desde la relación actividad–edad, más el riesgo de fulguraciones.
- S_geo: calentamiento de marea y fotoevaporación energy-limited, esta última como bandera de erosión para sub-Neptunos.
- S_mag: prior Beta amplio, actualizado solo con restricciones indirectas.
- Tiers T0–T3 tomados de la §4 del documento UHI v2.0 y del Anexo A.

## B. Supuestos propios (los documentos no los especifican)
**Datos y muestreo**
1. Tabla `pscomppars` completa: 6470 filas. Se puntúa un planeta solo si tiene T_eff, R*, M*, (P o a), (Rp o Mp) y al menos una magnitud J o K. Con eso quedan 6112 puntuados y 358 excluidos. Los planetas marcados como controvertidos se incluyen, con la columna `controversial`.
2. Cada parámetro se muestrea con una normal partida, usando err1 y err2. Si falta la incertidumbre, se asume: P 0.01 %, Rp 10 % (15 % si es estimado), Mp 20 %, T_eff 2 %, R* 7 %, M* 7 %, a 5 %.
3. Excentricidad: si falta, e = 0. Si hay valor pero no error, σ = 0.05. Siempre se recorta a [0, 0.95].
4. Edad estelar: si falta (1529 filas), se toma U(1, 10) Gyr.
5. Masa:
   - Se trata como medida si `pl_bmassprov` es "Mass" o "Msin(i)/sin(i)", o si es "Msini" en un planeta que transita (sin i ≈ 1).
   - Para "Msini" en un planeta que no transita: M = Msini/sin i, con inclinación isotrópica.
   - Si la masa viene de "M-R relationship", se recalcula desde el Rp muestreado con Chen & Kipping (2017), incluyendo la dispersión intrínseca. La rama joviana se toma log-uniforme entre 131.6 y 4000 M⊕.
   - Si el radio está marcado "Calculated Value" o el planeta no transita, el radio sale de la masa con Chen & Kipping (2017).
6. Si falta J se usa K, y si falta K se usa J (J≈K).
7. **Albedo de Bond (código v2, decidido por Effy: dependiente de la T_eff).** En v1 era A ~ U(0.2, 0.4) para todos los planetas; esa versión se conserva en `uhi_code_v1.py`. Ahora:
   - **Fuente:** Shields et al. 2013, *Astrobiology* 13:715–739, Tabla 2, columna "Earth-like gases, clouds, Rayleigh scattering on" (64 % de cobertura de nubes). Es el albedo planetario de banda ancha calculado con los espectros de HD128167 (F2V), el Sol (G2V), HD22049 (K2V) y AD Leo (M3V).
   - **Valores que tomo:** océano 0.32865 (F), 0.31948 (G), 0.30235 (K), 0.23372 (M); tierra firme 0.41428 (F), 0.41484 (G), 0.40133 (K), 0.33165 (M).
   - **Supuesto, T_eff de las estrellas ancla:** Shields et al. no dan T_eff en esa tabla. Asigné AD Leo 3477 K y σ Boo/HD128167 6700 K (dentro del rango publicado de 6594 a ~6800 K); HD22049 (ε Eri) 5020 K, tomado del NASA Exoplanet Archive; y el Sol 5772 K.
   - **Supuesto, interpolación:** lineal en T_eff entre las 4 anclas. Fuera de [3477, 6700] K uso el valor del extremo, de modo que las enanas M más frías que AD Leo (por ejemplo TRAPPIST-1, con 2566 K) reciben los valores M3V.
   - **Supuesto, mezcla de superficie:** A = A_océano + f_tierra·(A_tierra − A_océano), con f_tierra ~ U(0, 0.5) en cada draw (la Tierra tiene ~0.29). La T_eff también se muestrea con su incertidumbre en cada draw.
   - **Rangos resultantes:** M ≤3477 K: 0.234–0.283; K (5020 K): 0.302–0.352; G (5772 K): 0.319–0.367; F ≥6700 K: 0.329–0.372.
   - **Superficie sin hielo:** el EBM sigue sin retroalimentación hielo-albedo. Shields et al. 2013 y Joshi & Haberle 2012 (*Astrobiology* 12:3–8) muestran que esa retroalimentación es más débil alrededor de enanas M. No tomé valores numéricos de Joshi & Haberle; agregar el hielo queda como mejora pendiente.
   - **Usos del mismo albedo:** T_eq, el EBM (fHZ), la escala de altura de S_obs (v3: con el T_eq de la escala de altura; el TSM/ESM usa el T_eq de Kempton con A = 0) y la T_s del ESI (por eso también cambian algunos puestos del ESI).
   - **Redistribución:** f = 1.

**Subíndices**
8. **S_fis** = √(f_R·f_M), con trapecios: Rp (0.3, 0.5, 1.6, 2.5) R⊕ (**fR_d = 2.5: decidido por Effy**, 18 nov 2026) y Mp (0.03, 0.1, 5, 10) M⊕. El documento solo dice "logistic/trapezoid, physical thresholds".
9. **S_clima** = fHZ × p_HZ.
   - El EBM es 1-D de dos modos (North 1975; North, Cahalan & Coakley 1981): insolación media anual con oblicuidad terrestre (s₂ = −0.482), OLR = A + B·T con B = 2.09 W m⁻² °C⁻¹, D = 0.649, y albedo constante, sin retroalimentación hielo-albedo.
   - El intercepto A ~ U(150, 203.3) W m⁻² representa un invernadero incierto, desde el terrestre hasta uno reforzado por CO₂.
   - El EBM no modela rotación sincrónica, aunque esa es la situación probable en los planetas HZ de enanas M.
   - p_HZ vale 1 dentro de la HZ conservadora de Kopparapu et al. (2014, 1 M⊕), 0.5 solo dentro de la optimista y 0.1 fuera. Los valores 1/0.5/0.1 son propios.
10. **S_atm (v3).** Lo que dicen los documentos: Toolkit §3, "TSM/ESM per Kempton 2018; H = kT/(μg)"; tabla de normalización, "S_atm: TSM/ESM, H | Logistic | SNR targets | per instrument class". El documento UHI v2.0 solo dice "logística, tramos o kernels". No fijan pendiente ni centro. La pendiente k_log = 4 de v1/v2 era mía y, según el Sobol de la tarea 2, el 3.er supuesto más influyente. v3 la reemplaza por:
    - **S_atm = P(retener atmósfera) × P(detectar atmósfera).**
    - **P(detectar) = Φ(SNR − SNR_umbral)**, con Φ la CDF normal estándar. Esto es la probabilidad de que la SNR observada (SNR verdadera + ruido gaussiano unitario) supere el umbral. Es una probabilidad física, no una logística arbitraria.
    - **Transitantes:** P_det = 1 − (1 − Φ(SNR_t − umbral))·(1 − Φ(SNR_e − umbral)), es decir, se detecta en transmisión o en emisión.
      - SNR_t = TSM'·√(t_prog/10 h)·m₁·c.
        - **Kempton et al. 2018 (PASP 130:114401):** el TSM está normalizado para dar aproximadamente la S/N de un programa de 10 h con JWST/NIRISS (la mitad en tránsito), sin nubes. TSM' aplica la re-escala ×2.3/18 de Kempton (§3.1) a los planetas de 1.5–2.0 R⊕, de modo que su TSM supone una atmósfera de μ = 18 (como los de Rp < 1.5) en lugar de μ = 2.3. Los factores de escala del TSM (0.190 / 1.26 / 1.28 / 1.15) también son de Kempton.
        - c ~ U(cloud_min, 1) es la supresión por nubes o neblina (Kempton supone atmósfera sin nubes). **Supuesto:** cloud_min = 0.5.
      - SNR_e = ESM·√n_ecl·m₂. **Kempton 2018:** el ESM está normalizado a S/N = 7.5 para un eclipse de GJ 1132b en la banda MIRI LRS, y "will take at least two eclipses" para llegar a S/N > 10. Por eso **n_ecl = 2** (el valor exacto es mío).
      - m_i = exp(σ_mod·z_i), con z_i ~ N(0,1): incertidumbre de modelo de la SNR predicha. **Supuesto:** σ_mod = 0.5 (factor ~1.65).
      - **SNR_umbral: decidido por Effy (18 nov 2026), 5σ principal y 3σ variante.** La corrida principal (`uhi_results_v3.csv`) usa 5; la variante (`uhi_results_v3_3sigma.csv`) usa 3, con todo lo demás igual. Comparación en `compare_5v3sigma.png`. El mismo umbral se aplica a S_atm y S_obs. Con TSM = 10 (el corte de Kempton para Rp < 1.5) y sin nubes, P_det ≈ 1 con cualquiera de los dos umbrales.
    - **No transitantes (v4): modelo físico ELT/PCS de luz reflejada (ver 10b).** P_det = Φ(SNR_r − umbral), con SNR_r = SNR_HCI·m₁.
10b. **Observabilidad de los no transitantes (v4).** Reemplaza el proxy de v1–v3 (snr_adv·contraste/10^log_cref, A_g = 0.3, IWA = 2λ/D a 1 µm, d ≤ 20 pc).
    - **Contraste planeta–estrella:** C = A_g·Φ_L(α)·(Rp/a)².
      - Φ_L(α) = [sin α + (π − α)cos α]/π es la fase de Lambert. **Supuesto mío:** se observa en máxima elongación (α = 90°, Φ_L = 1/π, separación proyectada = a, para cualquier inclinación). Snellen et al. 2015 también suponen el planeta "half of the planet is seen to be illuminated by its host star, as expected at the greatest elongation"; la órbita es conocida por RV, así que la observación se puede programar. Se ignora la excentricidad.
      - A_g = (2/3)·A_Bond: relación de una esfera lambertiana. A_Bond es el albedo v2 de cada draw (Shields et al. 2013), así que A_g ≈ 0.16–0.25. **Supuesto:** albedo gris, el mismo en el óptico/NIR. Snellen et al. 2015 usan A_g gris = 0.3; la Tierra real no es lambertiana.
    - **Separación angular:** θ = a/d (a en AU y d = `sy_dist` en pc dan arcsec; se reporta en mas). Si `sy_dist` falta, SNR = 0.
    - **Ruido: halo estelar residual (Kasper et al. 2021, *The Messenger* 182:38, "PCS — A Roadmap for Exoearth Imaging with the ELT").** Kasper et al. explican que el ruido fotónico de la luz estelar dispersada a la posición del planeta es la fuente principal de ruido. Con un halo XAO residual de contraste h ≈ 5×10⁻⁵ a las menores separaciones, hacen falta 5×10³ fotones del planeta para 1σ con C = 10⁻⁸, y unos 10⁵ para 5σ. Eso equivale a SNR = C·√(N★/h), donde N★ son los fotones estelares recogidos. Uso esa fórmula (fotón-limitada, es decir, HDS + coronógrafo eliminan las motas cuasiestáticas; Kasper et al. suponen un contraste HDS mejor que 10⁻⁴).
    - **Fotones (Kasper et al. 2021):** un objeto de J = 26 da unos 10⁵ fotones en 10 h con el ELT (banda J, 50 nm, 10 % de throughput). Por eso N★ = 10⁵·10^(−0.4(J★ − 26))·(t/10 h), con J★ del catálogo (o K si falta J, como antes). Tiempo t = 10 h (el ejemplo de Kasper; Snellen et al. 2015 también simulan una noche, 10 h). **Supuesto:** se usa la banda J de Kasper aunque la banda A de O₂ esté en el óptico.
    - **Perfil del halo:** h(15 mas) = 5×10⁻⁵ (Kasper et al. 2021). **Supuesto mío:** el halo mejora un factor 10 de 15 a 100 mas (ley de potencias), igual que el requisito de contraste de PCS, que pasa de ~10⁻⁸ a 15 mas a 10⁻⁹ a 100 mas y más allá (Kasper et al. 2021). Es constante dentro de 15 mas y más allá de 100 mas.
    - **IWA:** λ/D = 0.75 µm/39 m ≈ 3.97 mas (λ = centro de 0.6–0.9 µm, el rango de Snellen et al. 2015; incluye la banda A de O₂). **Supuesto mío:** IWA ~ U(2, 4) λ/D en cada draw, es decir 7.9–15.9 mas. El límite superior es aproximadamente los 15 mas del requisito de PCS (~4 λ/D a 700 nm); el inferior, los ~2 λ/D que sugiere Effy. Dentro del IWA, SNR = 0. El número aleatorio del IWA sale de un generador propio por planeta (semilla crc32 del nombre), para no alterar la secuencia del resto del Monte Carlo.
    - **Verificación con valores nominales (A_Bond = 0.25):** Proxima b da θ = 37.3 mas, C = 4.7×10⁻⁸ y SNR ≈ 50 en 10 h. Del mismo orden, Snellen et al. 2015 obtienen S/N ≈ 10 en una noche para C = 6×10⁻⁷ a 25 mas con Strehl 0.3 (un sistema mucho menos eficiente que el XAO de PCS); Turbet et al. 2016 (A&A 596:A112) dan para Proxima b ~10⁻⁷ de contraste y 7 λ/D a 1 µm. Los tres coinciden en que Proxima b es accesible.
    - **Limitaciones (supuestos míos):** se ignoran el fondo de cielo (Snellen et al. 2015: más de 23 mag por debajo de la estrella por píxel en el óptico), la eficiencia real del HDS y los residuos de motas. Snellen et al. 2015 y Kasper et al. 2021 advierten que eliminar la luz estelar por debajo de 10⁻⁴ con HDS+HCI no está demostrado. Por eso el modelo es optimista.
    - La incertidumbre de modelo m₁ = exp(σ_mod·z) (σ_mod = 0.5) y el umbral (5σ principal, 3σ variante) son los mismos de v3.
    - **P(retener) = Φ(ln(I_costa/I)/σ_costa)**, la "línea de costa cósmica" de **Zahnle & Catling 2017 (ApJ 843:122)**: frontera empírica I ∝ v_esc⁴ entre cuerpos con y sin atmósfera.
      - Zahnle & Catling la trazan a ojo y no dan constante numérica. **Supuesto mío:** la normalizo para que pase por Marte (I = 1/1.524² = 0.431 S⊕, v_esc = 5.03 km/s), por un factor 10^shore_logc con shore_logc = 0.
      - v_esc = 11.186 km/s·√(Mp/Rp).
      - **Supuesto:** ancho σ_costa = 1.0 en ln I.
    - El tier ya no entra en S_atm.
11. **S_obs (v3).** Lo que dicen los documentos: Toolkit §3, "SNR_tran-spec ∝ (H/Rp)·δ·√N_tr/σ_phot; S_obs = logistic(SNR_eff)·Availability(instrument)"; tabla, "S_obs: SNR_eff, availability | Logistic × mask | SNR*, A_inst | Tier T0–T3 mapping". v3 conserva la estructura (SNR de rasgos espectrales × Availability) pero cambia la logística por **S_obs = Φ(SNR_eff − umbral) × Availability(tier)**, con Availability T0 = 1, T1 = 0.75, T2 = 0.5, T3 = 0.25 (sin cambios).
    - **Transitantes:**
      - Amplitud de la señal A = n_H·2H·Rp/R*², con H = kT_eq/(μg). μ = 28 o 44 (N₂/CO₂, 50/50) si Rp < 1.6, y 2.3 si no. **Supuesto:** n_H = 2.
      - Ruido fotónico tomado del propio TSM: S/N₁₀ = TSM·μ_K/μ, donde μ_K = 18 (Rp < 1.5) o 2.3 es el μ que supone Kempton/Louie 2018. Así σ_fot = (A/S/N₁₀)·√(5 h/(N_tr·T₁₄)).
      - T₁₄ = P/π·asin(R*/a), y N_tr = min(365.25/P, ntr_cap). **Supuesto:** ntr_cap = 10.
      - Piso sistemático: **Greene et al. 2016 (ApJ 817:17)** dan pisos de 20 ppm (NIRISS SOSS), 30 ppm (NIRCam) y 50 ppm (MIRI LRS) por tránsito o eclipse, que "no bajan al sumar más datos". Tomo **20 ppm**.
      - **Supuesto:** el piso es por canal y se promedia sobre n_bins = 10 canales independientes: σ_piso = 20 ppm/√10. Sin esto (n_bins = 1), el piso de 20 ppm frente a rasgos de ~40 ppm dejaba S_obs ≈ 0 en TRAPPIST-1.
      - SNR_eff = A/√(σ_fot² + σ_piso²)·m₃·c.
    - **No transitantes:** SNR_eff = SNR_r (el mismo valor que en S_atm, del modelo v4 ELT/PCS, punto 10b).
    - **Nota:** S_atm (¿hay y se detecta alguna atmósfera?) y S_obs (¿se miden rasgos de una atmósfera secundaria de μ alto con tránsitos acumulados y piso de ruido?) comparten los factores de nubes, de modo que no son independientes. La doble cuenta parcial ya estaba en el diseño de los documentos.
12. **S_estrella** = (1 − estrés_XUV) × (1 − penalización por fulguraciones).
    - L_XUV/L_bol = 10^−3.5 durante la saturación y luego ∝ (t/t_sat)^−1.23 (pendiente de Ribas et al. 2005).
    - t_sat = 0.1 Gyr si M* ≥ 0.6 y 1 Gyr si no.
    - F_XUV se normaliza a la Tierra actual.
    - Estrés = logística con un 50 % a 10× el valor terrestre.
    - Fulguraciones: 0.5 si M* < 0.35, 0.25 si 0.35 ≤ M* < 0.6, y 0 en otro caso. El catálogo no tiene datos de fulguraciones.
13. **S_geo** = trapecio sobre log10(q_rad + q_marea), con óptimo entre 0.03 y 1 W m⁻², y 0 por debajo de 0.01 o por encima de 10.
    - q_rad = 0.087·(Mp/Rp²)·2^{−(edad−4.5)/3}.
    - Marea: Ė = (21/2)(k₂/Q)·G·M*²·Rp⁵·**n**·e²/a⁶, con k₂ = 0.3 y Q = 100. **Nota:** la fórmula del Toolkit omite el movimiento medio n, y sin él no da unidades de potencia. Aquí se incluyó la n de la forma estándar.
    - Erosión, aplicada solo a sub-Neptunos (1.6 ≤ Rp ≤ 4): energy-limited con ε = 0.1 y K = 1, sobre una envoltura del 1 % de Mp. El riesgo es la logística de t_pérdida frente a 1 Gyr, y el factor aplicado es (1 − 0.5·riesgo).
14. **S_mag** ~ Beta(2, 2). Pasa a Beta(2, 3) si es probable la rotación sincrónica (M* < 0.6 y a < 0.2 AU), y β sube en 1 si Mp < 0.3 M⊕.

**Tier por planeta** (sobre valores nominales; el documento no da una regla operativa)
- **T0:** transita, masa medida con σ/M ≤ 0.5, y TSM ≥ umbral o ESM ≥ 7.5.
- **T1:** transita y max(TSM/umbral, ESM/7.5) ≥ 0.5.
- **T2:** transita con Rp ≥ 1.6 (ARIEL) y métrica ≥ 0.1; o no transita y (v4) su SNR nominal ELT/PCS (A_Bond = 0.3, IWA = 2λ/D) es ≥ 0.1 × umbral (0.5 con 5σ; 0.3 en la variante 3σ). Es lo análogo de «métrica ≥ 0.1» de los transitantes y no tiene corte de distancia. Con esta regla, los no transitantes T2 pasan de 95 (v3) a 662 (v4); muchos son gigantes alrededor de estrellas brillantes, con UHI = 0 por S_fis.
- **T3:** el resto.
- Vector por subíndice, constante y tomado del Anexo A: `fis:T0|clima:T1|atm:T1-T2|estrella:T0|geo:T3*|mag:T3*|obs:T0` (* = proxy).

**Ranking**
15. **v5 (hurdle):** se ordena por UHI_rank = P_alive × mediana(UHI | UHI > 0), con P_alive = P(UHI > 0); los empates se rompen por la media del posterior. Se reportan también la mediana y el RIC condicionales, P(rango ≤ 10) con ranking conjunto por draw y el intervalo creíble del 90 % del rango. (v1–v4 ordenaban por mediana; en v4 había 6057 planetas con mediana 0.)

## C. ESI (comparación)
Schulze-Makuch D. et al. (2011), *A Two-Tiered Approach to Assessing the Habitability of Exoplanets*, Astrobiology 11(10):1041–1052, doi:10.1089/ast.2010.0592.
- ESI = √(ESI_I·ESI_S), con ESI_I y ESI_S = Π(1 − |(x − x₀)/(x + x₀)|)^{w/n} y n = 2 por grupo. Equivale a la forma de 4 propiedades con exponente w/4.
- **Errata v5.2:** v1–v5.1 usaban w/4 dentro de cada grupo *y además* la raíz, es decir un exponente total w/8. Eso comprimía el ESI hacia 1 (mediana del catálogo 0.547 → 0.299 corregido). Comprobación con la fórmula corregida: Marte 0.70 y Venus 0.44, como en Schulze-Makuch+2011. El UHI no cambia porque el ESI no entra en él. Los valores de ESI de los Docs v4/v5 usan la fórmula vieja.
- Pesos: radio 0.57, densidad 1.07, velocidad de escape 0.70, temperatura superficial 5.58. Referencias terrestres: 1, 1, 1 y 288 K.
- **Supuesto:** T_s = T_eq + 33 K. El ESI se calcula en el mismo Monte Carlo y se reporta su mediana.

## D. Qué no se hizo (el documento lo pide, queda pendiente)
VoI, agregación Borda/Kemeny y validación cruzada de fHZ entre familias de modelos. Sí se hicieron Sobol y PRCC (con LHS) sobre los supuestos (`thresholds_for_effy.md`; la versión v3 está en `sens_out_v3/`). La corrida principal usa MC simple.


## B-v5. Cambios de la auditoría v4→v5 (spec de Neon/R-Aetherios)

Las marcas **[elección nuestra]** indican valores o formas que no vienen de un artículo.

16. **Desacople ortogonal S_atm / S_obs.**
    - S_atm = P_ret(v_esc, I) × f_flare, con f_flare = 1 porque el catálogo no trae un trazador de actividad por planeta. P_ret es la línea de costa cósmica de Zahnle & Catling 2017 en su versión de insolación bolométrica, calibrada en Marte (igual que v3–v4). S_atm ya no contiene ningún término instrumental (TSM/ESM/SNR). La (Mp, Rp) de cada draw es la misma para v_esc y para H.
    - La versión con fluencia XUV acumulada de Zahnle & Catling no se usa. El estrés XUV sigue en S_estrella, como en v1–v4. Es un solapamiento parcial entre subíndices, no un doble conteo dentro de S_atm.
    - S_obs = P_obs × Availability(tier).
      - Transitantes: P_obs = 1 − (1 − Φ(SNR_trans − u))(1 − Φ(SNR_emis − u)). SNR_trans usa μ log-sigmoide, piso de ruido y nubes (v3); SNR_emis = ESM·√n_ecl. La combinación «O» es la misma que v3 usaba en S_atm **[elección nuestra]**.
      - No transitantes: Φ(SNR_reflejada − u), con ELT/PCS igual que v4.
    - Piloto α_obs 2 vs 2.5: mismo top 10 (10/10), top 20 18/20, Spearman 0.987 sobre los 613 planetas con P_alive > 0. Se mantiene α = 2 del documento.
17. **μ(Rp) log-sigmoide.**
    - w_H2 = 1/(1 + exp(−(Rp − R0)/ΔR)), con R0 ~ N(1.65, 0.08) por draw y ΔR = 0.12 (valores del spec).
    - ln μ = (1 − w) ln μ_sec + w ln 2.3 (interpolación geométrica), con μ_sec = 28 o 44 con probabilidad 1/2.
    - Override: si la masa es medida con M/σ_M > 3 y ρ > 5.5 g cm⁻³ (por draw), entonces w_H2 = 0.
    - Cita verificada: Fulton et al. 2017 (AJ 154, 109) detectan un déficit ≥ 2× en 1.5–2.0 R⊕. R0 = 1.65 cae dentro de ese intervalo, pero no es su centro (1.75); el valor es del spec y R0 entra al Sobol en [1.5, 1.8]. Van Eylen et al. 2018 no se usa.
    - Se conserva μ_k (el μ que Kempton/Louie suponen en cada intervalo del TSM), porque es la normalización del TSM y no física del planeta. El salto residual de tsm_scale × μ_k en 1.5 R⊕ es de ~15 % (3.42 vs 2.90).
    - Se eliminó el reescalado ×2.3/18 de v3–v4 (TSM entre 1.5 y 2.0), que solo existía dentro del término TSM de S_atm.
18. **Clima de enanas M sincrónicas.** Aplica si el planeta es sincrónico (regla de v1–v4: M* < 0.6 y a < 0.2 AU) y T_eff < 3800 K en el draw.
    - Borde interno, albedo: Yang, Cowan & Abbot 2013 (ApJL 771, L45), Tabla 1, grupo 3 (CAM3, 1:1, estrella de 3400 K). Albedo planetario 0.42, 0.45, 0.47, 0.49, 0.50, 0.51 y 0.53 para S0 = 1000, 1200, …, 2200 W m⁻².
      - A = A_Shields + eff·(A_Yang(S0) − A_Shields), con eff ~ U(0.5, 1) **[elección nuestra]**. Representa la reducción del efecto de nubes por rotación más rápida en M tardías que reporta Kopparapu et al. 2016 (ApJ 819, 84).
      - Solo dentro del rango tabulado (S0 ≥ 1000 W m⁻²); por encima de 2200 se mantiene 0.53. No se extrapola por debajo de 1000 **[elección nuestra]**. Esto deja un salto en S = 0.735 S⊕.
      - El albedo nuevo también entra en T_eq, H y el contraste en luz reflejada.
    - Borde interno, prior: los bordes rg y rv de Kopparapu+2014 se multiplican por k_in ~ U(1.5, 2.0). Yang+2013 dice «habitable at twice the stellar flux found by previous studies»; Kopparapu+2016 dice «up to ~25% lower» que los GCM previos para M tardías, y 2 × 0.75 = 1.5. Intervalo uniforme **[elección nuestra]**.
    - Borde externo (S_eff < 0.8): si la presión de CO2 disponible es P_CO2 ≥ P_req(S), entonces fHZ = max(fHZ_EBM, f_eye).
      - Anclas de Turbet et al. 2018 (A&A 612, A86), para TRAPPIST-1 y planetas sincrónicos:
        - S ≥ 0.662: agua líquida subestelar «whatever its atmosphere», si el planeta es rico en agua. Proxima b, ~890 W m⁻², es análogo.
        - S = 0.382 (TRAPPIST-1f): CO2 > ~1 bar; con menos, la atmósfera colapsa.
        - S = 0.258 (TRAPPIST-1g): CO2 > ~5 bar.
        - Borde externo con 70 bar: S = 0.225.
      - Interpolación lineal en P entre 0.382 y 0.662 y log-lineal entre 0.225 y 0.382 **[elección nuestra]**.
      - Prior de presión log10(P_CO2/bar) ~ U(−3, 2) **[elección nuestra]**.
      - Fracción habitable f_eye ~ U(0.15, 0.45) **[elección nuestra]**.
      - Supuesto implícito: el planeta es rico en agua (Turbet lo exige) y las anclas de TRAPPIST-1 (2559 K) se aplican a toda M < 3800 K.
    - Los números aleatorios nuevos salen de un generador con semilla crc32(nombre) + 1, así la secuencia principal no cambia.
19. **Hurdle y ε.** UHI = Π S_i^{w_i} con piso ε = 1e-12 antes del log. Si algún S_i ≤ ε, UHI = 0 exacto (veto).
    - Como Φ(SNR − 5) ≥ Φ(−5) ≈ 2.9e-7 > ε, la parte instrumental nunca veta: P_alive mide la viabilidad física (S_fis, S_clima, S_geo…).
    - Ranking conjunto por draw: en cada draw k se ordenan todos los planetas con UHI_k > 0; los planetas con UHI_k = 0 reciben rango = (vivos en k) + 1.
    - Los draws de planetas distintos (incluso del mismo sistema) son independientes **[elección nuestra / limitación]**.
20. **Controles.** Tierra, Venus y Marte sintéticos alrededor de un análogo solar a 10 pc, en geometría transitante y no transitante. Magnitudes del Sol: M_J = 3.67 y M_Ks = 3.27 (Willmer 2018, ApJS 236, 47). Errores de 1 % en R y M **[elección nuestra]**. Script: `anchors_v5.py`.
21. **Estado de los controles (18 nov 2026) y diagnóstico de Venus.** Resultados v5 (`anchors_v5.csv`): Tierra sintética percentil 99.21 (no transita) / 99.12 (transita) → pasa; Marte UHI = 0 (S_clima = 0: fHZ = 0 con S = 0.43; P_ret ≈ 0.49) → pasa; **Venus percentil 98.74 / 98.58 → NO pasa**.
    - Venus sí recibe castigo climático, pero no veto: S_clima mediana 0.015 (fHZ del EBM 0.15 × prior Kopparapu fuera de la HZ optimista 0.1, porque S = 1.91 > rv = 1.776 para 5772 K).
    - Causa 1, EBM: OLR lineal A + B·T sin límite de radiación del vapor de agua. Con A ~ U(150, 203.3) y albedo ~0.34, los polos de Venus quedan entre 0 y 100 °C en parte de los draws, y nada dispara el invernadero desbocado.
    - Causa 2, agregación: UHI = Π S_i^{w_i} con pesos de Dirichlet (α_clima = 2 de 12). Un S_clima de 0.015 baja el UHI solo por un factor ~0.015^{w}, y no veta porque 0.015 > ε. El prior de Kopparapu es blando (0.1) y no un límite físico.
    - Causa 3, contexto del control: a 10 pc S_obs ≈ 1e-6 para Tierra y Venus. Φ(SNR−5) nunca veta (≥ 2.9e-7 > ε), así que ese factor domina a ambos por igual; y 5499 de los 6112 planetas del catálogo tienen UHI = 0, de modo que el percentil 99 equivale a quedar ~#61 entre 613 vivos. Tierra queda #49, Venus #78; la razón Venus/Tierra en UHI_rank es 0.57 (nt) y 0.48 (tr).
    - No es un error de geometría: el refuerzo de nubes de Yang (punto 18) no se aplica, porque la estrella es G (5772 K) y no es sincrónica; el albedo mediano de Venus es 0.344, igual al de la Tierra.
    - **Corrección propuesta (NO aplicada al catálogo; prueba local `venus_fix_test.py`, mismos draws y semilla):** límite de invernadero desbocado como veto de fHZ cuando el flujo absorbido S·S0/4·(1−A) > 282 W m⁻² (límite de radiación saliente del vapor de agua; Goldblatt et al. 2013, Nature Geosci 6, 661; cifra verificada el 18 nov 2026, ver punto 23). Variante "híbrida": solo para planetas no sincrónicos o T_eff ≥ 3800 K; las M sincrónicas conservan el borde interno 3D de v5 (Yang+2013 × k_in), porque el límite 1-D con redistribución /4 contradice ese resultado.
      - Híbrida: Venus → UHI = 0 (P_alive 0); Tierra percentil 99.30/99.26; vivos 613 → 515; top 10 idéntico; el primer cambio fuerte es Kepler-61 b (#17 → #50).
      - Global (también M sincrónicas): Venus 0; vivos 463; cambia el top 10 (TOI-4336 A b #1 → #14, Proxima d #7 → #15, TOI-2094 b y K2-125 b salen). Se descarta por chocar con Yang+2013.
      - Kopparapu "recent Venus" como veto (p_HZ = 0 si S > rv·k_in): Venus P_alive 0.0006, top 10 idéntico, vivos 572. Es circular, porque rv se calibró con Venus.
      - Decidido por Effy (18 nov 2026): variante híbrida adoptada en v5.1, ver punto 23.
22. **Sobol post-v5:** terminado (punto 5 del spec de Neon). Resumen en `thresholds_for_effy.md` (sección v5) y `sobol_total_v5.png`.
23. **v5.1: veto híbrido de invernadero desbocado (decisión de Effy, 18 nov 2026).**
    - Valor verificado: límite de radiación térmica saliente de una atmósfera húmeda **F_rg = 282 W m⁻²**. Fuente: Goldblatt, C., Robinson, T. D., Zahnle, K. J. & Crisp, D. (2013), «Low simulated radiation limit for runaway greenhouse climates», *Nature Geoscience* 6, 661–667, doi:10.1038/ngeo1892. Abstract (nature.com, VPL/UW, repositorio UVic): «We find a thermal radiation limit of 282 W m−2 (lower than previously reported)».
    - Corrigendum (Nature Geoscience, publicado nov/dic 2014): corrige la radiación solar *absorbida* del estado desbocado (294 → 286 W m⁻²) y Fig. 4; **no cambia el límite de 282 W m⁻²**, que es el que usamos.
    - Regla: si el flujo absorbido S·S0/4·(1 − A) > 282 W m⁻² (S0 = 1361 W m⁻², redistribución global /4, A = albedo del draw), entonces fHZ = 0 → S_clima = 0 → UHI = 0 (veto exacto). El criterio del paper es justamente «absorbido > límite térmico ⇒ desbocado».
    - Alcance híbrido: solo planetas **no sincrónicos o con T_eff ≥ 3800 K**. Las M sincrónicas (T_eff < 3800 K) conservan el borde interno 3D de v5 (albedo de Yang+2013 × k_in), porque el límite 1-D con redistribución /4 contradice los GCM de rotación 1:1 **[elección nuestra, aprobada por Effy]**.
    - Interruptor en código: `V3['F_rg'] = 282.0`, `V3['rg_veto'] = True` (en `climate_v5`). No consume números aleatorios nuevos (mismos draws y semilla que v5).
    - Limitaciones: es un límite 1-D de cielo despejado y atmósfera saturada; las nubes y la subsaturación (que el propio paper señala como necesarias para evitar el desbocado en la Tierra) se representan solo a través del albedo del draw. No se calibró con Venus.
    - Resultado: Venus sintético UHI = 0 (P_alive 0) en ambas geometrías → el control pasa. Detalle en `CHANGELOG_v4_v5.md` (sección v5.1).

24. **v5.2-final (spec de cierre de Neon/R-Aetherios; versión final aprobada por Effy, 18 nov 2026).**
    - **Rampa de Yang+2013.** En v5.1 el albedo de nubes subestelares se aplicaba con un escalón: `alb = where(md & (s0 >= 1000), alb_sh + eff·(A_Yang − alb_sh), alb_sh)` en `climate_v5`. En v5.2 es `alb = alb_sh + w_eff·eff·(A_Yang(S0) − alb_sh)` con w_eff = clip((S0 − 800)/200, 0, 1). Entre 800 y 1000 W m⁻², A_Yang queda en el primer valor tabulado (0.42) y no se extrapola. Parámetros: `V3['yang_ramp_lo'] = 800`, `V3['yang_ramp_hi'] = 1000`.
    - **Frontera Teff del veto separable.** El veto usa `V3['teff_rg']`, que por defecto es igual a `teff_md` (3800 K). Sirve para el barrido 1-D; los resultados de la base son idénticos.
    - **Barrido 1-D del veto** (`sweep_v52.py`, una variable a la vez, 10 000 draws, semilla 42). F_rg ∈ {275, 280, 282, 286, 290, 295} W m⁻² pasa en todos los valores (incluye 282 W m⁻² de Goldblatt et al. 2013 y 286 W m⁻² del corrigendum 2014, que es el solar absorbido y no el límite térmico). T_rg ∈ {3700, 3750, 3800, 3850, 3900} K **falla en 3700 y 3750 K**. Ver `tables/tables_v52.md` §5 y §7.
    - **Limitación aceptada, sin ajuste:** en T_rg = 3700 K, K2-125 b (T_eff = 3654 K) pasa de #10 (UHI_rank 0.333, P_alive 0.841) a #18 (0.265, 0.652); en T_rg = 3750 K queda #13 (0.306, 0.766). En ambos casos entra Ross 128 b (UHI_rank 0.323 en la base). El empate en la frontera por defecto (3800 K) es marginal: 0.333 frente a 0.323. Los otros 9 del top 10 no cambian en el barrido primario de T_rg. Venus sintético mantiene P_alive = 0 y Spearman > 0.99 en todos los valores primarios. **No se retocó T_rg** (sigue en 3800 K = teff_md). Es una limitación de la frontera térmica del veto, no un fallo de las anclas.
    - **ESI corregido** (ver C). La errata (exponente total w/8 en vez de w/4) inflaba los valores de ESI de los documentos v4 y v5; el UHI no cambia porque el ESI no entra en él.
    - **Sobol:** el análisis post-v5 (51 supuestos) no incluye F_rg (añadido en v5.1) ni la rampa w_eff (añadida en v5.2). No se volvió a correr, porque el top 10 no cambió de v5 a v5.2.
    - **Licencias (decisión de Effy):** código MIT (© Effy Viridis, 2026); datos, resultados y figuras CC-BY 4.0. El catálogo PHL no se redistribuye (se descarga al correr, bajo sus propios términos). El DOI de Zenodo queda pendiente.
    - Esta es la **versión final** del modelo (`v5.2-final`).
