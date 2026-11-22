# Sensibilidad al prior de clase - rocky, 5 sigma

## Global (vs principal)

| escenario | top10_overlap | top20_overlap | kendall_top20_union | spearman_top20_union | spearman_vivos | n_vivos | nuevo_n1 |
|---|---|---|---|---|---|---|---|
| A1 prior rocoso (π=2:1:1, M_v 5-10) | 8 | 16 | 0.522 | 0.719 | 0.982 | 262 | GJ 273 b |
| A2 prior volátil (π=0.5:1:1, M_v 1.5-3) | 4 | 13 | 0.214 | 0.256 | 0.947 | 240 | Proxima Cen b |
| A3 transición CK17 (M_v=2.04) | 5 | 13 | 0.288 | 0.346 | 0.959 | 249 | Proxima Cen b |
| A4 empírico archivo M (M_v=4.10, s=0.40) | 9 | 19 | 0.952 | 0.992 | 0.998 | 256 | GJ 273 b |
| A5 empírico M, P>5 d (M_v=2.83, s=0.70) | 9 | 18 | 0.835 | 0.941 | 0.997 | 255 | GJ 273 b |
| B1 M = M sin i | 9 | 18 | 0.766 | 0.913 | 0.993 | 257 | GJ 273 b |
| B2 M sin i + prior Howard+2010 | 10 | 19 | 0.933 | 0.986 | 0.993 | 253 | GJ 273 b |
| C1 M-R Chen & Kipping 2017 | 10 | 20 | 0.958 | 0.992 | 0.962 | 274 | GJ 273 b |
| D1 frontera línea Earth-like | 10 | 18 | 0.913 | 0.958 | 0.934 | 248 | GJ 273 b |
| D2 frontera 100% H2O | 10 | 20 | 1.000 | 1.000 | 0.983 | 256 | GJ 273 b |

## Top 20 principal

| planet | evidencia | mass_source | rank_main | rank_min | rank_max | peor_escenario | mejor_escenario | UHI_main | UHI_min | UHI_max | P_R_main | P_R_min | P_R_max | P_top10_main | P_top10_min | P_top10_max | veredicto |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GJ 273 b | prior | msini | 1 | 1 | 9 | A2 | A1 | 0.287 | 0.056 | 0.422 | 0.609 | 0.118 | 0.896 | 0.527 | 0.107 | 0.720 | robusto (top 10 en todos) |
| Wolf 1061 c | prior | msini | 2 | 2 | 14 | A2 | A1 | 0.249 | 0.039 | 0.410 | 0.507 | 0.080 | 0.843 | 0.473 | 0.075 | 0.725 | DEPENDE DEL PRIOR |
| GJ 1061 d | prior | measured | 3 | 2 | 8 | A1 | A2 | 0.222 | 0.131 | 0.230 | 0.954 | 0.561 | 0.995 | 0.491 | 0.311 | 0.499 | robusto (top 10 en todos) |
| Proxima Cen b | prior | msini | 4 | 1 | 10 | A1 | A2 | 0.202 | 0.139 | 0.215 | 0.938 | 0.661 | 0.990 | 0.505 | 0.387 | 0.519 | robusto (top 10 en todos) |
| GJ 3378 b | prior | msini | 5 | 5 | 11 | A2 | B2 | 0.164 | 0.049 | 0.217 | 0.718 | 0.221 | 0.926 | 0.336 | 0.122 | 0.372 | estable (top 20 en todos) |
| Ross 128 b | prior | msini | 6 | 3 | 11 | A1 | A2 | 0.163 | 0.093 | 0.177 | 0.892 | 0.505 | 0.975 | 0.420 | 0.282 | 0.447 | estable (top 20 en todos) |
| GJ 251 b | prior | msini | 7 | 3 | 36 | A2 | A1 | 0.146 | 0.017 | 0.271 | 0.432 | 0.052 | 0.811 | 0.298 | 0.039 | 0.456 | DEPENDE DEL PRIOR |
| GJ 682 b | prior | msini | 8 | 7 | 15 | A2 | A1 | 0.137 | 0.038 | 0.239 | 0.394 | 0.149 | 0.690 | 0.269 | 0.087 | 0.426 | estable (top 20 en todos) |
| GJ 433 d | prior | msini | 9 | 4 | 51 | A2 | A1 | 0.105 | 0.009 | 0.260 | 0.261 | 0.025 | 0.675 | 0.180 | 0.018 | 0.413 | DEPENDE DEL PRIOR |
| Gl 725 B c | prior | msini | 10 | 10 | 38 | A3 | B2 | 0.099 | 0.017 | 0.165 | 0.512 | 0.091 | 0.849 | 0.185 | 0.034 | 0.292 | DEPENDE DEL PRIOR |
| GJ 887 d | prior | msini | 11 | 5 | 49 | A2 | A1 | 0.092 | 0.010 | 0.256 | 0.199 | 0.021 | 0.593 | 0.158 | 0.017 | 0.413 | DEPENDE DEL PRIOR |
| Teegarden's Star c | prior | msini | 12 | 6 | 18 | A1 | A2 | 0.092 | 0.063 | 0.097 | 0.937 | 0.667 | 0.991 | 0.261 | 0.189 | 0.278 | DEPENDE DEL PRIOR |
| GJ 625 b | prior | msini | 13 | 12 | 34 | A2 | B2 | 0.089 | 0.018 | 0.127 | 0.623 | 0.140 | 0.890 | 0.217 | 0.061 | 0.255 | DEPENDE DEL PRIOR |
| TRAPPIST-1 e | densidad | measured | 14 | 4 | 77 | D1 | A2 | 0.082 | 0.006 | 0.082 | 0.996 | 0.079 | 0.997 | 0.217 | 0.018 | 0.493 | DEPENDE DEL PRIOR |
| GJ 367 d | prior | msini | 15 | 6 | 61 | A2 | A1 | 0.078 | 0.006 | 0.246 | 0.175 | 0.012 | 0.589 | 0.134 | 0.010 | 0.397 | DEPENDE DEL PRIOR |
| Proxima Cen d | prior | msini | 16 | 5 | 22 | A1 | A2 | 0.074 | 0.058 | 0.075 | 0.997 | 0.943 | 1.000 | 0.215 | 0.090 | 0.440 | DEPENDE DEL PRIOR |
| Kepler-442 b | radio | mr | 17 | 7 | 24 | A1 | A2 | 0.067 | 0.063 | 0.070 | 0.911 | 0.846 | 0.951 | 0.242 | 0.141 | 0.383 | DEPENDE DEL PRIOR |
| GJ 1002 c | prior | msini | 18 | 16 | 23 | A1 | A3 | 0.065 | 0.035 | 0.072 | 0.893 | 0.513 | 0.976 | 0.141 | 0.094 | 0.146 | DEPENDE DEL PRIOR |
| TRAPPIST-1 f | densidad | measured | 19 | 8 | 124 | D1 | A2 | 0.061 | 0.001 | 0.061 | 0.989 | 0.011 | 0.992 | 0.165 | 0.002 | 0.289 | DEPENDE DEL PRIOR |
| L 98-59 f | prior | msini | 20 | 17 | 45 | A3 | B1 | 0.057 | 0.013 | 0.082 | 0.624 | 0.130 | 0.894 | 0.221 | 0.066 | 0.247 | DEPENDE DEL PRIOR |
