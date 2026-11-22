# Sensibilidad al prior de clase - combined, 5 sigma

## Global (vs principal)

| escenario | top10_overlap | top20_overlap | kendall_top20_union | spearman_top20_union | spearman_vivos | n_vivos | nuevo_n1 |
|---|---|---|---|---|---|---|---|
| A1 prior rocoso (π=2:1:1, M_v 5-10) | 8 | 16 | 0.522 | 0.719 | 0.978 | 289 | GJ 273 b |
| A2 prior volátil (π=0.5:1:1, M_v 1.5-3) | 4 | 13 | 0.219 | 0.261 | 0.945 | 281 | Proxima Cen b |
| A3 transición CK17 (M_v=2.04) | 5 | 13 | 0.293 | 0.349 | 0.958 | 288 | Proxima Cen b |
| A4 empírico archivo M (M_v=4.10, s=0.40) | 10 | 19 | 0.962 | 0.994 | 0.999 | 288 | GJ 273 b |
| A5 empírico M, P>5 d (M_v=2.83, s=0.70) | 9 | 18 | 0.853 | 0.950 | 0.996 | 290 | GJ 273 b |
| B1 M = M sin i | 9 | 18 | 0.775 | 0.914 | 0.991 | 291 | GJ 273 b |
| B2 M sin i + prior Howard+2010 | 10 | 19 | 0.943 | 0.988 | 0.993 | 290 | GJ 273 b |
| C1 M-R Chen & Kipping 2017 | 10 | 20 | 0.958 | 0.992 | 0.963 | 306 | GJ 273 b |
| D1 frontera línea Earth-like | 10 | 18 | 0.913 | 0.958 | 0.965 | 282 | GJ 273 b |
| D2 frontera 100% H2O | 10 | 20 | 1.000 | 1.000 | 0.984 | 288 | GJ 273 b |

## Top 20 principal

| planet | evidencia | mass_source | rank_main | rank_min | rank_max | peor_escenario | mejor_escenario | UHI_main | UHI_min | UHI_max | P_R_main | P_R_min | P_R_max | P_top10_main | P_top10_min | P_top10_max | veredicto |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GJ 273 b | prior | msini | 1 | 1 | 9 | A2 | A1 | 0.287 | 0.056 | 0.422 | 0.609 | 0.118 | 0.896 | 0.525 | 0.107 | 0.717 | robusto (top 10 en todos) |
| Wolf 1061 c | prior | msini | 2 | 2 | 14 | A2 | A1 | 0.249 | 0.039 | 0.410 | 0.507 | 0.080 | 0.843 | 0.470 | 0.075 | 0.723 | DEPENDE DEL PRIOR |
| GJ 1061 d | prior | measured | 3 | 2 | 8 | A1 | A2 | 0.222 | 0.131 | 0.230 | 0.954 | 0.561 | 0.995 | 0.486 | 0.309 | 0.494 | robusto (top 10 en todos) |
| Proxima Cen b | prior | msini | 4 | 1 | 10 | A1 | A2 | 0.202 | 0.140 | 0.215 | 0.938 | 0.661 | 0.990 | 0.500 | 0.389 | 0.513 | robusto (top 10 en todos) |
| GJ 3378 b | prior | msini | 5 | 5 | 11 | A2 | B2 | 0.164 | 0.049 | 0.217 | 0.718 | 0.221 | 0.926 | 0.331 | 0.120 | 0.365 | estable (top 20 en todos) |
| Ross 128 b | prior | msini | 6 | 3 | 11 | A1 | A2 | 0.163 | 0.093 | 0.177 | 0.892 | 0.505 | 0.975 | 0.409 | 0.276 | 0.436 | estable (top 20 en todos) |
| GJ 251 b | prior | msini | 7 | 3 | 41 | A2 | A1 | 0.146 | 0.017 | 0.271 | 0.432 | 0.052 | 0.811 | 0.294 | 0.039 | 0.452 | DEPENDE DEL PRIOR |
| GJ 682 b | prior | msini | 8 | 7 | 15 | A2 | A1 | 0.137 | 0.038 | 0.239 | 0.394 | 0.149 | 0.690 | 0.267 | 0.086 | 0.424 | estable (top 20 en todos) |
| GJ 433 d | prior | msini | 9 | 4 | 57 | A2 | A1 | 0.105 | 0.009 | 0.260 | 0.261 | 0.025 | 0.675 | 0.178 | 0.018 | 0.412 | DEPENDE DEL PRIOR |
| Gl 725 B c | prior | msini | 10 | 10 | 34 | A2 | A4 | 0.101 | 0.020 | 0.165 | 0.512 | 0.091 | 0.849 | 0.187 | 0.041 | 0.292 | DEPENDE DEL PRIOR |
| GJ 887 d | prior | msini | 11 | 5 | 55 | A2 | A1 | 0.092 | 0.010 | 0.256 | 0.199 | 0.021 | 0.593 | 0.158 | 0.017 | 0.410 | DEPENDE DEL PRIOR |
| Teegarden's Star c | prior | msini | 12 | 6 | 18 | A1 | A2 | 0.091 | 0.063 | 0.097 | 0.937 | 0.667 | 0.991 | 0.254 | 0.186 | 0.268 | DEPENDE DEL PRIOR |
| GJ 625 b | prior | msini | 13 | 12 | 40 | A2 | B2 | 0.089 | 0.018 | 0.127 | 0.623 | 0.140 | 0.890 | 0.212 | 0.059 | 0.248 | DEPENDE DEL PRIOR |
| TRAPPIST-1 e | densidad | measured | 14 | 4 | 83 | D1 | A2 | 0.082 | 0.006 | 0.082 | 0.996 | 0.079 | 0.997 | 0.202 | 0.017 | 0.450 | DEPENDE DEL PRIOR |
| GJ 367 d | prior | msini | 15 | 6 | 68 | A2 | A1 | 0.078 | 0.006 | 0.246 | 0.175 | 0.012 | 0.589 | 0.134 | 0.010 | 0.396 | DEPENDE DEL PRIOR |
| Proxima Cen d | prior | msini | 16 | 5 | 24 | B1 | A2 | 0.074 | 0.058 | 0.075 | 0.997 | 0.943 | 1.000 | 0.203 | 0.086 | 0.407 | DEPENDE DEL PRIOR |
| Kepler-442 b | radio | mr | 17 | 7 | 24 | A1 | A2 | 0.067 | 0.063 | 0.070 | 0.911 | 0.846 | 0.951 | 0.231 | 0.138 | 0.359 | DEPENDE DEL PRIOR |
| GJ 1002 c | prior | msini | 18 | 15 | 23 | A1 | A3 | 0.066 | 0.038 | 0.072 | 0.893 | 0.513 | 0.976 | 0.140 | 0.105 | 0.146 | DEPENDE DEL PRIOR |
| TRAPPIST-1 f | densidad | measured | 19 | 8 | 146 | D1 | A2 | 0.061 | 0.001 | 0.061 | 0.989 | 0.011 | 0.992 | 0.155 | 0.002 | 0.273 | DEPENDE DEL PRIOR |
| L 98-59 f | prior | msini | 20 | 17 | 49 | A2 | B1 | 0.057 | 0.013 | 0.082 | 0.624 | 0.130 | 0.894 | 0.216 | 0.063 | 0.238 | DEPENDE DEL PRIOR |
