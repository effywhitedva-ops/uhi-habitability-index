# Sensibilidad al prior de clase - combined, 3 sigma

## Global (vs principal)

| escenario | top10_overlap | top20_overlap | kendall_top20_union | spearman_top20_union | spearman_vivos | n_vivos | nuevo_n1 |
|---|---|---|---|---|---|---|---|
| A1 prior rocoso (π=2:1:1, M_v 5-10) | 6 | 14 | 0.317 | 0.400 | 0.979 | 289 | GJ 273 b |
| A2 prior volátil (π=0.5:1:1, M_v 1.5-3) | 4 | 13 | 0.117 | 0.071 | 0.948 | 281 | Kepler-442 b |
| A3 transición CK17 (M_v=2.04) | 4 | 13 | 0.225 | 0.146 | 0.958 | 288 | Kepler-442 b |
| A4 empírico archivo M (M_v=4.10, s=0.40) | 8 | 19 | 0.905 | 0.973 | 0.999 | 288 | GJ 273 b |
| A5 empírico M, P>5 d (M_v=2.83, s=0.70) | 7 | 17 | 0.573 | 0.698 | 0.997 | 290 | Kepler-442 b |
| B1 M = M sin i | 9 | 16 | 0.609 | 0.746 | 0.994 | 291 | GJ 273 b |
| B2 M sin i + prior Howard+2010 | 9 | 18 | 0.844 | 0.939 | 0.997 | 290 | GJ 273 b |
| C1 M-R Chen & Kipping 2017 | 8 | 20 | 0.884 | 0.965 | 0.964 | 306 | GJ 273 b |
| D1 frontera línea Earth-like | 9 | 17 | 0.787 | 0.848 | 0.962 | 282 | GJ 273 b |
| D2 frontera 100% H2O | 10 | 20 | 0.989 | 0.998 | 0.985 | 288 | GJ 273 b |

## Top 20 principal

| planet | evidencia | mass_source | rank_main | rank_min | rank_max | peor_escenario | mejor_escenario | UHI_main | UHI_min | UHI_max | P_R_main | P_R_min | P_R_max | P_top10_main | P_top10_min | P_top10_max | veredicto |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GJ 273 b | prior | msini | 1 | 1 | 30 | A2 | A1 | 0.287 | 0.056 | 0.423 | 0.609 | 0.118 | 0.896 | 0.489 | 0.105 | 0.633 | DEPENDE DEL PRIOR |
| Wolf 1061 c | prior | msini | 2 | 2 | 44 | A2 | A1 | 0.249 | 0.039 | 0.410 | 0.507 | 0.080 | 0.843 | 0.434 | 0.073 | 0.629 | DEPENDE DEL PRIOR |
| Kepler-442 b | radio | mr | 3 | 1 | 8 | A1 | A2 | 0.230 | 0.213 | 0.240 | 0.911 | 0.846 | 0.951 | 0.374 | 0.272 | 0.499 | robusto (top 10 en todos) |
| GJ 1061 d | prior | measured | 4 | 2 | 10 | A1 | A3 | 0.226 | 0.133 | 0.235 | 0.954 | 0.561 | 0.995 | 0.399 | 0.271 | 0.423 | robusto (top 10 en todos) |
| Proxima Cen b | prior | msini | 5 | 3 | 11 | A1 | A2 | 0.202 | 0.140 | 0.215 | 0.938 | 0.661 | 0.990 | 0.371 | 0.247 | 0.398 | estable (top 20 en todos) |
| GJ 3378 b | prior | msini | 6 | 3 | 27 | A2 | B1 | 0.199 | 0.061 | 0.258 | 0.718 | 0.221 | 0.926 | 0.316 | 0.120 | 0.337 | DEPENDE DEL PRIOR |
| Ross 128 b | prior | msini | 7 | 7 | 14 | A1 | A4 | 0.165 | 0.093 | 0.178 | 0.892 | 0.505 | 0.975 | 0.240 | 0.142 | 0.273 | estable (top 20 en todos) |
| GJ 251 b | prior | msini | 8 | 4 | 62 | A2 | A1 | 0.147 | 0.017 | 0.273 | 0.432 | 0.052 | 0.811 | 0.242 | 0.035 | 0.341 | DEPENDE DEL PRIOR |
| L 98-59 f | prior | msini | 9 | 8 | 50 | A3 | B1 | 0.146 | 0.031 | 0.208 | 0.624 | 0.130 | 0.894 | 0.233 | 0.071 | 0.253 | DEPENDE DEL PRIOR |
| TRAPPIST-1 e | densidad | measured | 10 | 2 | 89 | D1 | A2 | 0.145 | 0.011 | 0.145 | 0.996 | 0.079 | 0.997 | 0.158 | 0.014 | 0.353 | DEPENDE DEL PRIOR |
| Kepler-452 b | radio | mr | 11 | 5 | 18 | A1 | A3 | 0.139 | 0.119 | 0.153 | 0.731 | 0.621 | 0.812 | 0.223 | 0.175 | 0.299 | DEPENDE DEL PRIOR |
| GJ 682 b | prior | msini | 12 | 9 | 44 | A3 | A1 | 0.137 | 0.039 | 0.239 | 0.394 | 0.149 | 0.690 | 0.241 | 0.079 | 0.364 | DEPENDE DEL PRIOR |
| Kepler-1544 b | radio | mr | 13 | 6 | 18 | B1 | A3 | 0.135 | 0.092 | 0.178 | 0.459 | 0.311 | 0.616 | 0.218 | 0.193 | 0.282 | DEPENDE DEL PRIOR |
| Teegarden's Star c | prior | msini | 14 | 13 | 23 | A1 | C1 | 0.128 | 0.088 | 0.135 | 0.937 | 0.667 | 0.991 | 0.196 | 0.120 | 0.214 | DEPENDE DEL PRIOR |
| Gliese 12 b | densidad | measured | 15 | 5 | 26 | A1 | A2 | 0.125 | 0.102 | 0.128 | 0.884 | 0.678 | 0.911 | 0.082 | 0.030 | 0.265 | DEPENDE DEL PRIOR |
| TOI-1452 b | densidad | measured | 16 | 8 | 71 | D1 | A2 | 0.121 | 0.024 | 0.137 | 0.579 | 0.116 | 0.657 | 0.131 | 0.032 | 0.274 | DEPENDE DEL PRIOR |
| TOI-700 d | radio | mr | 17 | 7 | 29 | A1 | A2 | 0.118 | 0.114 | 0.120 | 0.958 | 0.921 | 0.977 | 0.109 | 0.049 | 0.265 | DEPENDE DEL PRIOR |
| K2-72 e | radio | mr | 18 | 10 | 27 | A1 | A3 | 0.115 | 0.097 | 0.127 | 0.769 | 0.652 | 0.857 | 0.124 | 0.061 | 0.235 | DEPENDE DEL PRIOR |
| GJ 1002 b | prior | msini | 19 | 16 | 28 | A1 | D1 | 0.114 | 0.078 | 0.122 | 0.933 | 0.653 | 0.989 | 0.130 | 0.070 | 0.189 | DEPENDE DEL PRIOR |
| GJ 433 d | prior | msini | 20 | 3 | 78 | A2 | A1 | 0.111 | 0.010 | 0.275 | 0.261 | 0.025 | 0.675 | 0.178 | 0.018 | 0.408 | DEPENDE DEL PRIOR |
