#!/bin/bash
# v6.1: Sobol N=512 a 5s y 3s (4 procesos cada uno) en directorios nuevos; luego reporte S_T corregido y marcador DONE
cd "$(dirname "$0")"
PY="./.venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "No está el venv de UHI en $PY. Recréalo con: bash /workspace/env/bootstrap.sh" >&2
  exit 1
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
"$PY" src/sens_v6.py 5 512 4 sens_out_v6_5s_n512 > logs/v6.1/sobol_5s_n512.log 2>&1 &
P5=$!
"$PY" src/sens_v6.py 3 512 4 sens_out_v6_3s_n512 > logs/v6.1/sobol_3s_n512.log 2>&1 &
P3=$!
echo "sobol lanzado 5s pid=$P5 3s pid=$P3 $(date '+%Y-%m-%d %H:%M')"
wait $P5; R5=$?; wait $P3; R3=$?
echo "sobol fin rc5=$R5 rc3=$R3 $(date '+%Y-%m-%d %H:%M')"
"$PY" src/sobol_report_v6.py sens_out_v6_5s_n512 sens_out_v6_3s_n512 v6.1_n512 > logs/v6.1/sobol_report_n512.log 2>&1
echo "reporte rc=$? $(date '+%Y-%m-%d %H:%M')"
echo "UHI_V61_SOBOL_N512_DONE $(date '+%Y-%m-%d %H:%M')"
