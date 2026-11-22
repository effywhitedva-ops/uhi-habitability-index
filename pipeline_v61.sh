#!/bin/bash
# Cadena v6.1 desatendida: MC -> tablas de prior -> test de compatibilidad -> Sobol N=512 -> reporte -> DONE (logs/v6.1/pipeline.log)
cd "$(dirname "$0")"
PY="./.venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "No está el venv de UHI en $PY. Recréalo con: bash /workspace/env/bootstrap.sh" >&2
  exit 1
fi
./run_v61_mc.sh
UHI_RES=results/v6.1 UHI_TAB=tables/v6.1 UHI_FIG=figures/v6.1 "$PY" src/prior_sensitivity_v6.py > logs/v6.1/prior_sens_v61.log 2>&1; echo "prior_sens rc=$?"
"$PY" tests/test_v6.py > logs/v6.1/tests_v6_compat.log 2>&1; echo "test_v6 (compat v5.2) rc=$?"
echo SOBOL_START $(date +%H:%M)
./run_v61_sobol.sh
