#!/bin/bash
# v6.1 Monte Carlo batch: main, madhu and ten class-prior variants at 5 and 3 sigma (24 runs, seed 42,
# 10,000 draws) into results/v6.1, plus one legacy run (v6.1 fixes off) into results/v6.1_legacy_check,
# which must reproduce v6.0. Runs 8 jobs in parallel with single-threaded BLAS. Logs: logs/v6.1/mc_*.log.
cd "$(dirname "$0")"
PY="./.venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "UHI virtual environment not found at $PY. Create it with: python3.13 -m venv .venv && .venv/bin/pip install -r requirements.txt" >&2
  exit 1
fi
export PY
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
{ echo "legacy 5"; for t in 5 3; do for v in main madhu pA1_rocoso pA2_volatil pA3_ck17 pA4_empM pA5_empM_P5 pB1_msini_min pB2_msini_howard pC1_mr_ck17 pD1_linea_earthlike pD2_linea_h2o100; do echo "$v $t"; done; done; } |
  xargs -P 8 -L 1 bash -c 'if [ "$0" = legacy ]; then "$PY" run_uhi_v6.py --variant main --thr $1 --legacy --outdir results/v6.1_legacy_check > logs/v6.1/mc_legacy_main_$1s.log 2>&1; else "$PY" run_uhi_v6.py --variant $0 --thr $1 --outdir results/v6.1 > logs/v6.1/mc_$0_$1s.log 2>&1; fi; echo "fin $0 $1 rc=$? $(date +%H:%M)"'
echo MC_V61_LISTO $(date +%H:%M)
