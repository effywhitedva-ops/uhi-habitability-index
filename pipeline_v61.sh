#!/bin/bash
# Unattended v6.1 chain: MC -> prior-sensitivity tables -> v5.2 compatibility test -> Sobol N=512 -> report.
# Run from anywhere; it changes to the repository root. Requires a virtual environment at ./.venv with
# requirements.txt installed (see README, "Reproduce").
# Expected runtime is dominated by the Sobol step (many hours to days). Logs go to logs/v6.1/; that folder must exist.
# Outputs: results/v6.1, tables/v6.1, figures/v6.1, sens_out_v6_{5,3}s_n512 and the report tables/figures.
cd "$(dirname "$0")"
PY="./.venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "UHI virtual environment not found at $PY. Create it with: python3.13 -m venv .venv && .venv/bin/pip install -r requirements.txt" >&2
  exit 1
fi
./run_v61_mc.sh
UHI_RES=results/v6.1 UHI_TAB=tables/v6.1 UHI_FIG=figures/v6.1 "$PY" src/prior_sensitivity_v6.py > logs/v6.1/prior_sens_v61.log 2>&1; echo "prior_sens rc=$?"
"$PY" tests/test_v6.py > logs/v6.1/tests_v6_compat.log 2>&1; echo "test_v6 (compat v5.2) rc=$?"
echo SOBOL_START $(date +%H:%M)
./run_v61_sobol.sh
