#!/usr/bin/env bash
# Run the whole analysis from raw public summary statistics to figures.
# Needs curl, gzip, awk and Rscript on PATH.
set -euo pipefail
cd "$(dirname "$0")"
bash scripts/fetch_exposure.sh
Rscript R/01_select_instruments.R
bash scripts/fetch_outcome.sh
Rscript R/02_run_mr.R
