#!/usr/bin/env bash
# Download the harmonised type 1 diabetes summary statistics (GCST90014023,
# Chiou 2021, 2.5 GB) once, then keep only the rows for the selected instruments.
# The download resumes if interrupted; the raw file stays in data/raw/ (not tracked).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw
URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/harmonised/34012112-GCST90014023-EFO_0001359.h.tsv.gz"
RAW="data/raw/GCST90014023.h.tsv.gz"
IDS="data/instruments.tsv"
OUT="data/outcome_instruments.tsv"
curl -fSL -C - -o "$RAW" "$URL"
# Pass 1 reads the instrument list (column 1, header skipped); pass 2 filters the stream.
gzip -dc "$RAW" | awk -F'\t' -v OFS='\t' '
  FILENAME == ARGV[1] { if (FNR > 1) want[$1] = 1; next }
  FNR == 1 { print "SNP","chr","pos","effect_allele","other_allele","beta","se","eaf","pval"; next }
  $2 in want && $7 != "NA" { print $2, $3, $4, $6, $5, $7, $21, $11, $14 }
' "$IDS" - > "$OUT"
echo "outcome rows kept: $(($(wc -l < "$OUT") - 1))"
