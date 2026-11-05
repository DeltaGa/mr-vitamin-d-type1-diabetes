#!/usr/bin/env bash
# Stream the harmonised 25OHD summary statistics (GCST90000618) and keep only
# genome-wide significant variants. The 460 MB file is never stored in full.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data
URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000618/harmonised/32242144-GCST90000618-EFO_0004631.h.tsv.gz"
OUT="data/exposure_p5e-8.tsv"
# Columns kept: rsid, chr, pos, effect allele, other allele, beta, se, eaf, p, n
curl -fsSL "$URL" | gzip -dc | awk -F'\t' -v OFS='\t' '
  NR == 1 { print "SNP","chr","pos","effect_allele","other_allele","beta","se","eaf","pval","n"; next }
  $2 ~ /^rs[0-9]+$/ && $7 != "NA" && $19 != "NA" && $20 != "NA" && ($20 + 0) < 5e-8 {
    print $2, $3, $4, $6, $5, $7, $19, $11, $20, $16
  }' > "$OUT"
echo "variants kept: $(($(wc -l < "$OUT") - 1))"
