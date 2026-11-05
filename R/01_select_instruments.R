# Select independent instruments for serum 25OHD from genome-wide significant variants.
# Independence is approximated by greedy positional pruning (see README, Limitations).

suppressMessages(library(data.table))

window_bp <- 500000L
expo <- fread("data/exposure_p5e-8.tsv")
expo <- expo[chr %in% as.character(1:22)]
expo[, chr := as.integer(chr)]
expo[, f_stat := (beta / se)^2]

# Greedy pruning: take the smallest p-value, drop every variant within the window, repeat.
setorder(expo, pval)
kept <- list()
remaining <- copy(expo)
while (nrow(remaining) > 0) {
  lead <- remaining[1]
  kept[[length(kept) + 1L]] <- lead
  remaining <- remaining[!(chr == lead$chr & abs(pos - lead$pos) <= window_bp)]
}
instruments <- rbindlist(kept)
setorder(instruments, chr, pos)

# Drop weak instruments (F < 10) and variants with no allele frequency (needed for palindromes).
instruments <- instruments[f_stat >= 10 & !is.na(eaf)]

fwrite(instruments, "data/instruments.tsv", sep = "\t")
cat("significant variants:", nrow(expo), "\n")
cat("instruments after pruning:", nrow(instruments), "\n")
cat("median F statistic:", round(median(instruments$f_stat), 1), "\n")
