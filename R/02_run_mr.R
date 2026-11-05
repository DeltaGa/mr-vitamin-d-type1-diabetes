# Two-sample Mendelian randomization: serum 25OHD -> type 1 diabetes.
# Exposure: GCST90000618 (n = 496,946, European). Outcome: GCST90014023 (Chiou 2021).

suppressMessages({
  library(TwoSampleMR)
  library(ggplot2)
  library(data.table)
})

dir.create("results", showWarnings = FALSE)

expo_raw <- fread("data/instruments.tsv")
out_raw <- fread("data/outcome_instruments.tsv")

exposure <- format_data(
  as.data.frame(expo_raw), type = "exposure",
  snp_col = "SNP", beta_col = "beta", se_col = "se",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  eaf_col = "eaf", pval_col = "pval", samplesize_col = "n", chr_col = "chr", pos_col = "pos"
)
exposure$exposure <- "Serum 25OHD (SD)"

outcome <- format_data(
  as.data.frame(out_raw), type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  eaf_col = "eaf", pval_col = "pval"
)
outcome$outcome <- "Type 1 diabetes (log OR)"

dat <- harmonise_data(exposure, outcome, action = 2)
dat <- dat[dat$mr_keep, ]
cat("instruments after harmonisation:", nrow(dat), "\n")

# Main estimates: IVW (multiplicative random effects), weighted median, MR-Egger.
res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
res$or <- exp(res$b)
res$or_lci95 <- exp(res$b - 1.96 * res$se)
res$or_uci95 <- exp(res$b + 1.96 * res$se)
fwrite(res, "results/mr_estimates.csv")

het <- mr_heterogeneity(dat)
ple <- mr_pleiotropy_test(dat)
fwrite(het, "results/heterogeneity.csv")
fwrite(ple, "results/egger_intercept.csv")

# Leave-one-out and single-SNP analyses.
loo <- mr_leaveoneout(dat)
single <- mr_singlesnp(dat)
fwrite(loo, "results/leave_one_out.csv")
fwrite(single, "results/single_snp.csv")

ggsave("results/scatter.png", mr_scatter_plot(res, dat)[[1]], width = 7, height = 6, dpi = 150)
ggsave("results/forest.png", mr_forest_plot(single)[[1]], width = 7, height = 9, dpi = 150)
ggsave("results/leave_one_out.png", mr_leaveoneout_plot(loo)[[1]], width = 7, height = 9, dpi = 150)
ggsave("results/funnel.png", mr_funnel_plot(single)[[1]], width = 7, height = 6, dpi = 150)

print(res[, c("method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95")])
print(het[, c("method", "Q", "Q_df", "Q_pval")])
print(ple[, c("egger_intercept", "se", "pval")])
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
