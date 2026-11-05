# Sensitivity analysis without the HLA region (chr 6, 25-34 Mb, covers GRCh37 and GRCh38).
# HLA is the dominant type 1 diabetes locus; a 25OHD variant there tags HLA risk, not vitamin D.
# Same inputs and estimators as 02_run_mr.R.

suppressMessages({
  library(TwoSampleMR)
  library(ggplot2)
  library(data.table)
})

set.seed(2026)  # weighted median uses a bootstrap; fix the seed for reproducibility
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

hla <- dat$chr.exposure == "6" & dat$pos.exposure >= 25e6 & dat$pos.exposure <= 34e6
cat("instruments in HLA region removed:", sum(hla), "\n")
print(dat[hla, c("SNP", "chr.exposure", "pos.exposure")])
dat2 <- dat[!hla, ]

res <- mr(dat2, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
res$or <- exp(res$b)
res$or_lci95 <- exp(res$b - 1.96 * res$se)
res$or_uci95 <- exp(res$b + 1.96 * res$se)
fwrite(res, "results/mr_estimates_no_hla.csv")
het <- mr_heterogeneity(dat2)
ple <- mr_pleiotropy_test(dat2)
fwrite(het, "results/heterogeneity_no_hla.csv")
fwrite(ple, "results/egger_intercept_no_hla.csv")
single <- mr_singlesnp(dat2)
loo <- mr_leaveoneout(dat2)
fwrite(loo, "results/leave_one_out_no_hla.csv")
fwrite(single, "results/single_snp_no_hla.csv")

ggsave("results/scatter_no_hla.png", mr_scatter_plot(res, dat2)[[1]], width = 7, height = 6, dpi = 150)
ggsave("results/forest_no_hla.png", mr_forest_plot(single)[[1]], width = 8, height = 26, dpi = 150)
ggsave("results/funnel_no_hla.png", mr_funnel_plot(single)[[1]], width = 7, height = 6, dpi = 150)

print(res[, c("method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95")])
print(het[, c("method", "Q", "Q_df", "Q_pval")])
print(ple[, c("egger_intercept", "se", "pval")])
