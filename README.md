[![R](https://img.shields.io/badge/R-4.5+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# mr-vitamin-d-type1-diabetes v1.0.0

## Two-sample Mendelian randomization of serum 25-hydroxyvitamin D on type 1 diabetes, from public GWAS summary statistics

**Version:** 1.0.0
**Release Date:** November 4, 2026
**R:** 4.5 or later (developed on 4.5.3, TwoSampleMR 0.7.12)

---

## Table of Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Architecture](#architecture)
- [Troubleshooting](#troubleshooting)
- [Limitations](#limitations)
- [Support](#support)
- [Changelog](#changelog)
- [Citations](#citations)
- [Author](#author)

---

## Overview

This repository is a learning exercise: a complete two-sample Mendelian randomization (MR) analysis written in R by a first-year nutrition student with no prior R experience. It asks whether genetically predicted serum 25-hydroxyvitamin D (25OHD) changes the risk of type 1 diabetes, a question already studied in the literature, using the same estimators a published analysis would use. It is not a research result and makes no causal claim.

### Key Capabilities

- **Tokenless data access** - both datasets are streamed from the GWAS Catalog FTP; no account or API token is needed
- **Instrument selection** - genome-wide significant variants (p < 5e-8), greedy positional pruning, F-statistic filter (F >= 10)
- **Three estimators** - inverse-variance weighted (IVW), weighted median, MR-Egger
- **Sensitivity analyses** - Cochran's Q heterogeneity, MR-Egger intercept, leave-one-out, single-variant Wald ratios
- **Figures and tables** - scatter, forest, leave-one-out and funnel plots, plus CSV tables of every estimate

---

## System Requirements

- **R:** 4.5 or later, with the packages `TwoSampleMR`, `data.table`, `ggplot2`
- **Disk:** about 2.5 GB for the type 1 diabetes download (kept in `data/raw/`, not tracked by git)
- **Memory:** 4 GB or more

### External Tools

- `bash`, `curl`, `gzip` and `awk` on `PATH` (Git Bash provides all four on Windows)
- Network access to `ftp.ebi.ac.uk`

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/DeltaGa/mr-vitamin-d-type1-diabetes.git
cd mr-vitamin-d-type1-diabetes

# 2. Install the R packages (TwoSampleMR is distributed through the MRC IEU r-universe)
Rscript -e 'options(repos = c(mrcieu = "https://mrcieu.r-universe.dev", CRAN = "https://cloud.r-project.org")); install.packages(c("TwoSampleMR", "data.table", "ggplot2"))'
```

---

## Usage

### Basic Usage

```bash
# Run the full pipeline: download, select instruments, run MR, write results/
bash run_all.sh
```

The type 1 diabetes download is large and the server is sometimes slow. The download resumes if interrupted; run the script again.

### Advanced Usage

```bash
# Run one stage at a time
bash scripts/fetch_exposure.sh        # 25OHD statistics, filtered to p < 5e-8
Rscript R/01_select_instruments.R     # pruning and F-statistic filter
bash scripts/fetch_outcome.sh         # type 1 diabetes rows for the instruments
Rscript R/02_run_mr.R                 # estimators, sensitivity analyses, figures
```

### Results

Results are written to `results/`. See the Changelog for the run that produced the files committed here.

---

## Architecture

| Component | Responsibility |
|-----------|----------------|
| `scripts/fetch_exposure.sh` | Streams the harmonised 25OHD file and keeps variants with p < 5e-8 |
| `R/01_select_instruments.R` | Prunes to one variant per 500 kb window, drops weak instruments |
| `scripts/fetch_outcome.sh` | Downloads the type 1 diabetes file once, extracts the instrument rows |
| `R/02_run_mr.R` | Harmonises alleles, runs the estimators and sensitivity analyses, writes figures |
| `run_all.sh` | Runs the four stages in order |

```
mr-vitamin-d-type1-diabetes/
+-- README.md
+-- LICENSE
+-- run_all.sh
+-- scripts/    fetch_exposure.sh, fetch_outcome.sh
+-- R/          01_select_instruments.R, 02_run_mr.R
+-- data/       instruments.tsv, outcome_instruments.tsv  (raw downloads are not tracked)
+-- results/    estimates, tables, figures
```

---

## Troubleshooting

### The type 1 diabetes download is very slow

**Cause:** the GWAS Catalog FTP serves large files at a few hundred kilobytes per second at busy times, and range requests into the file are not usable because it is a plain gzip, not an indexed one.
**Solution:** leave `scripts/fetch_outcome.sh` running; `curl -C -` resumes from the bytes already on disk.

### Rscript fails to start on Windows with a conda-installed R

**Cause:** the conda environment's `Library/bin` and `lib/R/bin/x64` folders are not on `PATH`.
**Solution:** prepend both folders to `PATH` before calling `Rscript`.

---

## Limitations

- Independence between instruments is approximated by positional pruning (500 kb windows), not by linkage disequilibrium clumping against a reference panel. Residual correlation between instruments is possible.
- The exposure and outcome studies both draw on UK Biobank participants, so sample overlap may bias the estimate toward the observational association.
- No Steiger directionality filter, MR-PRESSO or multivariable analysis is run. These are the natural next steps.
- Single analysis on European-ancestry data; no replication.
- The output is a demonstration of method, not evidence about vitamin D and type 1 diabetes.

---

## Support

- **Issues:** open an issue in this repository's GitHub issue tracker
- **Contact:** see the Author section

---

## Changelog

### Version 1.0.0 (November 4, 2026)

- Initial release: instrument selection, IVW, weighted median, MR-Egger, heterogeneity, intercept, leave-one-out, single-variant analyses and plots

---

## Citations

- [**TwoSampleMR**](https://mrcieu.github.io/TwoSampleMR/): Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal inference across the human phenome. eLife 2018.
- [**GWAS Catalog GCST90000618**](https://www.ebi.ac.uk/gwas/studies/GCST90000618): serum 25-hydroxyvitamin D summary statistics, European ancestry, n = 496,946 (Revez JA et al., Nature Communications 2020).
- [**GWAS Catalog GCST90014023**](https://www.ebi.ac.uk/gwas/studies/GCST90014023): type 1 diabetes summary statistics (Chiou J et al., Nature 2021).
- [**GWAS Catalog**](https://www.ebi.ac.uk/gwas/): EMBL-EBI, source and harmonisation of both files.

---

## Author

**Tchicdje Kouojip Joram Smith (DeltaGa)**
Email: dev.github.tkjoramsmith@outlook.com
GitHub: [https://github.com/DeltaGa](https://github.com/DeltaGa)

---

© 2026 DeltaGa. All rights reserved.
