# GPP-cleaner-and-detrend

This repo contains an R workflow that cleans a very wide GPP CSV (rows = pixels, columns = dates) and produces three detrended variants plus a clean “no-NA” baseline—without blowing up memory on multi-GB files.
It is built for files shaped like:
pxID, x, y, [optional meta...], 2000-03-05, 2000-03-21, 2000-04-06, ...
Key ideas:
•	Streaming in chunks (no massive long-table in RAM)
•	Strict missing-data policy (normalize sentinels; optional zero→NA)
•	Row-wise detrend across time: STL remainder, SSA remainder, Climatology anomaly
•	Append-mode outputs so you can process 2–10 GB CSVs safely
•	Optional removal of three meta columns from outputs:
TW_area_per_PX(m2), woody_percent, herb_percent
Outputs
All outputs preserve the original date grid (same set/order of date columns), but only include rows that pass the “complete series” filter.
1.	gpp_noNA.csv — raw values for rows with no missing after normalization
2.	gpp_detrended_STL.csv — STL remainder: y − (trend + seasonal)
3.	gpp_detrended_climatology.csv — value − mean_by_MMDD (seasonality removed, slow drift preserved)
4.	gpp_detrended_SSA.csv — SSA remainder (data-driven smooth trend removed). If Rssa isn’t installed and DO_SSA=TRUE, the script installs it automatically.
What “missing after normalization” means
In date columns we first turn sentinel codes (e.g., -9999) into NA. Optionally, we also treat zeros (within a tolerance) as NA if zeros really mean “no data” for your product.
Any row with any NA after this step is dropped from all outputs.
Quick start
1.	Clone this repo and open R/RStudio.
2.	Open gpp_clean_detrend_streaming.R.
3.	Set the path to your wide CSV:
4.	csv_path <- "/path/to/TW_GPP_2025_all_pixels.csv"
5.	(Optional) Adjust configuration at the top (see “Configuration” below).
6.	Source the script:
7.	source("gpp_clean_detrend_streaming.R")
The script creates an output folder next to your CSV (e.g., 3claire_outputs/) and writes the four CSVs incrementally.
Requirements
•	R ≥ 4.1 recommended
•	CRAN packages: data.table, stringr, zoo
(the script installs anything missing)
•	Optional for SSA: Rssa
If DO_SSA=TRUE and Rssa is missing, the script will install it:
•	if (!requireNamespace("Rssa", quietly = TRUE)) install.packages("Rssa")
Tip: For reproducible environments consider renv::init() and snapshotting the library.
Configuration (edit at the top of the script)
# Input / output
csv_path <- "/path/to/TW_GPP_2025_all_pixels.csv"
out_dir  <- file.path(dirname(csv_path), "3claire_outputs")

# Streaming
CHUNK_ROWS <- 10000L                             # rows per fread() chunk (10k–25k typical)
N_THREADS  <- max(1L, parallel::detectCores() - 1L)

# STL settings
FREQ_GUESS <- 23L                                # ~16-day cadence (≈23 obs/year)

# SSA settings
DO_SSA     <- TRUE                               # compute SSA remainder
SSA_MAX_N  <- Inf                                # cap series length (optional)
SSA_L_FRAC <- 1/3                                # SSA window L ≈ n/3

# Missing-data normalization
TREAT_ZERO_AS_MISSING <- TRUE                    # treat ~0.0 as NA if zeros mean “no data”
ZERO_TOL              <- 1e-12                   # tolerance for numeric zero
MISSING_SENTINELS     <- c(-9999, -3000, -999.0) # domain-specific fill values → NA

# Remove these meta columns from outputs
REMOVE_META_FROM_OUTPUTS <- TRUE
META_CANDIDATES <- c("TW_area_per_PX(m2)", "woody_percent", "herb_percent")
Examples:
•	Keep zeros as valid data:
•	TREAT_ZERO_AS_MISSING <- FALSE
•	Keep meta columns in outputs:
•	REMOVE_META_FROM_OUTPUTS <- FALSE
•	Skip SSA:
•	DO_SSA <- FALSE
•	Limit SSA to the most recent 800 columns:
•	SSA_MAX_N <- 800

How detrending works:
•	STL remainder: we convert each row to a ts with frequency = FREQ_GUESS and compute trend + seasonal via stats::stl. The output is the remainder (y − (trend + seasonal)). If STL can’t fit, we fall back to a robust rolling-median trend.
•	Climatology anomaly: we derive MM-DD from each date column name and compute the mean for each calendar day within the row (across years). Each value is y − mean[MM-DD]. This removes seasonality but preserves slow multi-year drift.
•	SSA remainder: we fit SSA with L ≈ n/3, reconstruct a low-rank “trend” (first components), and output y − trend. If SSA is disabled or fails, we reuse the robust rolling-median fallback.
All three methods run row-wise across the time axis of each pixel and write results in append mode per chunk.
Input assumptions
•	The first three columns are pxID, x, y.
•	All remaining columns are dates formatted as YYYY-MM-DD in the header.
•	Date headers are treated as plain strings (not parsed)—safer and much faster.

Directory structure (after a run)
your_data/
  TW_GPP_2025_all_pixels.csv
  3claire_outputs/
    gpp_noNA.csv
    gpp_detrended_STL.csv
    gpp_detrended_climatology.csv
    gpp_detrended_SSA.csv
    
Troubleshooting
•	No rows in outputs: Your “complete series” filter is strict. Consider setting TREAT_ZERO_AS_MISSING <- FALSE or removing sentinel values from MISSING_SENTINELS.
•	Too many rows dropped: Check whether zeros truly are “no data” for your source. If not, disable the zero-as-missing rule.
•	SSA slow: Set SSA_MAX_N to cap per-row length, or DO_SSA <- FALSE.


