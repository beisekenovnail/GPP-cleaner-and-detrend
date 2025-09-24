##############################################
# GPP cleaner and detrend (row = pixel, cols = dates)
# 
# Last update: 24 Sep 2025 
# Author: Naıl Beisekenov
# 
# Input (WIDE; columns may be in ANY order):
#   pxID, x, y, [optional meta… can be anywhere],
#   2000-03-05, 2000-03-21, 2000-04-06, ...
#
# What this script does:
#   • Reads only the header first; header tokens are treated as PLAIN STRINGS
#     and NEVER parsed as dates (prevents charToDate errors).
#   • Detects date columns anywhere in the header via strict YYYY-MM-DD regex.
#   • Streams the CSV by row chunks (memory-safe for multi-GB files).
#   • Normalizes “no-data” sentinels to NA (and optional zero→NA).
#   • Drops any row (pixel) that has NA in ANY date column after normalization.
#   • Computes three row-wise detrended outputs across the time axis (columns):
#       1) STL remainder  (trend+seasonality removed; robust fallback)
#       2) Climatology    (value − mean_by_MMDD across timeline)
#       3) SSA remainder  (trend via SSA; robust fallback if SSA unavailable)
#   • Appends to four CSVs incrementally (never holding full results in RAM).
#
# Outputs (append mode):
#   1) gpp_noNA.csv
#   2) gpp_detrended_STL.csv
#   3) gpp_detrended_climatology.csv
#   4) gpp_detrended_SSA.csv
##############################################

## 0) Packages (incl. optional Rssa) ------------------------------------------
need <- c("data.table","stringr","zoo")
inst <- setdiff(need, rownames(installed.packages()))
if (length(inst)) install.packages(inst, repos = "https://cran.rstudio.com/")
invisible(lapply(need, require, character.only = TRUE))

## 1) Paths & configuration ---------------------------------------------------
# >>> EDIT THIS ONE LINE <<<
csv_path <- "/Users/beisekenovnail/Documents/2025 PHD/Dr. Feagin/tippy/R_codes/interesting_pixels/TW_GPP_2025_all_pixels.csv"

out_dir  <- file.path(dirname(csv_path), "6claire_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Streaming & compute
CHUNK_ROWS <- 10000L                                  # rows per chunk
N_THREADS  <- max(1L, parallel::detectCores() - 1L)   # fread threads
FREQ_GUESS <- 23L                                     # ~16-day cadence (~23 per year) for STL
SSA_MAX_N  <- Inf                                     # cap SSA length if needed
SSA_L_FRAC <- 1/3                                     # SSA window length L ≈ n/3

# Missing-data normalization
TREAT_ZERO_AS_MISSING <- TRUE          # zeros → NA?
ZERO_TOL              <- 1e-12         # near-zero tolerance
MISSING_SENTINELS     <- c(-9999, -3000, -999.0)  # domain “no-data” codes

# Meta toggle — drop these 3 from ALL outputs (no matter where they are)
REMOVE_META_FROM_OUTPUTS <- TRUE
META_CANDIDATES <- c("TW_area_per_PX(m^2)", "Woody_percent", "Herb_percent")

# SSA toggle (+ auto-install if missing)
DO_SSA <- TRUE

## 2) Header-only read & column classification --------------------------------
hdr <- names(data.table::fread(csv_path, nrows = 0, showProgress = FALSE))
if (length(hdr) < 4L) stop("Header must have ≥4 columns (pxID, x, y, and date columns).")

# Strict date header detector (regex only; NEVER parse)
is_date_like <- function(s) grepl("^\\d{4}-\\d{2}-\\d{2}$", s)

# Detect date columns ANYWHERE in the header
date_cols <- hdr[is_date_like(hdr)]
if (length(date_cols) < 2L) {
  stop("No date columns detected (need at least 2 headers like YYYY-MM-DD).")
}

# Non-date columns are ID/META (order preserved as in the source file)
non_date_cols <- setdiff(hdr, date_cols)

# We require base ID columns somewhere among non-dates (front, middle, or end)
base_id_cols <- c("pxID","x","y")
if (!all(base_id_cols %in% non_date_cols)) {
  stop("Missing required ID columns pxID, x, y (they can be anywhere in the header).")
}

# Optional meta columns (may be anywhere)
meta_cols <- intersect(META_CANDIDATES, non_date_cols)

# All ID/META columns in the order they appear in the source (kept in memory)
idmeta_cols_all <- non_date_cols

# Columns written to outputs (drop meta if requested)
output_id_cols <- if (isTRUE(REMOVE_META_FROM_OUTPUTS)) {
  # Keep only base IDs in outputs
  base_id_cols
} else {
  idmeta_cols_all
}

# Final output header = chosen ID/META + all date columns (dates follow IDs)
output_headers <- c(output_id_cols, date_cols)

## 2b) Optional: auto-install Rssa if SSA is requested ------------------------
HAS_RSSA <- requireNamespace("Rssa", quietly = TRUE)
if (isTRUE(DO_SSA) && !HAS_RSSA) {
  message("Rssa not found. Attempting to install Rssa ...")
  try({
    install.packages("Rssa", repos = "https://cran.rstudio.com/")
  }, silent = TRUE)
  HAS_RSSA <- requireNamespace("Rssa", quietly = TRUE)
  if (!HAS_RSSA) {
    message("WARNING: Rssa could not be installed automatically. ",
            "SSA output will use a robust fallback (no failure).")
  } else {
    message("Rssa installed successfully.")
  }
}

## 3) Detrending helpers (row-wise across time axis) --------------------------

# STL remainder with robust fallback
stl_anomaly_row <- function(y, freq_guess = FREQ_GUESS) {
  y <- as.numeric(y)
  tsx <- tryCatch(stats::ts(y, frequency = freq_guess),
                  error = function(e) stats::ts(y, frequency = 1))
  fit <- tryCatch(stats::stl(tsx,
                             s.window = ifelse(freq_guess >= 12, "periodic", "periodic"),
                             robust   = TRUE),
                  error = function(e) NULL)
  if (is.null(fit)) {
    n <- length(y)
    k <- max(7, 2*floor(max(7, round(n*0.10))/2)+1)    # odd window, ~10% of series
    trend <- tryCatch(zoo::rollmedian(y, k = k, fill = "extend", align = "center"),
                      error = function(e) rep(stats::median(y, na.rm = TRUE), n))
    return(y - trend)
  }
  as.numeric(fit$time.series[, "remainder"])
}

# Climatology anomaly: value − mean_by_MMDD across the timeline
climatology_anomaly_row <- function(y, date_keys) {
  y <- as.numeric(y)
  if (length(y) != length(date_keys)) {
    stop(sprintf("Length mismatch: y=%d, keys=%d", length(y), length(date_keys)))
  }
  mmdd <- stringr::str_sub(date_keys, 6, 10)          # "MM-DD"
  climo_by_mmdd <- tapply(y, mmdd, mean, na.rm = TRUE)
  climo_vec     <- unname(climo_by_mmdd[mmdd])        # back to original order
  y - climo_vec
}

# SSA remainder: trend via SSA; robust fallback if unavailable
ssa_anomaly_row <- function(y) {
  y <- as.numeric(y); n <- length(y)
  # If SSA not available/short/degenerate → robust fallback (never fail)
  if (!HAS_RSSA || !isTRUE(DO_SSA) || n < 20 || !any(is.finite(y))) {
    k <- max(7, 2*floor(max(7, round(n*0.10))/2)+1)
    trend <- tryCatch(zoo::rollmedian(y, k = k, fill = "extend", align = "center"),
                      error = function(e) rep(stats::median(y, na.rm = TRUE), n))
    return(y - trend)
  }
  if (n > SSA_MAX_N) {
    y <- tail(y, SSA_MAX_N); n <- length(y)
  }
  L <- max(10L, min(n - 1L, floor(n * SSA_L_FRAC)))   # window length ~ n/3
  ssa_fit <- tryCatch(Rssa::ssa(y, L = L), error = function(e) NULL)
  if (is.null(ssa_fit)) {
    k <- max(7, 2*floor(max(7, round(n*0.10))/2)+1)
    trend <- tryCatch(zoo::rollmedian(y, k = k, fill = "extend", align = "center"),
                      error = function(e) rep(stats::median(y, na.rm = TRUE), n))
    return(y - trend)
  }
  comps_max <- min(6L, length(ssa_fit$values))
  trend <- tryCatch({
    rec <- Rssa::reconstruct(ssa_fit, groups = list(trend = 1L:min(2L, comps_max)))
    as.numeric(rec$trend)
  }, error = function(e) rep(NA_real_, n))
  if (!any(is.finite(trend))) {
    k <- max(7, 2*floor(max(7, round(n*0.10))/2)+1)
    trend <- tryCatch(zoo::rollmedian(y, k = k, fill = "extend", align = "center"),
                      error = function(e) rep(stats::median(y, na.rm = TRUE), n))
  }
  y - trend
}

## 4) Create output files (headers only, 0 rows) ------------------------------
out_noNA <- file.path(out_dir, "gpp_noNA.csv")
out_stl  <- file.path(out_dir, "gpp_detrended_STL.csv")
out_clm  <- file.path(out_dir, "gpp_detrended_climatology.csv")
out_ssa  <- file.path(out_dir, "gpp_detrended_SSA.csv")

mk_empty_dt <- function(headers, id_cols_force_char) {
  dt <- data.table::as.data.table(setNames(rep(list(numeric(0)), length(headers)), headers))
  common <- intersect(id_cols_force_char, names(dt))
  if (length(common)) dt[, (common) := lapply(.SD, as.character), .SDcols = common]
  dt
}

data.table::fwrite(mk_empty_dt(output_headers, output_id_cols), out_noNA, append = FALSE)
data.table::fwrite(mk_empty_dt(output_headers, output_id_cols), out_stl,  append = FALSE)
data.table::fwrite(mk_empty_dt(output_headers, output_id_cols), out_clm,  append = FALSE)
data.table::fwrite(mk_empty_dt(output_headers, output_id_cols), out_ssa,  append = FALSE)

if (isTRUE(DO_SSA) && !HAS_RSSA) {
  message("NOTE: SSA fallbacks will be used (Rssa not available). Install 'Rssa' for true SSA.")
}

## 5) Streaming loop: read → normalize → filter → detrend → append ------------
rows_read <- 0L
chunk_idx <- 0L
date_keys <- date_cols

pixels_kept_total    <- 0L
pixels_dropped_total <- 0L

repeat {
  # fread in chunks; guard EOF “skip=… only has … lines” with tryCatch
  DT <- tryCatch(
    data.table::fread(
      csv_path,
      skip         = 1L + rows_read,
      nrows        = CHUNK_ROWS,
      header       = FALSE,
      col.names    = hdr,                     # use original header for column names
      showProgress = (chunk_idx %% 10L == 0L),
      na.strings   = c("NA","NaN","","null","NULL"),
      nThread      = N_THREADS
    ),
    error = function(e) data.table::data.table()
  )
  if (nrow(DT) == 0L) break  # EOF
  
  chunk_idx <- chunk_idx + 1L
  
  # Ensure all columns exist and keep the canonical order:
  # [all ID/META in source order] + [all DATE columns]
  # (This preserves any non-date metas even if they were at the end originally.)
  data.table::setcolorder(DT, c(idmeta_cols_all, date_cols))
  
  # Make date columns numeric
  suppressWarnings({
    for (nm in date_cols) if (!is.numeric(DT[[nm]])) DT[[nm]] <- as.numeric(DT[[nm]])
  })
  
  # Normalize sentinels and zeros to NA (DATE columns only)
  if (length(MISSING_SENTINELS)) {
    for (nm in date_cols) {
      v <- DT[[nm]]
      bad <- is.finite(v) & (v %in% MISSING_SENTINELS)
      if (any(bad)) v[bad] <- NA_real_
      DT[[nm]] <- v
    }
  }
  if (isTRUE(TREAT_ZERO_AS_MISSING)) {
    for (nm in date_cols) {
      v <- DT[[nm]]
      bad0 <- is.finite(v) & (abs(v) <= ZERO_TOL)
      if (any(bad0)) v[bad0] <- NA_real_
      DT[[nm]] <- v
    }
  }
  
  # Keep only rows with NO NA across all date columns
  ok <- DT[, stats::complete.cases(.SD), .SDcols = date_cols]
  kept    <- sum(ok); dropped <- nrow(DT) - kept
  pixels_kept_total    <- pixels_kept_total    + kept
  pixels_dropped_total <- pixels_dropped_total + dropped
  if (!any(ok)) { rows_read <- rows_read + nrow(DT); next }
  
  DT_ok <- DT[ok]
  
  # --- WRITE 1) Original (no-NA) rows — with ID columns filtered per config
  # Select exactly the columns we want to appear in outputs: output_id_cols + date_cols
  data.table::fwrite(DT_ok[, c(output_id_cols, date_cols), with = FALSE], out_noNA, append = TRUE)
  
  # Build numeric matrix for detrending (rows = pixels, cols = dates)
  M <- as.matrix(DT_ok[, ..date_cols]); storage.mode(M) <- "double"
  
  # Row-wise detrends
  stl_mat <- t(apply(M, 1L, function(y) stl_anomaly_row(y)))
  clm_mat <- t(apply(M, 1L, function(y) climatology_anomaly_row(y, date_keys)))
  ssa_mat <- if (isTRUE(DO_SSA)) t(apply(M, 1L, function(y) ssa_anomaly_row(y))) else NULL
  
  # Bind back only the ID columns configured for outputs
  stl_out <- data.table::data.table(DT_ok[, ..output_id_cols], data.table::as.data.table(stl_mat))
  data.table::setnames(stl_out, old = names(stl_out)[-(seq_along(output_id_cols))], new = date_cols)
  
  clm_out <- data.table::data.table(DT_ok[, ..output_id_cols], data.table::as.data.table(clm_mat))
  data.table::setnames(clm_out, old = names(clm_out)[-(seq_along(output_id_cols))], new = date_cols)
  
  data.table::fwrite(stl_out, out_stl, append = TRUE)
  data.table::fwrite(clm_out, out_clm, append = TRUE)
  
  if (!is.null(ssa_mat)) {
    ssa_out <- data.table::data.table(DT_ok[, ..output_id_cols], data.table::as.data.table(ssa_mat))
    data.table::setnames(ssa_out, old = names(ssa_out)[-(seq_along(output_id_cols))], new = date_cols)
    data.table::fwrite(ssa_out, out_ssa, append = TRUE)
  }
  
  rows_read <- rows_read + nrow(DT)
  if (chunk_idx %% 10L == 0L) {
    message(sprintf("Processed ~%d rows (chunk %d) — kept=%d, dropped=%d",
                    rows_read, chunk_idx, kept, dropped))
  }
}

## 6) Final messages ----------------------------------------------------------
message("Done.")
message(" - Saved: ", out_noNA)
message(" - Saved: ", out_stl)
message(" - Saved: ", out_clm)
message(" - Saved: ", out_ssa)
message(sprintf("Pixels kept (no NA after normalization): %d", pixels_kept_total))
message(sprintf("Pixels dropped due to missing/sentinel/zero: %d", pixels_dropped_total))
if (isTRUE(DO_SSA) && !HAS_RSSA) {
  message("SSA ran with robust fallbacks (Rssa not available). Install 'Rssa' for true SSA.")
}
