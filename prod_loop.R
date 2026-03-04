library(dplyr)
library(lubridate)
library(readr)

# ============================================================
# Batch runner: convert forecast CSVs -> dated outputs + latest summary
# ============================================================

# -----------------------------
# Config
# -----------------------------
forecast_dir <- "forecast_temps"
accu_prefix  <- "accuweather_dc_"         # expects accuweather_dc_YYYY-MM-DD.csv
accu_regex   <- "^accuweather_dc_\\d{4}-\\d{2}-\\d{2}\\.csv$"

# Historic inputs
hist_temp_file <- "historic_temps/dc_temps.csv"
bloom_file     <- "data/washingtondc.csv"

# Residuals from forward-chaining backtest
residuals_rds  <- "model_residuals_dc.rds"

# Output root
out_root <- "outputs"
runs_dir <- file.path(out_root, "runs")   # per-forecast-date outputs go here
summary_file <- file.path(out_root, "latest_summary.csv")
probs_file_latest <- file.path(out_root, "bloom_probabilities.csv")
temp_file_latest <- file.path(out_root, "temp_paths_sample50.csv")

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(runs_dir, showWarnings = FALSE, recursive = TRUE)

# Phenology params (your best combo)
params <- list(
  start_month = 2,
  start_day   = 1,
  Ts_c        = 20,
  Ea_kj       = 80
)

# Simulation knobs
N_paths <- 10000
phi     <- 0.85

# Mean blending (Accu -> climo)
tau_blend_days <- 30
w_min <- 0.10

# Uncertainty schedule (simple exponential ramp)
sigma0    <- 0.5        # °C near-term SD
tau_sigma <- 14         # days to approach long-range target

# Climatology window
clim_start_year <- 1990
clim_spread_scale <- 1.0

# Optional anomaly-based shrink of long-range sigma target
use_anom_shrink <- TRUE
min_scale <- 0.40
C_anom    <- 2.5

# Clamp simulated temps to climo q01/q99 (safety net)
do_clamp <- TRUE

# Simulation horizon (should go beyond plausible bloom date)
end_date <- as.Date("2026-05-31")

# Export a sample of temperature paths
n_paths_export <- 50

# Residual scaling (optional): set to 1.0 to keep as-is
res_scale <- 1.0
# If you want to keep DOY integer while scaling residuals, we round at the end:
# bloom_doy_adj <- round(bloom_doy_sim + res_scale * err)

# Deterministic seeds per run
seed_base <- 1234

# ============================================================
# Helpers
# ============================================================

extract_date_from_filename <- function(path) {
  fname <- basename(path)
  m <- regmatches(fname, regexpr("\\d{4}-\\d{2}-\\d{2}", fname))
  if (length(m) == 0 || is.na(m) || m == "") stop("Could not extract date from: ", fname)
  as.Date(m)
}

read_accu_file <- function(path) {
  x <- read_csv(path, show_col_types = FALSE)
  if (all(c("DATE", "TAVG") %in% names(x))) {
    x <- x %>% transmute(date = as.Date(DATE), tavg = as.numeric(TAVG))
  } else if (all(c("date", "tavg") %in% names(x))) {
    x <- x %>% transmute(date = as.Date(date), tavg = as.numeric(tavg))
  } else {
    stop("Accu file must contain either DATE/TAVG or date/tavg columns: ", path)
  }
  x %>% filter(!is.na(date), is.finite(tavg)) %>% arrange(date)
}

make_dts_fun <- function(Ts_c, Ea_kj) {
  R <- 8.314462618
  TsK <- Ts_c + 273.15
  Ea  <- Ea_kj * 1000
  function(tmean_c) {
    Tk <- tmean_c + 273.15
    exp((Ea / R) * (1 / TsK - 1 / Tk))
  }
}

predict_bloom_from_series <- function(dates, tmean_c, Fstar, params, dts_fun) {
  dfy <- tibble(DATE = dates, tmean = tmean_c) %>% arrange(DATE)
  yr <- year(dfy$DATE[1])
  start_date <- make_date(yr, params$start_month, params$start_day)
  dfy <- dfy %>% filter(DATE >= start_date)
  if (nrow(dfy) == 0) return(NA_real_)
  cum <- cumsum(dts_fun(dfy$tmean))
  idx <- which(cum >= Fstar)[1]
  if (is.na(idx)) return(NA_real_)
  yday(dfy$DATE[idx])
}

doy_to_date <- function(doy, year) {
  as.Date(floor(doy) - 1, origin = as.Date(sprintf("%d-01-01", year)))
}

sim_ar1 <- function(n, phi) {
  e <- numeric(n)
  innov_sd <- sqrt(1 - phi^2)
  for (i in 1:n) {
    eta <- rnorm(1, 0, innov_sd)
    e[i] <- if (i == 1) eta else phi * e[i - 1] + eta
  }
  e
}

# ============================================================
# 1) Load historic temps + build climatology once
# ============================================================

temps_raw <- read_csv(hist_temp_file, show_col_types = FALSE) %>%
  mutate(
    DATE = as.Date(DATE),
    year = year(DATE),
    TAVG = as.numeric(TAVG),
    TMAX = as.numeric(TMAX),
    TMIN = as.numeric(TMIN)
  )

temps_hist <- temps_raw %>%
  mutate(
    tmean = case_when(
      is.finite(TAVG) ~ TAVG,
      is.finite(TMAX) & is.finite(TMIN) ~ (TMAX + TMIN) / 2,
      TRUE ~ NA_real_
    ),
    doy = yday(DATE)
  ) %>%
  select(DATE, year, doy, tmean) %>%
  filter(!is.na(DATE), is.finite(tmean)) %>%
  arrange(year, DATE)

clim <- temps_hist %>%
  filter(year >= clim_start_year) %>%
  group_by(doy) %>%
  summarise(
    clim_mean = mean(tmean, na.rm = TRUE),
    clim_q10  = quantile(tmean, 0.10, na.rm = TRUE, type = 7),
    clim_q90  = quantile(tmean, 0.90, na.rm = TRUE, type = 7),
    clim_q01  = quantile(tmean, 0.01, na.rm = TRUE, type = 7),
    clim_q99  = quantile(tmean, 0.99, na.rm = TRUE, type = 7),
    .groups = "drop"
  ) %>%
  arrange(doy)

# ============================================================
# 2) Fit F* once from history (median cumDTS at observed bloom)
#    (If you prefer: replace with readRDS of pre-fit F*.)
# ============================================================

bloom <- read_csv(bloom_file, show_col_types = FALSE) %>%
  mutate(
    bloom_date = as.Date(bloom_date),
    year = as.integer(year),
    bloom_doy = as.integer(bloom_doy)
  ) %>%
  select(year, bloom_date, bloom_doy) %>%
  arrange(year)

years_common <- intersect(unique(temps_hist$year), unique(bloom$year))
temps_fit <- temps_hist %>% filter(year %in% years_common)
bloom_fit <- bloom %>% filter(year %in% years_common)

dts_fun <- make_dts_fun(params$Ts_c, params$Ea_kj)

cum_at_bloom <- temps_fit %>%
  group_by(year) %>%
  group_modify(~{
    yr <- .y$year
    dfy <- .x %>%
      arrange(DATE) %>%
      filter(DATE >= make_date(yr, params$start_month, params$start_day)) %>%
      mutate(cum = cumsum(dts_fun(tmean)))
    bd <- bloom_fit %>% filter(year == yr) %>% pull(bloom_date)
    val <- dfy %>% filter(DATE == bd) %>% pull(cum)
    tibble(cum_at_bloom = ifelse(length(val) == 0, NA_real_, val))
  }) %>%
  ungroup()

Fstar <- median(cum_at_bloom$cum_at_bloom, na.rm = TRUE)
cat("Using F* =", Fstar, "\n")

# ============================================================
# 3) Load residuals once
# ============================================================

residuals_df <- readRDS(residuals_rds)
residuals_vec <- residuals_df$residual
stopifnot(is.numeric(residuals_vec), length(residuals_vec) > 10)

# ============================================================
# 4) List forecast files & process missing outputs
# ============================================================

files <- list.files(forecast_dir, pattern = accu_regex, full.names = TRUE)
if (length(files) == 0) stop("No matching forecast files found in: ", forecast_dir)

# Sort by forecast date
forecast_dates_chr <- vapply(files, function(f) as.character(extract_date_from_filename(f)), character(1))
forecast_dates <- as.Date(forecast_dates_chr)

ord <- order(forecast_dates)
files <- files[ord]
forecast_dates <- forecast_dates[ord]

mt <- file.info(files)$mtime
latest_idx <- tail(order(forecast_dates, mt), 1)
latest_file <- files[latest_idx]
latest_date <- forecast_dates[latest_idx]

cat("Latest forecast detected:", basename(latest_file), "(", as.character(latest_date), ")\n")

summary_rows <- list()

for (i in seq_along(files)) {
  accu_path <- files[i]
  forecast_date <- forecast_dates[i]  # this will be "today" for this run
  
  run_dir <- file.path(runs_dir, format(forecast_date, "%Y-%m-%d"))
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_prob_csv  <- file.path(run_dir, "bloom_probabilities.csv")
  out_q_txt     <- file.path(run_dir, "bloom_quantiles.txt")
  out_paths_csv <- file.path(run_dir, "temp_paths_sample50.csv")
  
  outputs_exist <- file.exists(out_prob_csv) && 
    file.exists(out_q_txt) && 
    file.exists(out_paths_csv)
  
  is_latest <- identical(normalizePath(accu_path), normalizePath(latest_file))
  
  need_run <- is_latest || !outputs_exist
  
  cat("\n----", as.character(forecast_date), "----\n")
  cat("Input:", accu_path, "\n")
  cat("Output dir:", run_dir, "\n")
  cat("Need run?", need_run, "\n")
  
  if (!need_run) {
    # Load quantiles from existing txt for summary row (simple parse)
    # Expect lines: Q10: <doy>  (<date>)
    txt <- readLines(out_q_txt, warn = FALSE)
    get_q <- function(prefix) {
      line <- txt[grepl(paste0("^", prefix, ":"), txt)]
      if (length(line) == 0) return(NA_real_)
      as.numeric(sub(paste0("^", prefix, ":\\s*([0-9\\.]+).*"), "\\1", line[1]))
    }
    q10 <- get_q("Q10"); q50 <- get_q("Q50"); q90 <- get_q("Q90")
    
    # Mode from probabilities
    prob_tbl <- read_csv(out_prob_csv, show_col_types = FALSE)
    mode_date <- prob_tbl$DATE[which.max(prob_tbl$prob_percent)]
    mode_doy <- yday(as.Date(mode_date))
    
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      forecast_date = forecast_date,
      q10 = q10, q50 = q50, q90 = q90,
      mode = mode_doy
    )
    next
  }
  
  # Seed per run (reproducible)
  set.seed(seed_base + as.integer(forecast_date))
  
  # ---- Read Accu file and build calendar df ----
  accu <- read_accu_file(accu_path)
  
  dates_all <- seq(min(accu$date), end_date, by = "1 day")
  df <- tibble(DATE = dates_all, doy = yday(dates_all)) %>%
    left_join(clim, by = "doy") %>%
    left_join(accu %>% rename(DATE = date, accu_tavg = tavg), by = "DATE") %>%
    arrange(DATE)
  
  # Fill any missing climatology (DOY 366 edge cases)
  if (any(!is.finite(df$clim_mean))) {
    miss <- which(!is.finite(df$clim_mean))
    for (j in miss) {
      d <- df$doy[j]
      nearest <- clim$doy[which.min(abs(clim$doy - d))]
      df[j, c("clim_mean","clim_q10","clim_q90","clim_q01","clim_q99")] <- clim %>%
        filter(doy == nearest) %>%
        select(clim_mean, clim_q10, clim_q90, clim_q01, clim_q99)
    }
  }
  
  today <- forecast_date  # IMPORTANT: "today" comes from filename
  
  lead <- as.integer(df$DATE - today)
  k <- pmax(0, lead)
  has_accu <- is.finite(df$accu_tavg)
  
  # ---- Mean path mu ----
  w_decay <- exp(-k / tau_blend_days)
  w <- ifelse(df$DATE < today, 1,
              ifelse(has_accu, pmax(w_min, w_decay), 0))
  
  mu <- ifelse(
    df$DATE < today & has_accu,
    df$accu_tavg,
    ifelse(has_accu, w * df$accu_tavg + (1 - w) * df$clim_mean, df$clim_mean)
  )
  
  # ---- sigma_k (simple exponential ramp) ----
  z <- 1.2815515655446004
  sigma_clim <- (df$clim_q90 - df$clim_q10) / (2 * z)
  sigma_clim <- sigma_clim * clim_spread_scale
  
  if (use_anom_shrink) {
    mu_anom <- mu - df$clim_mean
    anom_scale <- pmax(min_scale, exp(-abs(mu_anom) / C_anom))
  } else {
    anom_scale <- rep(1.0, nrow(df))
  }
  
  sigma_target <- sigma_clim * anom_scale
  grow <- 1 - exp(-k / tau_sigma)
  
  sigma_k <- ifelse(df$DATE < today, 0,
                    sigma0 + (sigma_target - sigma0) * grow)
  
  sigma_k <- pmin(sigma_k, sigma_clim)
  
  # ---- Simulate temps + bloom ----
  temp_matrix <- matrix(NA_real_, nrow = nrow(df), ncol = N_paths)
  bloom_doy_sim <- numeric(N_paths)
  
  for (j in seq_len(N_paths)) {
    e <- sim_ar1(nrow(df), phi)
    t_path <- mu + sigma_k * e
    t_path[df$DATE < today] <- mu[df$DATE < today]
    
    if (do_clamp) {
      t_path <- pmin(t_path, df$clim_q99)
      t_path <- pmax(t_path, df$clim_q01)
    }
    
    temp_matrix[, j] <- t_path
    bloom_doy_sim[j] <- predict_bloom_from_series(df$DATE, t_path, Fstar, params, dts_fun)
  }
  
  ok <- is.finite(bloom_doy_sim)
  bloom_doy_sim <- bloom_doy_sim[ok]
  temp_matrix <- temp_matrix[, ok, drop = FALSE]
  
  # ---- Add phenology model error (bootstrap residuals) ----
  err <- sample(residuals_vec, size = length(bloom_doy_sim), replace = TRUE)
  bloom_doy_adj <- round(bloom_doy_sim + res_scale * err)  # integer day output
  
  year0 <- year(min(df$DATE))
  bloom_date_adj <- doy_to_date(bloom_doy_adj, year0)
  
  # ---- Output 1: probabilities CSV ----
  prob_tbl <- tibble(DATE = bloom_date_adj) %>%
    count(DATE, name = "n") %>%
    mutate(prob_percent = 100 * n / sum(n)) %>%
    arrange(DATE) %>%
    select(DATE, prob_percent)
  
  prob_tbl <- prob_tbl %>%
    mutate(DATE = format(as.Date(DATE), "%Y-%m-%d"))
  
  write_csv(prob_tbl, out_prob_csv)
  write_csv(prob_tbl, probs_file_latest)
  
  # ---- Output 2: quantiles TXT ----
  qs <- quantile(bloom_doy_adj, probs = c(0.10, 0.50, 0.90), na.rm = TRUE, type = 7)
  q_dates <- doy_to_date(as.numeric(qs), year0)
  
  mode_date <- prob_tbl$DATE[which.max(prob_tbl$prob_percent)]
  mode_doy <- yday(as.Date(mode_date))
  
  q_lines <- c(
    sprintf("Run timestamp: %s", as.character(Sys.time())),
    sprintf("Forecast date (today): %s", as.character(today)),
    sprintf("Forecast year: %d", year0),
    sprintf("N_paths used: %d", length(bloom_doy_adj)),
    sprintf("Params: start=%02d-%02d, Ts_c=%.1f, Ea_kj=%.1f", params$start_month, params$start_day, params$Ts_c, params$Ea_kj),
    sprintf("F*: %.6f", Fstar),
    sprintf("Residuals: file=%s | res_scale=%.2f", residuals_rds, res_scale),
    "",
    "Bloom DOY quantiles (after residual bootstrap):",
    sprintf("Q10: %.2f  (%s)", qs[[1]], as.character(q_dates[[1]])),
    sprintf("Q50: %.2f  (%s)", qs[[2]], as.character(q_dates[[2]])),
    sprintf("Q90: %.2f  (%s)", qs[[3]], as.character(q_dates[[3]])),
    sprintf("MODE: %d  (%s)  [%.2f%%]",
            mode_doy, as.character(mode_date), max(prob_tbl$prob_percent))
  )
  writeLines(q_lines, out_q_txt)
  
  # ---- Output 3: sample 50 temperature paths ----
  set.seed(42 + as.integer(today))
  cols <- sample(seq_len(ncol(temp_matrix)), size = min(n_paths_export, ncol(temp_matrix)), replace = FALSE)
  
  paths_sample <- tibble(
    DATE = rep(df$DATE, times = length(cols)),
    path_id = rep(seq_along(cols), each = nrow(df)),
    tmean_c = as.vector(temp_matrix[, cols])
  ) %>%
    mutate(observed = DATE < today)
  
  paths_sample <- paths_sample %>%
    mutate(DATE = format(as.Date(DATE), "%Y-%m-%d"))
  
  write_csv(paths_sample, out_paths_csv)
  write_csv(paths_sample, temp_file_latest)
  
  cat("Wrote outputs:\n",
      " - ", out_prob_csv, "\n",
      " - ", out_q_txt, "\n",
      " - ", out_paths_csv, "\n", sep = "")
  
  summary_rows[[length(summary_rows) + 1]] <- tibble(
    forecast_date = forecast_date,
    q10 = as.numeric(qs[[1]]),
    q50 = as.numeric(qs[[2]]),
    q90 = as.numeric(qs[[3]]),
    mode = as.numeric(mode_doy)
  )
}

# ============================================================
# 5) Write ONE summary file: forecast_date,q10,q50,q90,mode
# ============================================================

summary_df <- bind_rows(summary_rows) %>%
  arrange(forecast_date)

summary_df <- summary_df %>%
  mutate(forecast_date = format(as.Date(forecast_date), "%Y-%m-%d"))

write_csv(summary_df, summary_file)
cat("\nWrote summary:", summary_file, "\n")
cat("Rows:", nrow(summary_df), "| from", as.character(min(summary_df$forecast_date)),
    "to", as.character(max(summary_df$forecast_date)), "\n")