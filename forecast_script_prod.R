library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)

# ============================================================
# PRODUCTION RUNNER (DC) — AccuWeather + Climatology + DTS + Model Error
#
# Outputs:
# 1) CSV  : outputs/bloom_probabilities.csv      (DATE, prob_percent)
# 2) TXT  : outputs/bloom_quantiles.txt          (q10/q50/q90 DOY+DATE + metadata)
# 3) CSV  : outputs/temp_paths_sample50.csv      (DATE, path_id, tmean_c)
# 4) PNG  : outputs/bloom_histogram.png          (histogram of bloom DOY after model error)
# 5) PNG  : outputs/temp_uncertainty_ribbon.png  (simulated daily q10–q90 band + median + climo + Accu)
# 6) PNG  : outputs/temp_width_over_time.png     (width (q90-q10) sim vs climo)
# ============================================================

# -----------------------------
# 0) Files + knobs
# -----------------------------

hist_temp_file <- "historic_temps/dc_temps.csv"
bloom_file     <- "data/washingtondc.csv"
accu_file      <- "forecast_temps/accuweather_dc_2026-02-26.csv"

# Residuals from your forward-chaining backtest
residuals_rds  <- "model_residuals_dc.rds"

# Output directory + files
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_prob_csv   <- file.path(out_dir, "bloom_probabilities.csv")
out_q_txt      <- file.path(out_dir, "bloom_quantiles.txt")
out_paths_csv  <- file.path(out_dir, "temp_paths_sample50.csv")

out_hist_png   <- file.path(out_dir, "bloom_histogram.png")
out_temp_png   <- file.path(out_dir, "temp_uncertainty_ribbon.png")
out_width_png  <- file.path(out_dir, "temp_width_over_time.png")

# Phenology parameters (use your best combo)
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
sigma0    <- 0.5        # °C (near-term SD)
tau_sigma <- 14         # days (how fast uncertainty approaches long-range target)

# Long-range sigma target derived from climo q10/q90, optionally scaled
clim_spread_scale <- 1.0

# Optional anomaly-based shrink of long-range sigma target
use_anom_shrink <- TRUE
min_scale <- 0.40
C_anom    <- 2.5

# Clamp simulated daily temps to climo q01/q99 (safety net)
do_clamp <- TRUE

# Forecast simulation end date
end_date <- as.Date("2026-05-31")

# Sample paths to export
n_paths_export <- 50

# Plot window (x-axis) for temperature diagnostics
plot_start <- as.Date("2026-02-01")
plot_end   <- as.Date("2026-04-30")

set.seed(1)

# -----------------------------
# 1) Load historical temps + climatology (1990+)
# -----------------------------

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
  filter(year >= 1990) %>%
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

# -----------------------------
# 2) Load bloom + fit F* (median cumDTS at bloom) from history
#    (Later: replace with readRDS of saved F* + params.)
# -----------------------------

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

make_dts_fun <- function(Ts_c, Ea_kj) {
  R <- 8.314462618
  TsK <- Ts_c + 273.15
  Ea  <- Ea_kj * 1000
  function(tmean_c) {
    Tk <- tmean_c + 273.15
    exp((Ea / R) * (1 / TsK - 1 / Tk))
  }
}
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
cat("F*:", Fstar, "\n")

predict_bloom_from_series <- function(dates, tmean_c, Fstar, params, dts_fun) {
  dfy <- tibble(DATE = dates, tmean = tmean_c) %>% arrange(DATE)
  yr <- year(dfy$DATE[1])
  start_date <- make_date(yr, params$start_month, params$start_day)
  dfy <- dfy %>% filter(DATE >= start_date)
  if (nrow(dfy) == 0) return(NA_integer_)
  cum <- cumsum(dts_fun(dfy$tmean))
  idx <- which(cum >= Fstar)[1]
  if (is.na(idx)) return(NA_integer_)
  yday(dfy$DATE[idx])
}

# -----------------------------
# 3) Load AccuWeather CSV
# -----------------------------

accu_raw <- read_csv(accu_file, show_col_types = FALSE)

if (all(c("DATE", "TAVG") %in% names(accu_raw))) {
  accu <- accu_raw %>% transmute(date = as.Date(DATE), tavg = as.numeric(TAVG))
} else if (all(c("date", "tavg") %in% names(accu_raw))) {
  accu <- accu_raw %>% transmute(date = as.Date(date), tavg = as.numeric(tavg))
} else {
  stop("Accu file must contain either DATE/TAVG or date/tavg columns.")
}

accu <- accu %>% filter(!is.na(date), is.finite(tavg)) %>% arrange(date)
cat("Accu rows:", nrow(accu), "|", as.character(min(accu$date)), "to", as.character(max(accu$date)), "\n")

# -----------------------------
# 4) Build calendar, join climo + Accu
# -----------------------------

fname <- basename(accu_file)

# extract first YYYY-MM-DD pattern
date_str <- sub(".*?(\\d{4}-\\d{2}-\\d{2}).*", "\\1", fname)

if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", date_str)) {
  stop("Could not extract date from filename: ", fname)
}

today <- as.Date(date_str)

cat("Using forecast date derived from filename:", as.character(today), "\n")

dates_all <- seq(min(accu$date), end_date, by = "1 day")

df <- tibble(DATE = dates_all, doy = yday(dates_all)) %>%
  left_join(clim, by = "doy") %>%
  left_join(accu %>% rename(DATE = date, accu_tavg = tavg), by = "DATE") %>%
  arrange(DATE)

# fill missing climatology (DOY edge cases)
if (any(!is.finite(df$clim_mean))) {
  miss <- which(!is.finite(df$clim_mean))
  for (i in miss) {
    d <- df$doy[i]
    nearest <- clim$doy[which.min(abs(clim$doy - d))]
    df[i, c("clim_mean","clim_q10","clim_q90","clim_q01","clim_q99")] <- clim %>%
      filter(doy == nearest) %>%
      select(clim_mean, clim_q10, clim_q90, clim_q01, clim_q99)
  }
}

# -----------------------------
# 5) Mean path mu (Accu blended to climo, with minimum weight)
# -----------------------------

lead <- as.integer(df$DATE - today)
k <- pmax(0, lead)
has_accu <- is.finite(df$accu_tavg)

w_decay <- exp(-k / tau_blend_days)
w <- ifelse(df$DATE < today, 1,
            ifelse(has_accu, pmax(w_min, w_decay), 0))

mu <- ifelse(
  df$DATE < today & has_accu,
  df$accu_tavg,  # treat past as observed
  ifelse(has_accu, w * df$accu_tavg + (1 - w) * df$clim_mean, df$clim_mean)
)

cat("Mean Accu weight (future, where exists):",
    mean(w[df$DATE >= today & has_accu], na.rm = TRUE), "\n")

# -----------------------------
# 6) Uncertainty sigma_k (simple exponential ramp)
# -----------------------------

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

i0 <- which(df$DATE >= today)[1]
cat("sigma_k day0:", sigma_k[i0], " (target sigma0=", sigma0, ")\n")

# -----------------------------
# 7) Load model residuals for phenology error
# -----------------------------

residuals_df <- readRDS(residuals_rds)
residuals_vec <- residuals_df$residual
stopifnot(is.numeric(residuals_vec), length(residuals_vec) > 10)

# -----------------------------
# 8) Simulate temperature paths + bloom DOYs
# -----------------------------

sim_ar1 <- function(n, phi) {
  e <- numeric(n)
  innov_sd <- sqrt(1 - phi^2)
  for (i in 1:n) {
    eta <- rnorm(1, 0, innov_sd)
    e[i] <- if (i == 1) eta else phi * e[i - 1] + eta
  }
  e
}

temp_matrix <- matrix(NA_real_, nrow = nrow(df), ncol = N_paths)
bloom_doy_sim <- numeric(N_paths)

for (j in seq_len(N_paths)) {
  e <- sim_ar1(nrow(df), phi)
  t_path <- mu + sigma_k * e
  
  # past locked to observed mean
  t_path[df$DATE < today] <- mu[df$DATE < today]
  
  # optional clamp to avoid wild tails
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

# -----------------------------
# 9) Add phenology model error (bootstrap residuals)
# -----------------------------

err <- sample(residuals_vec, size = length(bloom_doy_sim), replace = TRUE)
bloom_doy_adj <- round(bloom_doy_sim + err*0.8)

year0 <- year(min(df$DATE))
doy_to_date <- function(doy, year) {
  as.Date(doy - 1, origin = as.Date(sprintf("%d-01-01", year)))
}
bloom_date_adj <- doy_to_date(bloom_doy_adj, year0)

# -----------------------------
# 10) Output 1: DATE + probability of peak bloom (%)
# -----------------------------

prob_tbl <- tibble(DATE = bloom_date_adj) %>%
  count(DATE, name = "n") %>%
  mutate(prob_percent = 100 * n / sum(n)) %>%
  arrange(DATE) %>%
  select(DATE, prob_percent)

write_csv(prob_tbl, out_prob_csv)
cat("Wrote:", out_prob_csv, "\n")

# -----------------------------
# 11) Output 2: q10 / q50 / q90 of predicted peak bloom
# -----------------------------

qs <- quantile(bloom_doy_adj, probs = c(0.10, 0.50, 0.90), na.rm = TRUE, type = 7)
q_dates <- doy_to_date(as.numeric(qs), year0)

q_lines <- c(
  sprintf("Run date: %s", as.character(Sys.time())),
  sprintf("Forecast year: %d", year0),
  sprintf("N_paths used: %d", length(bloom_doy_adj)),
  sprintf("Phenology params: start=%02d-%02d, Ts_c=%.1f, Ea_kj=%.1f", params$start_month, params$start_day, params$Ts_c, params$Ea_kj),
  sprintf("Forcing threshold F*: %.6f", Fstar),
  "",
  "Bloom DOY quantiles (after model-error bootstrap):",
  sprintf("Q10: %.2f  (%s)", qs[[1]], as.character(q_dates[[1]])),
  sprintf("Q50: %.2f  (%s)", qs[[2]], as.character(q_dates[[2]])),
  sprintf("Q90: %.2f  (%s)", qs[[3]], as.character(q_dates[[3]]))
)

writeLines(q_lines, out_q_txt)
cat("Wrote:", out_q_txt, "\n")

# -----------------------------
# 12) Output 3: subsample 50 temperature paths to CSV
# -----------------------------

set.seed(42)
cols <- sample(seq_len(ncol(temp_matrix)), size = min(n_paths_export, ncol(temp_matrix)), replace = FALSE)

paths_sample <- tibble(
  DATE = rep(df$DATE, times = length(cols)),
  path_id = rep(seq_along(cols), each = nrow(df)),
  tmean_c = as.vector(temp_matrix[, cols])
)

write_csv(paths_sample, out_paths_csv)
cat("Wrote:", out_paths_csv, "\n")

# -----------------------------
# 13) PLOTS
# -----------------------------

# 13a) Bloom histogram (after model error)
hist_df <- tibble(bloom_doy = bloom_doy_adj)

p_hist <- ggplot(hist_df, aes(x = bloom_doy)) +
  geom_histogram(binwidth = 1, boundary = 0) +
  geom_vline(xintercept = as.numeric(qs), linetype = "dashed") +
  coord_cartesian(xlim = c(70, 115)) +
  labs(
    title = "Ensemble peak bloom forecast (with model error)",
    subtitle = sprintf("N=%d | Q10=%.1f | Q50=%.1f | Q90=%.1f", nrow(hist_df), qs[[1]], qs[[2]], qs[[3]]),
    x = "Predicted peak bloom day-of-year",
    y = "Count"
  )

print(p_hist)

# 13b) Daily temperature uncertainty ribbon from ALL sims (q10–q90) + median
sim_daily_q <- tibble(
  DATE = df$DATE,
  sim_q10 = apply(temp_matrix, 1, quantile, probs = 0.10, na.rm = TRUE, type = 7),
  sim_q50 = apply(temp_matrix, 1, quantile, probs = 0.50, na.rm = TRUE, type = 7),
  sim_q90 = apply(temp_matrix, 1, quantile, probs = 0.90, na.rm = TRUE, type = 7)
) %>%
  left_join(df %>% select(DATE, clim_mean, clim_q10, clim_q90, accu_tavg), by = "DATE") %>%
  mutate(
    sim_width  = sim_q90 - sim_q10,
    clim_width = clim_q90 - clim_q10
  )

p_temp <- ggplot() +
  geom_ribbon(data = sim_daily_q,
              aes(x = DATE, ymin = sim_q10, ymax = sim_q90),
              fill = "grey60", alpha = 0.35) +
  geom_line(data = sim_daily_q, aes(x = DATE, y = sim_q50),
            color = "black", linewidth = 1) +
  geom_ribbon(data = sim_daily_q,
              aes(x = DATE, ymin = clim_q10, ymax = clim_q90),
              fill = "blue", alpha = 0.12) +
  geom_line(data = sim_daily_q, aes(x = DATE, y = clim_mean),
            color = "blue", linewidth = 1) +
  geom_line(data = sim_daily_q %>% filter(is.finite(accu_tavg)),
            aes(x = DATE, y = accu_tavg),
            color = "red", linewidth = 1) +
  geom_vline(xintercept = today, linetype = "dashed") +
  coord_cartesian(xlim = c(plot_start, plot_end)) +
  labs(
    title = "Projected daily temperatures (ensemble q10–q90)",
    subtitle = "Grey band = simulated q10–q90 | Black = simulated median | Blue band/line = climatology q10–q90/mean | Red = Accu",
    x = NULL,
    y = "Daily mean temperature (°C)"
  )

print(p_temp)

# 13c) Width over time (sim vs climo)
p_width <- ggplot(sim_daily_q, aes(x = DATE)) +
  geom_line(aes(y = sim_width), color = "black", linewidth = 1) +
  geom_line(aes(y = clim_width), color = "blue", linewidth = 1) +
  geom_vline(xintercept = today, linetype = "dashed") +
  coord_cartesian(xlim = c(plot_start, plot_end)) +
  labs(
    title = "Temperature uncertainty width over time (q90 - q10)",
    subtitle = "Black = simulated | Blue = climatology reference",
    x = NULL,
    y = "Temperature spread (°C)"
  )

# ggsave(out_width_png, p_width, width = 10, height = 4.5, dpi = 150)
# cat("Wrote:", out_width_png, "\n")

# -----------------------------
# Final console summary
# -----------------------------
cat("\nBloom date probability mass check (should sum to ~100):",
    sum(prob_tbl$prob_percent), "\n")
cat("Mode (most likely date):",
    as.character(prob_tbl$DATE[which.max(prob_tbl$prob_percent)]),
    sprintf("(%.2f%%)", max(prob_tbl$prob_percent)), "\n")