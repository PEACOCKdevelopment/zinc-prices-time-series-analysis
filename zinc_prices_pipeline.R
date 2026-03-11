# Zinc Price Time Series Pipeline
#
# This script implements an end-to-end, reproducible workflow for monthly zinc
# price modeling and forecasting using ARIMA and benchmark methods.
#
# Outputs:
# - Clean data quality table
# - Transformation and stationarity diagnostics
# - Candidate model comparison and benchmark metrics
# - Holdout and production forecast plots
# - CSV artifacts for downstream reporting or versioning

# -----------------------------
# Configuration
# -----------------------------

get_project_dir <- function() {
  # Resolve project directory from --file when run via Rscript; otherwise use cwd.
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }

  getwd()
}

PROJECT_DIR <- get_project_dir()
INPUT_CSV <- file.path(PROJECT_DIR, "zinc_prices.csv")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
HOLDOUT_H <- 12

required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "forecast",
  "tseries", "urca", "purrr", "tibble", "scales"
)

# -----------------------------
# Utility Functions
# -----------------------------

ensure_packages <- function(packages) {
  missing <- setdiff(packages, rownames(installed.packages()))

  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing required packages: ",
        paste(missing, collapse = ", "),
        ". Install them before running this script."
      )
    )
  }

  invisible(lapply(packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

calc_metrics <- function(actual, predicted) {
  # Add small epsilons for numerical stability in percentage metrics.
  eps <- 1e-8

  tibble::tibble(
    MAE = mean(abs(actual - predicted)),
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAPE = mean(abs(actual - predicted) / pmax(abs(actual), eps)) * 100,
    sMAPE = mean(200 * abs(actual - predicted) / (abs(actual) + abs(predicted) + eps))
  )
}

extract_ur_row <- function(model, label) {
  test_name <- names(model@teststat)[1]

  tibble::tibble(
    series = label,
    test = test_name,
    statistic = unname(model@teststat[1]),
    critical_1pct = unname(model@cval[1, "1pct"]),
    critical_5pct = unname(model@cval[1, "5pct"]),
    critical_10pct = unname(model@cval[1, "10pct"]),
    selected_lags = model@lags
  )
}

save_acf_pacf_plot <- function(series, path_png) {
  # Save ACF/PACF diagnostics in a single image using base plotting.
  png(path_png, width = 1400, height = 600)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(mfrow = c(1, 2))
  acf(series, lag.max = 48, main = "ACF: diff(log(price))")
  pacf(series, lag.max = 48, main = "PACF: diff(log(price))")
}

save_stl_components_plot <- function(stl_fit, path_png) {
  # Save STL decomposition components using base plotting.
  png(path_png, width = 1400, height = 800)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  plot(stl_fit, main = "STL decomposition components")
}

save_boxcox_profile_plot <- function(price_series, lambda_guerrero, path_png) {
  # Save profile log-likelihood curve for Box-Cox transformation.
  price_clean <- stats::na.omit(as.numeric(price_series))
  time_index <- seq_along(price_clean)
  boxcox_df <- data.frame(price_clean = price_clean, time_index = time_index)

  png(path_png, width = 1400, height = 800)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  MASS::boxcox(
    price_clean ~ time_index,
    data = boxcox_df,
    lambda = seq(-1, 1, 0.1)
  )
  graphics::title(main = "Figure 3. Box-Cox profile for raw zinc price")
  graphics::abline(v = lambda_guerrero, lty = 2, col = "firebrick")
}

# -----------------------------
# Pipeline Steps
# -----------------------------

load_and_validate_data <- function(csv_path) {
  raw_data <- readr::read_csv(csv_path, show_col_types = FALSE)

  data <- raw_data %>%
    dplyr::transmute(
      date = as.Date(paste0("01-", Date), format = "%d-%b-%y"),
      price = as.numeric(gsub(",", ".", gsub("\\s+", "", Price)))
    ) %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(log_price = log(price))

  if (anyNA(data$date) || anyNA(data$price)) {
    stop("Date or price parsing generated missing values.")
  }
  if (any(data$price <= 0)) {
    stop("Non-positive prices detected. Log transformation is invalid.")
  }
  if (dplyr::n_distinct(data$date) != nrow(data)) {
    stop("Duplicate monthly timestamps detected.")
  }

  data
}

build_time_series <- function(data) {
  start_year <- as.integer(format(min(data$date), "%Y"))
  start_month <- as.integer(format(min(data$date), "%m"))

  price_ts <- ts(data$price, start = c(start_year, start_month), frequency = 12)
  log_price_ts <- ts(data$log_price, start = c(start_year, start_month), frequency = 12)
  diff_log_price_ts <- diff(log_price_ts)

  list(
    price_ts = price_ts,
    log_price_ts = log_price_ts,
    diff_log_price_ts = diff_log_price_ts
  )
}

compute_qa_table <- function(data) {
  tibble::tibble(
    metric = c(
      "Observations",
      "Start date",
      "End date",
      "Missing dates",
      "Missing prices",
      "Duplicate dates"
    ),
    value = c(
      nrow(data),
      as.character(min(data$date)),
      as.character(max(data$date)),
      sum(is.na(data$date)),
      sum(is.na(data$price)),
      nrow(data) - dplyr::n_distinct(data$date)
    )
  )
}

compute_decomposition <- function(log_price_ts, data) {
  stl_fit <- stl(log_price_ts, s.window = "periodic", robust = TRUE)

  sa_log <- forecast::seasadj(stl_fit)
  sa_price <- exp(as.numeric(sa_log))

  decomp_df <- tibble::tibble(
    date = data$date,
    original_price = data$price,
    sa_price = sa_price
  )

  seasonality_strength <- max(
    0,
    1 - stats::var(stl_fit$time.series[, "remainder"], na.rm = TRUE) /
      stats::var(stl_fit$time.series[, "remainder"] + stl_fit$time.series[, "seasonal"], na.rm = TRUE)
  )

  decomp_summary <- tibble::tibble(
    metric = c(
      "Seasonality strength (0 to 1)",
      "Correlation original vs seasonally adjusted"
    ),
    value = c(
      round(seasonality_strength, 4),
      round(stats::cor(decomp_df$original_price, decomp_df$sa_price), 4)
    )
  )

  list(
    stl_fit = stl_fit,
    decomp_df = decomp_df,
    decomp_summary = decomp_summary,
    seasonality_strength = seasonality_strength
  )
}

compute_stationarity_tables <- function(price_ts, log_price_ts, diff_log_price_ts) {
  lambda_guerrero <- forecast::BoxCox.lambda(price_ts, method = "guerrero")

  boxcox_table <- tibble::tibble(
    diagnostic = "Guerrero Box-Cox lambda",
    value = round(lambda_guerrero, 3)
  )

  adf_level <- tseries::adf.test(log_price_ts, k = 12)
  adf_diff <- tseries::adf.test(diff_log_price_ts, k = 12)

  adf_table <- tibble::tibble(
    series = c("log(price)", "diff(log(price))"),
    statistic = c(unname(adf_level$statistic), unname(adf_diff$statistic)),
    p_value = c(adf_level$p.value, adf_diff$p.value)
  )

  ur_level <- urca::ur.df(log_price_ts, type = "trend", selectlags = "AIC")
  ur_diff <- urca::ur.df(diff_log_price_ts, type = "drift", selectlags = "AIC")

  ur_table <- dplyr::bind_rows(
    extract_ur_row(ur_level, "log(price)"),
    extract_ur_row(ur_diff, "diff(log(price))")
  )

  list(
    lambda_guerrero = lambda_guerrero,
    boxcox_table = boxcox_table,
    adf_table = adf_table,
    ur_table = ur_table,
    adf_level = adf_level,
    adf_diff = adf_diff,
    ur_level = ur_level,
    ur_diff = ur_diff
  )
}

build_train_test <- function(log_price_ts, full_dates, h) {
  n <- length(log_price_ts)
  if (n <= h) {
    stop("Not enough observations for holdout split.")
  }

  train_log <- ts(
    head(as.numeric(log_price_ts), n - h),
    start = start(log_price_ts),
    frequency = 12
  )

  test_log <- tail(as.numeric(log_price_ts), h)
  test_dates <- tail(full_dates, h)

  list(
    train_log = train_log,
    test_log = test_log,
    test_dates = test_dates
  )
}

fit_candidate_models <- function(train_log, test_log, h) {
  candidate_orders <- tibble::tibble(
    model = c("ARIMA(1,1,1)", "ARIMA(1,1,2)", "ARIMA(4,1,1)", "ARIMA(4,1,2)"),
    p = c(1, 1, 4, 4),
    d = c(1, 1, 1, 1),
    q = c(1, 2, 1, 2)
  )

  safe_arima <- purrr::possibly(
    function(p, d, q) {
      forecast::Arima(train_log, order = c(p, d, q), include.drift = TRUE)
    },
    otherwise = NULL
  )

  model_results <- candidate_orders %>%
    dplyr::mutate(fit = purrr::pmap(list(p, d, q), safe_arima)) %>%
    dplyr::filter(!purrr::map_lgl(fit, is.null)) %>%
    dplyr::mutate(
      fc = purrr::map(fit, forecast::forecast, h = h),
      metrics = purrr::map(fc, ~calc_metrics(test_log, as.numeric(.x$mean))),
      AIC = purrr::map_dbl(fit, AIC),
      BIC = purrr::map_dbl(fit, BIC),
      logLik = purrr::map_dbl(fit, ~as.numeric(logLik(.x))),
      Ljung_Box_p = purrr::map_dbl(
        fit,
        ~Box.test(
          residuals(.x),
          lag = 24,
          type = "Ljung-Box",
          fitdf = length(coef(.x))
        )$p.value
      )
    ) %>%
    tidyr::unnest(metrics) %>%
    dplyr::arrange(RMSE, AIC)

  if (nrow(model_results) == 0) {
    stop("No candidate ARIMA model converged.")
  }

  best_row <- dplyr::slice(model_results, 1)

  list(
    candidate_orders = candidate_orders,
    model_results = model_results,
    best_model_name = best_row$model[[1]],
    best_fit = best_row$fit[[1]],
    best_fc = best_row$fc[[1]]
  )
}

compute_benchmarks <- function(train_log, test_log, best_model_name, best_fc, h) {
  auto_fit <- forecast::auto.arima(train_log, seasonal = FALSE)
  auto_fc <- forecast::forecast(auto_fit, h = h)

  holt_fit <- forecast::holt(train_log, h = h, damped = TRUE, initial = "optimal")

  benchmark_log <- dplyr::bind_rows(
    calc_metrics(test_log, as.numeric(best_fc$mean)) %>% dplyr::mutate(model = best_model_name),
    calc_metrics(test_log, as.numeric(auto_fc$mean)) %>% dplyr::mutate(model = "auto.arima"),
    calc_metrics(test_log, as.numeric(holt_fit$mean)) %>% dplyr::mutate(model = "Holt (damped)")
  ) %>%
    dplyr::select(model, MAE, RMSE, MAPE, sMAPE) %>%
    dplyr::arrange(RMSE)

  test_level <- exp(test_log)
  benchmark_level <- dplyr::bind_rows(
    calc_metrics(test_level, exp(as.numeric(best_fc$mean))) %>% dplyr::mutate(model = best_model_name),
    calc_metrics(test_level, exp(as.numeric(auto_fc$mean))) %>% dplyr::mutate(model = "auto.arima"),
    calc_metrics(test_level, exp(as.numeric(holt_fit$mean))) %>% dplyr::mutate(model = "Holt (damped)")
  ) %>%
    dplyr::select(model, MAE, RMSE, MAPE, sMAPE) %>%
    dplyr::arrange(RMSE)

  list(
    auto_fit = auto_fit,
    auto_fc = auto_fc,
    holt_fit = holt_fit,
    benchmark_log = benchmark_log,
    benchmark_level = benchmark_level
  )
}

fit_final_forecast <- function(log_price_ts, candidate_orders, best_model_name, h) {
  best_order <- candidate_orders %>%
    dplyr::filter(model == best_model_name) %>%
    dplyr::slice(1)

  final_fit <- forecast::Arima(
    log_price_ts,
    order = c(best_order$p, best_order$d, best_order$q),
    include.drift = TRUE
  )

  future_fc <- forecast::forecast(final_fit, h = h, level = c(80, 95))

  list(
    best_order = best_order,
    final_fit = final_fit,
    future_fc = future_fc
  )
}

build_forecast_table <- function(data, future_fc) {
  last_month <- as.Date(format(max(data$date), "%Y-%m-01"))
  future_dates <- seq(from = last_month, by = "month", length.out = 13)[-1]

  tibble::tibble(
    date = future_dates,
    point_forecast = exp(as.numeric(future_fc$mean)),
    lo80 = exp(as.numeric(future_fc$lower[, "80%"])),
    hi80 = exp(as.numeric(future_fc$upper[, "80%"])),
    lo95 = exp(as.numeric(future_fc$lower[, "95%"])),
    hi95 = exp(as.numeric(future_fc$upper[, "95%"]))
  )
}

save_plots <- function(data, diff_log_price_ts, decomp_df, stl_fit, holdout_df, forecast_table, best_model_name, lambda_guerrero) {
  plot_df <- data %>%
    dplyr::select(date, price, log_price) %>%
    tidyr::pivot_longer(cols = c(price, log_price), names_to = "series", values_to = "value") %>%
    dplyr::mutate(
      series = dplyr::recode(
        series,
        price = "Price (USD/ton)",
        log_price = "log(Price)"
      )
    )

  fig1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = date, y = value)) +
    ggplot2::geom_line(color = "#1f78b4", linewidth = 0.8) +
    ggplot2::facet_wrap(~series, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = "Figure 1. Zinc price and log-price over time",
      x = "Date",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)

  fig2 <- ggplot2::ggplot(decomp_df, ggplot2::aes(x = date)) +
    ggplot2::geom_line(ggplot2::aes(y = original_price, color = "Original"), linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = sa_price, color = "Seasonally adjusted"), linewidth = 0.8, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("Original" = "#1f78b4", "Seasonally adjusted" = "#e31a1c")) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title = "Figure 2. Original vs seasonally adjusted series",
      x = "Date",
      y = "Price (USD)",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)

  fig4 <- ggplot2::ggplot(holdout_df, ggplot2::aes(x = date, y = value, color = series)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title = "Figure 4. Holdout performance in USD levels",
      x = "Date",
      y = "Price (USD)",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)

  history_tail <- data %>%
    dplyr::filter(date >= min(holdout_df$date)) %>%
    dplyr::select(date, price)

  fig5 <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = history_tail,
      ggplot2::aes(x = date, y = price),
      color = "#1f78b4",
      linewidth = 0.9
    ) +
    ggplot2::geom_ribbon(
      data = forecast_table,
      ggplot2::aes(x = date, ymin = lo95, ymax = hi95),
      fill = "#a6cee3",
      alpha = 0.25
    ) +
    ggplot2::geom_ribbon(
      data = forecast_table,
      ggplot2::aes(x = date, ymin = lo80, ymax = hi80),
      fill = "#1f78b4",
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      data = forecast_table,
      ggplot2::aes(x = date, y = point_forecast),
      color = "#e31a1c",
      linewidth = 1
    ) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title = paste("Figure 5. 12-month forecast using", best_model_name),
      x = "Date",
      y = "Price (USD)"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  ggplot2::ggsave(file.path(OUTPUT_DIR, "figure_1_price_log.png"), fig1, width = 10, height = 6, dpi = 120)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "figure_2_decomposition_overlay.png"), fig2, width = 10, height = 5, dpi = 120)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "figure_4_holdout_levels.png"), fig4, width = 10, height = 5, dpi = 120)
  ggplot2::ggsave(file.path(OUTPUT_DIR, "figure_5_final_forecast.png"), fig5, width = 10, height = 5, dpi = 120)

  save_boxcox_profile_plot(data$price, lambda_guerrero, file.path(OUTPUT_DIR, "figure_3_boxcox_profile.png"))
  save_stl_components_plot(stl_fit, file.path(OUTPUT_DIR, "figure_2b_stl_components.png"))
  save_acf_pacf_plot(diff_log_price_ts, file.path(OUTPUT_DIR, "figure_3_acf_pacf.png"))
}

# -----------------------------
# Main
# -----------------------------

main <- function() {
  ensure_packages(required_packages)

  if (!file.exists(INPUT_CSV)) {
    stop("Input file not found: ", INPUT_CSV)
  }

  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }

  message("[1/8] Loading and validating data...")
  data <- load_and_validate_data(INPUT_CSV)

  message("[2/8] Building time series and QA artifacts...")
  ts_bundle <- build_time_series(data)
  qa_table <- compute_qa_table(data)

  message("[3/8] Running decomposition diagnostics...")
  decomp <- compute_decomposition(ts_bundle$log_price_ts, data)

  message("[4/8] Running stationarity diagnostics...")
  stationarity <- compute_stationarity_tables(
    ts_bundle$price_ts,
    ts_bundle$log_price_ts,
    ts_bundle$diff_log_price_ts
  )

  message("[5/8] Fitting ARIMA candidate models...")
  split_bundle <- build_train_test(ts_bundle$log_price_ts, data$date, HOLDOUT_H)
  candidates <- fit_candidate_models(
    split_bundle$train_log,
    split_bundle$test_log,
    HOLDOUT_H
  )

  message("[6/8] Computing benchmark models...")
  benchmarks <- compute_benchmarks(
    split_bundle$train_log,
    split_bundle$test_log,
    candidates$best_model_name,
    candidates$best_fc,
    HOLDOUT_H
  )

  message("[7/8] Fitting production forecast...")
  final_model <- fit_final_forecast(
    ts_bundle$log_price_ts,
    candidates$candidate_orders,
    candidates$best_model_name,
    HOLDOUT_H
  )
  forecast_table <- build_forecast_table(data, final_model$future_fc)

  plot_start_date <- if (nrow(data) > 120) {
    data$date[nrow(data) - 119]
  } else {
    min(data$date)
  }

  holdout_df <- dplyr::bind_rows(
    data %>%
      dplyr::slice(1:(nrow(data) - HOLDOUT_H)) %>%
      dplyr::transmute(date, value = price, series = "Actual (train)"),
    data %>%
      dplyr::slice((nrow(data) - HOLDOUT_H + 1):n()) %>%
      dplyr::transmute(date, value = price, series = "Actual (test)"),
    tibble::tibble(
      date = split_bundle$test_dates,
      value = exp(as.numeric(candidates$best_fc$mean)),
      series = paste0("Forecast ", candidates$best_model_name)
    )
  ) %>%
    dplyr::filter(date >= plot_start_date)

  message("[8/8] Writing CSV outputs and saving figures...")

  readr::write_csv(qa_table, file.path(OUTPUT_DIR, "table_1_data_quality.csv"))
  readr::write_csv(decomp$decomp_summary, file.path(OUTPUT_DIR, "table_2_decomposition_diagnostics.csv"))
  readr::write_csv(stationarity$boxcox_table, file.path(OUTPUT_DIR, "table_3_boxcox.csv"))
  readr::write_csv(stationarity$adf_table, file.path(OUTPUT_DIR, "table_4_adf.csv"))
  readr::write_csv(stationarity$ur_table, file.path(OUTPUT_DIR, "table_5_urdf.csv"))

  readr::write_csv(
    candidates$model_results %>% dplyr::select(model, p, d, q, AIC, BIC, logLik, MAE, RMSE, MAPE, sMAPE, Ljung_Box_p),
    file.path(OUTPUT_DIR, "table_6_arima_candidates.csv")
  )

  best_lb <- Box.test(
    residuals(candidates$best_fit),
    lag = 24,
    type = "Ljung-Box",
    fitdf = length(coef(candidates$best_fit))
  )

  best_model_table <- tibble::tibble(
    selected_model = candidates$best_model_name,
    holdout_rmse_log = round(candidates$model_results$RMSE[1], 4),
    holdout_mape_log = round(candidates$model_results$MAPE[1], 4),
    residual_ljung_box_p = round(best_lb$p.value, 4)
  )

  readr::write_csv(best_model_table, file.path(OUTPUT_DIR, "table_7_selected_model.csv"))
  readr::write_csv(forecast_table, file.path(OUTPUT_DIR, "table_8_final_forecast.csv"))
  readr::write_csv(benchmarks$benchmark_log, file.path(OUTPUT_DIR, "table_9a_benchmark_log.csv"))
  readr::write_csv(benchmarks$benchmark_level, file.path(OUTPUT_DIR, "table_9b_benchmark_level.csv"))

  save_plots(
    data = data,
    diff_log_price_ts = ts_bundle$diff_log_price_ts,
    decomp_df = decomp$decomp_df,
    stl_fit = decomp$stl_fit,
    holdout_df = holdout_df,
    forecast_table = forecast_table,
    best_model_name = candidates$best_model_name,
    lambda_guerrero = stationarity$lambda_guerrero
  )

  message("Done. Outputs saved in: ", OUTPUT_DIR)
  message("Selected model: ", candidates$best_model_name)
}

if (sys.nframe() == 0) {
  main()
}
