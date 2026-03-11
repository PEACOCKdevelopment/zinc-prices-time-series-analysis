# Zinc Prices Time-Series Analysis

End-to-end econometric analysis and forecasting of monthly zinc prices (1990-2025) in R, with a reproducible pipeline and publication-ready report.

This project combines data quality checks, transformation diagnostics, stationarity testing, ARIMA model selection, and holdout forecast evaluation against an extrapolative Holt benchmark.

## Why This Project Matters

Commodity prices are noisy, cyclical, and difficult to forecast. This repository demonstrates a transparent statistical workflow for:

- validating assumptions before model fitting,
- comparing multiple ARIMA specifications with consistent metrics,
- benchmarking against a simpler extrapolative method,
- generating reusable artifacts (tables and figures) for reporting.

## RPubs Report

Full published report:

- https://rpubs.com/timothypawelczyk/1408744

## Analytical Scope

- Frequency: monthly
- Horizon: January 1990 to March 2025
- Target: zinc price (USD per metric ton)
- Main model family: ARIMA
- Benchmark model: Holt exponential smoothing (damped)

## Methodology Overview

1. Data ingestion and validation
2. Log transformation and decomposition diagnostics
3. Stationarity checks (ADF, ur.df) and Box-Cox assessment
4. ARIMA candidate fitting and ranking (AIC, BIC, holdout errors)
5. Residual diagnostics (including Ljung-Box)
6. Benchmark comparison and 12-month forward forecast generation

## Repository Structure

- zinc_prices.csv: source time-series data
- zinc_prices_pipeline.R: reproducible modeling pipeline that writes outputs to outputs/
- zinc_prices_rpubs_report.Rmd: full analytical narrative and diagnostics
- zinc_prices_rpubs_report.html: knitted HTML report
- outputs/: generated CSV tables and figures used for interpretation and reporting

## Reproducibility

### Requirements

- R (recommended: 4.5.x)
- Installed packages:
	readr, dplyr, tidyr, ggplot2, forecast, tseries, urca, purrr, tibble, scales, MASS, lmtest, knitr, rmarkdown

### Run The Pipeline

```bash
Rscript zinc_prices_pipeline.R
```

### Knit The Report

```bash
Rscript -e "rmarkdown::render('zinc_prices_rpubs_report.Rmd', output_format = 'html_document')"
```

## Main Outputs

The pipeline exports a complete artifact set into outputs/, including:

- quality checks and decomposition diagnostics,
- Box-Cox and stationarity tables,
- ARIMA candidate ranking and selected model diagnostics,
- final 12-month forecast table,
- benchmark comparison tables in both log and level scales,
- diagnostic and forecast figures (ACF/PACF, decomposition, holdout fit, forecast intervals).

## Key Findings


- The series shows strong cyclical behavior with weak seasonal influence.
- After transformation/differencing, ARIMA assumptions are better supported.
- Among tested candidates, ARIMA(4,1,1) delivers the strongest overall forecast performance in the report.
- ARIMA outperforms the Holt benchmark on holdout accuracy metrics for this dataset.

## Author

Tymoteusz Pawelczyk
