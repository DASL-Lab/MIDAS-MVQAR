
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(MASS)
library(tidyr)
library(scales)
library(fBasics)
library(urca)
library(tseries)

# Load dataset
df <- read.csv(file.choose())

# Convert 'Date' to proper Date format and set as row names
df$Date <- parse_date_time(df$Date, orders = "b-y")
df <- df %>% arrange(Date)
# Check the correlation among the variables
var_cols <- vapply(df, is.numeric, logical(1))
var_df <- df[, var_cols, drop = FALSE]
corr_pearson <- cor(var_df, method = "pearson")

# Re-create event period data
event_periods <- data.frame(
  event = c("Global Financial Crisis", "Oil Price Crash", 
            "COVID-19 Pandemic", "Post-COVID Inflation Surge"),
  start = as.Date(c("2007-12-01", "2014-06-01", "2020-03-01", "2021-01-01")),
  end = as.Date(c("2009-06-30", "2016-01-31", "2021-12-31", "2023-12-31"))
)

# Reshape data into long format
df_long <- df %>%
  tidyr::pivot_longer(
    cols = -Date,
    names_to = "Variable",
    values_to = "Index"
  )

# Ensure df_long$Date is Date
df_long <- df_long %>%
  mutate(Date = as.Date(Date))  

desired_order <- c("CPI", "IR", "NHPI", "GSCPI", "IPPI", "EPI")
df_long$Variable <- factor(df_long$Variable, levels = desired_order)

# Ensure event_periods$start/end are Date
event_periods <- event_periods %>%
  mutate(
    start = as.Date(start), 
    end   = as.Date(end)
  )

# Plot
event_plot <- ggplot(df_long, aes(x = Date, y = Index)) +
  # Add event shading
  geom_rect(data = event_periods, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = event),
            inherit.aes = FALSE, alpha = 0.28) +
  geom_line(color = "blue") +
  facet_wrap(~ Variable, ncol = 2, dir = "v", scales = "free_y") +  # <<< key line
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years")) +
  labs(title = "Economic Indicators Over Time with Global Volatility Events Highlighted",
       x = "Date", y = "Index", fill = "Event") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(face = "bold", size = 10),
    strip.text   = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave("time_series_event_plot_no_stand.png", plot = event_plot, width = 10, 
       height = 8, dpi = 600, limitsize = FALSE)


## Descriptive statistics (1998-2025)
des_stat <- df %>% 
  dplyr::select(-Date)

#Test for stationarity using ADF and KPSS tests
for (col in colnames(des_stat)) {
  cat("ADF test for", col, ":\n")
  print(summary(ur.df(des_stat[,col], type = "trend", lags = 12, selectlags = "BIC")))
  cat("KPSS test for", col, ":\n")
  print(summary(ur.kpss(des_stat[,col], type = "tau", lags = "short")))
}

# Taking differencing (Low frequency variables) for stationarity
cpi_lf_diff   <- diff(des_stat$CPI, 1)
nhpi_lf_diff  <- diff(des_stat$NHPI, 1)
ir_lf_diff  <- diff(des_stat$IR, 1)   

# Create the low frequency timeseries using the stationary data
lf_stationary <- na.omit(data.frame(
  CPI  = cpi_lf_diff,
  IR   = ir_lf_diff,
  NHPI = nhpi_lf_diff
))

# # Taking differencing (high frequency variables) for stationarity 
gscpi_hf_diff <- diff(des_stat$GSCPI, 1)
ippi_hf_diff <- diff(des_stat$IPPI, 1)  
epi_hf_diff  <- diff(des_stat$EPI, 1) 


# Create the high frequency timeseries using the stationary data
hf_stationary <- na.omit(data.frame(
  GSCPI = gscpi_hf_diff,
  IPPI  = ippi_hf_diff,
  EPI   = epi_hf_diff
))

#Test for stationarity after diferencing using ADF and KPSS tests (low frequency variables)
for (col in colnames(lf_stationary)) {
  cat("ADF test for", col, ":\n")
  print(summary(ur.df(lf_stationary[,col], type = "trend", lags = 12, selectlags = "BIC")))
  cat("KPSS test for", col, ":\n")
  print(summary(ur.kpss(lf_stationary[,col], type = "tau", lags = "short")))
}

#Test for stationarity after diferencing using ADF and KPSS tests (high frequency variables)
for (col in colnames(hf_stationary)) {
  cat("ADF test for", col, ":\n")
  print(summary(ur.df(hf_stationary[,col], type = "trend", lags = 12, selectlags = "BIC")))
  cat("KPSS test for", col, ":\n")
  print(summary(ur.kpss(hf_stationary[,col], type = "tau", lags = "short")))
}


## Data after differencing for stationarity
Date <- df$Date[-1]
df_diff <- cbind(Date, lf_stationary, hf_stationary)

# Reshape data into long format
df_long_diff <- df_diff %>%
  tidyr::pivot_longer(
    cols = -Date,
    names_to = "Variable",
    values_to = "Index"
  )

# Ensure df_long$Date is Date
df_long_diff <- df_long_diff %>%
  mutate(Date = as.Date(Date))  

df_long_diff$Variable <- factor(df_long_diff$Variable, levels = desired_order)

# Plot after differencing
event_plot_diff <- ggplot(df_long_diff, aes(x = Date, y = Index)) +
  # Add event shading
  geom_rect(data = event_periods, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = event),
            inherit.aes = FALSE, alpha = 0.28) +
  geom_line(color = "blue") +
  facet_wrap(~ Variable, ncol = 2, dir = "v", scales = "free_y") +  # <<< key line
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("3 years")) +
  labs(title = "Differenced Economic Indicators Over Time with Global Volatility Events Highlighted",
       x = "Date", y = "Index", fill = "Event") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", size = 10, angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(face = "bold", size = 10),
    strip.text   = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave("diff_time_series_event_plot_new.png", plot = event_plot_diff, width = 10, 
       height = 8, dpi = 600, limitsize = FALSE)


# Create time series objects for each variable starting from Jan 1998
ts_cpi <- ts(lf_stationary$CPI, start = c(1998, 2), frequency = 12)
ts_ir <- ts(lf_stationary$IR, start = c(1998, 2), frequency = 12)
ts_nhpi <- ts(lf_stationary$NHPI, start = c(1998, 2), frequency = 12)
ts_gscpi <- ts(hf_stationary$GSCPI, start = c(1998, 2), frequency = 12)
ts_ippi <- ts(hf_stationary$IPPI, start = c(1998, 2), frequency = 12)
ts_epi <- ts(hf_stationary$EPI, start = c(1998, 2), frequency = 12)

# Create a list of all time series for easier processing
ts_list <- list(
  diff_CPI = ts_cpi,
  diff_IR = ts_ir,
  diff_NHPI = ts_nhpi,
  diff_GSCPI = ts_gscpi,
  diff_IPPI = ts_ippi,
  diff_EPI = ts_epi
)


# Function to analyze a each time series
analyze_ts <- function(ts_data, var_name, out_dir = "figs") {
  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("ANALYSIS FOR:", var_name, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  file_3panel <- file.path(out_dir, paste0("ts_acf_pacf_", var_name, ".png"))
  png(
    filename = file_3panel,
    width    = 8,      # inches
    height   = 6,      # inches
    units    = "in",
    res      = 600,   
    type     = "cairo-png"
  )
  
  par(mfrow = c(2, 2))
  
  # Original time series
  plot(ts_data,
       main = paste("Time Series Plot of", var_name),
       ylab = var_name,
       xlab = "Time")
  grid()
  
  # ACF
  acf(as.numeric(ts_data),
      main    = paste("ACF of", var_name),
      lag.max = 36,
      na.action = na.pass,
      xlab    = "Lag (Months)")
  
  # PACF
  pacf(as.numeric(ts_data),
       main    = paste("PACF of", var_name),
       lag.max = 36,
       na.action = na.pass,
       xlab    = "Lag (Months)")
  
  plot.new()
  title(main = "")
  
  dev.off() 
  par(mfrow = c(1, 1))  # reset layout
  
 # SUMMARY STATISTICS
  cat("\nSUMMARY STATISTICS:\n")
  cat("Mean:", mean(ts_data, na.rm = TRUE), "\n")
  cat("Standard Deviation:", sd(ts_data, na.rm = TRUE), "\n")
  cat("Minimum:", min(ts_data, na.rm = TRUE), "\n")
  cat("Maximum:", max(ts_data, na.rm = TRUE), "\n")
}

# Analyze each time series
for (var_name in names(ts_list)) {
  analyze_ts(ts_list[[var_name]], var_name)
}

