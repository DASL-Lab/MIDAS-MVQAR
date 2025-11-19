library(xts)
library(zoo)
library(quantreg)
library(midasr)

source("functions.R")

##Set Parameters in an easy place to change
m <- 3
K <- 4
p <- 1
taus <- c(0.25, 0.50, 0.75, 0.90)


start_time <- Sys.time()

## Simulate Data
model_data_sim_1 <- simulate(m = m, K = K, p = p, taus = taus)
n_obs <- nrow(model_data_sim_1)

invisible(lapply(model_data_sim_1[c("cpi","ir","nhpi")], check_mid_gaps))       # Mid-series NA gap check (like midas_qr warning)

true_par <- create_true_par_vector(p)

cat("=== DATA SUMMARY ===\n")
cat("Sample size:", n_obs, "quarters\n")
cat("Total parameters:", length(true_par), "\n")
cat("Parameters per observation:", round(n_obs / length(true_par), 2), "\n")


##Run Estimation with Analytic Gradients
# First verify gradients work correctly
#cat("=== VERIFYING ANALYTIC GRADIENTS ===\n")
#par_test <- true_par + rnorm(length(true_par), sd = 0.1)
#grad_check <- verify_gradients(par_test, model_data_sim_1[1:30, ], tau = 0.5, K, p)

cat("\nStarting estimation with analytic gradients...\n")
fits_sim_1 <- estimate_joint_qarmidas(
  taus = taus, data = model_data_sim_1, K = K, p = p,
  kappa = 0.01,
  theta_shrink = 2e-3,     
  start_par = true_par,
  optimizer = "BFGS",
  use_bounds = FALSE,
  n_restarts = 2,
  jitter_sd = 0.01,
  verbose = TRUE,
  dry_run = FALSE,
  use_analytic_grad = TRUE  # Use our analytic gradients!
)

params <- extract_params(fits_sim_1, p)

print_true_parameters(true_par, p)

cat("\n=== ESTIMATED PARAMETERS ===\n")
for (tau_name in names(params)) {
  cat("\nResults for", tau_name, ":\n")
  cat("Convergence:", params[[tau_name]]$convergence, " Loss:", round(params[[tau_name]]$loss, 6), "\n")
  cat("Equation 1 (CPI):\n"); print(round(params[[tau_name]]$eq1_params, 4))
  cat("\nEquation 2 (Interest Rate):\n"); print(round(params[[tau_name]]$eq2_params, 4))
  cat("\nEquation 3 (NHPI):\n"); print(round(params[[tau_name]]$eq3_params, 4))
  cat("Time taken:", round(params[[tau_name]]$elapsed_minutes, 2), "minutes\n")
}

## PLOTS
# Create fitted series per tau from joint fits_sim_1
tau_keys <- paste0("tau_", format(taus, nsmall = 2))
fitted_by_tau_cpi  <- setNames(lapply(tau_keys, function(k) fits_sim_1[[k]]$fitted_cpi), tau_keys)
fitted_by_tau_ir   <- setNames(lapply(tau_keys, function(k) fits_sim_1[[k]]$fitted_ir),  tau_keys)
fitted_by_tau_nhpi <- setNames(lapply(tau_keys, function(k) fits_sim_1[[k]]$fitted_nhpi),tau_keys)


par(mfrow = c(1, 1))#, mar = c(4, 4, 3, 1) + 0.1)
plot_actual_vs_fitted(model_data_sim_1$cpi,  fitted_by_tau_cpi,  ylab = "CPI",           main_prefix = "CPI - ")
plot_actual_vs_fitted(model_data_sim_1$ir,   fitted_by_tau_ir,   ylab = "Interest Rate", main_prefix = "Interest Rate - ")
plot_actual_vs_fitted(model_data_sim_1$nhpi, fitted_by_tau_nhpi, ylab = "NHPI",          main_prefix = "NHPI - ")
par(mfrow = c(1, 1))


diagnostic_results <- diagnose_parameter_recovery(params, true_par, p, abs_floor = 5e-3)
create_variable_plots(params, model_data_sim_1)

cat("\n=== SUMMARY TABLE (All Metrics in %) ===\n")
summary_table <- data.frame()
for (tau_name in names(params)) {
  est <- c(params[[tau_name]]$eq1_params, params[[tau_name]]$eq2_params, params[[tau_name]]$eq3_params)
  tru <- true_par
  ape <- abs(est - tru) / pmax(abs(tru), 5e-3) * 100
  summary_table <- rbind(summary_table, data.frame(
    Quantile = tau_name,
    Convergence = params[[tau_name]]$convergence,
    Loss = round(params[[tau_name]]$loss, 6),
    RMSE = sprintf("%.3f", sqrt(mean((est - tru)^2))),
    MAPE = sprintf("%.2f%%", mean(ape)),
    MdAPE = sprintf("%.2f%%", median(ape)),
    Time_Minutes = round(params[[tau_name]]$elapsed_minutes, 2)
  ))
}
print(summary_table, row.names = FALSE)

end_time <- Sys.time()
cat("\n=== SINGLE-FUNCTION JOINT ESTIMATION WITH ANALYTIC GRADIENTS COMPLETE ===\n")
cat("Total estimation time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
cat("Quantiles estimated:", paste(taus, collapse = ", "), "\n")


