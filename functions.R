

# Almon polynomial weight function 
almon_weight_function <- function(j, theta1, theta2, K, eps2 = 5e-3) {
  s_j <- if (K > 0) (j / K) else j                    # s_j = j/K, j=0,1,..,K. 
  theta2_eff <- -pmax(eps2, -theta2)                   # theta2_eff = -max (epsilon_2, -theta2)
  un_norm <- exp(theta1 * s_j + theta2_eff * s_j^2)    # Quadratic exponential kernel (Almon polynomial) (unnormalized_weight = exp(theta1*s_j + theta2_eff*s_j^2))
  u_norm <- pmax(un_norm, .Machine$double.eps)         # unnormalized_weight = max(u,.Machine$double.eps)
  norm_weight <- un_norm / sum(un_norm)               # normalized_weight = unnormalized_weight/sum(unnormalized_weight)
  norm_weight
}

check_mid_gaps <- function(x) {
  indices <- seq_along(x); ok <- !is.na(x); valid_indices <- indices[ok]
  if (length(valid_indices) > 1 && any(abs(diff(valid_indices) - 1) > 0))
    warning("There are NAs in the middle of the time series")
  invisible(NULL)
}

# Derivative of Almon weights with respect to theta parameters
almon_weight_gradient <- function(j, theta1, theta2, K, eps2 = 5e-3) {
  s_j <- if (K > 0) (j / K) else j
  theta2_eff <- -pmax(eps2, -theta2)
  
  # Unnormalized weight and its derivatives
  u_j <- exp(theta1 * s_j + theta2_eff * s_j^2)
  du_dtheta1 <- u_j * s_j
  
  # the theta2_eff transformation in derivative
  dtheta2_eff_dtheta2 <- ifelse(-theta2 > eps2, 1, 0)
  du_dtheta2 <- u_j * s_j^2 * dtheta2_eff_dtheta2
  
  # Sum for normalization
  U_sum <- sum(u_j)
  dU_sum_dtheta1 <- sum(du_dtheta1)
  dU_sum_dtheta2 <- sum(du_dtheta2)
  
  # Derivatives of normalized weights (quotient rule)
  dw_dtheta1 <- (du_dtheta1 * U_sum - u_j * dU_sum_dtheta1) / U_sum^2
  dw_dtheta2 <- (du_dtheta2 * U_sum - u_j * dU_sum_dtheta2) / U_sum^2
  
  return(list(weights = u_j / U_sum, 
              grad_theta1 = dw_dtheta1, 
              grad_theta2 = dw_dtheta2))
}

create_true_par_vector <- function(p) {
  n_params_eq1 <- 1 + 3*p + 9
  n_params_eq2 <- 2 + 3*p + 9
  n_params_eq3 <- 3 + 3*p + 9
  
  # Equation 1: CPI  
  par_eq1 <- c(
    eq1_alpha = 0.3,    #intercept 
    if (p > 0) {
      c(eq1_phi1 = 0.6, eq1_phi2 = 0.08, eq1_phi3 = -0.06) 
    } else {
      numeric(0)          # Cross-lags: [cpi_lag1, ir_lag1, nhpi_lag1] for p=1
    },
    
    # MIDAS scales
    eq1_midas_scale_gscpi = 0.45,
    eq1_midas_scale_ippi  = 0.15,
    eq1_midas_scale_epi   = 0.08,
    
    # theta parameters for each predictor
    eq1_theta1_gscpi = -0.25, eq1_theta2_gscpi = 0.035,
    eq1_theta1_ippi  = -0.12, eq1_theta2_ippi  = 0.018,
    eq1_theta1_epi   = -0.18, eq1_theta2_epi   = 0.025
  )
  
  # Equation 2: IR
  par_eq2 <- c(
    eq2_alpha = 0.2,                                        # intercept
    eq2_phi_contemp = 0.35,                               # contemporaneous CPI effect
    if (p > 0) {
      c(eq2_phi1 = 0.05, eq2_phi2 = 0.55, eq2_phi3 = 0.04)
    } else {
      numeric(0)    # Cross-lags: [cpi_lag1, ir_lag1, nhpi_lag1] for p=1
    },

    # MIDAS scales
    eq2_midas_scale_gscpi = 0.12,
    eq2_midas_scale_ippi  = -0.38,
    eq2_midas_scale_epi   = 0.08,
    
    # theta parameters for each predictor
    eq2_theta1_gscpi = 0.15, eq2_theta2_gscpi = -0.022,
    eq2_theta1_ippi  = 0.22, eq2_theta2_ippi  = -0.032,
    eq2_theta1_epi   = 0.10, eq2_theta2_epi   = -0.015
  )
  
  #Equation 3: NHPI
  par_eq3 <- c(
    eq3_alpha = 0.25,                                      #intercept
    eq3_phi_contemp_cpi = 0.28, eq3_phi_contemp_ir = -0.22,     # contemporaneous CPI and IR effect
    if (p > 0) {
      c(eq3_phi1 = 0.06, eq3_phi2 = -0.07, eq3_phi3 = 0.48) 
    } else {
      numeric(0)    # Cross-lags: [cpi_lag1, ir_lag1, nhpi_lag1] for p=1
    },
    
    # MIDAS scales
    eq3_midas_scale_gscpi = 0.10,
    eq3_midas_scale_ippi  = 0.18,
    eq3_midas_scale_epi   = -0.32,
    
    # theta parameters for each predictor
    eq3_theta1_gscpi = 0.20, eq3_theta2_gscpi = -0.028,
    eq3_theta1_ippi  = 0.14, eq3_theta2_ippi  = -0.020,
    eq3_theta1_epi   = 0.25, eq3_theta2_epi   = -0.038
  )
  
  true_par <- c(par_eq1, par_eq2, par_eq3)
  stopifnot(length(true_par) == (n_params_eq1 + n_params_eq2 + n_params_eq3))
  true_par
}

# Simulate data using the true parameter vector (to verify DGP)
simulate_data_from_par <- function(true_par, n, m, K, p, dates_lf, dates_hf, gscpi_hf, ippi_hf, epi_hf) {
  # Unpack parameter vector
  n_params_eq1 <- 1 + 3*p + 9
  n_params_eq2 <- 2 + 3*p + 9
  n_params_eq3 <- 3 + 3*p + 9
  
  par_eq1 <- true_par[1:n_params_eq1]
  par_eq2 <- true_par[(n_params_eq1 + 1):(n_params_eq1 + n_params_eq2)]
  par_eq3 <- true_par[(n_params_eq1 + n_params_eq2 + 1):length(true_par)]
  
  # Extract parameters for equation 1 (CPI)
  alpha1  <- par_eq1[1]
  phi1 <- if (p > 0) par_eq1[2:(1 + 3*p)] else numeric(0)
  midas_scales_1 <- par_eq1[(1 + 3*p + 1):(1 + 3*p + 3)]
  theta_gscpi1 <- par_eq1[(1 + 3*p + 3 + 1):(1 + 3*p + 3 + 2)]
  theta_ippi1  <- par_eq1[(1 + 3*p + 3 + 3):(1 + 3*p + 3 + 4)]
  theta_epi1   <- par_eq1[(1 + 3*p + 3 + 5):(1 + 3*p + 3 + 6)]
  
  # Extract parameters for equation 2 (IR)
  alpha2  <- par_eq2[1] 
  phi_contemp2 <- par_eq2[2]
  phi2 <- if (p > 0) par_eq2[3:(2 + 3*p)] else numeric(0)
  midas_scales_2 <- par_eq2[(2 + 3*p + 1):(2 + 3*p + 3)]
  theta_gscpi2 <- par_eq2[(2 + 3*p + 3 + 1):(2 + 3*p + 3 + 2)]
  theta_ippi2  <- par_eq2[(2 + 3*p + 3 + 3):(2 + 3*p + 3 + 4)]
  theta_epi2   <- par_eq2[(2 + 3*p + 3 + 5):(2 + 3*p + 3 + 6)]
  
  # Extract parameters for equation 3 (NHPI)
  alpha3  <- par_eq3[1] 
  phi_contemp_cpi3 <- par_eq3[2] 
  phi_contemp_ir3  <- par_eq3[3]
  phi3 <- if (p > 0) par_eq3[4:(3 + 3*p)] else numeric(0)
  midas_scales_3 <- par_eq3[(3 + 3*p + 1):(3 + 3*p + 3)]
  theta_gscpi3 <- par_eq3[(3 + 3*p + 3 + 1):(3 + 3*p + 3 + 2)]
  theta_ippi3  <- par_eq3[(3 + 3*p + 3 + 3):(3 + 3*p + 3 + 4)]
  theta_epi3   <- par_eq3[(3 + 3*p + 3 + 5):(3 + 3*p + 3 + 6)]
  
  hf_block <- function(hf_series, n, m, K) {
    t(vapply(1:n, function(t) {
      hf_index_end <- (t - 1) * m 
      hf_index_start <- hf_index_end - K
      if (hf_index_start < 1) return(rep(NA, K + 1))
      hf_series[hf_index_start:hf_index_end]
    }, numeric(K + 1)))
  }
  standardize_cols <- function(X_matrix) {
    col_mean <- colMeans(X_matrix, na.rm = TRUE)
    col_sd   <- apply(X_matrix, 2, sd, na.rm = TRUE)
    col_sd[col_sd < 1e-8] <- 1
    scaled <- scale(X_matrix, center = col_mean, scale = col_sd)
    scaled
  }
  
  gscpi_std <- standardize_cols(hf_block(gscpi_hf, n, m, K))
  ippi_std  <- standardize_cols(hf_block(ippi_hf,  n, m, K))
  epi_std   <- standardize_cols(hf_block(epi_hf,   n, m, K))
  
  j <- 0:K
  w_gscpi1 <- almon_weight_function(j, theta_gscpi1[1], theta_gscpi1[2], K)
  w_ippi1  <- almon_weight_function(j, theta_ippi1[1],  theta_ippi1[2],  K)
  w_epi1   <- almon_weight_function(j, theta_epi1[1],   theta_epi1[2],   K)
  w_gscpi2 <- almon_weight_function(j, theta_gscpi2[1], theta_gscpi2[2], K)
  w_ippi2  <- almon_weight_function(j, theta_ippi2[1],  theta_ippi2[2],  K)
  w_epi2   <- almon_weight_function(j, theta_epi2[1],   theta_epi2[2],   K)
  w_gscpi3 <- almon_weight_function(j, theta_gscpi3[1], theta_gscpi3[2], K)
  w_ippi3  <- almon_weight_function(j, theta_ippi3[1],  theta_ippi3[2],  K)
  w_epi3   <- almon_weight_function(j, theta_epi3[1],   theta_epi3[2],   K)
  
  cpi_lf  <- rep(NA, n)
  ir_lf   <- rep(NA, n)
  nhpi_lf <- rep(NA, n)
  cpi_lf[1] <- 0.05; ir_lf[1] <- 0.03; nhpi_lf[1] <- 0.04
  
  for (t in 2:n) {
    if (any(is.na(gscpi_std[t, ])) || any(is.na(ippi_std[t, ])) || any(is.na(epi_std[t, ]))) {
      gscpi_term1 <- ippi_term1 <- epi_term1 <- 0
      gscpi_term2 <- ippi_term2 <- epi_term2 <- 0
      gscpi_term3 <- ippi_term3 <- epi_term3 <- 0
    } else {
      gscpi_row <- as.numeric(gscpi_std[t, ])
      ippi_row  <- as.numeric(ippi_std[t, ])
      epi_row   <- as.numeric(epi_std[t, ])
      gscpi_term1 <- sum(w_gscpi1 * gscpi_row);  ippi_term1 <- sum(w_ippi1 * ippi_row)
      epi_term1 <- sum(w_epi1 * epi_row)
      gscpi_term2 <- sum(w_gscpi2 * gscpi_row);  ippi_term2 <- sum(w_ippi2 * ippi_row)  
      epi_term2 <- sum(w_epi2 * epi_row)
      gscpi_term3 <- sum(w_gscpi3 * gscpi_row);  ippi_term3 <- sum(w_ippi3 * ippi_row)  
      epi_term3 <- sum(w_epi3 * epi_row)
    }
    
    # Equation 1: CPI
    ar_terms1 <- if (p > 0) (phi1[1]*cpi_lf[t-1] + phi1[2]*ir_lf[t-1] + phi1[3]*nhpi_lf[t-1]) else 0
    cpi_lf[t] <- alpha1 + ar_terms1 +
      midas_scales_1[1]*gscpi_term1 + midas_scales_1[2]*ippi_term1 + midas_scales_1[3]*epi_term1 +
      rnorm(1, sd = 0.15)
    
    # Equation 2: IR
    ar_terms2 <- if (p > 0) (phi2[1]*cpi_lf[t-1] + phi2[2]*ir_lf[t-1] + phi2[3]*nhpi_lf[t-1]) else 0
    ir_lf[t] <- alpha2 + phi_contemp2*cpi_lf[t] + ar_terms2 +
      midas_scales_2[1]*gscpi_term2 + midas_scales_2[2]*ippi_term2 + midas_scales_2[3]*epi_term2 +
      rnorm(1, sd = 0.12)
    
    # Equation 3: NHPI
    ar_terms3 <- if (p > 0) (phi3[1]*cpi_lf[t-1] + phi3[2]*ir_lf[t-1] + phi3[3]*nhpi_lf[t-1]) else 0
    nhpi_lf[t] <- alpha3 + phi_contemp_cpi3*cpi_lf[t] + phi_contemp_ir3*ir_lf[t] + ar_terms3 +
      midas_scales_3[1]*gscpi_term3 + midas_scales_3[2]*ippi_term3 + midas_scales_3[3]*epi_term3 +
      rnorm(1, sd = 0.10)
  }
  list(cpi = cpi_lf, ir = ir_lf, nhpi = nhpi_lf)
}


create_midas_blocks <- function(hf_series, n, m, K, prefix) {
  blocks <- t(vapply(1:n, function(t) {
    hf_index_end <- (t - 1) * m 
    hf_index_start <- hf_index_end - K
    if (hf_index_start < 1) return(rep(NA, K + 1))
    hf_series[hf_index_start:hf_index_end]
  }, numeric(K + 1)))
  colnames(blocks) <- paste0(prefix, "_lag", 0:K)
  blocks
}


standardize_block <- function(X_matrix) {
  col_mean  <- colMeans(X_matrix, na.rm = TRUE)
  col_sd <- apply(X_matrix, 2, sd, na.rm = TRUE)
  col_sd[col_sd < 1e-8] <- 1
  hf_block_matrix <- scale(X_matrix, center = col_mean, scale = col_sd)
  colnames(hf_block_matrix) <- colnames(X_matrix)
  hf_block_matrix
}

create_lag_matrix <- function(x, p, var_name) {
  if (p == 0) return(NULL)
  lag_matrix <- sapply(1:p, function(lag_order) c(rep(NA, lag_order), head(x, -lag_order)))
  colnames(lag_matrix) <- paste0(var_name, "_lag", 1:p)
  lag_matrix
}


simulate <- function(m = 3, K = 4, p = 1, taus = c(0.25, 0.50, 0.75, 0.90)) {
# Simulate data for three response variables (CPI, IR, NHPI) and three predictors (GSCPI, IPPI, epi)
  start_date <- as.Date("1998-01-01")
  end_date   <- as.Date("2025-02-01")

  dates_lf <- seq(start_date, end_date, by = "quarter")
  n        <- length(dates_lf)       # number of low-frequency observations (quarters)

  dates_hf <- seq(start_date, end_date, by = "month")
  n_months       <- length(dates_hf)    # total monthly observations

  # Simulate high-frequency predictors
  gscpi_hf <- arima.sim(list(ar = c(0.8, 0.1)), n = n_months, sd = 0.15)
  ippi_hf  <- arima.sim(list(ar = c(0.7, 0.15)), n = n_months, sd = 0.12)
  epi_hf   <- arima.sim(list(ar = c(0.6, 0.1)), n = n_months, sd = 0.10)
  epi_hf_ts   <- xts(epi_hf,   order.by = dates_hf)

  gscpi_blocks <- create_midas_blocks(gscpi_hf, n, m, K, "gscpi")
  ippi_blocks  <- create_midas_blocks(ippi_hf,  n, m, K, "ippi")
  epi_blocks   <- create_midas_blocks(epi_hf,   n, m, K, "epi")

  gscpi_blocks <- standardize_block(gscpi_blocks)
  ippi_blocks  <- standardize_block(ippi_blocks)
  epi_blocks   <- standardize_block(epi_blocks)

  # Create true parameter vector
  true_par <- create_true_par_vector(p)

  # Simulate data using the true parameter vector
  simulated_data <- simulate_data_from_par(true_par, n, m, K, p, dates_lf, dates_hf, gscpi_blocks, ippi_blocks, epi_blocks)

  # Create time series objects with the exact dates
  cpi_ts  <- xts(simulated_data$cpi,  order.by = dates_lf)
  ir_ts   <- xts(simulated_data$ir,   order.by = dates_lf)
  nhpi_ts <- xts(simulated_data$nhpi, order.by = dates_lf)

  gscpi_hf_ts <- xts(gscpi_hf, order.by = dates_hf)
  ippi_hf_ts  <- xts(ippi_hf,  order.by = dates_hf)

  cpi_lags  <- create_lag_matrix(cpi_ts,  p, "cpi")
  ir_lags   <- create_lag_matrix(ir_ts,   p, "ir")
  nhpi_lags <- create_lag_matrix(nhpi_ts, p, "nhpi")

  model_data_res_var <- data.frame(
    cpi  = as.numeric(cpi_ts),
    ir   = as.numeric(ir_ts),
    nhpi = as.numeric(nhpi_ts)
  )
  model_data_res <- if (p > 0) cbind(model_data_res_var, cpi_lags, ir_lags, nhpi_lags) else model_data_res_var
  model_data <- cbind(model_data_res, gscpi_blocks, ippi_blocks, epi_blocks)
  model_data <- na.omit(model_data)
  model_data
}

## LOSS and FIT FUNCTIONS
# Huberized pinball loss 
qhuber_loss_vec <- function(u, tau, kappa) {
  # Define transition thresholds
  lower_bound <- -(1 - tau) * kappa
  upper_bound <-  tau * kappa
  loss <- numeric(length(u))                               # Initialize vector for losses
  in_smooth_zone <- (u > lower_bound) & (u < upper_bound)  # Inside quadratic smoothing zone
  loss[in_smooth_zone] <- 0.5 * (u[in_smooth_zone]^2) / kappa # Quadratic region (Huberized center)
  # Right linear region (u >= upper_bound): slope = tau
  loss[u >= upper_bound] <- tau * u[u >= upper_bound] - 0.5 * (tau^2) * kappa
  # Left linear region (u <= lower_bound): slope = 1 - tau
  loss[u <= lower_bound] <- (1 - tau) * (-u[u <= lower_bound]) - 0.5 * ((1 - tau)^2) * kappa
  loss
}

# Derivative of Huberized pinball loss
dqhuber_loss_du <- function(u, tau, kappa) {
  lower_bound <- -(1 - tau) * kappa
  upper_bound <- tau * kappa
  dL_du <- numeric(length(u))
  # Smooth zone: u in (lower_bound, upper_bound)
  in_smooth_zone <- (u > lower_bound) & (u < upper_bound)
  dL_du[in_smooth_zone] <- u[in_smooth_zone] / kappa
  dL_du[u >= upper_bound] <- tau    # Upper tail: u >= upper_bound
  dL_du[u <= lower_bound] <- -(1 - tau)    # Lower tail: u <= lower_bound 
  return(dL_du)
}

joint_loss_matrix_smooth <- function(par, data, tau, K, p, kappa = 0.01, theta_shrink = 2e-3) {
  n_params_eq1 <- 1 + 3*p + 9
  n_params_eq2 <- 2 + 3*p + 9
  n_params_eq3 <- 3 + 3*p + 9
  n  <- nrow(data)
  
  par_eq1 <- par[1:n_params_eq1]
  par_eq2 <- par[(n_params_eq1 + 1):(n_params_eq1 + n_params_eq2)]
  par_eq3 <- par[(n_params_eq1 + n_params_eq2 + 1):length(par)]
  
  j  <- 0:K
  X_gscpi <- as.matrix(data[, paste0("gscpi_lag", 0:K)])
  X_ippi  <- as.matrix(data[, paste0("ippi_lag",  0:K)])
  X_epi   <- as.matrix(data[, paste0("epi_lag",   0:K)])
  
  ## Eq 1: CPI
  alpha1 <- par_eq1[1]
  phi1   <- if (p > 0) par_eq1[2:(1 + 3*p)] else numeric(0)
  phi_end1 <- 1 + 3*p
  midas_scales_1 <- par_eq1[(phi_end1 + 1):(phi_end1 + 3)]
  theta_gscpi1   <- par_eq1[(phi_end1 + 4):(phi_end1 + 5)]
  theta_ippi1    <- par_eq1[(phi_end1 + 6):(phi_end1 + 7)]
  theta_epi1     <- par_eq1[(phi_end1 + 8):(phi_end1 + 9)]
  w_gscpi1 <- almon_weight_function(j, theta_gscpi1[1], theta_gscpi1[2], K)
  w_ippi1  <- almon_weight_function(j, theta_ippi1[1],  theta_ippi1[2],  K)
  w_epi1   <- almon_weight_function(j, theta_epi1[1],   theta_epi1[2],   K)
  Z_cpi <- if (p > 0) cbind(1,
                            as.matrix(data[, paste0("cpi_lag",  1:p)]),
                            as.matrix(data[, paste0("ir_lag",   1:p)]),
                            as.matrix(data[, paste0("nhpi_lag", 1:p)])) else cbind(rep(1, n))
  beta_cpi <- c(alpha1, phi1)
  midas_cpi <- cbind(X_gscpi %*% w_gscpi1, X_ippi %*% w_ippi1, X_epi %*% w_epi1)
  y_hat1 <- as.vector(Z_cpi %*% beta_cpi + midas_cpi %*% matrix(midas_scales_1, ncol = 1))
  u1 <- data$cpi - y_hat1
  
  ## Eq 2: IR
  alpha2 <- par_eq2[1]; phi_contemp2 <- par_eq2[2]
  phi2   <- if (p > 0) par_eq2[3:(2 + 3*p)] else numeric(0)
  phi_end2 <- 2 + 3*p
  midas_scales_2 <- par_eq2[(phi_end2 + 1):(phi_end2 + 3)]
  theta_gscpi2   <- par_eq2[(phi_end2 + 4):(phi_end2 + 5)]
  theta_ippi2    <- par_eq2[(phi_end2 + 6):(phi_end2 + 7)]
  theta_epi2     <- par_eq2[(phi_end2 + 8):(phi_end2 + 9)]
  w_gscpi2 <- almon_weight_function(j, theta_gscpi2[1], theta_gscpi2[2], K)
  w_ippi2  <- almon_weight_function(j, theta_ippi2[1],  theta_ippi2[2],  K)
  w_epi2   <- almon_weight_function(j, theta_epi2[1],   theta_epi2[2],   K)
  Z_ir <- if (p > 0) cbind(1, data$cpi,
                           as.matrix(data[, paste0("cpi_lag",  1:p)]),
                           as.matrix(data[, paste0("ir_lag",   1:p)]),
                           as.matrix(data[, paste0("nhpi_lag", 1:p)])) else cbind(rep(1, n), data$cpi)
  beta_ir <- c(alpha2, phi_contemp2, phi2)
  midas_ir <- cbind(X_gscpi %*% w_gscpi2, X_ippi %*% w_ippi2, X_epi %*% w_epi2)
  y_hat2 <- as.vector(Z_ir %*% beta_ir + midas_ir %*% matrix(midas_scales_2, ncol = 1))
  u2 <- data$ir - y_hat2
  
  ## Eq 3: NHPI
  alpha3 <- par_eq3[1]; phi_contemp_cpi3 <- par_eq3[2]; phi_contemp_ir3 <- par_eq3[3]
  phi3   <- if (p > 0) par_eq3[4:(3 + 3*p)] else numeric(0)
  phi_end3 <- 3 + 3*p
  midas_scales_3 <- par_eq3[(phi_end3 + 1):(phi_end3 + 3)]
  theta_gscpi3   <- par_eq3[(phi_end3 + 4):(phi_end3 + 5)]
  theta_ippi3    <- par_eq3[(phi_end3 + 6):(phi_end3 + 7)]
  theta_epi3     <- par_eq3[(phi_end3 + 8):(phi_end3 + 9)]
  w_gscpi3 <- almon_weight_function(j, theta_gscpi3[1], theta_gscpi3[2], K)
  w_ippi3  <- almon_weight_function(j, theta_ippi3[1],  theta_ippi3[2],  K)
  w_epi3   <- almon_weight_function(j, theta_epi3[1],   theta_epi3[2],   K)
  Z_nhpi <- if (p > 0) cbind(1, data$cpi, data$ir,
                             as.matrix(data[, paste0("cpi_lag",  1:p)]),
                             as.matrix(data[, paste0("ir_lag",   1:p)]),
                             as.matrix(data[, paste0("nhpi_lag", 1:p)])) else cbind(rep(1, n), data$cpi, data$ir)
  beta_nhpi <- c(alpha3, phi_contemp_cpi3, phi_contemp_ir3, phi3)
  midas_nhpi <- cbind(X_gscpi %*% w_gscpi3, X_ippi %*% w_ippi3, X_epi %*% w_epi3)
  y_hat3 <- as.vector(Z_nhpi %*% beta_nhpi + midas_nhpi %*% matrix(midas_scales_3, ncol = 1))
  u3 <- data$nhpi - y_hat3
  
  loss <- sum(qhuber_loss_vec(u1, tau, kappa)) +
    sum(qhuber_loss_vec(u2, tau, kappa)) +
    sum(qhuber_loss_vec(u3, tau, kappa))
  
  if (theta_shrink > 0) {
    theta_param_indices <- c(
      (1 + 3*p + 4):(1 + 3*p + 9),
      (n_params_eq1 + 2 + 3*p + 4):(n_params_eq1 + 2 + 3*p + 9),
      (n_params_eq1 + n_params_eq2 + 3 + 3*p + 4):(n_params_eq1 + n_params_eq2 + 3 + 3*p + 9)
    )
    loss <- loss + theta_shrink * sum(par[theta_param_indices]^2)
  }
  loss
}


## ANALYTIC GRADIENT FUNCTION
joint_gradient_matrix_smooth <- function(par, data, tau, K, p, kappa = 0.01, theta_shrink = 2e-3) {
  n_params_eq1 <- 1 + 3*p + 9
  n_params_eq2 <- 2 + 3*p + 9
  n_params_eq3 <- 3 + 3*p + 9
  n <- nrow(data)
  
  # Initialize gradient vector
  grad <- numeric(length(par))
  
  # Unpack parameters
  par_eq1 <- par[1:n_params_eq1]
  par_eq2 <- par[(n_params_eq1 + 1):(n_params_eq1 + n_params_eq2)]
  par_eq3 <- par[(n_params_eq1 + n_params_eq2 + 1):length(par)]
  
  j <- 0:K
  X_gscpi <- as.matrix(data[, paste0("gscpi_lag", 0:K)])
  X_ippi  <- as.matrix(data[, paste0("ippi_lag",  0:K)])
  X_epi   <- as.matrix(data[, paste0("epi_lag",   0:K)])
  
  ## Equation 1: CPI
  alpha1 <- par_eq1[1]
  phi1   <- if (p > 0) par_eq1[2:(1 + 3*p)] else numeric(0)
  phi_end1 <- 1 + 3*p
  midas_scales_1 <- par_eq1[(phi_end1 + 1):(phi_end1 + 3)]
  theta_gscpi1   <- par_eq1[(phi_end1 + 4):(phi_end1 + 5)]
  theta_ippi1    <- par_eq1[(phi_end1 + 6):(phi_end1 + 7)]
  theta_epi1     <- par_eq1[(phi_end1 + 8):(phi_end1 + 9)]
  
  # Compute weights and their derivatives
  w_gscpi1 <- almon_weight_function(j, theta_gscpi1[1], theta_gscpi1[2], K)
  w_ippi1  <- almon_weight_function(j, theta_ippi1[1],  theta_ippi1[2],  K)
  w_epi1   <- almon_weight_function(j, theta_epi1[1],   theta_epi1[2],   K)
  
  # Compute weight gradients
  dw_gscpi1_dtheta <- almon_weight_gradient(j, theta_gscpi1[1], theta_gscpi1[2], K)
  dw_ippi1_dtheta  <- almon_weight_gradient(j, theta_ippi1[1],  theta_ippi1[2],  K)
  dw_epi1_dtheta   <- almon_weight_gradient(j, theta_epi1[1],   theta_epi1[2],   K)
  
  # Design matrix and predictions
  Z_cpi <- if (p > 0) cbind(1,
                            as.matrix(data[, paste0("cpi_lag",  1:p)]),
                            as.matrix(data[, paste0("ir_lag",   1:p)]),
                            as.matrix(data[, paste0("nhpi_lag", 1:p)])) else cbind(rep(1, n))
  beta_cpi <- c(alpha1, phi1)
  
  midas_cpi <- cbind(X_gscpi %*% w_gscpi1, X_ippi %*% w_ippi1, X_epi %*% w_epi1)
  y_hat1 <- as.vector(Z_cpi %*% beta_cpi + midas_cpi %*% matrix(midas_scales_1, ncol = 1))
  u1 <- data$cpi - y_hat1
  
  ## Equation 2: IR  
  alpha2 <- par_eq2[1]; phi_contemp2 <- par_eq2[2]
  phi2   <- if (p > 0) par_eq2[3:(2 + 3*p)] else numeric(0)
  phi_end2 <- 2 + 3*p
  midas_scales_2 <- par_eq2[(phi_end2 + 1):(phi_end2 + 3)]
  theta_gscpi2   <- par_eq2[(phi_end2 + 4):(phi_end2 + 5)]
  theta_ippi2    <- par_eq2[(phi_end2 + 6):(phi_end2 + 7)]
  theta_epi2     <- par_eq2[(phi_end2 + 8):(phi_end2 + 9)]
  
  w_gscpi2 <- almon_weight_function(j, theta_gscpi2[1], theta_gscpi2[2], K)
  w_ippi2  <- almon_weight_function(j, theta_ippi2[1],  theta_ippi2[2],  K)
  w_epi2   <- almon_weight_function(j, theta_epi2[1],   theta_epi2[2],   K)
  
  dw_gscpi2_dtheta <- almon_weight_gradient(j, theta_gscpi2[1], theta_gscpi2[2], K)
  dw_ippi2_dtheta  <- almon_weight_gradient(j, theta_ippi2[1],  theta_ippi2[2],  K)
  dw_epi2_dtheta   <- almon_weight_gradient(j, theta_epi2[1],   theta_epi2[2],   K)
  
  Z_ir <- if (p > 0) cbind(1, data$cpi,
                           as.matrix(data[, paste0("cpi_lag",  1:p)]),
                           as.matrix(data[, paste0("ir_lag",   1:p)]),
                           as.matrix(data[, paste0("nhpi_lag", 1:p)])) else cbind(rep(1, n), data$cpi)
  beta_ir <- c(alpha2, phi_contemp2, phi2)
  
  midas_ir <- cbind(X_gscpi %*% w_gscpi2, X_ippi %*% w_ippi2, X_epi %*% w_epi2)
  y_hat2 <- as.vector(Z_ir %*% beta_ir + midas_ir %*% matrix(midas_scales_2, ncol = 1))
  u2 <- data$ir - y_hat2
  
  ## Equation 3: NHPI
  alpha3 <- par_eq3[1]; phi_contemp_cpi3 <- par_eq3[2]; phi_contemp_ir3 <- par_eq3[3]
  phi3   <- if (p > 0) par_eq3[4:(3 + 3*p)] else numeric(0)
  phi_end3 <- 3 + 3*p
  midas_scales_3 <- par_eq3[(phi_end3 + 1):(phi_end3 + 3)]
  theta_gscpi3   <- par_eq3[(phi_end3 + 4):(phi_end3 + 5)]
  theta_ippi3    <- par_eq3[(phi_end3 + 6):(phi_end3 + 7)]
  theta_epi3     <- par_eq3[(phi_end3 + 8):(phi_end3 + 9)]
  
  w_gscpi3 <- almon_weight_function(j, theta_gscpi3[1], theta_gscpi3[2], K)
  w_ippi3  <- almon_weight_function(j, theta_ippi3[1],  theta_ippi3[2],  K)
  w_epi3   <- almon_weight_function(j, theta_epi3[1],   theta_epi3[2],   K)
  
  dw_gscpi3_dtheta <- almon_weight_gradient(j, theta_gscpi3[1], theta_gscpi3[2], K)
  dw_ippi3_dtheta  <- almon_weight_gradient(j, theta_ippi3[1],  theta_ippi3[2],  K)
  dw_epi3_dtheta   <- almon_weight_gradient(j, theta_epi3[1],   theta_epi3[2],   K)
  
  Z_nhpi <- if (p > 0) cbind(1, data$cpi, data$ir,
                             as.matrix(data[, paste0("cpi_lag",  1:p)]),
                             as.matrix(data[, paste0("ir_lag",   1:p)]),
                             as.matrix(data[, paste0("nhpi_lag", 1:p)])) else cbind(rep(1, n), data$cpi, data$ir)
  beta_nhpi <- c(alpha3, phi_contemp_cpi3, phi_contemp_ir3, phi3)
  
  midas_nhpi <- cbind(X_gscpi %*% w_gscpi3, X_ippi %*% w_ippi3, X_epi %*% w_epi3)
  y_hat3 <- as.vector(Z_nhpi %*% beta_nhpi + midas_nhpi %*% matrix(midas_scales_3, ncol = 1))
  u3 <- data$nhpi - y_hat3
  
  ## Compute Huberized pinball loss derivatives
  dL_du1 <- dqhuber_loss_du(u1, tau, kappa)
  dL_du2 <- dqhuber_loss_du(u2, tau, kappa) 
  dL_du3 <- dqhuber_loss_du(u3, tau, kappa)
  
  ## Equation 1 gradients
  # Gradient for alpha and phi parameters (Eq1)
  grad[1:(1 + 3*p)] <- -t(Z_cpi) %*% dL_du1
  
  # Gradient for MIDAS scales (Eq1)
  grad[(phi_end1 + 1):(phi_end1 + 3)] <- -t(midas_cpi) %*% dL_du1
  
  # Gradient for theta parameters (Eq1)
  # theta_gscpi1
  dmidas_dtheta1_gscpi1 <- X_gscpi %*% dw_gscpi1_dtheta$grad_theta1
  dmidas_dtheta2_gscpi1 <- X_gscpi %*% dw_gscpi1_dtheta$grad_theta2
  grad[phi_end1 + 4] <- -sum(dL_du1 * (midas_scales_1[1] * dmidas_dtheta1_gscpi1))
  grad[phi_end1 + 5] <- -sum(dL_du1 * (midas_scales_1[1] * dmidas_dtheta2_gscpi1))
  
  # theta_ippi1
  dmidas_dtheta1_ippi1 <- X_ippi %*% dw_ippi1_dtheta$grad_theta1
  dmidas_dtheta2_ippi1 <- X_ippi %*% dw_ippi1_dtheta$grad_theta2
  grad[phi_end1 + 6] <- -sum(dL_du1 * (midas_scales_1[2] * dmidas_dtheta1_ippi1))
  grad[phi_end1 + 7] <- -sum(dL_du1 * (midas_scales_1[2] * dmidas_dtheta2_ippi1))
  
  # theta_epi1
  dmidas_dtheta1_epi1 <- X_epi %*% dw_epi1_dtheta$grad_theta1
  dmidas_dtheta2_epi1 <- X_epi %*% dw_epi1_dtheta$grad_theta2
  grad[phi_end1 + 8] <- -sum(dL_du1 * (midas_scales_1[3] * dmidas_dtheta1_epi1))
  grad[phi_end1 + 9] <- -sum(dL_du1 * (midas_scales_1[3] * dmidas_dtheta2_epi1))
  
  ## Equation 2 gradients
  idx_start_eq2 <- n_params_eq1 + 1
  
  # Gradient for alpha, phi_contemp, phi (Eq2)
  grad[idx_start_eq2:(idx_start_eq2 + 1 + 3*p)] <- -t(Z_ir) %*% dL_du2
  
  # Gradient for MIDAS scales (Eq2)
  grad[idx_start_eq2 + 2 + 3*p + (0:2)] <- -t(midas_ir) %*% dL_du2
  
  # Gradient for theta parameters (Eq2)
  # theta_gscpi2
  dmidas_dtheta1_gscpi2 <- X_gscpi %*% dw_gscpi2_dtheta$grad_theta1
  dmidas_dtheta2_gscpi2 <- X_gscpi %*% dw_gscpi2_dtheta$grad_theta2
  grad[idx_start_eq2 + 2 + 3*p + 3] <- -sum(dL_du2 * (midas_scales_2[1] * dmidas_dtheta1_gscpi2))
  grad[idx_start_eq2 + 2 + 3*p + 4] <- -sum(dL_du2 * (midas_scales_2[1] * dmidas_dtheta2_gscpi2))
  
  # theta_ippi2
  dmidas_dtheta1_ippi2 <- X_ippi %*% dw_ippi2_dtheta$grad_theta1
  dmidas_dtheta2_ippi2 <- X_ippi %*% dw_ippi2_dtheta$grad_theta2
  grad[idx_start_eq2 + 2 + 3*p + 5] <- -sum(dL_du2 * (midas_scales_2[2] * dmidas_dtheta1_ippi2))
  grad[idx_start_eq2 + 2 + 3*p + 6] <- -sum(dL_du2 * (midas_scales_2[2] * dmidas_dtheta2_ippi2))
  
  # theta_epi2
  dmidas_dtheta1_epi2 <- X_epi %*% dw_epi2_dtheta$grad_theta1
  dmidas_dtheta2_epi2 <- X_epi %*% dw_epi2_dtheta$grad_theta2
  grad[idx_start_eq2 + 2 + 3*p + 7] <- -sum(dL_du2 * (midas_scales_2[3] * dmidas_dtheta1_epi2))
  grad[idx_start_eq2 + 2 + 3*p + 8] <- -sum(dL_du2 * (midas_scales_2[3] * dmidas_dtheta2_epi2))
  
  ## Equation 3 gradients
  idx_start_eq3 <- n_params_eq1 + n_params_eq2 + 1
  
  # Gradient for alpha, phi_contemp_cpi, phi_contemp_ir, phi (Eq3)
  grad[idx_start_eq3:(idx_start_eq3 + 2 + 3*p)] <- -t(Z_nhpi) %*% dL_du3
  
  # Gradient for MIDAS scales (Eq3)
  grad[idx_start_eq3 + 3 + 3*p + (0:2)] <- -t(midas_nhpi) %*% dL_du3
  
  # Gradient for theta parameters (Eq3)
  # theta_gscpi3
  dmidas_dtheta1_gscpi3 <- X_gscpi %*% dw_gscpi3_dtheta$grad_theta1
  dmidas_dtheta2_gscpi3 <- X_gscpi %*% dw_gscpi3_dtheta$grad_theta2
  grad[idx_start_eq3 + 3 + 3*p + 3] <- -sum(dL_du3 * (midas_scales_3[1] * dmidas_dtheta1_gscpi3))
  grad[idx_start_eq3 + 3 + 3*p + 4] <- -sum(dL_du3 * (midas_scales_3[1] * dmidas_dtheta2_gscpi3))
  
  # theta_ippi3
  dmidas_dtheta1_ippi3 <- X_ippi %*% dw_ippi3_dtheta$grad_theta1
  dmidas_dtheta2_ippi3 <- X_ippi %*% dw_ippi3_dtheta$grad_theta2
  grad[idx_start_eq3 + 3 + 3*p + 5] <- -sum(dL_du3 * (midas_scales_3[2] * dmidas_dtheta1_ippi3))
  grad[idx_start_eq3 + 3 + 3*p + 6] <- -sum(dL_du3 * (midas_scales_3[2] * dmidas_dtheta2_ippi3))
  
  # theta_epi3
  dmidas_dtheta1_epi3 <- X_epi %*% dw_epi3_dtheta$grad_theta1
  dmidas_dtheta2_epi3 <- X_epi %*% dw_epi3_dtheta$grad_theta2
  grad[idx_start_eq3 + 3 + 3*p + 7] <- -sum(dL_du3 * (midas_scales_3[3] * dmidas_dtheta1_epi3))
  grad[idx_start_eq3 + 3 + 3*p + 8] <- -sum(dL_du3 * (midas_scales_3[3] * dmidas_dtheta2_epi3))
  
  ## Add regularization for theta parameters
  if (theta_shrink > 0) {
    theta_param_indices <- c(
      (1 + 3*p + 4):(1 + 3*p + 9),
      (n_params_eq1 + 2 + 3*p + 4):(n_params_eq1 + 2 + 3*p + 9),
      (n_params_eq1 + n_params_eq2 + 3 + 3*p + 4):(n_params_eq1 + n_params_eq2 + 3 + 3*p + 9)
    )
    grad[theta_param_indices] <- grad[theta_param_indices] + 2 * theta_shrink * par[theta_param_indices]
  }
  
  return(grad)
}

calculate_fitted_values <- function(par, data, p, K) {
  n_params_eq1 <- 1 + 3*p + 9
  n_params_eq2 <- 2 + 3*p + 9
  n_params_eq3 <- 3 + 3*p + 9
  
  par_eq1 <- par[1:n_params_eq1]
  par_eq2 <- par[(n_params_eq1 + 1):(n_params_eq1 + n_params_eq2)]
  par_eq3 <- par[(n_params_eq1 + n_params_eq2 + 1):length(par)]
  
  n <- nrow(data); j <- 0:K
  X_gscpi <- as.matrix(data[, paste0("gscpi_lag", 0:K)])
  X_ippi  <- as.matrix(data[, paste0("ippi_lag",  0:K)])
  X_epi   <- as.matrix(data[, paste0("epi_lag",   0:K)])
  
  # Eq 1
  alpha1 <- par_eq1[1]
  phi1   <- if (p > 0) par_eq1[2:(1 + 3*p)] else numeric(0)
  phi_end1 <- 1 + 3*p
  midas_scales_1 <- par_eq1[(phi_end1+1):(phi_end1+3)]
  theta_gscpi1   <- par_eq1[(phi_end1+4):(phi_end1+5)]
  theta_ippi1    <- par_eq1[(phi_end1+6):(phi_end1+7)]
  theta_epi1     <- par_eq1[(phi_end1+8):(phi_end1+9)]
  w_gscpi1 <- almon_weight_function(j, theta_gscpi1[1], theta_gscpi1[2], K)
  w_ippi1  <- almon_weight_function(j, theta_ippi1[1],  theta_ippi1[2],  K)
  w_epi1   <- almon_weight_function(j, theta_epi1[1],   theta_epi1[2],   K)
  Z_cpi <- if (p > 0) cbind(1,
                            as.matrix(data[, paste0("cpi_lag",  1:p)]),
                            as.matrix(data[, paste0("ir_lag",   1:p)]),
                            as.matrix(data[, paste0("nhpi_lag", 1:p)]))
  else cbind(rep(1, n))
  beta_cpi <- c(alpha1, phi1)
  midas_cpi <- cbind(X_gscpi %*% w_gscpi1, X_ippi %*% w_ippi1, X_epi %*% w_epi1)
  yhat_cpi <- as.vector(Z_cpi %*% beta_cpi + midas_cpi %*% matrix(midas_scales_1, ncol = 1))
  
  # Eq 2
  alpha2 <- par_eq2[1]; phi_contemp2 <- par_eq2[2]
  phi2   <- if (p > 0) par_eq2[3:(2 + 3*p)] else numeric(0)
  phi_end2 <- 2 + 3*p
  midas_scales_2 <- par_eq2[(phi_end2+1):(phi_end2+3)]
  theta_gscpi2   <- par_eq2[(phi_end2+4):(phi_end2+5)]
  theta_ippi2    <- par_eq2[(phi_end2+6):(phi_end2+7)]
  theta_epi2     <- par_eq2[(phi_end2+8):(phi_end2+9)]
  w_gscpi2 <- almon_weight_function(j, theta_gscpi2[1], theta_gscpi2[2], K)
  w_ippi2  <- almon_weight_function(j, theta_ippi2[1],  theta_ippi2[2],  K)
  w_epi2   <- almon_weight_function(j, theta_epi2[1],   theta_epi2[2],   K)
  Z_ir <- if (p > 0) cbind(1, data$cpi,
                           as.matrix(data[, paste0("cpi_lag",  1:p)]),
                           as.matrix(data[, paste0("ir_lag",   1:p)]),
                           as.matrix(data[, paste0("nhpi_lag", 1:p)]))
  else cbind(rep(1, n), data$cpi)
  beta_ir <- c(alpha2, phi_contemp2, phi2)
  midas_ir <- cbind(X_gscpi %*% w_gscpi2, X_ippi %*% w_ippi2, X_epi %*% w_epi2)
  yhat_ir <- as.vector(Z_ir %*% beta_ir + midas_ir %*% matrix(midas_scales_2, ncol = 1))
  
  # Eq 3
  alpha3 <- par_eq3[1]; phi_contemp_cpi3 <- par_eq3[2]; phi_contemp_ir3 <- par_eq3[3]
  phi3   <- if (p > 0) par_eq3[4:(3 + 3*p)] else numeric(0)
  phi_end3 <- 3 + 3*p
  midas_scales_3 <- par_eq3[(phi_end3+1):(phi_end3+3)]
  theta_gscpi3   <- par_eq3[(phi_end3+4):(phi_end3+5)]
  theta_ippi3    <- par_eq3[(phi_end3+6):(phi_end3+7)]
  theta_epi3     <- par_eq3[(phi_end3+8):(phi_end3+9)]
  w_gscpi3 <- almon_weight_function(j, theta_gscpi3[1], theta_gscpi3[2], K)
  w_ippi3  <- almon_weight_function(j, theta_ippi3[1],  theta_ippi3[2],  K)
  w_epi3   <- almon_weight_function(j, theta_epi3[1],   theta_epi3[2],   K)
  Z_nhpi <- if (p > 0) cbind(1, data$cpi, data$ir,
                             as.matrix(data[, paste0("cpi_lag",  1:p)]),
                             as.matrix(data[, paste0("ir_lag",   1:p)]),
                             as.matrix(data[, paste0("nhpi_lag", 1:p)]))
  else cbind(rep(1, n), data$cpi, data$ir)
  beta_nhpi <- c(alpha3, phi_contemp_cpi3, phi_contemp_ir3, phi3)
  midas_nhpi <- cbind(X_gscpi %*% w_gscpi3, X_ippi %*% w_ippi3, X_epi %*% w_epi3)
  yhat_nhpi <- as.vector(Z_nhpi %*% beta_nhpi + midas_nhpi %*% matrix(midas_scales_3, ncol = 1))
  
  list(cpi_fitted = yhat_cpi, ir_fitted = yhat_ir, nhpi_fitted = yhat_nhpi)
}

# Calculate the confidence interval using the Hessian from the optim function
calc_ci_from_optim <- function(opt, p, level = 0.95) {
  if (is.null(opt$hessian)) {
    stop("Hessian not available in optim result. Make sure optim(..., hessian = TRUE).")
  }
  H <- opt$hessian
  if (any(!is.finite(H))) stop("Non-finite entries in Hessian.")
  
  # Invert Hessian to get covariance matrix
  vcov_mat <- tryCatch(
    solve(H),
    error = function(e) stop("Hessian not invertible (nearly singular): ", e$message)
  )
  
  est <- opt$par
  se  <- sqrt(diag(vcov_mat))
  
  alpha <- 1 - level
  z <- qnorm(1 - alpha/2)
  lower <- est - z * se
  upper <- est + z * se
  
  ## Parameter names to match your layout
  n_params_eq1 <- 1 + 3*p + 9
  n_params_eq2 <- 2 + 3*p + 9
  n_params_eq3 <- 3 + 3*p + 9
  
  phi_names_3p <- if (p > 0) paste0("phi", 1:(3*p)) else character(0)
  param_names_eq1 <- c("alpha", phi_names_3p,
                       "midas_scale_gscpi","midas_scale_ippi","midas_scale_epi",
                       "theta1_gscpi","theta2_gscpi",
                       "theta1_ippi","theta2_ippi",
                       "theta1_epi","theta2_epi")
  param_names_eq2 <- c("alpha","phi_contemp", phi_names_3p,
                       "midas_scale_gscpi","midas_scale_ippi","midas_scale_epi",
                       "theta1_gscpi","theta2_gscpi",
                       "theta1_ippi","theta2_ippi",
                       "theta1_epi","theta2_epi")
  param_names_eq3 <- c("alpha","phi_contemp_cpi","phi_contemp_ir", phi_names_3p,
                       "midas_scale_gscpi","midas_scale_ippi","midas_scale_epi",
                       "theta1_gscpi","theta2_gscpi",
                       "theta1_ippi","theta2_ippi",
                       "theta1_epi","theta2_epi")
  
  full_names <- c(param_names_eq1, param_names_eq2, param_names_eq3)
  stopifnot(length(full_names) == length(est))
  
  data.frame(
    param    = full_names,
    estimate = est,
    se       = se,
    ci_lower = lower,
    ci_upper = upper,
    row.names = NULL
  )
}

## Estimator function with analytic gradients - FIXED
estimate_joint_qarmidas <- function(
    taus, data, K, p,
    kappa = 0.01,
    theta_shrink = 2e-3,
    start_par = NULL,
    optimizer = c("BFGS", "L-BFGS-B"),
    use_bounds = FALSE, bounds_range = 5,
    n_restarts = 2, jitter_sd = 0.01,
    verbose = TRUE,
    dry_run = FALSE,
    use_analytic_grad = TRUE,  # New option for analytic gradients
    ci_level = 0.95            # Confidence level added     
) {
  optimizer <- match.arg(optimizer)
  n_params_eq1 <- 1 + 3*p + 9
  n_params_eq2 <- 2 + 3*p + 9
  n_params_eq3 <- 3 + 3*p + 9
  par  <- n_params_eq1 + n_params_eq2 + n_params_eq3
  
  if (is.null(start_par)) start_par <- create_true_par_vector(p)
  stopifnot(length(start_par) == par)
  
  result <- list()
  current_par <- start_par
  
  for (tau in sort(taus)) {
    key <- paste0("tau_", formatC(tau, digits = 2, format = "f"))
    if (verbose) cat("Estimating tau =", tau, "...\n")
    t0 <- Sys.time()
    
    kappa_tau <- if (abs(tau - 1) < 1e-12) max(0.02, kappa) else kappa
    
    ci_tab <- NULL   # reset CI holder for this tau
    
    if (dry_run) {
      best <- list(par = current_par, value = NA_real_, convergence = NA_integer_, hessian = NULL)
    } else {
      starts <- c(
        list(current_par),
        if (n_restarts > 0) replicate(n_restarts, current_par + rnorm(par, sd = jitter_sd), simplify = FALSE)
      )
      best <- NULL
      for (st in starts) {
        ctrl <- list(maxit = 4000, reltol = 1e-8)
        
        # Use analytic gradients if requested
        if (use_analytic_grad) {
          gr_func <- function(par) {
            joint_gradient_matrix_smooth(
              par = par, 
              data = data, 
              tau = tau, 
              K = K, 
              p = p,
              kappa = kappa_tau, 
              theta_shrink = theta_shrink
            )
          }
        } else {
          gr_func <- NULL
        }
        
        # Define objective function
        fn_obj <- function(par) {
          joint_loss_matrix_smooth(
            par = par,
            data = data,
            tau = tau,
            K = K,
            p = p,
            kappa = kappa_tau,
            theta_shrink = theta_shrink
          )
        }
        
        if (use_bounds || optimizer == "L-BFGS-B") {
          lo <- rep(-bounds_range, par); up <- rep(bounds_range, par)
          opt <- stats::optim(
            par = st, 
            fn = fn_obj,
            gr = gr_func,
            method = "L-BFGS-B", 
            lower = lo, 
            upper = up, 
            control = ctrl,
            hessian = TRUE                  # hessian added
          )
        } else {
          opt <- stats::optim(
            par = st, 
            fn = fn_obj,
            gr = gr_func,
            method = "BFGS", 
            control = ctrl,
            hessian = TRUE                 # hessian added
          )
        }
        if (is.null(best) || opt$value < best$value) best <- opt
      }
      
      # Compute CI table for this tau if Hessian is available
      ci_tab <- NULL
      if (!is.null(best$hessian)) {
        ci_tab <- calc_ci_from_optim(best, p, level = ci_level)
      }
    }
    
    fitted_values <- calculate_fitted_values(best$par, data, p, K)
    t1 <- Sys.time()
    result[[key]] <- list(
      opt_result    = best,
      ci_table      = ci_tab,                      # CIs stored added
      fitted_cpi    = fitted_values$cpi_fitted,
      fitted_ir     = fitted_values$ir_fitted,
      fitted_nhpi   = fitted_values$nhpi_fitted,
      elapsed_minutes = as.numeric(difftime(t1, t0, units = "mins")),
      convergence   = best$convergence,
      loss          = best$value
    )
    current_par <- best$par
    if (verbose && !dry_run)
      cat(" Done τ =", tau, "conv =", best$convergence, "loss =", round(best$value, 6), "\n")
    
  }
  
  # Turn result into an S3 object and store p and taus
  attr(result, "p")    <- p
  attr(result, "taus") <- taus
  class(result) <- "mvqarmidas_fit"
  
  result
}

# A summary() method for objects of class "mvqarmidas_fit"
summary.mvqarmidas_fit <- function(object, ...) {
  p    <- attr(object, "p")
  taus <- attr(object, "taus")
  if (is.null(p) || is.null(taus)) {
    stop("Object is missing 'p' or 'taus' attributes.")
  }
  
  n_params_eq1 <- 1 + 3*p + 9
  n_params_eq2 <- 2 + 3*p + 9
  n_params_eq3 <- 3 + 3*p + 9
  
  tau_keys <- paste0("tau_", format(taus, nsmall = 2))
  out <- vector("list", length(tau_keys))
  names(out) <- tau_keys
  
  enhance_table <- function(tab) {
    est <- tab$estimate
    se  <- tab$se
    z   <- est / se
    z[!is.finite(z)] <- NA
    pval <- 2 * pnorm(-abs(z))
    
    mat <- cbind(
      Estimate    = est,
      "Std.Error"= se,
      "z value"   = z,
      "Pr(>|z|)"  = pval,
      "2.5 %"     = tab$ci_lower,
      "97.5 %"    = tab$ci_upper
    )
    rownames(mat) <- tab$param
    as.data.frame(mat, check.names = FALSE)
  }
  
  for (i in seq_along(tau_keys)) {
    key <- tau_keys[i]
    res_k <- object[[key]]
    ci_tab <- res_k$ci_table
    if (is.null(ci_tab)) {
      out[[i]] <- NULL
      next
    }
    
    idx1 <- 1:n_params_eq1
    idx2 <- (n_params_eq1 + 1):(n_params_eq1 + n_params_eq2)
    idx3 <- (n_params_eq1 + n_params_eq2 + 1):(n_params_eq1 + n_params_eq2 + n_params_eq3)
    
    tab1 <- enhance_table(ci_tab[idx1, ])
    tab2 <- enhance_table(ci_tab[idx2, ])
    tab3 <- enhance_table(ci_tab[idx3, ])
    
    out[[i]] <- list(
      tau         = taus[i],
      eq1         = tab1,
      eq2         = tab2,
      eq3         = tab3,
      loss        = res_k$loss,
      convergence = res_k$convergence
    )
  }
  
  class(out) <- "summary.mvqarmidas_fit"
  out
}

print.summary.mvqarmidas_fit <- function(x, digits = 4, ...) {
  for (i in seq_along(x)) {
    tau_summary <- x[[i]]
    if (is.null(tau_summary)) next
    
    cat("Quantile tau =",tau_summary$tau, "\n")
    cat("Loss:", formatC(tau_summary$loss, digits = digits),
        "Convergence code:",tau_summary$convergence, "\n")
    
    cat("\nEquation 1 (CPI):\n")
    print(round(tau_summary$eq1, digits), ...)
    
    cat("\nEquation 2 (IR):\n")
    print(round(tau_summary$eq2, digits), ...)
    
    cat("\nEquation 3 (NHPI):\n")
    print(round(tau_summary$eq3, digits), ...)
    cat("\n")
  }
  invisible(x)
}


## 5) GRADIENT VERIFICATION
verify_gradients <- function(par_test, data, tau, K, p, kappa = 0.01, theta_shrink = 2e-3) {
  # Analytic gradient
  analytic_grad <- joint_gradient_matrix_smooth(par_test, data, tau, K, p, kappa, theta_shrink)
  
  # Numerical gradient
  numeric_grad <- grad(
    function(par) joint_loss_matrix_smooth(par, data, tau, K, p, kappa, theta_shrink),
    par_test
  )
  
  # Compare
  diff <- analytic_grad - numeric_grad
  max_diff <- max(abs(diff))
  rmse_diff <- sqrt(mean(diff^2))
  rel_rmse <- rmse_diff / (sqrt(mean(numeric_grad^2)) + 1e-10)
  
  cat("=== GRADIENT VERIFICATION ===\n")
  cat("Max difference:", max_diff, "\n")
  cat("RMSE difference:", rmse_diff, "\n")
  cat("Relative RMSE:", rel_rmse, "\n")
  
  if (max_diff < 1e-4 && rel_rmse < 0.01) {
    cat(" Gradients are accurate!\n")
  } else {
    cat("Gradient discrepancies detected\n")
  }
  
  # Plot comparison
  plot(analytic_grad, numeric_grad, pch = 19, col = "blue",
       xlab = "Analytic Gradient", ylab = "Numeric Gradient",
       main = "Gradient Verification")
  abline(0, 1, col = "red", lwd = 2)
  
  return(list(analytic = analytic_grad, numeric = numeric_grad, 
              max_diff = max_diff, rmse_diff = rmse_diff, rel_rmse = rel_rmse))
}

## Parameter Utilities and Diagnostics
extract_params <- function(fits_sim_1, p) {
  n_params_eq1 <- 1 + 3*p + 9; n_params_eq2 <- 2 + 3*p + 9; n_params_eq3 <- 3 + 3*p + 9
  phi_names_3p <- if (p > 0) paste0("phi", 1:(3*p)) else character(0)
  param_names_eq1 <- paste0(
    "eq1_", 
    c("alpha", phi_names_3p, "midas_scale_gscpi","midas_scale_ippi", 
      "midas_scale_epi", "theta1_gscpi","theta2_gscpi","theta1_ippi",
      "theta2_ippi","theta1_epi","theta2_epi")
    )
  param_names_eq2 <- paste0(
    "eq2_", 
    c("alpha","phi_contemp", phi_names_3p, "midas_scale_gscpi",
      "midas_scale_ippi","midas_scale_epi","theta1_gscpi","theta2_gscpi",
      "theta1_ippi","theta2_ippi","theta1_epi","theta2_epi")
    )
  param_names_eq3 <- paste0(
    "eq3_", 
    c("alpha","phi_contemp_cpi","phi_contemp_ir", phi_names_3p,
      "midas_scale_gscpi","midas_scale_ippi","midas_scale_epi",
      "theta1_gscpi","theta2_gscpi","theta1_ippi","theta2_ippi",
      "theta1_epi","theta2_epi")
    )
  
  out <- list()
  for (tau_name in names(fits_sim_1)) {
    par <- fits_sim_1[[tau_name]]$opt_result$par
    out[[tau_name]] <- list(
      eq1_params = setNames(par[1:n_params_eq1], param_names_eq1),
      eq2_params = setNames(par[(n_params_eq1 + 1):(n_params_eq1 + n_params_eq2)], param_names_eq2),
      eq3_params = setNames(par[(n_params_eq1 + n_params_eq2 + 1):length(par)], param_names_eq3),
      elapsed_minutes = fits_sim_1[[tau_name]]$elapsed_minutes,
      convergence  = fits_sim_1[[tau_name]]$convergence,
      loss         = fits_sim_1[[tau_name]]$loss
    )
  }
  out
}

print_true_parameters <- function(true_par, p) {
  n_params_eq1 <- 1 + 3*p + 9; n_params_eq2 <- 2 + 3*p + 9
  phi_names_3p <- if (p > 0) paste0("phi", 1:(3*p)) else character(0)
  param_names_eq1 <- c("alpha", phi_names_3p, "midas_scale_gscpi","midas_scale_ippi","midas_scale_epi",
                       "theta1_gscpi","theta2_gscpi","theta1_ippi","theta2_ippi","theta1_epi","theta2_epi")
  param_names_eq2 <- c("alpha","phi_contemp", phi_names_3p, "midas_scale_gscpi","midas_scale_ippi","midas_scale_epi",
                       "theta1_gscpi","theta2_gscpi","theta1_ippi","theta2_ippi","theta1_epi","theta2_epi")
  param_names_eq3 <- c("alpha","phi_contemp_cpi","phi_contemp_ir", phi_names_3p,
                       "midas_scale_gscpi","midas_scale_ippi","midas_scale_epi",
                       "theta1_gscpi","theta2_gscpi","theta1_ippi","theta2_ippi","theta1_epi","theta2_epi")
  
  cat("Equation 1 (CPI):\n"); print(round(setNames(true_par[1:n_params_eq1], param_names_eq1), 4))
  cat("\nEquation 2 (IR):\n"); print(round(setNames(true_par[(n_params_eq1 + 1):(n_params_eq1 + n_params_eq2)], param_names_eq2), 4))
  cat("\nEquation 3 (NHPI):\n"); print(round(setNames(true_par[(n_params_eq1 + n_params_eq2 + 1):length(true_par)], param_names_eq3), 4))
}

diagnose_parameter_recovery <- function(params, true_par, p, abs_floor = 5e-3) {
  n_params_eq1 <- 1 + 3*p + 9; n_params_eq2 <- 2 + 3*p + 9; n_params_eq3 <- 3 + 3*p + 9
  ar_terms1 <- 1:(1 + 3*p);  midas_scales_1 <- (1 + 3*p + 1):(1 + 3*p + 3)  
  theta_gscpi1 <- (1 + 3*p + 4):(1 + 3*p + 9)
  ar_terms2 <- (n_params_eq1 + 1):(n_params_eq1 + 2 + 3*p)
  midas_scales_2 <- (n_params_eq1 + 2 + 3*p + 1):(n_params_eq1 + 2 + 3*p + 3) 
  theta_ippi1 <- (n_params_eq1 + 2 + 3*p + 4):(n_params_eq1 + 2 + 3*p + 9)
  ar_terms3 <- (n_params_eq1 + n_params_eq2 + 1):(n_params_eq1 + n_params_eq2 + 3 + 3*p) 
  midas_scales_3 <- (n_params_eq1 + n_params_eq2 + 3 + 3*p + 1):(n_params_eq1 + n_params_eq2 + 3 + 3*p + 3) 
  theta_epi1 <- (n_params_eq1 + n_params_eq2 + 3 + 3*p + 4):(n_params_eq1 + n_params_eq2 + 3 + 3*p + 9)
  parameter_sets <- list(AR = c(ar_terms1, ar_terms2, ar_terms3), 
                         SCALES = c(midas_scales_1, midas_scales_2, midas_scales_3), 
                         theta = c(theta_gscpi1, theta_ippi1, theta_epi1))
  
  cat("\n=== PARAMETER RECOVERY ===\n")
  for (tau_name in names(params)) {
    est <- c(params[[tau_name]]$eq1_params, params[[tau_name]]$eq2_params, params[[tau_name]]$eq3_params)
    tru <- true_par
    ae  <- abs(est - tru)
    rmse <- sqrt(mean((est - tru)^2, na.rm = TRUE))
    ape  <- ae / pmax(abs(tru), abs_floor) * 100
    mape <- mean(ape, na.rm = TRUE)
    mdap <- median(ape, na.rm = TRUE)
    
    cat("\n", tau_name, ":\n", sep = "")
    cat("Convergence:", params[[tau_name]]$convergence,
        " Loss:", round(params[[tau_name]]$loss, 6),
        " | RMSE:", round(rmse, 5),
        " | MAPE:", round(mape, 2),
        " | MdAPE:", round(mdap, 2), "\n")
    cat("(APE uses abs_floor =", abs_floor, ")\n")
    
    for (set_name in names(parameter_sets)) {
      idx <- parameter_sets[[set_name]]
      set_rmse <- sqrt(mean((est[idx] - tru[idx])^2, na.rm = TRUE))
      set_mae  <- mean(abs(est[idx] - tru[idx]), na.rm = TRUE)
      set_mape <- mean((abs(est[idx] - tru[idx]) / pmax(abs(tru[idx]), abs_floor)) * 100, na.rm = TRUE)
      cat("  ", sprintf("%-6s", set_name), ": RMSE=", round(set_rmse, 5),
          ", MAE=", round(set_mae, 5),
          ", MAPE=", round(set_mape, 2), "\n", sep = "")
    }
  }
}

plot_actual_vs_fitted <- function(actual_series, fitted_by_tau_list,
                                  ylab = "Value", main_prefix = "", n_show = 50) {
  n_show <- min(n_show, length(actual_series))
  actual_vals <- tail(actual_series, n_show)
  tails <- lapply(fitted_by_tau_list, function(x) as.numeric(tail(x, n_show)))
  valid_idx <- vapply(tails, function(v) is.numeric(v) && length(v) == n_show && any(is.finite(v)),
                      logical(1))
  tails <- tails[valid_idx]
  
  if (length(tails) == 0) {
    x <- seq_len(n_show)
    plot(x, actual_vals, type = "l", lwd = 2, col = "black",
         ylab = ylab, xlab = "Quarter (most recent)",
         main = sprintf("%sActual vs Fitted (no valid fitted series)", main_prefix),
         xaxs = "i", yaxs = "i")
    grid(); return(invisible(NULL))
  }
  
  fit_mat <- do.call(cbind, tails)
  colnames(fit_mat) <- names(tails)
  taus_selected <- as.numeric(sub("^tau_", "", names(tails)))
  
  all_vals <- c(actual_vals, as.vector(fit_mat))
  rng <- range(all_vals, na.rm = TRUE)
  pad <- 0.05 * diff(rng); if (!is.finite(pad) || pad == 0) pad <- 0.1
  ylim <- c(rng[1] - pad, rng[2] + pad)
  
  base_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#66a61e",
                 "#e7298a", "#e6ab02", "#a6761d", "#666666")
  cols_use <- base_cols[seq_len(ncol(fit_mat))]
  
  x <- seq_len(n_show)
  plot(x, actual_vals, type = "l", lwd = 3, col = "black",
       ylab = ylab, xlab = "Quarter (most recent)",
       main = sprintf("%sActual vs Fitted Quantiles\n(Last %d Quarters)", main_prefix, n_show),
       ylim = ylim, xaxs = "i", yaxs = "i")
  grid()
  for (i in seq_len(ncol(fit_mat))) lines(x, fit_mat[, i], lwd = 2, col = cols_use[i], lty = 1)
  legend("topright",
         legend = c("Actual", paste0("\u03C4 = ", formatC(taus_selected, format = "f", digits = 2))),
         lty = 1, lwd = c(3, rep(2, length(taus_selected))), 
         col = c("black", cols_use), bty = "n", inset = 0.02, cex = 0.9)
}

create_variable_plots <- function(params, model_data_sim_1) {
  colors <- setNames(c("blue","red","darkgreen","purple")[seq_along(taus)],
                     paste0("tau_", formatC(sort(taus), digits = 2, format = "f")))
  fitted_values <- list()
  for (tau_name in names(params)) {
    par_eq1 <- unname(params[[tau_name]]$eq1_params)
    par_eq2 <- unname(params[[tau_name]]$eq2_params)
    par_eq3 <- unname(params[[tau_name]]$eq3_params)
    par <- c(par_eq1, par_eq2, par_eq3)
    fitted_values[[tau_name]] <- calculate_fitted_values(par, model_data_sim_1, p, K)
  }
  plot(1:nrow(model_data_sim_1), model_data_sim_1$cpi, type = "l", lwd = 2, col = "black",
       xlab = "Time", ylab = "CPI", main = "CPI - Actual vs Fitted by Quantile")
  for (tn in names(params)) lines(1:nrow(model_data_sim_1), fitted_values[[tn]]$cpi_fitted, col = colors[tn], lwd = 1.5, lty = 2)
  legend("topright", legend = c("Actual", names(colors)), 
         col = c("black", colors), lty = c(1, rep(2, length(colors))), lwd = c(2, rep(1.5, length(colors))), bg = "white")
  
  plot(1:nrow(model_data_sim_1), model_data_sim_1$ir, type = "l", lwd = 2, col = "black",
       xlab = "Time", ylab = "Interest Rate", main = "Interest Rate - Actual vs Fitted by Quantile")
  for (tn in names(params)) lines(1:nrow(model_data_sim_1), fitted_values[[tn]]$ir_fitted, col = colors[tn], lwd = 1.5, lty = 2)
  legend("topright", legend = c("Actual", names(colors)), 
         col = c("black", colors), lty = c(1, rep(2, length(colors))), lwd = c(2, rep(1.5, length(colors))), bg = "white")
  
  plot(1:nrow(model_data_sim_1), model_data_sim_1$nhpi, type = "l", lwd = 2, col = "black",
       xlab = "Time", ylab = "NHPI", main = "NHPI - Actual vs Fitted by Quantile")
  for (tn in names(params)) lines(1:nrow(model_data_sim_1), fitted_values[[tn]]$nhpi_fitted, col = colors[tn], lwd = 1.5, lty = 2)
  legend("topright", legend = c("Actual", names(colors)), 
         col = c("black", colors), lty = c(1, rep(2, length(colors))), lwd = c(2, rep(1.5, length(colors))), bg = "white")
}
