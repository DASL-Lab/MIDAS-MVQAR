suppressPackageStartupMessages({
  library(xts)
  library(zoo)
  library(quantreg)
  library(midasr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  theme_set(theme_bw())
  library(ggridges)
  library(stringr)
})

source("functions.R")

##Set Parameters in an easy place to change
m <- 3
K <- 4
p <- 1
taus <- 0.5
stopifnot(length(taus) == 1) # Only works for one tau at a time.
N <- 1000
sim_file_name <- paste0("sim", "-m", m, "-K", K, "-p", p, "-tau", taus, "-N", N, ".csv")
force_rerun <- FALSE

true_par <- create_true_par_vector(p)

res <- matrix(0.0, ncol = length(true_par) + 1, nrow = N + 1)
colnames(res) <- c("sim", names(true_par))
res[1, ] <- c("true", true_par)

if (file.exists(sim_file_name) && !force_rerun) {
  res <- read.csv(sim_file_name)
} else {
  cat("Starting Simulation.\n")

  loop_times <- c()
  t0 <- Sys.time()
  for (i in seq_len(N)) {
    t1 <- Sys.time()
    model_data_sim_1 <- simulate(m = m, K = K, p = p, taus = taus)
    n_obs <- nrow(model_data_sim_1)

    fits_sim <- tryCatch(
      estimate_joint_qarmidas(
        taus = taus, data = model_data_sim_1, K = K, p = p,
        kappa = 0.01,
        theta_shrink = 2e-3,     
        start_par = true_par,
        optimizer = "BFGS",
        use_bounds = FALSE,
        n_restarts = 2,
        jitter_sd = 0.01,
        verbose = FALSE,
        dry_run = FALSE,
        use_analytic_grad = TRUE  # Use our analytic gradients!
      ),
      error = function(e) e
    )
    if (inherits(fits_sim, "error")) next

    params <- extract_params(fits_sim, p)

    res[i + 1, ] <- c(
      paste0("sim", i),
      with(params[[1]], c(eq1_params, eq2_params, eq3_params))
    )

    loop_times[i] <- difftime(Sys.time(), t1, units = "mins")
    mean_loop_time <- mean(loop_times, na.rm = TRUE)
    expected_time <- mean_loop_time * (N - i)
    cat(
      "\rCompleted loop", i, "of", N,
      "in", round(loop_times[i] * 60, 2), "seconds.",
      "Expected time remaining:", round(expected_time, 2), "mins.   "
    )
  }
  cat(
    "\nCompleted in",
    round(difftime(Sys.time(), t0, units = "mins"), 2),
    "mins.\n"
  )
  write.csv(res, file = sim_file_name)
}

res_df <- res |> as.data.frame() |>
  pivot_longer(-sim) |>
  mutate(
    value = as.numeric(value),
    eq = str_extract(name, "eq(\\d)", group = TRUE),
    var = case_when(
      str_detect(name, "ippi") ~ "IPPI",
      str_detect(name, "gcspi") ~ "GSCPI",
      str_detect(name, "epi") ~ "EPI",
      str_detect(name, "cpi") ~ "CPI",
      str_detect(name, "ir") ~ "IR",
      .default = "none"
    ),
    is_theta = str_detect(name, "theta"),
    param = str_split_i(
      str_remove(name, "_midas"), i = 2, pattern = "_"
    ),
    param2 = str_remove(name, "eq\\d_"),
    sim = if_else(sim == "true", "truth", "simulated")
  ) 

gg_bees <- res_df |>
  ggplot() +
  aes(x = value, y = param2,
    colour = sim, shape = sim, size = sim) +
  geom_point() +
  scale_shape_manual(values = c("o", "|")) +
  scale_size_manual(values = c(2, 6)) +
  scale_colour_manual(values = c("black", "red")) +
  facet_wrap(is_theta ~ eq, scales = "free_y", labeller = label_both)
print(gg_bees)


# A different view
res_sim <- filter(res_df, sim != "truth")
res_tru <- filter(res_df, sim == "truth")

gg_box <- ggplot() +
  geom_boxplot( # also try geom_density_ridges(
    data = res_sim,
    mapping = aes(x = value, y = param2)
  ) +
  facet_wrap(is_theta ~ eq, scales = "free", labeller = label_both) +
  geom_point(
    data = res_tru,
    mapping = aes(x = value, y = param2),
    colour = "red", shape = "|", size = 10
  ) +
  coord_cartesian(xlim = c(-1.5, 1.5))
print(gg_box)

