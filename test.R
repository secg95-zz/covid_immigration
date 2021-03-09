"
Simulate several epidemics with the same parameters and gathers goodness of fit
metrics for the available models.
"
library(RJSONIO)
sim = new.env(); source("simulation.R", local=sim)
est = new.env(); source("estimation.R", local=est)
source("discretize_dist.R")

OUT_DIR = "results"

config = list(
  n_epidemics = 10,
  si_dist="weibull",
  si_shape = 3,
  si_scale = 3,
  initial_cases = 100,
  immigration_rate = 15,
  R = rep(1.2, 50)
)
config$time_span = length(config$R)
# R estimations to be compared by Mean Absolute Percentage Error
mape = function(truth, estimator) {
  mean(abs((truth - estimator) / truth), na.rm=TRUE)
}
# Simulate and estimate R
discrete_si = discretize_dist(
  si_dist = config$si_dist,
  si_shape = config$si_shape,
  si_scale = config$si_scale,
  n_vals = config$time_span + 3
)
results = list(
  mean_mape_epiestim = 0,
  mean_mape_epiestim_immigration = 0,
  mean_mape_ekf = 0
)
for (epidemic in 1:config$n_epidemics) {
  # Simulate epidemic
  simulation = sim$simulate_with_immigrants(
    discrete_si = discrete_si,
    initial_cases = config$initial_cases,
    time_span = config$time_span,
    immigration_rate = config$immigration_rate,
    R = config$R
  )
  # Estimate R using each of the available algorithms
  R_epiestim = est$epiestim(
    I = simulation$I,
    discrete_si = discrete_si
  )
  R_epiestim_immigration = est$epiestim_immigration(
    I = simulation$I,
    XI = simulation$XI,
    discrete_si = discrete_si
  )
  ekf_result = est$ekf(
    I = simulation$I,
    discrete_si = discrete_si
  )
  ekf2_result = est$ekf2(
    I = simulation$I,
    discrete_si = discrete_si
  )
  # Collect results
  simulation_results = list(
    I = simulation$I,
    XI = simulation$XI,
    epiestim = list(
      R_hat = R_epiestim,
      mape = mape(config$R, R_epiestim)
    ),
    epiestim_immigration = list(
      R_hat = R_epiestim_immigration,
      mape = mape(config$R, R_epiestim_immigration)
    ),
    ekf = list(
      R_hat = ekf_result$R_hat,
      mape = mape(config$R, ekf_result$R_hat),
      a = ekf_result$a,
      P = ekf_result$P,
      alpha = ekf_result$alpha
    ),
    ekf2 = list(
      R_hat = ekf2_result$R_hat,
      mape = mape(config$R, ekf2_result$R_hat),
      a = ekf2_result$a,
      P = ekf2_result$P,
      alpha_hat = ekf2_result$alpha_hat
    )
  )
  results$simulations[[epidemic]] = simulation_results
  results$mean_mape_epiestim =
    results$mean_mape_epiestim + simulation_results$epiestim$mape
  results$mean_mape_epiestim_immigration =
    results$mean_mape_epiestim_immigration +
    simulation_results$epiestim_immigration$mape
  results$mean_mape_ekf =
    results$mean_mape_ekf + simulation_results$ekf$mape
  results$mean_mape_ekf2 =
    results$mean_mape_ekf2 + simulation_results$ekf2$mape
}
results$mean_mape_epiestim = results$mean_mape_epiestim / config$n_epidemics
results$mean_mape_epiestim_immigration = results$mean_mape_epiestim_immigration / config$n_epidemics
results$mean_mape_ekf = results$mean_mape_ekf / config$n_epidemics
results$mean_mape_ekf2 = results$mean_mape_ekf2 / config$n_epidemics
# Store results in disk
summary = list(
  parameters = config,
  results = results,
  discrete_si = discrete_si
)
dir.create(OUT_DIR, recursive=TRUE)
timestamp = format(Sys.time(), "%Y%m%d%H%M")
out_path = paste0(OUT_DIR, "/", timestamp, ".json")
write(toJSON(summary, pretty=TRUE), out_path)
print(paste0("Test results written at ./", out_path))