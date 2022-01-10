"
Simulate several epidemics under different scenarios, and gather goodness of fit
metrics for the available models
"
library(RJSONIO)
source("discretize_dist.R")
SIM = new.env(); source("simulation.R", local=SIM)
EST = new.env(); source("estimation.R", local=EST)
ESTIMATORS = list(
  epiestim = EST$epiestim,
  epiestim_immigration = EST$epiestim_immigration,
  ekf = EST$ekf,
  mukf = EST$mukf
)

OUT_DIR = "simulation-study/results"
CONFIG_FILE = "simulation-study/config.json"

config = fromJSON(CONFIG_FILE)

# Comparison metric: Mean Absolute Percentage Error
mape = function(truth, estimator, mape_starts_at) {
  n = length(truth)
  mean(
    abs((truth[mape_starts_at:n] - estimator[mape_starts_at:n]) / truth[mape_starts_at:n]),
    na.rm=FALSE
  )
}
# Discretize serial interval
discrete_si = discretize_dist(
  si_dist = config$si_dist,
  si_shape = config$si_shape,
  si_scale = config$si_scale,
  n_vals = config$T + 3
)
# Simulate multiple epidemics per scenario
results = list()
for (scenario in config$scenarios) {
  print(paste("Simulating epidemics for scenario", scenario$id))
  scenario_start_time = Sys.time()
  scenario_results = list(
    mean_mape = list(),
    epidemic_results = list()
  )
  for (estimator_name in names(ESTIMATORS)) {
    scenario_results$mean_mape[[estimator_name]] = 0
  }
  for (epidemic in 1:config$simulations_per_scenario) {
    # Simulate epidemic
    simulation = SIM$simulate_with_immigrants(
      discrete_si = discrete_si,
      time_span = config$T,
      index_cases = config$index_cases,
      immigration_rate = scenario$immigration_rate,
      R = scenario$R
    )
    infectivity = overall_infectivity(simulation$I, c(0, discrete_si))
    # Estimate R using each of the available algorithms
    simulation_results = list(
      I = simulation$I,
      XI = simulation$XI,
      R_hat = list(),
      mape = list()
    )
    for (estimator_name in names(ESTIMATORS)) {
      # Estimate R_hat
      estimator = ESTIMATORS[[estimator_name]]
      if (estimator_name %in% c("ekf", "mukf")) {
        skip_initial = 2
      } else {
        skip_initial = 0
      }
      estimation_results = estimator(
        I = simulation$I,
        XI = simulation$XI,
        discrete_si = discrete_si,
        infectivity = infectivity,
        skip_initial = skip_initial
      )
      simulation_results$R_hat[[estimator_name]] = estimation_results$R_hat
      simulation_results$raw_output[[estimator_name]] = estimation_results
      simulation_results$mape[[estimator_name]] = mape(
        scenario$R,
        simulation_results$R_hat[[estimator_name]],
        config$mape_starts_at
      )
      # Accumulate MAPE
      scenario_results$mean_mape[[estimator_name]] =
        scenario_results$mean_mape[[estimator_name]] +
        simulation_results$mape[[estimator_name]]
    }
    scenario_results$epidemic_results[[epidemic]] = simulation_results
  }
  scenario_results$mean_mape = lapply(
    scenario_results$mean_mape,
    function(x) {x / config$simulations_per_scenario}
  )
  results[[scenario$id]] = scenario_results
  results[[scenario$id]][["scenario_id"]] = scenario$id
  scenario_total_time = Sys.time() - scenario_start_time
  print(paste("Processing scenario", scenario$id, "took",
        round(scenario_total_time / 3600, digits=2), "hours"))
}
# Write results to disk
dir.create(OUT_DIR, recursive=TRUE)
timestamp = format(Sys.time(), "%Y%m%d%H%M")
out_path = paste0(OUT_DIR, "/", timestamp, ".json")
write(toJSON(results, pretty=TRUE), out_path)
print(paste0("Test results written at ./", out_path))