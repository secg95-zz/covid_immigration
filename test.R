"
Simulates several epidemics with the same parameters and gathers goodness of fit
metrics for the available models.
"
library(RJSONIO)
sim = new.env(); source("simulation.R", local=sim)
est = new.env(); source("estimation.R", local=est)
source("discretize_dist.R")

OUT_DIR = "results"

config = list(
  n_epidemics = 100,
  si_dist="weibull",
  si_shape = 3,
  si_scale = 3,
  initial_cases = 10,
  time_span = 50,
  immigration_rate = 2,
  R = rep(1.2, 50)
)
# Simulate and estimate R
discrete_si = discretize_dist(
  si_dist = config$si_dist,
  si_shape = config$si_shape,
  si_scale = config$si_scale,
  n_vals = config$n_epidemics
)
results = NULL
for (. in 1:config$n_epidemics) {
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
  )$R$`Mean(R)`
  # Collect estimations and goodness of fit measures
  results = c(
    results,
    list(list(
      I = simulation$I,
      XI = simulation$XI,
      R_epiestim = R_epiestim
    ))
  )
}
# Store results in disk
summary = list(parameters = config, results = results)
dir.create(OUT_DIR, recursive=TRUE)
timestamp = format(Sys.time(), "%Y%m%d%H%M")
write(toJSON(summary), paste0(OUT_DIR, "/", timestamp, ".json"))