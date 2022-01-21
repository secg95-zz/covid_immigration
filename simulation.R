library(EpiEstim)

source("discretize_dist.R")

simulate_weibull = function(initial_cases, time_span, si_scale, si_shape, R) {
  discrete_weibull = diff(pweibull(0:time_span, shape=si_shape, scale=si_scale))
  I = initial_cases
  for (t in 1:(time_span - 1)) {
    next_incidence_mean = R[t] * sum(rev(I) * discrete_weibull[1:t])
    I = c(I, rpois(1, next_incidence_mean))
  }
  return(I)
}

simulate_gamma = function(initial_cases, time_span, si_scale, si_shape, R) {
  discrete_gamma = diff(pgamma(0:time_span, shape=si_shape, scale=si_scale))
  I = initial_cases
  for (t in 1:(time_span - 1)) {
    next_incidence_mean = R[t] * sum(rev(I) * discrete_gamma[1:t])
    I = c(I, rpois(1, next_incidence_mean))
  }
  return(I)
}

simulate_with_immigrants = function(
  discrete_si, T, R, I0, global_R, global_I0, immig_rate
) {
  "
  Simulates an epidemic which starts on a global scale and develops locally due
  to immigration. Its hypotheses are:
  
  1. I[t] ~ Poisson(R[t] * overall_infectivity[t])
  2. global_I[t] ~ Poisson(global_R[t] * global_infectivity[t])
  3. XI[t] ~ Poisson(global_I[t] * immig_rate[t])
  
  Parameters
  ----------
  discrete_si : double
    Discretized distribution of the serial interval. discretete_si[i + 1] :=
    probability of serial interval = i.
  T : integer
    Number of infection timesteps to be simulated.
  global_I0 : integer
    Number of incidences in the first period of the global epidemic.
  R : double
    Series of local reproduction numbers. Must be of length == T.
  I0 : double
    Number of local index cases at the start of the epidemic.
  global_R : double
    Series of global reproduction numbers. Must be of length == T.
  global_I0 : double
    Number of global index cases at the start of the epidemic.
  immig_rate : double
    Series of the rate of the global population that immigrates daily. 

  Returns
  -------
  I : double
    Series of incidences.
  XI : double
    Incidence series of immigrants.
  "
  # Simulate global epidemic
  global_I = global_I0
  for (t in 1:(T - 1)) {
    global_infectivity = overall_infectivity(global_I, c(0, discrete_si))
    next_global_I_mean = global_infectivity[t] * global_R[t]
    global_I = c(global_I, rpois(1, next_global_I_mean))
  }
  # Simulate series of incoming infected immigrants
  XI = NULL
  for (t in 1:T) {
    next_XI_mean = global_I[t] * immig_rate[t]
    XI = c(XI, rpois(1, next_XI_mean))
  }
  # Simulate local epidemic
  I = I0
  for (t in 1:(T - 1)) {
    incid_t = data.frame(list(
      local = I[1:t],
      imported = XI[1:t]
    ))
    local_infectivity = overall_infectivity(incid_t, c(0, discrete_si))
    next_I_mean = R[t] * local_infectivity[t]
    I = c(I, rpois(1, next_I_mean))
  }
  return(list("I"=I, "XI"=XI))
}