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
  discrete_si, time_span, index_cases, immigration_rate, R
) {
  "
  Simulates an incidence series with immigrants in which
  1. I[t], the incidences at t, is a Poisson variable.
  2. X[t], the arriving infectious immigrants at time t, is a Poisson variable.
  3. An infectious immigrant arriving at t may have become infectious within a
     window of length `ìmm_infectious_period` previous to their arrival, and
     their time of infection is uniformly distributed within that window.
     
  Parameters
  ----------
  discrete_si : double
    Discretized distribution of the serial interval. discretete_si[i + 1] :=
    probability of serial interval = i.
  time_span : integer
    Number of infection timesteps to be simulated.
  index_cases : integer
    Number of incidences in the first period of the epidemic.
  immigration_rate : double
    Series of E[X[t]].
  R : double
    Series of reproduction numbers. Must be of length >= time_span.
    
  Returns
  -------
  I : double
    Series of incidences.
  X : double
    Series of arriving infectious immigrants.
  XI : double
    Incidence series of immigrants. This one is longer than the other series,
    but matches with them on the right end. In other words, the last elements of
    all series represent the same point in time.
  "
  imm_infection_window = 5
  I = index_cases
  X = rpois(1, immigration_rate[1]) # Immigrant arrival series
  XI = rep(0, imm_infection_window) # Immigrant incidence series
  for (i in 1:X) {
    onset = length(XI) - sample(1:imm_infection_window, 1) + 1
    XI[onset] = XI[onset] + 1
  }
  for (t in 1:(time_span - 1)) {
    infected_by_locals = R[t] * sum(rev(I) * discrete_si[1:t])
    infected_by_imms = R[t] * sum(
      rev(XI) * discrete_si[1:(t + imm_infection_window - 1)]
      )
    next_incidence_mean = infected_by_locals + infected_by_imms
    I = c(I, rpois(1, next_incidence_mean))
    X = c(X, rpois(1, immigration_rate[t + 1]))
    XI = c(XI, 0)
    for (i in 1:tail(X, 1)) {
      onset = length(XI) - sample(1:imm_infection_window, 1) + 1
      XI[onset] = XI[onset] + 1
    }
  }
  return(
    list(
      "I"=I, "XI"=XI, "X"=X
    )
  )
}