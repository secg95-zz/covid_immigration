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
  si_dist="weibull", si_shape, si_scale, initial_cases, time_span,
  immigration_rate, R
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
  si_dist : character
    Distribution of the serial interval. `weibull` and `gamma` are available.
  si_shape : double
    Shape parameter of the serial interval.
  si_scale : double
    Scale parameter of the serial interval.
  initial_cases : integer
    Number of infectious immigrants arriving at the begining of the epidemic.
  time_span : integer
    Number of epidemic periods to be simulated.
  immigration_rate : double
    Mean number of infectious immigrants arriving at any period.
  R : double
    Series of reproduction numbers. Must be of length >= time_span.
    
  Returns
  -------
  All received parameters as well as
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
  discrete_si = discretize_dist(
    si_dist = si_dist,
    si_shape = si_shape,
    si_scale = si_scale,
    n_vals = time_span + imm_infection_window
  )
  I = 0
  X = initial_cases # Immigrant arrival series
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
    X = c(X, rpois(1, immigration_rate))
    XI = c(XI, 0)
    for (i in 1:tail(X, 1)) {
      onset = length(XI) - sample(1:imm_infection_window, 1) + 1
      XI[onset] = XI[onset] + 1
    }
  }
  return(
    list(
      "I"=I, "XI"=XI, "X"=X,
      "si_dist"=si_dist, "si_shape"=si_shape, "si_scale"=si_scale,
      "initial_cases"=initial_cases, "time_span"=time_span,
      "immigration_rate"=immigration_rate, "R"=R
    )
  )
}