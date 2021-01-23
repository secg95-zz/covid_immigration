'
Estimators for R given an incidence series. All estimators accept the same
parameters.

I : double
    Series of incidences.
discrete_si : double
  Discretized distribution of the serial interval. discretete_si[i + 1] :=
  probability of serial interval = i.
'
library(EpiEstim)

epiestim = function(I, discrete_si) {
  R_hat = estimate_R(
    incid = I,
    method = "non_parametric_si", # user-specified discrete distribution
    config = make_config(list(si_distr = c(0, discrete_si)))
  )
  return(R_hat)
}

ekf = function(I, si_dist, si_scale, si_shape) {

  
}