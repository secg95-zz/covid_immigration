source("discretize_dist.R")

epiestim = function(I, si_dist, si_scale, si_shape) {
  discrete_si = discretize_dist(
    si_dist = si_dist,
    si_shape = si_shape,
    si_scale = si_scale,
    n_vals = length(I)
  )
  R_hat = estimate_R(
    incid = I,
    method = "non_parametric_si", # user-specified discrete distribution
    config = make_config(list(si_distr = c(0, discrete_si)))
  )
  return(R_hat)
}