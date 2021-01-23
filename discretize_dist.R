discretize_dist = function(si_dist, si_shape, si_scale, n_vals) {
  if (si_dist == "weibull") {
    discrete_si = diff(pweibull(0:n_vals, shape=si_shape, scale=si_scale))
  } else if (si_dist == "gamma") {
    discrete_si = diff(pgamma(0:n_vals, shape=si_shape, scale=si_scale))
  } else {
    stop(paste("SI distribution `", si_dist, "` unavailable."))
  }
}
