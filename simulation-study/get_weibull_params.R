# Find the parameters of a Weibull distribution with a given mean and sd
library(nloptr)

objective_mean = 8.4
objective_sd = 3.8
objective = function(weibull_params) {
  scale = weibull_params[1]
  shape = weibull_params[2]
  weibull_mean = scale * gamma(1 + 1/shape)
  weibull_sd = sqrt((scale^2) * (gamma(1 + 2/shape) - gamma(1 + 1/shape)^2))
  return((objective_mean - weibull_mean)^2 + (objective_sd - weibull_sd)^2)
}
opt = nlminb(c(1, 1), objective, lower = c(0, 0), upper = c(1000, 1000))