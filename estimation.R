'
Estimators for R given an incidence series.
parameters.

I : double
    Series of incidences.
discrete_si : double
  Discretized distribution of the serial interval. discretete_si[i + 1] :=
  probability of serial interval = i.
XI : double
  Incidence series of immigrants. This one is longer than I, but matches it on
  the right end. In other words, the last elements of both series represent the
  same points in time.
R_hat : double
  Estimated reproduction number series. Equal in length to I.
'
library(EpiEstim)

source("raw_code/EKF.R")

epiestim = function(I, discrete_si) {
  R_hat = estimate_R(
    incid = I,
    method = "non_parametric_si", # user-specified discrete distribution
    config = make_config(list(si_distr = c(0, discrete_si)))
  )$R$`Mean(R)`
  # Pad estimation to length of incidence series
  R_hat = c(rep(NA, length(I) - length(R_hat)), R_hat)
  return(R_hat)
}
epiestim_immigration = function(I, XI, discrete_si) {
  R_hat = estimate_R(
    incid = data.frame(list(
      local = c(rep(0, length(XI) - length(I)), I),
      imported = XI
    )),
    method = "non_parametric_si", # user-specified discrete distribution
    config = make_config(list(si_distr = c(0, discrete_si)))
  )$R$`Mean(R)`
  # Adjust estimation length to that of the incidence series
  if (length(I) > length(R_hat)) {
    R_hat = c(rep(NA, length(I) - length(R_hat)), R_hat)
  } else if (length(I) < length(R_hat)) {
    R_hat = R_hat[(length(R_hat) - length(I)):length(R_hat)]
  }
  return(R_hat)
}
ekf = function(I, discrete_si) {
  lags = 7
  lags2 = 5
  cases = I
  casesMatrix_T=matrix(0,length(cases)-lags2+1,lags2)
  for(i in 1:dim(casesMatrix_T)[1]){
    casesMatrix_T[i,]=cases[i:(i+lags2-1)]
  }
  # x_nicolas[i] = x[i + 5]
  # x_nicolas=t(casesMatrix_T %*% discrete_si[5:1])
  x = overall_infectivity(I, c(0, discrete_si))
  par1=c(0,0,0)
  par2=10^6*c(1,1,1)
  par3=c(1,1,1)
  EKF=EKF_Complete(I[2:50], x[2:50], par1,par2,par3)
  R_hat = exp(EKF$a[1,])
  for (i in 1:length(R_hat)) {
    if (R_hat[i] == Inf | R_hat[i] == -Inf) {
      R_hat[i] = NaN
    }
  }
  return(R_hat)
}