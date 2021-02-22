'
Estimators for R given an incidence series.

Parameters
----------
I : double
    Series of incidences.
discrete_si : double
  Discretized distribution of the serial interval. discretete_si[i + 1] :=
  probability of serial interval = i.
XI : double
  Incidence series of immigrants. This one is longer than I, but matches it on
  the right end. In other words, the last elements of both series represent the
  same points in time.

Returns
-------
R_hat : double
  Estimated reproduction number series. Equal in length to I.
a : list
  Filter state expected values. a[[t]] is the expected value of the state vector
  at step t.
P : list
  Filter state covariance matrices. P[[t]] is the covariance matrix of the state
  vector at step t.
alpha : list
  Smoothed state values. alpha[[t]] is the smoothed value of the state vector at
  step t.
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
  # R_hat = exp(EKF$a[1,])
  R_hat = c(rep(NaN, length(I) - dim(EKF$alpha)[2]), exp(EKF$alpha[1,]))
  for (i in 1:length(R_hat)) {
    if (!is.na(R_hat[i]) & (R_hat[i] == Inf | R_hat[i] == -Inf)) {
      R_hat[i] = NaN
    }
  }
  return(list(
    R_hat=R_hat,
    a=lapply(seq_len(ncol(EKF$a)), function(i) EKF$a[,i]),
    P=EKF$P,
    alpha=lapply(seq_len(ncol(EKF$alpha)), function(i) EKF$alpha[,i])
  ))
}

ekf2 = function(I, discrete_si) {
  a = list()
  P = list()
  infectivity = overall_infectivity(I, c(0, discrete_si))
  infectivity[1] = 0 # EpiEstim makes this NA by default
  # Distribution of initial hidden state - might change later
  a[[1]] = c(log(1), log(10), log(1))
  P[[1]] = matrix(
    c(
      c(1, -1, 0), # R and M are negatively correlated
      c(-1, 1, 0),
      c(0, 0, 1)
    ),
    nrow = 3
  )
  for (t in 1:(lenght(I) - 1)) {
    # Apply Kalman filter recursion
    Z_t = exp(a[[t]][1]) * (infectivity[t] + exp(a[[t]][2])) # = Z_t(a_t)
    Z_prime_t = c(
        exp(a[[t]][1]) * (infectivity[t] + exp(a[[t]][2])),
        exp(a[[t]][1]) * exp(a[[t]][2]),
        0
    )
    d_t = Z_t - Z_prime_t %*% a[[t]]
    c_t = c(0, 0, 0)
    v_t = I[t] - Z_t
    F_t = Z_prime_t %*% P[[t]] %*% Z_prime_t
  }
}