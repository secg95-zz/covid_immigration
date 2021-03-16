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
library(matlib)
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
  # lags = 7
  # lags2 = 5
  # cases = I
  # casesMatrix_T=matrix(0,length(cases)-lags2+1,lags2)
  # for(i in 1:dim(casesMatrix_T)[1]){
  #   casesMatrix_T[i,]=cases[i:(i+lags2-1)]
  # }
  # x_nicolas[i] = x[i + 5]
  # x_nicolas=t(casesMatrix_T %*% discrete_si[5:1])
  n = length(I)
  x = overall_infectivity(I, c(0, discrete_si))
  par1=c(0,0,0) # a_1
  par2=10^6*c(1,1,1) # P_1
  par3=c(1,1,1) # Q
  EKF=EKF_Complete(I[2:n], x[2:n], par1,par2,par3)
  # Correct dimensions: We skipped t=1 and ignore t=n+1
  a = matrix(NA, nrow=3, ncol=n)
  a[2:n] = EKF$a[1:(n - 1)]
  alpha = matrix(NA, nrow=3, ncol=n)
  alpha[2:n] = EKF$alpha[1:(n - 1)]
  R_hat = exp(alpha[1,])
  # Change infinity to NaN to prevent json storage bugs
  for (i in 1:length(R_hat)) {
    if (!is.na(R_hat[i]) & (R_hat[i] == Inf | R_hat[i] == -Inf)) {
      R_hat[i] = NaN
    }
  }
  return(list(
    R_hat=R_hat,
    a=lapply(seq_len(ncol(EKF$a)), function(i) EKF$a[,i]),
    P=EKF$P,
    alpha=lapply(seq_len(ncol(alpha)), function(i) alpha[,i])
  ))
}

ekf2 = function(I, discrete_si) {
  # Extract ground-truth information
  n = length(I)
  theta = I
  infectivity = overall_infectivity(I, c(0, discrete_si))
  # Skip t=1 because infectivity isn't defined
  theta = theta[2:n]
  infectivity = infectivity[2:n]
  # Model matrices
  T_prime = matrix(
    c(
      c(1, 0, 1), # R and M are negatively correlated
      c(0, 1, 0),
      c(0, 0, 1)
    ),
    nrow = 3
  )
  R = diag(3)
  Q = diag(3)
  # Initialize state variables
  a = list(NULL)
  P = list(NULL)
  a[[2]] = matrix(c(0, log(10), 0), ncol=1)
  P[[2]] = matrix(
    c(
      c(1, -.5, 0), # R_t and M_t are
      c(-.5, 1, 0),
      c(0, 0, 1)
    ),
    nrow = 3
  )
  # Store matrices that we'll need for smoothing
  Z = list() # = Z_t(a_t)
  Z_prime = list() # = Z'_t %*% a_t
  F_inv = list()
  K = list()
  v = list()
  for (t in 2:n) {
    # Calculate variables for Kalman filter
    Z[[t]] = exp(a[[t]][1]) * (infectivity[t] + exp(a[[t]][2])) # = Z_t(a_t)
    Z_prime[[t]] = matrix(c(
        exp(a[[t]][1]) * (infectivity[t] + exp(a[[t]][2])),
        exp(a[[t]][1]) * exp(a[[t]][2]),
        0
    ), nrow=1) # = Z'_t %*% a_t
    d_t = Z[[t]] - Z_prime[[t]] %*% a[[t]]
    c_t = matrix(c(0, 0, 0), ncol=1)
    v[[t]] = theta[t] - Z[[t]]
    F_t = Z_prime[[t]] %*% P[[t]] %*% t(Z_prime[[t]])
    F_inv[[t]] = solve(F_t)
    # Kalman Filter recursion
    a_t_t = a[[t]] + (P[[t]] %*% t(Z_prime[[t]]) %*% F_inv[[t]] %*% v[[t]])
    a[[t + 1]] = T_prime %*% a_t_t
    P_t_t = P[[t]] - (P[[t]] %*% t(Z_prime[[t]]) %*% F_inv[[t]] %*% Z_prime[[t]] %*% P[[t]])
    P[[t + 1]] = (T_prime %*% P_t_t %*% t(T_prime)) + (R %*% Q %*% t(R))
  }
  # Kalman Smoothing
  alpha_hat = list(NULL)
  r_t_prev = matrix(c(0, 0, 0), ncol=1)
  for (t in n:2) {
    # Update r_t_prev
    K_t = T_prime %*% P[[t]] %*% t(Z_prime[[t]]) %*% F_inv[[t]]
    L_t = T_prime - (K_t %*% Z_prime[[t]])
    r_t_prev = (t(Z_prime[[t]]) %*% F_inv[[t]] %*% v[[t]]) + (t(L_t) %*% r_t_prev)
    # Calculate next (or previous) smoothed state
    alpha_hat[[t]] = a[[t]] + (P[[t]] %*% r_t_prev)
  }
  return(list(
    R_hat=unlist(lapply(a, function(x) exp(x[1]))),
    a=a[1:n],
    P=P[1:n],
    alpha_hat=alpha_hat
  ))
}