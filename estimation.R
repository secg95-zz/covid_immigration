'
Estimators for R given an incidence series.

Parameters
----------
I : double
    Series of incidences.
discrete_si : double
  Discretized distribution of the serial interval. discretete_si[i + 1] :=
  probability of serial interval = i.
infectivity : double
  Infectivity series, i.e. the rolling sum of `I` weighted by `discrete_si`.
XI : double
  Incidence series of immigrants. This one is longer than I, but matches it on
  the right end. In other words, the last elements of both series represent the
  same points in time.
skip_initial : integer
  Number of initial timesteps to skip, e.g. I[t] is discarded for all t <=
  skip_initial.

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
library(geometry)
source("raw_code/EKF.R")

create_sample=function(sqrt_P, a_t, m, k){
  "Creates the sampling for the UKF"
    sigmas = matrix(0, m, 2*m + 1)
    sigmas[,1] = a_t
    w = rep(0, 2*m + 1) 
    for(j in 1:m){
      P_ti = sqrt_P[,j]
      x_ti = a_t + sqrt(m + k)*P_ti
      x_im = a_t - sqrt(m + k)*P_ti
      sigmas[,j + 1] = x_ti
      sigmas[,(m + j) + 1] = x_im
      w[1] = k/(k + m)
      w[j + 1] = 1/(2*(m + k))
      w[m + j + 1] =  1/(2*(m + k))
    }
   avg_at = 0 
   for(j in 1:(2*m+1)){
      xtj = sigmas[,j]
      avg_at = avg_at + w[j]*c(exp(xtj[1]),exp(xtj[2]), 0)
    }
  return(list(sigmas=sigmas, w=w, avg_at=avg_at))
}

create_msample=function(sqrt_P, a_t, m, k){
  "Create the sampling for the MUKF"
    sample_points = c(-0.7, -0.5, 0.5, 0.7)
    num_sam = length(sample_points) 
    sigmas = matrix(0, m, num_sam*m + 1)
    w = rep(0, num_sam*m + 1) 
    sigmas[,1] = a_t
    phis = dnorm(sample_points, mean = 0, sd = 1, log = FALSE)
    aux = dot(phis, sample_points^4)
    aux_2 = dot(phis, sample_points^2)

    w_0 = 1 - ((sum(phis)*aux)/(3*aux_2)^2)
    lambda_square = sum(phis)/((1 - w_0)*aux_2)
    w[1] = w_0
    for(j in 1:m){
      P_ti = sqrt_P[,j]
      for (k in 1:(num_sam/2)){
        w[((j-1)*num_sam) + k + 1] = (1-w_0)*phis[k]/(2*m*sum(phis))
        w[((j-1)*num_sam) + k + num_sam/2 + 1] = (1-w_0)*phis[k]/(2*m*sum(phis))
        x_ti = a_t + sqrt(lambda_square)*sample_points[k]*P_ti
        x_im = a_t - sqrt(lambda_square)*sample_points[k]*P_ti
        sigmas[,((j-1)*num_sam) + k + 1] = x_ti
        sigmas[,((j-1)*num_sam) + k + num_sam/2 + 1] = x_im
      } 
    }
   avg_at = 0 
   for(j in 1:(num_sam*m + 1)){
      xtj = sigmas[,j]
      avg_at = avg_at + w[j]*c(exp(xtj[1]),exp(xtj[2]), 0)
    }
  return(list(sigmas=sigmas, w=w, avg_at=avg_at))
}


epiestim = function(I, discrete_si, ...) {
  R_hat = estimate_R(
    incid = I,
    method = "non_parametric_si", # user-specified discrete distribution
    config = make_config(list(si_distr = c(0, discrete_si)))
  )$R$`Mean(R)`
  # Pad estimation to length of incidence series
  R_hat = c(rep(NA, length(I) - length(R_hat)), R_hat)
  return(list(
    R_hat = R_hat
  ))
}
epiestim_immigration = function(I, XI, discrete_si, ...) {
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
  return(list(
    R_hat = R_hat
  ))
}
ekf = function(I, infectivity, skip_initial = 1, ...) {
  n = length(I)
  EKF = EKF_Complete(
    I[(skip_initial + 1):n],
    infectivity[(skip_initial + 1):n],
    par1 = c(0,0,0), # a_1,
    par2 = 10^6*c(1,1,1), # P_1
    par3 = c(1,1,1) # Q
  )
  # Correct dimensions given skipped timesteps
  n_not_skipped = n - skip_initial
  a = list()
  P = list()
  alpha = list()
  for (t in 1:n) {
    if (t <= skip_initial) {
      a[[t]] = NA
      P[[t]] = NA
      alpha[[t]] = NA
    }
    else {
      a[[t]] = EKF$a[[t - skip_initial]]
      P[[t]] = EKF$P[[t - skip_initial]]
      alpha[[t]] = EKF$alpha[[t - skip_initial]]
    }
  }
  R_hat = unlist(lapply(alpha, function(x) if (is.null(x)) NA else exp(x[1])))
  # Change infinity to NaN to prevent json storage bugs
  for (i in 1:length(R_hat)) {
    if (!is.na(R_hat[i]) & (R_hat[i] == Inf | R_hat[i] == -Inf)) {
      R_hat[i] = NaN
    }
  }
  return(list(
    R_hat=R_hat,
    a=a,
    P=P,
    alpha=alpha,
    Q=EKF$Q
  ))
}

compute_likelihood = function(y, x, avg_a_t, Tmat, P, Q){
  logLikelihood=sum(dpois(y,avg_a_t[1,]*(x + avg_a_t[2,]), log=TRUE))
  for(i in 1:length(y)){
    if(i>0)
      {error = dmvnorm(avg_a_t[,i], mean=avg_a_t[,i], sigma=P[[i]], log=TRUE)}
    else
      {error = dmvnorm(avg_a_t[,i], mean=Tmat%*%avg_a_t[,i-1], sigma=Q, log=TRUE)}
    logLikelihood=logLikelihood + error
  }
  return(logLikelihood)
}

ukf_computation = function(I, infectivity, par1, par2, par3, skip_initial = 4){
  #y,x, par1,par2,par3, m=3
  "Unscentend Kalman filter
  y : list
    Series of observed cases.
  x : double
    True cases, sums of lagged casses and migration.
  par1 : double
    Initial expected value of alpha given y_t (suggested value:0) .
  par2 : double
    Initial variance of alpha given y_t.
  par3 : double
    eta's initial variance, as these errors are uncorralted the suggested value for this paramter is: 1"
  m=3
  aux = length(infectivity)
  x = infectivity[skip_initial:aux]
  y = I[skip_initial:aux]
  a=as.matrix(par1)
  P=list()
  v=list()
  a_tts = list()
  samples_x = list()
  P[[1]]=diag(par2)
  Tmat=diag(3)
  Tmat[1,3]=1
  Q=diag(par3)
  #Heuristic 
  k = 5 - m
  sizeY=length(y)
  avg_a_t = matrix(0, m, sizeY)
  for(i in 1:sizeY){
    a_t=a[,i]
    P_t=P[[i]]
    P_aux = chol(P_t)
    samples_a_t = create_sample(P_aux, a_t, m, k)
    sigmas = samples_a_t$sigmas
    samples_x[[i]] = sigmas
    w = samples_a_t$w
    avg_a_t[,i] = samples_a_t$avg_at
    y_t_expected = 0
    z_t = rep(0, 2*m +1) 
    for (j in 1:(2*m+1)){
      xtj = sigmas[,j]
      z_tj = exp(xtj[1])*(x[i]+ exp(xtj[2]))
      z_t[j] = z_tj
      y_t_expected = y_t_expected + w[j]*z_tj
    }

    P_vvt = 0
    P_alphavt = 0
    for(j in 1:(2*m+1)){
      xtj = sigmas[,j]
      z_tj  = z_t[j]
      P_alphavt = P_alphavt + w[j]*((xtj - a_t)*(z_tj - y_t_expected))
      P_vvt = P_vvt + w[j]*((z_tj - y_t_expected)*(z_tj - y_t_expected))
    }
    v_t = y[i] - y_t_expected
    a_tt = a_t + P_alphavt*(1/P_vvt)*v_t
    a_tts[[i]] = a_tt
    P_tt=P_t-P_alphavt%*%((1/P_vvt)*t(P_alphavt))

    P_aux_n = chol(P_tt)
    samples_a_t1 = create_sample(P_aux_n, a_tt, m, k)
    sigmas_n = samples_a_t1$sigmas
    w_n = samples_a_t1$w

    a_t1=0
    for(j in 1:(2*m+1)){
      xtj_n = sigmas_n[,j]
      a_t1 = a_t1 + w_n[j]*(Tmat%*%xtj_n)
    }
    P_t1 = 0
    for(j in 1:(2*m+1)){
      xtj_n = sigmas_n[,j]
      P_t1 =  P_t1 + w_n[j]*((Tmat%*%xtj_n - a_t1)%*%t(Tmat%*%xtj_n - a_t1)) + w_n[j]*Q
    }
    a=cbind(a,a_t1)
    P[[i+1]]=P_t1
    v[[i]]=v_t
  }
  # Pad R_hat to expected length
  n = length(I) - skip_initial
  R_hat = c(rep(NA, skip_initial), avg_a_t[1, 1:n])
  logLikelihood = compute_likelihood(y, x, avg_a_t, Tmat, P, Q)  
  print(logLikelihood)
  return(list(
    R_hat=R_hat,
    a=a,
    P=P,
    alpha_hat=a,
    logLikelihood=logLikelihood
  ))
} 

inicializadorLogLik=function(I, infectivity, par1, par2, par3, kfilter, skip_initial=4){
  "Bayesian search for the best hyperparameters.
  I : double
      Series of incidences.
  infectivity : double
    Infectivity series, i.e. the rolling sum of `I` weighted by `discrete_si`.
  par1 : double
    Initial expected value of alpha given y_t (suggested value:0) .
  par2 : double
    Initial variance of alpha given y_t.
  par3 : double
    eta's initial variance, as these errors are uncorralted the suggested value for this paramter is: 1

  Returns
  -------
  par1 : double
    Initial expected value of alpha given y_t (suggested value:0) .
  par2 : double
    Initial variance of alpha given y_t.
  par3 : double
    eta's initial variance, as these errors are uncorralted the suggested value for this paramter is: 1.
  "
  pars=c(par1, par2, par3)
  obj = function(par) {
    par1=par[1:3]
    par2=par[4:6]
    par3=par[7:9]
    kfilter_result=kfilter(I, infectivity, par1, par2, par3, skip_initial=skip_initial)
    objectiveFunction=-kfilter_result$logLikelihood       
    return(objectiveFunction)
  }
  #TODO : a sin restriccion, para P y Q solo positivos.
  lb=c(0, 0, 0, 0.09, 0.09, 0.09, 0.8, 0.8, 0.8)
  ub=c(rep(Inf, 3), rep( 0.11, 3), rep( Inf, 3))
  initial_condintions=pars
  opt = nlminb(initial_condintions, obj , lower = lb, upper = ub)
  pars=opt$par
  par1=pars[1:3]
  par2=pars[4:6]
  par2=pars[7:9]
  return(list(par1=par1, par2=par2, par3=par3))
}

ukf=function(I, infectivity, skip_initial = 4){
  par1 = c(  0,   0,   0)
  par2 = c(0.1, 0.1, 0.1)
  par3 = c(  1,   1,   1)
  parm=inicializadorLogLik(I, infectivity, par1, par2, par3, ukf_computation, skip_initial)
  ukf_result=ukf_computation(I, infectivity, parm$par1, parm$par2, parm$par3, skip_initial)
  return(ukf_result)
}

mukf_computation = function(I, infectivity, par1, par2, par3, skip_initial = 4){
  #y,x, par1,par2,par3, m=3
  "Unscentend Kalman filter
  y : list
    Series of observed cases.
  x : double
    True cases, sums of lagged casses and migration.
  par1 : double
    Initial expected value of alpha given y_t (suggested value:0) .
  par2 : double
    Initial variance of alpha given y_t.
  par3 : double
    eta's initial variance, as these errors are uncorralted the suggested value for this paramter is: 1"
  m=3
  aux = length(infectivity)
  x = infectivity[skip_initial:aux]
  y = I[skip_initial:aux]
  a=as.matrix(par1)
  P=list()
  v=list()
  a_tts = list()
  samples_x = list()
  P[[1]]=diag(par2)
  Tmat=diag(3)
  Tmat[1,3]=1
  Q=diag(par3)
  #Heuristic 
  k = 5 - m
  sizeY=length(y)
  avg_a_t = matrix(0, m, sizeY)
  for(i in 1:sizeY){
    a_t=a[,i]
    P_t=P[[i]]
    P_aux = chol(P_t)
    samples_a_t = create_msample(P_aux, a_t, m, k)
    sigmas = samples_a_t$sigmas
    samples_x[[i]] = sigmas
    w = samples_a_t$w
    avg_a_t[,i] = samples_a_t$avg_at
    y_t_expected = 0
    z_t = rep(0, length(w)) 
    for (j in 1:length(w)){
      xtj = sigmas[,j]
      z_tj = exp(xtj[1])*(x[i]+ exp(xtj[2]))
      z_t[j] = z_tj
      y_t_expected = y_t_expected + w[j]*z_tj
    }

    P_vvt = 0
    P_alphavt = 0
    for(j in 1:length(w)){
      xtj = sigmas[,j]
      z_tj  = z_t[j]
      P_alphavt = P_alphavt + w[j]*((xtj - a_t)*(z_tj - y_t_expected))
      P_vvt = P_vvt + w[j]*((z_tj - y_t_expected)*(z_tj - y_t_expected))
    }
    v_t = y[i] - y_t_expected
    a_tt = a_t + P_alphavt*(1/P_vvt)*v_t
    a_tts[[i]] = a_tt
    P_tt=P_t-P_alphavt%*%((1/P_vvt)*t(P_alphavt))
 
    P_aux_n = chol(P_tt)
    samples_a_t1 = create_msample(P_aux_n, a_tt, m, k)
    sigmas_n = samples_a_t1$sigmas
    w_n = samples_a_t1$w

    a_t1=0
    for(j in 1:length(w)){
      xtj_n = sigmas_n[,j]
      a_t1 = a_t1 + w_n[j]*(Tmat%*%xtj_n)
    }
    P_t1 = 0
    for(j in 1:length(w)){
      xtj_n = sigmas_n[,j]
      P_t1 =  P_t1 + w_n[j]*((Tmat%*%xtj_n - a_t1)%*%t(Tmat%*%xtj_n - a_t1)) + w_n[j]*Q
    }
    a=cbind(a,a_t1)
    P[[i+1]]=P_t1
    v[[i]]=v_t
  }
  # Pad R_hat to expected length
  n = length(I) - skip_initial
  R_hat = c(rep(NA, skip_initial), avg_a_t[1, 1:n])

  logLikelihood = compute_likelihood(y, x, avg_a_t, Tmat, P, Q)

  return(list(
    R_hat=R_hat,
    a=a,
    P=P,
    alpha_hat=a,
    logLikelihood=logLikelihood
  ))
}

mukf=function(I, infectivity, skip_initial = 4, ...){
  par1 = c(  0,   0,   0)
  par2 = c(0.1, 0.1, 0.1)
  par3 = c(  1,   1,   1)
  parm=inicializadorLogLik(I, infectivity, par1, par2, par3, mukf_computation, skip_initial)
  mukf_result=mukf_computation(I, infectivity, parm$par1, parm$par2, parm$par3, skip_initial)
  return(mukf_result)
}

# library(matlib)
# ekf2 = function(I, infectivity, skip_initial = 1, ...) {
#   # Extract ground-truth information
#   n = length(I)
#   theta = I
#   # Model matrices
#   T_prime = matrix(
#     c(
#       c(1, 0, 1), # R and M are negatively correlated
#       c(0, 1, 0),
#       c(0, 0, 1)
#     ),
#     nrow = 3,
#     byrow=TRUE
#   )
#   R = diag(3)
#   Q = diag(3)
#   # Initialize state variables
#   a = list()
#   P = list()
#   alpha_hat = list()
#   for (t in 1:skip_initial) {
#     a[[t]] = NA
#     P[[t]] = NA
#     alpha_hat[[t]] = NA
#   }
#   a[[skip_initial + 1]] = matrix(c(0, 0, 0), ncol=1)
#   P[[skip_initial + 1]] = matrix(
#     c(
#       c(1, 0, -.5), # R_t and M_t are
#       c(-.5, 1, 0),
#       c(0, 0, 1)
#     ),
#     nrow = 3,
#     byrow = TRUE
#   )
#   # Store matrices that we'll need for smoothing
#   Z = list() # = Z_t(a_t)
#   Z_prime = list() # = Z'_t %*% a_t
#   F_inv = list()
#   K = list()
#   v = list()
#   for (t in (skip_initial + 1):n) {
#     # Calculate variables for Kalman filter
#     Z[[t]] = exp(a[[t]][1]) * (infectivity[t] + exp(a[[t]][2])) # = Z_t(a_t)
#     Z_prime[[t]] = matrix(c(
#         exp(a[[t]][1]) * (infectivity[t] + exp(a[[t]][2])),
#         exp(a[[t]][1]) * exp(a[[t]][2]),
#         0
#     ), nrow=1) # = Z'_t %*% a_t
#     d_t = Z[[t]] - Z_prime[[t]] %*% a[[t]]
#     c_t = matrix(c(0, 0, 0), ncol=1)
#     v[[t]] = theta[t] - Z[[t]]
#     F_t = Z_prime[[t]] %*% P[[t]] %*% t(Z_prime[[t]])
#     # Matrix inversion might fail
#     F_inv[[t]] = tryCatch(
#       {solve(F_t)},
#       error = function(e) {
#         print(paste0("F_t = ", F_t, "; system is singular (ekf2)"))
#       },
#       finally = {
#         return(list(R_hat=NULL, a=a, P=P, alpha_hat=NULL))
#       }
#     )
#     # Kalman Filter recursion
#     a_t_t = a[[t]] + (P[[t]] %*% t(Z_prime[[t]]) %*% F_inv[[t]] %*% v[[t]])
#     a[[t + 1]] = T_prime %*% a_t_t
#     P_t_t = P[[t]] - (P[[t]] %*% t(Z_prime[[t]]) %*% F_inv[[t]] %*% Z_prime[[t]] %*% P[[t]])
#     P[[t + 1]] = (T_prime %*% P_t_t %*% t(T_prime)) + (R %*% Q %*% t(R))
#   }
#   # Kalman Smoothing
#   r_t_prev = matrix(c(0, 0, 0), ncol=1)
#   for (t in n:(skip_initial + 1)) {
#     # Update r_t_prev
#     K_t = T_prime %*% P[[t]] %*% t(Z_prime[[t]]) %*% F_inv[[t]]
#     L_t = T_prime - (K_t %*% Z_prime[[t]])
#     r_t_prev = (t(Z_prime[[t]]) %*% F_inv[[t]] %*% v[[t]]) + (t(L_t) %*% r_t_prev)
#     # Calculate next (or previous) smoothed state
#     alpha_hat[[t]] = a[[t]] + (P[[t]] %*% r_t_prev)
#   }
#   return(list(
#     R_hat=unlist(lapply(a, function(x) if (is.numeric(x)) exp(x[1]) else NA)),
#     a=a[1:n],
#     P=P[1:n],
#     alpha_hat=alpha_hat[1:n]
#   ))
# }