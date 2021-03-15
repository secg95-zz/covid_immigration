library(readxl)
library(expm)
library(mvtnorm)
library(nloptr)
library(reshape2)
library(np)

create_sample=function(sqrt_P, a_t, m, k){
  "The is aproblem here, the ones abowe weights more, exponential function! ponerle logs????"
    sigmas = list()
    sigmas[[1]] = a_t
    w = rep(0, 2*m + 1) 
    for(j in 1:m){
      P_ti = sqrt_P[,j]
      x_ti = a_t + sqrt(m + k)*P_ti
      x_im = a_t - sqrt(m + k)*P_ti
      sigmas[[j + 1]] = x_ti
      sigmas[[(m + j) + 1]] = x_im
      w[1] = k/(k + m)
      w[j + 1] = 1/(2*(m + k))
      w[m + j + 1] =  1/(2*(m + k))
    }
  return(list(sigmas=sigmas,w=w))
}

UKF_FS=function(y,x, par1,par2,par3, m=3){
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
  a=as.matrix(par1)
  P=list()
  v=list()
  a_tts = list()
  P[[1]]=diag(par2)
  Tmat=diag(3)
  Tmat[1,3]=1
  Q=diag(par3)
  #Heuristic 
  k = 5 - m
  sizeY=length(y)
  samples_x = list()
  for(i in 1:sizeY){
    a_t=a[,i]
    P_t=P[[i]]
    P_aux = chol(P_t)
    samples_a_t = create_sample(P_aux, a_t, m, k)
    samples_x[[i]] = samples_a_t$sigmas
    sigmas = samples_a_t$sigmas
    w = samples_a_t$w
    y_t_expected = 0
    z_t = rep(0, 2*m +1) 
    for (j in 1:(2*m+1)){
      xtj = sigmas[[j]]
      z_tj = exp(xtj[1])*(x[i]+ exp(xtj[2]))
      print("z_j")
      print(z_tj)
      z_t[j] = z_tj
      y_t_expected = y_t_expected + w[j]*z_tj
    }

    P_vvt = 0
    P_alphavt = 0
    for(j in 1:(2*m+1)){
      xtj = sigmas[[j]]
      z_tj  = z_t[j]
      P_alphavt = P_alphavt + w[j]*((xtj - a_t)*(z_tj - y_t_expected))
      P_vvt = P_vvt + w[j]*((z_tj - y_t_expected)*(z_tj - y_t_expected))
    }
    
    print("predictions")
    print(y[i])
    print(y_t_expected)
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
      xtj_n = sigmas_n[[j]]
      a_t1 = a_t1 + w_n[j]*(Tmat%*%xtj_n)
    }
    P_t1 = 0
    for(j in 1:(2*m+1)){
      xtj_n = sigmas_n[[j]]
      P_t1 =  P_t1 + w_n[j]*((Tmat%*%xtj_n - a_t1)%*%t(Tmat%*%xtj_n - a_t1)) + w_n[j]*Q
    }
    a=cbind(a,a_t1)
    print("diferences")
    print(a_t1 - a_t)
    P[[i+1]]=P_t1
    v[[i]]=v_t
  }
  print("expected before smoothing")
  print(a)
  print(exp(a[1,]))
  print(x+exp(a[2,]))
  print(t(y))
  #smoothing

  alpha=matrix(0,3,sizeY)
  for(i in rev(2:sizeY)){
    C_t1 = 0
    a_t1 = a[i]
    a_t = a[i - 1]
    sigmas = samples_x[i - 1]
    for(j in 1:(2*m+1)){
      xti = sigmas[[1]][[j]]
      C_t1 = C_t1 + w[j]*((xti - a_t)%*%t(Tmat%*%xti - a_t1))
    }
    alpha[,i - 1] = a_tts[[i - 1]] + C_t1%*%(P[[i]]^-1)%*%(alpha[,i] - a_t1)
  }
  #TODO: compute Poisson Maximun likelihood
  logLikelihood=sum(dpois(y,exp(alpha[1,])*(x+exp(alpha[2,])), log=TRUE))
  for(i in 1:sizeY){}
    #logLikelihood=logLikelihood+if(i>0){dmvnorm(alpha[,i], mean=a[,i], sigma=P[[i]], log=TRUE)}else{dmvnorm(alpha[,i], mean=Tmat%*%alpha[,i-1], sigma=Q, log=TRUE)}
    #}
  return(list(alpha=alpha,a=a,P=P,logL=logLikelihood))
}

inicializadorLogLik=function(y,x, par1,par2,par3){
  "
  Function that starts the hyperparameter tunning for par1, par2 and par3
  "
  pars=c(par1,par2,par3)
  obj = function(par) {
    par1=par[1:3]
    par2=exp(par[4:6])
    par3=par[7:9]
    UKF=UKF_FS(y,x, par1,par2,par3)
    objectiveFunction=-UKF$logL        
    return(objectiveFunction)
  }
  
  lb=rep(-5,9)
  ub=rep(5,9)
  ini=pars
  opt = nlminb(ini, obj, lower = lb, upper = ub)
  pars=opt$par
  par1=pars[1:3]
  par2=exp(pars[4:6])
  par2=pars[7:9]
  return(list(par1=par1,par2=par2,par3=par3))
}

UKF_Complete=function(y,x, par1,par2,par3){
  "
  Run the optimizer
  "
  parm=inicializadorLogLik(y,x, par1,par2,par3)
  UKF=UKF_FS(y,x, parm$par1,parm$par2,parm$par3)
  return(UKF)
}




