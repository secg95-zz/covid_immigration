library(readxl)
library(expm)
library(mvtnorm)
library(nloptr)
library(reshape2)
library(np)
library(geometry) 


#Funciones: 
weightConstructionIncubacionInfeccion=function(lags=7, par1_inf=24.206087, par2_inf=2.984198, distribucion_inf ){
  texto=paste0("p",distribucion_inf,"(1:lags, par1_inf, par2_inf, lower.tail = FALSE)")
  omegaLag=rev(eval(parse(text = texto)))
  omegaLag=omegaLag*valorEsperado(par1=par1_inf, par2=par2_inf, distribucion=distribucion_inf)/sum(omegaLag)
  return(omegaLag)
}
weightConstructionIncubacion=function(lags=7, par1_inc=3.169434, par2_inc=5.163921, distribucion_inc){
  texto=paste0("p",distribucion_inc,"(0:lags,par1_inc, par2_inc)")
  omegaLag=diff(eval(parse(text = texto)))
  omegaLag=omegaLag/sum(omegaLag)
  return(omegaLag)
}
valorEsperado=function(par1=3.169434, par2=5.163921, distribucion="gamma"){
  texto=paste0("d",distribucion,"(x,par1, par2)")
  f_t=function(x){
    return(x*eval(parse(text = texto)))
  }
  t=integrate(f_t,0,Inf)$value
  return(t)
}

EKF_FS=function(y,x, par1,par2,par3){
  a=as.matrix(par1)
  P=list()
  dz=list()
  Fmat=list()
  v=list()
  C = list()
  
  P[[1]]=diag(par2)
  
  Tmat=diag(3)
  Tmat[1,3]=1
  Q=diag(par3)
  
  #Filtering
  sizeY=length(y)
  for(i in 1:sizeY){
    a_t=a[,i]
    P_t=P[[i]]
    v_t=y[i]-exp(a_t[1])*(x[i]+exp(a_t[2]))
    dz_t=t(as.matrix(c(exp(a_t[1])*(x[i]+exp(a_t[2])),exp(a_t[1])*exp(a_t[2]),0)))
    F_t=dz_t%*%P_t%*%t(dz_t)
    a_tt=a_t+P_t%*%t(dz_t)%*%(F_t^-1)%*%v_t
    P_tt=P_t-P_t%*%t(dz_t)%*%(F_t^-1)%*%dz_t%*%P_t
    a_t1=Tmat%*%a_tt
    P_t1=Tmat%*%P_tt%*%t(Tmat)+Q
    a=cbind(a,a_t1)
    P[[i+1]]=P_t1
    dz[[i]]=dz_t
    Fmat[[i]]=F_t
    v[[i]]=v_t
  }
  
  #Smoothing
  r=matrix(0,3,1)
  alpha=matrix(0,3,sizeY)
  for(i in rev(1:sizeY)){
    K_t=Tmat%*%P[[i]]%*%t(dz[[i]])%*%( Fmat[[i]]^-1)
    L_t=Tmat-K_t%*%dz[[i]]
    r_tm=t(dz[[i]])%*%( Fmat[[i]]^-1)%*% v[[i]]+t(L_t)%*%r
    alpha[,i]=a[,i]+P[[i]]%*%r_tm
    r=r_tm
  }
  
  logLikelihood=sum(dpois(y,exp(alpha[1,])*(x+exp(alpha[2,])), log=TRUE))
  for(i in 1:sizeY){
    logLikelihood=logLikelihood+if(i>0){dmvnorm(alpha[,i], mean=a[,i], sigma=P[[i]], log=TRUE)}else{dmvnorm(alpha[,i], mean=Tmat%*%alpha[,i-1], sigma=Q, log=TRUE)}
  }
  return(list(alpha=alpha,a=a,P=P,logL=logLikelihood))
}

UKF_FS=function(y,x, par1,par2,par3, m=3){
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
    # revisar ?? 
    P_t=P[[i]]
    print(P_t)
    P_aux = chol(P_t)
    sigmas = list()
    x_t0 = a_t
    sigmas[[1]] = x_t0
    #arreglar x, Z y ybarra son univariados
    sample_points = c(-0.7, -0.5, 0.5, 0.7) 
    phis = dnorm(sample_points, mean = 0, sd = 1, log = FALSE)
    aux = dot(phis, sample_points^4)
    aux_2 = dot(phis, sample_points^2)
    w_0 = 1 - ((sum(phis)*aux)/(3*aux_2)^2)
    lambda_square = sum(phis)/((1 - w_0)*aux_2)
    w = list()
    w[[1]] = w_0
    m = 4/2
    for(j in 1:3){
      P_ti = P_aux[[j]]
      x_ti = rep(0, length(sample_points))
      for (k in 1:(4/2)){
        print(j*k + 1)
        print((j*k + 4/2) + 1)
        w[[j*k + 1]] = (1-w_0)*phis[k]/(2*m*sum(phis))
        w[[(j*k + 4/2) + 1]] = (1-w_0)*phis[k]/(2*m*sum(phis))
        x_ti = a_t + sqrt(lambda_square)*sample_points[k]*P_ti
        x_im = a_t - sqrt(lambda_square)*sample_points[k]*P_ti
        sigmas[[j*k + 1]] = x_ti
        sigmas[[(j*k + 4/2) + 1]] = x_im
      }
      
    }
    print(sigmas)
    print(w)
    samples_x[[i]] = sigmas
    y_t_expected = c(0,0,0)
    for (j in 1:(2*m)+1){
      xti = sigmas[j]
      z_ti  =  c( exp(x_ti[1])*(x[i]+ exp(x_ti[2])), exp(x_ti[1])*exp(x_ti[2]),0)
      y_t_expected = y_t_expected + w[j]*z_ti
    }
    
    P_alphavt = 0
    for(j in 1:(2*m)+1){
      xti = sigmas[j]
      z_ti  =  c(exp(x_ti[1])*(x[i] + exp(x_ti[2])),exp(x_ti[1])*exp(x_ti[2]),0)
      P_alphavt = P_alphavt + w[i]*((x_ti - a_t)%*%t(z_ti - y_t_expected))
    }
    P_vvt = matrix(0, nrow = 3, ncol = 3)
    for(j in 1:(2*m)+1){
      xti = sigmas[j]
      z_ti  =  c(exp(x_ti[1])*(x[i] + exp(x_ti[2])),exp(x_ti[1])*exp(x_ti[2]),0)
      P_vvt = P_vvt + w[i]*((z_ti - y_t_expected)%*%t(z_ti - y_t_expected))
    }
    v_t = y[i] - y_t_expected
    print(a_t)
    print(P_alphavt)
    print(P_vvt)
    print(P_vvt^-1)
    print(v_t)
    a_tt=a_t+P_alphavt%*%(P_vvt^-1)%*%v_t
    a_tts[[i]] = a_tt
    P_tt=P_t-P_alphavt%*%t(P_vvt^-1)%*%t(P_alphavt)
    
    sigmas_n = list()
    x_t0 = a_tt
    sigmas_n[[1]] = x_t0
   
    a_t1=0
    for(j in 1:(2*m)+1){
      xti = sigmas_n[j]
      a_t1 = a_t1 + w[j]*(Tmat%*%x_ti)
    }
    P_t1 = 0
    for(j in 1:(2*m)+1){
      xti = sigmas_n[j]
      P_t1 =  P_t1 + w[j]*((Tmat%*%x_ti - a_t1)%*%t(Tmat%*%x_ti - a_t1)) + w[j]*Q
    }
    a=cbind(a,a_t1)
    P[[i+1]]=P_t1
    v[[i]]=v_t
  }
  
  #smoothing
  alpha=matrix(0,3,sizeY)
  for(i in rev(2:sizeY)){
    C_t1 = 0
    sigmas = samples_x[i]
    for(j in 1:(2*m)+1){
      xti = sigmas[j]
      C_t1 = C_t1 + w[i]*((x_ti - a_t)%*%t(Tmat%*%x_ti - a[i]))
    }
    alpha[,i - 1] = a_tts[[i - 1]] + C_t1%*%(P[[i]]^-1)%*%(alpha[,i] - a[i])
  }
  
  logLikelihood=sum(dpois(y,exp(alpha[1,])*(x+exp(alpha[2,])), log=TRUE))
  for(i in 1:sizeY){
    logLikelihood=logLikelihood+if(i>0){dmvnorm(alpha[,i], mean=a[,i], sigma=P[[i]], log=TRUE)}else{dmvnorm(alpha[,i], mean=Tmat%*%alpha[,i-1], sigma=Q, log=TRUE)}
  }
  return(list(alpha=alpha,a=a,P=P,logL=logLikelihood))
}


inicializadorLogLik=function(y,x, par1,par2,par3){
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
  parm=inicializadorLogLik(y,x, par1,par2,par3)
  UKF=UKF_FS(y,x, parm$par1,parm$par2,parm$par3)
  return(UKF)
}




betaStateSpace=function(base, initialCases, distribucion_inc, distribucion_inf, lags=7, lags2=5, par1_inc=3.169434, par2_inc=5.163921, par1_inf=24.206087, par2_inf=2.984198, CI=0.95){
  omega=weightConstructionIncubacionInfeccion(lags, par1_inf, par2_inf, distribucion_inf)
  f=weightConstructionIncubacion(lags2, par1_inc, par2_inc, distribucion_inf)
  
  priorMean=3/sum(omega)
  
  cases=base
  initialVector=initialCases
  
  #Dependiente
  casesMatrix_T=matrix(0,length(cases)-lags2+1,lags2)
  for(i in 1:dim(casesMatrix_T)[1]){
    casesMatrix_T[i,]=cases[i:(i+lags2-1)]
  }
  incubPotential_T=casesMatrix_T%*%f
  
  
  #Independiente
  longCases=c(initialVector,cases)
  casesMatrix_T=matrix(0,length(longCases)-lags,lags)
  for(i in 1:dim(casesMatrix_T)[1]){
    casesMatrix_T[i,]=longCases[i:(i+lags-1)]
  }
  infectPotential_T=casesMatrix_T%*%omega
  infectPotential_T=infectPotential_T[1:length(incubPotential_T)]
  
  cases=incubPotential_T
  y=cases
  x=infectPotential_T
  
  
  par1=c(0,0,0)
  par2=10^2*c(1,1,1)
  par3=c(1,1,1)
  UKF=UKF_Complete(y,x, par1,par2,par3)
  P=UKF$P
  a=UKF$a
  alpha=UKF$alpha
  
  alphaLB=matrix(0,0,3)
  alphaUB=matrix(0,0,3)
  for(i in 1:dim(alpha)[2]){
    t=rmvnorm(10000, mean=a[,i], sigma=P[[i]])
    a1=quantile(t[,1],c((1-CI)/2,1-(1-CI)/2))  
    a2=quantile(t[,2],c((1-CI)/2,1-(1-CI)/2))
    a3=quantile(t[,3],c((1-CI)/2,1-(1-CI)/2))
    alphaLB=rbind(alphaLB,c(a1[1],a2[1],a3[1]))
    alphaUB=rbind(alphaUB,c(a1[2],a2[2],a3[2]))
  }
  
  betaLB=exp(alphaLB[,1])
  beta=exp(alpha[1,])
  betaUB=exp(alphaUB[,1])
  migrationLB=exp(alphaLB[,2])
  migration=exp(alpha[2,])
  migrationUB=exp(alphaUB[,2])
  expectedCases=exp(alpha[1,])*(x+exp(alpha[2,]))
  expectedCasesLB=qpois((1-CI)/2,exp(alpha[1,])*(x+exp(alpha[2,])))
  expectedCasesUB=qpois((1-(1-CI)/2),exp(alpha[1,])*(x+exp(alpha[2,])))
  
  data=as.data.frame(cbind(betaLB,beta,betaUB,migrationLB,migration,migrationUB,expectedCasesLB,expectedCases,expectedCasesUB,y))
  
  data$R=data$beta*sum(omega)
  data$Rlb=data$lb*sum(omega)
  data$Rub=data$ub*sum(omega)
  
  
  #R_Caso
  lb=data$betaLB
  ub=data$betaUB
  beta=data$beta
  longCaseslb=c(rep(lb[1], length=lags-1),lb)
  longCasesub=c(rep(ub[1], length=lags-1),ub)
  longCasesbeta=c(rep(beta[1], length=lags-1),beta)
  
  data$R_c=data$beta*sum(omega)
  #data$R_clb=data$lb*sum(omega)
  #data$R_cub=data$ub*sum(omega)    
  
  #for(i in 1:dim(betaT)[1]){
  #  data$R_c[i]=sum(longCasesbeta[i:(i+lags-1)]*omega)
  #  data$R_clb[i]=sum(longCaseslb[i:(i+lags-1)]*omega)      
  #  data$R_cub[i]=sum(longCasesub[i:(i+lags-1)]*omega)      
  
  #}
  
  
  return(list(data=data))
}

