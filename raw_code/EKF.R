library(readxl)
library(expm)
library(mvtnorm)
library(nloptr)
library(reshape2)
library(np)

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
inicializadorLogLik=function(y,x, par1,par2,par3){
  pars=c(par1,par2,par3)
  obj = function(par) {
    par1=par[1:3]
    par2=exp(par[4:6])
    par3=par[7:9]
    EKF=EKF_FS(y,x, par1,par2,par3)
    objectiveFunction=-EKF$logL        
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
EKF_Complete=function(y,x, par1,par2,par3){
  parm=inicializadorLogLik(y,x, par1,par2,par3)
  EKF=EKF_FS(y,x, parm$par1,parm$par2,parm$par3)
  # Format response
  EKF$Q = parm$par3
  EKF$a = lapply(seq_len(ncol(EKF$a)), function(i) EKF$a[,i]) # Matrix to list
  EKF$alpha = lapply(seq_len(ncol(EKF$alpha)), function(i) EKF$alpha[,i]) # Matrix to list
  return(EKF)
}
betaStateSpace=function(base, initialCases, distribucion_inc, distribucion_inf, lags=7, lags2=5, par1_inc=3.169434, par2_inc=5.163921, par1_inf=24.206087, par2_inf=2.984198, CI=0.95){
  print('1')
  omega=weightConstructionIncubacionInfeccion(lags, par1_inf, par2_inf, distribucion_inf)
  f=weightConstructionIncubacion(lags2, par1_inc, par2_inc, distribucion_inf)
  
  priorMean=3/sum(omega)
  
  cases=base
  initialVector=initialCases
  print(length(cases))
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
  par2=10^6*c(1,1,1)
  par3=c(1,1,1)
  EKF=EKF_Complete(y,x, par1,par2,par3)
  P=EKF$P
  a=EKF$a
  alpha=EKF$alpha
  
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
  #data$Rlb=data$lb*sum(omega)
  #data$Rub=data$ub*sum(omega)
  
  print('2')
  
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
  #why not betaT?
  #for(i in 1:dim(beta)[1]){
   # data$R_c[i]=sum(longCasesbeta[i:(i+lags-1)]*omega)
    #data$R_clb[i]=sum(longCaseslb[i:(i+lags-1)]*omega)      
    #data$R_cub[i]=sum(longCasesub[i:(i+lags-1)]*omega)      
    
  #}
  
  
  return(list(data=data))
}

