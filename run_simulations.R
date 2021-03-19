prueba = simulate_with_immigrants('weibull', 2.8 ,3.5 ,5 , 120, rep(1.1,5), rep(1.1,120))
incidences = prueba$I[8:length(prueba$I)]
inital_cases = prueba$I[1:7]
result = betaStateSpace(incidences, inital_cases, 'weibull', 'weibull')

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
  "
  1.000000    0.000000    0.654472    5.088467   14.530000   24.271440
 [7]   33.229114   41.794700   48.577486   53.906599   58.368293   62.996933
[13]   68.637997   74.420266   80.976262   90.128817  100.714679  111.131145
[19]  121.243183  128.085423  132.039996  140.131932  152.895572  165.415540
[25]  176.705081  190.549712  203.589459  219.146081  241.073498  257.609232
[31]  268.045601  282.099138  298.450864  319.025371  350.587436  387.807790
[37]  417.562571  438.360742  466.911986  512.049258  566.590726  618.243249
[43]  660.467622  703.744792  751.133435  797.094681  851.344023  917.037836
[49]  978.978582 1036.318638
  " 
  par1=c(0,0,0)
  par2= c(0.1,0.1,0.1) #10^6*c(1,1,1)
  par3=c(1,1,1)
  UKF=UKF_FS(y,x, par1,par2,par3)
  #UKF=UKF_Complete(y,x, par1,par2,par3)
  P=UKF$P
  a=UKF$a
  print(x*a[1,] + (a[2,]))
  print(a[1,]*sum(omega))
  alpha=UKF$alpha
  print(x)
  
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
  print(data$beta*sum(omega))
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

