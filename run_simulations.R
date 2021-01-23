prueba = simulate_with_immigrants('weibull', 2.8 ,3.5 ,5 , 120, rep(1.1,5), rep(1.1,120))
incidences = prueba$I[8:length(prueba$I)]
inital_cases = prueba$I[1:7]
result = betaStateSpace(incidences, inital_cases, 'weibull', 'weibull')
