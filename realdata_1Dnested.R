#'' 
#' Spatiotemporal modelling
#' Underreported COVID-19 cases
#' 
#' @author Dinah Lope/ Haydar Demirhan
#' @created 17/09/2021
#''

# model specification
modelString = "
  model {
  for(i in 1:N){
  
      pi[i]  <- ifelse(m == 1, (ifelse((i <= 2409 && i >=  2337)|| (i <= 1022 && i >= 950) || (i <= 5329 && i >= 5257) || (i <= 1898 && i >= 1826) || (i <= 5475 && i >= 5403)
                                            || (i <= 730 && i >= 658) || (i <= 3212 && i >= 3140) || (i <= 3723 && i >= 3651),

                                          (guess*(1/2) + (1.0-guess)* ilogit(b[1] + tested[i,1]*b[2] +
                                                    mintemp[i,1]*b[3] + wind_N[i,1]*b[4] + b[5]*daycount[i,1] + b[6]*daycount2[i,1] + b[7]*daycount3[i,1])),

                                   (guess*(1/2) + (1.0-guess)* ilogit(b[1] + tested[i,1]*b[2] + mintemp[i,1]*b[3] + wind_N[i,1]*b[4])))),
                        (guess*(1/2)  + (1.0-guess)*ilogit(b[1] + tested[i,1]*b[2] + mintemp[i,1]*b[3] + wind_N[i,1]*b[4])))
                     
    log.lambda[i] <- ifelse(m == 1, log(pop[i,1]+ a[6] + vac2[i,1]*a[1] + vacn[i,1]*a[2] + theta[index[i]] + phi[index[i]]),
                              (ifelse((i <= 2409 && i >=  2337) || (i <= 1022 && i >= 950) || (i <= 5329 && i >= 5257) || (i <= 1898 && i >= 1826) || (i <= 5475 && i >= 5403)
                                            || (i <= 730 && i >= 658) || (i <= 3212 && i >= 3140) || (i <= 3723 && i >= 3651), 
                                           
                                           log(pop[i,1]+ a[6] + vac2[i,1]*a[1] + vacn[i,1]*a[2] + 
                                              theta[index[i]] + phi[index[i]] + a[4]*daycount[i,1] + a[3]*daycount2[i,1] + a[5]*daycount3[i,1]),
                                       
                                           log(pop[i,1]+ a[6] + vac2[i,1]*a[1] + vacn[i,1]*a[2] +  theta[index[i]] + phi[index[i]]))))
    
    z[i] ~ dpois(pi[i]*exp(log.lambda[i]))
  }  
  
  sigma.theta ~ dnorm(0,1)T(0,)
  for(j in 1:R){
    theta[j] ~ dnorm(0, 1/sigma.theta)
  }
    
  phi[1:R] ~ dmnorm(zrs[1:R], precMatrixCAR(adj[1:R, 1:R], rho.car, tau))
  
  a[1] ~ dnorm(0, 1/sigmaA)
  a[2] ~ dnorm(0, 1/sigmaA) 
  a[3] ~ dnorm(0, 1/sigmaA)
  a[4] ~ dnorm(0, 1/sigmaA)
  a[5] ~ dnorm(0, 1/sigmaA)
  a[6] ~ dnorm(0, 1/sigmaA) 
  sigmaA ~ dnorm(20,0.001)T(0,)
  
  b[1] ~ dnorm(0, 1/sigmaB)
  for(i in 2:7){
    b[i] ~ dnorm(0, 1/sigmaB)
  }
  
  sigmaB ~ dnorm(30,0.001)T(0,)
  nu ~ dnorm(0,1)T(0,) 
  tau <- 1/(nu*nu) 
  rho.car <- 0.999

  m ~ dcat(mPriorProb[])
  mPriorProb[1]  <- 0.5
  mPriorProb[2]  <- 0.5
  guess ~ dbeta(1,9)
}
"
# close quote for modelString
writeLines(modelString , con="TEMPmodel.txt")
