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
  for(t in 1:time){
  for(s in 1:space){
    
    pi[t,s]  <- ifelse(m == 1, (ifelse((s == 10 || s == 14 || s == 26 || s == 45 || s == 52 || s == 74 || s == 76 || s == 33),

                                          (guess*(1/2) + (1.0-guess)* ilogit(b[1] + tested[t,s]*b[2] +
                                                    mintemp[t,s]*b[3] + wind_N[t,s]*b[4] + b[5]*daycount[t,s] + b[6]*daycount2[t,s] + b[7]*daycount3[t,s])),

                                   (guess*(1/2) + (1.0-guess)* ilogit(b[1] + tested[t,s]*b[2] + mintemp[t,s]*b[3] + wind_N[t,s]*b[4])))),
                        (guess*(1/2)  + (1.0-guess)*ilogit(b[1] + tested[t,s]*b[2] + mintemp[t,s]*b[3] + wind_N[t,s]*b[4])))
                     
    log.lambda[t,s] <- ifelse(m == 1, log(pop[t,s]+ a[6] + vac2[t,s]*a[1] + vacn[t,s]*a[2] + theta[s] + phi[s]),
                              (ifelse((s == 10 || s == 14 || s == 26 || s == 45 || s == 52 || s == 74 || s == 76 || s == 33), 
                                           
                                           log(pop[t,s]+ a[6] + vac2[t,s]*a[1] + vacn[t,s]*a[2] + 
                                              theta[s] + phi[s] + a[4]*daycount[t,s] + a[3]*daycount2[t,s] + a[5]*daycount3[t,s]),
                                       
                                           log(pop[t,s]+ a[6] + vac2[t,s]*a[1] + vacn[t,s]*a[2] + theta[s] + phi[s]))))
    
    z[t,s] ~ dpois(pi[t,s]*exp(log.lambda[t,s]))  
    
  }
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
