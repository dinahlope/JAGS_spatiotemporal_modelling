#'
#' JAGS model specification for spatiotemporal epidemiological modelling
#' @author Dinah Lope, Haydar Demirhan
#' @version 1.0
#' @created 11/10/2023
#' @license GPLv3
#' 
#' 

library(rjags)
library(runjags)
library(coda)
library(nimble)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tictoc)

# source("0 settings.R")

#load JAGS utility file
source("DBDA2E-utilities.R")

#load the module GEOJAGS
load.module("GeoJAGS")

# load datafile from COVID study
# source: https://github.com/dinahlope/Spatiotemporal_Bayesian_COVID19
load(file="COVID19DATA.RData")


#load all other data required
load(file = "sup_data.RData")
load(file = "region_index.RData")

covdat <- datcov %>% 
  group_by(lganame) %>%
  mutate(daycount = n():1, 
         daycount2 = daycount^2, 
         daycount3 = daycount^3) %>%
  ungroup() %>% 
  rename(pct_second_tmp = seconddose_pct,
         vac_n = vaccentres_n) %>% 
  left_join(pop_dat)

cov_N = length(covdat$covidcase)
space = length(unique(covdat$lganame))
n_regions = space

TSmodel <- array(0,cov_N) 
TSmodel[2337:2409] <- 1
TSmodel[950:1022] <- 1
TSmodel[5257:5329] <- 1
TSmodel[1826:1898] <- 1
TSmodel[5403:5475] <- 1
TSmodel[658:730] <- 1
TSmodel[3140:3212] <- 1
TSmodel[3651:3723] <- 1

cov_constants = list(
  n = cov_N,
  n_adj = n_adj,
  adj = wm,

  vac2 = as.matrix(scale(covdat$pct_second_tmp)),
  vacn = as.matrix(scale(covdat$vac_n)),
  daycount3  = as.matrix(covdat$daycount3),
  daycount2  = as.matrix(covdat$daycount2),
  daycount  = as.matrix(covdat$daycount),
  
  tested = as.matrix(covdat$tested),
  mintemp = as.matrix(covdat$min_temp), 
  wind_N = as.matrix(covdat$wind_direc_N),
  pop = as.matrix(covdat$pop),
  
  R=n_regions,
  index=region_index,
  w=weights,
  l_adj=length(adjacency),
  
  z=covdat$covidcase,
  zrs = array(0,n_regions),
  TSmodel = TSmodel
)

# Specify the model
modelString = "
  model {
  for(i in 1:n){
  
      pi[i]  <- ifelse(m == 1, (ifelse(TSmodel[i] == 1, ( guess*(1/2) + (1.0-guess)* ilogit(b[1] + tested[i,1]*b[2] +
                                                    mintemp[i,1]*b[3] + wind_N[i,1]*b[4] + b[5]*daycount[i,1] + b[6]*daycount2[i,1] + b[7]*daycount3[i,1]) ),

                                   ( guess*(1/2) + (1.0-guess)* ilogit(b[1] + tested[i,1]*b[2] + mintemp[i,1]*b[3] + wind_N[i,1]*b[4]) ) )  ),
                        ( guess*(1/2)  + (1.0-guess)*ilogit(b[1] + tested[i,1]*b[2] + mintemp[i,1]*b[3] + wind_N[i,1]*b[4])) )
                     
    log.lambda[i] <- ifelse(m == 1, log(pop[i,1]+ a[6] + vac2[i,1]*a[1] + vacn[i,1]*a[2] + theta[index[i]] + phi[index[i]]),
                              (ifelse(TSmodel[i] == 1, log(pop[i,1]+ a[6] + vac2[i,1]*a[1] + vacn[i,1]*a[2] + 
                                              theta[index[i]] + phi[index[i]] + a[4]*daycount[i,1] + a[3]*daycount2[i,1] + a[5]*daycount3[i,1]),
                                           log(pop[i,1]+ a[6] + vac2[i,1]*a[1] + vacn[i,1]*a[2] +  theta[index[i]] + phi[index[i]]))))
    
    z[i] ~ dpois(pi[i]*exp(log.lambda[i]))
  }  
  
  sigma.theta ~ dnorm(0,1)T(0,)
  for(j in 1:R){
  theta[j] ~ dnorm(0, 1/sigma.theta) 
  }
  
  phi[1:R] ~ dmnorm(zrs[1:R], precMatrixCAR(adj[1:R, 1:R], rho.car, tau))
  
  for(i in 1:6){
  a[i] ~ dnorm(0, 1/sigmaA)
  }
  sigmaA ~ dnorm(20,0.001)T(0,) 
  
  for(i in 1:7){
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

#initialize chains
aInits = rnorm(6, 0, 1)
bInits = rnorm(7, 0, 1)

aInits[1] = 2.2
aInits[2] = 2.2
aInits[3] = 2.2
aInits[4] = 1.3
aInits[5] = 2.3
aInits[6] = 2.3

bInits[1] = -1.6
bInits[2] = -1.6
bInits[3] = -1.5
bInits[4] = 0.3
bInits[5] = -1.6
bInits[6] = -1.5
bInits[7] = 0.3

cov_inits=list(list(sigma.theta = 1.25, a = aInits, b = bInits,
                    phi = rep(0.55,n_regions),
                    theta =  rep(0.01,n_regions),
                    nu = 20, sigmaA = 20, sigmaB = 20, guess = 0.5,
                    m = 1,.RNG.name = "base::Mersenne-Twister",.RNG.seed = 22021),
               list(sigma.theta = 0.25,
                    phi =  rep(0.5,n_regions),
                    theta = rep(0.015,n_regions),
                    nu = 21,
                    m = 2,.RNG.name = "base::Mersenne-Twister",.RNG.seed = 32019),
               list(sigma.theta = 0.225,
                    phi =  rep(0.56,n_regions),
                    theta = rep(0.0125,n_regions),
                    nu = 19,
                    m = 1,.RNG.name = "base::Mersenne-Twister",.RNG.seed = 93019),
               list(sigma.theta = 0.325,
                    phi =  rep(0.575,n_regions),
                    theta = rep(0.012,n_regions),
                    nu = 18,
                    m = 2,.RNG.name = "base::Mersenne-Twister",.RNG.seed = 13209))

#generate chains
parameters <- c("a[1]", "a[2]", "a[3]", "a[4]", "a[5]", "a[6]",
                "b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "b[6]", "b[7]",
                "nu", "sigma.theta", "m", "sigmaA", "sigmaB", "guess")

for (i in 1:nrow(covdat)){
  parameters <- c(parameters , 
                  paste0("pi[",i,"]"),
                  paste0("log.lambda[",i,"]"))
}

for (r in 1:n_regions){
  parameters <- c(parameters, 
                  paste0("phi[",r,"]"),
                  paste0("theta[",r,"]"))
}

adaptSteps = 100
burnInSteps = 5000
nChains = 4
thinSteps = 11
numSavedSteps = 43000

# adaptSteps = 10
# burnInSteps = 10
# nChains = 4
# thinSteps = 1
# numSavedSteps = 10
nIter = ceiling((numSavedSteps*thinSteps)/nChains)

set.seed(2021)
tic()
jmodel <- jags.model("TEMPmodel.txt",
                     data=cov_constants,
                     n.chains=nChains,
                     inits=cov_inits ,
                     n.adapt=adaptSteps)

update(jmodel, n.iter = burnInSteps)
samples <- coda.samples(jmodel, variable.names=parameters, n.iter=nIter, thin = thinSteps)
toc()
