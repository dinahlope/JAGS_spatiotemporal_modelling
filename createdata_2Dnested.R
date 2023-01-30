#'
#' JAGS indexing study
#' @author Dinah Lope, Haydar Demirhan
#' @version 1.0
#' @created 07/06/2022
#'

library(rjags)
library(runjags)
library(coda)
library(nimble)
library(magrittr)
library(tidyverse)
library(reshape2)

source("0 settings.R")
source("DBDA2E-utilities.R")

count <- 0
tmp1 <- list()
tmp2 <- list()
tmp3 <- list()

set.seed(2021)
for (num_records in c(100, 500, 1000, 2000, 4000, 8000, 10000)){
  for (numSavedSteps in c(100, 500, 1000, 2000, 4000, 5000, 10000)){
    
    count <- count + 1
    num_area <- 1:10
    num_time <- 1:(num_records/(num_area %>% length()))
    
    s <- num_area %>% length() 
    t <- num_time %>% length() 
    if((s*t == num_records) %in% TRUE){"ok proceed"}else{"wrong matrix count"}
    
    zdata <- rpois(num_records, 10)
    covariate1 <- runif(num_records)
    
    region_index <- NULL
    for(i in 1:10){
      region_index[i] <- num_area[i]
    }
    if((region_index %>% length() == s) %in% TRUE){"ok proceed"}else{"wrong index count"}
    
    data_all <- tibble(cov = covariate1, 
                       flu = zdata, 
                       datetime = data.frame(tmp = rep(num_time, s)) %>% arrange(tmp) %>% pull(tmp),
                       lga = rep(num_area, t)) %>% as.data.frame()
    
    flu_constants = list(x_var = acast(data_all, lga~datetime, value.var="cov"),
                         # index = region_index,
                         z =  acast(data_all, lga~datetime, value.var="flu"),
                         space = s,
                         time = t)
                         # n = num_records)
    
    # Specify the model
    modelString = "
      model {
        for(s in 1:space){
          for(t in 1:time){
            lambda[s,t] <- exp(a[1] + x_var[s,t]*a[2] + theta[s])
            z[s,t] ~ dpois(lambda[s,t])
          }
        }
        for(j in 1:space){
          theta[j] ~ dnorm(0, 1/10)
        }
        for(i in 1:2){
          a[i] ~ dnorm(0, 1/10)
        }
      }
    " # close quote for modelString
    writeLines(modelString , con="TEMPmodel.txt")
    
    # Generate chains
    parameters <- c("a[1]", "a[2]")
    for (s in 1:s){
      for(t in 1:t){
        parameters <- c(parameters , 
                        paste0("lambda[",s,",", t,"]"))
      }
    }
    for (i in 1:s){
      parameters <- c(parameters , 
                      paste0("theta[",i,"]"))
    }

    options(scipen = 999)
    adaptSteps = 500
    burnInSteps = 100
    nChains = 1
    thinSteps = 1
    # numSavedSteps = 1000
    nIter = ceiling((numSavedSteps*thinSteps)/nChains)
    # nIter
    
    set.seed(2021)
    startTime = proc.time()
    jmodel <- jags.model("TEMPmodel.txt",
                         data=flu_constants,
                         n.chains=nChains,
                         #inits=flu_inits,
                         n.adapt=adaptSteps)
    update(jmodel, n.iter = burnInSteps)
    samples <- coda.samples(jmodel, variable.names=parameters, n.iter=nIter, thin = thinSteps)
    stopTime = proc.time()
    tmp_b = stopTime - startTime
    
    tmp_a <- paste(numSavedSteps, "numSavedSteps", nIter, "nIter", num_records, "num_records")
    # tmp_b <- system.time(update(jmodel, n.iter = burnInSteps))
    # samples <- coda.samples(jmodel, variable.names=parameters, n.iter=nIter, thin = thinSteps)
    
    tmp1[[count]] <- list(tmp_a, tmp_b[1], tmp_b[2], tmp_b[3], samples)
    tmp2[[count]] <- list(tmp_a, tmp_b[1], tmp_b[2], tmp_b[3], samples)
    tmp3[[count]] <- list(tmp_a, tmp_b[1], tmp_b[2], tmp_b[3])
  }
}

save(tmp1, tmp2, file = "data_all_nested.RData")


