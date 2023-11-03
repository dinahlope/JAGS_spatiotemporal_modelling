#' JAGS model specification for spatiotemporal epidemiological modelling
#' @author Dinah Lope/ Haydar Demirhan
#' @version 2.0
#' @created 12/10/2023
#' @license GPLv3
#' 

options(scipen = 999)
graphics.off()
rm(list=ls())
library(rjags)
library(runjags)
library(coda)
library(nimble)
library(magrittr)
library(tidyverse)
library(reshape2)
load.module("GeoJAGS")

source("DBDA2E-utilities.R")
count <- 0
tmp1 <- list()
tmp2 <- list()
tmp3 <- list()
set.seed(2021)
# Create a random adjacency matrix
num_area <- 10
num_neighbours <- 5
adj <- matrix(0, nrow = num_area, ncol = num_area)
upperIndicies <- which(upper.tri(adj, diag = FALSE), arr.ind = TRUE)
indices <- (as.matrix(upperIndicies[sample(1:nrow(upperIndicies), num_neighbours , replace = F),]))
adj[indices] <- 1
adj[(indices[,c(2,1)])] <- 1

for (num_records in c(100, 500, 1000, 2000, 4000, 8000, 10000)){
  for (numSavedSteps in c(100, 500, 1000, 2000, 4000, 5000, 10000)){
    cat("number of records = ",num_records,"\n")
    cat("number of saved steps = ",numSavedSteps,"\n" )
    
    count <- count + 1
    
    areas <- 1:10
    num_time <- 1:(num_records/(areas %>% length()))
    
    s <- areas %>% length()
    t <- num_time %>% length() 
    if((s*t == num_records) %in% TRUE){"ok proceed"}else{"wrong matrix count"}

    zdata <- rpois(num_records, 10)
    covariate1 <-  runif(num_records)
    
    region_index <- NULL
    for(i in 1:10){
      region_index[i] <- areas[i]
    }
    
    if((region_index %>% length() == s) %in% TRUE){"ok proceed"}else{"wrong index count"}
    
    data_all <- tibble(cov = covariate1, 
                       flu = zdata, 
                       datetime = data.frame(tmp = rep(num_time, s)) %>% arrange(tmp) %>% pull(tmp),
                       lga = rep(areas, t)) %>% as.data.frame()
    
    flu_constants = list(x_var = acast(data_all, lga~datetime, value.var="cov"),
                         # index = region_index,
                         z =  acast(data_all, lga~datetime, value.var="flu"),
                         space = s,
                         time = t,
                         mu = array(0,num_records),
                         adj = adj)
                         # n = num_records)
    
    # Specify the model
    modelString = "
      model {
        for(s in 1:space){
          for(t in 1:time){
            lambda[s,t] <- exp(a[1] + x_var[s,t]*a[2] + theta[s] + phi[s])
            z[s,t] ~ dpois(lambda[s,t])
          }
        }
        for(j in 1:space){
          theta[j] ~ dnorm(0, 0.10)
        }
        for(i in 1:2){
          a[i] ~ dnorm(0, 0.10)
        }
        phi[1:space] ~ dmnorm(mu[1:space], precMatrixCAR(adj[1:space, 1:space], 0.999, 10))
      }
    " 
    writeLines(modelString , con="TEMPmodel2D.txt")
    
    # Generate chains
    parameters <- c("a[1]", "a[2]")
    for (s in 1:s){
      for(t in 1:t){
        parameters <- c(parameters , 
                        paste0("lambda[",s,",", t,"]"))
      }
    }
    for (i in 1:s){
      parameters <- c(parameters, 
                      paste0("theta[",i,"]"), paste0("phi[",i,"]"))
    }
    # parameters
    
    options(scipen = 999)
    adaptSteps = 500
    burnInSteps = 100
    nChains = 1
    thinSteps = 1
    nIter = ceiling((numSavedSteps*thinSteps)/nChains)

    startTime = proc.time()
    jmodel <- jags.model("TEMPmodel2D.txt",
                         data=flu_constants,
                         n.chains=nChains,
                         # inits=inits,
                         n.adapt=adaptSteps)

    update(jmodel, n.iter = burnInSteps)
    samples <- coda.samples(jmodel, variable.names=parameters, n.iter=nIter, thin = thinSteps)
    stopTime = proc.time()
    tmp_b = stopTime - startTime
    tmp_a <- paste(numSavedSteps, "; numSavedSteps; ", nIter, "; nIter ;", num_records, "; num_records ;")
    tmp3 <- rbind(tmp3,data.frame(tmp_a, tmp_b[1], tmp_b[2], tmp_b[3]))
  }
}

write.csv(tmp3, file = "data_100Iter_2D_HD_tmp3_GEOJags.csv")
