#'
#' JAGS model specification for spatiotemporal epidemiological modelling
#' @author Dinah Lope, Haydar Demirhan
#' @version 1.0
#' @created 11/10/2023
#' @license GPLv3    
#' 

options(scipen = 999)
graphics.off()
rm(list=ls())
library(rjags)
library(runjags)
library(coda)
library(dplyr)
library(reshape2)
load.module("GeoJAGS")

# General functions by Kruschke, J., 2014. Doing Bayesian data analysis: A tutorial with R, JAGS.
source("DBDA2E-utilities.R")
set.seed(2021)

indexing <- "2D" # or "1D"
num_records <- 6 
cat("number of records = ", num_records, "\n")
num_area <- 2

areas <- 1:num_area
num_time <- 1:(num_records/(areas %>% length()))

s <- areas %>% length() 
t <- num_time %>% length() 
# Check
if((s*t == num_records) %in% TRUE){"ok proceed"}else{"wrong matrix count"}

# Create a random adjacency matrix
num_neighbours <- 1
adj <- matrix(0, nrow = num_area, ncol = num_area)
upperIndicies <- which(upper.tri(adj, diag = FALSE), arr.ind = TRUE)
indices <- t(as.matrix(upperIndicies[sample(1:nrow(upperIndicies), num_neighbours , replace = F),]))
adj[indices] <- 1
adj[t(indices[,c(2,1)])] <- 1

# Generate random data arbitrarily
zdata <- rpois(num_records, num_area)
covariate1 <-  runif(num_records)

# Create inputs for JAGS
if (indexing == "1D"){
  region_index <- NULL
  for(i in 1:num_area){
    region_index[i] <- areas[i]
  }
  region_index <- rep(region_index, num_records/num_area)
  
  # Check
  if((region_index %>% length() == s) %in% TRUE){"ok proceed"}else{"wrong index count"}
  
  input_data = list(x_var = as.matrix(covariate1),
                    index = region_index,
                    z = zdata,
                    space = s,
                    mu = array(0,num_records),
                    adj = adj,
                    n = num_records)
  
  # Specify the model with 1D
  modelString = "
        model {
          for(i in 1:n){
            lambda[i] <- exp(a[1] + x_var[i,1]*a[2] + theta[index[i]] + phi[index[i]])
            z[i] ~ dpois(lambda[i])
          }
          for(j in 1:space){
            theta[j] ~ dlogis(0, 0.1) #dnorm(0, 0.1)
          }
          for(i in 1:2){
            a[i] ~ dlogis(0, 0.1) #dnorm(0, 0.1)
          }
          phi[1:space] ~ dmnorm(mu[1:space], precMatrixCAR(adj[1:space, 1:space], 0.999, 1))
        }
      " 
  writeLines(modelString , con="TEMPmodel.txt")
  # Set the parameters
  parameters <- c("a[1]", "a[2]")
  for (i in 1:num_records){
    parameters <- c(parameters , 
                    paste0("lambda[",i,"]"))
  }
} else if (indexing == "2D"){
  data_all <- tibble(cov = covariate1, 
                     flu = zdata, 
                     datetime = data.frame(tmp = rep(num_time, s)) %>% arrange(tmp) %>% pull(tmp),
                     lga = rep(areas, t)) %>% as.data.frame()
  
  input_data = list(x_var = acast(data_all, lga~datetime, value.var="cov"),
                    z =  acast(data_all, lga~datetime, value.var="flu"),
                    space = s,
                    time = t,
                    mu = array(0,num_records),
                    adj = adj)
  
  # Specify the model with 2D
  modelString = "
      model {
        for(s in 1:space){
          for(t in 1:time){
            lambda[s,t] <- exp(a[1] + x_var[s,t]*a[2] + theta[s] + phi[s])
            z[s,t] ~ dpois(lambda[s,t])
          }
        }
        for(j in 1:space){
          theta[j] ~ dt(0,0.1,30)#dlogis(0, 0.1) #dnorm(0, 0.1)
        }
        for(i in 1:2){
          a[i] ~ dt(0,0.1,30)#dlogis(0, 0.1) #dnorm(0, 0.1)
        }
        phi[1:space] ~ dmnorm(mu[1:space], precMatrixCAR(adj[1:space, 1:space], 0.999, 1))
      }
      " 
  writeLines(modelString , con="TEMPmodel.txt")
  # Set the parameters
  parameters <- c("a[1]", "a[2]")
  for (s in 1:s){
    for(t in 1:t){
      parameters <- c(parameters, 
                      paste0("lambda[",s,",", t,"]"))
    }
  }
}

for (i in 1:s){
  parameters <- c(parameters, 
                  paste0("theta[",i,"]"), paste0("phi[",i,"]"))
}

# Set initials if needed
# flu_inits=list(
#   a = rnorm(2, 0, 1),
#   theta = rnorm(s, 0, 0.25))

# MCMC settings
adaptSteps = 5000
burnInSteps = 100000
nChains = 3
thinSteps = 251
numSavedSteps = 10000
nIter = ceiling((numSavedSteps*thinSteps)/nChains)
cat("Number of saved steps =", numSavedSteps, "\n" )
cat("Number of iterations =", nIter, "\n" )

# Generate chains
startTime = proc.time()

jmodel <- jags.model("TEMPmodel.txt",
                     data = input_data,
                     n.chains = nChains,
                     n.adapt = adaptSteps)

update(jmodel, n.iter = burnInSteps)
samples <- coda.samples(jmodel, variable.names=parameters, n.iter=nIter, thin = thinSteps)
stopTime = proc.time()
tmp_b = stopTime - startTime
cat("Elapsed time with", indexing, "indexing =", tmp_b[3], "\n" )


# Examine chains and plot posteriors
graphics.off()
diagMCMC(samples, parName="a[1]")
diagMCMC(samples, parName="a[2]")

diagMCMC(samples, parName=paste0("theta[1]"))
diagMCMC(samples, parName=paste0("theta[2]"))

diagMCMC(samples, parName=paste0("phi[1]"))
diagMCMC(samples, parName=paste0("phi[2]"))
