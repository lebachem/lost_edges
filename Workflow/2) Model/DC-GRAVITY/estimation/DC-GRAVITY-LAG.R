#-----------------------------------------------------------#
# File that implements the Model by Cimini et al            #
#-----------------------------------------------------------#


# Settings ----
# remove the old stuff
rm(list=ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
setwd("../../..")

# set seed
set.seed(42)

# Load the packges
library(GA)
library(PRROC)
library(networkTomography)

# load the data
load("Data/reduced_Data.RData")

# Estimate z for all networks ----
Z_hat<-c()

for (len in 2:N_t) {
  
  
  # Define the fitnes variables 
  fit<-log(1.1+X_inform[,len-1])

  dens<-sum(X_inform[,len]>0)/N
  
  # Run a genetic algorithm in order to find z:
  z_hat<-my_ga(tol_accept = 0.0001,lower=0,upper=10,popsize=100,maxiter=100,fac_pop=1.01,fac_iter=1.01,fac_upper=0.8)

  Z_hat<-c(Z_hat,z_hat)
  
  # save the "z"
  save(z_hat,file=paste("2) Model/DC-GRAVITY/results/z_small_cov_lag_",t[len],".RData",sep=""))
  print(len)
}

