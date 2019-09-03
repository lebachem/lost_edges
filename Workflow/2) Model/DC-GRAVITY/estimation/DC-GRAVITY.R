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

# load the full data
load("Data/full_data.RData")


# Estimate z for all networks ----

Z_hat<-c()

for (len in 178:N_t) {


Y<-A%*%X_inform[,len]

# Define the fitnes variables 
fit<-fitness(Y)

# real density
dens<-sum(X_inform[,len]>0)/N

# Run a genetic algorithm in order to find z:
z_hat<-my_ga(tol_accept = 0.0001)
Z_hat<-c(Z_hat,z_hat)

# save the "z"
save(z_hat,file=paste("2) Model/DC-GRAVITY/results/z_",t[len],".RData",sep=""))
print(len)
}

