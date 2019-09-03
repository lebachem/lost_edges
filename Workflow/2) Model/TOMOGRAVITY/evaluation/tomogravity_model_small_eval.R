#-----------------------------------------------------------#
# File that evaluates the TOMOGRAVITY Model                 #
#-----------------------------------------------------------#

# Settings ----
# remove the old stuff
rm(list=ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
setwd("../../..")

# Load the packges
library(ipfp)
library(stargazer)

# load the data
load("Data/reduced_Data.RData")

## load data----
X_hat_tomo  <-c()
for (tau in t){
  load(paste("2) Model/TOMOGRAVITY/results/fit_",tau,".RData",sep=""))
  X_hat_tomo <-cbind(X_hat_tomo,t(fit)) 
  print(tau)
}


res_tomo<-rowbuilder(X_hat_tomo,X_inform[,1:dim(X_hat_tomo)[2]])

stargazer(res_tomo$mean_row,type="text")



