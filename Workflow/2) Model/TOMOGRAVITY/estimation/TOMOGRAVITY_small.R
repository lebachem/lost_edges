#-----------------------------------------------------------#
# File that implements the Tomogravity model                #
#-----------------------------------------------------------#

# Settings ----
# remove the old stuff
rm(list=ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
setwd("../../..")

# Load the packges
library(networkTomography)

# set seed
set.seed(42)

# load the data
load("Data/reduced_data.RData")

# find the optimal threshold value for the tomogravity model----

for (tau in 1:N_t){
  
  fit<-tomogravity(t(as.matrix(Y[,tau])),A,lambda=0.01)$Xhat
  save(fit,file=paste("2) Model/TOMOGRAVITY/results/fit_",t[tau],".RData",sep=""))
  print(tau)
}
# save results----
rm(list=ls()[-which(ls()=="X_hat_tomogravity")])
save.image("2) Model/TOMOGRAVITY/results/Tomogravity.RData")


