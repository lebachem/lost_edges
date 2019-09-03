#-----------------------------------------------------------#
# File that implements the IPFP Model                       #
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
library(ipfp)
library(stargazer)

# load the data
load("Data/full_data.RData")


## Fit a model----
X_hat_entrop<-c()
for (tau in 1:N_t){
  X_hat_entrop<-rbind(X_hat_entrop,ipfp(y=Y[tau,],A=A,x0=rep(1,n*(n-1))))
  print(tau)
}

X_hat_ipfp<-t(X_hat_entrop)

res_ipfp<-rowbuilder(X_hat_ipfp,X_inform)

stargazer(res_ipfp$mean_row)

## Save the results----
rm(list=ls()[-which(ls()=="X_hat_ipfp")])
save.image("2) Model/IPFP/results/ipfp.RData")
