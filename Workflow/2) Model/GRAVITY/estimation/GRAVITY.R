#-----------------------------------------------------------#
# File that implements the Gravity Model                    #
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
library(networkTomography)
library(stargazer)

# load the data
load("Data/full_Data.RData")

## Fit the model----
# Calculate the 
X_hat_gravity<-t(gravity(Y,getSrcDstIndices(A)))

res_gravity<-rowbuilder(X_hat_gravity,X_inform)
stargazer(res_gravity$mean_row)

## Save the results----
rm(list=ls()[-which(ls()=="X_hat_gravity")])
save.image("2) Model/GRAVITY/results/gravity.RData")
