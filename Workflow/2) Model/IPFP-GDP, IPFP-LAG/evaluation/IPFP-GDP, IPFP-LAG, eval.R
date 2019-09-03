#------------------------------------------------------------------#
# File that evaluates the approach by Lebacher and Kauermann 2019  #
#------------------------------------------------------------------#

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


## load data with GDP ----
X_hat_LK  <-c()
for (tau in 1:N_t){
  load(paste("2) Model/IPFP-GDP, IPFP-LAG/results/res_",tau,".RData",sep=""))
X_hat_LK <-cbind(X_hat_LK,vec_E) 
print(tau)
}


res_lk<-rowbuilder(X_hat_LK,X_inform[,1:dim(X_hat_LK)[2]])
stargazer(res_lk$mean_row)


## load data with LAG ----
X_hat_LK2  <-c()
for (tau in 2:N_t){
  load(paste("2) Model/IPFP-GDP, IPFP-LAG/results/res_lag_",tau,".RData",sep=""))
  X_hat_LK2 <-cbind(X_hat_LK2,vec_E) 
  print(tau)
}

X_hat_LK2<-cbind(ipfp(y=t(as.matrix(Y[,1])),A=A,x0=rep(1,n*(n-1))),X_hat_LK2)
res_lk2<-rowbuilder(X_hat_LK2,X_inform[,1:dim(X_hat_LK2)[2]])
stargazer(res_lk2$mean_row)

