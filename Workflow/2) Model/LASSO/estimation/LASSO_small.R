#-----------------------------------------------------------#
# File that implements the non-negative LASSO               #
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
library(glmnet)
library(stargazer)

# load the data
load("Data/reduced_data.RData")

Y<-t(Y)

X_hat_lasso<-c()
for (tau in 1:N_t){
  print(tau)
  real_dens<-mean(X_inform[,tau]>0)
  
  root<-function(l){return(sum(glmnet(x=A,y=Y[tau,],lower.limits=0,lambda=l,standardize = FALSE)$beta>0)-N*real_dens)}
  
  # check the sign
  b=sign(root(0)*root(10))
  
  
  if (b==-1){
    l=uniroot(root,c(0,10))$root
    
  }
  
  if (b==1){

    red=100
    while(b==1){
      red=red*0.9
      b=sign(root(0)*root(red))
      print("find threshold value")
      
    print(b)
      
    }
    
    l=uniroot(root,c(0,red))$root
    
    
  }
  



  
  
  X_hat_lasso<-cbind(X_hat_lasso,glmnet(x=A,y=Y[tau,],lower.limits=0,lambda=l,standardize = FALSE)$beta)
}



X_hat_lasso<-as.matrix(X_hat_lasso)
res_lasso<-rowbuilder(X_hat_lasso,X_inform)
stargazer(res_lasso$mean_row)


### Save the results ----
rm(list=ls()[-which(ls()=="X_hat_lasso")])
save.image("2) Model/LASSO/results/lasso_small.RData")
