#--------------------------------------------------------------#
# File that implements the Hiearchical Models (EmpiricalBayes) #
#--------------------------------------------------------------#

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
library(rlist)
library(nloptr)
library(pracma)
library(truncdist)
library(systemicrisk)

# load the full data
Data<-load("Data/full_data.RData")



rep<-100

# Run 100 sequences of the ER and fitness model
X_hat_fitness<-list()
X_hat_ER<-list()

Y<-Y/1000

for (tau in 1:length(t)){
  print(tau)
  
  l<-Y[tau,1:n]
  l_index<-which(l==0)
  
  a<-Y[tau,(n+1):(2*n)]
  a_index<-which(a==0)
  
  dens<-sum(as.numeric(X_inform[,tau]>0))/N
  
  
  
  L_fixed<-matrix(NA,ncol=n,nrow=n)
  
  if (length(a_index)>0&length(l_index)>0){
    L_fixed[l_index,a_index]<-0
  }
  
  if (length(a_index)>0&length(l_index)==0){
    L_fixed[,a_index]<-0
  }
  
  if (length(a_index)==0&length(l_index)>0){
    L_fixed[l_index,]<-0
  }
  diag(L_fixed)<-0
  
  ## Fitness Model----
  model1<-calibrate_FitnessEmp(l,a,L_fixed = L_fixed,targetdensity = dens)
  Lsamp1 <- sample_HierarchicalModel(l=l,a=a,L_fixed = L_fixed,
                                     model=model1,nsamples=rep,thin=1e2)
  
  ## ER Model----
  model2<-calibrate_ER(l,a,L_fixed = L_fixed,targetdensity = dens)
  Lsamp2 <- sample_HierarchicalModel(l=l,a=a,L_fixed = L_fixed,
                                     model=model2,nsamples=rep,thin=1e2)
  

  X_hat_fitness[[tau]]<-Lsamp1
  X_hat_ER[[tau]]<-Lsamp2
  save(Lsamp1,file=paste("2) Model/H-FIT, H-ER/results/fit_",t[tau],".RData",sep=""))
  save(Lsamp2,file=paste("2) Model/H-FIT, H-ER/results/ER_",t[tau],".RData",sep=""))
  
}



## save the results----
rm(list=ls()[-c(which(ls()=="meanX_hat_ER"),which(ls()=="meanX_hat_fitness"))])
save.image("2) Model/H-FIT, H-ER/results/hierarch.RData")


