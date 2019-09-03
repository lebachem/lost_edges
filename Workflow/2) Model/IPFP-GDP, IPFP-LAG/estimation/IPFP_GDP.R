#------------------------------------------------------------------#
# File that implements the approach by Lebacher and Kauermann 2019 #
#------------------------------------------------------------------#

# Settings ----
# remove the old stuff
rm(list=ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
setwd("../../..")

# Load the packges
library(ipfp)
library(rlist)
library(nloptr)
library(pracma)
library(truncdist)
library(stargazer)

set.seed(42)

# load data
load("Data/reduced_Data.RData")

sender<-as.factor(iS[,1])
receiver<-as.factor(iS[,2])
u<-matrix(0,nrow=dim(iS)[1],ncol=n)
v<-matrix(0,nrow=dim(iS)[1],ncol=n)

for (i in 1:dim(iS)[1]){
  sender<-iS[i,1]
  receiver<-iS[i,2]
  u[i,sender]<-1
  v[i,receiver]<-1


} 


# Functions and settings ----

ll<-function(theta){
  logl<-sum(vec_E*Z%*%theta-exp(Z%*%theta))
  return(-logl)
}


llgradient<-function(theta){
  ob<-grad(f=ll,x0=theta)
  return(-ob)
}

constraint<-function(theta){
  
  return(A%*%exp(Z%*%theta)-Y)
  
}


cgradient<-function(theta){
  A%*%diag(c(exp(Z%*%theta)))%*%Z
}

local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG_EQ",
                    "xtol_rel" = 1.0e-20 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG_EQ",
              "xtol_rel" = 1.0e-20,
              "maxeval" = 100,"local_opts" = local_opts )



## Start the loop ----


for (tau in 1:N_t){

Z<-cbind(u,v,line_gdp[,tau])
comb<-X_inform[,tau]
Y=A%*%comb
X_hat<-ipfp(y=Y,A=A,x0=rep(1,n*(n-1)))
vec_E<-X_hat

nl<-nloptr(x0=rep(1,2*n+1), eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
old<-nl$solution

ob<-c()
ob2<-c()
crit=100
while (crit >0){
  vec_E<-exp(Z%*%nl$solution)
  nl<-nloptr(x0=old, eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
  crit<-t(nl$solution-old)%*%(nl$solution-old)
  old<-nl$solution

}

save.image(paste("2) Model/IPFP-GDP, IPFP-LAG/results/res_",tau,".RData",sep=""))

}


