#--------------------------------------------------------------#
# File to evaluate the degree reconstruction of the models     #
#--------------------------------------------------------------#

# Settings ----
# remove the old stuff
rm(list=ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
setwd("../")

# set seed
set.seed(42)

# Load the packges
library(ipfp)
library(PRROC)
library(GA)
library(zoo)
library(igraph)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(reshape2)
library(plyr)
library(gridExtra)


# load the data
load("2) Model/IPFP/results/ipfp_small.RData")
load("2) Model/GRAVITY/results/gravity_small.RData")
load("2) Model/LASSO/results/lasso_small.RData")
load("Data/reduced_Data.RData")

# Extract the binary structure ----

d1<-c()
d2<-c()
d3<-c()
d4<-c()

for (tau in 2:N_t){

print(tau)

# Extract the real density  
dens_real<-mean(X_inform[,tau]>0)
  
#1 IPFP
est1<-X_hat_ipfp[,tau]
survive1<-which(est1%in%est1[rev(order(est1))][1:(dens_real*N)])
est1[est1>0]<-0
est1[survive1]<-1

#2 GRAVITY
est2<-X_hat_gravity[,tau]
survive2<-which(est2%in%est2[rev(order(est2))][1:(dens_real*N)])
est2[est2>0]<-0
est2[survive2]<-1

#3 DC-GRAVITY
load(file=paste("2) Model/DC-GRAVITY/results/z_small_",t[tau],".RData",sep=""))
Y<-A%*%X_inform[,tau]
fit<-fitness(Y)
est3<-pij(z_hat)
survive3<-which(est3%in%est3[rev(order(est2))][1:(dens_real*N)])
est3[est3>0]<-0
est3[survive3]<-1
  
#4 H-FIT  
load(paste("2) Model/H-FIT, H-ER/results/small_fit_",t[tau],".RData",sep=""))
l1<-Reduce("+",lapply(Lsamp1$L,fun))
est4<-mat_to_vec(l1,dim(l1)[1])/100
survive4<-which(est4%in%est4[rev(order(est3))][1:(dens_real*N)])
est4[est4>0]<-0
est4[survive4]<-1  

#5 H-ER  
load(paste("2) Model/H-FIT, H-ER/results/small_ER_",t[tau],".RData",sep=""))
l2<-Reduce("+",lapply(Lsamp2$L,fun))
est5<-mat_to_vec(l2,dim(l2)[1])/100
survive5<-which(est5%in%est5[rev(order(est5))][1:(dens_real*N)])
est5[est5>0]<-0
est5[survive5]<-1

#6 DC-GRAVITY-GDP
load(file=paste("2) Model/DC-GRAVITY/results/z_small_cov_",t[tau],".RData",sep=""))
Y<-A%*%X_inform[,tau]
fit<-line_gdp[,tau]
est6<-pij(z_hat)
survive6<-which(est6%in%est6[rev(order(est6))][1:(dens_real*N)])
est6[est6>0]<-0
est6[survive6]<-1
  
#7 DC-GRAVITY-LAG
load(file=paste("2) Model/DC-GRAVITY/results/z_small_cov_lag_",t[tau],".RData",sep=""))
Y<-A%*%X_inform[,tau]
fit<-log(1.1+X_inform[,tau-1])
est7<-pij(z_hat)
survive7<-which(est7>0.5)
est7[est7>0]<-0
est7[survive7]<-1
  
#8 IPFP-LAG  
load(file=paste("2) Model/IPFP-GDP, IPFP-LAG/results/res_lag_",tau,".RData",sep=""))
est8=vec_E
survive8<-which(vec_E%in%vec_E[rev(order(vec_E))][1:(dens_real*N)])
est8[est8>0]<-0
est8[survive8]<-1

#9 IPFP-GDP
load(file=paste("2) Model/IPFP-GDP, IPFP-LAG/results/res_",tau,".RData",sep=""))
est9=vec_E
survive9<-which(vec_E%in%vec_E[rev(order(vec_E))][1:(dens_real*N)])
est9[est9>0]<-0
est9[survive9]<-1
  
#10 LASSO  
est10<-X_hat_lasso[,tau]
est10[est10>0]<-1
  
#11 TOMOGRAVITY
load(file=paste("2) Model/TOMOGRAVITY/results/fit_",t[tau],".RData",sep=""))
est11=fit
survive11<-which(est11%in%est11[rev(order(est11))][1:(dens_real*N)])
est11[est11>0]<-0
est11[survive11]<-1
  

  
  element1<-stat_plotter(X_inform[,tau],est1,return_diff = T,plot=F)
  element2<-stat_plotter(X_inform[,tau],est2,return_diff = T,plot=F)
  element3<-stat_plotter(X_inform[,tau],est3,return_diff = T,plot=F)
  element4<-stat_plotter(X_inform[,tau],est4,return_diff = T,plot=F)
  element5<-stat_plotter(X_inform[,tau],est5,return_diff = T,plot=F)
  element6<-stat_plotter(X_inform[,tau],est6,return_diff = T,plot=F)
  element7<-stat_plotter(X_inform[,tau],est7,return_diff = T,plot=F)
  element8<-stat_plotter(X_inform[,tau],est8,return_diff = T,plot=F)
  element9<-stat_plotter(X_inform[,tau],est9,return_diff = T,plot=F)
  element10<-stat_plotter(X_inform[,tau],est10,return_diff = T,plot=F)
  element11<-stat_plotter(X_inform[,tau],est11,return_diff = T,plot=F)

  
  
  d1<-cbind(d1,c(element1[[1]],element2[[1]],element3[[1]],element4[[1]],element5[[1]],element6[[1]],element7[[1]],element8[[1]],element9[[1]],element10[[1]],element11[[1]]))
  d2<-cbind(d2,c(element1[[2]],element2[[2]],element3[[2]],element4[[2]],element5[[2]],element6[[2]],element7[[2]],element8[[2]],element9[[2]],element10[[2]],element11[[2]]))
  d3<-cbind(d3,c(element1[[3]],element2[[3]],element3[[3]],element4[[3]],element5[[3]],element6[[3]],element7[[3]],element8[[3]],element9[[3]],element10[[3]],element11[[3]]))
  d4<-cbind(d4,c(element1[[4]],element2[[4]],element3[[4]],element4[[4]],element5[[4]],element6[[4]],element7[[4]],element8[[4]],element9[[4]],element10[[4]],element11[[4]]))
  }


# Plot the RMSE ----

# Define a variable for the time
tmon <- as.yearmon(2003 + seq(0, 180)/12)
tmon<-as.Date(tmon)

data <- data.frame( CIMI_gdp=d1[6,]/sqrt(N),
                    LASSO=d1[10,]/sqrt(N),
                    regression_gdp=d1[9,]/sqrt(N),
                    CIMI=d1[3,]/sqrt(N),
                    ER=d1[5,]/sqrt(N),
                    TOMO=d1[11,]/sqrt(N),
                    IPFP=d1[1,]/sqrt(N),
                    Gravity=d1[2,]/sqrt(N),
                    FIT=d1[4,]/sqrt(N),
                    regression_lag=d1[8,]/sqrt(N),
                    CIMI_lag=d1[7,]/sqrt(N), date = tmon )
data_long <- melt(data, id="date")  # convert to long format

# Outdegree
p1<-ggplot(data=data_long, aes(x=date, y=value,color=factor(variable,labels =c("DC-GRAVITY-GDP",
                                                                               "LASSO",
                                                                               "IPFP-GDP",
                                                                               "DC-GRAVITY",
                                                                               "H-ER",
                                                                               "TOMOGRAVITY",
                                                                               "IPFP",
                                                                               "GRAVITY",
                                                                               "H-FIT",
                                                                               "IPFP-LAG",
                                                                               "DC-GRAVITY-LAG")))) +
  geom_line()  + theme_bw()+labs(color = "Model") +
  ylab("RMSE Outdegree")+xlab("Time")+ggpubr::rotate_x_text()+
  scale_color_brewer(palette="Paired")

data <- data.frame( CIMI_gdp=d2[6,]/sqrt(N),
                    LASSO=d2[10,]/sqrt(N),
                    regression_gdp=d2[9,]/sqrt(N),
                    CIMI=d2[3,]/sqrt(N),
                    ER=d2[5,]/sqrt(N),
                    TOMO=d2[11,]/sqrt(N),
                    IPFP=d2[1,]/sqrt(N),
                    Gravity=d2[2,]/sqrt(N),
                    FIT=d2[4,]/sqrt(N),
                    regression_lag=d2[8,]/sqrt(N),
                    CIMI_lag=d2[7,]/sqrt(N), date = tmon )

data_long <- melt(data, id="date")  # convert to long format

p2<-ggplot(data=data_long, aes(x=date, y=value,color=factor(variable,labels =c("LASSO",
                                                                               "H-ER",
                                                                               "DC-GRAVITY-GDP",
                                                                               "GRAVITY",
                                                                               "IPFP",
                                                                               "TOMOGRAVITY",
                                                                               "H-FIT",
                                                                               "IPFP-GDP",
                                                                               "DC-GRAVITY",
                                                                               "IPFP-LAG",
                                                                               "DC-GRAVITY-LAG")))) +
  geom_line()  + theme_bw()+labs(color = "Model") +
  ylab("RMSE Indegree")+xlab("Time")+ggpubr::rotate_x_text()+
theme(legend.position="none")+
  scale_color_brewer(palette="Paired")


# joint legend
tmp <- ggplot_gtable(ggplot_build(p1))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
p1 <- p1 + theme(legend.position="none")

# plot the stuff
pdf("3) Evaluation/Figures/degree_recon_small.pdf",width = 9,height = 5)
grid.arrange(p1, legend,p2, ncol=3, widths=c(2.3,1, 2.2))
dev.off()


# Plot the binary network structure ----

#1 IPFP
pdf("3) Evaluation/Figures/degree_recon_ipfp_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est1,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_ipfp_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est1)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#2 GRAVITY
pdf("3) Evaluation/Figures/degree_recon_gravity_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est2,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_gravity_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est2)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#3 DC-GRAVITY
pdf("3) Evaluation/Figures/degree_recon_cimi_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est3,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_cimi_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est2)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#4 H-FIT
pdf("3) Evaluation/Figures/degree_recon_fit_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est4,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_fit_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est3)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#5 H-ER
pdf("3) Evaluation/Figures/degree_recon_ER_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est5,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_ER_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est5)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#6 DC-GRAVITY-GDP
pdf("3) Evaluation/Figures/degree_recon_cimi_gpd_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est6,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_cimi_gdp_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est6)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#7 DC-GRAVITY-LAG
pdf("3) Evaluation/Figures/degree_recon_cimi_lag_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est7,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_cimi_lag_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est7)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#8 IPFP-LAG
pdf("3) Evaluation/Figures/degree_recon_regression_lag_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est8,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_regression_lag_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est8)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#9 IPFP-GDP
pdf("3) Evaluation/Figures/degree_recon_regression_gdp_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est9,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_regression_gdp_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est9)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#10 LASSO
pdf("3) Evaluation/Figures/degree_recon_lasso_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est10,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_lasso_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est10)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#11 TOMOGRAVITY
pdf("3) Evaluation/Figures/degree_recon_tomo_small.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est11,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_tomo_small.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est11)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()


#save ----
save.image("3) Evaluation/Figures/data_degree_eval_small.RData")

