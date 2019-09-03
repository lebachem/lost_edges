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
load("Data/full_data.RData")
load("2) Model/IPFP/results/ipfp.RData")
load("2) Model/GRAVITY/results/gravity.RData")
load("2) Model/LASSO/results/lasso.RData")


# Extract the binary structure ----
d1<-c()
d2<-c()
d3<-c()
d4<-c()

for (tau in 1:N_t){

print(tau)

# extract the real density
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
load(file=paste("2) Model/DC-GRAVITY/results/z_",t[tau],".RData",sep=""))
Y<-A%*%X_inform[,tau]
fit<-fitness(Y)
est3<-pij(z_hat)
survive3<-which(est3%in%est3[rev(order(est2))][1:(dens_real*N)])
est3[est3>0]<-0
est3[survive3]<-1

#4 H-FIT
load(paste("2) Model/H-FIT, H-ER/results/fit_",t[tau],".RData",sep=""))
l1<-Reduce("+",lapply(Lsamp1$L,fun))
est4<-mat_to_vec(l1,dim(l1)[1])/100
survive4<-which(est4%in%est4[rev(order(est3))][1:(dens_real*N)])
est4[est4>0]<-0
est4[survive4]<-1

#5 H-ER
load(paste("2) Model/H-FIT, H-ER/results/ER_",t[tau],".RData",sep=""))
l2<-Reduce("+",lapply(Lsamp2$L,fun))
est5<-mat_to_vec(l2,dim(l2)[1])/100
survive5<-which(est5%in%est5[rev(order(est5))][1:(dens_real*N)])
est5[est5>0]<-0
est5[survive5]<-1

#7 LASSO
est6<-X_hat_lasso[,tau]
est6[est6>0]<-1

element1<-stat_plotter(X_inform[,tau],est1,return_diff = T,plot=F)
element2<-stat_plotter(X_inform[,tau],est2,return_diff = T,plot=F)
element3<-stat_plotter(X_inform[,tau],est3,return_diff = T,plot=F)
element4<-stat_plotter(X_inform[,tau],est4,return_diff = T,plot=F)
element5<-stat_plotter(X_inform[,tau],est5,return_diff = T,plot=F)
element6<-stat_plotter(X_inform[,tau],est6,return_diff = T,plot=F)


d1<-cbind(d1,c(element1[[1]],element2[[1]],element3[[1]],element4[[1]],element5[[1]],element6[[1]]))
d2<-cbind(d2,c(element1[[2]],element2[[2]],element3[[2]],element4[[2]],element5[[2]],element6[[2]]))
d3<-cbind(d3,c(element1[[3]],element2[[3]],element3[[3]],element4[[3]],element5[[3]],element6[[3]]))
d4<-cbind(d4,c(element1[[4]],element2[[4]],element3[[4]],element4[[4]],element5[[4]],element6[[4]]))
}


# Plot the RMSE ----

# Define a variable for the time
tmon <- as.yearmon(2003 + seq(0, 181)/12)
tmon<-as.Date(tmon)


data <- data.frame(
                   LASSO=d1[6,]/sqrt(N),
                   H_ER=d1[5,]/sqrt(N),
                   DC_GRAVITY=d1[3,]/sqrt(N),
                   H_FIT=d1[4,]/sqrt(N),
                   GRAVITY=d1[2,]/sqrt(N),
                   IPFP=d1[1,]/sqrt(N), date = tmon )
data_long <- melt(data, id="date")  # convert to long format

# Outdegree
p1<-ggplot(data=data_long, aes(x=date, y=value,color=factor(variable,labels =c("LASSO",
                                                                               "H-ER",
                                                                               "DC-GRAVITY",
                                                                               "H-FIT",
                                                                               "GRAVITY",
                                                                               "IPFP")))) +
                                                                               
  geom_line()  + theme_bw()+labs(color = "Model")+
  ylab("RMSE Outdegree")+xlab("Time")+ggpubr::rotate_x_text()+
  scale_color_brewer(palette="Paired")

data <- data.frame(LASSO=d2[6,]/sqrt(N),
                   H_ER=d2[5,]/sqrt(N),
                   DC_GRAVITY=d2[3,]/sqrt(N),
                   H_FIT=d2[4,]/sqrt(N),
                   GRAVITY=d2[2,]/sqrt(N),
                   IPFP=d2[1,]/sqrt(N), date = tmon )

data_long <- melt(data, id="date")  # convert to long format

# Indegree
p2<-ggplot(data=data_long, aes(x=date, y=value,color=factor(variable,labels =c("LASSO",
                                                                               "H-ER",
                                                                               "DC-GRAVITY",
                                                                               "H-FIT",
                                                                               "GRAVITY",
                                                                               "IPFP") ))) +
  geom_line()  + theme_bw()+labs(color = "Model")+
  ylab("RMSE Indegree")+xlab("Time")+ggpubr::rotate_x_text()+
  theme(legend.position="none")+
  scale_color_brewer(palette="Paired")

  
# joint legend
tmp <- ggplot_gtable(ggplot_build(p1))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
p1 <- p1 + theme(legend.position="none")

# Plot the stuff
pdf("3) Evaluation/Figures/degree_recon.pdf",width = 9,height = 5)
grid.arrange(p1, legend,p2, ncol=3, widths=c(2.3,0.8, 2.2))
dev.off()


# Plot the binary network structure ----

#1 IPFP
pdf("3) Evaluation/Figures/degree_recon_ipfp.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est1,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_ipfp.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est1)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#2 GRAVITY
pdf("3) Evaluation/Figures/degree_recon_gravity.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est2,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_gravity.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est2)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#3 DC-GRAVITY
pdf("3) Evaluation/Figures/degree_recon_cimi.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est3,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_cimi.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est2)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#4 H-FIT
pdf("3) Evaluation/Figures/degree_recon_fit.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est4,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_fit.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est3)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#5 H-ER
pdf("3) Evaluation/Figures/degree_recon_ER.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est5,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_ER.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est5)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()

#7 LASSO
pdf("3) Evaluation/Figures/degree_recon_lasso.pdf",width = 6,height = 2)
stat_plotter(X_inform[,tau],est6,return_diff = T)
dev.off()
pdf("3) Evaluation/Figures/recon_lasso.pdf",width = 10,height = 5)
par(mfrow=c(1,2))
matplotter(as.matrix(vec_to_mat(n,est6)$mat_bin))
matplotter(as.matrix(vec_to_mat(n,X_inform[,tau])$mat_bin))
dev.off()


# save ----
save.image("3) Evaluation/Figures/data_degree_eval.RData")
  
