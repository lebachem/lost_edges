#----------------------------------------------------------------#
# File to evaluate the probabilisitc reconstruction of the models#
#----------------------------------------------------------------#

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
library(SpecsVerification)

# load the data
load("Data/full_data.RData")
load("2) Model/IPFP/results/ipfp.RData")

# Extract the probabilistic structure ----
auc_roc1<-c()
auc_pr1<-c()
auc_roc2<-c()
auc_pr2<-c()
auc_roc_fit<-c()
auc_roc_ER<-c()
auc_pr_fit<-c()
auc_pr_ER<-c()

br1<-0
br2<-0
br3<-0
br4<-0
UNC<-0
SAVE=X_inform
for (tau in 1:N_t){
  
  p<-c()
  print(tau)
  
  # DC-GRAVITY
  load(file=paste("2) Model/DC-GRAVITY/results/z_",t[tau],".RData",sep=""))
  save=SAVE[,tau]
  Y<-A%*%save
  X_inform<-save
  fit<-fitness(Y)
  p_hat<-pij(z_hat)
  p<-rbind(p,p_hat)
  object<-roc_pr(p_hat,X_inform,ret_val=T)
  auc_roc2<-c(auc_roc2,object[[1]])
  auc_pr2<-c(auc_pr2,object[[2]])
  
  
  # IPFP
  # if the values of IPFP become to high, then the predictions
  # are exclusively zero
  p_hat<- 1-exp(-X_hat_ipfp[,tau])
  if (sum(p_hat==1)==length(p_hat)){
    auc_roc1<-c(auc_roc1,0.5)
    auc_pr1<-c(auc_pr1,0)  
    p<-rbind(p,p_hat)
  }
  
  if (sum(p_hat==1)<length(p_hat)){
  object<-roc_pr(p_hat,X_inform,ret_val = T)
  auc_roc1<-c(auc_roc1,object[[1]])
  auc_pr1<-c(auc_pr1,object[[2]])
  }
  # H-FIT
  load(paste("2) Model/H-FIT, H-ER/results/fit_",t[tau],".RData",sep=""))
  l1<-Reduce("+",lapply(Lsamp1$L,fun))
  vec1<-mat_to_vec(l1,dim(l1)[1])/100
  p<-rbind(p,vec1)
  object1<-roc_pr(vec1,X_inform,ret_val = T)
  auc_roc_fit<-c(auc_roc_fit,object1[[1]])
  auc_pr_fit<-c(auc_pr_fit,object1[[2]])
  
  # H-ER
  load(paste("2) Model/H-FIT, H-ER/results/ER_",t[tau],".RData",sep=""))
  l2<-Reduce("+",lapply(Lsamp2$L,fun))
  vec2<-mat_to_vec(l2,dim(l2)[1])/100
  p<-rbind(p,vec2)
  object2<-roc_pr(vec2,X_inform,ret_val=T)
  auc_roc_ER<-c(auc_roc_ER,object2[[1]])
  auc_pr_ER<-c(auc_pr_ER,object2[[2]])
  
  # aggregated Brier score
  real_bin<-as.numeric(X_inform>0)
  b1<-BrierDecomp(p[1,],real_bin)
  UNC<-UNC+b1[1,3]
  b1[,3]=b1[,3]-b1[,2]
  b2<-BrierDecomp(p[2,],real_bin)
  b2[,3]=b2[,3]-b2[,2]
  b3<-BrierDecomp(p[3,],real_bin)
  b3[,3]=b3[,3]-b3[,2]
  b4<-BrierDecomp(p[4,],real_bin)
  b4[,3]=b4[,3]-b4[,2]
  br1<-b1+br1
  br2<-b2+br2
  br3<-b3+br3
  br4<-b4+br4
  
}


# Plot the Brier score decomposition
data<-data.frame(c("REL","UNC-RES"),br2[1,c(1,2)],br4[1,c(1,2)],br3[1,c(1,2)],br1[1,c(1,2)])
colnames(data)=c("id","IPFP","H-ER","H-FIT","DC-GRAVITY")
data_long <- melt(data,id="id")  # convert to long format

pdf("3) Evaluation/Figures/decomp.pdf",width = 9,height = 2)
p<-ggplot(data_long, aes(fill=id, y=value, x=variable))+
  geom_bar( stat="identity")+labs(fill = "Decompositon")+xlab("Model")+ggpubr::rotate_x_text()+ coord_flip()
p +  geom_errorbar(aes(ymax=UNC, ymin=UNC), linetype="dashed") + theme_bw()
dev.off()

# Plot ROC and PR
tmon <- as.yearmon(2003 + seq(0, 181)/12)
tmon<-as.Date(tmon)

data <- data.frame( H_FIT=auc_roc_fit,DC_GRAVITY=auc_roc2, IPFP=auc_roc1,H_ER=auc_roc_ER, date = tmon )
data_long <- melt(data, id="date")  # convert to long format

p1<-ggplot(data=data_long, aes(x=date, y=value,color=factor(variable,labels =c("H-FIT","DC-GRAVITY","IPFP","H-ER") ))) +
  geom_line()  + theme_bw()+labs(color = "Model")+
  ylab("Receiver-Operating-Characteristic, AUC values")+xlab("Time")+ggpubr::rotate_x_text()+
  scale_color_brewer(palette="Paired")

data <- data.frame( H_FIT=auc_pr_fit,DC_GRAVITY=auc_pr2, IPFP=auc_pr1,H_ER=auc_pr_ER, date = tmon )
data_long <- melt(data, id="date")  # convert to long format

p2<-ggplot(data=data_long, aes(x=date, y=value,color=factor(variable,labels =c("H-FIT","DC-GRAVITY","IPFP","H-ER") ))) +
  geom_line()  + theme_bw()+labs(color = "Model")+
  ylab("Precision-Recall, AUC values")+xlab("Time")+ggpubr::rotate_x_text()+
  theme(legend.position="none")+
  scale_color_brewer(palette="Paired")


# joint legend
tmp <- ggplot_gtable(ggplot_build(p1))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
p1 <- p1 + theme(legend.position="none")
  
pdf("3) Evaluation/Figures/PR_ROC.pdf",width = 9,height = 5)
grid.arrange(p1, legend,p2, ncol=3, widths=c(2.3,0.8, 2.2))
dev.off()


# save ----
save.image("3) Evaluation/Figures/data_roc_pre_eval.RData")
