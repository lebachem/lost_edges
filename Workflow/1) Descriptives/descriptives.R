#-----------------------------------------------------------#
# File that creates some descriptives                       #
#-----------------------------------------------------------#

# Settings ----
# remove the old stuff
rm(list=ls())

# load packages
library(devtools)
install_github("thomasp85/patchwork")
library(patchwork)
library(zoo)
library(igraph)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(plyr)
library(MazamaSpatialUtils)
library(countrycode)
library(stargazer)


# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
setwd("..")

# load the reduced data
load("Data/reduced_Data.RData")
countries_small=countries
combinations_small=combinations
n_small<-length(countries_small)
A_small<-routing_mat(n_small)$A
X_reduced<-X_inform[which(rownames(X_inform)%in%combinations_small),]
Y_small<-t(A_small%*%X_reduced)
N_small<-n_small*(n_small-1)


# load the full data
load("Data/full_data.RData")


## Plot the density and aggregated quantities--

# Density

d<-c()
d_small<-c()
for (tau in 1:N_t){
  d<-c(d,sum(X_inform[,tau]>0))
  d_small<-c(d_small,sum(X_reduced[,tau]>0))
}

# Number of zero-valued marginals
m<-c()
m_small<-c()
for (tau in 1:N_t){
m<-c(m,sum(Y[tau,]==0))
m_small<-c(m_small,sum(Y_small[tau,]==0))
}

## Plot the correlation between the marginals and the GDP for the reduce dataset

# start with the outgoing margins
Y_out = t(Y_small[,1:59])
Y_in = t(Y_small[,60:118])
cor_out = c()
cor_in = c()
cor_dyad = c()
cor_lag = c()
for (tau in 1:length(t)){
  cor_out = c(cor_out, cor(Y_out[,tau],node_gdp[,tau]))
  cor_in = c(cor_in, cor(Y_in[,tau], node_gdp[,tau]))
  cor_dyad = c(cor_dyad, cor(X_reduced[,tau], exp(line_gdp[,tau])))
}


tmon <- as.yearmon(2003 + seq(0, 181)/11)
tmon<-as.Date(tmon)

total<-rowSums(Y[,1:n])
total_small<-rowSums(Y_small[,1:n_small])

message<-data.frame(cbind(d/N,1-m/dim(Y)[1],total,tmon))
colnames(message)<-c("density","marginals","total","time")

message_small<-data.frame(cbind(d_small/N_small,1-m_small/dim(Y_small)[1],total_small,tmon,cor_in, cor_out, cor_dyad))
colnames(message_small)<-c("density","marginals","total","time")

# Density
p1<-ggplot(data=message,aes(x=tmon,y=density))+
  theme_bw()+geom_line()+
  scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+
  ylab("Density of the Network")+
  xlab("Time")+ggpubr::rotate_x_text()

p1_small<-ggplot(data=message_small,aes(x=tmon,y=density))+
  theme_bw()+geom_line()+
  scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+
  ylab("Density of the Network")+xlab("Time")+ggpubr::rotate_x_text()

# Share of valued marginals
p2<-ggplot(data=message,aes(x=tmon,y=marginals))+
  theme_bw()+geom_line()+ 
  scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+
  ylab("Share of valued marginals")+xlab("Time")+ggpubr::rotate_x_text()

p2_small<-ggplot(data=message_small,aes(x=tmon,y=marginals))+ 
  theme_bw()+geom_line()+ scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+ 
  ylab("Share of valued marginals")+xlab("Time")+ggpubr::rotate_x_text()

# Development of total Messages
p3<-ggplot(data=message,aes(x=tmon,y=total))+
  theme_bw()+geom_line()+
  scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+
    ylab("Total Messages")+xlab("Time")+ggpubr::rotate_x_text()

p3_small<-ggplot(data=message_small,aes(x=tmon,y=total))+
  theme_bw()+geom_line()+
  scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+ 
  ylab("Total Messages")+xlab("Time")+ggpubr::rotate_x_text()

# Correlation of marginals and gdps
p4_in<<-ggplot(data=message,aes(x=tmon,y=cor_in))+
  theme_bw()+geom_line()+
  scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+
  ylab("Correlation: GDP and valued indegree")+xlab("Time")+ggpubr::rotate_x_text()

p4_out<-ggplot(data=message_small,aes(x=tmon,y=cor_out))+
  theme_bw()+geom_line()+
  scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+
  ylab("Correlation: GDP and valued outdegree")+xlab("Time")+ggpubr::rotate_x_text()



p4_dyad<-ggplot(data=message_small,aes(x=tmon,y=cor_dyad))+
  theme_bw()+geom_line()+
  scale_x_date(date_labels = "%m/%Y",breaks = pretty(tmon, n = 10))+
  ylab("Correlation: (gdp_i+gdp_j) and valued edges")+xlab("Time")+ggpubr::rotate_x_text()

  pdf("1) Descriptives/Figures/gdp_correlation.pdf",width = 9,height = 4.5)
  p4_in + p4_out + p4_dyad
  dev.off()

pdf("1) Descriptives/Figures/density_and_total.pdf",width = 9,height = 4.5)
p1+p2+p3
dev.off()

pdf("1) Descriptives/Figures/density_and_total_small.pdf",width = 9,height = 4.5)
p1_small+p2_small+p3_small
dev.off()



## Plot all edges ----

dataf <-data.frame(t(X_inform),date = tmon)

dataf_long <- melt(dataf, id="date")  # convert to long format

p4<-ggplot(data=dataf_long, aes(x=date, y=value,colour=variable)) +
  geom_line()  + theme_bw()+
  scale_colour_grey()+ ylab("Messages sent per Dyad")+xlab("Time")+ggpubr::rotate_x_text()+
  theme(legend.position="none")  #  +  scale_y_continuous(trans='log10')

pdf("1) Descriptives/Figures/messages.pdf",width = 9,height = 5.5)
p4
dev.off()


  ## Degrees and degree distribution ----

dout<-c()
din<-c()
deg<-c()

dout_small<-c()
din_small<-c()
deg_small<-c()

for (tau in 1:N_t){
  
  new<-graph_from_adjacency_matrix(vec_to_mat(n,X_inform[,tau])$mat_bin)
  dout<-cbind(dout,degree(new,mode="out",normalized = F))
  din<-cbind(din,degree(new,mode="in",normalized = F))
  deg<-cbind(deg,degree(new,mode="total",normalized = F))
  
  new_small<-graph_from_adjacency_matrix(vec_to_mat(n_small,X_reduced[,tau])$mat_bin)
  dout_small<-cbind(dout_small,degree(new_small,mode="out",normalized = F))
  din_small<-cbind(din_small,degree(new_small,mode="in",normalized = F))
  deg_small<-cbind(deg_small,degree(new_small,mode="total",normalized = F))  
}



doutCS_small<-c()
doutCS<-c()
for (tau in 1:N_t){
  doutCS<-cbind(doutCS,cumsum(din[,tau]/sum(din[,tau])))
  doutCS_small<-cbind(doutCS_small,cumsum(din_small[,tau]/sum(din_small[,tau])))
  
}

mout<-apply(doutCS,1,median)
minout<-apply(doutCS,1,min)
maxout<-apply(doutCS,1,max)

mout_small<-apply(doutCS_small,1,median)
minout_small<-apply(doutCS_small,1,min)
maxout_small<-apply(doutCS_small,1,max)

df<-data.frame(mout,minout,maxout,idx=1:n)
p1<-ggplot(data=df,aes(x=idx,y=mout))+
  theme_bw()+
  geom_line()+geom_ribbon(aes(x = idx, ymax = maxout, ymin = minout), alpha = 0.9, fill = "gray")+
  ylab("Cumulative Degree Distribution")+xlab("Indegree")+geom_abline(slope = 1/n,intercept = 0)


df_small<-data.frame(mout_small,minout_small,maxout_small,idx=1:n_small)
p1_small<-ggplot(data=df_small,aes(x=idx,y=mout_small))+
  theme_bw()+geom_line()+
  geom_ribbon(aes(x = idx, ymax = maxout_small, ymin = minout_small), alpha = 0.9, fill = "gray")+
  ylab("Cumulative Degree Distribution")+xlab("Indegree")+geom_abline(slope = 1/n_small,intercept = 0)


doutCS<-c()
doutCS_small<-c()
for (tau in 1:N_t){
  doutCS<-cbind(doutCS,cumsum(dout[,tau]/sum(dout[,tau])))
  doutCS_small<-cbind(doutCS_small,cumsum(dout_small[,tau]/sum(dout_small[,tau])))
  
}

mout<-apply(doutCS,1,median)
minout<-apply(doutCS,1,min)
maxout<-apply(doutCS,1,max)

mout_small<-apply(doutCS_small,1,median)
minout_small<-apply(doutCS_small,1,min)
maxout_small<-apply(doutCS_small,1,max)


df<-data.frame(mout,minout,maxout,idx=1:n)
p2<-ggplot(data=df,aes(x=idx,y=mout))+
  theme_bw()+geom_line()+
  geom_ribbon(aes(x = idx, ymax = maxout, ymin = minout), alpha = 0.9, fill = "gray")+
  ylab("Cumulative Degree Distribution")+xlab("Outdegree")+geom_abline(slope = 1/n,intercept = 0)

df_small<-data.frame(mout_small,minout_small,maxout_small,idx=1:n_small)
p2_small<-ggplot(data=df_small,aes(x=idx,y=mout_small))+
  theme_bw()+geom_line()+
  geom_ribbon(aes(x = idx, ymax = maxout_small, ymin = minout_small), alpha = 0.9, fill = "gray")+
  ylab("Cumulative Degree Distribution")+xlab("Outdegree")+geom_abline(slope = 1/n_small,intercept = 0)


df<-data.frame(inlist=CS_plots(din,return_deg = T,indeg = F,n=n),outlist=CS_plots(dout,return_deg = T,indeg = F,n=n))

df_small<-data.frame(inlist=CS_plots(din_small,return_deg = T,indeg = F,n=n_small),
               outlist=CS_plots(dout_small,return_deg = T,indeg = F,n=n_small))
p3  <- ggplot(df, aes(x=inlist,y=outlist))+ geom_point(alpha=0.05)+
  theme_bw()+ xlab("Indegree")+
  ylab("Outdegree")+geom_abline(slope = 1,intercept = 0)

p3_small  <- ggplot(df_small, aes(x=inlist,y=outlist))+ 
  geom_point(alpha=0.05)+ theme_bw()+ xlab("Indegree")+
  ylab("Outdegree")+geom_abline(slope = 1,intercept = 0)


pdf("1) Descriptives/Figures/deg_dist.pdf",width = 9,height = 4.5)
p1+p2+p3
dev.off()

pdf("1) Descriptives/Figures/deg_dist_small.pdf",width = 9,height = 4.5)
p1_small+p2_small+p3_small
dev.off()

#rm(list=ls())
