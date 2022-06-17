library(RPANDA)
library(Morpho)
library(pspline)
library(BB)
library(numDeriv)

# Retrieve value from slurm job
id_job<-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
path_id <- paste("./results/tree1_",id_job, sep = "")

source('clim_pl2.R') 
source('GIC-bis.R')  # GIC for the new function
source('clim_curves.R') # function to plot the reconstructed curve

## Load data
load(file="/mnt/shared/projects/nhm/goswamilab/Eutheria/climate/shape.data.322.R")#Y## procrustes coordinates 
eco_dat <- read.csv("/mnt/shared/projects/nhm/goswamilab/Eutheria/climate/full_species_data.csv", row.names = 1)
sample_trees <- read.nexus("/mnt/shared/projects/nhm/goswamilab/Eutheria/rates/tree1/modules/t1_list_of_trees.nex")
phy<-sample_trees[[id_job]]

temperatures <- read.table("/mnt/shared/projects/nhm/goswamilab/Eutheria/climate/Cramer2011trend.csv", sep=";", dec=".", header=T)

## Sub-selection of landmarks only (1:66); curves(67:754), if subsetting desired
#subech<-c(1:66)
#Y.lm<-shape.data[c(subech),,]
Y<-shape.data
Y = vecx(Y, byrow = TRUE)
Y <- Y[phy$tip.label, ] # check that phylo and data order match

#If want to separate analyses by eco categories:
# Biogeography
#grp <- eco_dat$Biogeographic_Regions
#names(grp)=rownames(eco_dat)

# Making the simmap tree with mapped states (here we use a reconstruction under ARD, but more fancy models can be used! I think that now corHMM allows reconstructing simmap trees too)
#tree_ecol <- make.simmap(phy, grp, model="ARD", nsim=1)

## Prepare data
cur1 = scale(temperatures[,3],min(temperatures[,3],na.rm=T),max(temperatures[,3], na.rm=T)-min(temperatures[,3],na.rm=T))
env_data = splinefun(temperatures[,1], cur1) # scaled (raw) temperature curve

# Let's try to remove the global cooling trend from the climatic curve
mtot = max(nodeHeights(phy))
indtime = (temperatures[,1] <= mtot+1)
mod=lm(temperatures[indtime,3]~temperatures[indtime,1])

# We make a second function
cc = temperatures[indtime,3]-mod$fitted.values
#will produce an error on object length, it's fine

cur1=scale(cc,min(cc,na.rm=T),max(cc, na.rm=T)-min(cc,na.rm=T))
env_data2=splinefun( temperatures[indtime,1],cur1) # scaled (detrended) temperature curve

# Let's define a function for the trend that we removed - an affine function [y = ax + b]
trend = function(x) mod$coefficients[1] + mod$coefficients[2]*x 

#looking at the derivative of the temperature curve
# Now perform the differentiation:
gr = grad(func = env_data, temperatures[indtime,1])
#plot(gr~temperatures[indtime,1], type="l", ylim=c(-0.5,0.5))

#plot(diff(env_data2( temperatures[indtime,1]))~ temperatures[indtime,1][-739], ylim=c(-0.05,0.05), type="l")
require(pspline)

## 10. Loading required package: pspline
spline_result <- sm.spline(x=temperatures[indtime,1], y=gr, df=50) # can play with df to smooth the curve
env_func <- function(t){predict(spline_result,t)}
t<-unique(temperatures[indtime,1])

# We build the interpolated smoothing spline function
env_data3<-splinefun(t,env_func(t))
#curve(env_data3, 0, mtot)

#Prepare the models and fit them:
#clim_1: beta×Temperatures(t)
#clim_2: beta×Detrended_temperature(t)
#clim_3: beta1×trend(t) + beta2×Detrended_temperature(t)
#clim_4: beta1×trend(t) + beta2×Detrended_temperature(t) + beta3 × Detrended_temperature(t) × trend(t) - this one includes an interaction btw temp and trend
#clim_5-9: repeats but with derivative climate curve and interactions with either trend or temperature
#clim_1b: beta_australasian×Temperatures(t) | beta_Neotropical×Temperatures(t)

## Climatic models
clim_1 <- function(x, beta) exp(beta*env_data(mtot-x))
clim_2 <- function(x, beta) exp(beta*env_data2(mtot-x))
clim_3 <- function(x, beta) exp(beta[1]*trend(mtot-x) + beta[2]*env_data2(mtot-x))
clim_4 <- function(x, beta) exp(beta[1]*trend(mtot-x) + beta[2]*env_data2(mtot-x) + beta[3]*env_data2(mtot-x)*trend(mtot-x))
clim_5 <- function(x, beta) exp(beta*env_data3(mtot-x))
clim_6 <- function(x, beta) exp(beta[1]*trend(mtot-x) + beta[2]*env_data3(mtot-x))
clim_7 <- function(x, beta) exp(beta[1]*trend(mtot-x) + beta[2]*env_data3(mtot-x) + beta[3]*env_data3(mtot-x)*trend(mtot-x))
clim_8 <- function(x, beta) exp(beta[1]*env_data2(mtot-x) + beta[2]*env_data3(mtot-x))
clim_9 <- function(x, beta) exp(beta[1]*env_data2(mtot-x) + beta[2]*env_data3(mtot-x) + beta[3]*env_data2(mtot-x)*env_data3(mtot-x))
clim_10 <- function(x, beta) exp(beta[1]*trend(mtot-x) + beta[2]*env_data2(mtot-x)+ beta[3]*env_data3(mtot-x))
clim_11 <- function(x, beta) exp(beta[1]*trend(mtot-x) + beta[2]*env_data2(mtot-x) + beta[3]*env_data3(mtot-x)+ beta[4]*trend(mtot-x)*env_data2(mtot-x))
clim_12 <- function(x, beta) exp(beta[1]*trend(mtot-x) + beta[2]*env_data2(mtot-x) + beta[3]*env_data3(mtot-x)+ beta[4]*trend(mtot-x)*env_data3(mtot-x))
clim_13 <- function(x, beta) exp(beta[1]*trend(mtot-x) + beta[2]*env_data2(mtot-x) + beta[3]*env_data3(mtot-x)+ beta[4]*env_data2(mtot-x)*env_data3(mtot-x))

## Fit models
#Standard models:
fit_bm <- mvgls(Y~1, tree=phy, model="BM", method = "H&L", error=TRUE)
names_id<-paste(path_id,"_fit_bm.Rdata", sep="")
save(fit_bm, file = names_id) 
fit_eb <- mvgls(Y~1, tree=phy, model="EB", method = "H&L", error=TRUE)
names_id<-paste(path_id,"_fit_eb.Rdata", sep="")
save(fit_eb, file = names_id)  
#fit_ou <- mvgls(Y~1, tree=phy, model="OU", method = "PL", error=TRUE) 

#Climate models
fit_clim1 <- fit_mvenv(Y, phy, model=clim_1, nparam=1, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim1.Rdata", sep="")
save(fit_clim1, file = names_id)  
fit_clim2 <- fit_mvenv(Y, phy, model=clim_2, nparam=1, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim2.Rdata", sep="")
save(fit_clim2, file = names_id)   
fit_clim3 <- fit_mvenv(Y, phy, model=clim_3, nparam=2, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim3.Rdata", sep="")
save(fit_clim3, file = names_id)   
fit_clim4 <- fit_mvenv(Y, phy, model=clim_4, nparam=3, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim4.Rdata", sep="")
save(fit_clim4, file = names_id)  
fit_clim5 <- fit_mvenv(Y, phy, model=clim_5, nparam=1, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim5.Rdata", sep="")
save(fit_clim5, file = names_id)  
fit_clim6 <- fit_mvenv(Y, phy, model=clim_6, nparam=2, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim6.Rdata", sep="")
save(fit_clim6, file = names_id)  
fit_clim7 <- fit_mvenv(Y, phy, model=clim_7, nparam=3, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim7.Rdata", sep="")
save(fit_clim7, file = names_id)   
fit_clim8 <- fit_mvenv(Y, phy, model=clim_8, nparam=2, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim8.Rdata", sep="")
save(fit_clim8, file = names_id)  
fit_clim9 <- fit_mvenv(Y, phy, model=clim_9, nparam=3, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim9.Rdata", sep="")
save(fit_clim9, file = names_id)   
fit_clim10 <- fit_mvenv(Y, phy, model=clim_10, nparam=3, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim10.Rdata", sep="")
save(fit_clim10, file = names_id)  
fit_clim11 <- fit_mvenv(Y, phy, model=clim_11, nparam=4, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim11.Rdata", sep="")
save(fit_clim11, file = names_id)   
fit_clim12 <- fit_mvenv(Y, phy, model=clim_12, nparam=4, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
snames_id<-paste(path_id,"_fit_clim12.Rdata", sep="")
save(fit_clim12, file = names_id)  
fit_clim13 <- fit_mvenv(Y, phy, model=clim_13, nparam=4, method=c("RidgeArch"), 
                       targM=c("unitVariance"), SE=TRUE)
names_id<-paste(path_id,"_fit_clim13.Rdata", sep="")
save(fit_clim13, file = names_id)  


###climate model for separate eco categories, use simm tree instead of normal tree, shown here just for the first model for example:
#fit_clim1b <- fit_mvenv(Y, tree_ecol, model=clim_1, nparam=1, method=c("RidgeArch"), 
#                        targM=c("unitVariance"), SE=TRUE)

GIC_scores <- list( 
  fit_bm = GIC(fit_bm),
  fit_eb = GIC(fit_eb),
  fit_clim1 = GIC3(fit_clim1),
  fit_clim2 = GIC3(fit_clim2),
  fit_clim3 = GIC3(fit_clim3),
  fit_clim4 = GIC3(fit_clim4),
  fit_clim5 = GIC3(fit_clim5),
  fit_clim6 = GIC3(fit_clim6),
  fit_clim7 = GIC3(fit_clim7),
  fit_clim8 = GIC3(fit_clim8),
  fit_clim9 = GIC3(fit_clim9),
  fit_clim10 = GIC3(fit_clim10),
  fit_clim11 = GIC3(fit_clim11),
  fit_clim12 = GIC3(fit_clim12),
  fit_clim13 = GIC3(fit_clim13)
  )

# compare the fit (ranks, weights...)
aic<-aicw(sapply(GIC_scores, function(x) x$GIC))
aic_table<-paste(path_id,"aic.csv",sep="")
write.csv(aic, file = aic_table)
