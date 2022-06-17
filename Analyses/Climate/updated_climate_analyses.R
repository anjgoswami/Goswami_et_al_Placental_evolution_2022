library(RPANDA)
library(Morpho)
library(pspline)
library(BB)
library(numDeriv)

source('./Analyses/Climate/clim_pl2.R')
source('./Analyses/Climate/GIC-bis.R')  # GIC for the new function
source('./Analyses/Climate/clim_curves.R') # function to plot the reconstructed curve

## Load data
load(file="./Data/shape.data.322.R")#Y## procrustes coordinates
eco_dat <- read.csv("./Raw_Data/full_species_data.csv", row.names = 1)
phy <- read.nexus("./Raw_Data/trees/tree1_random_set/tree75_80_95.nex")
temperatures <- read.table("./Analyses/Climate/Cramer2011trend.csv", sep=";", dec=".", header=T)

## Sub-selection of landmarks only (1:66); curves(67:754), if subsetting desired
#subech<-c(1:66)
#Y.lm<-shape.data[c(subech),,]
#Y = vecx(Y.lm, byrow = TRUE)
#Y <- Y[phy$tip.label, ] # check that phylo and data order match

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
fit_eb <- mvgls(Y~1, tree=phy, model="EB", method = "H&L", error=TRUE)
fit_ou <- mvgls(Y~1, tree=phy, model="OU", method = "PL", error=TRUE)

#Climate models
fit_clim1 <- fit_mvenv(Y, phy, model=clim_1, nparam=1, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim2 <- fit_mvenv(Y, phy, model=clim_2, nparam=1, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim3 <- fit_mvenv(Y, phy, model=clim_3, nparam=2, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim4 <- fit_mvenv(Y, phy, model=clim_4, nparam=3, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim5 <- fit_mvenv(Y, phy, model=clim_5, nparam=1, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim6 <- fit_mvenv(Y, phy, model=clim_6, nparam=2, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim7 <- fit_mvenv(Y, phy, model=clim_7, nparam=3, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim8 <- fit_mvenv(Y, phy, model=clim_8, nparam=2, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim9 <- fit_mvenv(Y, phy, model=clim_9, nparam=3, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim10 <- fit_mvenv(Y, phy, model=clim_10, nparam=3, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim11 <- fit_mvenv(Y, phy, model=clim_11, nparam=4, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim12 <- fit_mvenv(Y, phy, model=clim_12, nparam=4, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)

fit_clim13 <- fit_mvenv(Y, phy, model=clim_13, nparam=4, method=c("RidgeArch"),
                       targM=c("unitVariance"), SE=TRUE)



###climate model for separate eco categories, use simm tree instead of normal tree, shown here just for the first model for example:
#fit_clim1b <- fit_mvenv(Y, tree_ecol, model=clim_1, nparam=1, method=c("RidgeArch"),
#                        targM=c("unitVariance"), SE=TRUE)

GIC_scores <- list(
  fit_bm = GIC(fit_bm),
  fit_eb = GIC(fit_eb),
  fit_ou = GIC(fit_ou),
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
res=aicw(sapply(GIC_scores, function(x) x$GIC))
res=aicw(results)
aic_t<-data.frame(sapply(res,c))
write.csv(aic_t, file = "./Analyses/Climate/results/aic.csv")


#######################still working on this one#########################
##################Plot rates on Time bins

totage <- max(nodeHeights(phy))
time_bins <- c(0,seq(0,totage, by = 5))
treemap <- make.era.map(phy, time_bins)
treemap = reorderSimmap(treemap, "postorder")


# BMM model like in mvBMM (with constraint="proportional") - WITH ERROR
fit=mvgls(Y~1, tree=treemap, model="BMM", method="H&L", error=TRUE)

# retrieve the results. We add "1" because the others rates are proportional changes to the first time-bin
x=c(1,fit$param)
names(x) = colnames(treemap$mapped.edge) # names according to the time-bin

# plot (using staircase lines) the rates through time
par(mfrow=c(2,1))
plot(time_bins, x[order(as.numeric(as.character(names(x))))], type="s", lwd=2, xaxt = "n", xlab='Times MA', ylab='Rates (relative)' )
axis(1, at=time_bins, labels=-round(totage-time_bins, digits = 3))

# same but without estimating an error term (variance at the tips)
# fit1=mvgls(Y~1, tree=treemap, model="BMM", method="H&L")
# x1=c(1,fit1$param)
# names(x1) = colnames(treemap$mapped.edge)

# 10Ma steps
time_bins10 <- c(0,seq(20,totage, by = 10)) # I start the time binning at 20Ma after the root as there is low information near the root anyway
treemap10 <- make.era.map(phy, time_bins10)
treemap10 = reorderSimmap(treemap10, "postorder")

fit2=mvgls(Y~1, tree=treemap10, model="BMM", method="H&L", error=TRUE)

# retrieve the results. We add "1" because the others rates are proportional changes to the first time-bin
x2=c(1,fit2$param)
names(x2) = colnames(treemap10$mapped.edge) # names according to the time-bin


# SAME with thinner time-bin
# 2.5Ma steps
time_bins2 <- c(0, seq(20, totage, by = 2.5))
treemap2 <- make.era.map(phy, time_bins2)
treemap2 = reorderSimmap(treemap2, "postorder")

fit3=mvgls(Y~1, tree=treemap2, model="BMM", method="H&L", error=TRUE)

# retrieve the results. We add "1" because the others rates are proportional changes to the first time-bin
x3=c(1,fit3$param)
names(x3) = colnames(treemap2$mapped.edge) # names according to the time-bin

# plot (using staircase lines) the rates through time for 10Ma bins
par(mfrow=c(2,1))
plot(time_bins10, x2[order(as.numeric(as.character(names(x2))))], type="s", lwd=2, xaxt = "n", xlab='Times MA', ylab='Rates (relative)' )
axis(1, at=time_bins10, labels=-round(totage-time_bins10, digits = 3))

# plot (using staircase lines) the rates (2.5Ma bin) through time
plot(time_bins2, x3[order(as.numeric(as.character(names(x3))))], type="s", lwd=2, xaxt = "n", xlab='Times MA', ylab='Rates (relative)' )
axis(1, at=time_bins2, labels=-round(totage-time_bins2, digits = 3))

# simulate tree and data
phy <- pbtree(n=322, scale=80)
Y <- mvSIM(phy, model="BM1", param=list(sigma=crossprod(matrix(runif(50*50),50))))
# check that the matrix of data and the tree names matches.
Y <- Y[phy$tip.label, ]
totage <- max(nodeHeights(phy))
time_points <- seq(20,totage, by = 5)
time_bins <- c(0,time_points) # we start the first bin 20Ma after the root
treemap <- make.era.map(phy, time_bins)
treemap = reorderSimmap(treemap, "postorder")
#plot(treemap)
fit=mvgls(Y~1, tree=treemap, model="BMM", method="H&L", error=TRUE)
x=fit$param
names(x) = colnames(treemap$mapped.edge)
# plot
plot(time_bins, c(x[1], x[order(as.numeric(as.character(names(x))))]), type="S", lwd=2, xaxt = "n", xlab='Times MA', ylab='Rates' )
axis(1, at=time_bins, labels=-round(totage-time_bins, digits = 3))

#COMPARE TO CLIMATE
par(mfrow=c(2,1))
#beta5.all=fitenv5.all$model.par
#f5<-function(x, beta.all){exp(beta.all[1]*env_data2((mtot)-x) + beta.all[2]*trend(x) + beta.all[3]*env_data((mtot)-x))} # The function used to link the rates to the climatic curve
#t <- seq(0, mtot, length.out=100) # a vector of times
#rates5.all <- f5(t, beta5.all) # Retrieve the rates estimates for each times steps
#plot(-t, rev(rates5.all), type='l', xlab="Times", ylab=bquote(paste("Evolutionary rates ", sigma)), main="reconstructed rates under the ClimExp5 model")


beta_ac_clim2.all=fitac_clim2.all$model.par
f_ac_clim2.all<-function(x, beta.all){exp(beta.all[1]*env_data2((mtot)-x) + beta.all[2]*trend_3(x) + beta.all[3]*env_data2((mtot)-x)*trend_3(x))} # The function used to link the rates to the climatic curve
t <- seq(0, mtot, length.out=100)
rates_ac_clim2.all <- f_ac_clim2.all(t, beta_ac_clim2.all) # Retrieve the rates estimates for each times steps
plot(-t, rev(rates_ac_clim2.all), type='l', xlab="Times", ylab=bquote(paste("Evolutionary rates ", sigma)), main="reconstructed rates under the ClimExp model")

#plot using staircase rates through time
plot(time_bins, x[order(as.numeric(as.character(names(x))))],type = "s", lwd=2, xaxt = "n", xlab = "Times MA", ylab = "Rates (relatve)")
axis(1, at = time_bins, labels = -round(totage-time_bins, digits = 0))
