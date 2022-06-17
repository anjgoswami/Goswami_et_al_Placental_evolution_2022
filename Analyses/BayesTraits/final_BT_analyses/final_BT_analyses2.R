# Prepping files for BayesTraits Rates Analysis ----------------------------------------------
library(here)
library(ape)
library(phytools)
library(paleotree)
library(tidyverse)
library(SURGE)
library(geomorph)

load("./Data/shape.data.R")

#starter files for module-specific lm+curve lists
curve_table <- read_csv("./Raw_Data/placental_curves.csv")
my_curves <- create_curve_info(curve_table, n_fixed = 66)

lm_premax_d <- my_curves$Curves[which(curve_table$Module%in%c("premax_d"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_max_d <- my_curves$Curves[which(curve_table$Module%in%c("max_d"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_frontal <- my_curves$Curves[which(curve_table$Module%in%c("frontal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_nasal <- my_curves$Curves[which(curve_table$Module%in%c("nasal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_jugal <- my_curves$Curves[which(curve_table$Module%in%c("jugal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_parietal <- my_curves$Curves[which(curve_table$Module%in%c("parietal"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_premax_v <- my_curves$Curves[which(curve_table$Module%in%c("premax_v"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_max_v <- my_curves$Curves[which(curve_table$Module%in%c("max_v"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_squam_v <- my_curves$Curves[which(curve_table$Module%in%c("squam_v"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_glenoid <- my_curves$Curves[which(curve_table$Module%in%c("glenoid"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_squam_z <- my_curves$Curves[which(curve_table$Module%in%c("squam_z"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_supraocc <- my_curves$Curves[which(curve_table$Module%in%c("supraocc"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_occ_cond <- my_curves$Curves[which(curve_table$Module%in%c("occ_cond", "supraocc_cond"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_basiocc <- my_curves$Curves[which(curve_table$Module%in%c("basiocc"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_basisph <- my_curves$Curves[which(curve_table$Module%in%c("basisph"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_pterygoid <- my_curves$Curves[which(curve_table$Module%in%c("pterygoid"))]%>%unlist(.)%>%unique(.)%>%sort(.)
lm_palatine <- my_curves$Curves[which(curve_table$Module%in%c("palatine"))]%>%unlist(.)%>%unique(.)%>%sort(.)

lm_list <- list(
  lm_premax_d,
  lm_max_d,
  lm_frontal,
  lm_nasal,
  lm_jugal,
  lm_parietal,
  lm_premax_v,
  lm_max_v,
  lm_squam_v,
  lm_glenoid,
  lm_squam_z,
  lm_supraocc,
  lm_occ_cond,
  lm_basiocc,
  lm_basisph,
  lm_pterygoid,
  lm_palatine
)

module_names <- list(
  "premax_d",
  "max_d",
  "frontal",
  "nasal", 
  "jugal",
  "parietal",
  "premax_v",
  "max_v",
  "squam_v",
  "glenoid",
  "squam_z",
  "supraocc",
  "occ_cond",
  "basiocc",
  "basisph",
  "pterygoid",
  "palatine"
  )
#First, run same tree from 75_80 time bin, tree #95 (randomly generated) for all modules and for whole skull

#trees75_80_95
tree75_80<-read.nexus(file = "./Analyses/Bayes_Traits/final_BT_analyses/tree1_75_80_t95_all/tree75_80_95.nex")

btfolder <- "./Analyses/Bayes_Traits/final_BT_analyses/tree1_75_80_t95_all/pcScores2/" 

#this is a little function to write the bayestraits control files with informed priors. 
#x is the pc scores, module name is a character string for naming the file, folder is the destination folder
#for all of them, we are using the following options:
#7- independent constrast
#2- MCMC
#stones 500 5000
#iterations 750,000,000
#cores 8
#burnin 350,000,000
#sample 2500
#PriorAll gamma 5 12 <<<<------this sets the prior for the rate (sigma) and ancestral state (alpha) for each trait to a gamma distribution with k=5 and mu=12

alpha = 5
beta = 12
mean = alpha * beta
std = sqrt(alpha*beta^2)
print(mean); print(std) 
range = seq(0, mean + 4*std, 0.02)
gy1 = dgamma(range, alpha, rate = 1/beta)
plot(range, gy1, type ="l")
#this is a good prior for rate and a bad prior for anc state
#the function then sets the prior for alpha to be a uniform distribution from the min to the max of that PC score (trait)

list_of_models <- c("", "OU", "lambda", "kappa", "delta")

write_cmd_file<-function(x, module_name, varrate=TRUE, model = "", folder){
  filename_cmd<-paste0(folder,module_name,"_",model,"_control.cmd") 
  sink(filename_cmd)
  cat("7","2","stones 500 5000","iterations 750000000","cores 8","burnin 350000000","sample 2500","PriorAll gamma 5 12",sep="\n")
  cat("\n")
  for (j in 1:ncol(x)){
    cat(paste0("Prior Alpha-",j," uniform ",signif(min(x[,j]),5)," ",signif(max(x[,j]),5)))
    cat("\n")
  }
  cat(model)
  cat("\n")
  if(varrate==TRUE){
    cat("varrates")
    cat("\n")
  }
  cat("run")
  sink()
}

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end


phypc <- gm.prcomp(A = shape.data, phy = tree75_80, GLS=TRUE)
for (i in 1:length(runs)){
      #keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
      scores <- phypc$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.90)[1])]*1000 
      write.table(scores, file = paste0(btfolder,"/Phylo_PC_SCORES_","tree75_80_95","_",list_of_models[ii],"_",runs[i],"_run_",".txt"), quote = FALSE, col.names=FALSE)
      write_cmd_file(x=scores, module_name = paste0("tree75_80_95_",runs[i]), varrate=TRUE, model = list_of_models[ii], folder=btfolder)
      write_cmd_file(x=scores, module_name = paste0("tree75_80_95_",runs[i]), varrate=FALSE, model = "", folder=btfolder)
      }

#nasal
for (i in 1:length(lm_list)){
  phypc <- gm.prcomp(A = shape.data[lm_list[[i]],,], phy = tree75_80, GLS=TRUE)
  #keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
  scores <- phypc$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.90)[1])]*1000 
  for (ii in 1:length(runs)){
    lapply(1:length(list_of_models), function(j) write.table(scores, 
                                                             file = paste0(btfolder,"/Phylo_PC_SCORES_","tree75_80_95",module_names[i], "_", list_of_models[j],"_",runs[ii],"_run",".txt"),
                                                             quote = FALSE, col.names=FALSE))
    lapply(1:length(list_of_models), function(j) write_cmd_file(x=scores, module_name = paste0(module_names[[i]], "_",runs[ii]), model = list_of_models[j], varrate=TRUE,folder=btfolder))
    write_cmd_file(x=scores, module_name = paste0(module_names[[i]],"_single_rate_",runs[ii]), model = "", varrate=FALSE, folder=btfolder)
    write.table(scores, 
                file = paste0(btfolder,"/Phylo_PC_SCORES_","tree75_80_95",module_names[i], "_", "single_rate","_",runs[ii],"_run",".txt"),
                quote = FALSE, col.names=FALSE)
  }
}

#Second, run one tree from each of the 6 time bins for the 3 different topologies for the whole skull dataset.  If there are differences, repeat for each module.  Use a different random tree for each time bin.

btfolder <- "./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees/" 

#trees70_75_14
tree70_75<-read.nexus(file = "./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees/tree70_75_14.nex")

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end

for (i in 1:length(runs)){
  phypc <- gm.prcomp(A = shape.data, phy = tree70_75, GLS=TRUE)
  #keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
  scores <- phypc$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.90)[1])]*1000 
  write.table(scores, file = paste0(btfolder,"/Phylo_PC_SCORES_","tree70_75_14","_", list_of_models[ii],"_",runs[i],"_run",".txt"), quote = FALSE, col.names=FALSE)
  write_cmd_file(x=scores, module_name = paste0("tree70_75_14_",runs[i]), model = list_of_models[ii], folder=btfolder)
  
  }

#trees80_85_85
tree80_85<-read.nexus(file = "./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees/tree80_85_85.nex")

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end

for (i in 1:length(runs)){
  phypc <- gm.prcomp(A = shape.data, phy = tree80_85, GLS=TRUE)
  #keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
  scores <- phypc$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.90)[1])]*1000 
  write.table(scores, file = paste0(btfolder,"/Phylo_PC_SCORES_","tree80_85_85","_run_",runs[i],".txt"), quote = FALSE, col.names=FALSE)
  write_cmd_file(x=scores, module_name = paste0("tree80_85_85_",runs[i]), folder=btfolder)
  }


#trees85_90_83
tree85_90<-read.nexus(file = "./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees/tree85_90_83.nex")

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end

for (i in 1:length(runs)){
  phypc <- gm.prcomp(A = shape.data, phy = tree85_90, GLS=TRUE)
  #keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
  scores <- phypc$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.90)[1])]*1000 
  #write.table(scores, file = paste0(btfolder,"/Phylo_PC_SCORES_","tree85_90_83","_run_",runs[i],".txt"), quote = FALSE, col.names=FALSE)
  write_cmd_file(x=scores, module_name = paste0("tree85_90_83_",runs[i]), folder=btfolder)
}



#trees90_95_38
tree90_95<-read.nexus(file = "./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees/tree90_95_38.nex")

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end

for (i in 1:length(runs)){
  phypc <- gm.prcomp(A = shape.data, phy = tree90_95, GLS=TRUE)
  #keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
  scores <- phypc$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.90)[1])]*1000 
  #write.table(scores, file = paste0(btfolder,"/Phylo_PC_SCORES_","tree90_95_38","_run_",runs[i],".txt"), quote = FALSE, col.names=FALSE)
  write_cmd_file(x=scores, module_name = paste0("tree90_95_38_",runs[i]), folder=btfolder)
  
  }


#trees95_100_49
tree95_100<-read.nexus(file = "./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees/tree95_100_49.nex")

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end

for (i in 1:length(runs)){
  phypc <- gm.prcomp(A = shape.data, phy = tree95_100, GLS=TRUE)
  #keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
  scores <- phypc$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.90)[1])]*1000 
  #write.table(scores, file = paste0(btfolder,"/Phylo_PC_SCORES_","tree95_100_49","_run_",runs[i],".txt"), quote = FALSE, col.names=FALSE)
  write_cmd_file(x=scores, module_name = paste0("tree95_100_49_",runs[i]), folder=btfolder)
  }
 
