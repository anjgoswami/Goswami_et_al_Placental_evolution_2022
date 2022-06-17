library(devtools)
library(tidyverse)
install_github("JClavel/mvMORPH", ref="devel_1.0.4")
library(mvMORPH)
library(SURGE)


# phylogenetic tree
tree1_75_80_95 <- read.nexus("./Raw_Data/trees/tree1_random_set/tree75_80_95.nex")
#All of my species data is already rearranged to the same order as the tree
species.data<-read.csv("./Raw_Data/full_species_data.csv")
#Pull out the factor to paint the tree with below 

# Diet factor
fact_diet <- read.csv("./Raw_Data/fine_diet.csv", row.names = 1)

# locomotion factor
fact_loc <- read.csv("./Raw_Data/broad_locomotion.csv", row.names = 1)

# Body size
fact_size <- read.csv("./Data/csize.csv", row.names = 1)

# Morphology > a 2D matrix [species x coordinates]
data <- read.csv("./Data/shape.data.csv", row.names = 1)

# #Add the function to reconstruct from a sample of trees 
# paintAllTree <- function(tree, ancestral, tips){  
#   if(inherits(ancestral, "describe.simmap")){
#     names_grps <- colnames(ancestral$ace)
#     statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
#   }else{
#     names_grps <- colnames(ancestral$lik.anc)
#     statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
#   }  
#   combined = as.character(c(tips, statesNodes))
#   treebis=tree
#   for(i in sort(tree$edge[,2])){
#     treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
#   }  
#   return(treebis)
# }

indices <- matrix(1:ncol(data), ncol = 3, byrow = TRUE)

###########split out regions
curve_table <- read.csv("./Raw_Data/placental_curves.csv")
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

# list of samples
samples <- list(nasal = lm_nasal, premax_d = lm_premax_d, max_d = lm_max_d, jugal = lm_jugal, frontal = lm_frontal, parietal = lm_parietal, squam_v = lm_squam_v, squam_z = lm_squam_z, glenoid = lm_glenoid, suprocc = lm_supraocc, occ_cond = lm_occ_cond, basiocc = lm_basiocc, basisph = lm_basisph, pterygoid = lm_pterygoid, palatine = lm_palatine, max_v = lm_max_v, premax_v = lm_premax_v)



# build a SIMMAP tree > the idea is to reconstruct the history of the discrete states for which we test differences in rates. Ideally this should be performed on a sample of stochastic maps. Or by using a mapping constructed from the ML reconstruction for instance. Here I use a single map for illustration purpose
cat_diet = as.factor(fact_diet[tree1_75_80_95$tip.label,])
names(cat_diet) = tree1_75_80_95$tip.label
tree_simmap <- make.simmap(tree1_75_80_95, cat_diet , model="ARD", nsim=1) 



#Ideally this should be performed on a sample of stochastic maps.
#Or by using a mapping constructed from the ML reconstruction for instance
# simm_tr<- make.simmap(treeWHOLE, model="ARD", cat_feeding, nsim=100)
# a_single_simmap_tree <- paintAllTree(treeWHOLE, describe.simmap(simm_tr), as.character(cat_feeding))
# subsets


list_of_results<-list()
for(i in 1:length(samples)){
  subsamples = as.numeric(t(indices[samples[[i]], ]))
  subdat <- data[,subsamples]
  # Loop over the distribution of phylogenetic trees
  # Prepare data for analyses
  dat <- list(Y=as.matrix(subdat[tree1_75_80_95$tip.label,]),
              size=as.numeric(fact_size[tree1_75_80_95$tip.label,]))
  # fit linear model > here I include size in the model as a covariate. This should be equivalent to estimating rates on size corrected dat
  fit <- mvgls(Y~size, data=dat, tree=tree_simmap, model="BMM", method = "PL", error=TRUE)
  list_of_results[[i]]<-fit
}
names(list_of_results)<-names(samples)
save(list_of_results, file = "list_of_results_diet_abs.Rdata")

###repeating for locomotion

# build a SIMMAP tree > the idea is to reconstruct the history of the discrete states for which we test differences in rates. Ideally this should be performed on a sample of stochastic maps. Or by using a mapping constructed from the ML reconstruction for instance. Here I use a single map for illustration purpose
cat_loc = as.factor(fact_loc[tree1_75_80_95$tip.label,])
names(cat_loc) = tree1_75_80_95$tip.label
tree_simmap_loc <- make.simmap(tree1_75_80_95, cat_loc , model="ARD", nsim=1) 

list_of_results<-list()
for(i in 1:length(samples)){
  subsamples = as.numeric(t(indices[samples[[i]], ]))
  subdat <- data[,subsamples]
  # Loop over the distribution of phylogenetic trees
  # Prepare data for analyses
  dat <- list(Y=as.matrix(subdat[tree1_75_80_95$tip.label,]),
              size=as.numeric(fact_size[tree1_75_80_95$tip.label,]))
  # fit linear model > here I include size in the model as a covariate. This should be equivalent to estimating rates on size corrected dat
  fit <- mvgls(Y~size, data=dat, tree=tree_simmap_loc, model="BMM", method = "PL", error=TRUE)
  list_of_results[[i]]<-fit
}

names(list_of_results)<-names(samples)
save(list_of_results, file = "./Analyses/mvMorph/rates/list_of_results_loc.Rdata")


#for a list of trees
# trees1_75_80 <- list()
# filenames = dir(path = "C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/BayesTraits/BT_input/tree1/tree1_75_80/", pattern = "*.nexus")
# for (i in filenames){
#   dat <- read.nexus(i)  dat <- read.nexus(i)
#   trees1_75_80[[i]]<-dat
# }

trees1_75_80<-read.tree(file = "./Raw_Data/trees/dated_tree1_subsets/trees70_75subset.tre")

my_list_of_trees <- lapply(trees1_75_80, function(a_tree_to_transform) make.simmap(a_tree_to_transform, model="ARD", as.factor(cat_diet), nsim=1))
results_diet_75_80_all<-lapply(my_list_of_trees, function(trees_sim) mvgls(Y~size, data = dat, tree=trees_sim, model="BMM", error=TRUE)$param) # if you want to save only the rates# results_diet_75_80_all<-lapply(my_list_of_trees, function(trees_sim) mvgls(Y~size, data = dat, tree=trees_sim, model="BMM", error=TRUE)$param) # if you want to save only the rates# results_diet_75_80_all<-lapply(my_list_of_trees, function(trees_sim) mvgls(Y~size, data = dat, tree=trees_sim, model="BMM", error=TRUE)$param) # if you want to save only the rates# results_diet_75_80_all<-lapply(my_list_of_trees, function(trees_sim) mvgls(Y~size, data = dat, tree=trees_sim, model="BMM", error=TRUE)$param) # if you want to save only the rates
save(results_diet_75_80_all, file = "./Analyses/mvMorph/rates/list_of_results_diet_75_80_all.Rdata")
