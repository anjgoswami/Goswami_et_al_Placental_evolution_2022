library(dplyr)
library(mvMORPH)
library(geiger)

# Retrieve value from slurm job
id_job<-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


## -------------------------------------- ##
## Data and trees                         ##
## -------------------------------------- ##

## In this part, you should load the trees and data, and start the reconstruction of the
## ancestral states for the "regimes" on which you want to estimate the rates

# Sample of trees
sample_trees <- read.tree("./Raw_Data/trees/updated_tree1_subsets/trees70_75subset.tre")

#Pull out the factor to paint the tree with below 
#Eco
species.data <- read.csv("./Raw_Data/extant_species_data.csv", row.names = "Tip_Label")
fact_eco <- species.data %>% dplyr::select(Development)

# Body size?
fact_size <- read.csv("./Data/extant.size.csv", row.names = 1)


# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
dataset <- read.csv("./Data/extant.data.csv", row.names = 1)

nc1<-name.check(sample_trees[[1]],dataset)
sample_trees[[1]]<-drop.tip(sample_trees[[1]], nc1$tree_not_data)

## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##

cat_eco = as.factor(fact_eco[sample_trees[[1]]$tip.label,])
names(cat_eco) = sample_trees[[1]]$tip.label

tree_simmap <- make.simmap(sample_trees[[1]], cat_eco , model="ARD", nsim=1) # replace ARD by SYM?


## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##

## !! check that the object names match the ones you used above

# Prepare data for analyses in a list

dat <- list(data=as.matrix(dataset[sample_trees[[1]]$tip.label,]), 
            size=as.numeric(fact_size[sample_trees[[1]]$tip.label,]))

# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)


# Save the results.

# save a bit more of data in a list, like the reconstructed simmap tree and rates, or the complete object

results <- list(fit=fit$param, tree=tree_simmap)
names_id <- paste("./results/dev/tree1_70_75_",id_job,"dev_size_BMM_.Rdata", sep = "")
save(results, file=names_id)

treefolder <- "./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees_322/" 
treefiles <- dir(treefolder)
mytrees <- lapply(1:6, function(x) read.nexus(paste0(treefolder,treefiles[x])))
