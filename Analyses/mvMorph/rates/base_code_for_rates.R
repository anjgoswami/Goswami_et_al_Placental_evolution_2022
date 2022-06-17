## -------------------------------------- ##
##  Rates analysis (BMM model under PL)   ##
## -------------------------------------- ##

library(mvMORPH)

# Retrieve value from slurm job
id_job<-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


## -------------------------------------- ##
## Data and trees                         ##
## -------------------------------------- ##

## In this part, you should load the trees and data, and start the reconstruction of the
## ancestral states for the "regimes" on which you want to estimate the rates

# Sample of trees
sample_trees <- read.tree("trees75_80subset.tre")

#Pull out the factor to paint the tree with below 
#Diet
# Ecological factor

species_data <- read.csv("./Raw_Data/full_species_data.csv", row.names = "Tip_Label")
fact_diet <- species_data %>% dplyr::select(Fine_diet)

fact_loc <- species_data %>% dplyr::select(Broad_Locomotion)

# Body size?
fact_size <- read.csv("./Data/csize.322.csv", row.names = 1)

dat <- list(data=as.matrix(dataset[sample_trees[[1]]$tip.label,]),
            diet=as.factor(fact_diet[sample_trees[[1]]$tip.label,]),
            loc=as.factor(fact_loc[sample_trees[[1]]$tip.label,]),
            size=as.numeric(log(fact_size[sample_trees[[1]]$tip.label,])))

table(dat$diet,dat$loc)


# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
# we can use "Morpho" package to do that from an array of superimposed data: data = vecx(gpa$rotated, byrow = TRUE)
dataset <- read.csv("./Data/shape.data.322.csv", row.names = 1)


## -------------------------------------- ##
## Perform the ancestral reconstruction   ##
## -------------------------------------- ##

# We need to reconstruct first the "regimes" on which evolutionary rates will be estimated
# build a tree of class "SIMMAP" from phytools => three options here:
# i) if you run your analyses on a set of trees, maybe we can simply use a stochastic map each time. 
# Averaged over the trees this should represent uncertainty in both the trees and the reconstructions
# ii) instead, you can use the ML reconstruction for the ancestral states and build a tree to the simmap format.
# that is, we're only using the "format" of simmap trees but not actually using the stochastic mapping technique to reconstruct the ancestral states.
# iii) use the averaged reconstruction across multiple stochastic mapping (=> this should be pretty similar to ML)

cat_eco = as.factor(fact_diet[sample_trees[[1]]$tip.label,])
names(cat_eco) = sample_trees[[1]]$tip.label

# i)
tree_simmap <- make.simmap(sample_trees[[1]], cat_eco , model="ARD", nsim=1) # replace ARD by SYM?

# # or ii) use the ML ancestral state reconstruction for the ecological factor.
# # you need the paintAllTree function:
# # tree = the phylogenetic tree
# # ancestral = either a "describe.simmap" object from phytools or an "ace" object from ape
# # tips = states at the tips. Check that the tree and data are in same order!
# 
# paintAllTree <- function(tree, ancestral, tips){  
#   
#   if(inherits(ancestral, "describe.simmap")){
#     names_grps <- colnames(ancestral$ace)
#     statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
#   }else{
#     names_grps <- colnames(ancestral$lik.anc)
#     statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
#   }  
#   
#   combined = as.character(c(tips, statesNodes))
#   treebis=tree
#   for(i in sort(tree$edge[,2])){
#     treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
#   }  
#   return(treebis)
# }
# 
# # First estimate ancestral state using ML in "ape"
# ace_habitat <- ace(cat_eco, sample_trees[[id_job]], type = "discrete", model="ARD")
# # convert the ML reconstruction to a "simmap" tree object
# tree_simmap <- paintAllTree(sample_trees[[id_job]], ace_habitat, as.character(cat_eco))

# # iii) marginal reconstructions by stochastic mapping
# # i.e. do the same using stochastic mapping instead of ML
# nsim = 100 # Number of stochastic mappings
# my_trees<-make.simmap(sample_trees[[id_job]], cat_eco , model="ARD", nsim=nsim)
# 
# # Use the collection of stochastic mappings to reconstruct the marginal ancestral states
# distrib <-summary(my_trees, plot=TRUE)
# # Use these reconstruction to build a simmap tree
# tree_simmap  <- paintAllTree(sample_trees[[id_job]], distrib, as.character(cat_eco))




## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##

## !! check that the object names match the ones you used above

# Prepare data for analyses in a list
#dat <- list(data=as.matrix(dataset[sample_trees[[id_job]]$tip.label,]))

# If you want to account for size, rather than removing size first, you can use the following instead
dat <- list(data=as.matrix(dataset[sample_trees[[1]]$tip.label,]), 
            size=as.numeric(fact_size[sample_trees[[1]]$tip.label,]))

# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
fit <- mvgls(data~1, data=dat, tree=tree_simmap, model="BMM", method = "H&L", error=TRUE)

# if you want to add a covariate, e.g. size, then use the following instead:
# fit <- mvgls(data~size, data=dat, tree=tree_simmap, model="BMM", method = "PL", error=TRUE)


# Save the results.
# Multiple options here again. If you just want the rates (and not say, the evolutionary covariances), you can save a lot of memory space by saving only that:

#name_rates_file <- paste("estimated_rates_tree_",id_job,".txt",sep="")
#write.table(fit$param, file=name_rates_file)

# # or save a bit more of data in a list, like the reconstructed simmap tree and rates, or the complete object!
results <- list(fit=fit$param, tree=tree_simmap)
names_id <- paste("tree1_75_80_",id_job,"diet_size_BMM_.Rdata", sep = "")
save(results, file=names_id)
