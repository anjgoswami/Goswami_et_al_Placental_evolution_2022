## -------------------------------------- ##
##  Rates analysis (BMM model under PL)   ##
## -------------------------------------- ##

library(mvMORPH)
library(dplyr)
library(geiger)

# Retrieve value from slurm job
id_job<-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## -------------------------------------- ##
## Data and trees                         ##
## -------------------------------------- ##

## In this part, you should load the trees and data, and start the reconstruction of the
## ancestral states for the "regimes" on which you want to estimate the rates

# Sample of trees
sample_trees <- read.nexus("./Analyses/mvMorph/t1_list_of_trees.nex")


species_data <- read.csv("./Raw_Data/full_species_data.csv", row.names = "Tip_Label")
fact_eco <- species_data %>% dplyr::select(Fine_Diet)

# Body size?
fact_size <- read.csv("./Data/csize.322.csv", row.names = 1)


# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
# we can use "Morpho" package to do that from an array of superimposed data: data = vecx(gpa$rotated, byrow = TRUE)
dataset <- read.csv("./Data/shape.data.322.csv", row.names = 1)
load("./Data/shape.data.322.R")

moduleLMs <- read.csv('./Raw_Data/lm_ids.csv', row.names = 1)



#starter files for module-specific lm+curve lists
# curve_table <- read_csv("./Raw_Data/placental_curves.csv")
# my_curves <- create_curve_info(curve_table, n_fixed = 66)

lm_premax_d <- which(moduleLMs$Module=="premax_d")
lm_max_d <- which(moduleLMs$Module=="max_d")
lm_frontal <- which(moduleLMs$Module=="frontal")
lm_nasal <- which(moduleLMs$Module=="nasal")
lm_jugal <- which(moduleLMs$Module=="jugal")
lm_parietal <- which(moduleLMs$Module=="parietal")
lm_premax_v <- which(moduleLMs$Module=="premax_v")
lm_max_v <- which(moduleLMs$Module=="max_v")
lm_squam_v <- which(moduleLMs$Module=="squam_v")
lm_glenoid <- which(moduleLMs$Module=="glenoid")
lm_squam_z <- which(moduleLMs$Module=="squam_z")
lm_supraocc <- which(moduleLMs$Module=="supraocc")
lm_occ_cond <- which(moduleLMs$Module=="occ_cond")
lm_basiocc <- which(moduleLMs$Module=="basiocc")
lm_basisph <- which(moduleLMs$Module=="basisph")
lm_pterygoid <- which(moduleLMs$Module=="pterygoid")
lm_palatine <- which(moduleLMs$Module=="palatine")

list_of_modules <- list(
  lm_nasal,
  lm_premax_d,
  lm_max_d,
  lm_jugal,
  lm_frontal,
  lm_parietal,
  lm_squam_v,
  lm_squam_z,
  lm_glenoid,
  lm_supraocc,
  lm_occ_cond,
  lm_basiocc,
  lm_basisph,
  lm_pterygoid,
  lm_palatine,
  lm_max_v,
  lm_premax_v
)

list_of_module_names <- list(
  "nasal",
  "premax_d",
  "max_d",
  "jugal",
  "frontal",
  "parietal",
  "squam_v",
  "squam_z",
  "glenoid",
  "supraocc",
  "occ_cond",
  "basiocc",
  "basisph",
  "pterygoid",
  "palatine",
  "max_v",
  "premax_v"
)

cat_eco = as.factor(fact_eco[sample_trees[[1]]$tip.label,])
names(cat_eco) = sample_trees[[1]]$tip.label

tree_simmap <- make.simmap(sample_trees[[1]], cat_eco , model="ARD", nsim=1) # replace ARD by SYM?


for(j in 1:length(list_of_modules)){
  module.shape.data<-two.d.array(shape.data[list_of_modules[[j]],,])
  
  
  ## -------------------------------------- ##
  ## Perform the analyses                   ##
  ## -------------------------------------- ##
  dat <- list(data=as.matrix(module.shape.data[sample_trees[[1]]$tip.label,]), 
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
  names_id <- paste0("mod_rates_tree1_bin_",id_job,"_",list_of_module_names[j],"_diet_size_BMM_.Rdata")
  save(results, file=names_id)
}
