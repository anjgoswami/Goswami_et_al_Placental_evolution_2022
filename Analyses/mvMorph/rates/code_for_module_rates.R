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
fact_eco <- read.csv("./Raw_Data/fine_diet.csv", row.names = 1)

# Body size?
fact_size <- read.csv("./Data/csize.322.csv", row.names = 1)


# Morphology => assumed to be of matrix format [n*p] where n is the number of species and p is the number of dimensions
# we can use "Morpho" package to do that from an array of superimposed data: data = vecx(gpa$rotated, byrow = TRUE)
dataset <- read.csv("./Data/shape.data.322.csv", row.names = 1)

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

list_of_modules <- list(
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

list_of_module_names <- list(
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
names_id <- paste0("tree1_75_80_",id_job,"_",list_of_module_names[j],"_diet_size_BMM_.Rdata")
save(results, file=names_id)
}
