library(mvMORPH)
library(parallel)

# parametrizations
nperm = 999L # number of permutations. Probably enough...
nbcores = Sys.getenv('SLURM_CPUS_PER_TASK') # retrieve the number of cores specified in the SLURM batch file. This is used in the MANOVA to parallelize the computations

# phylogenetic tree
tree <- read.nexus("subtree_eco_devt.nexus")

# Ecological factor
fact_eco <- read.csv("eco.csv", row.names = 1)

# Development factor
fact_dev <- read.csv("devt.csv", row.names = 1)

# Body size
fact_size <- read.csv("log_size.csv", row.names = 1)

# Morphology > a 2D matrix [species x coordinates]
data <- read.table("data.txt")


# Prepare data for analyses > check that data is in the same order as tips in the tree.
dat <- list(data=as.matrix(data[tree$tip.label,]),
eco=as.factor(fact_eco[tree$tip.label,]),
dev=as.factor(fact_dev[tree$tip.label,]),
size=as.numeric(fact_size[tree$tip.label,]))

# fit linear model > Penalized likelihood with Pagel's lambda
# blas_set_num_threads(6) # for specifying to use multithreading through library(RhpcBLASctl)
fit <- mvgls(data~size*eco*dev, data=dat, tree=tree, model="lambda")

# Multivariate test > Pillai trace
multivariate_test <- manova.gls(fit, nperm=nperm, test="Pillai", nbcores=nbcores, verbose=FALSE, type="II")
print(multivariate_test)

# save results
results <- list(fit=fit, test=multivariate_test)
save(results, file=names_i)
    
# end

