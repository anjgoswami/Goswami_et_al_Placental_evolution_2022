library(mvMORPH)

# phylogenetic tree
tree <- read.nexus("tree75_80_1.nex")

# Diet factor
fact_diet <- read.csv("broad_diet.csv", row.names = 1)

# Morphology
data <- read.csv("shape.lm.csv")

# Body size
fact_size <- read.csv("csize.csv", row.names = 1)


# Prepare data for analyses and check that tips match the data labels
dat <- list(Y=as.matrix(data[tree$tip.label,]),
            size=as.numeric(fact_size[tree$tip.label,]))


# build a SIMMAP tree > the idea is to reconstruct the history of the discrete states for which we test differences in rates. Ideally this should be performed on a sample of stochastic maps. Or by using a mapping constructed from the ML reconstruction for instance. Here I use a single map for illustration purpose
cat_diet = as.factor(fact_diet[tree$tip.label,])
names(cat_diet) = tree$tip.label
tree_simmap <- make.simmap(tree, cat_diet , model="SYM", nsim=1) # replace SYM by ARD?

# fit linear model > here I include size in the model as a covariate. This should be equivalent to estimating rates on size corrected data
fit <- mvgls(Y~size, data=dat, tree=tree_simmap, model="BMM", method = "PL", error=TRUE) # if you remove the error argument, then it is assumed that there's no ME nor intraspecific variance that may blur the results.

# Save results
results <- list(fit=fit$param, tree=tree_simmap)
save(results, file="lm_diet_BMM.Rdata")


