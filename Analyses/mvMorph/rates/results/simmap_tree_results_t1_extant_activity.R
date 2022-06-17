library(tidyverse)

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_75_80/full/results/diet"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)

tree_extractor_func<-function(x){
  load(x)
  results_tree<-results$tree
  return(results_tree)
}

raw_tree_results <- lapply(resultslist, function(x) tree_extractor_func(x))
tree_results<- raw_tree_results %>% unlist() %>% enframe()

full_tree=read.nexus("./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees_322/tree75_80_95.nex")

plotSimmap(raw_tree_results[[1]], colors=NULL, fsize=0.5, ftype="reg", lwd=2, pts=FALSE, 
           node.numbers=FALSE, mar=NULL, add=FALSE, offset=NULL,
           direction="rightwards", type="fan", setEnv=TRUE, part=1.0, 
           xlim=NULL, ylim=NULL, nodes="intermediate", tips=NULL, maxY=NULL, 
           hold=TRUE, split.vertical=FALSE,
           lend=2, asp=NA, outline=FALSE, plot=TRUE)
XX<-describe.simmap(raw_tree_results[[1]],plot=FALSE)
plot(full_tree,no.margin=TRUE,label.offset=.02,edge.width=2)
nodelabels(pie=XX$states,piecol=c("mediumblue", "green4", "sienna4", "deepskyblue", "darkolivegreen3", "tan3", "red"), cex=0.5)


##### estimating ancestral reconstruction for ecology and devt for frogs
## need to run couple of scripts before this I think to load in/subset data- please see
## Loading_and_subsetting_data

library(ape) #for ancestral states
library(coda) #this is for the mcmc
library(geiger)
library(phytools) #for plotTree
library(mvMORPH)
library(evobiR) 
#install.packages("evobiR")

#setwd("C:/Users/carlb/Box Sync/Carla Bardua")
full_tree=read.nexus("./Analyses/Bayes_Traits/final_BT_analyses/whole_skull_all_trees_322/tree75_80_95.nex")
species.data<-read.csv("./Raw_Data/full_species_data.csv")

tree_rotated=rotateNodes(full_tree, nodes="all") ##
tree_rotated_di <- multi2di(tree_rotated, random = TRUE)
tree_rotated_di$edge.length[which(tree_rotated_di$edge.length == 0)] = 1e-5

##reordering tree and data to be the same names/data match
Reordered <- ReorderData(full_tree, species.data, taxa.names=4)
Ecology=Reordered$Fine_Diet


fitSYMTEST<-ace(Ecology, full_tree, model="ER", type="discrete")

cols<-setNames(c("pink", "red", "green4", "sienna4", "yellow", "deepskyblue", "tan3"),levels(as.factor(Ecology)))

# plot the tree states

plotTree(tree_1_w_data@phylo,type="fan",fsize=0.3,ftype="i",lwd=1)
nodelabels(node=1:full_tree$Nnode+Ntip(full_tree),
           pie=fitSYMTEST$lik.anc, cex=0.2, piecol = cols)+
tiplabels(pie = to.matrix(Ecology, sort(unique(Ecology))), piecol = cols, cex = 0.2)
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1], y = -max(nodeHeights(full_tree)))


##reordering tree and data to be the same names/data match
Reordered <- ReorderData(tree_1_w_data@phylo, species.data, taxa.names="Tip_Label")
Ecology=Reordered$Fine_Diet
fitSYMTEST<-ace(tree_1_w_data@extraInfo$Fine_Diet, tree_1_w_data@phylo, model="ER", type="discrete")
cols<-setNames(c("pink", "red", "green4", "sienna4", "yellow", "deepskyblue", "tan3"),levels(as.factor(Ecology)))
# plot the tree states
plotTree(tree_1_w_data@phylo,type="fan",fsize=0.3,ftype="i",lwd=1)
nodelabels(node=1:tree_1_w_data@phylo$Nnode+Ntip(tree_1_w_data@phylo),
           pie=fitSYMTEST$lik.anc, cex=0.2, piecol = cols)+
  tiplabels(pie = to.matrix(Ecology, sort(unique(Ecology))), piecol = cols, cex = 0.2)
add.simmap.legend(colors = cols, prompt = FALSE, x = 0.9 * par()$usr[1], y = -max(nodeHeights(full_tree)))
