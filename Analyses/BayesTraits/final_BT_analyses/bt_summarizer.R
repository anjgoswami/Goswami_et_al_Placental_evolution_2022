library(BTRTools)
library(tidyverse)
library(treeio)

setwd("/home/rfelice/projects/nhm/goswamilab/BT_placental/whole_skull_all_trees/")

add_rjpp_to_tree <- function(rjpp_out){
    rjpp_data <- as_tibble(rjpp_out$data)
    timetree <- rjpp_out$meantree
    timetree$edge.length <- rjpp_data$orgBL[-1]
    timetree$root.time <- max(nodeHeights(timetree))
    rjpp_data_nodes <- rjpp_data %>% rename(., node=descNode) %>% mutate(., iters = rjpp_out$niter) %>% mutate(., ppRate = round(nOrgnNRate/iters,2))
    timetree <- treeio::as.treedata(timetree)
    treedata <- treeio::full_join(timetree, y = rjpp_data_nodes, by = "node")
    return(treedata)
}



tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/Phylo_PC_SCORES_tree75_80_95_skull_run_a.txt.VarRates.txt",
                               rjtrees = "./pPCscores/Phylo_PC_SCORES_tree75_80_95_skull_run_a.txt.Output.trees",
                               tree = "./trees/tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(tree_1_BTraits)

save(tree_results, "./pPCscores/tree_results_75_80.Rdata")

rm(tree_1_BTraits)

tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/Phylo_PC_SCORES_tree70_75_14_skull_run_a.txt.VarRates.txt",
                               rjtrees = "./pPCscores/Phylo_PC_SCORES_tree70_75_14_skull_run_a.txt.Output.trees",
                               tree = "./trees/tree70_75_14.nex")

tree_results<-add_rjpp_to_tree(tree_1_BTraits)

save(tree_results, "./pPCscores/tree_results_70_75.Rdata")

rm(tree_1_BTraits)

tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/Phylo_PC_SCORES_tree80_85_85_skull_run_a.txt.VarRates.txt",
                               rjtrees = "./pPCscores/Phylo_PC_SCORES_tree80_85_85_skull_run_a.txt.Output.trees",
                               tree = "./trees/tree80_85_85.nex")

tree_results<-add_rjpp_to_tree(tree_1_BTraits)

save(tree_results, "./pPCscores/tree_results_80_85.Rdata")

rm(tree_1_BTraits)

tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/Phylo_PC_SCORES_tree85_90_83_skull_run_a.txt.VarRates.txt",
                               rjtrees = "./pPCscores/Phylo_PC_SCORES_tree85_90_83_skull_run_a.txt.Output.trees",
                               tree = "./trees/tree85_90_83.nex")

tree_results<-add_rjpp_to_tree(tree_1_BTraits)

save(tree_results, "./pPCscores/tree_results_85_90.Rdata")

rm(tree_1_BTraits)

tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/Phylo_PC_SCORES_tree90_95_38_skull_run_a.txt.VarRates.txt",
                               rjtrees = "./pPCscores/Phylo_PC_SCORES_tree90_95_38_skull_run_a.txt.Output.trees",
                               tree = "./trees/tree90_95_38.nex")

tree_results<-add_rjpp_to_tree(tree_1_BTraits)

save(tree_results, "./pPCscores/tree_results_90_95.Rdata")

rm(tree_1_BTraits)

tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/Phylo_PC_SCORES_tree95_100_49_skull_run_a.txt.VarRates.txt",
                               rjtrees = "./pPCscores/Phylo_PC_SCORES_tree95_100_49_skull_run_a.txt.Output.trees",
                               tree = "./trees/tree95_100_49.nex")

tree_results<-add_rjpp_to_tree(tree_1_BTraits)

save(tree_results, "./pPCscores/tree_results_95_100.Rdata")

rm(tree_1_BTraits)
