library(BTRTools)
library(ape)
library(dplyr)
library(treeio)


setwd("/home/agoswami/projects/nhm/goswamilab/BT_placental/tree1_75_80_95_all/")

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


#basioccipital
bo_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95basiocc__a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95basiocc__a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(bo_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/basiocc_tree_results_75_80.Rdata")

rm(tree_results)

#basisphenoid
bs_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95basisph__a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95basisph__a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(bs_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/basisph_tree_results_75_80.Rdata")

rm(tree_results)

#frontal
fr_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95frontal_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95frontal_single_rate_a_run.txt.Output",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(fr_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/frontal_tree_results_75_80.Rdata")

rm(tree_results)

#glenoid
gl_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95glenoid_single_rate_a_run.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95glenoid_single_rate_a_run.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(gl_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/glenoid_tree_results_75_80.Rdata")

rm(tree_results)

#jugal
jg_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95jugal_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95jugal_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(jg_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/jugal_tree_results_75_80.Rdata")

rm(tree_results)

#max_d
md_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95max_d_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95max_d_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(md_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/max_d_tree_results_75_80.Rdata")

rm(tree_results)

#max_v
mv_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95max_v_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95max_v_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(mv_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/max_v_tree_results_75_80.Rdata")

rm(tree_results)

#nasal
nasal_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95nasal_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95nasal_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(nasal_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/nasal_tree_results_75_80.Rdata")

rm(tree_results)

#occ_cond
occ_cond_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95occ_cond_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95occ_cond_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(occ_cond_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/occ_cond_tree_results_75_80.Rdata")

rm(tree_results)

#palatine
pal_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95palatine_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95palatine_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(pal_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/palatine_tree_results_75_80.Rdata")

rm(tree_results)

#parietal
parietal_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95parietal_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95parietal_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(parietal_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/parietal_tree_results_75_80.Rdata")

rm(tree_results)

#premax_d
md_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95premax_d_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95premax_d_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(premax_d_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/premax_d_tree_results_75_80.Rdata")

rm(tree_results)

#premax_v
premax_v_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95premax_v_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95premax_v_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(premax_v_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/premax_v_tree_results_75_80.Rdata")

rm(tree_results)

#pterygoid
pterygoid_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95pterygoid_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95pterygoid_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(pterygoid_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/pterygoid_tree_results_75_80.Rdata")

rm(tree_results)

#squam_v
squam_v_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95squam_v_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95squam_v_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(squam_v_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/squam_v_tree_results_75_80.Rdata")

rm(tree_results)

#squam_z
squam_z_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_9squam_z_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95squam_z_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(squam_z_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/squam_z_tree_results_75_80.Rdata")

rm(tree_results)

#supraoccipital
supraocc_tree_1_BTraits<-BTRTools::rjpp(rjlog = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95supraocc_single_rate_a_run.txt.VarRates.txt",
                               rjtrees = "./pPCscores/modules_jan/Phylo_PC_SCORES_tree75_80_95supraocc_single_rate_a_run.txt.Output.trees",
                               tree = "./tree75_80_95.nex")

tree_results<-add_rjpp_to_tree(supraocc_tree_1_BTraits)

save(tree_results, "./pPCscores/module_results_BM_var/supraocc_tree_results_75_80.Rdata")

rm(tree_results)

