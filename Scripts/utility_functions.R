####Load posterior to test convergence####
tracePlots <- function(file, burnin=0, thinning=1, plot=TRUE, display=c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")){
  require(BTRTools)
  require(coda)
  
  rjout <- loadRJ(file, burnin = burnin, thinning = thinning)
  chain_out <- type.convert(rjout$rj_output)
  rownames(chain_out) = chain_out[,"It"]
  chain_out = chain_out[,-1]
  # Retrieve numerical
  index <- sapply(chain_out,function(x) is.numeric(x))
  chain <- mcmc(chain_out[,index])
  
  # plot the trace
  if(plot){
    plot(chain[,display])
  }
  
  # Just compute some statistics (autocorrelation...)
  cat("Effective sample size:","\n")
  print(effectiveSize(chain[,display]))
  
  # return results
  invisible(chain)
  
}


####function to obtain post probabilies of rate shifts on node####
return_pprob <- function(PP, threshold=0.5, cl = "nOrgnNRate"){
  nodes <- PP$data$descNode[which((PP$data[ , cl] / PP$niter) >= threshold)]
  pprobs <- PP$data[which((PP$data[ , cl] / PP$niter) >= threshold) , cl] / PP$niter
  
  result = list(nodes=nodes, pprobs=pprobs)
  return(result)
}

#uses treeio to combine tree topology and bayestraits posterior into single S4 object for plotting.
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


##############rate distributions by clade##########

rate_summary <- function(rjpp_out, time_tree, lookup.table, group_column, taxa, min_clade_size=2){
  #select the species from  the lookup table that are part of the current dataset
  lookup.table <- rename(lookup.table, GROUP = group_column)
  lookup.table <- rename(lookup.table, SPECIES = taxa)
  current.species <- filter(lookup.table, SPECIES %in% rjpp_out$meantree$tip.label)
  clades <- unique(lookup.table$GROUP)
  #select the sub tree for each group
  time_trees <- lapply(1:length(clades), function(x) keep.tip(time_tree, (filter(current.species, GROUP == clades[x]) %>% pull(SPECIES))))
  names(time_trees)<-clades
  taxa_per_clade <- unlist(lapply(1:length(time_trees), function(x) Ntip(time_trees[[x]])))
  time_trees <- time_trees[which(taxa_per_clade>=min_clade_size)]
  #get the time represented by each sub tree
  time_per_tree <- unlist(lapply(1:length(time_trees), function(x) max(nodeHeights(time_trees[[x]]))))
  #the first column is the rate on the root and our time tree has no root edge so this is removed with the "[-1,]"
  rate_table <- rjpp_out$scalars$rates[-1,]
  #make copy the mean tree x times where x is the number of trees in the posterior distribution
  treelist<-rep(list(rjpp_out$meantree),dim(rate_table)[2])
  #now convert the edge lengths of each of those trees to be the rate scalar for that edge
  #this allows us to keep rates associated with the appropriate edges as we go to the next step
  for (k in 1:length(treelist)){
    treelist[[k]]$edge.length<-rate_table[,k]
  }
  #now subset those trees based on group identity
  clades <- names(time_trees)
  results_table <- as_tibble(t(setNames(as.numeric(rep("", length(clades))), clades))[0, ])
  results_table_scaled <- results_table 
  for (i in 1:length(treelist)){
    scaled_trees <- lapply(1:length(clades), function(x) keep.tip(treelist[[i]], time_trees[[x]]$tip.label))
    mean_rate_list <- unlist(lapply(1:length(scaled_trees), function(x) mean(scaled_trees[[x]]$edge.length)))
    mean_rate_list_scaled <- mean_rate_list/time_per_tree
    names(mean_rate_list) <- names(mean_rate_list_scaled) <- colnames(results_table)
    results_table <- bind_rows(results_table, mean_rate_list)
    results_table_scaled <- bind_rows(results_table_scaled, mean_rate_list_scaled)
  }
  #comnbine into two tables, one for unscaled one for scaled 
  results <- results_table %>% pivot_longer(cols = clades, 
                                            names_to = "Group",
                                            values_to = "MeanRate")
  results_scaled <- results_table_scaled %>% pivot_longer(cols = clades, 
                                                          names_to = "Group",
                                                          values_to = "MeanRate")
  list1<-list(results, results_scaled)
  names(list1)<-c("results","results_scaled")
  return(list1)
}