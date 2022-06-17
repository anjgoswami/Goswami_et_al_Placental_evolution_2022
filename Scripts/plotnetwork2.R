plotNetwork2<-function (rhos, module_names = NULL, linecolour = "#56B4E9", 
          title = NULL, layout = NULL, minimum = 0) 
{
  withins <- grep("Module*", colnames(rhos))
  nmodule <- length((grep("Module*", colnames(rhos))))
  rholist <- t(rhos)
  words <- strsplit(rownames(rholist), " ")
  plotcorr <- matrix(data = NA, nrow = nmodule, ncol = nmodule)
  modnums <- unlist(lapply(withins, function(x) strsplit(colnames(rhos)[x], 
                                                         " ")[[1]][2]))
  mods <- sapply(words[withins], function(x) paste(x, collapse = " "))
  mods <- gsub("Module", "M", mods)
  mods <- mods[order(sapply(mods, function(x) as.numeric(strsplit(x, 
                                                                  "M ")[[1]][[2]])))]
  colnames(plotcorr) <- rownames(plotcorr) <- mods
  for (i in 1:(length(words) - 1)) {
    if (length(words[[i]]) == 2) {
      module <- paste("M", words[[i]][2])
      plotcorr[module, module] <- rholist[i, "MaxL_p"]
    }
    if (length(words[[i]]) == 3) {
      from_module <- paste("M", (words[[i]][1]))
      to_module <- paste("M", words[[i]][3])
      plotcorr[from_module, to_module] <- rholist[i, "MaxL_p"]
      plotcorr[to_module, from_module] <- rholist[i, "MaxL_p"]
    }
  }
  within <- diag(plotcorr)
  between <- plotcorr
  if (is.null(module_names)) {
    mod.names <- mods
  }
  else if (module_names == "rhos") {
    mod.names <- within
  }
  if (linecolour == "viridis") {
    cls <- viridis::viridis(100, direction = -1)
    vcols <- apply(between * 100, 1, function(x) cls[x])
    linecolour <- NULL
  }
  else {
    vcols <- NULL
  }
  qgraph::qgraph(between, shape = "circle", posCol = linecolour, 
                 edge.color = vcols, labels = mod.names, vsize = within * 
                   10, diag = FALSE, title = title, layout = layout, minimum = minimum)
}