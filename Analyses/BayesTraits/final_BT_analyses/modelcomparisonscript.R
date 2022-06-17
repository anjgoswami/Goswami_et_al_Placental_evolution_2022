library(tidyverse)


setwd("/mnt/shared/home/agoswami/projects/nhm/goswamilab/BT_placental/whole_skull_tree1_all_trees_322/pPCscores")

#get list of output folders
folderlist <- dir()


stones_list <- lapply(folderlist, function(x) dir(x, pattern="*a_run_.txt.Stones.txt$"))
stones_list <- unlist(stones_list)
####### NEED A LINE HERE TO FILTER SsINGLE RATE ONLY

modellist <- c(
  "kappa",
  "lambda",
  "delta",
  "OU",
  "__")

var_or_single <- c("single_rate", "var_rate")

stonestable <- tibble(file1 = stones_list)
stonestable <- stonestable %>% mutate(Tree = str_extract(file1, paste0(folderlist, collapse = "|"))) %>%
  mutate(Model = str_extract(file1, paste0(modellist, collapse = "|"))) %>%
  mutate(VarOrSingle = str_extract(file1, paste0(var_or_single, collapse = "|"))) %>%
  mutate(Model=replace(Model, Model=="__", "BM")) %>%
  mutate(file1 = paste0("./", Tree, "/", file1))
stonestable <- stonestable %>% mutate(cleanratelabel = recode(stonestable$VarOrSingle, single_rate = "single rate", var_rate = "variable rates"))
stonestable <- stonestable %>% mutate(Model = recode(stonestable$Model, delta = "Delta", lambda = "Lambda", kappa = "Kappa"))



getStones <- function(stones, labels = NULL) {
  res <- matrix(ncol = 2, nrow = length(stones))
  if (!is.null(labels)) {
    if (length(labels) != length(stones)) {
      stop("The number of labels must equal the number of stones files.")
    }
  }
  colnames(res) <- c("logfile", "marginalLh")
  for (i in 1:length(stones)) {
    raw <- readLines(stones[[i]])
    res[i, 1] <- stones[[i]]
    res[i, 2] <- as.numeric(strsplit(raw[length(raw)], "\t")[[1]][2])
  }
  res <- data.frame(res, stringsAsFactors = FALSE)
  if (!is.null(labels)) {
    res$logfile <- labels
  }
  res$marginalLh <- as.numeric(as.character(res$marginalLh))
  res <- tibble::as_tibble(res)
  class(res) <- append("bt_stones", class(res))
  return(res)
}


modellik <- getStones(stonestable$file1)
stonestable <- stonestable %>% mutate(marginalLh = modellik$marginalLh) %>% group_by(Tree)
bestmodels <- stonestable %>% filter(marginalLh == max(marginalLh))

write.csv(stonestable, "/mnt/shared/home/rfelice/projects/nhm/goswamilab/BT_placental/whole_skull_tree1_all_trees_322/full_stones_table.csv")
write.csv(bestmodels, "/mnt/shared/home/rfelice/projects/nhm/goswamilab/BT_placental/whole_skull_tree1_all_trees_322/best_stones_table.csv")



library(reshape2)
#plot:
autoplot.bt_stones <- function(x) {
  d <- matrix(ncol = nrow(x), nrow = nrow(x))
  colnames(d) <- rownames(d) <- paste(x$Model, x$cleanratelabel, sep=", ")

  for (i in seq_len(nrow(d))) {
    for (j in seq_len(ncol(d))) {
      d[i, j] <- round(2 * (x$marginalLh[i] - x$marginalLh[j]), 2)
    }
  }

  d <- reshape2::melt(d[nrow(d):1, ])
  d$bf_cat <- cut(d$value,
    breaks = c(min(d$value) - 1, 0, 2, 5, 10, max(d$value)),
    labels = c("<0", "0-2", "2-5", "5-10", ">10"))
  d$value[d$value == 0] <- d$bf_cat[d$value == 0] <- NA
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Var2, y = Var1, fill = bf_cat)) +
    ggplot2::geom_tile(colour = "white", size = 0.25) +
    ggplot2::geom_text(data = d, ggplot2::aes(Var2, Var1, label = value),
      na.rm = TRUE) +
    ggplot2::scale_y_discrete(expand=c(0,0))+
    ggplot2::scale_x_discrete(expand=c(0,0), position = "top") +
    ggplot2::scale_fill_manual(values = viridis::viridis(8)[3:7],
      na.value = "grey90", drop = FALSE) +
    ggplot2::labs(x = "", y = "", fill = "Bayes factor") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_family = "Helvetica") +
    ggplot2::ggtitle(paste0("Tree Topology 1: Root Age = ", gsub("_","-",x$Tree[1]), " Ma"))
  p
}

plotlist <- stonestable %>% do(plots = autoplot.bt_stones(.))

setwd("/mnt/shared/home/agoswami/projects/nhm/goswamilab/BT_placental/whole_skull_tree1_all_trees_322/plots/")

for (i in 1:nrow(plotlist)){
  filename <- paste0(plotlist$Tree[i], ".pdf")
  pdf(file=filename, width=21, height=27)
  print(plotlist$plots[i])
  dev.off()
}
