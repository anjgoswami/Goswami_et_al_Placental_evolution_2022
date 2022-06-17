#make a pretty table
#load packages
library(tidyverse)
library(gt)

########tree1#############################
#load data:
setwd("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/pmanovas/pmanovas_fullset_322/tree1/results/extant/")
filenames<-filenames<-dir(pattern=".Rdata")
treenames<-substr(filenames, start = 1, stop = 7)
data<-tibble() 
for (i in 1:length(filenames)) {
  load(filenames[[i]])
  tree <- paste0(treenames[[i]])
  pillai <- results$test$stat
  terms <- results$test$terms
    p_val <- results$test$pvalue
  ses <- (log(results$test$stat) - colMeans(log(results$test$nullstat))) / apply(log(results$test$nullstat), 2, sd)
  datmat<-rbind(pillai,ses,p_val)
  colnames(datmat)<- c("Size", "Diet", "Locomotion", "Habitat", "Development", "Activity", "Social", "Size:Diet", "Size:Loc", "Size:Hab","Size:Act", "Size:Dev", "Size:Soc")
  data_tmp <- as_tibble(datmat)
  data_tmp <- data_tmp %>% mutate(tree=tree) %>% mutate(result = c("Pillai's Test Statistic","SES", "P value"))
  data <- bind_rows(data,data_tmp)
}

data <- data %>% mutate(tree = recode(data$tree, tree1_1 = "Tree Topology 1: Root Age 70-75",
                                      tree1_2 = "Tree Topology 1: Root Age 75-80",
                                      tree1_3 = "Tree Topology 1: Root Age 80-85",
                                      tree1_4 = "Tree Topology 1: Root Age 85-90",
                                      tree1_5 = "Tree Topology 1: Root Age 90-95",
                                      tree1_6 = "Tree Topology 1: Root Age 95-100"))


#data <- rename(data, Size = size, Diet = diet, Loc = loc, Hab = hab, Dev = dev, Act = act, Soc = soc)
#Size:Diet = size:diet, Size:Loc = size:loc, Size:Hab = size:hab, Size:Act = size:act, Size:Dev = size:dev, Size:Soc = = size:soc)

gttest<-data %>%
  gt(rowname_col = "result", groupname_col = "tree")%>%
  fmt_number(columns=everything(),decimals = 3)%>%
  tab_spanner(label = "Factors",columns=everything()) %>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_row_groups(everything()))%>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_body(rows = "P value" < 0.005 ))%>%
  tab_header(title = "Phylogenetic Regressions: Tree Topology 1, Extant Taxa only")%>%
  cols_align(align = "center")%>%
  cols_width(c("Size", "Diet", "Locomotion", "Habitat", "Development", "Activity", "Social", "Size:Diet", "Size:Loc", "Size:Hab","Size:Act", "Size:Dev", "Size:Soc")~px(80))
#  cols_label(size:diet = "Size:Diet", size:loc = "Size:Loc", size:hab = "Size:Hab" , size:act = "Size:Act", size:dev = "Size:Dev", size:soc = "Size:Soc")
gtsave(gttest, file="tree1_extant_table.pdf")

#load data:
setwd("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/pmanovas/pmanovas_fullset_322/tree1/results/full/")
filenames<-filenames<-dir(pattern=".Rdata")
treenames<-substr(filenames, start = 1, stop = 7)
data<-tibble() 
for (i in 1:length(filenames)) {
  load(filenames[[i]])
  tree <- paste0(treenames[[i]])
  pillai <- results$test$stat
  terms <- results$test$terms
  p_val <- results$test$pvalue
  ses <- (log(results$test$stat) - colMeans(log(results$test$nullstat))) / apply(log(results$test$nullstat), 2, sd)
  datmat<-rbind(pillai,ses,p_val)
  colnames(datmat)<-c("Size", "Diet", "Locomotion", "Size:Diet", "Size:Loc")
  data_tmp <- as_tibble(datmat)
  data_tmp <- data_tmp %>% mutate(tree=tree) %>% mutate(result = c("Pillai's Test Statistic","SES", "P value"))
  data <- bind_rows(data,data_tmp)
}

data <- data %>% mutate(tree = recode(data$tree, tree1_1 = "Tree Topology 1: Root Age 70-75",
                                      tree1_2 = "Tree Topology 1: Root Age 75-80",
                                      tree1_3 = "Tree Topology 1: Root Age 80-85",
                                      tree1_4 = "Tree Topology 1: Root Age 85-90",
                                      tree1_5 = "Tree Topology 1: Root Age 90-95",
                                      tree1_6 = "Tree Topology 1: Root Age 95-100"))


gttest<-data %>%
  gt(rowname_col = "result", groupname_col = "tree")%>%
  fmt_number(columns=everything(),decimals = 3)%>%
  tab_spanner(label = "Factors",columns=everything()) %>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_row_groups(everything()))%>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_body(rows = "P value" < 0.005 ))%>%
  tab_header(title = "Phylogenetic Regressions: Tree Topology 1, Full Dataset")%>%
  cols_align(align = "center")%>%
  cols_width(c("Size", "Diet", "Locomotion", "Size:Diet", "Size:Loc")~px(80))
gtsave(gttest, file="tree1_full_table.pdf")

########tree2#############################
#load data:
setwd("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/pmanovas/pmanovas_fullset_322/tree2/results/extant/")
filenames<-filenames<-dir(pattern=".Rdata")
treenames<-substr(filenames, start = 1, stop = 7)
data<-tibble() 
for (i in 1:length(filenames)) {
  load(filenames[[i]])
  tree <- paste0(treenames[[i]])
  pillai <- results$test$stat
  terms <- results$test$terms
  p_val <- results$test$pvalue
  ses <- (log(results$test$stat) - colMeans(log(results$test$nullstat))) / apply(log(results$test$nullstat), 2, sd)
  datmat<-rbind(pillai,ses,p_val)
  colnames(datmat)<-c("Size", "Diet", "Locomotion", "Habitat", "Development", "Activity", "Social", "Size:Diet", "Size:Loc", "Size:Hab","Size:Act", "Size:Dev", "Size:Soc")
  data_tmp <- as_tibble(datmat)
  data_tmp <- data_tmp %>% mutate(tree=tree) %>% mutate(result = c("Pillai's Test Statistic","SES", "P value"))
  data <- bind_rows(data,data_tmp)
}

data <- data %>% mutate(tree = recode(data$tree, tree2_1 = "Tree Topology 2: Root Age 70-75",
                                      tree2_2 = "Tree Topology 2: Root Age 75-80",
                                      tree2_3 = "Tree Topology 2: Root Age 80-85",
                                      tree2_4 = "Tree Topology 2: Root Age 85-90",
                                      tree2_5 = "Tree Topology 2: Root Age 90-95",
                                      tree2_6 = "Tree Topology 2: Root Age 95-100"))


gttest<-data %>%
  gt(rowname_col = "result", groupname_col = "tree")%>%
  fmt_number(columns=everything(),decimals = 3)%>%
  tab_spanner(label = "Factors",columns=everything()) %>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_row_groups(everything()))%>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_body(rows = "P value" < 0.005 ))%>%
  tab_header(title = "Phylogenetic Regressions: Tree Topology 2, Extant Taxa only")%>%
  cols_align(align = "center")%>%
  cols_width(c("Size", "Diet", "Locomotion", "Habitat", "Development", "Activity", "Social", "Size:Diet", "Size:Loc", "Size:Hab","Size:Act", "Size:Dev", "Size:Soc")~px(80))
gtsave(gttest, file="tree2_extant_table.pdf")

#load data:
setwd("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/pmanovas/pmanovas_fullset_322/tree2/results/full/")
filenames<-filenames<-dir(pattern=".Rdata")
treenames<-substr(filenames, start = 1, stop = 7)
data<-tibble() 
for (i in 1:length(filenames)) {
  load(filenames[[i]])
  tree <- paste0(treenames[[i]])
  pillai <- results$test$stat
  terms <- results$test$terms
  p_val <- results$test$pvalue
  ses <- (log(results$test$stat) - colMeans(log(results$test$nullstat))) / apply(log(results$test$nullstat), 2, sd)
  datmat<-rbind(pillai,ses,p_val)
  colnames(datmat)<-c("Size", "Diet", "Locomotion", "Size:Diet", "Size:Loc")
  data_tmp <- as_tibble(datmat)
  data_tmp <- data_tmp %>% mutate(tree=tree) %>% mutate(result = c("Pillai's Test Statistic","SES", "P value"))
  data <- bind_rows(data,data_tmp)
}

data <- data %>% mutate(tree = recode(data$tree, tree2_1 = "Tree Topology 2: Root Age 70-75",
                                      tree2_2 = "Tree Topology 2: Root Age 75-80",
                                      tree2_3 = "Tree Topology 2: Root Age 80-85",
                                      tree2_4 = "Tree Topology 2: Root Age 85-90",
                                      tree2_5 = "Tree Topology 2: Root Age 90-95",
                                      tree2_6 = "Tree Topology 2: Root Age 95-100"))


gttest<-data %>%
  gt(rowname_col = "result", groupname_col = "tree")%>%
  fmt_number(columns=everything(),decimals = 3)%>%
  tab_spanner(label = "Factors",columns=everything()) %>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_row_groups(everything()))%>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_body(rows = "P value" < 0.005 ))%>%
  tab_header(title = "Phylogenetic Regressions: Tree Topology 2, Full Dataset")%>%
  cols_align(align = "center")%>%
  cols_width(c("Size", "Diet", "Locomotion", "Size:Diet", "Size:Loc")~px(80))

gtsave(gttest, file="tree2_full_table.pdf")

##############################tree3###############################
#load data:
setwd("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/pmanovas/pmanovas_fullset_322/tree3/results/extant/")
filenames<-filenames<-dir(pattern=".Rdata")
treenames<-substr(filenames, start = 1, stop = 7)
data<-tibble() 
for (i in 1:length(filenames)) {
  load(filenames[[i]])
  tree <- paste0(treenames[[i]])
  pillai <- results$test$stat
  terms <- results$test$terms
  p_val <- results$test$pvalue
  ses <- (log(results$test$stat) - colMeans(log(results$test$nullstat))) / apply(log(results$test$nullstat), 2, sd)
  datmat<-rbind(pillai,ses,p_val)
  colnames(datmat)<-c("Size", "Diet", "Locomotion", "Habitat", "Development", "Activity", "Social", "Size:Diet", "Size:Loc", "Size:Hab","Size:Act", "Size:Dev", "Size:Soc")
  data_tmp <- as_tibble(datmat)
  data_tmp <- data_tmp %>% mutate(tree=tree) %>% mutate(result = c("Pillai's Test Statistic","SES", "P value"))
  data <- bind_rows(data,data_tmp)
}

data <- data %>% mutate(tree = recode(data$tree, tree3_1 = "Tree Topology 3: Root Age 70-75",
                                      tree3_2 = "Tree Topology 3: Root Age 75-80",
                                      tree3_3 = "Tree Topology 3: Root Age 80-85",
                                      tree3_4 = "Tree Topology 3: Root Age 85-90",
                                      tree3_5 = "Tree Topology 3: Root Age 90-95",
                                      tree3_6 = "Tree Topology 3: Root Age 95-100"))


gttest<-data %>%
  gt(rowname_col = "result", groupname_col = "tree")%>%
  fmt_number(columns=everything(),decimals = 3)%>%
  tab_spanner(label = "Factors",columns=everything()) %>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_row_groups(everything()))%>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_body(rows = "P value" < 0.005 ))%>%
  tab_header(title = "Phylogenetic Regressions: Tree Topology 3, Extant Taxa only")%>%
  cols_align(align = "center")%>%
  cols_width(c("Size", "Diet", "Locomotion", "Habitat", "Development", "Activity", "Social", "Size:Diet", "Size:Loc", "Size:Hab","Size:Act", "Size:Dev", "Size:Soc")~px(80))

gtsave(gttest, file="tree3_extant_table.pdf")

#load data:
setwd("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/pmanovas/pmanovas_fullset_322/tree3/results/full/")
filenames<-filenames<-dir(pattern=".Rdata")
treenames<-substr(filenames, start = 1, stop = 7)
data<-tibble() 
for (i in 1:length(filenames)) {
  load(filenames[[i]])
  tree <- paste0(treenames[[i]])
  pillai <- results$test$stat
  terms <- results$test$terms
  p_val <- results$test$pvalue
  ses <- (log(results$test$stat) - colMeans(log(results$test$nullstat))) / apply(log(results$test$nullstat), 2, sd)
  datmat<-rbind(pillai,ses,p_val)
  colnames(datmat)<-c("Size", "Diet", "Locomotion", "Size:Diet", "Size:Loc")
  data_tmp <- as_tibble(datmat)
  data_tmp <- data_tmp %>% mutate(tree=tree) %>% mutate(result = c("Pillai's Test Statistic","SES", "P value"))
  data <- bind_rows(data,data_tmp)
}

data <- data %>% mutate(tree = recode(data$tree, tree3_1 = "Tree Topology 3: Root Age 70-75",
                                      tree3_2 = "Tree Topology 3: Root Age 75-80",
                                      tree3_3 = "Tree Topology 3: Root Age 80-85",
                                      tree3_4 = "Tree Topology 3: Root Age 85-90",
                                      tree3_5 = "Tree Topology 3: Root Age 90-95",
                                      tree3_6 = "Tree Topology 3: Root Age 95-100"))


gttest<-data %>%
  gt(rowname_col = "result", groupname_col = "tree")%>%
  fmt_number(columns=everything(),decimals = 3)%>%
  tab_spanner(label = "Factors",columns=everything()) %>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_row_groups(everything()))%>%
  tab_style(style = list(cell_fill(color = 'grey'), 
                         cell_text(weight = 'bold')),
            locations = cells_body(rows = "P value" < 0.005 ))%>%
  tab_header(title = "Phylogenetic Regressions: Tree Topology 3, Full Dataset")%>%
  cols_align(align = "center")%>%
  cols_width(c("Size", "Diet", "Locomotion", "Size:Diet", "Size:Loc")~px(80))

gtsave(gttest, file="tree3_full_table.pdf")