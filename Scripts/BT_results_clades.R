
#requires input files with clade associations, here from full_species_data, and output from BayesTraits analyses, here myrateslist.R, code for generating this from BayesTraits output included below
species_data<-read.csv("./Raw_Data/full_species_data.csv")
load(file="myrateslist.R")


tree_labels<-c("70_75","75_80","80_85","85_90","90_95","95_100")
names(my_rates_list)<-tree_labels
for (ii in 1:length(my_rates_list)){
results2<-my_rates_list[[ii]]$results
phyloseq <- c("Zalambdalestidae", "Cimolesta", "Leptictida", "Cingulata", "Pilosa", "Afrosoricida", "Macroscelidea","Tubulidentata", "Hyracoidea", "Proboscidea", "Desmostylia", "Sirenia", "Embrithopoda", "Dermoptera", "Scandentia", "Primates", "Rodentia", "Lagomorpha", "Chiroptera", "Eulipotyphla", "Acreodi", "Amblypoda", "Artiodactyla", "Cetacea", "Litopterna", "Notoungulata", "Astrapotheria","Perissodactyla", "Creodonta", "Carnivora", "Pholidota")
phyloseq <- phyloseq[which((phyloseq %in% unique(results2$Group)))]
results3 <- results2
results3$Group <- factor(results3$Group, levels = phyloseq)
p2 <- ggplot(data=results3, aes(x=log(MeanRate), group=Group, fill=Group)) +
  geom_density(adjust=1.5, alpha=.9)+
  scale_fill_viridis_d(breaks=phyloseq)

species_data2 <- species_data %>% filter(Order %in% unique(results2$Group))
species_data2 <- species_data2 %>%
  filter(!duplicated(paste0(pmax(Superorder, Order), pmin(Superorder, Order))))
species_data2 <- species_data2 %>% select(c(Superorder, Order))

 x<-list("Carnivora"=  "Laurasiatheria" ,
 "Pilosa"=  "Xenarthra"      ,
 "Primates"=  "Euarchontaglires",
 "Notoungulata"  = "Laurasiatheria"  ,
"Artiodactyla"= "Laurasiatheria"  ,
  "Cetacea" = "Laurasiatheria"  ,
  "Afrosoricida"= "Afrotheria"    ,
  "Perissodactyla"= "Laurasiatheria" ,
"Rodentia"= "Euarchontaglires",
  "Chiroptera"=  "Laurasiatheria" ,
  "Amblypoda"=  "Laurasiatheria"  ,
  "Cingulata"     =  "Xenarthra"       ,
"Sirenia"=  "Afrotheria"    ,
  "Eulipotyphla"=  "Laurasiatheria" ,
  "Proboscidea"=  "Afrotheria"      ,
  "Creodonta"     =  "Laurasiatheria"  ,
"Leptictida"=  "Stem" ,
  "Lagomorpha"=  "Euarchontaglires",
  "Litopterna"=  "Laurasiatheria"  ,
  "Macroscelidea" =  "Afrotheria" ,
"Scandentia" =  "Euarchontaglires")

 results3<-results3 %>% mutate(Superorder = recode(results3$Group, !!!x))
 results4 <- split(results3,f = results3$Superorder)
 library(gridExtra)
 library(RColorBrewer)
 library(patchwork)
 p1 <- ggplot(data=results4$Stem, aes(x=log(MeanRate), group=Group, fill=Group)) +
   geom_density(adjust=1.5, alpha=.9)+
   scale_fill_brewer(palette = "Set3", breaks=phyloseq)+
   facet_grid(Superorder~.)+
   theme_bw()+
   guides(fill=guide_legend(ncol=2))+
   xlim(-4.5,1)+
   ylim(0,6.7)+
   theme(strip.text.y = element_text(size = 7))

 p2 <- p1 %+% results4$Xenarthra 
 p3 <- p2 %+% results4$Afrotheria
 p4 <- p3 %+% results4$Euarchontaglires
 p5 <- p4 %+% results4$Laurasiatheria
 
 p1 <- p1 + theme( axis.text.x = element_blank(),  axis.title.x = element_blank(),  axis.ticks.x = element_blank())
 p2 <- p2 + theme( axis.text.x = element_blank(),  axis.title.x = element_blank(),  axis.ticks.x = element_blank())
 p3 <- p3 + theme( axis.text.x = element_blank(),  axis.title.x = element_blank(),  axis.ticks.x = element_blank())
 p4 <- p4 + theme( axis.text.x = element_blank(),  axis.title.x = element_blank(),  axis.ticks.x = element_blank())
 
 p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 1) + plot_annotation(title = paste0("Tree Topology 1: Root Age ", names(my_rates_list)[ii]) )
 
 ggsave(filename=paste0("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/BayesTraits/Results/whole_skull_322/", "Tree_1b_", names(my_rates_list)[ii], ".pdf"), device="pdf", width=7, height=7, units="in")
}

############code for generating myrateslist from BayesTraits output#################

####load my rate Summary function
rate_summary <- function(rjpp_out,lookup.table,  group_column, taxa, min_clade_size=2){
  lookup.table <- rename(lookup.table, GROUP = group_column)
  lookup.table <- rename(lookup.table, SPECIES = taxa)
  current.species <- filter(lookup.table, SPECIES %in% rjpp_out$tree_summary$tree_summaries$original_tree$tip.label)
  clades <- unique(lookup.table$GROUP)
  #select the sub tree for each group
  time_tree<-rjpp_out$tree_summary$tree_summaries$original_tree
  time_trees <- lapply(1:length(clades), function(x) keep.tip(time_tree, (filter(current.species, GROUP == clades[x]) %>% pull(SPECIES))))
  names(time_trees)<-clades
  taxa_per_clade <- unlist(lapply(1:length(time_trees), function(x) Ntip(time_trees[[x]])))
  time_trees <- time_trees[which(taxa_per_clade>=min_clade_size)]
  #get the time represented by each sub tree
  time_per_tree <- unlist(lapply(1:length(time_trees), function(x) max(phytools::nodeHeights(time_trees[[x]]))))
  #the first column is the rate on the root and our time tree has no root edge so this is removed with the "[-1,]"

  rate_tibble <- rjpp_out$rates
  mean_rate_tibble<-rjpp_out$scalars %>% filter(nTips==1)
  species_key<-as_tibble(matrix(unlist(rjpp_out$species_key[mean_rate_tibble$node_id]), ncol=3,byrow=TRUE)) %>%
    rename(ancNode = V1, descNode=V2, SPECIES=V3)
  mean_rate_tibble <- full_join(mean_rate_tibble,species_key)
  mean_rate_tibble <- mean_rate_tibble %>% full_join(lookup.table, by = "SPECIES")
  mean_rate_tibble <- mean_rate_tibble %>% select(node_id, branch, ancNode, descNode, meanRate, SPECIES, GROUP)
  rate_tibble <- left_join(mean_rate_tibble, rate_tibble)

  rate_tibble_long <- rate_tibble %>% pivot_longer(cols = -c(node_id,branch,ancNode,descNode,meanRate,SPECIES,GROUP), names_to = "iter", values_to = "rate")
  rate_tibble_split<-split(rate_tibble_long,f = rate_tibble_long$iter)

  f1<-function(x){
    x %>%
      group_by(GROUP)%>%
      summarize(mean_rate = mean(rate, na.rm = TRUE))
  }
  clade_summaries <- lapply(rate_tibble_split, f1)
  clade_summaries_joined <- map_dfr(clade_summaries, bind_rows)
  return(clade_summaries_joined )
  }



library(tidyverse)
library(ape)
library(BTprocessR)
library(pbapply)
library(phytools)
treelist <- dir()

varratesA <- list.files(pattern = "*a_run_.txt.VarRates.txt$", recursive = TRUE)
#in this example, varrates with lambda transformation is the best supported model, so extracting those results for run a (from multiple runs for convergence check)
varratesA <- varratesA[grepl("var_rate", varratesA)]
varratesA <- varratesA[grepl("lambda", varratesA)]

outtreesA <- list.files(pattern = "*a_run_.txt.Output.trees$", recursive = TRUE)
outtreesA <- outtreesA[grepl("var_rate", outtreesA)]
outtreesA <- outtreesA[grepl("lambda", outtreesA)]

treefolder <- ("insert folder with trees")
timetrees <- paste0(treefolder,dir(treefolder))

filetable <- tibble(varrate = varratesA, outtrees = outtreesA, timetrees = timetrees)

summaries <- lapply(1:nrow(filetable), function(x) rjpp(rjlog = filetable$varrate[x], rjtrees = filetable$outtrees[x], tree = read.nexus(filetable$timetrees[x])))

save(summaries, file="summaries.R")

species_data<-read_csv("full_species_data.csv")

my_rates_list <- lapply(1:length(summaries), function(x) rate_summary(rjpp_out = summaries[[x]],
                                                                  lookup.table = species_data,
                                                                  group_column = "Order",
                                                                  taxa = "Tip_Label",
                                                                  min_clade_size = 2))

save(my_rates_list, file="myrateslist.R")
