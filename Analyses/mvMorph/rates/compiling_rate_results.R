figfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/Figures/violins/"

resultsfolder<-"./Analyses/mvMorph/rates/results/fine_diet_322/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
ggplot(rate_results, aes(x = Diet, y=log(Rate), fill=Diet)) +
  geom_violin()+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6))

#install.packages(ggstatsplot)
library(ggstatsplot)
rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))
rate_results_2$Diet <- recode(rate_results$Diet, `Bulk invertivore` = "Bulk/ninvertivore",`Social insectivore` = "Social/ninsectivore")

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = FALSE, type = "np",
  x = Diet,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 45,hjust=0.5,vjust=0.2))

#########Tree 2#################################
##########Bin 3 diet############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/full/results/diet/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
# ggplot(rate_results, aes(x = Diet, y=log(Rate), fill=Diet)) +
#   geom_violin()
# 
# ggplot(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, fill=Diet)) +
#   geom_boxplot()+
#   theme_bw()+
#   scale_fill_brewer(palette = "Dark2")

rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))
rate_results_2$Diet <- recode(rate_results$Diet, `Bulk invertivore` = "Bulk\ninvertivore",`Social insectivore` = "Social\ninsectivore")

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = FALSE, type = "np",
  x = Diet,   y = LogRate, xlab = "", ylab= "Log(Rate)", centrality.plotting = FALSE, results.subtitle=FALSE)

violins+theme(axis.text.x = element_text(angle = 45,hjust=0.5,vjust=0.2))+
  scale_color_manual(values=c("royalblue2","red2","springgreen3","tan2","gold2","deeppink","tan4"))
ggsave(filename = paste0(figfolder, "rates_diet_full2.pdf"))

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = TRUE, type = "np",
  x = Diet,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_diet_full_stats.pdf"))

ggplot(data = rate_results_2,aes(x = Diet,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Diet))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+
  scale_color_manual(values=c("royalblue2","red2","springgreen3","tan2","gold2","deeppink","tan4"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none")

ggsave(filename = paste0(figfolder, "rates_diet_full_3.pdf"))

kruskal.test(rate_results)
kruskal.test(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")))

##########Bin 3 locomotion############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/full/results/loc/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Locomotion=name, Rate=value)
#ggplot(rate_results, aes(x = Locomotion, y=log(Rate), fill=Locomotion)) +
  geom_violin()

rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = FALSE, type = "np",
  x = Locomotion,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_loc_full.pdf"))

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = TRUE, type = "np",
  x = Locomotion,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_loc_full_stats.pdf"))

ggplot(data = rate_results_2,aes(x = Locomotion,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Locomotion))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+
  scale_color_manual(values=c("royalblue2","springgreen4","tan4","grey","turquoise","springgreen2","tan2", "brown","skyblue"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1), legend.position = "none")

ggsave(filename = paste0(figfolder, "rates_loc_full_2.pdf"))

 
kruskal.test(rate_results)
#kruskal.test(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")))

##########Bin 3 extant diet############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/extant/results/diet/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
# ggplot(rate_results, aes(x = Diet, y=log(Rate), fill=Diet)) +
#   geom_violin()
# 
# ggplot(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, fill=Diet)) +
#   geom_boxplot()+
#   theme_bw()+
#   scale_fill_brewer(palette = "Dark2")

rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = FALSE, type = "np",
  x = Diet,   y = LogRate,     centrality.plotting = FALSE)
violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_diet_extant.pdf"))

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = TRUE, type = "np",
  x = Diet,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_diet_extant_stats.pdf"))



kruskal.test(rate_results)
kruskal.test(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")))


##########Bin 3 extant loc############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/extant/results/loc/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Locomotion=name, Rate=value)
# ggplot(rate_results, aes(x = Locomotion, y=log(Rate), fill=Locomotion)) +
#   geom_violin()
# 
# ggplot(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, fill=Diet)) +
#   geom_boxplot()+
#   theme_bw()+
#   scale_fill_brewer(palette = "Dark2")

rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = FALSE, type = "np",
  x = Locomotion,   y = LogRate,     centrality.plotting = FALSE)
violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_loc_extant.pdf"))

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = TRUE, type = "np",
  x = Locomotion,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_loc_extant_stats.pdf"))

kruskal.test(rate_results)
kruskal.test(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")))


##########Bin 3 extant social############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/extant/results/soc/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Social=name, Rate=value)
rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

violins<-ggbetweenstats(
  data = subset(rate_results_2, Social !="unknown"),
  pairwise.comparisons = FALSE, type = "np",
  x = Social,   y = LogRate,     centrality.plotting = FALSE)
violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_soc_extant.pdf"))

violins<-ggbetweenstats(
  data = subset(rate_results_2, Social !="unknown"),
  pairwise.comparisons = TRUE, type = "np",
  x = Social,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_social_extant_stats.pdf"))

ggplot(data = subset(rate_results_2, Social !="unknown"),aes(x = Social,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Social))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+
  #ylim(-5,2.5)+
  scale_color_manual(values=c("springgreen2","red2"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none")

ggsave(filename = paste0(figfolder, "rates_soc_all.pdf"))


# ggplot(subset(rate_results, Social !="unknown"), aes(x = Social, y=Rate, fill=Social)) +
#   geom_violin()
# 
# ggplot(subset(rate_results, Social !="unknown"), aes(x = Social, y=log(Rate), fill=Social)) +
#   geom_violin()
# 
# ggplot(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, fill=Diet)) +
#   geom_boxplot()+
#   theme_bw()+
#   scale_fill_brewer(palette = "Dark2")
# 
# 
# kruskal.test(rate_results)
# kruskal.test(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")))


##########Bin 3 extant development############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/extant/results/dev/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Development=name, Rate=value)
rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

violins<-ggbetweenstats(
  data = subset(rate_results_2, Development !="Unknown"),
  pairwise.comparisons = FALSE, type = "np",
  x = Development,   y = LogRate,     centrality.plotting = FALSE)
violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_dev_extant.pdf"))

violins<-ggbetweenstats(
  data = subset(rate_results_2, Development !="Unknown"),
  pairwise.comparisons = TRUE, type = "np",
  x = Social,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_dev_extant_stats.pdf"))

ggplot(data = subset(rate_results_2, Development !="Unknown"),aes(x = Development,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Development))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+
  scale_color_manual(values=c("salmon2","turquoise2"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none")

ggsave(filename = paste0(figfolder, "rates_dev_2.pdf"))

# ggplot(subset(rate_results, Development != "Unknown"), aes(x = Development, y=log(Rate), fill=Development)) +
#   geom_violin()
# 
# ggplot(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, fill=Diet)) +
#   geom_boxplot()+
#   theme_bw()+
#   scale_fill_brewer(palette = "Dark2")
# 
# 
# kruskal.test(rate_results)
# kruskal.test(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")))


##########Bin 3 Activity extant############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/extant/results/act/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Activity=name, Rate=value)
rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

violins<-ggbetweenstats(
  data = subset(rate_results_2, Activity !="unknown"),
  pairwise.comparisons = FALSE, type = "np",
  x = Activity,   y = LogRate,     centrality.plotting = FALSE)
violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_act_extant.pdf"))

violins<-ggbetweenstats(
  data = subset(rate_results_2, Activity !="unknown"),
  pairwise.comparisons = TRUE, type = "np",
  x = Activity,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_act_extant_stats.pdf"))

ggplot(data = rate_results_2,aes(x = Activity,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Activity))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+
  scale_color_manual(values=c("springgreen2","royalblue2","orange2", "blue4"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none")

ggsave(filename = paste0(figfolder, "rates_act_2.pdf"))

# ggplot(rate_results, aes(x = Activity, y=log(Rate), fill=Activity)) +
#   geom_violin()
# 
# ggplot(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, fill=Diet)) +
#   geom_boxplot()+
#   theme_bw()+
#   scale_fill_brewer(palette = "Dark2")
# 
# 
# kruskal.test(rate_results)
# kruskal.test(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")))




##########Bin 3 Habitat extant############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/extant/results/hab/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Habitat=name, Rate=value)

rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))
rate_results_2$Habitat <- recode(rate_results$Habitat, `Mixed_terrestrial` = "Mixed terrestrial")

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = FALSE, type = "np",
  x = Habitat,   y = LogRate,     centrality.plotting = FALSE)
violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_hab_extant.pdf"))

violins<-ggbetweenstats(
  data = rate_results_2,
  pairwise.comparisons = TRUE, type = "np",
  x = Habitat,   y = LogRate,     centrality.plotting = FALSE)

violins+theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust=0.2))
ggsave(filename = paste0(figfolder, "rates_hab_extant_stats.pdf"))

ggplot(data = rate_results_2,aes(x = Habitat,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Habitat))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+
  scale_color_manual(values=c("turquoise4", "darkgreen", "tan1", "springgreen2","gold2"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none")

ggsave(filename = paste0(figfolder, "rates_hab_2.pdf"))


# ggplot(rate_results, aes(x = Habitat, y=log(Rate), fill=Habitat)) +
#   geom_violin()
# 
# ggplot(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, fill=Diet)) +
#   geom_boxplot()+
#   theme_bw()+
#   scale_fill_brewer(palette = "Dark2")
# 
# 
# kruskal.test(rate_results)
# kruskal.test(subset(rate_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")))

