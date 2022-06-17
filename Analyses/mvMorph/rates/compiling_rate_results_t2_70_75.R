figfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/Figures/violins/T2_70_75/"

#########Tree 2#################################
##########Bin 3 diet############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/full/results/diet/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))
rate_results_2$Diet <- recode(rate_results$Diet, `Bulk invertivore` = "Bulk\ninvertivore",`Social insectivore` = "Social\ninsectivore")

pdiet<-ggplot(data = rate_results_2,aes(x = Diet,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Diet))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+ggtitle("Diet - Fossil + Extant")+
  scale_color_manual(values=c("royalblue2","red2","springgreen3","tan2","gold2","deeppink","tan4"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none", plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

##########Bin 3 locomotion############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/full/results/loc/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Locomotion=name, Rate=value)
rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

ploc<-ggplot(data = rate_results_2,aes(x = Locomotion,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Locomotion))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+ggtitle("Locomotion - Fossil + Extant")+
  scale_color_manual(values=c("royalblue2","springgreen4","tan4","grey","turquoise","springgreen2","tan2", "brown","skyblue"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1), legend.position = "none", plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

##########Bin 3 extant diet############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/extant/results/diet/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Diet=name, Rate=value)
rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

pdietex<-ggplot(data = rate_results_2,aes(x = Diet,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Diet))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+ggtitle("Diet - Extant")+
  scale_color_manual(values=c("royalblue2","red2","springgreen3","tan2","gold2","deeppink","tan4"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none", plot.title = element_text(size = 10, face = "bold", hjust = 0.5))


##########Bin 3 extant loc############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/extant/results/loc/"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(Locomotion=name, Rate=value)

rate_results_2<-rate_results%>%mutate(LogRate=log(Rate))

plocex<-ggplot(data = rate_results_2,aes(x = Locomotion,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Locomotion))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+ggtitle("Locomotion - Extant")+
  scale_color_manual(values=c("royalblue2","springgreen4","tan4","grey","turquoise","springgreen2","tan2", "brown","skyblue"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1), legend.position = "none", plot.title = element_text(size = 10, face = "bold", hjust = 0.5))



##########Bin 3 extant social############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/extant/results/soc/"
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

psoc<-ggplot(data = subset(rate_results_2, Social !="unknown"),aes(x = Social,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Social))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+ggtitle("Sociality")+
  #ylim(-5,2.5)+
  scale_color_manual(values=c("springgreen2","red2"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none", plot.title = element_text(size = 10, face = "bold", hjust = 0.5))



##########Bin 3 extant development############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/extant/results/dev/"
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

pdev<-ggplot(data = subset(rate_results_2, Development !="Unknown"),aes(x = Development,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Development))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+ggtitle("Development")+
  scale_color_manual(values=c("salmon2","turquoise2"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none", plot.title = element_text(size = 10, face = "bold", hjust = 0.5))


##########Bin 3 Activity extant############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/extant/results/act/"
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

pact<-ggplot(data = rate_results_2,aes(x = Activity,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Activity))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+ggtitle("Activity")+
  scale_color_manual(values=c("springgreen2","royalblue2","orange2", "blue4"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none", plot.title = element_text(size = 10, face = "bold", hjust = 0.5))



##########Bin 3 Habitat extant############################

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/extant/results/hab/"
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

phab<-ggplot(data = rate_results_2,aes(x = Habitat,   y = LogRate))+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=.4, aes(color=Habitat))+
  geom_boxplot(fill=NA,width=.4)+
  geom_violin(fill=NA,width=.7)+
  xlab("") + ylab("Log(Rate)")+ggtitle("Habitat")+
  scale_color_manual(values=c("turquoise4", "darkgreen", "tan1", "springgreen2","gold2"))+
  ggstatsplot::theme_ggstatsplot()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),legend.position = "none", plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

p.combined <- grid.arrange(pdiet, ploc, pdietex, plocex, phab, pdev, psoc, pact, nrow =2, top = text_grob("Tree Topology 2: Root Age = 75-80 Ma", size = 16, face = "bold", just = "centre"))

ggsave(filename = paste0(figfolder, "t2_75_violins.pdf"),
       plot = p.combined, device = cairo_pdf,
       width = 30, height = 22, units = "cm")

