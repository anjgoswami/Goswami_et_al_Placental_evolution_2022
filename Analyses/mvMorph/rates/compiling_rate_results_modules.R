library(tidyverse)
library(fs)
library(stringr)
devtools::install_github("kassambara/ggpubr")
library(ggpubr)

########Tree 1 extant########################
###########Activity####################

data_dir<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/act/"
resultsfiles <- dir(data_dir)
resultsfiles <- paste("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/act/",resultsfiles, sep = "")
my_results <- tibble(filename=resultsfiles)
moduleslist <- c("basiocc","basisph","frontal","glenoid","jugal","max_d","max_v","nasal","occ_cond","palatine","parietal","premax_d","premax_v","pterygoid","squam_v","squam_z","supraocc")
my_results<- my_results %>%
  mutate(., Tree = str_extract(filename, "(?<=tree).+(?=_bin)")) %>%
  mutate(., Module = str_extract(filename, paste0(moduleslist,collapse = "|")))
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultsfiles, function(x) extractor_func(x))
rate_results<- do.call(rbind, raw_rate_results) %>% as_tibble()
my_results <- bind_cols(my_results,rate_results)
my_results <- my_results %>% pivot_longer(
  cols = `All`:`Nocturnal`,
  names_to = "Activity",
  values_to = "Rate"
)


ggplot(my_results, aes(y=Rate, x = Module, fill=Activity))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x= element_text(angle=45, hjust=1.3, vjust=1.2))

ggplot(subset(my_results, Activity != "All"), aes(y=Rate, x = Module, fill=Activity))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x= element_text(angle=45, hjust=1.3, vjust=1.2))


a<-ggplot(my_results, aes(y=Rate, x = Module , fill=Activity))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


a2<-ggplot(subset(my_results, Activity != "All"), aes(y=Rate, x = Module , fill=Activity))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_act.pdf",
       plot = a, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_act_wo_all.pdf",
       plot = a2, device = cairo_pdf,
       width = 11, height = 11, units = "cm")


###########Social####################

data_dir<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/soc/"
resultsfiles <- dir(data_dir)
resultsfiles <- paste("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/soc/",resultsfiles, sep = "")
my_results <- tibble(filename=resultsfiles)
moduleslist <- c("nasal","premax_d","premax_v","max_d","max_v","palatine","pterygoid","jugal","squam_z","glenoid","squam_v","frontal","parietal" ,"supraocc","basisph","basiocc","occ_cond")
my_results<- my_results %>%
  mutate(., Tree = str_extract(filename, "(?<=tree).+(?=_bin)")) %>%
  mutate(., Module = str_extract(filename, paste0(moduleslist,collapse = "|")))
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultsfiles, function(x) extractor_func(x))
rate_results<- do.call(rbind, raw_rate_results) %>% as_tibble()
my_results <- bind_cols(my_results,rate_results)
my_results <- my_results %>% pivot_longer(
  cols = `Social`:`Solitary`,
  names_to = "Sociality",
  values_to = "Rate"
)

my_results$Module <- factor(my_results$Module, levels = moduleslist)

s<-ggplot(my_results, aes(y=Rate, x = Module , fill=Sociality))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_soc.pdf",
       plot = s, device = cairo_pdf,
       width = 11, height = 11, units = "cm")


###########Development####################

data_dir<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/dev/"
resultsfiles <- dir(data_dir)
resultsfiles <- paste("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/dev/",resultsfiles, sep = "")
my_results <- tibble(filename=resultsfiles)
moduleslist <- c("nasal","premax_d","premax_v","max_d","max_v","palatine","pterygoid","jugal","squam_z","glenoid","squam_v","frontal","parietal" ,"supraocc","basisph","basiocc","occ_cond")
my_results<- my_results %>%
  mutate(., Tree = str_extract(filename, "(?<=tree).+(?=_bin)")) %>%
  mutate(., Module = str_extract(filename, paste0(moduleslist,collapse = "|")))
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultsfiles, function(x) extractor_func(x))
rate_results<- do.call(rbind, raw_rate_results) %>% as_tibble()
my_results <- bind_cols(my_results,rate_results)
my_results <- my_results %>% pivot_longer(
  cols = `Altricial`:`Precocial`,
  names_to = "Development",
  values_to = "Rate"
)

my_results$Module <- factor(my_results$Module, levels = moduleslist)

d<-ggplot(my_results, aes(y=Rate, x = Module , fill=Development))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_dev.pdf",
       plot = d, device = cairo_pdf,
       width = 11, height = 11, units = "cm")



###########Habitat####################

data_dir<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/hab/"
resultsfiles <- dir(data_dir)
resultsfiles <- paste("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/hab/",resultsfiles, sep = "")
my_results <- tibble(filename=resultsfiles)
moduleslist <- c("nasal","premax_d","premax_v","max_d","max_v","palatine","pterygoid","jugal","squam_z","glenoid","squam_v","frontal","parietal" ,"supraocc","basisph","basiocc","occ_cond")
my_results<- my_results %>%
  mutate(., Tree = str_extract(filename, "(?<=tree).+(?=_bin)")) %>%
  mutate(., Module = str_extract(filename, paste0(moduleslist,collapse = "|")))
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultsfiles, function(x) extractor_func(x))
rate_results<- do.call(rbind, raw_rate_results) %>% as_tibble()
my_results <- bind_cols(my_results,rate_results)
my_results <- my_results %>% pivot_longer(
  cols = `Aquatic`:`Open`,
  names_to = "Habitat",
  values_to = "Rate"
)

my_results$Module <- factor(my_results$Module, levels = moduleslist)

h<-ggplot(my_results, aes(y=Rate, x = Module , fill=Habitat))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_hab.pdf",
       plot = h, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

h2<-ggplot(subset(my_results, Habitat != "Aquatic"), aes(y=Rate, x = Module , fill=Habitat))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_hab_wo_aquatic.pdf",
       plot = h2, device = cairo_pdf,
       width = 11, height = 11, units = "cm")



###########Locomotion####################

data_dir<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/loc/"
resultsfiles <- dir(data_dir)
resultsfiles <- paste("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/loc/",resultsfiles, sep = "")
my_results <- tibble(filename=resultsfiles)
moduleslist <- c("nasal","premax_d","premax_v","max_d","max_v","palatine","pterygoid","jugal","squam_z","glenoid","squam_v","frontal","parietal" ,"supraocc","basisph","basiocc","occ_cond")
my_results<- my_results %>%
  mutate(., Tree = str_extract(filename, "(?<=tree).+(?=_bin)")) %>%
  mutate(., Module = str_extract(filename, paste0(moduleslist,collapse = "|")))
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultsfiles, function(x) extractor_func(x))
rate_results<- do.call(rbind, raw_rate_results) %>% as_tibble()
my_results <- bind_cols(my_results,rate_results)
my_results <- my_results %>% pivot_longer(
  cols = `Aquatic`:`Volant`,
  names_to = "Locomotion",
  values_to = "Rate"
)

my_results$Module <- factor(my_results$Module, levels = moduleslist)

l<-ggplot(my_results, aes(y=Rate, x = Module , fill=Locomotion))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_loc.pdf",
       plot = l, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

l2<-ggplot(subset(my_results, Locomotion != "Aquatic"), aes(y=Rate, x = Module , fill=Locomotion))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_loc_wo_aquatic.pdf",
       plot = l2, device = cairo_pdf,
       width = 11, height = 11, units = "cm")


###########Diet####################

data_dir<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/diet/"
resultsfiles <- dir(data_dir)
resultsfiles <- paste("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/extant/results/diet/",resultsfiles, sep = "")
my_results <- tibble(filename=resultsfiles)
moduleslist <- c("nasal","premax_d","premax_v","max_d","max_v","palatine","pterygoid","jugal","squam_z","glenoid","squam_v","frontal","parietal" ,"supraocc","basisph","basiocc","occ_cond")
my_results<- my_results %>%
  mutate(., Tree = str_extract(filename, "(?<=tree).+(?=_bin)")) %>%
  mutate(., Module = str_extract(filename, paste0(moduleslist,collapse = "|")))
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultsfiles, function(x) extractor_func(x))
rate_results<- do.call(rbind, raw_rate_results) %>% as_tibble()
my_results <- bind_cols(my_results,rate_results)
my_results <- my_results %>% pivot_longer(
  cols = `Bulk invertivore`:`Social insectivore`,
  names_to = "Diet",
  values_to = "Rate"
)

my_results$Module <- factor(my_results$Module, levels = moduleslist)

f<-ggplot(my_results, aes(y=Rate, x = Module , fill=Diet))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_diet.pdf",
       plot = f, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

f2<-ggplot(subset(my_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, x = Module , fill=Diet))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/extant_mods_tree1_diet_wo_whales.pdf",
       plot = f2, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

#################Full dataset################################################
###########Locomotion####################

data_dir<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/full/results/loc/"
resultsfiles <- dir(data_dir)
resultsfiles <- paste("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/full/results/loc/",resultsfiles, sep = "")
my_results <- tibble(filename=resultsfiles)
moduleslist <- c("nasal","premax_d","premax_v","max_d","max_v","palatine","pterygoid","jugal","squam_z","glenoid","squam_v","frontal","parietal" ,"supraocc","basisph","basiocc","occ_cond")
my_results<- my_results %>%
  mutate(., Tree = str_extract(filename, "(?<=tree).+(?=_bin)")) %>%
  mutate(., Module = str_extract(filename, paste0(moduleslist,collapse = "|")))
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultsfiles, function(x) extractor_func(x))
rate_results<- do.call(rbind, raw_rate_results) %>% as_tibble()
my_results <- bind_cols(my_results,rate_results)
my_results <- my_results %>% pivot_longer(
  cols = `Aquatic`:`Volant`,
  names_to = "Locomotion",
  values_to = "Rate"
)

my_results$Module <- factor(my_results$Module, levels = moduleslist)

l<-ggplot(my_results, aes(y=Rate, x = Module , fill=Locomotion))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/full_mods_tree1_loc.pdf",
       plot = l, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

l2<-ggplot(subset(my_results, Locomotion != "Aquatic"), aes(y=Rate, x = Module , fill=Locomotion))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/full_mods_tree1_loc_wo_aquatic.pdf",
       plot = l2, device = cairo_pdf,
       width = 11, height = 11, units = "cm")


###########Diet####################

data_dir<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/full/results/diet/"
resultsfiles <- dir(data_dir)
resultsfiles <- paste("C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/modules/full/results/diet/",resultsfiles, sep = "")
my_results <- tibble(filename=resultsfiles)
moduleslist <- c("nasal","premax_d","premax_v","max_d","max_v","palatine","pterygoid","jugal","squam_z","glenoid","squam_v","frontal","parietal" ,"supraocc","basisph","basiocc","occ_cond")
my_results<- my_results %>%
  mutate(., Tree = str_extract(filename, "(?<=tree).+(?=_bin)")) %>%
  mutate(., Module = str_extract(filename, paste0(moduleslist,collapse = "|")))
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultsfiles, function(x) extractor_func(x))
rate_results<- do.call(rbind, raw_rate_results) %>% as_tibble()
my_results <- bind_cols(my_results,rate_results)
my_results <- my_results %>% pivot_longer(
  cols = `Bulk invertivore`:`Social insectivore`,
  names_to = "Diet",
  values_to = "Rate"
)

my_results$Module <- factor(my_results$Module, levels = moduleslist)

f<-ggplot(my_results, aes(y=Rate, x = Module , fill=Diet))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/full_mods_tree1_diet.pdf",
       plot = f, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

f2<-ggplot(subset(my_results, Diet %in% c("Carnivore", "Herbivore","Omnivore", "Invertivore", "Social insectivore")), aes(y=Rate, x = Module , fill=Diet))+
  geom_boxplot(lwd=.2, outlier.size = .5)+
  theme(axis.text.x= element_text(angle=45, hjust=1, vjust=1, size =6), panel.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))	+
  scale_x_discrete(labels=c("Nasal","Premaxilla (d)","Premaxilla (v)","Maxilla (d)","Maxilla (v)","Palatine","Pterygoid","Jugal","Squamosal (z)","Glenoid","Squamosal (v)","Frontal","Parietal","Supraoccipital","Basisphenoid","Basioccipital","Occipital Condyle"))+
  geom_vline(xintercept=seq(1.5, length(unique(my_results$Module))-0.5, 1), lwd=.2, colour="grey")+
  theme(panel.border = element_rect(size=.5, fill=NA))


ggsave(filename = "./Analyses/Figures/Rates/full_mods_tree1_diet_wo_whales.pdf",
       plot = f2, device = cairo_pdf,
       width = 11, height = 11, units = "cm")
