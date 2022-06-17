library(tidyverse)

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_70_75/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()


resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_75_80/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_80_85/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_85_90/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_90_95/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree1/322_fullset_95_100/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

#####################tree2################################################
resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_70_75/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()


resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_75_80/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_80_85/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_85_90/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_90_95/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree2/322_fullset_95_100/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

####################################tree3#######################################
resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree3/322_fullset_70_75/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()


resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree3/322_fullset_75_80/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree3/322_fullset_80_85/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree3/322_fullset_85_90/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree3/322_fullset_90_95/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()

resultsfolder<-"C:/Users/Anjali Goswami/Dropbox/Work/Projects/ERC/Synmammals/all_landmarked_placentals/analyses/rates/rates_fullset_322/tree3/322_fullset_95_100/extant/results/act"
resultslist<-dir(resultsfolder, pattern = "*.Rdata$", full.names = TRUE)
extractor_func<-function(x){
  load(x)
  results_fit<-results$fit
  return(results_fit)
}
raw_rate_results <- lapply(resultslist, function(x) extractor_func(x))
rate_results<- raw_rate_results %>% unlist() %>% enframe()

rate_results <- rate_results %>% rename(act=name, Rate=value)
ggplot(rate_results, aes(y=Rate, fill=act)) +
  geom_boxplot()
  